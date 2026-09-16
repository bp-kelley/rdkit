//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Search/RotorRefine.h"
#include "Search/ThompsonSamplingSearch.h"

#include <RDGeneral/RDLog.h>

#include <cstdio>

#include <set>

#include "Joiner/JoinerProfiling.h"
#include "Search/InterFragScore.h"
#include "Search/RotorDriver.h"
#include "Utils/SymmetricRmsd.h"
#include "Sampler/TorsionSampler.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace RDKit {

using detail::nowNs;
using detail::prof;
using detail::profiling;

namespace {
inline double circDiff(double a, double b) {
  double d = std::fmod(a - b + 540.0, 360.0) - 180.0;
  return std::fabs(d);
}
}  // namespace

// XXX FIX ME: I'm not fully versed on thompson sampling in general.
//  I have a few candidates to look over this.
//  It doesn't help as much as I would have hoped with larger nRotors
//   where it really should shine.

std::string ThompsonSamplingSearch::validateParams(
    const RigidRotorSearchParams &sp, const std::string &ffVariant) const {
  const std::string err = RigidRotorSearch::validateParams(sp, ffVariant);
  if (!err.empty()) {
    return err;
  }
  if (sp.thompson.minBudget > sp.thompson.maxBudget) {
    return "search.thompson.minBudget must be <= search.thompson.maxBudget";
  }
  return {};
}

std::vector<SearchResult> ThompsonSamplingSearch::search(
    const FragmentJoinerContext &ctx, const RigidRotorSearchParams &sp) {
  const auto &params = sp;
  std::vector<SearchResult> out;
  unsigned int nNaN = 0;  //!< failed evaluations; see the NaN handling below

  const std::vector<unsigned int> extraRotors = ctx.getIntraRotorBonds(sp);
  std::vector<unsigned int> driveBonds = ctx.rotorBonds;
  driveBonds.insert(driveBonds.end(), extraRotors.begin(), extraRotors.end());
  
  RotorDriver drv(ctx.mol, driveBonds, -1, ctx.scorer);
  const size_t nr = drv.numRotors();
  
  // hash rotors to central bond
  std::vector<char> rotorIsIntra(nr, 0), rotorIsAtrop(nr, 0);
  if (!extraRotors.empty()) {
    const std::set<unsigned int> intraSet(extraRotors.begin(),
                                          extraRotors.end());
    for (size_t r = 0; r < nr; ++r) {
      const auto q = drv.torsion(static_cast<unsigned int>(r));
      const Bond *cb = ctx.mol.getBondBetweenAtoms(q[1], q[2]);
      if (cb && intraSet.count(cb->getIdx())) {
        rotorIsIntra[r] = 1;
        if (ctx.getIntraType(cb->getIdx()) == IntraRotorType::Atropisomer) {
          rotorIsAtrop[r] = 1;
        }
      }
    }
  }
  const size_t nf = ctx.frags.size();

  // Assume rotor angles are informed priors
  // angles start at Beta alpha = priorStrength; a uniform-grid is simply alpha
  // ROTOR WEIGHT:  Add a size weight so a kept conformer credits rotors that move more atoms
  //  THIS is s novel approach, but may be too clever for it's own good
  const double priorS = std::max(1.0, params.thompson.priorStrength);
  const double backstop = params.thompson.backstopStepDeg;
  const double sizeExp = std::max(0.0, params.thompson.sizePriorExp);
  std::vector<std::vector<double>> rotorArms(nr), rotorPrior(nr);
  std::vector<double> rotorW(nr, 1.0);
  double maxMove = 1.0;
  for (size_t r = 0; r < nr; ++r) {
    maxMove = std::max(maxMove, static_cast<double>(drv.movingAtoms(r).size()));
  }
  for (size_t r = 0; r < nr; ++r) {
    // ROTOR WEIGHT: a rotor's affect on the pose scales with atoms moved.
    //  so scale the prior based on # atoms, small # atoms have smaller effects
    //  on diversity
    const double effect =
        static_cast<double>(drv.movingAtoms(r).size()) / maxMove;
    const double prefPrior = 1.0 + (priorS - 1.0) * std::pow(effect, sizeExp);
    std::vector<double> pref;
    if (params.torsionSampler) {
      try {
	// remember we have intra and inter rotors: intra
	//  can optionally tweek ring puckers/etc
        auto t = drv.torsion(static_cast<unsigned int>(r));
        const bool basin = rotorIsIntra[r]
                               ? true
                               : useBasinAnglesForRotor(params.junctionBasinAngles,
							ctx.mol, t[1], t[2]);
        pref = params.torsionSampler->getAngles(ctx.mol, t[0], t[1], t[2], t[3],
                                                basin);
      } catch (...) {
      }
    }
    for (double ang : pref) {
      rotorArms[r].push_back(ang);
      rotorPrior[r].push_back(prefPrior);  // add strong prior
    }
    if (backstop > 0.0) {
      for (double g = -180.0; g < 180.0; g += backstop) {
        bool near = false;
        for (double p : pref) {
          if (circDiff(g, p) < backstop * 0.5) {
            near = true;
            break;
          }
        }
        if (!near) {
          rotorArms[r].push_back(g);
          rotorPrior[r].push_back(1.0);  // add weak prior
        }
      }
    }
    if (rotorArms[r].empty()) {
      rotorArms[r] = params.defaultAngles;
      rotorPrior[r].assign(rotorArms[r].size(), 1.0);
    }
    if (rotorIsAtrop[r]) {
      // The far well is the OTHER enantiomer: keep the candidates on the side
      // the declared axis is already on.
      const auto kept = basinLimitAngles(
          rotorArms[r], drv.dihedralDeg(static_cast<unsigned int>(r)));
      std::vector<double> pri(kept.size(), 1.0);
      for (size_t k = 0; k < kept.size(); ++k) {
        for (size_t j = 0; j < rotorArms[r].size(); ++j) {
          if (rotorArms[r][j] == kept[k]) { pri[k] = rotorPrior[r][j]; break; }
        }
      }
      rotorArms[r] = kept;
      rotorPrior[r] = pri;
    }
    rotorW[r] = static_cast<double>(drv.movingAtoms(r).size()) / maxMove;
  }

  // rotorW affects two unrelated params.
  //  we need to be able to twiddle both independently for A/B testing
  double minMove = std::numeric_limits<double>::max();
  for (size_t r = 0; r < nr; ++r) {
    minMove = std::min(minMove, static_cast<double>(drv.movingAtoms(r).size()));
  }
  minMove = std::max(1.0, minMove);
  auto weightFor = [&](RotorWeighting w, size_t r, double exp) {
    switch (w) {
      case RotorWeighting::MovingAtoms:
        return exp == 1.0 ? rotorW[r] : std::pow(rotorW[r], exp);
      case RotorWeighting::Inverted:
        // mirror of MovingAtoms: the SMALLEST rotor gets 1.0 and larger ones
        // less, so this is the same shape reflected rather than a new scale.
	//  XXX FIX ME -> this was just to see if the opposite
	//                guess from what I thought would be better would
	//                actually be worse.  It is.
        return minMove /
               std::max(1.0, static_cast<double>(drv.movingAtoms(r).size()));
      case RotorWeighting::Uniform:
        break;
    }
    return 1.0;
  };
  std::vector<double> novW(nr, 1.0), postW(nr, 1.0);
  for (size_t r = 0; r < nr; ++r) {
    novW[r] = weightFor(params.thompson.noveltyWeighting, r,
                        std::max(0.0, params.thompson.noveltyWeightExp));
    postW[r] = weightFor(params.thompson.posteriorWeighting, r, 1.0);
  }
  
  // sample fragment confs
  std::vector<unsigned int> fragK(nf);
  for (size_t f = 0; f < nf; ++f) {
    fragK[f] = std::min<unsigned int>(
        std::max(1u, params.fragConfBranch),
        static_cast<unsigned int>(ctx.frags[f].confs.size()));
    if (fragK[f] == 0) fragK[f] = 1;
  }

  std::vector<std::vector<double>> fA(nf), fB(nf), rA(nr), rB(nr);
  for (size_t f = 0; f < nf; ++f) {
    fA[f].assign(fragK[f], 1.0);
    fB[f].assign(fragK[f], 1.0);
  }
  for (size_t r = 0; r < nr; ++r) {
    rA[r] = rotorPrior[r];  // rotor priors start higher
    rB[r].assign(rotorArms[r].size(), 1.0);
  }

  std::mt19937 rng(params.randomSeed);
  auto betaSample = [&](double a, double b) {
    std::gamma_distribution<double> ga(a, 1.0), gb(b, 1.0);
    double x = ga(rng), y = gb(rng);
    return (x + y > 0.0) ? x / (x + y) : 0.5;
  };

  // Reward:  true if novel and < ewindow
  const ThompsonParams &tp = params.thompson;
  const double window = params.energyWindow;
  const double divThr = params.diversityRmsThresh;

  // use novel rotor state OR RMSD for pruning
  const double angThr = params.thompson.noveltyAngleDeg;
  const bool useFp = angThr > 0.0 && divThr > 0.0;
  const double angThrSq = angThr * angThr;
  const auto order = drv.numRotorAtoms();

  std::vector<SearchResult> kept;
  std::vector<std::vector<double>> keptAng;  // weighted-torsion-fingerprint angles
  double best = std::numeric_limits<double>::infinity();
  std::vector<unsigned int> confChoice(nf, 0);
  std::vector<int> rotArm(nr, 0);

  unsigned int budget = params.thompsonBudget;
  if (budget == 0) {
    size_t fragArms = 0;
    for (size_t f = 0; f < nf; ++f) fragArms += fragK[f];
    long long b =
        static_cast<long long>(params.thompson.perRotor) * nr +
        static_cast<long long>(params.thompson.perFragConf) * fragArms;
    b = std::max<long long>(params.thompson.minBudget,
                            std::min<long long>(params.thompson.maxBudget, b));
    budget = static_cast<unsigned int>(b);
  }
  d_lastBudget = budget;

  std::map<std::vector<unsigned int>, std::vector<double>> conformerCache;

  // Use the context Rotor driving hierarchical
  const unsigned int nAtoms = ctx.mol.getNumAtoms();
  std::vector<char> varies(nAtoms, 0);
  for (size_t r = 0; r < nr; ++r) {
    for (unsigned int a : drv.movingAtoms(static_cast<unsigned int>(r))) {
      varies[a] = 1;
    }
  }
  std::vector<char> dirtyFrag(nf, 0);
  for (size_t f = 0; f < nf; ++f) {
    if (fragK[f] > 1) dirtyFrag[f] = 1;
  }

  // ctx.edges is in BFS so we know fixed and moving atom
  for (const auto &ed : ctx.edges) {
    if (dirtyFrag[ed.parentFrag]) dirtyFrag[ed.childFrag] = 1;
  }
  for (size_t f = 0; f < nf; ++f) {
    if (dirtyFrag[f]) {
      for (unsigned int a : ctx.frags[f].atoms) varies[a] = 1;
    }
  }
  std::vector<unsigned int> cmpIdx;
  cmpIdx.reserve(static_cast<size_t>(nAtoms) * 3);
  for (unsigned int a = 0; a < nAtoms; ++a) {
    if (varies[a]) {
      cmpIdx.push_back(3 * a);
      cmpIdx.push_back(3 * a + 1);
      cmpIdx.push_back(3 * a + 2);
    }
  }
  const double divLimit =
      divThr * divThr * static_cast<double>(nAtoms ? nAtoms : 1);

  const auto deadline = timeOut(params);
  for (unsigned int s = 0; s < budget; ++s) {
    // Stop drawing on the budget; `kept` already holds every conformer accepted
    // so far.
    //  XXX FIX ME -> why check every 64 draws?
    //   This is an odd magic number
    if ((s & 0x3F) == 0 && timedOut(deadline)) break;
    // Thompson-pick a conformer for each fragment
    for (size_t f = 0; f < nf; ++f) {
      int ba = 0;
      double bt = -1.0;
      for (unsigned int k = 0; k < fragK[f]; ++k) {
        double t = betaSample(fA[f][k], fB[f][k]);
        if (t > bt) {
          bt = t;
          ba = static_cast<int>(k);
        }
      }
      confChoice[f] = ba;
    }
    long long tp0 = profiling() ? nowNs() : 0;
    auto pit = conformerCache.find(confChoice);
    if (pit == conformerCache.end()) {
      pit = conformerCache.emplace(confChoice, ctx.placeAll(confChoice)).first;
    }
    drv.positions() = pit->second;  // need to copy copy placement; driving mutates it
    long long tp1 = profiling() ? nowNs() : 0;
    
    // Thompson-pick + rotor angle, coarse to fine
    for (unsigned int r : order) {
      const auto &ar = rotorArms[r];
      if (ar.empty()) continue;
      int ba = 0;
      double bt = -1.0;
      for (size_t k = 0; k < ar.size(); ++k) {
        double t = betaSample(rA[r][k], rB[r][k]);
        if (t > bt) {
          bt = t;
          ba = static_cast<int>(k);
        }
      }
      rotArm[r] = ba;
      drv.setDihedral(r, ar[ba]);
    }
    
    long long tp2 = profiling() ? nowNs() : 0;
    double sc = drv.score();
    long long tp3 = profiling() ? nowNs() : 0;
    if (profiling()) {
      prof().tPlace += tp1 - tp0;
      prof().tDrive += tp2 - tp1;
      prof().tScore += tp3 - tp2;
      prof().nSamples += 1;
      prof().nDrives += order.size();
      prof().nScores += 1;
    }
    if (!std::isnan(sc) && sc < best) best = sc;
    
    // Drop any score which is noncomputable
    if (std::isnan(sc)) {
      ++nNaN;
    }
    const bool inWindow = !std::isnan(sc) &&
                          (!std::isfinite(best) || sc <= best + window);
    bool novel = true;
    std::vector<double> curAng;
    if (useFp) {  // use novelty fp for pruning
      curAng.resize(nr);
      for (size_t r = 0; r < nr; ++r)
        curAng[r] = rotorArms[r].empty() ? 0.0 : rotorArms[r][rotArm[r]];
      for (const auto &ka : keptAng) {
        double num = 0.0, den = 0.0;
        for (size_t r = 0; r < nr; ++r) {
          const double d = circDiff(curAng[r], ka[r]);
          num += novW[r] * d * d;
          den += novW[r];
        }
        if (den > 0.0 && num / den < angThrSq) {
          novel = false;
          break;
        }
      }
    } else if (divThr > 0.0) {
      // Use as-is RMSD for a quick check, non trans/rot optimized
      const std::vector<double> &pos = drv.positions();
      for (const auto &k : kept) {
        double sq = 0.0;
        bool within = true;
        for (unsigned int ci : cmpIdx) {
          const double d = pos[ci] - k.coords[ci];
          sq += d * d;
          if (sq >= divLimit) {
            within = false;
            break;
          }
        }
        if (within) {
          novel = false;
          break;
        }
      }
    }
    if (profiling()) prof().tDiv += nowNs() - tp3;
    const bool isKeeper = inWindow && novel;
    if (isKeeper) {
      kept.push_back({drv.positions(), sc});
      if (useFp) keptAng.push_back(std::move(curAng));
    }

    const double reward = isKeeper ? 1.0 : 0.0;

    // reward the chosen arms this draw chose
    for (size_t f = 0; f < nf; ++f) {
      fA[f][confChoice[f]] += reward;
      fB[f][confChoice[f]] += 1.0 - reward;
    }
    
    for (size_t r = 0; r < nr; ++r) {
      if (rotorArms[r].empty()) continue;
      rA[r][rotArm[r]] += reward * postW[r];
      rB[r][rotArm[r]] += (1.0 - reward) * postW[r];
    }
  }

  // TS_ARMSTATS: (NOTE: CLAUDE REVIEW)
  // is the posterior actually LEARNING per rotor, and does that
  // depend on how many atoms the rotor moves?  rotorW scales every update by
  // movingAtoms/maxMoving, so a terminal rotor can finish a run still close to
  // its prior no matter how often it was pulled.  Bucketing by moving-atom
  // quartile separates "never pulled" from "pulled but never learned".
  //
  // evidence = total Beta mass added beyond the prior (a+b - prior), i.e. how
  // much this arm was actually updated.  entropy is over the arm-choice
  // distribution mean(a/(a+b)) -- high entropy at the end of a run means the
  // posterior never committed.
  if (params.diagnostics.TS_ARMSTATS && nr) {
    std::vector<size_t> byMove(nr);
    for (size_t r = 0; r < nr; ++r) {
      byMove[r] = r;
    }
    std::sort(byMove.begin(), byMove.end(), [&](size_t x, size_t y) {
      return drv.movingAtoms(x).size() < drv.movingAtoms(y).size();
    });
    std::printf("[ts_armstats] rotors=%zu budget=%u\n", nr, budget);
    std::printf(
        "[ts_armstats] %-6s %6s %8s %8s %9s %9s\n", "quart", "rotors",
        "movAtoms", "rotorW", "evidence", "entropy");
    for (int q = 0; q < 4; ++q) {
      const size_t lo = nr * q / 4, hi = nr * (q + 1) / 4;
      if (lo >= hi) {
        continue;
      }
      double mov = 0.0, w = 0.0, ev = 0.0, ent = 0.0;
      size_t cnt = 0;
      for (size_t i = lo; i < hi; ++i) {
        const size_t r = byMove[i];
        if (rotorArms[r].empty()) {
          continue;
        }
        mov += static_cast<double>(drv.movingAtoms(r).size());
        w += postW[r];
        double tot = 0.0;
        std::vector<double> p(rotorArms[r].size(), 0.0);
        for (size_t k = 0; k < rotorArms[r].size(); ++k) {
          // subtract the prior so this reports EVIDENCE, not initial mass
          ev += (rA[r][k] + rB[r][k]) - (rotorPrior[r][k] + 1.0);
          p[k] = rA[r][k] / std::max(1e-9, rA[r][k] + rB[r][k]);
          tot += p[k];
        }
        double e = 0.0;
        for (double v : p) {
          const double pp = v / std::max(1e-9, tot);
          if (pp > 0.0) {
            e -= pp * std::log(pp);
          }
        }
        // normalise so 1.0 == uniform over this rotor's arms, i.e. "learned
        // nothing"; comparable across rotors with different arm counts
        ent += rotorArms[r].size() > 1
                   ? e / std::log(static_cast<double>(rotorArms[r].size()))
                   : 0.0;
        ++cnt;
      }
      if (!cnt) {
        continue;
      }
      std::printf("[ts_armstats] Q%-5d %6zu %8.1f %8.3f %9.1f %9.4f\n", q + 1,
                  cnt, mov / cnt, w / cnt, ev / cnt, ent / cnt);
    }
  }

  for (auto &k : kept) {
    if (!std::isnan(k.score) &&
        (!std::isfinite(best) || k.score <= best + window)) {
      out.push_back(std::move(k));
    }
  }
  if (nNaN) {
    BOOST_LOG(rdWarningLog)
        << "ThompsonSamplingSearch: " << nNaN << " of " << budget
        << " draws scored NaN and were discarded" << std::endl;
  }
  
  std::sort(out.begin(), out.end(),
            [](const SearchResult &a, const SearchResult &b) {
              if (std::isnan(a.score)) return false;
              if (std::isnan(b.score)) return true;
              return a.score < b.score;
            });

  // Local refinement search around the MMFF basin (optional)
  if (params.thompson.refineSteps > 0) {
    // refine in inter-frag scores
    //  XXX FIX ME why make ScoreFn here?
    RotorDriver::ScoreFn refineFn;
    refineRotors(drv, out, params.thompson.refineSteps,
                        params.thompson.refineStepDeg,
                        params.thompson.refinePasses, refineFn);
  }


  if (params.finalSymmetryDedup) {
    ctx.symmetryDedup(out, divThr);
  }

  // Keep maxConfs by the selected out mode
  {
    const size_t maxConfs = tp.maxConfs;
    const bool diverse = tp.outMode == OutputSelection::Diverse;
    const bool stratify = tp.outMode == OutputSelection::Stratify;
    // Note: diversityRmsThresh == 0 (DISABLED) makes Diverse selection
    // degenerate to energy order ony
    const double outLimit =
        divThr * divThr * static_cast<double>(nAtoms ? nAtoms : 1);
    if (maxConfs > 0 && out.size() > maxConfs) {
      if (diverse) {
        std::vector<SearchResult> sel;  // energy-ordered greedy RMSD-diverse pick
        sel.reserve(maxConfs);
        for (auto &r : out) {
          bool distinct = true;
          for (const auto &s : sel) {
            double sq = 0.0;
            for (unsigned int ci : cmpIdx) {
              const double d = r.coords[ci] - s.coords[ci];
              sq += d * d;
              if (sq >= outLimit) break;
            }
            if (sq < outLimit) {
              distinct = false;
              break;
            }
          }
          if (distinct) {
            sel.push_back(std::move(r));
            if (sel.size() >= maxConfs) break;
          }
        }
        out.swap(sel);
      } else if (stratify) {
        // energy-stratified stride: N confs evenly spaced across the
        // energy-sorted pool, attempt to bind high-strain buried
	// ligands
	//  XXX FIX ME -> apparently doesn't work.
        std::vector<SearchResult> sel;
        sel.reserve(maxConfs);
        const double step =
            static_cast<double>(out.size()) / static_cast<double>(maxConfs);
        for (size_t i = 0; i < maxConfs; ++i) {
          size_t idx = static_cast<size_t>(i * step);
          if (idx >= out.size()) idx = out.size() - 1;
          sel.push_back(std::move(out[idx]));
        }
        out.swap(sel);
      } else {
        out.resize(maxConfs);  // energy-lowest
      }
    }
  }
  if (profiling()) {
    prof().nKept += static_cast<long long>(kept.size());
    prof().nOut += static_cast<long long>(out.size());
  }
  return out;
}

}  // namespace RDKit
