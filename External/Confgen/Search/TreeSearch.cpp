//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Search/TreeSearch.h"

#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/Nonbonded.h>
#include <ForceField/MMFF/TorsionAngle.h>
#include <ForceField/MMFF/Params.h>
#include "Utils/NonbondedLookup.h"
#include "Search/TreeSearch.h"
#include "Joiner/JoinerProfiling.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <set>

namespace RDKit {

using detail::nowNs;
using detail::prof;
using detail::profiling;

RotorTree::RotorTree(RotorDriver &driver, RigidRotorSearchParams params)
    : d_driver(driver), d_params(std::move(params)) {
  d_angles.resize(driver.numRotors());
}

void RotorTree::setAngles(unsigned int rotor, std::vector<double> angles) {
  if (rotor < d_angles.size()) {
    d_angles[rotor] = std::move(angles);
  }
}

void RotorTree::enableIncremental(ForceFields::MMFF::InterFragVdWContrib *vdw,
                                  ForceFields::MMFF::TorsionAngleContrib *tor) {
  d_incremental = IncrementalInterFragScore(d_driver, InterFragScoreHandles{vdw, tor});
}

const std::vector<double> &RotorTree::anglesForRotor(unsigned int rotor) const {
  if (rotor < d_angles.size() && !d_angles[rotor].empty()) {
    return d_angles[rotor];
  }
  return d_params.defaultAngles;
}

namespace {
// smallest unsigned circular difference between two angles (degrees)
inline double circDiff(double a, double b) {
  double d = std::fmod(a - b + 540.0, 360.0) - 180.0;
  return std::fabs(d);
}
// True iff the coordinate RMSD between two same-length flat buffers is < thr
// (shared frame -> no alignment needed; all partials keep the fixed root-side
// atoms in place).  RMSD = sqrt(sum_sqdist / n), so RMSD < thr <=> sum_sqdist <
// thr^2 * n; we accumulate the squared distance and bail the instant it reaches
// that limit.  Most candidate pairs differ far more than thr, so the scan
// aborts after a handful of coordinates instead of touching every atom --
// bit-for-bit the same boolean as the old coordRms(a,b) < thr, just without the
// wasted tail.
inline bool withinRms(const std::vector<double> &a,
                      const std::vector<double> &b, double thr) {
  const size_t n = a.size() / 3;
  const double limit = thr * thr * static_cast<double>(n ? n : 1);
  double s = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    const double d = a[i] - b[i];
    s += d * d;
    if (s >= limit) return false;
  }
  return true;
}
}  // namespace

std::vector<SearchResult> RotorTree::search() {
  if (d_params.thompsonBudget > 0) {
    return searchThompson();
  }
  struct State {
    std::vector<double> coords;
    double score;
    std::vector<double> setAngles;  //!< chosen dihedral per rotor set so far
  };

  std::vector<State> beam;
  beam.push_back({d_driver.positions(), d_driver.score(), {}});

  const double tol = d_params.tree.angleTolerance;
  const double window = d_params.energyWindow;
  const double divThr = d_params.diversityRmsThresh;

  auto selectBeam = [&](std::vector<State> &&cand) -> std::vector<State> {
    std::sort(cand.begin(), cand.end(), [](const State &a, const State &b) {
      // NaN scores (no scorer) sort as equal-ish; push them to the back
      if (std::isnan(a.score)) return false;
      if (std::isnan(b.score)) return true;
      return a.score < b.score;
    });
    double best = cand.empty() ? 0.0 : cand.front().score;
    bool haveWindow = !std::isnan(best);
    std::vector<State> kept;
    for (auto &st : cand) {
      if (haveWindow && !std::isnan(st.score) && st.score > best + window) {
        break;  // sorted, so nothing further is in-window
      }
      // DIVERSITY-preserving retention: within the energy window, drop a
      // candidate only if it is geometrically too close to one already kept.
      // This keeps higher-energy-but-distinct conformers (the bioactive pose)
      // that pure energy-greedy top-K would prune.  When off (divThr == 0),
      // fall back to the cheap angle-fingerprint dedup and let the energy order
      // + beamWidth decide.
      bool drop = false;
      if (divThr > 0.0) {
        for (const auto &k : kept) {
          if (withinRms(st.coords, k.coords, divThr)) {
            drop = true;
            break;
          }
        }
      } else if (tol > 0.0) {
        for (const auto &k : kept) {
          bool same = true;
          for (size_t r = 0; r < st.setAngles.size(); ++r) {
            if (circDiff(st.setAngles[r], k.setAngles[r]) > tol) {
              same = false;
              break;
            }
          }
          if (same) {
            drop = true;
            break;
          }
        }
      }
      if (drop) {
        continue;
      }
      kept.push_back(std::move(st));
      if (kept.size() >= d_params.tree.beamWidth) {
        break;
      }
    }
    return kept;
  };

  // INCREMENTAL junction-local rescoring.  Driving rotor r moves only the atoms
  // in movingAtoms(r), so of the cross-fragment vdW pairs only those with
  // EXACTLY ONE endpoint in that moving set change distance (a pair with BOTH
  // endpoints moved rotates rigidly -> its distance, hence energy, is
  // unchanged; XOR, not OR).  We precompute those "changing" pair indices per
  // rotor once, then rescore a single- rotor move as  newTotal = oldScore -
  // oldChangingVdW + newChangingVdW  plus the full junction-torsion delta
  // (O(junction bonds), cheap).  Bit-for-bit equal to the full rescore;
  // validated under ASM_SCOREVALIDATE.  Needs the inter-fragment vdW contrib
  // handle (RotorTree::enableIncremental); otherwise fall back to full drive().
  // The packing is built ONCE in enableIncremental() and reused across every
  // search(), so a RotorTree hoisted out of the per-seed loop never re-packs.
  // isValid() is false for the CHNOPS lookup variant -> full drive() below.
  const bool usePacked = d_incremental.isValid();
  const bool validate = usePacked && d_diag.ASM_SCOREVALIDATE;
  size_t mismatches = 0;
  if (usePacked && d_diag.ASM_INCRSTATS) {
    size_t nPairs = d_incremental.numPairs();
    size_t sumCp = 0, maxCp = 0;
    for (unsigned int r : d_driver.numRotorAtoms()) {
      sumCp += d_incremental.numChangingPairs(r);
      maxCp = std::max(maxCp, d_incremental.numChangingPairs(r));
    }
    size_t nr = d_driver.numRotorAtoms().size();
    BOOST_LOG(rdWarningLog)
        << "[INCRSTATS] nRotors=" << nr << " nPairs=" << nPairs
        << " avgChanging=" << (nr ? (double)sumCp / nr : 0.0)
        << " maxChanging=" << maxCp << " avgRatio="
        << ((nr && nPairs) ? (double)sumCp / (nr * nPairs) : 0.0) << "\n";
  }

  for (unsigned int r : d_driver.numRotorAtoms()) {
    const auto &angles = anglesForRotor(r);
    std::vector<State> next;
    next.reserve(beam.size() * angles.size());
    for (const auto &st : beam) {
      // The "old" contributions at st.coords are INVARIANT across this rotor's
      // candidate angles, so hoist them out of the angle loop (compute once per
      // parent state, not once per candidate).  We also restore st.coords ONCE
      // here: setDihedral sets the ABSOLUTE angle of only movingAtoms(r),
      // leaving every other atom == st.coords, so successive candidate angles
      // need no re-restore (rotor r is not nested in itself and this level's
      // rotor is fixed).
      double oldVdw = 0.0, oldTor = 0.0;
      if (usePacked) {
        d_driver.positions() = st.coords;
        const double *buf = d_driver.positions().data();
        oldVdw = d_incremental.changingVdw(buf, r);
        oldTor = d_incremental.torsionEnergy(buf);
      }
      for (double a : angles) {
        double sc;
        if (usePacked) {
          d_driver.setDihedral(r, a);  // rotate moving set to absolute angle a
          const double *nbuf = d_driver.positions().data();
          double newVdw = d_incremental.changingVdw(nbuf, r);
          double newTor = d_incremental.torsionEnergy(nbuf);
          sc = st.score - oldVdw + newVdw - oldTor + newTor;
          if (validate) {
            double full =
                d_driver.score();  // full inter-fragment sum at new coords
            double allowed = 1.0e-4 * (1.0 + std::fabs(full));
            if (std::fabs(sc - full) > allowed) {
              ++mismatches;
              // The incremental rescore is meant to be bit-equal to the full
              // one; a mismatch is a bug in the packing, hence rdErrorLog.
              BOOST_LOG(rdErrorLog)
                  << "[ASM_SCOREVALIDATE] MISMATCH rotor=" << r
                  << " angle=" << a << " incr=" << sc << " full=" << full
                  << " diff=" << (sc - full)
                  << " (dVdw=" << (newVdw - oldVdw)
                  << " dTor=" << (newTor - oldTor) << ")\n";
            }
          }
        } else {
          d_driver.positions() = st.coords;  // restore this partial
          sc = d_driver.drive(r, a);  // set rotor r (coords only), rescore
        }
        State ns;
        ns.coords = d_driver.positions();
        ns.score = sc;
        ns.setAngles = st.setAngles;
        ns.setAngles.push_back(a);
        next.push_back(std::move(ns));
      }
    }
    if (!next.empty()) {
      beam = selectBeam(std::move(next));
    }
  }
  if (validate) {
    auto &log = mismatches ? rdErrorLog : rdWarningLog;
    BOOST_LOG(log) << "[ASM_SCOREVALIDATE] " << mismatches << " mismatch(es)\n";
  }

  std::vector<SearchResult> out;
  out.reserve(beam.size());
  for (auto &st : beam) {
    out.push_back({std::move(st.coords), st.score});
  }
  // leave the driver holding the best result
  if (!out.empty()) {
    d_driver.positions() = out.front().coords;
  }
  return out;
}

std::vector<SearchResult> RotorTree::searchThompson() {
  const auto &order = d_driver.numRotorAtoms();
  const size_t nr = d_driver.numRotors();
  const std::vector<double> base = d_driver.positions();  // seed geometry

  // bandit arms = each rotor's candidate angles, with Beta(a,b) reward params
  std::vector<std::vector<double>> arms(nr);
  std::vector<std::vector<double>> a(nr), b(nr);
  for (size_t r = 0; r < nr; ++r) {
    arms[r] = anglesForRotor(static_cast<unsigned int>(r));
    a[r].assign(arms[r].size(), 1.0);
    b[r].assign(arms[r].size(), 1.0);
  }

  std::mt19937 rng(d_params.randomSeed);
  auto betaSample = [&](double alpha, double beta) {
    std::gamma_distribution<double> ga(alpha, 1.0), gb(beta, 1.0);
    double x = ga(rng), y = gb(rng);
    return (x + y > 0.0) ? x / (x + y) : 0.5;
  };

  const double window = d_params.energyWindow;
  const double divThr = d_params.diversityRmsThresh;

  std::vector<SearchResult> kept;
  double best = std::numeric_limits<double>::infinity();
  std::vector<int> chosen(nr, 0);

  for (unsigned int s = 0; s < d_params.thompsonBudget; ++s) {
    d_driver.positions() = base;
    // Thompson-pick + drive each rotor (coarse to fine, so a parent is set
    // before its children -- the driven dihedrals stay independent)
    for (unsigned int r : order) {
      const auto &ar = arms[r];
      if (ar.empty()) {
        continue;
      }
      int bestArm = 0;
      double bestT = -1.0;
      for (size_t k = 0; k < ar.size(); ++k) {
        double t = betaSample(a[r][k], b[r][k]);
        if (t > bestT) {
          bestT = t;
          bestArm = static_cast<int>(k);
        }
      }
      chosen[r] = bestArm;
      d_driver.drive(r, ar[bestArm]);
    }
    double sc = d_driver.score();
    if (!std::isnan(sc) && sc < best) {
      best = sc;
    }
    // keep if in the (running) energy window and geometrically novel
    bool inWindow =
        std::isnan(sc) || !std::isfinite(best) || sc <= best + window;
    bool novel = true;
    if (divThr > 0.0) {
      for (const auto &k : kept) {
        if (withinRms(d_driver.positions(), k.coords, divThr)) {
          novel = false;
          break;
        }
      }
    }
    const bool keptIt = inWindow && novel;
    if (keptIt) {
      kept.push_back({d_driver.positions(), sc});
    }
    // reward the arms this sample chose
    for (size_t r = 0; r < nr; ++r) {
      if (arms[r].empty()) {
        continue;
      }
      if (keptIt) {
        a[r][chosen[r]] += 1.0;
      } else {
        b[r][chosen[r]] += 1.0;
      }
    }
  }

  // prune to the energy window of the final best, sort best-first
  std::vector<SearchResult> out;
  for (auto &k : kept) {
    if (std::isnan(k.score) || !std::isfinite(best) ||
        k.score <= best + window) {
      out.push_back(std::move(k));
    }
  }
  std::sort(out.begin(), out.end(),
            [](const SearchResult &x, const SearchResult &y) {
              if (std::isnan(x.score)) return false;
              if (std::isnan(y.score)) return true;
              return x.score < y.score;
            });
  if (!out.empty()) {
    d_driver.positions() = out.front().coords;
  }
  return out;
}

// ---------------------------------------------------------------------------
// TreeSearch: the deterministic beam RigidRotorSearch (declared in
// TreeSearch.h). This is the old FragmentJoiner::assemble() deterministic path,
// now a strategy driven by the shared FragmentJoinerContext.  It lives here,
// next to the RotorTree it uses.
// ---------------------------------------------------------------------------
std::string TreeSearch::validateParams(const RigidRotorSearchParams &sp,
                                       const std::string &ffVariant) const {
  const std::string err = RigidRotorSearch::validateParams(sp, ffVariant);
  if (!err.empty()) {
    return err;
  }
  if (sp.tree.beamWidth == 0) {
    return "search.tree.beamWidth must be > 0";
  }
  return {};
}

std::vector<SearchResult> TreeSearch::search(const FragmentJoinerContext &ctx,
                                             const RigidRotorSearchParams &sp) {
  const auto &params = sp;
  std::vector<SearchResult> all;

  // Deterministic ROTOR Tree search: a few fragment-conformer seeds, each
  // beam-searched.
  const unsigned int branch = std::max(1u, params.fragConfBranch);
  std::vector<std::vector<unsigned int>> seeds;
  std::set<std::vector<unsigned int>> seen;
  for (unsigned int s = 0; s < params.rootSeeds; ++s) {
    std::vector<unsigned int> choice(ctx.frags.size(), 0);
    for (size_t f = 0; f < ctx.frags.size(); ++f) {
      unsigned int nUse = std::min<unsigned int>(
          branch, static_cast<unsigned int>(ctx.frags[f].confs.size()));
      if (s != 0 && nUse > 1) {
        unsigned int h =
            (s * 2654435761u) ^ (static_cast<unsigned int>(f) * 40503u);
        choice[f] = h % nUse;
      }
    }
    if (seen.insert(choice).second) {
      seeds.push_back(std::move(choice));
    }
  }

  // Build the driver + tree once; enable incremental junction-local scoring if
  // possible.
  long long tSetup0 = profiling() ? nowNs() : 0;
  RotorDriver drv(ctx.mol, ctx.rotorBonds, -1, ctx.scorer);
  RotorTree tree(drv, params);
  tree.setDiagnostics(sp.diagnostics);
  // disableIncrementalScoring is an A/B against the full per-move rescore; the
  // incremental path is bit-identical.  No-op when handles are null
  // (useFullFFScorer).
  if (ctx.scoreHandles.vdw && !sp.disableIncrementalScoring) {
    tree.enableIncremental(ctx.scoreHandles.vdw, ctx.scoreHandles.tor);
  }

  // Prep the per-rotor candidate angles from the torsion sampler.
  if (params.torsionSampler) {
    for (size_t e = 0; e < ctx.edges.size(); ++e) {
      auto t = drv.torsion(static_cast<unsigned int>(e));
      try {
        auto angles = params.torsionSampler->getAngles(
            ctx.mol, t[0], t[1], t[2], t[3],
            useBasinAnglesForRotor(params.junctionBasinAngles, ctx.mol, t[1], t[2]));
        if (!angles.empty()) {
          tree.setAngles(static_cast<unsigned int>(e), angles);
        }
      } catch (...) {
      }
    }
  }
  if (profiling())
    prof().tDrive += nowNs() - tSetup0;  // one-time driver+angles setup

  for (const auto &choice : seeds) {
    long long ts0 = profiling() ? nowNs() : 0;
    std::vector<double> buffer = ctx.placeAll(choice);
    drv.positions() = std::move(buffer);
    long long ts1 = profiling() ? nowNs() : 0;

    // search the tree
    auto r = tree.search();
    if (profiling()) {
      long long ts2 = nowNs();
      prof().tPlace += ts1 - ts0;  // per-seed placeAll
      prof().tScore += ts2 - ts1;  // per-seed beam search
      prof().nSamples += 1;
    }
    all.insert(all.end(), std::make_move_iterator(r.begin()),
               std::make_move_iterator(r.end()));
  }

  std::sort(all.begin(), all.end(),
            [](const SearchResult &a, const SearchResult &b) {
              if (std::isnan(a.score)) return false;
              if (std::isnan(b.score)) return true;
              return a.score < b.score;
            });

  // diversity-prune across seeds so we do not return near-duplicate conformers
  double thr = params.diversityRmsThresh;
  if (thr <= 0.0) {
    if (profiling()) prof().nOut += static_cast<long long>(all.size());
    return all;
  }

  long long td0 = profiling() ? nowNs() : 0;
  std::vector<SearchResult> kept;
  if (params.finalSymmetryDedup) {
    kept = std::move(all);
    ctx.symmetryDedupInPlace(
        kept, thr);  // symmetry-aware optimal-superposition (QCP)
  } else {
    for (auto &r : all) {
      bool dup = false;
      for (const auto &k : kept) {
        double sq = 0.0;
        for (size_t i = 0; i < r.coords.size(); ++i) {
          double d = r.coords[i] - k.coords[i];
          sq += d * d;
        }
        if (std::sqrt(sq / (r.coords.size() / 3)) < thr) {
          dup = true;
          break;
        }
      }
      if (!dup) {
        kept.push_back(std::move(r));
      }
    }
  }
  if (profiling()) {
    prof().tDiv += nowNs() - td0;
    prof().nOut += static_cast<long long>(kept.size());
  }
  return kept;
}

}  // namespace RDKit
