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

#include <RDGeneral/RDLog.h>

#include <algorithm>
#include <mutex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <memory>
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
// fast as-is rms check to throw away confs in the same basin.
// fails as soon as possible.
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
      // Fast check RMS for removal fallback to the angle fingerprint
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

  // INCREMENTAL junction-local rescoring.
  // fragments are already scored so energy is InterFrag Energy + Fragments Energy
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
      // score the incremental changes, oldVdw and oldTor are invariant
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
          if (validate) { // debugging loop
            double full = d_driver.score();  // full inter-fragment sum at new coords
            double allowed = 1.0e-4 * (1.0 + std::fabs(full));
            if (std::fabs(sc - full) > allowed) {
              ++mismatches;
              // Validate the incremental versus full score
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

std::string TreeSearch::validateParams(const RigidRotorSearchParams &sp,
                                       const std::string &ffVariant) const {
  const std::string err = RigidRotorSearch::validateParams(sp, ffVariant);
  if (!err.empty()) {
    return err;
  }
  if (sp.tree.beamWidth == 0) {
    return "search.tree.beamWidth must be > 0";
  }
  if (sp.thompsonBudget > 0) {
    // Warn that no thompson will be done
    static std::once_flag warned;
    std::call_once(warned, [] {
      BOOST_LOG(rdWarningLog)
          << "search.thompsonBudget is ignored by RigidRotorSearchMode::Tree; "
             "use Thompson or Auto for a sampled search"
          << std::endl;
    });
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

  // diversity-prune across seeds
  double thr = params.diversityRmsThresh;
  if (thr <= 0.0) {
    if (profiling()) prof().nOut += static_cast<long long>(all.size());
    return all;
  }

  long long td0 = profiling() ? nowNs() : 0;
  std::vector<SearchResult> kept;
  if (params.finalSymmetryDedup) {
    kept = std::move(all);
    ctx.symmetryDedup(kept, thr); 
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
