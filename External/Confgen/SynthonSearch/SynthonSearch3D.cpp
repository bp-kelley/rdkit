//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SynthonSearch/SynthonSearch3D.h"

#include <GraphMol/GaussianShape/GaussianShape.h>
#include <RDGeneral/RDLog.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <RDGeneral/Exceptions.h>

#include <numeric>
#include <random>
#include <sstream>
#include <thread>

namespace RDKit {

// ---------------------------------------------------------------- scorer ---

ShapeScorer::ShapeScorer(const ROMol &query, int queryConfId,
                                       bool allCarbonRadii,
                                       double colorWeight)
    : d_allCarbonRadii(allCarbonRadii), d_colorWeight(colorWeight) {
  if (colorWeight < 0.0 || colorWeight > 1.0) {
    throw ValueErrorException(
        "ShapeScorer: colorWeight must be in [0, 1], got " +
        std::to_string(colorWeight));
  }
  RDLog::LogStateSetter blocker;
  // NOTE: confId -1 means ALL conformers of the query, not just the first --
  // the query is a multi-conformer shape.
  GaussianShape::ShapeInputOptions opts;
  opts.allCarbonRadii = allCarbonRadii;
  d_queryShape = std::make_unique<GaussianShape::ShapeInput>(
      query, queryConfId, opts);
}

ShapeScorer::~ShapeScorer() = default;

std::optional<double> ShapeScorer::score(ROMol &product) const {
  const auto s = scoreComponents(product);
  return s ? std::optional<double>((*s)[0]) : std::nullopt;
}

std::optional<std::array<double, 3>> ShapeScorer::scoreComponents(
    ROMol &product) const {
  if (!d_queryShape || !product.getNumConformers()) {
    return std::nullopt;
  }
  RDLog::LogStateSetter blocker;
  std::array<double, 3> best{-std::numeric_limits<double>::infinity(), 0.0,
                             0.0};
  // Iterate the conformers themselves: getConformer() takes an ID, not an
  // index, so indexing works only while the IDs happen to run 0..n-1.  A
  // molecule whose conformers were filtered, or one read from a file with
  // arbitrary IDs, would throw ConformerException here.
  for (auto ci = product.beginConformers(); ci != product.endConformers();
       ++ci) {
    const int confId = (*ci)->getId();
    try {
      // [0] combined shape+color, [1] shape, [2] color.  We rank on the
      // combination (modified tanimoto combo)
      GaussianShape::ShapeInputOptions opts;
      opts.allCarbonRadii = d_allCarbonRadii;
      const auto s = GaussianShape::AlignMolecule(
          *d_queryShape, product, opts, nullptr,
          GaussianShape::ShapeOverlayOptions(), confId);
      // The winner's parts travel with it
      auto w = s;
      if (d_colorWeight != 0.5) {
        w[0] = (1.0 - d_colorWeight) * s[1] + d_colorWeight * s[2];
      }
      if (w[0] > best[0]) {
        best = w;
      }
    } catch (...) {
      // one conformer failing to align is not the product failing to score
    }
  }
  if (!std::isfinite(best[0])) {
    return std::nullopt;
  }
  return best;
}

// ---------------------------------------- research-only: even-pair seeder ---
// Outside the production baseline-restart recipe. Retained to reproduce the
// alternate-seeder experiments and as scaffolding for possible pair models.


namespace {
//! Scores from one a particular reagent position, i.e. R1
/*!
  This is a cheap way to fill a posterior by examining one
  reagent at at time.
   `score[p][r]` is the best score seen for reagent `r` at
   position `p`; -inf means never scored (nothing buildable).
 
 Note that this is a guess, some reagents don't add shape
  linearly as they allow the whole manifold to move.
  In general, you can't take the highest score from R1, R2 and R3
   and expect it to make the best matching product, but it sure
   helps!
*/
struct SynthonArmScores {
  std::vector<std::vector<double>> score;
  std::vector<unsigned int> best;
  double bestScore = -std::numeric_limits<double>::infinity();
  unsigned int evaluations = 0;
  unsigned int unscorable = 0;
  double confgenMs = 0.0;
  double scoreMs = 0.0;
  bool empty() const { return score.empty(); }
};



using Reagents = std::vector<unsigned int>;

struct Timing {
  std::atomic<long long> confgenNs{0};
  std::atomic<long long> scoreNs{0};
};

long long nowNs() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

unsigned int resolveThreads(unsigned int requested) {
  if (requested) {
    return requested;
  }
  const auto available = std::thread::hardware_concurrency();
  return available ? available : 1u;
}

struct Request {
  unsigned int trajectory = 0;
  Reagents reagents;
};

struct BatchResult {
  std::vector<std::optional<double>> scores;
  std::vector<bool> cacheHit;
  std::vector<bool> crossTrajectoryCacheHit;
};

//! Persistent scalar cache plus one worker pool per combined batch.
class SharedTupleScorer {
 public:
  SharedTupleScorer(const EnumerateSynthons3D &lib,
                    const SynthonProductScorer &scorer, unsigned int numThreads)
      : d_lib(lib),
        d_scorer(scorer),
        d_numThreads(resolveThreads(numThreads)) {}

  BatchResult score(const std::vector<Request> &requests) {
    BatchResult result;
    result.scores.resize(requests.size());
    result.cacheHit.assign(requests.size(), false);
    result.crossTrajectoryCacheHit.assign(requests.size(), false);

    struct Work {
      Reagents reagents;
      unsigned int trajectory = 0;
      std::optional<double> score;
    };
    std::vector<Work> work;
    std::map<Reagents, size_t> pending;
    std::vector<std::optional<size_t>> workOf(requests.size());

    for (size_t i = 0; i < requests.size(); ++i) {
      const auto cached = d_cache.find(requests[i].reagents);
      if (cached != d_cache.end()) {
        result.scores[i] = cached->second.score;
        result.cacheHit[i] = true;
        result.crossTrajectoryCacheHit[i] =
            cached->second.trajectory != requests[i].trajectory;
        continue;
      }
      const auto queued = pending.find(requests[i].reagents);
      if (queued != pending.end()) {
        workOf[i] = queued->second;
        result.cacheHit[i] = true;  // same-batch in-flight deduplication
        result.crossTrajectoryCacheHit[i] =
            work[queued->second].trajectory != requests[i].trajectory;
        continue;
      }
      const size_t slot = work.size();
      pending.emplace(requests[i].reagents, slot);
      work.push_back(
          {requests[i].reagents, requests[i].trajectory, std::nullopt});
      workOf[i] = slot;
    }

    const unsigned int nThreads = std::min<unsigned int>(
        d_numThreads, static_cast<unsigned int>(work.size()));
    std::atomic<size_t> next{0};
    auto worker = [&]() {
      for (;;) {
        const size_t i = next++;
        if (i >= work.size()) {
          return;
        }
        const long long t0 = nowNs();
        long long t1 = t0;
        // A malformed synthon can throw out of assembly or scoring -- one
        // Freedom Space library contains an aromatic flag on a non-ring atom,
        // which raises AtomKekulizeException. An exception escaping a
        // std::thread calls std::terminate and kills the PROCESS, so a single
        // bad molecule takes down an entire multi-library sweep. Treat it the
        // way an unbuildable product is treated: no score, search continues.
        try {
          SynthonProduct product = d_lib.getProduct(work[i].reagents);
          t1 = nowNs();
          if (product) {
            work[i].score = d_scorer.score(*product.mol);
          }
        } catch (...) {
          if (t1 == t0) {
            t1 = nowNs();
          }
          work[i].score = std::nullopt;
        }
        d_timing.confgenNs += t1 - t0;
        d_timing.scoreNs += nowNs() - t1;
      }
    };
    if (nThreads <= 1) {
      worker();
    } else {
      std::vector<std::thread> pool;
      pool.reserve(nThreads);
      for (unsigned int i = 0; i < nThreads; ++i) {
        pool.emplace_back(worker);
      }
      for (auto &thread : pool) {
        thread.join();
      }
    }

    for (const auto &item : work) {
      d_cache.emplace(item.reagents, CacheEntry{item.score, item.trajectory});
      ++d_uniqueAttempts;
      if (item.score) {
        ++d_uniqueScored;
      } else {
        ++d_uniqueUnscorable;
      }
    }
    for (size_t i = 0; i < requests.size(); ++i) {
      if (workOf[i]) {
        result.scores[i] = work[*workOf[i]].score;
      }
    }
    return result;
  }

  unsigned int uniqueAttempts() const { return d_uniqueAttempts; }
  unsigned int uniqueScored() const { return d_uniqueScored; }
  unsigned int uniqueUnscorable() const { return d_uniqueUnscorable; }
  double confgenMs() const { return d_timing.confgenNs / 1.0e6; }
  double scoreMs() const { return d_timing.scoreNs / 1.0e6; }

 private:
  const EnumerateSynthons3D &d_lib;
  const SynthonProductScorer &d_scorer;
  unsigned int d_numThreads = 1;
  struct CacheEntry {
    std::optional<double> score;
    unsigned int trajectory = 0;
  };
  std::map<Reagents, CacheEntry> d_cache;
  Timing d_timing;
  unsigned int d_uniqueAttempts = 0;
  unsigned int d_uniqueScored = 0;
  unsigned int d_uniqueUnscorable = 0;
};

struct Trajectory {
  explicit Trajectory(unsigned int index, int seed, const EnumerateSynthons3D &lib)
      : rng(static_cast<unsigned int>(seed)), fixed(lib.arity(), -1) {
    summary.trajectory = index;
    summary.randomSeed = seed;
    arms.score.resize(lib.arity());
    for (unsigned int p = 0; p < lib.arity(); ++p) {
      arms.score[p].assign(lib.numReagents(p),
                           -std::numeric_limits<double>::infinity());
    }
  }

  std::mt19937 rng;
  std::vector<int> fixed;
  SynthonArmScores arms;
  SynthonTrajectorySummary summary;
  bool active = true;
  bool changed = false;
};

std::vector<unsigned int> sweepOrder(const EnumerateSynthons3D &lib,
                                     const SynthonSearch3DParams &params) {
  std::vector<unsigned int> order(lib.arity());
  std::iota(order.begin(), order.end(), 0u);
  if (params.largestFirst) {
    std::vector<double> meanAtoms(lib.arity(), 0.0);
    for (unsigned int p = 0; p < lib.arity(); ++p) {
      double total = 0.0;
      size_t count = 0;
      for (const auto &mol : lib.getReagents()[p]) {
        if (mol) {
          total += mol->getNumHeavyAtoms();
          ++count;
        }
      }
      meanAtoms[p] = count ? total / static_cast<double>(count) : 0.0;
    }
    std::stable_sort(order.begin(), order.end(),
                     [&](unsigned int a, unsigned int b) {
                       return meanAtoms[a] > meanAtoms[b];
                     });
  }
  if (params.firstPosition >= 0 &&
      static_cast<unsigned int>(params.firstPosition) < lib.arity()) {
    const auto first = static_cast<unsigned int>(params.firstPosition);
    const auto at = std::find(order.begin(), order.end(), first);
    if (at != order.end()) {
      std::rotate(order.begin(), at, at + 1);
    }
  }
  return order;
}

void noteObserved(std::map<Reagents, double> &observed,
                  const Reagents &reagents, double score) {
  const auto old = observed.find(reagents);
  if (old == observed.end()) {
    observed.emplace(reagents, score);
  } else if (score > old->second) {
    old->second = score;
  }
}

void assembleStats(const std::vector<Request> &requests, const BatchResult &batch,
             std::vector<Trajectory> &trajectories, std::map<Reagents, double> &observed,
             MultipleTrajectoryStats &stats) {
  for (size_t i = 0; i < requests.size(); ++i) {
    auto &summary = trajectories[requests[i].trajectory].summary;
    ++summary.logicalRequests;
    ++stats.logicalRequests;
    if (batch.cacheHit[i]) {
      ++summary.cacheHits;
      ++stats.cacheHits;
      if (batch.crossTrajectoryCacheHit[i]) {
        ++summary.crossTrajectoryCacheHits;
        ++stats.crossTrajectoryCacheHits;
      } else {
        ++stats.intraTrajectoryCacheHits;
      }
    }
    if (batch.scores[i]) {
      ++summary.logicalScored;
      ++stats.logicalScored;
      noteObserved(observed, requests[i].reagents, *batch.scores[i]);
    } else {
      ++summary.logicalUnscorable;
      ++stats.logicalUnscorable;
    }
  }
}

std::vector<unsigned int> topArms(const Trajectory &trajectory, unsigned int position,
                                  int requested) {
  std::vector<unsigned int> result;
  const auto &row = trajectory.arms.score[position];
  for (unsigned int r = 0; r < row.size(); ++r) {
    if (std::isfinite(row[r])) {
      result.push_back(r);
    }
  }
  const size_t keep = std::min<size_t>(
      result.size(), static_cast<size_t>(std::max(0, requested)));
  std::partial_sort(
      result.begin(), result.begin() + keep, result.end(),
      [&](unsigned int a, unsigned int b) { return row[a] > row[b]; });
  result.resize(keep);
  return result;
}

void materializeFinalists(const EnumerateSynthons3D &lib,
                          const SynthonProductScorer &scorer,
                          unsigned int numThreads,
                          const std::vector<std::pair<Reagents, double>> &top,
                          SynthonSearchResult &result, Timing &timing,
                          MultipleTrajectoryStats &stats) {
  struct Finalist {
    Reagents reagents;
    std::optional<double> score;
    ROMOL_SPTR mol;
  };
  std::vector<Finalist> finalists;
  finalists.reserve(top.size());
  for (const auto &item : top) {
    finalists.push_back({item.first, std::nullopt, ROMOL_SPTR()});
  }

  std::atomic<size_t> next{0};
  auto worker = [&]() {
    for (;;) {
      const size_t i = next++;
      if (i >= finalists.size()) {
        return;
      }
      const long long t0 = nowNs();
      long long t1 = t0;
      try {  // see the note in SharedTupleScorer: never let a thread throw
        SynthonProduct product = lib.getProduct(finalists[i].reagents);
        t1 = nowNs();
        if (product) {
          finalists[i].score = scorer.score(*product.mol);
          if (finalists[i].score) {
            finalists[i].mol = product.mol;
          }
        }
      } catch (...) {
        if (t1 == t0) {
          t1 = nowNs();
        }
        finalists[i].score = std::nullopt;
      }
      timing.confgenNs += t1 - t0;
      timing.scoreNs += nowNs() - t1;
    }
  };
  const unsigned int workers = std::min<unsigned int>(
      resolveThreads(numThreads), static_cast<unsigned int>(finalists.size()));
  if (workers <= 1) {
    worker();
  } else {
    std::vector<std::thread> pool;
    pool.reserve(workers);
    for (unsigned int i = 0; i < workers; ++i) {
      pool.emplace_back(worker);
    }
    for (auto &thread : pool) {
      thread.join();
    }
  }

  for (auto &finalist : finalists) {
    if (!finalist.score) {
      ++result.unscorable;
      continue;
    }
    ++result.evaluations;
    ++stats.finalistEvaluations;
    result.best.push_back({std::move(finalist.reagents), *finalist.score,
                           std::move(finalist.mol)});
  }
  std::sort(result.best.begin(), result.best.end(),
            [](const SynthonHit &a, const SynthonHit &b) {
              if (a.score != b.score) {
                return a.score > b.score;
              }
              return a.reagents < b.reagents;
            });
}

}  // namespace

std::string describeParams(const SynthonSearch3DParams &params,
                           const EnumerateSynthons3D &lib) {
  std::ostringstream os;
  os << std::boolalpha;
  const unsigned int threads =
      params.numThreads ? params.numThreads
                        : std::max(1u, std::thread::hardware_concurrency());
  os << "search.numTrajectories=" << params.numTrajectories << "\n"
     << "search.randomSeed=" << params.randomSeed << "\n"
     << "search.seedStride=" << params.seedStride << "\n"
     << "search.samplesPerReagent=" << params.samplesPerReagent << "\n"
     << "search.largestFirst=" << params.largestFirst << "\n"
     << "search.firstPosition=" << params.firstPosition << "\n"
     << "search.refineIters=" << params.refineIters << "\n"
     << "search.pairRefineTopK=" << params.pairRefineTopK << "\n"
     << "search.numThreads=" << threads
     << (params.numThreads ? "" : " (resolved from hardware concurrency)")
     << "\n"
     << "search.numBestProducts=" << params.numBestProducts << "\n"
     << "search.keepScoreThreshold=" << params.keepScoreThreshold
     << (params.keepScoreThreshold > Disabled ? "" : " (disabled, top-N only)")
     << "\n"
     << "search.maxKeptProducts=" << params.maxKeptProducts
     << (params.maxKeptProducts > 0 ? "" : " (no ceiling)") << "\n";

  const bool on = params.queryHeavyAtoms && params.pruneMinimum > 0.0;
  os << "prune.enabled=" << on << "\n"
     << "prune.queryHeavyAtoms=" << params.queryHeavyAtoms << "\n"
     << "prune.minimum=" << params.pruneMinimum << "\n"
     << "prune.maximum=" << params.pruneMaximum << "\n";
  if (on) {
    const double cutoff =
        std::min(params.pruneMinimum, params.pruneMaximum);
    os << "prune.appliesTo=refine,pairCrossProduct (exact product size)\n"
       << "prune.reagentBoundAppliesTo=sweep,refine (exact reagent bound)\n"
       << "prune.dynamicSource=bestFound\n"
       << "prune.effectiveInitial=" << cutoff << "\n"
       << "prune.boundType=heavyAtomHeuristic (false-negative rate unmeasured)"
       << "\n";
    const auto range = lib.productSizeRange();
    os << "prune.libraryProductHeavyAtoms=[" << range.first << ", "
       << range.second << "]\n";
    // A query inside the library's size band can always reach ratio 1.0, so
    // no product is filtered on size alone until a real score raises the
    // cutoff.  Outside it, the whole library is bounded below 1.0 and the
    // ceiling is worth stating: it is the best score the filter will admit.
    if (params.queryHeavyAtoms < range.first ||
        params.queryHeavyAtoms > range.second) {
      const double nearest = params.queryHeavyAtoms < range.first
                                 ? static_cast<double>(range.first)
                                 : static_cast<double>(range.second);
      const double q = params.queryHeavyAtoms;
      os << "prune.libraryBound=" << std::min(q, nearest) / std::max(q, nearest)
         << " (query outside library size range)\n";
    } else {
      os << "prune.libraryBound=1 (query inside library size range)\n";
    }
  }
  return os.str();
}

SynthonSearchResult synthonSearch3D(
    const EnumerateSynthons3D &lib, const SynthonProductScorer &scorer,
    const SynthonSearch3DParams &params, MultipleTrajectoryStats *outStats) {
  BOOST_LOG(rdInfoLog) << describeParams(params, lib);
  // The size filter bounds SHAPE only -- it is a heavy-atom ratio, we
  //  can't do the same for color
  if (params.queryHeavyAtoms && params.pruneMinimum > 0.0) {
    if (const auto *shape = dynamic_cast<const ShapeScorer *>(&scorer)) {
      if (shape->colorWeight() > 0.5) {
        BOOST_LOG(rdWarningLog)
            << "prune.boundVsColorWeight=MISMATCH: colorWeight="
            << shape->colorWeight()
            << " weights color above shape, but the size filter bounds SHAPE "
               "only (heavy-atom ratio); it cannot bound color, so pruned "
               "candidates may be false negatives.  Disable pruning "
               "(pruneMinimum=Disabled) for a color-weighted search you want "
               "to trust.\n";
      }
    }
  }

  MultipleTrajectoryStats stats;
  SynthonSearchResult result;
  if (!lib.isValid() || !lib.arity() || !params.numTrajectories) {
    if (outStats) {
      *outStats = std::move(stats);
    }
    return result;
  }

  std::vector<Trajectory> trajectories;
  trajectories.reserve(params.numTrajectories);
  for (unsigned int i = 0; i < params.numTrajectories; ++i) {
    trajectories.emplace_back(
        i, params.randomSeed + static_cast<int>(i) * params.seedStride, lib);
  }

  SharedTupleScorer sharedScorer(lib, scorer, params.numThreads);
  // Skip candidates the heavy-atom size filter rejects BEFORE assembling
  // them: product size is a sum of cached per-synthon counts, which is free
  // against the ~5 ms an assembly costs.  The cutoff is the running best
  // clamped into [pruneMinimum, pruneMaximum] -- the minimum discards hopeless
  // candidates before any hit exists, the maximum stops one strong hit from
  // suppressing everything else still worth reporting.
  //
  // Applied ONLY to refinement and the pair cross product, never to the
  // position sweeps.  The distinction is not about cost, it is about whether
  // dropping a candidate can change the answer:
  //
  //   refine / pair cross product -- every candidate is a COMPLETE assignment
  //     competing against the incumbent.  One the bound puts below the cutoff
  //     cannot beat the incumbent, so skipping it is outcome-neutral, and the
  //     topK x topK cross product is the largest batch in the search.
  //
  //   position sweep -- candidates are PARTIAL assignments with random
  //     completions, and the sweep takes an argmax over them.  A low score
  //     there is a stepping stone, not a reject; removing entries moves the
  //     argmax and sends the trajectory into a different basin.  The argument
  //     is structural -- it has NOT been demonstrated experimentally, because
  //     search-to-search score noise on this library (~0.01, from wall-clock
  //     confgen timeouts under varying thread load) is larger than any prune
  //     effect measured so far.
  // Best ratio a reagent could reach under ANY completion of the other
  // positions.  productHeavyCount is an exact sum over positions, so fixing
  // one synthon leaves the product size in a known interval and this bound is
  // exact, not heuristic: if it falls under the cutoff then NO combination of
  // the other reagents passes, and the reagent cannot appear in any admissible
  // product.  That makes it safe in the sweep, where per-candidate filtering
  // is not -- we drop a reagent only when the whole column is impossible,
  // never merely because one random completion looked small.
  auto reagentBound = [&](unsigned int position, unsigned int reagent) {
    const unsigned int h = lib.fragmentHeavyCount(position, reagent);
    if (!h) {
      return 1.0;  // a null synthon: no size information to bound with
    }
    double lo = h, hi = h;
    for (unsigned int q = 0; q < lib.arity(); ++q) {
      if (q == position) {
        continue;
      }
      const auto r = lib.positionSizeRange(q);
      lo += r.first;
      hi += r.second;
    }
    const double nq = params.queryHeavyAtoms;
    if (nq < lo) {
      return nq / lo;
    }
    if (nq > hi) {
      return hi / nq;
    }
    return 1.0;  // the query size is reachable, so the ratio can reach 1
  };

  // Cutoff shared by both filters: the running best clamped into the
  // [minimum, maximum] band.
  auto cutoff = [&]() {
    const double best =
        result.score > -std::numeric_limits<double>::max() ? result.score : 0.0;
    return std::min(std::max(best, params.pruneMinimum), params.pruneMaximum);
  };

  const bool pruning = params.queryHeavyAtoms && params.pruneMinimum > 0.0;

  //! True when this reagent cannot contribute to any admissible product.
  auto reagentImpossible = [&](unsigned int position, unsigned int reagent) {
    // Dynamic programming - skip work whos result is known
    if (lib.synthonUnusable(position, reagent)) {
      ++stats.deadReagentsSkipped;
      return true;
    }
    if (!pruning) {
      return false;
    }
    if (reagentBound(position, reagent) < cutoff()) {
      ++stats.reagentsFiltered;
      return true;
    }
    return false;
  };

  //! Psuedo shape-filter based on number of hvy atoms
  auto sizeFilterOne = [&](const Reagents &reagents) {
    if (!pruning) {
      return false;
    }
    const unsigned int np = lib.productHeavyCount(reagents);
    if (!np) {
      return false;
    }
    const double lo = std::min<double>(params.queryHeavyAtoms, np);
    const double hi = std::max<double>(params.queryHeavyAtoms, np);
    if (hi > 0.0 && lo / hi < cutoff()) {
      ++stats.sizeFiltered;
      return true;
    }
    return false;
  };

  auto sizeFilter = [&](std::vector<Request> &requests) {
    if (!pruning || requests.empty()) {
      return;
    }
    const double bar = cutoff();
    std::vector<Request> kept;
    kept.reserve(requests.size());
    for (auto &r : requests) {
      const unsigned int np = lib.productHeavyCount(r.reagents);
      const double lo = std::min<double>(params.queryHeavyAtoms, np);
      const double hi = std::max<double>(params.queryHeavyAtoms, np);
      if (np && hi > 0.0 && lo / hi < bar) {
        ++stats.sizeFiltered;
        continue;
      }
      kept.push_back(std::move(r));
    }
    requests.swap(kept);
  };


  std::map<Reagents, double> observed;
  const auto order = sweepOrder(lib, params);
  const unsigned int samples = std::max(1u, params.samplesPerReagent);

  // Greedy trajectories advance in lockstep, but every trajectory owns its RNG and fixed
  // choices. Interleaving changes only how physical work is scheduled.
  for (unsigned int step = 0; step < lib.arity(); ++step) {
    const unsigned int position = order[step];
    // Cost of ONE position sweep -- the unit greedy is made of.  Snapshot the
    // shared scorer's counters so the split is this position's, not cumulative.
    const long long sweepT0 = nowNs();
    const double confgenBefore = sharedScorer.confgenMs();
    const double scoreBefore = sharedScorer.scoreMs();
    const unsigned int uniqueBefore = sharedScorer.uniqueAttempts();
    const unsigned int cacheBefore = stats.cacheHits;
    std::vector<Request> requests;
    for (unsigned int trajectoryIndex = 0; trajectoryIndex < trajectories.size(); ++trajectoryIndex) {
      auto &trajectory = trajectories[trajectoryIndex];
      for (unsigned int reagent = 0; reagent < lib.numReagents(position);
           ++reagent) {
        if (reagentImpossible(position, reagent)) {
          continue;
        }
        // Note:  don't filter away every sample, otherwise we can't
	//  properly evaluate this window or, worse, entirely rely on the
	//  size estimate which isn't really shape.
        size_t kept = 0;
        Reagents first;
        for (unsigned int sample = 0; sample < samples; ++sample) {
          Reagents candidate(lib.arity(), 0);
          for (unsigned int p = 0; p < lib.arity(); ++p) {
            if (p == position) {
              candidate[p] = reagent;
            } else if (trajectory.fixed[p] >= 0) {
              candidate[p] = static_cast<unsigned int>(trajectory.fixed[p]);
            } else {
              candidate[p] = std::uniform_int_distribution<unsigned int>(
                  0, lib.numReagents(p) - 1)(trajectory.rng);
            }
          }
          if (sample == 0) {
            first = candidate;
          }
          if (sizeFilterOne(candidate)) {
            continue;  // cannot reach the cutoff: skip assembly AND scoring
          }
          requests.push_back({trajectoryIndex, std::move(candidate)});
          ++kept;
        }
        if (!kept) {
          requests.push_back({trajectoryIndex, std::move(first)});
        }
      }
    }

    // Deliberately NOT size-filtered: see sizeFilter's comment.
    //  We still want an argmax of the samples to keep a realistic
    //  estimate
    const auto batch = sharedScorer.score(requests);
    assembleStats(requests, batch, trajectories, observed, stats);
    {
      SynthonPositionSweep sweep;
      sweep.position = position;
      sweep.reagents = lib.numReagents(position);
      sweep.evaluations = sharedScorer.uniqueAttempts() - uniqueBefore;
      sweep.cacheHits = stats.cacheHits - cacheBefore;
      sweep.wallMs = (nowNs() - sweepT0) / 1.0e6;
      sweep.confgenMs = sharedScorer.confgenMs() - confgenBefore;
      sweep.scoreMs = sharedScorer.scoreMs() - scoreBefore;
      stats.positionSweeps.push_back(sweep);
    }
    std::vector<std::vector<std::vector<double>>> samplesByTrajectory(
        trajectories.size(),
        std::vector<std::vector<double>>(lib.numReagents(position)));
    for (size_t i = 0; i < requests.size(); ++i) {
      if (batch.scores[i]) {
        samplesByTrajectory[requests[i].trajectory][requests[i].reagents[position]]
            .push_back(*batch.scores[i]);
      }
    }
    for (unsigned int trajectoryIndex = 0; trajectoryIndex < trajectories.size(); ++trajectoryIndex) {
      auto &trajectory = trajectories[trajectoryIndex];
      double best = -std::numeric_limits<double>::infinity();
      int bestReagent = -1;
      for (unsigned int reagent = 0; reagent < lib.numReagents(position);
           ++reagent) {
        const auto &values = samplesByTrajectory[trajectoryIndex][reagent];
        if (values.empty()) {
          continue;
        }
        const double value = *std::max_element(values.begin(), values.end());
        trajectory.arms.score[position][reagent] = value;
        if (value > best) {
          best = value;
          bestReagent = static_cast<int>(reagent);
        }
      }
      if (bestReagent >= 0) {
        trajectory.fixed[position] = bestReagent;
      }
    }
  }

  // Use the same per-position argmax reconstruction as the measured baseline.
  // Under Max it has the score of an actual best sampled product; preserving
  // this detail makes each trajectory comparable with an independent legacy call.
  std::vector<Request> seeds;
  for (unsigned int trajectoryIndex = 0; trajectoryIndex < trajectories.size(); ++trajectoryIndex) {
    Reagents seed(lib.arity(), 0);
    for (unsigned int p = 0; p < lib.arity(); ++p) {
      const auto &row = trajectories[trajectoryIndex].arms.score[p];
      seed[p] = static_cast<unsigned int>(
          std::max_element(row.begin(), row.end()) - row.begin());
    }
    trajectories[trajectoryIndex].summary.reagents = seed;
    seeds.push_back({trajectoryIndex, std::move(seed)});
  }
  {
    const auto batch = sharedScorer.score(seeds);
    assembleStats(seeds, batch, trajectories, observed, stats);
    for (size_t i = 0; i < seeds.size(); ++i) {
      if (batch.scores[i]) {
        trajectories[i].summary.score = *batch.scores[i];
      }
    }
  }

  for (unsigned int iteration = 0; iteration < params.refineIters;
       ++iteration) {
    for (auto &trajectory : trajectories) {
      trajectory.changed = false;
    }

    // Coordinate refinement remains sequential within a trajectory, while the same
    // coordinate from all active trajectories shares one physical batch.
    for (unsigned int position = 0; position < lib.arity(); ++position) {
      std::vector<Request> requests;
      for (unsigned int trajectoryIndex = 0; trajectoryIndex < trajectories.size(); ++trajectoryIndex) {
        const auto &trajectory = trajectories[trajectoryIndex];
        if (!trajectory.active) {
          continue;
        }
        for (unsigned int reagent = 0; reagent < lib.numReagents(position);
             ++reagent) {
          if (reagentImpossible(position, reagent)) {
            continue;
          }
          auto candidate = trajectory.summary.reagents;
          candidate[position] = reagent;
          requests.push_back({trajectoryIndex, std::move(candidate)});
        }
      }
      sizeFilter(requests);
      const auto batch = sharedScorer.score(requests);
      assembleStats(requests, batch, trajectories, observed, stats);
      std::vector<double> bestScore(trajectories.size(),
                                    -std::numeric_limits<double>::infinity());
      std::vector<Reagents> bestCandidate(trajectories.size());
      for (size_t i = 0; i < requests.size(); ++i) {
        if (batch.scores[i] &&
            *batch.scores[i] > bestScore[requests[i].trajectory]) {
          bestScore[requests[i].trajectory] = *batch.scores[i];
          bestCandidate[requests[i].trajectory] = requests[i].reagents;
        }
      }
      for (unsigned int trajectoryIndex = 0; trajectoryIndex < trajectories.size(); ++trajectoryIndex) {
        auto &trajectory = trajectories[trajectoryIndex];
        if (trajectory.active && bestScore[trajectoryIndex] > trajectory.summary.score) {
          trajectory.summary.score = bestScore[trajectoryIndex];
          trajectory.summary.reagents = std::move(bestCandidate[trajectoryIndex]);
          trajectory.changed = true;
        }
      }
    }

    // Preserve the production KxK pair-refine after coordinate descent stalls.
    if (params.pairRefineTopK > 0) {
      for (unsigned int first = 0; first < lib.arity(); ++first) {
        for (unsigned int second = first + 1; second < lib.arity(); ++second) {
          std::vector<Request> requests;
          for (unsigned int trajectoryIndex = 0; trajectoryIndex < trajectories.size();
               ++trajectoryIndex) {
            const auto &trajectory = trajectories[trajectoryIndex];
            if (!trajectory.active || trajectory.changed) {
              continue;
            }
            const auto a = topArms(trajectory, first, params.pairRefineTopK);
            const auto b = topArms(trajectory, second, params.pairRefineTopK);
            for (const unsigned int ra : a) {
              for (const unsigned int rb : b) {
                if (ra == trajectory.summary.reagents[first] &&
                    rb == trajectory.summary.reagents[second]) {
                  continue;
                }
                auto candidate = trajectory.summary.reagents;
                candidate[first] = ra;
                candidate[second] = rb;
                requests.push_back({trajectoryIndex, std::move(candidate)});
              }
            }
          }
          sizeFilter(requests);
      const auto batch = sharedScorer.score(requests);
          assembleStats(requests, batch, trajectories, observed, stats);
          std::vector<double> bestScore(
              trajectories.size(), -std::numeric_limits<double>::infinity());
          std::vector<Reagents> bestCandidate(trajectories.size());
          for (size_t i = 0; i < requests.size(); ++i) {
            if (batch.scores[i] &&
                *batch.scores[i] > bestScore[requests[i].trajectory]) {
              bestScore[requests[i].trajectory] = *batch.scores[i];
              bestCandidate[requests[i].trajectory] = requests[i].reagents;
            }
          }
          for (unsigned int trajectoryIndex = 0; trajectoryIndex < trajectories.size();
               ++trajectoryIndex) {
            auto &trajectory = trajectories[trajectoryIndex];
            if (trajectory.active && !trajectory.changed &&
                bestScore[trajectoryIndex] > trajectory.summary.score) {
              trajectory.summary.score = bestScore[trajectoryIndex];
              trajectory.summary.reagents = std::move(bestCandidate[trajectoryIndex]);
              trajectory.changed = true;
            }
          }
        }
      }
    }

    bool anyActive = false;
    for (auto &trajectory : trajectories) {
      if (trajectory.active && !trajectory.changed) {
        trajectory.active = false;
      }
      anyActive = anyActive || trajectory.active;
    }
    if (!anyActive) {
      break;
    }
  }

  stats.uniqueSearchAttempts = sharedScorer.uniqueAttempts();
  stats.uniqueSearchScored = sharedScorer.uniqueScored();
  stats.uniqueSearchUnscorable = sharedScorer.uniqueUnscorable();
  stats.trajectories.reserve(trajectories.size());
  for (auto &trajectory : trajectories) {
    stats.trajectories.push_back(std::move(trajectory.summary));
  }

  std::vector<std::pair<Reagents, double>> ranked(observed.begin(),
                                                  observed.end());
  std::sort(ranked.begin(), ranked.end(), [](const auto &a, const auto &b) {
    if (a.second != b.second) {
      return a.second > b.second;
    }
    return a.first < b.first;
  });
  if (!ranked.empty()) {
    result.reagents = ranked.front().first;
    result.score = ranked.front().second;
  }

  result.evaluations = stats.uniqueSearchScored;
  result.unscorable = stats.uniqueSearchUnscorable;
  result.cacheHits = stats.cacheHits;
  result.confgenMs = sharedScorer.confgenMs();
  result.scoreMs = sharedScorer.scoreMs();

  Timing finalistTiming;
  // How many of `ranked` to report/save.  `ranked` is sorted best-first, so
  // both limits are prefix lengths and the threshold scan is a single walk.
  size_t keep = params.numBestProducts > 0
                    ? static_cast<size_t>(params.numBestProducts)
                    : 0;
  if (params.keepScoreThreshold > Disabled) {
    size_t above = 0;
    while (above < ranked.size() &&
           ranked[above].second >= params.keepScoreThreshold) {
      ++above;
    }
    // max, not assignment: the threshold only ever adds, so a run that clears
    // it nowhere still reports its best products.
    keep = std::max(keep, above);
  }
  if (params.maxKeptProducts > 0) {
    keep = std::min(keep, static_cast<size_t>(params.maxKeptProducts));
  }
  if (keep > 0 && !ranked.empty()) {
    ranked.resize(std::min<size_t>(ranked.size(), keep));
    materializeFinalists(lib, scorer, params.numThreads, ranked, result,
                         finalistTiming, stats);
    result.confgenMs += finalistTiming.confgenNs / 1.0e6;
    result.scoreMs += finalistTiming.scoreNs / 1.0e6;
  }

  if (outStats) {
    *outStats = std::move(stats);
  }
  return result;
}

}  // namespace RDKit
