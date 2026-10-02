//
//  Copyright (c) 2026, Glysade Inc.
//
//  Did the search miss a better product than the one it returned?
//
//  A search result is only interpretable against the library it searched.
//  synthonSearch3D touches a few percent of a large library, so a returned
//  score answers "what did the greedy trajectories find" and not "what is in
//  here".  This tool attacks the second question two ways, because they fail
//  differently:
//
//    * COLUMN SCAN -- exhaustively score every reagent at each position with
//      the other positions pinned to the returned answer.  Cheap (sum of
//      reagents, not their product) and decisive about one thing: whether the
//      answer is a coordinate-wise local maximum.  If a column scan beats it,
//      refinement simply failed and the search has a bug.
//
//    * RANDOM SCAN -- score uniformly random products.  A column scan cannot
//      see a better product in a different basin; a random sample can, and it
//      also calibrates the returned score against the library's own score
//      distribution.  Reported as a count of samples beating the reference,
//      which extrapolates to an expected number of better products library
//      wide.
//
//  Both report the best product found, so a genuine miss can be inspected.
//
#include <RDGeneral/RDLog.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <Confgen/SynthonSearch/EnumerateSynthons3D.h>
#include <Confgen/SynthonSearch/SynthonSearch3D.h>

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <set>
#include <string>
#include <thread>
#include <vector>

using namespace RDKit;

namespace {

void usage(const char *argv0) {
  std::cerr
      << "usage: " << argv0 << " --load LIB --query SMILES [options]\n\n"
         "  --load FILE      a saved 3D synthon library\n"
         "  --query SMILES   the query molecule\n"
         "  --anchor A,B,..  reagent indices to scan the columns through;\n"
         "                   omit to skip the column scan\n"
         "  --random N       score N uniformly random products (default 0)\n"
         "  --dead-synthons  report synthons the cache knows are unembeddable\n"
         "  --reference S    report how many products beat score S\n"
         "  --threads N      worker threads (default hardware concurrency)\n"
         "  --seed N         RNG seed for the random scan (default 0xf00d)\n"
         "  --query-confs N  query conformers to embed (default 10)\n"
         "  -h/--help\n";
}

std::vector<unsigned int> parseIdx(const std::string &s) {
  std::vector<unsigned int> v;
  std::stringstream ss(s);
  std::string cell;
  while (std::getline(ss, cell, ',')) {
    if (!cell.empty()) {
      v.push_back(static_cast<unsigned int>(std::stoul(cell)));
    }
  }
  return v;
}

std::string idxStr(const std::vector<unsigned int> &v) {
  std::string s;
  for (size_t i = 0; i < v.size(); ++i) {
    s += (i ? "," : "") + std::to_string(v[i]);
  }
  return s;
}

//! Score a batch of products across threads, returning the best seen.
struct ScanResult {
  double best = -1.0;
  std::vector<unsigned int> bestReagents;
  std::string bestSmiles;
  unsigned int scored = 0;
  unsigned int unscorable = 0;
  std::vector<double> scores;
  //! Why the unscorable ones failed, indexed by SynthonBuildStatus, plus a
  //! trailing bucket for products that BUILT but the scorer still rejected.
  std::vector<unsigned int> status{std::vector<unsigned int>(6, 0)};
  unsigned int fellBackToFull = 0;
  std::vector<std::vector<unsigned int>> failures;
};

ScanResult scanBatch(const EnumerateSynthons3D &lib,
                     const SynthonProductScorer &scorer,
                     const std::vector<std::vector<unsigned int>> &batch,
                     unsigned int numThreads, const char *label) {
  ScanResult out;
  out.scores.reserve(batch.size());
  std::mutex mtx;
  std::atomic<size_t> next{0};
  std::atomic<size_t> done{0};

  auto worker = [&]() {
    try {
      for (;;) {
        const size_t i = next.fetch_add(1);
        if (i >= batch.size()) {
          return;
        }
        std::optional<double> s;
        std::string smi;
        int st = 5;  // "built but the scorer declined it"
        bool fallback = false;
        try {
          auto prod = lib.getProduct(batch[i]);
          if (prod.status != SynthonBuildStatus::Ok) {
            st = static_cast<int>(prod.status);
          }
          fallback = prod.usedCoarseFallback;
          if (prod.mol) {
            s = scorer.score(*prod.mol);
            if (s) {
              smi = MolToSmiles(*prod.mol);
            }
          }
        } catch (const std::exception &) {
          st = 2;  // an assembly throw is a zip failure by another name
        }
        {
          std::lock_guard<std::mutex> lock(mtx);
          if (fallback) {
            ++out.fellBackToFull;
          }
          if (s) {
            out.scores.push_back(*s);
            ++out.scored;
            if (*s > out.best) {
              out.best = *s;
              out.bestReagents = batch[i];
              out.bestSmiles = smi;
            }
          } else {
            ++out.unscorable;
            ++out.status[st];
            if (out.failures.size() < 20) {
              out.failures.push_back(batch[i]);
            }
          }
        }
        const size_t n = done.fetch_add(1) + 1;
        if (n % 2000 == 0) {
          std::lock_guard<std::mutex> lock(mtx);
          std::cerr << "  [" << label << "] " << n << "/" << batch.size()
                    << " best=" << out.best << "\r" << std::flush;
        }
      }
    } catch (const std::exception &e) {
      std::cerr << "\n  worker died: " << e.what() << "\n";
    }
  };

  std::vector<std::thread> pool;
  for (unsigned int t = 0; t < numThreads; ++t) {
    pool.emplace_back(worker);
  }
  for (auto &t : pool) {
    t.join();
  }
  std::cerr << "                                                        \r";
  return out;
}

void report(const char *label, const ScanResult &r, double reference) {
  std::printf("\n  %s\n", label);
  std::printf("    scored %u, unscorable %u   (coarse->full fallbacks: %u)\n",
              r.scored, r.unscorable, r.fellBackToFull);
  if (r.unscorable) {
    static const char *names[6] = {"Ok",
                                   "BadReagentIndex",
                                   "ZipFailed",
                                   "NoReagentConfs (MMFF cannot type it)",
                                   "ConfGenFailed (coarse AND full found none)",
                                   "built, scorer declined"};
    std::printf("    why unscorable:\n");
    for (int i = 0; i < 6; ++i) {
      if (r.status[i]) {
        std::printf("      %-44s %u\n", names[i], r.status[i]);
      }
    }
    std::printf("    example failing reagents:");
    for (size_t i = 0; i < r.failures.size() && i < 8; ++i) {
      std::printf(" [%s]", idxStr(r.failures[i]).c_str());
    }
    std::printf("\n");
  }
  if (!r.scored) {
    return;
  }
  auto s = r.scores;
  std::sort(s.begin(), s.end());
  auto pct = [&](double p) {
    return s[std::min(s.size() - 1,
                      static_cast<size_t>(p * (double)s.size()))];
  };
  std::printf("    p50 %.4f  p90 %.4f  p99 %.4f  p99.9 %.4f  max %.4f\n",
              pct(0.50), pct(0.90), pct(0.99), pct(0.999), r.best);
  if (reference > 0.0) {
    const size_t better =
        s.end() - std::upper_bound(s.begin(), s.end(), reference);
    std::printf("    beating %.4f: %zu of %u (%.4f%%)\n", reference, better,
                r.scored, 100.0 * (double)better / (double)r.scored);
  }
  std::printf("    best %.4f  [%s]\n", r.best, idxStr(r.bestReagents).c_str());
  std::printf("         %s\n", r.bestSmiles.c_str());
}

}  // namespace

int main(int argc, char *argv[]) {
  std::string loadFile, query, anchorStr;
  size_t randomN = 0;
  double reference = 0.0;
  unsigned int numThreads = std::max(1u, std::thread::hardware_concurrency());
  int seed = 0xf00d;
  int queryConfs = 10;
  bool deadSynthons = false;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto val = [&](const char *flag) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "missing value for " << flag << "\n";
        std::exit(1);
      }
      return argv[++i];
    };
    if (a == "--load") {
      loadFile = val("--load");
    } else if (a == "--query") {
      query = val("--query");
    } else if (a == "--anchor") {
      anchorStr = val("--anchor");
    } else if (a == "--dead-synthons") {
      deadSynthons = true;
    } else if (a == "--random") {
      randomN = static_cast<size_t>(std::stoul(val("--random")));
    } else if (a == "--reference") {
      reference = std::stod(val("--reference"));
    } else if (a == "--threads") {
      numThreads = static_cast<unsigned int>(std::stoul(val("--threads")));
    } else if (a == "--seed") {
      seed = std::stoi(val("--seed"));
    } else if (a == "--query-confs") {
      queryConfs = std::stoi(val("--query-confs"));
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else {
      std::cerr << "unknown argument: " << a << "\n";
      usage(argv[0]);
      return 1;
    }
  }
  if (loadFile.empty() || query.empty()) {
    usage(argv[0]);
    return 1;
  }

  EnumerateSynthons3D lib;
  {
    std::ifstream in(loadFile, std::ios_base::binary);
    if (!in) {
      std::cerr << "could not open " << loadFile << "\n";
      return 1;
    }
    lib.initFromStream(in);
  }
  if (!lib.isValid()) {
    std::cerr << "library is not valid\n";
    return 1;
  }

  std::vector<unsigned int> counts;
  double total = 1.0;
  for (unsigned int p = 0; p < lib.arity(); ++p) {
    counts.push_back(lib.numReagents(p));
    total *= (double)lib.numReagents(p);
  }
  std::printf("  library: arity %u, reagents [%s], %.4g products\n",
              lib.arity(), idxStr(counts).c_str(), total);
  const auto range = lib.productSizeRange();
  std::printf("  product heavy atoms: [%u, %u]\n", range.first, range.second);

  if (deadSynthons) {
    // A synthon whose cache entry is a TOMBSTONE can never yield a buildable
    // product, and the library already knows it: prefill attempted the embed
    // at build time and recorded the failure.  Nothing new is stored here, we
    // are only reading back a fact that was computed once and then ignored.
    const auto &fl = lib.embedder();
    if (!fl) {
      std::printf("\n  DEAD SYNTHONS: this library carries no fragment "
                  "cache, so there is nothing to read back\n");
    } else {
      std::printf("\n  DEAD SYNTHONS (cache: %zu entries, %zu unembeddable)\n",
                  fl->size(), fl->numUnembeddable());
      const bool coarse =
          isCoarseAssembly(lib.params3D().embedStyle);
      const std::set<std::string> exitSymbols(
          lib.molzipParams().atomSymbols.begin(),
          lib.molzipParams().atomSymbols.end());
      unsigned int grandDead = 0, grandUnknown = 0;
      for (unsigned int p = 0; p < lib.arity(); ++p) {
        unsigned int nDead = 0, nUnknown = 0, nLive = 0;
        std::vector<unsigned int> deadIdx;
        for (unsigned int r = 0; r < lib.numReagents(p); ++r) {
          const auto &syn = lib.getReagents()[p][r];
          if (!syn) {
            continue;
          }
          // Match prefill's key exactly: coarse embeds the WHOLE synthon with
          // its exits reduced to plain dummies.
          RWMol whole(*syn);
          if (coarse) {
            for (auto atom : whole.atoms()) {
              if (exitSymbols.count(atom->getSymbol())) {
                atom->setAtomicNum(0);
                atom->setIsotope(0);
                atom->setNoImplicit(true);
                atom->setNumExplicitHs(0);
              }
            }
          }
          try {
            MolOps::sanitizeMol(whole);
          } catch (...) {
            ++nUnknown;
            continue;
          }
          const auto n = fl->numFragmentConfs(whole);
          if (!n) {
            ++nUnknown;
          } else if (!*n) {
            ++nDead;
            if (deadIdx.size() < 12) {
              deadIdx.push_back(r);
            }
          } else {
            ++nLive;
          }
        }
        std::printf("    pos%u: %u live, %u DEAD, %u not in cache (of %u)\n", p,
                    nLive, nDead, nUnknown, lib.numReagents(p));
        if (!deadIdx.empty()) {
          std::printf("      dead indices:");
          for (const unsigned int r : deadIdx) {
            std::printf(" %u", r);
          }
          std::printf("\n");
          for (size_t i = 0; i < deadIdx.size() && i < 3; ++i) {
            std::printf("      e.g. [%u] %s\n", deadIdx[i],
                        MolToSmiles(*lib.getReagents()[p][deadIdx[i]]).c_str());
          }
        }
        grandDead += nDead;
        grandUnknown += nUnknown;
      }
      // What the dead synthons cost: every product containing one is built
      // twice (coarse then full) and discarded.
      double reachable = 1.0;
      for (unsigned int p = 0; p < lib.arity(); ++p) {
        unsigned int nDead = 0;
        for (unsigned int r = 0; r < lib.numReagents(p); ++r) {
          const auto &syn = lib.getReagents()[p][r];
          if (!syn) {
            continue;
          }
          RWMol whole(*syn);
          if (coarse) {
            for (auto atom : whole.atoms()) {
              if (exitSymbols.count(atom->getSymbol())) {
                atom->setAtomicNum(0);
                atom->setIsotope(0);
                atom->setNoImplicit(true);
                atom->setNumExplicitHs(0);
              }
            }
          }
          try {
            MolOps::sanitizeMol(whole);
          } catch (...) {
            continue;
          }
          const auto n = fl->numFragmentConfs(whole);
          if (n && !*n) {
            ++nDead;
          }
        }
        reachable *= (double)(lib.numReagents(p) - nDead) /
                     (double)lib.numReagents(p);
      }
      std::printf("    %u dead synthons, %u not in cache\n", grandDead,
                  grandUnknown);
      std::printf("    unbuildable products: %.2f%% of the library\n",
                  100.0 * (1.0 - reachable));
    }
  }

  std::unique_ptr<ROMol> q(SmilesToMol(query));
  if (!q) {
    std::cerr << "could not parse query\n";
    return 1;
  }
  RWMol qh(*q);
  MolOps::addHs(qh);
  DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
  ps.randomSeed = seed;
  DGeomHelpers::EmbedMultipleConfs(qh, queryConfs, ps);
  std::vector<std::pair<int, double>> mmffRes;
  MMFF::MMFFOptimizeMoleculeConfs(qh, mmffRes, 1, 500);
  std::printf("  query: %u heavy atoms, %u conformers\n", q->getNumHeavyAtoms(),
              qh.getNumConformers());
  ShapeScorer scorer(qh, qh.getNumConformers() ? 0 : -1);

  std::printf("  threads: %u   seed: %d   reference: %.4f\n", numThreads, seed,
              reference);

  if (!anchorStr.empty()) {
    const auto anchor = parseIdx(anchorStr);
    if (anchor.size() != lib.arity()) {
      std::cerr << "--anchor needs " << lib.arity() << " indices\n";
      return 1;
    }
    for (unsigned int p = 0; p < lib.arity(); ++p) {
      std::vector<std::vector<unsigned int>> batch;
      batch.reserve(lib.numReagents(p));
      for (unsigned int r = 0; r < lib.numReagents(p); ++r) {
        auto cand = anchor;
        cand[p] = r;
        batch.push_back(std::move(cand));
      }
      const std::string label = "column r" + std::to_string(p);
      const auto res = scanBatch(lib, scorer, batch, numThreads, label.c_str());
      report(("EXHAUSTIVE " + label + " through [" + idxStr(anchor) + "]")
                 .c_str(),
             res, reference);
    }
  }

  if (randomN) {
    std::mt19937 rng(static_cast<unsigned int>(seed));
    std::vector<std::vector<unsigned int>> batch;
    batch.reserve(randomN);
    for (size_t i = 0; i < randomN; ++i) {
      std::vector<unsigned int> cand(lib.arity());
      for (unsigned int p = 0; p < lib.arity(); ++p) {
        cand[p] = std::uniform_int_distribution<unsigned int>(
            0, lib.numReagents(p) - 1)(rng);
      }
      batch.push_back(std::move(cand));
    }
    const auto res = scanBatch(lib, scorer, batch, numThreads, "random");
    report("RANDOM SCAN", res, reference);
    if (reference > 0.0 && res.scored) {
      auto s = res.scores;
      std::sort(s.begin(), s.end());
      const size_t better =
          s.end() - std::upper_bound(s.begin(), s.end(), reference);
      if (better) {
        std::printf(
            "\n    extrapolated to the full library: ~%.0f of %.4g products "
            "beat %.4f\n",
            total * (double)better / (double)res.scored, total, reference);
      } else {
        // Zero hits does not mean zero exist.  The rule of three puts the
        // 95% upper bound on an unobserved rate at 3/n, which is the only
        // defensible number to quote from an empty sample.
        std::printf(
            "\n    none of %u samples beat %.4f; 95%% upper bound is "
            "3/%u, so AT MOST ~%.0f of %.4g products do (<%.3f%%)\n",
            res.scored, reference, res.scored,
            total * 3.0 / (double)res.scored, total,
            100.0 * 3.0 / (double)res.scored);
      }
    }
  }
  return 0;
}
