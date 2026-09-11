//
//  Copyright (c) 2026, Glysade Inc.
//
//  How reachable is a MISSED product from the one the search actually
//  returned?
//
//  A recovery failure is only interesting if the target was findable. Two very
//  different situations look identical in a recovery count:
//
//    * the target sits a couple of reagent substitutions away, along a path on
//      which the shape score rises -- a local search SHOULD have walked there,
//      and the failure is in the search;
//    * the target is separated by a valley: every single-reagent step toward
//      it makes the shape score worse, so no hill-climbing search can reach it
//      regardless of budget, and the failure is in the objective's landscape.
//
//  This walks from the returned product to the planted one through
//  single-reagent substitutions, choosing at each step the neighbour most
//  similar (2D fingerprint) to the target, and reports the SHAPE score along
//  the way.  The product space is far too large to enumerate, but the walk
//  only ever expands sum(reagents) neighbours per step, which is cheap.
//
#include <RDGeneral/RDLog.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <DataStructs/BitOps.h>
#include <Confgen/SynthonSearch/EnumerateSynthons3D.h>
#include <Confgen/SynthonSearch/SynthonSearch3D.h>

#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace RDKit;

namespace {

void usage(const char *argv0) {
  std::cerr
      << "usage: " << argv0 << " --load LIB --pairs FILE [options]\n\n"
         "  --load FILE     a saved 3D synthon library\n"
         "  --pairs FILE    lines of 'planted found [score]', reagent indices\n"
         "                  comma separated -- e.g. '61,94,140 39,97,51'.\n"
         "                  The bench prints exactly these two columns.\n"
         "  --max-steps N   abandon a walk after N substitutions (default 6)\n"
         "  --query-pose N  which planted conformer is the query (default 0)\n"
         "  -h/--help\n\n"
         "  Reports, per pair: Hamming distance, whether an FP-guided walk\n"
         "  reaches the target, and whether the SHAPE score rises along it.\n";
}

std::vector<unsigned int> parseIdx(const std::string &s) {
  std::vector<unsigned int> v;
  std::stringstream ss(s);
  std::string cell;
  while (std::getline(ss, cell, ',')) {
    v.push_back(static_cast<unsigned int>(std::stoul(cell)));
  }
  return v;
}

std::unique_ptr<ExplicitBitVect> fingerprint(
    const EnumerateSynthons3D &lib, const std::vector<unsigned int> &idx,
    FingerprintGenerator<std::uint32_t> &gen) {
  ROMOL_SPTR m = lib.get2D(idx);
  if (!m) {
    return nullptr;
  }
  try {
    RWMol tmp(*m);
    MolOps::sanitizeMol(tmp);  // get2D returns a raw graph
    return std::unique_ptr<ExplicitBitVect>(gen.getFingerprint(tmp));
  } catch (...) {
    return nullptr;
  }
}

}  // namespace

int main(int argc, char *argv[]) {
  std::string libPath, pairPath;
  unsigned int maxSteps = 6, queryPose = 0;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto val = [&](const char *what) {
      if (i + 1 >= argc) {
        std::cerr << "missing value for " << what << "\n";
        std::exit(1);
      }
      return std::string(argv[++i]);
    };
    if (a == "--load") {
      libPath = val("--load");
    } else if (a == "--pairs") {
      pairPath = val("--pairs");
    } else if (a == "--max-steps") {
      maxSteps = static_cast<unsigned int>(std::stoul(val("--max-steps")));
    } else if (a == "--query-pose") {
      queryPose = static_cast<unsigned int>(std::stoul(val("--query-pose")));
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else {
      std::cerr << "unknown option " << a << "\n";
      usage(argv[0]);
      return 1;
    }
  }
  if (libPath.empty() || pairPath.empty()) {
    usage(argv[0]);
    return 1;
  }

  EnumerateSynthons3D lib;
  {
    std::ifstream in(libPath, std::ios_base::binary);
    if (!in) {
      std::cerr << "cannot open " << libPath << "\n";
      return 1;
    }
    lib.initFromStream(in);
  }
  if (!lib.isValid()) {
    std::cerr << "library did not load\n";
    return 1;
  }
  std::printf("[reach] %s: arity=%u reagents:", libPath.c_str(), lib.arity());
  for (unsigned int p = 0; p < lib.arity(); ++p) {
    std::printf(" %u", lib.numReagents(p));
  }
  std::printf("\n");

  std::unique_ptr<FingerprintGenerator<std::uint32_t>> gen(
      MorganFingerprint::getMorganGenerator<std::uint32_t>(2));

  std::ifstream pairs(pairPath);
  if (!pairs) {
    std::cerr << "cannot open " << pairPath << "\n";
    return 1;
  }

  unsigned int nPairs = 0, reached = 0, monotone = 0, blocked = 0;
  std::string line;
  std::printf(
      "\n  # ham  reach steps | shape(found) shape(target) minShapeOnPath  "
      "verdict\n");
  while (std::getline(pairs, line)) {
    if (line.empty()) {
      continue;
    }
    std::stringstream ls(line);
    std::string a, b;
    if (!(ls >> a >> b)) {
      continue;
    }
    const auto planted = parseIdx(a), found = parseIdx(b);
    if (planted.size() != lib.arity() || found.size() != lib.arity()) {
      continue;
    }
    ++nPairs;

    unsigned int hamming = 0;
    for (unsigned int p = 0; p < lib.arity(); ++p) {
      hamming += planted[p] != found[p];
    }

    // The query is one pose of the planted product, matching the benchmark.
    SynthonProduct target = lib.getProduct(planted);
    if (!target || !target.mol->getNumConformers()) {
      std::printf("  %2u  %u   (planted product would not build)\n", nPairs,
                  hamming);
      continue;
    }
    const unsigned int nc = target.mol->getNumConformers();
    const int qid =
        target.mol->getConformer(std::min(queryPose, nc - 1)).getId();
    ShapeScorer scorer(*target.mol, qid);

    auto shapeOf = [&](const std::vector<unsigned int> &idx) {
      SynthonProduct p = lib.getProduct(idx);
      if (!p) {
        return -1.0;
      }
      const auto v = scorer.score(*p.mol);
      return v ? *v : -1.0;
    };

    auto targetFp = fingerprint(lib, planted, *gen);
    if (!targetFp) {
      continue;
    }

    const double shapeFound = shapeOf(found);
    const double shapeTarget = shapeOf(planted);

    // Greedy walk: at each step take the single-reagent substitution whose
    // product is most similar to the target.  Knowing the destination is the
    // point -- this asks whether a path EXISTS that a local search could
    // follow, not whether the search would have guessed the direction.
    std::vector<unsigned int> cur = found;
    double minShape = shapeFound;
    unsigned int steps = 0;
    bool arrived = (cur == planted);
    bool dipped = false;
    while (!arrived && steps < maxSteps) {
      double bestSim = -1.0;
      std::vector<unsigned int> bestIdx;
      for (unsigned int p = 0; p < lib.arity(); ++p) {
        if (cur[p] == planted[p]) {
          continue;  // already correct at this position
        }
        for (unsigned int r = 0; r < lib.numReagents(p); ++r) {
          if (r == cur[p]) {
            continue;
          }
          auto cand = cur;
          cand[p] = r;
          auto fp = fingerprint(lib, cand, *gen);
          if (!fp) {
            continue;
          }
          const double sim = TanimotoSimilarity(*fp, *targetFp);
          if (sim > bestSim) {
            bestSim = sim;
            bestIdx = cand;
          }
        }
      }
      if (bestIdx.empty()) {
        break;
      }
      cur = bestIdx;
      ++steps;
      const double sh = shapeOf(cur);
      if (sh >= 0 && sh < minShape) {
        minShape = sh;
        dipped = true;
      }
      arrived = (cur == planted);
    }

    if (arrived) {
      ++reached;
    }
    // "Monotone" means a hill-climber on SHAPE could follow this path: the
    // score never drops below where it started.
    const bool mono = arrived && !dipped;
    if (mono) {
      ++monotone;
    } else if (arrived) {
      ++blocked;
    }
    std::printf(
        "  %2u  %u   %-5s %2u  |   %6.4f      %6.4f      %6.4f      %s\n",
        nPairs, hamming, arrived ? "yes" : "no", steps, shapeFound, shapeTarget,
        minShape,
        !arrived ? "unreached"
                 : (mono ? "REACHABLE (shape never dips)"
                         : "valley: shape drops en route"));
  }

  std::printf(
      "\n[reach] %u pairs: %u reached by an FP walk, of which %u are shape-"
      "monotone and %u cross a valley\n",
      nPairs, reached, monotone, blocked);
  std::printf(
      "  shape-monotone => a hill-climber COULD have walked there; the failure "
      "is in the search.\n"
      "  valley         => every step toward the target scores worse first; no "
      "local search reaches it.\n");
  return 0;
}
