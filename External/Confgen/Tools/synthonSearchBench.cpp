//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  synthonSearchBench -- recovery benchmark for the 3D synthon search.
//
//  THE EXPERIMENT: plant a product, then try to find it again.
//    1. build (or load) a 3D synthon library for one reaction
//    2. pick a product at random -- that IS the ground truth
//    3. generate its conformers and use them as the shape QUERY
//    4. run the selected synthon search against that query
//    5. did the search recover the planted reagent indices?
//
//  A planted product is the only query whose right answer is known, so
//  recovery rate is measurable rather than a matter of opinion.  Note the
//  search can legitimately return a DIFFERENT product that overlays the query
//  as well or better -- that is reported separately from a miss, because a
//  shape search is not obliged to prefer the molecule we happened to plant.
//
//  No data ships with this tool: point it at a synthon CSV of the form
//    SMILES,synton_id,synton_role,reaction_id
//  with [U] exit vectors (see External/Confgen/README.md).
//
#include <Confgen/SynthonSearch/EnumerateSynthons3D.h>
#include <Confgen/SynthonSearch/SynthonSearch3D.h>
#include <Confgen/Utils/ParamsIO.h>

#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <functional>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace RDKit;

namespace {

double nowMs() {
  return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

void usage(const char *argv0) {
  std::cerr
      << "usage: " << argv0
      << " [options]\n\n"
         "  Library (one of):\n"
         "    --synthons FILE --rxn ID  build from a synthon table (CSV or\n"
         "                          TSV; columns found from the header:\n"
         "                          SMILES, synton role/#, reaction id)\n"
         "    --load FILE           load a previously saved 3D library\n\n"
         "  Options:\n"
         "    --save FILE           save the built library (with its fragment\n"
         "                          cache) for reuse\n"
         "    --max-per-pos N       cap reagents per position (0 = all)\n"
         "    --no-prefill          do NOT pre-embed every synthon fragment\n"
         "    --coarse              cut ONLY at the synthon junctions: each\n"
         "                          synthon stays whole, so far fewer rotors\n"
         "    --product-confs N     conformers generated PER PRODUCT.  This is\n"
         "                          the dominant cost: assembly is 94-99% of\n"
         "                          the search.  Low = a coarse sample.\n"
         "    --trials N            planted products to search for (default 5)\n"
         "    --budget N            evaluations per search (-1 AUTO, 0 off)\n"
         "    --noise-var F         observation noise; LARGER explores longer\n"
         "    --prior-var F         arm prior spread; larger explores longer\n"
         "    --threads N           workers (0 = hardware concurrency)\n"
         "    --samples-per-reagent N  completions scored per reagent (greedy)\n"
         "    --trajectories N      independent search trajectories (k)\n"
         "    --trajectory-seed-stride N  seed gap between trajectories\n"
         "    --batch-size N        Thompson candidates per feedback round\n"
         "    --duplicate-retries N fresh TS proposals before mutation\n"
         "    --mutation-retries N  one-coordinate proposals before uniform\n"
         "    --pair-tail N         redundant experiment: joint tail draws\n"
         "    --pair-escape-draws N diagnostic two-coordinate moves\n"
         "    --pair-escape-batch N scoring batch size for pair escape\n"
         "    --pair-escape-rounds N improving jumps allowed (default 1)\n"
         "    --draws N             candidates drawn by evenpairs\n"
         "    --topk N              evenpairs aggregate: mean of a synthon's\n"
         "                          best N scores (1 = plain max)\n"
         "    --refine              coordinate-descent polish after sampling\n"
         "    --cross-lib FILE      a SECOND library (the other embed style) --\n"
         "                          reports how the planted product scores\n"
         "                          ACROSS fidelities: does a coarse product\n"
         "                          reproduce the shape a fine query defines?\n"

         "    --query-pose N        which conformer of the planted product is\n"
         "                          the query pose (default 0 = lowest energy)\n"
         "    --sweep-trace         per-step report of where the target\n"
         "                          synthons rank during a greedy sweep\n"
         "    --determinism-check   build the same molecule twice and compare:\n"
         "                          is confgen reproducible at fixed seed?\n"
         "    --failure-census      enumerate the whole space and report WHY\n"
         "                          products fail to build, and which synthons\n"
         "                          are implicated.  Small spaces only.\n"
         "    --exhaustive          ALSO score every product, to report where\n"
         "                          the planted one actually ranks.  Only\n"
         "                          sane for small spaces.\n"
         "    --legacy-ts-replay    replay duplicate TS draws (A/B only)\n"
         "    --seed N              RNG seed for the planted picks\n"
         "    --search-seed-offset N add only to search RNGs; targets unchanged\n"
         "    -h/--help\n";
}

//! Read one reaction out of a synthon table, preserving position order.
/*!
  Delimiter and column order come from the HEADER, so both layouts in the wild
  work unchanged:
      SMILES,synton_id,synton_role,reaction_id     (comma, "synton_2")
      SMILES\tsynton_id\tsynton#\treaction_id       (tab,   "2")
  The role column is whatever names a synthon position; roles are sorted as
  strings, which orders 1/2/3 and synton_1/synton_2/synton_3 alike.
*/
EnumerationTypes::BBS loadSynthonReaction(const std::string &path,
                                          const std::string &rxnId,
                                          size_t maxPerPosition) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "cannot open " << path << "\n";
    return {};
  }
  std::string line;
  if (!std::getline(in, line)) {
    return {};
  }
  const char delim = line.find('\t') != std::string::npos ? '\t' : ',';
  auto split = [delim](const std::string &l) {
    std::vector<std::string> f;
    std::stringstream ss(l);
    std::string cell;
    while (std::getline(ss, cell, delim)) {
      if (!cell.empty() && cell.back() == '\r') {
        cell.pop_back();
      }
      f.push_back(cell);
    }
    return f;
  };

  int cSmiles = -1, cRole = -1, cRxn = -1;
  const auto hdr = split(line);
  for (size_t i = 0; i < hdr.size(); ++i) {
    const std::string &h = hdr[i];
    if (h == "SMILES" || h == "smiles") {
      cSmiles = static_cast<int>(i);
    } else if (h.find("reaction") != std::string::npos) {
      cRxn = static_cast<int>(i);
    } else if (h.find("role") != std::string::npos ||
               h.find('#') != std::string::npos) {
      cRole = static_cast<int>(i);
    }
  }
  if (cSmiles < 0 || cRole < 0 || cRxn < 0) {
    std::cerr << "unrecognised header in " << path
              << " -- need a SMILES "
                 "column, a synthon role/# column and a reaction id column\n";
    return {};
  }

  std::map<std::string, MOL_SPTR_VECT> byRole;
  const int need = std::max(cSmiles, std::max(cRole, cRxn));
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    const auto f = split(line);
    if (static_cast<int>(f.size()) <= need || f[cRxn] != rxnId) {
      continue;
    }
    if (maxPerPosition && byRole[f[cRole]].size() >= maxPerPosition) {
      continue;
    }
    ROMOL_SPTR m(SmilesToMol(f[cSmiles]));
    if (m) {
      byRole[f[cRole]].push_back(m);
    }
  }
  EnumerationTypes::BBS bbs;
  for (auto &kv : byRole) {  // string order == position order
    bbs.push_back(kv.second);
  }
  return bbs;
}

}  // namespace

int main(int argc, char **argv) {
  std::string csv, rxn, loadFile, saveFile;
  size_t maxPerPos = 0;
  int productConfs = -1;  // -1 = leave the library default
  int trials = 5, budget = 8000, seed = 0xf00d;
  // Added to the SEARCH RNGs only, never to the planted-product draw, so a
  // RESTART can be measured against the identical set of targets.  The
  // measured complementarity between two searches came from them landing in
  // different basins, not from their mechanisms composing -- this isolates
  // that.
  int searchSeedOffset = 0;
  //! Independent search trajectories run SEQUENTIALLY in this one process,
  //! reduced at the end by keeping the best product.  Restarts share nothing
  //! -- no model, no cache, no coordination -- so running them in one process
  //! is equivalent to running them in k containers, and needs no orchestration.
  int restarts = 1;
  double noiseVar = -1.0, priorVar = -1.0;  // <0 = leave the default
  unsigned int threads = 0;
  bool prefill = true, refine = false, rejectDuplicates = true, coarse = false;
  bool exhaustive = false;
  bool sweepTrace = false;
  bool failureCensus = false;
  bool determinismCheck = false;
  std::vector<std::string> zipTriples;
  long redundancy = 0;
  long slice = -1;
  int topN = 10;
  std::string greedyAggregate = "max";
  bool fraglibStats = false;
  long decoyScores = 0;
  int pairTopK = 0;
  int pairTailDraws = 0;
  int pairEscapeDraws = 4000;
  unsigned int pairEscapeBatch = 256;
  unsigned int pairEscapeRounds = 1;
  unsigned int trajectories = 2;
  int trajectorySeedStride = 1000;
  std::string partners = "random";
  bool centralFirst = false;
  unsigned int queryPose = 0;
  std::string crossFile;
  unsigned int samplesPerReagent = 1;
  int evenDraws = 4000;
  int topK = 3;
  unsigned int batchSize = 64;
  unsigned int duplicateRetries = 16;
  unsigned int mutationRetries = 64;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto val = [&](const char *what) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "missing value for " << what << "\n";
        std::exit(1);
      }
      return argv[++i];
    };
    if (a == "--synthons" || a == "--csv") {
      csv = val("--synthons");
    } else if (a == "--rxn") {
      rxn = val("--rxn");
    } else if (a == "--load") {
      loadFile = val("--load");
    } else if (a == "--save") {
      saveFile = val("--save");
    } else if (a == "--max-per-pos") {
      maxPerPos = static_cast<size_t>(std::stoul(val("--max-per-pos")));
    } else if (a == "--product-confs") {
      productConfs = std::stoi(val("--product-confs"));
    } else if (a == "--topk") {
      topK = std::stoi(val("--topk"));
    } else if (a == "--draws") {
      evenDraws = std::stoi(val("--draws"));
    } else if (a == "--samples-per-reagent") {
      samplesPerReagent =
          static_cast<unsigned int>(std::stoul(val("--samples-per-reagent")));
    } else if (a == "--trajectories") {
      trajectories =
          static_cast<unsigned int>(std::stoul(val("--trajectories")));
    } else if (a == "--trajectory-seed-stride") {
      trajectorySeedStride = std::stoi(val("--trajectory-seed-stride"));
    } else if (a == "--batch-size") {
      batchSize = static_cast<unsigned int>(std::stoul(val("--batch-size")));
    } else if (a == "--duplicate-retries") {
      duplicateRetries =
          static_cast<unsigned int>(std::stoul(val("--duplicate-retries")));
    } else if (a == "--mutation-retries") {
      mutationRetries =
          static_cast<unsigned int>(std::stoul(val("--mutation-retries")));
    } else if (a == "--cross-lib") {
      crossFile = val("--cross-lib");
    } else if (a == "--library-query") {
      /* accepted and ignored: the query is always one pose of the built
         product now -- see the comment at the query construction */
    } else if (a == "--partners") {
      partners = val("--partners");
    } else if (a == "--pair-tail") {
      pairTailDraws = std::stoi(val("--pair-tail"));
    } else if (a == "--pair-escape-draws") {
      pairEscapeDraws = std::stoi(val("--pair-escape-draws"));
    } else if (a == "--pair-escape-batch") {
      pairEscapeBatch =
          static_cast<unsigned int>(std::stoul(val("--pair-escape-batch")));
    } else if (a == "--pair-escape-rounds") {
      pairEscapeRounds =
          static_cast<unsigned int>(std::stoul(val("--pair-escape-rounds")));
    } else if (a == "--central-first") {
      centralFirst = true;
    } else if (a == "--pair-refine") {
      pairTopK = std::stoi(val("--pair-refine"));
    } else if (a == "--decoy-scores") {
      decoyScores = std::stol(val("--decoy-scores"));
    } else if (a == "--fraglib-stats") {
      fraglibStats = true;
    } else if (a == "--greedy-aggregate") {
      greedyAggregate = val("--greedy-aggregate");
    } else if (a == "--top-n") {
      topN = std::stoi(val("--top-n"));
    } else if (a == "--slice") {
      slice = std::stol(val("--slice"));
    } else if (a == "--redundancy") {
      redundancy = std::stol(val("--redundancy"));
    } else if (a == "--zip") {
      zipTriples.push_back(val("--zip"));
    } else if (a == "--query-pose") {
      queryPose = static_cast<unsigned int>(std::stoul(val("--query-pose")));
    } else if (a == "--determinism-check") {
      determinismCheck = true;
    } else if (a == "--failure-census") {
      failureCensus = true;
    } else if (a == "--sweep-trace") {
      sweepTrace = true;
    } else if (a == "--exhaustive") {
      exhaustive = true;
    } else if (a == "--coarse") {
      coarse = true;
    } else if (a == "--no-prefill") {
      prefill = false;
    } else if (a == "--trials") {
      trials = std::stoi(val("--trials"));
    } else if (a == "--noise-var") {
      noiseVar = std::stod(val("--noise-var"));
    } else if (a == "--prior-var") {
      priorVar = std::stod(val("--prior-var"));
    } else if (a == "--budget") {
      budget = std::stoi(val("--budget"));
    } else if (a == "--threads") {
      threads = static_cast<unsigned int>(std::stoul(val("--threads")));
    } else if (a == "--refine") {
      refine = true;
    } else if (a == "--legacy-ts-replay" || a == "--no-cache") {
      // --no-cache is retained as a compatibility alias for old A/B scripts.
      rejectDuplicates = false;
    } else if (a == "--restarts") {
      restarts = std::stoi(val("--restarts"));
    } else if (a == "--search-seed-offset") {
      searchSeedOffset = std::stoi(val("--search-seed-offset"));
    } else if (a == "--seed") {
      seed = std::stoi(val("--seed"));
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else {
      std::cerr << "unknown argument: " << a << "\n";
      usage(argv[0]);
      return 1;
    }
  }
  if (loadFile.empty() && (csv.empty() || rxn.empty())) {
    usage(argv[0]);
    return 1;
  }

  RDLog::InitLogs();
  EnumerateSynthons3DParams params;
  params.prefillFraglib = prefill;
  params.storeFraglib = true;
  params.embedStyle =
      coarse ? SynthonEmbedStyle::Coarse : SynthonEmbedStyle::Full;
  if (productConfs > 0) {
    params.confgen.numOutputConfs = productConfs;
  }

  EnumerateSynthons3D lib;
  double buildMs = 0.0;
  if (!loadFile.empty()) {
    std::ifstream in(loadFile, std::ios::binary);
    if (!in) {
      std::cerr << "cannot open " << loadFile << "\n";
      return 1;
    }
    const double t0 = nowMs();
    try {
      lib.initFromStream(in);
    } catch (const std::exception &e) {
      std::cerr << "failed to load " << loadFile << ": " << e.what() << "\n";
      return 1;
    }
    buildMs = nowMs() - t0;
    std::printf("[synthonbench] loaded %s in %.0f ms\n", loadFile.c_str(),
                buildMs);
  } else {
    const auto bbs = loadSynthonReaction(csv, rxn, maxPerPos);
    if (bbs.empty()) {
      std::cerr << "no synthons for reaction '" << rxn << "' in " << csv
                << "\n";
      return 1;
    }
    const double t0 = nowMs();
    lib = EnumerateSynthons3D(bbs, params);
    buildMs = nowMs() - t0;
    std::printf("[synthonbench] built %s/%s in %.0f ms%s\n", csv.c_str(),
                rxn.c_str(), buildMs,
                prefill ? " (fragments pre-embedded)" : "");
    std::printf("[synthonbench] embed style: %s\n",
                coarse ? "COARSE (synthon junctions only)" : "FULL");
  }

  if (!lib.isValid()) {
    std::cerr << "library is not valid (exit-vector labels inconsistent?)\n";
    return 1;
  }

  // Report the size of the space we are searching.
  double space = 1.0;
  const size_t fraglibAtStart = lib.fraglib() ? lib.fraglib()->size() : 0;
  std::printf("[synthonbench] arity=%u  reagents:", lib.arity());
  for (unsigned int p = 0; p < lib.arity(); ++p) {
    std::printf(" %u", lib.numReagents(p));
    space *= lib.numReagents(p);
  }
  std::printf("  -> %.3g products\n", space);

  if (!saveFile.empty()) {
    std::ofstream out(saveFile, std::ios::binary);
    const double t0 = nowMs();
    lib.toStream(out);
    out.close();
    std::printf("[synthonbench] saved -> %s (%.0f ms)\n", saveFile.c_str(),
                nowMs() - t0);
  }

  // Optional second library, the other embed style, for cross-fidelity scoring.
  EnumerateSynthons3D crossLib;
  bool haveCross = false;
  if (!crossFile.empty()) {
    std::ifstream in(crossFile, std::ios::binary);
    if (!in) {
      std::cerr << "cannot open " << crossFile << "\n";
      return 1;
    }
    try {
      crossLib.initFromStream(in);
      haveCross = crossLib.isValid();
    } catch (const std::exception &e) {
      std::cerr << "failed to load " << crossFile << ": " << e.what() << "\n";
      return 1;
    }
    if (haveCross && crossLib.arity() != lib.arity()) {
      std::cerr << "cross library arity does not match\n";
      return 1;
    }
    std::printf("[synthonbench] cross library: %s (%s)\n", crossFile.c_str(),
                crossLib.params3D().embedStyle == SynthonEmbedStyle::Coarse
                    ? "COARSE"
                    : "FULL");
  }

  // ALWAYS emit the complete effective experiment configuration.  A/B output
  // without this block is not a reproducible measurement.
  std::cout
      << "[params.begin]\n"
      << "search=multipleTrajectory\n"
      << "library.source=" << (loadFile.empty() ? csv : loadFile) << "\n"
      << "library.reaction=" << rxn << "\n"
      << "library.maxPerPosition=" << maxPerPos << "\n"
      << "library.embedStyle="
      << (lib.params3D().embedStyle == SynthonEmbedStyle::Coarse ? "Coarse"
                                                                 : "Full")
      << "\n"
      << "library.prefill=" << lib.params3D().prefillFraglib << "\n"
      << "trials=" << trials << "\n"
      << "seed=" << seed << "\n"
      << "trialSeed=seed+trialIndex\n"
      << "searchSeedOffset=" << searchSeedOffset << "\n"
      << "searchSeed=seed+trialIndex+searchSeedOffset\n"
      << "threads=" << threads << "\n"
      << "queryPose=" << queryPose << "\n"
      << "refine=" << refine << "\n"
      << "topN=" << topN << "\n"
      << "greedy.samplesPerReagent=" << samplesPerReagent << "\n"
      << "greedy.centralFirst=" << centralFirst << "\n"
      << "greedy.pairRefineTopK=" << pairTopK << "\n"
      << "greedy.internalCoordinateRefine=true\n"
      << "multipleTrajectories.count=" << trajectories << "\n"
      << "multipleTrajectories.seedStride=" << trajectorySeedStride << "\n"
      << "multipleTrajectories.seed="
         "seed+trialIndex+searchSeedOffset+trajectory*seedStride\n"
      << "multipleTrajectories.workerPool=shared\n"
      << "multipleTrajectories.tupleCache=per-query\n"
      << "multipleTrajectories.cacheSemantics=trajectory-preserving\n"
      << "multipleTrajectories.finalists=rebuilt-and-rescored\n"
      << "timing.wallS=search-only; excludes library load, prefill, and query setup\n"
      << "timing.componentMs=sum across worker threads, not wall time\n"
      << "library.fraglibEntriesAtStart="
      << (lib.fraglib() ? lib.fraglib()->size() : 0) << "\n"
      << "conformer.params.begin\n"
      << fragmentConfGenParamsToString(lib.params3D().confgen)
      << "conformer.params.end\n"
      << "[params.end]\n";

  // How many conformers does each SYNTHON actually get?  The per-class
  // embedding recipes assume FINE fragmentation, where an acyclic fragment has
  // no internal rotors and one conformer is right.  A COARSE library caches
  // whole synthons, which keep their internal rotors -- so if they land in the
  // same class they are frozen to a single internal geometry.
  if (fraglibStats) {
    const auto &fl = lib.fraglib();
    if (!fl) {
      std::printf("[fraglib] library has no stored fragment library\n");
      return 1;
    }
    std::map<unsigned int, unsigned int> hist, histRing, histAcyclic;
    unsigned int missing = 0, dead = 0, total = 0;
    for (const auto &position : lib.getReagents()) {
      for (const auto &synthon : position) {
        if (!synthon) {
          continue;
        }
        ++total;
        RWMol whole(*synthon);
        for (auto atom : whole.atoms()) {
          if (atom->getAtomicNum() == 92 || atom->getAtomicNum() == 93) {
            atom->setAtomicNum(0);
            atom->setIsotope(0);
            atom->setNoImplicit(true);
            atom->setNumExplicitHs(0);
          }
        }
        try {
          MolOps::sanitizeMol(whole);
        } catch (...) {
          continue;
        }
        const auto nConfs = fl->numFragmentConfs(whole);
        if (!nConfs) {
          ++missing;  // never attempted
          continue;
        }
        if (!*nConfs) {
          ++dead;  // tombstoned: attempted and unembeddable
          continue;
        }
        const unsigned int n = *nConfs;
        ++hist[n];
        if (whole.getRingInfo()->numRings()) {
          ++histRing[n];
        } else {
          ++histAcyclic[n];
        }
      }
    }
    auto dump = [](const char *what,
                   const std::map<unsigned int, unsigned int> &h) {
      unsigned int tot = 0, confs = 0;
      for (const auto &kv : h) {
        tot += kv.second;
        confs += kv.first * kv.second;
      }
      std::printf("[fraglib] %-9s %4u synthons, mean %.2f confs:", what, tot,
                  tot ? static_cast<double>(confs) / tot : 0.0);
      for (const auto &kv : h) {
        std::printf("  %ux%u", kv.second, kv.first);
      }
      std::printf("\n");
    };
    std::printf("[fraglib] %u synthons, %u not in cache, %u tombstoned\n",
                total, missing, dead);
    dump("all", hist);
    dump("ring", histRing);
    dump("acyclic", histAcyclic);
    return 0;
  }

  // Exact route multiplicity over a full slice: fix position 0 and enumerate
  // every remaining combination.  Unlike the random sample this counts routes
  // per product directly, with no birthday extrapolation.
  if (slice >= 0 && lib.arity() >= 2) {
    std::map<std::string, unsigned int> routes;
    long built = 0, total = 0;
    std::vector<unsigned int> idx(lib.arity(), 0);
    idx[0] = static_cast<unsigned int>(slice);
    std::function<void(unsigned int)> rec = [&](unsigned int pos) {
      if (pos == lib.arity()) {
        ++total;
        ROMOL_SPTR m = lib.get2D(idx);
        if (m) {
          ++built;
          ++routes[MolToSmiles(*m)];
        }
        return;
      }
      for (unsigned int r = 0; r < lib.numReagents(pos); ++r) {
        idx[pos] = r;
        rec(pos + 1);
      }
    };
    rec(1);
    std::map<unsigned int, unsigned int> hist;
    for (const auto &kv : routes) {
      ++hist[kv.second];
    }
    std::printf("[slice %ld] %ld combinations -> %ld built, %zu distinct\n",
                slice, total, built, routes.size());
    for (const auto &kv : hist) {
      std::printf("           %u route(s): %u products (%.2f%%)\n", kv.first,
                  kv.second,
                  100.0 * kv.second / static_cast<double>(routes.size()));
    }
    std::printf("[slice %ld] mean routes/product = %.3f\n", slice,
                routes.empty() ? 0.0
                               : static_cast<double>(built) /
                                     static_cast<double>(routes.size()));
    return 0;
  }

  // How much of the nominal product space is actually DISTINCT?  A synthon
  // library reaches the same molecule from different reagent triples, so the
  // true space is smaller than the index product -- and any benchmark that
  // scores recovery by reagent indices under-counts by roughly this rate.
  if (redundancy > 0) {
    auto join = [](const std::vector<unsigned int> &v) {
      std::string r;
      for (size_t i = 0; i < v.size(); ++i) {
        r += (i ? "," : "") + std::to_string(v[i]);
      }
      return r;
    };
    std::mt19937 rng(1234u);
    std::map<std::string, std::vector<unsigned int>> byProduct;
    long built = 0;
    for (long i = 0; i < redundancy; ++i) {
      std::vector<unsigned int> idx(lib.arity());
      for (unsigned int p = 0; p < lib.arity(); ++p) {
        idx[p] = std::uniform_int_distribution<unsigned int>(
            0, lib.numReagents(p) - 1)(rng);
      }
      ROMOL_SPTR m = lib.get2D(idx);
      if (!m) {
        continue;
      }
      ++built;
      const std::string smi = MolToSmiles(*m);
      auto it = byProduct.find(smi);
      if (it == byProduct.end()) {
        byProduct.emplace(smi, idx);
      } else if (it->second != idx) {
        std::printf("[redundancy] %s\n  also from %s (first %s)\n", smi.c_str(),
                    join(idx).c_str(), join(it->second).c_str());
      }
    }
    std::printf(
        "[redundancy] %ld random triples -> %ld built, %zu distinct products "
        "(%.2f%% collide)\n",
        redundancy, built, byProduct.size(),
        built ? 100.0 * (built - static_cast<long>(byProduct.size())) / built
              : 0.0);
    return 0;
  }

  // Zip named reagent triples and print their canonical SMILES.  Distinct
  // indices can build the SAME product, which is the difference between a
  // real ranking failure and a bookkeeping artifact.
  if (!zipTriples.empty()) {
    std::string first;
    for (const auto &spec : zipTriples) {
      std::vector<unsigned int> idx;
      std::stringstream ss(spec);
      std::string cell;
      while (std::getline(ss, cell, ',')) {
        idx.push_back(static_cast<unsigned int>(std::stoul(cell)));
      }
      if (idx.size() != lib.arity()) {
        std::printf("[zip] %s: need %u indices, got %zu\n", spec.c_str(),
                    lib.arity(), idx.size());
        continue;
      }
      ROMOL_SPTR m = lib.get2D(idx);
      const std::string smi =
          m ? MolToSmiles(*m) : std::string("<unbuildable>");
      if (first.empty()) {
        first = smi;
      }
      std::printf("[zip] %-22s %s%s\n", spec.c_str(), smi.c_str(),
                  (!first.empty() && smi == first && &spec != &zipTriples[0])
                      ? "   <-- SAME as first"
                      : "");
    }
    return 0;
  }

  if (determinismCheck) {
    // Is FragmentConfGen reproducible?  Build the SAME molecule twice with two
    // fresh generators, same params, same seed, both with cold private
    // fraglibs -- then again through the library's warm shared cache.  If any
    // of these disagree, conformer generation depends on cache state or
    // ordering rather than on (molecule, seed), and every warm-vs-cold
    // comparison we have run measures different OUTPUTS, not just speed.
    std::mt19937 drng(1234);
    for (int t = 0; t < 3; ++t) {
      std::vector<unsigned int> idx(lib.arity(), 0);
      for (unsigned int p = 0; p < lib.arity(); ++p) {
        idx[p] = std::uniform_int_distribution<unsigned int>(
            0, lib.numReagents(p) - 1)(drng);
      }
      ROMOL_SPTR graph = lib.get2D(idx);
      if (!graph) {
        continue;
      }
      auto buildOnce = [&](void) -> ROMOL_SPTR {
        FragmentConfGenParams qp;
        qp.numOutputConfs = 10;
        FragmentConfGen g(qp);
        FragmentConfGenResult r = g.build(*graph);
        if (r.conformers.empty()) {
          return {};
        }
        auto m = boost::make_shared<RWMol>(*r.conformers.front());
        for (size_t i = 1; i < r.conformers.size(); ++i) {
          if (r.conformers[i]->getNumConformers()) {
            m->addConformer(new Conformer(r.conformers[i]->getConformer()),
                            true);
          }
        }
        return m;
      };
      ROMOL_SPTR a = buildOnce();
      ROMOL_SPTR b = buildOnce();
      SynthonProduct viaLib = lib.getProduct(idx);
      if (!a || !b) {
        continue;
      }
      unsigned int nrot = 0;
      {
        RWMol tmp(*a);
        nrot = Descriptors::calcNumRotatableBonds(tmp);
      }
      // compare A's FIRST conformer as the query against each set
      ShapeScorer sc(*a, a->getConformer(0).getId());
      const auto sa = sc.score(*a);
      const auto sb = sc.score(*b);
      // Best RMSD between the two ensembles: for each conformer of A, the
      // closest conformer of B.  Symmetry-aware, heavy atoms only.  Shape
      // overlay conflates alignment quality with how different the geometries
      // are; this says it in Angstrom.
      // NB conformer IDs are not indices -- collect the real IDs first.
      auto confIds = [](const ROMol &m) {
        std::vector<int> ids;
        for (auto c = m.beginConformers(); c != m.endConformers(); ++c) {
          ids.push_back((*c)->getId());
        }
        return ids;
      };
      double bestPair = 1e30, firstConfBest = 1e30;
      int bestIdA = -1, bestIdB = -1;
      const std::vector<int> aIds = confIds(*a), bIds = confIds(*b);
      {
        ROMOL_SPTR ah(MolOps::removeHs(*a));
        ROMOL_SPTR bh(MolOps::removeHs(*b));
        for (size_t i = 0; i < aIds.size(); ++i) {
          for (size_t j = 0; j < bIds.size(); ++j) {
            RWMol probe(*ah), ref(*bh);
            double r = 1e30;
            try {
              r = MolAlign::getBestRMS(probe, ref, aIds[i], bIds[j]);
            } catch (...) {
              continue;
            }
            if (r < bestPair) {
              bestPair = r;
              bestIdA = aIds[i];
              bestIdB = bIds[j];
            }
            if (i == 0) {
              firstConfBest = std::min(firstConfBest, r);
            }
          }
        }
      }
      // Score the CLOSEST pair specifically: query = A conf bestI, candidate =
      // B conf bestJ alone.  If a 0.5A pair does not score near 1, the scorer
      // is over-sensitive; if it does, the earlier low numbers were simply the
      // query conformer having no close partner (A conf 0 was 1.74A away).
      double closestPairScore = 0.0;
      if (bestIdA >= 0 && bestIdB >= 0) {
        ShapeScorer pairSc(*a, bestIdA);
        auto bOne = boost::make_shared<RWMol>(*b);
        std::vector<int> drop;
        for (auto ci = bOne->beginConformers(); ci != bOne->endConformers();
             ++ci) {
          if ((int)(*ci)->getId() != bestIdB) {
            drop.push_back((*ci)->getId());
          }
        }
        for (int id : drop) {
          bOne->removeConformer(id);
        }
        const auto ps = pairSc.score(*bOne);
        closestPairScore = ps ? *ps : 0.0;
      }
      std::printf(
          "[determinism] rotors=%u  A confs=%u B confs=%u  "
          "A-vs-A=%.4f  A-vs-B(fresh, same seed)=%.4f  "
          "bestRMSD(A,B)=%.2fA (confIds %d/%d -> combo %.4f confA0->nearestB=%.2fA",
          nrot, a->getNumConformers(), b->getNumConformers(), sa ? *sa : 0.0,
          sb ? *sb : 0.0, bestPair, bestIdA, bestIdB, closestPairScore,
          firstConfBest);
      if (viaLib) {
        const auto sl = sc.score(*viaLib.mol);
        std::printf("  A-vs-LIBRARY(warm cache)=%.4f", sl ? *sl : 0.0);
      }
      std::printf("\n");
      std::fflush(stdout);
    }
  }

  if (failureCensus) {
    // Walk the whole space, tally failures by status, and count how often each
    // reagent takes part in one.  A synthon that cannot be embedded poisons
    // EVERY product containing it, so a per-reagent count separates "this one
    // building block is broken" from "these particular combinations are".
    std::map<int, size_t> byStatus;
    std::vector<std::map<unsigned int, size_t>> byReagent(lib.arity());
    std::vector<std::map<unsigned int, size_t>> totalReagent(lib.arity());
    std::vector<unsigned int> idx(lib.arity(), 0);
    size_t total = 0, failed = 0;
    for (;;) {
      const SynthonProduct p = lib.getProduct(idx);
      ++total;
      for (unsigned int q = 0; q < lib.arity(); ++q) {
        totalReagent[q][idx[q]]++;
      }
      if (!p) {
        ++failed;
        byStatus[static_cast<int>(p.status)]++;
        for (unsigned int q = 0; q < lib.arity(); ++q) {
          byReagent[q][idx[q]]++;
        }
      }
      size_t pos = 0;
      for (; pos < idx.size(); ++pos) {
        if (++idx[pos] < lib.numReagents(static_cast<unsigned int>(pos))) {
          break;
        }
        idx[pos] = 0;
      }
      if (pos == idx.size()) {
        break;
      }
    }
    std::printf("\n[failure census] %zu products, %zu failed (%.1f%%)\n", total,
                failed, 100.0 * failed / total);
    for (const auto &kv : byStatus) {
      std::printf(
          "  status %-28s %zu\n",
          synthonBuildStatusMessage(static_cast<SynthonBuildStatus>(kv.first)),
          kv.second);
    }
    for (unsigned int q = 0; q < lib.arity(); ++q) {
      std::printf("  position %u: reagents implicated in a failure:\n", q);
      for (const auto &kv : byReagent[q]) {
        const size_t tot = totalReagent[q][kv.first];
        if (kv.second * 4 >= tot) {  // fails in >=25% of its products
          std::printf("    reagent %-5u fails %zu/%zu (%.0f%%)\n", kv.first,
                      kv.second, tot, 100.0 * kv.second / tot);
        }
      }
    }
    std::fflush(stdout);
  }

  std::mt19937 rng(static_cast<unsigned int>(seed));
  unsigned int recovered = 0, tied = 0, missed = 0, unbuildable = 0;
  unsigned int inTopN = 0;

  std::printf("\n%-6s %-22s %-22s %8s %8s %9s %9s %9s %9s %9s %7s\n", "trial",
              "planted", "found", "score", "target", "draws", "cached",
              "rejected", "confCpuMs", "shapeCpu", "wallS");

  auto join = [](const std::vector<unsigned int> &v) {
    std::string s;
    for (size_t i = 0; i < v.size(); ++i) {
      s += (i ? "," : "") + std::to_string(v[i]);
    }
    return s;
  };

  for (int t = 0; t < trials; ++t) {
    const size_t cacheAtTrialStart = lib.fraglib() ? lib.fraglib()->size() : 0;
    // 1. plant a product at random
    std::vector<unsigned int> planted(lib.arity());
    for (unsigned int p = 0; p < lib.arity(); ++p) {
      planted[p] = std::uniform_int_distribution<unsigned int>(
          0, lib.numReagents(p) - 1)(rng);
    }
    // THE QUERY IS ONE POSE the planted product actually has.
    //
    // It comes from the SAME generation as the candidates and is then pruned to
    // a single conformer.  Both halves matter:
    //  - one pose, not the ensemble: using all 10 conformers let ANY of ten
    //    geometries match, which made the target far too easy to find;
    //  - from the same generation: a separately-built query (its own cold
    //    fraglib, fragments embedded in a different order -- ETKDG is
    //    order-sensitive) produces conformers that are in NOBODY's ensemble, so
    //    nothing can match it and the correct answer scored 0.26-0.49 while
    //    unrelated products scored higher.  That measured conformer
    //    reproducibility, not search quality.
    // With one achievable pose the planted product scores 1.0 and nothing else
    // can, so recovery means what it says.
    SynthonProduct target = lib.getProduct(planted);
    int queryConfId = -1;
    if (target && target.mol->getNumConformers()) {
      const unsigned int n = target.mol->getNumConformers();
      const unsigned int pick = std::min<unsigned int>(queryPose, n - 1);
      queryConfId = target.mol->getConformer(pick).getId();
    }
    if (!target) {
      ++unbuildable;
      std::printf("%-6d (planted product could not be built: %s)\n", t,
                  synthonBuildStatusMessage(target.status));
      continue;
    }

    // 2. its conformers ARE the query
    ShapeScorer scorer(*target.mol, queryConfId);
    const auto selfScore = scorer.score(*target.mol);

    // Noise floor: how well do RANDOM products score against this query?  More
    // conformers per synthon gives every candidate more chances to fit, which
    // raises this floor and compresses the contrast the greedy marginal needs
    // -- a discrimination cost that does not show up in the planted product's
    // own score (it stays 1.0).
    if (decoyScores > 0) {
      std::mt19937 drng(static_cast<unsigned int>(seed + t));
      std::vector<double> ds;
      for (long i = 0; i < decoyScores; ++i) {
        std::vector<unsigned int> idx(lib.arity());
        for (unsigned int q = 0; q < lib.arity(); ++q) {
          idx[q] = std::uniform_int_distribution<unsigned int>(
              0, lib.numReagents(q) - 1)(drng);
        }
        if (idx == planted) {
          continue;
        }
        SynthonProduct sp = lib.getProduct(idx);
        if (!sp) {
          continue;
        }
        if (const auto v = scorer.score(*sp.mol)) {
          ds.push_back(*v);
        }
      }
      if (!ds.empty()) {
        std::sort(ds.begin(), ds.end());
        const double mean =
            std::accumulate(ds.begin(), ds.end(), 0.0) / ds.size();
        std::printf(
            "       decoys n=%zu  mean %.4f  median %.4f  p90 %.4f  max %.4f\n",
            ds.size(), mean, ds[ds.size() / 2], ds[(ds.size() * 9) / 10],
            ds.back());
      }
    }
    // The ACHIEVABLE ceiling is the planted product as the LIBRARY builds it,
    // scored against this query -- not the query against itself.  With an
    // independent query nothing can reproduce the query's own conformers, so
    // the self-overlay 1.0 is unreachable and comparing to it is meaningless.
    {
      SynthonProduct built = lib.getProduct(planted);
      if (built) {
        const auto bs = scorer.score(*built.mol);
        unsigned int nrot = 0;
        {
          RWMol tmp(*built.mol);
          nrot = Descriptors::calcNumRotatableBonds(tmp);
        }
        std::printf(
            "       planted: rotors=%u  query=1 pose (of %u)  BUILT confs=%u  "
            "as-built scores %.4f\n",
            nrot, target.mol->getNumConformers(), built.mol->getNumConformers(),
            bs ? *bs : 0.0);
      }
    }
    // A self-overlay is the ceiling: shape 1.0 AND colour 1.0, i.e. combo 2.0.
    // Anything less means one of the two terms is not contributing.
    if (t == 0) {
      std::printf(
          "[synthonbench] self-overlay: combo=%.4f\n",
          selfScore ? *selfScore : 0.0);
    }

    // How does the SAME product, built at the other fidelity, score against
    // this query?  If a coarse product cannot reproduce the shape a fine query
    // defines, then a coarse first pass is filtering on a different objective
    // than the fine pass will judge -- which is the whole premise of staging.
    if (haveCross) {
      SynthonProduct other = crossLib.getProduct(planted);
      if (other) {
        const auto xs = scorer.score(*other.mol);
        std::printf(
            "       planted: this-style confs=%u self combo=%.4f "
            "other-style confs=%u cross combo=%.4f "
            "\n",
            target.mol->getNumConformers(), selfScore ? *selfScore : 0.0,
            other.mol->getNumConformers(),
            xs ? *xs : 0.0);
      } else {
        std::printf("       planted: other style could not build it (%s)\n",
                    synthonBuildStatusMessage(other.status));
      }
    }

    // 3. search for it
    const size_t cacheAtSearchStart = lib.fraglib() ? lib.fraglib()->size() : 0;
    const double w0 = nowMs();

    // ---- the search.  One entry point: k independent greedy+refine
    // trajectories through a shared scorer, reduced to the best products.
    // Every alternative that was measured against a budget-matched control --
    // Thompson sampling, even-pair sweeps, a beam, partner panels, pair-tail
    // and pair-escape -- lost or was equalled by simply running more
    // trajectories, and has been retired.  See SynthonSearch/restarts.md.
    SynthonSearch3DParams mp;
    mp.numTrajectories = std::max(1u, trajectories);
    mp.randomSeed = seed + t + searchSeedOffset;
    mp.seedStride = trajectorySeedStride;
    mp.samplesPerReagent = samplesPerReagent;
    mp.refineIters = 3;
    mp.pairRefineTopK = pairTopK;
    mp.numThreads = threads;
    mp.numBestProducts = topN;
    MultipleTrajectoryStats multipleStats;
    SynthonSearchResult res =
        synthonSearch3D(lib, scorer, mp, &multipleStats);
    const double wallS = (nowMs() - w0) / 1000.0;
    const size_t cacheAfterSearch = lib.fraglib() ? lib.fraglib()->size() : 0;

    // Recovery is a question about the MOLECULE, not about which synthons
    // were used to build it.  A synthon library is redundant: the same product
    // is reachable from different reagent triples (a linker atom assigned to
    // one position or the next), so index equality under-counts recovery and
    // reports a correct answer as a miss.  Compare canonical SMILES.
    bool exact = (res.reagents == planted);
    bool viaOtherRoute = false;
    if (!exact && !res.reagents.empty()) {
      ROMOL_SPTR pm = lib.get2D(planted);
      ROMOL_SPTR fm = lib.get2D(res.reagents);
      if (pm && fm && MolToSmiles(*pm) == MolToSmiles(*fm)) {
        exact = true;
        viaOtherRoute = true;
      }
    }
    // A different product that overlays at least as well is NOT a failure of
    // the search -- shape is the objective, the planted molecule is not.
    const bool asGood = !exact && selfScore && res.score >= *selfScore - 1e-6;
    if (exact) {
      ++recovered;
    } else if (asGood) {
      ++tied;
    } else {
      ++missed;
    }

    // Where does the planted product ACTUALLY rank?  "miss" conflates two very
    // different failures: never drawn (an exploration problem) and drawn but
    // out-ranked (a scoring problem, or a planted product that simply is not
    // the best overlay of itself).
    if (exhaustive) {
      std::vector<unsigned int> idx(lib.arity(), 0);
      size_t better = 0, scored = 0, unbuilt = 0;
      double plantedScore = selfScore ? *selfScore : 0.0;
      double bestSeen = -1e30;
      std::vector<unsigned int> bestIdx;
      for (;;) {
        SynthonProduct p = lib.getProduct(idx);
        if (p) {
          const auto sc = scorer.score(*p.mol);
          if (sc) {
            ++scored;
            if (*sc > plantedScore + 1e-9) {
              ++better;
            }
            if (*sc > bestSeen) {
              bestSeen = *sc;
              bestIdx = idx;
            }
          } else {
            ++unbuilt;
          }
        } else {
          ++unbuilt;
        }
        size_t pos = 0;  // odometer over the reagent indices
        for (; pos < idx.size(); ++pos) {
          if (++idx[pos] < lib.numReagents(static_cast<unsigned int>(pos))) {
            break;
          }
          idx[pos] = 0;
        }
        if (pos == idx.size()) {
          break;
        }
      }
      std::printf(
          "       exhaustive: planted scores %.4f, rank %zu of %zu scored "
          "(%zu unbuildable); best is %s at %.4f\n",
          plantedScore, better + 1, scored, unbuilt, join(bestIdx).c_str(),
          bestSeen);
    }

    std::printf(
        "%-6d %-22s %-22s %8.4f %8.4f %9u %9u %9u %9.0f %9.0f %7.2f%s\n", t,
        join(planted).c_str(), join(res.reagents).c_str(), res.score,
        selfScore ? *selfScore : 0.0, res.evaluations, res.cacheHits,
        res.duplicatesRejected, res.confgenMs, res.scoreMs, wallS,
        exact ? (viaOtherRoute ? "  HIT*" : "  HIT")
              : (asGood ? "  tie" : "  miss"));
    std::printf("       fragment-cache growth: setup +%zu, search +%zu\n",
                cacheAtSearchStart - cacheAtTrialStart,
                cacheAfterSearch - cacheAtSearchStart);
    std::printf(
        "       search accounting: scored=%u unscorable=%u attempted=%u\n",
        res.evaluations, res.unscorable, res.evaluations + res.unscorable);
    
    
    // The deliverable of a shape search is the TOP HITS, and the planted
    // product is a PROXY for "did we surface good overlays" -- not the
    // definition of correct.  A different molecule that overlays as well is a
    // legitimate hit (rare, but shape+colour 1.0 does not imply identity), so
    // report where the planted product LANDS in the returned list rather than
    // only whether it came first.
    if (!res.best.empty()) {
      ROMOL_SPTR pm = lib.get2D(planted);
      const std::string ps = pm ? MolToSmiles(*pm) : std::string();
      size_t rank = 0;
      double plantedInSearch = 0.0;
      for (size_t i = 0; i < res.best.size(); ++i) {
        ROMOL_SPTR hm = lib.get2D(res.best[i].reagents);
        if (hm && !ps.empty() && MolToSmiles(*hm) == ps) {
          rank = i + 1;
          plantedInSearch = res.best[i].score;
          break;
        }
      }
      if (rank) {
        ++inTopN;
      }
      // `as-built` is what the planted product scores when built standalone;
      // `in-search` is what it scored when the search actually evaluated it.
      // They differ if the rebuild does not reproduce the query pose.
      std::printf(
          "       planted ranks %s of %zu returned  (top %.4f, planted "
          "as-built %.4f, in-search %s)\n",
          rank ? std::to_string(rank).c_str() : "NOT IN TOP-N", res.best.size(),
          res.best.front().score, selfScore ? *selfScore : 0.0,
          rank ? (std::to_string(plantedInSearch)).c_str() : "never scored");
    }

    // A non-exact winner is only interesting if it is a DIFFERENT molecule.
    // Distinct reagent indices can still zip to the same product, and that
    // shows up as a spurious "tie" at combo 1.0 -- print both so the two
    // cases are never confused.
    if (!exact && !res.reagents.empty()) {
      ROMOL_SPTR pm = lib.get2D(planted);
      ROMOL_SPTR fm = lib.get2D(res.reagents);
      if (pm && fm) {
        std::printf("       planted: %s\n         found: %s\n",
                    MolToSmiles(*pm).c_str(), MolToSmiles(*fm).c_str());
      }
    }
    std::fflush(stdout);
  }

  // Did anything embed ON DEMAND during the search?  A prefilled coarse
  // library should already hold every fragment the assembler asks for, since
  // the junction bonds are always cut and each fragment is therefore a whole
  // synthon.  Growth here means the prefill is incomplete and production would
  // pay that cost on every cold query.
  std::printf(
      "[synthonbench] fraglib entries: %zu at start -> %zu at end (%+d "
      "embedded on demand)\n",
      fraglibAtStart, lib.fraglib() ? lib.fraglib()->size() : 0,
      static_cast<int>((lib.fraglib() ? lib.fraglib()->size() : 0) -
                       fraglibAtStart));

  const int scored = trials - static_cast<int>(unbuildable);
  std::printf(
      "\n[synthonbench] %d/%d recovered (HIT* = same product, other route), "
      "%d tied (>= planted score), %d missed%s\n"
      "[synthonbench] planted product in the returned top-N: %d/%d\n",
      recovered, scored, tied, missed,
      unbuildable
          ? (", " + std::to_string(unbuildable) + " unbuildable").c_str()
          : "",
      inTopN, scored);
  return 0;
}
