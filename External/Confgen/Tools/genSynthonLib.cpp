//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  genSynthonLib -- build a 3D synthon library fragments ONCE and save it.
//
#include <Confgen/SynthonSearch/EnumerateSynthons3D.h>
#include <Confgen/SynthonSearch/SynthonZipUtils.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
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
  std::cerr << "usage: " << argv0
            << " --synthons FILE --rxn ID --output FILE [options]\n\n"
               "  -i/--synthons FILE   synthon table, CSV or TSV; columns are\n"
               "                       found from the header (SMILES, synthon\n"
               "                       role/#, reaction id)\n"
               "  -r/--rxn ID          reaction to build\n"
               "  -o/--output FILE     output library path\n"
               "      --acyclic-confs N  conformers kept per ACYCLIC fragment.  [default 1]\n"
      "      --frag-rmsd F    RMSD radius (A) for the fragment pool [Default 0.1]\n"
      "      --style S        coarse | full (default coarse).  coarse embeds the\n"
               "                       entire synthon to reduce the number of rotatable bonds\n"
               "                       falling back to full per product when coarse fails\n"
               "                       n.b. coarse leaves a synthon on the ACYCLIC fragment class\n"
               "                       of 1 conformer -- see --acyclic-confs\n"
               "      --max-per-pos N  limit reagents per position (0 = all)\n"
               "      --product-confs N  conformers per assembled product\n"
               "      --no-prefill     skip embedding entirely\n"
               "      --threads N      prefill workers (0 = all cores)\n"
               "      --fraglib FILE   fragment file to create or append to\n"
               "                       Note: keep the same embedder for each reaction.\n"
               "  -h/--help\n";
}

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
    std::cerr << "unrecognised header in " << path << "\n";
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
  for (auto &kv : byRole) {
    bbs.push_back(kv.second);
  }
  return bbs;
}

}  // namespace

int main(int argc, char **argv) {
  std::string synthons, rxn, output, style = "coarse";
  size_t maxPerPos = 0;
  int acyclicConfs = -1;   // -1 = leave the class recipe alone
  double fragRms = -1.0;   // -1 = leave the prune radius alone
  int productConfs = -1;
  bool prefill = true;
  unsigned int threads = 1;
  std::string embedderPath;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto val = [&](const char *what) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "missing value for " << what << "\n";
        std::exit(1);
      }
      return argv[++i];
    };
    if (a == "-i" || a == "--synthons") {
      synthons = val("--synthons");
    } else if (a == "-r" || a == "--rxn") {
      rxn = val("--rxn");
    } else if (a == "-o" || a == "--output") {
      output = val("--output");
    } else if (a == "--acyclic-confs") {
      acyclicConfs = std::stoi(val("--acyclic-confs"));
    } else if (a == "--frag-rmsd") {
      fragRms = std::stod(val("--frag-rmsd"));
    } else if (a == "--style") {
      style = val("--style");
      SynthonEmbedStyle parsed;
      if (synthonEmbedStyleFromName(style, parsed) &&
          parsed == SynthonEmbedStyle::CoarseSampled) {
        std::cerr << "--style coarsesampled: the uncapped full-molecule "
                     "sampling that builds the fraglib is NOT IMPLEMENTED "
                     "yet.";

        return 1;
      }
      if (!synthonEmbedStyleFromName(style, parsed)) {
        std::cerr << "--style must be full, coarse or coarsesampled\n";
        return 1;
      }
    } else if (a == "--max-per-pos") {
      maxPerPos = static_cast<size_t>(std::stoul(val("--max-per-pos")));
    } else if (a == "--product-confs") {
      productConfs = std::stoi(val("--product-confs"));
    } else if (a == "--threads") {
      threads = static_cast<unsigned int>(std::stoul(val("--threads")));
    } else if (a == "--fraglib") {
      embedderPath = val("--fraglib");
    } else if (a == "--no-prefill") {
      prefill = false;
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else {
      std::cerr << "unknown argument: " << a << "\n";
      usage(argv[0]);
      return 1;
    }
  }
  if (synthons.empty() || rxn.empty() || output.empty()) {
    usage(argv[0]);
    return 1;
  }

  RDLog::InitLogs();
  const auto bbs = loadSynthonReaction(synthons, rxn, maxPerPos);
  if (bbs.empty()) {
    std::cerr << "no synthons for reaction '" << rxn << "' in " << synthons
              << "\n";
    return 1;
  }

  //  Hard fought fix, if the synthon reactions make ring closures, we
  //   either fail embeddings or make terrible conformers.
  //  Detect this so that the ring forms it's own fragment during embedding
  //   and anything hanging off we can treat embed seperately.
  std::vector<unsigned int> exitCounts;
  for (const auto &position : bbs) {
    unsigned int most = 0;
    for (const auto &syn : position) {
      if (!syn) {
        continue;
      }
      unsigned int n = 0;
      for (const auto atom : syn->atoms()) {
        if (atom->getAtomicNum() == 0 || atom->getAtomicNum() == 92 ||
            atom->getAtomicNum() == 93 || atom->getAtomicNum() == 94) {
          ++n;
        }
      }
      most = std::max(most, n);
    }
    exitCounts.push_back(most);
  }
  const unsigned int rings = SynthonZip::numRingClosures(exitCounts);

  EnumerateSynthons3DParams params;
  if (!synthonEmbedStyleFromName(style, params.embedStyle)) {
    std::cerr << "unknown --style '" << style
              << "'; expected full, coarse or coarsesampled\n";
    return 1;
  }
  
  // Do the ring detection and set the correct synthon cut bonds
  //  so we don't have to do Full conf gen whenever these fragments
  //  are used.
  auto ringAttachmentBonds = [](const ROMol &syn) {
    std::vector<unsigned int> cuts;
    std::vector<unsigned int> exits;
    for (const auto atom : syn.atoms()) {
      const int z = atom->getAtomicNum();
      if (z == 0 || z == 92 || z == 93 || z == 94) {
        exits.push_back(atom->getIdx());
      }
    }
    if (exits.size() < 2) {
      return cuts;  // only a ring-forming synthon has two or more
    }
    // atoms on a path between two exits: they end up in the assembled ring
    std::set<unsigned int> onPath;
    for (size_t i = 0; i < exits.size(); ++i) {
      for (size_t j = i + 1; j < exits.size(); ++j) {
        const auto path = MolOps::getShortestPath(syn, static_cast<int>(exits[i]),
                                                  static_cast<int>(exits[j]));
        for (int a : path) {
          if (a >= 0) {
            onPath.insert(static_cast<unsigned int>(a));
          }
        }
      }
    }
    for (const auto bond : syn.bonds()) {
      const bool a = onPath.count(bond->getBeginAtomIdx()) > 0;
      const bool b = onPath.count(bond->getEndAtomIdx()) > 0;
      // Found the spinach, add the cut bond here
      if (a != b && bond->getBeginAtom()->getAtomicNum() > 1 &&
          bond->getEndAtom()->getAtomicNum() > 1) {
        cuts.push_back(bond->getIdx());
      }
    }
    return cuts;
  };

  if (rings && isCoarseAssembly(params.embedStyle)) {
    // Log that we found ring formation
    std::printf(
        "[gensynthonlib] this synthon reaction closes %u ring(s): pre-creating final rigid fragments\n",
        rings);
  }
  params.prefillEmbedder = prefill;
  
  // Embedding overrides, mostly heuristic
  // If we miss some ring forms, this might be where it happens.
  //  acylicy count - uses expensive embedding for a whole synthon
  //  RMSD prune radius [0.1 A]
  if (acyclicConfs > 0 || fragRms > 0) {
    auto cp = params.confgen.embedding.classParams;
    for (auto &kv : cp) {
      if (fragRms > 0) {
        kv.second.rmsd = fragRms;
      }
      if (acyclicConfs > 0 && kv.first == FragmentClass::Acyclic) {
        kv.second.maxConfs = acyclicConfs;
      }
    }
    params.confgen.embedding.classParams = cp;
  }
  params.storeEmbedder = true;
  if (productConfs > 0) {
    params.confgen.numOutputConfs = productConfs;
  }

  size_t nSynthons = 0;
  double space = 1.0;
  for (const auto &pos : bbs) {
    nSynthons += pos.size();
    space *= pos.size();
  }
  std::printf("[gensynthonlib] %s/%s: %zu positions, %zu synthons -> %.3g products\n",
              synthons.c_str(), rxn.c_str(), bbs.size(), nSynthons, space);
  std::printf("[gensynthonlib] style=%s prefill=%s\n",
              synthonEmbedStyleName(params.embedStyle),
              prefill ? "yes" : "no");
  std::fflush(stdout);

  // We always add to the shared embedder
  std::shared_ptr<Embedder> shared;
  if (!embedderPath.empty()) {
    EmbedderParams flp = params.confgen.embedding;
    flp.randomSeed = params.confgen.randomSeed;
    flp.ffVariant = params.confgen.zipper.ffVariant;
    shared = std::make_shared<Embedder>(flp);
    std::ifstream in(embedderPath, std::ios_base::binary);
    if (in) {
      shared->initFromStream(in);
      std::printf("[gensynthonlib] shared embedder: %zu entries on entry\n",
                  shared->size());
    }
    params.confgen.embedder = shared;
    params.prefillEmbedder = false;  
    params.storeEmbedder = false;    
  }

  const double t0 = nowMs();
  EnumerateSynthons3D lib(bbs, params);
  const double buildMs = nowMs() - t0;
  if (!lib.isValid()) {
    std::cerr << "library is not valid: the synthons do not share one "
                 "bonding label schemes\n";
    return 1;
  }
  std::printf("[gensynthonlib] built in %.1f s\n", buildMs / 1000.0);
  if (rings && isCoarseAssembly(lib.params3D().embedStyle)) {
    size_t tagged = 0, cutTotal = 0;
    for (unsigned int p = 0; p < lib.arity(); ++p) {
      const auto &row = lib.getReagents()[p];
      for (unsigned int i = 0; i < row.size(); ++i) {
        if (!row[i]) {
          continue;
        }
        auto cuts = ringAttachmentBonds(*row[i]);
        if (!cuts.empty()) {
          cutTotal += cuts.size();
          ++tagged;
          lib.setSynthonCutBonds(p, i, std::move(cuts));
        }
      }
    }
    std::printf(
        "[gensynthonlib] precomputed cut bonds for %zu ring forming synthons (%zu bonds)\n",
        tagged, cutTotal);
  }
  // n.b. we need to track unembeddable synthons so the search doesn't keep
  //  trying to make them.  Report them per position and price them in
  //  products, because "15 dead synthons" and "111k dead products" are the
  //  same fact and only the second one tells you whether to care.
  if (lib.embedder()) {
    if (const size_t bad = lib.embedder()->numUnembeddable()) {
      std::printf(
		  "[gensynthonlib] WARNING: %zu synthon(s) marked as unusable\n",
		  bad);

      double reachable = 1.0;
      for (unsigned int p = 0; p < lib.arity(); ++p) {
        unsigned int nDead = 0;
        std::string examples;
        for (unsigned int r = 0; r < lib.numReagents(p); ++r) {
          if (!lib.synthonUnusable(p, r)) {
            continue;
          }
          ++nDead;
          if (nDead <= 3) {
            examples += (examples.empty() ? "" : ", ") +
                        MolToSmiles(*lib.getReagents()[p][r]);
          }
        }
        if (nDead) {
          std::printf("[gensynthonlib]   position %u: %u of %u marked unusable  (%s%s)\n",
                      p, nDead, lib.numReagents(p), examples.c_str(),
                      nDead > 3 ? ", ..." : "");
        }
        reachable *= (double)(lib.numReagents(p) - nDead) /
                     (double)std::max(1u, lib.numReagents(p));
      }
      std::printf(
          "[gensynthonlib]   %.2f%% of products are unbuildable and are now "
          "will be skipped\n",
          100.0 * (1.0 - reachable));
    }
  }
  if (shared && prefill) {
    const double p0 = nowMs();
    const unsigned int added = lib.prefill(threads);
    std::printf(
        "[gensynthonlib] prefill (%u threads): +%u new entries, %zu total, "
        "%.1f s\n",
        threads, added, shared->size(), (nowMs() - p0) / 1000.0);

    std::ofstream out(embedderPath, std::ios_base::binary);
    if (out) {
      shared->serialize(out);
    }
  }

  const double t1 = nowMs();
  std::ofstream out(output, std::ios::binary);
  if (!out) {
    std::cerr << "cannot write " << output << "\n";
    return 1;
  }
  lib.toStream(out);
  out.close();
  std::printf("[gensynthonlib] wrote %s in %.1f s\n", output.c_str(),
              (nowMs() - t1) / 1000.0);
  return 0;
}
