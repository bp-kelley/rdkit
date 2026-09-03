//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  genSynthonLib -- build a 3D synthon library ONCE and save it.
//
//  Production searches never pay for embedding: the library is prepared
//  offline, with every synthon fragment already in its cache, and searches
//  load it.  This is that offline step, and it is also what makes algorithm
//  comparisons honest -- prefill is a fixed cost that would otherwise sit
//  inside whichever arm ran first.
//
//  The saved file carries its embed STYLE.  A coarse library caches whole
//  synthons and a full one caches their rotatable-bond pieces, so the two are
//  not interchangeable: loading one as the other misses every lookup and
//  silently re-embeds on demand.
//
#include <Confgen/SynthonSearch/EnumerateSynthons3D.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
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
               "  -o/--output FILE     where to write the library\n"
               "      --acyclic-confs N  conformers kept per ACYCLIC fragment.  The\n"
      "                       default recipe keeps 1, which is right for FINE\n"
      "                       fragments (no internal rotors) but freezes a\n"
      "                       whole COARSE synthon to one internal geometry\n"
      "      --frag-rmsd F    RMSD radius (A) for the pool -> keepN diverse\n"
      "                       prune, all classes.  Default recipe is 0.1\n"
      "      --style S        coarse | full (default full).  coarse\n"
               "                       cuts ONLY at the synthon junctions, so\n"
               "                       each synthon is embedded whole\n"
               "      --max-per-pos N  cap reagents per position (0 = all)\n"
               "      --product-confs N  conformers per assembled product\n"
               "      --no-prefill     skip embedding (search pays instead)\n"
               "  -h/--help\n";
}

//! Header-driven reader; see synthonSearchBench for the same loader.
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
  std::string synthons, rxn, output, style = "full";
  size_t maxPerPos = 0;
  int acyclicConfs = -1;   // -1 = leave the class recipe alone
  double fragRms = -1.0;   // -1 = leave the prune radius alone
  int productConfs = -1;
  bool prefill = true;

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
      if (style != "coarse" && style != "full") {
        std::cerr << "--style must be coarse or full\n";
        return 1;
      }
    } else if (a == "--max-per-pos") {
      maxPerPos = static_cast<size_t>(std::stoul(val("--max-per-pos")));
    } else if (a == "--product-confs") {
      productConfs = std::stoi(val("--product-confs"));
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

  EnumerateSynthons3DParams params;
  params.embedStyle = style == "coarse" ? SynthonEmbedStyle::Coarse
                                        : SynthonEmbedStyle::Full;
  params.prefillFraglib = prefill;
  // Class-recipe overrides.  Both are hypotheses under test: the acyclic
  // conformer count (fine-fragment recipe applied to whole synthons) and the
  // prune radius (0.1 A retains near-duplicates, so keepN fills with one
  // basin's variations).
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
  params.storeFraglib = true;
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
  std::printf("[gensynthonlib] style=%s prefill=%s\n", style.c_str(),
              prefill ? "yes" : "no");
  std::fflush(stdout);

  const double t0 = nowMs();
  EnumerateSynthons3D lib(bbs, params);
  const double buildMs = nowMs() - t0;
  if (!lib.isValid()) {
    std::cerr << "library is not valid: the synthons do not share one "
                 "exit-vector label scheme\n";
    return 1;
  }
  std::printf("[gensynthonlib] built in %.1f s\n", buildMs / 1000.0);

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
