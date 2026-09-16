//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  genFragLib -- build a serialized Embedder (cache of embedded rigid fragments)
//  from a SMILES file of molecules.  Each molecule is cut into its rigid
//  fragments exactly as FragmentConfGen does at run time, and each distinct
//  fragment is embedded once and cached (keyed by canonical enhanced-stereo
//  CXSMILES).  Multi-threaded over molecules.
//
//  Usage:
//    genFragLib --input mols.smi --output frags.embedder [options]
//    genFragLib --input more.smi --add frags.embedder --output frags2.embedder
//
//  Options (embedding options default to the Embedder defaults):
//    -i, --input FILE      SMILES file, one molecule per line (first token
//    used) -o, --output FILE     output fragment-library file (required) -a,
//    --add FILE        existing fragment library to load and extend; its
//                          embedding params MUST match the requested ones
//                          (error otherwise).  If no embedding options are
//                          given on the command line, the loaded library's
//                          params are adopted.
//        --num-confs N     conformers embedded per fragment
//        --embed dg|etkdg  fragment embedding method
//        --seed N          embedder random seed
//        --ff VARIANT      force field (default MMFF94).  MUST match the
//                          consumer's joiner.ffVariant: it sets the
//                          minimised geometry AND the stored energies,
//                          and sameEmbeddingType() now enforces it.
//        --minimize        MMFF-minimize each fragment (electrostatics off)
//        --min-iters N     max iterations for --minimize
//        --sample-trivial-rotors  also cut 3-fold symmetric-top rotors
//    -t, --threads N       worker threads (default: hardware concurrency)
//    -h, --help
//
#include <Confgen/Embedder/Embedder.h>
#include <Confgen/FragmentConfGen.h>

#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/RDLog.h>

#include <atomic>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using namespace RDKit;

namespace {

void usage(const char *argv0) {
  std::cerr
      << "usage: " << argv0
      << " --input mols.smi --output frags.embedder [options]\n\n"
         "  -i/--input FILE     SMILES file (one molecule per line)\n"
         "  -o/--output FILE    output fragment-library file (required)\n"
         "  -a/--add FILE       existing library to extend (params must match)\n"
         "      --num-confs N    conformers per fragment\n"
         "      --embed dg|etkdg fragment embedding method\n"
         "      --seed N         embedder random seed\n"
         "      --ff VARIANT     force field: MMFF94, MMFF94s, MMFF94_TOR,\n"
         "                       MMFF94s_TOR (default MMFF94).  MUST match the\n"
         "                       consumer's joiner.ffVariant -- it sets the\n"
         "                       minimised geometry and the stored energies\n"
         "      --minimize       MMFF-minimize fragments (electrostatics off)\n"
         "      --min-iters N    max iters for --minimize\n"
         "      --perclass       per-fragment-class embedding recipe\n"
         "      --whole-acyclic  keep acyclic groups whole (cut ring-adjacent only)\n"
         "      --sample-trivial-rotors  also cut 3-fold symmetric tops\n"
         "      --sdf            write a readable SDF (one record per conformer,\n"
         "                       SD props fragment_smiles/conf_id/mmff_energy)\n"
         "                       instead of the binary fragment library\n"
         "  -t/--threads N       worker threads (default: hw concurrency)\n"
         "  -h/--help\n";
}

std::string describe(const EmbedderParams &p) {
  std::string s;
  s += "ff=" + p.ffVariant;
  s += " numConfs=" + std::to_string(p.numConfsPerFragment);
  s += " embed=" + std::string(p.fragmentEmbedMode == FragmentEmbedMode::ETKDG
                                   ? "etkdg"
                                   : "dg");
  s += " seed=" + std::to_string(p.randomSeed);
  s += " minimize=" +
       std::string(p.minimizeMode != FragmentMinimize::None ? "1" : "0");
  s += " minIters=" + std::to_string(p.minimizeMaxIters);
  s += " acyclicPool=" +
       std::to_string(p.classParams.at(FragmentClass::Acyclic).maxSamples);
  return s;
}

}  // namespace

int main(int argc, char **argv) {
  std::string inFile, outFile, addFile;
  EmbedderParams params;  // defaults = current Embedder defaults
  bool anyEmbedFlag = false;
  bool sampleTrivial = false;
  bool wholeAcyclic = false;  // --whole-acyclic: keep acyclic groups whole
  bool sdfOut =
      false;  // --sdf: write a readable SDF instead of the binary embedder
  unsigned int nThreads = std::thread::hardware_concurrency();
  if (nThreads == 0) {
    nThreads = 1;
  }

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto val = [&](const char *name) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << name << " requires a value\n";
        std::exit(2);
      }
      return argv[++i];
    };
    if (a == "-i" || a == "--input") {
      inFile = val("--input");
    } else if (a == "-o" || a == "--output") {
      outFile = val("--output");
    } else if (a == "-a" || a == "--add") {
      addFile = val("--add");
    } else if (a == "--num-confs") {
      params.numConfsPerFragment =
          static_cast<unsigned int>(std::stoul(val("--num-confs")));
      anyEmbedFlag = true;
    } else if (a == "--embed") {
      const std::string e = val("--embed");
      params.fragmentEmbedMode = (e == "dg" || e == "DG")
                                     ? FragmentEmbedMode::DG
                                     : FragmentEmbedMode::ETKDG;
      anyEmbedFlag = true;
    } else if (a == "--seed") {
      params.randomSeed = std::stoi(val("--seed"));
      anyEmbedFlag = true;
    } else if (a == "--ff") {
      params.ffVariant = val("--ff");
      if (!isValidFF(params.ffVariant)) {
        std::cerr << "genFragLib: unknown --ff \"" << params.ffVariant
                  << "\" (MMFF94, MMFF94s, MMFF94_TOR, MMFF94s_TOR)\n";
        return 1;
      }
      anyEmbedFlag = true;
    } else if (a == "--minimize") {
      params.minimizeMode = FragmentMinimize::Full;
      anyEmbedFlag = true;
    } else if (a == "--min-iters") {
      params.minimizeMaxIters =
          static_cast<unsigned int>(std::stoul(val("--min-iters")));
      anyEmbedFlag = true;
    } else if (a == "--pool") {
      //! Set a default pool size for ALL fragment classes
      params.setFlatPool(static_cast<int>(std::stoul(val("--pool"))));
      anyEmbedFlag = true;
    } else if (a == "--perclass") {
      params.perClassEmbedding = true;
      anyEmbedFlag = true;
    } else if (a == "--sdf") {
      sdfOut = true;
    } else if (a == "--whole-acyclic") {
      wholeAcyclic = true;
    } else if (a == "--sample-trivial-rotors") {
      sampleTrivial = true;
    } else if (a == "-t" || a == "--threads") {
      nThreads = static_cast<unsigned int>(std::stoul(val("--threads")));
      if (nThreads == 0) {
        nThreads = 1;
      }
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else if (!a.empty() && a[0] != '-' && inFile.empty()) {
      inFile = a;  // positional: input
    } else if (!a.empty() && a[0] != '-' && outFile.empty()) {
      outFile = a;  // positional: output
    } else {
      std::cerr << "unknown argument: " << a << "\n";
      usage(argv[0]);
      return 2;
    }
  }
  if (inFile.empty() || outFile.empty()) {
    usage(argv[0]);
    return 2;
  }

  auto lib = std::make_shared<Embedder>(params);

  // Extend an existing library: load it, enforce matching params.
  if (!addFile.empty()) {
    std::ifstream in(addFile);
    if (!in) {
      std::cerr << "error: cannot open --add library: " << addFile << "\n";
      return 1;
    }
    try {
      lib->initFromStream(in);  // sets lib params + loads entries
    } catch (const std::exception &e) {
      std::cerr << "error: failed to read --add library: " << e.what() << "\n";
      return 1;
    }
    if (anyEmbedFlag && lib->params() != params) {
      std::cerr << "error: --add library params differ from the requested "
                   "params; refusing to mix.\n"
                << "  requested: " << describe(params) << "\n"
                << "  library:   " << describe(lib->params()) << "\n";
      return 1;
    }
    params = lib->params();  // adopt (they now govern new embeds)
    std::cerr << "loaded " << lib->size() << " fragments from " << addFile
              << " (" << describe(params) << ")\n";
  }

  // Read input molecules.
  std::vector<std::string> smis;
  {
    std::ifstream in(inFile);
    if (!in) {
      std::cerr << "error: cannot open input: " << inFile << "\n";
      return 1;
    }
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty()) {
        continue;
      }
      const std::string smi = line.substr(0, line.find_first_of(" \t"));
      if (!smi.empty()) {
        smis.push_back(smi);
      }
    }
  }
  std::cerr << "read " << smis.size() << " molecules from " << inFile << "\n";
  std::cerr << "embedding with " << describe(params) << " on " << nThreads
            << " thread(s)\n";

  FragmentConfGenParams gp;
  gp.embedder = lib;
  gp.sampleTrivialRotors = sampleTrivial;
  gp.wholeAcyclicFragments = wholeAcyclic;
  const FragmentConfGen gen(gp);

  // Warm up any the generator when threading.
  if (nThreads > 1) {
    std::unique_ptr<ROMol> m(SmilesToMol("CCOc1ccccc1"));
    if (m) {
      try {
        gen.fragmentAndEmbed(*m);
      } catch (...) {
      }
    }
  }

  std::atomic<size_t> next{0};
  std::atomic<size_t> done{0};
  std::atomic<size_t> failed{0};
  auto worker = [&]() {
    size_t i;
    while ((i = next.fetch_add(1)) < smis.size()) {
      std::unique_ptr<ROMol> mol;
      try {
        mol.reset(SmilesToMol(smis[i]));
      } catch (...) {
      }
      if (!mol) {
        failed.fetch_add(1);
      } else {
        try {
          gen.fragmentAndEmbed(*mol);  // populates the shared lib
        } catch (...) {
          failed.fetch_add(1);
        }
      }
      const size_t d = done.fetch_add(1) + 1;
      if (d % 2000 == 0) {
        std::cerr << "  " << d << "/" << smis.size() << " molecules, "
                  << lib->size() << " fragments\n";
      }
    }
  };

  std::vector<std::thread> pool;
  for (unsigned int t = 0; t < nThreads; ++t) {
    pool.emplace_back(worker);
  }
  for (auto &th : pool) {
    th.join();
  }
  std::cerr << "processed " << done.load() << " molecules (" << failed.load()
            << " failed to parse); " << lib->size() << " unique fragments\n";

  // Write the library out.
  std::ofstream out(outFile);
  if (!out) {
    std::cerr << "error: cannot open output: " << outFile << "\n";
    return 1;
  }
  if (sdfOut) {
    lib->writeSDF(out);
  } else {
    lib->serialize(out);
  }
  out.close();
  if (!out) {
    std::cerr << "error: failed while writing " << outFile << "\n";
    return 1;
  }
  std::cerr << "wrote " << lib->size() << " fragments to " << outFile << "\n";
  return 0;
}
