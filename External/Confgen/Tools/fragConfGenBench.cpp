//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Timing / measurement harness for the fragment conformer
//  generator.  This is NOT a correctness test (those live in
//  testFragmentConfGen.cpp); it is a standalone benchmark that sweeps real
//  drug-like molecules from canonSmiles.long.smi, bucketed by rotatable-bond
//  count, and reports each method's conformer count, single-point MMFF strain,
//  coverage of the ETKDG reference ensembles and wall-clock time.
//
//  It is built but never run by the automated suite (no add_test); run a
//  profile by name:
//
//      fragConfGenBench breakdown     # full rotor-bucket method comparison
//      fragConfGenBench thompson2     # Thompson single-shot at high rotors
//      fragConfGenBench thompson      # energy-window sweep, tree, 6-8 rotors
//      fragConfGenBench strain        # crowded-junction strain sweep
//      fragConfGenBench census        # symmetric-top rotor census + skip
//      timing fragConfGenBench diag          # reference-scored bucket
//      diagnostic
//
//  Needs RDBASE set (for canonSmiles.long.smi).
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Confgen/Utils/TheobaldRmsd.h>
#include <Confgen/Utils/SymmetricRmsd.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <Confgen/FragmentConfGen.h>
#include <Confgen/Search/InterFragScore.h>
#include <Confgen/Sampler/TorsionSampler.h>
#include <Confgen/Embedder/Fraglib.h>
#include <Confgen/Joiner/FragmentJoiner.h>
#include <Confgen/Search/RigidRotorSearch.h>  // runRigidRotorSearch (xtalrecon)
#include <Confgen/Utils/ParamsIO.h>  // FRAGCG_PARAMS: serialized-params knob overrides
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/StretchBend.h>  // StretchBendContrib (inter-fragment SB, termtable)
#include <ForceField/MMFF/Params.h>  // MMFFStbn/MMFFBond/MMFFAngle
#include <RDGeneral/RDLog.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <unistd.h>

using namespace RDKit;

namespace {

//! The Hamburg torsion library (few chemically-preferred angles per bond) from
//! an external data dir; env-overridable (FCG_HAMBURG_LIB).  maxAngles
//! caps the arms per bond (0 =
//! all).
std::shared_ptr<TorsionSampler> hamburgSampler(unsigned int maxAngles = 0) {
  // No default path: the Hamburg XML is external data we do not ship, so it
  // must be pointed at explicitly.  An unset variable yields an invalid
  // sampler, which is reported below rather than silently degrading to a
  // uniform grid.
  const char *e = std::getenv("FCG_HAMBURG_LIB");
  const std::string lib = e ? e : "";
  auto s = std::make_shared<TorsionLibrarySampler>(lib, 30.0, maxAngles);
  // A mistyped path used to be invisible (silent uniform-grid degradation);
  // say so loudly instead of reporting a TorLib arm that never ran.
  if (!s->isValid()) {
    std::fprintf(stderr, "[bench] TorsionLibrarySampler(\"%s\"): %s\n",
                 lib.c_str(), torsionSamplerStatusMessage(s->status()));
    std::fflush(stderr);
  }
  if (std::getenv("TORLIB_NOTOL"))
    s->options.toleranceRanges = false;  // A/B: sparse wells
  if (std::getenv("TORLIB_TOL2"))
    s->options.wideToleranceRanges = true;  // A/B: +tolerance2 wells
  if (std::getenv("TORLIB_DUMP"))
    s->options.dump = true;  // per-bond angle breakdown
  return s;
}

//! SMIRNOFF (OpenFF .offxml) sampler; path from FCG_SMIRNOFF (CC-BY data,
//! external).
std::shared_ptr<TorsionSampler> smirnoffSampler(unsigned int maxAngles = 0) {
  std::string path;
  if (const char *e = std::getenv("FCG_SMIRNOFF")) path = e;
  return std::make_shared<SmirnoffTorsionSampler>(path, 60.0, maxAngles);
}

std::string testDataPath(const std::string &relPath) {
  if (const char *rdBase = std::getenv("RDBASE")) {
    return std::string(rdBase) + "/" + relPath;
  }
  return "../" + relPath;
}

//! The recommended default torsion sampler: composite ETKDG (RDKit CrystalFF
//! preferred angles) + a coarse uniform-grid backstop.  Self-contained (no
//! external XML), clean geometry, no crowded-junction outliers.  Matches the
//! FragmentConfGenParams::torsionSampler default.
std::shared_ptr<TorsionSampler> defaultSampler() {
  return std::make_shared<CompositeTorsionSampler>(
      std::vector<std::shared_ptr<TorsionSampler>>{
          std::make_shared<ETKDGTorsionSampler>(),
          std::make_shared<UniformTorsionSampler>(60.0)});
}

//! The current best-default configuration, kept in ONE place so every
//! experiment is measured against the same reference (and nobody panics at a
//! scary number from a deliberately-weak experimental config).  This is what
//! "DEFAULT" rows in the profiles run.
FragmentConfGenParams bestDefaultParams(int seed = 42) {
  FragmentConfGenParams
      p;  // library defaults (composite-ETKDG sampler, CoordsOnly)
  p.numOutputConfs = 50;
  p.energyWindow = 10.0;
  p.randomSeed = seed;
  // Match the PRODUCTION library default: bare ETKDG preferred angles (the
  // joiner adds its own low-weight 60-deg backstop).  The old defaultSampler()
  // here was the Composite(ETKDG+Uniform60) union, which flattens the informed
  // prior and measurably hurts Thompson sampling -- so this benchmark was not
  // measuring the shipped default.
  p.search.torsionSampler = std::make_shared<ETKDGTorsionSampler>();
  return p;
}

double nowMs() {
  using namespace std::chrono;
  return duration<double, std::milli>(steady_clock::now().time_since_epoch())
      .count();
}

//! Validate params and STOP the run if they fail
void requireValidParams(const FragmentConfGenParams &p, const char *where) {
  const std::string err = p.validate();
  if (!err.empty()) {
    std::fprintf(stderr, "\n[%s] FATAL parameter error: %s\n", where,
                 err.c_str());
    std::fflush(stderr);
    std::exit(2);
  }
}

//! ASM_SEARCH env -> RigidRotorSearchMode override (unset -> the runner's own
//! default, i.e. RigidRotorSearchMode::Auto).  Lets any bench pick a specific
//! search:
//!   tree | thompson(ts) | systematic(sys) | merged(union).
RigidRotorSearchMode searchModeFromEnv(
    RigidRotorSearchMode fallback = RigidRotorSearchMode::Auto) {
  const char *e = std::getenv("ASM_SEARCH");
  if (!e) return fallback;
  const std::string s(e);
  if (s == "tree") return RigidRotorSearchMode::Tree;
  if (s == "merged" || s == "union") return RigidRotorSearchMode::Merged;
  if (s == "thompson" || s == "ts") return RigidRotorSearchMode::Thompson;
  if (s == "systematic" || s == "sys") return RigidRotorSearchMode::Systematic;
  return fallback;
}

//! Reference ensemble.  This is the eventual comparison
//! target.

// ---------------------------------------------------------------------------
// Profiles
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Platinum bioactive-conformer reproduction benchmark
// ---------------------------------------------------------------------------

//! Largest connected fragment by heavy-atom count, conformers preserved.
//! Mirrors what FragmentConfGen.build() keeps (largest component), so the
//! reference pose and our generated conformers share the same atom set for
//! RMSD.
ROMOL_SPTR largestFrag(const ROMol &m) {
  std::vector<ROMOL_SPTR> frags =
      MolOps::getMolFrags(m, /*sanitizeFrags=*/false);
  if (frags.empty()) {
    return nullptr;
  }
  ROMOL_SPTR best = frags.front();
  unsigned int bestHeavy = 0;
  for (const auto &f : frags) {
    unsigned int h = 0;
    for (const auto *a : f->atoms()) {
      if (a->getAtomicNum() > 1) {
        ++h;
      }
    }
    if (h > bestHeavy) {
      bestHeavy = h;
      best = f;
    }
  }
  return best;
}

//! Minimum symmetry-aware heavy-atom RMSD from any conformer in `gen` (mols
//! with Hs, one conformer each) to the heavy-atom reference pose `refHeavy`.
//! This is the Platinum reproduction metric: getBestRMS enumerates the graph
//! automorphisms so symmetry-equivalent poses (flipped phenyl, equivalent
//! termini) match, and it aligns by rotation+translation only -- a wrong
//! enantiomer (mirror image) cannot be superposed and so is correctly
//! penalised.  NaN if nothing is comparable. min & median full-molecule MMFF
//! energy (no electrostatics, given variant) over an ensemble of
//! single-conformer molecules -- for comparing the geometric quality of our
//! assembled conformers (a blown-up bond/angle shows as a huge
//! energy). Each molecule is typed independently, so runs are
//! directly comparable.
std::pair<double, double> mmffEnergyStats(const std::vector<ROMOL_SPTR> &confs,
                                          const std::string &variant) {
  std::vector<double> es;
  for (const auto &m : confs) {
    if (!m || m->getNumConformers() == 0) continue;
    auto fn = makeFullFFScoreFn(*m, /*electrostatics=*/false, variant);
    if (!fn) continue;
    const Conformer &conf = m->getConformer();
    const unsigned int n = m->getNumAtoms();
    std::vector<double> buf(3 * static_cast<size_t>(n));
    for (unsigned int a = 0; a < n; ++a) {
      const RDGeom::Point3D &p = conf.getAtomPos(a);
      buf[3 * a] = p.x;
      buf[3 * a + 1] = p.y;
      buf[3 * a + 2] = p.z;
    }
    double e = fn(buf.data(), n);
    if (std::isfinite(e)) es.push_back(e);
  }
  if (es.empty())
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
  std::sort(es.begin(), es.end());
  return {es.front(), es[es.size() / 2]};
}

//! Per-term MMFF94(s) energy decomposition of a single molecule/conformer (no
//! estat): {bondStretch, angleBend, stretchBend, oopBend, torsion, vdW}.  Each
//! term is isolated by enabling ONLY that term on the MMFFMolProperties before
//! building the force field -- so we can see WHERE our assembled geometry is
//! strained (bad bonds vs bad angles vs clashes) rather than just the inflated
//! total.
struct MmffTerms {
  double bond = 0, angle = 0, sb = 0, oop = 0, tors = 0, vdw = 0;
};
MmffTerms mmffDecompose(const ROMol &mol, const std::string &variant) {
  MmffTerms t;
  const bool s = (variant == "MMFF94s" || variant == "MMFF94S");
  auto termE = [&](void (MMFF::MMFFMolProperties::*setter)(bool)) -> double {
    RWMol rw(mol);
    MMFF::MMFFMolProperties props(rw, s ? "MMFF94s" : "MMFF94");
    if (!props.isValid()) return std::numeric_limits<double>::quiet_NaN();
    props.setMMFFBondTerm(false);
    props.setMMFFAngleTerm(false);
    props.setMMFFStretchBendTerm(false);
    props.setMMFFOopTerm(false);
    props.setMMFFTorsionTerm(false);
    props.setMMFFVdWTerm(false);
    props.setMMFFEleTerm(false);
    (props.*setter)(true);
    std::unique_ptr<ForceFields::ForceField> ff(
        MMFF::constructForceField(rw, &props, 1.0e8, 0));
    if (!ff) return std::numeric_limits<double>::quiet_NaN();
    ff->initialize();
    return ff->calcEnergy();
  };
  t.bond = termE(&MMFF::MMFFMolProperties::setMMFFBondTerm);
  t.angle = termE(&MMFF::MMFFMolProperties::setMMFFAngleTerm);
  t.sb = termE(&MMFF::MMFFMolProperties::setMMFFStretchBendTerm);
  t.oop = termE(&MMFF::MMFFMolProperties::setMMFFOopTerm);
  t.tors = termE(&MMFF::MMFFMolProperties::setMMFFTorsionTerm);
  t.vdw = termE(&MMFF::MMFFMolProperties::setMMFFVdWTerm);
  return t;
}

//! Decompose the LOWEST-energy conformer of an ensemble (by full no-estat
//! MMFF).
MmffTerms decomposeMinConf(const std::vector<ROMOL_SPTR> &confs,
                           const std::string &variant) {
  double best = std::numeric_limits<double>::infinity();
  const ROMol *bm = nullptr;
  for (const auto &m : confs) {
    if (!m || m->getNumConformers() == 0) continue;
    auto fn = makeFullFFScoreFn(*m, false, variant);
    if (!fn) continue;
    const Conformer &c = m->getConformer();
    const unsigned int n = m->getNumAtoms();
    std::vector<double> buf(3 * static_cast<size_t>(n));
    for (unsigned int a = 0; a < n; ++a) {
      const auto &p = c.getAtomPos(a);
      buf[3 * a] = p.x;
      buf[3 * a + 1] = p.y;
      buf[3 * a + 2] = p.z;
    }
    double e = fn(buf.data(), n);
    if (std::isfinite(e) && e < best) {
      best = e;
      bm = m.get();
    }
  }
  if (!bm) return {};
  return mmffDecompose(*bm, variant);
}

double minHeavyRms(const std::vector<ROMOL_SPTR> &gen, const ROMol &refHeavy) {
  RWMol ref(refHeavy);
  double best = std::numeric_limits<double>::max();
  for (const auto &g : gen) {
    if (!g || g->getNumConformers() == 0) {
      continue;
    }
    ROMOL_SPTR probe;
    try {
      probe.reset(MolOps::removeHs(*g));
    } catch (...) {
      continue;
    }
    if (probe->getNumAtoms() != ref.getNumAtoms()) {
      continue;  // e.g. we dropped a component the reference still carries
    }
    try {
      best = std::min(best, MolAlign::getBestRMS(*probe, ref));
    } catch (...) {
      // isomorphic-but-unmatched; skip
    }
  }
  return std::isfinite(best) ? best : std::numeric_limits<double>::quiet_NaN();
}

//! Diagnostic: over the energy-sorted output ensemble `gen`, return the ENERGY
//! RANK (0-based index) of the conformer closest to the crystal reference, plus
//! that best RMSD.  Because gen is sorted lowest-energy-first, the rank tells
//! us whether the bioactive pose is a low-energy conformer (rank near 0 ->
//! energy-greedy top-N keeps it) or buried deep (energy is a poor predictor of
//! the bioactive pose).
std::pair<int, double> winnerEnergyRank(const std::vector<ROMOL_SPTR> &gen,
                                        const ROMol &refHeavy) {
  RWMol ref(refHeavy);
  double best = std::numeric_limits<double>::max();
  int bestRank = -1;
  for (size_t i = 0; i < gen.size(); ++i) {
    const auto &g = gen[i];
    if (!g || g->getNumConformers() == 0) continue;
    ROMOL_SPTR probe;
    try {
      probe.reset(MolOps::removeHs(*g));
    } catch (...) {
      continue;
    }
    if (probe->getNumAtoms() != ref.getNumAtoms()) continue;
    try {
      double r = MolAlign::getBestRMS(*probe, ref);
      if (r < best) {
        best = r;
        bestRank = static_cast<int>(i);
      }
    } catch (...) {
    }
  }
  return {bestRank, best};
}

//! Distance-geometry ensemble baseline: plain DG (no torsion knowledge) or
//! ETKDGv3 (CrystalFF torsion preferences), `nConf` confs, one single-conformer
//! mol (with Hs) each, in embed order.  No minimisation, no RMS pruning.
std::vector<ROMOL_SPTR> embedEns(const ROMol &input, bool etkdg,
                                 unsigned int nConf, int seed, double &outMs) {
  std::vector<ROMOL_SPTR> out;
  std::unique_ptr<ROMol> molH(MolOps::addHs(input));
  auto rw = boost::make_shared<RWMol>(*molH);
  DGeomHelpers::EmbedParameters ps =
      etkdg ? DGeomHelpers::EmbedParameters(DGeomHelpers::ETKDGv3)
            : DGeomHelpers::EmbedParameters();  // plain distance geometry
  ps.randomSeed = seed;
  ps.pruneRmsThresh = -1.0;
  RDKit::INT_VECT cids;
  const double t0 = nowMs();
  {
    RDLog::LogStateSetter blocker;
    try {
      DGeomHelpers::EmbedMultipleConfs(*rw, cids, nConf, ps);
    } catch (...) {
    }
  }
  outMs = nowMs() - t0;
  for (const auto cid : cids) {
    auto single = boost::make_shared<RWMol>(*rw);
    single->clearConformers();
    auto *c = new Conformer(rw->getConformer(cid));
    c->setId(0);
    single->addConformer(c, false);
    out.push_back(single);
  }
  return out;
}

//! MMFF-minimise (electrostatics off) every conformer of an ensemble, returning
//! minimised copies.  Untypeable molecules pass through unminimised.
std::vector<ROMOL_SPTR> mmffMin(const std::vector<ROMOL_SPTR> &ens,
                                double &outMs) {
  std::vector<ROMOL_SPTR> out;
  const double t0 = nowMs();
  for (const auto &m : ens) {
    if (!m || m->getNumConformers() == 0) {
      continue;
    }
    auto rw = boost::make_shared<RWMol>(*m);
    try {
      MMFF::MMFFMolProperties props(*rw);
      if (props.isValid()) {
        props.setMMFFEleTerm(false);
        std::unique_ptr<ForceFields::ForceField> ff(
            MMFF::constructForceField(*rw, &props, 1.0e8, 0));
        if (ff) {
          ff->initialize();
          ff->minimize(1000);
        }
      }
    } catch (...) {
    }
    out.push_back(rw);
  }
  outMs = nowMs() - t0;
  return out;
}

//! One method's (rms, count, time) for a single molecule.
struct Cell {
  double rms = std::numeric_limits<double>::quiet_NaN();
  size_t n = 0;
  double ms = 0.0;
};

//! Per-method, per-bucket tally for the Platinum reproduction rates.
struct PlatAcc {
  unsigned int scored = 0;  //!< mols with a finite min-RMSD
  unsigned int failed = 0;  //!< genN==0 or non-comparable graph
  unsigned int lt1 = 0;  //!< min-RMSD < 1.0 A (the headline reproduction rate)
  unsigned int lt2 = 0;  //!< min-RMSD < 2.0 A
  double sumRms = 0.0;
  double genN = 0.0;
  double ms = 0.0;
  std::vector<double> rmss;  //!< for the median
  void add(double rms, size_t nConf, double t) {
    genN += static_cast<double>(nConf);
    ms += t;
    if (std::isfinite(rms)) {
      // Clamp pathological blown-up geometries (huge but finite getBestRMS) so
      // the mean can't overflow; thresholds/median are unaffected (clamp only
      // bites values already far above 2 A).
      const double rr = std::min(rms, 25.0);
      ++scored;
      if (rr < 1.0) ++lt1;
      if (rr < 2.0) ++lt2;
      sumRms += rr;
      rmss.push_back(rr);
    } else {
      ++failed;
    }
  }
  double median() {
    if (rmss.empty()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(rmss.begin(), rmss.end());
    const size_t n = rmss.size();
    return (n % 2) ? rmss[n / 2] : 0.5 * (rmss[n / 2 - 1] + rmss[n / 2]);
  }
};

//! `msOverride` >= 0 replaces the per-mol time (used to show the same method's
//! warm-cache timing against its cold-cache RMSD stats).
void printPlatRow(const char *name, PlatAcc &a, unsigned int nMols,
                  double msOverride = -1.0) {
  const double denom = a.scored ? a.scored : 1;
  const double msPerMol =
      msOverride >= 0.0 ? msOverride : (nMols ? a.ms / nMols : 0.0);
  std::printf(
      "  %-12s %5u %5u  %6.1f%% %6.1f%%  %6.2f  %6.2f  %6.1f  %7.0f\n", name,
      nMols, a.failed, 100.0 * a.lt1 / denom, 100.0 * a.lt2 / denom, a.median(),
      a.scored ? a.sumRms / a.scored : std::numeric_limits<double>::quiet_NaN(),
      nMols ? a.genN / nMols : 0.0, msPerMol);
  std::fflush(stdout);
}

//! Bioactive-conformer reproduction on the Platinum Diverse Dataset: for each
//! crystal ligand pose, generate an ensemble and report the fraction whose best
//! member is within 1.0 / 2.0 A heavy-atom RMSD of the crystal pose, bucketed
//! by rotatable-bond count.  vs plain ETKDG (in-process) and, if
  //! by rotatable-bond count, against plain ETKDG (in-process).
void runPlatinum() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) {
    sdf = e;
  }
  int seed =
      42;  // PLATINUM_SEED overrides -- the ETKDG embed + Thompson search are
  if (const char *e = std::getenv("PLATINUM_SEED"))
    seed = std::atoi(e);  // stochastic
  int maxConfs = 50;  // SIGNED: -1 = AUTO (rotor-driven), 0 = DISABLED
  if (const char *e = std::getenv("PLATINUM_MAXCONFS")) {
    maxConfs = std::atoi(e);
  }
  size_t maxMols = 0;  // 0 = all
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) {
    maxMols = static_cast<size_t>(std::atoi(e));
  }
  // Diagnostic knobs: override per-fragment conf count, turn on fragment
  // MMFF-min (cold builds only -- the warm lib is fixed at build time), and
  // skip the heavy DG/ETKDG baseline matrix for fast ours-only
  // iteration.
  unsigned int nconfFrag = 0;  // 0 = keep default
  if (const char *e = std::getenv("PLATINUM_NCONF")) {
    nconfFrag = static_cast<unsigned int>(std::atoi(e));
  }
  const bool fragMin = std::getenv("PLATINUM_FRAGMIN") != nullptr;
  unsigned int poolFrag = 0;
  if (const char *e = std::getenv("PLATINUM_POOL")) {
    poolFrag = static_cast<unsigned int>(std::atoi(e));
  }
  const bool acyclic = std::getenv("PLATINUM_ACYCLIC") != nullptr;
  const bool doBase = std::getenv("PLATINUM_NOBASE") == nullptr;

  FragmentConfGenParams pCold =
      bestDefaultParams(seed);  // shipped default config
  pCold.numOutputConfs = maxConfs;
  if (std::getenv("PLATINUM_COORDSONLY")) {
    // route "ours" through the new coordinate-only FragmentJoiner backend
    if (const char *e = std::getenv("PLATINUM_ENERGYWINDOW"))
      pCold.energyWindow = std::atof(e);
  }

  // Optional warm fragment cache (production scenario): identical geometry to
  // the cold build (same seed/embed), just embedding served from the cache ->
  // lets us report warm-cache timing next to the cold-cache RMSD.
  std::shared_ptr<Fraglib> lib;
  if (const char *libEnv = std::getenv("FRAGLIB")) {
    lib = std::make_shared<Fraglib>();
    std::ifstream in(libEnv);
    if (!in) {
      std::cerr << "platinum: cannot open FRAGLIB " << libEnv << "\n";
      return;
    }
    try {
      lib->initFromStream(in);
    } catch (const std::exception &e) {
      std::cerr << "platinum: failed to load FRAGLIB: " << e.what() << "\n";
      return;
    }
  }
  FragmentConfGenParams pWarm = pCold;
  if (lib) {
    pWarm.fraglib = lib;
    // PLATINUM_FRAGLIB_GEOM: use the loaded library for the ACCURACY (cold)
    // path too -- substitute its fragment GEOMETRY into the measured ensemble
    // (fragment isolation experiment).  Cache HIT -> loaded conformers;
    // MISS -> embed+add (so the final lib->size() vs the initial count reveals
    // coverage).
    if (std::getenv("PLATINUM_FRAGLIB_GEOM")) {
      std::cerr << "platinum: FRAGLIB_GEOM on -- " << lib->size()
                << " loaded fragments feed the accuracy path\n";
      pCold.fraglib = lib;
    }
  }

  // torlib+backstop variant: Hamburg TorsionLibrary preferred angles unioned
  // with a coarse uniform-grid backstop (kills crowded-junction outliers).
  // Needs the external TorsionLibrary.xml; skipped if it can't be loaded.  Same
  // fraglib works (embedding is sampler-independent; only the assembly torsion
  // sweep changes).
  bool doTlib = false;
  FragmentConfGenParams pTlibCold = pCold;
  FragmentConfGenParams pTlibWarm = pWarm;
  try {
    auto tlib = std::make_shared<CompositeTorsionSampler>(
        std::vector<std::shared_ptr<TorsionSampler>>{
            hamburgSampler(), std::make_shared<UniformTorsionSampler>(60.0)});
    pTlibCold.search.torsionSampler = tlib;
    pTlibWarm.search.torsionSampler = tlib;
    doTlib = true;
  } catch (const std::exception &e) {
    std::cerr << "platinum: torlib+backstop unavailable (" << e.what()
              << "); skipping that variant\n";
  }

  std::cout
      << "\n[platinum] bioactive-conformer reproduction on the Platinum "
         "Diverse Dataset (2017_01)\n"
      << "  config: composite-ETKDG sampler, maxConfs=" << maxConfs
      << ", seed=" << seed
      << (lib ? (", warm FRAGLIB=" + std::to_string(lib->size()) + " fragments")
              : std::string(", cold cache only"))
      << ", ETKDG baseline" << "\n"
      << "  metric: min symmetry-aware heavy-atom RMSD of the ensemble to "
         "the crystal pose\n";

  // The 8 embedding-baseline cells, in report order.  DG/ETKDG x raw/+MMFF x
  // {best of first 50, best of all 100} conformers from a single 100-conf
  // embed.
  static const char *const kBaseLabels[8] = {
      "dg@50",    "dg@100",    "dg+mmff@50",    "dg+mmff@100",
      "etkdg@50", "etkdg@100", "etkdg+mmff@50", "etkdg+mmff@100"};
  struct Bucket {
    const char *label;
    unsigned int lo, hi;
    unsigned int nMols = 0;
    double oursWarmMs =
        0.0;  //!< warm-cache build time (RMSD stats live in `ours`)
    double tlibWarmMs = 0.0;
    PlatAcc ours = {}, tlib = {};
    PlatAcc base[8] = {};
  };
  std::vector<Bucket> buckets = {{"0", 0, 0},     {"1", 1, 1},     {"2", 2, 2},
                                 {"3", 3, 3},     {"4", 4, 4},     {"5", 5, 5},
                                 {"6", 6, 6},     {"7", 7, 7},     {"8", 8, 8},
                                 {"9-11", 9, 11}, {"12+", 12, 999}};
  auto bucketIdx = [&](unsigned int nrot) -> int {
    for (size_t i = 0; i < buckets.size(); ++i) {
      if (nrot >= buckets[i].lo && nrot <= buckets[i].hi) {
        return static_cast<int>(i);
      }
    }
    return static_cast<int>(buckets.size()) - 1;
  };

  // --- Prep every molecule ONCE, single-threaded (SDMolSupplier is not
  // thread-safe; stereo perception + largest-component reduction done here).
  // ---
  struct Job {
    ROMOL_SPTR refHeavy;  //!< heavy-atom crystal pose (the RMSD target)
    ROMOL_SPTR input;     //!< graph only (no conformer) -> generator input
    int bucket = -1;
    std::string smi;
  };
  std::vector<Job> jobs;
  size_t skipped = 0;
  {
    std::unique_ptr<SDMolSupplier> supplier;
    try {
      supplier = std::make_unique<SDMolSupplier>(sdf, /*sanitize=*/true,
                                                 /*removeHs=*/false);
    } catch (...) {
      std::cerr << "platinum: cannot open " << sdf << "\n";
      return;
    }
    RDLog::LogStateSetter blocker;
    while (!supplier->atEnd()) {
      if (maxMols && jobs.size() >= maxMols) {
        break;
      }
      std::unique_ptr<ROMol> raw;
      try {
        raw.reset(supplier->next());
      } catch (...) {
        ++skipped;
        continue;
      }
      if (!raw || raw->getNumConformers() == 0) {
        ++skipped;
        continue;
      }
      try {
        MolOps::assignStereochemistryFrom3D(*raw);
        std::unique_ptr<ROMol> heavy(MolOps::removeHs(*raw));
        ROMOL_SPTR refHeavy = largestFrag(*heavy);
        if (!refHeavy || refHeavy->getNumAtoms() < 2) {
          ++skipped;
          continue;
        }
        auto input = boost::make_shared<RWMol>(*refHeavy);
        input->clearConformers();
        unsigned int nrot = 0;
        try {
          nrot = Descriptors::calcNumRotatableBonds(*refHeavy);
        } catch (...) {
        }
        std::string smi;
                jobs.push_back({refHeavy, input, bucketIdx(nrot), smi});
      } catch (...) {
        ++skipped;
      }
    }
  }
  for (const auto &j : jobs) {
    ++buckets[j.bucket].nMols;
  }
  const unsigned int allMols = static_cast<unsigned int>(jobs.size());

  unsigned int nThreads = std::thread::hardware_concurrency();
  if (const char *e = std::getenv("PLATINUM_THREADS")) {
    nThreads = static_cast<unsigned int>(std::atoi(e));
  }
  nThreads =
      std::max(1u, std::min<unsigned int>(nThreads, allMols ? allMols : 1));
  std::cout << "  parallel: " << nThreads << " threads over " << allMols
            << " molecules (" << skipped << " skipped)"
            << (doTlib ? ", +torlib+backstop variant" : "") << "\n";

  // Per-molecule work -> a self-contained result (no shared writes).
  struct MolResult {
    int bucket = -1;
    double oursRms = std::numeric_limits<double>::quiet_NaN();
    size_t oursN = 0;
    double oursColdMs = 0, oursWarmMs = 0;
    double tlibRms = std::numeric_limits<double>::quiet_NaN();
    size_t tlibN = 0;
    double tlibColdMs = 0, tlibWarmMs = 0;
    Cell base[8];
  };
  auto timedBuild = [&](const FragmentConfGenParams &pp, const ROMol &in,
                        double &ms) -> std::vector<ROMOL_SPTR> {
    const double t0 = nowMs();
    FragmentConfGenResult res;
    try {
      requireValidParams(pp, "runPlatinum");
      res = FragmentConfGen(pp).build(in);
    } catch (...) {
    }
    ms = nowMs() - t0;
    return res.conformers;
  };
  auto processJob = [&](const Job &j, const FragmentConfGenParams &pc,
                        const FragmentConfGenParams &pw,
                        const FragmentConfGenParams &ptc,
                        const FragmentConfGenParams &ptw) -> MolResult {
    MolResult r;
    r.bucket = j.bucket;
    auto oursConfs = timedBuild(pc, *j.input, r.oursColdMs);
    r.oursRms = minHeavyRms(oursConfs, *j.refHeavy);
    r.oursN = oursConfs.size();
    if (lib) {
      timedBuild(pw, *j.input, r.oursWarmMs);  // timing only; same geometry
    }
    if (doTlib) {
      auto confs = timedBuild(ptc, *j.input, r.tlibColdMs);
      r.tlibRms = minHeavyRms(confs, *j.refHeavy);
      r.tlibN = confs.size();
      if (lib) {
        timedBuild(ptw, *j.input, r.tlibWarmMs);
      }
    }
    // Embedding-baseline matrix: DG/ETKDG embedded at 100 confs, then the same
    // ensemble MMFF-minimised; each scored at cutoffs 50 and 100
    // (best-of-first-N).
    if (doBase) {
      double tDg = 0, tEtk = 0, tDgm = 0, tEtkm = 0;
      auto dg = embedEns(*j.input, /*etkdg=*/false, 100, seed, tDg);
      auto etk = embedEns(*j.input, /*etkdg=*/true, 100, seed, tEtk);
      auto dgm = mmffMin(dg, tDgm);
      auto etkm = mmffMin(etk, tEtkm);
      auto cut = [](const std::vector<ROMOL_SPTR> &v, size_t k) {
        return std::vector<ROMOL_SPTR>(v.begin(),
                                       v.begin() + std::min(k, v.size()));
      };
      auto fill = [&](Cell &c, const std::vector<ROMOL_SPTR> &e, double ms) {
        auto s = e;
        c.rms = minHeavyRms(s, *j.refHeavy);
        c.n = s.size();
        c.ms = ms;
      };
      fill(r.base[0], cut(dg, 50), tDg);
      fill(r.base[1], cut(dg, 100), tDg);
      fill(r.base[2], cut(dgm, 50), tDg + tDgm);
      fill(r.base[3], cut(dgm, 100), tDg + tDgm);
      fill(r.base[4], cut(etk, 50), tEtk);
      fill(r.base[5], cut(etk, 100), tEtk);
      fill(r.base[6], cut(etkm, 50), tEtk + tEtkm);
      fill(r.base[7], cut(etkm, 100), tEtk + tEtkm);
    }
        return r;
  };
  auto accumulate = [&](const MolResult &r) {
    Bucket &b = buckets[r.bucket];
    b.ours.add(r.oursRms, r.oursN, r.oursColdMs);
    b.oursWarmMs += r.oursWarmMs;
    if (doTlib) {
      b.tlib.add(r.tlibRms, r.tlibN, r.tlibColdMs);
      b.tlibWarmMs += r.tlibWarmMs;
    }
    for (int k = 0; k < 8; ++k) {
      b.base[k].add(r.base[k].rms, r.base[k].n, r.base[k].ms);
    }
      };
  // Fresh per-thread params: each worker owns its samplers so concurrent
  // build()s never share sampler state (the fraglib IS shared -- it is
  // thread-safe).
  auto makeParams = [&](FragmentConfGenParams &pc, FragmentConfGenParams &pw,
                        FragmentConfGenParams &ptc,
                        FragmentConfGenParams &ptw) {
    pc = bestDefaultParams(seed);
    pc.numOutputConfs = maxConfs;
    // Reproduction benchmark scores with MMFF94s: the crystal pose is
    // time-averaged, which MMFF94s reproduces best (NB not necessarily the best
    // variant for a single BIOACTIVE conformer -- the library/joiner default
    // stays MMFF94).  ASM_MMFF94 reverts.
    pc.joiner.ffVariant = std::getenv("ASM_MMFF94") ? "MMFF94" : "MMFF94s";
    if (std::getenv("PLATINUM_COORDSONLY")) {
      if (const char *e = std::getenv("ASM_VDWCUTOFF"))
        pc.joiner.interFragVdwCutoff = std::atof(e);
      if (const char *e = std::getenv("PLATINUM_ENERGYWINDOW"))
        pc.energyWindow = std::atof(e);
    }
    if (nconfFrag) {
      pc.embedding.numConfsPerFragment = nconfFrag;
    }
    // Fragment minimiser mode comes from params (default ShrugScore, perClass
    // ON).  DG embedding needs Full (keep the minimised geometry, no torsions
    // to trust); PLATINUM_NOMIN forces None (single-point).  PLATINUM_FRAGMIN
    // is retained only to force Full explicitly.
    if (std::getenv("PLATINUM_DG")) {
      pc.embedding.fragmentEmbedMode =
          FragmentEmbedMode::DG;  // else ETKDG (bestDefaultParams)
      pc.embedding.minimizeMode = FragmentMinimize::Full;
    }
    if (fragMin) pc.embedding.minimizeMode = FragmentMinimize::Full;
    if (std::getenv("PLATINUM_NOMIN"))
      pc.embedding.minimizeMode = FragmentMinimize::None;
    if (poolFrag) {
      pc.embedding.setFlatPool(
          poolFrag);  // diversity-only unless PLATINUM_FRAGMIN also set
    }
    pc.wholeAcyclicFragments = acyclic;
    if (const char *e = std::getenv("ASM_FPNOVELTY")) {
      pc.search.thompson.noveltyAngleDeg = std::atof(e);
    }
    if (const char *e = std::getenv("ASM_PERCLASS")) {
      pc.embedding.perClassEmbedding = std::atoi(e) != 0;
    }
    if (const char *e = std::getenv("ASM_SYMDEDUP")) {
      pc.search.finalSymmetryDedup = std::atoi(e) != 0;
    }
    // Ceiling / search-width knobs (measure the best the current fraglib can
    // do).
    if (const char *e = std::getenv("ASM_FRAGBRANCH"))
      pc.search.fragConfBranch = std::atoi(e);
    if (const char *e = std::getenv("ASM_ROOTSEEDS"))
      pc.search.rootSeeds = std::atoi(e);
    if (const char *e = std::getenv("ASM_BEAM"))
      pc.search.tree.beamWidth = std::atoi(e);
    if (const char *e =
            std::getenv("ASM_TSAUTO"))  // 0 + budget 0 -> exhaustive beam
      pc.search.thompson.autoBudget = std::atoi(e) != 0;
    if (const char *e = std::getenv("ASM_TSBUDGET"))
      pc.search.thompsonBudget = std::atoi(e);
    if (const char *e = std::getenv("ASM_TSMAX"))
      pc.search.thompson.maxBudget = std::atoi(e);
    if (const char *e = std::getenv("ASM_TSPERROTOR"))
      pc.search.thompson.perRotor = std::atoi(e);
    if (const char *e = std::getenv("ASM_TSPERFRAG"))
      pc.search.thompson.perFragConf = std::atoi(e);
    pw = pc;
    if (lib) {
      pw.fraglib = lib;
    }
    ptc = pc;
    ptw = pw;
    if (doTlib) {
      auto ts = std::make_shared<CompositeTorsionSampler>(
          std::vector<std::shared_ptr<TorsionSampler>>{
              hamburgSampler(), std::make_shared<UniformTorsionSampler>(60.0)});
      ptc.search.torsionSampler = ts;
      ptw.search.torsionSampler = ts;
    }
  };

  const double runStart = nowMs();
  // Warmup mol 0 single-threaded so every lazy global init fires before
  // threading.
  if (!jobs.empty()) {
    FragmentConfGenParams pc, pw, ptc, ptw;
    makeParams(pc, pw, ptc, ptw);
    accumulate(processJob(jobs[0], pc, pw, ptc, ptw));
  }
  std::atomic<size_t> nextJob{1};
  std::atomic<size_t> nDone{jobs.empty() ? size_t(0) : size_t(1)};
  std::mutex accMtx;
  auto worker = [&]() {
    FragmentConfGenParams pc, pw, ptc, ptw;
    makeParams(pc, pw, ptc, ptw);
    for (;;) {
      const size_t i = nextJob.fetch_add(1);
      if (i >= jobs.size()) {
        break;
      }
      if (const char *e = std::getenv("PLATINUM_ONLY")) {
        if (i != static_cast<size_t>(std::atol(e))) {
          continue;
        }
      }
      if (const char *e = std::getenv("PLATINUM_FROM")) {
        if (i < static_cast<size_t>(std::atol(e))) {
          continue;
        }
      }
      if (std::getenv("PLATINUM_TRACE")) {
        std::string s;
        try {
          s = MolToSmiles(*jobs[i].refHeavy);
        } catch (...) {
          s = "<smiles-failed>";
        }
        std::fprintf(stderr, "[trace] mol %zu natoms=%u %s\n", i,
                     jobs[i].refHeavy->getNumAtoms(), s.c_str());
        std::fflush(stderr);
      }
      MolResult r = processJob(jobs[i], pc, pw, ptc, ptw);
      if (std::getenv("PLATINUM_TRACE")) {
        unsigned int nrot = 0;
        try {
          nrot = Descriptors::calcNumRotatableBonds(*jobs[i].refHeavy);
        } catch (...) {
        }
        std::string nm = jobs[i].refHeavy->hasProp("_Name")
                             ? jobs[i].refHeavy->getProp<std::string>("_Name")
                             : "?";
        std::fprintf(
            stderr,
            "[trace] mol %2zu rot=%2u nAt=%2u  ours=%.3f  tlib=%.3f  name=%s\n",
            i, nrot, jobs[i].refHeavy->getNumAtoms(), r.oursRms, r.tlibRms,
            nm.c_str());
        std::fflush(stderr);
      }
      {
        std::lock_guard<std::mutex> lk(accMtx);
        accumulate(r);
      }
      const size_t d = ++nDone;
      if (d % 100 == 0) {
        std::fprintf(stderr, "  ...%zu/%u mols, %.0fs elapsed\n", d, allMols,
                     (nowMs() - runStart) / 1000.0);
      }
    }
  };
  {
    std::vector<std::thread> pool;
    for (unsigned int t = 0; t < nThreads; ++t) {
      pool.emplace_back(worker);
    }
    for (auto &t : pool) {
      t.join();
    }
  }

  // Report.
  std::cout << "\n  columns: nMols failed  %<1.0A  %<2.0A  medRMS  meanRMS  "
               "genN  ms/mol   (rates over comparable/scored mols)\n"
            << "  ours=composite-ETKDG sampler; ours-tlib=torlib+backstop; "
               "+lib=warm fragment cache (timing only)\n";
  PlatAcc allOurs, allTlib, allCdp, allBase[8];
  double allOursWarmMs = 0.0, allTlibWarmMs = 0.0;
  auto accum = [](PlatAcc &dst, const PlatAcc &s) {
    dst.scored += s.scored;
    dst.failed += s.failed;
    dst.lt1 += s.lt1;
    dst.lt2 += s.lt2;
    dst.sumRms += s.sumRms;
    dst.genN += s.genN;
    dst.ms += s.ms;
    dst.rmss.insert(dst.rmss.end(), s.rmss.begin(), s.rmss.end());
  };
  auto emitBucket = [&](const char *label, Bucket &b) {
    std::printf("\n=== %s rotatable bonds (%u mols) ===\n", label, b.nMols);
    std::printf(
        "  method       nMols fail   %%<1.0   %%<2.0   medRMS meanRMS   genN   ms/mol\n");
    printPlatRow("ours", b.ours, b.nMols);
    if (lib) {
      printPlatRow("ours+lib", b.ours, b.nMols,
                   b.nMols ? b.oursWarmMs / b.nMols : 0.0);
    }
    if (doTlib) {
      printPlatRow("ours-tlib", b.tlib, b.nMols);
      if (lib) {
        printPlatRow("ours-tlib+lib", b.tlib, b.nMols,
                     b.nMols ? b.tlibWarmMs / b.nMols : 0.0);
      }
    }
    if (doBase) {
      for (int k = 0; k < 8; ++k) {
        printPlatRow(kBaseLabels[k], b.base[k], b.nMols);
      }
    }
      };
  for (auto &b : buckets) {
    emitBucket(b.label, b);
    accum(allOurs, b.ours);
    accum(allTlib, b.tlib);
    for (int k = 0; k < 8; ++k) {
      accum(allBase[k], b.base[k]);
    }
    allOursWarmMs += b.oursWarmMs;
    allTlibWarmMs += b.tlibWarmMs;
  }
  std::printf("\n=== OVERALL (%u mols, %zu skipped, %.0fs) ===\n", allMols,
              skipped, (nowMs() - runStart) / 1000.0);
  std::printf(
      "  method       nMols fail   %%<1.0   %%<2.0   medRMS meanRMS   genN   ms/mol\n");
  printPlatRow("ours", allOurs, allMols);
  if (lib) {
    printPlatRow("ours+lib", allOurs, allMols,
                 allMols ? allOursWarmMs / allMols : 0.0);
  }
  if (doTlib) {
    printPlatRow("ours-tlib", allTlib, allMols);
    if (lib) {
      printPlatRow("ours-tlib+lib", allTlib, allMols,
                   allMols ? allTlibWarmMs / allMols : 0.0);
    }
  }
  if (doBase) {
    for (int k = 0; k < 8; ++k) {
      printPlatRow(kBaseLabels[k], allBase[k], allMols);
    }
  }
    if (lib && std::getenv("PLATINUM_FRAGLIB_GEOM")) {
    // Coverage: the loaded lib grows by one master per CACHE MISS (fragment not
    // not in the set -> embedded fresh).  final == initial means every
    // request hit the cache.
    std::cerr << "platinum: FRAGLIB_GEOM coverage -- lib grew to "
              << lib->size()
              << " (started 195; growth = misses embedded from our ETKDG)\n";
  }
}

//! Embed one crystal fragment under a given mode (index-preserving), align its
//! heavy atoms onto the crystal fragment, and report the shape RMSD plus the
//! angle error of each exit vector vs the crystal's proper exit direction.
//! `exits` are (dummyIdx, nbrIdx) pairs in the crystal `piece` (dummies =
//! atomic number 0).
struct FragDiag {
  double shapeRms = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> exitAngles;  // degrees (best conf)
  // distribution over the pool of the per-conformer WORST exit-angle error:
  double minExit = std::numeric_limits<double>::quiet_NaN();  // best conf
  double medExit = std::numeric_limits<double>::quiet_NaN();  // typical conf
  double frac20 =
      std::numeric_limits<double>::quiet_NaN();  // fraction < 20 deg
  unsigned int nConf = 0;
  double minShapeRms =
      std::numeric_limits<double>::quiet_NaN();  // best shape in pool
};
//! Exit cap used when embedding the isolated fragment.
enum class ExitCap {
  H,
  C,
  Cyclopropyl
};
//! Embed a POOL of `nPool` conformers under the mode, MMFF-min each (if mmff),
//! and return the diag of the conformer that BEST matches the crystal exits
//! (min sum of exit-angle errors).  nPool=1 = single draw; nPool large = "does
//! the pool cover the bioactive pucker/N-inversion state" oracle.
FragDiag embedFragDiag(
    const ROMol &piece, bool etkdg, bool mmff, int seed,
    const std::vector<std::pair<unsigned int, unsigned int>> &exits,
    unsigned int nPool = 1, ExitCap cap = ExitCap::H,
    const Conformer *refConf = nullptr, bool invN = false) {
  FragDiag best;
  auto f = boost::make_shared<RWMol>(piece);
  for (const auto &e : exits) {  // cap each exit dummy (H / C / cyclopropyl)
    auto *a = f->getAtomWithIdx(e.first);
    a->setIsotope(0);
    if (cap == ExitCap::H) {
      a->setAtomicNum(1);
      a->setNoImplicit(true);
    } else {
      a->setAtomicNum(6);  // C (methyl after addHs)
      a->setNoImplicit(false);
      if (cap == ExitCap::Cyclopropyl) {  // exit-C + 2 more C -> cyclopropane
        unsigned int ca = f->addAtom(new Atom(6), true, true);
        unsigned int cb = f->addAtom(new Atom(6), true, true);
        f->addBond(e.first, ca, Bond::SINGLE);
        f->addBond(e.first, cb, Bond::SINGLE);
        f->addBond(ca, cb, Bond::SINGLE);
      }
    }
  }
  f->clearConformers();
  try {
    MolOps::sanitizeMol(*f);
  } catch (...) {
    return best;
  }
  auto emb = boost::make_shared<RWMol>(*f);
  MolOps::addHs(*emb);  // in-place; heavy indices preserved, Hs appended
  DGeomHelpers::EmbedParameters ps =
      etkdg ? DGeomHelpers::EmbedParameters(DGeomHelpers::ETKDGv3)
            : DGeomHelpers::EmbedParameters();
  ps.randomSeed = seed;
  RDKit::INT_VECT cids;
  {
    RDLog::LogStateSetter blk;
    try {
      DGeomHelpers::EmbedMultipleConfs(*emb, cids, nPool, ps);
    } catch (...) {
    }
  }
  if (cids.empty()) {
    return best;
  }
  std::unique_ptr<MMFF::MMFFMolProperties> props;
  if (mmff) {
    props.reset(new MMFF::MMFFMolProperties(*emb));
    if (!props->isValid()) {
      props.reset();
    } else {
      props->setMMFFEleTerm(false);
    }
  }
  MatchVectType
      hmap;  // ring atoms only (piece dummies are atomicNum 0 -> excluded)
  for (const auto *a : piece.atoms()) {
    if (a->getAtomicNum() > 1) {
      hmap.emplace_back(static_cast<int>(a->getIdx()),
                        static_cast<int>(a->getIdx()));
    }
  }
  // reference conformer: refConf (e.g. a good final-product fragment) if given,
  // else the piece's own conformer (the crystal).
  RWMol ref(piece);
  if (refConf) {
    ref.clearConformers();
    auto *rc = new Conformer(*refConf);
    rc->setId(0);
    ref.addConformer(rc, false);
  }
  const int refCid = ref.getConformer().getId();
  const Conformer &cc = ref.getConformer();
  std::vector<double> maxErrs;  // per-conf WORST exit-angle error
  double bestScore = std::numeric_limits<double>::max();
  double minShape =
      std::numeric_limits<double>::max();  // best (lowest) shape RMS in pool
  for (const auto cid : cids) {
    if (props) {
      try {
        std::unique_ptr<ForceFields::ForceField> ff(
            MMFF::constructForceField(*emb, props.get(), 1.0e8, cid));
        if (ff) {
          ff->initialize();
          ff->minimize(1000);
        }
      } catch (...) {
      }
    }
    FragDiag d;
    try {
      d.shapeRms = MolAlign::alignMol(*emb, ref, cid, refCid, &hmap);
    } catch (...) {
      continue;
    }
    minShape = std::min(minShape, d.shapeRms);
    const Conformer &ec = emb->getConformer(cid);
    double worst = 0.0;
    auto angDeg = [](RDGeom::Point3D a, RDGeom::Point3D b) {
      a.normalize();
      b.normalize();
      return std::acos(std::max(-1.0, std::min(1.0, a.dotProduct(b)))) * 180.0 /
             M_PI;
    };
    for (const auto &e : exits) {
      RDGeom::Point3D cDir = cc.getAtomPos(e.first) - cc.getAtomPos(e.second);
      RDGeom::Point3D oDir = ec.getAtomPos(e.first) - ec.getAtomPos(e.second);
      if (cDir.lengthSq() < 1e-9 || oDir.lengthSq() < 1e-9) {
        continue;
      }
      double ang = angDeg(cDir, oDir);
      // Enumerate the OTHER nitrogen-inversion state: reflect the exit through
      // the plane of the ring-N's two other neighbours (the umbrella flip) and
      // take the better-matching face.
      if (invN) {
        const Atom *na = piece.getAtomWithIdx(e.second);
        if (na->getAtomicNum() == 7) {
          std::vector<unsigned int> oth;
          for (const auto nb : piece.atomNeighbors(na)) {
            if (nb->getIdx() != e.first && nb->getAtomicNum() > 1) {
              oth.push_back(nb->getIdx());
            }
          }
          if (oth.size() == 2) {
            RDGeom::Point3D N = ec.getAtomPos(e.second);
            RDGeom::Point3D P = ec.getAtomPos(e.first);
            RDGeom::Point3D nrm = (ec.getAtomPos(oth[0]) - N)
                                      .crossProduct(ec.getAtomPos(oth[1]) - N);
            if (nrm.lengthSq() > 1e-9) {
              nrm.normalize();
              RDGeom::Point3D Pr = P - nrm * (2.0 * (P - N).dotProduct(nrm));
              RDGeom::Point3D o2 = Pr - N;
              if (o2.lengthSq() > 1e-9) {
                ang = std::min(ang, angDeg(cDir, o2));
              }
            }
          }
        }
      }
      d.exitAngles.push_back(ang);
      worst = std::max(worst, ang);
    }
    maxErrs.push_back(worst);
    if (worst < bestScore) {
      bestScore = worst;
      best = d;
    }
  }
  if (!maxErrs.empty()) {
    std::sort(maxErrs.begin(), maxErrs.end());
    best.minExit = maxErrs.front();
    best.medExit = maxErrs[maxErrs.size() / 2];
    unsigned int below = 0;
    for (double v : maxErrs) {
      if (v < 20.0) ++below;
    }
    best.frac20 = static_cast<double>(below) / maxErrs.size();
    best.nConf = static_cast<unsigned int>(maxErrs.size());
  }
  if (minShape != std::numeric_limits<double>::max())
    best.minShapeRms = minShape;
  return best;
}

//! Fragment a conformer (heavy + its coords) and return each rigid fragment's
//! canonical SMILES + ele-off MMFF single-point energy (exits capped as H, Hs
//! added at ideal positions).  Used to compare per-fragment strain between our
//! output and the crystal.
std::vector<std::pair<std::string, double>> fragEnergies(
    const ROMol &heavyConf) {
  std::vector<std::pair<std::string, double>> out;
  std::vector<unsigned int> links;
  try {
    links = FragmentConfGen::findLinkBonds(heavyConf);
  } catch (...) {
    return out;
  }
  std::unique_ptr<ROMol> fragged;
  if (links.empty()) {
    fragged.reset(new ROMol(heavyConf));
  } else {
    std::vector<std::pair<unsigned int, unsigned int>> lab;
    for (size_t i = 0; i < links.size(); ++i) {
      unsigned int L = 1000 + static_cast<unsigned int>(i);
      lab.emplace_back(L, L);
    }
    try {
      fragged.reset(
          MolFragmenter::fragmentOnBonds(heavyConf, links, true, &lab));
    } catch (...) {
      return out;
    }
  }
  std::vector<ROMOL_SPTR> pieces = MolOps::getMolFrags(*fragged, true);
  for (const auto &piece : pieces) {
    auto f = boost::static_pointer_cast<RWMol>(
        piece);  // fresh frag; mutate in place
    for (auto *a : f->atoms()) {
      if (a->getAtomicNum() == 0) {
        a->setAtomicNum(1);
        a->setIsotope(0);
        a->setNoImplicit(true);
      }
    }
    try {
      MolOps::sanitizeMol(*f);
      MolOps::addHs(*f, false, true);  // add remaining Hs at ideal positions
    } catch (...) {
      continue;
    }
    double e = std::numeric_limits<double>::quiet_NaN();
    try {
      MMFF::MMFFMolProperties props(*f);
      if (props.isValid()) {
        props.setMMFFEleTerm(false);
        std::unique_ptr<ForceFields::ForceField> ff(
            MMFF::constructForceField(*f, &props, 1.0e8, 0));
        if (ff) {
          ff->initialize();
          e = ff->calcEnergy();
        }
      }
    } catch (...) {
    }
    std::string s;
    try {
      RWMol key(*piece);
      for (auto *a : key.atoms()) {
        if (a->getAtomicNum() == 0)
          a->setIsotope(0);  // [*] so it matches across mols
      }
      s = MolToSmiles(key);
    } catch (...) {
    }
    out.emplace_back(s, e);
  }
  return out;
}

//! Case study on ONE ligand (PLATINUM_NAME): whole-molecule RMSD,
//! and per-fragment exit-vector fidelity + shape RMSD across embedding modes
//! (ETKDG / DG / ETKDG+MMFF / DG+MMFF).  Dumps crystal/ours to
//! PLATINUM_SDFOUT if set.
void runCaseStudy() {
  const char *nm = std::getenv("PLATINUM_NAME");
  if (!nm) {
    std::cerr << "casestudy: set PLATINUM_NAME=<ligand id from `worst`>\n";
    return;
  }
  const std::string want = nm;
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;

  std::unique_ptr<ROMol> raw;
  {
    SDMolSupplier sup(sdf, true, false);
    RDLog::LogStateSetter blk;
    while (!sup.atEnd()) {
      std::unique_ptr<ROMol> m;
      try {
        m.reset(sup.next());
      } catch (...) {
        continue;
      }
      if (m && m->hasProp("_Name") &&
          m->getProp<std::string>("_Name") == want) {
        raw = std::move(m);
        break;
      }
    }
  }
  if (!raw) {
    std::cerr << "casestudy: ligand " << want << " not found\n";
    return;
  }
  MolOps::assignStereochemistryFrom3D(*raw);
  std::unique_ptr<ROMol> h(MolOps::removeHs(*raw));
  ROMOL_SPTR heavy = largestFrag(*h);
  const std::string smi = MolToSmiles(*heavy);
  unsigned int nrot = Descriptors::calcNumRotatableBonds(*heavy);

  auto in = boost::make_shared<RWMol>(*heavy);
  in->clearConformers();
  FragmentConfGenParams op = bestDefaultParams(42);
  if (const char *e = std::getenv("PLATINUM_NCONF"))
    op.embedding.numConfsPerFragment = std::atoi(e);
  if (const char *e = std::getenv("PLATINUM_POOL")) {
    op.embedding.setFlatPool(
        std::atoi(e));  // diversity-only unless PLATINUM_FRAGMIN set
  }
  if (std::getenv("PLATINUM_FRAGMIN"))
    op.embedding.minimizeMode = FragmentMinimize::Full;
  if (std::getenv("PLATINUM_ACYCLIC")) op.wholeAcyclicFragments = true;
  if (const char *e = std::getenv("PLATINUM_OUTCONFS"))
    op.numOutputConfs = std::atoi(e);
  std::printf("  [ours cfg: nConf=%u out=%u fragMin=%d]\n",
              op.embedding.numConfsPerFragment, op.numOutputConfs,
              (int)op.embedding.minimizeMode);
  requireValidParams(op, "runCaseStudy");
  auto res = FragmentConfGen(op).build(*in);
  double oursRms = minHeavyRms(res.conformers, *heavy);
    std::printf("\n[casestudy] %s   nrot=%u\n  %s\n", want.c_str(), nrot,
              smi.c_str());
  std::printf("  ours  minRMS=%.2f (genN=%zu)\n", oursRms,
              res.conformers.size());

  // Per-fragment exit-vector fidelity across embedding modes.
  auto links = FragmentConfGen::findLinkBonds(*heavy);
  std::vector<std::pair<unsigned int, unsigned int>> lab;
  for (size_t i = 0; i < links.size(); ++i) {
    unsigned int L = 1000 + static_cast<unsigned int>(i);
    lab.emplace_back(L, L);
  }
  std::unique_ptr<ROMol> fragged(
      MolFragmenter::fragmentOnBonds(*heavy, links, true, &lab));
  std::vector<ROMOL_SPTR> pieces = MolOps::getMolFrags(*fragged, true);

  struct Mode {
    const char *label;
    bool etkdg, mmff;
    unsigned int nPool;
  };
  const std::vector<Mode> modes = {
      {"ETKDG", true, false, 1},
      {"DG", false, false, 1},
      {"ETKDG+MMFF", true, true, 1},
      {"DG+MMFF", false, true, 1},
      {"ETKDG+MMFF x1000(best)", true, true, 1000}};
  std::printf(
      "\n  per-fragment: shape RMSD (heavy) and exit-vector angle error(s) "
      "vs crystal\n  (x1000 = best-of-1000-pool = does overkill cover the "
      "crystal state?)\n");
  for (size_t pi = 0; pi < pieces.size(); ++pi) {
    const ROMol &piece = *pieces[pi];
    unsigned int nHeavy = 0, nExit = 0;
    std::vector<std::pair<unsigned int, unsigned int>> exits;
    for (const auto *a : piece.atoms()) {
      if (a->getAtomicNum() > 1) ++nHeavy;
      if (a->getAtomicNum() == 0) {
        ++nExit;
        const Atom *nbr = nullptr;
        for (const auto nb : piece.atomNeighbors(a)) {
          nbr = nb;
          break;
        }
        if (nbr) exits.emplace_back(a->getIdx(), nbr->getIdx());
      }
    }
    std::string psmi;
    try {
      psmi = MolToSmiles(piece);
    } catch (...) {
    }
    std::printf("  frag %zu: %u heavy, %u exit(s)  %s\n", pi, nHeavy, nExit,
                psmi.c_str());
    for (const auto &m : modes) {
      FragDiag d = embedFragDiag(piece, m.etkdg, m.mmff, 42, exits, m.nPool);
      std::printf("     %-22s shapeRMS=%5.2f (poolMin=%5.2f)  exitErr=",
                  m.label, d.shapeRms, d.minShapeRms);
      if (d.exitAngles.empty()) {
        std::printf("(none)");
      } else {
        for (double a : d.exitAngles) std::printf("%5.1f ", a);
      }
      std::printf("deg\n");
    }
  }

  // --- per-fragment MMFF strain: crystal vs our best output ---
  auto bestHeavy = [&](const std::vector<ROMOL_SPTR> &confs) -> ROMOL_SPTR {
    double best = std::numeric_limits<double>::max();
    ROMOL_SPTR bc;
    RWMol ref(*heavy);
    for (const auto &c : confs) {
      if (!c || c->getNumConformers() == 0) continue;
      try {
        ROMOL_SPTR p(MolOps::removeHs(*c));
        if (p->getNumAtoms() != heavy->getNumAtoms()) continue;
        RWMol probe(*p);
        double r = MolAlign::getBestRMS(probe, ref);
        if (r < best) {
          best = r;
          bc = p;
        }
      } catch (...) {
      }
    }
    return bc;
  };
  ROMOL_SPTR oursBest = bestHeavy(res.conformers);
  auto emap = [](const std::vector<std::pair<std::string, double>> &v) {
    std::map<std::string, double> m;
    for (const auto &p : v) m[p.first] = p.second;
    return m;
  };
  auto mCry = emap(fragEnergies(*heavy));
  auto mOur = oursBest ? emap(fragEnergies(*oursBest))
                       : std::map<std::string, double>{};
  std::printf(
      "\n  per-fragment ele-off MMFF strain of the best-RMS output conformer "
      "(kcal/mol):\n    %-38s %8s %8s %8s\n",
      "fragment", "crystal", "ours");
  for (const auto &kv : mCry) {
    double o = mOur.count(kv.first) ? mOur[kv.first]
                                    : std::numeric_limits<double>::quiet_NaN();
    std::string lab =
        kv.first.size() > 38 ? kv.first.substr(0, 35) + "..." : kv.first;
    std::printf("    %-38s %8.1f %8.1f\n", lab.c_str(), kv.second, o);
  }

  if (const char *out = std::getenv("PLATINUM_SDFOUT")) {
    try {
      SDWriter w(out);
      auto tag = [&](ROMol &m, const std::string &t) {
        m.setProp("_Name", want + "_" + t);
        w.write(m);
      };
      RWMol cry(*heavy);
      tag(cry, "crystal");
      if (!res.conformers.empty()) {
        RWMol o(*res.conformers.front());
        tag(o, "ours");
      }
      w.close();
      std::printf("\n  wrote crystal/ours -> %s\n", out);
    } catch (...) {
    }
  }
}

//! Find low-rotor (<=2) Platinum ligands where OURS reproduces the crystal
//! poorly -- i.e. our geometry is the problem, not an
//! intrinsically hard molecule.  Prints them worst-ours first with names for
//! the exit-vector deep dive (`casestudy`, set PLATINUM_NAME).
void runWorst() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  unsigned int maxRot = 2;
  if (const char *e = std::getenv("PLATINUM_MAXROT")) maxRot = std::atoi(e);
  unsigned int nThreads = std::thread::hardware_concurrency();
  if (const char *e = std::getenv("PLATINUM_THREADS")) nThreads = std::atoi(e);

  struct Job {
    ROMOL_SPTR heavy;
    std::string name, smi;
    unsigned int nrot;
  };
  std::vector<Job> jobs;
  {
    std::unique_ptr<SDMolSupplier> sup;
    try {
      sup = std::make_unique<SDMolSupplier>(sdf, true, false);
    } catch (...) {
      std::cerr << "worst: cannot open\n";
      return;
    }
    RDLog::LogStateSetter blk;
    while (!sup->atEnd()) {
      std::unique_ptr<ROMol> raw;
      try {
        raw.reset(sup->next());
      } catch (...) {
        continue;
      }
      if (!raw || raw->getNumConformers() == 0) continue;
      std::string name =
          raw->hasProp("_Name") ? raw->getProp<std::string>("_Name") : "?";
      try {
        MolOps::assignStereochemistryFrom3D(*raw);
        std::unique_ptr<ROMol> h(MolOps::removeHs(*raw));
        ROMOL_SPTR heavy = largestFrag(*h);
        if (!heavy || heavy->getNumAtoms() < 2) continue;
        unsigned int nrot = Descriptors::calcNumRotatableBonds(*heavy);
        if (nrot > maxRot) continue;
        jobs.push_back({heavy, name, MolToSmiles(*heavy), nrot});
      } catch (...) {
      }
    }
  }
  std::cout << "[worst] " << jobs.size() << " ligands with <=" << maxRot
            << " rotatable bonds; ours vs "
            << "ours" << "\n";

  struct Rec {
    std::string name, smi;
    unsigned int nrot;
    double ours, cdp;
    size_t genN;
  };
  std::vector<Rec> recs(jobs.size());
  std::atomic<size_t> next{0};
  auto worker = [&]() {
    FragmentConfGenParams p = bestDefaultParams(42);
    for (;;) {
      size_t i = next.fetch_add(1);
      if (i >= jobs.size()) break;
      auto in = boost::make_shared<RWMol>(*jobs[i].heavy);
      in->clearConformers();
      FragmentConfGenResult res;
      try {
        requireValidParams(p, "runWorst");
        res = FragmentConfGen(p).build(*in);
      } catch (...) {
      }
      double ours = minHeavyRms(res.conformers, *jobs[i].heavy);
      double cdp = std::numeric_limits<double>::quiet_NaN();
            recs[i] = {jobs[i].name, jobs[i].smi, jobs[i].nrot,
                 ours,         cdp,         res.conformers.size()};
    }
  };
  std::vector<std::thread> pool;
  for (unsigned int t = 0; t < std::max(1u, nThreads); ++t)
    pool.emplace_back(worker);
  for (auto &t : pool) t.join();

  std::sort(recs.begin(), recs.end(), [](const Rec &a, const Rec &b) {
    double ra = std::isfinite(a.ours) ? a.ours : -1,
           rb = std::isfinite(b.ours) ? b.ours : -1;
    return ra > rb;
  });
  std::printf("\n  %-16s rot  genN  oursRMS  cdpRMS   smiles\n", "name");
  for (size_t i = 0; i < recs.size() && i < 30; ++i) {
    const auto &r = recs[i];
    std::printf("  %-16s %3u  %4zu  %7.2f  %6.2f   %s\n", r.name.c_str(),
                r.nrot, r.genN, r.ours, r.cdp, r.smi.c_str());
  }
}

//! Sweep junction-angle samplers (ETKDG / TorsionLib / Uniform, with/without a
//! uniform-grid backstop) through the CoordsOnly joiner on the Platinum set (50
//! by default).  Env: PLATINUM_MAXMOLS, PLATINUM_THREADS, ASM_NCONF.
void runAsmSampler() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = 50;
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);
  unsigned int threads = 8;
  if (const char *e = std::getenv("PLATINUM_THREADS")) threads = std::atoi(e);
  unsigned int nConf = 16;
  if (const char *e = std::getenv("ASM_NCONF")) nConf = std::atoi(e);
  unsigned int thompson =
      0;  // >0 = FIXED budget; unset -> auto-scaled (default)
  if (const char *e = std::getenv("ASM_THOMPSON")) thompson = std::atoi(e);
  const bool exhaustive = std::getenv("ASM_EXHAUSTIVE") != nullptr;
  // fragment embedding: "etkdg" (ETKDG, no MMFF-min) or "dgmmff" (DG +
  // MMFF-min). default (unset) = current shipped behaviour (ETKDG + MMFF-min).
  bool fragUseDG = false, fragMin = true;
  std::string fragEmbed = "etkdg+min";
  if (const char *e = std::getenv("ASM_FRAGEMBED")) {
    fragEmbed = e;
    if (fragEmbed == "etkdg") {
      fragUseDG = false;
      fragMin = false;
    } else if (fragEmbed == "dgmmff") {
      fragUseDG = true;
      fragMin = true;
    }
  }
  int seed = 0xf00d;
  if (const char *e = std::getenv("ASM_SEED"))
    seed = std::atoi(e);  // match platinum's 42
  // auto-budget coefficients (tunable): budget = clamp(perRotor*nRot +
  // perFrag*fragArms)
  int tsPerRotor = -1, tsPerFrag = -1, tsMin = -1, tsMax = -1;
  if (const char *e = std::getenv("ASM_TS_PERROTOR")) tsPerRotor = std::atoi(e);
  if (const char *e = std::getenv("ASM_TS_PERFRAG")) tsPerFrag = std::atoi(e);
  if (const char *e = std::getenv("ASM_TS_MIN")) tsMin = std::atoi(e);
  if (const char *e = std::getenv("ASM_TS_MAX")) tsMax = std::atoi(e);

  std::vector<ROMOL_SPTR> jobs;
  {
    std::unique_ptr<SDMolSupplier> sup;
    try {
      sup = std::make_unique<SDMolSupplier>(sdf, true, false);
    } catch (...) {
      std::cerr << "asmsampler: cannot open " << sdf << "\n";
      return;
    }
    RDLog::LogStateSetter blk;
    while (!sup->atEnd() && jobs.size() < maxMols) {
      std::unique_ptr<ROMol> raw;
      try {
        raw.reset(sup->next());
      } catch (...) {
        continue;
      }
      if (!raw || raw->getNumConformers() == 0) continue;
      try {
        MolOps::assignStereochemistryFrom3D(*raw);
        std::unique_ptr<ROMol> heavy(MolOps::removeHs(*raw));
        ROMOL_SPTR refHeavy = largestFrag(*heavy);
        if (refHeavy && refHeavy->getNumAtoms() >= 2) jobs.push_back(refHeavy);
      } catch (...) {
      }
    }
  }

  auto grid = [](double s) {
    return std::make_shared<UniformTorsionSampler>(s);
  };
  auto comp = [](std::shared_ptr<TorsionSampler> a,
                 std::shared_ptr<TorsionSampler> b) {
    return std::make_shared<CompositeTorsionSampler>(
        std::vector<std::shared_ptr<TorsionSampler>>{std::move(a),
                                                     std::move(b)});
  };
  std::vector<std::pair<std::string, std::shared_ptr<TorsionSampler>>> cfgs;
  cfgs.emplace_back("Uniform60", grid(60.0));
  cfgs.emplace_back("Uniform30", grid(30.0));
  cfgs.emplace_back("Uniform15", grid(15.0));
  // ASM_ETKDG_MAXANG caps the ETKDG sampler's arms/rotor (default 12; it pads
  // with a 60-deg grid up to this).  Fewer arms -> proportionally fewer
  // Thompson Beta draws (the drive cost); this sweeps whether 12 is overkill
  // for accuracy.
  unsigned int etkdgMaxAng =
      6;  // matches the ETKDGTorsionSampler class default
  if (const char *e = std::getenv("ASM_ETKDG_MAXANG")) {
    etkdgMaxAng = static_cast<unsigned int>(std::atoi(e));
  }
  cfgs.emplace_back(
      "ETKDG", std::make_shared<ETKDGTorsionSampler>(10.0, 60.0, etkdgMaxAng));
  cfgs.emplace_back("ETKDG+grid", comp(std::make_shared<ETKDGTorsionSampler>(
                                           10.0, 60.0, etkdgMaxAng),
                                       grid(60.0)));
  try {
    cfgs.emplace_back("TorLib", hamburgSampler());
    cfgs.emplace_back("TorLib+grid", comp(hamburgSampler(), grid(60.0)));
  } catch (const std::exception &e) {
    std::cerr << "asmsampler: torsion library unavailable (" << e.what()
              << "); skipping TorLib configs\n";
  }
  // ASM_ONLY=<name>[,<name>...]: keep only samplers whose EXACT name is listed
  // (e.g. ASM_ONLY=TorLib,ETKDG,Uniform30 to run just those three).
  if (const char *only = std::getenv("ASM_ONLY")) {
    std::string spec(only);
    auto listed = [&](const std::string &n) {
      size_t p = 0;
      while (p < spec.size()) {
        size_t c = spec.find(',', p);
        if (c == std::string::npos) c = spec.size();
        if (spec.compare(p, c - p, n) == 0) return true;
        p = c + 1;
      }
      return false;
    };
    std::vector<std::pair<std::string, std::shared_ptr<TorsionSampler>>> keep;
    for (auto &c : cfgs) {
      if (listed(c.first)) keep.push_back(c);
    }
    cfgs.swap(keep);
  }

  std::printf(
      "[asmsampler] junction-angle samplers through CoordsOnly joiner, %zu "
      "mols, nConf=%u, frag=%s%s%s\n\n  %-14s %-7s %-7s %-6s %-6s %-7s %-8s\n",
      jobs.size(), nConf, fragEmbed.c_str(), exhaustive ? ", EXHAUSTIVE" : "",
      thompson ? (", THOMPSON=" + std::to_string(thompson)).c_str() : "",
      "sampler", "%<1", "%<2", "med", "genN", "warmMs", "avgBudg");

  // WARM (assembly-only, cached-fragment) ms/mol -- the production timing.
  // Enable the profiler from code so the reported time excludes the cold
  // one-off embed cost.
  setJoinerProfiling(true);
  takeJoinerWarmMsPerMol(0);  // clear any warmup accumulation

  for (auto &cfg : cfgs) {
    PlatAcc acc;
    long long budgetSum = 0;
    std::mutex mtx;
    std::atomic<size_t> next{0};
    auto worker = [&]() {
      RDLog::LogStateSetter blk;
      // per-thread sampler copy: TorsionLibrarySampler lazily loads mutable
      // rules, so a shared instance races across threads
      auto localSampler =
          cfg.second ? cfg.second->copy() : std::shared_ptr<TorsionSampler>();
      while (true) {
        size_t i = next++;
        if (i >= jobs.size()) break;
        if (const char *e =
                std::getenv("ASM_MOL")) {  // probe ONE molecule by index
          if (i != static_cast<size_t>(std::atol(e))) continue;
        }
        double rms = std::numeric_limits<double>::quiet_NaN();
        size_t nc = 0;
        unsigned int bud = 0;
        const double t0 = nowMs();
        try {
          FragmentConfGenParams pp;
          pp.numOutputConfs = 2000;
          if (const char *e = std::getenv("ASM_OUTCONFS"))
            pp.numOutputConfs = std::atoi(e);
          pp.randomSeed = seed;
          pp.search.torsionSampler = localSampler;
          pp.search.thompsonBudget =
              thompson;  // 0 -> auto (unless exhaustive)
          pp.search.thompson.autoBudget = !exhaustive && thompson == 0;
          pp.search.searchMode = searchModeFromEnv();  // ASM_SEARCH override
          // ASM_BEAM: widen (or effectively disable) the exhaustive-beam
          // pruning so a uniform sweep is a TRUE ceiling -- Uniform15 (superset
          // of Uniform30's angles) must then dominate Uniform30.  A small beam
          // prunes good partial assemblies.
          if (const char *e = std::getenv("ASM_BEAM"))
            pp.search.tree.beamWidth = std::atoi(e);
          if (const char *e = std::getenv("ASM_ROOTSEEDS"))
            pp.search.rootSeeds = std::atoi(e);
          if (tsPerRotor >= 0) pp.search.thompson.perRotor = tsPerRotor;
          if (tsPerFrag >= 0) pp.search.thompson.perFragConf = tsPerFrag;
          if (tsMin >= 0) pp.search.thompson.minBudget = tsMin;
          if (tsMax >= 0) pp.search.thompson.maxBudget = tsMax;
          // crystal-structure reproduction -> MMFF94s; override w/ env
          pp.joiner.ffVariant =
              std::getenv("ASM_MMFF94") ? "MMFF94" : "MMFF94s";
          if (const char *e = std::getenv("ASM_VDWCUTOFF"))
            pp.joiner.interFragVdwCutoff = std::atof(e);
          if (const char *e = std::getenv("ASM_PRIOR"))
            pp.search.thompson.priorStrength = std::atof(e);
          if (const char *e = std::getenv("ASM_BACKSTOP"))
            pp.search.thompson.backstopStepDeg = std::atof(e);
          if (const char *e = std::getenv("ASM_SIZEPRIOR"))
            pp.search.thompson.sizePriorExp = std::atof(e);
          // respect the shipped default (on) unless ASM_SYMDEDUP is set (=0
          // disables)
          if (const char *e = std::getenv("ASM_SYMDEDUP"))
            pp.search.finalSymmetryDedup = std::atoi(e) != 0;
          if (const char *e = std::getenv("ASM_FPNOVELTY"))
            pp.search.thompson.noveltyAngleDeg = std::atof(e);
          if (const char *e = std::getenv("ASM_DIVRMS"))
            pp.search.diversityRmsThresh = std::atof(e);
          if (const char *e = std::getenv("ASM_PERCLASS"))
            pp.embedding.perClassEmbedding = std::atoi(e) != 0;
          pp.embedding.fragmentEmbedMode =
              fragUseDG ? FragmentEmbedMode::DG : FragmentEmbedMode::ETKDG;
          pp.embedding.minimizeMode =
              fragMin ? FragmentMinimize::Full : FragmentMinimize::None;
          requireValidParams(pp, "runAsmSampler");
          auto rr = FragmentConfGen(pp).build(*jobs[i]);
          nc = rr.conformers.size();
          bud = rr.joinerBudget;
          rms = minHeavyRms(rr.conformers, *jobs[i]);
          if (std::getenv("ASM_RANKDIAG")) {
            auto wr = winnerEnergyRank(rr.conformers, *jobs[i]);
            std::fprintf(stderr,
                         "[rankdiag] mol=%zu nconf=%zu winnerRank=%d "
                         "winnerRms=%.3f\n",
                         i, nc, wr.first, wr.second);
            std::fflush(stderr);
          }
        } catch (...) {
        }
        const double ms = nowMs() - t0;
        if (std::getenv("ASM_TRACE")) {
          unsigned int nrot = 0;
          try {
            nrot = Descriptors::calcNumRotatableBonds(*jobs[i]);
          } catch (...) {
          }
          std::fprintf(stderr,
                       "[asmTRACE] %-12s mol %2zu rot=%2u nAt=%2u nc=%4zu "
                       "minRMS=%.3f\n",
                       cfg.first.c_str(), i, nrot, jobs[i]->getNumAtoms(), nc,
                       rms);
          std::fflush(stderr);
        }
        std::lock_guard<std::mutex> lk(mtx);
        acc.add(rms, nc, ms);
        budgetSum += bud;
      }
    };
    std::vector<std::thread> pool;
    for (unsigned int t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto &t : pool) t.join();
    const double denom = acc.scored ? acc.scored : 1;
    if (std::getenv("ASM_PROFILE"))
      printJoinerProfile(cfg.first.c_str(), jobs.size());  // prints + resets
    const double warmMs = takeJoinerWarmMsPerMol(jobs.size());
    std::printf(
        "  %-14s %6.1f%% %6.1f%% %6.2f %6.1f %7.1f %8.0f\n", cfg.first.c_str(),
        100.0 * acc.lt1 / denom, 100.0 * acc.lt2 / denom, acc.median(),
        jobs.empty() ? 0.0 : acc.genN / jobs.size(), warmMs,
        jobs.empty() ? 0.0 : static_cast<double>(budgetSum) / jobs.size());
    std::fflush(stdout);
  }
}

//! Timing of RMSD methods (heavy atoms): raw same-index RMSD, QCP (optimal
//! superposition), symmetry-aware min-QCP over automorphisms, and RDKit's
//! symmetry-aware getBestRMS.  Also reports heavy count + automorphism count +
//! one-time automorphism build time.
void runRmsdBench() {
  using clk = std::chrono::steady_clock;
  auto nsPer = [](clk::time_point a, clk::time_point b, size_t k) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() /
           static_cast<double>(k ? k : 1);
  };
  struct Case {
    const char *name, *smi;
  };
  const std::vector<Case> cases = {
      {"asym-drug", "CC(=O)Nc1ccc(O)cc1CCN"},
      {"benzene", "c1ccccc1"},
      {"naphthalene", "c1ccc2ccccc2c1"},
      {"p-diCl-benzene", "Clc1ccc(Cl)cc1"},
      {"benzophenone", "O=C(c1ccccc1)c1ccccc1"},
      {"biphenyl", "c1ccc(-c2ccccc2)cc1"},
      {"caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"},
      {"adamantane", "C1C2CC3CC1CC(C2)C3"},
      {"steroid", "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"},
  };
  const unsigned int nconf = 24;
  RDLog::LogStateSetter blk;
  std::printf(
      "[rmsdbench] heavy-atom RMSD timing, %u confs, all-pairs (ns per comparison)\n",
      nconf);
  std::printf("  %-15s %5s %5s %6s  %8s %8s %10s %12s   %8s\n", "molecule",
              "heavy", "autos", "pairs", "rawIndex", "qcp", "qcp+symm",
              "getBestRMS", "autoBuild");
  for (const auto &cs : cases) {
    std::unique_ptr<RWMol> molH(SmilesToMol(cs.smi));
    if (!molH) continue;
    MolOps::addHs(*molH);
    DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
    ps.randomSeed = 42;
    INT_VECT cids;
    DGeomHelpers::EmbedMultipleConfs(*molH, cids, nconf, ps);
    // heavy-only molecule (keeps conformers) so every method is heavy-atom only
    RWMol mol(*molH);
    MolOps::removeHs(mol);
    const unsigned int C = mol.getNumConformers();
    if (C < 2) continue;
    const unsigned int H = mol.getNumAtoms();

    auto tb0 = clk::now();
    std::vector<unsigned int> heavy;
    auto autos = heavyAtomAutomorphisms(mol, heavy);  // heavy == all atoms here
    auto tb1 = clk::now();
    std::vector<std::vector<unsigned int>> permB;
    for (const auto &pm : autos) {
      std::vector<unsigned int> ib(H);
      for (unsigned int p = 0; p < H; ++p) ib[p] = heavy[pm[p]];
      permB.push_back(std::move(ib));
    }

    std::vector<std::vector<double>> co(C, std::vector<double>(3 * H));
    std::vector<int> id(C);
    unsigned int ci = 0;
    for (auto cit = mol.beginConformers(); cit != mol.endConformers();
         ++cit, ++ci) {
      id[ci] = (*cit)->getId();
      for (unsigned int i = 0; i < H; ++i) {
        const auto p = (*cit)->getAtomPos(i);
        co[ci][3 * i] = p.x;
        co[ci][3 * i + 1] = p.y;
        co[ci][3 * i + 2] = p.z;
      }
    }
    const size_t pairs = static_cast<size_t>(C) * (C - 1) / 2;
    volatile double sink = 0.0;

    auto t0 = clk::now();
    for (unsigned int i = 0; i < C; ++i)
      for (unsigned int j = i + 1; j < C; ++j) {
        double s = 0.0;
        for (unsigned int a = 0; a < 3 * H; ++a) {
          const double q = co[i][a] - co[j][a];
          s += q * q;
        }
        sink = sink + std::sqrt(s / H);
      }
    auto t1 = clk::now();
    for (unsigned int i = 0; i < C; ++i)
      for (unsigned int j = i + 1; j < C; ++j)
        sink = sink + qcpRmsd(co[i], co[j], heavy, heavy);
    auto t2 = clk::now();
    for (unsigned int i = 0; i < C; ++i)
      for (unsigned int j = i + 1; j < C; ++j) {
        double m = 1e18;
        for (const auto &ib : permB) {
          const double r = qcpRmsd(co[i], co[j], heavy, ib);
          if (r < m) m = r;
        }
        sink = sink + m;
      }
    auto t3 = clk::now();
    double best = 0.0;
    for (unsigned int i = 0; i < C; ++i)
      for (unsigned int j = i + 1; j < C; ++j) {
        try {
          best += MolAlign::getBestRMS(mol, mol, id[i], id[j]);
        } catch (...) {
        }
      }
    auto t4 = clk::now();
    sink = sink + best;

    std::printf("  %-15s %5u %5zu %6zu  %8.0f %8.0f %10.0f %12.0f   %6.1fus\n",
                cs.name, H, autos.size(), pairs, nsPer(t0, t1, pairs),
                nsPer(t1, t2, pairs), nsPer(t2, t3, pairs),
                nsPer(t3, t4, pairs), nsPer(tb0, tb1, 1) / 1000.0);
    (void)sink;
  }
}

// ---- matrix: embeddings x torsion-profiles x search cells -------------------
// Full factorial performance matrix over the three embedding modes, three
// torsion profiles and three search cells, at fixed MMFF94s + carbon caps +
// default budget/ ring-state counts.  Each cell runs the whole set TWICE
// against a per-cell in-memory shared library: pass 1 fills the cache
// (cold/first-run) and pass 2 reuses it (warm), giving both ms/mol figures.
// Single-thread by default for clean timing (set PLATINUM_THREADS>1 only for a
// quick, contention-inflated accuracy pass).  Rows are printed as each cell
// finishes so a long (overnight) run yields partial results.
void runMatrix() {
  std::string sdf =
      std::string(std::getenv("RDBASE") ? std::getenv("RDBASE") : ".") +
      "/Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf";
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = 0;
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);

  struct Job {
    ROMOL_SPTR refHeavy;  // heavy-atom crystal pose (RMSD target)
    ROMOL_SPTR input;     // graph only (generator input)
    std::string smi;
  };
  std::vector<Job> jobs;
  {
    std::unique_ptr<SDMolSupplier> supplier;
    try {
      supplier = std::make_unique<SDMolSupplier>(sdf, true, false);
    } catch (...) {
      std::cerr << "matrix: cannot open " << sdf << "\n";
      return;
    }
    RDLog::LogStateSetter blocker;
    while (!supplier->atEnd()) {
      if (maxMols && jobs.size() >= maxMols) break;
      std::unique_ptr<ROMol> raw;
      try {
        raw.reset(supplier->next());
      } catch (...) {
        continue;
      }
      if (!raw || raw->getNumConformers() == 0) continue;
      try {
        MolOps::assignStereochemistryFrom3D(*raw);
        std::unique_ptr<ROMol> heavy(MolOps::removeHs(*raw));
        ROMOL_SPTR refHeavy = largestFrag(*heavy);
        if (!refHeavy || refHeavy->getNumAtoms() < 2) continue;
        auto input = boost::make_shared<RWMol>(*refHeavy);
        input->clearConformers();
        std::string smi;
                jobs.push_back({refHeavy, input, smi});
      } catch (...) {
      }
    }
  }

  struct Embed {
    const char *name;
    FragmentEmbedMode mode;
    bool minimize;
  };
  const Embed embeds[] = {{"ETKDG", FragmentEmbedMode::ETKDG, false},
                          {"ETKDG+MMFF", FragmentEmbedMode::ETKDG, true},
                          {"DG+MMFF", FragmentEmbedMode::DG, true}};
  const char *torsionNames[] = {"torlib", "uniform30", "etkdg"};
  struct Search {
    const char *name;
    bool thompson;
    double fpDeg;  // Thompson novelty angle (0 = coord-RMSD); ignored for
                   // rotortree
  };
  const Search searches[] = {{"thompson+fp", true, 30.0},
                             {"thompson+rmsd", true, 0.0},
                             {"rotortree", false, 0.0}};

  // Resume: skip any cell whose "embed torsion search" row already appears in
  // the file named by MATRIX_DONE (a partial run's output).  Point MATRIX_DONE
  // at the partial table and redirect stdout to a new file; concatenate the two
  // for the full 27-cell matrix.
  std::set<std::string> doneCells;
  if (const char *e = std::getenv("MATRIX_DONE")) {
    std::ifstream f(e);
    std::string line;
    while (std::getline(f, line)) {
      std::istringstream iss(line);
      std::string a, b, c;
      if (iss >> a >> b >> c) doneCells.insert(a + "|" + b + "|" + c);
    }
    std::fprintf(stderr, "[matrix] resume: %zu cells already done in %s\n",
                 doneCells.size(), e);
  }

  const int seed = 0xf00d;
  std::printf(
      "[matrix] %zu mols, single-thread, MMFF94s, carbon caps, "
      "numConfsPerFragment=6, numOutputConfs=50\n",
      jobs.size());
  std::printf("  %-11s %-10s %-14s %6s %6s %7s %9s %9s\n", "embed", "torsion",
              "search", "%<1", "%<2", "genN", "msCold", "msWarm");
  std::fflush(stdout);

  for (const auto &em : embeds) {
    for (const char *tn : torsionNames) {
      for (const auto &se : searches) {
        const std::string cellKey =
            std::string(em.name) + "|" + tn + "|" + se.name;
        if (doneCells.count(cellKey)) continue;  // resume: already computed
        FragmentConfGenParams base;
        base.embedding.numConfsPerFragment = 6;
        base.numOutputConfs = 50;
        base.energyWindow = 10.0;
        base.randomSeed = seed;
        base.joiner.ffVariant = "MMFF94s";
        base.embedding.fragmentEmbedMode = em.mode;
        base.embedding.minimizeMode =
            em.minimize ? FragmentMinimize::Full : FragmentMinimize::None;
        if (std::string(tn) == "torlib") {
          base.search.torsionSampler = hamburgSampler();
        } else if (std::string(tn) == "uniform30") {
          base.search.torsionSampler =
              std::make_shared<UniformTorsionSampler>(30.0);
        } else {
          base.search.torsionSampler = std::make_shared<ETKDGTorsionSampler>();
        }
        if (se.thompson) {
          base.search.thompson.autoBudget = true;
          base.search.thompson.noveltyAngleDeg = se.fpDeg;
        } else {
          base.search.thompson.autoBudget = false;
          base.search.thompsonBudget = 0;  // -> exhaustive tree beam
        }

        // Per-cell in-memory shared library, embedding-matched to `base` so the
        // load guard is satisfied and pass 2 is genuinely warm.
        FraglibParams flp;
        flp.numConfsPerFragment = base.embedding.numConfsPerFragment;
        flp.fragmentEmbedMode = base.embedding.fragmentEmbedMode;
        flp.randomSeed = seed;
        flp.minimizeMode = base.embedding.minimizeMode;
        flp.classParams = base.embedding.classParams;
        base.fraglib = std::make_shared<Fraglib>(flp);

        double coldMs = 0.0, warmMs = 0.0, genN = 0.0;
        size_t lt1 = 0, lt2 = 0, scored = 0;
        // pass 1: cold (library fills as fragments are first seen)
        for (const auto &j : jobs) {
          requireValidParams(base, "runMatrix");
          const double t0 = nowMs();
          FragmentConfGen(base).build(*j.input);
          coldMs += nowMs() - t0;
        }
        // pass 2: warm (every fragment now cached)
        for (const auto &j : jobs) {
          requireValidParams(base, "runMatrix");
          const double t0 = nowMs();
          auto res = FragmentConfGen(base).build(*j.input);
          warmMs += nowMs() - t0;
          const double rms = minHeavyRms(res.conformers, *j.refHeavy);
          genN += static_cast<double>(res.conformers.size());
          if (std::isfinite(rms)) {
            ++scored;
            if (rms < 1.0) ++lt1;
            if (rms < 2.0) ++lt2;
          }
        }
        const double denom = scored ? static_cast<double>(scored) : 1.0;
        const double nm = jobs.empty() ? 1.0 : static_cast<double>(jobs.size());
        std::printf("  %-11s %-10s %-14s %5.1f%% %5.1f%% %7.1f %9.0f %9.0f\n",
                    em.name, tn, se.name, 100.0 * lt1 / denom,
                    100.0 * lt2 / denom, genN / nm, coldMs / nm, warmMs / nm);
        std::fflush(stdout);
      }
    }
  }

    // Timing is generation only -- process startup and I/O are
  // not a cell of our factorial.  Timing is the bridge's internal gen_ms
  // (prepare+generate only) -- Python startup, imports, and SDF I/O are
  // excluded, and any built-in fragment library loads at
  // ConformerGenerator construction (before that timer), so gen_ms is a warm,
  // startup-free per-molecule figure.  Reported in the warm column.
  }

//! Build a submol of `parent` over parent-atom indices `keep`, copying the
//! crystal conformer.  parentToSub maps parent atom idx -> submol atom idx.
boost::shared_ptr<RWMol> submolWithCoords(const ROMol &parent,
                                          const std::vector<int> &keep,
                                          std::map<int, int> &parentToSub) {
  auto m = boost::make_shared<RWMol>();
  parentToSub.clear();
  for (int pi : keep) {
    parentToSub[pi] =
        m->addAtom(new Atom(*parent.getAtomWithIdx(pi)), false, true);
  }
  for (const auto b : parent.bonds()) {
    int a = b->getBeginAtomIdx(), c = b->getEndAtomIdx();
    if (parentToSub.count(a) && parentToSub.count(c)) {
      m->addBond(parentToSub[a], parentToSub[c], b->getBondType());
    }
  }
  auto *conf = new Conformer(m->getNumAtoms());
  const Conformer &pc = parent.getConformer();
  for (auto &kv : parentToSub) {
    conf->setAtomPos(kv.second, pc.getAtomPos(kv.first));
  }
  m->addConformer(conf, false);
  return m;
}

//! Embed `probe` fresh (ETKDG, single draw unless nPool>1), optionally MMFF-min
//! each, align the `probeCoreSub` atoms onto `refCore` (built over the same
//! core atoms in `refCoreIdx` order) and return the best heavy-atom core RMSD
//! vs the crystal. One aromatic exit to score: origin core atom (probe
//! sub-idx), the neighbour that defines the exit direction (probe sub-idx: real
//! aromatic atom in context, carbon cap in isolated), and the crystal's true
//! unit exit vector.
struct ExitProbe {
  int coreSub;
  int nbrSub;
  RDGeom::Point3D crystalVec;
};
double embedAlignCore(RWMol probe, const ROMol &refCore,
                      const std::vector<int> &probeCoreSub,
                      const std::vector<int> &refCoreIdx,
                      const std::vector<ExitProbe> &exits, bool mmff, int seed,
                      unsigned int nPool, double *exitErrOut) {
  if (exitErrOut) *exitErrOut = std::numeric_limits<double>::quiet_NaN();
  probe.clearConformers();
  try {
    MolOps::sanitizeMol(probe);
  } catch (...) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  MolOps::addHs(probe);
  DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
  ps.randomSeed = seed;
  RDKit::INT_VECT cids;
  {
    RDLog::LogStateSetter blk;
    try {
      DGeomHelpers::EmbedMultipleConfs(probe, cids, nPool, ps);
    } catch (...) {
    }
  }
  if (cids.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::unique_ptr<MMFF::MMFFMolProperties> props;
  if (mmff) {
    props.reset(new MMFF::MMFFMolProperties(probe));
    if (!props->isValid()) {
      props.reset();
    } else {
      props->setMMFFEleTerm(false);
    }
  }
  MatchVectType amap;  // (probeIdx, refIdx)
  for (size_t i = 0; i < probeCoreSub.size(); ++i) {
    amap.emplace_back(probeCoreSub[i], refCoreIdx[i]);
  }
  const int refCid = refCore.getConformer().getId();
  double best = std::numeric_limits<double>::quiet_NaN();
  for (int cid : cids) {
    if (props) {
      try {
        std::unique_ptr<ForceFields::ForceField> ff(
            MMFF::constructForceField(probe, props.get(), 1.0e8, cid));
        if (ff) {
          ff->initialize();
          ff->minimize(1000);
        }
      } catch (...) {
      }
    }
    double r;
    try {
      // alignMol transforms `probe`'s conformer into the crystal frame (core
      // aligned), so exit vectors can be compared directly against the
      // crystal's.
      r = MolAlign::alignMol(probe, refCore, cid, refCid, &amap);
    } catch (...) {
      continue;
    }
    if (std::isnan(best) || r < best) best = r;
    // Coverage metric: the BEST exit match anywhere in the pool (does SOME
    // conformer reach the bioactive exit vector?), tracked independently of
    // core RMSD.
    if (exitErrOut && !exits.empty()) {
      const Conformer &pc = probe.getConformer(cid);
      double sum = 0;
      for (const auto &e : exits) {
        RDGeom::Point3D v = pc.getAtomPos(e.nbrSub) - pc.getAtomPos(e.coreSub);
        v.normalize();
        double d = std::max(-1.0, std::min(1.0, v.dotProduct(e.crystalVec)));
        sum += std::acos(d) * 180.0 / M_PI;
      }
      const double meanErr = sum / exits.size();
      if (std::isnan(*exitErrOut) || meanErr < *exitErrOut)
        *exitErrOut = meanErr;
    }
  }
  return best;
}

//! PoC: for each AROMATIC-adjacent fragment of each crystal ligand, compare the
//! crystal core-RMSD of the isolated carbon-capped embed (current pipeline)
//! against embedding the fragment WITH its real aromatic ring system (the
//! neighbour's whole aromatic system, kept from the crystal graph), then
//! stripping it for the RMSD.  Tests whether real aromatic context gives a
//! better SINGLE conformer (option 2).  Env: PLATINUM_SDF, PLATINUM_MAXMOLS
//! (default 40), POC_MMFF (MMFF-min the draw).
void runPocCtx() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = 40;
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);
  const bool mmff = std::getenv("POC_MMFF") != nullptr;
  unsigned int poolN = mmff ? 4 : 1;  // pool size: coverage sweep via POC_POOL
  if (const char *e = std::getenv("POC_POOL")) poolN = std::atoi(e);
  const int seed = 42;
  std::unique_ptr<SDMolSupplier> sup;
  try {
    sup = std::make_unique<SDMolSupplier>(sdf, true, false);
  } catch (...) {
    std::cerr << "pocctx: cannot open " << sdf << "\n";
    return;
  }
  RDLog::LogStateSetter blocker;
  int nFrag = 0, betterCtx = 0, worseCtx = 0;
  double sumIso = 0, sumCtx = 0;
  int nBucket[2] = {0, 0}, rBetterB[2] = {0, 0};
  double rIsoB[2] = {0, 0}, rCtxB[2] = {0, 0};  // core RMSD by [rigid, flex]
  size_t nMol = 0;
  while (!sup->atEnd() && nMol < maxMols) {
    std::unique_ptr<ROMol> raw;
    try {
      raw.reset(sup->next());
    } catch (...) {
      continue;
    }
    if (!raw || raw->getNumConformers() == 0) continue;
    boost::shared_ptr<RWMol> heavy;
    try {
      MolOps::assignStereochemistryFrom3D(*raw);
      std::unique_ptr<ROMol> h(MolOps::removeHs(*raw));
      ROMOL_SPTR frag = largestFrag(*h);
      if (!frag || frag->getNumAtoms() < 3) continue;
      heavy = boost::make_shared<RWMol>(*frag);
      MolOps::fastFindRings(*heavy);
    } catch (...) {
      continue;
    }
    ++nMol;
    const ROMol &H = *heavy;
    std::vector<unsigned int> links;
    try {
      links = FragmentConfGen::findLinkBonds(H);
    } catch (...) {
      continue;
    }
    if (links.empty()) continue;
    std::set<std::pair<int, int>> linkbonds;
    for (unsigned int bi : links) {
      const Bond *b = H.getBondWithIdx(bi);
      int x = b->getBeginAtomIdx(), y = b->getEndAtomIdx();
      linkbonds.insert({std::min(x, y), std::max(x, y)});
    }
    const int n = static_cast<int>(H.getNumAtoms());
    // fragment components: BFS not crossing link bonds
    std::vector<int> comp(n, -1);
    int nc = 0;
    for (int s = 0; s < n; ++s) {
      if (comp[s] >= 0) continue;
      std::vector<int> stack{s};
      comp[s] = nc;
      while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        for (const auto nb : H.atomNeighbors(H.getAtomWithIdx(u))) {
          int v = nb->getIdx();
          if (linkbonds.count({std::min(u, v), std::max(u, v)})) continue;
          if (comp[v] < 0) {
            comp[v] = nc;
            stack.push_back(v);
          }
        }
      }
      ++nc;
    }
    for (int ci = 0; ci < nc; ++ci) {
      std::vector<int> core;
      for (int a = 0; a < n; ++a) {
        if (comp[a] == ci) core.push_back(a);
      }
      if (core.size() < 2) continue;
      std::vector<std::pair<int, int>> exits;  // (coreAtom, neighbour)
      bool hasAro = false;
      std::set<int> aroCtx;
      for (const auto &lb : linkbonds) {
        bool ain = comp[lb.first] == ci, bin = comp[lb.second] == ci;
        if (ain == bin) continue;
        int c = ain ? lb.first : lb.second, nbr = ain ? lb.second : lb.first;
        exits.push_back({c, nbr});
        if (H.getAtomWithIdx(nbr)->getIsAromatic()) {
          hasAro = true;
          std::vector<int> st{nbr};
          while (!st.empty()) {
            int u = st.back();
            st.pop_back();
            if (aroCtx.count(u)) continue;
            aroCtx.insert(u);
            for (const auto nb : H.atomNeighbors(H.getAtomWithIdx(u))) {
              if (nb->getIsAromatic()) st.push_back(nb->getIdx());
            }
          }
        }
      }
      if (!hasAro) continue;  // only aromatic-adjacent fragments
      // reference: crystal core (core atoms only)
      std::map<int, int> r2s;
      auto refCore = submolWithCoords(H, core, r2s);
      try {
        MolOps::sanitizeMol(*refCore);
      } catch (...) {
        continue;
      }
      std::vector<int> refCoreIdx;
      for (int pi : core) refCoreIdx.push_back(r2s[pi]);
      const Conformer &Hc = H.getConformer();
      // ISOLATED probe: core + a carbon cap per exit (the aromatic exit's
      // direction is defined by the carbon cap embedded in isolation).
      std::map<int, int> p2s;
      auto iso = submolWithCoords(H, core, p2s);
      std::vector<int> isoCoreSub;
      for (int pi : core) isoCoreSub.push_back(p2s[pi]);
      std::vector<ExitProbe> isoExits;
      for (const auto &e : exits) {
        unsigned int cc = iso->addAtom(new Atom(6), false, true);
        iso->addBond(p2s[e.first], cc, Bond::SINGLE);
        if (H.getAtomWithIdx(e.second)->getIsAromatic()) {
          RDGeom::Point3D cv = Hc.getAtomPos(e.second) - Hc.getAtomPos(e.first);
          cv.normalize();
          isoExits.push_back({p2s[e.first], static_cast<int>(cc), cv});
        }
      }
      // CONTEXT probe: core + real aromatic system; carbon-cap only aliphatic
      // exits.
      std::vector<int> ctxAtoms = core;
      for (int a : aroCtx) {
        if (comp[a] != ci) ctxAtoms.push_back(a);
      }
      std::map<int, int> c2s;
      auto ctx = submolWithCoords(H, ctxAtoms, c2s);
      std::vector<int> ctxCoreSub;
      for (int pi : core) ctxCoreSub.push_back(c2s[pi]);
      std::vector<ExitProbe> ctxExits;
      for (const auto &e : exits) {
        if (!H.getAtomWithIdx(e.second)->getIsAromatic()) {
          unsigned int cc = ctx->addAtom(new Atom(6), false, true);
          ctx->addBond(c2s[e.first], cc, Bond::SINGLE);
        } else if (c2s.count(e.second)) {
          RDGeom::Point3D cv = Hc.getAtomPos(e.second) - Hc.getAtomPos(e.first);
          cv.normalize();
          ctxExits.push_back({c2s[e.first], c2s[e.second], cv});
        }
      }
      // Is the core flexible (a non-aromatic ring that can pucker/N-invert,
      // where the aromatic substituent context could reshape the conformation)?
      bool flexCore = false;
      for (int pi : core) {
        const Atom *a = H.getAtomWithIdx(pi);
        if (!a->getIsAromatic() && H.getRingInfo()->numAtomRings(pi) > 0) {
          flexCore = true;
          break;
        }
      }
      double eIso, eCtx;
      double rIso = embedAlignCore(*iso, *refCore, isoCoreSub, refCoreIdx,
                                   isoExits, mmff, seed, poolN, &eIso);
      double rCtx = embedAlignCore(*ctx, *refCore, ctxCoreSub, refCoreIdx,
                                   ctxExits, mmff, seed, poolN, &eCtx);
      if (std::isnan(rIso) || std::isnan(rCtx) || std::isnan(eIso) ||
          std::isnan(eCtx)) {
        continue;
      }
      ++nFrag;
      sumIso += eIso;
      sumCtx += eCtx;
      if (eCtx < eIso - 2.0) ++betterCtx;  // >2 deg better exit vector
      if (eCtx > eIso + 2.0) ++worseCtx;
      const int b = flexCore ? 1 : 0;
      nBucket[b]++;
      rIsoB[b] += rIso;
      rCtxB[b] += rCtx;
      if (rCtx < rIso - 0.05) rBetterB[b]++;
    }
  }
  std::printf(
      "\nPoC: aromatic-context vs isolated carbon-cap embed (%s)\n"
      "  %d aromatic exits over %zu mols\n",
      mmff ? "ETKDG+MMFF" : "ETKDG", nFrag, nMol);
  std::printf(
      "  EXIT-VECTOR angle err: isolated %.1f deg vs context %.1f deg  "
      "(ctx better %d / worse %d / tie %d)\n",
      nFrag ? sumIso / nFrag : 0.0, nFrag ? sumCtx / nFrag : 0.0, betterCtx,
      worseCtx, nFrag - betterCtx - worseCtx);
  const char *bl[2] = {"rigid core ", "FLEX core  "};
  for (int b = 0; b < 2; ++b) {
    if (!nBucket[b]) continue;
    std::printf(
        "  CORE RMSD [%s n=%3d]: isolated %.3f vs context %.3f A  "
        "(ctx better %d)\n",
        bl[b], nBucket[b], rIsoB[b] / nBucket[b], rCtxB[b] / nBucket[b],
        rBetterB[b]);
  }
  std::fflush(stdout);
}

// ---------------------------------------------------------------------------
// fragrank: SAMPLING vs RANKING diagnostic for flexible-ring fragments
// ---------------------------------------------------------------------------

//! Replicates Fraglib.cpp's classifyFragment (which lives in an anonymous
//! namespace and cannot be linked here) for a fragment CORE submol.  Returns
//! the flexible-ring class as a small code: 1=SmallRing (5-9-membered
//! saturated), 2=LargeRing (>=10), 0=other (Rigid: all-aromatic / <=4-membered
//! / fused; or Acyclic).  maxRingOut gets the largest ring size.
int classifyCoreRing(RWMol &frag, int &maxRingOut) {
  maxRingOut = 0;
  if (!frag.getRingInfo()->isInitialized()) MolOps::fastFindRings(frag);
  const RingInfo *ri = frag.getRingInfo();
  if (ri->numRings() == 0) return 0;  // Acyclic
  std::size_t maxRing = 0;
  bool allAromatic = true, fused = false;
  for (const auto &ring : ri->atomRings()) {
    maxRing = std::max(maxRing, ring.size());
    for (int idx : ring) {
      const auto u = static_cast<unsigned int>(idx);
      if (!frag.getAtomWithIdx(u)->getIsAromatic()) allAromatic = false;
      if (ri->numAtomRings(u) >= 2) fused = true;
    }
  }
  maxRingOut = static_cast<int>(maxRing);
  if (allAromatic || maxRing <= 4 || fused) return 0;  // Rigid
  return maxRing >= 10 ? 2 : 1;                        // LargeRing : SmallRing
}

//! Single-point MMFF94s (no electrostatics) energy of conformer `cid` of `mol`,
//! evaluated with an already-built topology scorer `scoreFn`.  No minimisation.
double singlePointE(
    const std::function<double(const double *, unsigned int)> &scoreFn,
    const ROMol &mol, int cid) {
  const Conformer &c = mol.getConformer(cid);
  const unsigned int n = mol.getNumAtoms();
  std::vector<double> buf(3 * static_cast<size_t>(n));
  for (unsigned int a = 0; a < n; ++a) {
    const RDGeom::Point3D &p = c.getAtomPos(a);
    buf[3 * a] = p.x;
    buf[3 * a + 1] = p.y;
    buf[3 * a + 2] = p.z;
  }
  return scoreFn(buf.data(), n);
}

//! Answers ONE question for flexible-ring fragments: is our library's failure a
//! SAMPLING issue (the bioactive ring pucker is never generated) or a RANKING
//! issue (it IS generated but ranks high-energy so it isn't selected)?
//!
//! For each crystal ligand, each flexible-ring fragment (SmallRing/LargeRing)
//! is cut out (BFS not crossing link bonds, exactly like pocctx), embedded as
//! an ISOLATED fragment with carbon-capped exits via ETKDGv3 (matching
//! embedFragment
//! -- raw pool, NO MMFF minimisation, FRAGRANK_POOL confs).  For every pool
//! conformer we take (a) its single-point MMFF94s no-estat energy (the score
//! our pipeline ranks by) and (b) its heavy-atom core RMSD to the crystal
//! fragment geometry (aligned on the core atoms).  We then locate the pool
//! conformer CLOSEST to the crystal pucker and report where its ENERGY ranks in
//! the pool.
//!   minRMS small  -> SAMPLING ok (pucker is in the pool)
//!   rank/pctile high with minRMS small -> RANKING problem (present but ranked
//!   out)
//! Env: PLATINUM_SDF, PLATINUM_MAXMOLS (default 60), FRAGRANK_POOL (default
  //! 50).
void runFragRank() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = 60;
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);
  unsigned int poolN = 50;
  if (const char *e = std::getenv("FRAGRANK_POOL")) poolN = std::atoi(e);
  const int seed = 42;
  const std::string var = "MMFF94s";

  std::unique_ptr<SDMolSupplier> sup;
  try {
    sup = std::make_unique<SDMolSupplier>(sdf, true, false);
  } catch (...) {
    std::cerr << "fragrank: cannot open " << sdf << "\n";
    return;
  }
  RDLog::LogStateSetter blocker;

  // aggregates over all flexible-ring fragments
  int nRingFrag = 0, nSampled = 0, nSmall = 0, nLarge = 0;
  double sumPctileSampled = 0.0, sumDEsampled = 0.0;
  int nTop20 = 0, nBottom50 = 0;

  size_t nMol = 0;
  while (!sup->atEnd() && nMol < maxMols) {
    std::unique_ptr<ROMol> raw;
    try {
      raw.reset(sup->next());
    } catch (...) {
      continue;
    }
    if (!raw || raw->getNumConformers() == 0) continue;
    boost::shared_ptr<RWMol> heavy;
    try {
      MolOps::assignStereochemistryFrom3D(*raw);
      std::unique_ptr<ROMol> h(MolOps::removeHs(*raw));
      ROMOL_SPTR frag = largestFrag(*h);
      if (!frag || frag->getNumAtoms() < 3) continue;
      heavy = boost::make_shared<RWMol>(*frag);
      MolOps::fastFindRings(*heavy);
    } catch (...) {
      continue;
    }
    const size_t molIdx = nMol;  // 0-based index of accepted molecules
    ++nMol;
    const ROMol &H = *heavy;

    std::vector<unsigned int> links;
    try {
      links = FragmentConfGen::findLinkBonds(H);
    } catch (...) {
      continue;
    }
    if (links.empty()) continue;
    std::set<std::pair<int, int>> linkbonds;
    for (unsigned int bi : links) {
      const Bond *b = H.getBondWithIdx(bi);
      int x = b->getBeginAtomIdx(), y = b->getEndAtomIdx();
      linkbonds.insert({std::min(x, y), std::max(x, y)});
    }
    const int n = static_cast<int>(H.getNumAtoms());
    // fragment components: BFS not crossing link bonds (same recipe as pocctx)
    std::vector<int> comp(n, -1);
    int nc = 0;
    for (int s = 0; s < n; ++s) {
      if (comp[s] >= 0) continue;
      std::vector<int> stack{s};
      comp[s] = nc;
      while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        for (const auto nb : H.atomNeighbors(H.getAtomWithIdx(u))) {
          int v = nb->getIdx();
          if (linkbonds.count({std::min(u, v), std::max(u, v)})) continue;
          if (comp[v] < 0) {
            comp[v] = nc;
            stack.push_back(v);
          }
        }
      }
      ++nc;
    }

    // --- best-to-crystal conformer (whole molecule), computed once ---
    // pose.
    std::map<int, int> cry2cdp;
    
    for (int ci = 0; ci < nc; ++ci) {
      std::vector<int> core;
      for (int a = 0; a < n; ++a) {
        if (comp[a] == ci) core.push_back(a);
      }
      if (core.size() < 3) continue;

      // classify the fragment core; keep only flexible rings
      std::map<int, int> c2class;
      auto classMol = submolWithCoords(H, core, c2class);
      try {
        MolOps::sanitizeMol(*classMol);
      } catch (...) {
        continue;
      }
      int ringSize = 0;
      const int cls = classifyCoreRing(*classMol, ringSize);
      if (cls == 0) continue;  // not a flexible ring
      const char *clsName = cls == 2 ? "LargeRing" : "SmallRing";

      // exits crossing into this component (coreAtom, neighbour)
      std::vector<std::pair<int, int>> exits;
      for (const auto &lb : linkbonds) {
        bool ain = comp[lb.first] == ci, bin = comp[lb.second] == ci;
        if (ain == bin) continue;
        int c = ain ? lb.first : lb.second, nbr = ain ? lb.second : lb.first;
        exits.push_back({c, nbr});
      }

      // crystal reference core (core atoms only, crystal geometry)
      std::map<int, int> r2s;
      auto refCore = submolWithCoords(H, core, r2s);
      try {
        MolOps::sanitizeMol(*refCore);
      } catch (...) {
        continue;
      }
      std::vector<int> refCoreIdx;
      for (int pi : core) refCoreIdx.push_back(r2s[pi]);

      // isolated probe = core + a CARBON cap per exit (matches
      // embedFragment's dummy->carbon capping); this is exactly the
      // fragment our pipeline embeds.
      std::map<int, int> p2s;
      auto probe = submolWithCoords(H, core, p2s);
      std::vector<int> probeCoreSub;
      for (int pi : core) probeCoreSub.push_back(p2s[pi]);
      for (const auto &e : exits) {
        unsigned int cc = probe->addAtom(new Atom(6), false, true);
        probe->addBond(p2s[e.first], cc, Bond::SINGLE);
      }
      try {
        MolOps::sanitizeMol(*probe);
      } catch (...) {
        continue;
      }
      MolOps::addHs(*probe);
      DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
      ps.randomSeed = seed;
      probe->clearConformers();
      INT_VECT cids;
      try {
        DGeomHelpers::EmbedMultipleConfs(*probe, cids, poolN, ps);
      } catch (...) {
        continue;
      }
      if (cids.empty()) continue;

      auto scoreFn = makeFullFFScoreFn(*probe, /*electrostatics=*/false, var);
      if (!scoreFn) continue;  // cannot MMFF-type -> cannot rank, skip

      MatchVectType amap;  // (probeIdx, refIdx) over core atoms
      for (size_t i = 0; i < probeCoreSub.size(); ++i) {
        amap.emplace_back(probeCoreSub[i], refCoreIdx[i]);
      }
      const int refCid = refCore->getConformer().getId();

      std::vector<double> energies, rmsds;
      energies.reserve(cids.size());
      rmsds.reserve(cids.size());
      for (int cid : cids) {
        double e = singlePointE(scoreFn, *probe, cid);
        double r;
        try {
          // aligns pool conformer cid into the crystal frame over the core
          // atoms
          r = MolAlign::alignMol(*probe, *refCore, cid, refCid, &amap);
        } catch (...) {
          continue;
        }
        if (!std::isfinite(e) || !std::isfinite(r)) continue;
        energies.push_back(e);
        rmsds.push_back(r);
      }
      const int k = static_cast<int>(energies.size());
      if (k == 0) continue;

      // pool conformer closest to the crystal pucker
      int minIdx = 0;
      for (int i = 1; i < k; ++i)
        if (rmsds[i] < rmsds[minIdx]) minIdx = i;
      const double minRMS = rmsds[minIdx];
      const double eClosest = energies[minIdx];
      const double eMin = *std::min_element(energies.begin(), energies.end());
      int rank =
          0;  // 0-based energy rank: how many pool confs are strictly lower
      for (int i = 0; i < k; ++i)
        if (energies[i] < eClosest) ++rank;
      const double pctile = 100.0 * rank / k;
      const double dE = eClosest - eMin;

      ++nRingFrag;
      if (cls == 2)
        ++nLarge;
      else
        ++nSmall;
      const bool sampled = minRMS < 0.5;
      if (sampled) {
        ++nSampled;
        sumPctileSampled += pctile;
        sumDEsampled += dE;
        if (pctile < 20.0) ++nTop20;
        if (pctile >= 50.0) ++nBottom50;
      }

      std::printf(
          "[fragrank] mol=%zu class=%s ringSize=%d poolN=%d minRMS=%.3f "
          "rankOfClosest=%d/%d energyPctile=%.0f dEfromMin=%.2f",
          molIdx, clsName, ringSize, k, minRMS, rank, k, pctile, dE);

      // --- fragment fields ---
            std::printf("\n");
      std::fflush(stdout);
    }
  }

  std::printf(
      "\n[fragrank] AGGREGATE over %d flexible-ring fragments "
      "(%d SmallRing, %d LargeRing) from %zu mols, poolN=%u\n",
      nRingFrag, nSmall, nLarge, nMol, poolN);
  const double sampFrac = nRingFrag ? 100.0 * nSampled / nRingFrag : 0.0;
  std::printf(
      "  SAMPLING: %d/%d = %.1f%% have minRMS<0.5A "
      "(crystal pucker present in the pool)\n",
      nSampled, nRingFrag, sampFrac);
  if (nSampled) {
    std::printf(
        "  RANKING (among the %d sampled): mean energyPctile of the "
        "crystal-closest conf = %.1f%%\n"
        "    in top-20%% energy: %d/%d = %.1f%%   in bottom-50%%: %d/%d = %.1f%%\n"
        "    mean dEfromMin = %.2f kcal/mol\n",
        nSampled, sumPctileSampled / nSampled, nTop20, nSampled,
        100.0 * nTop20 / nSampled, nBottom50, nSampled,
        100.0 * nBottom50 / nSampled, sumDEsampled / nSampled);
  }
    std::fflush(stdout);
}

// ---------------------------------------------------------------------------
// termtable: per-conformer FULL MMFF94s (no-estat) energy-term breakdown for
// our pipeline's ensemble, so we can see WHICH energy term
// makes the crystal-closest ("bioactive") conformer rank the way it does.
// ---------------------------------------------------------------------------

//! Exact INTER-fragment stretch-bend energy callable, built by mirroring the
//! MMFF Builder::addStretchBend enumeration (central atom j, neighbour pair
//! i,k) but keeping ONLY angles whose atom set {i,j,k} spans >1 fragment (i.e.
//! crosses a link bond).  Same buffered coords contract as
//! makeInterFragmentScoreFn: the FF is typed/built once and evaluated on
//! whatever coord buffer (matching `mol`'s atom order) it is handed.  Null
//! callable if MMFF typing fails.
RotorDriver::ScoreFn makeInterStretchBendScoreFn(const ROMol &mol,
                                                 const std::vector<int> &fragOf,
                                                 const std::string &variant) {
  auto owned = boost::make_shared<ROMol>(mol);
  if (owned->getNumConformers() == 0) {
    owned->addConformer(new Conformer(owned->getNumAtoms()), true);
  }
  MMFF::MMFFMolProperties props(*owned, variant);
  if (!props.isValid()) {
    return {};
  }
  props.setMMFFEleTerm(false);
  const unsigned int nAtoms = owned->getNumAtoms();
  Conformer &conf = owned->getConformer();
  auto ff = boost::make_shared<ForceFields::ForceField>();
  for (unsigned int i = 0; i < nAtoms; ++i) {
    ff->positions().push_back(&conf.getAtomPos(i));
  }
  auto contrib =
      std::make_unique<ForceFields::MMFF::StretchBendContrib>(ff.get());
  bool any = false;
  for (unsigned int j = 0; j < nAtoms; ++j) {
    const Atom *jA = owned->getAtomWithIdx(j);
    if (jA->getDegree() < 2) {
      continue;
    }
    std::vector<unsigned int> nbrs;
    for (const auto nb : owned->atomNeighbors(jA)) {
      nbrs.push_back(nb->getIdx());
    }
    for (size_t a = 0; a < nbrs.size(); ++a) {
      for (size_t b = a + 1; b < nbrs.size(); ++b) {
        const unsigned int i = nbrs[a], k = nbrs[b];
        // inter iff the i-j-k atom set is NOT all one fragment
        if (fragOf[i] == fragOf[j] && fragOf[k] == fragOf[j]) {
          continue;
        }
        unsigned int sbType;
        ForceFields::MMFF::MMFFStbn stbn;
        ForceFields::MMFF::MMFFBond bond[2];
        ForceFields::MMFF::MMFFAngle ang;
        if (props.getMMFFStretchBendParams(*owned, i, j, k, sbType, stbn, bond,
                                           ang)) {
          contrib->addTerm(i, j, k, &stbn, &ang, &bond[0], &bond[1]);
          any = true;
        }
      }
    }
  }
  if (!any) {
    return [](const double *, unsigned int) { return 0.0; };  // no inter-SB
  }
  ff->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
  ff->initialize();
  return [owned, ff](const double *pos, unsigned int) -> double {
    try {
      double *p = const_cast<double *>(pos);
      double e = 0.0;
      for (const auto &c : ff->contribs()) {
        e += c->getEnergy(p);
      }
      return e;
    } catch (...) {
      return 0.0;
    }
  };
}

void runTermTable() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = 40;  // # SDF records to SCAN
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);
  size_t maxProc = 5;  // # QUALIFYING (>1 fragment) mols to actually process
  if (const char *e = std::getenv("TERMTABLE_MOLS")) maxProc = std::atoi(e);
  unsigned int ourConfs = 40;  // cap OUR ensemble (readable tables)
  if (const char *e = std::getenv("TERMTABLE_OURCONFS"))
    ourConfs = std::atoi(e);
  const char *wantName =
      std::getenv("PLATINUM_NAME");  // restrict to one ligand
  const int seed = 42;
  const std::string var = "MMFF94s";

  std::unique_ptr<SDMolSupplier> sup;
  try {
    sup = std::make_unique<SDMolSupplier>(sdf, true, false);
  } catch (...) {
    std::cerr << "termtable: cannot open " << sdf << "\n";
    return;
  }
  RDLog::LogStateSetter blocker;

  std::printf(
        "\n[termtable] per-conformer MMFF94s (no-estat) term breakdown\n"
        "  ourConfs<=%u\n"
      "  total = bond+angle+stretchBend+oop+torsion+vdw (sum of the six decomposed "
      "terms)\n"
      "  interVdw = cross-fragment vdW (makeInterFragmentScoreFn, no junctions); "
      "interSB = stretch-bend terms whose i-j-k span a link bond\n",
      ourConfs);

  // A single output-conformer row.
  struct Row {
    double crystalRms = std::numeric_limits<double>::quiet_NaN();
    double total = 0, bond = 0, angle = 0, sb = 0, oop = 0, tors = 0, vdw = 0;
    double interVdw = 0, interSB = 0;
    bool closest = false;
  };

  auto coordBuf = [](const ROMol &m) {
    const Conformer &c = m.getConformer();
    const unsigned int n = m.getNumAtoms();
    std::vector<double> buf(3 * static_cast<size_t>(n));
    for (unsigned int a = 0; a < n; ++a) {
      const auto &p = c.getAtomPos(a);
      buf[3 * a] = p.x;
      buf[3 * a + 1] = p.y;
      buf[3 * a + 2] = p.z;
    }
    return buf;
  };

  // Build fragOf (BFS not crossing link bonds) and the junction bond list for a
  // (Hs-bearing) representative conformer -- topology is constant within an
  // ensemble, so this is computed once per ensemble.  Returns false if the mol
  // has no link bonds (single fragment).
  auto buildFragOf =
      [](const ROMol &m, std::vector<int> &fragOf,
         std::vector<std::pair<unsigned int, unsigned int>> &junctions)
      -> bool {
    std::vector<unsigned int> links;
    try {
      links = FragmentConfGen::findLinkBonds(m);
    } catch (...) {
      return false;
    }
    if (links.empty()) return false;
    std::set<std::pair<int, int>> linkbonds;
    junctions.clear();
    for (unsigned int bi : links) {
      const Bond *b = m.getBondWithIdx(bi);
      int x = b->getBeginAtomIdx(), y = b->getEndAtomIdx();
      linkbonds.insert({std::min(x, y), std::max(x, y)});
      junctions.emplace_back(b->getBeginAtomIdx(), b->getEndAtomIdx());
    }
    const int n = static_cast<int>(m.getNumAtoms());
    fragOf.assign(n, -1);
    int nc = 0;
    for (int s = 0; s < n; ++s) {
      if (fragOf[s] >= 0) continue;
      std::vector<int> stack{s};
      fragOf[s] = nc;
      while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        for (const auto nb : m.atomNeighbors(m.getAtomWithIdx(u))) {
          int v = nb->getIdx();
          if (linkbonds.count({std::min(u, v), std::max(u, v)})) continue;
          if (fragOf[v] < 0) {
            fragOf[v] = nc;
            stack.push_back(v);
          }
        }
      }
      ++nc;
    }
    return true;
  };

  // Score one ensemble: build rows, sort ascending by total MMFF, print the
  // table
  // + the crystal-closest summary line.  `refHeavy` is the crystal heavy-atom
  // pose.
  auto emitEnsemble = [&](const char *label,
                          const std::vector<ROMOL_SPTR> &confs,
                          const ROMol &refHeavy) {
    RWMol ref(refHeavy);
    std::vector<int> fragOf;
    std::vector<std::pair<unsigned int, unsigned int>> junctions;
    RotorDriver::ScoreFn interVdwFn, interSbFn;
    if (!confs.empty()) {
      if (buildFragOf(*confs.front(), fragOf, junctions)) {
        // interVdw: cross-fragment vdW ONLY (pass NO junction bonds -> no
        // torsion terms, so the callable returns pure cross-fragment vdW).
        interVdwFn = makeInterFragmentScoreFn(
            *confs.front(), fragOf, /*junctionBonds=*/{},
            /*electrostatics=*/false, var, /*vdwCutoff=*/0.0, nullptr);
        interSbFn = makeInterStretchBendScoreFn(*confs.front(), fragOf, var);
      }
    }

    std::vector<Row> rows;
    for (const auto &m : confs) {
      if (!m || m->getNumConformers() == 0) continue;
      MmffTerms t = mmffDecompose(*m, var);
      const double total = t.bond + t.angle + t.sb + t.oop + t.tors + t.vdw;
      if (!std::isfinite(total)) continue;  // untypeable conformer -> skip
      Row r;
      r.bond = t.bond;
      r.angle = t.angle;
      r.sb = t.sb;
      r.oop = t.oop;
      r.tors = t.tors;
      r.vdw = t.vdw;
      r.total = total;
      const auto buf = coordBuf(*m);
      if (interVdwFn) r.interVdw = interVdwFn(buf.data(), m->getNumAtoms());
      if (interSbFn) r.interSB = interSbFn(buf.data(), m->getNumAtoms());
      // crystal heavy-atom RMSD (symmetry-aware, per single conformer)
      try {
        std::unique_ptr<ROMol> probe(MolOps::removeHs(*m));
        if (probe->getNumAtoms() == ref.getNumAtoms()) {
          r.crystalRms = MolAlign::getBestRMS(*probe, ref);
        }
      } catch (...) {
      }
      rows.push_back(r);
    }

    std::printf("\n  --- %s (%zu conformers) ---\n", label, rows.size());
    if (rows.empty()) {
      std::printf("    (no MMFF-typeable conformers)\n");
      return;
    }
    std::sort(rows.begin(), rows.end(),
              [](const Row &a, const Row &b) { return a.total < b.total; });
    // mark the crystal-closest row (min finite crystalRMS)
    int closestRank = -1;
    double bestRms = std::numeric_limits<double>::max();
    for (size_t i = 0; i < rows.size(); ++i) {
      if (std::isfinite(rows[i].crystalRms) && rows[i].crystalRms < bestRms) {
        bestRms = rows[i].crystalRms;
        closestRank = static_cast<int>(i);
      }
    }
    if (closestRank >= 0) rows[closestRank].closest = true;
    // best-FINDABLE hit: the LOWEST-ENERGY conf that is still <1A -- what
    // energy selection actually surfaces -- vs the crystal-CLOSEST, which may
    // be a strained outlier that merely nudged nearest in RMSD. (rows are
    // energy-sorted ascending.)
    int bestHitRank = -1;
    for (size_t i = 0; i < rows.size(); ++i)
      if (std::isfinite(rows[i].crystalRms) && rows[i].crystalRms < 1.0) {
        bestHitRank = static_cast<int>(i);
        break;
      }

    std::printf(
        "    rank  crysRMS     total     bond    angle  strBend      oop  "
        "torsion      vdw  interVdw   interSB\n");
    for (size_t i = 0; i < rows.size(); ++i) {
      const Row &r = rows[i];
      char rms[16];
      if (std::isfinite(r.crystalRms))
        std::snprintf(rms, sizeof(rms), "%7.2f", r.crystalRms);
      else
        std::snprintf(rms, sizeof(rms), "%7s", "n/a");
      std::printf(
          "    %3zu%c %s %9.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %9.1f %9.1f\n",
          i + 1, r.closest ? '*' : ' ', rms, r.total, r.bond, r.angle, r.sb,
          r.oop, r.tors, r.vdw, r.interVdw, r.interSB);
    }

    if (closestRank < 0) {
      std::printf(
          "    summary: no conformer is graph-comparable to the crystal pose "
          "(RMSD n/a)\n");
      return;
    }
    const Row &cl = rows[closestRank];
    const Row &r1 = rows.front();
    // Which of the six primary terms pushes the bioactive conf's total above
    // rank-1?
    struct TD {
      const char *name;
      double delta;
    };
    const TD cand[6] = {
        {"bond", cl.bond - r1.bond},    {"angle", cl.angle - r1.angle},
        {"stretchBend", cl.sb - r1.sb}, {"oop", cl.oop - r1.oop},
        {"torsion", cl.tors - r1.tors}, {"vdw", cl.vdw - r1.vdw}};
    const TD *worst = &cand[0];
    for (int i = 1; i < 6; ++i)
      if (cand[i].delta > worst->delta) worst = &cand[i];
    if (closestRank == 0) {
      std::printf(
          "    summary: crystal-closest conf is rank 1/%zu (RMSD=%.2f), total=%.1f "
          "-- the bioactive pose IS the lowest-energy conformer\n",
          rows.size(), cl.crystalRms, cl.total);
    } else {
      std::printf(
          "    summary: crystal-closest conf is rank %d/%zu (RMSD=%.2f), "
          "total=%.1f (rank-1 total=%.1f, +%.1f above min); its biggest term above "
          "the rank-1 conf is %s (+%.1f kcal)\n",
          closestRank + 1, rows.size(), cl.crystalRms, cl.total, r1.total,
          cl.total - r1.total, worst->name, worst->delta);
    }
    if (bestHitRank >= 0) {
      const Row &bh = rows[bestHitRank];
      std::printf(
          "    besthit: lowest-energy <1A conf is rank %d/%zu (RMSD=%.2f), +%.1f "
          "above min\n",
          bestHitRank + 1, rows.size(), bh.crystalRms, bh.total - r1.total);
    } else {
      std::printf("    besthit: no conf <1A in ensemble\n");
    }
  };

  size_t nScanned = 0, nProc = 0;
  while (!sup->atEnd() && nScanned < maxMols && nProc < maxProc) {
    std::unique_ptr<ROMol> raw;
    try {
      raw.reset(sup->next());
    } catch (...) {
      continue;
    }
    ++nScanned;
    if (!raw || raw->getNumConformers() == 0) continue;
    ROMOL_SPTR heavy;
    std::string name;
    try {
      if (raw->hasProp("_Name")) name = raw->getProp<std::string>("_Name");
      MolOps::assignStereochemistryFrom3D(*raw);
      std::unique_ptr<ROMol> h(MolOps::removeHs(*raw));
      heavy = largestFrag(*h);
      if (!heavy || heavy->getNumAtoms() < 3) continue;
    } catch (...) {
      continue;
    }
    if (wantName && name != std::string(wantName)) continue;
    // qualify: must actually have >1 fragment (a link bond)
    std::vector<unsigned int> links;
    try {
      links = FragmentConfGen::findLinkBonds(*heavy);
    } catch (...) {
      continue;
    }
    if (links.empty()) continue;  // single rigid fragment -> nothing to split

    std::string smi;
    try {
      smi = MolToSmiles(*heavy);
    } catch (...) {
      continue;
    }
    unsigned int nrot = 0;
    try {
      nrot = Descriptors::calcNumRotatableBonds(*heavy);
    } catch (...) {
    }
    ++nProc;
    std::printf(
        "\n============================================================\n"
        "[mol %zu] %s   nAtoms(heavy)=%u  nrot=%u  linkBonds=%zu\n  %s\n",
        nProc, name.empty() ? "?" : name.c_str(), heavy->getNumAtoms(), nrot,
        links.size(), smi.c_str());

    // OUR ensemble: raw-ETKDG TS+TorLib default (the shipped pipeline default).
    auto in = boost::make_shared<RWMol>(*heavy);
    in->clearConformers();
    FragmentConfGenParams op = bestDefaultParams(seed);
    op.numOutputConfs = ourConfs;
    op.search.searchMode =
        searchModeFromEnv();         // ASM_SEARCH override (e.g. thompson)
    if (std::getenv("TT_TORLIB")) {  // TT_TORLIB=1 -> Hamburg TorLib sampler
      try {
        op.search.torsionSampler = hamburgSampler();
      } catch (...) {
      }
    }
    std::vector<ROMOL_SPTR> oursConfs;
    try {
      requireValidParams(op, "runTermTable");
      oursConfs = FragmentConfGen(op).build(*in).conformers;
    } catch (...) {
    }
    emitEnsemble("OURS", oursConfs, *heavy);

    
    // Link/rotor bonds the joiner drives (findLinkBonds), reported on the HEAVY
    // molecule -- NO explicit Hs are added here.  Every link bond joins two
    // heavy atoms, and addHs appends Hs AFTER all heavy atoms, so these
    // heavy-atom indices (and hence these bonds) are stable across a subsequent
    // addHs.
    {
      std::vector<int> fragOf;
      std::vector<std::pair<unsigned int, unsigned int>> junctions;
      buildFragOf(*heavy, fragOf, junctions);
      std::printf(
          "\n  [linkbonds] mol %zu %s  (HEAVY mol, no explicit Hs; all link atoms "
          "are heavy -> indices stable across addHs)\n",
          nProc, name.empty() ? "?" : name.c_str());
      std::printf("    used SMILES        : %s\n", smi.c_str());
      std::printf(
          "    canonical isomeric : %s   (identical -- the used SMILES IS "
          "Chem::MolToSmiles(heavy))\n",
          smi.c_str());
      std::printf("    %zu link/rotor bond(s) [bondIdx on heavy mol]:\n",
                  links.size());
      std::printf(
          "      bondIdx  atomI(idx:elem)  atomJ(idx:elem)  fragI->fragJ\n");
      for (unsigned int bi : links) {
        const Bond *b = heavy->getBondWithIdx(bi);
        const unsigned int ai = b->getBeginAtomIdx(), aj = b->getEndAtomIdx();
        std::printf("      %7u   %5u:%-2s         %5u:%-2s        %d->%d\n", bi,
                    ai, heavy->getAtomWithIdx(ai)->getSymbol().c_str(), aj,
                    heavy->getAtomWithIdx(aj)->getSymbol().c_str(), fragOf[ai],
                    fragOf[aj]);
      }
    }
    std::fflush(stdout);
  }
  std::printf(
      "\n[termtable] processed %zu qualifying molecule(s) (%zu scanned)\n",
      nProc, nScanned);
}

//! FRAGMENT GEOMETRY: how close is our EMBEDDED fragment pool to the crystal's
//! own geometry for that fragment, ignoring placement/orientation?
//!
//! This is the question three search-based probes failed to answer (see
//! each attempt to inject the crystal
//! fragment into the search perturbed the ensemble it was meant to hold fixed.
//! Measuring the geometry DIRECTLY sidesteps that entirely -- no search, no
//! assembly, no scoring, nothing to confound.  For every fragment we take the
//! best superposed (QCP) heavy-atom RMSD between the crystal's geometry for
//! that fragment and ANY conformer in the pool our embedder produced.
//!
//!   small  -> our fragments are right; the tail gap is downstream
//!   (assembly/selection) large  -> the embedder never produces the needed
//!   fragment conformer, and fragment-library
//!             work is the lever
void runFragRms() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = 150;
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);
  int minRot = 0, maxRot = std::numeric_limits<int>::max();
  if (const char *e = std::getenv("PLATINUM_MINROT")) minRot = std::atoi(e);
  if (const char *e = std::getenv("PLATINUM_MAXROT")) maxRot = std::atoi(e);
  const std::string var = "MMFF94s";
  const unsigned int nConfs = std::getenv("FRAGRMS_CONFS")
                                  ? std::atoi(std::getenv("FRAGRMS_CONFS"))
                                  : 16;

  std::unique_ptr<SDMolSupplier> sup;
  try {
    sup = std::make_unique<SDMolSupplier>(sdf, true, false);
  } catch (...) {
    std::cerr << "fragrms: cannot open " << sdf << "\n";
    return;
  }
  RDLog::LogStateSetter blk;

  std::vector<double> best;         // per fragment: best RMSD over the pool
  std::vector<double> worstPerMol;  // per molecule: its worst fragment
  std::vector<size_t> poolSize;
  size_t done = 0, nBig = 0;
  while (!sup->atEnd() && done < maxMols) {
    std::unique_ptr<ROMol> raw;
    try {
      raw.reset(sup->next());
    } catch (...) {
      continue;
    }
    if (!raw || raw->getNumConformers() == 0) continue;
    try {
      MolOps::assignStereochemistryFrom3D(*raw);
      ROMOL_SPTR big = largestFrag(*raw);
      if (!big || big->getNumAtoms() < 4) continue;
      int nr = 0;
      try {
        nr = static_cast<int>(Descriptors::calcNumRotatableBonds(*big));
      } catch (...) {
      }
      if (nr < minRot || nr > maxRot) continue;
      auto m = boost::make_shared<RWMol>(*big);
      MolOps::addHs(*m, false, true);

      // normal embedded pool
      unsetenv("ASM_EXACT_GEOM");
      FragmentJoinerInput emb = buildFragmentJoinerInput(
          *m, nConfs, 0xf00d, var, false, true, nullptr);
      if (!emb.mol || emb.fragments.empty()) continue;
      const Conformer &xc = m->getConformer();

      double worst = 0.0;
      bool any = false;
      for (const auto &f : emb.fragments) {
        if (f.confs.empty() || f.atoms.empty()) continue;
        // heavy atoms of this fragment, index-aligned between crystal and pool
        // conformers
        std::vector<unsigned int> idx;
        for (unsigned int a : f.atoms)
          if (m->getAtomWithIdx(a)->getAtomicNum() > 1) idx.push_back(a);
        if (idx.size() < 3) continue;  // RMSD is meaningless for 1-2 atoms
        std::vector<RDGeom::Point3D> xtal(m->getNumAtoms());
        for (unsigned int a : idx) xtal[a] = xc.getAtomPos(a);
        double bestR = std::numeric_limits<double>::max();
        for (const auto &c : f.confs) {
          try {
            bestR = std::min(bestR, qcpRmsd(xtal, c.pos, idx));
          } catch (...) {
          }
        }
        if (!std::isfinite(bestR) ||
            bestR == std::numeric_limits<double>::max())
          continue;
        best.push_back(bestR);
        poolSize.push_back(f.confs.size());
        if (bestR > 0.5) ++nBig;
        worst = std::max(worst, bestR);
        any = true;
      }
      if (any) worstPerMol.push_back(worst);
      ++done;
    } catch (...) {
      continue;
    }
  }

  auto q = [](std::vector<double> v, double p) {
    if (v.empty()) return std::nan("");
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, static_cast<size_t>(p * (v.size() - 1)))];
  };
  std::printf(
      "\n[fragrms] %zu molecules (rot %d..%s), %zu fragments, pool<=%u confs\n",
      done, minRot,
      maxRot == std::numeric_limits<int>::max()
          ? "inf"
          : std::to_string(maxRot).c_str(),
      best.size(), nConfs);
  std::printf(
      "  BEST superposed heavy-atom RMSD from the crystal fragment to our pool:\n");
  std::printf("    median=%.3f  p75=%.3f  p90=%.3f  p99=%.3f  max=%.3f\n",
              q(best, 0.5), q(best, 0.75), q(best, 0.90), q(best, 0.99),
              q(best, 1.0));
  std::printf(
      "    fragments worse than 0.25A: %.1f%%   worse than 0.5A: %.1f%%\n",
      100.0 *
          std::count_if(best.begin(), best.end(),
                        [](double x) { return x > 0.25; }) /
          std::max<size_t>(1, best.size()),
      100.0 * nBig / std::max<size_t>(1, best.size()));
  std::printf(
      "  per-MOLECULE worst fragment: median=%.3f  p90=%.3f  max=%.3f\n",
      q(worstPerMol, 0.5), q(worstPerMol, 0.90), q(worstPerMol, 1.0));
  std::printf(
      "  (a molecule can only reach sub-1A overall if EVERY fragment is close)\n");
}

//! TORSION COVERAGE: for every INTER-FRAGMENT junction, how far is the
//! CRYSTAL's actual dihedral from the nearest angle our samplers would ever
//! propose?  This is the direct test of "are we choosing the right junction
//! angles" -- it needs no search, no assembly and no scoring, so it cannot be
//! confounded by budget, window, ranking or output caps.
//!   deviation small  -> the right angles ARE on the table; any failure is the
//!   SEARCH deviation large  -> we can never build the crystal pose, whatever
//!   the search does
//! Compares ETKDG vs Hamburg TorLib, base angles vs full basin/tolerance sets,
//! plus a 30deg uniform grid as a reference ceiling.
//! PLATINUM_MINROT/MAXROT/MAXMOLS apply.
void runTorCheck() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = 200;
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);
  int minRot = 0, maxRot = std::numeric_limits<int>::max();
  if (const char *e = std::getenv("PLATINUM_MINROT")) minRot = std::atoi(e);
  if (const char *e = std::getenv("PLATINUM_MAXROT")) maxRot = std::atoi(e);
  const std::string var = "MMFF94s";

  // gapFill > 0: after taking the sampler's angles, fill any circular GAP wider
  // than `gapFill` degrees with evenly spaced midpoints.  This buys coverage
  // only where the library is sparse/badly spread, instead of unioning a full
  // uniform grid everywhere (which doubles arms per rotor and starves a fixed
  // stochastic budget).
  struct Samp {
    const char *name;
    std::shared_ptr<TorsionSampler> s;
    bool basin;
    double gapFill;
  };
  std::vector<Samp> samps = {
      {"ETKDG  base ", std::make_shared<ETKDGTorsionSampler>(), false, 0},
      {"ETKDG  basin", std::make_shared<ETKDGTorsionSampler>(), true, 0},
      {"TorLib base ", hamburgSampler(), false, 0},
      {"TorLib basin", hamburgSampler(), true, 0},
      {"TorLib g120 ", hamburgSampler(), true, 120},
      {"TorLib g90  ", hamburgSampler(), true, 90},
      {"TorLib g60  ", hamburgSampler(), true, 60},
      {"ETKDG  g60  ", std::make_shared<ETKDGTorsionSampler>(), true, 60},
  };
  auto fillGaps = [](std::vector<double> a, double maxGap) {
    if (maxGap <= 0 || a.size() < 2) return a;
    std::sort(a.begin(), a.end());
    std::vector<double> out = a;
    for (size_t i = 0; i < a.size(); ++i) {
      const double lo = a[i];
      const double hi = (i + 1 < a.size()) ? a[i + 1] : a[0] + 360.0;
      const double gap = hi - lo;
      if (gap <= maxGap) continue;
      const int nAdd = static_cast<int>(std::ceil(gap / maxGap)) - 1;
      for (int k = 1; k <= nAdd; ++k) {
        double v = lo + gap * k / (nAdd + 1);
        while (v > 180.0) v -= 360.0;
        out.push_back(v);
      }
    }
    return out;
  };
  const std::vector<double> grid30 = {-180, -150, -120, -90, -60, -30,
                                      0,    30,   60,   90,  120, 150};

  // per sampler: every rotor's deviation, and every molecule's WORST rotor
  // deviation
  std::vector<std::vector<double>> dev(samps.size() + 1),
      worst(samps.size() + 1);
  std::vector<size_t> nAngles(samps.size() + 1, 0), nRotor(samps.size() + 1, 0);

  auto circDev = [](double a, double b) {
    double d = std::fmod(a - b + 540.0, 360.0) - 180.0;
    return std::fabs(d);
  };

  std::unique_ptr<SDMolSupplier> sup;
  try {
    sup = std::make_unique<SDMolSupplier>(sdf, true, false);
  } catch (...) {
    std::cerr << "torcheck: cannot open " << sdf << "\n";
    return;
  }
  RDLog::LogStateSetter blk;

  size_t done = 0;
  while (!sup->atEnd() && done < maxMols) {
    std::unique_ptr<ROMol> raw;
    try {
      raw.reset(sup->next());
    } catch (...) {
      continue;
    }
    if (!raw || raw->getNumConformers() == 0) continue;
    try {
      MolOps::assignStereochemistryFrom3D(*raw);
      ROMOL_SPTR big = largestFrag(*raw);
      if (!big || big->getNumAtoms() < 4) continue;
      int nr = 0;
      try {
        nr = static_cast<int>(Descriptors::calcNumRotatableBonds(*big));
      } catch (...) {
      }
      if (nr < minRot || nr > maxRot) continue;

      auto m = boost::make_shared<RWMol>(*big);
      MolOps::addHs(*m, /*explicitOnly=*/false, /*addCoords=*/true);
      // exact crystal geometry, 1 conf/frag -- we only need the TOPOLOGY
      // (junctions + torsion quartets); the crystal's own coordinates supply
      // every dihedral.
      setenv("ASM_EXACT_GEOM", "1", 1);
      FragmentJoinerInput in =
          buildFragmentJoinerInput(*m, 1, 0xf00d, var, false, false, nullptr);
      if (!in.mol || in.fragments.size() < 2 || in.junctions.empty()) continue;

      FragmentJoinerParams gp;
      gp.ffVariant = var;
      // This diagnostic reads only the rotor topology and the crystal coords,
      // so it does not need a searchable context -- an empty rotor set is the
      // only thing that stops it.
      const auto ctx = joinFragments(in, gp);
      if (ctx.rotorBonds.empty()) continue;
      RotorDriver drv(*in.mol, ctx.rotorBonds, -1, nullptr);
      const Conformer &xc =
          in.mol->getConformer();  // == crystal coords (exact-geom path)

      std::vector<double> wr(samps.size() + 1, 0.0);
      size_t nr2 = ctx.rotorBonds.size();
      for (size_t r = 0; r < nr2; ++r) {
        const auto q = drv.torsion(static_cast<unsigned int>(r));
        double xtalAng = 0.0;
        try {
          xtalAng = MolTransforms::getDihedralDeg(xc, q[0], q[1], q[2], q[3]);
        } catch (...) {
          continue;
        }
        for (size_t si = 0; si < samps.size(); ++si) {
          std::vector<double> cand;
          try {
            cand = samps[si].s->getAngles(*in.mol, q[0], q[1], q[2], q[3],
                                          samps[si].basin);
          } catch (...) {
          }
          if (cand.empty()) continue;
          cand = fillGaps(std::move(cand), samps[si].gapFill);
          double best = 1e9;
          for (double a : cand) best = std::min(best, circDev(xtalAng, a));
          dev[si].push_back(best);
          wr[si] = std::max(wr[si], best);
          nAngles[si] += cand.size();
          ++nRotor[si];
        }
        double bg = 1e9;
        for (double a : grid30) bg = std::min(bg, circDev(xtalAng, a));
        dev[samps.size()].push_back(bg);
        wr[samps.size()] = std::max(wr[samps.size()], bg);
        nAngles[samps.size()] += grid30.size();
        ++nRotor[samps.size()];
      }
      for (size_t si = 0; si <= samps.size(); ++si)
        if (!dev[si].empty()) worst[si].push_back(wr[si]);
      ++done;
    } catch (...) {
      continue;
    }
  }

  auto pct = [](std::vector<double> v, double q) {
    if (v.empty()) return std::nan("");
    std::sort(v.begin(), v.end());
    size_t i = static_cast<size_t>(q * (v.size() - 1));
    return v[i];
  };
  auto frac = [](const std::vector<double> &v, double thr) {
    if (v.empty()) return std::nan("");
    return 100.0 *
           std::count_if(v.begin(), v.end(),
                         [&](double x) { return x <= thr; }) /
           v.size();
  };

  std::printf(
      "\n[torcheck] %zu molecules (rot %d..%s), inter-fragment junctions only\n",
      done, minRot,
      maxRot == std::numeric_limits<int>::max()
          ? "inf"
          : std::to_string(maxRot).c_str());
  std::printf(
      "  How far is the CRYSTAL dihedral from the NEAREST angle each sampler proposes?\n\n");
  std::printf("  %-13s %6s %7s %7s | %8s %8s %8s | %7s %8s\n", "sampler",
              "nrot", "med", "p90", "<=15deg", "<=30deg", "<=45deg", "ang/rot",
              "molOK30");
  for (size_t si = 0; si <= samps.size(); ++si) {
    if (dev[si].empty()) continue;
    const char *nm = si < samps.size() ? samps[si].name : "30deg grid ";
    std::printf(
        "  %-13s %6zu %7.1f %7.1f | %7.1f%% %7.1f%% %7.1f%% | %7.1f %7.1f%%\n",
        nm, dev[si].size(), pct(dev[si], 0.5), pct(dev[si], 0.9),
        frac(dev[si], 15.0), frac(dev[si], 30.0), frac(dev[si], 45.0),
        nRotor[si] ? double(nAngles[si]) / nRotor[si] : 0.0,
        frac(worst[si], 30.0));
  }
  std::printf(
      "\n  med/p90 = per-ROTOR deviation (deg).  molOK30 = %% of molecules whose EVERY\n"
      "  rotor is within 30deg of a proposed angle -- i.e. the crystal pose is\n"
      "  reachable at all.  ang/rot = candidate angles offered per rotor (cost).\n");
}

//! XTAL RECONSTRUCTION: feed each crystal ligand's EXACT fragment geometry (no
//! embed, no minimise -- ASM_EXACT_GEOM) through our real assembly pathway, and
//! ask two things:
//!  (1) REPRODUCTION -- can the search rebuild the crystal pose? (best-of-N
//!  heavy RMSD) (2) RANKING -- where does the TRUE crystal pose land in our
//!  scorers?  We score the
//!      crystal pose itself with the joiner's inter-fragment scorer AND full
//!      MMFF, and report what fraction of the generated ensemble our scorer
//!      prefers over it.  Since every fragment is EXACT crystal geometry, a
//!      poor crystal rank here is INTRINSIC (unrelaxed-bioactive strain), not a
//!      fragment-embedding artifact.  Per-term MMFF attribution (crystal vs our
//!      top pick) shows WHICH term buries the crystal.
//! The crystal is unrelaxed, so its absolute energy is high by construction --
//! the signal is the crystal's RANK vs the ensemble, not its absolute energy.
//! XTAL_RELAX=<kcal/A^2> does a position-restrained MMFF min first (shrug off
//! local strain, keep the pose) for cleaner fragment energies, per the
//! "flat-bottom well" idea.
void runXtalRecon() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = 40;
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);
  // rotor-count filter, same semantics as `bank` -- lets xtalrecon be pointed
  // straight at the high-rotor tail
  // instead of a rotor-mixed sample.
  int minRot = 0, maxRot = std::numeric_limits<int>::max();
  if (const char *e = std::getenv("PLATINUM_MINROT")) minRot = std::atoi(e);
  if (const char *e = std::getenv("PLATINUM_MAXROT")) maxRot = std::atoi(e);
  const std::string var = std::getenv("XTAL_MMFF94") ? "MMFF94" : "MMFF94s";
  const int seed = 0xf00d;
  const double relaxK =
      std::getenv("XTAL_RELAX") ? std::atof(std::getenv("XTAL_RELAX")) : 0.0;
  // XTAL_EMBED: apples-to-apples CONTROL -- run the SAME search/scorer/RMSD
  // code but with normal EMBEDDED fragments (ETKDG+MMFF) instead of exact
  // crystal geometry, to isolate the fragment-geometry contribution (expect
  // ~81% vs exact ~89%).
  const bool embed = std::getenv("XTAL_EMBED") != nullptr;
  if (!embed)
    setenv("ASM_EXACT_GEOM", "1", 1);  // else buildFragmentJoinerInput embeds

  // XTAL_OVERSAMPLE: run the warm-path fraglib with an oversampled pool +
  // RMSD-diverse selection.  ETKDG -> keep ETKDG coords, RANK by shrugged-MMFF
  // energy (FRAGLIB_SCORE_SHRUG); DG (XTAL_DG) -> MMFF-min each (loose 0.25
  // gradient) and keep the minimised coords.
  const bool oversample = embed && std::getenv("XTAL_OVERSAMPLE");
  const bool useDGtop = std::getenv("XTAL_DG") != nullptr;
  // XTAL_ETKDG_MIN: the ETKDG *minimizer* path -- ETKDG (knowledge-based) embed
  // THEN loose MMFF-min, keeping the minimised coords (vs the default keep-raw
  // + shrug-score).
  const bool etkdgMin = std::getenv("XTAL_ETKDG_MIN") != nullptr;
  const bool doMin = useDGtop || etkdgMin;
  std::shared_ptr<Fraglib> lib;
  if (oversample) {
    FraglibParams flp;
    flp.perClassEmbedding = false;  // xtalrecon uses an explicit flat pool/keep
    flp.fragmentEmbedMode =
        useDGtop ? FragmentEmbedMode::DG : FragmentEmbedMode::ETKDG;
    // DG or ETKDG-min -> Full (keep minimised coords); ETKDG keep-coords ->
    // ShrugScore.
    flp.minimizeMode =
        doMin ? FragmentMinimize::Full : FragmentMinimize::ShrugScore;
    flp.minimizeGradTol = 0.25;  // loose (no pucker collapse)
    flp.shrugDisplacement = 0.1;
    flp.numConfsPerFragment =
        std::getenv("XTAL_KEEP") ? std::atoi(std::getenv("XTAL_KEEP")) : 24;
    // perClassEmbedding is off -> everything embeds as the Exhaustive class.
    // Flat pool = XTAL_POOL (was flp.poolSize); keepN (maxConfs) stays the
    // Exhaustive default, matching the previous behaviour (the old
    // flp.numConfsPerFragment never affected the kept count).
    const int xtalPool =
        std::getenv("XTAL_POOL") ? std::atoi(std::getenv("XTAL_POOL")) : 1000;
    flp.setFlatPool(xtalPool);
    flp.randomSeed = seed;
    lib = std::make_shared<Fraglib>(flp);
    const char *label = useDGtop   ? "DG+MMFF(loose)"
                        : etkdgMin ? "ETKDG+MMFF(loose)"
                                   : "ETKDG+shrug-score";
    std::printf("  [oversample] %s pool=%d keep=%u %s\n", label, xtalPool,
                flp.numConfsPerFragment,
                doMin ? "(minimised coords)" : "(ETKDG coords kept)");
  }

  std::printf(
      "[xtalrecon] EXACT crystal fragments through our assembly pathway, "
      "%s, up to %zu mols%s\n",
      var.c_str(), maxMols, relaxK > 0 ? " (position-restrained pre-min)" : "");
  std::printf(
      "  reproduction = best-of-N heavy RMSD; ranking = crystal pose's place in "
      "our scorer over the generated ensemble\n\n");
  std::printf(
      "  %-24s %5s %5s | repro | interfrag-rank fullMMFF-rank | dVdW  dTors\n",
      "name", "rot", "nfrag");

  std::unique_ptr<SDMolSupplier> sup;
  try {
    sup = std::make_unique<SDMolSupplier>(sdf, true, false);
  } catch (...) {
    std::cerr << "xtalrecon: cannot open " << sdf << "\n";
    return;
  }
  RDLog::LogStateSetter blk;

  size_t done = 0, reproLt1 = 0, reproLt2 = 0, buriedByScore = 0, buriedByE = 0;
  double sumMinRms = 0, sumScorePct = 0, sumEpct = 0, sumDvdw = 0, sumDtors = 0;
  // capped best-of-K reproduction: keep only the top-K output confs by a
  // ranking, then best RMSD among them.  IF = inter-frag score ranking; FULL =
  // full-MMFF ranking.
  size_t capIF50 = 0, capFULL50 = 0, capIF25 = 0, capFULL25 = 0, capDIV50 = 0,
         capDIV25 = 0;
  while (!sup->atEnd() && done < maxMols) {
    std::unique_ptr<ROMol> raw;
    try {
      raw.reset(sup->next());
    } catch (...) {
      continue;
    }
    if (!raw || raw->getNumConformers() == 0) continue;
    try {
      MolOps::assignStereochemistryFrom3D(*raw);
      ROMOL_SPTR big = largestFrag(*raw);
      if (!big || big->getNumAtoms() < 4) continue;
      {  // rotor-count filter (PLATINUM_MINROT/MAXROT), applied before any
         // expensive work
        int nr = 0;
        try {
          nr = static_cast<int>(Descriptors::calcNumRotatableBonds(*big));
        } catch (...) {
        }
        if (nr < minRot || nr > maxRot) continue;
      }
      auto m = boost::make_shared<RWMol>(*big);
      MolOps::addHs(*m, /*explicitOnly=*/false, /*addCoords=*/true);
      if (relaxK > 0.0) {
        // XTAL_RELAX = #iterations of a light unrestrained MMFF min: shrug off
        // the worst unrelaxed strain while staying in the crystal's local well
        // (few steps
        // ~= a flat-bottom stay), for cleaner fragment energies -- per the
        // user's idea.
        MMFF::MMFFMolProperties props(*m, var);
        if (props.isValid()) {
          std::unique_ptr<ForceFields::ForceField> ff(
              MMFF::constructForceField(*m, &props, 1.0e8, 0));
          if (ff) {
            ff->initialize();
            ff->minimize(static_cast<int>(relaxK));
          }
        }
      }
      const unsigned int n = m->getNumAtoms();
      const std::string name =
          m->hasProp("_Name") ? m->getProp<std::string>("_Name") : "?";

      auto fullE = makeFullFFScoreFn(*m, false, var);
      if (!fullE) continue;
      std::vector<double> xtalBuf(3 * static_cast<size_t>(n));
      const Conformer &xc = m->getConformer();
      for (unsigned int a = 0; a < n; ++a) {
        const auto &p = xc.getAtomPos(a);
        xtalBuf[3 * a] = p.x;
        xtalBuf[3 * a + 1] = p.y;
        xtalBuf[3 * a + 2] = p.z;
      }
      const double Extal = fullE(xtalBuf.data(), n);
      const MmffTerms xt = mmffDecompose(*m, var);

      // XTAL_EMBED_RAW: embedded but NO fragment minimisation (raw ETKDG
      // geometry, whose own DG strain is closer to the crystal's -> tests the
      // relaxed-vs-strained mismatch). XTAL_DG: DG + MMFF-min fragments
      // instead of ETKDG.
      const bool useDG = std::getenv("XTAL_DG") != nullptr;
      const bool fragMin = embed && (useDG || !std::getenv("XTAL_EMBED_RAW"));
      // warm path (oversample lib) ignores the useDG/minMMFF args -- its
      // FraglibParams govern.
      FragmentJoinerInput in = buildFragmentJoinerInput(
          *m, embed ? 16 : 1, seed, var, useDG, /*minMMFF=*/fragMin, lib.get());
      if (!in.mol || in.fragments.size() < 2 || in.junctions.empty()) continue;

      // crystal pose's inter-fragment score (the joiner's ranking criterion)
      std::vector<std::pair<unsigned int, unsigned int>> jb;
      for (const auto &J : in.junctions) jb.emplace_back(J.atomA, J.atomB);
      auto iff = makeInterFragmentScoreFn(*in.mol, in.atomFragment, jb, false, var);
      const double Sxtal = iff ? iff(xtalBuf.data(), n)
                               : std::numeric_limits<double>::quiet_NaN();

      FragmentJoinerParams gp;
      gp.ffVariant = var;
      RigidRotorSearchParams sp;
      sp.fragConfBranch =
          embed ? 4 : 1;  // exact -> 1 conf/frag; embed -> sample pool
      sp.diversityRmsThresh =
          0.5;  // keep diverse poses incl. near-crystal
      sp.tree.beamWidth = 100;
      sp.defaultAngles = {-180, -150, -120, -90, -60, -30,
                               0,    30,   60,   90,  120, 150};  // 30deg
      sp.thompson.autoBudget = true;
      sp.thompson.maxBudget = 6000;
      auto res = runRigidRotorSearch(joinFragments(in, gp), sp).results;
      if (res.empty()) continue;

      ROMOL_SPTR xHeavy(MolOps::removeHs(static_cast<const ROMol &>(*m)));
      double minRms = std::numeric_limits<double>::max();
      int nearIdx = -1;
      double minEnsE = std::numeric_limits<double>::max();
      size_t betterScore = 0, betterE = 0;
      const MmffTerms bt = [&]() {  // best-SCORE conf term breakdown (res[0])
        auto pm = boost::make_shared<RWMol>(*in.mol);
        auto *c = new Conformer(n);
        for (unsigned int a = 0; a < n; ++a)
          c->setAtomPos(a, {res[0].coords[3 * a], res[0].coords[3 * a + 1],
                            res[0].coords[3 * a + 2]});
        pm->clearConformers();
        pm->addConformer(c, true);
        return mmffDecompose(*pm, var);
      }();
      // per-conf full-MMFF energy + heavy RMSD (res is already sorted by
      // inter-frag score)
      std::vector<double> rmsV(res.size(),
                               std::numeric_limits<double>::quiet_NaN());
      std::vector<double> eV(res.size(),
                             std::numeric_limits<double>::quiet_NaN());
      for (size_t k = 0; k < res.size(); ++k) {
        const double e = fullE(res[k].coords.data(), n);
        eV[k] = e;
        if (e < Extal) ++betterE;
        if (std::isfinite(res[k].score) && std::isfinite(Sxtal) &&
            res[k].score < Sxtal)
          ++betterScore;
        minEnsE = std::min(minEnsE, e);
        try {
          auto pm = boost::make_shared<RWMol>(*in.mol);
          auto *c = new Conformer(n);
          for (unsigned int a = 0; a < n; ++a)
            c->setAtomPos(a, {res[k].coords[3 * a], res[k].coords[3 * a + 1],
                              res[k].coords[3 * a + 2]});
          pm->clearConformers();
          pm->addConformer(c, true);
          ROMOL_SPTR ph(MolOps::removeHs(static_cast<const ROMol &>(*pm)));
          double r = MolAlign::getBestRMS(*ph, *xHeavy);
          rmsV[k] = r;
          if (r < minRms) {
            minRms = r;
            nearIdx = static_cast<int>(k);
          }
        } catch (...) {
        }
      }
      // capped best-of-K: top-K by inter-frag score (res order) vs top-K by
      // full MMFF
      std::vector<size_t> byE(res.size());
      std::iota(byE.begin(), byE.end(), 0);
      std::sort(byE.begin(), byE.end(),
                [&](size_t a, size_t b) { return eV[a] < eV[b]; });
      auto bestOf = [&](const std::vector<size_t> &ord, size_t K) {
        double b = std::numeric_limits<double>::max();
        for (size_t i = 0; i < std::min(K, ord.size()); ++i)
          if (std::isfinite(rmsV[ord[i]])) b = std::min(b, rmsV[ord[i]]);
        return b;
      };
      std::vector<size_t> byScore(res.size());
      std::iota(byScore.begin(), byScore.end(), 0);  // already inter-frag order
      const double if50 = bestOf(byScore, 50), full50 = bestOf(byE, 50);
      const double if25 = bestOf(byScore, 25), full25 = bestOf(byE, 25);
      if (if50 < 1.0) ++capIF50;
      if (full50 < 1.0) ++capFULL50;
      if (if25 < 1.0) ++capIF25;
      if (full25 < 1.0) ++capFULL25;
      // DIVERSITY-PRESERVING cap (no minimisation): walk poses in
      // inter-frag-score order, KEEP one only if it is > divThr heavy-atom
      // raw-RMSD from every already-kept pose, so the geometrically distinct
      // near-crystal pose survives even if it is not low-energy.
      const double divThr = std::getenv("XTAL_DIVTHR")
                                ? std::atof(std::getenv("XTAL_DIVTHR"))
                                : 1.0;
      auto rawRms = [&](size_t a, size_t b) {
        double s = 0;
        size_t nh = 0;
        for (unsigned int h = 0; h < n; ++h) {
          if (in.mol->getAtomWithIdx(h)->getAtomicNum() == 1) continue;
          double dx = res[a].coords[3 * h] - res[b].coords[3 * h];
          double dy = res[a].coords[3 * h + 1] - res[b].coords[3 * h + 1];
          double dz = res[a].coords[3 * h + 2] - res[b].coords[3 * h + 2];
          s += dx * dx + dy * dy + dz * dz;
          ++nh;
        }
        return nh ? std::sqrt(s / nh) : 0.0;
      };
      auto bestOfDiverse = [&](size_t K) {
        std::vector<size_t> kept;
        for (size_t idx : byScore) {  // energy(inter-frag)-ordered candidates
          bool distinct = true;
          for (size_t j : kept)
            if (rawRms(idx, j) < divThr) {
              distinct = false;
              break;
            }
          if (distinct) kept.push_back(idx);
          if (kept.size() >= K) break;
        }
        double b = std::numeric_limits<double>::max();
        for (size_t idx : kept)
          if (std::isfinite(rmsV[idx])) b = std::min(b, rmsV[idx]);
        return b;
      };
      const double div50 = bestOfDiverse(50), div25 = bestOfDiverse(25);
      if (div50 < 1.0) ++capDIV50;
      if (div25 < 1.0) ++capDIV25;

      unsigned int nrot = 0;
      try {
        nrot = Descriptors::calcNumRotatableBonds(*xHeavy);
      } catch (...) {
      }
      const double scorePct = 100.0 * betterScore / res.size();
      const double ePct = 100.0 * betterE / res.size();
      const double dVdw = bt.vdw - xt.vdw, dTors = bt.tors - xt.tors;
      std::printf(
          "  %-24s %5u %5zu | %.2f%s | %5.0f%% (n=%-4zu) %5.0f%%       | "
          "%+5.1f %+5.1f\n",
          name.substr(0, 24).c_str(), nrot, in.fragments.size(), minRms,
          minRms < 1.0 ? "*" : " ", scorePct, res.size(), ePct, dVdw, dTors);
      std::fflush(stdout);

      ++done;
      if (minRms < 1.0) ++reproLt1;
      if (minRms < 2.0) ++reproLt2;
      if (scorePct > 50.0) ++buriedByScore;
      if (ePct > 50.0) ++buriedByE;
      sumMinRms += minRms;
      sumScorePct += scorePct;
      sumEpct += ePct;
      sumDvdw += dVdw;
      sumDtors += dTors;
    } catch (...) {
      continue;
    }
  }

  if (!done) {
    std::printf("\n  no molecules scored\n");
    return;
  }
  std::printf("\n=== summary over %zu mols ===\n", done);
  std::printf(
      "  REPRODUCTION (exact frags): best-of-N RMSD mean=%.2f  %%<1A=%.0f  %%<2A=%.0f\n",
      sumMinRms / done, 100.0 * reproLt1 / done, 100.0 * reproLt2 / done);
  std::printf(
      "  RANKING of the true crystal pose (fraction of ensemble our scorer prefers):\n");
  std::printf(
      "    inter-frag scorer: mean=%.0f%%  (crystal buried >50%% in %.0f%% of mols)\n",
      sumScorePct / done, 100.0 * buriedByScore / done);
  std::printf(
      "    full MMFF (%s):    mean=%.0f%%  (crystal buried >50%% in %.0f%% of mols)\n",
      var.c_str(), sumEpct / done, 100.0 * buriedByE / done);
  std::printf(
      "  ENERGY ATTRIBUTION (our top pick - crystal, kcal): dVdW mean=%+.1f  dTors mean=%+.1f\n",
      sumDvdw / done, sumDtors / done);
  std::printf(
      "  (negative dVdW/dTors => our top-ranked pose has LOWER vdW/torsion than the\n"
      "   crystal, i.e. the scorer prefers a non-crystal pose -- ranking, not geometry)\n");
  std::printf(
      "\n  CAPPED best-of-K reproduction (%%<1A) -- keep top-K output confs, then best RMSD:\n");
  std::printf(
      "    top-25:  inter-frag rank=%.0f%%   full-MMFF rank=%.0f%%   DIVERSE=%.0f%%\n",
      100.0 * capIF25 / done, 100.0 * capFULL25 / done,
      100.0 * capDIV25 / done);
  std::printf(
      "    top-50:  inter-frag rank=%.0f%%   full-MMFF rank=%.0f%%   DIVERSE=%.0f%%\n",
      100.0 * capIF50 / done, 100.0 * capFULL50 / done,
      100.0 * capDIV50 / done);
  std::printf(
      "    (uncapped best-of-N was %.0f%%; the gap to it is what OUTPUT RANKING costs)\n",
      100.0 * reproLt1 / done);
}

//! Per-rotor-count timing crossover: run ETKDG and TorLib samplers, each with
//! BOTH the deterministic tree/beam search AND Thompson sampling, over the full
//! Platinum set, and bucket ASSEMBLY-ONLY wall time (embed excluded,
//! per-molecule) by the FRAGCG_PARAMS=<file>: load a serialized
//! FragmentConfGenParams text file (key = value; see Utils/ParamsIO) and layer
//! it over each run's params, so any knob can be driven from a config file
//! instead of a bespoke env var.  Returns the file text (empty if unset);
//! validates once up front and hard-exits on a parse error so a typo can't
//! silently no-op.
std::string loadParamsOverrideText() {
  const char *f = std::getenv("FRAGCG_PARAMS");
  if (!f) return {};
  std::ifstream in(f);
  if (!in) {
    std::fprintf(stderr, "FRAGCG_PARAMS: cannot open '%s'\n", f);
    std::exit(1);
  }
  std::stringstream ss;
  ss << in.rdbuf();
  std::string text = ss.str();
  RDKit::FragmentConfGenParams probe;
  std::string err = RDKit::fragmentConfGenParamsFromString(text, probe);
  if (!err.empty()) {
    std::fprintf(stderr, "FRAGCG_PARAMS('%s'): %s\n", f, err.c_str());
    std::exit(1);
  }
  // Parsing only rejects UNKNOWN KEYS; a known key with an out-of-range or
  // misspelled VALUE parses fine and would otherwise run as the default -- i.e.
  // silently as the control arm of whatever experiment this file encodes.
  err = probe.validate();
  if (!err.empty()) {
    std::fprintf(stderr, "FRAGCG_PARAMS('%s'): %s\n", f, err.c_str());
    std::exit(1);
  }
  return text;
}

//! molecule's rotatable-bond count.  The point is to see whether there is a
//! rotor count at which TS overtakes the exhaustive beam (or vice versa) on
//! speed, and what it costs in accuracy.  Env: PLATINUM_MAXMOLS (default all),
//! PLATINUM_THREADS (8), PLATINUM_SDF, ASM_SEED, ASM_MMFF94 (else MMFF94s),
//! ROTORCROSS_MAXBIN (tail bin, 14), FRAGCG_PARAMS (serialized-params override
//! file).
void runRotorCross() {
  const std::string paramsOverride = loadParamsOverrideText();
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = std::numeric_limits<size_t>::max();
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);
  unsigned int threads = 8;
  if (const char *e = std::getenv("PLATINUM_THREADS")) threads = std::atoi(e);
  int seed = 0xf00d;
  if (const char *e = std::getenv("ASM_SEED")) seed = std::atoi(e);
  int maxBin = 14;  // rotor counts >= this fold into one tail bin
  if (const char *e = std::getenv("ROTORCROSS_MAXBIN")) maxBin = std::atoi(e);
  const std::string mmff = std::getenv("ASM_MMFF94") ? "MMFF94" : "MMFF94s";

  // Load mols + precompute each one's rotatable-bond count ONCE (shared across
  // configs).
  std::vector<ROMOL_SPTR> jobs;
  std::vector<int> rotOf;
  {
    std::unique_ptr<SDMolSupplier> sup;
    try {
      sup = std::make_unique<SDMolSupplier>(sdf, true, false);
    } catch (...) {
      std::cerr << "rotorcross: cannot open " << sdf << "\n";
      return;
    }
    RDLog::LogStateSetter blk;
    while (!sup->atEnd() && jobs.size() < maxMols) {
      std::unique_ptr<ROMol> raw;
      try {
        raw.reset(sup->next());
      } catch (...) {
        continue;
      }
      if (!raw || raw->getNumConformers() == 0) continue;
      try {
        MolOps::assignStereochemistryFrom3D(*raw);
        std::unique_ptr<ROMol> heavy(MolOps::removeHs(*raw));
        ROMOL_SPTR refHeavy = largestFrag(*heavy);
        if (!refHeavy || refHeavy->getNumAtoms() < 2) continue;
        int nrot = 0;
        try {
          nrot =
              static_cast<int>(Descriptors::calcNumRotatableBonds(*refHeavy));
        } catch (...) {
        }
        jobs.push_back(refHeavy);
        rotOf.push_back(nrot);
      } catch (...) {
      }
    }
  }

  // The four configs: {sampler} x {search mode}.  TorLib may be unavailable
  // (needs the Hamburg XML) -- skip it gracefully.
  struct Cfg {
    std::string name;
    std::shared_ptr<TorsionSampler> sampler;
    bool exhaustive;  // true = deterministic tree/beam; false = Thompson (auto
                      // budget)
  };
  std::vector<Cfg> cfgs;
  cfgs.push_back({"ETKDG-tree", std::make_shared<ETKDGTorsionSampler>(), true});
  cfgs.push_back({"ETKDG-TS", std::make_shared<ETKDGTorsionSampler>(), false});
  try {
    auto s1 = hamburgSampler();
    auto s2 = hamburgSampler();
    cfgs.push_back({"TorLib-tree", s1, true});
    cfgs.push_back({"TorLib-TS", s2, false});
  } catch (const std::exception &e) {
    std::cerr << "rotorcross: torsion library unavailable (" << e.what()
              << "); running ETKDG configs only\n";
  }
  if (std::getenv("FCG_SMIRNOFF")) {  // external CC-BY OpenFF .offxml
    cfgs.push_back({"SMIRNOFF-tree", smirnoffSampler(), true});
    cfgs.push_back({"SMIRNOFF-TS", smirnoffSampler(), false});
  }
  // ASM_ONLY=<name>[,<name>...] keeps only the listed configs (e.g. ETKDG-tree)
  // -- lets a sweep run one config instead of all four (each ASM_SEARCH
  // override applies to all).
  if (const char *only = std::getenv("ASM_ONLY")) {
    const std::string spec(only);
    std::vector<Cfg> keep;
    for (auto &c : cfgs) {
      if (spec.find(c.name) != std::string::npos) keep.push_back(c);
    }
    if (!keep.empty()) cfgs.swap(keep);
  }

  // Per-rotor-bin accumulator (bin = min(rot, maxBin)).
  struct Bin {
    size_t n = 0, lt1 = 0, lt2 = 0;
    long long asmNs = 0;  // summed assembly-only wall ns
    double genN = 0.0;
    long long budget = 0;
  };

  std::printf(
      "[rotorcross] assembly-only wall time by rotatable-bond count, %zu mols, "
      "%s, %u threads\n"
      "  (embed excluded; ms = per-mol assembly wall time; tail bin = rot>=%d)\n\n",
      jobs.size(), mmff.c_str(), threads, maxBin);

  std::map<std::string, std::map<int, Bin>>
      results;  // cfg name -> bin -> stats

  for (auto &cfg : cfgs) {
    std::map<int, Bin> bins;
    std::mutex mtx;
    std::atomic<size_t> next{0};
    auto worker = [&]() {
      RDLog::LogStateSetter blk;
      auto localSampler =
          cfg.sampler ? cfg.sampler->copy() : std::shared_ptr<TorsionSampler>();
      while (true) {
        size_t i = next++;
        if (i >= jobs.size()) break;
        double rms = std::numeric_limits<double>::quiet_NaN();
        size_t nc = 0;
        long long asmNs = 0;
        unsigned int bud = 0;
        try {
          FragmentConfGenParams pp;
          pp.numOutputConfs = 2000;
          pp.randomSeed = seed;
          pp.search.torsionSampler = localSampler;
          pp.search.thompsonBudget = 0;
          pp.search.thompson.autoBudget = !cfg.exhaustive;
          pp.search.searchMode = searchModeFromEnv();  // ASM_SEARCH override
          pp.joiner.ffVariant = mmff;
          // fragment-conformer coverage sweep (fragment-state-starvation test)
          if (const char *e = std::getenv("ASM_BRANCH"))
            pp.search.fragConfBranch = std::atoi(e);
          if (const char *e = std::getenv("ASM_ROOTSEEDS"))
            pp.search.rootSeeds = std::atoi(e);
          // output-retention granularity: finer diversityRms keeps MORE
          // distinct confs (tests whether the 0.5A dedup discards good poses vs
          // the pool lacking them).
          if (const char *e = std::getenv("ASM_DIVRMS"))
            pp.search.diversityRmsThresh = std::atof(e);
          // FRAGCG_PARAMS overrides win over the bench defaults above (already
          // validated; sampler pointer is not serialized, so localSampler
          // survives).
          if (!paramsOverride.empty()) {
            RDKit::fragmentConfGenParamsFromString(paramsOverride, pp);
          }
          requireValidParams(pp, "runRotorCross");
          auto rr = FragmentConfGen(pp).build(*jobs[i]);
          nc = rr.conformers.size();
          bud = rr.joinerBudget;
          asmNs = rr.joinerAssemblyNs;
          rms = minHeavyRms(rr.conformers, *jobs[i]);
        } catch (...) {
        }
        const int bin = std::min(rotOf[i], maxBin);
        std::lock_guard<std::mutex> lk(mtx);
        Bin &b = bins[bin];
        ++b.n;
        b.asmNs += asmNs;
        b.genN += static_cast<double>(nc);
        b.budget += bud;
        if (std::isfinite(rms)) {
          if (rms < 1.0) ++b.lt1;
          if (rms < 2.0) ++b.lt2;
        }
      }
    };
    std::vector<std::thread> pool;
    for (unsigned int t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto &t : pool) t.join();
    results[cfg.name] = bins;

    // per-config table
    std::printf("  == %-12s ==\n    %4s %6s %10s %7s %7s %8s %9s\n",
                cfg.name.c_str(), "rot", "nMol", "msAsm/mol", "%<1", "%<2",
                "genN", "avgBudg");
    long long totNs = 0;
    size_t totN = 0, totLt1 = 0, totLt2 = 0;
    for (auto &kv : bins) {
      const Bin &b = kv.second;
      const double denom = b.n ? b.n : 1;
      std::printf("    %3d%s %6zu %10.3f %6.1f%% %6.1f%% %8.1f %9.0f\n",
                  kv.first, kv.first >= maxBin ? "+" : " ", b.n,
                  (b.asmNs / 1.0e6) / denom, 100.0 * b.lt1 / denom,
                  100.0 * b.lt2 / denom, b.genN / denom, b.budget / denom);
      totNs += b.asmNs;
      totN += b.n;
      totLt1 += b.lt1;
      totLt2 += b.lt2;
    }
    const double dN = totN ? totN : 1;
    std::printf("    %-4s %6zu %10.3f %6.1f%% %6.1f%%\n\n", "ALL", totN,
                (totNs / 1.0e6) / dN, 100.0 * totLt1 / dN, 100.0 * totLt2 / dN);
    std::fflush(stdout);
  }

  // Crossover summary: tree-vs-TS ms per rotor bin for each sampler family,
  // marking where TS becomes the faster of the two.
  auto crossover = [&](const std::string &treeName, const std::string &tsName) {
    if (!results.count(treeName) || !results.count(tsName)) return;
    std::printf(
        "  == crossover %s vs %s (assembly ms/mol; * = TS faster) ==\n"
        "    %4s %10s %10s %9s\n",
        treeName.c_str(), tsName.c_str(), "rot", "tree ms", "TS ms", "TS/tree");
    const auto &tb = results[treeName];
    const auto &sb = results[tsName];
    std::set<int> allBins;
    for (auto &kv : tb) allBins.insert(kv.first);
    for (auto &kv : sb) allBins.insert(kv.first);
    for (int bin : allBins) {
      double treeMs = std::numeric_limits<double>::quiet_NaN();
      double tsMs = std::numeric_limits<double>::quiet_NaN();
      if (tb.count(bin) && tb.at(bin).n)
        treeMs = (tb.at(bin).asmNs / 1.0e6) / tb.at(bin).n;
      if (sb.count(bin) && sb.at(bin).n)
        tsMs = (sb.at(bin).asmNs / 1.0e6) / sb.at(bin).n;
      const double ratio = (std::isfinite(treeMs) && treeMs > 0)
                               ? tsMs / treeMs
                               : std::numeric_limits<double>::quiet_NaN();
      const bool tsFaster = std::isfinite(ratio) && ratio < 1.0;
      std::printf("    %3d%s %10.3f %10.3f %8.2fx%s\n", bin,
                  bin >= maxBin ? "+" : " ", treeMs, tsMs, ratio,
                  tsFaster ? " *" : "");
    }
    std::printf("\n");
    std::fflush(stdout);
  };
  crossover("ETKDG-tree", "ETKDG-TS");
  crossover("TorLib-tree", "TorLib-TS");
}

//! Single-point MMFF energy (electrostatics OFF), variant-aware, for a
//! conformer that already carries explicit Hs + 3D coords.  NaN if untypeable.
double mmffSinglePoint(const ROMol &m, const std::string &variant) {
  try {
    RWMol rw(m);
    MMFF::MMFFMolProperties props(rw, variant);
    if (!props.isValid()) return std::numeric_limits<double>::quiet_NaN();
    props.setMMFFEleTerm(false);
    std::unique_ptr<ForceFields::ForceField> ff(
        MMFF::constructForceField(rw, &props, 1.0e8, 0));
    if (!ff) return std::numeric_limits<double>::quiet_NaN();
    ff->initialize();
    return ff->calcEnergy();
  } catch (...) {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

//! Per-molecule placement of the bioactive (crystal) pose inside an
//! energy-ranked ensemble: for each conformer its energy + heavy-atom RMSD to
//! the crystal, an energy rank (0 = lowest energy), and the ranks at which a <1
//! A / <2 A pose first appears.  rankLt* = -1 when no conformer meets that
//! threshold.
struct BankStat {
  int nConfs = 0;
  double bestRms = std::numeric_limits<double>::quiet_NaN();
  int bestRank = -1;  //!< energy rank of the closest-to-crystal conformer
  int rankLt1 = -1;   //!< energy rank of the lowest-energy conf with RMSD < 1 A
  int rankLt2 = -1;   //!< energy rank of the lowest-energy conf with RMSD < 2 A
  std::vector<double> energy;  //!< per input conf (NaN = untypeable)
  std::vector<double> rms;     //!< per input conf, heavy RMSD to crystal
  std::vector<int> rank;       //!< per input conf, energy rank
};

BankStat computeBankStat(const std::vector<ROMOL_SPTR> &ens,
                         const ROMol &refHeavy,
                         const std::vector<double> &energies) {
  BankStat s;
  s.nConfs = static_cast<int>(ens.size());
  s.energy = energies;
  s.rms.assign(ens.size(), std::numeric_limits<double>::quiet_NaN());
  s.rank.assign(ens.size(), -1);
  RWMol ref(refHeavy);
  for (size_t i = 0; i < ens.size(); ++i) {
    const auto &g = ens[i];
    if (!g || g->getNumConformers() == 0) continue;
    ROMOL_SPTR probe;
    try {
      probe.reset(MolOps::removeHs(*g));
    } catch (...) {
      continue;
    }
    if (probe->getNumAtoms() != ref.getNumAtoms()) continue;
    try {
      s.rms[i] = MolAlign::getBestRMS(*probe, ref);
    } catch (...) {
    }
  }
  // energy rank (ascending energy; NaN energies sort last, keeping input order)
  std::vector<size_t> order(ens.size());
  for (size_t i = 0; i < ens.size(); ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    const double ea = s.energy[a], eb = s.energy[b];
    if (std::isnan(ea)) return false;
    if (std::isnan(eb)) return true;
    return ea < eb;
  });
  for (size_t k = 0; k < order.size(); ++k)
    s.rank[order[k]] = static_cast<int>(k);
  // best RMSD + threshold ranks
  double best = std::numeric_limits<double>::max();
  for (size_t i = 0; i < ens.size(); ++i) {
    if (!std::isfinite(s.rms[i])) continue;
    if (s.rms[i] < best) {
      best = s.rms[i];
      s.bestRank = s.rank[i];
    }
    if (s.rms[i] < 1.0 && (s.rankLt1 < 0 || s.rank[i] < s.rankLt1))
      s.rankLt1 = s.rank[i];
    if (s.rms[i] < 2.0 && (s.rankLt2 < 0 || s.rank[i] < s.rankLt2))
      s.rankLt2 = s.rank[i];
  }
  if (std::isfinite(best)) s.bestRms = best;
  return s;
}

//! Bank an OURS(tree search) run over the full Platinum
//! set: writes per-method SDF ensembles (each conformer tagged with MMFF energy
//! + energy rank; each molecule tagged with #rotors, gen time, best RMSD and
//! the crystal-pose placement at <1/<2 A) plus a per-molecule stats TSV (->
//! xlsx offline). Env: PLATINUM_SDF, PLATINUM_MAXMOLS (all), PLATINUM_THREADS
  //! (8), ASM_SEED, BANK_OUTDIR (scratch), BANK_MAXCONFS (250),
  //! ASM_MMFF94 (else MMFF94s), BANK_TS (also run our
//! Thompson config as a 3rd method).
void runBank() {
  std::string sdf = testDataPath(
      "Code/GraphMol/test_data/platinum_diverse_dataset_2017_01.sdf");
  if (const char *e = std::getenv("PLATINUM_SDF")) sdf = e;
  size_t maxMols = std::numeric_limits<size_t>::max();
  if (const char *e = std::getenv("PLATINUM_MAXMOLS")) maxMols = std::atoi(e);
  //! PLATINUM_FROM: skip the first N ACCEPTED mols, so a second,
  //! non-overlapping set can be carved out (e.g. FROM=300 MAXMOLS=300 -> the
  //! 2nd 300-mol set, disjoint from the 1st).
  size_t from = 0;
  if (const char *e = std::getenv("PLATINUM_FROM")) from = std::atoi(e);
  //! Rotor-count filter: iterate ONLY on the regime under investigation (e.g.
    //! the rot>=7 tail, ~622 of 2859 mols) instead of paying for a full-set
    //! run to get a handful of
  //! relevant molecules.
  int minRot = 0, maxRot = std::numeric_limits<int>::max();
  if (const char *e = std::getenv("PLATINUM_MINROT")) minRot = std::atoi(e);
  if (const char *e = std::getenv("PLATINUM_MAXROT")) maxRot = std::atoi(e);
  unsigned int threads = 8;
  if (const char *e = std::getenv("PLATINUM_THREADS")) threads = std::atoi(e);
  int seed = 0xf00d;
  if (const char *e = std::getenv("ASM_SEED")) seed = std::atoi(e);
  // SIGNED: -1 = AUTO (rotor-driven, see autoOutputConfs), 0 = DISABLED.
  // Was unsigned, which turned BANK_MAXCONFS=-1 into 4294967295 -- AUTO could
  // not be benchmarked at all.
  int maxConfs = 250;
  if (const char *e = std::getenv("BANK_MAXCONFS")) maxConfs = std::atoi(e);
  const std::string mmff = std::getenv("ASM_MMFF94") ? "MMFF94" : "MMFF94s";
  const bool doTS = std::getenv("BANK_TS") != nullptr;
  //! BANK_TORLIB: use the Hamburg TorLib sampler for OURS instead of ETKDG
  //! (embedding is still ETKDG+MMFF; this only changes the junction-angle
  //! sampler).  ASM_SEARCH selects the search (e.g. systematic).
  const bool useTorlib = std::getenv("BANK_TORLIB") != nullptr;
  const std::string paramsOverride =
      loadParamsOverrideText();  // FRAGCG_PARAMS override file
  // ONE shared fragment cache for all worker threads (Fraglib is thread-safe:
  // mutex-guarded map, masters immutable once inserted).  Without this each of
  // the N threads builds its own private cache -> N copies of every embedded
  // fragment -> memory blows up under parallel runs (and no cross-thread
  // reuse).  Default ETKDG+MMFF embedding at the run seed. The shared lib MUST
  // embed with the SAME params the per-molecule pp uses, or a mismatch (e.g.
  // FRAGCG_PARAMS bumps pp.embedding.poolSize but the shared lib stays at the
  // default) makes the pipeline serve 0 confs.  Mirror the FRAGCG_PARAMS
  // override into sharedFlp by routing it through a FragmentConfGenParams and
  // taking its .embedding.
  FragmentConfGenParams sharedTmpl;
  sharedTmpl.randomSeed = seed;
  sharedTmpl.joiner.ffVariant = mmff;  // as the per-molecule pp below sets it
  if (!paramsOverride.empty())
    RDKit::fragmentConfGenParamsFromString(paramsOverride, sharedTmpl);
  FraglibParams sharedFlp = sharedTmpl.embedding;
  // getFraglibParams() drives the embedder's variant from the joiner's (one
  // force field for the whole pipeline).  This path builds the Fraglib
  // directly, so apply the same rule -- otherwise the shared cache embeds under
  // a different force field than every per-molecule run assumes.
  sharedFlp.ffVariant = sharedTmpl.joiner.ffVariant;
  sharedFlp.randomSeed = seed;
  auto sharedLib = std::make_shared<Fraglib>(sharedFlp);
  // BANK_FRAGLIB: preload a CANNED library built offline by genFragLib -- the
  // shipping warm-cache scenario, as opposed to BANK_SHARE's within-run
  // warming.  Implies BANK_SHARE (a preloaded library is useless unpassed).
  // The embedding must match or FragmentConfGenParams::validate() rejects it,
  // which is the point: a library built under different settings serves
  // fragments this run does not mean.
  bool haveCannedLib = false;
  if (const char *e = std::getenv("BANK_FRAGLIB")) {
    std::ifstream in(e, std::ios::binary);
    if (!in) {
      std::cerr << "[bank] cannot open BANK_FRAGLIB " << e << "\n";
      return;
    }
    try {
      sharedLib->initFromStream(in);
      haveCannedLib = true;
      std::printf("[bank] BANK_FRAGLIB %s -> %zu fragments preloaded\n", e,
                  sharedLib->size());
    } catch (const std::exception &ex) {
      std::cerr << "[bank] failed to load BANK_FRAGLIB: " << ex.what() << "\n";
      return;
    }
  }
  std::string outdir = "./bank_results";
  if (const char *e = std::getenv("BANK_OUTDIR")) outdir = e;
  std::system(("mkdir -p '" + outdir + "'").c_str());

  // Load mols + crystal refs + rotor counts + SMILES.
  struct Job {
    ROMOL_SPTR ref;  // largest heavy component, crystal conformer
    int nrot = 0;
    std::string name, smi;
  };
  std::vector<Job> jobs;
  {
    std::unique_ptr<SDMolSupplier> sup;
    try {
      sup = std::make_unique<SDMolSupplier>(sdf, true, false);
    } catch (...) {
      std::cerr << "bank: cannot open " << sdf << "\n";
      return;
    }
    RDLog::LogStateSetter blk;
    size_t idx = 0;
    size_t accepted =
        0;  // filter-passing mols seen so far (for the PLATINUM_FROM offset)
    while (!sup->atEnd() && jobs.size() < maxMols) {
      std::unique_ptr<ROMol> raw;
      try {
        raw.reset(sup->next());
      } catch (...) {
        continue;
      }
      ++idx;
      if (!raw || raw->getNumConformers() == 0) continue;
      try {
        MolOps::assignStereochemistryFrom3D(*raw);
        std::string nm = raw->hasProp("_Name")
                             ? raw->getProp<std::string>("_Name")
                             : ("mol" + std::to_string(idx));
        std::unique_ptr<ROMol> heavy(MolOps::removeHs(*raw));
        ROMOL_SPTR refHeavy = largestFrag(*heavy);
        if (!refHeavy || refHeavy->getNumAtoms() < 2) continue;
        Job j;
        j.ref = refHeavy;
        try {
          j.nrot =
              static_cast<int>(Descriptors::calcNumRotatableBonds(*refHeavy));
        } catch (...) {
        }
        try {
          j.smi = MolToSmiles(*refHeavy);
        } catch (...) {
        }
        if (nm.empty()) nm = "mol" + std::to_string(idx);
        j.name = nm;
        if (j.nrot < minRot || j.nrot > maxRot) continue;  // rotor-count filter
        if (accepted++ < from) continue;  // skip the first `from` accepted mols
        jobs.push_back(std::move(j));
      } catch (...) {
      }
    }
  }

  {
    const char *sm = std::getenv("ASM_SEARCH");
    std::printf(
        "[bank] %zu mols (from %zu, rot %d-%d), %s, ours=%s %s%s, up to %d confs, %u threads\n"
        "  outdir=%s\n\n",
        jobs.size(), from, minRot,
        maxRot == std::numeric_limits<int>::max() ? 99 : maxRot, mmff.c_str(),
        useTorlib ? "TorLib" : "ETKDG", sm ? sm : "auto",
        doTS ? "+Thompson" : "",
        maxConfs,
        threads, outdir.c_str());
  }

  // Method-tagged SDF writers + a stats accumulator, all mutex-guarded.
  std::map<std::string, std::unique_ptr<SDWriter>> writers;
  auto mkWriter = [&](const std::string &m) {
    writers[m] = std::make_unique<SDWriter>(outdir + "/" + m + ".sdf");
  };
  mkWriter("ours_tree");
  if (doTS) mkWriter("ours_ts");
  std::mutex ioMtx;
  std::vector<std::string> statRows;
  std::atomic<size_t> next{0}, done{0};

  auto emit = [&](const std::string &method, const Job &j,
                  const std::vector<ROMOL_SPTR> &ens,
                  const std::vector<double> &energies, double genMs) {
    BankStat s = computeBankStat(ens, *j.ref, energies);
    std::lock_guard<std::mutex> lk(ioMtx);
    auto &w = writers[method];
    for (size_t i = 0; i < ens.size(); ++i) {
      if (!ens[i] || ens[i]->getNumConformers() == 0) continue;
      RWMol c(*ens[i]);
      c.setProp("_Name", j.name + "_" + method + "_" + std::to_string(i));
      if (std::isfinite(s.energy[i])) c.setProp("mmff_energy", s.energy[i]);
      c.setProp("energy_rank", s.rank[i]);
      if (std::isfinite(s.rms[i])) c.setProp("rms_to_crystal", s.rms[i]);
      c.setProp("n_rotors", j.nrot);
      c.setProp("gen_ms", genMs);
      w->write(c);
    }
    w->flush();
    char row[512];
    std::snprintf(row, sizeof(row),
                  "%s\t%s\t%d\t%d\t%.2f\t%.3f\t%d\t%d\t%d\t%d\t%d",
                  method.c_str(), j.name.c_str(), j.nrot, s.nConfs, genMs,
                  s.bestRms, s.bestRank, s.rankLt1, s.rankLt2,
                  (std::isfinite(s.bestRms) && s.bestRms < 1.0) ? 1 : 0,
                  (std::isfinite(s.bestRms) && s.bestRms < 2.0) ? 1 : 0);
    statRows.emplace_back(row);
  };

  auto worker = [&]() {
    RDLog::LogStateSetter blk;
    // Per-worker sampler (ETKDG default, or Hamburg TorLib under BANK_TORLIB).
    std::shared_ptr<TorsionSampler> sampler;
    if (useTorlib) {
      sampler = hamburgSampler();
    } else {
      auto es = std::make_shared<ETKDGTorsionSampler>();
      // A/B: symmetry-fold the ETKDG angle set (TorLib folds; this sampler
      // never did)
      // ETKDG_FOLDSYM=0 disables, =1 enables; unset keeps the class default
      // (ON).  It must be able to turn folding OFF now that ON is the default,
      // or the A/B cannot be re-run.
      if (const char *e = std::getenv("ETKDG_FOLDSYM"))
        es->foldBySymmetry = (std::atoi(e) != 0);
      sampler = es;
    }
    const bool traceBank = std::getenv("BANK_TRACE") != nullptr;
    while (true) {
      size_t i = next++;
      if (i >= jobs.size()) break;
      const Job &j = jobs[i];
      if (traceBank) {
        std::fprintf(stderr, "[bank] >>> %zu %s  nrot=%d  smi=%s\n", i,
                     j.name.c_str(), j.nrot, j.smi.c_str());
        std::fflush(stderr);
      }

      // OURS -- tree search (deterministic beam), full pipeline wall time.
      auto runOurs = [&](bool ts, const std::string &method) {
        std::vector<ROMOL_SPTR> ens;
        std::vector<double> en;
        double ms = 0.0;
        try {
          FragmentConfGenParams pp;
          pp.numOutputConfs = maxConfs;
          pp.randomSeed = seed;
          if (haveCannedLib || std::getenv("BANK_SHARE"))
            pp.fraglib = sharedLib;  // opt-in: cross-molecule reuse (SPEED), at
                                     // the cost of
          // accumulating every unique fragment for the whole run (MORE memory,
          // not less -- the default per-build cache is per-molecule and tiny).
          pp.search.torsionSampler = sampler;
          pp.search.thompsonBudget = 0;
          pp.search.thompson.autoBudget = ts;          // false = tree/beam
          pp.search.searchMode = searchModeFromEnv();  // ASM_SEARCH override
          pp.joiner.ffVariant = mmff;
          if (!paramsOverride.empty())
            RDKit::fragmentConfGenParamsFromString(paramsOverride,
                                                   pp);  // FRAGCG_PARAMS
          // classParams is a map-of-structs (not in the FRAGCG_PARAMS text
          // format): inject the Acyclic recipe from env so the
          // whole-acyclic/driving experiments can size the pool and kept-conf
          // count for fragments that now carry real internal rotors.
          // ACYCLIC_KEEP = maxConfs (keepN), ACYCLIC_PERROTOR = pool as
          // confsPerRotor*rotors. ACYCLIC_EMBED=dg|etkdg,
          // ACYCLIC_MINMODE=Full|Score|ShrugScore|None, ACYCLIC_ITERS=N
          {
            const char *keep = std::getenv("ACYCLIC_KEEP");
            const char *perRot = std::getenv("ACYCLIC_PERROTOR");
            const char *emb = std::getenv("ACYCLIC_EMBED");
            const char *mm = std::getenv("ACYCLIC_MINMODE");
            const char *iters = std::getenv("ACYCLIC_ITERS");
            if (keep || perRot || emb || mm || iters) {
              FragmentParams ac =
                  pp.embedding.classParams[FragmentClass::Acyclic];
              if (keep) ac.maxConfs = std::atoi(keep);
              if (perRot) {
                ac.op = FragmentSamplesOperator::MULTIPLY;
                ac.confsPerRotor = std::atoi(perRot);
                ac.maxSamples = std::max(ac.maxSamples, ac.confsPerRotor * 12);
              }
              if (emb) {
                ac.embedMode =
                    (std::string(emb) == "dg" || std::string(emb) == "DG")
                        ? FragmentEmbedMode::DG
                        : FragmentEmbedMode::ETKDG;
              }
              if (mm) {
                const std::string v(mm);
                ac.minimizeMode = v == "None"    ? FragmentMinimize::None
                                  : v == "Score" ? FragmentMinimize::Score
                                  : v == "ShrugScore"
                                      ? FragmentMinimize::ShrugScore
                                      : FragmentMinimize::Full;
              }
              if (iters) ac.minimizeMaxIters = std::atoi(iters);
              pp.embedding.setClassParam(FragmentClass::Acyclic, ac);
            }
          }
          // RIGID_POOL / RIGID_KEEP: rigid (aromatic/fused) fragments dominate
          // embed cost -- the default recipe embeds maxSamples=100 ETKDG confs
          // and keeps ~2 after RMSD dedup.
          {
            const char *rp = std::getenv("RIGID_POOL");
            const char *rk = std::getenv("RIGID_KEEP");
            if (rp || rk) {
              FragmentParams rg =
                  pp.embedding.classParams[FragmentClass::Rigid];
              if (rp) {
                rg.maxSamples = std::atoi(rp);
                rg.minSamples = std::min(rg.minSamples, rg.maxSamples);
              }
              if (rk) rg.maxConfs = std::atoi(rk);
              pp.embedding.setClassParam(FragmentClass::Rigid, rg);
            }
          }
          if (std::getenv("FRAGLIB_ENERGYLOG"))
            pp.embedding.logFragmentEnergies = true;
          requireValidParams(pp, "runBank");
          const double t0 = nowMs();
          auto rr = FragmentConfGen(pp).build(*j.ref);
          ms = nowMs() - t0;
          ens = rr.conformers;
          en.reserve(ens.size());
          for (auto &c : ens) {
            double e = std::numeric_limits<double>::quiet_NaN();
            if (c && c->hasProp(kEnergyProp)) {
              try {
                e = c->getProp<double>(kEnergyProp);
              } catch (...) {
              }
            }
            if (std::isnan(e) && c) e = mmffSinglePoint(*c, mmff);
            en.push_back(e);
          }
        } catch (const std::exception &ex) {
          if (std::getenv("BANK_SHOWERR"))
            std::fprintf(stderr, "[bank] EXCEPTION %s: %s\n", j.name.c_str(),
                         ex.what());
        } catch (...) {
          if (std::getenv("BANK_SHOWERR"))
            std::fprintf(stderr, "[bank] EXCEPTION %s: (unknown)\n",
                         j.name.c_str());
        }
        emit(method, j, ens, en, ms);
      };
      runOurs(false, "ours_tree");
      if (doTS) runOurs(true, "ours_ts");
      if (traceBank) {
        std::fprintf(stderr, "[bank] <<< %zu %s ours done\n", i,
                     j.name.c_str());
        std::fflush(stderr);
      }

      
      const size_t d = ++done;
      if (d % 100 == 0) {
        std::fprintf(stderr, "[bank] %zu/%zu mols\n", d, jobs.size());
        std::fflush(stderr);
      }
    }
  };
  std::vector<std::thread> pool;
  for (unsigned int t = 0; t < threads; ++t) pool.emplace_back(worker);
  for (auto &t : pool) t.join();
  for (auto &kv : writers) kv.second->close();

  // per-molecule stats TSV (analysis / xlsx offline)
  const std::string tsv = outdir + "/stats.tsv";
  {
    std::ofstream f(tsv);
    f << "method\tname\tn_rotors\tn_confs\tgen_ms\tbest_rms\tbest_rank\t"
         "rank_lt1\trank_lt2\thit_lt1\thit_lt2\n";
    for (const auto &r : statRows) f << r << "\n";
  }
  // Speed breakdown (diagnostics.ASM_PROFILE): per-class embed cost, interFrag
    // vs full-MMFF score counts, combine pool pressure, prune/dedup time.
  if (!paramsOverride.empty() || std::getenv("BANK_PROFILE")) {
    printJoinerProfile("bank", jobs.size());
  }
  std::printf("\n[bank] wrote %zu stat rows -> %s\n", statRows.size(),
              tsv.c_str());
  std::printf("[bank] SDF ensembles -> %s/{ours_tree%s}.sdf\n",
              outdir.c_str(), doTS ? ",ours_ts" : "");
  std::fflush(stdout);
}

//! Dump OUR Hamburg (TorLib) per-bond torsion assignment for one molecule, to
//! matched rule / symmetry / angles.
//! Validates our rule matcher + angle expansion against the torsion library's use
//! of the same library. Env: TORDUMP_SMI (required), FCG_HAMBURG_LIB/FALLBACK
//! (the XML).  Sets TORLIB_DUMP.
void runTorsionDump() {
  const char *smiEnv = std::getenv("TORDUMP_SMI");
  if (!smiEnv) {
    std::cerr << "torsiondump: set TORDUMP_SMI=<smiles>\n";
    return;
  }
  // TORDUMP_CMP: print BOTH samplers' RETURNED (driver-frame) angles per bond,
  // so ETKDG's divergence from the TorLib assignment is directly
  // visible.  Default (unset) sets TORLIB_DUMP and prints only TorLib's
  // lib-frame assignment (diff vs confgen -v DEBUG via torcompare.py).
  const bool cmp = std::getenv("TORDUMP_CMP") != nullptr;
  if (!cmp) {
#ifdef _WIN32
    _putenv_s("TORLIB_DUMP", "1");
#else
    setenv("TORLIB_DUMP", "1", 1);
#endif
  }
  RDLog::LogStateSetter blk;
  std::unique_ptr<RWMol> mol(SmilesToMol(smiEnv));
  if (!mol) {
    std::cerr << "torsiondump: bad SMILES\n";
    return;
  }
  MolOps::addHs(*mol);
  DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
  ps.randomSeed = 0xf00d;
  if (DGeomHelpers::EmbedMolecule(*mol, ps) < 0) {
    std::cerr << "torsiondump: embed failed\n";
    return;
  }
  auto torlib =
      hamburgSampler();  // faithful by default (TORLIB_LEAN etc. still honored)
  auto etkdg = std::make_shared<ETKDGTorsionSampler>();
  std::printf("[torsiondump] %s%s\n", smiEnv,
              cmp ? "  (TorLib vs ETKDG, driver-frame angles)"
                  : "  (TorLib; compare to confgen -v DEBUG)");
  auto fmt = [](std::vector<double> a) {
    std::sort(a.begin(), a.end());
    std::string s;
    char buf[16];
    for (double x : a) {
      std::snprintf(buf, sizeof(buf), "%.0f ", x);
      s += buf;
    }
    return s;
  };
  for (const auto b : mol->bonds()) {
    if (b->getBondType() != Bond::SINGLE || b->getIsAromatic()) continue;
    if (mol->getRingInfo()->numBondRings(b->getIdx())) continue;
    const unsigned int j = b->getBeginAtomIdx(), k = b->getEndAtomIdx();
    auto heavyNbr = [&](unsigned int c, unsigned int other) -> int {
      for (const auto nbr : mol->atomNeighbors(mol->getAtomWithIdx(c))) {
        if (nbr->getIdx() != other && nbr->getAtomicNum() > 1)
          return static_cast<int>(nbr->getIdx());
      }
      return -1;
    };
    const int i = heavyNbr(j, k), l = heavyNbr(k, j);
    if (i < 0 || l < 0) continue;  // terminal bond, not a real rotor
    if (mol->getAtomWithIdx(j)->getDegree() < 2 ||
        mol->getAtomWithIdx(k)->getDegree() < 2)
      continue;
    auto tl =
        torlib->getAngles(*mol, static_cast<unsigned int>(i), j, k,
                          static_cast<unsigned int>(l));  // TORLIB_DUMP if !cmp
    if (cmp) {
      auto et = etkdg->getAngles(*mol, static_cast<unsigned int>(i), j, k,
                                 static_cast<unsigned int>(l));
      std::printf("  bond %u-%u  torlib[%zu]: %-42s etkdg[%zu]: %s\n", j, k,
                  tl.size(), fmt(tl).c_str(), et.size(), fmt(et).c_str());
    }
  }
}

//! Dump a full FragmentConfGenParams template (every knob at its default) as
//! editable text, to FRAGCG_PARAMS_OUT if set, else stdout.  Feed the edited
//! file back via FRAGCG_PARAMS.
void runParamsTemplate() {
  RDKit::FragmentConfGenParams defaults;
  if (const char *out = std::getenv("FRAGCG_PARAMS_OUT")) {
    std::ofstream os(out);
    if (!os) {
      std::fprintf(stderr, "FRAGCG_PARAMS_OUT: cannot open '%s'\n", out);
      std::exit(1);
    }
    RDKit::writeFragmentConfGenParams(os, defaults);
    std::fprintf(stderr, "wrote default params template to %s\n", out);
  } else {
    RDKit::writeFragmentConfGenParams(std::cout, defaults);
  }
}

struct Profile {
  const char *name;
  void (*fn)();
  const char *help;
};

const std::vector<Profile> &profiles() {
  static const std::vector<Profile> ps = {
      {"platinum", runPlatinum,
       "bioactive-conformer reproduction on the Platinum Diverse Dataset vs ETKDG "
       ""},
      {"worst", runWorst,
       "low-rotor ligands where ours reproduces poorly (worst-first, names) "
       "(PLATINUM_MAXROT)"},
      {"casestudy", runCaseStudy,
       "one ligand (PLATINUM_NAME): per-fragment exit-vector fidelity + shape "
       "across ETKDG/DG/+MMFF; SDF dump (PLATINUM_SDFOUT)"},
      {"asmsampler", runAsmSampler,
       "sweep junction-angle samplers (ETKDG/TorLib/Uniform +/- backstop) through "
       "the CoordsOnly joiner (PLATINUM_MAXMOLS default 50)"},
      {"rotorcross", runRotorCross,
       "ETKDG & TorLib x tree-search & Thompson: assembly-only time bucketed by "
       "rotor count to find the TS-vs-tree speed crossover (PLATINUM_MAXMOLS, "
       "FRAGCG_PARAMS override file)"},
      {"params", runParamsTemplate,
       "dump a full default FragmentConfGenParams template (editable key=value text) "
       "to FRAGCG_PARAMS_OUT else stdout; feed back via FRAGCG_PARAMS"},
      {"torsiondump", runTorsionDump,
       "dump OUR Hamburg per-bond torsion assignment (rule/symmetry/angles) for "
       "TORDUMP_SMI"},
      {"bank", runBank,
       "OURS(tree) over Platinum: per-method SDF "
       "ensembles (MMFF energy, #rotors, gen ms, crystal-pose placement <1/<2A) + "
       "stats.tsv (BANK_OUTDIR, BANK_TS adds Thompson, PLATINUM_MAXMOLS)"},
      {"rmsdbench", runRmsdBench,
       "heavy-atom RMSD timing: raw-index vs QCP vs symmetry-QCP vs RDKit getBestRMS"},
      {"matrix", runMatrix,
       "full factorial: 3 embeddings x 3 torsion profiles x 3 search cells at "
       "MMFF94s; %<1/%<2/genN + cold&warm ms/mol (PLATINUM_MAXMOLS to subset)"},
      {"pocctx", runPocCtx,
       "PoC: aromatic-context vs isolated carbon-cap fragment embed, crystal core "
       "RMSD (PLATINUM_MAXMOLS default 40, POC_MMFF)"},
      {"fragrank", runFragRank,
       "SAMPLING vs RANKING for flexible-ring fragments: is the crystal pucker in "
       "our pool, and where does its energy rank (FRAGRANK_POOL)"},
      {"termtable", runTermTable,
       "per-conformer FULL MMFF94s (no-estat) term breakdown (bond/angle/sb/oop/"
       "tors/vdw + interVdw/interSB) for OURS, sorted by total, crystal-"
       "closest marked (PLATINUM_MAXMOLS, TERMTABLE_MOLS/OURCONFS, "
       "PLATINUM_NAME)"},
      {"fragrms", runFragRms,
       "best superposed (QCP) RMSD from each CRYSTAL fragment to our embedded pool -- search-free fragment-geometry check (PLATINUM_MINROT/MAXROT/MAXMOLS, FRAGRMS_CONFS)"},
      {"torcheck", runTorCheck,
       "how far is each CRYSTAL junction dihedral from the nearest angle our samplers "
       "propose (ETKDG vs TorLib, base vs basin) -- search-free test of angle CHOICE "
       "(PLATINUM_MINROT/MAXROT/MAXMOLS)"},
      {"xtalrecon", runXtalRecon,
       "EXACT crystal fragments through our assembly pathway: reproduction (best-of-N "
       "RMSD) + where the true crystal pose RANKS in our scorer (PLATINUM_MAXMOLS, "
       "XTAL_MMFF94, XTAL_RELAX)"},
  };
  return ps;
}

void usage(const char *argv0) {
  std::cerr << "usage: " << argv0 << " <profile>\n\nprofiles:\n";
  for (const auto &p : profiles()) {
    std::fprintf(stderr, "  %-11s %s\n", p.name, p.help);
  }
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    usage(argv[0]);
    return 2;
  }
  const std::string want = argv[1];
  for (const auto &p : profiles()) {
    if (want == p.name) {
      p.fn();
      return 0;
    }
  }
  std::cerr << "unknown profile: " << want << "\n\n";
  usage(argv[0]);
  return 2;
}
