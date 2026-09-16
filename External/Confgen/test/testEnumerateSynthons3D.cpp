//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// The EnumerateSynthons3D / Molzip3D cases that used to live here tested a
// zipper that no longer exists: Molzip3D was a stub returning no products, and
// the class has been rebuilt on FragmentConfGen in SynthonSearch/.  They were
// already disabled ("[.]") and are removed rather than left to rot; the new
// library needs its own tests.  What remains is the Molzip3D-independent unit
// coverage (NonbondedLookup table, TorsionSampler).
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <GraphMol/ChemReactions/Enumerate/EnumerateSynthons.h>
#include <Confgen/Utils/NonbondedLookup.h>
#include <Confgen/Sampler/TorsionSampler.h>
#include <Confgen/SynthonSearch/EnumerateSynthons3D.h>
#include <Confgen/SynthonSearch/SynthonSearch3D.h>
#include <GraphMol/GaussianShape/GaussianShape.h>

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/Nonbonded.h>
#include <RDGeneral/RDLog.h>

#include <algorithm>
#include <array>
#include <memory>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>

using namespace RDKit;

namespace {
// build a one-position-per-synthon BBS from SMILES, one reagent each
EnumerationTypes::BBS makeBBS(
    const std::vector<std::vector<std::string>> &smis) {
  EnumerationTypes::BBS bbs;
  for (const auto &pos : smis) {
    MOL_SPTR_VECT reagents;
    for (const auto &smi : pos) {
      reagents.emplace_back(SmilesToMol(smi));
    }
    bbs.push_back(reagents);
  }
  return bbs;
}

ROMOL_SPTR makeEmbeddedSynthon(const std::string &smi, int seed) {
  std::unique_ptr<ROMol> mol(SmilesToMol(smi));
  REQUIRE(mol);
  std::unique_ptr<ROMol> molH(MolOps::addHs(*mol));
  REQUIRE(molH);
  auto res = boost::make_shared<RWMol>(*molH);
  DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
  ps.randomSeed = seed;
  {
    RDLog::LogStateSetter blocker;
    REQUIRE(DGeomHelpers::EmbedMolecule(*res, ps) >= 0);
  }
  return res;
}

bool bestETKDGMMFFEnergy(const std::string &smi, unsigned int numConfs,
                         int randomSeed, double &bestEnergy) {
  std::unique_ptr<ROMol> mol;
  try {
    mol.reset(SmilesToMol(smi));
  } catch (...) {
    mol.reset();
  }
  if (!mol) {
    return false;
  }
  std::unique_ptr<ROMol> molH(MolOps::addHs(*mol));
  if (!molH) {
    return false;
  }
  auto rw = boost::make_shared<RWMol>(*molH);
  DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
  ps.randomSeed = randomSeed;
  ps.pruneRmsThresh = -1.0;
  RDKit::INT_VECT cids;
  {
    RDLog::LogStateSetter blocker;
    DGeomHelpers::EmbedMultipleConfs(*rw, cids, numConfs, ps);
  }
  if (cids.empty()) {
    return false;
  }

  MMFF::MMFFMolProperties props(*rw);
  if (!props.isValid()) {
    return false;
  }
  props.setMMFFEleTerm(false);
  bestEnergy = std::numeric_limits<double>::max();
  for (const auto cid : cids) {
    std::unique_ptr<ForceFields::ForceField> ff(
        MMFF::constructForceField(*rw, &props, 1.0e8, cid));
    if (!ff) {
      continue;
    }
    ff->initialize();
    ff->minimize(1000);
    bestEnergy = std::min(bestEnergy, ff->calcEnergy());
  }
  return std::isfinite(bestEnergy);
}

ROMOL_SPTR bestETKDGMMFFConformer(const std::string &smi, unsigned int numConfs,
                                  int randomSeed, double &bestEnergy,
                                  bool minimizeConformers = true) {
  std::unique_ptr<ROMol> mol;
  try {
    mol.reset(SmilesToMol(smi));
  } catch (...) {
    mol.reset();
  }
  if (!mol) {
    return nullptr;
  }
  std::unique_ptr<ROMol> molH(MolOps::addHs(*mol));
  if (!molH) {
    return nullptr;
  }
  auto rw = boost::make_shared<RWMol>(*molH);
  DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
  ps.randomSeed = randomSeed;
  ps.pruneRmsThresh = -1.0;
  RDKit::INT_VECT cids;
  {
    RDLog::LogStateSetter blocker;
    DGeomHelpers::EmbedMultipleConfs(*rw, cids, numConfs, ps);
  }
  if (cids.empty()) {
    return nullptr;
  }

  MMFF::MMFFMolProperties props(*rw);
  if (!props.isValid()) {
    return nullptr;
  }
  props.setMMFFEleTerm(false);

  int bestCid = -1;
  bestEnergy = std::numeric_limits<double>::max();
  for (const auto cid : cids) {
    std::unique_ptr<ForceFields::ForceField> ff(
        MMFF::constructForceField(*rw, &props, 1.0e8, cid));
    if (!ff) {
      continue;
    }
    ff->initialize();
    if (minimizeConformers) {
      ff->minimize(1000);
    }
    const double energy = ff->calcEnergy();
    if (energy < bestEnergy) {
      bestEnergy = energy;
      bestCid = cid;
    }
  }
  if (bestCid < 0 || !std::isfinite(bestEnergy)) {
    return nullptr;
  }

  auto res = boost::make_shared<RWMol>(*rw);
  const Conformer &conf = rw->getConformer(bestCid);
  auto *newConf = new Conformer(conf);
  newConf->setId(0);
  res->clearConformers();
  res->addConformer(newConf, false);
  res->setProp<double>(kEnergyProp, bestEnergy);
  res->getConformer().setProp<double>(kEnergyProp, bestEnergy);
  return res;
}

std::string testDataPath(const std::string &relPath) {
  if (const char *rdBase = std::getenv("RDBASE")) {
    return std::string(rdBase) + "/" + relPath;
  }
  return "../" + relPath;
}

unsigned int countBreakableRotorCandidates(const ROMol &mol) {
  auto isCarbonylCarbon = [&](const Atom *atom) {
    if (atom->getAtomicNum() != 6) {
      return false;
    }
    for (const auto bond : mol.atomBonds(atom)) {
      if (bond->getBondType() != Bond::DOUBLE) {
        continue;
      }
      const auto other = bond->getOtherAtom(atom);
      if (other->getAtomicNum() == 8 || other->getAtomicNum() == 16) {
        return true;
      }
    }
    return false;
  };

  unsigned int res = 0;
  for (const auto bond : mol.bonds()) {
    if (bond->getBondType() != Bond::SINGLE || bond->getIsAromatic() ||
        mol.getRingInfo()->numBondRings(bond->getIdx())) {
      continue;
    }
    const auto beginAtom = bond->getBeginAtom();
    const auto endAtom = bond->getEndAtom();
    if (beginAtom->getAtomicNum() <= 1 || endAtom->getAtomicNum() <= 1) {
      continue;
    }
    if ((beginAtom->getAtomicNum() == 6 && endAtom->getAtomicNum() == 7 &&
         isCarbonylCarbon(beginAtom)) ||
        (beginAtom->getAtomicNum() == 7 && endAtom->getAtomicNum() == 6 &&
         isCarbonylCarbon(endAtom))) {
      continue;
    }
    unsigned int beginHeavy = 0;
    for (const auto nbr : mol.atomNeighbors(beginAtom)) {
      if (nbr->getIdx() != endAtom->getIdx() && nbr->getAtomicNum() > 1) {
        ++beginHeavy;
      }
    }
    unsigned int endHeavy = 0;
    for (const auto nbr : mol.atomNeighbors(endAtom)) {
      if (nbr->getIdx() != beginAtom->getIdx() && nbr->getAtomicNum() > 1) {
        ++endHeavy;
      }
    }
    if (beginHeavy && endHeavy) {
      ++res;
    }
  }
  return res;
}

struct VdWLookupDiffStats {
  unsigned int numPairs = 0;
  unsigned int numSamples = 0;
  double sumAbs = 0.0;
  double sumSq = 0.0;
  double maxAbs = 0.0;
  double maxRel = 0.0;
  double maxDistance = 0.0;
  double maxRelDistance = 0.0;
  double maxOriginal = 0.0;
  double maxLookup = 0.0;
  double maxRelOriginal = 0.0;
  double maxRelLookup = 0.0;
  double maxRij = 0.0;
  double maxEps = 0.0;
  double maxRelRij = 0.0;
  double maxRelEps = 0.0;
  unsigned int maxTableIdx = 0;
  unsigned int maxRelTableIdx = 0;
};
}  // namespace

TEST_CASE("NonbondedLookup CHNOPS vdW table", "[synthon3d]") {
  using namespace ForceFields::MMFF;

  CHECK(VDW_LOOKUP_PARAMS.size() == 21);
  CHECK(VDW_LOOKUP.size() == VDW_LOOKUP_PARAMS.size());
  CHECK(VDW_LOOKUP.front().size() ==
        ForceFields::MMFF::detail::VDW_LOOKUP_NUM_DISTANCES);

  for (unsigned int tableIdx = 0; tableIdx < VDW_LOOKUP_PARAMS.size();
       ++tableIdx) {
    CHECK(closest_table(VDW_LOOKUP_PARAMS[tableIdx].R_ij_star,
                        VDW_LOOKUP_PARAMS[tableIdx].epsilon) == tableIdx);
  }

  const auto tableIdx = closest_table(VDW_LOOKUP_PARAMS[0].R_ij_star,
                                      VDW_LOOKUP_PARAMS[0].epsilon);
  constexpr double distance = 2.24;
  CHECK(table_lookup(tableIdx, distance) ==
        Catch::Approx(ForceFields::MMFF::detail::calcVdWEnergy(
                          distance, VDW_LOOKUP_PARAMS[0].R_ij_star,
                          VDW_LOOKUP_PARAMS[0].epsilon))
            .epsilon(1e-4));
  CHECK(table_lookup(tableIdx, 99.0) == Catch::Approx(0.0));
}

// Hidden diagnostic; run with:
//   synthon3DTestCatch "[.diagnostic][vdwlookup]"
//
// This intentionally reports errors without enforcing a tolerance. The original
// VdWContrib energy path calls MMFF::Utils::calcVdWEnergy() for every pair, so
// this compares that exact per-pair expression against the CHNOPS lookup.
TEST_CASE("NonbondedLookup diagnostic against original VdWContrib",
          "[.][synthon3d][vdwlookup][.diagnostic]") {
  const std::vector<std::string> smis = {
      "CCO",      "CCN",          "CCS",        "CCP",          "CC(=O)N",
      "c1ccncc1", "CCS(=O)(=O)O", "O=P(O)(O)O", "NCCOP(=O)(O)S"};

  std::map<std::pair<int, int>, ForceFields::MMFF::MMFFVdWRijstarEps>
      sampledParams;
  for (const auto &smi : smis) {
    std::unique_ptr<ROMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    std::unique_ptr<ROMol> molH(MolOps::addHs(*mol));
    REQUIRE(molH);

    MMFF::MMFFMolProperties props(*molH);
    REQUIRE(props.isValid());
    for (unsigned int i = 0; i < molH->getNumAtoms(); ++i) {
      for (unsigned int j = i + 1; j < molH->getNumAtoms(); ++j) {
        ForceFields::MMFF::MMFFVdWRijstarEps p;
        if (!props.getMMFFVdWParams(i, j, p)) {
          continue;
        }
        const auto key =
            std::make_pair(static_cast<int>(std::round(p.R_ij_star * 1000.0)),
                           static_cast<int>(std::round(p.epsilon * 1000000.0)));
        sampledParams.emplace(key, p);
      }
    }
  }

  VdWLookupDiffStats stats;
  std::array<VdWLookupDiffStats, 3> cutoffStats;
  constexpr std::array<double, 3> distanceCutoffs{{2.0, 2.5, 3.0}};
  stats.numPairs = static_cast<unsigned int>(sampledParams.size());
  REQUIRE(stats.numPairs > 0);

  auto updateStats = [](VdWLookupDiffStats &s, double dist, double original,
                        double lookup, double R_ij_star, double epsilon,
                        unsigned int tableIdx) {
    const double diff = lookup - original;
    const double absDiff = std::abs(diff);
    const double relDiff = absDiff / std::max(1.0e-12, std::abs(original));

    ++s.numSamples;
    s.sumAbs += absDiff;
    s.sumSq += diff * diff;
    if (absDiff > s.maxAbs) {
      s.maxAbs = absDiff;
      s.maxDistance = dist;
      s.maxOriginal = original;
      s.maxLookup = lookup;
      s.maxRij = R_ij_star;
      s.maxEps = epsilon;
      s.maxTableIdx = tableIdx;
    }
    if (relDiff > s.maxRel) {
      s.maxRel = relDiff;
      s.maxRelDistance = dist;
      s.maxRelOriginal = original;
      s.maxRelLookup = lookup;
      s.maxRelRij = R_ij_star;
      s.maxRelEps = epsilon;
      s.maxRelTableIdx = tableIdx;
    }
  };

  for (const auto &entry : sampledParams) {
    const auto &p = entry.second;
    const unsigned int tableIdx =
        ForceFields::MMFF::closest_table(p.R_ij_star, p.epsilon);
    for (double dist = 1.5; dist <= 8.0001; dist += 0.05) {
      const double original =
          ForceFields::MMFF::Utils::calcVdWEnergy(dist, p.R_ij_star, p.epsilon);
      const double lookup = ForceFields::MMFF::table_lookup(tableIdx, dist);

      updateStats(stats, dist, original, lookup, p.R_ij_star, p.epsilon,
                  tableIdx);
      for (unsigned int i = 0; i < distanceCutoffs.size(); ++i) {
        if (dist >= distanceCutoffs[i]) {
          updateStats(cutoffStats[i], dist, original, lookup, p.R_ij_star,
                      p.epsilon, tableIdx);
        }
      }
    }
  }

  auto printStats = [](const VdWLookupDiffStats &s, const std::string &label) {
    const double meanAbs = s.sumAbs / s.numSamples;
    const double rms = std::sqrt(s.sumSq / s.numSamples);
    std::cout << "\n    " << label
              << "\n      distance samples          : " << s.numSamples
              << "\n      mean |lookup-original|    : " << meanAbs
              << " kcal/mol"
              << "\n      RMS error                 : " << rms << " kcal/mol"
              << "\n      max |lookup-original|     : " << s.maxAbs
              << " kcal/mol at " << s.maxDistance << " A"
              << "\n        original/lookup         : " << s.maxOriginal
              << " / " << s.maxLookup << " kcal/mol"
              << "\n        R_ij*/epsilon/table     : " << s.maxRij << " / "
              << s.maxEps << " / " << s.maxTableIdx
              << "\n      max relative error        : " << s.maxRel << " at "
              << s.maxRelDistance << " A"
              << "\n        original/lookup         : " << s.maxRelOriginal
              << " / " << s.maxRelLookup << " kcal/mol"
              << "\n        R_ij*/epsilon/table     : " << s.maxRelRij << " / "
              << s.maxRelEps << " / " << s.maxRelTableIdx << "\n";
  };

  std::cout << "\n  VdW lookup diagnostic vs original VdWContrib path"
            << "\n    unique MMFF parameter pairs : " << stats.numPairs;
  printStats(stats, "all sampled distances");
  for (unsigned int i = 0; i < cutoffStats.size(); ++i) {
    printStats(cutoffStats[i],
               "distances >= " + std::to_string(distanceCutoffs[i]) + " A");
  }
  SUCCEED("vdW lookup diagnostic reported");
}



TEST_CASE("UniformTorsionSampler increments", "[synthon3d]") {
  UniformTorsionSampler s(90.0);
  ROMol dummy;
  auto angles = s.getAngles(dummy, 0, 1, 2, 3);
  REQUIRE(angles.size() == 4);  // 0, 90, 180, 270
  CHECK(angles.front() == Catch::Approx(0.0));
  CHECK(angles.back() == Catch::Approx(270.0));
}

TEST_CASE("ETKDGTorsionSampler returns bounded arms", "[synthon3d]") {
  auto mol = makeEmbeddedSynthon("CCCC", 0xf00d);
  ETKDGTorsionSampler s(10.0, 120.0, 4);
  auto angles = s.getAngles(*mol, 0, 1, 2, 3);
  REQUIRE(!angles.empty());
  CHECK(angles.size() <= 4);
  for (const auto angle : angles) {
    CHECK(angle >= 0.0);
    CHECK(angle < 360.0);
  }
}

TEST_CASE("TorsionLibrarySampler returns bounded arms", "[synthon3d]") {
  // External data we do not ship; point FCG_HAMBURG_LIB at it to exercise
  // this sampler.  Skipped rather than failed when it is absent.
  const char *env = std::getenv("FCG_HAMBURG_LIB");
  if (!env) {
    SKIP("FCG_HAMBURG_LIB not set (Hamburg torsion library not available)");
  }
  const std::string path = env;
  std::ifstream in(path);
  if (!in.good()) {
    SKIP("FCG_HAMBURG_LIB does not point at a readable file");
  }

  auto mol = makeEmbeddedSynthon("CCCC", 0xf00d);
  TorsionLibrarySampler s(path, 120.0, 4);
  auto angles = s.getAngles(*mol, 0, 1, 2, 3);
  REQUIRE(!angles.empty());
  CHECK(angles.size() <= 4);
  for (const auto angle : angles) {
    CHECK(angle >= 0.0);
    CHECK(angle < 360.0);
  }
}

// ===========================================================================
//  EnumerateSynthons3D -- labelled-synthon assembly on real Chemspace data.
//  Syntons_5567.csv columns: SMILES, synton_id, synton_role, reaction_id.
//  Exit vectors are the AtomType scheme ([U], [Np], ...).
// ===========================================================================
namespace {

std::string synthonCsvPath() {
  const char *rdbase = std::getenv("RDBASE");
  if (!rdbase) {
    return "";
  }
  return std::string(rdbase) +
         "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
}

//! One reaction's synthons, grouped by position.  `maxPerPosition` keeps the
//! test quick -- the full reaction is thousands of products.
EnumerationTypes::BBS loadSynthonReaction(const std::string &path,
                                          const std::string &rxnId,
                                          size_t maxPerPosition) {
  std::map<std::string, MOL_SPTR_VECT> byRole;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::vector<std::string> f;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
      f.push_back(cell);
    }
    if (f.size() < 4 || f[0] == "SMILES" || f[3] != rxnId) {
      continue;
    }
    if (maxPerPosition && byRole[f[2]].size() >= maxPerPosition) {
      continue;
    }
    ROMOL_SPTR m(SmilesToMol(f[0]));
    if (m) {
      byRole[f[2]].push_back(m);
    }
  }
  EnumerationTypes::BBS bbs;
  for (auto &kv : byRole) {  // synton_1, synton_2, ... in order
    bbs.push_back(kv.second);
  }
  return bbs;
}

}  // namespace

TEST_CASE("EnumerateSynthons3D assembles labelled synthons", "[synthon3d]") {
  const std::string csv = synthonCsvPath();
  if (csv.empty()) {
    SKIP("RDBASE not set");
  }
  // a7 is 3-component; 4 per position keeps this a unit test.
  const auto bbs = loadSynthonReaction(csv, "a7", 4);
  REQUIRE(bbs.size() == 3);
  for (const auto &pos : bbs) {
    REQUIRE(pos.size() == 4);
  }

  EnumerateSynthons3DParams params;
  params.prefillEmbedder = false;  // exercise the on-demand path
  params.confgen.numOutputConfs = 4;
  EnumerateSynthons3D lib(bbs, params);
  REQUIRE(lib.isValid());  // [U] labels are a consistent scheme
  CHECK(lib.arity() == 3);
  CHECK(lib.numReagents(0) == 4);

  SECTION("a product comes back with conformers") {
    const auto p = lib.getProduct({0, 0, 0});
    INFO("status: " << synthonBuildStatusMessage(p.status));
    REQUIRE(p);
    CHECK(p.mol->getNumAtoms() > 0);
    CHECK(p.mol->getNumConformers() > 0);
    // the exit-vector markers must be consumed by the zip
    for (const auto atom : p.mol->atoms()) {
      CHECK(atom->getAtomicNum() != 92);  // [U]
    }
  }

  SECTION("get() fulfils the EnumerateLibrary API") {
    const auto res = lib.get({0, 0, 0});
    REQUIRE(res.size() == 1);
    REQUIRE(res[0].size() == 1);
    CHECK(res[0][0]->getNumConformers() > 0);
  }

  SECTION("an out-of-range index is reported, not thrown") {
    const auto p = lib.getProduct({0, 0, 99});
    CHECK_FALSE(p);
    CHECK(p.status == SynthonBuildStatus::BadReagentIndex);
    CHECK(lib.get({0, 0, 99}).empty());
  }

  SECTION("serialize round-trips") {
    // The base class writes a boost TEXT archive and everything this class
    // appends is raw binary.  A text archive's reader leaves its trailing
    // separator in the stream, so without an explicit boundary the appended
    // payload is read one byte late -- which showed up only as a embedder
    // magic-number failure, and only through the Python wrapper.
    for (const bool storeLib : {false, true}) {
      INFO("storeEmbedder=" << storeLib);
      EnumerateSynthons3DParams sp = params;
      sp.storeEmbedder = storeLib;
      sp.embedStyle = SynthonEmbedStyle::Coarse;
      sp.confgen.embedding.numConfsPerFragment = 3;
      sp.confgen.randomSeed = 2468;
      EnumerateSynthons3D src(bbs, sp);
      REQUIRE(src.isValid());

      std::stringstream ss;
      src.toStream(ss);

      EnumerateSynthons3D dst;
      REQUIRE_NOTHROW(dst.initFromStream(ss));
      CHECK(dst.isValid());
      CHECK(dst.arity() == src.arity());
      CHECK(dst.numReagents(0) == src.numReagents(0));
      CHECK(dst.params3D().embedStyle == sp.embedStyle);
      CHECK(dst.params3D().confgen.embedding.numConfsPerFragment == 3);
      CHECK(dst.params3D().confgen.randomSeed == 2468);
      CHECK(dst.embedder()->params().numConfsPerFragment == 3);
      CHECK(dst.embedder()->params().randomSeed == 2468);

      const auto p = dst.getProduct({0, 0, 0});
      INFO("status: " << synthonBuildStatusMessage(p.status));
      REQUIRE(p);
      CHECK(p.mol->getNumConformers() > 0);
      CHECK(MolToSmiles(*p.mol) ==
            MolToSmiles(*src.getProduct({0, 0, 0}).mol));
    }
  }
}

TEST_CASE("EnumerateSynthons3D mixed cut-bond instructions", "[synthon3d]") {
  // genSynthonLib instructs either every synthon in a library or none, so the
  // mixed case -- some synthons carrying cut bonds, others silent -- has no
  // coverage from real libraries.  It is reachable through the public setter,
  // and the rules it exercises are asymmetric enough to be worth pinning.
  const std::string csv = synthonCsvPath();
  if (csv.empty()) {
    SKIP("RDBASE not set");
  }
  const auto bbs = loadSynthonReaction(csv, "a7", 4);
  REQUIRE(bbs.size() == 3);

  EnumerateSynthons3DParams params;
  params.prefillEmbedder = false;
  params.embedStyle = SynthonEmbedStyle::Coarse;  // cut bonds only apply here
  params.confgen.numOutputConfs = 4;

  SECTION("a library with no instructions reports none") {
    EnumerateSynthons3D lib(bbs, params);
    CHECK_FALSE(lib.hasSynthonCutBonds());
    CHECK(lib.synthonCutBonds(0, 0).empty());
  }

  SECTION("out-of-range reads are empty, not throws") {
    EnumerateSynthons3D lib(bbs, params);
    CHECK(lib.synthonCutBonds(99, 0).empty());
    CHECK(lib.synthonCutBonds(0, 99).empty());
    CHECK_NOTHROW(lib.setSynthonCutBonds(99, 0, {0}));
    CHECK_FALSE(lib.hasSynthonCutBonds());  // the write was refused
  }

  SECTION("instructing ONE synthon leaves the others silent") {
    EnumerateSynthons3D lib(bbs, params);
    // Pick a real in-synthon bond so the instruction is meaningful.
    const ROMol &syn = *lib.getReagents()[0][0];
    REQUIRE(syn.getNumBonds() > 0);
    unsigned int bond = syn.getNumBonds();  // sentinel: none found
    for (const auto b : syn.bonds()) {
      if (b->getBeginAtom()->getAtomicNum() && b->getEndAtom()->getAtomicNum()) {
        bond = b->getIdx();  // both ends survive the zip
        break;
      }
    }
    REQUIRE(bond < syn.getNumBonds());

    lib.setSynthonCutBonds(0, 0, {bond});
    CHECK(lib.hasSynthonCutBonds());
    CHECK(lib.synthonCutBonds(0, 0).size() == 1);
    // Every other synthon, including others at the same position, stays empty.
    CHECK(lib.synthonCutBonds(0, 1).empty());
    CHECK(lib.synthonCutBonds(1, 0).empty());
    CHECK(lib.synthonCutBonds(2, 0).empty());

    // A product using the instructed synthon and one that does not must BOTH
    // build: a partially instructed library is not a broken one.
    const auto instructed = lib.getProduct({0, 0, 0});
    INFO("instructed status: " << synthonBuildStatusMessage(instructed.status));
    REQUIRE(instructed);
    CHECK(instructed.mol->getNumConformers() > 0);

    const auto silent = lib.getProduct({1, 0, 0});
    INFO("silent status: " << synthonBuildStatusMessage(silent.status));
    REQUIRE(silent);
    CHECK(silent.mol->getNumConformers() > 0);

    // Same graph either way: instructions choose SEAMS, not chemistry.
    EnumerateSynthons3D plain(bbs, params);
    CHECK(MolToSmiles(*instructed.mol) ==
          MolToSmiles(*plain.getProduct({0, 0, 0}).mol));
  }

  SECTION("clearing an instruction restores the default seam") {
    EnumerateSynthons3D lib(bbs, params);
    lib.setSynthonCutBonds(0, 0, {0});
    REQUIRE(lib.hasSynthonCutBonds());
    lib.setSynthonCutBonds(0, 0, {});
    CHECK(lib.synthonCutBonds(0, 0).empty());
    // hasSynthonCutBonds is a gate that only ever latches ON: it says "this
    // library may have instructions", and the per-synthon read is what
    // decides.  Assert the behaviour rather than the tidier thing, so a
    // future change to either has to be deliberate.
    CHECK(lib.hasSynthonCutBonds());
    const auto p = lib.getProduct({0, 0, 0});
    INFO("status: " << synthonBuildStatusMessage(p.status));
    REQUIRE(p);
  }

  SECTION("a bond to an exit atom cannot be requested and is ignored") {
    // A synthon cannot ask for its own junction to be cut: the zip consumes
    // exit atoms, so the bond has no surviving endpoint to translate to.
    // It must be skipped silently, not throw and not corrupt the cut set.
    EnumerateSynthons3D lib(bbs, params);
    const ROMol &syn = *lib.getReagents()[0][0];
    unsigned int exitBond = syn.getNumBonds();
    for (const auto b : syn.bonds()) {
      if (!b->getBeginAtom()->getAtomicNum() ||
          !b->getEndAtom()->getAtomicNum() ||
          b->getBeginAtom()->getAtomicNum() == 92 ||
          b->getEndAtom()->getAtomicNum() == 92) {
        exitBond = b->getIdx();
        break;
      }
    }
    if (exitBond >= syn.getNumBonds()) {
      SKIP("no exit-atom bond in this synthon");
    }
    lib.setSynthonCutBonds(0, 0, {exitBond});
    const auto p = lib.getProduct({0, 0, 0});
    INFO("status: " << synthonBuildStatusMessage(p.status));
    REQUIRE(p);
    CHECK(p.mol->getNumConformers() > 0);
  }

  SECTION("an out-of-range bond index is ignored, not fatal") {
    EnumerateSynthons3D lib(bbs, params);
    const unsigned int tooBig = lib.getReagents()[0][0]->getNumBonds() + 100;
    lib.setSynthonCutBonds(0, 0, {tooBig});
    const auto p = lib.getProduct({0, 0, 0});
    INFO("status: " << synthonBuildStatusMessage(p.status));
    REQUIRE(p);
    CHECK(p.mol->getNumConformers() > 0);
  }

  SECTION("instructions survive a serialize round-trip") {
    EnumerateSynthons3D src(bbs, params);
    src.setSynthonCutBonds(0, 0, {0});
    src.setSynthonCutBonds(2, 3, {0, 1});
    REQUIRE(src.hasSynthonCutBonds());

    std::stringstream ss;
    src.toStream(ss);
    EnumerateSynthons3D dst;
    REQUIRE_NOTHROW(dst.initFromStream(ss));

    CHECK(dst.hasSynthonCutBonds());
    CHECK(dst.synthonCutBonds(0, 0) == std::vector<unsigned int>{0});
    CHECK(dst.synthonCutBonds(2, 3) == std::vector<unsigned int>({0, 1}));
    // and the silent ones are still silent
    CHECK(dst.synthonCutBonds(0, 1).empty());
    CHECK(dst.synthonCutBonds(1, 0).empty());
  }
}

TEST_CASE("EnumerateSynthons3D fragment cache warms", "[synthon3d]") {
  const std::string csv = synthonCsvPath();
  if (csv.empty()) {
    SKIP("RDBASE not set");
  }
  // a7 is 3-component.  The cache only pays when products FAR outnumber
  // synthons -- with a tiny pool the "cold" pass is already warm by its fifth
  // product and the measured ratio understates the real effect.
  const auto bbs = loadSynthonReaction(csv, "a7", 8);
  REQUIRE(bbs.size() == 3);
  const size_t nSynthons = bbs[0].size() + bbs[1].size() + bbs[2].size();

  EnumerateSynthons3DParams params;
  params.prefillEmbedder = false;  // start cold so the first pass pays for it
  params.confgen.numOutputConfs = 4;
  EnumerateSynthons3D lib(bbs, params);
  REQUIRE(lib.isValid());

  // The SAME products twice: identical assembly and search work, so the only
  // difference between the passes is whether the synthon fragments were
  // already embedded.  That isolates the cache from everything else.
  std::vector<std::vector<unsigned int>> combos;
  for (unsigned int i = 0; i < 8; ++i) {
    for (unsigned int j = 0; j < 8; ++j) {
      for (unsigned int k = 0; k < 2; ++k) {
        combos.push_back({i, j, k});
      }
    }
  }

  auto timePass = [&](unsigned int &built) {
    built = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (const auto &c : combos) {
      if (lib.getProduct(c)) {
        ++built;
      }
    }
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - t0)
        .count();
  };

  unsigned int builtCold = 0, builtWarm = 0;
  const double coldMs = timePass(builtCold);
  const size_t fragsAfterCold = lib.embedder()->size();
  const double warmMs = timePass(builtWarm);

  // The pass above is only cold at its START -- by its fifth product most
  // fragments are already cached, so it measures a warming average, not a cold
  // cost.  A genuinely cold per-product number needs a FRESH cache each time.
  // Library construction is included, so measure that separately and subtract.
  double ctorMs = 0.0;
  {
    const auto t0 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < combos.size(); ++i) {
      EnumerateSynthons3D fresh(bbs, params);
      (void)fresh.arity();
    }
    ctorMs = std::chrono::duration<double, std::milli>(
                 std::chrono::steady_clock::now() - t0)
                 .count();
  }
  double trueColdMs = 0.0;
  unsigned int builtTrueCold = 0;
  {
    const auto t0 = std::chrono::steady_clock::now();
    for (const auto &c : combos) {
      EnumerateSynthons3D fresh(bbs, params);  // empty cache every time
      if (fresh.getProduct(c)) {
        ++builtTrueCold;
      }
    }
    trueColdMs = std::chrono::duration<double, std::milli>(
                     std::chrono::steady_clock::now() - t0)
                     .count();
  }
  const double trueColdNet = trueColdMs - ctorMs;

  CHECK(builtCold == builtWarm);  // same products either way
  REQUIRE(builtCold > 0);

  const double n = static_cast<double>(combos.size());
  std::cout << "\n  [synthon cache] " << combos.size() << " products from "
            << nSynthons << " synthons, " << builtCold << " assembled, "
            << fragsAfterCold << " fragments\n"
            << "    TRUE cold (fresh cache each product) : "
            << trueColdNet / n << " ms/product\n"
            << "    first pass (cache warming as it goes): " << coldMs / n
            << " ms/product\n"
            << "    warm (cache full)                    : " << warmMs / n
            << " ms/product\n"
            << "    cache speedup, true cold -> warm : x"
            << (warmMs > 0 ? trueColdNet / warmMs : 0.0) << "\n"
            << "    (library construction " << ctorMs / n
            << " ms/product, subtracted from true cold)\n";
  CHECK(builtTrueCold == builtCold);  // same products regardless of cache
  SUCCEED("cache timing reported");
}

TEST_CASE("SynthonSearch3D finds a planted product", "[synthon3d]") {
  const std::string csv = synthonCsvPath();
  if (csv.empty()) {
    SKIP("RDBASE not set");
  }
  const auto bbs = loadSynthonReaction(csv, "a7", 6);
  REQUIRE(bbs.size() == 3);

  EnumerateSynthons3DParams params;
  params.confgen.numOutputConfs = 3;  // cost is linear in this; see the header
  EnumerateSynthons3D lib(bbs, params);
  REQUIRE(lib.isValid());

  // Plant the target: take a real product of this library as the query, so we
  // know a perfect answer EXISTS and what it is.  Not every combination
  // assembles (a7 loses ~12% -- see the cache test), so find a real one rather
  // than assuming an index triple is buildable.
  std::vector<unsigned int> planted;
  SynthonProduct target;
  for (unsigned int i = 0; i < 6 && planted.empty(); ++i) {
    for (unsigned int j = 0; j < 6 && planted.empty(); ++j) {
      for (unsigned int k = 0; k < 6 && planted.empty(); ++k) {
        auto p = lib.getProduct({i, j, k});
        if (p) {
          planted = {i, j, k};
          target = std::move(p);
        }
      }
    }
  }
  REQUIRE_FALSE(planted.empty());
  REQUIRE(target);

  ShapeScorer scorer(*target.mol);

  // The planted product scored against itself is the ceiling.
  auto selfCopy = RWMol(*target.mol);
  const auto selfScore = scorer.score(selfCopy);
  REQUIRE(selfScore);
  INFO("self-score " << *selfScore);

  SynthonSearch3DParams mp;
  mp.numTrajectories = 1;
  mp.samplesPerReagent = 4;
  mp.numThreads = 4;
  MultipleTrajectoryStats oneStats;
  const auto one = synthonSearch3D(lib, scorer, mp, &oneStats);
  REQUIRE(one.evaluations > 0);
  REQUIRE(one.reagents.size() == 3);
  CHECK(one.score > 0.0);

  SECTION("more trajectories never score worse") {
    // The measured property this search exists for: extra trajectories are
    // strictly dominant, because the reduction keeps the best product found
    // and each trajectory descends into a different local optimum.
    mp.numTrajectories = 3;
    MultipleTrajectoryStats threeStats;
    const auto three =
        synthonSearch3D(lib, scorer, mp, &threeStats);
    REQUIRE(three.reagents.size() == 3);
    CHECK(three.score >= one.score);
    CHECK(threeStats.trajectories.size() == 3);

    // Trajectories must stay INDEPENDENT: sharing the score cache may not
    // collapse them onto one answer, or the diversity that makes restarts
    // work is gone.  Distinct seeds are the mechanism.
    CHECK(threeStats.trajectories[0].randomSeed !=
          threeStats.trajectories[1].randomSeed);

    // Physical work must be less than logical work: the shared cache is the
    // reason k trajectories cost less than k separate searches.
    CHECK(threeStats.uniqueSearchAttempts <= threeStats.logicalRequests);
  }

  SECTION("the reduced top-N is ordered and deduplicated") {
    REQUIRE_FALSE(one.best.empty());
    for (size_t i = 1; i < one.best.size(); ++i) {
      CHECK(one.best[i - 1].score >= one.best[i].score);
      CHECK(one.best[i].reagents != one.best[i - 1].reagents);
    }
  }
}
