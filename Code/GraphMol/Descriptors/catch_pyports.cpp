//
//  Copyright (C) 2026 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Tests for the C++ ports of descriptors that used to be python-only.
//  The reference values in test_data/py_descriptor_reference.tsv were
//  generated with test_data/make_py_descriptor_reference.py

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

#include <cmath>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>

using namespace RDKit;
using namespace RDKit::Descriptors;

namespace {
std::vector<std::string> splitTabs(const std::string &line) {
  std::vector<std::string> res;
  std::stringstream ss(line);
  std::string tok;
  while (std::getline(ss, tok, '\t')) {
    res.push_back(tok);
  }
  return res;
}

std::map<std::string, std::function<double(const ROMol &)>> getFunctions() {
  std::map<std::string, std::function<double(const ROMol &)>> fns = {
      {"MaxAbsEStateIndex", calcMaxAbsEStateIndex},
      {"MaxEStateIndex", calcMaxEStateIndex},
      {"MinAbsEStateIndex", calcMinAbsEStateIndex},
      {"MinEStateIndex", calcMinEStateIndex},
      {"qed", [](const ROMol &m) { return calcQED(m); }},
      {"SPS", [](const ROMol &m) { return calcSPS(m); }},
      {"HeavyAtomMolWt", calcHeavyAtomMolWt},
      {"NumValenceElectrons", calcNumValenceElectrons},
      {"NumRadicalElectrons", calcNumRadicalElectrons},
      {"MaxPartialCharge", calcMaxPartialCharge},
      {"MinPartialCharge", calcMinPartialCharge},
      {"MaxAbsPartialCharge", calcMaxAbsPartialCharge},
      {"MinAbsPartialCharge", calcMinAbsPartialCharge},
      {"FpDensityMorgan1",
       [](const ROMol &m) { return calcFpDensityMorgan(m, 1); }},
      {"FpDensityMorgan2",
       [](const ROMol &m) { return calcFpDensityMorgan(m, 2); }},
      {"FpDensityMorgan3",
       [](const ROMol &m) { return calcFpDensityMorgan(m, 3); }},
      {"AvgIpc", calcAvgIpc},
      {"BalabanJ", calcBalabanJ},
      {"BertzCT", [](const ROMol &m) { return calcBertzCT(m); }},
      {"Chi0", calcChi0},
      {"Chi1", calcChi1},
      {"Ipc", [](const ROMol &m) { return calcIpc(m); }},
  };
  for (unsigned int i = 0; i < estateBins.size() + 1; ++i) {
    fns["EState_VSA" + std::to_string(i + 1)] = [i](const ROMol &m) {
      return calcEState_VSA(m)[i];
    };
  }
  for (unsigned int i = 0; i < vsaEStateBins.size() + 1; ++i) {
    fns["VSA_EState" + std::to_string(i + 1)] = [i](const ROMol &m) {
      return calcVSA_EState(m)[i];
    };
  }
  const auto &fragNames = getFragmentDescriptorNames();
  for (const auto &nm : fragNames) {
    fns[nm] = [nm](const ROMol &m) {
      return static_cast<double>(calcFragmentDescriptor(m, nm));
    };
  }
  return fns;
}
}  // namespace

TEST_CASE("python ported descriptors match the python reference values") {
  std::string rdbase = std::getenv("RDBASE");
  std::string fname =
      rdbase +
      "/Code/GraphMol/Descriptors/test_data/py_descriptor_reference.tsv";
  std::ifstream inf(fname);
  REQUIRE(inf.good());
  std::string line;
  std::getline(inf, line);
  auto names = splitTabs(line);
  REQUIRE(names.size() > 1);
  auto fns = getFunctions();
  for (size_t i = 1; i < names.size(); ++i) {
    INFO(names[i]);
    CHECK(fns.count(names[i]) == 1);
  }
  unsigned int nMols = 0;
  while (std::getline(inf, line)) {
    auto fields = splitTabs(line);
    REQUIRE(fields.size() == names.size());
    std::unique_ptr<RWMol> m(SmilesToMol(fields[0]));
    REQUIRE(m);
    ++nMols;
    for (size_t i = 1; i < names.size(); ++i) {
      INFO(fields[0] << " " << names[i]);
      auto ref = std::stod(fields[i]);
      auto val = fns[names[i]](*m);
      if (std::isnan(ref)) {
        // Gasteiger charges are NaN for some species
        CHECK(std::isnan(val));
        continue;
      }
      CHECK_THAT(val, Catch::Matchers::WithinRel(ref, 1e-8) ||
                          Catch::Matchers::WithinAbs(ref, 1e-10));
    }
  }
  CHECK(nMols > 50);
}

TEST_CASE("simple python-ported descriptors") {
  SECTION("NumValenceElectrons") {
    std::vector<std::pair<std::string, int>> data = {
        {"CC", 14}, {"C(=O)O", 18}, {"C(=O)[O-]", 18}, {"C(=O)", 12}};
    for (const auto &[smi, val] : data) {
      std::unique_ptr<RWMol> m(SmilesToMol(smi));
      REQUIRE(m);
      CHECK(calcNumValenceElectrons(*m) == val);
    }
  }
  SECTION("NumRadicalElectrons") {
    std::vector<std::pair<std::string, unsigned int>> data = {
        {"CC", 0}, {"C[CH3]", 0}, {"C[CH2]", 1}, {"C[CH]", 2}, {"C[C]", 3}};
    for (const auto &[smi, val] : data) {
      std::unique_ptr<RWMol> m(SmilesToMol(smi));
      REQUIRE(m);
      CHECK(calcNumRadicalElectrons(*m) == val);
    }
  }
  SECTION("HeavyAtomMolWt") {
    std::unique_ptr<RWMol> m(SmilesToMol("[NH4+].[Cl-]"));
    REQUIRE(m);
    CHECK_THAT(calcHeavyAtomMolWt(*m), Catch::Matchers::WithinAbs(49.46, 0.01));
  }
  SECTION("empty molecules") {
    RWMol m;
    CHECK(calcChi0(m) == 0.0);
    CHECK(calcBalabanJ(m) == 0.0);
    CHECK(calcBertzCT(m) == 0.0);
    CHECK(calcIpc(m) == 0.0);
    CHECK(calcMaxEStateIndex(m) == 0.0);
    CHECK(calcMaxPartialCharge(m) == -500.0);
    CHECK(calcMinPartialCharge(m) == 500.0);
    CHECK(calcFpDensityMorgan(m, 2) == 0.0);
    CHECK(calcSPS(m) == 0.0);
  }
}

TEST_CASE("QED") {
  // examples from the QED paper, reference values from Peter G's original
  // implementation (also used in the python doctests)
  std::vector<std::pair<std::string, double>> data = {
      {"N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O", 0.253},
      {"CNC(=NCCSCc1nc[nH]c1C)NC#N", 0.234},
      {"CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12", 0.234}};
  for (const auto &[smi, val] : data) {
    std::unique_ptr<RWMol> m(SmilesToMol(smi));
    REQUIRE(m);
    CHECK_THAT(calcQED(*m), Catch::Matchers::WithinAbs(val, 0.001));
  }
  SECTION("weights") {
    std::unique_ptr<RWMol> m(SmilesToMol("CNC(=NCCSCc1nc[nH]c1C)NC#N"));
    REQUIRE(m);
    auto props = calcQEDProperties(*m);
    CHECK(calcQED(props) == calcQED(*m));
    CHECK(calcQED(props, QEDWeightsMax) != calcQED(props, QEDWeightsNone));
  }
}

TEST_CASE("fragment descriptors") {
  const auto &names = getFragmentDescriptorNames();
  const auto &smarts = getFragmentDescriptorSmarts();
  REQUIRE(names.size() == 85);
  REQUIRE(smarts.size() == names.size());
  CHECK(names.front() == "fr_C_O");

  SECTION("table matches Data/FragmentDescriptors.csv") {
    std::string rdbase = std::getenv("RDBASE");
    std::ifstream inf(rdbase + "/Data/FragmentDescriptors.csv");
    REQUIRE(inf.good());
    std::string line;
    size_t idx = 0;
    while (std::getline(inf, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }
      if (line.back() == '\r') {
        line.pop_back();
      }
      auto fields = splitTabs(line);
      if (fields.size() < 3) {
        continue;
      }
      auto name = fields[0];
      std::replace(name.begin(), name.end(), '=', '_');
      std::replace(name.begin(), name.end(), '-', '_');
      REQUIRE(idx < names.size());
      CHECK(names[idx] == name);
      CHECK(smarts[idx] == fields[2]);
      ++idx;
    }
    CHECK(idx == names.size());
  }
  SECTION("counts") {
    std::unique_ptr<RWMol> m(SmilesToMol("OC(=O)c1ccccc1OC(=O)C"));
    REQUIRE(m);
    CHECK(calcFragmentDescriptor(*m, "fr_benzene") == 1);
    CHECK(calcFragmentDescriptor(*m, "fr_ester") == 1);
    CHECK(calcFragmentDescriptor(*m, "fr_COO") == 1);
    auto all = calcFragmentDescriptors(*m);
    REQUIRE(all.size() == names.size());
    CHECK(all[0] == calcFragmentDescriptor(*m, "fr_C_O"));
    CHECK_THROWS_AS(calcFragmentDescriptor(*m, "fr_not_there"),
                    KeyErrorException);
  }
}
