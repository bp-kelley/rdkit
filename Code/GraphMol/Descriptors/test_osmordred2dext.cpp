//
//  Copyright (C) 2026 Osmo Labs
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Osmordred2DExt.h"
#include <catch2/catch_all.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace RDKit;
namespace E = RDKit::Descriptors::Osmordred2DExt;

namespace {
double valueOf(const std::vector<double> &values, const std::string &name) {
  static const auto names = E::getOsmordred2DExtDescriptorNames();
  for (size_t i = 0; i < names.size(); ++i) {
    if (names[i] == name) {
      return values[i];
    }
  }
  FAIL("no descriptor named " << name);
  return 0.0;
}
}  // namespace

TEST_CASE("Osmordred extension: sizes and names") {
  const auto names = E::getOsmordred2DExtDescriptorNames();
  REQUIRE(names.size() == E::OSMORDRED_2D_EXT_NUM_DESCRIPTORS);
  REQUIRE(names.size() == 2533);
  std::set<std::string> unique(names.begin(), names.end());
  CHECK(unique.size() == names.size());

  std::unique_ptr<ROMol> mol(SmilesToMol("CC(=O)Nc1ccc(O)cc1"));
  REQUIRE(mol);
  CHECK(E::calcOsmordred2DExt(*mol).size() == names.size());
  CHECK(E::calcNarumi(*mol).size() == 3);
  CHECK(E::calcEdgeAdjacency(*mol).size() == 60);
  CHECK(E::calcBurdenEigenvalues(*mol).size() == 96);
  CHECK(E::calcCATS2D(*mol).size() == 150);
  CHECK(E::calcMatrix2D(*mol).size() == 608);
  CHECK(E::calcMDE(*mol).size() == 19);
  CHECK(E::calcAtomPairs2D(*mol).size() == 1596);
}

TEST_CASE("Osmordred extension: single heavy atom") {
  std::unique_ptr<ROMol> mol(SmilesToMol("C"));
  REQUIRE(mol);
  CHECK(E::calcOsmordred2DExt(*mol).size() == E::OSMORDRED_2D_EXT_NUM_DESCRIPTORS);
}

TEST_CASE("Osmordred extension: nCIR counts connected cycles, not SSSR rings") {
  // A disconnected union of two rings is not a circuit: diphenyl ether has 2, not 3.
  const std::vector<std::pair<std::string, unsigned int>> cases = {
      {"c1ccccc1", 1}, {"c1ccc2ccccc2c1", 3}, {"c1ccc(Oc2ccccc2)cc1", 2}};
  for (const auto &[smi, expected] : cases) {
    std::unique_ptr<ROMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    CHECK(E::calcNumCircuits(*mol) == expected);
  }
}

TEST_CASE("Osmordred extension: reference values") {
  struct Ref {
    std::string smiles;
    double gnar, sm02bo, sm02ed, spmax1m, spmin1m, mdec22, mdeo11;
  };
  // Values from the build validated against Dragon/alvaDesc reference output.
  const std::vector<Ref> refs = {
      {"c1ccccc1", 2.0, 3.2771447329921766, 3.6109179126442239, 3.7830235262153682,
       2.0222262431597509, 9.125465128398087, 0.0},
      {"c1ccc2ccccc2c1", 2.1689435423953971, 3.9843436670077725, 4.6539603501575231,
       4.0176605051591734, 2.1597684430848623, 11.678631656953222, 0.0},
      {"c1ccc(Oc2ccccc2)cc1", 2.1287318524121646, 4.1588830833596724, 4.7957905455967387,
       3.8789979008367657, 2.0250907205016095, 13.995261637820169, 0.0},
      {"CC(=O)Nc1ccc(O)cc1", 1.849080398343679, 3.9219733362813125, 4.6347289882296385,
       3.8550708588900937, 1.9879530643983032, 3.3019272488946267, 0.14285714285714288},
  };
  for (const auto &r : refs) {
    INFO(r.smiles);
    std::unique_ptr<ROMol> mol(SmilesToMol(r.smiles));
    REQUIRE(mol);
    const auto v = E::calcOsmordred2DExt(*mol);
    CHECK(valueOf(v, "GNar") == Catch::Approx(r.gnar).epsilon(1e-10));
    CHECK(valueOf(v, "SM02_AEA(bo)") == Catch::Approx(r.sm02bo).epsilon(1e-10));
    CHECK(valueOf(v, "SM02_AEA(ed)") == Catch::Approx(r.sm02ed).epsilon(1e-10));
    CHECK(valueOf(v, "SpMax1_Bh(m)") == Catch::Approx(r.spmax1m).epsilon(1e-10));
    CHECK(valueOf(v, "SpMin1_Bh(m)") == Catch::Approx(r.spmin1m).epsilon(1e-10));
    CHECK(valueOf(v, "MDEC-22") == Catch::Approx(r.mdec22).epsilon(1e-10));
    CHECK(valueOf(v, "MDEO-11") == Catch::Approx(r.mdeo11).margin(1e-12));
  }
}
