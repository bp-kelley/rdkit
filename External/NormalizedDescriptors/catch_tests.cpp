//
//  Copyright (C) 2026 The RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>
#include "NormalizedDescriptors.h"

using namespace RDKit::NormalizedDescriptors;
using Catch::Matchers::WithinAbs;

TEST_CASE("CDFTable interpolation") {
  CDFTable table(0.0, 10.0, {0.0, 2.0, 4.0, 8.0}, {0.1, 0.5, 0.7, 0.9});
  SECTION("table points") {
    CHECK_THAT(table.normalize(0.0), WithinAbs(0.1, 1e-12));
    CHECK_THAT(table.normalize(2.0), WithinAbs(0.5, 1e-12));
    CHECK_THAT(table.normalize(8.0), WithinAbs(0.9, 1e-12));
  }
  SECTION("linear interpolation") {
    CHECK_THAT(table.normalize(1.0), WithinAbs(0.3, 1e-12));
    CHECK_THAT(table.normalize(3.0), WithinAbs(0.6, 1e-12));
    CHECK_THAT(table(6.0), WithinAbs(0.8, 1e-12));
  }
  SECTION("clipping to min/max and outside the table") {
    CHECK_THAT(table.normalize(-5.0), WithinAbs(0.1, 1e-12));
    CHECK_THAT(table.normalize(9.0), WithinAbs(0.9, 1e-12));
    CHECK_THAT(table.normalize(1e6), WithinAbs(0.9, 1e-12));
  }
  SECTION("non-finite values") {
    CHECK(table.normalize(std::numeric_limits<double>::quiet_NaN()) == 0.0);
    CHECK(table.normalize(std::numeric_limits<double>::infinity()) == 0.0);
  }
  SECTION("results are clipped to [0, 1]") {
    CDFTable bad(0.0, 1.0, {0.0, 1.0}, {-0.5, 1.5});
    CHECK(bad.normalize(0.0) == 0.0);
    CHECK(bad.normalize(1.0) == 1.0);
    CHECK_THAT(bad.normalize(0.5), WithinAbs(0.5, 1e-12));
  }
  SECTION("min clipping applies before the table lookup") {
    CDFTable clipped(1.0, 3.0, {0.0, 4.0}, {0.0, 1.0});
    CHECK_THAT(clipped.normalize(0.0), WithinAbs(0.25, 1e-12));
    CHECK_THAT(clipped.normalize(4.0), WithinAbs(0.75, 1e-12));
  }
  SECTION("single point table") {
    CDFTable single(0.0, 0.0, {0.0}, {0.42});
    CHECK_THAT(single.normalize(-1.0), WithinAbs(0.42, 1e-12));
    CHECK_THAT(single.normalize(12.0), WithinAbs(0.42, 1e-12));
  }
  SECTION("invalid tables") {
    CHECK_THROWS_AS(CDFTable(0.0, 1.0, {}, {}), ValueErrorException);
    CHECK_THROWS_AS(CDFTable(0.0, 1.0, {0.0, 1.0}, {0.0}), ValueErrorException);
    CHECK_THROWS_AS(CDFTable(0.0, 1.0, {1.0, 0.0}, {0.0, 1.0}),
                    ValueErrorException);
    CHECK_THROWS_AS(CDFTable(1.0, 0.0, {0.0, 1.0}, {0.0, 1.0}),
                    ValueErrorException);
  }
}

TEST_CASE("CDFTableSet") {
  std::string text = R"TXT(# a comment
descriptor foo norm 0 10 3
0 0
# comments are allowed between points
5 0.5

10 1
descriptor bar tabulated -1 1 2
-1 0.25
1 0.75
)TXT";
  CDFTableSet tables;
  std::istringstream inStream(text);
  tables.loadFromStream(inStream);
  REQUIRE(tables.size() == 2);
  CHECK(tables.getNames() == std::vector<std::string>{"bar", "foo"});
  CHECK(tables.hasTable("foo"));
  CHECK_THAT(tables.normalize("foo", 2.5), WithinAbs(0.25, 1e-12));
  CHECK_THAT(tables.normalize("bar", 0.0), WithinAbs(0.5, 1e-12));
  CHECK(tables.getTable("foo").getDistribution() == "norm");

  SECTION("missing descriptors normalize to 0") {
    CHECK(!tables.hasTable("baz"));
    CHECK(tables.normalize("baz", 1.0) == 0.0);
    CHECK_THROWS_AS(tables.getTable("baz"), KeyErrorException);
  }
  SECTION("round trip") {
    std::ostringstream outStream;
    tables.writeToStream(outStream);
    CDFTableSet tables2;
    std::istringstream inStream2(outStream.str());
    tables2.loadFromStream(inStream2);
    REQUIRE(tables2.getNames() == tables.getNames());
    for (const auto &name : tables.getNames()) {
      CHECK(tables2.getTable(name).getXs() == tables.getTable(name).getXs());
      CHECK(tables2.getTable(name).getCDF() == tables.getTable(name).getCDF());
      CHECK(tables2.getTable(name).getDistribution() ==
            tables.getTable(name).getDistribution());
    }
  }
  SECTION("malformed input") {
    for (const auto &bad : {
             "foo norm 0 1 1\n0 0",               // missing keyword
             "descriptor foo norm 0 1",           // short header
             "descriptor foo norm 0 1 0",         // no points
             "descriptor foo norm 0 1 2\n0 0",    // too few points
             "descriptor foo norm 0 1 1\n0 0 1",  // extra field
             "descriptor foo norm a 1 1\n0 0",    // bad number
             "descriptor foo norm 0 1 1 7\n0 0",  // extra header field
         }) {
      INFO(bad);
      CDFTableSet badTables;
      std::istringstream badStream(bad);
      CHECK_THROWS_AS(badTables.loadFromStream(badStream), ValueErrorException);
    }
  }
}

TEST_CASE("default tables match descriptastorus") {
  const auto &tables = getDefaultTables();
  CHECK(tables.size() >= 201);
  // reference values from descriptastorus' RDKit2DNormalized, i.e.
  // scipy.stats.<dist>.cdf(clip(v, minV, maxV), *params)
  struct {
    const char *name;
    double value;
    double expected;
  } refs[] = {
      {"MolLogP", 2.5, 0.28147344406620506},
      {"MolLogP", -100.0, 1.737987092370867e-06},
      {"MolLogP", 100.0, 0.9999973063987906},
      {"NumHDonors", 2.0, 0.7316602057680455},
      {"ExactMolWt", 350.1, 0.34378524188563014},
      {"TPSA", 75.3, 0.5290342513421387},
      {"fr_benzene", 1.0, 0.35822936499573854},
      {"qed", 0.6, 0.506022609149951},
  };
  for (const auto &ref : refs) {
    INFO(ref.name << " " << ref.value);
    CHECK_THAT(tables.normalize(ref.name, ref.value),
               WithinAbs(ref.expected, 1e-4));
  }
}

namespace {
std::vector<std::string> splitTabs(const std::string &line) {
  std::vector<std::string> res;
  std::istringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) {
    res.push_back(field);
  }
  return res;
}
}  // namespace

TEST_CASE("calcNormalizedDescriptors matches rdkit.Chem.Descriptors") {
  const auto &names = getNormalizedDescriptorNames();
  REQUIRE(names.size() == 217);
  CHECK(names.front() == "MaxAbsEStateIndex");
  CHECK(names.back() == "fr_urea");

  const char *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase);
  std::ifstream inStream(
      std::string(rdbase) +
      "/External/NormalizedDescriptors/test_data/reference_descriptors.tsv");
  REQUIRE(inStream);
  std::string line;
  REQUIRE(std::getline(inStream, line));
  auto header = splitTabs(line);
  REQUIRE(header.size() == names.size() + 1);
  for (size_t i = 0; i < names.size(); ++i) {
    CHECK(header[i + 1] == names[i]);
  }

  const auto &tables = getDefaultTables();
  unsigned int nMols = 0;
  while (std::getline(inStream, line)) {
    auto fields = splitTabs(line);
    REQUIRE(fields.size() == names.size() + 1);
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(fields[0]));
    REQUIRE(mol);
    ++nMols;
    auto raw = calcDescriptorValues(*mol);
    auto normalized = calcNormalizedDescriptors(*mol);
    REQUIRE(raw.size() == names.size());
    REQUIRE(normalized.size() == names.size());
    for (size_t i = 0; i < names.size(); ++i) {
      INFO(fields[0] << " " << names[i]);
      double expected = std::stod(fields[i + 1]);
      double tol = 1e-4 * std::max(1.0, std::fabs(expected));
      CHECK_THAT(raw[i], WithinAbs(expected, tol));
      CHECK_THAT(normalized[i],
                 WithinAbs(tables.normalize(names[i], expected), 1e-4));
      CHECK(normalized[i] >= 0.0);
      CHECK(normalized[i] <= 1.0);
    }
  }
  CHECK(nMols == 10);
}
