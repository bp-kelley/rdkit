//
//  Copyright (C) 2026 The RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "NormalizedDescriptors.h"

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>

namespace RDKit {
namespace NormalizedDescriptors {

CDFTable::CDFTable(double minV, double maxV, std::vector<double> xs,
                   std::vector<double> cdf)
    : d_minV(minV), d_maxV(maxV), d_xs(std::move(xs)), d_cdf(std::move(cdf)) {
  if (d_xs.empty()) {
    throw ValueErrorException("CDFTable requires at least one point");
  }
  if (d_xs.size() != d_cdf.size()) {
    throw ValueErrorException("CDFTable x and cdf sizes differ");
  }
  if (d_maxV < d_minV) {
    throw ValueErrorException("CDFTable maxV is smaller than minV");
  }
  if (!std::is_sorted(d_xs.begin(), d_xs.end())) {
    throw ValueErrorException("CDFTable x values must be sorted");
  }
}

double CDFTable::normalize(double value) const {
  if (d_xs.empty() || !std::isfinite(value)) {
    return 0.0;
  }
  double v = std::clamp(value, d_minV, d_maxV);
  double res;
  if (v <= d_xs.front()) {
    res = d_cdf.front();
  } else if (v >= d_xs.back()) {
    res = d_cdf.back();
  } else {
    // first point with x > v; v is strictly inside the table so 0 < hi < size
    auto hi = std::upper_bound(d_xs.begin(), d_xs.end(), v) - d_xs.begin();
    auto lo = hi - 1;
    double dx = d_xs[hi] - d_xs[lo];
    double frac = dx > 0 ? (v - d_xs[lo]) / dx : 0.0;
    res = d_cdf[lo] + frac * (d_cdf[hi] - d_cdf[lo]);
  }
  if (!std::isfinite(res)) {
    return 0.0;
  }
  return std::clamp(res, 0.0, 1.0);
}

void CDFTableSet::addTable(const std::string &name, CDFTable table) {
  d_tables[name] = std::move(table);
}

bool CDFTableSet::hasTable(const std::string &name) const {
  return d_tables.find(name) != d_tables.end();
}

const CDFTable &CDFTableSet::getTable(const std::string &name) const {
  auto it = d_tables.find(name);
  if (it == d_tables.end()) {
    throw KeyErrorException(name);
  }
  return it->second;
}

std::vector<std::string> CDFTableSet::getNames() const {
  std::vector<std::string> res;
  res.reserve(d_tables.size());
  for (const auto &[name, _] : d_tables) {
    res.push_back(name);
  }
  return res;
}

double CDFTableSet::normalize(const std::string &name, double value) const {
  auto it = d_tables.find(name);
  if (it == d_tables.end()) {
    return 0.0;
  }
  return it->second.normalize(value);
}

void CDFTableSet::loadFromStream(std::istream &inStream) {
  std::string line;
  unsigned int lineNum = 0;
  while (std::getline(inStream, line)) {
    ++lineNum;
    auto start = line.find_first_not_of(" \t\r\n");
    if (start == std::string::npos || line[start] == '#') {
      continue;
    }
    std::istringstream ls(line);
    ls.imbue(std::locale::classic());
    std::string name;
    double minV, maxV;
    size_t npts;
    if (!(ls >> name >> minV >> maxV >> npts) || npts == 0) {
      std::ostringstream errout;
      errout << "bad CDF table header on line " << lineNum;
      throw ValueErrorException(errout.str());
    }
    std::vector<double> xs(npts), cdf(npts);
    for (size_t i = 0; i < npts; ++i) {
      if (!(ls >> xs[i] >> cdf[i])) {
        std::ostringstream errout;
        errout << "CDF table for " << name << " on line " << lineNum
               << " has fewer than " << npts << " points";
        throw ValueErrorException(errout.str());
      }
    }
    std::string extra;
    if (ls >> extra) {
      std::ostringstream errout;
      errout << "CDF table for " << name << " on line " << lineNum
             << " has more than " << npts << " points";
      throw ValueErrorException(errout.str());
    }
    addTable(name, CDFTable(minV, maxV, std::move(xs), std::move(cdf)));
  }
}

void CDFTableSet::loadFromFile(const std::string &fileName) {
  std::ifstream inStream(fileName);
  if (!inStream) {
    throw BadFileException("could not open CDF table file " + fileName);
  }
  loadFromStream(inStream);
}

void CDFTableSet::writeToStream(std::ostream &outStream) const {
  auto oldPrecision =
      outStream.precision(std::numeric_limits<double>::max_digits10);
  for (const auto &[name, table] : d_tables) {
    outStream << name << " " << table.getMin() << " " << table.getMax() << " "
              << table.getXs().size();
    for (size_t i = 0; i < table.getXs().size(); ++i) {
      outStream << " " << table.getXs()[i] << " " << table.getCDF()[i];
    }
    outStream << "\n";
  }
  outStream.precision(oldPrecision);
}

std::string getDefaultTablePath() {
  const char *rdbase = std::getenv("RDBASE");
  std::string base = rdbase ? rdbase : "";
  return base + "/External/NormalizedDescriptors/data/rdkit2d_cdf_v1.txt";
}

const CDFTableSet &getDefaultTables() {
  // thread-safe one-time initialization
  static const CDFTableSet tables = [] {
    CDFTableSet res;
    res.loadFromFile(getDefaultTablePath());
    return res;
  }();
  return tables;
}

}  // namespace NormalizedDescriptors
}  // namespace RDKit
