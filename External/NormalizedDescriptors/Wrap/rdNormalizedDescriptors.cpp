//
//  Copyright (C) 2026 The RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <sstream>
#include <string>
#include <vector>

#include "../NormalizedDescriptors.h"

namespace python = boost::python;
using namespace RDKit::NormalizedDescriptors;

namespace {
CDFTable *makeCDFTable(double minV, double maxV, const python::object &xs,
                       const python::object &cdf) {
  std::vector<double> xsv, cdfv;
  pythonObjectToVect<double>(xs, xsv);
  pythonObjectToVect<double>(cdf, cdfv);
  return new CDFTable(minV, maxV, std::move(xsv), std::move(cdfv));
}

python::tuple getXs(const CDFTable &table) {
  python::list res;
  for (auto v : table.getXs()) {
    res.append(v);
  }
  return python::tuple(res);
}

python::tuple getCDF(const CDFTable &table) {
  python::list res;
  for (auto v : table.getCDF()) {
    res.append(v);
  }
  return python::tuple(res);
}

python::list getNames(const CDFTableSet &tables) {
  python::list res;
  for (const auto &name : tables.getNames()) {
    res.append(name);
  }
  return res;
}

void loadFromString(CDFTableSet &tables, const std::string &text) {
  std::istringstream inStream(text);
  tables.loadFromStream(inStream);
}

std::string toString(const CDFTableSet &tables) {
  std::ostringstream outStream;
  tables.writeToStream(outStream);
  return outStream.str();
}

double normalizeWithDefaults(const std::string &name, double value) {
  return getDefaultTables().normalize(name, value);
}
}  // namespace

BOOST_PYTHON_MODULE(rdNormalizedDescriptors) {
  python::scope().attr("__doc__") =
      "Module containing CDF-table normalization of molecular descriptors, a "
      "port of the RDKit2DNormalized descriptors from descriptastorus";

  python::class_<CDFTable>(
      "CDFTable",
      "A tabulated cumulative distribution function for one descriptor.\n"
      "Values are clipped to [minV, maxV], linearly interpolated in the "
      "table and clipped to [0, 1]. Non-finite values normalize to 0.0.",
      python::no_init)
      .def("__init__",
           python::make_constructor(
               makeCDFTable, python::default_call_policies(),
               (python::arg("minV"), python::arg("maxV"), python::arg("xs"),
                python::arg("cdf"))),
           "Constructor. xs must be sorted and the same length as cdf")
      .def("Normalize", &CDFTable::normalize,
           (python::arg("self"), python::arg("value")),
           "returns the normalized value, in [0, 1]")
      .def("__call__", &CDFTable::normalize,
           (python::arg("self"), python::arg("value")))
      .def("GetMin", &CDFTable::getMin, python::arg("self"))
      .def("GetMax", &CDFTable::getMax, python::arg("self"))
      .def("GetXs", getXs, python::arg("self"))
      .def("GetCDF", getCDF, python::arg("self"));

  python::class_<CDFTableSet>("CDFTableSet",
                              "A collection of CDF tables keyed by "
                              "descriptor name",
                              python::init<>(python::args("self")))
      .def("AddTable", &CDFTableSet::addTable,
           (python::arg("self"), python::arg("name"), python::arg("table")),
           "adds (or replaces) the table for a descriptor")
      .def("HasTable", &CDFTableSet::hasTable,
           (python::arg("self"), python::arg("name")))
      .def("GetTable", &CDFTableSet::getTable,
           (python::arg("self"), python::arg("name")),
           python::return_internal_reference<1>(),
           "returns the table for a descriptor, raises KeyError if missing")
      .def("GetNames", getNames, python::arg("self"),
           "returns the sorted descriptor names")
      .def("Normalize", &CDFTableSet::normalize,
           (python::arg("self"), python::arg("name"), python::arg("value")),
           "normalizes a descriptor value, returns 0.0 for descriptors "
           "without a table")
      .def("LoadFromFile", &CDFTableSet::loadFromFile,
           (python::arg("self"), python::arg("fileName")),
           "adds the tables in a file to this set")
      .def("LoadFromString", loadFromString,
           (python::arg("self"), python::arg("text")),
           "adds the tables in a string to this set")
      .def("ToString", toString, python::arg("self"),
           "returns the tables in the text format read by LoadFromString")
      .def("__len__", &CDFTableSet::size, python::arg("self"));

  python::def("GetDefaultTablePath", getDefaultTablePath,
              "returns the path of the default CDF tables");
  python::def("GetDefaultTables", getDefaultTables,
              python::return_value_policy<python::reference_existing_object>(),
              "returns the default CDF tables (fitted to descriptastorus' "
              "RDKit2DNormalized distributions)");
  python::def("NormalizeDescriptor", normalizeWithDefaults,
              (python::arg("name"), python::arg("value")),
              "normalizes a descriptor value with the default tables, returns "
              "0.0 for descriptors without a table");
}
