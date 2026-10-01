//
//  Copyright (C) 2026 The RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  CDF-table normalization of molecular descriptors, ported from the
//  RDKit2DNormalized descriptors of descriptastorus
//  (https://github.com/bp-kelley/descriptastorus).
//
#include <RDGeneral/export.h>
#ifndef RDKIT_NORMALIZEDDESCRIPTORS_H
#define RDKIT_NORMALIZEDDESCRIPTORS_H

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

namespace RDKit {
namespace NormalizedDescriptors {

//! A tabulated cumulative distribution function for a single descriptor.
/*!
  The table is a set of (x, cdf(x)) points sorted by x. A descriptor value v is
  normalized by:
    1. clipping v to [minV, maxV],
    2. linearly interpolating the CDF between the bracketing table points
       (values outside the tabulated x range take the first/last cdf value),
    3. clipping the result to [0, 1].
  Non-finite values normalize to 0.0.

  This reproduces descriptastorus' RDKit2DNormalized, which evaluates a fitted
  scipy.stats distribution's CDF on the clipped value; the table is that CDF
  sampled on a grid (see tools/fit_normalized_descriptors.py).
*/
class RDKIT_NORMALIZEDDESCRIPTORS_EXPORT CDFTable {
 public:
  CDFTable() = default;
  //! \c xs must be non-empty, sorted ascending, and the same size as \c cdf
  CDFTable(double minV, double maxV, std::vector<double> xs,
           std::vector<double> cdf, std::string distribution = "");

  //! returns the normalized value, in [0, 1], for \c value
  double normalize(double value) const;
  double operator()(double value) const { return normalize(value); }

  double getMin() const { return d_minV; }
  double getMax() const { return d_maxV; }
  const std::vector<double> &getXs() const { return d_xs; }
  const std::vector<double> &getCDF() const { return d_cdf; }
  //! name of the fitted scipy.stats distribution the table was sampled
  //! from (informational only)
  const std::string &getDistribution() const { return d_distribution; }

 private:
  double d_minV = 0.0;
  double d_maxV = 0.0;
  std::vector<double> d_xs;
  std::vector<double> d_cdf;
  std::string d_distribution;
};

//! A named collection of CDF tables, one per descriptor.
/*!
  The text format read and written by this class (the one produced by
  tools/fit_normalized_descriptors.py) has, for each descriptor, a header line
  followed by npoints lines of points sorted by x:

    descriptor <name> <scipy distribution> <min> <max> <npoints>
    <x> <cdf>
    ...

  Fields are whitespace separated. Blank lines and lines starting with '#'
  are ignored.
*/
class RDKIT_NORMALIZEDDESCRIPTORS_EXPORT CDFTableSet {
 public:
  CDFTableSet() = default;

  //! adds (or replaces) the table for descriptor \c name
  void addTable(const std::string &name, CDFTable table);
  bool hasTable(const std::string &name) const;
  //! throws a KeyErrorException if \c name has no table
  const CDFTable &getTable(const std::string &name) const;
  std::vector<std::string> getNames() const;
  size_t size() const { return d_tables.size(); }

  //! normalizes \c value with the table for \c name.
  /*!
    Following descriptastorus, a descriptor without a table normalizes to 0.0.
  */
  double normalize(const std::string &name, double value) const;

  //! reads tables from \c inStream, adding them to this set
  /*! throws a ValueErrorException on malformed input */
  void loadFromStream(std::istream &inStream);
  //! reads tables from the file \c fileName, adding them to this set
  void loadFromFile(const std::string &fileName);
  //! writes the tables in this set to \c outStream (sorted by name)
  void writeToStream(std::ostream &outStream) const;

 private:
  std::map<std::string, CDFTable> d_tables;
};

//! returns the path of the CDF tables fitted to descriptastorus'
//! RDKit2DNormalized distributions:
//!   $RDBASE/External/NormalizedDescriptors/data/normalized_descriptor_cdfs.txt
RDKIT_NORMALIZEDDESCRIPTORS_EXPORT std::string getDefaultTablePath();

//! returns the default CDF tables, loading them from getDefaultTablePath()
//! on first use
RDKIT_NORMALIZEDDESCRIPTORS_EXPORT const CDFTableSet &getDefaultTables();

}  // namespace NormalizedDescriptors
}  // namespace RDKit

#endif
