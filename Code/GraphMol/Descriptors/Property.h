//
//  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_PROPERTIES_H
#define RDKIT_PROPERTIES_H

#include <GraphMol/RDKitBase.h>
#include <limits>
#include <string>
#include <utility>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/shared_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <Query/Query.h>
#include <RDGeneral/Exceptions.h>

namespace RDKit {
namespace Descriptors {
struct RDKIT_DESCRIPTORS_EXPORT PropertyFunctor {
  // Registry of property functions
  //  See REGISTER_DESCRIPTOR
  std::string propName;
  std::string propVersion;
  double (*d_dataFunc)(const ROMol &);

  PropertyFunctor(std::string name, std::string version,
                  double (*func)(const ROMol &) = nullptr)
      : propName(std::move(name)),
        propVersion(std::move(version)),
        d_dataFunc(func) {}
  virtual ~PropertyFunctor() {}

  //! Compute the value of the property
  virtual double operator()(const RDKit::ROMol &mol) const {
    if (d_dataFunc == nullptr) {
      throw ValueErrorException("PropertyFunctor has no data function");
    }
    return (*d_dataFunc)(mol);
  }

  //! Return the name of the property
  const std::string getName() const { return propName; }
  //! Return the properties version
  const std::string getVersion() const { return propVersion; }
};

//! A property that is one element of a descriptor vector (e.g. one bin of
//! PEOE_VSA, or one BCUT2D value).
/*!
  Used on its own it computes the whole vector and returns its element.
  Properties::computeProperties() computes each vector once per molecule and
  shares it between the properties that are elements of it.
*/
struct RDKIT_DESCRIPTORS_EXPORT VectorElementPropertyFunctor
    : public PropertyFunctor {
  using VectorFunc = std::vector<double> (*)(const ROMol &);
  VectorFunc d_vectorFunc;
  unsigned int d_index;

  //! \c elementFunc must return element \c index of \c vectorFunc; it is
  //! the plain function used for property queries.
  VectorElementPropertyFunctor(std::string name, std::string version,
                               VectorFunc vectorFunc, unsigned int index,
                               double (*elementFunc)(const ROMol &))
      : PropertyFunctor(std::move(name), std::move(version), elementFunc),
        d_vectorFunc(vectorFunc),
        d_index(index) {}
};

//! Holds a collection of properties for computation purposes
class RDKIT_DESCRIPTORS_EXPORT Properties {
 protected:
  std::vector<boost::shared_ptr<PropertyFunctor>> m_properties;
  double d_failureValue = std::numeric_limits<double>::quiet_NaN();

  //! computes the properties for \c mol, giving \c failureValue to those
  //! that fail
  std::vector<double> computeValues(const RDKit::ROMol &mol,
                                    double failureValue) const;

 public:
  Properties();
  Properties(const std::vector<std::string> &propNames);
  virtual ~Properties() = default;

  std::vector<std::string> getPropertyNames() const;
  //! computes the properties for \c mol
  /*!
    A property whose calculation throws a std::exception gets the failure
    value (see setFailureValue()) instead.
  */
  virtual std::vector<double> computeProperties(const RDKit::ROMol &mol,
                                                bool annotate = false) const;
  void annotateProperties(RDKit::ROMol &mol) const;

  //! sets the value given to a property that fails to compute (default NaN)
  void setFailureValue(double val) { d_failureValue = val; }
  double getFailureValue() const { return d_failureValue; }

  //! Register a property function - takes ownership
  static int registerProperty(PropertyFunctor *ptr);
  static int registerProperty(boost::shared_ptr<PropertyFunctor> prop);
  static boost::shared_ptr<PropertyFunctor> getProperty(
      const std::string &name);
  static std::vector<std::string> getAvailableProperties();
  static std::vector<boost::shared_ptr<PropertyFunctor>> registry;
};

typedef Queries::Query<bool, const ROMol &, true> PROP_BOOL_QUERY;
typedef Queries::AndQuery<int, const ROMol &, true> PROP_AND_QUERY;
typedef Queries::OrQuery<int, const ROMol &, true> PROP_OR_QUERY;
typedef Queries::XOrQuery<int, const ROMol &, true> PROP_XOR_QUERY;

typedef Queries::EqualityQuery<double, const ROMol &, true> PROP_EQUALS_QUERY;

typedef Queries::GreaterQuery<double, const ROMol &, true> PROP_GREATER_QUERY;

typedef Queries::GreaterEqualQuery<double, const ROMol &, true>
    PROP_GREATEREQUAL_QUERY;

typedef Queries::LessQuery<double, const ROMol &, true> PROP_LESS_QUERY;

typedef Queries::LessEqualQuery<double, const ROMol &, true>
    PROP_LESSEQUAL_QUERY;

typedef Queries::RangeQuery<double, const ROMol &, true> PROP_RANGE_QUERY;

template <class T>
T *makePropertyQuery(const std::string &name, double what) {
  T *t = new T(what);
  t->setDataFunc(Properties::getProperty(name)->d_dataFunc);
  return t;
}

RDKIT_DESCRIPTORS_EXPORT PROP_RANGE_QUERY *makePropertyRangeQuery(
    const std::string &name, double min, double max);

}  // namespace Descriptors
}  // namespace RDKit
#endif
