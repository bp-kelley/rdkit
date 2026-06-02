//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Enumerate.h"
#include "EnumerateSynthons.h"
#include "CartesianProduct.h"
#include "RandomSample.h"
#include "RandomSampleAllBBs.h"
#include "EvenSamplePairs.h"
#include "../ReactionPickler.h"
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/multiprecision/cpp_int.hpp>
#ifdef RDK_USE_BOOST_SERIALIZATION
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>
#endif
#include <RDGeneral/BoostEndInclude.h>

#ifdef RDK_USE_BOOST_SERIALIZATION
// Since we are exporting the classes for serialization,
//  we should declare the archives types used here
BOOST_CLASS_EXPORT(RDKit::EnumerateSynthons);
#endif

namespace RDKit {
using namespace EnumerationTypes;
namespace {

// This little bit of nonesense code tries to sniff the synthons to see
// how we can recombine them with molzip
enum class SynthonLabelScheme {
  None,
  Multiple,          // Probably throw an error
  IsotopeDummy,      // [1*], [2*]
  AtomMappedDummy,  // [*:1], [*:2]
  AtomType          // [U], [Np], [Pu], [Am], etc.
};

struct SynthonZipperInfo {
  SynthonLabelScheme scheme = SynthonLabelScheme::None;
  std::set<std::string> atomSymbols;
};

static std::string atomSymbol(const RDKit::Atom *atom) {
  return RDKit::PeriodicTable::getTable()->getElementSymbol(atom->getAtomicNum());
}

static bool isIsotopeDummy(const RDKit::Atom *atom) {
  return atom->getAtomicNum() == 0 && atom->getIsotope() != 0;
}

static bool isAtomMappedDummy(const RDKit::Atom *atom) {
  return atom->getAtomicNum() == 0 && atom->getAtomMapNum() != 0;
}

static bool isSynthonAtomTypeDesignator(const RDKit::Atom *atom) {
  // We could/should be more clever about this...
  switch (atom->getAtomicNum()) {
    case 92: // U
    case 93: // Np
    case 94: // Pu
    case 95: // Am
      return true;
    default:
      return false;
  }
}

SynthonZipperInfo sniffSynthons(const RDKit::ROMol &mol) {
  SynthonZipperInfo info;

  bool hasIsotopeDummies = false;
  bool hasAtomMappedDummies = false;
  bool hasAtomTypes = false;

  for (const auto atom : mol.atoms()) {
    if (isIsotopeDummy(atom)) {
      hasIsotopeDummies = true;
    } else if (isAtomMappedDummy(atom)) {
      hasAtomMappedDummies = true;
    } else if (isSynthonAtomTypeDesignator(atom)) {
      hasAtomTypes = true;
      info.atomSymbols.insert(atomSymbol(atom));
    }
  }

  const int numJoinTypes =
      static_cast<int>(hasIsotopeDummies) +
      static_cast<int>(hasAtomMappedDummies) +
      static_cast<int>(hasAtomTypes);

  if (numJoinTypes == 0) {
    info.scheme = SynthonLabelScheme::None;
    return info;
  }

  if (numJoinTypes > 1) {
    info.scheme = SynthonLabelScheme::Multiple;
    return info;
  }

  if (hasIsotopeDummies) {
    info.scheme = SynthonLabelScheme::IsotopeDummy;
  } else if (hasAtomMappedDummies) {
    info.scheme = SynthonLabelScheme::AtomMappedDummy;
  } else {
    info.scheme = SynthonLabelScheme::AtomType;
  }

  return info;
}

bool sniffMolzipParams(
  MolzipParams &params,
  const std::vector<SynthonZipperInfo> &zippers) {

  params.enforceValenceRules = true; // maybe not...
  if(zippers.size() == 0)
    return false;

  std::set<std::string> symbols;
  for(auto &zipper:zippers) {
    if(zipper.scheme != zippers[0].scheme)
      return false;
    symbols.insert(zipper.atomSymbols.begin(), zipper.atomSymbols.end());
  }
  auto &info = zippers.front();
  
  switch (info.scheme) {
  case SynthonLabelScheme::IsotopeDummy:
    params.label = RDKit::MolzipLabel::Isotope;
    break;
    
  case SynthonLabelScheme::AtomMappedDummy:
    params.label = RDKit::MolzipLabel::AtomMapNumber;
    break;
    
  case SynthonLabelScheme::AtomType: {
    params.label = RDKit::MolzipLabel::AtomType;
    std::vector<std::string> atomSymbols(symbols.begin(), symbols.end());
    params.atomSymbols = atomSymbols;
    break;
  }
  case SynthonLabelScheme::None:
    return false;
  case SynthonLabelScheme::Multiple:      
    return false;
  }
  
  return true;
}

bool sniffMolzipParams(MolzipParams &params,
		       const MOL_SPTR_VECT &synthons) {
  std::vector<SynthonZipperInfo> zippers;
  for(auto &mol:synthons) {
    if(!mol)
      return false;
    zippers.push_back(sniffSynthons(*mol));
  }
  return sniffMolzipParams(params, zippers);
}
}

EnumerateSynthons::EnumerateSynthons(
     const EnumerationTypes::BBS &reagents,
     const EnumerationParams &params)
  : EnumerateLibrary(ChemicalReaction(), reagents, params) {
  // make fake bbs to initialize enumerator  
  m_valid = sniffMolzipParams(m_molzipParams, m_bbs[0]);
}

EnumerateSynthons::EnumerateSynthons(
     const EnumerationTypes::BBS &reagents,
     const EnumerationStrategyBase &enumerator,
     const EnumerationParams &params)
  : EnumerateLibrary(ChemicalReaction(), reagents, enumerator, params) {
  m_valid = sniffMolzipParams(m_molzipParams, m_bbs[0]);
}

EnumerateSynthons::EnumerateSynthons(const EnumerateSynthons &rhs)
  : EnumerateLibrary(rhs),
    m_molzipParams(rhs.m_molzipParams),
    m_valid(rhs.m_valid) {
}

std::vector<MOL_SPTR_VECT> EnumerateSynthons::next() {
  PRECONDITION(static_cast<bool>(*this), "No more enumerations");
  const RGROUPS &reactantIndices = m_enumerator->next();
  if(!m_valid || m_bbs.size() == 0 || reactantIndices.size() == 0)
    return std::vector<MOL_SPTR_VECT>();

  RWMol m(*m_bbs[0][reactantIndices[0]]);
  for(size_t i=1; i<m_bbs.size(); ++i) {
    m.insertMol(*m_bbs[i][reactantIndices[i]]);
  }
  
  auto zipped = molzip(m, m_molzipParams);

  std::vector<MOL_SPTR_VECT> res(1);
  res[0].emplace_back(zipped.release());
  return res;
}

void EnumerateSynthons::toStream(std::ostream &ss) const {
#ifdef RDK_USE_BOOST_SERIALIZATION
  boost::archive::text_oarchive ar(ss);
  ar << *this;
#else
  PRECONDITION(0, "BOOST SERIALIZATION NOT INSTALLED");
#endif
}

void EnumerateSynthons::initFromStream(std::istream &ss) {
#ifdef RDK_USE_BOOST_SERIALIZATION
  boost::archive::text_iarchive ar(ss);
  ar >> *this;
#else
  PRECONDITION(0, "BOOST SERIALIZATION NOT INSTALLED");
#endif
}
  
}  // namespace RDKit
