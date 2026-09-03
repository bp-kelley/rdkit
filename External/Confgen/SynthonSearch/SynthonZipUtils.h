//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_SYNTHON_ZIP_UTILS_H
#define RDKIT_SYNTHON_ZIP_UTILS_H

//! \file SynthonZipUtils.h
//! Internal helpers shared by the synthon enumerators for recognizing the
//! exit-vector ("dummy") labeling scheme used to zip synthons together and for
//! deriving the corresponding MolzipParams.  These are header-only/inline so
//! both the 2D (EnumerateSynthons) and 3D (EnumerateSynthons3D) paths use one
//! implementation.

#include <GraphMol/RDKitBase.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

#include <set>
#include <string>
#include <vector>

namespace RDKit {
namespace SynthonZip {

// How the exit vectors on a synthon are labeled so molzip can rejoin them.
enum class SynthonLabelScheme {
  None,
  Multiple,         // mixed schemes - we can't safely zip these
  IsotopeDummy,     // [1*], [2*]
  AtomMappedDummy,  // [*:1], [*:2]
  AtomType          // [U], [Np], [Pu], [Am], ...
};

struct SynthonZipperInfo {
  SynthonLabelScheme scheme = SynthonLabelScheme::None;
  std::set<std::string> atomSymbols;
};

inline std::string atomSymbol(const Atom *atom) {
  return PeriodicTable::getTable()->getElementSymbol(atom->getAtomicNum());
}

inline bool isIsotopeDummy(const Atom *atom) {
  return atom->getAtomicNum() == 0 && atom->getIsotope() != 0;
}

inline bool isAtomMappedDummy(const Atom *atom) {
  return atom->getAtomicNum() == 0 && atom->getAtomMapNum() != 0;
}

inline bool isSynthonAtomTypeDesignator(const Atom *atom) {
  switch (atom->getAtomicNum()) {
    case 92:  // U
    case 93:  // Np
    case 94:  // Pu
    case 95:  // Am
      return true;
    default:
      return false;
  }
}

//! True if the atom is any kind of synthon exit-vector marker.
inline bool isExitDummy(const Atom *atom) {
  return isIsotopeDummy(atom) || isAtomMappedDummy(atom) ||
         isSynthonAtomTypeDesignator(atom);
}

//! The value used to pair an exit vector across two synthons (same key == same
//! junction).  Returns "" if the atom is not an exit-vector marker.
inline std::string dummyLabelKey(const Atom *atom, SynthonLabelScheme scheme) {
  switch (scheme) {
    case SynthonLabelScheme::IsotopeDummy:
      return isIsotopeDummy(atom) ? std::to_string(atom->getIsotope()) : "";
    case SynthonLabelScheme::AtomMappedDummy:
      return isAtomMappedDummy(atom) ? std::to_string(atom->getAtomMapNum())
                                     : "";
    case SynthonLabelScheme::AtomType:
      return isSynthonAtomTypeDesignator(atom) ? atomSymbol(atom) : "";
    default:
      return "";
  }
}

inline SynthonZipperInfo sniffSynthons(const ROMol &mol) {
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

  const int numJoinTypes = static_cast<int>(hasIsotopeDummies) +
                           static_cast<int>(hasAtomMappedDummies) +
                           static_cast<int>(hasAtomTypes);

  if (numJoinTypes == 0) {
    info.scheme = SynthonLabelScheme::None;
  } else if (numJoinTypes > 1) {
    info.scheme = SynthonLabelScheme::Multiple;
  } else if (hasIsotopeDummies) {
    info.scheme = SynthonLabelScheme::IsotopeDummy;
  } else if (hasAtomMappedDummies) {
    info.scheme = SynthonLabelScheme::AtomMappedDummy;
  } else {
    info.scheme = SynthonLabelScheme::AtomType;
  }

  return info;
}

inline bool sniffMolzipParams(MolzipParams &params,
                              const std::vector<SynthonZipperInfo> &zippers) {
  params.enforceValenceRules = true;
  if (zippers.empty()) {
    return false;
  }

  std::set<std::string> symbols;
  for (const auto &zipper : zippers) {
    if (zipper.scheme != zippers[0].scheme) {
      return false;
    }
    symbols.insert(zipper.atomSymbols.begin(), zipper.atomSymbols.end());
  }

  switch (zippers.front().scheme) {
    case SynthonLabelScheme::IsotopeDummy:
      params.label = MolzipLabel::Isotope;
      return true;
    case SynthonLabelScheme::AtomMappedDummy:
      params.label = MolzipLabel::AtomMapNumber;
      return true;
    case SynthonLabelScheme::AtomType:
      params.label = MolzipLabel::AtomType;
      params.atomSymbols =
          std::vector<std::string>(symbols.begin(), symbols.end());
      return true;
    case SynthonLabelScheme::None:
    case SynthonLabelScheme::Multiple:
      return false;
  }
  return false;
}

inline bool sniffMolzipParams(MolzipParams &params,
                              const MOL_SPTR_VECT &synthons) {
  std::vector<SynthonZipperInfo> zippers;
  for (const auto &mol : synthons) {
    if (!mol) {
      return false;
    }
    zippers.push_back(sniffSynthons(*mol));
  }
  return sniffMolzipParams(params, zippers);
}

//! Return the scheme common to a set of synthons, or None/Multiple if they are
//! inconsistent.
inline SynthonLabelScheme getSynthonScheme(const MOL_SPTR_VECT &synthons) {
  SynthonLabelScheme scheme = SynthonLabelScheme::None;
  bool first = true;
  for (const auto &mol : synthons) {
    if (!mol) {
      return SynthonLabelScheme::Multiple;
    }
    auto s = sniffSynthons(*mol).scheme;
    if (s == SynthonLabelScheme::None) {
      continue;  // a synthon with no exit vectors doesn't constrain the scheme
    }
    if (first) {
      scheme = s;
      first = false;
    } else if (s != scheme) {
      return SynthonLabelScheme::Multiple;
    }
  }
  return scheme;
}

}  // namespace SynthonZip
}  // namespace RDKit

#endif
