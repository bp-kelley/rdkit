//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_ENUMERATE_SYNTHON_H
#define RDKIT_ENUMERATE_SYNTHON_H

#include "Enumerate.h"
#include <GraphMol/ChemTransforms/ChemTransforms.h>

namespace RDKit {

  class RDKIT_CHEMREACTIONS_EXPORT EnumerateSynthons
    : public EnumerateLibrary {
  MolzipParams m_molzipParams;
  bool m_valid = false;

 public:
  EnumerateSynthons() : EnumerateLibrary() {}
  EnumerateSynthons(const std::string &s) : EnumerateLibrary() {
    initFromString(s);
  }

  EnumerateSynthons(const EnumerationTypes::BBS &reagents,
                   const EnumerationParams &params = EnumerationParams());

  EnumerateSynthons(const EnumerationTypes::BBS &reagents,
		    const EnumerationStrategyBase &enumerator,
		    const EnumerationParams &params = EnumerationParams());

  EnumerateSynthons(const EnumerateSynthons &rhs);

  std::vector<MOL_SPTR_VECT> next() override;

  void toStream(std::ostream &ss) const override;
  void initFromStream(std::istream &ss) override;


private:
#ifdef RDK_USE_BOOST_SERIALIZATION  
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int /*version*/) {
   ar &boost::serialization::base_object<EnumerateLibrary>(*this);
   ar &m_valid;
   ar &m_molzipParams.label;
   ar &m_molzipParams.atomSymbols;
   ar &m_molzipParams.atomProperty;
   ar &m_molzipParams.enforceValenceRules;
   ar &m_molzipParams.generateCoordinates;
   ar &m_molzipParams.alignCoordinates;
  }
#endif
};
  
}

#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_CLASS_VERSION(RDKit::EnumerateSynthons, 1)
#endif

#endif
