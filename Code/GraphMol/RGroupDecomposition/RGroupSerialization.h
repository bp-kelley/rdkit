//
//  Copyright (C) 2024 Brian P. Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_RGROUP_SERIALIZATION_H
#define RDKIT_RGROUP_SERIALIZATION_H

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/archive/archive_exception.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "../RDKitBase.h"
#include <GraphMol/MolPickler.h>

namespace RDKit {

template<class Archive>
void save_mol(Archive &ar, const RWMOL_SPTR &mol) {
  std::string res;
  if (mol.get() != nullptr) {
    MolPickler::pickleMol(mol.get(), res, PicklerOps::AtomProps | PicklerOps::MolProps);
  }
  ar << res;
}

template<class Archive>
void restore_mol(Archive &ar, RWMOL_SPTR &mol) {
  std::string res;
  ar >> res;
  if (res.size()) {
    mol = RWMOL_SPTR(new RWMol(res));
  }
}

template<class Archive>
void save_mols(Archive &ar, const std::vector<ROMOL_SPTR> &mols) {
  std::vector<std::string> molstr;

  for(auto &mol : mols) {
    std::string res;
    MolPickler::pickleMol(mol.get(), res, PicklerOps::AtomProps | PicklerOps::MolProps);
    molstr.push_back(res);
  }
  ar << molstr;
}

template<class Archive>
void restore_mols(Archive &ar, std::vector<ROMOL_SPTR> &mols) {
  std::vector<std::string> molstr;
  ar >> molstr;
  mols.clear();
  for(auto &str : molstr) {
    mols.emplace_back(new RWMol(str));
  }
}

template<class Archive>
void save_fingerprint(Archive &ar, const std::unique_ptr<ExplicitBitVect> &bv) {
  std::string res;
  ar << res;
}

template<class Archive>
void restore_fingerprint(Archive &ar, std::unique_ptr<ExplicitBitVect> &bv) {
  std::string res;
  ar >> res;
}

}


#endif
#endif
