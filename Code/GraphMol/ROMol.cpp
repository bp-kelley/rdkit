//
//  Copyright (C) 2003-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <iostream>

// our stuff
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "ROMol.h"
#include "Atom.h"
#include "QueryAtom.h"
#include "Bond.h"
#include "QueryBond.h"
#include "MolPickler.h"
#include "Conformer.h"
#include "SubstanceGroup.h"

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

namespace RDKit {
class QueryAtom;
class QueryBond;

const int ci_RIGHTMOST_ATOM = -0xBADBEEF;
const int ci_LEADING_BOND = -0xBADBEEF + 1;
const int ci_ATOM_HOLDER = -0xDEADD06;

void ROMol::destroy() {
  d_atomBookmarks.clear();
  d_bondBookmarks.clear();

    for(auto atom: _atoms) {
        delete atom;
    }
    for(auto bond: _bonds) {
        delete bond;
    }
    _atoms.clear();
    _bonds.clear();
    
  delete dp_ringInfo;

  d_sgroups.clear();
  d_stereo_groups.clear();
}

ROMol::ROMol(const std::string &pickle) : RDProps() {
  initMol();
  numBonds = 0;
  MolPickler::molFromPickle(pickle, *this);
    numBonds = static_cast<unsigned int>(_bonds.size());
}

ROMol::ROMol(const std::string &pickle, unsigned int propertyFlags)
    : RDProps() {
  initMol();
  numBonds = 0;
  MolPickler::molFromPickle(pickle, *this, propertyFlags);
        numBonds = static_cast<unsigned int>(_bonds.size());
}

void ROMol::initFromOther(const ROMol &other, bool quickCopy, int confId) {
  if (this == &other) {
    return;
  }
  numBonds = 0;
  // std::cerr<<"    init from other: "<<this<<" "<<&other<<std::endl;
  // copy over the atoms
  for (const auto oatom : other.atoms()) {
    constexpr bool updateLabel = false;
    constexpr bool takeOwnership = true;
    addAtom(oatom->copy(), updateLabel, takeOwnership);
  }

  // and the bonds:
  for (const auto obond : other.bonds()) {
    addBond(obond->copy(), true);
  }

  // ring information
  delete dp_ringInfo;
  if (other.dp_ringInfo) {
    dp_ringInfo = new RingInfo(*(other.dp_ringInfo));
  } else {
    dp_ringInfo = new RingInfo();
  }

  // enhanced stereochemical information
  d_stereo_groups.clear();
  for (auto &otherGroup : other.d_stereo_groups) {
    std::vector<Atom *> atoms;
    for (auto &otherAtom : otherGroup.getAtoms()) {
      atoms.push_back(getAtomWithIdx(otherAtom->getIdx()));
    }
    std::vector<Bond *> bonds;
    for (auto &otherBond : otherGroup.getBonds()) {
      bonds.push_back(getBondWithIdx(otherBond->getIdx()));
    }
    d_stereo_groups.emplace_back(otherGroup.getGroupType(), std::move(atoms),
                                 std::move(bonds), otherGroup.getReadId());
    d_stereo_groups.back().setWriteId(otherGroup.getWriteId());
  }

  if (other.dp_delAtoms) {
    dp_delAtoms.reset(new boost::dynamic_bitset<>(*other.dp_delAtoms));
  } else {
    dp_delAtoms.reset(nullptr);
  }
  if (other.dp_delBonds) {
    dp_delBonds.reset(new boost::dynamic_bitset<>(*other.dp_delBonds));
  } else {
    dp_delBonds.reset(nullptr);
  }

  if (!quickCopy) {
    // copy conformations
    for (const auto &conf : other.d_confs) {
      if (confId < 0 || rdcast<int>(conf->getId()) == confId) {
        this->addConformer(new Conformer(*conf));
      }
    }

    // Copy sgroups
    for (const auto &sg : getSubstanceGroups(other)) {
      addSubstanceGroup(*this, sg);
    }

    d_props = other.d_props;

    // Bookmarks should be copied as well:
    for (auto abmI : other.d_atomBookmarks) {
      for (const auto *aptr : abmI.second) {
        setAtomBookmark(getAtomWithIdx(aptr->getIdx()), abmI.first);
      }
    }
    for (auto bbmI : other.d_bondBookmarks) {
      for (const auto *bptr : bbmI.second) {
        setBondBookmark(getBondWithIdx(bptr->getIdx()), bbmI.first);
      }
    }
  } else {
    d_props.reset();
    STR_VECT computed;
    d_props.setVal(RDKit::detail::computedPropName, computed);
  }

  // std::cerr<<"---------    done init from other: "<<this<<"
  // "<<&other<<std::endl;
}

void ROMol::initMol() {
  d_props.reset();
  dp_ringInfo = new RingInfo();
  // ok every molecule contains a property entry called
  // RDKit::detail::computedPropName
  // which provides
  //  list of property keys that correspond to value that have been computed
  // this can used to blow out all computed properties while leaving the rest
  // along
  // initialize this list to an empty vector of strings
  STR_VECT computed;
  d_props.setVal(RDKit::detail::computedPropName, computed);
}

unsigned int ROMol::getAtomDegree(const Atom *at) const {
  PRECONDITION(at, "no atom");
  PRECONDITION(&at->getOwningMol() == this,
               "atom not associated with this molecule");
    return rdcast<unsigned int>(at->bonds().size());
};

unsigned int ROMol::getNumAtoms(bool onlyExplicit) const {
    unsigned int res = static_cast<unsigned int>(_atoms.size());
  if (!onlyExplicit) {
    // if we are interested in hydrogens as well add them up from
    // each
    for (const auto atom : atoms()) {
      res += atom->getTotalNumHs();
    }
  }
  return res;
};
unsigned int ROMol::getNumHeavyAtoms() const {
  unsigned int res = 0;
  for (const auto atom : atoms()) {
    if (atom->getAtomicNum() > 1) {
      ++res;
    }
  }
  return res;
};

// returns the first inserted atom with the given bookmark
Atom *ROMol::getAtomWithBookmark(int mark) {
  auto lu = d_atomBookmarks.find(mark);
  PRECONDITION((lu != d_atomBookmarks.end() && !lu->second.empty()),
               "atom bookmark not found");
  return lu->second.front();
};

// returns all atoms with the given bookmark
ROMol::ATOM_PTR_LIST &ROMol::getAllAtomsWithBookmark(int mark) {
  auto lu = d_atomBookmarks.find(mark);
  PRECONDITION(lu != d_atomBookmarks.end(), "atom bookmark not found");
  return lu->second;
};

// returns the unique atom with the given bookmark
Atom *ROMol::getUniqueAtomWithBookmark(int mark) {
  auto lu = d_atomBookmarks.find(mark);
  PRECONDITION((lu != d_atomBookmarks.end()), "bookmark not found");
  return lu->second.front();
}

// returns the first inserted bond with the given bookmark
Bond *ROMol::getBondWithBookmark(int mark) {
  auto lu = d_bondBookmarks.find(mark);
  PRECONDITION((lu != d_bondBookmarks.end() && !lu->second.empty()),
               "bond bookmark not found");
  return lu->second.front();
};

// returns all bonds with the given bookmark
ROMol::BOND_PTR_LIST &ROMol::getAllBondsWithBookmark(int mark) {
  auto lu = d_bondBookmarks.find(mark);
  PRECONDITION(lu != d_bondBookmarks.end(), "bond bookmark not found");
  return lu->second;
};

// returns the unique bond with the given bookmark
Bond *ROMol::getUniqueBondWithBookmark(int mark) {
  auto lu = d_bondBookmarks.find(mark);
  PRECONDITION((lu != d_bondBookmarks.end()), "bookmark not found");
  return lu->second.front();
}

void ROMol::clearAtomBookmark(const int mark) { d_atomBookmarks.erase(mark); }

void ROMol::clearAtomBookmark(int mark, const Atom *atom) {
  PRECONDITION(atom, "no atom");
  auto lu = d_atomBookmarks.find(mark);
  if (lu != d_atomBookmarks.end()) {
    auto &marks = lu->second;
    unsigned int tgtIdx = atom->getIdx();
    auto entry = std::find_if(marks.begin(), marks.end(), [&tgtIdx](auto ptr) {
      return ptr->getIdx() == tgtIdx;
    });
    if (entry != marks.end()) {
      marks.erase(entry);
    }
    if (marks.empty()) {
      d_atomBookmarks.erase(mark);
    }
  }
}

void ROMol::clearBondBookmark(int mark) { d_bondBookmarks.erase(mark); }
void ROMol::clearBondBookmark(int mark, const Bond *bond) {
  PRECONDITION(bond, "no bond");
  auto lu = d_bondBookmarks.find(mark);
  if (lu != d_bondBookmarks.end()) {
    auto &marks = lu->second;
    unsigned int tgtIdx = bond->getIdx();
    auto entry = std::find_if(marks.begin(), marks.end(), [&tgtIdx](auto ptr) {
      return ptr->getIdx() == tgtIdx;
    });
    if (entry != marks.end()) {
      marks.erase(entry);
    }
    if (marks.empty()) {
      d_bondBookmarks.erase(mark);
    }
  }
}

unsigned int ROMol::getNumBonds(bool onlyHeavy) const {
  // By default return the bonds that connect only the heavy atoms
  // hydrogen connecting bonds are ignores
  auto res = numBonds;
  if (!onlyHeavy) {
    // If we need hydrogen connecting bonds add them up
    for (const auto atom : atoms()) {
      res += atom->getTotalNumHs();
    }
  }
  return res;
}

Bond *ROMol::getBondBetweenAtoms(unsigned int idx1, unsigned int idx2) {
  return const_cast<Bond *>(
      static_cast<const ROMol *>(this)->getBondBetweenAtoms(
          idx1, idx2));  // avoid code duplication
}

const Bond *ROMol::getBondBetweenAtoms(unsigned int idx1,
                                       unsigned int idx2) const {
  URANGE_CHECK(idx1, getNumAtoms());
  URANGE_CHECK(idx2, getNumAtoms());
  
  return _atoms[idx1]->getBondTo(idx2);
}

ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
  PRECONDITION(at, "no atom");
  PRECONDITION(&at->getOwningMol() == this,
               "atom not associated with this molecule");
    ADJ_ITER begin = at->nbrs().begin();
    ADJ_ITER end = at->nbrs().end();
    return std::make_pair(begin, end);
};

ROMol::OBOND_ITER_PAIR ROMol::getAtomBonds(Atom const *at) const {
  PRECONDITION(at, "no atom");
  PRECONDITION(&at->getOwningMol() == this,
               "atom not associated with this molecule");
    return std::make_pair(at->bonds().begin(), at->bonds().end());
}

ROMol::ATOM_ITER_PAIR ROMol::getVertices() {
    return std::make_pair(_atoms.begin(), _atoms.end()); }
ROMol::BOND_ITER_PAIR ROMol::getEdges() {         
    EDGE_ITER begin = _bonds.begin();
    EDGE_ITER  end = _bonds.end();
    return std::make_pair(begin, end);
 }
ROMol::ATOM_ITER_PAIR ROMol::getVertices() const {
    return std::make_pair(_atoms.begin(), _atoms.end());
}
ROMol::CONST_BOND_ITER_PAIR ROMol::getEdges() const {         
    CONST_EDGE_ITER begin = _bonds.begin();
    CONST_EDGE_ITER  end = _bonds.end();
    return std::make_pair(begin, end);
 }

unsigned int ROMol::addAtom(Atom *atom_pin, bool updateLabel,
                            bool takeOwnership) {
  PRECONDITION(atom_pin, "null atom passed in");
  PRECONDITION(!takeOwnership || !atom_pin->hasOwningMol() ||
                   &atom_pin->getOwningMol() == this,
               "cannot take ownership of an atom which already has an owner");
  Atom *atom_p;
  if (!takeOwnership) {
    atom_p = atom_pin->copy();
  } else {
    atom_p = atom_pin;
  }

  atom_p->setOwningMol(this);
    auto which = static_cast<unsigned int>(_atoms.size());
    _atoms.push_back(atom_p);
  atom_p->setIdx(which);
  if (updateLabel) {
    replaceAtomBookmark(atom_p, ci_RIGHTMOST_ATOM);
  }
  for (auto &conf : d_confs) {
    conf->setAtomPos(which, RDGeom::Point3D(0.0, 0.0, 0.0));
  }
  return rdcast<unsigned int>(which);
};

unsigned int ROMol::addBond(Bond *bond_pin, bool takeOwnership) {
  PRECONDITION(bond_pin, "null bond passed in");
  PRECONDITION(!takeOwnership || !bond_pin->hasOwningMol() ||
                   &bond_pin->getOwningMol() == this,
               "cannot take ownership of an bond which already has an owner");
  URANGE_CHECK(bond_pin->getBeginAtomIdx(), getNumAtoms());
  URANGE_CHECK(bond_pin->getEndAtomIdx(), getNumAtoms());
  PRECONDITION(bond_pin->getBeginAtomIdx() != bond_pin->getEndAtomIdx(),
               "attempt to add self-bond");
  PRECONDITION(getBondBetweenAtoms(bond_pin->getBeginAtomIdx(),
                                   bond_pin->getEndAtomIdx())==nullptr,
               "bond already exists");

  Bond *bond_p;
  if (!takeOwnership) {
    bond_p = bond_pin->copy();
  } else {
    bond_p = bond_pin;
  }
    auto a1 = _atoms[bond_p->getBeginAtomIdx()];
    auto a2 = _atoms[bond_p->getEndAtomIdx()];
    assert(a1->getIdx() != a2->getIdx());
    a1->_bonds.push_back(bond_p);
    a2->_bonds.push_back(bond_p);
    assert(a1->getIdx() != a2->getIdx());
    a1->_oatoms.push_back(a2);
    a2->_oatoms.push_back(a1);
    _bonds.push_back(bond_p);
  bond_p->setOwningMol(this);
  bond_p->setIdx(numBonds);
  numBonds++;
  return numBonds;
}

void ROMol::setStereoGroups(std::vector<StereoGroup> stereo_groups) {
  d_stereo_groups = std::move(stereo_groups);
}

void ROMol::debugMol(std::ostream &str) const {
  str << "Atoms:" << std::endl;
  for (const auto atom : atoms()) {
    str << "\t" << *atom << std::endl;
  }

  str << "Bonds:" << std::endl;
  for (const auto bond : bonds()) {
    str << "\t" << *bond << std::endl;
  }

  const auto &sgs = getSubstanceGroups(*this);
  if (!sgs.empty()) {
    str << "Substance Groups:" << std::endl;
    for (const auto &sg : sgs) {
      str << "\t" << sg << std::endl;
    }
  }

  const auto &stgs = getStereoGroups();
  if (!stgs.empty()) {
    unsigned idx = 0;
    str << "Stereo Groups:" << std::endl;
    for (const auto &stg : stgs) {
      str << "\t" << idx << ' ' << stg << std::endl;
      ++idx;
    }
  }
}

// --------------------------------------------
//
//  Iterators
//
// --------------------------------------------
ROMol::AtomIterator ROMol::beginAtoms() { return AtomIterator(this); }
ROMol::ConstAtomIterator ROMol::beginAtoms() const {
  return ConstAtomIterator(this);
}
ROMol::AtomIterator ROMol::endAtoms() {
  return AtomIterator(this, getNumAtoms());
}
ROMol::ConstAtomIterator ROMol::endAtoms() const {
  return ConstAtomIterator(this, getNumAtoms());
}

ROMol::AromaticAtomIterator ROMol::beginAromaticAtoms() {
  return AromaticAtomIterator(this);
}
ROMol::ConstAromaticAtomIterator ROMol::beginAromaticAtoms() const {
  return ConstAromaticAtomIterator(this);
}
ROMol::AromaticAtomIterator ROMol::endAromaticAtoms() {
  return AromaticAtomIterator(this, getNumAtoms());
}
ROMol::ConstAromaticAtomIterator ROMol::endAromaticAtoms() const {
  return ConstAromaticAtomIterator(this, getNumAtoms());
}

ROMol::HeteroatomIterator ROMol::beginHeteros() {
  return HeteroatomIterator(this);
}
ROMol::ConstHeteroatomIterator ROMol::beginHeteros() const {
  return ConstHeteroatomIterator(this);
}
ROMol::HeteroatomIterator ROMol::endHeteros() {
  return HeteroatomIterator(this, getNumAtoms());
}
ROMol::ConstHeteroatomIterator ROMol::endHeteros() const {
  return ConstHeteroatomIterator(this, getNumAtoms());
}

bool ROMol::hasQuery() const {
  for (auto atom : atoms()) {
    if (atom->hasQuery()) {
      return true;
    }
  }
  for (auto bond : bonds()) {
    if (bond->hasQuery()) {
      return true;
    }
  }
  return false;
}

ROMol::QueryAtomIterator ROMol::beginQueryAtoms(QueryAtom const *what) {
  return QueryAtomIterator(this, what);
}
ROMol::ConstQueryAtomIterator ROMol::beginQueryAtoms(
    QueryAtom const *what) const {
  return ConstQueryAtomIterator(this, what);
}
ROMol::QueryAtomIterator ROMol::endQueryAtoms() {
  return QueryAtomIterator(this, getNumAtoms());
}
ROMol::ConstQueryAtomIterator ROMol::endQueryAtoms() const {
  return ConstQueryAtomIterator(this, getNumAtoms());
}
ROMol::MatchingAtomIterator ROMol::beginMatchingAtoms(bool (*what)(Atom *)) {
  return MatchingAtomIterator(this, what);
}
ROMol::ConstMatchingAtomIterator ROMol::beginMatchingAtoms(
    bool (*what)(const Atom *)) const {
  return ConstMatchingAtomIterator(this, what);
}
ROMol::MatchingAtomIterator ROMol::endMatchingAtoms() {
  return MatchingAtomIterator(this, getNumAtoms());
}
ROMol::ConstMatchingAtomIterator ROMol::endMatchingAtoms() const {
  return ConstMatchingAtomIterator(this, getNumAtoms());
}

ROMol::BondIterator ROMol::beginBonds() { return BondIterator(this); }
ROMol::ConstBondIterator ROMol::beginBonds() const {
  return ConstBondIterator(this);
}
ROMol::BondIterator ROMol::endBonds() {
  auto [beg, end] = getEdges();
  return BondIterator(this, end);
}
ROMol::ConstBondIterator ROMol::endBonds() const {
  auto [beg, end] = getEdges();
  return ConstBondIterator(this, end);
}

void ROMol::clearComputedProps(bool includeRings) const {
  // the SSSR information:
  if (includeRings) {
    this->dp_ringInfo->reset();
  }

  RDProps::clearComputedProps();

  for (auto atom : atoms()) {
    atom->clearComputedProps();
  }

  for (auto bond : bonds()) {
    bond->clearComputedProps();
  }
}

void ROMol::updatePropertyCache(bool strict) {
  for (auto atom : atoms()) {
    atom->updatePropertyCache(strict);
  }
  for (auto bond : bonds()) {
    bond->updatePropertyCache(strict);
  }
}

bool ROMol::needsUpdatePropertyCache() const {
  for (const auto atom : atoms()) {
    if (atom->needsUpdatePropertyCache()) {
      return true;
    }
  }
  // there is no test for bonds yet since they do not obtain a valence property
  return false;
}

const Conformer &ROMol::getConformer(int id) const {
  // make sure we have more than one conformation
  if (d_confs.size() == 0) {
    throw ConformerException("No conformations available on the molecule");
  }

  if (id < 0) {
    return *(d_confs.front());
  }
  auto cid = (unsigned int)id;
  for (auto conf : d_confs) {
    if (conf->getId() == cid) {
      return *conf;
    }
  }
  // we did not find a conformation with the specified ID
  std::string mesg = "Can't find conformation with ID: ";
  mesg += id;
  throw ConformerException(mesg);
}

Conformer &ROMol::getConformer(int id) {
  return const_cast<Conformer &>(
      static_cast<const ROMol *>(this)->getConformer(id));
}

void ROMol::removeConformer(unsigned int id) {
  for (auto ci = d_confs.begin(); ci != d_confs.end(); ++ci) {
    if ((*ci)->getId() == id) {
      d_confs.erase(ci);
      return;
    }
  }
}

unsigned int ROMol::addConformer(Conformer *conf, bool assignId) {
  PRECONDITION(conf, "bad conformer");
  PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
               "Number of atom mismatch");
  if (assignId) {
    int maxId = -1;
    for (auto cptr : d_confs) {
      maxId = std::max((int)(cptr->getId()), maxId);
    }
    maxId++;
    conf->setId((unsigned int)maxId);
  }
  conf->setOwningMol(this);
  CONFORMER_SPTR nConf(conf);
  d_confs.push_back(nConf);
  return conf->getId();
}

#ifdef RDK_USE_BOOST_SERIALIZATION
template <class Archive>
void ROMol::save(Archive &ar, const unsigned int) const {
  std::string pkl;
  MolPickler::pickleMol(*this, pkl, PicklerOps::AllProps);
  ar << pkl;
}

template <class Archive>
void ROMol::load(Archive &ar, const unsigned int) {
  std::string pkl;
  ar >> pkl;

  delete dp_ringInfo;
  initMol();

  numBonds = 0;
  MolPickler::molFromPickle(pkl, *this, PicklerOps::AllProps);
    numBonds = static_cast<unsigned int>(_bonds.size());
}

template RDKIT_GRAPHMOL_EXPORT void ROMol::save<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, const unsigned int) const;
template RDKIT_GRAPHMOL_EXPORT void ROMol::load<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, const unsigned int);
#endif

}  // namespace RDKit
