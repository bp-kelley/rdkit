//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Embedder/Fraglib.h"
#include "Search/InterFragScore.h"

#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/Nonbonded.h>
#include <ForceField/MMFF/TorsionAngle.h>
#include <ForceField/MMFF/Params.h>

#include "Utils/NonbondedLookup.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <vector>

namespace RDKit {

namespace {
//! Build a scoring function BETWEEN rotatable fragments
/*!
  Note: only MMFF is supported, any other force field will just
  return a full scorer for now.

  XXX FIX ME -> not a full abstraction yet, several searchers
    build their own FF for convenience and only workk with MMFF
*/
enum class FFFamily { MMFF, Unsupported };

FFFamily ffFamily(const std::string &ffVariant) {
  if (isValidFF(ffVariant)) {
    return FFFamily::MMFF;
  }
  return FFFamily::Unsupported;
}
}  // namespace

IncrementalInterFragScore::IncrementalInterFragScore(
    const RotorDriver &driver, const InterFragScoreHandles &handles)
    : d_vdw(handles.vdw), d_tor(handles.tor) {
  // XXX FIX ME: The lookup variant can't currently be packed so fall back to full
  d_valid = (d_vdw != nullptr) && !d_vdw->usesLookup();
  if (!d_valid) {
    return;
  }
  d_cut2 = d_vdw->cut2();
  const size_t nr = driver.numRotors();
  const size_t nPairs = d_vdw->numPairs();
  d_packed.resize(nr);
  std::vector<char> moved(driver.numAtoms(), 0);
  for (size_t r = 0; r < nr; ++r) {
    const auto &mv = driver.movingAtoms(static_cast<unsigned int>(r));
    for (unsigned int at : mv) {
      moved[at] = 1;
    }
    auto &pr = d_packed[r];
    for (size_t p = 0; p < nPairs; ++p) {
      // XOR: exactly one endpoint moves -> the pair distance changes.
      // Both-moved pairs rotate rigidly (distance fixed) and MUST be excluded.
      if (moved[d_vdw->at1(p)] != moved[d_vdw->at2(p)]) {
        pr.a1.push_back(d_vdw->at1(p));
        pr.a2.push_back(d_vdw->at2(p));
        pr.Rs.push_back(d_vdw->rStar(p));
        pr.wd.push_back(d_vdw->wellDepth(p));
      }
    }
    for (unsigned int at : mv) {
      moved[at] = 0;  // reset for the next rotor
    }
  }
}

double IncrementalInterFragScore::changingVdw(const double *pos,
                                              unsigned int r) const {
  if (!d_valid || r >= d_packed.size()) {
    return 0.0;
  }
  const auto &pr = d_packed[r];
  return ForceFields::MMFF::InterFragVdWContrib::energyOfPacked(
      const_cast<double *>(pos), pr.a1.size(), pr.a1.data(), pr.a2.data(),
      pr.Rs.data(), pr.wd.data(), d_cut2);
}

double IncrementalInterFragScore::torsionEnergy(const double *pos) const {
  return d_tor ? d_tor->getEnergy(const_cast<double *>(pos)) : 0.0;
}

size_t IncrementalInterFragScore::numPairs() const {
  return d_vdw ? d_vdw->numPairs() : 0;
}

std::vector<std::vector<JunctionTorsionTerm>> junctionTorsionTerms(
    const ROMol &mol, MMFF::MMFFMolProperties &props,
    const std::vector<std::pair<unsigned int, unsigned int>> &junctionBonds) {
  std::vector<std::vector<JunctionTorsionTerm>> out(junctionBonds.size());
  if (!props.isValid()) {
    return out;
  }
  // the junction bonds may not actually be MMFF torsion bonds, so we still
  // need to check against the MMFF rules
  std::vector<MatchVectType> matchVect;
  if (const ROMol *query = MMFF::Tools::DefaultTorsionBondSmarts::query()) {
    SubstructMatch(mol, *query, matchVect);
  }
  std::set<std::pair<unsigned int, unsigned int>> torBonds;
  for (const auto &m : matchVect) {
    const unsigned int a = m[0].second, b = m[1].second;
    torBonds.insert({std::min(a, b), std::max(a, b)});
  }
  auto isSp23 = [](const Atom *a) {
    return a->getHybridization() == Atom::SP2 ||
           a->getHybridization() == Atom::SP3;
  };

  for (size_t bi = 0; bi < junctionBonds.size(); ++bi) {
    const unsigned int j = junctionBonds[bi].first;
    const unsigned int k = junctionBonds[bi].second;
    if (!torBonds.count({std::min(j, k), std::max(j, k)})) {
      continue;
    }
    const Atom *jA = mol.getAtomWithIdx(j);
    const Atom *kA = mol.getAtomWithIdx(k);
    if (!isSp23(jA) || !isSp23(kA)) {
      continue;
    }
    for (const auto n1 : mol.atomNeighbors(jA)) {
      const unsigned int a1 = n1->getIdx();
      if (a1 == k) {
        continue;
      }
      for (const auto n2 : mol.atomNeighbors(kA)) {
        const unsigned int a4 = n2->getIdx();
        if (a4 == j || a4 == a1) {
          continue;
        }
        unsigned int torType;
        ForceFields::MMFF::MMFFTor params;
        if (props.getMMFFTorsionParams(mol, a1, j, k, a4, torType, params)) {
          out[bi].push_back({a1, j, k, a4, params.V1, params.V2, params.V3});
        }
      }
    }
  }
  return out;
}

RotorDriver::ScoreFn makeInterFragmentScoreFn(
    const ROMol &mol, const std::vector<int> &atomFragments,
    const std::vector<std::pair<unsigned int, unsigned int>> &junctionBonds,
    bool electrostatics, const std::string &ffVariant, double vdwCutoff,
    InterFragScoreHandles *handles) {
  if (handles) {
    *handles = InterFragScoreHandles{};
  }
  if (ffFamily(ffVariant) != FFFamily::MMFF) {
    BOOST_LOG(rdWarningLog)
        << "makeInterFragmentScoreFn: unsupported force field \"" << ffVariant
        << "\"; returning empty scoring function";
    return {};
  }
  auto scoreFxnMol = boost::make_shared<ROMol>(mol);
  if (scoreFxnMol->getNumConformers() == 0) {
    const bool ownsConf = true;
    scoreFxnMol->addConformer(new Conformer(scoreFxnMol->getNumAtoms()), ownsConf);
  }
  MMFF::MMFFMolProperties props(*scoreFxnMol, ffVariant);
  if (!props.isValid()) {
    return {};
  }
  props.setMMFFEleTerm(electrostatics);


  const unsigned int nAtoms = scoreFxnMol->getNumAtoms();
  Conformer &conf = scoreFxnMol->getConformer();
  auto ff = boost::make_shared<ForceFields::ForceField>();
  for (unsigned int i = 0; i < nAtoms; ++i) {
    ff->positions().push_back(&conf.getAtomPos(i));
  }

  // build the neighbor matrix we need in the contribs
  auto nbrMat = MMFF::Tools::buildNeighborMatrix(*scoreFxnMol);
  auto vdw = std::make_unique<ForceFields::MMFF::InterFragVdWContrib>(
      ff.get(), vdwCutoff, /*useLookup=*/false);
  bool anyVdw = false;
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = i + 1; j < nAtoms; ++j) {
      if (atomFragments[i] == atomFragments[j]) {
        continue;
      }
      if (MMFF::Tools::getTwoBitCell(
              nbrMat, MMFF::Tools::twoBitCellPos(nAtoms, i, j)) >=
          MMFF::Tools::RELATION_1_4) {
        ForceFields::MMFF::MMFFVdWRijstarEps c;
        if (props.getMMFFVdWParams(i, j, c)) {
          vdw->addTerm(i, j, &c);
          anyVdw = true;
        }
      }
    }
  }
  if (anyVdw) {
    if (handles) {
      handles->vdw = vdw.get();  // stays valid: ff owns it
    }
    ff->contribs().push_back(ForceFields::ContribPtr(vdw.release()));
  }

  auto tor = std::make_unique<ForceFields::MMFF::TorsionAngleContrib>(ff.get());
  bool hasTor = false;
  for (const auto &perBond :
       junctionTorsionTerms(*scoreFxnMol, props, junctionBonds)) {
    for (const auto &t : perBond) {
      ForceFields::MMFF::MMFFTor params;
      params.V1 = t.V1;
      params.V2 = t.V2;
      params.V3 = t.V3;
      tor->addTerm(t.a1, t.j, t.k, t.a4, &params);
      hasTor = true;
    }
  }


  if (hasTor) {
    if (handles) {
      handles->tor = tor.get();
    }
    ff->contribs().push_back(ForceFields::ContribPtr(tor.release()));
  }

  ff->initialize();

  // Return the dynmically created ScoreFxn.
  //  lambda are neat.
  // note: the closure here owns scoreFxnMol and ff
  return [scoreFxnMol, ff](const double *pos, unsigned int) -> double {
    // We make our own function here to avoid recreating the distance cache
    //  RDKit's ForceField::calcEnergy() is a bit too generic for this purpose
    //  so we inline it here.
    // Also, we have a sentinel, MMFF calcEnergy can throw for degenerate placements
    //  (identical coords) we catch this and return an absurd energy.
    try {
      double *p = const_cast<double *>(pos);
      double e = 0.0;
      for (const auto &c : ff->contribs()) {
        e += c->getEnergy(p);
      }
      return e;
    } catch (...) {
      return 1.0e12;
    }
  };
}

RotorDriver::ScoreFn makeFullFFScoreFn(const ROMol &mol, bool electrostatics,
                                       const std::string &ffVariant,
                                       double nonBondedThresh) {
  if (ffFamily(ffVariant) != FFFamily::MMFF) {
    BOOST_LOG(rdWarningLog)
        << "makeFullFFScoreFn: unsupported force field \"" << ffVariant
        << "\"; returning empty scoring function";
    return {};
  }
  auto scoreFxnMol = boost::make_shared<ROMol>(mol);
  if (scoreFxnMol->getNumConformers() == 0) {
    // We need a conf to satisfy the constructor, our internal call uses a double * buffer
    const bool ownsConf = true;
    scoreFxnMol->addConformer(new Conformer(scoreFxnMol->getNumAtoms()), ownsConf);
  }
  MMFF::MMFFMolProperties props(*scoreFxnMol, ffVariant);
  if (!props.isValid()) {
    BOOST_LOG(rdWarningLog) << "Could not prepare ForceField: " << ffVariant << " returning empty scoring function";
    return {};
  }
  props.setMMFFEleTerm(electrostatics);

  std::shared_ptr<ForceFields::ForceField> ff(
      MMFF::constructForceField(*scoreFxnMol, &props, nonBondedThresh));
  if (!ff) {
    return {};
  }
  ff->initialize();
  // return an absurd energy when mmff finds degenerate coords
  //  XXX FIX ME ->
  //  n.b. this is slower than it should be due to the calcEnergy recreating the
  //   nbr matrix, we could probably sum energies like we do above
  // note: the closure here owns scoreFxnMol and ff  
  return [scoreFxnMol, ff](const double *pos, unsigned int) -> double {
    try {
      return ff->calcEnergy(const_cast<double *>(pos));
    } catch (...) {
      return 1.0e12;
    }
  };
}

}  // namespace RDKit
