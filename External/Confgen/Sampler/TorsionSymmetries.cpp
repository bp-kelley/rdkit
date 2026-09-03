//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Sampler/TorsionSymmetries.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/RDLog.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace RDKit {

const std::string SYMMETRY_CLASSES = "__torsionSymmetryClasses";

namespace {

double normalizeAngle(double angle) {
  angle = std::fmod(angle, 360.0);
  if (angle < 0.0) {
    angle += 360.0;
  }
  return angle;
}

double angleDistance(double a, double b) {
  const double diff = std::fabs(normalizeAngle(a) - normalizeAngle(b));
  return std::min(diff, 360.0 - diff);
}

bool addUniqueAngle(std::vector<double> &angles, double angle,
                    double minSeparation) {
  angle = normalizeAngle(angle);
  for (const double existing : angles) {
    if (angleDistance(existing, angle) < minSeparation) {
      return false;
    }
  }
  angles.push_back(angle);
  return true;
}

bool isExplicitHydrogen(const Atom *atom) {
  return atom && atom->getAtomicNum() == 1;
}

bool hasLinearSubstituent(const ROMol &mol, unsigned int center,
                          unsigned int axisOther) {
  const Atom *atom = mol.getAtomWithIdx(center);
  if (atom->getHybridization() != Atom::SP) {
    return false;
  }
  for (const auto nbr : mol.atomNeighbors(atom)) {
    const unsigned int n = nbr->getIdx();
    if (n == axisOther) {
      continue;
    }
    const Bond *bond = mol.getBondBetweenAtoms(center, n);
    if (bond && bond->getBondType() == Bond::TRIPLE) {
      return true;
    }
  }
  return false;
}

//! Rotational-symmetry order of `center`'s substituents about the j-k axis.
//!
//!  Find symmetry via the topological ranks and then check geometry for
//!   rotational symmetry.
//!
unsigned int sideSymmetryFromRanks(const ROMol &mol,
                                   const std::vector<SymmetryClass> &ranks,
                                   unsigned int center,
                                   unsigned int axisOther) {
  const Atom *atom = mol.getAtomWithIdx(center);
  unsigned int nH = 0;
  std::map<SymmetryClass, unsigned int> heavy;
  for (const auto nbr : mol.atomNeighbors(atom)) {
    const unsigned int n = nbr->getIdx();
    if (n == axisOther) {
      continue;
    }
    if (isExplicitHydrogen(nbr)) {
      ++nH;
    } else {
      ++heavy[ranks[n]];
    }
  }
  // GEOMETRIC guard: topological (graph) equivalence is necessary but not
  // sufficient --
  //  the topological ranks can over fold, need to lok at sp3 and sp2 centers
  //  differently
  // ALL-HYDROGEN substituents.  Heavy-atom rotors (CF3, CCl3, C(CH3)3) are NOT
  // handled here -- they have no hydrogens to count and fall to the branch
  // below, which is what gives them their C3.
  if (heavy.empty()) {
    // three fold Hs or just one don't need rotation sampling.
    //  XXX FIX ME - I'm not really sure about the one h...
    return (nH == 3 && atom->getHybridization() == Atom::SP3) ? nH : 1;
  }
  // Ignore mirrored symmetry (NH2) and the like
  if (heavy.size() == 1 && nH == 0) {
    const unsigned int n = heavy.begin()->second;
    if (n >= 3 && atom->getHybridization() == Atom::SP3)
      return n;                                     // C(CH3)3 / CX3
    if (n == 2 && atom->getIsAromatic()) return 2;  // aromatic ring flip
    return 1;
  }
  return 1;  // no rotational redundancy
}

//! Rotational-symmetry fold rule:
//!  driven bond has atom mapped to :1 with a given symmetry redudancy
struct SymPattern {
  std::shared_ptr<ROMol> query;
  std::array<int, 2> torsionBondAtoms{
      {-1, -1}};  // query atom idxs carrying map number 1
  unsigned int symmetry = 1;
};

std::string foldRulesPath() {
  if (const char *e = std::getenv("FCG_FOLD_RULES")) {
    return e;
  }
  return "";
}

//! Parse the fold-rules file from a smarts like file
//! Format is:
//!  "<order> <SMARTS>  [# comment]".
//! If the file is missing, no rules load and
//! no folding is applied
const std::vector<SymPattern> &symmetryPatterns() {
  static const std::vector<SymPattern> pats = []() {
    std::vector<SymPattern> out;
    const auto &filePath = foldRulesPath();
    std::ifstream in(filePath);
    if (!in.good()) {
      BOOST_LOG(rdWarningLog)
          << std::string("Smarts Torsion Fold rules file '") << filePath
          << "' could not be read, check environment FCG_FOLD_RULES";
      return out;
    }
    RDLog::LogStateSetter blocker;
    std::string line;
    while (std::getline(in, line)) {
      std::istringstream ss(line);
      unsigned int order = 0;
      std::string smarts;
      // read the simple smarts patterns "smarts" n-fold
      if (!(ss >> order >> smarts) || order == 0) {
        continue;
      }
      std::shared_ptr<ROMol> q(SmartsToMol(smarts));
      if (!q) {
        continue;
      }
      SymPattern sp;
      sp.query = q;
      sp.symmetry = order;
      int n = 0;
      for (const auto atom : q->atoms()) {
        if (atom->getAtomMapNum() == 1 && n < 2) {
          sp.torsionBondAtoms[n++] = static_cast<int>(atom->getIdx());
        }
      }
      if (n == 2) out.push_back(std::move(sp));
    }
    return out;
  }();
  return pats;
}

}  // namespace

std::vector<SymmetryClass> getSymmetryClasses(const ROMol &mol) {
  // Since this will be called a lot, cache ranks on the molecule
  std::vector<SymmetryClass> ranks;
  if (mol.getPropIfPresent(SYMMETRY_CLASSES, ranks)) {
    return ranks;
  }
  Canon::rankMolAtoms(mol, ranks, /*breakTies=*/false);
  mol.setProp<std::vector<SymmetryClass>>(SYMMETRY_CLASSES, ranks);
  return ranks;
}

unsigned int torsionRotationalSymmetryByRanks(
    const ROMol &mol, const std::vector<SymmetryClass> &ranks, unsigned int j,
    unsigned int k) {
  if (j >= mol.getNumAtoms() || k >= mol.getNumAtoms() ||
      !mol.getBondBetweenAtoms(j, k)) {
    return 1;
  }
  if (hasLinearSubstituent(mol, j, k) || hasLinearSubstituent(mol, k, j)) {
    return 360;
  }
  const unsigned int symJ = sideSymmetryFromRanks(mol, ranks, j, k);
  const unsigned int symK = sideSymmetryFromRanks(mol, ranks, k, j);
  return std::max(symJ, symK);
}

unsigned int torsionRotationalSymmetryByRanks(const ROMol &mol, unsigned int j,
                                              unsigned int k) {
  return torsionRotationalSymmetryByRanks(mol, getSymmetryClasses(mol), j, k);
}

unsigned int torsionRotationalSymmetryBySmarts(const ROMol &mol, unsigned int j,
                                               unsigned int k) {
  // optimization: anchor the match to just the junction bond, no need for a
  // global search
  for (const auto &p : symmetryPatterns()) {
    if (p.torsionBondAtoms[0] < 0 || p.torsionBondAtoms[1] < 0) {
      continue;
    }
    const unsigned int qa = static_cast<unsigned int>(p.torsionBondAtoms[0]);
    const unsigned int qb = static_cast<unsigned int>(p.torsionBondAtoms[1]);
    const ROMol &query = *p.query;
    SubstructMatchParameters ps;
    ps.uniquify = false;
    ps.maxMatches = 8;
    ps.extraAtomCheck = [qa, qb, j, k, &query](const Atom &q, const Atom &m) {
      // XXX FIX ME - do we really need this check?
      if (&q.getOwningMol() != &query) {
        return true;
      }
      const unsigned int qi = q.getIdx();
      if (qi == qa || qi == qb) {
        const unsigned int mi = m.getIdx();
        return mi == j || mi == k;
      }
      return true;
    };

    for (const auto &match : SubstructMatch(mol, query, ps)) {
      int a = -1, b = -1;
      for (const auto &pr : match) {
        if (pr.first == p.torsionBondAtoms[0])
          a = pr.second;
        else if (pr.first == p.torsionBondAtoms[1])
          b = pr.second;
      }
      if ((a == static_cast<int>(j) && b == static_cast<int>(k)) ||
          (a == static_cast<int>(k) && b == static_cast<int>(j))) {
        return p.symmetry;
      }
    }
  }
  return 1;
}

std::vector<double> foldAnglesByTorsionSymmetry(
    unsigned int sym, const std::vector<double> &angles, double minSeparation) {
  if (sym <= 1) {
    return angles;
  }
  std::vector<double> folded;
  if (sym == 360) {
    addUniqueAngle(folded, 0.0, minSeparation);
    return folded;
  }
  const double period = 360.0 / static_cast<double>(sym);
  for (const double angle : angles) {
    addUniqueAngle(folded, std::fmod(normalizeAngle(angle), period),
                   minSeparation);
  }
  return folded;
}

}  // namespace RDKit
