//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Search/RotorDriver.h"
#include "FragmentConfGen.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <utility>

namespace RDKit {

namespace {
// First neighbor of atom `a` other than atom `avoid`, prefers a heavy atom so the
// driven dihedral is defined by real substituents where possible.  -1 if none.
int pickNbr(const ROMol &mol, unsigned int a, unsigned int avoid) {
  int fallback = -1;
  for (const auto nbr : mol.atomNeighbors(mol.getAtomWithIdx(a))) {
    unsigned int idx = nbr->getIdx();
    if (idx == avoid) {
      continue;
    }
    if (nbr->getAtomicNum() > 1) {
      return static_cast<int>(idx);
    }
    if (fallback < 0) {
      fallback = static_cast<int>(idx);
    }
  }
  return fallback;
}

// Find the downstream tree that we will be rotating
std::vector<unsigned int> sideAtoms(const ROMol &mol, unsigned int start,
                                    unsigned int blocked) {
  std::vector<char> seen(mol.getNumAtoms(), 0);
  seen[blocked] = 1;
  seen[start] = 1;
  std::vector<unsigned int> stack{start};
  std::vector<unsigned int> out;
  while (!stack.empty()) {
    unsigned int a = stack.back();
    stack.pop_back();
    out.push_back(a);
    for (const auto nbr : mol.atomNeighbors(mol.getAtomWithIdx(a))) {
      unsigned int idx = nbr->getIdx();
      if (!seen[idx]) {
        seen[idx] = 1;
        stack.push_back(idx);
      }
    }
  }
  return out;
}

// Pick the root of the rotor hierarchy at the largest rigid fragment: cut the rotor bonds,
//  the thinking is this frag is fixed saving some computations
unsigned int chooseRoot(const ROMol &mol,
                        const std::vector<unsigned int> &rotorBonds) {
  unsigned int n = mol.getNumAtoms();
  std::vector<unsigned int> parent(n);
  for (unsigned int a = 0; a < n; ++a) {
    parent[a] = a;
  }
  std::function<unsigned int(unsigned int)> find =
      [&](unsigned int x) -> unsigned int {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  std::vector<char> isRotor(mol.getNumBonds(), 0);
  for (unsigned int bi : rotorBonds) {
    isRotor[bi] = 1;
  }
  for (const auto b : mol.bonds()) {
    if (!isRotor[b->getIdx()]) {
      unsigned int ra = find(b->getBeginAtomIdx());
      unsigned int rb = find(b->getEndAtomIdx());
      if (ra != rb) {
        parent[ra] = rb;
      }
    }
  }
  std::vector<unsigned int> fragSize(n, 0);
  for (unsigned int a = 0; a < n; ++a) {
    ++fragSize[find(a)];
  }
  unsigned int bestRep = 0;
  unsigned int bestSize = 0;
  for (unsigned int a = 0; a < n; ++a) {
    unsigned int rep = find(a);
    if (fragSize[rep] > bestSize) {
      bestSize = fragSize[rep];
      bestRep = rep;
    }
  }
  for (unsigned int a = 0; a < n; ++a) {
    if (find(a) == bestRep) {
      return a;  // lowest-index atom of the largest fragment
    }
  }
  return 0;
}
}  // namespace

RotorDriver::RotorDriver(const ROMol &mol, int confId, ScoreFn fullFF,
                         ScoreFn fragmentFF)
    : d_fullFF(std::move(fullFF)), d_fragFF(std::move(fragmentFF)) {
  init(mol, FragmentConfGen::findLinkBonds(mol), confId);
}

RotorDriver::RotorDriver(const ROMol &mol,
                         const std::vector<unsigned int> &rotorBonds,
                         int confId, ScoreFn fullFF, ScoreFn fragmentFF)
    : d_fullFF(std::move(fullFF)), d_fragFF(std::move(fragmentFF)) {
  init(mol, rotorBonds, confId);
}

void RotorDriver::init(const ROMol &mol,
                       const std::vector<unsigned int> &rotorBonds,
                       int confId) {
  d_nAtoms = mol.getNumAtoms();
  const Conformer &conf = mol.getConformer(confId);
  d_pos.resize(3 * static_cast<size_t>(d_nAtoms));
  for (unsigned int a = 0; a < d_nAtoms; ++a) {
    const auto &p = conf.getAtomPos(a);
    d_pos[3 * a] = p.x;
    d_pos[3 * a + 1] = p.y;
    d_pos[3 * a + 2] = p.z;
  }

  d_rootAtom = chooseRoot(mol, rotorBonds);

  d_rotors.reserve(rotorBonds.size());
  for (unsigned int bi : rotorBonds) {
    const Bond *b = mol.getBondWithIdx(bi);
    unsigned int a = b->getBeginAtomIdx();
    unsigned int c = b->getEndAtomIdx();
    // In a j-k junction, orient the rotor so k is the end AWAY from the root.
    //  everything on side k will be rotated
    std::vector<unsigned int> aSide = sideAtoms(mol, a, c);
    bool rootOnA =
        std::find(aSide.begin(), aSide.end(), d_rootAtom) != aSide.end();
    unsigned int j = rootOnA ? a : c;  // fixed (root) side
    unsigned int k = rootOnA ? c : a;  // moving side
    int i = pickNbr(mol, j, k);
    int l = pickNbr(mol, k, j);
    if (i < 0 || l < 0) {
      // a terminal bond has no dihedral to drive; skip it
      continue;
    }
    Rotor rt;
    rt.i = static_cast<unsigned int>(i);
    rt.j = j;
    rt.k = k;
    rt.l = static_cast<unsigned int>(l);
    // moving = the k-side minus k itself (j-k) will remain fixed
    if (rootOnA) {
      rt.moving = sideAtoms(mol, c, a);
    } else {
      rt.moving = std::move(aSide);
    }
    rt.moving.erase(std::remove(rt.moving.begin(), rt.moving.end(), k),
                    rt.moving.end());
    d_rotors.push_back(std::move(rt));
  }

  // XXX Fix ME: Some searchers don't need the hierarchy, so maybe optionally build it
  buildRotorHierarchy();
}

void RotorDriver::buildRotorHierarchy() {
  size_t nr = d_rotors.size();
  // per-rotor membership mask over atoms, for O(1) "is atom in r's moving set"
  std::vector<std::vector<char>> mask(nr, std::vector<char>(d_nAtoms, 0));
  for (size_t r = 0; r < nr; ++r) {
    for (unsigned int atom : d_rotors[r].moving) {
      mask[r][atom] = 1;
    }
  }
  // Rotor r2 is contained in r iff r2's away-anchor (r2.k) lies in r's moving set:
  // then all of r2 is inside r's subtree
  d_contained.assign(nr, {});
  for (size_t r = 0; r < nr; ++r) {
    for (size_t r2 = 0; r2 < nr; ++r2) {
      if (r2 != r && mask[r][d_rotors[r2].k]) {
        d_contained[r].push_back(static_cast<unsigned int>(r2));
      }
    }
  }
  auto byScaleDesc = [this](unsigned int x, unsigned int y) {
    size_t mx = d_rotors[x].moving.size();
    size_t my = d_rotors[y].moving.size();
    if (mx != my) {
      return mx > my;  // largest movement first
    }
    return x < y;  // stable, deterministic tie-break
  };
  for (auto &lst : d_contained) {
    std::sort(lst.begin(), lst.end(), byScaleDesc);
  }
  d_scaleOrder.resize(nr);
  for (size_t r = 0; r < nr; ++r) {
    d_scaleOrder[r] = static_cast<unsigned int>(r);
  }
  std::sort(d_scaleOrder.begin(), d_scaleOrder.end(), byScaleDesc);
}

RDGeom::Point3D RotorDriver::pos(unsigned int a) const {
  return RDGeom::Point3D(d_pos[3 * a], d_pos[3 * a + 1], d_pos[3 * a + 2]);
}

void RotorDriver::setPos(unsigned int a, const RDGeom::Point3D &p) {
  d_pos[3 * a] = p.x;
  d_pos[3 * a + 1] = p.y;
  d_pos[3 * a + 2] = p.z;
}

std::array<unsigned int, 4> RotorDriver::torsion(unsigned int r) const {
  const Rotor &rt = d_rotors[r];
  return {rt.i, rt.j, rt.k, rt.l};
}

const std::vector<unsigned int> &RotorDriver::movingAtoms(
    unsigned int r) const {
  return d_rotors[r].moving;
}

const std::vector<unsigned int> &RotorDriver::rotorSubsets(
    unsigned int r) const {
  return d_contained[r];
}

double RotorDriver::dihedralDeg(unsigned int r) const {
  const Rotor &rt = d_rotors[r];
  RDGeom::Point3D p0 = pos(rt.i);
  RDGeom::Point3D p1 = pos(rt.j);
  RDGeom::Point3D p2 = pos(rt.k);
  RDGeom::Point3D p3 = pos(rt.l);
  RDGeom::Point3D b0 = p0 - p1;
  RDGeom::Point3D b1 = p2 - p1;
  RDGeom::Point3D b2 = p3 - p2;
  b1.normalize();
  // components of b0/b2 perpendicular to the axis b1
  RDGeom::Point3D v = b0 - b1 * b0.dotProduct(b1);
  RDGeom::Point3D w = b2 - b1 * b2.dotProduct(b1);
  double x = v.dotProduct(w);
  double y = b1.crossProduct(v).dotProduct(w);
  return atan2(y, x) * 180.0 / M_PI;
}

void RotorDriver::setDihedral(unsigned int r, double angleDeg) {
  const Rotor &rt = d_rotors[r];
  double delta = (angleDeg - dihedralDeg(r)) * M_PI / 180.0;
  RDGeom::Point3D pj = pos(rt.j);
  RDGeom::Point3D pk = pos(rt.k);
  RDGeom::Point3D axis = pk - pj;
  double alen = axis.length();
  if (alen < 1e-9) {
    return;
  }
  axis /= alen;
  double c = std::cos(delta);
  double s = std::sin(delta);

  for (unsigned int a : rt.moving) {
    RDGeom::Point3D p = pos(a) - pk;
    RDGeom::Point3D cross = axis.crossProduct(p);
    double dot = axis.dotProduct(p);
    RDGeom::Point3D pr = p * c + cross * s + axis * (dot * (1.0 - c));
    setPos(a, pr + pk);
  }
}

double RotorDriver::drive(unsigned int r, double angleDeg) {
  setDihedral(r, angleDeg);
  return score();
}

double RotorDriver::scoreFull() const {
  return d_fullFF ? d_fullFF(d_pos.data(), d_nAtoms)
                  : std::numeric_limits<double>::quiet_NaN();
}

double RotorDriver::scoreFragment() const {
  return d_fragFF ? d_fragFF(d_pos.data(), d_nAtoms)
                  : std::numeric_limits<double>::quiet_NaN();
}

double RotorDriver::score() const {
  if (d_fragFF) {
    return d_fragFF(d_pos.data(), d_nAtoms);
  }
  if (d_fullFF) {
    return d_fullFF(d_pos.data(), d_nAtoms);
  }
  return std::numeric_limits<double>::quiet_NaN();
}


}  // namespace RDKit
