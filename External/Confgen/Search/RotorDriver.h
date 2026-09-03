//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_ROTORDRIVER_H
#define RDKIT_ROTORDRIVER_H

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <Geometry/point.h>

#include <array>
#include <functional>
#include <vector>

namespace RDKit {

//! Coordinate-only rotor driver.
/*!
  Base class for rotor driving a molecule.

  The rotors are rooted at the largest rigid fragment: every rotor moves the
  atoms on the side of its bond AWAY from that root.

  This side effect is that earlier rotors drive more atoms, see
    rotorSubsets and numRotorAtoms

  Scoring is delegated to optional callables.

  Main function is score = drive(rotor, angle)
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT RotorDriver {
 public:
  using ScoreFn = std::function<double(const double *pos, unsigned int nAtoms)>;

  explicit RotorDriver(const ROMol &mol, int confId = -1, ScoreFn fullFF = {},
                       ScoreFn fragmentFF = {});
  //! Drive an explicit set of rotor bond indices.
  RotorDriver(const ROMol &mol, const std::vector<unsigned int> &rotorBonds,
              int confId = -1, ScoreFn fullFF = {}, ScoreFn fragmentFF = {});

  unsigned int numAtoms() const { return d_nAtoms; }
  size_t numRotors() const { return d_rotors.size(); }

  //! Return The driven i-j-k-l torsion quartet for rotor `r` 
  std::array<unsigned int, 4> torsion(unsigned int r) const;
  //! Return atom indices for rotor r
  const std::vector<unsigned int> &movingAtoms(unsigned int r) const;

  //! Return the rotors contained in the set that torsion r rotates
  const std::vector<unsigned int> &rotorSubsets(unsigned int r) const;
  //! For each rotor, return the number of atoms moved
  const std::vector<unsigned int> &numRotorAtoms() const {
    return d_scaleOrder;
  }
  //! The atom the hierarchy is rooted at (in the largest rigid fragment); its
  //! side of every rotor bond is held fixed.  0 when there are no rotors.
  unsigned int rootAtom() const { return d_rootAtom; }

  //! Return the current torsion angle for rotor r
  double dihedralDeg(unsigned int r) const;

  //! Set the torsion angle for rotor r
  void setDihedral(unsigned int r, double angleDeg);

  //! Set the torsion for rotor r AND score
  double drive(unsigned int r, double angleDeg);

  //! Score the current working coordinates
  double score() const;
  double scoreFull() const;
  double scoreFragment() const;

  //! Return the working coordinate buffer
  const std::vector<double> &positions() const { return d_pos; }
  std::vector<double> &positions() { return d_pos; }

 private:
  void init(const ROMol &mol, const std::vector<unsigned int> &rotorBonds,
            int confId);
  void buildRotorHierarchy();
  RDGeom::Point3D pos(unsigned int a) const;
  void setPos(unsigned int a, const RDGeom::Point3D &p);

  struct Rotor {
    unsigned int i, j, k, l;           //!< driven dihedral quartet (axis = j-k)
    std::vector<unsigned int> moving;  //!< atoms rotated (the k-side, minus k)
  };

  unsigned int d_nAtoms = 0;
  unsigned int d_rootAtom = 0;
  std::vector<double> d_pos;  //!< working coordinates, 3*numAtoms
  std::vector<Rotor> d_rotors;
  //! per rotor: contained (rigidly-carried) rotor indices, largest movement
  //! first
  std::vector<std::vector<unsigned int>> d_contained;
  //! all rotor indices, largest movement first
  std::vector<unsigned int> d_scaleOrder;
  ScoreFn d_fullFF, d_fragFF;
};

}  // namespace RDKit

#endif
