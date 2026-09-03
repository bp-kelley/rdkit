//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Rotational symmetry of a driven torsion bond j-k, used to FOLD equivalent
//  library-frame angles together (a CH3 cap has 3-fold symmetry, so its
//  0/120/240 junction angles are redundant).  Two interchangeable
//  implementations:
//
//    * BY RANKS  (torsionRotationalSymmetryByRanks) -- the general, data-free
//    path.
//      Canonical atom ranks WITHOUT tie-breaking ARE symmetry classes, so
//      equivalent terminal substituents are found by grouping each side's
//      non-axis neighbours by rank. [default]
//    * BY SMARTS (torsionRotationalSymmetryBySmarts) -- hand-curated
//       smarts defining torsion rules, requires external file
//
#ifndef RDKIT_TORSIONSYMMETRIES_H
#define RDKIT_TORSIONSYMMETRIES_H

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>

#include <string>
#include <vector>

namespace RDKit {

using SymmetryClass = unsigned int;

//! Molecule property key under which the per-atom symmetry classes are cached.
RDKIT_FRAGMENTCONFGEN_EXPORT extern const std::string SYMMETRY_CLASSES;

//! Return the toplogical (not rotational) symmetry class in order of atoms.
//! By VALUE: the classes are cached on the molecule, but getProp hands back a
//! copy, so a reference return would dangle.
RDKIT_FRAGMENTCONFGEN_EXPORT std::vector<SymmetryClass> getSymmetryClasses(
    const ROMol &mol);

//! The goal here is to find rotational symmetries for the given torsion
//!   this helps remove redundant sampling.
//! Result:  1 == no folding, 3 = fold module 120 deg, 360 = nada
RDKIT_FRAGMENTCONFGEN_EXPORT unsigned int torsionRotationalSymmetryByRanks(
    const ROMol &mol, const std::vector<SymmetryClass> &ranks, unsigned int j,
    unsigned int k);

//! Return the torsion rotation classes around atoms i,j
RDKIT_FRAGMENTCONFGEN_EXPORT unsigned int torsionRotationalSymmetryByRanks(
    const ROMol &mol, unsigned int j, unsigned int k);

//! Return the torsion rotation classes around atoms i,j using
//!  canned smarts patterns
RDKIT_FRAGMENTCONFGEN_EXPORT unsigned int torsionRotationalSymmetryBySmarts(
    const ROMol &mol, unsigned int j, unsigned int k);

//! Fold a library-frame angle set into one period of the given rotational
//! symmetry. sym<=1 returns the angles unchanged; sym==360 collapses to a
//! single angle.
RDKIT_FRAGMENTCONFGEN_EXPORT std::vector<double> foldAnglesByTorsionSymmetry(
    unsigned int sym, const std::vector<double> &angles, double minSeparation);

}  // namespace RDKit

#endif
