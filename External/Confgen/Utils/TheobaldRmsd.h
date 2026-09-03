//  Copyright (C) 2026 Pande Group, Glysade Inc
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

//  NOTE: code ported by claude.ai
//   the decision was made to drop explicit SSE and let the compiler do the work
//   using the original code, we could perhaps squeeze some more performance
//    out of this using the original code on x64 architectures
//
//  Portable scalar port of the Pande group's Theobald QCP RMSD (IRMSD),
//  atom-major. The original (theobald_rmsd.c) vectorises the 3x3 inner product
//  with SSE; here that is a plain loop (fragments are small and each conf-gen
//  runs single-threaded), while the clever part -- the closed-form largest
//  eigenvalue of the QCP key matrix via a Newton solve on the characteristic
//  quartic (msdFromMandG) -- is carried over as-is. Reference: Liu, Agrawal &
//  Theobald, J Comput Chem 31:1561 (2010).
//
#ifndef RDKIT_THEOBALD_RMSD_H
#define RDKIT_THEOBALD_RMSD_H

#include <RDGeneral/export.h>
#include <Geometry/point.h>
#include <vector>

namespace RDKit {

//! QCP mean-square deviation between two structures already CENTERED on their
//! centroids, stored atom-major: a = {x0,y0,z0, x1,y1,z1, ...} (length
//! 3*nAtoms). `Ga`/`Gb` are the sums of squared coordinates of a/b (i.e. tr(A^T
//! A), tr(B^T B)). Returns the optimal-superposition MSD (rotation removed)
//! without ever forming the rotation matrix.  Negative round-off results are
//! clamped to 0.
RDKIT_FRAGMENTCONFGEN_EXPORT double theobaldMSD(unsigned int nAtoms,
                                                const double *a,
                                                const double *b, double Ga,
                                                double Gb);

//! Optimal-superposition (rotation + translation removed) RMSD between
//! coordinate sets `a` and `b`, over the atoms listed in `idx` (indices into
//! both).  Centres the selected atoms, builds the QCP inputs and calls
//! theobaldMSD.  `a` and `b` must index-align (a[idx[k]] corresponds to
//! b[idx[k]]).
RDKIT_FRAGMENTCONFGEN_EXPORT double qcpRmsd(
    const std::vector<RDGeom::Point3D> &a,
    const std::vector<RDGeom::Point3D> &b,
    const std::vector<unsigned int> &idx);

//! QCP RMSD between two atom-major coordinate buffers (a,b hold x,y,z per
//! atom), matching a's atom idxA[k] to b's atom idxB[k] (k = 0..idxA.size()-1).
//! The two index lists let the caller apply an atom permutation to b (e.g. a
//! symmetry automorphism) without copying coordinates.  idxA and idxB must have
//! equal length.
RDKIT_FRAGMENTCONFGEN_EXPORT double qcpRmsd(
    const std::vector<double> &a, const std::vector<double> &b,
    const std::vector<unsigned int> &idxA,
    const std::vector<unsigned int> &idxB);

}  // namespace RDKit

#endif
