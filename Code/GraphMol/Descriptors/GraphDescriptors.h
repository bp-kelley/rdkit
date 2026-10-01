//
//  Copyright (C) 2026 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file GraphDescriptors.h

  \brief C++ ports of the topological descriptors from
  rdkit.Chem.GraphDescriptors that were previously only available from
  Python.

  The implementations follow the Python code closely so that the values
  match those returned by rdkit.Chem.Descriptors.

*/
#include <RDGeneral/export.h>
#ifndef RD_GRAPHDESCRIPTORS_H
#define RD_GRAPHDESCRIPTORS_H

#include <string>

namespace RDKit {
class ROMol;
namespace Descriptors {

//! Chi0: From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422,
//! (1991)
/*!
  Uses the heavy-atom degree of each atom; atoms with degree 0 are skipped.
*/
RDKIT_DESCRIPTORS_EXPORT double calcChi0(const ROMol &mol);
const std::string chi0Version = "1.0.0";

//! Chi1: From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422,
//! (1991)
RDKIT_DESCRIPTORS_EXPORT double calcChi1(const ROMol &mol);
const std::string chi1Version = "1.0.0";

//! Balaban's J value for a molecule
/*!
  Chem. Phys. Lett. vol 89, 399-404, (1982)
*/
RDKIT_DESCRIPTORS_EXPORT double calcBalabanJ(const ROMol &mol);
const std::string balabanJVersion = "1.0.0";

//! Bertz's complexity index
/*!
  From S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981)

  \param mol    the molecule of interest
  \param cutoff vertices are considered topologically identical if their
                sorted distance vectors are equal out to the cutoff-th
                nearest neighbor
*/
RDKIT_DESCRIPTORS_EXPORT double calcBertzCT(const ROMol &mol,
                                            unsigned int cutoff = 100);
const std::string bertzCTVersion = "2.0.0";

//! Information content of the coefficients of the characteristic polynomial
//! of the adjacency matrix of a hydrogen-suppressed graph of a molecule.
/*!
  From Eq 6 of D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67, 4517-4533
  (1977)

  \param mol  the molecule of interest
  \param avg  if set, the information content divided by the total
              population (AvgIpc) is returned

  <b>Note</b> like the python implementation, this uses the Le
  Verrier-Faddeev-Frame method in double precision. For large molecules
  (roughly 80+ heavy atoms) the polynomial coefficients become too large to
  be represented exactly and the results are numerically unstable; in those
  cases values may differ from the python implementation, which is equally
  affected.
*/
RDKIT_DESCRIPTORS_EXPORT double calcIpc(const ROMol &mol, bool avg = false);
const std::string ipcVersion = "1.0.0";

//! Average information content (Eq 7 of Bonchev & Trinajstic), see calcIpc()
RDKIT_DESCRIPTORS_EXPORT double calcAvgIpc(const ROMol &mol);
const std::string avgIpcVersion = "1.0.0";

}  // namespace Descriptors
}  // namespace RDKit
#endif
