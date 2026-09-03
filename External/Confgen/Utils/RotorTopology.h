//  Copyright (C) 2026 Glysade Inc
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_CONFGEN_ROTORTOPOLOGY_H
#define RDKIT_CONFGEN_ROTORTOPOLOGY_H

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <vector>

namespace RDKit {

//! How a molecule's rotatable bonds are arrange
/*!
  Statistics of a rotatable bond
    How many rotors
    How many chains in the rotated group
    Number of rotors in the largest chain
    longestChainFraction - 1.0 all chain, 0.2 lots of fiddly bits in between
*/
struct RDKIT_FRAGMENTCONFGEN_EXPORT RotorTopology {
  unsigned int nRotors = 0;
  unsigned int nChains = 0;
  unsigned int maxChain = 0;  //!< XXX FIX ME: bad name rotors in the longest chain
  double longestChainFraction = 0.0;
};

RDKIT_FRAGMENTCONFGEN_EXPORT RotorTopology
rotorTopology(const ROMol &mol, const std::vector<unsigned int> &rotorBonds);

//! AUTO diversity radius: how far apart two output conformers must be to count
//! as distinct.
/*!
  Measured on Platinum + PDBbind (1049 molecules, rot>=7), widening 0.5 -> 0.75
  helps most where the rotors are concentrated in ONE long chain:

  | longestChainFraction | n   | delta %<1 | p     |
  | <= 0.40              | 301 | +2.33     | 0.143 |
  | 0.40-0.60            | 443 | +2.03     | 0.150 |
  | > 0.60               | 305 | +5.25     | 0.004 |

  Below rot 7 a wider radius only SHRINKS the ensemble (rot 0-4: 38 -> 15
  conformers) with no accuracy change, so the rule is gated on both.
*/
RDKIT_FRAGMENTCONFGEN_EXPORT double autoDiversityRms(const RotorTopology &topo);

//! Given a rotor set, return the rotor diversty RMS to use for uniqueness
RDKIT_FRAGMENTCONFGEN_EXPORT double diversityRmsForRotors(
    const std::vector<double> &byRotor, unsigned int nRotors,
    double fallback = 0.5);

}  // namespace RDKit

#endif
