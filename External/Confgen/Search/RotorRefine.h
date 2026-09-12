//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_CONFGEN_ROTORREFINE_H
#define RDKIT_CONFGEN_ROTORREFINE_H

#include <RDGeneral/export.h>
#include <vector>

#include "Search/RotorDriver.h"
#include "Search/SearchResult.h"

namespace RDKit {

//! Coordinate descent over the rotors of already-found conformers.
/*!
  This is a local refinement to optimize an existing set of
  rotors, i.e. a basin search or small "shrug"

  In general, this is not useful, but is kept for posterity.

  \param drv        driver over the rotors to walk (its positions are clobbered)
  \param out        results to refine, IN PLACE; re-sorted by score on return
  \param nRefine    how many of the leading results to refine (0 = none)
  \param stepDeg    initial descent step in degrees
  \param nPasses    number of halving passes
  \param objective  optional alternative score; when null, drv.score() is used
*/
RDKIT_FRAGMENTCONFGEN_EXPORT void refineRotors(
    RotorDriver &drv, std::vector<SearchResult> &out, unsigned int nRefine,
    double stepDeg, unsigned int nPasses,
    const RotorDriver::ScoreFn &objective = {});

}  // namespace RDKit

#endif
