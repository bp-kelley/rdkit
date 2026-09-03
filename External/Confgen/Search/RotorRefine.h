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
  Takes the best `nRefine` results and walks each rotor by +/- step, keeping any
  move that lowers the score; the step halves each pass and the descent stops
  when a pass improves nothing.  Results are re-sorted by score afterwards.

  This is a LOCAL move: it improves a conformer within the basin the search
  already put it in, and cannot cross into another one.  It exists because a
  search can find the right basin and sit a few degrees off the bottom of it.

  Stays in ROTOR SPACE by construction, so refined coordinates remain
  expressible as rotor offsets -- free-atom minimisation would break that
  representation and is deliberately not done (see docs/thompson-sampling.md).

  Search-agnostic: it needs a driver positioned on a molecule and the results to
  polish, nothing about how those results were produced.  It began life inside
  the Thompson search, which is why only that search could use it.

  \param drv        driver over the rotors to walk (its positions are clobbered)
  \param out        results to refine, IN PLACE; re-sorted by score on return
  \param nRefine    how many of the leading results to refine (0 = none)
  \param stepDeg    initial descent step in degrees
  \param nPasses    number of halving passes
  \param objective  optional alternative score; when null, drv.score() is used
*/
RDKIT_FRAGMENTCONFGEN_EXPORT void refineRotorsInPlace(
    RotorDriver &drv, std::vector<SearchResult> &out, unsigned int nRefine,
    double stepDeg, unsigned int nPasses,
    const RotorDriver::ScoreFn &objective = {});

}  // namespace RDKit

#endif
