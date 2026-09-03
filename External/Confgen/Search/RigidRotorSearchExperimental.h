//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Experimental junction-angle search strategies, filled in one at a time (see
//  the checklist in the design plan).  Each is a drop-in RigidRotorSearch
//  selectable via FragmentJoinerParams::searchMode.  Until implemented, each
//  delegates to TreeSearch so the pipeline + bench selector work end-to-end.
//
#ifndef RDKIT_RIGIDROTORSEARCHEXPERIMENTAL_H
#define RDKIT_RIGIDROTORSEARCHEXPERIMENTAL_H

#include <RDGeneral/export.h>
#include <vector>

#include "Search/RigidRotorSearch.h"  // RigidRotorSearch base (pulls FragmentJoinerContext + RigidRotorSearchParams)

namespace RDKit {

//! EXPERIMENTAL -- run Systematic AND Thompson, merge their ensembles.
/*!
  Not a new search policy: the two existing searches reach largely DIFFERENT
  molecules.  Measured on rot>=11 (n=43): systematic 46.51%, thompson 41.86%,
  but they agree on only 27.9% and their union is 60.47% -- more headroom than
  any single-policy change tried so far.

  The union is an ORACLE number: it counts a molecule solved if EITHER run
  produced a sub-1A pose.  A merged ensemble still has to RANK that pose into a
  capped output, so this exists to measure how much of the oracle survives
  selection.  If little does, the limit is ranking, not search.
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT MergedSearch : public RigidRotorSearch {
 public:
  std::vector<SearchResult> search(
      const FragmentJoinerContext &ctx,
      const RigidRotorSearchParams &sp) override;
  unsigned int lastBudget() const override { return d_lastBudget; }

 private:
  unsigned int d_lastBudget = 0;
};

}  // namespace RDKit

#endif
