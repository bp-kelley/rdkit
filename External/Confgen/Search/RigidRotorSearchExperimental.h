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
//  selectable via FragmentZipperParams::searchMode.  Until implemented, each
//  delegates to TreeSearch so the pipeline + bench selector work end-to-end.
//
#ifndef RDKIT_RIGIDROTORSEARCHEXPERIMENTAL_H
#define RDKIT_RIGIDROTORSEARCHEXPERIMENTAL_H

#include <RDGeneral/export.h>
#include <vector>

#include "Search/RigidRotorSearch.h"  // RigidRotorSearch base (pulls FragmentZipperContext + RigidRotorSearchParams)

namespace RDKit {

//! EXPERIMENTAL -- run Systematic AND Thompson, merge their ensembles.
/*!
  Very expensive with small gains in practice.
  Still, it is interesting and worth review.
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT MergedSearch : public RigidRotorSearch {
 public:
  std::vector<SearchResult> search(
      const FragmentZipperContext &ctx,
      const RigidRotorSearchParams &sp) override;
  unsigned int lastBudget() const override { return d_lastBudget; }

 private:
  unsigned int d_lastBudget = 0;
};

}  // namespace RDKit

#endif
