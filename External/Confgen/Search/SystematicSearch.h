//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_SYSTEMATICSEARCH_H
#define RDKIT_SYSTEMATICSEARCH_H

#include <RDGeneral/export.h>
#include <vector>

#include "Search/RigidRotorSearch.h"

namespace RDKit {

//! Perform a systematic cartesian product pruning from the ground up
//!  When a rotor is proven to be > ewindow, drop it completely
class RDKIT_FRAGMENTCONFGEN_EXPORT SystematicSearch : public RigidRotorSearch {
 public:

  bool checkFF(const std::string &ff) const override { return ff.find("MMFF") != std::string::npos; }

  std::string validateParams(const RigidRotorSearchParams &sp,
                             const std::string &ffVariant) const override;

  std::vector<SearchResult> search(
      const FragmentZipperContext &ctx,
      const RigidRotorSearchParams &sp) override;
};

}  // namespace RDKit

#endif
