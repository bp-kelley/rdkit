//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_THOMPSONSAMPLINGSEARCH_H
#define RDKIT_THOMPSONSAMPLINGSEARCH_H

#include <RDGeneral/export.h>
#include <vector>

#include "Search/RigidRotorSearch.h"

namespace RDKit {

//! Joint Thompson-sampling assembly search.
//!   here the bandit arms are fragment conf AND rotor angle
//! Angles are assumed to be an informed prior based on the sampler angles
//!  again with a uniform grid for non matched angles
//! Reward is the energy, and budget should scale with num rotors (and # fragment confs)
class RDKIT_FRAGMENTCONFGEN_EXPORT ThompsonSamplingSearch
    : public RigidRotorSearch {
 public:
  //! Thompson adds its budget ordering to what the base already checks.
  std::string validateParams(const RigidRotorSearchParams &sp,
                             const std::string &ffVariant) const override;

  std::vector<SearchResult> search(
      const FragmentJoinerContext &ctx,
      const RigidRotorSearchParams &sp) override;
  unsigned int lastBudget() const override { return d_lastBudget; }

 private:
  unsigned int d_lastBudget = 0;
};

}  // namespace RDKit

#endif
