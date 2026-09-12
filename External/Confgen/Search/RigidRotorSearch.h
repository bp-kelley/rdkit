//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Rotor drives joined rigid fragments
//  returns the lowest energy ensemble found.
//
#ifndef RDKIT_CONFGEN_RIGIDROTORSEARCH_H
#define RDKIT_CONFGEN_RIGIDROTORSEARCH_H

#include <RDGeneral/export.h>
#include <memory>
#include <string>
#include <vector>

#include "Joiner/FragmentJoiner.h"
#include "Search/SearchParams.h"

#include <chrono>
#include "Search/SearchResult.h"

namespace RDKit {

//! A junction-angle search strategy via deterministic beam
class RDKIT_FRAGMENTCONFGEN_EXPORT RigidRotorSearch {
 public:
  virtual ~RigidRotorSearch() = default;
  virtual std::vector<SearchResult> search(
      const FragmentJoinerContext &ctx, const RigidRotorSearchParams &sp) = 0;
  //! If we have a sample budget, return it
  virtual unsigned int lastBudget() const { return 0; }

  //! Does this search REQUIRE a specific FF?
  virtual bool checkFF(const std::string &/*ff*/) const { return false; }

  //! Are these parameters viable for THIS search?
  /*!
    Returns an empty string when valid, otherwise a human readable error message

    \param sp         the search parameters
    \param ffVariant  the force field the joiner will score with -- separate
                      because it lives in FragmentJoinerParams, not here
  */
  virtual std::string validateParams(const RigidRotorSearchParams &sp,
                                     const std::string &ffVariant) const;

  //! True when validateParams() finds nothing wrong.
  bool isValid(const RigidRotorSearchParams &sp,
               const std::string &ffVariant) const {
    return validateParams(sp, ffVariant).empty();
  }
  //! Did the search time out?
  bool timedOut() const { return d_timedOut; }

 protected:
  //! Compute the timeout
  static std::chrono::steady_clock::time_point timeOut(
      const RigidRotorSearchParams &sp) {
    if (sp.timeBudgetMs <= 0)
      return std::chrono::steady_clock::time_point::max();
    return std::chrono::steady_clock::now() +
           std::chrono::milliseconds(sp.timeBudgetMs);
  }
  bool timedOut(const std::chrono::steady_clock::time_point &timeout) {
    if (timeout == std::chrono::steady_clock::time_point::max()) return false;
    if (std::chrono::steady_clock::now() < timeout) return false;
    d_timedOut = true;
    return true;
  }
  bool d_timedOut = false;
};

//! Construct the RigidRotorSearch
RDKIT_FRAGMENTCONFGEN_EXPORT std::unique_ptr<RigidRotorSearch>
makeRigidRotorSearch(RigidRotorSearchMode mode);

//! Which searches these parameters can actually run.
/*!
  Some Auto params switch between modes at different rotor counts,
  this takes that into account.
*/
RDKIT_FRAGMENTCONFGEN_EXPORT std::vector<RigidRotorSearchMode>
reachableSearchModes(const RigidRotorSearchParams &sp);

//! Validate the search parameters against the searches they can actually run.
RDKIT_FRAGMENTCONFGEN_EXPORT std::string validateSearchParams(
    const RigidRotorSearchParams &sp, const std::string &ffVariant);

//! results and timing for a rigid rotor search
struct RDKIT_FRAGMENTCONFGEN_EXPORT RigidRotorSearchResult {
  std::vector<SearchResult> results;
  unsigned int budget = 0;   //! Number of samples used from the torsion samplers
  bool timedOut = false;
};

//! Run a search on a joined context
RDKIT_FRAGMENTCONFGEN_EXPORT RigidRotorSearchResult runRigidRotorSearch(
    const FragmentJoinerContext &ctx, const RigidRotorSearchParams &sp);

}  // namespace RDKit

#endif
