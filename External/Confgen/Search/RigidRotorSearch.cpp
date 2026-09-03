//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Embedder/Fraglib.h"
#include "Search/RigidRotorSearch.h"

#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>

#include <iterator>
#include <string>

#include "Utils/ParamsIO.h"  // searchModeName
#include "Search/TreeSearch.h"
#include "Search/ThompsonSamplingSearch.h"
#include "Search/SystematicSearch.h"
#include "Utils/RotorTopology.h"
#include "Search/RigidRotorSearchExperimental.h"

namespace RDKit {

std::string RigidRotorSearch::validateParams(const RigidRotorSearchParams &sp,
                                             const std::string &ffVariant) const {
  if (!checkFF(ffVariant) && !isValidFF(ffVariant)) {
    return "this search builds its own MMFF terms and needs an MMFF force "
           "field (got \"" + ffVariant + "\")";
  }
  if (sp.energyWindow < 0.0) {
    return "search.energyWindow must be >= 0";
  }
  if (sp.defaultAngles.empty()) {
    return "search.defaultAngles must not be empty -- every search falls back "
           "to it when the sampler proposes nothing";
  }
  if (sp.torsionSampler && !sp.torsionSampler->isValid()) {
    return std::string("search.torsionSampler data did not load: ") +
           torsionSamplerStatusMessage(sp.torsionSampler->status());
  }
  return {};
}

std::vector<RigidRotorSearchMode> reachableSearchModes(
    const RigidRotorSearchParams &sp) {
  if (sp.searchMode != RigidRotorSearchMode::Auto) {
    return {sp.searchMode};
  }
  // Auto's pick depends on the molecule's rotor count, which we do not have --
  // but MOST of it is decided by the parameters alone.  Mirror the branch in
  // runRigidRotorSearch() and keep only the arms it can still take, so we
  // never fail a run over a knob belonging to a search that cannot be chosen.
  std::vector<RigidRotorSearchMode> modes;
  if (sp.autoSystematicMinRotors > 0) {
    modes.push_back(RigidRotorSearchMode::Systematic);
  }
  // The rotor count only selects BETWEEN Systematic and this arm; which of
  // Thompson/Tree the arm resolves to is fixed by the budget parameters.
  modes.push_back((sp.thompsonBudget > 0 || sp.thompson.autoBudget)
                      ? RigidRotorSearchMode::Thompson
                      : RigidRotorSearchMode::Tree);
  return modes;
}

std::string validateSearchParams(const RigidRotorSearchParams &sp,
                                 const std::string &ffVariant) {
  for (auto mode : reachableSearchModes(sp)) {
    const auto search = makeRigidRotorSearch(mode);
    const std::string err = search->validateParams(sp, ffVariant);
    if (!err.empty()) {
      return std::string(searchModeName(mode)) + " search: " + err;
    }
  }
  return {};
}

std::unique_ptr<RigidRotorSearch> makeRigidRotorSearch(
    RigidRotorSearchMode mode) {
  switch (mode) {
    case RigidRotorSearchMode::Thompson:
      return std::make_unique<ThompsonSamplingSearch>();
    case RigidRotorSearchMode::Merged:
      return std::make_unique<MergedSearch>();
    case RigidRotorSearchMode::Systematic:
      return std::make_unique<SystematicSearch>();
    case RigidRotorSearchMode::Tree:
    case RigidRotorSearchMode::Auto:  // resolved to Tree/Thompson by
                                      // runRigidRotorSearch
    default:
      return std::make_unique<TreeSearch>();
  }
}

RigidRotorSearchResult runRigidRotorSearch(const FragmentJoinerContext &ctx,
                                           const RigidRotorSearchParams &sp) {
  RigidRotorSearchResult out;
  if (!ctx.isValid()) {
    BOOST_LOG(rdWarningLog) << "Invalid context, skipping search..." << std::endl;;
    return out;
  }

  // Resolve any auto params here
  RigidRotorSearchParams resolved = sp;
  RigidRotorSearchMode mode = resolved.searchMode;
  if (mode == RigidRotorSearchMode::Auto) {
    const size_t nRot = ctx.rotorBonds.size();
    if (resolved.autoSystematicMinRotors > 0 && nRot >= resolved.autoSystematicMinRotors) {
      mode = RigidRotorSearchMode::Systematic;
    } else {
      mode = (resolved.thompsonBudget > 0 || resolved.thompson.autoBudget)
                 ? RigidRotorSearchMode::Thompson
                 : RigidRotorSearchMode::Tree;
    }
  }

  if (isAuto(resolved.diversityRmsThresh)) {
    resolved.diversityRmsThresh =
        mode == RigidRotorSearchMode::Tree
            ? 0.5
            : diversityRmsForRotors(resolved.autoDiversityRmsByRotor,
                                    static_cast<unsigned int>(ctx.rotorBonds.size()));
  }
  const RigidRotorSearchParams &spr = resolved;

  std::unique_ptr<RigidRotorSearch> searcher = makeRigidRotorSearch(mode);
  out.results = searcher->search(ctx, spr);
  out.budget = searcher->lastBudget();
  out.timedOut = searcher->timedOut();
  return out;
}

}  // namespace RDKit
