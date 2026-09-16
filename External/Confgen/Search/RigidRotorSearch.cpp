//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Embedder/Embedder.h"
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
  //  so add Systematic as a baseline
  // The table decides Systematic vs the Thompson/Tree arm by rotor count, which
  // we do not have here -- so offer every mode the table can produce.
  std::vector<RigidRotorSearchMode> modes;
  bool tableHasSystematic = false;
  for (const auto m : RigidRotorSearchParams::AutoModeAtRotor) {
    if (m == RigidRotorSearchMode::Systematic) {
      tableHasSystematic = true;
      break;
    }
  }
  if (tableHasSystematic) {
    modes.push_back(RigidRotorSearchMode::Systematic);
  }
  // If we have a thompson budget or are auto selecting the budget
  //  use thompson, otherwise use the tree search
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

std::vector<double> basinLimitAngles(const std::vector<double> &angles,
                                     double currentDeg) {
  std::vector<double> out;
  const bool positive = currentDeg > 0.0;
  for (const double a : angles) {
    if ((a > 0.0) == positive && std::abs(a) > 1e-9) {
      out.push_back(a);
    }
  }
  if (out.empty()) {
    out.push_back(currentDeg);  // never strand a rotor with no candidates
  }
  return out;
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
    mode = RigidRotorSearchParams::autoModeForRotors(nRot);
    if (mode == RigidRotorSearchMode::Thompson &&
        !(resolved.thompsonBudget > 0 || resolved.thompson.autoBudget)) {
      mode = RigidRotorSearchMode::Tree;  // Thompson is switched off
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
