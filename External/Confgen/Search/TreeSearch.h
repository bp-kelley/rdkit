//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_TREESEARCH_H
#define RDKIT_TREESEARCH_H

#include <RDGeneral/export.h>
#include <vector>

#include "Search/RigidRotorSearch.h"  // RigidRotorSearch base (pulls FragmentJoinerContext + RigidRotorSearchParams)
#include "Search/RotorDriver.h"
#include "Search/InterFragScore.h"
#include "Search/SearchResult.h"
#include "Utils/DiagnosticsParams.h"

namespace ForceFields {
namespace MMFF {
class InterFragVdWContrib;
class TorsionAngleContrib;
}  // namespace MMFF
}  // namespace ForceFields

namespace RDKit {

//! The level-synchronous angle beam that TreeSearch drives.
/*!
  A beam search over the molecule's rotor angles.  Rotors are visited coarse to
  fine (RotorDriver::rotorsByScale): the widest-moving junction first, then the
  refinements nested inside it.  Because moving sets are rooted and nested,
  driving a rotor never disturbs a rotor already set on this path -- a parent
  was set before its children, and a child/sibling's moving set excludes the
  parent's torsion atoms
  -- so each level's choice is independent of the earlier ones.

  At each level every surviving partial conformer is branched over that rotor's
  candidate angles (RotorDriver::drive, coordinates only -- no molecule is ever
  rebuilt), scored on the working buffer, and the best `beamWidth` within
  `energyWindow` of the current best are carried forward.  The result is the
  final beam: a set of full-molecule coordinate buffers with their scores.

  This is the "search more angle combinations, cheaply" lever: the whole search
  is coordinate manipulation on one fixed topology with one force field built
  once.

  It is an implementation detail of TreeSearch -- exposed only because the tests
  drive it directly.
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT RotorTree {
 public:
  //! Takes the whole search params: the beam needs the SHARED settings
  //! (defaultAngles, diversityRmsThresh, randomSeed, energyWindow) as well as
  //! its own beamWidth/angleTolerance.
  explicit RotorTree(RotorDriver &driver,
                     RigidRotorSearchParams params = {});

  //! Override the candidate angles for one rotor (else defaultAngles is used).
  void setAngles(unsigned int rotor, std::vector<double> angles);

  //! Enable INCREMENTAL junction-local scoring for the deterministic beam.
  void enableIncremental(ForceFields::MMFF::InterFragVdWContrib *vdw,
                         ForceFields::MMFF::TorsionAngleContrib *tor);

  //! Run the search
  std::vector<SearchResult> search();

  //! DIAGNOSTIC toggles (incremental-score validation / stats); default
  //! all-off.
  void setDiagnostics(const DiagnosticsParams &d) { d_diag = d; }

 private:
  const std::vector<double> &anglesForRotor(unsigned int rotor) const;

  RotorDriver &d_driver;
  RigidRotorSearchParams d_params;
  std::vector<std::vector<double>> d_angles;  //!< per rotor; empty -> default
  //! Rescore only what each rotor changes.  Built ONCE by enableIncremental()
  //! (the packing is topology-invariant) and reused by every search().  Not
  //! tree-specific -- see IncrementalInterFragScore.
  IncrementalInterFragScore d_incremental;
  DiagnosticsParams d_diag;  //!< incremental-score validate/stats (default off)
};

//! Deterministic beam search over junction angles (the "tree" search).  Builds
//! a small set of fragment-conformer seeds (rootSeeds x fragConfBranch), and
//! for each seed runs a RotorTree level-synchronous beam over the rotor angles
//! (coarse->fine), carrying the best `beamWidth` partials within `energyWindow`
//! at each level.  Results are pooled across seeds, energy-sorted, and
//! diversity-pruned.
class RDKIT_FRAGMENTCONFGEN_EXPORT TreeSearch : public RigidRotorSearch {
 public:
  //! Tree adds the beam width to what the base already checks.
  std::string validateParams(const RigidRotorSearchParams &sp,
                             const std::string &ffVariant) const override;

  std::vector<SearchResult> search(const FragmentJoinerContext &ctx,
                                   const RigidRotorSearchParams &sp) override;
  // deterministic: lastBudget() stays 0 (base default)
};

}  // namespace RDKit

#endif
