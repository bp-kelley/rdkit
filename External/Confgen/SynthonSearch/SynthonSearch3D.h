//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//   3D search of a combinatorial synthon library.
//
//   The library is enumerable in principle and astronomically large in
//   practice, so the search never materialises it: it samples.
//
//   THE SEARCH IS A MULTI-TRAJECTORY GREEDY SWEEP + COORDINATE REFINEMENT.
//
//   Empirically, shape searches can get stuck in a local minimum.
//   Using multiple trajectories helps break out of this cage.
//
#ifndef RDKIT_CONFGEN_SYNTHONSEARCH3D_H
#define RDKIT_CONFGEN_SYNTHONSEARCH3D_H

#include <RDGeneral/export.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/GaussianShape/ShapeInput.h>

#include <limits>
#include <memory>
#include <string>
#include <array>
#include <optional>
#include <vector>

#include "SynthonSearch/EnumerateSynthons3D.h"
#include "Utils/ParamSentinels.h"  // AutoI / AutoR / Disabled

namespace RDKit {

//! Score one assembled product against whatever the caller is looking for.
/*!
  We use std::nullopt to indicate when we cannot score a product such
  as confgen or other failures.  This is a sentinel so we don't
  poison results.
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT SynthonProductScorer {
 public:
  virtual ~SynthonProductScorer() = default;
  //! n.b. we "own" the product during the scoring, it is up
  //!  to the caller to send in a fresh copy.
  //!  ASSUMPTION: higher scores are better
  virtual std::optional<double> score(ROMol &product) const = 0;

};

class RDKIT_FRAGMENTCONFGEN_EXPORT ShapeScorer
    : public SynthonProductScorer {
 public:
  //! ShapeScorer, return the tanimotoscore/2
  //! \param queryConfId    -1 uses ALL the query's conformers, not just one
  //! \param allCarbonRadii  true treats every atom as carbon. [default True]
  //! \param colorWeight   weight on color v shape; 0.5 is the plain mean [default]
  /*!
    `colorWeight` moves the objective off the 50/50 mean:
    `score = (1 - w) * shape + w * color`, with 0.5 reproducing the previous
  */
  explicit ShapeScorer(const ROMol &query, int queryConfId = -1,
                              bool allCarbonRadii = true,
                              double colorWeight = 0.5);
  ~ShapeScorer() override;

  //! Return the current weight on color
  double colorWeight() const { return d_colorWeight; }

  //! Note: score over all the conformers of product
  std::optional<double> score(ROMol &product) const override;

  //! The same best conformer as score(), split into its scoring components
  /*!
    \return {combined, shape, color}, or nullopt if unscorable.
    Combined is half of the combined score
  */
  std::optional<std::array<double, 3>> scoreComponents(ROMol &product) const;

 private:
  std::unique_ptr<GaussianShape::ShapeInput> d_queryShape;
  bool d_allCarbonRadii = true;
  double d_colorWeight = 0.5;
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT SynthonHit {
  std::vector<unsigned int> reagents;
  double score = -std::numeric_limits<double>::infinity();
  ROMOL_SPTR mol; // optional, null if requested
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT SynthonSearchResult {
  std::vector<unsigned int> reagents;
  double score = -std::numeric_limits<double>::infinity();
  unsigned int evaluations = 0;  //!< products actually assembled AND scored
  unsigned int unscorable = 0;   //!< assembly/scoring failures skipped

  //! Timing
  double confgenMs = 0.0;  //!< assembling products (zip + conformers)
  double scoreMs = 0.0;    //!< scoring them (e.g. shape overlay)

  //! Number of hits drawn from cache
  unsigned int cacheHits = 0;

  //! Duplicates in the draw stage from this batch.
  unsigned int duplicatesRejected = 0;

  //! Best hits currently found in best to worst order
  std::vector<SynthonHit> best;
};

//! Parameters for synthonSearch3D.
/*!
  Note: searches are done in random seeded "trajectories".  More
  trajectories increases sampling.
  
  A trajectory is independent in search but can share results between
  multile concurrent trajectories.

  In general increasing the number of trajectories increases the likelihood
  of returning the "best" match.
*/
struct RDKIT_FRAGMENTCONFGEN_EXPORT SynthonSearch3DParams {
  unsigned int numTrajectories = 2;
  int randomSeed = 0xf00d;
  int seedStride = 1000;

  unsigned int samplesPerReagent = 10;
  bool largestFirst = true;
  int firstPosition = AutoI;  //!< -1 leaves largest-first ordering unchanged

  unsigned int refineIters = 3;
  int pairRefineTopK = 16;  //!< 0 disables the joint pair-refine

  unsigned int numThreads = 0;  //!< one shared pool; 0 = hardware concurrency
  int numBestProducts = 10;

  //! Also retain EVERY product scoring at or above this, not just the top N.
  /*!
    A fixed count throws away real hits: a library rich in matches loses
    everything past `numBestProducts` even when those products score as well
    as the best in a library that was kept entire. A score threshold retains
    on merit instead, so what comes back does not depend on which library a
    product happened to live in.

    `Disabled` (0) = top `numBestProducts` only. Note 0 cannot mean "keep
    everything" here: combo scores live in [0, 1], so a 0.0 threshold would
    admit every product ever scored -- raise `numBestProducts` for that.

    Retention is never fewer than `numBestProducts`, so the threshold only
    ever ADDS hits; a run that finds nothing above it still reports its best.

    COSTS PER HIT, not per search: hits above the threshold are
    *materialized* (re-assembled with conformers and re-scored), unlike the
    reagents+score pairs the search keeps for free. Hence `maxKeptProducts`.
  */
  double keepScoreThreshold = Disabled;

  //! Ceiling on retained products, guarding `keepScoreThreshold`.
  /*!
    A permissive threshold on a large library can clear thousands of
    products, and each one materialized holds a molecule with conformers.
    The highest-scoring are kept. `Disabled` (0) = no ceiling.
  */
  int maxKeptProducts = 1000;

  //! Query heavy-atom count, enabling the size filter.  0 = filter off.
  unsigned int queryHeavyAtoms = Disabled;
  //! Size-filter thresholds, as [minimum, maximum].
  /*!
    Applied to coordinate refinement and the pair cross product ONLY, never
    to the position sweeps: refinement candidates are complete assignments
    competing against the incumbent, so a skip is outcome-neutral, while a
    sweep takes an argmax over partial assignments and drops would move it.

    A candidate is skipped without being assembled when
    `min(nq,np)/max(nq,np)` falls below the cutoff, where `np` comes from
    cached per-synthon atom counts and costs a couple of additions against the
    ~5 ms an assembly costs.

    The cutoff is the running Nth-best score CLAMPED into
    [pruneMinimum, pruneMaximum]. The minimum lets hopeless candidates go
    before any hit exists; the maximum stops one strong hit from suppressing
    every other candidate that could still be worth reporting (a combo of 0.7
    is already a good match).

    HEURISTIC, not a proven bound: Gaussian shape and color are not
    rigorously bounded by a heavy-atom ratio -- color depends on which
    pharmacophore features are present. Measure the false-negative rate
    before relying on it.

    FOR THE SAME REASON THE FILTER IS SHAPE-ONLY, and there is no equivalent
    cheap estimate for color. So it interacts badly with
    `ShapeScorer::colorWeight` above 0.5: the objective moves toward a
    quantity the filter cannot estimate, and its rejections become false
    negatives on precisely the candidates the raised weight was meant to
    surface. The search warns when both are on; prefer pruneMinimum=Disabled
    for a color-weighted search whose recall you need to trust.
  */
  double pruneMinimum = Disabled;
  double pruneMaximum = 1.0;
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT SynthonTrajectorySummary {
  unsigned int trajectory = 0;
  int randomSeed = 0;
  std::vector<unsigned int> reagents;
  double score = -std::numeric_limits<double>::infinity();
  unsigned int logicalRequests = 0;
  unsigned int logicalScored = 0;
  unsigned int logicalUnscorable = 0;
  unsigned int cacheHits = 0;
  unsigned int crossTrajectoryCacheHits = 0;
};

//! Cost of one position's sweep: the unit greedy search is actually made of.
/*!
  The sweep scores every reagent at a position against `samplesPerReagent`
  completions, so ONE position sweep is one full pass over that position. Total
  search cost is the sum of these, which is why cost is linear in
  sum(reagents) rather than in the product space.

  Reported per position because the positions are not alike: a 3-component
  library's central linker and its terminal caps differ in reagent count, in
  product size (hence assembly cost), and in how much the tuple cache absorbs.
  A single aggregate rate hides all of that.

  `confgenMs` and `scoreMs` are summed across worker threads, so on N threads
  they total roughly N x wallMs.
*/
struct RDKIT_FRAGMENTCONFGEN_EXPORT SynthonPositionSweep {
  unsigned int position = 0;
  unsigned int reagents = 0;     //!< reagents at this position
  unsigned int evaluations = 0;  //!< products physically built and scored
  unsigned int cacheHits = 0;    //!< requests answered without building
  double wallMs = 0.0;
  double confgenMs = 0.0;  //!< assembling products, summed over threads
  double scoreMs = 0.0;    //!< shape+color overlay, summed over threads
};

//! Accounting which separates logical search work from physical computation.
struct RDKIT_FRAGMENTCONFGEN_EXPORT MultipleTrajectoryStats {
  std::vector<SynthonTrajectorySummary> trajectories;

  //! Requests made by all independent search states, including cache hits.
  unsigned int logicalRequests = 0;
  unsigned int logicalScored = 0;
  unsigned int logicalUnscorable = 0;

  //! Unique tuple evaluations performed during search (scored + unscorable).
  unsigned int uniqueSearchAttempts = 0;
  unsigned int uniqueSearchScored = 0;
  unsigned int uniqueSearchUnscorable = 0;

  unsigned int cacheHits = 0;
  unsigned int crossTrajectoryCacheHits = 0;
  unsigned int intraTrajectoryCacheHits = 0;

  unsigned int finalistEvaluations = 0;

  //! Candidates skipped by the heavy-atom size filter, never assembled.
  unsigned int sizeFiltered = 0;

  //! Reagent-column skips: no completion of that reagent could reach the
  //! cutoff, so it was dropped from a sweep or refine pass without assembly.
  unsigned int reagentsFiltered = 0;

  //! Skips of synthons the library knows are unbuildable (tombstoned).  Not a
  //! prune: these products do not exist, whatever the cutoff.
  unsigned int deadReagentsSkipped = 0;

  //! Per-position sweep costs, in sweep order.  See SynthonPositionSweep.
  std::vector<SynthonPositionSweep> positionSweeps;
};

//! Dump the parameters used for this search
RDKIT_FRAGMENTCONFGEN_EXPORT std::string describeParams(
    const SynthonSearch3DParams &params, const EnumerateSynthons3D &lib);

//! Run independent greedy+refine trajectories through one shared scorer.
RDKIT_FRAGMENTCONFGEN_EXPORT SynthonSearchResult
synthonSearch3D(const EnumerateSynthons3D &lib,
		const SynthonProductScorer &scorer,
		const SynthonSearch3DParams &params = {},
		MultipleTrajectoryStats *stats = nullptr);

}  // namespace RDKit

#endif
