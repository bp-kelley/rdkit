//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Profiling tools to find hotspots in confgen
//
#ifndef RDKIT_JOINERPROFILING_H
#define RDKIT_JOINERPROFILING_H

#include <atomic>
#include <chrono>

namespace RDKit {
namespace detail {

//! accumulators for various parts of the assembly
struct JoinProf {
  std::atomic<long long> tEmbed{0}, tPlace{0}, tDrive{0}, tScore{0}, tDiv{0};
  std::atomic<long long> nEmbeds{0}, nSamples{0}, nDrives{0}, nScores{0},
      nKept{0}, nOut{0};

  //! --- speed breakdown (systematic path)
  //! ------------------------------------- Per-fragment-class embedding cost,
  //! indexed by FragmentClass (Rigid, SmallRing, LargeRing, Acyclic,
  //! Exhaustive, Fast).
  static constexpr int kNClass = 6;
  std::atomic<long long> tEmbedClass[kNClass];
  std::atomic<long long> nEmbedClass[kNClass];
  std::atomic<long long> nConfClass[kNClass];  //!< conformers KEPT per class
  //! Scoring work, split by kind: the decomposed inter-fragment terms the
  //! combine evaluates per candidate (vdW pair sums + junction torsion) vs
  //! whole-molecule MMFF evaluations.
  std::atomic<long long> nScoreInterFrag{
      0};  //!< inter-fragment (vdW + junction torsion) evals
  std::atomic<long long> nScoreFullMMFF{0};  //!< full-molecule MMFF evals
  //! Combine-step pool pressure: conformers materialised vs surviving the
  //! per-node prune.
  std::atomic<long long> nPoolMade{0}, nPoolKept{0};
  //! Rejection accounting: WHY does a candidate combination die?  (tests the
  //! "an early/deep node minimum tightens the energy window and kills
  //! conformers" hypothesis)
  std::atomic<long long> nRejQuick{
      0};  //!< live runningMin+eWindow+5 quick rejection
  std::atomic<long long> nRejWindow{
      0};  //!< final per-node window: E > nodeMin + eWindow
  std::atomic<long long> nRejRms{0};   //!< RMSD duplicate at the node
  std::atomic<long long> nRejCap{0};   //!< hit poolCap
  std::atomic<long long> nRejTrim{0};  //!< discarded by the in-loop memory trim
  std::atomic<long long> tCombine{0};  //!< time inside the bottom-up combine
  std::atomic<long long> tPrune{
      0};  //!< time in per-node prune (energy window + RMSD)

  JoinProf() {
    for (int i = 0; i < kNClass; ++i) {
      tEmbedClass[i] = 0;
      nEmbedClass[i] = 0;
      nConfClass[i] = 0;
    }
  }
};

//! Global assembly profiler
JoinProf &prof();
//! true if profiling is enabled
bool profiling();
//! Enable or disable profiling.
void setJoinerProfiling(bool on);

//! ns timer
inline long long nowNs() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

}  // namespace detail
}  // namespace RDKit

#endif
