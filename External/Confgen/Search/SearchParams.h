//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Plain parameter structs for the assembly search strategies.  These replace
//  the per-strategy environment-variable knobs; the top layer (a CLI/tool)
//  populates them, the library reads them.  Defaults reproduce the historical
//  (all-knobs-off) behavior.
//
#ifndef RDKIT_CONFGEN_SEARCHPARAMS_H
#define RDKIT_CONFGEN_SEARCHPARAMS_H

#include <RDGeneral/export.h>

#include <memory>
#include <vector>

#include "Sampler/TorsionSampler.h"   // TorsionSampler
#include "Utils/DiagnosticsParams.h"  // DiagnosticsParams
#include "Utils/ParamSentinels.h"     // kAutoD / kDisabled / resolveAuto

namespace RDKit {

//! Deterministic beam ("tree") search knobs.  See docs/tree-beam.md.
//! Parameters specific to the deterministic beam (TreeSearch) ONLY.
//! Anything read by more than one search lives on RigidRotorSearchParams --
//! it used to live here, which made shared settings look tree-specific.
struct RDKIT_FRAGMENTCONFGEN_EXPORT TreeSearchParams {
  unsigned int beamWidth = 50;  //!< partial conformers kept per level
  double angleTolerance = 15.0;
};

//! Thompson Search limit for output conformers.
enum class OutputSelection {
  Energy,   //!< lowest energy first
  Diverse,  //!< greedy RMSD-diverse over the energy-sorted pool
  Stratify  //!< energy-stratified stride across the pool
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT ThompsonParams {
  //! --- informed prior + budget scaling ---
  double priorStrength = 3.0;     //!< Beta alpha for sampler-preferred angles (>=1)
  double backstopStepDeg = 0.0;   //!< uniform-grid backstop step, deg (0 = none), alpha 1
  double sizePriorExp = 0.0;      //!< scale prior by (movingAtoms/max)^exp (0 = full prior)
  double noveltyAngleDeg = 0.0;  //!< in-sweep torsion-fingerprint novelty, deg (0 = coord RMSD)
  
  bool autoBudget = true;         //!< auto-scale #draws with rotors + frag-conf arms
  unsigned int perRotor = 400;    //!< draws per rotor
  unsigned int perFragConf = 30;  //!< draws per fragment-conf arm
  unsigned int minBudget = 300;   //!< minimum # samples (min)
  unsigned int maxBudget = 64000; //!< max samples

  unsigned int maxConfs = 0;      //!< max number of output conformers (0 = nolimit)
  OutputSelection outMode = OutputSelection::Energy;  //!< selection when capped
  bool flatContext = false;  //!< collapse the tree-descent context
  
  unsigned int refineSteps = 0; //! Sample K lower basins for better energyies
  double refineStepDeg = 8.0;     //!< initial coordinate-descent step
  unsigned int refinePasses = 3;  //!< step halves each pass
  
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT SystematicParams {
  double eWindow = 15.0;  //!< per-node energy window (kcal)
  long maxPoolConfs = 4000; //!< This is per pool, see combinePoolMultiplier as well
  int maxFragmentConfs = 16; //!< leaf fragment-conf limit

  double finalRms = kAutoD; //!< finalRms thresh, kAutoD is auto selectin
  double nodeRms = kAutoD;  //!< per node pruning
  bool energyOnly = false;  //!< energy-window retention only (no RMSD-diverse)

  double upperEnergyWindow = 5.0; //<! drop early candidates prior to stabilizing MMFF basin

  double dedupEnergyBand = 0.0; //!< Only RMS check within this energy window
  bool junctionAngleTerms = false; //!< use a full FF not vdw only for driving
  int combinePoolMultiplier = 4; //!< When to trim pools to prevent memory explosion
};

//! Which junction-angle search runs.
enum class JunctionBasinAngles {
  All,       //!< every junction gets the full (tolerance/basin) angle set
  None,      //!< every junction gets base angles only
  ChainOnly  //!< tolerance ONLY on chain-chain junctions
             //!<   XXX FIX ME -> I think this is wrong
};

inline bool useBasinAnglesForRotor(JunctionBasinAngles mode, const ROMol &mol,
				   unsigned int j, unsigned int k) {
  if (mode == JunctionBasinAngles::All) return true;
  if (mode == JunctionBasinAngles::None) return false;
  // ChainOnly: a junction is inter fragment when it is ring-adjacent, intra otherwise
  const RingInfo *ri = mol.getRingInfo();
  if (!ri || !ri->isInitialized()) return true;
  return ri->numAtomRings(j) == 0 && ri->numAtomRings(k) == 0;
}

enum class RigidRotorSearchMode {
  Auto,      //!< Thompson when budgeted, else the Tree beam XXX FIX ME -> bad auto
  Tree,      //!< greedy deterministic beam using per-edge minimum energies
  Thompson,  //!< one-armed-bandit sampling (maybe good for non-separable force
             //!< fields)
  Systematic,          //!< bottom-up systematic combine + prune
  Merged  //!< experimental: run Systematic + Thompson, merge ensembles
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT RigidRotorSearchParams {
  //! AUTO diversity radius per ROTOR COUNT: index by the molecule's rotatable
  //! bond count, clamped to the last entry.  Only consulted when
  //! diversityRmsThresh is AUTO (-1).
  /*!
    A floppy molecule's conformer pool is dominated by near-duplicates, so a
    tight "these are the same pose" radius fills a fixed output budget with one
    basin's worth of variations.  Widening it forces the same budget to span
    more of the space.  Measured on Platinum + PDBbind (1049 molecules, rot>=7):
    +2.0 to +7.1 %<1 at 0.75, with conformer counts unchanged.

    NOT applied to the Tree search -- see autoDiversityRmsForMode().  There the
    same number prunes the BEAM, so widening it explores less rather than
    selecting better: measured 0 wins / 6 losses, p=0.031.
  */
  std::vector<double> autoDiversityRmsByRotor = {0.5, 0.5, 0.5, 0.5, 0.5,
                                                 0.5, 0.5, 0.5, 0.75};

  //! --- SHARED across searches (formerly on TreeSearchParams) --------------
  //! How far apart two conformers must be to count as distinct.  Read by every
  //! search, though each applies it at a different stage: Thompson as its
  //! in-sweep novelty reward AND output dedup, systematic as per-node pool
  //! retention (finalRms/nodeRms derive from it), the tree as a beam prune.
  //! **-1 = AUTO** (autoDiversityRmsByRotor); 0 = no dedup.
  double diversityRmsThresh = 0.0;
  //! Fallback angles when the torsion sampler returns nothing for a junction.
  std::vector<double> defaultAngles = {-180.0, -120.0, -60.0, 0.0, 60.0, 120.0};
  unsigned int randomSeed = 0xf00d;  //!< RNG seed for the sampling searches
  double energyWindow = 25.0;        //!< retention window, kcal/mol
  //! >0: an explicit Thompson draw budget, overriding thompson.autoBudget.
  //! Also routes Auto: 0 with autoBudget off selects the tree.
  unsigned int thompsonBudget = 0;

  RigidRotorSearchMode searchMode = RigidRotorSearchMode::Auto;
  long timeBudgetMs = 0; //< Search Budget, 0 is no limit
  unsigned int autoSystematicMinRotors = 11; //!< switch to Systematic at this # rotors
  unsigned int rootSeeds = 6; //< max number of low fragments confs to seed
  unsigned int fragConfBranch = 4;
  std::shared_ptr<TorsionSampler> torsionSampler;
  JunctionBasinAngles junctionBasinAngles = JunctionBasinAngles::All;
  bool driveIntraFragmentTorsions = false; //< coarse sample non junction torsions

  TreeSearchParams tree;
  bool finalSymmetryDedup = true; //< XXX FIX ME->always do this
  bool disableIncrementalScoring = false; //< XXX FIX ME -> for testing FFs

  ThompsonParams thompson;
  SystematicParams systematic;
  DiagnosticsParams diagnostics;
};

}  // namespace RDKit

#endif
