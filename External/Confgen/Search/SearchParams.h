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
#include "Utils/ParamSentinels.h"     // AutoR / Disabled / resolveAuto

namespace RDKit {

//! Parameters for tree search, sometimes known as beam search
struct RDKIT_FRAGMENTCONFGEN_EXPORT TreeSearchParams {
  unsigned int beamWidth = 50;  //!< bucketed conformers kept per level
  double angleTolerance = 15.0;
};

//! Thompson Search limit for output conformers.
enum class OutputSelection {
  Energy,   //!< lowest energy first
  Diverse,  //!< greedy RMSD-diverse over the energy-sorted pool
  Stratify  //!< energy-stratified stride across the pool
};

//! ThompsonParams, so many parameters!
struct RDKIT_FRAGMENTCONFGEN_EXPORT ThompsonParams {
  //! --- informed prior + budget scaling ---
  double priorStrength = 3.0;     //!< Beta alpha for sampler-preferred angles (>=1)
  double backstopStepDeg = Disabled;   //!< uniform-grid backstop step, deg (0 = none), alpha 1
  double sizePriorExp = 0.0;      //!< scale prior by (movingAtoms/max)^exp (0 = full prior)
  double noveltyAngleDeg = 0.0;  //!< in-sweep torsion-fingerprint novelty, deg (0 = coord RMSD)
  
  bool autoBudget = true;         //!< auto-scale #draws with rotors + frag-conf arms
  unsigned int perRotor = 400;    //!< draws per rotor
  unsigned int perFragConf = 30;  //!< draws per fragment-conf arm
  unsigned int minBudget = 300;   //!< minimum # samples (min)
  unsigned int maxBudget = 64000; //!< max samples

  unsigned int maxConfs = Disabled;      //!< max number of output conformers (0 = nolimit)
  OutputSelection outMode = OutputSelection::Energy;  //!< selection when capped
  bool flatContext = false;  //!< collapse the tree-descent context
  
  unsigned int refineSteps = Disabled; //! Sample K lower basins for better energyies
  double refineStepDeg = 8.0;     //!< initial coordinate-descent step
  unsigned int refinePasses = 3;  //!< step halves each pass
  
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT SystematicParams {
  double eWindow = 15.0;  //!< per-node energy window (kcal)
  long maxPoolConfs = 4000; //!< This is per pool, see combinePoolMultiplier as well
  int maxFragmentConfs = 16; //!< leaf fragment-conf limit

  double finalRms = AutoR; //!< finalRms thresh, AutoR is auto selectin
  double nodeRms = AutoR;  //!< per node pruning
  bool energyOnly = false;  //!< energy-window retention only (no RMSD-diverse)

  double upperEnergyWindow = 5.0; //<! drop early candidates prior to stabilizing MMFF basin

  double dedupEnergyBand = Disabled; //!< Only RMS check within this energy window
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

//! Does the current mode add basin angles to the rotors for sampling
//!  Basin angles are small rotor additions to enhance sampling
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
  Auto,       //!< Thompson when a thompson rotor budget is specified, else the Tree beam
  Tree,       //!< Greedy deterministic tree/beam using per-edge minimum energies
  Thompson,   //!< one-armed-bandit sampling (maybe good for non-separable force
              //!< fields)
  Systematic, //!< bottom-up systematic combine + prune
  Merged      //!< experimental: run Systematic + Thompson, merge ensembles (expensive)
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT RigidRotorSearchParams {
  //! AUTO diversity radius per ROTOR COUNT: index by the molecule's rotatable
  //! bond count, clamped to the last entry.  Only consulted when
  //! diversityRmsThresh is AUTO (-1).
  /*!
    At a rotor count of 9, increase the RMSD so we can find more
    variety in conformers.
    A floppy molecule's conformer pool is dominated by near-duplicates, so a
    tight "these are the same pose" radius fills a fixed output budget with one
    basin's worth of variations.

    note: not applied to the TreeSearch as higher rms paradoxically explores
    fewer basins.
  */
  std::vector<double> autoDiversityRmsByRotor = {0.5, 0.5, 0.5, 0.5, 0.5,
                                                 0.5, 0.5, 0.5, 0.75};

  //! set the diversityRmsThreshold, 0 = no dedup, -1 = Auto pick RMS
  double diversityRmsThresh = 0.0;
  //! Fallback angles when the torsion sampler returns nothing for a junction.
  std::vector<double> defaultAngles = {-180.0, -120.0, -60.0, 0.0, 60.0, 120.0};
  unsigned int randomSeed = 0xf00d;  //!< RNG seed for the sampling searches
  double energyWindow = 25.0;        //!< energy retention window, kcal/mol
  
  unsigned int thompsonBudget = Disabled; //!< Thompson draw budget, >0 explicit override

  RigidRotorSearchMode searchMode = RigidRotorSearchMode::Auto;
  long timeBudgetMs = Disabled; //< Search time budget, 0 is no limit
  unsigned int autoSystematicMinRotors = 11; //!< switch to Systematic at this # rotors
  unsigned int rootSeeds = 6; //< max number of low energy fragments confs to seed in search
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
