//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Plain parameter structs for the search strategies. 
//
#ifndef RDKIT_CONFGEN_SEARCHPARAMS_H
#define RDKIT_CONFGEN_SEARCHPARAMS_H

#include <RDGeneral/export.h>

#include <memory>
#include <vector>

#include "Sampler/TorsionSampler.h"   // TorsionSampler
#include <array>

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
//! How a rotor's influence is weighted.
/*!
  MovingAtoms scales by movingAtoms/maxMovingAtoms -- a geometric importance
  heuristic, not Bayesian evidence strength.

  XXX FIX ME:  the uniform posterior works the best, everything else either
               biases away from small rotor exploration or large.
	       everything else should be dropped

	       This, however, works well for the novelty weighting for
	       pruning rotors.
*/
//! What to do with a bond the input DECLARES as an atropisomer.
/*!
  We can either freeze it or do mild basin sampling ensuring we
  don't change stereochemistry.
*/
enum class AtropisomerSampling {
  Basin,  //!< drive it, but only within the well its declared sign names
  Frozen  //!< neither cut nor driven; keep the embedded geometry
};

enum class RotorWeighting {
  MovingAtoms,  //!< weight by movingAtoms/maxMovingAtoms; largest rotor = 1
  Uniform,      //!< every rotor counts equally
  Inverted      //!< minMovingAtoms/movingAtoms; SMALLEST rotor = 1
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT ThompsonParams {
  //! --- informed prior + budget scaling ---
  double priorStrength = 3.0;     //!< Beta alpha for sampler-preferred angles (>=1)
  double backstopStepDeg = Disabled;   //!< uniform-grid backstop step, deg (0 = none), alpha 1
  double sizePriorExp = 0.0;      //!< scale prior by (movingAtoms/max)^exp (0 = full prior)
  double noveltyAngleDeg = 0.0;  //!< in-sweep torsion-fingerprint novelty, deg (0 = coord RMSD)
  
  //! Weight on the torsion-fingerprint novelty distance.
  RotorWeighting noveltyWeighting = RotorWeighting::MovingAtoms;
  //! Exponent on that moving atom weight
  double noveltyWeightExp = 1.0;
  //! Rotor Weight on the Beta posterior update.
  /*!
    See FIX ME above
  */
  RotorWeighting posteriorWeighting = RotorWeighting::Uniform;

  bool autoBudget = true;         //!< auto-scale #draws with rotors + frag-conf arms
  unsigned int perRotor = 400;    //!< draws per rotor
  unsigned int perFragConf = 30;  //!< draws per fragment-conf arm
  unsigned int minBudget = 300;   //!< minimum # samples (min)
  unsigned int maxBudget = 64000; //!< max samples

  unsigned int maxConfs = Disabled;      //!< max number of output conformers (0 = nolimit)
  OutputSelection outMode = OutputSelection::Energy;  //!< selection when capped
  
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
  All,       //!< every junction gets the full angle set + small basin samples
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
  
  //! What Auto resolves to, indexed by ROTATABLE-BOND COUNT.
  /*!
    This is just a heuristic based on PDBBind and Platinum datasets
    so it may not translate to new datasets
  */
  static constexpr std::array<RigidRotorSearchMode, 14> AutoModeAtRotor = {
      RigidRotorSearchMode::Thompson,    // 0
      RigidRotorSearchMode::Thompson,    // 1
      RigidRotorSearchMode::Thompson,    // 2
      RigidRotorSearchMode::Thompson,    // 3
      RigidRotorSearchMode::Thompson,    // 4
      RigidRotorSearchMode::Thompson,    // 5
      RigidRotorSearchMode::Thompson,    // 6
      RigidRotorSearchMode::Thompson,    // 7
      RigidRotorSearchMode::Thompson,    // 8
      RigidRotorSearchMode::Thompson,    // 9
      RigidRotorSearchMode::Thompson,    // 10
      RigidRotorSearchMode::Systematic,  // 11
      RigidRotorSearchMode::Systematic,  // 12
      RigidRotorSearchMode::Thompson,    // 13 and every higher count
  };

  //! Table lookup, clamped: every count past the end takes the last entry.
  static constexpr RigidRotorSearchMode autoModeForRotors(size_t nRotors) {
    return AutoModeAtRotor[nRotors < AutoModeAtRotor.size()
                               ? nRotors
                               : AutoModeAtRotor.size() - 1];
  }

  unsigned int rootSeeds = 6; //< max number of low energy fragments confs to seed in search
  unsigned int fragConfBranch = 4;
  std::shared_ptr<TorsionSampler> torsionSampler;
  JunctionBasinAngles junctionBasinAngles = JunctionBasinAngles::All;
  //! How a DECLARED atropisomer axis is handled.  Never cut either way.
  AtropisomerSampling atropisomerSampling = AtropisomerSampling::Basin;

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
