//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//   This is the core data structure for the conformer generator.
//   It holds the fragments from the embedder, puts them into
//    an initial geometric position and sets up the core for running
//    rotor driving and the various search methodologies
#ifndef RDKIT_FRAGMENTJOINER_H
#define RDKIT_FRAGMENTJOINER_H

#include <RDGeneral/export.h>
#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <Geometry/point.h>
#include <Geometry/Transform3D.h>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "Search/InterFragScore.h"
#include "Search/SearchResult.h"
#include "Search/SearchParams.h"
#include "Sampler/TorsionSampler.h"
#include "Utils/DiagnosticsParams.h"

namespace RDKit {

//! Fine-grained control over how torsions are sampled
enum class IntraRotorType {
  Free,         //!< ordinary rotor (wholeAcyclicFragments)
  PlanarAmide,  //!< conjugated C(=X)-Y: a real torsion, but 0/180 only
  Atropisomer   //!< declared axis: sample only WITHIN its own well, never
                //!< across the barrier -- the far well is the other enantiomer
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT IntraRotor {
  unsigned int bond = 0;
  IntraRotorType type = IntraRotorType::Free;
};


class Embedder;

//! Joins rigid fragments at exit vectors
/*!
  Algorithm:
     pick a parent - the largest fragment
     Everything else is considered a child, rigidly transform and place
     at appropriate bond angles and base torsion.

  Specifics: bonding atoms are called childBondAtoms, these are moved and placed
  along the fragment exit vector at the ideal bond length for the atoms types.

  exit vectors are solely determined by the embedding and are usually dummy
  atoms that are treated as carbon for the purposes of coordinate generation.

  These might not make ideal bond lengths, these are recalculated during the
  joining.

  \param childCoords  the child fragment's atoms (all of them), transformed in
                      place
  \param childBondAtom index into childCoords of the atom that bonds to the
                       parent
  \param childExit     index into childCoords of the child's exit marker
  \param parentBondAtom world position of the parent atom the bond forms to
  \param parentExit    world position of the parent's exit marker (direction
                       only)
  \param bondLen       MMFF-ideal length for the new junction bond
*/
RDKIT_FRAGMENTCONFGEN_EXPORT void placeChildCoords(
    std::vector<RDGeom::Point3D> &childCoords, unsigned int childBondAtom,
    unsigned int childExit, const RDGeom::Point3D &parentBondAtom,
    const RDGeom::Point3D &parentExit, double bondLen);

//! The rigid transform and associated data that statisfies the placeChildCoords
//! algorithm
//!  the origin of the rotation/translation is childNbr
struct ChildPlacement {
  RDGeom::Transform3D rot;  //!< identity unless a rotation is needed
  RDGeom::Point3D
      childNbr;  //!< subtract before rotating (child's bonding-atom position)
  RDGeom::Point3D site;  //!< add after rotating (the new junction site)

  RDGeom::Point3D apply(const RDGeom::Point3D &p) const {
    RDGeom::Point3D q = p - childNbr;
    rot.TransformPoint(q);
    return q + site;
  }
};

//! Compute the placement transform without touching any coordinates.
RDKIT_FRAGMENTCONFGEN_EXPORT ChildPlacement computeChildPlacement(
    const RDGeom::Point3D &childBondPos, const RDGeom::Point3D &childExitPos,
    const RDGeom::Point3D &parentBondAtom, const RDGeom::Point3D &parentExit,
    double bondLen);

//! One conformer of a single fragment awaiting joining.
//! Energy is whatever the fragment conf gen returned and can be used
//!  to optimize the rotor driving, i.e. calculated only once
struct RDKIT_FRAGMENTCONFGEN_EXPORT JoinFragmentConf {
  std::vector<RDGeom::Point3D> pos;
  double energy = 0.0;
};

//! Structure holding the indices of the fragment IN the parent mol
//!  conformers are kept lowest energy first
struct RDKIT_FRAGMENTCONFGEN_EXPORT JoinFragment {
  std::vector<unsigned int>
      atoms;  //!< conformers atom indices int the assembly
  std::vector<JoinFragmentConf> confs;  //!< sortest lowest to highest energy
};

//! An (undirected) inter-fragment rotatable bond.
//!  The joiner must orient the bond correctly and move the child's coords
struct RDKIT_FRAGMENTCONFGEN_EXPORT JoinJunction {
  unsigned int fragA, fragB;  //!< the two fragments this bond joins
  unsigned int atomA,
      atomB;  //!< the bond's atoms (atomA in fragA, atomB in fragB)
  double bondLen =
      1.5;  //!< MMFF-ideal length for the junction bond (default C-C)
};

//! Defines a rotatable bond between two fragments.
//! This edge is used when searching and rotating.
struct RDKIT_FRAGMENTCONFGEN_EXPORT JoinEdge {
  unsigned int parentFrag, childFrag, parentAtom, childAtom;
  double bondLen;
  unsigned int bondIdx;  //!< the rotatable bond's index in the mol
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT FragmentJoinerParams {
  //! Fragmentation flags forwarded to findLinkBonds for consistency
  bool sampleTrivialRotors = false;
  bool wholeAcyclicFragments = false;
  std::string ffVariant = "MMFF94";
  //! Optimization, don't bother with vdw calcs above this distance
  double interFragVdwCutoff = 10.0;
  //! Score by shrugging of conf energy but keeping the same basic geometry
  //!  this is useful when using non MMFF generated confs
  double fragShrugDisplacement = 0.0;
  double fragShrugForceConst = 100.0;
  //! When scoring, use the full MMFF energy, not the faster calculations
  //!  assuming constant fragments
  bool useFullFFScorer = false;
  //! Profiling diagnostics
  DiagnosticsParams diagnostics;

  //! Are these parameters viable?
  /*!
    \return empty when they are, otherwise a human readable reason.  Checked by
            joinFragments() (which reports BadParams) and by
            FragmentConfGenParams::validate().
  */
  std::string validate() const;

  //! True when validate() finds nothing wrong.
  bool isValid() const { return validate().empty(); }
};

//! Joiner Status
enum class FragmentJoinerStatus {
  Ok,
  BadParams,             //!< the FragmentJoinerParams are not viable
  NoFragments,           //!< nothing to join
  NoFragmentConformers,  //!< a fragment arrived with an empty conformer pool (embedding fail)
  NoScorer               //!< no force field could be built for the molecule
};

//! Human readable form of a FragmentJoinerStatus
RDKIT_FRAGMENTCONFGEN_EXPORT const char *fragmentJoinerStatusMessage(
    FragmentJoinerStatus s);

//! The search context used for rigid rotor searching, used by all search
//!  algorithms.
struct RDKIT_FRAGMENTCONFGEN_EXPORT FragmentJoinerContext {
  ROMol mol;  //!< the joined molecule, carrying the initial placement
  std::vector<JoinFragment> frags;
  std::vector<JoinEdge> edges;  //!< BFS sorted from parent out to the children
  std::vector<unsigned int> rotorBonds;  //!< bond idx per edge (INTER-fragment junctions)
  //! Left over rotor bonds internal to fragments.  Each has a designation
  //!  to help refine sampling
  std::vector<IntraRotor> intraRotorBonds;

  //! The uncut rotors a search may actually drive, with every gate applied.
  /*!
    One call, so a search never switches on IntraRotorType.  Adding a new kind
    changes this function and nothing else.

      Free         driven when driveIntraFragmentTorsions
      PlanarAmide  driven when driveIntraFragmentTorsions (a real 0/180 torsion)
      Atropisomer  driven when atropisomerSampling == Basin, and then the caller
                   must keep it inside its own well -- see basinLimitAngles()
  */
  std::vector<unsigned int> getIntraRotorBonds(
      const RigidRotorSearchParams &sp) const;

  //! The rule for one bond, or nullopt when it is not an uncut rotor.
  std::optional<IntraRotorType> getIntraType(unsigned int bond) const;
  unsigned int root = 0;
  RotorDriver::ScoreFn scorer;
  //! raw contrib handles into `scorer`'s force field for INCREMENTAL scoring
  InterFragScoreHandles scoreHandles;
  std::string ffVariant;

  FragmentJoinerStatus status = FragmentJoinerStatus::Ok;

  //! True when this context can be searched.
  bool isValid() const { return status == FragmentJoinerStatus::Ok; }

  //! Generate initial placements
  std::vector<double> placeAll(
      const std::vector<unsigned int> &confChoice) const;

  //! symettrically RMSD prune
  void symmetryDedup(std::vector<SearchResult> &v, double thr) const;
};

//! Coordinate-only greedy fragment join.
/*!
  Roots the fragment tree at the largest fragment, orients the junctions into
  parent->child edges by BFS, places every fragment at its lowest-energy
  conformer to seed the molecule, and builds the inter-fragment scorer.

  \param mol             the assembled molecule
  \param atomFragment          fragment index per atom
  \param fragments       the rigid fragments and their conformer pools
  \param junctions       the inter-fragment rotatable bonds
  \param params          joiner parameters
  \param intraRotorBonds rotatable bonds NOT cut (they live inside a fragment)

  return the search context; check isValid() -- failures are reported through
          FragmentJoinerContext::status, not thrown.
*/
RDKIT_FRAGMENTCONFGEN_EXPORT FragmentJoinerContext joinFragments(
    const ROMol &mol, std::vector<int> atomFragment,
    std::vector<JoinFragment> fragments, std::vector<JoinJunction> junctions,
    FragmentJoinerParams params = {},
    std::vector<IntraRotor> intraRotorBonds = {});

//! Prepared input for the FragmentJoiner (full topology + fragment pools +
//! junctions), all in full-molecule atom indices.
struct RDKIT_FRAGMENTCONFGEN_EXPORT FragmentJoinerInput {
  ROMOL_SPTR mol;
  std::vector<int> atomFragment;
  std::vector<JoinFragment> fragments;
  std::vector<JoinJunction> junctions;
  //! Rotatable bonds NOT cut (they live inside a fragment) -- see
  //! FragmentJoinerContext.
  std::vector<IntraRotor> intraRotorBonds;
};

//! Join a prepared FragmentJoinerInput.  See joinFragyments() above.
RDKIT_FRAGMENTCONFGEN_EXPORT FragmentJoinerContext joinFragments(
    FragmentJoinerInput in, FragmentJoinerParams params = {});

//! Build the joiner input using existing params
//!  see the individual parameter details for descriptions
//! \param linkBondsOverride  if non-null, cut ONLY these bonds instead of every
//!        inter-fragment rotatable bond.  Rotatable bonds not in the set become
//!        intra-fragment rotors.  Lets a caller that knows the real seams (the
//!        synthon junctions of a zipped product, say) keep the pieces whole.
RDKIT_FRAGMENTCONFGEN_EXPORT FragmentJoinerInput buildFragmentJoinerInput(
    const ROMol &mol, unsigned int nConfs = 16, int seed = 0xf00d,
    const std::string &ffVariant = "MMFF94", bool fragUseDG = false,
    bool fragMinimizeFF = true, const Embedder *lib = nullptr,
    const FragmentJoinerParams *asmParams = nullptr,
    const std::vector<unsigned int> *linkBondsOverride = nullptr);

//! Print (to stderr) and reset the accumulated joiner diagnostics
RDKIT_FRAGMENTCONFGEN_EXPORT void printJoinerProfile(const char *label = "",
                                                     size_t nMols = 0);

//! Turn the joiner profiler on/off
RDKIT_FRAGMENTCONFGEN_EXPORT void setJoinerProfiling(bool on);

//! Tracks the joiner only timings in the profile
//! Requires profiling enabled.
RDKIT_FRAGMENTCONFGEN_EXPORT double takeJoinerWarmMsPerMol(size_t nMols);

}  // namespace RDKit

#endif
