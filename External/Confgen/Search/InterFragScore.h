//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//   Force-field scoring for scoring torsion driven molecules
//    E score is just bewteen fragments, fragments are assumed
//    to already be scored
#ifndef RDKIT_CONFGEN_INTERFRAGSCORE_H
#define RDKIT_CONFGEN_INTERFRAGSCORE_H

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>

#include <string>
#include <utility>
#include <vector>

#include "Search/RotorDriver.h"

namespace ForceFields {
namespace MMFF {
class InterFragVdWContrib;
class TorsionAngleContrib;
}  // namespace MMFF
}  // namespace ForceFields

namespace RDKit {
namespace MMFF {
class MMFFMolProperties;
}  // namespace MMFF

//! Contributions for the current FF only scoring between fragments.
struct RDKIT_FRAGMENTCONFGEN_EXPORT InterFragScoreHandles {
  ForceFields::MMFF::InterFragVdWContrib *vdw = nullptr;
  ForceFields::MMFF::TorsionAngleContrib *tor = nullptr;
};

//! Rescore only what a rotor actually changes.
/*!
  Allow changes between fixed and non-fixed rotors.

  I.e. A-r1-B-r2-C
   Allows r2 to be fixed OR r1 to be fixed and adjust the scoring
   appropriately.  IntraFragScore assumes r1 and r2 are changeable.

  This is used when incrementally building inside - out for instance.
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT IncrementalInterFragScore {
 public:
  IncrementalInterFragScore() = default;
  IncrementalInterFragScore(const RotorDriver &driver,
                            const InterFragScoreHandles &handles);

  //! False when incremental scoring cannot apply; rescore in full instead.
  bool isValid() const { return d_valid; }

  //! vdW energy over just the pairs rotor `r` changes.
  double changingVdw(const double *pos, unsigned int r) const;

  //! Junction-torsion energy (whole contrib; it is already junction-local).
  double torsionEnergy(const double *pos) const;

  size_t numChangingPairs(unsigned int r) const {
    return r < d_packed.size() ? d_packed[r].a1.size() : 0;
  }
  size_t numPairs() const;

 private:

  struct PackedRotor {
    std::vector<int> a1, a2;
    std::vector<double> Rs, wd;
  };
  std::vector<PackedRotor> d_packed;
  ForceFields::MMFF::InterFragVdWContrib *d_vdw = nullptr;
  ForceFields::MMFF::TorsionAngleContrib *d_tor = nullptr;
  double d_cut2 = 0.0;
  bool d_valid = false;
};

//! One MMFF torsion quartet spanning a junction bond, with its parameters.
/*!
  a1-j-k-a4, where j-k is the junction bond as the CALLER ordered it.  V1/V2/V3
  are the MMFF torsion constants; the energy is
  ForceFields::MMFF::Utils::calcTorsionEnergy(V1, V2, V3, cosPhi).
*/
struct RDKIT_FRAGMENTCONFGEN_EXPORT JunctionTorsionTerm {
  unsigned int a1, j, k, a4;
  double V1, V2, V3;
};

//! MMFF torsion quartets spanning each junction bond.
/*!
  Compute the torsion terms for a molecule
*/
RDKIT_FRAGMENTCONFGEN_EXPORT std::vector<std::vector<JunctionTorsionTerm>>
junctionTorsionTerms(
    const ROMol &mol, MMFF::MMFFMolProperties &props,
    const std::vector<std::pair<unsigned int, unsigned int>> &junctionBonds);

//! Build a WHOLE-MOLECULE force-field scorer.
/*!
  \param mol              Molecule to score
  \param electrostatics   if true turn electrostatics on
  \param ffVariant        selects the force field; see kFFVariants.  Only the
                          MMFF family is implemented -- anything else logs and
                          returns an empty scoring function.
  \param nonBondedThresh  ignore non-bonded pairs beyond this distance
*/
RDKIT_FRAGMENTCONFGEN_EXPORT RotorDriver::ScoreFn makeFullFFScoreFn(
    const ROMol &mol, bool electrostatics = false,
    const std::string &ffVariant = "MMFF94", double nonBondedThresh = 1.0e8);

//! Build an optimized FF to score between fragments
/*!
  \param mol  Molecule to score
  \param atomFragments fragment beloning to the atomIdx
  \param juntionBonds  bonds seperating the fragments
  \param electrostatics if true turn electrostatics on
  \param ffVariant      selects the force field; see kFFVariants.  Only the
                        MMFF family is implemented -- anything else logs and
                        returns an empty scoring function.
  \param vdwCuttof      If supplied, don't check VdW energies about this distance
  \param handles        optional vdw handles sent in by caller
 */
  
RDKIT_FRAGMENTCONFGEN_EXPORT RotorDriver::ScoreFn makeInterFragmentScoreFn(
    const ROMol &mol, const std::vector<int> &atomFragments,
    const std::vector<std::pair<unsigned int, unsigned int>> &junctionBonds,
    bool electrostatics = false, const std::string &ffVariant = "MMFF94",
    double vdwCutoff = 0.0, InterFragScoreHandles *handles = nullptr);

}  // namespace RDKit

#endif
