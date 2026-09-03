//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#include "Joiner/FragmentJoiner.h"  // RigidRotorSearchMode (joinerSearchMode passthrough)
#ifndef RDKIT_FRAGMENT_CONFGEN_H
#define RDKIT_FRAGMENT_CONFGEN_H

#include "Sampler/TorsionSampler.h"
#include "Embedder/Fraglib.h"  // FraglibParams + FragmentEmbedMode/FragmentMinimize enums
#include <GraphMol/RDKitBase.h>
#include <memory>
#include <vector>

namespace RDKit {

//! version of getMolFrags that returns RWMols
inline RWMOL_SPTR_VECT getRWMolFrags(
    const ROMol &mol, bool sanitizeFrags = true,
    std::vector<int> *frags = nullptr,
    std::vector<std::vector<int>> *fragsMolAtomMapping = nullptr,
    bool copyConformers = true) {
  auto v = MolOps::getMolFrags(mol, sanitizeFrags, frags, fragsMolAtomMapping,
                               copyConformers);
  RWMOL_SPTR_VECT out;
  out.reserve(v.size());
  for (auto &m : v) {
    out.push_back(boost::static_pointer_cast<RWMol>(m));
  }
  return out;
}

//! Fragment-based conformer generator
/*!
  This conformer generated was implemented as a test bed for using different
  techniques for rotor-driven approaches.  This is a tried and true approach
  for many different conformation generators:

    References:
      https://pubs.acs.org/doi/10.1021/ci100031x (omega)
      https://pmc.ncbi.nlm.nih.gov/articles/PMC2896087/ (frog2)
      https://pubs.acs.org/doi/10.1021/acs.jcim.3c00563 (conforge)

   The basic algorithm used here is as follows:

    Note: MMFF WITHOUT electrostatics is preferred for finding biologically
        relevant conformations in a protein context.  With electrostatics the
        field favours gas-phase minima -- intramolecular hydrogen bonds and
        folded, charge-paired geometries that a solvated or bound ligand does
        not adopt.

    1. SPLIT: perceive fragment-link (rotatable) bonds and split the molecule into its
       maximal rigid fragments,
    2. EMBED: generate 3D geometry for each fragment with RDKit distance geometry
       (plain DG) or ETKDG.
    3. JOIN: re-assemble the fragments using preffered torsions
    4. SAMPLE+SEARCH: search the rotor space between fragments (and optionally sample intra fragment
        rotors as well)
    4. score, energy-window filter and RMSD-diversity select the ensemble.

    SAMPLE+TorsionLibrary:

    (A) ETKDG torsion sampling (default)

    (B) Uniform torsion sampling with various angles.  This is also used as a back-stop
        for missing torsions

    (C) various torsion libaries can be used:
        Note: these libraries are not shipped with the RDKit due to licensing (non free for
	  commercial use)
	  
     References:
     1. Schärfer, C., Schulz-Gasch, T., Ehrlich, H.-C., Guba, W.,
        Rarey, M., & Stahl, M. (2013). Torsion Angle Preferences in
        Druglike Chemical Space: A Comprehensive Guide. Journal of
        Medicinal Chemistry, 56(5), 2016–2028.
        https://doi.org/10.1021/jm3016816
     2. Guba, W., Meyder, A., Rarey, M., & Hert, J. (2016).
        Torsion Library Reloaded: A New Version of Expert-Derived
        SMARTS Rules for Assessing Conformations of Small Molecules.
        Journal of Chemical Information and Modeling, 56(1), 1–5.
        https://doi.org/10.1021/acs.jcim.5b00522
     3. Penner, P., Guba, W., Schmidt, R., Meyder, A., Stahl, M.,
        & Rarey, M. (2022). The Torsion Library: Semi-automated
        Improvement of Torsion Rules with SMARTScompare. Journal of
        Chemical Information and Modeling, 62(7), 1644–1653.
        https://doi-org/10.1021/acs.jcim.2c00043

  SEARCH: for experimental purposes, multiple search mechanisms are included

  At each step MMFF94 or MMFF94s is used for ranking candidates.  Default is MMFF94s

*/

//! Which score the output ensemble is ranked / capped by.
enum class OutputRanking {
  FFEnergy,     //!< full single-point (or basin) MMFF energy label (default)
  InterFragScore  //!< junction-only inter-fragment term (vdW + torsion + stretchbend, etc; no
                  //!< intra strain)
};

//! The molecule's rotatable bonds
//!   inter are between fragments
//!   intra are internal to fragments (only used when wholeAcyclicFragments is true)
struct RDKIT_FRAGMENTCONFGEN_EXPORT RotatableBonds {
  std::vector<unsigned int> inter;  
  std::vector<unsigned int> intra;
};

struct RDKIT_FRAGMENTCONFGEN_EXPORT FragmentConfGenParams {
  int randomSeed = -1;

  //! Number of conformers to output, lowest energy first.
  //!  **AUTO** -- scale with the molecule's rotatable-bond count, floppier
  //!     means more confs generated
  //!
  int numOutputConfs = kAutoI;

  //! energy window (kcal/mol) above the best conformer to keep
  double energyWindow = 10.0;
  //! Maximum number of rotatable bonds to embed.   0 == no limit
  int maxRotatableBonds = 0;
  //! Wall-clock budget per molecule, milliseconds.  0 = no limit.
  //!  when hit, return current ensemble
  long timeBudgetMs = 0;

  //! sample "trivial" symmetric-top rotors too (e.g. -CF3, -C(CH3)3).
  bool sampleTrivialRotors = false;

  //! cut only rotatable bonds adjacent to a ring
  //!  See RotatableBonds
  bool wholeAcyclicFragments = false;

  FraglibParams embedding;
  FragmentJoinerParams joiner;
  RigidRotorSearchParams search;

  OutputRanking outputRanking = OutputRanking::FFEnergy;
  //! Use a fill minimization at the end to rank
  bool rankByBasinEnergy = false;
  //! Forcefield to use when ranking
  //!  XXX FIX ME -> should always be the same as the embedding
  std::string labelFFVariant;

  //! Shared thread-safe Fragment Library.  Prefer setFraglib(), which
  //! validates compatibility at the point of assignment.
  std::shared_ptr<Fraglib> fraglib;

  //! Set a shared fraglib
  /*!
    This can be useful for multi-threaded conf building for speed as
    new fragments are embedded and shared
    
    Throws std::invalid_argument if the fragib params are invalid
  */
  void setFraglib(std::shared_ptr<Fraglib> lib);

  //! diagnostics toggle
  //!  XXX FIX ME -> probably move into the result object.  The joiner profiler
  //!  it drives is a process-global flag, so under the shared-fraglib /
  //!  many-generators threading model the last generator constructed wins.
  DiagnosticsParams diagnostics;

  //! Current Defaults
  FragmentConfGenParams() {
    search.torsionSampler = std::make_shared<ETKDGTorsionSampler>();
    search.rootSeeds = 6;
    search.fragConfBranch = 8;
    search.tree.beamWidth = 80;
    // -1 AUTO: 0.5 normally, 0.75 for floppy molecules whose rotors are
    // concentrated in one long chain (see Utils/RotorTopology.h).
    search.diversityRmsThresh = kAutoD;
    search.thompson.backstopStepDeg = 60.0;
    search.thompson.noveltyAngleDeg = 30.0;
  }

  //! Validate parameters.  returns empty string on success
  //!  otherwise a human readable failure report.
  std::string validate() const;
};

// XXX FIX ME -> Eventually we will propagate more errors and diagnostics
//! Properties set on each output conformer.
/*!
  The energy and the force field that produced it are SEPARATE fields: the
  energy is only interpretable against the force field it came from, and
  burying the force field in the property NAME (as the old "MMFF_energy" did)
  makes it unreadable to a consumer and a lie the moment a non-MMFF field is
  added.
*/
//! Force-field energy of the conformer, kcal/mol (double).
inline constexpr const char *kEnergyProp = "energy";
//! Which force field kEnergyProp was computed with, e.g. "MMFF94s" (string).
inline constexpr const char *kForceFieldProp = "forcefield";
//! Inter-fragment score from the rigid-rotor search (double).
inline constexpr const char *kInterFragScoreProp = "interfrag_score";

enum class RDKIT_FRAGMENTCONFGEN_EXPORT FragConfGenResultType {
  OK,
  FF_FAIL,
  //! rejected before any work: more rotatable bonds than maxRotatableBonds
  //! allows
  TOO_MANY_ROTORS,
  //! the time budget expired; `conformers` holds the best ensemble found so
  //! far, which may be usable but is NOT the converged answer
  TIMED_OUT
};

//! Result of a conformer-generation run
struct RDKIT_FRAGMENTCONFGEN_EXPORT FragmentConfGenResult {
  std::vector<ROMOL_SPTR> conformers;

  //! number of fragment-link (rotatable) bonds that were cut
  unsigned int numLinkBonds = 0;
  //! number of rigid fragments the molecule was split into
  unsigned int numFragments = 0;
  //! lowest MMFF energy in the output ensemble (NaN if empty)
  double bestEnergy = std::numeric_limits<double>::quiet_NaN();
  //! number of Thompson angle draws the joiner spent
  //! 0 for the exhaustive-beam path or non-assembled)
  unsigned int joinerBudget = 0;
  //! wall time (ns) the rigid-rotor search spent (assembly/search only; embed
  //! excluded)
  long long joinerAssemblyNs = 0;

  //! set when the molecule was rejected or truncated by a guard rail; OK
  //! otherwise
  FragConfGenResultType status = FragConfGenResultType::OK;

  //! Returns true if the result is valid, false otherwise
  bool isValid() const { return resultType() == FragConfGenResultType::OK; }
  FragConfGenResultType resultType() const {
    if (status != FragConfGenResultType::OK) {
      return status;
    }
    if (conformers.size() == 0) {
      return FragConfGenResultType::FF_FAIL;
    }
    return FragConfGenResultType::OK;
  }
};

class RDKIT_FRAGMENTCONFGEN_EXPORT FragmentConfGen {
 public:
  explicit FragmentConfGen(FragmentConfGenParams params = FragmentConfGenParams());

  const FragmentConfGenParams &params() const { return d_params; }

  //! Build the conformer ensembles.  note: Hs are added if they don't exist.
  //!  \param mol   Molecule to generate conformers for
  //! \param linkBonds  if non-null, cut ONLY these bonds (COARSE assembly).
  //!        Everything between them stays whole and its rotors become
  //!        intra-fragment ones.  Used by callers that already know the seams,
  //!        e.g. a synthon library assembling at its own junctions.
  FragmentConfGenResult build(
      const ROMol &mol, const std::vector<unsigned int> *linkBonds = nullptr) const;

  //! Find rotatable bonds, see RotatableBonds class for more details
  //!  \param sampleTrivial           sample trivial rotors (CF3)
  //!  \param wholeAcyclicFragments   return whole linkers between rings
  static RotatableBonds findRotatableBonds(const ROMol &mol,
                                           bool sampleTrivial = false,
                                           bool wholeAcyclicFragments = false);

  //! Return interfrags only
  static std::vector<unsigned int> findLinkBonds(
      const ROMol &mol, bool sampleTrivial = false,
      bool wholeAcyclicFragments = false) {
    return findRotatableBonds(mol, sampleTrivial, wholeAcyclicFragments).inter;
  }

  //! Cut the molecule at its link bonds and embed each rigid fragment in 3D.
  /*!
    Each cut inserts a matched pair of isotope-labelled exit vectors
    so the fragments can be re-zipped by MolzipLabel::Isotope.

    \param mol       the molecule to fragment (added Hs are handled internally)
    \param linkBonds if non-null, receives the link-bond indices that were cut
    \return the embedded fragments; a single element (the whole molecule,
            embedded) when there are no link bonds.
  */
  //! \param linkBondsOut  receives the bonds that were cut
  //! \param linkBonds     if non-null, cut ONLY these (empty = cut nothing, so
  //!        the molecule is embedded and cached WHOLE).  Lets a caller prefill
  //!        the pieces a COARSE assembly will ask for.
  std::vector<ROMOL_SPTR> fragmentAndEmbed(
      const ROMol &mol, std::vector<unsigned int> *linkBondsOut = nullptr,
      const std::vector<unsigned int> *linkBonds = nullptr) const;

 private:
  //! Run the full single-molecule pipeline
  void buildEnsemble(const ROMol &mol, FragmentConfGenResult &result,
                     const std::vector<unsigned int> *linkBonds) const;

  FragmentConfGenParams d_params;
};

//! Given a set of parameters, build a FragmentConfGen
RDKIT_FRAGMENTCONFGEN_EXPORT FragmentConfGen
createFragmentConfGen(FragmentConfGenParams params = FragmentConfGenParams());

}  // namespace RDKit

#endif
