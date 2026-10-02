//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//   A combinatorial synthon library that hands back ASSEMBLED 3D products.
//
//   One reaction with `arity` positions and a reagent list per position; a
//   product is chosen by one reagent index per position.  The reaction, the
//   reagent sets and their serialization come from EnumerateLibrary; what is
//   added here is 3D: a fragment cache for this reaction's synthons and
//   conformer generation for the assembled product.
//
#ifndef RDKIT_CONFGEN_ENUMERATESYNTHONS3D_H
#define RDKIT_CONFGEN_ENUMERATESYNTHONS3D_H

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Enumerate/Enumerate.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <RDGeneral/RDLog.h>

#include <iosfwd>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "Embedder/Embedder.h"
#include "FragmentConfGen.h"

namespace RDKit {

//! Status of a Synthon build
enum class SynthonBuildStatus {
  Ok,
  BadReagentIndex,  
  ZipFailed,        
  NoReagentConfs,   
  ConfGenFailed     
};

RDKIT_FRAGMENTCONFGEN_EXPORT const char *synthonBuildStatusMessage(
    SynthonBuildStatus s);

struct RDKIT_FRAGMENTCONFGEN_EXPORT SynthonProduct {
  ROMOL_SPTR mol;  //!< null unless status == Ok; carries its conformers
  SynthonBuildStatus status = SynthonBuildStatus::Ok;
  bool usedCoarseFallback = false; // fell back to a full conformer build
  explicit operator bool() const {
    return status == SynthonBuildStatus::Ok && mol;
  }
};

//! Conformers per assembled product for a shape search.
/*!
  This controls the default output confs, generally for shape scoring.
  More confs is generally better but requires more computations
*/
constexpr int kDefaultSynthonProductConfs = 10;

//! Controls how the conformers are made for searching
/*!
  | style | fragment geometry from |
  |-------|------------------------|
  | `Full`          | RigidConfGenDefaults                      |
  | `Coarse`        | ETKDG on the CAPPED synthon, in situ      |
  | `CoarseSampled` | sampled from FULL UNCAPPED molecules      |

  The difference between CoarseSampled and Coarse is that Coarse
  samples fragments in the context of final products.  This is intended
  to sample better rotations, vdw with other side chains etc.
  However, it makes the generation of the embedded library much more
  expensive.
*/
enum class SynthonEmbedStyle {
  Full,          // Use RigidConfGenDefaults
  Coarse,        // Use ETKDG to make the reactants in-situ and rotor drive them
  CoarseSampled  // As Coarse, but fraglib geometries sampled from full molecules
};

//! returns true if we are using a coarse embedding method
inline bool isCoarseAssembly(SynthonEmbedStyle s) {
  return s == SynthonEmbedStyle::Coarse ||
         s == SynthonEmbedStyle::CoarseSampled;
}

//! Returns a name from a SynthonEmbedStyle
RDKIT_FRAGMENTCONFGEN_EXPORT const char *synthonEmbedStyleName(
    SynthonEmbedStyle s);
//! \return false if the given name is not a style.
RDKIT_FRAGMENTCONFGEN_EXPORT bool synthonEmbedStyleFromName(
    const std::string &name, SynthonEmbedStyle &out);

struct RDKIT_FRAGMENTCONFGEN_EXPORT EnumerateSynthons3DParams {
  FragmentConfGenParams confgen;

  SynthonEmbedStyle embedStyle = SynthonEmbedStyle::Full;

  bool prefillEmbedder = true; // pregen the synthon confs
  bool storeEmbedder = true;   // store the synthons

  EnumerateSynthons3DParams() {
    confgen.numOutputConfs = kDefaultSynthonProductConfs;
  }
};

//! A one-reaction combinatorial 3D synthon library.
/*!
  EnumerateSynthons3D thread safety:
    XXX FIX ME -> need to be able to pass in a embedder or clone an enumerate library

*/
class RDKIT_FRAGMENTCONFGEN_EXPORT EnumerateSynthons3D
    : public EnumerateLibrary {
 public:
  EnumerateSynthons3D() : EnumerateLibrary() {}
  explicit EnumerateSynthons3D(const std::string &s) : EnumerateLibrary() {
    initFromString(s);
  }

  //! n.b. reagents must be synthon compatiable
  explicit EnumerateSynthons3D(const EnumerationTypes::BBS &reagents,
                               EnumerateSynthons3DParams params = {});

  //! Is this synthon set viable
  bool isValid() const { return d_valid; }

  //! Number of reagent sets in the reaction, i.e. R1, R2, R3
  unsigned int arity() const {
    return static_cast<unsigned int>(getReagents().size());
  }
  
  unsigned int numReagents(unsigned int pos) const {
    return static_cast<unsigned int>(getReagents()[pos].size());
  }

  //! Get the molecule with no conformers
  ROMOL_SPTR get2D(const std::vector<unsigned int> &reagentIdx) const;

  //! get 3D conformations of the products
  std::vector<MOL_SPTR_VECT> get(
      const EnumerationTypes::RGROUPS &pos) const override;

  //! Get a single product with status
  SynthonProduct getProduct(const std::vector<unsigned int> &reagentIdx) const;

  //! Fragment heavy atom count
  unsigned int fragmentHeavyCount(unsigned int pos, unsigned int idx) const;

  //! Product heavy atom count
  unsigned int productHeavyCount(const std::vector<unsigned int> &reagentIdx) const;

  //! Smallest and largest product this library can make (heavy atoms).
  std::pair<unsigned int, unsigned int> productSizeRange() const;

  //! Smallest and largest synthon at one position (heavy atoms) ignoring
  //!  unembedable synthons
  std::pair<unsigned int, unsigned int> positionSizeRange(
      unsigned int pos) const;

  //! Bonds of one synthon that should be cut once the product is assembled.
  /*!
    n.b. after the product is formed, the cut bonds need to be remapped.
    if the cut bonds are empty, proceed as normal using the zip-junction
    bond to rotate.

    here "cut" really means how to split up the molecule for embedding and
    rotor driving.

    here is what happens when a synthon's cut bonds are set or empty
    
    synthon A   synthon B   A–B junction   A's set bonds   B's set bonds
    empty       empty       cut (standard) —               -
    set         empty       not cut        cut             —
    empty       set         not cut        —               cut
    set         set         not cut        cut             cut

    This functionality is required for ring-forming synthons.  The embedder
    cannot embed a partial ring, however it can embed the spinach.  This
    allows the ring to be formed and have a good chunk of the synthon
    pre-embedded.

    Note: this technique generally only works during ring formation
      as when noticed when generating synthon libraries.
  */
  const std::vector<unsigned int> &synthonCutBonds(unsigned int pos,
                                                   unsigned int idx) const;

  //! Set the cut bonds for one synthon (genSynthonLib does this at build).
  void setSynthonCutBonds(unsigned int pos, unsigned int idx,
                          std::vector<unsigned int> bonds);

  //! True if ANY synthon carries a cut-bond instruction.
  bool hasSynthonCutBonds() const { return d_haveCutBonds; }

  //! Embed all the fragments.  (n.b. thread-safe)
  /*!
    \param numThreads workers; 0 = hardware concurrency, 1 = serial.
  */
  unsigned int prefill(unsigned int numThreads = 1);

  const std::shared_ptr<Embedder> &embedder() const { return d_embedder; }

  //! Replace a new embedder
  void setEmbedder(std::shared_ptr<Embedder> lib) {
    d_unusable.reset();  // memoised from the OLD cache; no longer valid
    d_embedder = lib;
    d_params.confgen.embedder = lib;
  }
  const EnumerateSynthons3DParams &params3D() const { return d_params; }

  //! Exit-atom conventions this library was built with.
  const MolzipParams &molzipParams() const { return d_molzipParams; }

  //! Checks the synthon for usability
  bool synthonUnusable(unsigned int pos, unsigned int idx) const;

  //! Convert the synthon into the embedded fragment
  /*!
    \param synthon - the synthon to be embedded
    \returns the fragment that will be embedded.
    
    Note, if this returns nullptr, then the synthon CANNOT be
    embedded as-is.
  */
  std::unique_ptr<RWMol> cacheFragment(const ROMol &synthon) const;

  void toStream(std::ostream &ss) const override;
  void initFromStream(std::istream &ss) override;

 private:
  // helper to silence warnings about bad reactions
  struct LogBlocker {
    RDLog::LogStateSetter blocker;
  };
  EnumerateSynthons3D(LogBlocker, const EnumerationTypes::BBS &reagents,
                      EnumerateSynthons3DParams params);

  EnumerateSynthons3DParams d_params;
  MolzipParams d_molzipParams;
  bool d_valid = false;
  std::shared_ptr<Embedder> d_embedder;
  //! Records usable synthons for fast acceess
  struct UnusableCache {
    std::mutex mutex;
    std::vector<std::vector<signed char>> flags;
  };
  mutable std::shared_ptr<UnusableCache> d_unusable;
  //! heavy atoms per synthon, exit dummies excluded; [position][index]
  std::vector<std::vector<unsigned int>> d_fragmentHeavyCount;
  //! fill d_fragmentHeavyCount; called from every construction path
  void cacheSynthonSizes();
  //! [position][synthon] -> bonds of that synthon to cut after assembly
  std::vector<std::vector<std::vector<unsigned int>>> d_synthonCutBonds;
  bool d_haveCutBonds = false;
};

}  // namespace RDKit

#endif
