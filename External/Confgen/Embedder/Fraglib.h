//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_FRAGLIB_H
#define RDKIT_FRAGLIB_H

#include <GraphMol/RDKitBase.h>

#include "Sampler/TorsionSampler.h"
#include "Utils/ParamSentinels.h"

#include <Geometry/point.h>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <mutex>
#include <string>
#include <vector>

namespace RDKit {

inline constexpr const char *kFragConfEnergy = "fragConfEnergy";

//! Note: The "_TOR" variants are MMFF94(s) with the refined dihedral parameters
//! of Wahl et al., J. Cheminform. 2019, 11, 53 (doi:10.1186/s13321-019-0371-6),
//! which mainly correct the overestimated C(aromatic)-N(amide) rotation
//! barrier.
const char* const kFFVariants[] = {
    "MMFF94", "MMFF94s", "MMFF94_TOR", "MMFF94s_TOR"
};

inline bool isValidFF(const std::string &v) {
  auto it = std::find(std::begin(kFFVariants), std::end(kFFVariants), v);
  return it != std::end(kFFVariants);
}

enum class FragmentEmbedMode {
  DG,     //!< plain distance geometry
  ETKDG,  //!< experimental-torsion-knowledge distance geometry (ETKDGv3)
};

//! Minimize after embedding or just score to sort frags.
enum class FragmentMinimize {
  None,  //!< single-point energy at the embedded coords (no minimisation)
  Full,  //!< MMFF-minimise and KEEP the minimised coords -- for plain DG, whose
         //!< raw
         //!<   geometry must be relaxed (no knowledge-based torsions to trust)
  Score,      //!< MMFF-minimise for the ENERGY only, keep embedded coords
  ShrugScore  //!< as Score, but shrug off non MMFF energies a tad
};

//! Embedded Fragment class, this differentiates embedding parameters
enum class FragmentClass {
  Rigid,
  SmallRing,
  LargeRing,
  Acyclic,
  Exhaustive,
  Fast
};

//! Operators to modify # of rotors to # sampled conformers
enum class FragmentSamplesOperator {
  MAXIMUM,     //!< Sample the specified maxium
  MULTIPLY,    //!< max conformers * num rotors
  EXPONENTIAL  //!< min(maxSamples, expPower^num_rotors) i.e. 3**NROTORS
};

//! How to embed a particular fragment type
struct FragmentParams {
  int minSamples;
  int maxSamples;
  int maxConfs;
  double eWindow;
  double rmsd;
  FragmentSamplesOperator op;
  int confsPerRotor;  //!< 0 = don't scale the pool by rotor count
  double expPower = 2;

  //! Allow overrides for the defaults
  std::optional<FragmentEmbedMode> embedMode = std::nullopt;
  std::optional<FragmentMinimize> minimizeMode = std::nullopt;
  int minimizeMaxIters = AutoI;

  //! The number of conformers to sample to generate fragments
  int computeMaxSamples(int num_rotors) const {
    switch (op) {
      case FragmentSamplesOperator::MAXIMUM:
        return maxSamples;
      case FragmentSamplesOperator::MULTIPLY:
        return std::max(
            minSamples,
            static_cast<int>(std::min(maxSamples, num_rotors * confsPerRotor)));
      case FragmentSamplesOperator::EXPONENTIAL:
        return std::max(minSamples,
                        std::min(maxSamples, static_cast<int>(std::pow(
                                                 expPower, num_rotors))));
    }
    return maxSamples;
  }
};

inline bool operator==(const FragmentParams &a, const FragmentParams &b) {
  return a.minSamples == b.minSamples && a.maxSamples == b.maxSamples &&
         a.maxConfs == b.maxConfs && a.eWindow == b.eWindow &&
         a.rmsd == b.rmsd && a.op == b.op &&
         a.confsPerRotor == b.confsPerRotor && a.expPower == b.expPower &&
         a.embedMode == b.embedMode && a.minimizeMode == b.minimizeMode &&
         a.minimizeMaxIters == b.minimizeMaxIters;
}

//! Default embedding parameters for a class
RDKIT_FRAGMENTCONFGEN_EXPORT FragmentParams
getDefaultFragmentParams(FragmentClass cls);

//! Default embedding parameters
RDKIT_FRAGMENTCONFGEN_EXPORT std::map<FragmentClass, FragmentParams>
getDefaultFragmentParams();

//! Embedding parameters, controls DG type, post-minimization type
//!  and how to sample conformers
struct RDKIT_FRAGMENTCONFGEN_EXPORT FraglibParams {
  //! number of 3D conformers generated for each fragment
  unsigned int numConfsPerFragment = 10;
  //! DG vs ETKDG for the embedding
  FragmentEmbedMode fragmentEmbedMode = FragmentEmbedMode::ETKDG;
  //! random seed forwarded to the embedder (-1 == not set)
  int randomSeed = AutoI;
  //! How each embedded conformer is scored (some scores minimize as well)
  FragmentMinimize minimizeMode = FragmentMinimize::Full;

  //! Heavy-atom flat-bottom half-width (A) for ShrugScore.
  double shrugDisplacement = 0.1;
  //! FFFF variant used to score/minimise the embedded fragment conformers.
  std::string ffVariant = "MMFF94";
  //! MMFF minimisation gradient tolerance (loose default keeps distinct pucker
  //! basins).
  double minimizeGradTol = 0.25;
  //! max iterations for the per-conformer MMFF minimisation (signed: -1 = AUTO,
  //! 0 = off)
  int minimizeMaxIters = 1000;

  //! Fragment type parameterizations
  //! To modify, use:
  //!   params.setClassParam(FragmentClass::SmallRing, FragmentParams);
  std::map<FragmentClass, FragmentParams> classParams =
      getDefaultFragmentParams();

  //! Set (add or replace) the embedding recipe for one class.  Equivalent to
  //! `classParams[cls] = fp`, but named for intent.
  void setClassParam(FragmentClass cls, const FragmentParams &fp) {
    classParams[cls] = fp;
  }

  //! Helper to set a maximum sample size for ALL fragment types
  void setFlatPool(int n) {
    for (auto &kv : classParams) {
      kv.second.op = FragmentSamplesOperator::MAXIMUM;
      kv.second.minSamples = kv.second.maxSamples = n;
    }
  }
  //! Use a per class embedding using different params for Acyclic, Ring etc.
  bool perClassEmbedding = true;
  //! Energy window (kcal/mol) for pooled-conformer selection (perClassEmbedding
  //! sets this per class; 0 = keep all within the pool regardless of energy).
  double energyWindow = Disabled;

  //! XXX Dead Code - slated for removal
  bool FRAGLIB_TRACE = false;  //!< if true debug fraglib steps to stderr
  //! DIAGNOSTIC: per embedded fragment conformer, log class / rotor count /
  //! MMFF energy BEFORE and AFTER minimisation (and the delta).  Shows how much
  //! strain the minimiser actually removes per fragment class -- e.g. whether
  //! floppy Acyclic chains are still high after minimisation (i.e.
  //! under-converged and wanting more iterations) vs rings that relax cleanly.
  bool logFragmentEnergies = false;
};

//! Two FraglibParams are equal iff EVERY embedding option matches (seed
//! included).
inline bool operator==(const FraglibParams &a, const FraglibParams &b) {
  return a.numConfsPerFragment == b.numConfsPerFragment &&
         a.fragmentEmbedMode == b.fragmentEmbedMode &&
         a.randomSeed == b.randomSeed && a.minimizeMode == b.minimizeMode &&
         a.shrugDisplacement == b.shrugDisplacement &&
         a.minimizeGradTol == b.minimizeGradTol &&
         a.minimizeMaxIters == b.minimizeMaxIters &&
         a.perClassEmbedding == b.perClassEmbedding &&
         a.energyWindow == b.energyWindow && a.classParams == b.classParams;
}
inline bool operator!=(const FraglibParams &a, const FraglibParams &b) {
  return !(a == b);
}

//! Are two embedding parameters identical
inline bool sameEmbeddingType(const FraglibParams &a, const FraglibParams &b) {
  return a.ffVariant == b.ffVariant &&
         a.numConfsPerFragment == b.numConfsPerFragment &&
         a.fragmentEmbedMode == b.fragmentEmbedMode &&
         a.minimizeMode == b.minimizeMode &&
         a.shrugDisplacement == b.shrugDisplacement &&
         a.minimizeGradTol == b.minimizeGradTol &&
         a.minimizeMaxIters == b.minimizeMaxIters &&
         a.perClassEmbedding == b.perClassEmbedding &&
         a.energyWindow == b.energyWindow && a.classParams == b.classParams;
}

//! The fragment library
class RDKIT_FRAGMENTCONFGEN_EXPORT Fraglib {
 public:
  explicit Fraglib(FraglibParams params = FraglibParams())
      : d_params(std::move(params)) {}
  ~Fraglib();

  Fraglib(const Fraglib &) = delete;
  Fraglib &operator=(const Fraglib &) = delete;

  //! Return a the fraglib with generated conformations
  // \param frag - fragment to embed
  // \param cache - true to save int the fraglib, false otherwise
  ROMOL_SPTR get(const ROMol &frag, bool cache=true) const;

 private:
  //! Look `key` up in the cache; on a miss build the fragment with `make` and either retain it
  //!   returns nullptr on failure XXX FIX ME -> std::optional
  const RWMol *lookupOrEmbed(const std::string &key,
                             const std::function<RWMol *()> &make, bool cache,
                             std::unique_ptr<RWMol> &owned) const;

 public:
  bool getConformerCoords(RWMol &frag, unsigned int nMolAtoms,
                          const std::string &molIdxProp,
                          std::vector<std::vector<RDGeom::Point3D>> &out,
                          std::vector<double> *energiesOut = nullptr,
                          bool cache = true) const;

  static std::string cacheKey(const ROMol &frag);

  //! Generate a canoincal cache key for a molecule (here it is isomeric smiles
  //! with dummy
  //!  atoms.  Stereo dummy atoms are labeled with isotopes to preserve stereo
  //!   NOTE: if remap is true,
  //!     the fragment is always reordered in canonical order
  //!     this makes conformation copying and atom labeling trivial.
  //!     The cache key uses isotope labels to preserve stereo
  //!     While the molzip code uses isotope labels to preserve exit
  //!      vectors.  This allows us to trivial copy atom props when
  //!      necessary.
  //! If `outOrder` is non-null return the canonical atom order
  static std::string generateKey(RWMol &frag, bool remap = true,
                                 std::vector<unsigned int> *outOrder = nullptr);

  //! number of cached/embedded fragments (tombstones excluded)
  size_t size() const;

  //! Fragments that could NOT be embedded and are remembered as such.
  //!  These could be unphysical stereo or other effects
  size_t numUnembeddable() const;

  //! fragment count - number of conformers for fragment
  //! Record a fragment as permanently unusable.
  /*!
    This keeps the fragment in the fraglib but notes that we can't
    embed it for whatever reason to prevent future potentially costly
    attempts.

    It also can be used to log future efforts in forcefields/etc.

    Returns trues if this call marked the tombstone, or false
    if already marked or a live entry with coords exists
  */
  bool markUnembeddable(const ROMol &frag);

  //! Conformer count for a cached fragment.
  /*!
    returns
     nullopt if the fragment isn't in the library
     0 if the fragment has been marked as unembeddable
     >0 is the actual conformer count for the fragment.
  */
  std::optional<unsigned int> numFragmentConfs(const ROMol &frag) const;

  const FraglibParams &params() const { return d_params; }

  //! serialize the fragment library to a stream
  void serialize(std::ostream &os) const;
  //! unserialize the fragment library from a stream
  void initFromStream(std::istream &is);

  //! Write the fragments/conformers as an SDF complete with
  //!  MMFF energies and other diagnostics
  void writeSDF(std::ostream &os) const;

  //! Inject a pre-built fragment, can be used for testing/debugging
  std::string add(RWMol &frag) const;

 private:
  mutable std::map<std::string, RWMol *> d_fraglib;
  mutable std::mutex d_mutex;
  FraglibParams d_params;
};

}  // namespace RDKit

#endif
