//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_TORSION_SAMPLER_H
#define RDKIT_TORSION_SAMPLER_H

#include <GraphMol/RDKitBase.h>
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace RDKit {

//! Load state of a sampler backed by an external data file.
/*!
  The file-backed samplers (TorsionLibrarySampler, SmirnoffTorsionSampler) read
  data that is deliberately NOT bundled with this library. so an absent or

  When presented with no rules, the sampler will fall back to a uniform
  search.
*/
enum class TorsionSamplerStatus {
  Ok,              //!< data loaded; at least one rule is available
  NoPathGiven,     //!< no path was supplied
  FileUnreadable,  //!< a path was supplied but could not be opened
  NoRulesParsed    //!< the file was read but yielded no usable rules
};

//! Human-readable one-liner for a status, for error messages and logs.
RDKIT_FRAGMENTCONFGEN_EXPORT const char *torsionSamplerStatusMessage(
    TorsionSamplerStatus status);

//! Suppy dihedrals when rotor-driving a single rotatable bond.
/*!
  Generate dihedrals to sample given a torsion
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT TorsionSampler {
 public:
  virtual ~TorsionSampler() = default;

  //! Return the dihedral angles (degrees) to sample for the i-j-k-l torsion.
  //!  \param i
  //!  \param j
  //!  \param k
  //!  \param l
  //!  \param addBasinAngle - add samples around the basin angles [default True]
  virtual std::vector<double> getAngles(const ROMol &mol, unsigned int i,
                                        unsigned int j, unsigned int k,
                                        unsigned int l,
                                        bool addBasinAngle = true) const = 0;

  virtual std::shared_ptr<TorsionSampler> copy() const = 0;

  //! Whether this sampler's external data (if any) loaded.
  /*!
    Samplers that need no external data are always Ok.  A sampler that is not
    Ok will fall back to a uniform search.
  */
  virtual TorsionSamplerStatus status() const {
    return TorsionSamplerStatus::Ok;
  }

  bool isValid() const { return status() == TorsionSamplerStatus::Ok; }

 protected:
  //! Collapse rotationally-equivalent angles for a symmetric rotor (t-Bu, CF3,
  //! para-phenyl...).
  //! \param useSmarts  choose SMARTS-based symmetry perception over symmetry
  //! classes
  //!                   This requires an external smarts library (not supplied)
  static std::vector<double> foldSymmetricAngles(const ROMol &mol,
                                                 unsigned int j, unsigned int k,
                                                 std::vector<double> angles,
                                                 double minSeparation,
                                                 bool useSmarts = false);
};

//! Samples the full circle [0, 360) at a fixed angular increment.
class RDKIT_FRAGMENTCONFGEN_EXPORT UniformTorsionSampler
    : public TorsionSampler {
  double d_stepDeg{5.0};

 public:
  explicit UniformTorsionSampler(double stepDeg = 5.0) : d_stepDeg(stepDeg) {}

  std::vector<double> getAngles(const ROMol &mol, unsigned int i,
                                unsigned int j, unsigned int k, unsigned int l,
                                bool addBasinAngle = true) const override;

  std::shared_ptr<TorsionSampler> copy() const override {
    return std::make_shared<UniformTorsionSampler>(*this);
  }
};

//! Samples low-energy ETKDG experimental torsion-preference angles.
/*!
  The sampler uses the CrystalFF/ETKDG torsion preference rules to identify
  candidate torsions around the driven bond.

  If no experimental torsion rule matches, this falls back to a uniform rotation
  sampling.
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT ETKDGTorsionSampler : public TorsionSampler {
  double d_scanStepDeg{10.0};
  double d_fallbackStepDeg{60.0};
  //! Max candidate angles returned per rotor; **0 == return all**
  //!  note this limit is applied bost symmetry folding
  unsigned int d_maxAngles{6};
  unsigned int d_version{2};
  bool d_useSmallRingTorsions{false};
  bool d_useMacrocycleTorsions{false};

 public:
  ETKDGTorsionSampler(double scanStepDeg = 10.0, double fallbackStepDeg = 60.0,
                      unsigned int maxAngles = 6, unsigned int version = 2,
                      bool useSmallRingTorsions = false,
                      bool useMacrocycleTorsions = false)
      : d_scanStepDeg(scanStepDeg),
        d_fallbackStepDeg(fallbackStepDeg),
        d_maxAngles(maxAngles),
        d_version(version),
        d_useSmallRingTorsions(useSmallRingTorsions),
        d_useMacrocycleTorsions(useMacrocycleTorsions) {}

  //! Fold rotationally-equivalent angles for symmetric substituents (t-Bu, CF3,
  //! phenyl...).
  //!  3-4% faster than without with no accuracy loss
  bool foldBySymmetry{true};

  std::vector<double> getAngles(const ROMol &mol, unsigned int i,
                                unsigned int j, unsigned int k, unsigned int l,
                                bool addBasinAngle = true) const override;

  std::shared_ptr<TorsionSampler> copy() const override {
    return std::make_shared<ETKDGTorsionSampler>(*this);
  }
};

//! Samples torsions from a Hamburg style torsion library
/*!
  The current XML parser reads files in the the Hamburg `TorsionLibrary.xml`
  format
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT TorsionLibrarySampler
    : public TorsionSampler {
  struct Rule;
  struct RuleTable;

  //! Parse (or reuse) the rule table for this path.  Tables are immutable and
  //! shared process-wide, keyed on the path, so each distinct library is parsed
  //! ONCE however many samplers or threads use it.  Never returns null and
  //! never throws; a load failure comes back as a table whose status says why.
  //! NOTE: the cache is keyed on path only -- a file edited mid-process is not
  //! re-read.
  static std::shared_ptr<const RuleTable> acquireRuleTable(
      const std::string &torsionLibraryPath);

  std::string d_torsionLibraryPath;
  double d_fallbackStepDeg{30.0};  // for unknown torsions
  unsigned int d_maxAngles{0};  //!< max candidate angles per rotor; 0 == no cap
  //! Immutable, shared.  Acquired in the constructor, so there is no lazy load
  //! and no mutable state -- copies and concurrent readers are safe.
  std::shared_ptr<const RuleTable> d_table;

 public:
  struct Options {
    bool symmetryFolding =
        true;  //!< fold rotationally-symmetrically equivalent angles
    bool terminalOffsets = true;
    bool toleranceRanges = true;       //!< add samples to basin angles
    bool wideToleranceRanges = false;  //! sample basin angles more densely
    bool lean = false;                 //!< Force all over-basin sampling to off
    bool symmetryUseSmarts =
        false;  //!< Choose smarts symmetry folding over graph-based symmetry
    bool mostSpecificMatch =
        false;  //!< pick the most-specific matching rule vs first-in-file
    bool dump = false;  //!< DIAGNOSTIC: per-bond assignment to stderr
  };
  Options options;

  //! \param torsionLibraryPath  the ONE rule file to read.  To combine several
  //!        rule sets, concatenate them into a single file in precedence order
  //!        -- rules are matched first-in-file-wins, which a
  //!        CompositeTorsionSampler cannot express (it UNIONS its children's
  //!        angles, which flattens the informed prior rather than deferring).
  TorsionLibrarySampler(std::string torsionLibraryPath,
                        double fallbackStepDeg = 30.0,
                        unsigned int maxAngles = 0)
      : d_torsionLibraryPath(std::move(torsionLibraryPath)),
        d_fallbackStepDeg(fallbackStepDeg),
        d_maxAngles(maxAngles),
        d_table(acquireRuleTable(d_torsionLibraryPath)) {}

  std::vector<double> getAngles(const ROMol &mol, unsigned int i,
                                unsigned int j, unsigned int k, unsigned int l,
                                bool addBasinAngle = true) const override;

  TorsionSamplerStatus status() const override;

  std::shared_ptr<TorsionSampler> copy() const override {
    return std::make_shared<TorsionLibrarySampler>(*this);
  }
};

//! Samples torsions from a SMIRNOFF Open Force Field proper-torsion set
//! (.offxml).
/*!
  N.b. for smirnoff Last matching reaction smarts wins

  n.b. The OpenFF SMIRNOFF parameter files (.offxml, e.g.
  openff-2.x "Sage") are distributed by the Open Force Field Initiative under
  CC-BY-4.0 -- they are NOT bundled here.

  SMIRNOFF encodes an energy profile, angles are energy minima and are sampled
  from there. Falls back to a uniform grid if not pattern is found.
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT SmirnoffTorsionSampler
    : public TorsionSampler {
  struct Proper;
  struct ProperTable;

  //! See TorsionLibrarySampler::acquireRuleTable -- same shared, immutable,
  //! path-keyed, non-throwing contract.
  static std::shared_ptr<const ProperTable> acquireProperTable(
      const std::string &offxmlPath);

  std::string d_offxmlPath;
  double d_fallbackStepDeg{120.0};
  unsigned int d_maxAngles{0};  //!< max candidate angles per rotor; 0 == no cap
  double d_minimaScanStepDeg{2.0};
  std::shared_ptr<const ProperTable> d_table;

 public:
  explicit SmirnoffTorsionSampler(std::string offxmlPath,
                                  double fallbackStepDeg = 120.0,
                                  unsigned int maxAngles = 0)
      : d_offxmlPath(std::move(offxmlPath)),
        d_fallbackStepDeg(fallbackStepDeg),
        d_maxAngles(maxAngles),
        d_table(acquireProperTable(d_offxmlPath)) {}

  std::vector<double> getAngles(const ROMol &mol, unsigned int i,
                                unsigned int j, unsigned int k, unsigned int l,
                                bool addBasinAngle = true) const override;

  TorsionSamplerStatus status() const override;

  std::shared_ptr<TorsionSampler> copy() const override {
    return std::make_shared<SmirnoffTorsionSampler>(*this);
  }
};

//! Returns the deduplicated union of several samplers' angles.
/*!
  Composite Sampler provides the uniform grid when sampling fails.
  This helps to remove strain outliers occasionally.
*/
class RDKIT_FRAGMENTCONFGEN_EXPORT CompositeTorsionSampler
    : public TorsionSampler {
  std::vector<std::shared_ptr<TorsionSampler>> d_samplers;

 public:
  explicit CompositeTorsionSampler(
      std::vector<std::shared_ptr<TorsionSampler>> samplers)
      : d_samplers(std::move(samplers)) {}

  std::vector<double> getAngles(const ROMol &mol, unsigned int i,
                                unsigned int j, unsigned int k, unsigned int l,
                                bool addBasinAngle = true) const override;

  //! First non-Ok child status, else Ok -- a composite quietly holding a
  //! sampler whose library failed to load is exactly what this is here to
  //! catch.
  TorsionSamplerStatus status() const override;

  std::shared_ptr<TorsionSampler> copy() const override {
    return std::make_shared<CompositeTorsionSampler>(*this);
  }
};

}  // namespace RDKit

#endif
