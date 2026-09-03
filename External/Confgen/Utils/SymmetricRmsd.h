//
//  Symmetry-aware conformer deduplication: heavy-atom graph automorphisms (from
//  canonical ranks, no tie-breaking) + minimum QCP RMSD over those
//  automorphisms. Hydrogens are ignored throughout (they do not affect
//  heavy-atom RMSD and their symmetry follows the heavy skeleton).
//
#ifndef RDKIT_SYMMETRIC_RMSD_H
#define RDKIT_SYMMETRIC_RMSD_H

#include <RDGeneral/export.h>
#include <vector>

namespace RDKit {
class ROMol;

//! Return heavy atom automorphisms
RDKIT_FRAGMENTCONFGEN_EXPORT std::vector<std::vector<unsigned int>>
heavyAtomAutomorphisms(const ROMol &mol, std::vector<unsigned int> &heavyOut,
                       size_t maxAutos = 4096);

//! Prunes conformer RMSD favoring first added
//!  Symmetry aware
class RDKIT_FRAGMENTCONFGEN_EXPORT RMSDPruner {
 public:
  //! Precomputes the heavy-atom set and symmetry automorphisms of `mol`.
  //! `rmsThresh` is the diversity threshold (Angstrom); <= 0 keeps everything.
  RMSDPruner(const ROMol &mol, double rmsThresh);

  //! minimum RMSD of coords to current set
  double RMSD(const std::vector<double> &coords) const;

  //! Add the coords to the set if they are outside the prune rms
  //!  returns true if kept, false otherwise
  bool add(const std::vector<double> &coords);

  size_t size() const { return d_kept.size(); }
  size_t numAutomorphisms() const { return d_permB.size(); }
  const std::vector<std::vector<double>> &kept() const { return d_kept; }

 private:
  double d_thresh;
  std::vector<unsigned int> d_heavy;  //!< heavy atom indices (idxA)
  std::vector<std::vector<unsigned int>>
      d_permB;  //!< per-automorphism b-index list
  std::vector<std::vector<double>>
      d_kept;  //!< kept full-mol coordinate buffers
  //! Reusable scratch for minRmsdToKept, sized once in the ctor so the O(kept x
  //! automorphism) comparison loop never allocates: d_scratchA holds the
  //! current kept conformer's centred heavy coords; d_candB / d_candGb hold the
  //! candidate's centred heavy coords and inner product per automorphism
  //! (computed once per call, reused across every kept).  mutable:
  //! minRmsdToKept is logically const.
  mutable std::vector<double> d_scratchA;  //!< [3*H]
  mutable std::vector<double> d_candB;     //!< [P*3*H]
  mutable std::vector<double> d_candGb;    //!< [P]
};

}  // namespace RDKit

#endif
