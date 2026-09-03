//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_NONBONDEDLOOKUP_H
#define RDKIT_NONBONDEDLOOKUP_H

#include <RDGeneral/export.h>
#include <RDGeneral/Invariant.h>
#include <ForceField/Contrib.h>
#include <ForceField/ForceField.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace ForceFields {
namespace MMFF {

namespace detail {
constexpr unsigned int VDW_LOOKUP_NUM_ATOM_TYPES = 6;
constexpr unsigned int VDW_LOOKUP_NUM_TABLES =
    (VDW_LOOKUP_NUM_ATOM_TYPES * (VDW_LOOKUP_NUM_ATOM_TYPES + 1)) / 2;
constexpr unsigned int VDW_LOOKUP_NUM_DISTANCES = 401;
constexpr double VDW_LOOKUP_MIN_DISTANCE = 0.0;
constexpr double VDW_LOOKUP_DISTANCE_STEP = 0.02;
constexpr double VDW_LOOKUP_MAX_DISTANCE =
    VDW_LOOKUP_MIN_DISTANCE +
    VDW_LOOKUP_DISTANCE_STEP * (VDW_LOOKUP_NUM_DISTANCES - 1);

struct VdWTableParams {
  double R_ij_star;
  double epsilon;
};

constexpr double calcVdWEnergy(double dist, double R_star_ij,
                               double wellDepth) {
  constexpr double vdw1 = 1.07;
  constexpr double vdw1m1 = vdw1 - 1.0;
  constexpr double vdw2 = 1.12;
  constexpr double vdw2m1 = vdw2 - 1.0;
  const double dist2 = dist * dist;
  const double dist7 = dist2 * dist2 * dist2 * dist;
  const double aTerm = vdw1 * R_star_ij / (dist + vdw1m1 * R_star_ij);
  const double aTerm2 = aTerm * aTerm;
  const double aTerm7 = aTerm2 * aTerm2 * aTerm2 * aTerm;
  const double R_star_ij2 = R_star_ij * R_star_ij;
  const double R_star_ij7 = R_star_ij2 * R_star_ij2 * R_star_ij2 * R_star_ij;
  const double bTerm = vdw2 * R_star_ij7 / (dist7 + vdw2m1 * R_star_ij7) - 2.0;
  return wellDepth * aTerm7 * bTerm;
}

constexpr std::array<std::array<double, VDW_LOOKUP_NUM_DISTANCES>,
                     VDW_LOOKUP_NUM_TABLES>
makeVDWLookup(const std::array<VdWTableParams, VDW_LOOKUP_NUM_TABLES> &params) {
  std::array<std::array<double, VDW_LOOKUP_NUM_DISTANCES>,
             VDW_LOOKUP_NUM_TABLES>
      res{};
  for (unsigned int tableIdx = 0; tableIdx < params.size(); ++tableIdx) {
    for (unsigned int distIdx = 0; distIdx < res[tableIdx].size(); ++distIdx) {
      const double dist =
          VDW_LOOKUP_MIN_DISTANCE + distIdx * VDW_LOOKUP_DISTANCE_STEP;
      res[tableIdx][distIdx] = calcVdWEnergy(dist, params[tableIdx].R_ij_star,
                                             params[tableIdx].epsilon);
    }
  }
  return res;
}
}  // namespace detail

//  XXX FIX ME -> the lookup tables aren't worth the effort, remove
// Representative MMFF94 vdW pair parameters for H, C, N, O, P, and S using
// atom types HC, CR, NR, OR, P, and S.
inline constexpr std::array<detail::VdWTableParams,
                            detail::VDW_LOOKUP_NUM_TABLES>
    VDW_LOOKUP_PARAMS{{
        {2.9698484809835, 0.021574173607307283},     // HC-HC
        {3.5987879811530421, 0.02807765636226562},   // HC-CR
        {3.6670841353470527, 0.027718754793617923},  // HC-NR
        {3.3246176293946523, 0.035314214265179299},  // HC-OR
        {4.0871643086409239, 0.039131136790172737},  // HC-P
        {3.9291297295260761, 0.044493160327157992},  // HC-S
        {3.9377389919289634, 0.067796993042913714},  // CR-CR
        {3.9842650025999027, 0.069780057128632833},  // CR-NR
        {3.7706594524857651, 0.067936315743771972},  // CR-OR
        {4.3104969857175233, 0.11923642373602522},   // CR-P
        {4.1800332367337782, 0.12810300080946846},   // CR-S
        {4.0283209169689105, 0.072149830128808698},  // NR-NR
        {3.8274080367246333, 0.068688419743367798},  // NR-OR
        {4.3411549148254762, 0.1260356180528083},    // NR-P
        {4.2153004393535118, 0.13441224434621038},   // NR-S
        {3.5581488427996217, 0.076254597956916734},  // OR-OR
        {4.2043066898612116, 0.10435412898639047},   // OR-P
        {4.0573663924792962, 0.11651347774569312},   // OR-S
        {4.5731317025495741, 0.25957008521351627},   // P-P
        {4.4768030664416978, 0.26092978955926649},   // P-S
        {4.3693657230022742, 0.26808283125482379},   // S-S
    }};
inline constexpr auto VDW_LOOKUP = detail::makeVDWLookup(VDW_LOOKUP_PARAMS);

// Return the row in VDW_LOOKUP with the closest CHNOPS representative vdW
// parameters.
inline unsigned int closest_table(double R_ij_star, double epsilon) {
  PRECONDITION(R_ij_star > 0.0, "bad MMFF VdW R_ij_star");
  PRECONDITION(epsilon > 0.0, "bad MMFF VdW epsilon");

  unsigned int bestIdx = 0;
  double bestDistance = std::numeric_limits<double>::max();
  for (unsigned int tableIdx = 0; tableIdx < VDW_LOOKUP_PARAMS.size();
       ++tableIdx) {
    const double dR = (R_ij_star - VDW_LOOKUP_PARAMS[tableIdx].R_ij_star) /
                      VDW_LOOKUP_PARAMS[tableIdx].R_ij_star;
    const double dE = (epsilon - VDW_LOOKUP_PARAMS[tableIdx].epsilon) /
                      VDW_LOOKUP_PARAMS[tableIdx].epsilon;
    const double paramDistance = dR * dR + dE * dE;
    if (paramDistance < bestDistance) {
      bestDistance = paramDistance;
      bestIdx = tableIdx;
    }
  }
  return bestIdx;
}

// Return the interpolated energy from table_idx for distance.
inline double table_lookup(unsigned int table_idx, double distance) {
  URANGE_CHECK(table_idx, VDW_LOOKUP.size());
  if (distance <= detail::VDW_LOOKUP_MIN_DISTANCE) {
    return VDW_LOOKUP[table_idx].front();
  }
  if (distance >= detail::VDW_LOOKUP_MAX_DISTANCE) {
    return 0.0;
  }

  const double gridPosition = (distance - detail::VDW_LOOKUP_MIN_DISTANCE) /
                              detail::VDW_LOOKUP_DISTANCE_STEP;
  const auto lowerIdx = static_cast<unsigned int>(gridPosition);
  const double fraction = gridPosition - lowerIdx;
  return VDW_LOOKUP[table_idx][lowerIdx] * (1.0 - fraction) +
         VDW_LOOKUP[table_idx][lowerIdx + 1] * fraction;
}

//! the van der Waals term for MMFF using a precached CHNOPS lookup table
class RDKIT_FORCEFIELD_EXPORT VdWContribLookup : public ForceFieldContrib {
 public:
  VdWContribLookup() {}
  VdWContribLookup(ForceField *owner) {
    PRECONDITION(owner, "bad owner");
    dp_forceField = owner;
  }

  //! Track a new VdW pair
  void addTerm(unsigned int idx1, unsigned int idx2,
               const MMFFVdWRijstarEps *mmffVdWConstants) {
    PRECONDITION(mmffVdWConstants, "bad MMFFVdW parameters");
    URANGE_CHECK(idx1, dp_forceField->positions().size());
    URANGE_CHECK(idx2, dp_forceField->positions().size());
    d_at1Idxs.push_back(idx1);
    d_at2Idxs.push_back(idx2);
    d_tableIdxs.push_back(
        closest_table(mmffVdWConstants->R_ij_star, mmffVdWConstants->epsilon));
  }

  double getEnergy(double *pos) const override {
    PRECONDITION(dp_forceField, "no owner");
    PRECONDITION(pos, "bad vector");
    double energySum = 0.0;

    const auto numPairs = static_cast<unsigned int>(d_at1Idxs.size());
    for (unsigned int i = 0; i < numPairs; ++i) {
      const auto at1Idx = static_cast<unsigned int>(d_at1Idxs[i]);
      const auto at2Idx = static_cast<unsigned int>(d_at2Idxs[i]);
      const double dist = dp_forceField->distance(at1Idx, at2Idx, pos);
      energySum += table_lookup(d_tableIdxs[i], dist);
    }
    return energySum;
  }

  void getGrad(double *, double *) const override {}
  VdWContribLookup *copy() const override {
    return new VdWContribLookup(*this);
  }

 private:
  std::vector<int16_t> d_at1Idxs;
  std::vector<int16_t> d_at2Idxs;
  std::vector<std::uint8_t> d_tableIdxs;
};

//! Inter-fragment MMFF van der Waals with a distance cutoff (heuristic)
//!
class RDKIT_FORCEFIELD_EXPORT InterFragVdWContrib : public ForceFieldContrib {
 public:
  InterFragVdWContrib() {}
  InterFragVdWContrib(ForceField *owner, double cutoff, bool useLookup)
      : d_cut2(cutoff > 0.0 ? cutoff * cutoff : 0.0), d_useLookup(useLookup) {
    PRECONDITION(owner, "bad owner");
    dp_forceField = owner;
  }

  void addTerm(unsigned int idx1, unsigned int idx2,
               const MMFFVdWRijstarEps *mmffVdWConstants) {
    PRECONDITION(mmffVdWConstants, "bad MMFFVdW parameters");
    d_at1Idxs.push_back(idx1);
    d_at2Idxs.push_back(idx2);
    d_Rstars.push_back(mmffVdWConstants->R_ij_star);
    d_wellDepths.push_back(mmffVdWConstants->epsilon);
    d_tableIdxs.push_back(
        closest_table(mmffVdWConstants->R_ij_star, mmffVdWConstants->epsilon));
  }

  double getEnergy(double *pos) const override {
    PRECONDITION(dp_forceField, "no owner");
    PRECONDITION(pos, "bad vector");
    const size_t numPairs = d_at1Idxs.size();
    double energySum = 0.0;
    const double cut2 = d_cut2;

    if (d_useLookup) {
      // XXX the lookup isn't worth it in practice, left for posterity
      for (size_t i = 0; i < numPairs; ++i) {
        const unsigned int a = d_at1Idxs[i], b = d_at2Idxs[i];
        const double dx = pos[3 * a] - pos[3 * b];
        const double dy = pos[3 * a + 1] - pos[3 * b + 1];
        const double dz = pos[3 * a + 2] - pos[3 * b + 2];
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (cut2 > 0.0 && d2 > cut2) {
          continue;
        }
        energySum += table_lookup(d_tableIdxs[i], std::sqrt(d2));
      }
      return energySum;
    }

    // Exact buffered-14-7: IMPORTANT inlined so the loop VECTORISES and restricted
    //  so the compiler can do the right thing.
    // NOTE: this bit was suggested by claude.ai
    const int *__restrict a1 = d_at1Idxs.data();
    const int *__restrict a2 = d_at2Idxs.data();
    const double *__restrict Rs = d_Rstars.data();
    const double *__restrict wd = d_wellDepths.data();
#if defined(__clang__)
#pragma clang loop vectorize(enable)
#endif
    for (size_t i = 0; i < numPairs; ++i) {
      const int a = a1[i], b = a2[i];
      const double dx = pos[3 * a] - pos[3 * b];
      const double dy = pos[3 * a + 1] - pos[3 * b + 1];
      const double dz = pos[3 * a + 2] - pos[3 * b + 2];
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (cut2 > 0.0 && d2 > cut2) {
        continue;  // masked in the vectorised loop
      }
      energySum += detail::calcVdWEnergy(std::sqrt(d2), Rs[i], wd[i]);
    }
    return energySum;
  }

  //! XXX FIX ME: unused get the atoms on each side of a pair
  size_t numPairs() const { return d_at1Idxs.size(); }
  int at1(size_t i) const { return d_at1Idxs[i]; }
  int at2(size_t i) const { return d_at2Idxs[i]; }

  //! vdw params for pair i
  double rStar(size_t i) const { return d_Rstars[i]; }
  double wellDepth(size_t i) const { return d_wellDepths[i]; }
  //! Squared distance cutoff (0 = none).
  double cut2() const { return d_cut2; }
  bool usesLookup() const { return d_useLookup; }

  //! Suggested by claude.ai as an optimization
  //! Vectorized energy of a COMPACT, contiguous subset of pairs (SoA arrays
  //! a1/a2/Rstar/wd of length `n`).
  //! This is the fast
  //! incremental primitive: RotorTree packs each rotor's changing pairs once,
  //! then rescores a single-rotor move with two calls here.
  static double energyOfPacked(const double *pos, size_t n,
                               const int *__restrict a1,
                               const int *__restrict a2,
                               const double *__restrict Rs,
                               const double *__restrict wd, double cut2) {
    double energySum = 0.0;
#if defined(__clang__)
#pragma clang loop vectorize(enable)
#endif
    for (size_t i = 0; i < n; ++i) {
      const int a = a1[i], b = a2[i];
      const double dx = pos[3 * a] - pos[3 * b];
      const double dy = pos[3 * a + 1] - pos[3 * b + 1];
      const double dz = pos[3 * a + 2] - pos[3 * b + 2];
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (cut2 > 0.0 && d2 > cut2) {
        continue;
      }
      energySum += detail::calcVdWEnergy(std::sqrt(d2), Rs[i], wd[i]);
    }
    return energySum;
  }

  //! Energy of ONLY the listed pair indices (indices into the stored pair
  //! arrays)  This optimization allows the rotor hierarchy to adjust only the
  //!  moving bits
  double energyOfPairs(const double *pos,
                       const std::vector<unsigned int> &pairIdxs) const {
    PRECONDITION(pos, "bad vector");
    double energySum = 0.0;
    const double cut2 = d_cut2;
    if (d_useLookup) {
      for (unsigned int i : pairIdxs) {
        const unsigned int a = d_at1Idxs[i], b = d_at2Idxs[i];
        const double dx = pos[3 * a] - pos[3 * b];
        const double dy = pos[3 * a + 1] - pos[3 * b + 1];
        const double dz = pos[3 * a + 2] - pos[3 * b + 2];
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (cut2 > 0.0 && d2 > cut2) {
          continue;
        }
        energySum += table_lookup(d_tableIdxs[i], std::sqrt(d2));
      }
      return energySum;
    }
    for (unsigned int i : pairIdxs) {
      const int a = d_at1Idxs[i], b = d_at2Idxs[i];
      const double dx = pos[3 * a] - pos[3 * b];
      const double dy = pos[3 * a + 1] - pos[3 * b + 1];
      const double dz = pos[3 * a + 2] - pos[3 * b + 2];
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (cut2 > 0.0 && d2 > cut2) {
        continue;
      }
      energySum +=
          detail::calcVdWEnergy(std::sqrt(d2), d_Rstars[i], d_wellDepths[i]);
    }
    return energySum;
  }

  void getGrad(double *, double *) const override {}
  InterFragVdWContrib *copy() const override {
    return new InterFragVdWContrib(*this);
  }

 private:
  std::vector<int> d_at1Idxs;
  std::vector<int> d_at2Idxs;
  std::vector<double> d_Rstars;
  std::vector<double> d_wellDepths;
  std::vector<std::uint8_t> d_tableIdxs;
  double d_cut2 = 0.0;
  bool d_useLookup = false;
};

}  // namespace MMFF
}  // namespace ForceFields

#endif
