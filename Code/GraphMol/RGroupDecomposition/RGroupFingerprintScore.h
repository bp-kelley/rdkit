//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_RGROUPFINGERPRINTSCORE_H
#define RDKIT_RGROUPFINGERPRINTSCORE_H

#include "RGroupMatch.h"
#include <vector>

namespace RDKit {

// class to hold the bitcounts for an attachment point/rgroup label
struct VarianceDataForLabel {
  // rgroup label
  int label;
  // number of structures attached here
  int numberFingerprints;
  // bitcounts - size fingerprint size, each position is the count of bits set
  // over the fingerprints for all the structures
  std::vector<int> bitCounts;
  
  VarianceDataForLabel() : label(), numberFingerprints(), bitCounts() {
  }
  VarianceDataForLabel(const int &label, int numberFingerprints,
                       std::vector<int> bitCounts);
  VarianceDataForLabel(const int &label);
  VarianceDataForLabel(const VarianceDataForLabel &other) = default;
  VarianceDataForLabel &operator=(const VarianceDataForLabel &other) = delete;
  // add an rgroup structure to a bit counts array
  void addRgroupData(RGroupData *rgroupData);
  // remove an rgroup structure to a bit counts array
  void removeRgroupData(RGroupData *rgroupData);
  // calculate the mean variance for a bit counts array
  double variance() const;
private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int /*version*/) {
    ar &label;
    ar &numberFingerprints;
    ar &bitCounts;
  }
  #endif
};

struct FingerprintVarianceScoreData {
  size_t numberOfMissingUserRGroups = 0;
  size_t numberOfMolecules = 0;
  std::map<int, std::shared_ptr<VarianceDataForLabel>> labelsToVarianceData;

  // calculates fingerprint variance score from rgroup bit counts
  double fingerprintVarianceGroupScore();

  // Adds a molecule match to the rgroup fingerprint bit counts
  // vectors
  void addVarianceData(int matchNumber, int permutationNumber,
                       const std::vector<std::vector<RGroupMatch>> &matches,
                       const std::set<int> &labels);

  // Subtracts a molecule match from the rgroup fingerprint bit counts
  // vectors
  void removeVarianceData(int matchNumber, int permutationNumber,
                          const std::vector<std::vector<RGroupMatch>> &matches,
                          const std::set<int> &labels);

  void clear();

 private:
  void modifyVarianceData(int matchNumber, int permutationNumber,
                          const std::vector<std::vector<RGroupMatch>> &matches,
                          const std::set<int> &labels, bool add);
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void save(Archive &ar, const unsigned int /*version*/) const {
    ar &numberOfMissingUserRGroups;
    ar &numberOfMolecules;
    ar &labelsToVarianceData;
  }
  template <class Archive>
  void load(Archive &ar, const unsigned int /*version*/) {
    ar &numberOfMissingUserRGroups;
    ar &numberOfMolecules;
    ar &labelsToVarianceData;
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER();
#endif
};

// The arithmetic mean of the mean fingerprint bit variances for the
// fingerprints at each rgroup position.
RDKIT_RGROUPDECOMPOSITION_EXPORT double fingerprintVarianceScore(
    const std::vector<size_t> &permutation,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels,
    FingerprintVarianceScoreData *fingerprintVarianceScoreData = nullptr);

}  // namespace RDKit

#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_CLASS_VERSION(RDKit::VarianceDataForLabel, 1)
BOOST_CLASS_VERSION(RDKit::FingerprintVarianceScoreData, 1)
#endif

#endif  // RDKIT_RGROUPFINGERPRINTSCORE_H
