//
//  Copyright (C) 2026 The RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "NormalizedDescriptors.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/BCUT.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace RDKit {
namespace NormalizedDescriptors {
namespace {
// The non-fragment descriptors in rdkit.Chem.Descriptors._descList, in that
// order
const std::vector<std::string> descriptorNames = {
    "MaxAbsEStateIndex",
    "MaxEStateIndex",
    "MinAbsEStateIndex",
    "MinEStateIndex",
    "qed",
    "SPS",
    "MolWt",
    "HeavyAtomMolWt",
    "ExactMolWt",
    "NumValenceElectrons",
    "NumRadicalElectrons",
    "MaxPartialCharge",
    "MinPartialCharge",
    "MaxAbsPartialCharge",
    "MinAbsPartialCharge",
    "FpDensityMorgan1",
    "FpDensityMorgan2",
    "FpDensityMorgan3",
    "BCUT2D_MWHI",
    "BCUT2D_MWLOW",
    "BCUT2D_CHGHI",
    "BCUT2D_CHGLO",
    "BCUT2D_LOGPHI",
    "BCUT2D_LOGPLOW",
    "BCUT2D_MRHI",
    "BCUT2D_MRLOW",
    "AvgIpc",
    "BalabanJ",
    "BertzCT",
    "Chi0",
    "Chi0n",
    "Chi0v",
    "Chi1",
    "Chi1n",
    "Chi1v",
    "Chi2n",
    "Chi2v",
    "Chi3n",
    "Chi3v",
    "Chi4n",
    "Chi4v",
    "HallKierAlpha",
    "Ipc",
    "Kappa1",
    "Kappa2",
    "Kappa3",
    "LabuteASA",
    "PEOE_VSA1",
    "PEOE_VSA10",
    "PEOE_VSA11",
    "PEOE_VSA12",
    "PEOE_VSA13",
    "PEOE_VSA14",
    "PEOE_VSA2",
    "PEOE_VSA3",
    "PEOE_VSA4",
    "PEOE_VSA5",
    "PEOE_VSA6",
    "PEOE_VSA7",
    "PEOE_VSA8",
    "PEOE_VSA9",
    "SMR_VSA1",
    "SMR_VSA10",
    "SMR_VSA2",
    "SMR_VSA3",
    "SMR_VSA4",
    "SMR_VSA5",
    "SMR_VSA6",
    "SMR_VSA7",
    "SMR_VSA8",
    "SMR_VSA9",
    "SlogP_VSA1",
    "SlogP_VSA10",
    "SlogP_VSA11",
    "SlogP_VSA12",
    "SlogP_VSA2",
    "SlogP_VSA3",
    "SlogP_VSA4",
    "SlogP_VSA5",
    "SlogP_VSA6",
    "SlogP_VSA7",
    "SlogP_VSA8",
    "SlogP_VSA9",
    "TPSA",
    "EState_VSA1",
    "EState_VSA10",
    "EState_VSA11",
    "EState_VSA2",
    "EState_VSA3",
    "EState_VSA4",
    "EState_VSA5",
    "EState_VSA6",
    "EState_VSA7",
    "EState_VSA8",
    "EState_VSA9",
    "VSA_EState1",
    "VSA_EState10",
    "VSA_EState2",
    "VSA_EState3",
    "VSA_EState4",
    "VSA_EState5",
    "VSA_EState6",
    "VSA_EState7",
    "VSA_EState8",
    "VSA_EState9",
    "FractionCSP3",
    "HeavyAtomCount",
    "NHOHCount",
    "NOCount",
    "NumAliphaticCarbocycles",
    "NumAliphaticHeterocycles",
    "NumAliphaticRings",
    "NumAmideBonds",
    "NumAromaticCarbocycles",
    "NumAromaticHeterocycles",
    "NumAromaticRings",
    "NumAtomStereoCenters",
    "NumBridgeheadAtoms",
    "NumHAcceptors",
    "NumHDonors",
    "NumHeteroatoms",
    "NumHeterocycles",
    "NumRotatableBonds",
    "NumSaturatedCarbocycles",
    "NumSaturatedHeterocycles",
    "NumSaturatedRings",
    "NumSpiroAtoms",
    "NumUnspecifiedAtomStereoCenters",
    "Phi",
    "RingCount",
    "MolLogP",
    "MolMR",
};

const std::vector<std::string> &allNames() {
  static const std::vector<std::string> names = [] {
    auto res = descriptorNames;
    // _descList sorts the fragment descriptors by name
    auto frags = Descriptors::getFragmentDescriptorNames();
    std::sort(frags.begin(), frags.end());
    res.insert(res.end(), frags.begin(), frags.end());
    return res;
  }();
  return names;
}

using ValueMap = std::unordered_map<std::string, double>;

// Runs one group of descriptor calculations; if it throws, the descriptors
// it would have set stay missing and normalize to 0.0 (as in descriptastorus)
void tryGroup(const std::function<void()> &calc) {
  try {
    calc();
  } catch (...) {
  }
}

void setVector(ValueMap &vals, const std::string &prefix,
               const std::vector<double> &values, unsigned int firstIdx = 1) {
  for (unsigned int i = 0; i < values.size(); ++i) {
    vals[prefix + std::to_string(i + firstIdx)] = values[i];
  }
}

ValueMap calcAll(const ROMol &mol) {
  using namespace Descriptors;
  ValueMap v;
  tryGroup([&] {
    v["MaxAbsEStateIndex"] = calcMaxAbsEStateIndex(mol);
    v["MaxEStateIndex"] = calcMaxEStateIndex(mol);
    v["MinAbsEStateIndex"] = calcMinAbsEStateIndex(mol);
    v["MinEStateIndex"] = calcMinEStateIndex(mol);
  });
  tryGroup([&] { v["qed"] = calcQED(mol); });
  tryGroup([&] { v["SPS"] = calcSPS(mol); });
  tryGroup([&] {
    v["MolWt"] = calcAMW(mol);
    v["HeavyAtomMolWt"] = calcHeavyAtomMolWt(mol);
    v["ExactMolWt"] = calcExactMW(mol);
    v["NumValenceElectrons"] = calcNumValenceElectrons(mol);
    v["NumRadicalElectrons"] = calcNumRadicalElectrons(mol);
  });
  tryGroup([&] {
    v["MaxPartialCharge"] = calcMaxPartialCharge(mol);
    v["MinPartialCharge"] = calcMinPartialCharge(mol);
    v["MaxAbsPartialCharge"] = calcMaxAbsPartialCharge(mol);
    v["MinAbsPartialCharge"] = calcMinAbsPartialCharge(mol);
  });
  for (unsigned int radius = 1; radius <= 3; ++radius) {
    tryGroup([&] {
      v["FpDensityMorgan" + std::to_string(radius)] =
          calcFpDensityMorgan(mol, radius);
    });
  }
  tryGroup([&] {
    static const std::vector<std::string> bcutNames = {
        "BCUT2D_MWHI",   "BCUT2D_MWLOW",   "BCUT2D_CHGHI", "BCUT2D_CHGLO",
        "BCUT2D_LOGPHI", "BCUT2D_LOGPLOW", "BCUT2D_MRHI",  "BCUT2D_MRLOW"};
    auto bcut = BCUT2D(mol);
    for (size_t i = 0; i < bcut.size() && i < bcutNames.size(); ++i) {
      v[bcutNames[i]] = bcut[i];
    }
  });
  tryGroup([&] { v["AvgIpc"] = calcAvgIpc(mol); });
  tryGroup([&] { v["BalabanJ"] = calcBalabanJ(mol); });
  tryGroup([&] { v["BertzCT"] = calcBertzCT(mol); });
  tryGroup([&] {
    v["Chi0"] = calcChi0(mol);
    v["Chi1"] = calcChi1(mol);
    v["Chi0n"] = calcChi0n(mol);
    v["Chi0v"] = calcChi0v(mol);
    v["Chi1n"] = calcChi1n(mol);
    v["Chi1v"] = calcChi1v(mol);
    v["Chi2n"] = calcChi2n(mol);
    v["Chi2v"] = calcChi2v(mol);
    v["Chi3n"] = calcChi3n(mol);
    v["Chi3v"] = calcChi3v(mol);
    v["Chi4n"] = calcChi4n(mol);
    v["Chi4v"] = calcChi4v(mol);
    v["HallKierAlpha"] = calcHallKierAlpha(mol);
    v["Kappa1"] = calcKappa1(mol);
    v["Kappa2"] = calcKappa2(mol);
    v["Kappa3"] = calcKappa3(mol);
  });
  tryGroup([&] { v["Ipc"] = calcIpc(mol); });
  tryGroup([&] {
    v["LabuteASA"] = calcLabuteASA(mol);
    v["TPSA"] = calcTPSA(mol);
  });
  tryGroup([&] { setVector(v, "PEOE_VSA", calcPEOE_VSA(mol)); });
  tryGroup([&] { setVector(v, "SMR_VSA", calcSMR_VSA(mol)); });
  tryGroup([&] { setVector(v, "SlogP_VSA", calcSlogP_VSA(mol)); });
  tryGroup([&] { setVector(v, "EState_VSA", calcEState_VSA(mol)); });
  tryGroup([&] { setVector(v, "VSA_EState", calcVSA_EState(mol)); });
  tryGroup([&] {
    v["FractionCSP3"] = calcFractionCSP3(mol);
    v["HeavyAtomCount"] = mol.getNumHeavyAtoms();
    v["NHOHCount"] = calcLipinskiHBD(mol);
    v["NOCount"] = calcLipinskiHBA(mol);
    v["NumAliphaticCarbocycles"] = calcNumAliphaticCarbocycles(mol);
    v["NumAliphaticHeterocycles"] = calcNumAliphaticHeterocycles(mol);
    v["NumAliphaticRings"] = calcNumAliphaticRings(mol);
    v["NumAmideBonds"] = calcNumAmideBonds(mol);
    v["NumAromaticCarbocycles"] = calcNumAromaticCarbocycles(mol);
    v["NumAromaticHeterocycles"] = calcNumAromaticHeterocycles(mol);
    v["NumAromaticRings"] = calcNumAromaticRings(mol);
    v["NumBridgeheadAtoms"] = calcNumBridgeheadAtoms(mol);
    v["NumHAcceptors"] = calcNumHBA(mol);
    v["NumHDonors"] = calcNumHBD(mol);
    v["NumHeteroatoms"] = calcNumHeteroatoms(mol);
    v["NumHeterocycles"] = calcNumHeterocycles(mol);
    v["NumRotatableBonds"] = calcNumRotatableBonds(mol);
    v["NumSaturatedCarbocycles"] = calcNumSaturatedCarbocycles(mol);
    v["NumSaturatedHeterocycles"] = calcNumSaturatedHeterocycles(mol);
    v["NumSaturatedRings"] = calcNumSaturatedRings(mol);
    v["NumSpiroAtoms"] = calcNumSpiroAtoms(mol);
    v["Phi"] = calcPhi(mol);
    v["RingCount"] = calcNumRings(mol);
  });
  tryGroup([&] {
    v["NumAtomStereoCenters"] = numAtomStereoCenters(mol);
    v["NumUnspecifiedAtomStereoCenters"] = numUnspecifiedAtomStereoCenters(mol);
  });
  tryGroup([&] {
    double logp, mr;
    calcCrippenDescriptors(mol, logp, mr);
    v["MolLogP"] = logp;
    v["MolMR"] = mr;
  });
  tryGroup([&] {
    const auto &names = getFragmentDescriptorNames();
    auto counts = calcFragmentDescriptors(mol);
    for (size_t i = 0; i < names.size() && i < counts.size(); ++i) {
      v[names[i]] = counts[i];
    }
  });
  return v;
}
}  // namespace

const std::vector<std::string> &getNormalizedDescriptorNames() {
  return allNames();
}

std::vector<double> calcDescriptorValues(const ROMol &mol) {
  auto vals = calcAll(mol);
  const auto &names = allNames();
  std::vector<double> res(names.size(),
                          std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < names.size(); ++i) {
    auto it = vals.find(names[i]);
    if (it != vals.end()) {
      res[i] = it->second;
    }
  }
  return res;
}

int getNormalizedDescriptorIndex(const std::string &name) {
  static const std::unordered_map<std::string, int> positions = [] {
    std::unordered_map<std::string, int> res;
    const auto &names = allNames();
    for (size_t i = 0; i < names.size(); ++i) {
      res[names[i]] = static_cast<int>(i);
    }
    return res;
  }();
  auto it = positions.find(name);
  return it == positions.end() ? -1 : it->second;
}

std::vector<double> calcNormalizedDescriptors(const ROMol &mol,
                                              const CDFTableSet &tables) {
  return tables.normalizeDescriptors(calcDescriptorValues(mol));
}

std::vector<double> calcNormalizedDescriptors(const ROMol &mol) {
  return calcNormalizedDescriptors(mol, getDefaultTables());
}

}  // namespace NormalizedDescriptors
}  // namespace RDKit
