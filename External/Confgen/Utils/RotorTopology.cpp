//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Utils/RotorTopology.h"

#include <algorithm>
#include <map>


//*  This code is an attempt to return RMS values for different styles
//*  of rotatable groups.  I.e. long chains are floppy so prune at a higher
//*  RMS and so forth.  Most of this is heuristics and may not really scale
//*  beyond the test datasets.
namespace RDKit {

RotorTopology rotorTopology(const ROMol &mol,
                            const std::vector<unsigned int> &rotorBonds) {
  RotorTopology topo;
  const size_t nb = rotorBonds.size();
  topo.nRotors = static_cast<unsigned int>(nb);
  if (nb == 0 || !mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    return topo;
  }

  // Are two rotors in the same "chain", i.e. between ring systems
  std::vector<size_t> parent(nb);
  for (size_t i = 0; i < nb; ++i) {
    parent[i] = i;
  }
  const std::function<size_t(size_t)> find = [&](size_t x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  for (size_t i = 0; i < nb; ++i) {
    const Bond *bi = mol.getBondWithIdx(rotorBonds[i]);
    for (size_t j = i + 1; j < nb; ++j) {
      const Bond *bj = mol.getBondWithIdx(rotorBonds[j]);
      for (const unsigned int a : {bi->getBeginAtomIdx(), bi->getEndAtomIdx()}) {
        if ((a == bj->getBeginAtomIdx() || a == bj->getEndAtomIdx()) &&
            !mol.getRingInfo()->numAtomRings(a)) {
          parent[find(i)] = find(j);
        }
      }
    }
  }

  std::map<size_t, unsigned int> sizes;
  for (size_t i = 0; i < nb; ++i) {
    sizes[find(i)]++;
  }
  topo.nChains = static_cast<unsigned int>(sizes.size());
  for (const auto &kv : sizes) {
    topo.maxChain = std::max(topo.maxChain, kv.second);
  }
  topo.longestChainFraction =
      static_cast<double>(topo.maxChain) / static_cast<double>(nb);
  return topo;
}

double diversityRmsForRotors(const std::vector<double> &byRotor,
                             unsigned int nRotors, double fallback) {
  if (byRotor.empty()) {
    return fallback;
  }
  const size_t idx = std::min<size_t>(nRotors, byRotor.size() - 1);
  return byRotor[idx];
}

double autoDiversityRms(const RotorTopology &topo) {
  // HEURISTIC: Only floppy molecules whose flexibility is CONCENTRATED in one chain
  if (topo.nRotors > 7 && topo.longestChainFraction > 0.6) {
    return 0.75;
  }
  return 0.5;
}

}  // namespace RDKit
