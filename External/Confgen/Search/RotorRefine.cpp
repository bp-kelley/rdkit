//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Search/RotorRefine.h"

#include <algorithm>
#include <cmath>

namespace RDKit {

void refineRotorsInPlace(RotorDriver &drv, std::vector<SearchResult> &out,
                         unsigned int nRefine, double stepDeg,
                         unsigned int nPasses,
                         const RotorDriver::ScoreFn &objective) {
  const size_t nr = drv.numRotors();
  if (nRefine == 0 || out.empty() || nr == 0) {
    return;
  }
  const size_t nDo = std::min<size_t>(nRefine, out.size());
  auto refineScore = [&]() {
    return objective ? objective(drv.positions().data(), drv.numAtoms())
                     : drv.score();
  };

  for (size_t i = 0; i < nDo; ++i) {
    drv.positions() = out[i].coords;
    double bestSc = refineScore();
    if (std::isnan(bestSc)) {
      continue;
    }
    for (unsigned int pass = 0; pass < nPasses; ++pass) {
      const double step = stepDeg / static_cast<double>(1u << pass);
      bool improved = false;
      for (size_t r = 0; r < nr; ++r) {
        const double cur = drv.dihedralDeg(static_cast<unsigned int>(r));
        double bestAng = cur;
        for (const double d : {step, -step}) {
          drv.setDihedral(static_cast<unsigned int>(r), cur + d);
          const double sc = refineScore();
          if (!std::isnan(sc) && sc < bestSc) {
            bestSc = sc;
            bestAng = cur + d;
            improved = true;
          }
        }
        drv.setDihedral(static_cast<unsigned int>(r), bestAng);  // keep winner
      }
      if (!improved) {
        break;
      }
    }
    out[i].coords = drv.positions();
    // keep out[] on the SEARCH's scale even when the descent used another
    // objective, so downstream ranking stays comparable
    out[i].score = objective ? drv.score() : bestSc;
  }

  std::sort(out.begin(), out.end(),
            [](const SearchResult &a, const SearchResult &b) {
              if (std::isnan(a.score)) return false;
              if (std::isnan(b.score)) return true;
              return a.score < b.score;
            });
}

}  // namespace RDKit
