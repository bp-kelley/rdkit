//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Experimental RigidRotorSearch strategies -- selectable, but not the shipped default.
//
#include "Search/RigidRotorSearchExperimental.h"

#include "Search/RotorDriver.h"
#include "Search/SystematicSearch.h"
#include "Search/ThompsonSamplingSearch.h"
#include "Sampler/TorsionSampler.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <random>
#include <vector>

namespace RDKit {

namespace {
inline double circDiff(double a, double b) {
  double d = std::fmod(a - b + 540.0, 360.0) - 180.0;
  return std::fabs(d);
}
}  // namespace



std::vector<SearchResult> MergedSearch::search(const FragmentZipperContext &ctx,
                                               const RigidRotorSearchParams &sp) {
  SystematicSearch sys;
  ThompsonSamplingSearch ts;
  std::vector<SearchResult> out = sys.search(ctx, sp);
  std::vector<SearchResult> b = ts.search(ctx, sp);
  d_lastBudget = sys.lastBudget() + ts.lastBudget();
  out.insert(out.end(), std::make_move_iterator(b.begin()),
             std::make_move_iterator(b.end()));

  std::sort(out.begin(), out.end(),
            [](const SearchResult &x, const SearchResult &y) {
              if (std::isnan(x.score)) return false;
              if (std::isnan(y.score)) return true;
              return x.score < y.score;
            });
  // The two arms overlap on the easy poses; dedup before the output cap or the
  // merged ensemble is mostly duplicates and the extra coverage is squeezed out.
  const double finalRms =
      resolveAuto(sp.systematic.finalRms, sp.diversityRmsThresh > 0.0
                                              ? sp.diversityRmsThresh
                                              : 0.5);
  if (finalRms > 0.0) {
    ctx.symmetryDedup(out, finalRms);
  }
  return out;
}

}  // namespace RDKit
