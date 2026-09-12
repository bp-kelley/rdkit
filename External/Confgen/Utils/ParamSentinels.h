//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Sentinels for numeric parameters bounded at zero:
//
//      -1  = AUTO      derive a sensible value
//       0  = DISABLED  feature off
//      >0  = use it literally
//
//  Counts are declared SIGNED so -1 is expressible.  `0` must never mean
//  "derive" -- see docs/parameter-conventions.md for why, and for the
//  resolved-parameter dump that goes with it.
//
#ifndef RDKIT_CONFGEN_PARAMSENTINELS_H
#define RDKIT_CONFGEN_PARAMSENTINELS_H

#include <RDGeneral/export.h>

namespace RDKit {

constexpr int AutoI = -1;       //!< "derive it" for counts
constexpr double AutoR = -1.0;  //!< "derive it" for thresholds bounded at zero
constexpr int Disabled = 0;  //!< "feature off" (same for both; spelled out for intent)

//! Any negative value asks to be derived; these parameters are all bounded at
//! zero.
inline bool isAuto(double v) { return v < 0.0; }
inline bool isAuto(int v) { return v < 0; }

//! Resolve a possibly-AUTO parameter.  0 means DISABLED and is preserved as-is.
template <typename T>
inline T resolveAuto(T v, T derived) {
  return v < 0 ? derived : v;
}

}  // namespace RDKit

#endif
