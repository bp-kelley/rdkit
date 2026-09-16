//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  DEVELOPER Debugging Parameters for logging and fine grained diagnostics
//
#ifndef RDKIT_CONFGEN_DIAGNOSTICSPARAMS_H
#define RDKIT_CONFGEN_DIAGNOSTICSPARAMS_H

#include <RDGeneral/export.h>

namespace RDKit {

struct RDKIT_FRAGMENTCONFGEN_EXPORT DiagnosticsParams {
  bool FRAGLIB_TRACE = false;  //!< Embedder
  bool ASM_PROFILE = false;    //!< Joiner:  whole-process timing
  bool ASM_INCRSTATS = false;  //!< Joiner:  rigid FF score statistics
  bool ASM_SCOREVALIDATE =
      false;  //!< Joiner:  rigid vs full FF checking
  bool ASM_EXACT_GEOM =
      false;  //!< Joiner:  Allows seeding from exact input geometry (i.e. xtal)
              //!< sampler; slated for removal)
  bool SYS_VALIDATE =
      false;  //!< Search/sys: subtree-energy self-check vs full score
  bool FRAGCG_DUMP_PARAMS =
      false;  //!< top:        dump the resolved parameter tree
  bool TS_ARMSTATS =
      false;  //!< Search/ts:  per-rotor arm coverage + posterior entropy,
              //!< bucketed by moving-atom quartile
};

}  // namespace RDKit

#endif
