//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_CONFGEN_SEARCHRESULT_H
#define RDKIT_CONFGEN_SEARCHRESULT_H

#include <RDGeneral/export.h>

#include <vector>

namespace RDKit {

//! One conformer produced by a junction-angle search.
/*!
  Every search returns these -- tree, Thompson and systematic alike -- so the
  type is deliberately independent of any one of them.  It was previously nested
  inside the tree search, which made unrelated code include that header just to
  name its own return type.

  Coordinates are a flat, full-molecule buffer; nothing here is a molecule.  The
  searches are coordinate-only by design (no ROMol is rebuilt during a search),
  and these buffers are materialised into conformers once, at the end, by
  FragmentConfGen.
*/
struct RDKIT_FRAGMENTCONFGEN_EXPORT SearchResult {
  std::vector<double> coords;  //!< flat [x0,y0,z0,...] full-molecule buffer
  //! Inter-fragment score (kcal/mol) under the search's scorer.  NOT the full
  //! MMFF energy of the conformer -- the constant intra-fragment term is
  //! excluded, since it does not vary while driving junction angles.
  double score;
};

}  // namespace RDKit

#endif
