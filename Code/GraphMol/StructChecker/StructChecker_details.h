//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#ifndef RD_STRUCTCHECKER_DETAILS_H
#define RD_STRUCTCHECKER_DETAILS_H

#include <string>
#include <vector>
#include "../RDKitBase.h"

/* Example of Usage
1)  StructChecker chk;
    int flags = StructureCheck::checkMolStructure( mol ); // use defaults
 or
2)
    StructureCheck::StructCheckerOptions options;   // use defaults
    // To use external data
    StructureCheck::loadOptionsFromFiles(options, file1, file2, �);
    StructChecker chk(options);

    for( mol in mols ) {
        int flags = StructureCheck::checkMolStructure( mol, &options);
        if (0!=(flags & StructureCheck::StructureFlags::BAD_SET)) {
        // write to error file
        } else if (0!=(flags & StructureCheck::StructureFlags::TRANSFORMED_SET))
{
        // input molecule was transformed
        } else { // flag == NO_CHANGE
        // no change
        }
    }
*/

namespace RDKit {
namespace StructureCheck {

// TypeDefs for translating augmented atom pairs
static const int ANY_CHARGE = 8;
enum RadicalType {
  RT_NONE = 0,
  SINGLET = 1,
  DOUBLET = 2,
  TRIPLET = 3,
  ANY_RADICAL = 0xFF
};

enum AABondType {  // MDL CTFile bond types plus extensions
  BT_NONE = 0,     // means REMOVE Bond
  SINGLE = 1,
  DOUBLE = 2,
  TRIPLE = 3,
  AROMATIC = 4,
  SINGLE_DOUBLE = 5,
  SINGLE_AROMATIC = 6,
  DOUBLE_AROMATIC = 7,
  ANY_BOND = 8,
  ALL_BOND_TYPES = 0xF
};

enum AATopology {
  TP_NONE = 0,  // Don't care
  RING = 1,     // Ring
  CHAIN = 2     // Chain
};

struct Ligand {
  std::string AtomSymbol;
  int Charge;
  RadicalType Radical;
  unsigned SubstitutionCount;  // substitution count 0 = don't care
  AABondType BondType;
  Ligand()
      : Charge(ANY_CHARGE),
        Radical(ANY_RADICAL),
        SubstitutionCount(0),
        BondType(ANY_BOND) {}
};

struct AugmentedAtom {
  std::string AtomSymbol;
  std::string ShortName;
  int Charge;
  RadicalType Radical;
  AATopology Topology;
  std::vector<Ligand> Ligands;

  AugmentedAtom()
      : Charge(ANY_CHARGE), Radical(ANY_RADICAL), Topology(TP_NONE) {}

  AugmentedAtom(const std::string &symbol, const std::string &name, int charge,
                RadicalType radical, AATopology topology)
      : AtomSymbol(symbol),
        ShortName(name),
        Charge(charge),
        Radical(radical),
        Topology(topology) {}
};

struct IncEntry {
  std::string AtomSymbol;
  double LocalInc;
  double AlphaInc;
  double BetaInc;
  double MultInc;

  // Used for logging
  int local_inc_used;
  int alpha_inc_used;
  int beta_inc_used;
  int mult_inc_used;
};

struct PathEntry {
  AugmentedAtom Path;
  double Cond;
  // Used for logging
  int cond_used;
};
}
}
#endif
