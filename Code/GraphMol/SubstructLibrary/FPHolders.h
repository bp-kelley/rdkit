//  Copyright (c) 2017-2019, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef RDKIT_SUBSTRUCT_PATTERN_HOLDER_LIBRARY
#define RDKIT_SUBSTRUCT_PATTERN_HOLDER_LIBRARY
#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include "Holders.h"
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>

namespace RDKit {
  
//! Base FPI for the fingerprinter used to rule out impossible matches
class RDKIT_SUBSTRUCTLIBRARY_EXPORT FPHolderBase : public HolderBase<ExplicitBitVect*> {
 public:
  FPHolderBase() : HolderBase<ExplicitBitVect*>() {}
  FPHolderBase(const std::vector<ExplicitBitVect*> &vect) : HolderBase<ExplicitBitVect*>(vect) {}
  virtual ~FPHolderBase() {
    for (size_t i = 0; i < datavect.size(); ++i) delete datavect[i];
  }

  //! Adds a molecule to the fingerprinter
  unsigned int addMol(const ROMol &m) {
    datavect.push_back(makeFingerprint(m));
    return rdcast<unsigned int>(datavect.size() - 1);
  }

  //! Adds a raw bit vector to the fingerprinter
  unsigned int addFingerprint(const ExplicitBitVect &v) {
    datavect.push_back(new ExplicitBitVect(v));
    return rdcast<unsigned int>(datavect.size() - 1);
  }

  //! Return false if a substructure search can never match the molecule
  bool passesFilter(unsigned int idx, const ExplicitBitVect &query) const {
    if (idx >= datavect.size()) throw IndexErrorException(idx);

    return AllProbeBitsMatch(query, *datavect[idx]);
  }

  //! Get the bit vector at the specified index (throws IndexError if out of
  //! range)
  const ExplicitBitVect &getFingerprint(unsigned int idx) const {
    if (idx >= datavect.size()) throw IndexErrorException(idx);
    return *datavect[idx];
  }

  //! make the query vector
  //!  Caller owns the vector!
  virtual ExplicitBitVect *makeFingerprint(const ROMol &m) const = 0;

  std::vector<ExplicitBitVect *> &getFingerprints() { return getData(); }
  const std::vector<ExplicitBitVect *> &getFingerprints() const { return getData(); }

  virtual boost::shared_ptr<FPHolderBase> filter_holder(const std::vector<unsigned int> &indices) const = 0;

  virtual boost::shared_ptr<FPHolderBase> filter_remove_holder(const std::vector<unsigned int> &indices) const = 0;
};

//! Uses the pattern fingerprinter to rule out matches
class RDKIT_SUBSTRUCTLIBRARY_EXPORT PatternHolder : public FPHolderBase {
 public:
  //! Caller owns the vector!
  PatternHolder() : FPHolderBase() {}
  PatternHolder(const std::vector<ExplicitBitVect*> &vect) : FPHolderBase(vect) {}
  virtual ExplicitBitVect *makeFingerprint(const ROMol &m) const {
    return PatternFingerprintMol(m, 2048);
  }

  boost::shared_ptr<FPHolderBase> filter_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<FPHolderBase>(new PatternHolder(filter(indices)));
  }

  boost::shared_ptr<FPHolderBase> filter_remove_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<FPHolderBase>(new PatternHolder(remove(indices)));
  }
};

}

#endif
