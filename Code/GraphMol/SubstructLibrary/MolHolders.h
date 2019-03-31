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
#ifndef RDK_SUBSTRUCT_MOLHOLDERS_LIBRARY
#define RDK_SUBSTRUCT_MOLHOLDERS_LIBRARY
#include <RDGeneral/export.h>
#include "Holders.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {


//! Base class API for holding molecules to substructure search.
/*!
  This is an API that hides the implementation details used for
  indexing molecules for substructure searching.  It simply
  provides an API for adding and getting molecules from a set.
 */
class RDKIT_SUBSTRUCTLIBRARY_EXPORT MolHolderBase  {
 public:
  virtual ~MolHolderBase() {}

  //! Add a new molecule to the substructure search library
  //!  Returns the molecules index in the library
  virtual unsigned int addMol(const ROMol &m) = 0;

  // implementations should throw IndexError on out of range
  virtual boost::shared_ptr<ROMol> getMol(unsigned int) const = 0;

  //! Get the current library size
  virtual unsigned int size() const = 0;

  virtual boost::shared_ptr<MolHolderBase> filter_holder(
      const std::vector<unsigned int> &indices) const = 0;

  virtual boost::shared_ptr<MolHolderBase> filter_remove_holder(
      const std::vector<unsigned int> &indices) const = 0;
  
};

//! Concrete class that holds molecules in memory
/*!
    This is currently one of the faster implementations.
    However it is very memory intensive.
*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT MolHolder :
    public MolHolderBase,
    public HolderBase<boost::shared_ptr<ROMol>>
{
public:
  MolHolder() : MolHolderBase(), HolderBase<boost::shared_ptr<ROMol>>() {}
   MolHolder(const std::vector<boost::shared_ptr<ROMol>> &vect) : MolHolderBase(),
      HolderBase<boost::shared_ptr<ROMol>>(vect) {}

  virtual unsigned int addMol(const ROMol &m) {
    datavect.push_back(boost::make_shared<ROMol>(m));
    return size() - 1;
  }

  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    if (idx >= datavect.size()) throw IndexErrorException(idx);
    return datavect[idx];
  }

  virtual unsigned int size() const {
    return rdcast<unsigned int>(datavect.size());
  }

  boost::shared_ptr<MolHolderBase> filter_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<MolHolderBase>(new MolHolder(filter(indices)));
  }

  boost::shared_ptr<MolHolderBase> filter_remove_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<MolHolderBase>(new MolHolder(remove(indices)));
  }

        std::vector<boost::shared_ptr<ROMol>> & getMols() { return getData(); }
  const std::vector<boost::shared_ptr<ROMol>> & getMols() const { return getData(); }

};

//! Concrete class that holds binary cached molecules in memory
/*!
  This implementation uses quite a bit less memory than the
  non cached implementation.  However, due to the reduced speed
  it should be used in conjunction with a pattern fingerprinter.

  See RDKit::FPHolder
*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT CachedMolHolder :
    public MolHolderBase,
    public HolderBase<std::string>
{
 public:
  CachedMolHolder() : MolHolderBase(), HolderBase<std::string>() {}
  CachedMolHolder(const std::vector<std::string> &vect) : MolHolderBase(), HolderBase<std::string>(vect) {}
  
  virtual unsigned int addMol(const ROMol &m) {
    datavect.push_back(std::string());
    MolPickler::pickleMol(m, datavect.back());
    return size() - 1;
  }

  //! Adds a pickled binary molecule, no validity checking of the input
  //!  is done.
  unsigned int addBinary(const std::string &pickle) {
    datavect.push_back(pickle);
    return size() - 1;
  }

  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    if (idx >= datavect.size()) throw IndexErrorException(idx);
    boost::shared_ptr<ROMol> mol(new ROMol);
    MolPickler::molFromPickle(datavect[idx], mol.get());
    return mol;
  }

  virtual unsigned int size() const {
    return rdcast<unsigned int>(datavect.size());
  }

  boost::shared_ptr<MolHolderBase> filter_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<MolHolderBase>(new CachedMolHolder(filter(indices)));
  }

  boost::shared_ptr<MolHolderBase> filter_remove_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<MolHolderBase>(new CachedMolHolder(remove(indices)));
  }

        std::vector<std::string> &getMols() { return getData(); }
  const std::vector<std::string> &getMols() const { return getData(); }
};

//! Concrete class that holds smiles strings in memory
/*!
    This implementation uses quite a bit less memory than the
    cached binary or uncached implementation.  However, due to the
    reduced speed it should be used in conjunction with a pattern
    fingerprinter.

    See RDKit::FPHolder
*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT CachedSmilesMolHolder :
    public MolHolderBase,
    public HolderBase<std::string>
{
  std::vector<std::string> datavect;

 public:
  CachedSmilesMolHolder() : MolHolderBase(), HolderBase<std::string>() {}
  CachedSmilesMolHolder(const std::vector<std::string> &vect) : MolHolderBase(), HolderBase<std::string>(vect) {}
  
  virtual unsigned int addMol(const ROMol &m) {
    bool doIsomericSmiles = true;
    datavect.push_back(MolToSmiles(m, doIsomericSmiles));
    return size() - 1;
  }

  //! Add a smiles to the dataset, no validation is done
  //! to the inputs.
  unsigned int addSmiles(const std::string &smiles) {
    datavect.push_back(smiles);
    return size() - 1;
  }

  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    if (idx >= datavect.size()) throw IndexErrorException(idx);

    boost::shared_ptr<ROMol> mol(SmilesToMol(datavect[idx]));
    return mol;
  }

  virtual unsigned int size() const {
    return rdcast<unsigned int>(datavect.size());
  }

  boost::shared_ptr<MolHolderBase> filter_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<MolHolderBase>(new CachedSmilesMolHolder(filter(indices)));
  }

  boost::shared_ptr<MolHolderBase> filter_remove_holder(
      const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<MolHolderBase>(new CachedSmilesMolHolder(remove(indices)));
  }
  
        std::vector<std::string> &getMols() { return datavect; }
  const std::vector<std::string> &getMols() const { return datavect; }
};

//! Concrete class that holds trusted smiles strings in memory
/*!
    A trusted smiles is essentially a smiles string that
    RDKit has generated.  This indicates that fewer
    sanitization steps are required.  See
    http://rdkit.blogspot.com/2016/09/avoiding-unnecessary-work-and.html

    This implementation uses quite a bit less memory than the
    cached binary or uncached implementation.  However, due to the
    reduced speed it should be used in conjunction with a pattern
    fingerprinter.

    See RDKit::FPHolder
*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT CachedTrustedSmilesMolHolder :
    public MolHolderBase,
    public HolderBase<std::string>
{
 public:
  CachedTrustedSmilesMolHolder() : MolHolderBase(), HolderBase<std::string>() {}
  CachedTrustedSmilesMolHolder(const std::vector<std::string> &vect) : MolHolderBase(), HolderBase<std::string>(vect) {}

  virtual unsigned int addMol(const ROMol &m) {
    bool doIsomericSmiles = true;
    datavect.push_back(MolToSmiles(m, doIsomericSmiles));
    return size() - 1;
  }

  //! Add a smiles to the dataset, no validation is done
  //! to the inputs.
  unsigned int addSmiles(const std::string &smiles) {
    datavect.push_back(smiles);
    return size() - 1;
  }

  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    if (idx >= datavect.size()) throw IndexErrorException(idx);

    RWMol *m = SmilesToMol(datavect[idx], 0, false);
    m->updatePropertyCache();
    return boost::shared_ptr<ROMol>(m);
  }

  virtual unsigned int size() const {
    return rdcast<unsigned int>(datavect.size());
  }

  boost::shared_ptr<MolHolderBase> filter_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<MolHolderBase>(new CachedTrustedSmilesMolHolder(filter(indices)));
  }

  boost::shared_ptr<MolHolderBase> filter_remove_holder(
      const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<MolHolderBase>(new CachedTrustedSmilesMolHolder(remove(indices)));
  }

        std::vector<std::string> &getMols() { return datavect; }
  const std::vector<std::string> &getMols() const { return datavect; }
};

}

#endif
