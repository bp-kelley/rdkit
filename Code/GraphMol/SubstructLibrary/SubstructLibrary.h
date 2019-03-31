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
#ifndef RDK_SUBSTRUCT_LIBRARY
#define RDK_SUBSTRUCT_LIBRARY
#include <RDGeneral/export.h>
#include "MolHolders.h"
#include "FPHolders.h"
#include "KeyHolder.h"


namespace RDKit {

RDKIT_SUBSTRUCTLIBRARY_EXPORT bool SubstructLibraryCanSerialize();

//! Substructure Search a library of molecules
/*!  This class allows for multithreaded substructure searches os
     large datasets.

     The implementations can use fingerprints to speed up searches
     and have molecules cached as binary forms to reduce memory
     usage.

     basic usage:
     \code
     SubstructLibrary lib;
     lib.addMol(mol);
     std::vector<unsigned int> results = lib.getMatches(query);
     for(std::vector<unsigned int>::const_iterator matchIndex=results.begin();
             matchIndex != results.end();
             ++matchIndex) {
       boost::shared_ptr<ROMol> match = lib.getMol(*matchIndex);
     }
     \endcode
     
     Using different mol holders and pattern fingerprints.

     \code
     boost::shared_ptr<CachedTrustedSmilesMolHolder> molHolder = \
        boost::make_shared<CachedTrustedSmilesMolHolder>();
     boost::shared_ptr<PatternHolder> patternHolder = \
        boost::make_shared<PatternHolder>();

     SubstructLibrary lib(molHolder, patternHolder);
     lib.addMol(mol);
     \endcode

     Cached molecule holders create molecules on demand.  There are currently
     three styles of cached molecules.

       CachedMolHolder: stores molecules in the rdkit binary format.
       CachedSmilesMolHolder: stores molecules in smiles format.
       CachedTrustedSmilesMolHolder: stores molecules in smiles format.

     The CachedTrustedSmilesMolHolder is made to add molecules from
     a trusted source.  This makes the basic assumption that RDKit was
     used to sanitize and canonicalize the smiles string.  In practice
     this is considerably faster than using arbitrary smiles strings since
     certain assumptions can be made.

     When loading from external data, as opposed to using the "addMol" API,
     care must be taken to ensure that the pattern fingerprints and smiles
     are synchronized.

     Each pattern holder has an API point for making its fingerprint.  This
     is useful to ensure that the pattern stored in the database will be
     compatible with the patterns made when analyzing queries.
     
     \code
     boost::shared_ptr<CachedTrustedSmilesMolHolder> molHolder = \
         boost::make_shared<CachedTrustedSmilesMolHolder>();
     boost::shared_ptr<PatternHolder> patternHolder =    \
         boost::make_shared<PatternHolder>();
     
     // the PatternHolder instance is able to make fingerprints.
     //  These, of course, can be read from a file.  For demonstration
     //   purposes we construct them here.
     const std::string trustedSmiles = "c1ccccc1";
     ROMol *m = SmilesToMol(trustedSmiles);
     const ExplicitBitVect *bitVector = patternHolder->makeFingerprint(*m);
     
     // The trusted smiles and bitVector can be read from any source.
     //  This is the fastest way to load a substruct library.
     molHolder->addSmiles( trustedSmiles );
     patternHolder->addFingerprint( *bitVector );
     SubstructLibrary lib(molHolder, patternHolder);
     delete m;
     delete bitVector;
     \endcode
     
*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT SubstructLibrary {
  boost::shared_ptr<MolHolderBase> molholder;
  boost::shared_ptr<FPHolderBase> fpholder;
  boost::shared_ptr<KeyHolder> keyholder;

 public:
  SubstructLibrary()
      : molholder(new MolHolder),
        fpholder(),
        keyholder() {}

  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules)
      : molholder(molecules), fpholder(), keyholder() {}

  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules,
                   boost::shared_ptr<FPHolderBase>  fingerprints)
      : molholder(molecules),
        fpholder(fingerprints),
        keyholder() {}

  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules,
                   boost::shared_ptr<KeyHolder>     keys)
      : molholder(molecules),
        fpholder(),
        keyholder(keys) {}
  
  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules,
                   boost::shared_ptr<FPHolderBase>  fingerprints,
                   boost::shared_ptr<KeyHolder>     keys)
      : molholder(molecules),
        fpholder(fingerprints),
        keyholder(keys) {}
    
  SubstructLibrary(const std::string &pickle)
      : molholder(new MolHolder), fpholder(), keyholder() {
          initFromString(pickle);
  }

  //! Get the underlying molecule holder implementation
  boost::shared_ptr<MolHolderBase> &getMolHolder() {
    return molholder;
  }

  const boost::shared_ptr<MolHolderBase> &getMolHolder() const {
    return molholder;
  }

  //! Get the underlying fingerprint holder implementation
  boost::shared_ptr<FPHolderBase> &getFpHolder() {
    return fpholder;
  }

  //! Get the underlying fingerprint holder implementation
  const boost::shared_ptr<FPHolderBase> &getFpHolder() const {
    return fpholder;
  }

  //! Get the underlying fingerprint holder implementation
  boost::shared_ptr<KeyHolder> getKeyHolder() {
    return keyholder;
  }

  const MolHolderBase &getMolecules() const {
    PRECONDITION(molholder.get(), "Molecule holder NULL in SubstructLibrary");
    return *molholder.get();
  }

  //! Get the underlying fingerprint implementation.
  /*! Throws a value error if no fingerprints have been set */
  FPHolderBase &getFingerprints() {
    if (!fpholder.get())
      throw ValueErrorException("Substruct Library does not have fingerprints");
    return *fpholder.get();
  }

  const FPHolderBase &getFingerprints() const {
    if (!fpholder.get())
      throw ValueErrorException("Substruct Library does not have fingerprints");
    return *fpholder.get();
  }

  //! Add a molecule to the library
  /*!
    \param mol Molecule to add

    returns index for the molecule in the library
  */
  unsigned int addMol(const ROMol &mol);

  //! Add a molecule to the library
  /*!
    \param mol Molecule to add
    \param key key to assign to molecule (requires a KeyHolder)
            If not keyholder is present throws a ValueError

    returns index for the molecule in the library
  */
  unsigned int addMol(const ROMol &mol, RDValue key);
  
  //! Get the matching indices for the query
  /*!
    \param query       Query to match against molecules
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
    \param maxResults  Maximum results to return, -1 means return all [default
    -1]
  */
  std::vector<unsigned int> getMatches(const ROMol &query,
                                       bool recursionPossible = true,
                                       bool useChirality = true,
                                       bool useQueryQueryMatches = false,
                                       int numThreads = -1,
                                       int maxResults = -1);
  //! Get the matching indices for the query between the given indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (non-inclusive) of the search.
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
    \param maxResults  Maximum results to return, -1 means return all [default
    -1]
  */
  std::vector<unsigned int> getMatches(
      const ROMol &query, unsigned int startIdx, unsigned int endIdx,
      bool recursionPossible = true, bool useChirality = true,
      bool useQueryQueryMatches = false, int numThreads = -1,
      int maxResults = -1);

  //! Return the number of matches for the query
  /*!
    \param query       Query to match against molecules
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
  */
  unsigned int countMatches(const ROMol &query, bool recursionPossible = true,
                            bool useChirality = true,
                            bool useQueryQueryMatches = false,
                            int numThreads = -1);
  //! Return the number of matches for the query between the given indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (non-inclusive) of the search.
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
  */
  unsigned int countMatches(const ROMol &query, unsigned int startIdx,
                            unsigned int endIdx, bool recursionPossible = true,
                            bool useChirality = true,
                            bool useQueryQueryMatches = false,
                            int numThreads = -1);

  //! Returns true if any match exists for the query
  /*!
    \param query       Query to match against molecules
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
  */
  bool hasMatch(const ROMol &query, bool recursionPossible = true,
                bool useChirality = true, bool useQueryQueryMatches = false,
                int numThreads = -1);
  //! Returns true if any match exists for the query between the specified
  //! indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (inclusive) of the search.
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
  */
  bool hasMatch(const ROMol &query, unsigned int startIdx, unsigned int endIdx,
                bool recursionPossible = true, bool useChirality = true,
                bool useQueryQueryMatches = false, int numThreads = -1);

  //! Returns the molecule at the given index
  /*!
    \param idx       Index of the molecule in the library
  */
  boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    // expects implementation to throw IndexError if out of range
    PRECONDITION(molholder.get(), "molholder is null in SubstructLibrary");
    return molholder->getMol(idx);
  }

  //! Returns the molecule at the given index
  /*!
    \param idx       Index of the molecule in the library
  */
  boost::shared_ptr<ROMol> operator[](unsigned int idx) {
    // expects implementation to throw IndexError if out of range
    PRECONDITION(molholder.get(), "molholder is null in SubstructLibrary");
    return molholder->getMol(idx);
  }

  //! return the number of molecules in the library
  unsigned int size() const {
    PRECONDITION(molholder.get(), "molholder is null in SubstructLibrary");
    return rdcast<unsigned int>(molholder->size());
  }

  //! Returns the RDValue for key at the given index (assuming the library is indexed with keys)
  /*!
    \param idx       Index of the molecule in the library
  */  
  RDValue getKey(unsigned int idx) const {
    if (!keyholder.get()) throw ValueErrorException("Substruct Library made without keys");
    return keyholder->getKey(idx);
  }

  //! Returns the string value for the key at the given index (assuming the library is indexed with keys)
  //!  If the key can't be converted to a string, a ValueError is thrown
  /*!
    \param idx       Index of the molecule in the library
  */    
  std::string getKeyString(unsigned int idx) const {
    if (!keyholder.get()) throw ValueErrorException("Substruct Library made without keys");
    return keyholder->getKeyString(idx);
  }
  
  //! Return a new substruct library with only the specified indices
  SubstructLibrary filter(const std::vector<unsigned int> &indices) const;
  //! Return a new substructlibrary with the specified indices removed
  SubstructLibrary filter_remove(const std::vector<unsigned int> &indices) const;

  //! Keep the supplied indices and remove the rest
  //!  Note: this will affect any other library using the same Holders
  //!    see filter
  
  void keep_indices(const std::vector<unsigned int> &indices);
  //! Remove the specified indices from the library
  //!  Note: this will affect any other library using the same Holders
  //!     see filter_remove
  void remove_indices(const std::vector<unsigned int> &indices);

  //! serializes (pickles) to a stream
  void toStream(std::ostream &ss) const;
  //! returns a string with a serialized (pickled) representation
  std::string Serialize() const;
  //! initializes from a stream pickle
  void initFromStream(std::istream &ss);
  //! initializes from a string pickle
  void initFromString(const std::string &text);
};
}

#include "SubstructLibrarySerialization.h"
#endif
