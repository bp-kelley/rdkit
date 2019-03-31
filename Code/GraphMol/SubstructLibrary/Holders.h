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
#ifndef RDK_SUBSTRUCT_HOLDERS_LIBRARY
#define RDK_SUBSTRUCT_HOLDERS_LIBRARY
#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>

namespace RDKit {

template<class T>
class RDKIT_SUBSTRUCTLIBRARY_EXPORT HolderBase {
protected:
  std::vector<T> datavect;
  typedef std::vector<T> VT;
  
public:
        std::vector<T> &getData() { return datavect; }
  const std::vector<T> &getData() const { return datavect; }

  HolderBase() : datavect() {}
  HolderBase(VT vect) : datavect(vect) {}
  VT filter(const std::vector<unsigned int> &indices) const {
    VT vect;
    vect.reserve(indices.size());
    for(auto idx:indices) {
      if(idx >= datavect.size()) throw IndexErrorException(idx);
      vect.push_back(datavect[idx]);
    }
    return vect; 
  }

  //! remove the specified indices from the data vector
  VT remove(const std::vector<unsigned int> &indices) const {
    std::set<unsigned int> temp_indices(indices.begin(), indices.end());
    VT vect;
    for(size_t idx=0; idx<datavect.size(); ++idx) {
      if (temp_indices.find(idx) == temp_indices.end())
        vect.push_back( datavect[idx] );
    }
    return vect;
  }
  
  //! Keep the indices in indices, and throw out the rest
  void filter_indices(const std::vector<unsigned int> &indices) {
    datavect.swap(filter(indices));
  }

  //! remove the specified indices from the data vector
  void remove_indices(const std::vector<unsigned int> &indices) {
    datavect.swap(remove(indices));
  }
};
}
#endif
