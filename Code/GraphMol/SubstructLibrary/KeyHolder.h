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
#ifndef RDK_SUBSTRUCT_KEYHOLDER_LIBRARY
#define RDK_SUBSTRUCT_KEYHOLDER_LIBRARY
#include <RDGeneral/export.h>
#include "Holders.h"

namespace RDKit {

class RDKIT_SUBSTRUCTLIBRARY_EXPORT KeyHolder: public HolderBase<RDValue> {
public:

  KeyHolder(HolderBase<RDValue>::VT vect) : HolderBase<RDValue>(vect) {}
  
  unsigned int addKey(RDValue key) {
    datavect.push_back(key);
    return size() - 1;
  }

  const RDValue & getKey(unsigned int idx) const {
    if (idx >= datavect.size()) throw IndexErrorException(idx);
    return datavect[idx];
  }

  std::string getKeyString(unsigned int idx) const {
    if (idx >= datavect.size()) throw IndexErrorException(idx);
    const RDValue & val = datavect[idx];
    std::string res;
    if (rdvalue_tostring(val, res)) {
      return res;
    }
    throw ValueErrorException(
        std::string("Unable to convert Key to string at idx: ") +
        boost::lexical_cast<std::string>(idx));
  }

  unsigned int size() const {
    return rdcast<unsigned int>(datavect.size());
  }

  boost::shared_ptr<KeyHolder> filter_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<KeyHolder>(new KeyHolder(filter(indices)));
  }

  boost::shared_ptr<KeyHolder> filter_remove_holder(const std::vector<unsigned int> &indices) const {
    return boost::shared_ptr<KeyHolder>(new KeyHolder(remove(indices)));
  }
        std::vector<RDValue> & getKeys() { return datavect; }
  const std::vector<RDValue> & getKeys() const { return datavect; }
};

}

#endif
