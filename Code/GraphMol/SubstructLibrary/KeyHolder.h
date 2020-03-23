//
//  Copyright (C) 2020 Brian P. Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDK_KEYHOLDER
#define RDK_KEYHOLDER
#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <map>
#include <vector>

namespace RDKit {
  
class RDKIT_SUBSTRUCTLIBRARY_EXPORT KeyHolderBase {
public:
 virtual ~KeyHolderBase() {}
 virtual unsigned int size() const = 0;
 
 template<class T>
 unsigned int getIdx(const T&value) const {
   std::string val = boost::lexical_cast<std::string>(value);
   return getIdx(value);
 }
 virtual unsigned int addMol(const ROMol& mol) = 0;
 virtual unsigned int add(const std::string& key) = 0;

 virtual void apply(ROMol &mol, unsigned int idx) const;
 virtual unsigned int getIdx(const std::string &value) const = 0;
 virtual std::string  getKey(unsigned int idx) const = 0;
 virtual void remove(unsigned int idx) = 0;
};

class RDKIT_SUBSTRUCTLIBRARY_EXPORT StringKeyHolder : public KeyHolderBase {
  std::map<std::string, unsigned int> index;
  std::vector<std::string> index_key;
  std::string prop_name;
  
 public:
  StringKeyHolder(const std::string &prop_name="_Name") :
    KeyHolderBase(), index(), index_key(), prop_name(prop_name) {
  }
  virtual ~StringKeyHolder() {
  }
  
  virtual unsigned int size() const {
    return rdcast<unsigned int>(index.size());
   }

  virtual unsigned int addMol(const ROMol& mol) {
    unsigned int idx = index.size();
    std::string key;
    if(prop_name.size() == 0) {
      key = boost::lexical_cast<std::string>(idx);
    }
    else if(!mol.getPropIfPresent(prop_name, key)) {
      throw ValueErrorException(std::string("Property ") + prop_name + " isn't present in molecule");
    }

    return add(key);
  }

  virtual unsigned int add(const std::string &key) {
    if(index.find(key) != index.end()) {
      throw ValueErrorException(key + " already exists");
    }
    unsigned int idx = index.size();
    index[key] = idx;
    index_key.push_back(key);
    return idx;
  }


  virtual void apply(ROMol &m, unsigned int idx) {
    m.SetProp(prop_name, getKey(idx));
  }

  virtual unsigned int getIdx(const std::string &value) const {
    // do a scan of the values
    auto it = index.find(value);
    if (it != index.end()) {
      return it->second;
    }
    throw KeyErrorException(value);
  }
  
  virtual std::string getKey(unsigned int idx) const {
    if (idx > index_key.size()) {
      throw IndexErrorException(idx);
    }
    return index_key[idx];
  }

  void remove(unsigned int idx) {
    std::string key = getKey(idx); // throws
    index.erase(key);
     
    key = index_key.back();
    index[key] = idx;
    index_key[idx] = key;
    index_key.pop_back();
  }  
};

}

#endif
