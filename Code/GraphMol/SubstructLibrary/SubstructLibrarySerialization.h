//  Copyright (c) 2019, Novartis Institutes for BioMedical Research Inc.
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
// n.b. must be included at the END of SubstructLibrary.h
#ifndef RDK_SUBSTRUCT_LIBRARY_SERIALIZATION
#define RDK_SUBSTRUCT_LIBRARY_SERIALIZATION

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>



BOOST_SERIALIZATION_ASSUME_ABSTRACT(RDKit::MolHolderBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(RDKit::FPHolderBase)



namespace boost {
namespace serialization {

template <class Archive>
inline bool writeRDValue(Archive &ar, const RDKit::RDValue &value) {
  switch (value.getTag()) {
    case RDKit::RDTypeTag::StringTag:
      ar << RDKit::DTags::StringTag;
      ar << RDKit::rdvalue_cast<std::string>(value);
      break;
    case RDKit::RDTypeTag::IntTag:
      ar << RDKit::DTags::IntTag;
      ar << RDKit::rdvalue_cast<int>(value);
      break;
    case RDKit::RDTypeTag::UnsignedIntTag:
      ar << RDKit::DTags::UnsignedIntTag;
      ar << RDKit::rdvalue_cast<unsigned int>(value);
      break;
    case RDKit::RDTypeTag::BoolTag:
      ar << RDKit::DTags::BoolTag;
      ar << RDKit::rdvalue_cast<bool>(value);
      break;
    case RDKit::RDTypeTag::FloatTag:
      ar << RDKit::DTags::FloatTag;
      ar << RDKit::rdvalue_cast<float>(value);
      break;
    case RDKit::RDTypeTag::DoubleTag:
      ar << RDKit::DTags::DoubleTag;
      ar << RDKit::rdvalue_cast<double>(value);
      break;

    case RDKit::RDTypeTag::VecStringTag:
      ar << RDKit::DTags::VecStringTag;
      ar << RDKit::rdvalue_cast<std::vector<std::string>>(value);
      break;
    case RDKit::RDTypeTag::VecDoubleTag:
      ar << RDKit::DTags::VecDoubleTag;
      ar << RDKit::rdvalue_cast<std::vector<double>>(value);
      break;
    case RDKit::RDTypeTag::VecFloatTag:
      ar << RDKit::DTags::VecFloatTag;
      ar << RDKit::rdvalue_cast<std::vector<float>>(value);
      break;
    case RDKit::RDTypeTag::VecIntTag:
      ar << RDKit::DTags::VecIntTag;
      ar << RDKit::rdvalue_cast<std::vector<int>>(value);
      break;
    case RDKit::RDTypeTag::VecUnsignedIntTag:
      ar << RDKit::DTags::VecUIntTag;
      ar << RDKit::rdvalue_cast<std::vector<unsigned int>>(value);
      break;
    default:
      throw ValueErrorException("Unable to write RDValue to pickle");
      return false;
  }
  return true;
}

// This should be refactored into StreamOps.h, but that deals directly
//  with dict serialization not RDValue
template<class Archive, class T>
void _readRDValue(Archive &ar, RDKit::RDValue &value) {
  T v;
  ar >> v;
  value = v;
}

template<class Archive>
inline bool readRDValue(Archive &ar, RDKit::RDValue &value) {
  unsigned char type;
  ar >> type;
  switch(type) {
    case RDKit::DTags::IntTag: _readRDValue<Archive, int>(ar, value); break;
    case RDKit::DTags::UnsignedIntTag: _readRDValue<Archive, unsigned int>(ar, value); break;
    case RDKit::DTags::BoolTag: _readRDValue<Archive, bool>(ar, value); break;
    case RDKit::DTags::FloatTag: _readRDValue<Archive, float>(ar, value); break;
    case RDKit::DTags::DoubleTag: _readRDValue<Archive, double>(ar, value); break;

    case RDKit::DTags::StringTag:
      _readRDValue<Archive, std::string>(ar, value);
      break;
    case RDKit::DTags::VecStringTag:
      _readRDValue<Archive, std::vector<std::string>>(ar, value);
      break;
    case RDKit::DTags::VecIntTag:
      _readRDValue<Archive, std::vector<int>>(ar, value);
      break;
    case RDKit::DTags::VecUIntTag:
      _readRDValue<Archive, std::vector<unsigned int>>(ar, value);
      break;
    case RDKit::DTags::VecFloatTag:
      _readRDValue<Archive, std::vector<float>>(ar, value);
      break;
    case RDKit::DTags::VecDoubleTag:
      _readRDValue<Archive, std::vector<double>>(ar, value);
      break;

    default:
      throw ValueErrorException("Unable to read RDValue from pickle");
      return false;
  }
  return true;
}

template <class Archive>
void serialize(Archive &ar, RDKit::MolHolderBase &, const unsigned int version) {
  RDUNUSED_PARAM(version);
  RDUNUSED_PARAM(ar);
}

template <class Archive>
void save(Archive &ar, const RDKit::RDValue& value,
          const unsigned int version) {
  RDUNUSED_PARAM(version);
  writeRDValue(ar, value);
}
  
template <class Archive>
void load(Archive &ar, RDKit::RDValue& value,
          const unsigned int version) {
  RDUNUSED_PARAM(version);
  readRDValue(ar, value);
}

template <class Archive>
void save(Archive &ar, const RDKit::MolHolder& molholder,
          const unsigned int version) {
  RDUNUSED_PARAM(version);
  ar &boost::serialization::base_object<RDKit::MolHolderBase>(molholder);

  std::int64_t pkl_count = molholder.getMols().size();
  ar & pkl_count;
  
  for(auto &mol: molholder.getMols()) {
    std::string pkl;
    RDKit::MolPickler::pickleMol(*mol.get(), pkl);
    ar << pkl;
  }
}
  
template <class Archive>
void load(Archive &ar, RDKit::MolHolder& molholder,
          const unsigned int version) {
  RDUNUSED_PARAM(version);
  ar &boost::serialization::base_object<RDKit::MolHolderBase>(molholder);
  
  std::vector<boost::shared_ptr<RDKit::ROMol>> &mols = molholder.getMols();
  mols.clear();
  
  std::int64_t pkl_count = -1;
  ar & pkl_count;

  for(std::int64_t i = 0; i<pkl_count; ++i) {
    std::string pkl;
    ar >> pkl;
    mols.push_back(boost::make_shared<RDKit::ROMol>(pkl));
  }
}

template<class Archive, class MolHolder>
void serialize_strings(Archive &ar, MolHolder &molholder,
                       const unsigned int version) {
  RDUNUSED_PARAM(version);
  ar &boost::serialization::base_object<RDKit::MolHolderBase>(molholder);
  ar &molholder.getMols();
}

template <class Archive>
void serialize(Archive &ar, RDKit::CachedMolHolder &molholder,
               const unsigned int version) {
  serialize_strings(ar, molholder, version);
}

template <class Archive>
void serialize(Archive &ar, RDKit::CachedSmilesMolHolder &molholder,
               const unsigned int version) {
  serialize_strings(ar, molholder, version);  
}

template <class Archive>
void serialize(Archive &ar, RDKit::CachedTrustedSmilesMolHolder &molholder,
               const unsigned int version) {
  serialize_strings(ar, molholder, version);
}

template <class Archive>
void save(Archive &ar, const RDKit::FPHolderBase &fpholder,
          const unsigned int version) {
  RDUNUSED_PARAM(version);
  std::vector<std::string> pickles;
  for(auto &fp: fpholder.getFingerprints()) {
    pickles.push_back(fp->toString());
  }
  ar &pickles;
}
  
template <class Archive>
void load(Archive &ar, RDKit::FPHolderBase &fpholder,
          const unsigned int version) {
  RDUNUSED_PARAM(version);
  std::vector<std::string> pickles;
  std::vector<ExplicitBitVect *> &fps = fpholder.getFingerprints();
  
  ar &pickles;
  for (size_t i = 0; i < fps.size(); ++i) delete fps[i];
  fps.clear();
  
  for(auto &pkl: pickles) {
    fps.push_back( new ExplicitBitVect(pkl) );
  }
}

template <class Archive>
void serialize(Archive &ar, RDKit::PatternHolder &pattern_holder,
               const unsigned int version) {
  RDUNUSED_PARAM(version);
  ar &boost::serialization::base_object<RDKit::FPHolderBase>(pattern_holder);
}

template <class Archive>
void serialize(Archive &ar, RDKit::KeyHolder &key_holder,
               const unsigned int version) {
  RDUNUSED_PARAM(version);
  std::vector<RDKit::RDValue> &keys = key_holder.getKeys();
  ar & keys;
}

template <class Archive>
void registerSubstructLibraryTypes(Archive &ar) {
  ar.register_type(static_cast<RDKit::MolHolder *>(NULL));
  ar.register_type(static_cast<RDKit::CachedMolHolder *>(NULL));
  ar.register_type(static_cast<RDKit::CachedSmilesMolHolder *>(NULL));
  ar.register_type(static_cast<RDKit::CachedTrustedSmilesMolHolder *>(NULL));
  ar.register_type(static_cast<RDKit::PatternHolder *>(NULL));
}


template <class Archive>
void save(Archive &ar, const RDKit::SubstructLibrary &slib,
          const unsigned int version) {
  RDUNUSED_PARAM(version);
  registerSubstructLibraryTypes(ar);
  ar & slib.getMolHolder();
  ar & slib.getFpHolder();
}

template <class Archive>
void load(Archive &ar, RDKit::SubstructLibrary &slib,
          const unsigned int version) {
  RDUNUSED_PARAM(version);
  registerSubstructLibraryTypes(ar);
  ar & slib.getMolHolder();
  ar & slib.getFpHolder();
}

} // end namespace serialization
} // end namespace boost

BOOST_CLASS_VERSION(RDKit::MolHolder, 1);
BOOST_CLASS_VERSION(RDKit::CachedMolHolder, 1);
BOOST_CLASS_VERSION(RDKit::CachedSmilesMolHolder, 1);
BOOST_CLASS_VERSION(RDKit::CachedTrustedSmilesMolHolder, 1);
BOOST_CLASS_VERSION(RDKit::PatternHolder, 1);
BOOST_CLASS_VERSION(RDKit::SubstructLibrary, 1);

BOOST_SERIALIZATION_SPLIT_FREE(RDKit::MolHolder);
BOOST_SERIALIZATION_SPLIT_FREE(RDKit::FPHolderBase);
BOOST_SERIALIZATION_SPLIT_FREE(RDKit::SubstructLibrary);

#endif
#endif
