//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
#ifndef RDKIT_RDVALUE_PTRMAGIC_H
#define RDKIT_RDVALUE_PTRMAGIC_H

#include <stdint.h>
#include <cassert>
#include <boost/any.hpp>
#include "Invariant.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>

#include "LocaleSwitcher.h"

namespace RDKit {

  // Inspired by
  // https://nikic.github.io/2012/02/02/Pointer-magic-for-efficient-dynamic-value-representations.html
// 16 bit storage for value types using Quiet NaN spaces in
//  doubles
// Won't work on Solaris and some other os's as mmaping maps from
// top memory down
// Example check:
//     std::string *pointer = new std::string(v);
//      assert((reinterpret_cast<uint64_t>(pointer) & StringTag) == 0);

//  implementations, need a typedef at compile time to figure this out.

/*
  Encoding for storing other things as a double.  Use
  Quiet NaN
  Quiet NaN: // used to encode types
   F   F    F   1XXX < - X = type bits (first bit is set to one)

  seeeeeee|eeeemmmm|mmmmmmmm|mmmmmmmm|mmmmmmmm|mmmmmmmm|mmmmmmmm|mmmmmmmm
  s1111111|11111ppp|pppppppp|pppppppp|pppppppp|pppppppp|pppppppp|pppppppp
               ^- first mantissa bit 1    everything else is "payload" -^
   ^- exponent bits all 1                 and mustn't be all-zero (as it
  ^- any sign bit                         would be INF then)

  Available 
  8 = 1000 MaxDouble // Not really a tag, is a sentinel
  9 = 1001 Float
  b = 1010 Int32
  a = 1011 Uint32
  C = 1100 <none>
  D = 1101 <none>
  E = 1110 <none>
  F = 1111 PtrTag (look at lower 3 bits for type)
*/ 

      
class RDValue {
private:
  union {
    double doubleBits;
    uint64_t otherBits;
  };
public:
  static const uint64_t NaN       = 0xfff7FFFFFFFFFFFF; // signalling NaN
  static const uint64_t MaxDouble = 0xfff8000000000000; // 
  static const uint64_t FloatTag  = 0xfff9000000000000; // 
  static const uint64_t Int32Tag  = 0xfffa000000000000; // 
  static const uint64_t UInt32Tag = 0xfffb000000000000;

  // PTR Tags use the last 3 bits for typing info
  static const uint64_t PtrTag            = 0xffff000000000000;
  static const uint64_t StringTag         = 0xffff000000000001; //001
  static const uint64_t VecDoubleTag      = 0xffff000000000002; //010
  static const uint64_t VecFloatTag       = 0xffff000000000003; //011
  static const uint64_t VecIntTag         = 0xffff000000000004; //100
  static const uint64_t VecUnsignedIntTag = 0xffff000000000005; //101
  static const uint64_t VecStringTag      = 0xffff000000000006; //110
  static const uint64_t AnyTag            = 0xffff000000000007; //111
  
  static const uint64_t TagMask      = 0xFFFF000000000000;

  static const uint64_t PointerTagMask   = 0xFFFF000000000007;
  static const uint64_t ApplyMask        = 0x0000FFFFFFFFFFFF;
  static const uint64_t ApplyPtrMask     = 0x0000FFFFFFFFFFF8;

  inline RDValue() : doubleBits(0.0) {}
  
  inline RDValue(double number) {
    if (boost::math::isnan(number)) {
      otherBits = NaN;
      assert(boost::math::isnan(doubleBits));
    }
    else
      doubleBits = number;
  }

  inline RDValue(float number) {
    union { int i; float f; } fu;
    fu.f = number;
    otherBits = (((uint64_t)fu.i) & ApplyMask ) | FloatTag;
  }
    
  inline RDValue(int32_t number) {
    otherBits = (((uint64_t)number) & ApplyMask ) | Int32Tag;
    uint64_t tag = getTag();
    assert(tag   == Int32Tag);
    assert(tag  != StringTag);
    assert(tag  != AnyTag);
  }

  inline RDValue(unsigned int number) {
    otherBits = (((uint64_t)number) & ApplyMask ) | UInt32Tag;
  }
  
  inline RDValue(bool number) {
    otherBits = (static_cast<uint64_t>(number) & ApplyMask) | Int32Tag;
  }

  inline RDValue(boost::any *pointer) {
    // ensure that the pointer really is only 48 bit
    assert((reinterpret_cast<uint64_t>(pointer) & AnyTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | AnyTag;
  }
  
  inline RDValue(const boost::any &any) {
    // ensure that the pointer really is only 48 bit
    boost::any *pointer = new boost::any(any);
    assert((reinterpret_cast<uint64_t>(pointer) & AnyTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | AnyTag;
  }

  template <class T>
  inline RDValue(const T&v) {
    boost::any *pointer = new boost::any(v);
    assert((reinterpret_cast<uint64_t>(pointer) & AnyTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | AnyTag;
    uint64_t tag = getTag();
    assert(tag  != Int32Tag);
    assert(tag  != StringTag);
    assert(tag  == AnyTag);

  }
  
  inline RDValue(const std::string &v) {
    std::string *pointer = new std::string(v);
    assert((reinterpret_cast<uint64_t>(pointer) & StringTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | StringTag;
    uint64_t tag = getTag();
    assert(tag  != Int32Tag);
    assert(tag  == StringTag);
    assert(tag  != AnyTag);

  }

  inline RDValue(const std::vector<double> &v) {
    std::vector<double> *pointer = new std::vector<double>(v);
    assert((reinterpret_cast<uint64_t>(pointer) & VecDoubleTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | VecDoubleTag;
  }

  inline RDValue(const std::vector<float> &v) {
    std::vector<float> *pointer = new std::vector<float>(v);
    assert((reinterpret_cast<uint64_t>(pointer) & VecFloatTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | VecFloatTag;
  }
  
  inline RDValue(const std::vector<int> &v) {
    std::vector<int> *pointer = new std::vector<int>(v);
    assert((reinterpret_cast<uint64_t>(pointer) & VecIntTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | VecIntTag;    
  }

  inline RDValue(const std::vector<unsigned int> &v) {
    std::vector<unsigned int> *pointer = new std::vector<unsigned int>(v);
    assert((reinterpret_cast<uint64_t>(pointer) & VecIntTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | VecUnsignedIntTag;    
  }
  
  inline RDValue(const std::vector<std::string> &v) {
    std::vector<std::string> *pointer = new std::vector<std::string>(v);
    assert((reinterpret_cast<uint64_t>(pointer) & VecStringTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | VecStringTag;        
  }
  
  uint64_t getTag() const {
    uint64_t tag = otherBits & TagMask;
    if (tag == PtrTag)
      return otherBits & PointerTagMask;
    return tag;
  }
  
  inline bool isDouble() const {
    return otherBits < MaxDouble || (otherBits & NaN) == NaN;
  }

  inline bool isFloat() const { // check size?  rdcast?
    return getTag() == FloatTag;
  }
    
  inline bool isInt32() const {
    return getTag() == Int32Tag;
  }

  inline bool isUInt32() const {
    return getTag() == UInt32Tag;
  }
  
  inline bool isBool() const {
    return (getTag() == Int32Tag &&
            (getInt32() == 0 || getInt32() == 1));

  }
  
  inline bool isString() const {
    return getTag() == StringTag;
  }
  
  inline bool isVecDouble() const {
    return getTag() == VecDoubleTag;
  }

  inline bool isVecFloat() const {
    return getTag() == VecFloatTag;
  }
  
  inline bool isVecInt() const {
    return getTag() == VecIntTag;
  }

  inline bool isVecUnsignedInt() const {
    return getTag() == VecUnsignedIntTag;
  }
  
  inline bool isVecString() const {
    return getTag() == VecStringTag;
  }
  
  inline bool isAny() const {
    return getTag() == AnyTag;
  }

  /////////////////////////////////////////////////
  // POD getters
  inline double getDouble() const {
    assert(isDouble());

    return doubleBits;
  }

  inline float getFloat() const {
    assert(isDouble());

    return rdcast<float>(doubleBits);
  }
    
  inline int32_t getInt32() const {
    assert(isInt32());

    return static_cast<int32_t>(otherBits & ~Int32Tag);
  }

  inline int32_t getUInt32() const {
    assert(isUInt32());

    return static_cast<uint32_t>(otherBits & ~UInt32Tag);
  }
  
  inline bool getBool() const {
    assert(isBool());
    return getInt32() != 0;
  }
  //////////////////////////////////////////////////
  // variant const interface
  inline const double & asDouble() const {
    assert(isDouble());
    return doubleBits;
  }

  inline const float & asFloat() const {
    assert(isFloat());
    return reinterpret_cast<const float&>(otherBits);// little endian
  }
  
  inline const int32_t & asInt32() const {
    assert(isInt32());
    return reinterpret_cast<const int32_t&>(otherBits);// little endian
  }

  inline const uint32_t &asUInt32() const {
    assert(isUInt32());
    return reinterpret_cast<const uint32_t&>(otherBits);// little endian
  }

  
  inline const bool &asBool() const {
    assert(isBool());
    // does this even work?
    return reinterpret_cast<const bool&>(otherBits);// little endian
  }
  
  inline const std::string & asString() const {
    assert(isString());

    return *reinterpret_cast<std::string *>(otherBits & ~StringTag);
  }

  inline const std::vector<double> & asVecDouble() const {
    assert(isVecDouble());

    return *reinterpret_cast<std::vector<double> *>(otherBits & ~VecDoubleTag);
  }

  inline const std::vector<float> & asVecFloat() const {
    assert(isVecFloat());

    return *reinterpret_cast<std::vector<float> *>(otherBits & ~VecFloatTag);
  }
  
  inline const std::vector<int> & asVecInt() const {
    assert(isVecInt());
    
    return *reinterpret_cast<std::vector<int> *>(otherBits & ~VecIntTag);
  }

  inline const std::vector<unsigned int> & asVecUnsignedInt() const {
    assert(isVecUnsignedInt());
    
    return *reinterpret_cast<std::vector<unsigned int> *>(otherBits &
                                                          ~VecUnsignedIntTag);
  }
  
  inline const std::vector<std::string> & asVecString() const {
    assert(isVecString());
    
    return *reinterpret_cast<std::vector<std::string> *>(otherBits & ~VecStringTag);
  }

  inline boost::any & asAny() const {
    assert(isAny());
    return *reinterpret_cast<boost::any *>(otherBits & ~AnyTag);
  }

  template <class T>
  const T get() const {
    assert(isAny());
    
    return boost::any_cast<T>(reinterpret_cast<boost::any *>(otherBits & ~AnyTag));
  }

  //////////////////////////////////////////////////
  // variant const interface
  inline double & asDouble()  {
    assert(isDouble());
    return doubleBits;
  }

  inline float & asFloat()  {
    assert(isFloat());
    return *(float*)&otherBits; // little endian
  }

  inline int32_t & asInt32()  {
    assert(isInt32());
    return *(int32_t*)&otherBits;// little endian
  }

  inline uint32_t &asUInt32()  {
    assert(isUInt32());
    return *(uint32_t*)&otherBits;// little endian
  }

  
  inline  bool &asBool()  {
    assert(isBool());
    // does this even work?
    return *(bool*)&otherBits;// little endian
  }
  
  inline  std::string & asString()  {
    assert(isString());

    return *reinterpret_cast<std::string *>(otherBits & ~StringTag);
  }

  inline  std::vector<double> & asVecDouble()  {
    assert(isVecDouble());

    return *reinterpret_cast<std::vector<double> *>(otherBits & ~VecDoubleTag);
  }

  inline  std::vector<float> & asVecFloat()  {
    assert(isVecFloat());

    return *reinterpret_cast<std::vector<float> *>(otherBits & ~VecFloatTag);
  }

  inline  std::vector<int> & asVecInt()  {
    assert(isVecInt());
    
    return *reinterpret_cast<std::vector<int> *>(otherBits & ~VecIntTag);
  }

  inline  std::vector<unsigned int> & asVecUnsignedInt()  {
    assert(isVecUnsignedInt());
    
    return *reinterpret_cast<std::vector<unsigned int> *>(
        otherBits & ~VecUnsignedIntTag);
  }

  inline  std::vector<std::string> & asVecString()  {
    assert(isVecString());
    
    return *reinterpret_cast<std::vector<std::string> *>(otherBits & ~VecStringTag);
  }

  inline boost::any & asAny()  {
    assert(isAny());
    return *reinterpret_cast<boost::any *>(otherBits & ~AnyTag);
  }

  template <class T>
   T get()  {
    assert(isAny());
    
    return boost::any_cast<T>(reinterpret_cast<boost::any *>(otherBits & ~AnyTag));
  }
  
  static // Given a type and an RDAnyValue - delete the appropriate structure
  inline void cleanup_rdvalue(RDValue &val) {
    switch(val.getTag()) {
      case RDValue::StringTag:
        delete &val.asString();
        break;
      case RDValue::VecDoubleTag:
        delete &val.asVecDouble();
        break;
      case RDValue::VecFloatTag:
        delete &val.asVecFloat();
        break;
      case RDValue::VecIntTag:
        delete &val.asVecInt();
        break;
      case RDValue::VecUnsignedIntTag:
        delete &val.asVecUnsignedInt();
        break;
      case RDValue::VecStringTag:
        delete &val.asVecString();
        break;
      case RDValue::AnyTag:
        delete &val.asAny();
        break;
      default:
        break;
    }
  }

  RDValue& operator=(double v) {
    cleanup_rdvalue(*this);    
    if (boost::math::isnan(v)) otherBits = NaN;
    else doubleBits = v;
    return *this;
  }
    
  RDValue& operator=(float v) {
    cleanup_rdvalue(*this);
    union { int i; float f; } fu;
    fu.f = v;
    otherBits = (((uint64_t)fu.i) & ApplyMask ) | FloatTag;
    return *this;
  }
    
  RDValue& operator=(int v) {
    cleanup_rdvalue(*this);
    otherBits = (static_cast<int32_t>(v) & ApplyMask) | Int32Tag;
    return *this;
  }

  RDValue& operator=(unsigned int v) {
    cleanup_rdvalue(*this);
    otherBits = (static_cast<uint32_t>(v) & ApplyMask) | UInt32Tag;
    return *this;
  }
  
  RDValue& operator=(bool v) {
    cleanup_rdvalue(*this);
    otherBits = (static_cast<int32_t>(v) & ApplyMask) | Int32Tag;
    return *this;
  }
    
  RDValue& operator=(const std::string & v) {
    cleanup_rdvalue(*this);
    std::string * pointer = new std::string(v);
    otherBits = reinterpret_cast<uint64_t>(pointer) | StringTag;
    return *this;
  }
    
  RDValue& operator=(const std::vector<double> & v) {
    cleanup_rdvalue(*this);
    std::vector<double> * pointer = new std::vector<double>(v);
    otherBits = reinterpret_cast<uint64_t>(pointer) | VecDoubleTag;
    return *this;
  }

  RDValue& operator=(const std::vector<float> & v) {
    cleanup_rdvalue(*this);
    std::vector<float> * pointer = new std::vector<float>(v);
    otherBits = reinterpret_cast<uint64_t>(pointer) | VecFloatTag;
    return *this;
  }
  
  RDValue& operator=(const std::vector<int> & v) {
    cleanup_rdvalue(*this);
    std::vector<int> * pointer = new std::vector<int>(v);
    otherBits = reinterpret_cast<uint64_t>(pointer) | VecIntTag;
    return *this;
  }

  RDValue& operator=(const std::vector<unsigned int> & v) {
    cleanup_rdvalue(*this);
    std::vector<unsigned int> * pointer = new std::vector<unsigned int>(v);
    otherBits = reinterpret_cast<uint64_t>(pointer) | VecUnsignedIntTag;
    return *this;
  }
  
  RDValue& operator=(const std::vector<std::string> & v) {
    cleanup_rdvalue(*this);
    std::vector<std::string> * pointer = new std::vector<std::string>(v);
    otherBits = reinterpret_cast<uint64_t>(pointer) | VecStringTag;
    return *this;
  }

  RDValue& operator=(const boost::any & v) {
    cleanup_rdvalue(*this);
    boost::any * pointer = new boost::any(v);
    otherBits = reinterpret_cast<uint64_t>(pointer) | AnyTag;
    return *this;
  }
    
  template<class T>
  RDValue &operator=(const T&v) {
    cleanup_rdvalue(*this);
    boost::any *pointer = new boost::any(v);
    assert((reinterpret_cast<uint64_t>(pointer) & AnyTag) == 0);

    otherBits = reinterpret_cast<uint64_t>(pointer) | AnyTag;
    return *this;
  }
};

// Given two RDValue::Values - copy the appropriate structure
inline void copy_rdvalue(RDValue &dest,
                         const RDValue &src) {
  RDValue::cleanup_rdvalue(dest);
  switch(src.getTag()) {
    case RDValue::StringTag:
      dest = src.asString();
      break;
    case RDValue::VecDoubleTag:
      dest = src.asVecDouble();
      break;
    case RDValue::VecFloatTag:
      dest = src.asVecFloat();
      break;
    case RDValue::VecIntTag:
      dest = src.asVecInt();
      break;
    case RDValue::VecUnsignedIntTag:
      dest = src.asVecUnsignedInt();
      break;
    case RDValue::VecStringTag:
      dest = src.asVecString();
      break;
    case RDValue::AnyTag:
      dest = src.asAny();
      break;
    default:
      dest = src;
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Const access
template <class T>
const T &rdvalue_cast(const RDValue &v) {
  if (v.isAny())
    return boost::any_cast<const T &>(v.asAny());
  throw boost::bad_any_cast();
}


template<>
inline const boost::any &rdvalue_cast<boost::any>(const RDValue &v) {
  if (v.isAny()) return v.asAny();
  throw boost::bad_any_cast();
}

template<>
inline const double &rdvalue_cast<double>(const RDValue &v) {
  if (v.isDouble()) return v.asDouble();
  throw boost::bad_any_cast();
}

template<>
inline const float &rdvalue_cast<float>(const RDValue &v) {
  if (v.isFloat()) return v.asFloat();
  throw boost::bad_any_cast();
}

template<>
inline const int &rdvalue_cast<int>(const RDValue &v) {
  if (v.isInt32()) return v.asInt32();
  throw boost::bad_any_cast();
}
template<>
inline const unsigned int &rdvalue_cast<unsigned int>(const RDValue &v) {
  if (v.isUInt32()) return v.asUInt32();
  throw boost::bad_any_cast();
}

template<>
inline const bool &rdvalue_cast<bool>(const RDValue &v) {
  if (v.isBool()) return v.asBool();
  throw boost::bad_any_cast();
}

template<>
inline const std::string &rdvalue_cast<std::string>(const RDValue &v) {
  if (v.isString()) return v.asString();
  throw boost::bad_any_cast();
}

template<>
inline const std::vector<double> &rdvalue_cast<std::vector<double> >(const RDValue &v) {
  if (v.isVecDouble()) return v.asVecDouble();
  throw boost::bad_any_cast();
}

template<>
inline const std::vector<float> &rdvalue_cast<std::vector<float> >(const RDValue &v) {
  if (v.isVecFloat()) return v.asVecFloat();
  throw boost::bad_any_cast();
}

template<>
inline const std::vector<int> &rdvalue_cast<std::vector<int> >(const RDValue &v) {
  if (v.isVecInt()) return v.asVecInt();
  throw boost::bad_any_cast();
}

template<>
inline const std::vector<unsigned int> &rdvalue_cast<std::vector<unsigned int> >(
    const RDValue &v) {
  if (v.isVecUnsignedInt()) return v.asVecUnsignedInt();
  throw boost::bad_any_cast();
}

template<>
inline const std::vector<std::string> &rdvalue_cast<std::vector<std::string> >(
    const RDValue &v) {
  if (v.isVecString()) return v.asVecString();
  throw boost::bad_any_cast();
}

/////////////////////////////////////////////////////////////////////////////////////
// Direct access??
template <class T>
 T &rdvalue_cast( RDValue &v) {
  if (v.isAny())
    return boost::any_cast< T &>(v.asAny());
  throw boost::bad_any_cast();
}

template<>
inline boost::any &rdvalue_cast<boost::any>(RDValue &v) {
  if (v.isAny()) return v.asAny();
  throw boost::bad_any_cast();
}

template<>
inline  double &rdvalue_cast<double>( RDValue &v) {
  if (v.isDouble()) return v.asDouble();
  throw boost::bad_any_cast();
}

template<>
inline  float &rdvalue_cast<float>( RDValue &v) {
  if (v.isFloat()) return v.asFloat();
  throw boost::bad_any_cast();
}

template<>
inline  int &rdvalue_cast<int>( RDValue &v) {
  if (v.isInt32()) return v.asInt32();
  throw boost::bad_any_cast();
}
template<>
inline  unsigned int &rdvalue_cast<unsigned int>( RDValue &v) {
  if (v.isUInt32()) return v.asUInt32();
  throw boost::bad_any_cast();
}

template<>
inline  bool &rdvalue_cast<bool>( RDValue &v) {
  if (v.isBool()) return v.asBool();
  throw boost::bad_any_cast();
}

template<>
inline  std::string &rdvalue_cast<std::string>( RDValue &v) {
  if (v.isString()) return v.asString();
  throw boost::bad_any_cast();
}

template<>
inline  std::vector<double> &rdvalue_cast<std::vector<double> >( RDValue &v) {
  if (v.isVecDouble()) return v.asVecDouble();
  throw boost::bad_any_cast();
}

template<>
inline  std::vector<float> &rdvalue_cast<std::vector<float> >( RDValue &v) {
  if (v.isVecFloat()) return v.asVecFloat();
  throw boost::bad_any_cast();
}

template<>
inline  std::vector<int> &rdvalue_cast<std::vector<int> >( RDValue &v) {
  if (v.isVecInt()) return v.asVecInt();
  throw boost::bad_any_cast();
}

template<>
inline  std::vector<unsigned int> &rdvalue_cast<std::vector<unsigned int> >(
    RDValue &v) {
  if (v.isVecUnsignedInt()) return v.asVecUnsignedInt();
  throw boost::bad_any_cast();
}

template<>
inline  std::vector<std::string> &rdvalue_cast<std::vector<std::string> >(
     RDValue &v) {
  if (v.isVecString()) return v.asVecString();
  throw boost::bad_any_cast();
}

/////////////////////////////////////////////////////////////////////////////////////
// lexical casts...
template <class T>
std::string vectToString(const std::vector<T> &tv) {
  std::ostringstream sstr;
  sstr.imbue(std::locale("C"));
  sstr << std::setprecision(17);
  sstr << "[";
  std::copy(tv.begin(), tv.end(), std::ostream_iterator<T>(sstr, ","));
  sstr << "]";
  return sstr.str();
}

inline bool rdvalue_tostring(const RDValue &val, std::string &res) {
  Utils::LocaleSwitcher ls; // for lexical cast...
  switch (val.getTag() ) {
    case RDValue::Int32Tag:
      res = boost::lexical_cast<std::string>(val.asInt32());
      break;
    case RDValue::UInt32Tag: 
      res = boost::lexical_cast<std::string>(val.asUInt32());
      break;
    case RDValue::StringTag:
      res = val.asString();
      break;
    case RDValue::VecDoubleTag:
      res = vectToString<double>(val.asVecDouble());
      break;
    case RDValue::VecFloatTag:
      res = vectToString<float>(val.asVecFloat());
      break;
    case RDValue::VecIntTag:
      res = vectToString<int>(val.asVecInt());
      break;
    case RDValue::VecUnsignedIntTag:
      res = vectToString<unsigned int>(val.asVecUnsignedInt());
      break;
    case RDValue::VecStringTag:
      res = vectToString<std::string>(val.asVecString());
    case RDValue::AnyTag:
      try {
        res = boost::any_cast<std::string>(val.asAny());
      } catch (const boost::bad_any_cast &) {
        if (val.asAny().type() == typeid(long)) {
          res = boost::lexical_cast<std::string>(boost::any_cast<long>(val.asAny()));
        } else if (val.asAny().type() == typeid(unsigned long)) {
          res =
              boost::lexical_cast<std::string>(
                  boost::any_cast<unsigned long>(val.asAny()));
        } else {
          throw;
          return false;
        }
      }
      break;
    default:
      res = boost::lexical_cast<std::string>(val.asDouble());
  }
  return true;
}

// from_rdvalue -> converts string values to appropriate types
template <class T>
typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    const RDValue &arg) {
  T res;
  if (arg.getTag() == RDValue::StringTag) {
    Utils::LocaleSwitcher ls;
    try {
      res = rdvalue_cast<T>(arg);
    } catch (const boost::bad_any_cast &exc) {
      try {
        res = boost::lexical_cast<T>(rdvalue_cast<std::string>(arg));
      } catch (...) {
        throw exc;
      }
    }
  } else {
    res = rdvalue_cast<T>(arg);
  }
  return res;
}

template <class T>
typename boost::disable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    const RDValue &arg) {
  return rdvalue_cast<T>(arg);
}


} // namespace rdkit
#endif

