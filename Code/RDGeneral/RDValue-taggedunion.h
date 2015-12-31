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
#ifndef RDKIT_RDVALUE_TAGGED_UNION_H
#define RDKIT_RDVALUE_TAGGED_UNION_H

#include <stdint.h>
#include <cassert>
#include <boost/any.hpp>
#include "Invariant.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>
#include "LocaleSwitcher.h"


namespace RDKit {

// RDValue does not dynamically create POD types (kind of like
//  cdiggins::any)  However, it doesn't use RTTI type info
//  directly, it uses a companion short valued type
//  to determine what to do.
// For unregistered types, it falls back to boost::any.
//  The Size of an RDAny is (sizeof(double) + sizeof(short) == 10 bytes
//   (aligned to actually 16 so hard to pass as value type)
//
//   For the sake of compatibility, errors throw boost::bad_any_cast
//
// Examples:
//
//   RDAny v(2.);
//   v = 1;
//   std::vector<double> d;
//   v == d;
//   v.asDoubleVect().push_back(4.)
//   rdany_cast<std::vector<double>(v).push_back(4.)
//
//   Falls back to boost::any for non registered types
//   v = boost::shared_ptr<ROMol>(new ROMol(m));
//

// RDValue does not manange memory of non-pod data
//  this must be done externally (string, Any, vector...)
// Tagged union

struct RDValue {
  enum ValueType {EmptyTag = 0, BoolTag, DoubleTag, FloatTag, IntTag, UnsignedIntTag,
                  StringTag=100, AnyTag=101,
                  VectDoubleTag = 1001, VectFloatTag = 1002, VectIntTag = 1003,
                  VectUnsignedIntTag = 1004, VectStringTag = 1005};
  union Value {
    double d;
    float f;
    int i;
    unsigned u;
    bool b;
    std::string *s;
    boost::any *a;
    std::vector<double> *vd;
    std::vector<float> *vf;
    std::vector<int> *vi;
    std::vector<unsigned int> *vu;
    std::vector<std::string> *vs;

   inline Value() {}
   inline Value(double v) : d(v) {}
   inline Value(float v) : f(v) {}
   inline Value(int v) : i(v) {}
   inline Value(unsigned int v) : u(v) {}
   inline Value(bool  v) : b(v) {}
   inline Value(std::string *v) : s(v) {}
   inline Value(boost::any *v) : a(v) {}
   inline Value(std::vector<double> *v) : vd(v) {}
   inline Value(std::vector<float> *v) : vf(v) {}
   inline Value(std::vector<int> *v) : vi(v) {}
   inline Value(std::vector<unsigned int> *v) : vu(v) {}
   inline Value(std::vector<std::string> *v) : vs(v) {}
  } value;
  ValueType type;

 inline RDValue(): type(EmptyTag) {}
  // Pod Style (Direct storage)
 inline RDValue(double v)   : value(v), type(DoubleTag) {}
 inline RDValue(float v)    : value(v), type(FloatTag) {}
 inline RDValue(int v)      : value(v), type(IntTag) {}
 inline RDValue(unsigned v) : value(v), type(UnsignedIntTag) {}
 inline RDValue(bool v)     : value(v), type(BoolTag) {}

 inline RDValue(boost::any *v)  : value(v),type(AnyTag) {}

  // Copies passed in pointers
 inline RDValue(const boost::any &v)  : value(new boost::any(v)),type(AnyTag) {}
 inline RDValue(const std::string &v) : value(new std::string(v)),type(StringTag){};
 template <class T>
 inline RDValue(const T &v) : value(new boost::any(v)),type(AnyTag) {}

 inline RDValue(const std::vector<double> &v) : value(new std::vector<double>(v)),
    type(VectDoubleTag) {}
 inline RDValue(const std::vector<float> &v)  : value(new std::vector<float>(v)),
    type(VectFloatTag) {}
 inline RDValue(const std::vector<int> &v)    : value(new std::vector<int>(v)),
    type(VectIntTag) {}
 inline RDValue(const std::vector<unsigned int> &v) :
  value(new std::vector<unsigned int>(v)),
    type(VectUnsignedIntTag) {}
 inline RDValue(const std::vector<std::string> &v) :
  value(new std::vector<std::string>(v)),
    type(VectStringTag) {}

  ValueType getTag() const { return type; }
  static // Given a type and an RDAnyValue - delete the appropriate structure
  inline void cleanup_rdvalue(RDValue &rdvalue) {
    switch (rdvalue.type) {
      case RDValue::StringTag:
        delete rdvalue.value.s;
        break;
      case RDValue::AnyTag:
        delete rdvalue.value.a;
        break;
      case RDValue::VectDoubleTag:
        delete rdvalue.value.vd;
        break;
      case RDValue::VectFloatTag:
        delete rdvalue.value.vf;
        break;
      case RDValue::VectIntTag:
        delete rdvalue.value.vi;
        break;
      case RDValue::VectUnsignedIntTag:
        delete rdvalue.value.vu;
        break;
      case RDValue::VectStringTag:
        delete rdvalue.value.vs;
        break;
      default:
        break;
    }
    rdvalue.type = EmptyTag;
  }

  RDValue& operator=(double v) {
    cleanup_rdvalue(*this);
    value = v;
    type = DoubleTag;
    return *this;
  }
  RDValue& operator=(float v) {
    cleanup_rdvalue(*this);
    value = v;
    type = FloatTag;
    return *this;
  }
  RDValue& operator=(int  v) {
    cleanup_rdvalue(*this);
    value = v;
    type = IntTag;
    return *this;
  }
  RDValue& operator=(unsigned int  v) {
    cleanup_rdvalue(*this);
    value = v;
    type = UnsignedIntTag;
    return *this;
  }
  RDValue& operator=(bool v) {
    cleanup_rdvalue(*this);
    value = v;
    type = BoolTag;
    return *this;
  }
  RDValue& operator=(const std::string & v) {
    cleanup_rdvalue(*this);
    value = new std::string(v);
    type = StringTag;
    return *this;
  }
  RDValue& operator=(const std::vector<double> &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<double>(v);
    type = VectDoubleTag;
    return *this;
  }
  RDValue& operator=(const std::vector<float> &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<float>(v);
    type = VectFloatTag;
    return *this;
  }
  RDValue& operator=(const std::vector<int>  &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<int>(v);
    type = VectIntTag;
    return *this;
  }
  RDValue& operator=(const std::vector<unsigned int>  &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<unsigned int>(v);
    type = VectUnsignedIntTag;
    return *this;
  }
  RDValue& operator=(const std::vector<std::string>  &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<std::string>(v);
    type = VectStringTag;
    return *this;
  }
  template<class T>
  RDValue& operator=(const T&v) {
    cleanup_rdvalue(*this);
    value = new boost::any(v);
    type = AnyTag;
    return *this;
  }
    
};

// Given two RDValue::Values - copy the appropriate structure
inline void copy_rdvalue(RDValue &dest,
                         const RDValue &src) {
  RDValue::cleanup_rdvalue(dest);
  dest.type = src.type;
  switch (src.type) {
    case RDValue::StringTag:
      dest.value.s = new std::string(*src.value.s);
      break;
    case RDValue::AnyTag:
      dest.value.a = new boost::any(*src.value.a);
      break;
    case RDValue::VectDoubleTag:
      dest.value.vd = new std::vector<double>(*src.value.vd);
      break;
    case RDValue::VectFloatTag:
      dest.value.vf = new std::vector<float>(*src.value.vf);
      break;
    case RDValue::VectIntTag:
      dest.value.vi = new std::vector<int>(*src.value.vi);
      break;
    case RDValue::VectUnsignedIntTag:
      dest.value.vu = new std::vector<unsigned int>(*src.value.vu);
      break;
    case RDValue::VectStringTag:
      dest.value.vs = new std::vector<std::string>(*src.value.vs);
      break;
    default:
      dest = src;
  }
}


// Const access
template <class T>
const T &rdvalue_cast(const RDValue &v) {
  if (v.type == RDValue::AnyTag)
    return boost::any_cast<const T &>(*v.value.a);
  throw boost::bad_any_cast();
}

template<>
inline const double &rdvalue_cast<double>(const RDValue &v) {
  if (v.type == RDValue::DoubleTag) return v.value.d;
  throw boost::bad_any_cast();
}
template<>
inline const float &rdvalue_cast<float>(const RDValue &v) {
  if (v.type == RDValue::FloatTag) return v.value.f;
  throw boost::bad_any_cast();
}
template<>
inline const int &rdvalue_cast<int>(const RDValue &v) {
  if (v.type == RDValue::IntTag) return v.value.i;
  throw boost::bad_any_cast();
}
template<>
inline const unsigned int &rdvalue_cast<unsigned int>(const RDValue &v) {
  if (v.type == RDValue::UnsignedIntTag) return v.value.u;
  throw boost::bad_any_cast();
}
template<>
inline const bool &rdvalue_cast<bool>(const RDValue &v) {
  if (v.type == RDValue::BoolTag) return v.value.b;
  throw boost::bad_any_cast();
}
template<>
inline const std::string &rdvalue_cast<std::string>(const RDValue &v) {
  if (v.type == RDValue::StringTag) return *v.value.s;
  throw boost::bad_any_cast();
}
template<>
inline const boost::any &rdvalue_cast<boost::any>(const RDValue &v) {
  if (v.type == RDValue::AnyTag) return *v.value.a;
  throw boost::bad_any_cast();
}
template<>
inline const std::vector<double> &rdvalue_cast<std::vector<double> >(const RDValue &v) {
  if (v.type == RDValue::VectDoubleTag) return *v.value.vd;
  throw boost::bad_any_cast();
}
template<>  
inline const std::vector<float> &rdvalue_cast<std::vector<float> >(const RDValue &v) {
  if (v.type == RDValue::VectFloatTag) return *v.value.vf;
  throw boost::bad_any_cast();
}
template<>
inline const std::vector<int> &rdvalue_cast<std::vector<int> >(const RDValue &v) {
  if (v.type == RDValue::VectIntTag) {
    return *v.value.vi;
  }
  throw boost::bad_any_cast();
}
template<>
inline const std::vector<unsigned int> &rdvalue_cast<std::vector<unsigned int> >(const RDValue &v) {
  if (v.type == RDValue::VectUnsignedIntTag) return *v.value.vu;
  throw boost::bad_any_cast();
}
template<>  
inline const std::vector<std::string> &rdvalue_cast<std::vector<std::string> >(const RDValue &v) {
  if (v.type == RDValue::VectStringTag) return *v.value.vs;
  throw boost::bad_any_cast();
}
    
// Direct access (used in rdvalue_cast)
template <class T>
T &rdvalue_cast(RDValue &v) {
  if (v.type == RDValue::AnyTag)
    return boost::any_cast<T &>(*v.value.a);
  throw boost::bad_any_cast();
}

template<>
inline double &rdvalue_cast<double>(RDValue &v) {
  if (v.type == RDValue::DoubleTag) return v.value.d;
  throw boost::bad_any_cast();
}
template<>
inline float &rdvalue_cast<float>(RDValue &v) {
  if (v.type == RDValue::FloatTag) return v.value.f;
  throw boost::bad_any_cast();
}
template<>
inline int &rdvalue_cast<int>(RDValue &v) {
  if (v.type == RDValue::IntTag) return v.value.i;
  throw boost::bad_any_cast();
}
template<>
inline unsigned int &rdvalue_cast<unsigned int>(RDValue &v) {
  if (v.type == RDValue::UnsignedIntTag) return v.value.u;
  throw boost::bad_any_cast();
}
template<>
inline bool &rdvalue_cast<bool>(RDValue &v) {
  if (v.type == RDValue::BoolTag) return v.value.b;
  throw boost::bad_any_cast();
}
template<>
inline boost::any &rdvalue_cast<boost::any>(RDValue &v) {
  if (v.type == RDValue::AnyTag) return *v.value.a;
  throw boost::bad_any_cast();
}
template<>
inline std::vector<double> &rdvalue_cast<std::vector<double> >(RDValue &v) {
  if (v.type == RDValue::VectDoubleTag) return *v.value.vd;
  throw boost::bad_any_cast();
}
template<>
inline std::vector<float> &rdvalue_cast<std::vector<float> >(RDValue &v) {
  if (v.type == RDValue::VectFloatTag) return *v.value.vf;
  throw boost::bad_any_cast();
}
template<>
inline std::vector<int> &rdvalue_cast<std::vector<int> >(RDValue &v) {
  if (v.type == RDValue::VectIntTag) {
    return *v.value.vi;
  }
  throw boost::bad_any_cast();
}
template<>
inline std::vector<unsigned int> &rdvalue_cast<std::vector<unsigned int> >(RDValue &v) {
  if (v.type == RDValue::VectUnsignedIntTag) return *v.value.vu;
  throw boost::bad_any_cast();
}
template<>
inline std::vector<std::string> &rdvalue_cast<std::vector<std::string> >(RDValue &v) {
  if (v.type == RDValue::VectStringTag) return *v.value.vs;
  throw boost::bad_any_cast();
}

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

template<>
inline std::string &rdvalue_cast<std::string>(RDValue &v) {
  if (v.type == RDValue::StringTag) return *v.value.s;
  throw boost::bad_any_cast();
}


inline bool rdvalue_tostring(const RDValue &val, std::string &res) {
  Utils::LocaleSwitcher ls; // for lexical cast...
  switch (val.type) {
    case RDValue::StringTag:
      res = *val.value.s;
      break;
    case RDValue::IntTag:
      res = boost::lexical_cast<std::string>(val.value.i);
      break;
    case RDValue::DoubleTag:
      res = boost::lexical_cast<std::string>(val.value.d);
      break;
    case RDValue::UnsignedIntTag:
      res = boost::lexical_cast<std::string>(val.value.u);
      break;
    case RDValue::BoolTag:
      res = boost::lexical_cast<std::string>(val.value.b);
      break;
    case RDValue::FloatTag:
      res = boost::lexical_cast<std::string>(val.value.f);
      break;
    case RDValue::VectDoubleTag:
      res = vectToString<double>(*val.value.vd);
      break;
    case RDValue::VectFloatTag:
      res = vectToString<float>(*val.value.vf);
      break;
    case RDValue::VectIntTag:
      res = vectToString<int>(*val.value.vi);
      break;
    case RDValue::VectUnsignedIntTag:
      res = vectToString<unsigned int>(*val.value.vu);
      break;
    case RDValue::VectStringTag:
      res = vectToString<std::string>(*val.value.vs);
      break;
    case RDValue::EmptyTag:
      res = "";
      break;
    case RDValue::AnyTag:
      const boost::any &any = *val.value.a;
      try {
        res = boost::any_cast<std::string>(val);
      } catch (const boost::bad_any_cast &) {
        if (any.type() == typeid(long)) {
          res = boost::lexical_cast<std::string>(boost::any_cast<long>(val));
        } else if (any.type() == typeid(unsigned long)) {
          res =
              boost::lexical_cast<std::string>(boost::any_cast<unsigned long>(val));
        } else {
          throw;
          return false;
        }
      }
      break;
  }
  return true;
}

// from_rdvalue -> converts string values to appropriate types
template <class T>
typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    const RDValue &arg) {
  T res;
  if (arg.type == RDValue::StringTag) {
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

