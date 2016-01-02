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

namespace RDTypeTag {
  const short EmptyTag           = 0;
  const short IntTag             = 1;
  const short DoubleTag          = 2;
  const short StringTag          = 3;
  const short FloatTag           = 4;
  const short BoolTag            = 5;
  const short UnsignedIntTag     = 6;
  const short AnyTag             = 7;
  const short VectDoubleTag      = 8;
  const short VectFloatTag       = 9;
  const short VectIntTag         = 10;
  const short VectUnsignedIntTag = 11;
  const short VectStringTag      = 12;
}

struct RDValue {
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
  short type;
  short reserved_tag; // 16 bit alignment

 inline RDValue(): value(0.0), type(RDTypeTag::EmptyTag) {}
  // Pod Style (Direct storage)
 inline RDValue(double v)   : value(v), type(RDTypeTag::DoubleTag) {}
 inline RDValue(float v)    : value(v), type(RDTypeTag::FloatTag) {}
 inline RDValue(int v)      : value(v), type(RDTypeTag::IntTag) {}
 inline RDValue(unsigned v) : value(v), type(RDTypeTag::UnsignedIntTag) {}
 inline RDValue(bool v)     : value(v), type(RDTypeTag::BoolTag) {}

 inline RDValue(boost::any *v)  : value(v),type(RDTypeTag::AnyTag) {}

  // Copies passed in pointers
 inline RDValue(const boost::any &v)  : value(new boost::any(v)),type(RDTypeTag::AnyTag) {}
 inline RDValue(const std::string &v) : value(new std::string(v)),type(RDTypeTag::StringTag){};
 template <class T>
 inline RDValue(const T &v) : value(new boost::any(v)),type(RDTypeTag::AnyTag) {}

 inline RDValue(const std::vector<double> &v) : value(new std::vector<double>(v)),
    type(RDTypeTag::VectDoubleTag) {}
 inline RDValue(const std::vector<float> &v)  : value(new std::vector<float>(v)),
    type(RDTypeTag::VectFloatTag) {}
 inline RDValue(const std::vector<int> &v)    : value(new std::vector<int>(v)),
    type(RDTypeTag::VectIntTag) {}
 inline RDValue(const std::vector<unsigned int> &v) :
  value(new std::vector<unsigned int>(v)),
    type(RDTypeTag::VectUnsignedIntTag) {}
 inline RDValue(const std::vector<std::string> &v) :
  value(new std::vector<std::string>(v)),
    type(RDTypeTag::VectStringTag) {}

  short getTag() const { return type; }
  
  void destroy() {
    switch (type) {
      case RDTypeTag::StringTag:
        delete value.s;
        break;
      case RDTypeTag::AnyTag:
        delete value.a;
        break;
      case RDTypeTag::VectDoubleTag:
        delete value.vd;
        break;
      case RDTypeTag::VectFloatTag:
        delete value.vf;
        break;
      case RDTypeTag::VectIntTag:
        delete value.vi;
        break;
      case RDTypeTag::VectUnsignedIntTag:
        delete value.vu;
        break;
      case RDTypeTag::VectStringTag:
        delete value.vs;
        break;
      default:
        break;
    }
    type = RDTypeTag::EmptyTag;
  }

  static // Given a type and an RDAnyValue - delete the appropriate structure
  inline void cleanup_rdvalue(RDValue &rdvalue) {
    rdvalue.destroy();
  }

  RDValue& operator=(double v) {
    cleanup_rdvalue(*this);
    value = v;
    type = RDTypeTag::DoubleTag;
    return *this;
  }
  RDValue& operator=(float v) {
    cleanup_rdvalue(*this);
    value = v;
    type = RDTypeTag::FloatTag;
    return *this;
  }
  RDValue& operator=(int  v) {
    cleanup_rdvalue(*this);
    value = v;
    type = RDTypeTag::IntTag;
    return *this;
  }
  RDValue& operator=(unsigned int  v) {
    cleanup_rdvalue(*this);
    value = v;
    type = RDTypeTag::UnsignedIntTag;
    return *this;
  }
  RDValue& operator=(bool v) {
    cleanup_rdvalue(*this);
    value = v;
    type = RDTypeTag::BoolTag;
    return *this;
  }
  RDValue& operator=(const std::string & v) {
    cleanup_rdvalue(*this);
    value = new std::string(v);
    type = RDTypeTag::StringTag;
    return *this;
  }
  RDValue& operator=(const std::vector<double> &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<double>(v);
    type = RDTypeTag::VectDoubleTag;
    return *this;
  }
  RDValue& operator=(const std::vector<float> &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<float>(v);
    type = RDTypeTag::VectFloatTag;
    return *this;
  }
  RDValue& operator=(const std::vector<int>  &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<int>(v);
    type = RDTypeTag::VectIntTag;
    return *this;
  }
  RDValue& operator=(const std::vector<unsigned int>  &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<unsigned int>(v);
    type = RDTypeTag::VectUnsignedIntTag;
    return *this;
  }
  RDValue& operator=(const std::vector<std::string>  &v) {
    cleanup_rdvalue(*this);
    value = new std::vector<std::string>(v);
    type = RDTypeTag::VectStringTag;
    return *this;
  }
  template<class T>
  RDValue& operator=(const T&v) {
    cleanup_rdvalue(*this);
    value = new boost::any(v);
    type = RDTypeTag::AnyTag;
    return *this;
  }
    
};

// Given two RDValue::Values - copy the appropriate structure
inline void copy_rdvalue(RDValue &dest,
                         const RDValue &src) {
  dest.destroy();
  dest.type = src.type;
  switch (src.type) {
    case RDTypeTag::StringTag:
      dest.value.s = new std::string(*src.value.s);
      break;
    case RDTypeTag::AnyTag:
      dest.value.a = new boost::any(*src.value.a);
      break;
    case RDTypeTag::VectDoubleTag:
      dest.value.vd = new std::vector<double>(*src.value.vd);
      break;
    case RDTypeTag::VectFloatTag:
      dest.value.vf = new std::vector<float>(*src.value.vf);
      break;
    case RDTypeTag::VectIntTag:
      dest.value.vi = new std::vector<int>(*src.value.vi);
      break;
    case RDTypeTag::VectUnsignedIntTag:
      dest.value.vu = new std::vector<unsigned int>(*src.value.vu);
      break;
    case RDTypeTag::VectStringTag:
      dest.value.vs = new std::vector<std::string>(*src.value.vs);
      break;
    default:
      dest = src;
  }
}


// Const access
template <class T>
const T &rdvalue_cast(const RDValue &v) {
  if (v.type == RDTypeTag::AnyTag)
    return boost::any_cast<const T &>(*v.value.a);
  throw boost::bad_any_cast();
}

template<>
inline const double &rdvalue_cast<double>(const RDValue &v) {
  if (v.type == RDTypeTag::DoubleTag) return v.value.d;
  throw boost::bad_any_cast();
}
template<>
inline const float &rdvalue_cast<float>(const RDValue &v) {
  if (v.type == RDTypeTag::FloatTag) return v.value.f;
  throw boost::bad_any_cast();
}
template<>
inline const int &rdvalue_cast<int>(const RDValue &v) {
  if (v.type == RDTypeTag::IntTag) return v.value.i;
  throw boost::bad_any_cast();
}
template<>
inline const unsigned int &rdvalue_cast<unsigned int>(const RDValue &v) {
  if (v.type == RDTypeTag::UnsignedIntTag) return v.value.u;
  throw boost::bad_any_cast();
}
template<>
inline const bool &rdvalue_cast<bool>(const RDValue &v) {
  if (v.type == RDTypeTag::BoolTag) return v.value.b;
  throw boost::bad_any_cast();
}
template<>
inline const std::string &rdvalue_cast<std::string>(const RDValue &v) {
  if (v.type == RDTypeTag::StringTag) return *v.value.s;
  throw boost::bad_any_cast();
}
template<>
inline const boost::any &rdvalue_cast<boost::any>(const RDValue &v) {
  if (v.type == RDTypeTag::AnyTag) return *v.value.a;
  throw boost::bad_any_cast();
}
template<>
inline const std::vector<double> &rdvalue_cast<std::vector<double> >(const RDValue &v) {
  if (v.type == RDTypeTag::VectDoubleTag) return *v.value.vd;
  throw boost::bad_any_cast();
}
template<>  
inline const std::vector<float> &rdvalue_cast<std::vector<float> >(const RDValue &v) {
  if (v.type == RDTypeTag::VectFloatTag) return *v.value.vf;
  throw boost::bad_any_cast();
}
template<>
inline const std::vector<int> &rdvalue_cast<std::vector<int> >(const RDValue &v) {
  if (v.type == RDTypeTag::VectIntTag) {
    return *v.value.vi;
  }
  throw boost::bad_any_cast();
}
template<>
inline const std::vector<unsigned int> &rdvalue_cast<std::vector<unsigned int> >(const RDValue &v) {
  if (v.type == RDTypeTag::VectUnsignedIntTag) return *v.value.vu;
  throw boost::bad_any_cast();
}
template<>  
inline const std::vector<std::string> &rdvalue_cast<std::vector<std::string> >(const RDValue &v) {
  if (v.type == RDTypeTag::VectStringTag) return *v.value.vs;
  throw boost::bad_any_cast();
}
    
// Direct access (used in rdvalue_cast)
template <class T>
T &rdvalue_cast(RDValue &v) {
  if (v.type == RDTypeTag::AnyTag)
    return boost::any_cast<T &>(*v.value.a);
  throw boost::bad_any_cast();
}

template<>
inline double &rdvalue_cast<double>(RDValue &v) {
  if (v.type == RDTypeTag::DoubleTag) return v.value.d;
  throw boost::bad_any_cast();
}
template<>
inline float &rdvalue_cast<float>(RDValue &v) {
  if (v.type == RDTypeTag::FloatTag) return v.value.f;
  throw boost::bad_any_cast();
}
template<>
inline int &rdvalue_cast<int>(RDValue &v) {
  if (v.type == RDTypeTag::IntTag) return v.value.i;
  throw boost::bad_any_cast();
}
template<>
inline unsigned int &rdvalue_cast<unsigned int>(RDValue &v) {
  if (v.type == RDTypeTag::UnsignedIntTag) return v.value.u;
  throw boost::bad_any_cast();
}
template<>
inline bool &rdvalue_cast<bool>(RDValue &v) {
  if (v.type == RDTypeTag::BoolTag) return v.value.b;
  throw boost::bad_any_cast();
}
template<>
inline boost::any &rdvalue_cast<boost::any>(RDValue &v) {
  if (v.type == RDTypeTag::AnyTag) return *v.value.a;
  throw boost::bad_any_cast();
}
template<>
inline std::vector<double> &rdvalue_cast<std::vector<double> >(RDValue &v) {
  if (v.type == RDTypeTag::VectDoubleTag) return *v.value.vd;
  throw boost::bad_any_cast();
}
template<>
inline std::vector<float> &rdvalue_cast<std::vector<float> >(RDValue &v) {
  if (v.type == RDTypeTag::VectFloatTag) return *v.value.vf;
  throw boost::bad_any_cast();
}
template<>
inline std::vector<int> &rdvalue_cast<std::vector<int> >(RDValue &v) {
  if (v.type == RDTypeTag::VectIntTag) {
    return *v.value.vi;
  }
  throw boost::bad_any_cast();
}
template<>
inline std::vector<unsigned int> &rdvalue_cast<std::vector<unsigned int> >(RDValue &v) {
  if (v.type == RDTypeTag::VectUnsignedIntTag) return *v.value.vu;
  throw boost::bad_any_cast();
}
template<>
inline std::vector<std::string> &rdvalue_cast<std::vector<std::string> >(RDValue &v) {
  if (v.type == RDTypeTag::VectStringTag) return *v.value.vs;
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
  if (v.type == RDTypeTag::StringTag) return *v.value.s;
  throw boost::bad_any_cast();
}



inline bool rdvalue_tostring(const RDValue &val, std::string &res) {
  Utils::LocaleSwitcher ls; // for lexical cast...
  switch (val.type) {
    case RDTypeTag::StringTag:
      res = *val.value.s;
      break;
    case RDTypeTag::IntTag:
      res = boost::lexical_cast<std::string>(val.value.i);
      break;
    case RDTypeTag::DoubleTag:
      res = boost::lexical_cast<std::string>(val.value.d);
      break;
    case RDTypeTag::UnsignedIntTag:
      res = boost::lexical_cast<std::string>(val.value.u);
      break;
    case RDTypeTag::BoolTag:
      res = boost::lexical_cast<std::string>(val.value.b);
      break;
    case RDTypeTag::FloatTag:
      res = boost::lexical_cast<std::string>(val.value.f);
      break;
    case RDTypeTag::VectDoubleTag:
      res = vectToString<double>(*val.value.vd);
      break;
    case RDTypeTag::VectFloatTag:
      res = vectToString<float>(*val.value.vf);
      break;
    case RDTypeTag::VectIntTag:
      res = vectToString<int>(*val.value.vi);
      break;
    case RDTypeTag::VectUnsignedIntTag:
      res = vectToString<unsigned int>(*val.value.vu);
      break;
    case RDTypeTag::VectStringTag:
      res = vectToString<std::string>(*val.value.vs);
      break;
    case RDTypeTag::EmptyTag:
      res = "";
      break;
    case RDTypeTag::AnyTag:
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
  if (arg.type == RDTypeTag::StringTag) {
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

