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
#ifndef RDKIT_RDANY_H
#define RDKIT_RDANY_H
#include <boost/any.hpp>
#include <string>
#include <vector>

namespace RDKit {

// RDValue does not dynamically create POD types (kind of like
//  cdiggins::any)  However, it doesn't use RTTI type info
//  directly, it uses a companion short valued type
//  to determine what to do.
// For unregistered types, it falls back to boost::any.
//  The Size of an RDAny is (sizeof(double) + sizeof(short) == 10 bytes)
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
namespace RDValueTypes {
  const short Empty       = 0;
  // POD Style
  const short Bool        = 1;
  const short Double      = 2;
  const short Float       = 3;
  const short Int         = 4;
  const short UnsignedInt = 5;

  // Memory Managed
  const short String = 100;
  const short Any    = 101;

  const short VectDouble      = 1001;
  const short VectFloat       = 1002;
  const short VectInt         = 1003;
  const short VectUnsignedInt = 1004;
  const short VectString      = 1005;
}

// RDValue does not manange memory of non-pod data
//  this must be done externally (string, Any, vector...)
union RDValue {
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

  RDValue() {}
  // Pod Style (Direct storage)
  RDValue(double v)   : d(v) {}
  RDValue(float v)    : d(v) {}
  RDValue(int v)      : i(v) {}
  RDValue(unsigned v) : u(v) {}
  RDValue(bool v)     : b(v) {}

  RDValue(boost::any *v)  : a(v) {}

  // Stores passed in pointers -- does not destruct, use RDAny
  /* -- interface too dangerous ?
  RDValue(std::string *v) : s(v) {};

  RDValue(std::vector<double> *v)       : vd(v) {}
  RDValue(std::vector<float> *v)        : vf(v) {}
  RDValue(std::vector<int> *v)          : vi(v) {}
  RDValue(std::vector<unsigned int> *v) : vu(v) {}
  RDValue(std::vector<std::string> *v)  : vs(v) {}
  */
  
  // Copies passed in pointers
  RDValue(const boost::any &v)  : a(new boost::any(v)) {}
  RDValue(const std::string &v) : s(new std::string(v)){};

  RDValue(const std::vector<double> &v) : vd(new std::vector<double>(v)) {}
  RDValue(const std::vector<float> &v)  : vf(new std::vector<float>(v)) {}
  RDValue(const std::vector<int> &v)    : vi(new std::vector<int>(v)) {}
  RDValue(const std::vector<unsigned int> &v)
      : vu(new std::vector<unsigned int>(v)) {}
  RDValue(const std::vector<std::string> &v)
      : vs(new std::vector<std::string>(v)) {}
};

// Given a type and an RDAnyValue - delete the appropriate structure
void cleanup_rdany(const RDValue &any, short type) {
  switch (type) {
    case RDValueTypes::String:
      delete any.s;
      break;
    case RDValueTypes::Any:
      delete any.a;
      break;
    case RDValueTypes::VectDouble:
      delete any.vd;
      break;
    case RDValueTypes::VectFloat:
      delete any.vf;
      break;
    case RDValueTypes::VectInt:
      delete any.vi;
      break;
    case RDValueTypes::VectUnsignedInt:
      delete any.vu;
      break;
    case RDValueTypes::VectString:
      delete any.vs;
      break;
  }
}

// Given two RDValues - copy the appropriate structure
void copy_rdany(RDValue &dest, short desttype,
                const RDValue &src, short srctype) {
  cleanup_rdany(dest, desttype);
  switch (srctype) {
    case RDValueTypes::String:
      dest.s = new std::string(*src.s);
      break;
    case RDValueTypes::Any:
      dest.a = new boost::any(*src.a);
      break;
    case RDValueTypes::VectDouble:
      dest.vd = new std::vector<double>(*src.vd);
      break;
    case RDValueTypes::VectFloat:
      dest.vf = new std::vector<float>(*src.vf);
      break;
    case RDValueTypes::VectInt:
      dest.vi = new std::vector<int>(*src.vi);
      break;
    case RDValueTypes::VectUnsignedInt:
      dest.vu = new std::vector<unsigned int>(*src.vu);
      break;
    case RDValueTypes::VectString:
      dest.vs = new std::vector<std::string>(*src.vs);
      break;
    default:
      dest = src;
  }
}

struct RDAny {
  RDValue m_value;
  short type;

  RDAny() : m_value(), type(RDValueTypes::Empty) {}
  RDAny(double d) : m_value(d), type(RDValueTypes::Double) {}
  RDAny(float d) : m_value(d), type(RDValueTypes::Float) {}
  RDAny(int d) : m_value(d), type(RDValueTypes::Int) {}
  RDAny(unsigned int d) : m_value(d), type(RDValueTypes::UnsignedInt) {}
  RDAny(bool d) : m_value(d), type(RDValueTypes::Bool) {}

  RDAny(boost::any *d) : m_value(d), type(RDValueTypes::Any) {}

  // Takes ownership
  /* -- interface too dangerous ?
  RDAny(std::string *d) : m_value(d), type(RDValueTypes::String) {}
  RDAny(std::vector<double> *d) : m_value(d), type(RDValueTypes::VectDouble) {}
  RDAny(std::vector<float> *d) : m_value(d), type(RDValueTypes::VectFloat) {}
  RDAny(std::vector<int> *d) : m_value(d), type(RDValueTypes::VectInt) {}
  RDAny(std::vector<unsigned int> *d)
      : m_value(d), type(RDValueTypes::VectUnsignedInt) {}
  RDAny(std::vector<std::string> *d)
      : m_value(d), type(RDValueTypes::VectString) {}
  */
  
  // Copies
  RDAny(const std::string &d) : m_value(d), type(RDValueTypes::String) {}
  RDAny(const boost::any &d) : m_value(d), type(RDValueTypes::Any) {}
  RDAny(const std::vector<double> &d)
      : m_value(d), type(RDValueTypes::VectDouble) {}
  RDAny(const std::vector<float> &d)
      : m_value(d), type(RDValueTypes::VectFloat) {}
  RDAny(const std::vector<int> &d) : m_value(d), type(RDValueTypes::VectInt) {}
  RDAny(const std::vector<unsigned int> &d)
      : m_value(d), type(RDValueTypes::VectUnsignedInt) {}
  RDAny(const std::vector<std::string> &d)
      : m_value(d), type(RDValueTypes::VectString) {}

  // Falls back to boost::any for non-registered types
  template <class T>
  RDAny(const T &d)
      : m_value(new boost::any(d)), type(RDValueTypes::Any) {}

  RDAny(const RDAny &rhs) {
    copy_rdany(m_value, type, rhs.m_value, rhs.type);
    type = rhs.type; // can't set before cleaning up
  }

  ~RDAny() { cleanup_rdany(m_value, type); }

  // For easy of use:
  //   RDAny v;
  //   v = 2.0;
  //   v = std::string("foo...");

  RDAny &operator=(double d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::Double;
    return *this;
  }

  RDAny &operator=(float d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::Float;
    return *this;
  }

  RDAny &operator=(int d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::Int;
    return *this;
  }

  RDAny &operator=(unsigned int d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::UnsignedInt;
    return *this;
  }

  RDAny &operator=(bool d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::Bool;
    return *this;
  }

  RDAny &operator=(const std::string &d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::String;
    return *this;
  }

  RDAny &operator=(const std::vector<double> &d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::VectDouble;
    return *this;
  }

  RDAny &operator=(const std::vector<float> &d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::VectFloat;
    return *this;
  }

  RDAny &operator=(const std::vector<int> &d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::VectInt;
    return *this;
  }

  RDAny &operator=(const std::vector<unsigned int> &d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::VectDouble;
    return *this;
  }

  RDAny &operator=(const std::vector<std::string> &d) {
    cleanup_rdany(m_value, type);
    m_value = d;
    type = RDValueTypes::VectString;
    return *this;
  }

  RDAny &operator=(const boost::any &d) {
    cleanup_rdany(m_value, type);
    m_value.a = new boost::any(d);
    type = RDValueTypes::Any;
    return *this;
  }

  // Const access
  const double &asDouble() const {
    if (type == RDValueTypes::Double) return m_value.d;
    throw boost::bad_any_cast();
  }

  const float &asFloat() const {
    if (type == RDValueTypes::Float) return m_value.f;
    throw boost::bad_any_cast();
  }

  const int &asInt() const {
    if (type == RDValueTypes::Int) return m_value.i;
    throw boost::bad_any_cast();
  }

  const unsigned int &asUnsignedInt() const {
    if (type == RDValueTypes::UnsignedInt) return m_value.u;
    throw boost::bad_any_cast();
  }

  const bool &asBool() const {
    if (type == RDValueTypes::Bool) return m_value.b;
    throw boost::bad_any_cast();
  }

  const std::string &asString() const {
    if (type == RDValueTypes::String) return *m_value.s;
    throw boost::bad_any_cast();
  }

  const boost::any &asAny() const {
    if (type == RDValueTypes::Any) return *m_value.a;
    throw boost::bad_any_cast();
  }

  const std::vector<double> &asVectDouble() const {
    if (type == RDValueTypes::VectDouble) return *m_value.vd;
    throw boost::bad_any_cast();
  }
  
  const std::vector<float> &asVectFloat() const {
    if (type == RDValueTypes::VectFloat) return *m_value.vf;
    throw boost::bad_any_cast();
  }
  
  const std::vector<int> &asVectInt() const {
    if (type == RDValueTypes::VectInt) return *m_value.vi;
    throw boost::bad_any_cast();
  }
  
  const std::vector<unsigned int> &asVectUnsignedInt() const {
    if (type == RDValueTypes::VectUnsignedInt) return *m_value.vu;
    throw boost::bad_any_cast();
  }
  
  const std::vector<std::string> &asVectString() const {
    if (type == RDValueTypes::VectString) return *m_value.vs;
    throw boost::bad_any_cast();
  }
    
  // Direct access (used in rdany_cast)
  double &asDouble() {
    if (type == RDValueTypes::Double) return m_value.d;
    throw boost::bad_any_cast();
  }

  float &asFloat() {
    if (type == RDValueTypes::Float) return m_value.f;
    throw boost::bad_any_cast();
  }

  int &asInt() {
    if (type == RDValueTypes::Int) return m_value.i;
    throw boost::bad_any_cast();
  }

  unsigned int &asUnsignedInt() {
    if (type == RDValueTypes::UnsignedInt) return m_value.u;
    throw boost::bad_any_cast();
  }

  bool &asBool() {
    if (type == RDValueTypes::Bool) return m_value.b;
    throw boost::bad_any_cast();
  }

  std::string &asString() {
    if (type == RDValueTypes::String) return *m_value.s;
    throw boost::bad_any_cast();
  }

  boost::any &asAny() {
    if (type == RDValueTypes::Any) return *m_value.a;
    throw boost::bad_any_cast();
  }

  std::vector<double> &asVectDouble() {
    if (type == RDValueTypes::VectDouble) return *m_value.vd;
    throw boost::bad_any_cast();
  }
  
  std::vector<float> &asVectFloat() {
    if (type == RDValueTypes::VectFloat) return *m_value.vf;
    throw boost::bad_any_cast();
  }
  
  std::vector<int> &asVectInt() {
    if (type == RDValueTypes::VectInt) return *m_value.vi;
    throw boost::bad_any_cast();
  }
  
  std::vector<unsigned int> &asVectUnsignedInt() {
    if (type == RDValueTypes::VectUnsignedInt) return *m_value.vu;
    throw boost::bad_any_cast();
  }
  
  std::vector<std::string> &asVectString() {
    if (type == RDValueTypes::VectString) return *m_value.vs;
    throw boost::bad_any_cast();
  }
};

// rdany_cast
//  falls back to boost::any for non registered types
////////////////////////////////////////////////////////////////
// Const Access
template <class T>
const T &rdany_cast(const RDAny &d) {
  return boost::any_cast<const T &>(d.asAny());
}

template <>
const double &rdany_cast<double>(const RDAny &d) {
  return d.asDouble();
}

template <>
const float &rdany_cast<float>(const RDAny &d) {
  return d.asFloat();
}

template <>
const int &rdany_cast<int>(const RDAny &d) {
  return d.asInt();
}

template <>
const unsigned int &rdany_cast<unsigned int>(const RDAny &d) {
  return d.asUnsignedInt();
}

template <>
const std::string &rdany_cast<std::string>(const RDAny &d) {
  return d.asString();
}

template <>
const std::vector<double> &rdany_cast<std::vector<double> >(const RDAny &d) {
  return d.asVectDouble();
}

template <>
const std::vector<float> &rdany_cast<std::vector<float> >(const RDAny &d) {
  return d.asVectFloat();
}

template <>
const std::vector<int> &rdany_cast<std::vector<int> >(const RDAny &d) {
  return d.asVectInt();
}

template <>
const std::vector<unsigned int> &rdany_cast<std::vector<unsigned int> >(
    const RDAny &d) {
  return d.asVectUnsignedInt();
}

template <>
const std::vector<std::string> &rdany_cast<std::vector<std::string> >(
    const RDAny &d) {
  return d.asVectString();
}

/////////////////////////////////////////////////////////////////////////
// Direct
template <class T>
T &rdany_cast(RDAny &d) {
  return boost::any_cast<T &>(d.asAny());
}

template <>
double &rdany_cast<double>(RDAny &d) {
  return d.asDouble();
}

template <>
float &rdany_cast<float>(RDAny &d) {
  return d.asFloat();
}

template <>
int &rdany_cast<int>(RDAny &d) {
  return d.asInt();
}

template <>
unsigned int &rdany_cast<unsigned int>(RDAny &d) {
  return d.asUnsignedInt();
}

template <>
std::string &rdany_cast<std::string>(RDAny &d) {
  return d.asString();
}

template <>
std::vector<double> &rdany_cast<std::vector<double> >(RDAny &d) {
  return d.asVectDouble();
}

template <>
std::vector<float> &rdany_cast<std::vector<float> >(RDAny &d) {
  return d.asVectFloat();
}

template <>
std::vector<int> &rdany_cast<std::vector<int> >(RDAny &d) {
  return d.asVectInt();
}

template <>
std::vector<unsigned int> &rdany_cast<std::vector<unsigned int> >(
    RDAny &d) {
  return d.asVectUnsignedInt();
}

template <>
std::vector<std::string> &rdany_cast<std::vector<std::string> >(
    RDAny &d) {
  return d.asVectString();
}
}
#endif
