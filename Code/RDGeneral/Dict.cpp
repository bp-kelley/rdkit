// $Id$
//
// Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Dict.h"

#include <boost/shared_array.hpp>
#include <boost/cstdint.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace RDKit {

namespace {
/*
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
*/
bool rdany_tostring(const RDAny &val, std::string &res) {
  Utils::LocaleSwitcher ls; // for lexical cast...
  switch (val.m_value.type) {
    case RDValue::String:
      res = *val.m_value.value.s;
      break;
    case RDValue::Int:
      res = boost::lexical_cast<std::string>(val.m_value.value.i);
      break;
    case RDValue::Double:
      res = boost::lexical_cast<std::string>(val.m_value.value.d);
      break;
    case RDValue::UnsignedInt:
      res = boost::lexical_cast<std::string>(val.m_value.value.u);
      break;
    case RDValue::Bool:
      res = boost::lexical_cast<std::string>(val.m_value.value.b);
      break;
    case RDValue::Float:
      res = boost::lexical_cast<std::string>(val.m_value.value.f);
      break;
    case RDValue::VectDouble:
      res = vectToString<double>(*val.m_value.value.vd);
      break;
    case RDValue::VectFloat:
      res = vectToString<float>(*val.m_value.value.vf);
      break;
    case RDValue::VectInt:
      res = vectToString<int>(*val.m_value.value.vi);
      break;
    case RDValue::VectUnsignedInt:
      res = vectToString<unsigned int>(*val.m_value.value.vu);
      break;
    case RDValue::VectString:
      res = vectToString<std::string>(*val.m_value.value.vs);
      break;
    case RDValue::Empty:
      res = "";
      break;
    case RDValue::Any:
      const boost::any &any = *val.m_value.value.a;
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
}

void Dict::getVal(const std::string &what, std::string &res) const {
  //
  //  We're going to try and be somewhat crafty about this getVal stuff to make
  //  these
  //  containers a little bit more generic.  The normal behavior here is that
  //  the
  //  value being queried must be directly castable to type T.  We'll robustify
  //  a
  //  little bit by trying that and, if the cast fails, attempting a couple of
  //  other casts, which will then be lexically cast to type T.
  //
  DataType::const_iterator pos = _data.find(what);
  if (pos == _data.end())
    throw KeyErrorException(what);    
  const RDAny &val = pos->second;
  rdany_tostring(val.m_value, res);
}

bool Dict::getValIfPresent(const std::string &what, std::string &res) const {
  //
  //  We're going to try and be somewhat crafty about this getVal stuff to make
  //  these
  //  containers a little bit more generic.  The normal behavior here is that
  //  the
  //  value being queried must be directly castable to type T.  We'll robustify
  //  a
  //  little bit by trying that and, if the cast fails, attempting a couple of
  //  other casts, which will then be lexically cast to type T.
  //
  DataType::const_iterator pos = _data.find(what);
  if (pos == _data.end()) return false;
  const RDAny  &val = pos->second;
  return rdany_tostring(val, res);
}

}
