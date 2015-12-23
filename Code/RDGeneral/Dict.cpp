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

bool toString(const RDAny &val, std::string &res) {
  Utils::LocaleSwitcher ls; // for lexical cast...
  switch (val.type) {
    case RDValueTypes::String:
      res = *val.m_value.s;
      break;
    case RDValueTypes::Int:
      res = boost::lexical_cast<std::string>(val.m_value.i);
      break;
    case RDValueTypes::Double:
      res = boost::lexical_cast<std::string>(val.m_value.d);
      break;
    case RDValueTypes::UnsignedInt:
      res = boost::lexical_cast<std::string>(val.m_value.u);
      break;
    case RDValueTypes::Bool:
      res = boost::lexical_cast<std::string>(val.m_value.b);
      break;
    case RDValueTypes::Float:
      res = boost::lexical_cast<std::string>(val.m_value.f);
      break;
    case RDValueTypes::VectDouble:
      res = vectToString<double>(*val.m_value.vd);
      break;
    case RDValueTypes::VectFloat:
      res = vectToString<float>(*val.m_value.vf);
      break;
    case RDValueTypes::VectInt:
      res = vectToString<int>(*val.m_value.vi);
      break;
    case RDValueTypes::VectUnsignedInt:
      res = vectToString<unsigned int>(*val.m_value.vu);
      break;
    case RDValueTypes::Any:
      const boost::any &any = val.asAny();
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
  toString(val, res);
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
  return toString(val, res);
}

/*
namespace {
template <class T>
typename boost::enable_if<boost::is_arithmetic<T>, T>::type _fromany(
    const RDAny &arg) {
  std::cerr << "enable if is aritmetic " << arg.type << std::endl;
  T res;
  if (arg.type == RDValueTypes::String) {
    Utils::LocaleSwitcher ls;
    try {
      res = rdany_cast<T>(arg);
    } catch (const boost::bad_any_cast &exc) {
      try {
        res = boost::lexical_cast<T>(rdany_cast<std::string>(arg));
      } catch (...) {
        throw exc;
      }
    }
  } else {
    res = rdany_cast<T>(arg);
  }
  return res;
}
template <class T>
typename boost::disable_if<boost::is_arithmetic<T>, T>::type _fromany(
    const RDAny &arg) {
  std::cerr << "disable if is aritmetic " << arg.type << std::endl;
  return rdany_cast<T>(arg);
}
}

template <typename T>
T Dict::fromany(const RDAny &arg) const {
  std::cerr << "calling from any" << std::endl;
  return _fromany<T>(arg);
};
template <typename T>
RDAny Dict::toany(T arg) const {
  return RDAny(arg);
};

#define ANY_FORCE(T)                                        \
  template T Dict::fromany<T>(const RDAny &arg) const; \
  template RDAny Dict::toany<T>(T arg) const;

ANY_FORCE(bool);
ANY_FORCE(boost::shared_array<double>);
ANY_FORCE(boost::shared_array<int>);
ANY_FORCE(double);
ANY_FORCE(int);
ANY_FORCE(std::list<int>);
ANY_FORCE(std::string);
ANY_FORCE(std::vector<boost::shared_array<double> >);
ANY_FORCE(std::vector<boost::shared_array<int> >);
ANY_FORCE(std::vector<double>);
ANY_FORCE(std::vector<int>);
ANY_FORCE(std::vector<std::list<int> >);
ANY_FORCE(std::vector<std::string>);
ANY_FORCE(std::vector<std::vector<double> >);
ANY_FORCE(std::vector<std::vector<int> >);
ANY_FORCE(std::vector<unsigned int>);
ANY_FORCE(std::vector<unsigned long long>);
ANY_FORCE(unsigned int);

template const std::string &Dict::fromany<const std::string &>(
    const RDAny &arg) const;

typedef boost::tuples::tuple<boost::uint32_t, boost::uint32_t, boost::uint32_t>
    uint32_t_tuple;
typedef boost::tuples::tuple<double, double, double> double_tuple;

template uint32_t_tuple Dict::fromany<uint32_t_tuple>(
    const RDAny &arg) const;
template RDAny Dict::toany<uint32_t_tuple>(uint32_t_tuple arg) const;
template double_tuple Dict::fromany<double_tuple>(const RDAny &arg) const;
template RDAny Dict::toany<double_tuple>(double_tuple arg) const;
*/
}
