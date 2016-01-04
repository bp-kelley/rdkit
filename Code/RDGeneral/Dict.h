//
// Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Dict.h

  \brief Defines the Dict class

*/
#ifndef __RD_DICT_H__
#define __RD_DICT_H__

#include <map>
#include <string>
#include <vector>
#include "RDValue.h"
#include "Exceptions.h"

#include "BoostStartInclude.h"
#include <boost/lexical_cast.hpp>
#include "BoostEndInclude.h"

#include "tags.h"

namespace RDKit {
typedef std::vector<std::string> STR_VECT;

//! \brief The \c Dict class can be used to store objects of arbitrary
//!        type keyed by \c strings.
//!
//!  The actual storage is done using \c RDAny objects.
//!
class Dict {
  struct Pair {
    int     key;
    RDValue val;

   Pair() : key(), val() {}
   Pair(int k, RDValue_cast_t v) : key(k), val(v) {}
   Pair(const std::string &s, RDValue_cast_t v) : key(tagmap.get(s)), val(v) {}

  };
  
  typedef std::vector<Pair> DataType;
public:
  static RDTags tagmap;
  
  Dict() : _data(), _hasNonPodData(false) {  };

  Dict(const Dict &other) : _data(other._data) {
    _hasNonPodData = other._hasNonPodData;
    if (_hasNonPodData) {
      std::vector<Pair> data(other._data.size());
      _data.swap(data);
      for (size_t i=0; i< _data.size(); ++i) {
        _data[i].key = other._data[i].key;
        copy_rdvalue(_data[i].val, other._data[i].val);
      }
    }   
  }
  
  ~Dict() {
    reset(); // to clear pointers if necessary
  }
  
  Dict &operator=(const Dict &other) {
    _hasNonPodData = other._hasNonPodData;
    if (_hasNonPodData) {
      std::vector<Pair> data(other._data.size());
      _data.swap(data);
      for (size_t i=0; i< _data.size(); ++i) {
        _data[i].key = other._data[i].key;
        copy_rdvalue(_data[i].val, other._data[i].val);
      }
    } else {
      _data = other._data;      
    }    
    return *this;
  };

  //----------------------------------------------------------
  //! \brief Returns whether or not the dictionary contains a particular
  //!        key.
  bool hasVal(int tag) const {
    for(size_t i=0 ; i< _data.size(); ++i) {
      if (_data[i].key == tag) return true;
    }
    return false;
  };
  bool hasVal(const std::string & what) const {
    return hasVal(tagmap.get(what));
  }

  //----------------------------------------------------------
  //! Returns the set of keys in the dictionary
  /*!
     \return  a \c STR_VECT
  */
  STR_VECT keys() const {
    STR_VECT res;
    DataType::const_iterator item;
    for (item = _data.begin(); item != _data.end(); item++) {
      res.push_back(tagmap.get(item->key));
    }
    return res;
  }

  //----------------------------------------------------------
  //! \brief Gets the value associated with a particular key
  /*!
     \param what  the key to lookup
     \param res   a reference used to return the result

     <B>Notes:</b>
      - If \c res is a \c std::string, every effort will be made
        to convert the specified element to a string using the
        \c boost::lexical_cast machinery.
      - If the dictionary does not contain the key \c what,
        a KeyErrorException will be thrown.
  */
  template <typename T>
  void getVal(const std::string &what, T &res) const {
    res = getVal<T>(what);
  };
  
  template <typename T>
  void getVal(int tag, T &res) const {
    res = getVal<T>(tag);
  };

  //! \overload
  template <typename T>
  T getVal(int tag) const {
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].key == tag) {
        return from_rdvalue<T>(_data[i].val);
      }
    }
    throw KeyErrorException(common_properties::GetPropName(tag));
  }

  template <typename T>
  T getVal(const std::string &what) const {
    return getVal<T>(tagmap.get(what));
  }

  //! \overload
  void getVal(const std::string &what, std::string &res) const {
    return getVal(tagmap.get(what), res);
  }

  void getVal(int tag, std::string &res) const;

  //----------------------------------------------------------
  //! \brief Potentially gets the value associated with a particular key
  //!        returns true on success/false on failure.
  /*!
     \param what  the key to lookup
     \param res   a reference used to return the result

     <B>Notes:</b>
      - If \c res is a \c std::string, every effort will be made
        to convert the specified element to a string using the
        \c boost::lexical_cast machinery.
      - If the dictionary does not contain the key \c what,
        a KeyErrorException will be thrown.
  */

  template <typename T>
  bool getValIfPresent(const std::string &what, T &res) const {
    return getValIfPresent(tagmap.get(what), res);
  }
  
  template <typename T>
  bool getValIfPresent(int tag, T &res) const {
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].key == tag) {
        res = from_rdvalue<T>(_data[i].val);
        return true;
      }
    }
    return false;
  };


  //! \overload
  bool getValIfPresent(int tag, std::string &res) const;
  
  bool getValIfPresent(const std::string &what, std::string &res) const {
    return getValIfPresent(tagmap.get(what), res);
  }

  //----------------------------------------------------------
  //! \brief Sets the value associated with a key
  /*!

     \param what the key to set
     \param val  the value to store

     <b>Notes:</b>
        - If \c val is a <tt>const char *</tt>, it will be converted
           to a \c std::string for storage.
        - If the dictionary already contains the key \c what,
          the value will be replaced.
  */
  template <typename T>
  void setVal(const std::string &what, T &val) {
    return setVal(tagmap.get(what), val);
  }

  template <typename T>
  void setVal(int tag, T &val) {
    _hasNonPodData = true;
    
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].key == tag) {
        _data[i].val = val;
        return;
      }
    }
    _data.push_back(Pair(tag, val));
  };

  void setVal(const std::string &what, bool val) {
    setVal(tagmap.get(what), val);
  }
  
  void setVal(int tag, bool val) {
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].key == tag) {
        _data[i].val = val;
        return;
      }
    }
    _data.push_back(Pair(tag, val));
  }

  void setVal(const std::string &what, double val) {
    setVal(tagmap.get(what), val);
  }
  
  void setVal(int tag, double val) {
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].key == tag) {
        _data[i].val = val;
        return;
      }
    }
    _data.push_back(Pair(tag, val));
  }
  
  void setVal(const std::string &what, float val) {
    setVal(tagmap.get(what), val);
  }
  
  void setVal(int tag, float val) {        
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].key == tag) {
        _data[i].val = val;
        return;
      }
    }
    _data.push_back(Pair(tag, val));
  }
  
  void setVal(const std::string &what, int val) {
    setVal(tagmap.get(what), val);
  }
  
  void setVal(int tag, int val) {
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].key == tag) {
        _data[i].val = val;
        return;
      }
    }
    _data.push_back(Pair(tag, val));    
  }
  
  void setVal(const std::string &what, unsigned int val) {
    setVal(tagmap.get(what), val);
  }
  
  void setVal(int tag, unsigned int val) {
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].key == tag) {
        _data[i].val = val;
        return;
      }
    }
    _data.push_back(Pair(tag, val));
  }
  
  //! \overload
  void setVal(const std::string &what, const char *val) {
    std::string h(val);
    setVal(what, h);
  }

  void setVal(int tag, const char *val) {
    std::string h(val);
    setVal(tag, h);
  }

  //----------------------------------------------------------
  //! \brief Clears the value associated with a particular key,
  //!     removing the key from the dictionary.
  /*!

     \param what the key to clear

   <b>Notes:</b>
      - If the dictionary does not contain the key \c what,
        a KeyErrorException will be thrown.
  */
  void clearVal(const std::string &what) {
    clearVal(tagmap.get(what));
  }
  
  void clearVal(int tag) {
    for(DataType::iterator it = _data.begin(); it < _data.end() ; ++it) {
      if (it->key == tag) {
        _data.erase(it);
        return;
      }
    }
    throw KeyErrorException(common_properties::GetPropName(tag));
  };

  //----------------------------------------------------------
  //! \brief Clears all keys (and values) from the dictionary.
  //!
  void reset() {
    if (_hasNonPodData) {
      for (size_t i=0; i< _data.size(); ++i) {
        RDValue::cleanup_rdvalue(_data[i].val);
      }
    }
    DataType data;
    _data.swap(data);
  };

  //----------------------------------------------------------
  //! Converts a \c RDAny to type \c T
  /*!
     \param arg a \c RDAny reference

     \returns the converted object of type \c T
  */
  /*
  template <typename T>
      T fromany(const RDAny &arg) const {
    return from_rdany<T>(arg);
  }
  */
  //----------------------------------------------------------
  //! Converts an instance of type \c T to \c RDAny
  /*!
     \param arg the object to be converted

     \returns a \c RDAny instance
  */
  /*
  template <typename T>
      RDAny toany(T arg) const {
    return RDAny(arg);
  };
  */
 private:
  DataType _data;  //!< the actual dictionary
  bool     _hasNonPodData;
};
}
#endif
