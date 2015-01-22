// $Id$
//
//  Copyright (c) 2007-2014, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior written permission.
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

#ifndef __RD_PROPERTIES_H__
#define __RD_PROPERTIES_H__

#include "Dict.h"
#include <boost/foreach.hpp>

namespace RDKit
{
  class Properties
  {
    //! \brief The \c Properties class can be used to extend objects
    //!        with properties
    Dict dp_props;  // not a very good solution...
    STR_VECT *compLst;
    
  public:
  Properties() : dp_props(), compLst(0) {
    }

  Properties(const Properties &rhs) :
    dp_props(rhs.dp_props),
    compLst(rhs.compLst ? new STR_VECT(*rhs.compLst) : 0)
  {
  }

  ~Properties()
  {
    delete compLst;
  }

  void clear()
  {
    dp_props.reset();
  }

  void update(const Properties &other)
  {
    dp_props.update(other.dp_props);
    if(compLst)
    {
      if(other.compLst)
        *compLst = *other.compLst;
      else
        compLst->clear();
    }
    else
    {
      if(other.compLst)
        compLst = new STR_VECT(*other.compLst);
    }
  }
  
    // ------------------------------------
    //  Local Property Dict functionality
    // ------------------------------------
    //! returns a list with the names of our \c properties
    STR_VECT getPropList() const {
      return dp_props.keys();
    }

    //! sets a \c property value
    /*!
       \param key the name under which the \c property should be stored.
           If a \c property is already stored under this name, it will be
	   replaced.
       \param val the value to be stored
       \param computed (optional) allows the \c property to be flagged
           \c computed.
     */
    template <typename T>
    void setProp(const char *key, T val, bool computed=false) {
    
      //if(!dp_props) dp_props = new Dict();
      std::string what(key);
      setProp(what,val, computed);
    }

    //! \overload
    template <typename T>
    void setProp(const std::string &key, T val, bool computed=false)  {
      if (computed) {
        if (!compLst)
        {
          compLst = new STR_VECT;
          compLst->push_back(key);
        }
        else if (std::find(compLst->begin(), compLst->end(), key) == compLst->end()) {
            compLst->push_back(key);
        }
      }
      dp_props.setVal(key, val);
    }

    //! allows retrieval of a particular property value
    /*!

       \param key the name under which the \c property should be stored.
           If a \c property is already stored under this name, it will be
	   replaced.
       \param res a reference to the storage location for the value.

       <b>Notes:</b>
         - if no \c property with name \c key exists, a KeyErrorException will be thrown.
	 - the \c boost::lexical_cast machinery is used to attempt type conversions.
	   If this fails, a \c boost::bad_lexical_cast exception will be thrown.

    */
    template <typename T>
    void getProp(const char *key,T &res) const {
      dp_props.getVal(key,res);
    }
    //! \overload
    template <typename T>
    void getProp(const std::string &key,T &res) const {
      dp_props.getVal(key,res);
    }

    //! \overload
    template <typename T>
    T getProp(const char *key) const {
      return dp_props.getVal<T>(key);
    }
    //! \overload
    template <typename T>
    T getProp(const std::string &key) const {
      return dp_props.getVal<T>(key);
    }

    //! returns whether or not we have a \c property with name \c key
    //!  and assigns the value if we do
    template <typename T>
    bool getPropIfPresent(const char *key,T &res) const {
        return dp_props.getValIfPresent(key,res);
    }
    //! \overload
    template <typename T>
    bool getPropIfPresent(const std::string &key,T &res) const {
        return dp_props.getValIfPresent(key,res);
    }

    //! returns whether or not we have a \c property with name \c key
    bool hasProp(const char *key) const {
      return dp_props.hasVal(key);
    };
    //! \overload
    bool hasProp(const std::string &key) const {
      return dp_props.hasVal(key);
    };

    
    //! clears the value of a \c property
    /*!
       <b>Notes:</b>
         - if no \c property with name \c key exists, a KeyErrorException
	   will be thrown.
	 - if the \c property is marked as \c computed, it will also be removed
	   from our list of \c computedProperties
    */
    void clearProp(const char *key)  {
      std::string what(key);
      clearProp(what);
    };
    //! \overload
    void clearProp(const std::string &key)  {
      if(compLst)
      {
	STR_VECT_I svi = std::find(compLst->begin(), compLst->end(), key);
	if (svi != compLst->end()) {
	  compLst->erase(svi);
	}
      }
      dp_props.clearVal(key);
    };

    //! clears all of our \c computed \c properties
    void clearComputedProps()  {
      if (compLst)
      {
	BOOST_FOREACH(const std::string &sv,*compLst){
	  dp_props.clearVal(sv);
	}
	compLst->clear();
      }
    }
    
    bool isComputedProp(const std::string &key) const
    {
      if(compLst)
        return std::find(compLst->begin(), compLst->end(), key) != compLst->end();
      return false;
    }
    
    STR_VECT *getComputedProps() const {
      return compLst;
    }
  };
}

#endif
