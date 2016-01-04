#ifndef RDKIT_TAGS_H
#define RDKIT_TAGS_H

#include "BoostStartInclude.h"

#ifdef RDK_TEST_MULTITHREADED
#include <boost/thread/mutex.hpp>
#include <boost/thread/tss.hpp>
#endif

#include <boost/unordered_map.hpp>
#include "BoostEndInclude.h"
#include "types.h"

namespace RDKit
{

class RDTags
{
  typedef boost::unordered_map<std::string, int> MAP;
public:
  RDTags() {
    for(int i=0;i<=common_properties::MAX;++i) {
      m[common_properties::GetPropName(i)] = i;
      keys.push_back(common_properties::GetPropName(i));
    }
  };

  const std::string &get(int tag) const {
    return keys[tag]; // raises exception
  }
  
  int get(const std::string &k) const
  {
#ifdef RDK_TEST_MULTITHREADED
    static boost::thread_specific_ptr<MAP> instance;

    if ( !instance.get() ) {
      boost::unique_lock<boost::mutex> lock(_m);
      // enable TLS loaded with current map
      instance.reset(new MAP(m));
    }
    
    MAP &ptr = *instance;
    // check tls map
    {
      MAP::const_iterator it=ptr.find(k);
      if (it != ptr.end())
        return it->second;
    }
#endif
    
    {
      boost::unique_lock<boost::mutex> lock(_m);
      
      {
        MAP::const_iterator it=m.find(k);
        
        if (it != m.end()) {
#ifdef RDK_TEST_MULTITHREADED          
          // insert into local map        
          instance->insert(*it);
#endif
          return it->second;
        }
      }
      {
        int res = m.size();
        std::pair<std::string, int> p(k,res);
#ifdef RDK_TEST_MULTITHREADED        
        instance->insert(p);
#endif        
        m.insert(p);
        keys.push_back(k);
        return res;
      }
    }
  }
  
private:
    mutable boost::mutex _m;
public:
    mutable MAP m;
    mutable std::vector<std::string > keys;
};

}

#endif
