#include <iostream>
#include <string>
#include <unistd.h>

#include <boost/scope_exit.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/sync/interprocess_upgradable_mutex.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/interprocess/sync/sharable_lock.hpp>
#include <boost/interprocess/sync/upgradable_lock.hpp>

using namespace boost::interprocess;

namespace RDKit
{
typedef interprocess_upgradable_mutex upgradable_mutex_type;

struct RWSharedData {
  mutable upgradable_mutex_type mutex;
};

struct RWSharedLock {
  shared_memory_object shm;
  mapped_region *region;
  std::string name;
  RWSharedLock(const char * name) :
      shm(create_only, name, read_write),
      name(name)
  {
      shm.truncate(sizeof (RWSharedData));
      region = new mapped_region(shm, read_write);
      new (region->get_address()) RWSharedData;
  }

 ~RWSharedLock() {
   std::cerr << "Closing shared lock " << name << std::endl;
   delete region;
   shared_memory_object::remove(name.c_str());
  }
};
  
struct RWReadLock {
  shared_memory_object shm;
  mapped_region region;
  RWSharedData *shared_mutex;
  boost::interprocess::sharable_lock<upgradable_mutex_type> lock;
  RWReadLock(const char *name) :
      shm(open_only, name, read_write),
      region(shm, read_write),
      shared_mutex(static_cast<RWSharedData *>(region.get_address())),
      lock(shared_mutex->mutex)
  {
  }
};

struct RWWriteLock {
  shared_memory_object shm;
  mapped_region region;
  RWSharedData *shared_mutex;
  scoped_lock<upgradable_mutex_type> lock;
  RWWriteLock(const char *name) :
      shm(open_only, name, read_write),
      region(shm, read_write),
      shared_mutex(static_cast<RWSharedData *>(region.get_address())),
      lock(shared_mutex->mutex)
  {
  }
};
}
