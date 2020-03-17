#include <boost/thread.hpp>
#include <boost/thread/shared_mutex.hpp>
typedef boost::shared_mutex READWRITE_lock;

# define READERLOCK(lock) boost::upgrade_lock< boost::shared_mutex > lock(lock)
# define WRITERLOCK(lock) \
    boost::upgrade_lock< boost::shared_mutex > lock(lock);\
    boost::upgrade_to_unique_lock< boost::shared_mutex > uniqueLock(lock)


int main() {
  READWRITE_lock lock;
  {
    READERLOCK(lock);
  }
  {
    WRITERLOCK(lock);
  }
  return 0;
}


  
