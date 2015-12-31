#include "RDValue.h"
#include "Invariant.h"
#include <limits>
#include <vector>
#include <string>

using namespace RDKit;

template<class T>
void testLimits() {
  // check numeric limits
  {
    RDValue v(std::numeric_limits<T>::min());
    CHECK_INVARIANT(rdvalue_cast<T>(v) == std::numeric_limits<T>::min(), "bad min");
    CHECK_INVARIANT(rdvalue_cast<T>(RDValue(v)) == std::numeric_limits<T>::min(), "bad min");
    v = std::numeric_limits<T>::max();
    CHECK_INVARIANT(rdvalue_cast<T>(v) == std::numeric_limits<T>::max(), "bad max");
    CHECK_INVARIANT(rdvalue_cast<T>(RDValue(v)) == std::numeric_limits<T>::max(), "bad max");
    
  }
  {
    RDValue v(std::numeric_limits<T>::max());
    CHECK_INVARIANT(rdvalue_cast<T>(v) == std::numeric_limits<T>::max(), "bad max");
    RDValue vv(v);
    CHECK_INVARIANT(rdvalue_cast<T>(vv) == std::numeric_limits<T>::max(), "bad max");

    v = std::numeric_limits<T>::min();
    RDValue vvv(v);
    CHECK_INVARIANT(rdvalue_cast<T>(v) == std::numeric_limits<T>::min(), "bad min");
    CHECK_INVARIANT(rdvalue_cast<T>(vvv) == std::numeric_limits<T>::min(), "bad min");
  }
}

void testPOD() {
  testLimits<int>();
  testLimits<unsigned int>();
  testLimits<double>();
  testLimits<float>();
  testLimits<bool>();
}



template<class T>
void testVector() {
  T minv = std::numeric_limits<T>::min();
  T maxv = std::numeric_limits<T>::max();
  std::vector<T> data;
  data.push_back(minv);
  data.push_back(maxv);
  data.push_back(T());

  RDValue v(data);
  CHECK_INVARIANT(rdvalue_cast<std::vector<T> >(v) == data, "bad vec");
  RDValue vv; copy_rdvalue(vv,v);
  CHECK_INVARIANT(rdvalue_cast<std::vector<T> >(vv) == data, "bad copy constructor");
  RDValue::cleanup_rdvalue(v); // desctructor...
  RDValue::cleanup_rdvalue(vv);
}

void testPODVectors() {
  testVector<int>();
  testVector<unsigned int>();
  testVector<double>();
  testVector<float>();
  testVector<long double>(); // stored in anys
}

void testStringVect() {
  std::vector<std::string> vecs;
  vecs.push_back("my");
  vecs.push_back("dog");
  vecs.push_back("has");
  vecs.push_back("fleas");
  RDValue v(vecs);
  CHECK_INVARIANT(rdvalue_cast<std::vector<std::string> >(v) == vecs, "bad vect");
  RDValue vc; copy_rdvalue(vc, v);
  CHECK_INVARIANT(rdvalue_cast<std::vector<std::string> >(vc) == vecs, "bad vect");
  RDValue vv = vecs;
  RDValue vvc; copy_rdvalue(vvc, vv);
  CHECK_INVARIANT(rdvalue_cast<std::vector<std::string> >(vv) == vecs, "bad vect");          
  CHECK_INVARIANT(rdvalue_cast<std::vector<std::string> >(vvc) == vecs, "bad vect");

  RDValue::cleanup_rdvalue(v); // desctructor...
  RDValue::cleanup_rdvalue(vc); // desctructor...
  RDValue::cleanup_rdvalue(vv); // desctructor...
  RDValue::cleanup_rdvalue(vvc); // desctructor...
}

int main() {
  testPOD();
  testPODVectors();
  testStringVect();
}
