#include <iostream>
#include <thread>
#include <shared_mutex>

int value = 0;
std::shared_mutex mutex;

// Reads the value and sets v to that value
void readValue(int& v) {
  mutex.lock_shared();
  // Simulate some latency
  std::this_thread::sleep_for(std::chrono::seconds(1));
  v = value;
  mutex.unlock_shared();
}

// Sets value to v
void setValue(int v) {
  mutex.lock();
  // Simulate some latency
  std::this_thread::sleep_for(std::chrono::seconds(1));
  value = v;
  mutex.unlock();
}

int main() {
  return 0;
}
