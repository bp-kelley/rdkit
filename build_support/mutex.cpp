#include <iostream>
#include <thread>
#include <mutex>

int value = 0;
std::mutex mutex;

// Reads the value and sets v to that value
void readValue(int& v) {
  mutex.lock();
  // Simulate some latency
  std::this_thread::sleep_for(std::chrono::seconds(1));
  v = value;
  mutex.unlock();
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
}
