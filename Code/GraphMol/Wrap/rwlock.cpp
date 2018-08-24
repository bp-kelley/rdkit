//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior written
//       permission.
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
#ifndef _RD_WRAPPED_RWLOCK_H_
#define _RD_WRAPPED_RWLOCK_H_
#include <RDBoost/python.h>
#include <RDBoost/pyint_api.h>
#include <RDBoost/Wrap.h>
#include <RDGeneral/rwlock.h>

namespace RDKit {
RWSharedLock *SharedLock(const std::string &name) {
  shared_memory_object::remove(name.c_str());
  RWSharedLock *lock = new RWSharedLock(name.c_str());
  return lock;
}

RWReadLock *ScopedReadLock(const std::string &name) {
  return new RWReadLock(name.c_str());
}

RWWriteLock *ScopedWriteLock(const std::string &name) {
  return new RWWriteLock(name.c_str());
}
}

using namespace boost::python;

void wrap_rwlock() {
  class_<RDKit::RWSharedLock, boost::noncopyable>("RWSharedLock", init<const char*>());
  class_<RDKit::RWReadLock, boost::noncopyable>("RWReadLock", init<const char *>());
  class_<RDKit::RWWriteLock, boost::noncopyable>("RWWriteLock", init<const char *>());

  def("SharedLock", &RDKit::SharedLock,
      return_value_policy<manage_new_object>());
  def("ScopedWriteLock", &RDKit::ScopedWriteLock,
      return_value_policy<manage_new_object>());
  def("ScopedReadLock", &RDKit::ScopedReadLock,
      return_value_policy<manage_new_object>());  
}

#endif
