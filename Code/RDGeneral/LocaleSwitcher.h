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

#ifndef _RD_LOCALESWITCHER_H
#define _RD_LOCALESWITCHER_H

#include <locale>
#include <stdlib.h>

  // Creates a c_locale to use in boost::format and
  //  various locale aware functions
  //    Use the static PosixLocale below
  // example:
  //   boost::format("ffff", PosixLocale.cpp_locale)
  //   ostream::imbue(PosixLocale.cpp_locale)
  //   atof_l(..., PosixLocale.c_locale)
  //   sprintf_l(..., PosixLocale.c_locale, ...)
  //     n.b. atof is deprecated, use strold?

#ifndef _WIN32
# ifdef __APPLE__
#  include <xlocale.h>
# else
#  include <locale.h>
# endif

namespace RDKit {
namespace Utils {

struct PosixLocaleInitializer {
  locale_t c_locale;
  std::locale cpp_locale;
  PosixLocaleInitializer() : cpp_locale("C") {
    c_locale = newlocale(LC_ALL_MASK, "C", (locale_t)0);
    
  }
  ~PosixLocaleInitializer() {
    freelocale(c_locale);
  }

  double atof(const char *str) const {
    return atof_l(str, c_locale);
  }
};

static const PosixLocaleInitializer PosixLocale;
}
}

#else // Windows

namespace RDKit {
namespace Utils {

struct PosixLocaleInitializer {
  _locale_t c_locale;
  std::locale cpp_locale;
  PosixLocaleInitializer() : cpp_locale("C") {
    c_locale = _create_locale(LC_ALL_MASK, "C");
    
  }
  ~PosixLocaleInitializer() {
    _freelocale(c_locale);
  }

  double atof(const char *str) const {
    return _atof_l(str, c_locale);
  }

};

static const PosixLocaleInitializer PosixLocale;
}
}

#endif

#endif
