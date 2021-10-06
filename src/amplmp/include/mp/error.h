/*
 Exception classes and assertions

 Copyright (C) 2013 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef MP_ERROR_H_
#define MP_ERROR_H_

#include "mp/format.h"

namespace mp {

#ifndef MP_ASSERT
# define MP_ASSERT(condition, message) assert((condition) && message)
#endif

#define MP_RAISE(msg) throw std::runtime_error(msg)
#define MP_WARNING(msg) Print(msg)
// A general error.
class Error : public fmt::internal::RuntimeError {
 protected:
  Error() {}

  void SetMessage(const std::string &message) {
    std::runtime_error &base = *this;
    base = std::runtime_error(message);
  }

  void init(fmt::CStringRef format_str, fmt::ArgList args) {
    SetMessage(fmt::format(format_str, args));
  }

 public:
  FMT_VARIADIC_(char, , Error, init, fmt::CStringRef)
  ~Error() throw() {}
};

// The operation is not supported by the object.
class UnsupportedError : public Error {
 public:
  FMT_VARIADIC_(char, , UnsupportedError, init, fmt::CStringRef)
};

// Makes UnsupportedError with prefix "unsupported: ".
inline UnsupportedError MakeUnsupportedError(
    fmt::CStringRef format_str, fmt::ArgList args) {
  return UnsupportedError("unsupported: {}", fmt::format(format_str, args));
}
FMT_VARIADIC(UnsupportedError, MakeUnsupportedError, fmt::CStringRef)
}  // namespace mp

#endif  // MP_ERROR_H_
