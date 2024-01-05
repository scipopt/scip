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

#include <cstdlib>

#include "mp/format.h"

namespace mp {

#ifndef MP_ASSERT
  /// Debug assert
  #define MP_ASSERT(condition, message) \
    assert((condition) && message)
#endif

/// Assert even for Release
#define MP_ASSERT_ALWAYS(condition, message) \
  do { if (!(condition)) MP_RAISE(message); } while(0)

/// Use this to raise UnsupportedError
#define MP_UNSUPPORTED(name) \
  throw MakeUnsupportedError( name )

/// Raise error
#define MP_RAISE(msg) throw mp::Error(msg)

/// Raise with exit code
#define MP_RAISE_WITH_CODE(exit_code, msg) \
  throw mp::Error(msg, exit_code)

/// Raise infeasibility
#define MP_INFEAS(msg) \
  MP_RAISE_WITH_CODE(200, std::string("Model infeasible: ") + msg)

/// Silence unused parameter warnings
#define MP_UNUSED(x) (void)(x)
/// A general error.
class Error : public fmt::internal::RuntimeError {
  int exit_code_ = EXIT_FAILURE;
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
  /// Costruct from message formatting arguments?
  FMT_VARIADIC_(char, , Error, init, fmt::CStringRef)

  /// Construct from message and optional exit code
  Error(fmt::CStringRef msg, int c=-1) : exit_code_(c)
  { SetMessage(msg.c_str()); }

  /// The exit code
  int exit_code() const { return exit_code_; }
};

/// The operation is not supported by the object.
class UnsupportedError : public Error {
 public:
  FMT_VARIADIC_(char, , UnsupportedError, init, fmt::CStringRef)
};

/// Makes UnsupportedError with prefix "unsupported: ".
/// Use via macro MP_UNSUPPORTED
inline UnsupportedError MakeUnsupportedError(
    fmt::CStringRef format_str, fmt::ArgList args) {
  return UnsupportedError("unsupported: {}", fmt::format(format_str, args));
}
FMT_VARIADIC(UnsupportedError, MakeUnsupportedError, fmt::CStringRef)

/// An option error
class OptionError : public Error {
public:
  explicit OptionError(fmt::CStringRef message) : Error(message) {}
};

}  // namespace mp

#endif  // MP_ERROR_H_
