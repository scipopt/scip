/*
 NL reader

 Copyright (C) 2013 - 2016 AMPL Optimization Inc

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

#ifndef MP_NL_
#define MP_NL_

#include "mp/error.h"   // MP_ASSERT
#include "mp/format.h"  // fmt::CStringRef

#include <cstring>  // std::strlen, std::size_t
#include <string>   // std::string

namespace mp {

/** A reference to a null-terminated string with size. */
class NLStringRef {
 private:
  const char *data_;
  std::size_t size_;

 public:
  /**
    \rst
    Constructs a string reference object from a C string computing
    the size with ``std::strlen``.
    \endrst
   */
  NLStringRef(const char *s) : data_(s), size_(std::strlen(s)) {}

  /** Constructs a string reference object from a C string and a size. */
  NLStringRef(const char *s, std::size_t size) : data_(s), size_(size) {
    MP_ASSERT(!s[size], "string not null-terminated");
  }

  /**
    \rst
    Constructs a string reference from an ``std::string`` object.
    \endrst
   */
  NLStringRef(const std::string &s) : data_(s.c_str()), size_(s.size()) {}

  /** Returns the pointer to a C string. */
  const char *c_str() const { return data_; }

  /** Returns the string size. */
  std::size_t size() const { return size_; }
};

// Flags for ReadNLFile and ReadNLString.
enum {
  /** Read variable bounds before anything else. */
  READ_BOUNDS_FIRST = 1
};

/**
  \rst
  Reads an optimization problem in the NL format from the string *str*
  and sends notifications of the problem components to the *handler* object.
  The handler class can be one of the following

  * derived from `mp::NLHandler` or `mp::NullNLHandler`,
  * `mp::Problem`,
  * provide an interface compatible with one of the above.

  Both *str* and *name* can be C strings or ``std::string`` objects.
  The *name* argument is used as the name of the input when reporting errors.
  *flags* can be either 0, which is the default, to read all constructs in
  the order they appear in the input, or `mp::READ_BOUNDS_FIRST` to read
  variable bounds after the NL header and before other constructs such as
  nonlinear expressions.
  \endrst
 */
template <typename Handler>
void ReadNLString(NLStringRef str, Handler &handler,
                  fmt::CStringRef name = "(input)", int flags = 0);

/**
  \rst
  Reads an optimization problem in the NL format from the file *filename*
  and sends notifications of the problem components to the *handler* object.
  The handler class can be one of the following

  * derived from `mp::NLHandler` or `mp::NullNLHandler`,
  * `mp::Problem`,
  * provide an interface compatible with one of the above.

  The *filename* argument can be a C string or an ``std::string`` object.
  *flags* can be either 0, which is the default, to read all constructs in
  the order they appear in the input, or `mp::READ_BOUNDS_FIRST` to read
  variable bounds after the NL header and before other constructs such as
  nonlinear expressions.

  **Example**::

    // Count the number of variable references in all nonlinear expressions.
    struct VarCounter : mp::NullNLHandler<int> {
      int num_vars;
      VarCounter() : num_vars(0) {}
      Reference OnVariableRef(int) {
        ++num_vars;
        return Reference();
      }
    };
    VarCounter counter;
    mp::ReadNLFile("test.nl", counter);
    fmt::print("The number of variable references is {}.", counter.num_vars);
  \endrst
 */
template <typename Handler>
void ReadNLFile(fmt::CStringRef filename, Handler &handler, int flags = 0);
}  // namespace mp

#endif  // MP_NL_
