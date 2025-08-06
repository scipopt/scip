/*
 Copyright (C) 2024 AMPL Optimization Inc.

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

 Author: Gleb Belov
 */
#ifndef NLWRITER2MISC_H
#define NLWRITER2MISC_H

#include <string>
#include <algorithm>
#include <cstdio>
#include <cstdarg>
#include <cassert>

#include "mp/nl-header.h"
#include "mp/nl-utils.h"

namespace mp {

/// Text formatter
class TextFormatter {
public:
  /// Construct
  TextFormatter(NLUtils& u,
                bool nlc,
                int outprec=0)
    : utils_(u),
      nl_comments(nlc),
      output_prec(outprec)
  { }
  /// Mode query
  NLHeader::Format Mode() const { return NLHeader::TEXT; }
  /// Text printf
  int apr(File&, const char*, ...);
  /// Text nput
  void nput(File&, double);
  /// Retrieve utils
  NLUtils& Utils() { return utils_; }


private:
  NLUtils& utils_;
  const bool nl_comments{false};
  const int output_prec{0};  // 0 means full precision
};


/// Binary formatter
class BinaryFormatter {
public:
  /// Construct
  BinaryFormatter(NLUtils& u, bool, int)
    : utils_(u) { }
  /// Mode query
  NLHeader::Format Mode() const { return NLHeader::BINARY; }
  /// Binary printf
  int apr(File&, const char*, ...);
  /// Binary nput
  void nput(File&, double);
  /// Retrieve utils
  NLUtils& Utils() { return utils_; }


private:
  NLUtils& utils_;
};


class StringFileWriter;


/// NLWriter2 parameters.
/// @param Formatter: low-level writer: binary or text.
/// @param Feeder: a class implementing the NLFeeder
/// concept that feeds model information on request.
/// @param Utils: writer utilities
template <class Formatter, class Feeder>
struct NLWriter2Params {
  using FormatterType=Formatter;
  using FeederType=Feeder;
};

}  // namespace mp

#endif // NLWRITER2MISC_H
