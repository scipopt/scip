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
#include "mp/nl-utils.h"

/////////////////////// TECHNICAL /////////////////////////
namespace mp {

void File::Printf(const char *format, ...) {
  va_list args;
  va_start (args, format);
  std::vfprintf (f_, format, args);
  va_end (args);
}

File NLUtils::
openf(const std::string& fname, int Close, const char *mode)
{
  File f;
  if (Close) {
    std::remove(fname.c_str());
    return f;
  }
  f.Open(fname.c_str(), mode);
  if (!f) {
    log_warning("can't open %s", fname.c_str());
    return f;
  }
  if (if_show_filenames())
    log_message("File %s\n", fname.c_str());
  return f;
}

void NLUtils::log_message(const char *format, ...) {
  va_list args;
  va_start (args, format);
  std::vprintf (format, args);
  va_end (args);
}

void NLUtils::log_warning(const char *format, ...) {
  std::fprintf(stderr, "WARNING: ");
  va_list args;
  va_start (args, format);
  std::vfprintf (stderr, format, args);
  va_end (args);
  std::fprintf(stderr, "\n");
}

}  // namespace mp
