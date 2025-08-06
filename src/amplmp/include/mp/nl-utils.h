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
#ifndef NLUTILS_H
#define NLUTILS_H

#include <string>
#include <cstdio>
#include <cstdarg>
#include <cassert>


namespace mp {

/// A wrapper for file I/O.
class File {
public:
  /// Destructor
  ~File() { Close(); }
  /// Constructor.
  File() { }
  /// Constructor.
  File(const File& f) { assert(!f_ && !f.f_); }
  /// Constructor.
  File(File&& f) : f_(f.f_) { f.f_ = nullptr; }
  /// copy assign
  void operator=(const File& f) { assert(!f_ && !f.f_); }
  /// move assign
  void operator=(File&& f) noexcept { std::swap(f_, f.f_); }

  /// Check being open for now.
  operator bool() const { return f_; }
  /// Open file.
  bool Open(const char* nm, const char* mode) {
    f_ = std::fopen(nm, mode);
    return (*this);
  }
  /// Close file.
  void Close() {
    if (f_) {
      std::fclose(f_);
    }
    f_ = nullptr;
  }

  /// Printf into the file
  void Printf(const char* format, ...);

  /// Obtain handler
  FILE* GetHandle() const { return f_; }

private:
  FILE *f_{nullptr};
};


/// NL writer and SOL reader utilities.
/// It provides default facilities for logging
/// and error handling.
class NLUtils {
public:
  /// Destructor
  virtual ~NLUtils() { }
  /// File opener
  /// @param Close: remove file instead
  virtual File openf(const std::string& fname,
                     int Close, const char *mode);
  /// log message
  virtual void log_message(const char* format, ...);
  /// log warning (e.g., cannot write name file)
  virtual void log_warning(const char* format, ...);
  /// Whether to show file names
  /// in the default openf() via log_message().
  virtual bool if_show_filenames() { return show_fln_; }
  /// Set: whether to show file names
  virtual void set_show_filenames(bool v) { show_fln_ = v; }

  /// Override this to your error handler.
  /// Not using exceptions by default.
  /// Only called on a wrong output format string
  /// (internal error.)
  virtual void myexit(const std::string& msg) {
    // No help including <cstdio> after <stdio.h> on Mac
    using namespace std;
    fprintf(stderr, "%s\n", msg.c_str());
    exit(1);
  }


private:
  bool show_fln_{false};
};


}  // namespace mp

#endif // NLUTILS_H
