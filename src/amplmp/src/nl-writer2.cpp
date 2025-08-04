/*
 NL writer

 NL is a format for representing optimization problems such as linear,
 quadratic, nonlinear, complementarity and constraint programming problems
 in discrete or continuous variables. It is described in the technical report
 "Writing .nl Files" (http://www.cs.sandia.gov/~dmgay/nlwrite.pdf).

 This is a complete reusable C++ implementation of an NL writer.

 Usage:
   /// Write an NL file:
   WriteNLFile(filename, feeder);

   /// Write an NL string:
   WriteNLString(string, feeder);

 where feeder is an object that provides information on model
 components. See NLFeeder for an example of a feeder class.

 See also NLReader and NLHandler classes.

 Copyright (C) 2023 AMPL Optimization Inc.

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

#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>

#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wunused-function"
#endif

#include "mp/nl-writer2.h"
#include "mp/nl-writer2.hpp"

#ifndef ssize_t        // signed size_t (POSIX)
  #define ssize_t int
#endif

/// Check that the stub NLFeeder
/// compiles with NLWriter2
static void test_nlfeeder_with_nlwriter2() {
  using Expr = void*;
  class TestNLF : public mp::NLFeeder<TestNLF, Expr> {
  public:
  };
  TestNLF nlfeeder;
  mp::NLUtils utils;
  mp::WriteNLFile("test.nl", nlfeeder, utils);
}


namespace DAVID_GAY_GFMT {
char* gfmt(double x, int prec);
}

namespace NL_LIB_GFMT {
void gfmt(char*, size_t, double, int);
}


////////////////////////////////////////////////////////
////////////////////// NLW2 UTILS //////////////////////
////////////////////////////////////////////////////////
namespace mp {

extern "C"
NLHeader_C MakeNLHeader_C_Default(void) {
  NLHeader nlh;
  NLHeader_C nlh_c;
  nlh_c.pi = *(NLProblemInfo_C*)(&nlh);
  nlh_c.nli = *(NLInfo_C*)(&nlh);
  return nlh_c;
}


/// \brief For printing to .nl or auxiliary files.
int TextFormatter::apr(File& f, const char *fmt, ...)
{
  auto fd = f.GetHandle();
  char buf[32], *s;
  double x;
  int rc;
  ssize_t i, j;
  va_list ap;

  rc = 0;

  va_start(ap, fmt);
  if (*fmt != '%')
    rc++;
  for(;;) {
    for(;;) {
      switch(i = *fmt++) {
      case 0:	  goto done;
      case '%': break;
      case '\t':
        if (!nl_comments) {
          putc('\n', fd);
          va_end(ap);
          return rc;
        }
        /* no break */
      default:
        putc(i, fd);
        continue;
      }
      break;
    }
    rc++;
    switch(*fmt++) {
    case 'c':
      i = va_arg(ap, int);
      putc(i, fd);
      continue;
    case 'z':
      i = va_arg(ap, size_t);
      goto have_i;
    case 'd':
      i = va_arg(ap, int);
have_i:
      if (i < 0) {
        putc('-',fd);
        i = -i;
      }
      s = buf;
      do {
        j = i / 10;
        *s++ = i - 10*j + '0';
      }
      while((i = j));
      do {
        i = *--s;
        putc(i,fd);
      }	while (s > buf);
      continue;
    case '.':
      while (*fmt++ != 'g');
    case 'g':
      x = va_arg(ap, double);
#ifdef NL_LIB_USE_SPRINTF
      snprintf(s = buf,
               sizeof(buf), "%.*g", output_prec, x);
#elif NL_LIB_USE_OWN_GFMT
      NL_LIB_GFMT::gfmt(s = buf, sizeof(buf), x, output_prec);
#else
      s = DAVID_GAY_GFMT::gfmt(x, output_prec);
#endif
      goto have_s;
    case 's':
      s = va_arg(ap, char*);
have_s:
      while ((i = *s++))
        putc(i,fd);
      continue;
    default:
      Utils().myexit("aprintf bug: unexpected fmt: " +
                     std::string(fmt-1));
    }
  }
done:
  va_end(ap);
  return rc;
}

void TextFormatter::nput(File& nm, double r)
{
  apr(nm, "n%g\n", r);
  /* formerly fprintf(nm, "n%.16g\n", *r); */
}

void BinaryFormatter::nput(File& nm, double r)
{
  double x;
  long L;
  short sh;

  if ((x = r) <= 2147483647.
      && x >= -2147483648.
      && (L = (long)x, (double)x == L)) {
    sh = (short)L;
    if (sh == L)
      apr(nm, "s%h", sh);
    else
      apr(nm, "l%l", L);
  }
  else
    apr(nm, "n%g", x);
}

namespace {
typedef int Long;
typedef Long Int;
/// For printing to .nl/.mps. used in #outmps.
typedef union U { double x; short sh; Long L; Int i; char c; size_t z; Long LL[2];
                #ifndef NO_LONG_LONG
                  long long ll;
                #endif
                } U;
}

int BinaryFormatter::
apr(File& f, const char *fmt, ...)
{
  auto fd = f.GetHandle();
  U u;
  char *s;
  int i, rc;
  size_t len, x;
  va_list ap;

  rc = 0;
  va_start(ap, fmt);

  if ((i = *fmt) != '%') {
    fmt++;
#define DMG
#ifdef DMG
    if (i != 'i')
#endif
    {
      u.c = i;
      fwrite(&u.c, 1, 1, fd);
      rc++;
    }}

  for(;;) {
    while(*fmt == ' ')
      fmt++;
    if (*fmt++ != '%')
      break;
    switch(*fmt++) {
    case 'c':
      u.c = va_arg(ap, int);
      len = 1;
      break;
    case 'd':
      u.i = va_arg(ap, int);
      len = sizeof(Int);
      break;
    case '.':
      while(*fmt++ != 'g');
    case 'g':
      u.x = va_arg(ap, double);
      len = sizeof(double);
      break;
    case 'h':
      u.sh = va_arg(ap, int);
      len = sizeof(short);
      if (*fmt == 'd')
        fmt++;
      break;
    case 'l':
      u.L = (Long)va_arg(ap, long);
      len = sizeof(Long);
      if (*fmt == 'd')
        fmt++;
      break;
    case 's':
      s = va_arg(ap, char*);
      u.i = std::strlen(s);
      fwrite((char *)&u.i, sizeof(Int), 1, fd);
      fwrite(s, u.i, 1, fd);
      goto s_written;
    case 'z':
      x = va_arg(ap, size_t);
      u.L = (Long)x;
      len = sizeof(Long);
      break;
    default:
      Utils().myexit("bprintf bug: unexpected fmt: " +
                     std::string(fmt-1));
#define AVOID_BOGUS_WARNINGS
#ifdef AVOID_BOGUS_WARNINGS
      len = 0;
#endif
    }
    fwrite((char *)&u.L, len, 1, fd);
s_written:
    rc++;
    if (*fmt == ':') {
      /* special Hollerith */
      s = va_arg(ap, char*);
      fwrite(s, u.i, 1, fd);
      fmt += 3;
      rc++;
    }
  }
  va_end(ap);
  return rc;
}

} // namespace mp


extern "C" {
char *
dtoa_r_dmgay(double dd, int mode, int ndigits,
             int *decpt, int *sign, char **rve,
             char *buf, size_t blen);
}

namespace DAVID_GAY_GFMT {

typedef int Long;

/// Adjusts precision of floating point numbers.
int
g_fmt(char *b, double x, int prec)
{
  int i, k;
  char *s;
  int decpt, j, sign;
  char *b0, sbuf[400], *se;

  b0 = b;
  if (! x) {        // Report -0 as '0'
    *b++ = '0';
    *b = 0;
    goto done;
  }
  s = dtoa_r_dmgay(x, prec ? 2 : 0, prec,
                   &decpt, &sign, &se, sbuf, sizeof(sbuf));
  if (sign)
    *b++ = '-';
  if (decpt == 9999) /* Infinity or Nan */ {
    if (*s == 'N')	/* NaN */
      b = b0;
    while((*b = *s++))
      b++;
    goto done;
  }
  if (decpt <= -4 || decpt > se - s + (s[1] ? 5 : 4)) {
    *b++ = *s++;
    if (*s) {
      *b++ = '.';
      while((*b = *s++))
        b++;
    }
    *b++ = 'e';
    /* sprintf(b, "%+.2d", decpt - 1); */
    if (--decpt < 0) {
      *b++ = '-';
      decpt = -decpt;
    }
    else
      *b++ = '+';
    for(j = 2, k = 10; 10*k <= decpt; j++, k *= 10);
    for(;;) {
      i = decpt / k;
      *b++ = i + '0';
      if (--j <= 0)
        break;
      decpt -= i*k;
      decpt *= 10;
    }
    *b = 0;
  }
  else if (decpt <= 0) {
    *b++ = '0';
    *b++ = '.';
    for(; decpt < 0; decpt++)
      *b++ = '0';
    while((*b = *s++))
      b++;
  }
  else {
    while((*b = *s++)) {
      b++;
      if (--decpt == 0 && *s)
        *b++ = '.';
    }
    for(; decpt > 0; decpt--)
      *b++ = '0';
    *b = 0;
  }
done:
  return b - b0;
}

char Bf[400];

char *
gfmt(double x, int prec)
{
  g_fmt(Bf, x, prec);
  return Bf;
}

}  // namespace DAVID_GAY_GFMT


namespace NL_LIB_GFMT {

/// Adjusts precision of floating point numbers.
/// If too slow, should be replaced by fmt::format
/// or std::to_chars (C++ 17).
void gfmt(char *b, size_t sz, double x, int prec) {
  switch (std::fpclassify(x))
  {
  case FP_INFINITE:
    if (x>=0)
      std::strncpy(b, "Infinity", sz);
    else
      std::strncpy(b, "-Infinity", sz);
    return;
  case FP_NAN:
    std::strncpy(b, "NaN", sz);
    return;
    //			case FP_NORMAL:
    //					return "normal";
    //			case FP_SUBNORMAL:
    //					return "subnormal";
  case FP_ZERO:
    std::strncpy(b, "0", sz);
    return;
  default:
    if (x == std::floor(x)) {   // integer
      auto x_by_1e5 = x / 1e5;
      if (x_by_1e5 != std::round(x_by_1e5)) {
        snprintf(b, sz, "%.*g", prec, x);
        return;
      }
    }
    auto l10 = std::log10(std::fabs(x));
    auto p0 = int(-std::floor(l10));
    assert(std::fabs(x * pow(10.0, p0)) > 0.9999);
    assert(std::fabs(x * pow(10.0, p0)) < 10.0001);
    --p0;
    int i=1;         // we already have 1 digit
    for (; i<prec; ++i) {
      auto xl10 = x * pow(10.0, p0 + i);
      if (xl10 == std::floor(xl10))
        break;
    }
    do {
      if (b) std::snprintf(b, sz, "%.*g", i, x);
      if (i>=prec)
        break;
      auto x1 = std::stod(b);
      if (x==x1)
        return;
      ++i;
    } while (i<=prec);
  }
}

}  // namespace NL_LIB_GFMT
