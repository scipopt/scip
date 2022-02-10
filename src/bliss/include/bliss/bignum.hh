#ifndef BLISS_BIGNUM_HH
#define BLISS_BIGNUM_HH

/*
  Copyright (c) 2003-2015 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.

  This file is part of bliss.

  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

#if defined(BLISS_USE_GMP)
#include <gmp.h>
#endif

#include <cstdlib>
#include <cstdio>
#include "bliss/defs.hh"

namespace bliss {

/**
 * \brief A very simple class for big integers (or approximation of them).
 *
 * If the compile time flag BLISS_USE_GMP is set,
 * then the GNU Multiple Precision Arithmetic library (GMP) is used to
 * obtain arbitrary precision, otherwise "long double" is used to
 * approximate big integers.
 */

#if defined(BLISS_USE_GMP)

class BigNum
{
  mpz_t v;
public:
  /**
   * Create a new big number and set it to zero.
   */
  BigNum() {mpz_init(v); }

  /**
   * Destroy the number.
   */
  ~BigNum() {mpz_clear(v); }

  /**
   * Set the number to \a n.
   */
  void assign(const int n) {mpz_set_si(v, n); }

  /**
   * Multiply the number with \a n.
   */
  void multiply(const int n) {mpz_mul_si(v, v, n); }

  /**
   * Print the number in the file stream \a fp.
   */
  size_t print(FILE* const fp) const {return mpz_out_str(fp, 10, v); }
};

#else

class BigNum
{
  long double v;
public:
  /**
   * Create a new big number and set it to zero.
   */
  BigNum(): v(0.0) {}

  /**
   * Set the number to \a n.
   */
  void assign(const int n) {v = (long double)n; }

  /**
   * Multiply the number with \a n.
   */
  void multiply(const int n) {v *= (long double)n; }

  /**
   * Print the number in the file stream \a fp.
   */
  size_t print(FILE* const fp) const {return fprintf(fp, "%Lg", v); }
};

#endif

} //namespace bliss

#endif
