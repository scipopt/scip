#pragma once

/*
  Copyright (c) 2003-2021 Tommi Junttila
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
#else
#include <vector>
#include <string>
#endif

#include <cstdlib>
#include <cstdio>
#include "bliss/defs.hh"


namespace bliss {

/**
 * \brief A simple wrapper class for non-negative big integers (or approximation of them).
 *
 * If the compile time flag BLISS_USE_GMP is set,
 * then the GNU Multiple Precision Arithmetic library (GMP) is used to
 * obtain arbitrary precision.
 * Otherwise, if the compile time flag BLISS_BIGNUM_APPROX is set,
 * a "long double" is used to approximate a big integer.
 * Otherwise, by default, a big integer is represented by a product of
 * integer-sized factors.
 */

#if defined(BLISS_USE_GMP)

class BigNum
{
  mpz_t v;
public:
  /**
   * \brief Create a new big number and set it to zero.
   */
  BigNum() {mpz_init(v); }

  /**
   * \brief Destroy the number.
   */
  ~BigNum() {mpz_clear(v); }

  /**
   * \brief Set the number to \a n.
   */
  void assign(unsigned int n) {mpz_set_ui(v, n); }

  /**
   * \brief Multiply the number with \a n.
   */
  void multiply(unsigned int n) {mpz_mul_ui(v, v, n); }

  /**
   * \brief Print the number in the file stream \a fp.
   */
  size_t print(FILE* const fp) const {return mpz_out_str(fp, 10, v); }

  /**
   * Get a copy of the internal GNU GMP integer.
   * The caller is responsible for calling mpz_init before,
   * and mpz_clear afterwards on the \a result variable.
   */
  void get(mpz_t& result) const {mpz_set(result, v); }
};

#elif defined(BLISS_BIGNUM_APPROX)

class BigNum
{
  long double v;
public:
  /**
   * \brief Create a new big number and set it to zero.
   */
  BigNum(): v(0.0) {}

  /**
   * \brief Set the number to \a n.
   */
  void assign(unsigned int n) {v = (long double)n; }

  /**
   * \brief Multiply the number with \a n.
   */
  void multiply(unsigned int n) {v *= (long double)n; }

  /**
   * \brief Print the number in the file stream \a fp.
   */
  size_t print(FILE* const fp) const {return fprintf(fp, "%Lg", v); }
};

#else


class BigNum
{
  /* This is a version that does not actually compute the number
   * but rather only stores the factor integers.
   */
  std::vector<unsigned int> factors;
public:
  /**
   * Create a new big number and set it to zero.
   */
  BigNum() {
    factors.push_back(0);
  }

  /**
   * Destroy the number.
   */
  ~BigNum() {}

  /**
   * Set the number to \a n.
   */
  void assign(unsigned int n) {
    factors.clear();
    factors.push_back(n);
  }

  /**
   * Multiply the number with \a n.
   */
  void multiply(unsigned int n) {
    factors.push_back(n);
  }

  /**
   * Print the number in the file stream \a fp.
   * In the current version, the returned number of characters printed,
   * is incorrect (either -1 or 0).
   */
  size_t print(FILE* const fp) const {
    assert(not factors.empty());
    size_t r = 0;
    /*
    const char* sep = "";
    for(int v: factors) {
      r += fprintf(fp, "%s%d", sep, v);
      sep = "*";
    }
    */
    for(char d: to_string())
      r += fprintf(fp, "%c", d);
    return r;
  }

  /**
   * Get a reference to the factors vector.
   */
  const std::vector<unsigned int>& get_factors() const {
    return factors;
  }

  /**
   * Get the tring representation of the number.
   * Unoptimized, uses an elementary school algorithm to multiply the factors.
   */
  std::string to_string() const {
    // Base 100 result, in reverse order
    std::vector<unsigned int> result;
    result.push_back(1);
    for(unsigned int factor: factors) {
      std::vector<unsigned int> summand;
      unsigned int offset = 0;
      while(factor != 0) {
        const unsigned int multiplier = factor % 100;
        // Multiplication by a "digit"
        std::vector<unsigned int> product;
        for(unsigned int i = 0; i < offset; i++)
          product.push_back(0);
        unsigned int carry = 0;
        for(unsigned int digit: result) {
          unsigned int v = digit * multiplier + carry;
          product.push_back(v % 100);
          carry = v / 100;
        }
        if(carry > 0)
          product.push_back(carry);
        // Addition
        add(summand, product);
        // Next "digit" in factor
        factor = factor / 100;
        offset++;
      }
      result = summand;
    }
    return _string(result);
  }

protected:
  static void add(std::vector<unsigned int>& num, const std::vector<unsigned int>& summand) {
    unsigned int carry = 0;
    unsigned int i = 0;
    while(i < num.size() and i < summand.size()) {
      const unsigned int v = carry + num[i] + summand[i];
      num[i] = v % 100;
      carry = v / 100;
      i++;
    }
    while(i < summand.size()) {
      const unsigned int v = carry + summand[i];
      num.push_back(v % 100);
      carry = v / 100;
      i++;
    }
    while(i < num.size()) {
      const unsigned int v = carry + num[i];
      num[i] = v % 100;
      carry = v / 100;
      i++;
    }
    if(carry != 0)
      num.push_back(carry);
  }


  static std::string _string(const std::vector<unsigned int> n) {
    const char digits[] = {'0','1','2','3','4','5','6','7','8','9'};
    std::string r;
    bool first = true;
    for(auto it = n.crbegin(); it != n.crend(); it++) {
      unsigned int digit = *it;
      unsigned int high = digit / 10;
      if(not first or high > 0)
        r.push_back(digits[high]);
      first = false;
      r.push_back(digits[digit % 10]);
    }
    return r;
  }

};

#endif

} //namespace bliss
