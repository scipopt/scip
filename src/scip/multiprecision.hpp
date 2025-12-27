/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_rational.h
 * @brief  definition of wrapper class for rational numbers
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MULTIPRECISION_HPP_
#define __SCIP_MULTIPRECISION_HPP_

#include <numeric>
#include <vector>
#include <string>

#ifdef SCIP_WITH_BOOST
#include <boost/multiprecision/number.hpp>

#ifdef SCIP_WITH_GMP
#include <boost/multiprecision/gmp.hpp>
using namespace boost::multiprecision;
namespace scip
{
   using Rational = boost::multiprecision::number<gmp_rational, et_off>;
   using Integer = number<gmp_int, et_off>;
}
#else
#include <boost/multiprecision/cpp_int.hpp>
using namespace boost::multiprecision;
namespace scip
{
   using Rational = cpp_rational;
   using Integer = cpp_int;
}
#endif

namespace scip
{
   using sparsevec = std::vector<SCIP_RATIONAL>;
}
#else
namespace scip
{
   /**@todo replace this dummy class with a direct interface to GMP for systems where Boost is not available */
   class Rational
   {
   public:
      double val;
      Rational(){};
      /// copy constructor
      Rational(const Rational& r){};
      Rational(const char* s){val = atof(s);}
      Rational& operator=(const Rational& r){val = 0.0; return *this;};
      Rational(const long long num, const long long den){val = 0.0;};
      Rational(const Rational& num, const Rational& den){val = 0.0;};
      Rational(const std::string& s){val = 0.0;};
      Rational(const double d){val = 0.0;};
      inline int sign() const { return 0; }
      inline bool is_zero() const { return this->val == 0; }
      inline std::string str() const { return ""; }
      Rational& operator=(const double& r){val = r; return *this;};
      Rational& operator=(const int& i){val = i; return *this;};
      Rational& operator+(const Rational& r){return *this;};
      Rational& operator+(const double& r){return *this;};
      Rational& operator-(const Rational& r){return *this;};
      Rational& operator-(const double& r){return *this;};
      Rational& operator-(){return *this;};
      Rational& operator*(const Rational& r){return *this;};
      Rational& operator*(const double& r){return *this;};
      Rational& operator*=(const double& r){return *this;};
      Rational& operator+=(const Rational& r){return *this;};
      Rational& operator+=(const double& r){return *this;};
      Rational& operator-=(const Rational& r){return *this;};
      Rational& operator-=(const double& r){return *this;};
      Rational& operator/(const Rational& r){return *this;};
      Rational& operator/(const double& r){return *this;};
      friend std::ostream& operator<<(std::ostream& os, const Rational& q){return os;};
      friend double operator/(const double d, const Rational& r){return 0;};
      friend bool operator<(const Rational& r, const double& d){return true;};
      friend bool operator<(const Rational& r, const Rational& d){return true;};
      friend bool operator>(const Rational& r, const double& d){return true;};
      friend bool operator>(const Rational& r, const Rational& d){return true;};
      friend bool operator<=(const Rational& r, const double& d){return true;};
      friend bool operator<=(const Rational& r, const Rational& d){return true;};
      friend bool operator>=(const Rational& r, const double& d){return true;};
      friend bool operator>=(const Rational& r, const Rational& d){return true;};
      friend bool operator==(const Rational& r, const double& d){return true;};
      friend bool operator==(const Rational& r, const Rational& d){return true;};
   };
   using Integer = Rational;
   using sparsevec = std::vector<SCIP_RATIONAL>;
   INLINE Rational& abs(Rational& r){return r;}
   INLINE Rational& max(Rational& r1, Rational& r2){return r1;}
   INLINE Rational& min(Rational& r1, Rational& r2){return r1;}
   INLINE SCIP_Longint denominator(Rational& r){return 0L;}
   INLINE SCIP_Longint numerator(Rational& r){return 0L;}
}
#endif

#endif
