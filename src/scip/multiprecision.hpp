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
   Rational& abs(Rational& r){return r;}
   Rational& max(Rational& r1, Rational& r2){return r1;}
   Rational& min(Rational& r1, Rational& r2){return r1;}
   SCIP_Longint denominator(Rational& r){return 0L;}
   SCIP_Longint numerator(Rational& r){return 0L;}
}
#endif

#endif
