
/*  math.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__MATH_H)
#define __MATH_H

extern double HUGE_VAL;
extern double acos(double );
extern double asin(double );
extern double atan(double );
extern double atan2(double, double );
extern double ceil(double );
extern double cos(double );
extern double cosh(double );
extern double exp(double );
extern double fabs(double );
extern double floor(double );
extern double fmod(double , double );
extern double frexp(double , int *);
extern double ldexp(double , int );
extern double log(double );
extern double log10(double );
extern double modf(double , double *);
extern double nextafter(double , double );
extern double pow(double , double );
extern double sin(double );
extern double sinh(double );
extern double sqrt(double );
extern double tan(double );
extern double tanh(double );
extern long double sqrtl(long double );
extern long double fabsl(long double );
extern int isnan(double );
extern int isfinite(double );
extern int finite(double );

#endif
