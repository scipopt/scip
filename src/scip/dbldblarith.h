/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dbldblarith.h
 * @brief  defines macros for basic operations in double-double arithmetic giving roughly twice the precision of a double
 * @author Robert Lion Gottwald
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef _SCIP_DBLDBL_ARITH_
#define _SCIP_DBLDBL_ARITH_

#include "math.h"

#define __SCIPdbldblSplit(rhi, rlo, x) \
    do { \
       const double __tmp_split_dbl = 134217729.0 * (x); \
       (rhi) = __tmp_split_dbl - (__tmp_split_dbl - (x)); \
       (rlo) = (x) - (rhi);\
    } while(0)

/** multiply two floating point numbers, both given by one double, and return the result as two doubles. */
#define SCIPdbldblProd(rhi, rlo, a, b) \
    do { \
        double __tmp_dbldbl_prod_ahi; \
        double __tmp_dbldbl_prod_alo; \
        double __tmp_dbldbl_prod_bhi; \
        double __tmp_dbldbl_prod_blo; \
        __SCIPdbldblSplit(__tmp_dbldbl_prod_ahi, __tmp_dbldbl_prod_alo, a); \
        __SCIPdbldblSplit(__tmp_dbldbl_prod_bhi, __tmp_dbldbl_prod_blo, b); \
        (rhi) = (a) * (b); \
        (rlo) = __tmp_dbldbl_prod_alo * __tmp_dbldbl_prod_blo - \
           ((((rhi) - __tmp_dbldbl_prod_ahi * __tmp_dbldbl_prod_bhi) \
           - __tmp_dbldbl_prod_alo * __tmp_dbldbl_prod_bhi) \
           - __tmp_dbldbl_prod_ahi * __tmp_dbldbl_prod_blo); \
    } while(0)

/** square a floating point number given by one double and return the result as two doubles. */
#define SCIPdbldblSquare(rhi, rlo, a) \
    do { \
        double __tmp_dbldbl_square_ahi; \
        double __tmp_dbldbl_square_alo; \
        __SCIPdbldblSplit(__tmp_dbldbl_square_ahi, __tmp_dbldbl_square_alo, a); \
        (rhi) = (a) * (a); \
        (rlo) = __tmp_dbldbl_square_alo * __tmp_dbldbl_square_alo - \
           ((((rhi) - __tmp_dbldbl_square_ahi * __tmp_dbldbl_square_ahi) \
           - 2.0 * __tmp_dbldbl_square_alo * __tmp_dbldbl_square_ahi)); \
    } while(0)

/** add two floating point numbers, both given by one double, and return the result as two doubles. */
#define SCIPdbldblSum(rhi, rlo, a, b) \
    do { \
        double __tmp1_dbldbl_sum; \
        double __tmp2_dbldbl_sum; \
        __tmp2_dbldbl_sum = (a) + (b); \
        __tmp1_dbldbl_sum = __tmp2_dbldbl_sum - (a); \
        (rlo) = ((a) - (__tmp2_dbldbl_sum - __tmp1_dbldbl_sum)) + ((b) - __tmp1_dbldbl_sum); \
        (rhi) = __tmp2_dbldbl_sum; \
    } while(0)

/** divide two floating point numbers, both given by one double, and return the result as two doubles. */
#define SCIPdbldblDiv(rhi, rlo, a, b) \
    do { \
       double __tmp_dbldbl_div_hi; \
       double __tmp_dbldbl_div_lo; \
       double __estim_dbldbl_div = (a)/(b); \
       SCIPdbldblProd(__tmp_dbldbl_div_hi, __tmp_dbldbl_div_lo, b, __estim_dbldbl_div); \
       SCIPdbldblSum21(__tmp_dbldbl_div_hi, __tmp_dbldbl_div_lo, __tmp_dbldbl_div_hi, __tmp_dbldbl_div_lo, -(a)); \
       __tmp_dbldbl_div_hi /= (b); \
       __tmp_dbldbl_div_lo /= (b); \
       SCIPdbldblSum21(rhi, rlo, -__tmp_dbldbl_div_hi, -__tmp_dbldbl_div_lo, __estim_dbldbl_div); \
    } while(0)

/** add two floating point numbers, the first is given by two doubles, the second is given by one double,
 *  and return the result as two doubles.
 */
#define SCIPdbldblSum21(rhi, rlo, ahi, alo, b) \
   do { \
      double __tmp_dbldbl_sum21_hi; \
      double __tmp_dbldbl_sum21_lo; \
      SCIPdbldblSum(__tmp_dbldbl_sum21_hi, __tmp_dbldbl_sum21_lo, ahi, b); \
      (rlo) = __tmp_dbldbl_sum21_lo + (alo); \
      (rhi) = __tmp_dbldbl_sum21_hi; \
   } while(0)


/** multiply two floating point numbers, the first is given by two doubles, the second is given by one double,
 *  and return the result as two doubles.
 */
#define SCIPdbldblProd21(rhi, rlo, ahi, alo, b) \
    do { \
       double __tmp_dbldbl_prod21_hi; \
       double __tmp_dbldbl_prod21_lo; \
       SCIPdbldblProd(__tmp_dbldbl_prod21_hi, __tmp_dbldbl_prod21_lo, ahi, b); \
       (rlo) = (alo) * (b) + __tmp_dbldbl_prod21_lo; \
       (rhi) = __tmp_dbldbl_prod21_hi; \
    } while(0)

/** divide two floating point numbers, the first is given by one double, the second is given by two doubles,
 *  and return the result as two doubles.
 */
#define SCIPdbldblDiv12(rhi, rlo, a, bhi, blo) \
    do { \
       double __tmp_dbldbl_div12_hi; \
       double __tmp_dbldbl_div12_lo; \
       double __estim_dbldbl_div12 = (a)/(bhi); \
       SCIPdbldblProd21(__tmp_dbldbl_div12_hi, __tmp_dbldbl_div12_lo, bhi, blo, __estim_dbldbl_div12); \
       SCIPdbldblSum21(__tmp_dbldbl_div12_hi, __tmp_dbldbl_div12_lo, __tmp_dbldbl_div12_hi, __tmp_dbldbl_div12_lo, -(a)); \
       __tmp_dbldbl_div12_hi /= (bhi); \
       __tmp_dbldbl_div12_lo /= (bhi); \
       SCIPdbldblSum21(rhi, rlo, -__tmp_dbldbl_div12_hi, -__tmp_dbldbl_div12_lo, __estim_dbldbl_div12); \
    } while(0)


/** divide two floating point numbers, the first is given by two doubles, the second is given by one double,
 *  and return the result as two doubles.
 */
#define SCIPdbldblDiv21(rhi, rlo, ahi, alo, b) \
   do { \
      double __tmp_dbldbl_div21_hi; \
      double __tmp_dbldbl_div21_lo; \
      double __estim_dbldbl_div21_hi; \
      double __estim_dbldbl_div21_lo; \
      __estim_dbldbl_div21_hi = (ahi)/(b); \
      __estim_dbldbl_div21_lo = (alo)/(b); \
      SCIPdbldblProd21(__tmp_dbldbl_div21_hi, __tmp_dbldbl_div21_lo, __estim_dbldbl_div21_hi, __estim_dbldbl_div21_lo, b); \
      SCIPdbldblSum22(__tmp_dbldbl_div21_hi, __tmp_dbldbl_div21_lo, __tmp_dbldbl_div21_hi, __tmp_dbldbl_div21_lo, -(ahi), -(alo)); \
      __tmp_dbldbl_div21_hi /= (b); \
      __tmp_dbldbl_div21_lo /= (b); \
      SCIPdbldblSum22(rhi, rlo, __estim_dbldbl_div21_hi, __estim_dbldbl_div21_lo, -__tmp_dbldbl_div21_hi, -__tmp_dbldbl_div21_lo); \
   } while(0)

/** multiply two floating point numbers, both given by two doubles, and return the result as two doubles. */
#define SCIPdbldblProd22(rhi, rlo, ahi, alo, bhi, blo) \
   do { \
      double __tmp_dbldbl_prod22_hi; \
      double __tmp_dbldbl_prod22_lo; \
      SCIPdbldblProd(__tmp_dbldbl_prod22_hi, __tmp_dbldbl_prod22_lo, ahi, bhi); \
      SCIPdbldblSum21(__tmp_dbldbl_prod22_hi, __tmp_dbldbl_prod22_lo, \
                   __tmp_dbldbl_prod22_hi, __tmp_dbldbl_prod22_lo, (alo) * (bhi)); \
      SCIPdbldblSum21(rhi, rlo, \
                   __tmp_dbldbl_prod22_hi, __tmp_dbldbl_prod22_lo, (ahi) * (blo)); \
   } while(0)

/** add two floating point numbers, both given by two doubles, and return the result as two doubles. */
#define SCIPdbldblSum22(rhi, rlo, ahi, alo, bhi, blo) \
   do { \
      double __tmp_dbldbl_sum22_hi; \
      double __tmp_dbldbl_sum22_lo; \
      SCIPdbldblSum21(__tmp_dbldbl_sum22_hi, __tmp_dbldbl_sum22_lo, ahi, alo, bhi); \
      SCIPdbldblSum21(rhi, rlo, __tmp_dbldbl_sum22_hi, __tmp_dbldbl_sum22_lo, blo); \
   } while(0)

/** square a floating point number given by two doubles and return the result as two doubles. */
#define SCIPdbldblSquare2(rhi, rlo, ahi, alo) \
   do { \
      double __tmp_dbldbl_square2_hi; \
      double __tmp_dbldbl_square2_lo; \
      SCIPdbldblSquare(__tmp_dbldbl_square2_hi, __tmp_dbldbl_square2_lo, (ahi)); \
      SCIPdbldblSum21(rhi, rlo, __tmp_dbldbl_square2_hi, __tmp_dbldbl_square2_lo, 2 * (ahi) * (alo)); \
   } while(0)

/** divide two floating point numbers, both given by two doubles, and return the result as two doubles. */
#define SCIPdbldblDiv22(rhi, rlo, ahi, alo, bhi, blo) \
   do { \
      double __tmp_dbldbl_div22_hi; \
      double __tmp_dbldbl_div22_lo; \
      double __estim_dbldbl_div22_hi = (ahi) / (bhi); \
      double __estim_dbldbl_div22_lo = (alo) / (bhi); \
      SCIPdbldblProd22(__tmp_dbldbl_div22_hi, __tmp_dbldbl_div22_lo, \
                    bhi, blo, __estim_dbldbl_div22_hi, __estim_dbldbl_div22_lo); \
      SCIPdbldblSum22(__tmp_dbldbl_div22_hi, __tmp_dbldbl_div22_lo, \
                   __tmp_dbldbl_div22_hi, __tmp_dbldbl_div22_lo, -(ahi), -(alo)); \
      __tmp_dbldbl_div22_hi /= (bhi); \
      __tmp_dbldbl_div22_lo /= (bhi); \
      SCIPdbldblSum22(rhi, rlo, __estim_dbldbl_div22_hi, __estim_dbldbl_div22_lo, \
                                -__tmp_dbldbl_div22_hi, -__tmp_dbldbl_div22_lo); \
   } while(0)


/** take the square root of a floating point number given by one double and return the result as two doubles. */
#define SCIPdbldblSqrt(rhi, rlo, a) \
   do { \
      double __estim_dbldbl_sqrt = sqrt(a); \
      SCIPdbldblDiv(rhi, rlo, a, __estim_dbldbl_sqrt); \
      SCIPdbldblSum21(rhi, rlo, rhi, rlo, __estim_dbldbl_sqrt); \
      (rhi) *= 0.5; \
      (rlo) *= 0.5; \
   } while(0)


/** take the square root of a floating point number given by two doubles and return the result as two doubles. */
#define SCIPdbldblSqrt2(rhi, rlo, ahi, alo) \
   do { \
      double __estim_dbldbl_sqrt2 = sqrt(ahi); \
      SCIPdbldblDiv21(rhi, rlo, ahi, alo, __estim_dbldbl_sqrt2); \
      SCIPdbldblSum21(rhi, rlo, rhi, rlo, __estim_dbldbl_sqrt2); \
      (rhi) *= 0.5; \
      (rlo) *= 0.5; \
   } while(0)

/** compute the absolute value of the floating point number given by two doubles */
#define SCIPdbldblAbs2(rhi, rlo, ahi, alo) \
   do { \
      if( ahi < 0.0 ) \
      { \
         (rhi) = -(ahi); \
         (rlo) = -(alo); \
      } \
      else \
      { \
         (rhi) = (ahi); \
         (rlo) = (alo); \
      } \
   } while(0)

#endif
