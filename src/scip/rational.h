/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip/rational.h
 * @ingroup INTERNALAPI
 * @brief  rational wrapper
 * @author Leon Eifler
  */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RATIONAL_H__
#define __SCIP_RATIONAL_H__

#include <stdbool.h>
#include "scip/def.h"
#include "scip/intervalarith.h"
#include "scip/mem.h"
#include "scip/type_misc.h"
#ifdef SCIP_WITH_GMP
#include <gmp.h>
#ifdef SCIP_WITH_ZIMPL
#include "zimpl/numb.h"
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum SCIP_RoundModeR
{
   SCIP_ROUND_DOWNWARDS,
   SCIP_ROUND_UPWARDS,
   SCIP_ROUND_NEAREST
};

struct SCIP_RationalArray
{
   BMS_BLKMEM*           blkmem;             /**< block memory that stores the vals array */
   SCIP_Rational**       vals;               /**< array values */
   int                   valssize;           /**< size of vals array */
   int                   firstidx;           /**< index of first element in vals array */
   int                   minusedidx;         /**< index of first non zero element in vals array */
   int                   maxusedidx;         /**< index of last non zero element in vals array */
};

typedef enum SCIP_RoundModeR SCIP_ROUNDMODER;

/*
 * Creation methods
 */

/** Allocate and create a rational from nominator and denominator */
EXTERN
SCIP_Rational* RcreateInt(
   BMS_BLKMEM*           mem,
   int                   nom,                /**< the nominator */
   int                   denom               /**< the denominator */
   );

/*
 * Creation methods
 */

/** Allocate and create a rational from a scip_real */
EXTERN
SCIP_Rational* RcreateReal(
   BMS_BLKMEM*           mem,
   const SCIP_Real       num                 /**< the scip_real */
   );

/** Allocate and create a rational from a string in the format, e.g. "12/35" */
EXTERN
SCIP_Rational* RcreateString(
   BMS_BLKMEM*           mem,                /**< block memory */
   const char*           desc                /**< the String describing the rational */
   );

/** create an array of rationals */
EXTERN
SCIP_Rational** RcreateArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   int                   size                /**< the size of the array */
   );

/** copy an array of rationals */
void* RcopyArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      target,             /**< address to copy to */
   SCIP_Rational**       src,                /**< src array */
   int                   len                 /**< size of src array */
   );

/** create a copy of a rational */
EXTERN
SCIP_Rational* Rcopy(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational*        src                 /**< rational to copy */
   );

EXTERN
SCIP_Rational* RcreateTemp(
   BMS_BUFMEM*           buf
   );

EXTERN
SCIP_Rational* RcreateNoMem(void);

EXTERN
SCIP_Rational* Rcreate(
   BMS_BLKMEM*           buf
   );

#ifdef WITH_GMP
/** create a rational from an mpq_t */
EXTERN
SCIP_Rational* RcreateGMP(
   BMS_BLKMEM*           mem,                /**< block memory */
   const mpq_t           numb                /**< the mpq_rational */
   );

/** get the underlying mpq_t* */
EXTERN mpq_t* RgetGMP(
   SCIP_Rational*  r                   /**< the rational */
   );

EXTERN
void RsetGMP(
   SCIP_Rational*        r,
   const mpq_t           numb
   );
#endif

/** delete a rational and free the allocated memory */
EXTERN
void RdeleteNoMem(
   SCIP_Rational**       r                   /**< adress of the rational */
   );

/** delete a rational and free the allocated memory */
EXTERN
void Rdelete(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       r                   /**< adress of the rational */
   );

EXTERN
void RdeleteTemp(
   BMS_BUFMEM*           buf,
   SCIP_Rational**       r
   );

/** delete an array of rationals and free the allocated memory */
EXTERN
void RdeleteArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      array,              /**< address of rational array */
   int                   size                /**< size of the array */
   );

/** set a rational to the value of another rational */
EXTERN
void Rset(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*   src                 /**< the src */
   );

/** set a rational to a nom/denom value */
EXTERN
void RsetInt(
   SCIP_Rational*        res,                /**< the result */
   int                   nom,                /**< the nominator */
   int                   denom               /**< the denominator */
   );

/** set a rational to the value described by a string */
EXTERN
void RsetString(
   SCIP_Rational*        res,                /**< the result */
   const char*           desc                /**< the string describing the rational */
   );

/** set a rational to the value of another a real */
EXTERN
void RsetReal(
   SCIP_Rational*        r,
   SCIP_Real             real
   );

/*
 * Computing methods
 */

/** add two rationals and save the result in res*/
EXTERN
void Radd(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  op1,                /**< first operand */
   const SCIP_Rational*  op2                 /**< second operand */
   );

/** add a rational and a real and save the result in res*/
EXTERN
void RaddReal(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   );

/** subtract two rationals and save the result in res*/
EXTERN
void Rdiff(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  op1,                /**< first operand */
   const SCIP_Rational*  op2                 /**< second operand */
   );

/** subtract a rational and a real and save the result in res*/
EXTERN
void RdiffReal(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   );

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) of two rationals */
EXTERN
void RrelDiff(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  val1,               /**< first value to be compared */
   const SCIP_Rational*  val2                /**< second value to be compared */
   );

/** multiply two rationals and save the result in res*/
EXTERN
void Rmult(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  op1,                /**< first operand */
   const SCIP_Rational*  op2                 /**< second operand */
   );

/** multiply a rational and a real and save the result in res*/
EXTERN
void RmultReal(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/** divide two rationals and save the result in res*/
EXTERN
void Rdiv(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  op1,                /**< first operand */
   const SCIP_Rational*  op2                 /**< second operand */
   );

/** divide a rational and a real and save the result in res*/
EXTERN
void RdivReal(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/** set res to -op */
EXTERN
void Rneg(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  op                  /**< operand */
   );

/** set res to Abs(op) */
EXTERN
void Rabs(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  op                  /**< operand */
   );

/** set res to 1/op */
EXTERN
void Rinv(
   SCIP_Rational*        res,                /**< the result */
   const SCIP_Rational*  op                  /**< operand */
   );

/** compute the minimum of two rationals */
EXTERN
void Rmin(
   SCIP_Rational*        ret,                /**< the result */
   const SCIP_Rational*  r1,                 /**< the first rational */
   const SCIP_Rational*  r2                  /**< the second rational */
   );

/** compute the maximum of two rationals */
EXTERN
void Rmax(
   SCIP_Rational*        ret,                /**< the result */
   const SCIP_Rational*  r1,                 /**< the first rational */
   const SCIP_Rational*  r2                  /**< the second rational */
   );

/*
 * Comparisoon methods
 */

/** check if two rationals are equal */
EXTERN
SCIP_Bool RisEqual(
   const SCIP_Rational*  r1,                 /**< the first rational */
   const SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if a rational and a real are equal */
EXTERN
SCIP_Bool RisEqualReal(
   const SCIP_Rational*  r1,                 /**< the rational */
   SCIP_Real             r2                  /**< the real */
   );

/** check if the first rational is greater than the second*/
EXTERN
SCIP_Bool RisGT(
   const SCIP_Rational*  r1,                 /**< the first rational */
   const SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is smaller than the second*/
EXTERN
SCIP_Bool RisLT(
   const SCIP_Rational*  r1,                 /**< the first rational */
   const SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is smaller or equal than the second*/
EXTERN
SCIP_Bool RisLE(
   const SCIP_Rational*  r1,                 /**< the first rational */
   const SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is greater or equal than the second*/
EXTERN
SCIP_Bool RisGE(
   const SCIP_Rational*  r1,                 /**< the first rational */
   const SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the rational is zero */
EXTERN
SCIP_Bool RisZero(
   const SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is positive */
EXTERN
SCIP_Bool RisPositive(
   const SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is negative */
EXTERN
SCIP_Bool RisNegative(
   const SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is positive infinity */
EXTERN
SCIP_Bool RisInfinity(
   const SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is negative infinity */
EXTERN
SCIP_Bool RisNegInfinity(
   const SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is of infinite value */
EXTERN
SCIP_Bool RisAbsInfinity(
   const SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is of infinite value */
EXTERN
SCIP_Bool RisIntegral(
   const SCIP_Rational*  r                   /**< the rational to check */
   );

EXTERN
SCIP_Bool RisFpRepresentable(
   const SCIP_Rational*    r
   );

/*
 * Printing/Conversion methods
 */

/** print a Rational to std out */
void RtoString(
   const SCIP_Rational*  r,                  /**< the rational to print */
   char*                 str
   );

/** print a rational to command line (for debugging) */
void Rprint(
   const SCIP_Rational*  r                   /**< the rational to print */
   );

/** return approximation of Rational as SCIP_Real */
EXTERN
SCIP_Real RgetRealRelax(
   const SCIP_Rational*  r,                  /**< the rational to convert */
   SCIP_ROUNDMODE        roundmode           /**< rounding direction (not really working yet) */
   );

/** return approximation of Rational as SCIP_Real */
EXTERN
SCIP_Real RgetRealApprox(
   const SCIP_Rational*  r                   /**< the rational to convert */
   );

/*
 * Dynamic Arrays todo: use stl to do this
 */

/** creates a dynamic array of real values */
extern
SCIP_RETCODE SCIPrationalarrayCreate(
   SCIP_RATIONALARRAY**  rationalarray,          /**< pointer to store the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a copy of a dynamic array of real values */
extern
SCIP_RETCODE SCIPrationalarrayCopy(
   SCIP_RATIONALARRAY**  rationalarray,      /**< pointer to store the copied real array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_RATIONALARRAY*   sourcerationalarray /**< dynamic real array to copy */
   );

/** frees a dynamic array of real values */
extern
SCIP_RETCODE SCIPrationalarrayFree(
   SCIP_RATIONALARRAY**      rationalarray   /**< pointer to the real array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
SCIP_RETCODE SCIPrationalarrayExtend(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic real array */
extern
SCIP_RETCODE SCIPrationalarrayClear(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic real array */
   );

/** gets value of entry in dynamic array */
extern
void SCIPrationalarrayGetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   idx,                /**< array index to get value for */
   SCIP_Rational*        result              /**< store the result */
   );

/** sets value of entry in dynamic array */
extern
SCIP_RETCODE SCIPrationalarraySetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   const SCIP_Rational*  val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
SCIP_RETCODE SCIPrationalarrayIncVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
   const SCIP_Rational*  incval              /**< value to increase array index */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPrationalarrayGetMinIdx(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic real array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPrationalarrayGetMaxIdx(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic real array */
   );

#ifdef __cplusplus
}
#endif

#endif
