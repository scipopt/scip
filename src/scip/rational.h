/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving raint Integer Programs                      */
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
#include <stdlib.h>
#include "scip/def.h"
#include "scip/intervalarith.h"
#include "scip/mem.h"
#include "scip/type_misc.h"
#ifdef SCIP_WITH_GMP
#include <gmp.h>
#endif
#ifdef SCIP_WITH_ZIMPL
#include "zimpl/numb.h"
#endif
#ifdef SCIP_WITH_MPFR
#include <mpfr.h>
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

/** Allocate and create a rational from a string in the format, e.g. "12/35" */
EXTERN
SCIP_RETCODE RcreateString(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,           /**< pointer to the rational to create */
   char*                 desc                /**< the String describing the rational */
   );

/** create an array of rationals */
EXTERN
SCIP_RETCODE RcreateArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      rational,           /**< pointer to the array to create */
   int                   size                /**< the size of the array */
   );

/** create an array of rationals */
EXTERN
SCIP_RETCODE RcreateArrayTemp(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational***      rational,           /**< pointer to the array to create */
   int                   size                /** the size of the array */
   );

/** copy an array of rationals */
SCIP_RETCODE RcopyArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      target,             /**< address to copy to */
   SCIP_Rational**       src,                /**< src array */
   int                   len                 /**< size of src array */
   );

/** create a copy of a rational */
EXTERN
SCIP_RETCODE Rcopy(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,           /**< pointer to the rational to create */
   SCIP_Rational*        src                 /**< rational to copy */
   );

EXTERN
SCIP_RETCODE RcreateTemp(
   BMS_BUFMEM*           buf,
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   );

EXTERN
SCIP_RETCODE RcreateNoMem(
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   );

EXTERN
SCIP_RETCODE Rcreate(
   BMS_BLKMEM*           blkmem,
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   );

#ifdef SCIP_WITH_GMP
/** create a rational from an mpq_t */
EXTERN
SCIP_Rational* RcreateGMP(
   BMS_BLKMEM*           mem,                /**< block memory */
    mpq_t           numb                /**< the mpq_rational */
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

EXTERN
void RsetGMPArray(
   mpq_t*                res,
   SCIP_Rational**       src,
   int                   len 
   );

void RsetArrayGMP(
   SCIP_Rational**       res,
   mpq_t*                src,
   int                   len
   );

EXTERN
void RclearGMPArray(
   mpq_t*                ar,
   int                   len
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

/** free an array of rationals */
EXTERN
void RdeleteArrayTemp(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational***      array,              /**< pointer to the array */
   int                   size                /**< size of the array */
   );

/** set a rational to the value of another rational */
EXTERN
void Rset(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*   src                 /**< the src */
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
    char*           desc                /**< the string describing the rational */
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

/* transform rational into canonical form */
EXTERN
void Rcanonicalize(
   SCIP_Rational*        r                   /**< rational to put in canonical form */
   );

/** add two rationals and save the result in res*/
EXTERN
void Radd(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op1,                /**< first operand */
    SCIP_Rational*  op2                 /**< second operand */
   );

/** add a rational and a real and save the result in res*/
EXTERN
void RaddReal(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   );

/** subtract two rationals and save the result in res*/
EXTERN
void Rdiff(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op1,                /**< first operand */
    SCIP_Rational*  op2                 /**< second operand */
   );

/** subtract a rational and a real and save the result in res*/
EXTERN
void RdiffReal(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   );

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) of two rationals */
EXTERN
void RrelDiff(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  val1,               /**< first value to be compared */
    SCIP_Rational*  val2                /**< second value to be compared */
   );

/** multiply two rationals and save the result in res*/
EXTERN
void Rmult(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op1,                /**< first operand */
    SCIP_Rational*  op2                 /**< second operand */
   );

/** multiply a rational and a real and save the result in res*/
EXTERN
void RmultReal(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/** divide two rationals and save the result in res*/
EXTERN
void Rdiv(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op1,                /**< first operand */
    SCIP_Rational*  op2                 /**< second operand */
   );

/** divide a rational and a real and save the result in res*/
EXTERN
void RdivReal(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/* Computes res += op1 * op2 and saves the result in res */
void RaddProd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   );

/* Computes res -= op1 * op2 and saves the result in res */
void RdiffProd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   );

/** set res to -op */
EXTERN
void Rneg(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op                  /**< operand */
   );

/** set res to Abs(op) */
EXTERN
void Rabs(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op                  /**< operand */
   );

/** set res to 1/op */
EXTERN
void Rinv(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op                  /**< operand */
   );

/** compute the minimum of two rationals */
EXTERN
void Rmin(
   SCIP_Rational*        ret,                /**< the result */
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** compute the maximum of two rationals */
EXTERN
void Rmax(
   SCIP_Rational*        ret,                /**< the result */
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/*
 * Comparisoon methods
 */

/** check if two rationals are equal */
EXTERN
SCIP_Bool RisEqual(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if two rationals are equal */
EXTERN
SCIP_Bool RisAbsEqual(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );


/** check if a rational and a real are equal */
EXTERN
SCIP_Bool RisEqualReal(
    SCIP_Rational*  r1,                 /**< the rational */
   SCIP_Real             r2                  /**< the real */
   );

/** check if real approx of rational and a real are equal */
EXTERN
SCIP_Bool RApproxEqualReal(
   SCIP_Rational*        r1,                 /**< the rational */
   SCIP_Real             r2                  /**< the real */
   );

/** check if the first rational is greater than the second*/
EXTERN
SCIP_Bool RisGT(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is greater than the second*/
EXTERN
SCIP_Bool RisAbsGT(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is smaller than the second*/
EXTERN
SCIP_Bool RisLT(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is smaller or equal than the second*/
EXTERN
SCIP_Bool RisLE(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is greater or equal than the second*/
EXTERN
SCIP_Bool RisGE(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the rational is zero */
EXTERN
SCIP_Bool RisZero(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is positive */
EXTERN
SCIP_Bool RisPositive(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is negative */
EXTERN
SCIP_Bool RisNegative(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is positive infinity */
EXTERN
SCIP_Bool RisInfinity(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is negative infinity */
EXTERN
SCIP_Bool RisNegInfinity(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is of infinite value */
EXTERN
SCIP_Bool RisAbsInfinity(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is of infinite value */
EXTERN
SCIP_Bool RisIntegral(
    SCIP_Rational*  r                   /**< the rational to check */
   );

EXTERN
SCIP_Bool RisFpRepresentable(
    SCIP_Rational*    r
   );

/*
 * Printing/Conversion methods
 */

/** convert a Rational to a string for printing, returns the number of copied characters.
 * If return value is equal to strlen, it means the string was truncated.
 */
int RtoString(
   SCIP_Rational*        r,                  /**< the rational to print */
   char*                 str,                /**< the string to save the rational in */
   int                   strlen              /**< maximal length that can be copied to str */
   );

const char* RgetString(
   SCIP_Rational*        r
   );

/** return the strlen of a rational number */
SCIP_Longint Rstrlen(
   SCIP_Rational*        r                /** rational to consider */
   );

/** print a rational to command line (for debugging) */
void Rprint(
   SCIP_Rational*        r,                  /**< the rational to print */
   FILE*                 file                /**< file to print to (or NULL for std output) */
   );

/** print rational to file using message handler */
void Rmessage(
   SCIP_MESSAGEHDLR*     msg,                /**< message handler */
   FILE*                 file,               /**< file pointer */
   SCIP_Rational*        r                   /**< the rational to print */
   );

/** return approximation of Rational as SCIP_Real */
EXTERN
SCIP_Real RgetRealRelax(
   SCIP_Rational*  r,                  /**< the rational to convert */
   SCIP_ROUNDMODE        roundmode           /**< rounding direction (not really working yet) */
   );

/** return approximation of Rational as SCIP_Real */
EXTERN
SCIP_Real RgetRealApprox(
    SCIP_Rational*  r                   /**< the rational to convert */
   );

EXTERN
void Rround(
   SCIP_Rational*        retval,             /**< the resulting rounded integer */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   );

/** round rational to next integer in direction of roundmode */
EXTERN
SCIP_Bool RroundInteger(
   long int*                 retval,             /**< the resulting rounded lon int */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
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
    SCIP_Rational*  val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
SCIP_RETCODE SCIPrationalarrayIncVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
    SCIP_Rational*  incval              /**< value to increase array index */
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

extern
void testRuntimesRational(
   );

#ifdef __cplusplus
}
#endif

#endif
