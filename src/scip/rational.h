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


typedef enum SCIP_RoundModeR SCIP_ROUNDMODER;

/*
 * Creation methods
 */

/** allocates and creates a rational from a string in the format, e.g. "12/35" */
SCIP_EXPORT
SCIP_RETCODE RatCreateString(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,           /**< pointer to the rational to create */
   char*                 desc                /**< the string describing the rational */
   );

/** creates an array of rationals */
SCIP_EXPORT
SCIP_RETCODE RatCreateArray(
   SCIP_Rational***      rational,           /**< pointer to the array to create */
   int                   size                /**< the size of the array */
   );

/** creates an array of rationals */
SCIP_EXPORT
SCIP_RETCODE RatCreateBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      rational,           /**< pointer to the array to create */
   int                   size                /**< the size of the array */
   );

/** creates an array of rationals */
SCIP_EXPORT
SCIP_RETCODE RatCreateBufferArray(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational***      rational,           /**< pointer to the array to create */
   int                   size                /** the size of the array */
   );

/** copies an array of rationals */
SCIP_RETCODE RatCopyBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      target,             /**< address to copy to */
   SCIP_Rational**       src,                /**< src array */
   int                   len                 /**< size of src array */
   );

/** creates a copy of a rational */
SCIP_EXPORT
SCIP_RETCODE RatCopy(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,           /**< pointer to the rational to create */
   SCIP_Rational*        src                 /**< rational to copy */
   );

SCIP_EXPORT
SCIP_RETCODE RatCreateBuffer(
   BMS_BUFMEM*           buf,
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   );

SCIP_EXPORT
SCIP_RETCODE RatCreate(
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   );

SCIP_EXPORT
SCIP_RETCODE RatCreateBlock(
   BMS_BLKMEM*           blkmem,
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   );

#ifdef SCIP_WITH_GMP
/** creates a rational from an mpq_t */
SCIP_EXPORT
SCIP_RETCODE RatCreateGMP(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,            /**< pointer to the rational to create */
    mpq_t           numb                /**< the mpq_rational */
   );

/** gets the underlying mpq_t* */
SCIP_EXPORT mpq_t* RatGetGMP(
    SCIP_Rational*  r                   /**< the rational */
   );

SCIP_EXPORT
void RatSetGMP(
   SCIP_Rational*        r,
   const mpq_t           numb
   );

SCIP_EXPORT
void RatSetGMPArray(
   mpq_t*                res,
   SCIP_Rational**       src,
   int                   len
   );

void RatSetArrayGMP(
   SCIP_Rational**       res,
   mpq_t*                src,
   int                   len
   );

SCIP_EXPORT
void RatClearGMPArray(
   mpq_t*                ar,
   int                   len
   );
#endif

/** deletes a rational and frees the allocated memory */
SCIP_EXPORT
void RatFree(
   SCIP_Rational**       r                   /**< adress of the rational */
   );

/** deletes a rational and frees the allocated memory */
SCIP_EXPORT
void RatFreeBlock(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       r                   /**< adress of the rational */
   );

SCIP_EXPORT
void RatFreeBuffer(
   BMS_BUFMEM*           buf,
   SCIP_Rational**       r
   );

/** deletes an array of rationals and frees the allocated memory */
SCIP_EXPORT
void RatFreeArray(
   SCIP_Rational***      array,              /**< address of rational array */
   int                   size                /**< size of the array */
   );

/** deletes an array of rationals and frees the allocated memory */
SCIP_EXPORT
void RatFreeBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      array,              /**< address of rational array */
   int                   size                /**< size of the array */
   );

/** frees an array of rationals */
SCIP_EXPORT
void RatFreeBufferArray(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational***      array,              /**< pointer to the array */
   int                   size                /**< size of the array */
   );

/** sets a rational to the value of another rational */
SCIP_EXPORT
void RatSet(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*   src                 /**< the src */
   );

/** sets a rational to a nom/denom value */
SCIP_EXPORT
void RatSetInt(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Longint          nom,                /**< the nominator */
   SCIP_Longint          denom               /**< the denominator */
   );

/** sets a rational to the value described by a string */
SCIP_EXPORT
void RatSetString(
   SCIP_Rational*        res,                /**< the result */
   const char*           desc                /**< the string describing the rational */
   );

/** sets a rational to the value of another a real */
SCIP_EXPORT
void RatSetReal(
   SCIP_Rational*        r,
   SCIP_Real             real
   );

/*
 * Computing methods
 */

/* transforms rational into canonical form */
SCIP_EXPORT
void RatCanonicalize(
   SCIP_Rational*        r                   /**< rational to put in canonical form */
   );

/** adds two rationals and saves the result in res */
SCIP_EXPORT
void RatAdd(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op1,                /**< first operand */
    SCIP_Rational*  op2                 /**< second operand */
   );

/** adds a rational and a real and saves the result in res */
SCIP_EXPORT
void RatAddReal(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   );

/** subtracts two rationals and saves the result in res */
SCIP_EXPORT
void RatDiff(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*  op1,                /**< first operand */
    SCIP_Rational*  op2                 /**< second operand */
   );

/** subtracts a rational and a real and saves the result in res */
SCIP_EXPORT
void RatDiffReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   );

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) of two rationals */
SCIP_EXPORT
void RatRelDiff(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        val1,               /**< first value to be compared */
   SCIP_Rational*        val2                /**< second value to be compared */
   );

/** multiplies two rationals and saves the result in res */
SCIP_EXPORT
void RatMult(
   SCIP_Rational*        res,                /**< the result */
    SCIP_Rational*       op1,                /**< first operand */
    SCIP_Rational*       op2                 /**< second operand */
   );

/** multiply a rational and a real and save the result in res */
SCIP_EXPORT
void RatMultReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/** divide two rationals and save the result in res */
SCIP_EXPORT
void RatDiv(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   );

/** divide a rational by a real and save the result in res */
SCIP_EXPORT
void RatDivReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/* Computes res += op1 * op2 and saves the result in res */
SCIP_EXPORT
void RatAddProd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   );

/* Computes res -= op1 * op2 and saves the result in res */
SCIP_EXPORT
void RatDiffProd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   );

/** set res to -op */
SCIP_EXPORT
void RatNegate(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   );

/** set res to Abs(op) */
SCIP_EXPORT
void RatAbs(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   );

/** set res to 1/op */
SCIP_EXPORT
void RatInvert(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   );

/** compute the minimum of two rationals */
SCIP_EXPORT
void RatMIN(
   SCIP_Rational*        ret,                /**< the result */
   SCIP_Rational*        r1,                 /**< the first rational */
   SCIP_Rational*        r2                  /**< the second rational */
   );

/** compute the maximum of two rationals */
SCIP_EXPORT
void RatMAX(
   SCIP_Rational*        ret,                /**< the result */
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/*
 * Comparisoon methods
 */

/** check if two rationals are equal */
SCIP_EXPORT
SCIP_Bool RatIsEqual(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if two rationals are equal */
SCIP_EXPORT
SCIP_Bool RatIsAbsEqual(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );


/** check if a rational and a real are equal */
SCIP_EXPORT
SCIP_Bool RatIsEqualReal(
    SCIP_Rational*  r1,                 /**< the rational */
   SCIP_Real             r2                  /**< the real */
   );

/** check if real approx of rational and a real are equal */
SCIP_EXPORT
SCIP_Bool RatIsApproxEqualReal(
   SCIP_Rational*        r1,                 /**< the rational */
   SCIP_Real             r2                  /**< the real */
   );

/** check if the first rational is greater than the second*/
SCIP_EXPORT
SCIP_Bool RatIsGT(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is greater than the second*/
SCIP_EXPORT
SCIP_Bool RatIsAbsGT(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is smaller than the second*/
SCIP_EXPORT
SCIP_Bool RatIsLT(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is smaller or equal than the second*/
SCIP_EXPORT
SCIP_Bool RatIsLE(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the first rational is greater or equal than the second*/
SCIP_EXPORT
SCIP_Bool RatIsGE(
    SCIP_Rational*  r1,                 /**< the first rational */
    SCIP_Rational*  r2                  /**< the second rational */
   );

/** check if the rational is zero */
SCIP_EXPORT
SCIP_Bool RatIsZero(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is positive */
SCIP_EXPORT
SCIP_Bool RatIsPositive(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is negative */
SCIP_EXPORT
SCIP_Bool RatIsNegative(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is positive infinity */
SCIP_EXPORT
SCIP_Bool RatIsInfinity(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is negative infinity */
SCIP_EXPORT
SCIP_Bool RatIsNegInfinity(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is of infinite value */
SCIP_EXPORT
SCIP_Bool RatIsAbsInfinity(
    SCIP_Rational*  r                   /**< the rational to check */
   );

/** check if the rational is of infinite value */
SCIP_EXPORT
SCIP_Bool RatIsIntegral(
    SCIP_Rational*  r                   /**< the rational to check */
   );

SCIP_EXPORT
SCIP_Bool RatIsFpRepresentable(
    SCIP_Rational*    r
   );

/*
 * Printing/Conversion methods
 */

/** converts a rational to a string for printing, returns the number of copied characters.
 *
 *  @note If return value is equal to strlen, it means the string was truncated.
 */
SCIP_EXPORT
int RatToString(
   SCIP_Rational*        r,                  /**< the rational to print */
   char*                 str,                /**< the string to save the rational in */
   int                   strlen              /**< maximal length that can be copied to str */
   );

/** returns the strlen of a rational number */
SCIP_EXPORT
SCIP_Longint RatStrlen(
   SCIP_Rational*        r                /** rational to consider */
   );

/** prints a rational to command line (for debugging) */
SCIP_EXPORT
void RatPrint(
   SCIP_Rational*        r                   /**< the rational to print */
   );

/** printf extension for rationals (does not support all format options yet) */
SCIP_EXPORT
void RatPrintf(const char *format, ...);

/** rational extension for the SCIPdebugMsg */
#ifdef SCIP_DEBUG
#define RatDebugMessage           printf("[%s:%d] debug: ", __FILE__, __LINE__), RatPrintf
#else
#define RatDebugMessage           while( FALSE ) /*lint -e{530}*/ RatPrintf
#endif

/** prints rational to file using message handler */
SCIP_EXPORT
void RatMessage(
   SCIP_MESSAGEHDLR*     msg,                /**< message handler */
   FILE*                 file,               /**< file pointer */
   SCIP_Rational*        r                   /**< the rational to print */
   );

/** returns approximation of rational as SCIP_Real */
SCIP_EXPORT
SCIP_Real RatRoundReal(
   SCIP_Rational*        r,                  /**< the rational to convert */
   SCIP_ROUNDMODE        roundmode           /**< rounding direction (not really working yet) */
   );

/** returns approximation of rational as SCIP_Real */
SCIP_EXPORT
SCIP_Real RatApproxReal(
    SCIP_Rational*       r                   /**< the rational to convert */
   );

SCIP_EXPORT
void RatRound(
   SCIP_Rational*        retval,             /**< the resulting rounded integer */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   );

/** rounds rational to next integer in direction of roundmode */
SCIP_EXPORT
SCIP_Bool RatRoundInteger(
   SCIP_Longint*         retval,             /**< the resulting rounded lon int */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   );

/*
 * Dynamic Arrays
 */

/** creates a dynamic array of real values */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarrayCreate(
   SCIP_RATIONALARRAY**  rationalarray,      /**< pointer to store the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** resizes a dynamic array of real values */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarrayResize(
   SCIP_RATIONALARRAY*   rationalarray,      /**< pointer to store the real array */
   int                   newsize             /**< new size */
   );

/** creates a copy of a dynamic array of real values */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarrayCopy(
   SCIP_RATIONALARRAY**  rationalarray,      /**< pointer to store the copied real array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_RATIONALARRAY*   sourcerationalarray /**< dynamic real array to copy */
   );

/** frees a dynamic array of real values */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarrayFree(
   SCIP_RATIONALARRAY**  rationalarray,      /**< pointer to the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** clears a dynamic real array */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarrayClear(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic real array */
   );

/** gets value of entry in dynamic array */
SCIP_EXPORT
void SCIPrationalarrayGetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   idx,                /**< array index to get value for */
   SCIP_Rational*        result              /**< store the result */
   );

/** sets value of entry in dynamic array */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarraySetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   idx,                /**< array index to set value for */
   SCIP_Rational*        val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarrayIncVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Rational*        incval              /**< value to increase array index */
   );

/** prints a rationalarray to std out */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalArrayPrint(
   SCIP_RATIONALARRAY*   rationalarray      /**< dynamic rational array */
   );

/** returns the minimal index of all stored non-zero elements */
SCIP_EXPORT
int SCIPrationalarrayGetMinIdx(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic rational array */
   );

/** returns the maximal index of all stored non-zero elements */
SCIP_EXPORT
int SCIPrationalarrayGetMaxIdx(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic rational array */
   );

#ifdef __cplusplus
}
#endif

#endif
