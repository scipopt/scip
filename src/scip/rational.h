/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   rational.h
 * @ingroup PUBLICCOREAPI
 * @brief  wrapper for rational number arithmetic
 * @author Leon Eifler
 * @author Dominik Kamp
  */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RATIONAL_H__
#define __SCIP_RATIONAL_H__

#include <stdbool.h>
#include <stdlib.h>
#include "scip/def.h"
#include "scip/type_message.h"
#include "scip/type_rational.h"
#include "scip/type_set.h"
#include "blockmemshell/memory.h"
#ifdef SCIP_WITH_GMP
#include <gmp.h>
#endif

#ifdef SCIP_WITH_MPFR
#include <mpfr.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicRationalMethods
 *
 * @{
 */

/*
 * Creation methods
 */

/** creates a rational using standard memory allocation */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCreate(
   SCIP_RATIONAL**       rational            /**< pointer to the rational to create */
   );

/** creates a rational using buffer memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCreateBuffer(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_RATIONAL**       rational            /**< pointer to the rational to create */
   );

/** creates a rational using block memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCreateBlock(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_RATIONAL**       rational            /**< pointer to the rational to create */
   );

/** creates a copy of a rational using ordinary memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCopy(
   SCIP_RATIONAL**       result,             /**< pointer to the rational to create */
   SCIP_RATIONAL*        src                 /**< rational to copy */
   );

/** creates a copy of a rational using block memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCopyBlock(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_RATIONAL**       result,             /**< pointer to the rational to create */
   SCIP_RATIONAL*        src                 /**< rational to copy */
   );

/** creates a copy of a rational */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCopyBuffer(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_RATIONAL**       result,             /**< pointer to the rational to create */
   SCIP_RATIONAL*        src                 /**< rational to copy */
   );

/** creates an array of rationals using ordinary memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCreateArray(
   SCIP_RATIONAL***      rational,           /**< pointer to the array to create */
   int                   size                /**< the size of the array */
   );

/** creates an array of rationals using block memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCreateBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_RATIONAL***      rational,           /**< pointer to the array to create */
   int                   size                /**< the size of the array */
   );

/** creates an array of rationals using buffer memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCreateBufferArray(
   BMS_BUFMEM*           mem,                /**< buffer memory */
   SCIP_RATIONAL***      rational,           /**< pointer to the array to create */
   int                   size                /**< the size of the array */
   );

/** copies an array of rationals using ordinary memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCopyArray(
   SCIP_RATIONAL***      target,             /**< address to copy to */
   SCIP_RATIONAL**       src,                /**< src array */
   int                   len                 /**< size of src array */
   );

/** copies an array of rationals using block memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCopyBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_RATIONAL***      target,             /**< address to copy to */
   SCIP_RATIONAL**       src,                /**< src array */
   int                   len                 /**< size of src array */
   );

/** copy an array of rationals using buffer memory */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCopyBufferArray(
   BMS_BUFMEM*           mem,                /**< buffer memory */
   SCIP_RATIONAL***      result,             /**< address to copy to */
   SCIP_RATIONAL**       src,                /**< src array */
   int                   len                 /**< size of src array */
   );

/** realloc a rational ordinary array */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalReallocArray(
   SCIP_RATIONAL***      result,             /**< address to copy to */
   int                   oldlen,             /**< size of src array */
   int                   newlen              /**< size of src array */
   );

/** realloc a rational buffer arrray */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalReallocBufferArray(
   BMS_BUFMEM*           mem,                /**< buffer memory */
   SCIP_RATIONAL***      result,             /**< address to copy to */
   int                   oldlen,             /**< size of src array */
   int                   newlen              /**< size of src array */
   );

/** realloc a rational block arrray */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalReallocBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_RATIONAL***      result,             /**< address to copy to */
   int                   oldlen,             /**< size of src array */
   int                   newlen              /**< size of src array */
   );

#if defined(SCIP_WITH_BOOST) && defined(SCIP_WITH_GMP)
/** gets the underlying gmp rational pointer */
SCIP_EXPORT
mpq_t* SCIPrationalGetGMP(
   SCIP_RATIONAL*        rational            /**< rational to access */
   );

/** sets rational to gmp rational */
SCIP_EXPORT
void SCIPrationalSetGMP(
   SCIP_RATIONAL*        rational,           /**< rational to define */
   const mpq_t           numb                /**< gmp rational to set */
   );

/** creates rational from gmp rational */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCreateBlockGMP(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_RATIONAL**       rational,           /**< pointer to the rational to create */
   mpq_t                 numb                /**< gmp rational to set */
   );

/** sets gmp rational array to values of rational array */
SCIP_EXPORT
void SCIPrationalSetGMPArray(
   mpq_t*                mpqaaray,           /**< gmp rational array */
   SCIP_RATIONAL**       ratarrray,          /**< rational array */
   int                   len                 /**< array length */
   );

/** sets rational array to values of gmp rational array */
SCIP_EXPORT
void SCIPrationalSetArrayGMP(
   SCIP_RATIONAL**       ratarray,           /**< rational array */
   mpq_t*                mpqarray,           /**< gmp rational array */
   int                   len                 /**< array length */
   );

/** clears gmp rational array */
SCIP_EXPORT
void SCIPrationalClearArrayGMP(
   mpq_t*                mpqarray,           /**< gmp rational array */
   int                   len                 /**< array length */
   );
#endif

/** deletes a rational and frees the allocated ordinary memory */
SCIP_EXPORT
void SCIPrationalFree(
   SCIP_RATIONAL**       rational            /**< address of the rational */
   );

/** deletes a rational and frees the allocated block memory */
SCIP_EXPORT
void SCIPrationalFreeBlock(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_RATIONAL**       rational            /**< address of the rational */
   );

/** deletes a rational and frees the allocated buffer memory */
SCIP_EXPORT
void SCIPrationalFreeBuffer(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_RATIONAL**       rational            /**< address of the rational */
   );

/** deletes an array of rationals and frees the allocated ordinary memory */
SCIP_EXPORT
void SCIPrationalFreeArray(
   SCIP_RATIONAL***      ratarray,           /**< address of rational array */
   int                   size                /**< size of the array */
   );

/** deletes an array of rationals and frees the allocated block memory */
SCIP_EXPORT
void SCIPrationalFreeBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_RATIONAL***      ratblockarray,      /**< address of rational array */
   int                   size                /**< size of the array */
   );

/** deletes an array of rationals and frees the allocated buffer memory */
SCIP_EXPORT
void SCIPrationalFreeBufferArray(
   BMS_BUFMEM*           mem,                /**< buffer memory */
   SCIP_RATIONAL***      ratbufarray,        /**< pointer to the array */
   int                   size                /**< size of the array */
   );

/** sets a rational to the value of another rational */
SCIP_EXPORT
void SCIPrationalSetRational(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        src                 /**< the src */
   );

/** sets a rational to a nom/denom value */
SCIP_EXPORT
void SCIPrationalSetFraction(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_Longint          nom,                /**< the nominator */
   SCIP_Longint          denom               /**< the denominator */
   );

/** sets a rational to the value of another a real */
SCIP_EXPORT
void SCIPrationalSetReal(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_Real             real                /**< the value to set to */
   );

/** sets a rational to positive infinity */
SCIP_EXPORT
void SCIPrationalSetInfinity(
   SCIP_RATIONAL*        res                 /**< the result */
   );

/** sets a rational to negative infinity */
SCIP_EXPORT
void SCIPrationalSetNegInfinity(
   SCIP_RATIONAL*        res                 /**< the result */
   );

/** checks if a string describes a rational number */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsString(
   const char*           desc                /**< string to check */
   );

/** sets a rational to the value described by a string */
SCIP_EXPORT
void SCIPrationalSetString(
   SCIP_RATIONAL*        res,                /**< the result */
   const char*           desc                /**< the string describing the rational */
   );

/** allocates and creates a rational from a string if known, otherwise assigns a null pointer */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCreateString(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_RATIONAL**       rational,           /**< pointer to the rational to create */
   const char*           desc                /**< the string describing the rational */
   );

/** extract the next token as a rational value if it is one; in case no value is parsed the endptr is set to @p desc
 *
 *  @return Returns TRUE if a value could be extracted, otherwise FALSE
 */
SCIP_EXPORT
SCIP_Bool SCIPstrToRationalValue(
   char*                 desc,               /**< string to search */
   SCIP_RATIONAL*        value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p desc */
   );

/** resets the flag isfprepresentable to SCIP_ISFPREPRESENTABLE_UNKNOWN */
SCIP_EXPORT
void SCIPrationalResetFloatingPointRepresentable(
   SCIP_RATIONAL*        rat                 /**< the number to set flag for */
   );

/*
 * Computing methods
 */

/* transforms rational into canonical form */
SCIP_EXPORT
void SCIPrationalCanonicalize(
   SCIP_RATIONAL*        rational            /**< rational to put in canonical form */
   );

/* checks if the underlying Rational has a value >= infinity;
 * needed after underlying value was directly set, e.g. by exact lp solver
 */
SCIP_EXPORT
void SCIPrationalCheckInfByValue(
   SCIP_RATIONAL*        rational            /**< rational number */
   );

/** adds two rationals and saves the result in res */
SCIP_EXPORT
void SCIPrationalAdd(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_RATIONAL*        op2                 /**< second operand */
   );

/** adds a rational and a real and saves the result in res */
SCIP_EXPORT
void SCIPrationalAddReal(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   );

/** subtracts two rationals and saves the result in res */
SCIP_EXPORT
void SCIPrationalDiff(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_RATIONAL*        op2                 /**< second operand */
   );

/** subtracts a rational and a real and saves the result in res */
SCIP_EXPORT
void SCIPrationalDiffReal(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   );

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) of two rationals
 *
 *  @note this method handles infinity like finite numbers
 */
SCIP_EXPORT
void SCIPrationalRelDiff(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        val1,               /**< first value to be compared */
   SCIP_RATIONAL*        val2                /**< second value to be compared */
   );

/** multiplies two rationals and saves the result in res */
SCIP_EXPORT
void SCIPrationalMult(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_RATIONAL*        op2                 /**< second operand */
   );

/** multiply a rational and a real and save the result in res */
SCIP_EXPORT
void SCIPrationalMultReal(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/** divide two rationals and save the result in res */
SCIP_EXPORT
void SCIPrationalDiv(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_RATIONAL*        op2                 /**< second operand */
   );

/** divide a rational by a real and save the result in res */
SCIP_EXPORT
void SCIPrationalDivReal(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/* Computes res += op1 * op2 and saves the result in res */
SCIP_EXPORT
void SCIPrationalAddProd(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_RATIONAL*        op2                 /**< second operand */
   );

/* Computes res += op1 * op2 and saves the result in res */
SCIP_EXPORT
void SCIPrationalAddProdReal(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/* Computes res -= op1 * op2 and saves the result in res */
SCIP_EXPORT
void SCIPrationalDiffProd(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_RATIONAL*        op2                 /**< second operand */
   );

/* Computes res -= op1 * op2 and saves the result in res */
SCIP_EXPORT
void SCIPrationalDiffProdReal(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   );

/** set res to -op */
SCIP_EXPORT
void SCIPrationalNegate(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op                  /**< operand */
   );

/** set res to Abs(op) */
SCIP_EXPORT
void SCIPrationalAbs(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op                  /**< operand */
   );

/** set res to 1/op */
SCIP_EXPORT
void SCIPrationalInvert(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op                  /**< operand */
   );

/** compute the minimum of two rationals */
SCIP_EXPORT
void SCIPrationalMin(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< the first rational */
   SCIP_RATIONAL*        op2                 /**< the second rational */
   );

/** compute the maximum of two rationals */
SCIP_EXPORT
void SCIPrationalMax(
   SCIP_RATIONAL*        res,                /**< the result */
   SCIP_RATIONAL*        op1,                /**< the first rational */
   SCIP_RATIONAL*        op2                 /**< the second rational */
   );

/*
 * Comparison methods
 */

/** checks if two rationals are equal */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsEQ(
   SCIP_RATIONAL*        rat1,               /**< the first rational */
   SCIP_RATIONAL*        rat2                /**< the second rational */
   );

/** checks if two rationals are equal */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsAbsEQ(
   SCIP_RATIONAL*        rat1,               /**< the first rational */
   SCIP_RATIONAL*        rat2                /**< the second rational */
   );

/** checks if rational and real are equal */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsEQReal(
   SCIP_RATIONAL*        rat,                /**< the rational */
   SCIP_Real             real                /**< the real */
   );

/** checks if real approx of rational and real are equal */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsApproxEQReal(
   SCIP_SET*             set,                /**< SCIP set pointer */
   SCIP_RATIONAL*        rat,                /**< the rational */
   SCIP_Real             real,               /**< the real */
   SCIP_ROUNDMODE_RAT    roundmode           /**< the rounding mode to use */
   );

/** checks if first rational is greater than second rational */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsGT(
   SCIP_RATIONAL*        rat1,               /**< the first rational */
   SCIP_RATIONAL*        rat2                /**< the second rational */
   );

/** checks if first rational is smaller than second rational */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsLT(
   SCIP_RATIONAL*        rat1,               /**< the first rational */
   SCIP_RATIONAL*        rat2                /**< the second rational */
   );

/** checks if first rational is greater or equal than second rational */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsGE(
   SCIP_RATIONAL*        rat1,               /**< the first rational */
   SCIP_RATIONAL*        rat2                /**< the second rational */
   );

/** checks if first rational is less or equal than second rational */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsLE(
   SCIP_RATIONAL*        rat1,               /**< the first rational */
   SCIP_RATIONAL*        rat2                /**< the second rational */
   );

/** checks if first rational is greater than second rational */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsAbsGT(
   SCIP_RATIONAL*        rat1,               /**< the first rational */
   SCIP_RATIONAL*        rat2                /**< the second rational */
   );

/** checks if rational is greater than real */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsGTReal(
   SCIP_RATIONAL*        rat,                /**< the rational */
   SCIP_Real             real                /**< the real */
   );

/** checks if rational is less than real */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsLTReal(
   SCIP_RATIONAL*        rat,                /**< the rational */
   SCIP_Real             real                /**< the real */
   );

/** checks if rational is greater or equal than real */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsGEReal(
   SCIP_RATIONAL*        rat,                /**< the rational */
   SCIP_Real             real                /**< the real */
   );

/** checks if rational is less or equal than real */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsLEReal(
   SCIP_RATIONAL*        rat,                /**< the rational */
   SCIP_Real             real                /**< the real */
   );

/** checks if rational is zero */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsZero(
   SCIP_RATIONAL*        rational            /**< the rational to check */
   );

/** checks if rational is positive */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsPositive(
   SCIP_RATIONAL*        rational            /**< the rational to check */
   );

/** checks if rational is negative */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsNegative(
   SCIP_RATIONAL*        rational            /**< the rational to check */
   );

/** checks if rational is positive infinity */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsInfinity(
   SCIP_RATIONAL*        rational            /**< the rational to check */
   );

/** checks if rational is negative infinity */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsNegInfinity(
   SCIP_RATIONAL*        rational            /**< the rational to check */
   );

/** checks if rational is negative infinity */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsAbsInfinity(
   SCIP_RATIONAL*        rational            /**< the rational to check */
   );

/** checks if rational is integral */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsIntegral(
   SCIP_RATIONAL*        rational            /**< the rational to check */
   );

/** checks if rational is exactly representable as real */
SCIP_EXPORT
SCIP_Bool SCIPrationalIsFpRepresentable(
   SCIP_RATIONAL*        rational            /**< the rational to check */
   );

/*
 * Printing/Conversion methods
 */

/** converts a rational to a string for printing, returns the number of copied characters.
 *
 *  @return number of characters printed into string, see also SCIPstrncpy()
 *
 *  @note If return value is equal to strlen, it means the string was truncated.
 */
SCIP_EXPORT
int SCIPrationalToString(
   SCIP_RATIONAL*        rational,           /**< the rational to print */
   char*                 str,                /**< the string to save the rational in */
   int                   strlen              /**< maximal length that can be copied to str */
   );

/** returns the strlen of a rational number */
SCIP_EXPORT
int SCIPrationalStrLen(
   SCIP_RATIONAL*        rational           /** rational to consider */
   );

/* if we have a C99 compiler */
#ifdef SCIP_HAVE_VARIADIC_MACROS

/** rational extension for the SCIPdebugMsg */
/*lint -emacro(681,SCIPrationalDebugMessage) */
/*lint -emacro(506,SCIPrationalDebugMessage) */
/*lint -emacro(774,SCIPrationalDebugMessage) */
#ifdef SCIP_DEBUG
#define SCIPrationalDebugMessage(...)            SCIPrationalPrintDebugMessage(__FILE__, __LINE__, __VA_ARGS__)
#else
#define SCIPrationalDebugMessage(...)            while ( FALSE ) SCIPrationalPrintDebugMessage(__FILE__, __LINE__, __VA_ARGS__)
#endif

#else
/* if we do not have a C99 compiler, use a workaround that prints a message, but not the file and linenumber */

/** rational extension for the SCIPdebugMsg */
#ifdef SCIP_DEBUG
#define SCIPrationalDebugMessage                 printf("debug: "), SCIPrationalPrintf
#else
#define SCIPrationalDebugMessage                 while ( FALSE ) SCIPrationalPrintf
#endif

#endif

/** prints rational into a file using message handler */
SCIP_EXPORT
void SCIPrationalMessage(
   SCIP_MESSAGEHDLR*     msg,                /**< message handler */
   FILE*                 file,               /**< file pointer */
   SCIP_RATIONAL*        rational            /**< the rational to print */
   );

/** prints rational depending on the verbosity level */
SCIP_EXPORT
void SCIPrationalPrintVerbInfo(
   SCIP_MESSAGEHDLR*     msg,                /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   SCIP_RATIONAL*        rational            /**< the rational to print */
   );

/** prints a rational to command line (for debugging) */
SCIP_EXPORT
void SCIPrationalPrint(
   SCIP_RATIONAL*        rational            /**< the rational to print */
   );

/** printf extension for rationals (does not support all format options yet) */
SCIP_EXPORT
void SCIPrationalPrintf(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a debug message */
void SCIPrationalPrintDebugMessage(
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline,         /**< line in the source file where the function was called */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** returns the numerator of a rational as a long */
SCIP_EXPORT
SCIP_Longint SCIPrationalNumerator(
   SCIP_RATIONAL*        rational            /**< the rational */
   );

/** returns the denominator of a rational as a long */
SCIP_EXPORT
SCIP_Longint SCIPrationalDenominator(
   SCIP_RATIONAL*        rational            /**< the rational */
   );

/** compares denominator of a rational to a long */
SCIP_EXPORT
SCIP_Bool SCIPrationalDenominatorIsLE(
   SCIP_RATIONAL*        rational,           /**< the rational */
   SCIP_Longint          val                 /**< long value to compare to */
   );

/** returns the sign of the rational (1 if positive, -1 if negative, 0 if zero) */
SCIP_EXPORT
int SCIPrationalGetSign(
   const SCIP_RATIONAL*  rational            /**< the rational */
   );

/** computes fractional part of a rational number */
SCIP_EXPORT
void SCIPrationalGetFrac(
   SCIP_RATIONAL*        res,                /**< rational to save the frac */
   SCIP_RATIONAL*        src                 /**< src rational */
   );

/** returns approximation of rational as SCIP_Real */
SCIP_EXPORT
SCIP_Real SCIPrationalGetReal(
   SCIP_RATIONAL*        rational            /**< the rational to convert */
   );

/** gets the relaxation of a rational as a real
 *
 *  @note Requires MPFR if rational is not fp-representable and roundmode is different from SCIP_R_ROUND_NEAREST.
 */
SCIP_EXPORT
SCIP_Real SCIPrationalRoundReal(
   SCIP_RATIONAL*        rational,           /**< the rational to convert */
   SCIP_ROUNDMODE_RAT    roundmode           /**< the rounding direction */
   );

/** rounds a rational to an integer and saves it as a rational */
SCIP_EXPORT
void SCIPrationalRoundInteger(
   SCIP_RATIONAL*        res,                /**< the resulting rounded integer */
   SCIP_RATIONAL*        src,                /**< the rational to round */
   SCIP_ROUNDMODE_RAT    roundmode           /**< the rounding direction */
   );

/** rounds rational to next integer in direction of roundmode
 *
 *  @return FALSE if rational outside of long-range
 */
SCIP_EXPORT
SCIP_Bool SCIPrationalRoundLong(
   SCIP_Longint*         res,                /**< the resulting rounded long int */
   SCIP_RATIONAL*        src,                /**< the rational to round */
   SCIP_ROUNDMODE_RAT    roundmode           /**< the rounding direction */
   );

/** compute an approximate number with denominator <= maxdenom, closest to src and save it in res using continued fractions */
SCIP_EXPORT
void SCIPrationalComputeApproximation(
   SCIP_RATIONAL*        res,
   SCIP_RATIONAL*        src,
   SCIP_Longint          maxdenom,
   int                   forcegreater        /**< 1 if res >= src should be enforced, -1 if res <= src should be enforced, 0 else */
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

/** gets value of entry in dynamic array */
SCIP_EXPORT
void SCIPrationalarrayGetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   idx,                /**< array index to get value for */
   SCIP_RATIONAL*        result              /**< store the result */
   );

/** sets value of entry in dynamic array */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarraySetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   idx,                /**< array index to set value for */
   SCIP_RATIONAL*        val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarrayIncVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   idx,                /**< array index to increase value for */
   SCIP_RATIONAL*        incval              /**< value to increase array index */
   );

/** prints a rationalarray to std out */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalarrayPrint(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic rational array */
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

/** changes the infinity threshold to new value */
SCIP_EXPORT
void SCIPrationalChgInfinity(
   SCIP_Real             inf                 /**< new infinity value */
   );

/** return the infinity threshold for rationals */
SCIP_EXPORT
SCIP_Real SCIPrationalGetInfinity(
   void
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
