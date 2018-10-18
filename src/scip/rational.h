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


#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Rational SCIP_Rational;

/** Allocate and create a rational from nominator and denominator */
EXTERN
SCIP_Rational* RCreateInt(
   int                   nom,                /**< The nominator */
   int                   denom               /**< The denominator */
   );

/** Allocate and create a rational from a string in the format, e.g. "12/35" */
EXTERN
SCIP_Rational* RCreateString(
   const char*           desc                /**< The String describing the rational */
   );

EXTERN
SCIP_Rational** RCreateArray(
   int                   size
   );

/** Delete a rational and free the allocated memory */
EXTERN
void RDelete(
   SCIP_Rational**       r                   /**< Adress of the rational */
   );

EXTERN
void RDeleteArray(
   SCIP_Rational***      array,
   int                   size
   );

/** Set a rational to the value of another rational */
EXTERN
void RSet(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        src                 /**< The src */
   );

/** Set a rational to a nom/denom value */
EXTERN
void RSetInt(
   SCIP_Rational*        res,                /**< The result */
   int                   nom,                /**< The nominator */
   int                   denom               /**< The denominator */
   );

/** Add two rationals and save the result in res*/
EXTERN
void RAdd(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op1,                /**< First operand */
   SCIP_Rational*        op2                 /**< Second operand */
   );

/** Subtract two rationals and save the result in res*/
EXTERN
void RSub(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op1,                /**< First operand */
   SCIP_Rational*        op2                 /**< Second operand */
   );

/** Multiply two rationals and save the result in res*/
EXTERN
void RMult(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op1,                /**< First operand */
   SCIP_Rational*        op2                 /**< Second operand */
   );

/** Divide two rationals and save the result in res*/
EXTERN
void RDiv(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op1,                /**< First operand */
   SCIP_Rational*        op2                 /**< Second operand */
   );

/** Set res to -op */
EXTERN
void RNeg(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op                  /**< Operand */
   );

/** Set res to Abs(op) */
EXTERN
void RAbs(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op                  /**< Operand */
   );

/** Set res to 1/op */
EXTERN
void RInv(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op                  /**< Operand */
   );

/** Print a Rational to std out */
void RPrint(
   SCIP_Rational*        r                   /**< The rational to print */
   );


#ifdef __cplusplus
}
#endif

#endif
