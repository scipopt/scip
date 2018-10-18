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

/**@file    rational.cpp
 * @brief
 * @ingroup
 * @author  Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cmath>
#include <vector>
#include "blockmemshell/memory.h"
#include "scip/rational.h"
#include "scip/multiprecision.hpp"
#include "scip/type_message.h"
#include "scip/type_retcode.h"
#include "scip/pub_message.h"
#include <iostream>

extern "C"{

using std::vector;

struct SCIP_Rational{
   Rational r;
};


/** Allocate and create a rational from nominator and denominator */
EXTERN
SCIP_Rational* RCreateInt(
   int                   nom,                /**< The nominator */
   int                   denom               /**< The denominator */
   )
{
   SCIP_Rational* retrat;

   BMSallocMemory(& retrat);
   new (&retrat->r) Rational(nom, denom);

   return retrat;
}

/** Allocate and create a rational from a string in the format, e.g. "12/35" */
EXTERN
SCIP_Rational* RCreateString(
   const char*           desc                /**< The String describing the rational */
   )
{
   SCIP_Rational* retrat;

   BMSallocMemory(& retrat);
   new (&retrat->r) Rational(desc);

   return retrat;
}

EXTERN
SCIP_Rational** RCreateArray(
   int                   size
   )
{
   SCIP_Rational** retrat;

   BMSallocMemoryArray(&retrat, size);

   for( int i = 0; i < size; ++i )
   {
      BMSallocMemory(& retrat[i]);
      new (&retrat[i]->r) Rational();
   }

   return retrat;
}

EXTERN
void RDeleteArray(
   SCIP_Rational***      array,
   int                   size
   )
{
   assert(array != NULL);

   for( int i = 0; i < size; ++i )
   {
      RDelete(array[i]);
   }

   BMSfreeMemoryArray(array);
}

/** Delete a rational and free the allocated memory */
void RDelete(
   SCIP_Rational**       r                   /**< Adress of the rational */
   )
{
   if( r == NULL)
      return;

   (*r)->r.~Rational();
   BMSfreeMemory(r);
}

/** Set a rational to the value of another rational */
EXTERN
void RSet(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        src                 /**< The src */
   )
{
   res->r = src->r;
}

/** Set a rational to a nom/denom value */
EXTERN
void RSetInt(
   SCIP_Rational*        res,                /**< The result */
   int                   nom,                /**< The nominator */
   int                   denom               /**< The denominator */
   )
{
   res->r = (nom/denom);
}

/** Add two rationals and save the result in res*/
EXTERN
void RAdd(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op1,                /**< First operand */
   SCIP_Rational*        op2                 /**< Second operand */
   )
{
   res->r = op1->r + op2->r;
}

/** Subtract two rationals and save the result in res*/
EXTERN
void RSub(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op1,                /**< First operand */
   SCIP_Rational*        op2                 /**< Second operand */
   )
{
   res->r = op1->r - op2->r;
}

/** Multiply two rationals and save the result in res*/
EXTERN
void RMult(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op1,                /**< First operand */
   SCIP_Rational*        op2                 /**< Second operand */
   )
{
   res->r = op1->r * op2->r;
}

/** Divide two rationals and save the result in res*/
EXTERN
void RDiv(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op1,                /**< First operand */
   SCIP_Rational*        op2                 /**< Second operand */
   )
{
   res->r = op1->r / op2->r;
}

/** Set res to -op */
EXTERN
void RNeg(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op                  /**< Operand */
   )
{
   res->r = -op->r;
}


/** Set res to Abs(op) */
EXTERN
void RAbs(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op                  /**< Operand */
   )
{
   res->r = abs(op->r);
}


/** Set res to 1/op */
EXTERN
void RInv(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        op                  /**< Operand */
   )
{
   res->r = 1 / op->r;
}


/** Print a Rational to std out */
void RPrint(
   SCIP_Rational*        r                   /**< The rational to print */
   )
{
   std::cout << r->r << std::endl ;
}




}