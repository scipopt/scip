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
#include "scip/pub_misc.h"
#include "scip/intervalarith.h"
#include <iostream>

#ifdef WITH_GMP
#include <gmp.h>
#endif

extern "C"{

using std::vector;

struct SCIP_Rational{
   Rational r;
   SCIP_Bool isinf = FALSE;
};

/*
 * Creation methods
 */

/** allocate and create a rational from nominator and denominator */
SCIP_Rational* RcreateInt(
   BMS_BLKMEM*           mem,                /**< block memory */
   int                   nom,                /**< the nominator */
   int                   denom               /**< the denominator */
   )
{
   SCIP_Rational* retrat;

   BMSallocBlockMemory(mem, &retrat);

   retrat->isinf = FALSE;
   new (&retrat->r) Rational(nom, denom);

   return retrat;
}

/** allocate and create a rational from a string in the format, e.g. "12/35" */
SCIP_Rational* RcreateString(
   BMS_BLKMEM*           mem,                /**< block memory */
   const char*           desc                /**< the String describing the rational */
   )
{
   SCIP_Rational* retrat;

   BMSallocBlockMemory(mem, &retrat);
   if( 0 == strcmp(desc, "inf") )
   {
      new (&retrat->r) Rational(1);
      retrat->isinf = TRUE;
   }
   else if ( 0 == strcmp(desc, "-inf") )
   {
      new (&retrat->r) Rational(-1);
      retrat->isinf = TRUE;
   }
   else
   {
      new (&retrat->r) Rational(desc);
   }

   return retrat;
}

/** allocate and create a rational from a string in the format, e.g. "12/35" */
static
SCIP_Rational* RcreateReal(
   BMS_BLKMEM*           mem,                /**< block memory */
   const SCIP_Real       dbl                 /**< the string describing the rational */
   )
{
   SCIP_Rational* retrat;

   BMSallocBlockMemory(mem, &retrat);
   retrat->isinf = FALSE;
   new (&retrat->r) Rational(dbl);

   return retrat;
}

/** create an array of rationals */
SCIP_Rational** RcreateArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   int                   size                /** the size of the array */
   )
{
   SCIP_Rational** retrat;

   BMSallocBlockMemoryArray(mem, &retrat, size);

   for( int i = 0; i < size; ++i )
   {
      BMSallocBlockMemory(mem, &retrat[i]);
      new (&retrat[i]->r) Rational();
   }

   return retrat;
}

/* create a copy of a rational */
SCIP_Rational* Rcopy(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational*        src                 /**< rational to copy */
   )
{
   SCIP_Rational* ret;
   ret = RcreateInt(mem, 0,1);

   Rset(ret, src);
   return ret;
}


/** cretae a rational from an mpq_t */
#ifdef WITH_GMP
SCIP_Rational* RcreateGMP(
   BMS_BLKMEM*           mem,                /**< block memory */
   const mpq_t           numb                /**< the gmp rational */
   )
{
   SCIP_Rational* retrat;

   BMSallocBlockMemory(mem, &retrat);
   retrat->isinf = FALSE;
   new (&retrat->r) Rational(numb);

   return retrat;
}

/** get the underlying mpq_t* */
mpq_t* RgetMpq(
   SCIP_Rational*        r                   /**< the rational */
   )
{
   return &(r->r.backend().data());
}

#endif

/** free an array of rationals */
void RdeleteArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      array,              /**< pointer to the array */
   int                   size                /**< size of the array */
   )
{
   assert(array != NULL);

   for( int i = 0; i < size; ++i )
   {
      Rdelete(mem, array[i]);
   }

   BMSfreeBlockMemoryArray(mem, array, size);
}

/** delete a rational and free the allocated memory */
void Rdelete(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       r                   /**< adress of the rational */
   )
{
   if( r == NULL)
      return;

   (*r)->r.~Rational();
   BMSfreeBlockMemory(mem, r);
}

/** set a rational to the value of another rational */
void Rset(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        src                 /**< the src */
   )
{
   assert(res != NULL);

   res->r = src->r;
}

/** set a rational to a nom/denom value */
void RsetInt(
   SCIP_Rational*        res,                /**< the result */
   int                   nom,                /**< the nominator */
   int                   denom               /**< the denominator */
   )
{
   assert(res != NULL);

   res->r = (nom/denom);
}

/*
 * Computing methods
 */

/** add two rationals and save the result in res*/
void Radd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   res->r = op1->r + op2->r;
}

/** add a rational and a real and save the result in res*/
void RaddReal(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   )
{
   assert(res != NULL && rat != NULL);

   res->r = rat->r + real;
}

/** subtract two rationals and save the result in res*/
void Rdiff(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   res->r = op1->r - op2->r;
}

/** subtract a rational and a real and save the result in res*/
void RdiffReal(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   )
{
   assert(res != NULL && rat != NULL);

   res->r = rat->r - real;
}

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
void RrelDiff(
   SCIP_Rational*        res,
   SCIP_Rational*        val1,               /**< first value to be compared */
   SCIP_Rational*        val2                /**< second value to be compared */
   )
{
   Rational absval1;
   Rational absval2;
   Rational quot;

   assert(res != NULL && val1 != NULL && val2 != NULL);
   assert(!val1->isinf && !val2->isinf);

   absval1 = abs(val1->r);
   absval2 = abs(val2->r);
   quot = max(absval1, absval2);
   if( 1.0 > quot );
      quot = 1.0;

   res->r = ((val1->r)-(val2->r))/quot;
}

/** multiply two rationals and save the result in res*/
void Rmult(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   res->r = op1->r * op2->r;
}

/** multiply a rational and a real and save the result in res*/
void RmultReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   )
{
    assert(res != NULL && op1 != NULL);

    if( op2->isinf )
    {
       res->r = op2 > 0 ? op1->r : -op1->r;
       res->isinf = true;
    }
    else
    {
       res->r = op1->r * op2;
    }
}


/** divide two rationals and save the result in res*/
void Rdiv(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   res->r = op1->r / op2->r;
}

/** divide a rational and a real and save the result in res*/
void RdivReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL);

   assert(op2 != 0.0);

   RmultReal(res, op1, 1.0 / op2 );
}

/** set res to -op */
void Rneg(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   )
{
   assert(res != NULL && op != NULL);

   res->r = -op->r;
}


/** set res to Abs(op) */
void Rabs(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   )
{
   assert(res != NULL && op != NULL);

   res->r = abs(op->r);
}


/** set res to 1/op */
void Rinv(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   )
{
   assert(res != NULL && op != NULL);
   assert(!res->isinf);

   res->r = 1 / op->r;
}

/*
 * Comparisoon methods
 */

/** compute the minimum of two rationals */
void Rmin(
   SCIP_Rational*        ret,                /**< the result */
   SCIP_Rational*        r1,                 /**< the first rational */
   SCIP_Rational*        r2                  /**< the second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);
   assert(r1->r != NULL && r2->r != NULL);

   if( r2->isinf )
      ret->r = r2->r > 0 ? r1->r : r2->r;
   else if( r1->isinf )
      ret->r = r1->r > 0 ? r2->r : r1->r;
   else
      ret->r = min(r1->r, r2->r);
}

/** check if two rationals are equal */
SCIP_Bool RisEqual(
   SCIP_Rational*        r1,                 /**< the first rational */
   SCIP_Rational*        r2                  /**< the second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);

   if( r1->isinf )
   {
      if( r2->isinf && (r1->r) == (r2->r) )
         return TRUE;
      else
         return FALSE;
   }
   else if( r2->isinf )
      return FALSE;
   else
      return (r1->r) == (r2->r);
}

/** check if a rational and a real are equal */
SCIP_Bool RisEqualReal(
   SCIP_Rational*        r1,                 /**< the rational */
   SCIP_Real             r2                  /**< the real */
   )
{
   assert(r1 != NULL);

   if( r1->isinf )
      return FALSE;
   else
      return r1->r == r2;
}

/** check if the first rational is greater than the second*/
SCIP_Bool RisGT(
   SCIP_Rational*        r1,                 /**< The first rational */
   SCIP_Rational*        r2                  /**< The second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);

   if( r1->isinf )
   {
      if( r1->r < 0 )
         return FALSE;
      else if( r2->isinf && (r2->r) > 0 )
         return FALSE;
      else
         return TRUE;
   }
   else if( r2->isinf )
   {
      if( r2 > 0 )
         return FALSE;
      else
         return TRUE;
   }
   else
   {
      return r1->r > r2->r;
   }
}

/** check if the first rational is smaller than the second*/
SCIP_Bool RisLT(
   SCIP_Rational*        r1,                 /**< The first rational */
   SCIP_Rational*        r2                  /**< The second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);

   return RisGT(r2, r1);
}

/** check if the first rational is greater or equal than the second*/
EXTERN
SCIP_Bool RisGE(
   SCIP_Rational*        r1,                 /**< The first rational */
   SCIP_Rational*        r2                  /**< The second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);

   if( r1->isinf )
   {
      if( r1->r > 0 )
         return TRUE;
      else if( r2->isinf && (r2->r) < 0 )
         return TRUE;
      else
         return FALSE;
   }
   else if( r2->isinf )
   {
      if( r2 < 0 )
         return TRUE;
      else
         return FALSE;
   }
   else
   {
      return r1->r >= r2->r;
   }
}

/** check if the first rational is less or equal than the second*/
EXTERN
SCIP_Bool RisLE(
   SCIP_Rational*        r1,                 /**< The first rational */
   SCIP_Rational*        r2                  /**< The second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);

   return RisGE(r2, r1);
}


/** check if the rational is zero */
SCIP_Bool RisZero(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);
   assert(r->r != NULL);

   return r->r.is_zero();
}

/** check if the rational is positive */
SCIP_Bool RisPositive(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return r->r > 0;
}

/** check if the rational is negative */
SCIP_Bool RisNegative(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return r->r < 0;
}

/** check if the rational is positive infinity */
SCIP_Bool RisInfinity(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return r->isinf && r > 0;
}

/** check if the rational is negative infinity */
SCIP_Bool RisNegInfinity(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return r->isinf && r < 0;
}

/** check if the rational is negative infinity */
SCIP_Bool RisAbsInfinity(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return r->isinf;
}

/*
 * Printing/Conversion methods 
 */

/** convert a Rational to a string for printing */
void RtoString(
   SCIP_Rational*        r,                  /**< the rational to print */
   char*                 str                 /**< the string to save the rational in */
   )
{
   assert(r != NULL);

   if( r->isinf )
   {
      (void) SCIPstrncpy(str, "infinity value", SCIP_MAXSTRLEN);
   }
   else
   {
      std::string s = (r->r).str();
      (void) SCIPstrncpy(str, s.c_str(), SCIP_MAXSTRLEN);
   }
}

/** print a rational to command line (for debugging) */
void Rprint(
   SCIP_Rational*        r                   /**< the rational to print */
   )
{
   assert(r != NULL);
   if( r->isinf )
      std::cout << r->r << "inf" << std::endl;
   else
      std::cout << r->r << std::endl;
}

/** get the relaxation of a rational as a real, unfortunately you can't control the roundmode without using mpfr */
SCIP_Real RgetRealRelax(
   SCIP_Rational*        r,                  /**< the rational */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   )
{
   SCIP_Real realapprox;
   SCIP_Real nom, denom;
   SCIP_ROUNDMODE current;

   assert(r != NULL);

   if( r->isinf )
      return (r->r.sign() * SCIP_DEFAULT_INFINITY);

   current = SCIPintervalGetRoundingMode();
   if( current != roundmode )
   {
      switch(roundmode)
      {
      case SCIP_ROUND_DOWNWARDS:
         SCIPintervalSetRoundingModeDownwards();
         break;
      case SCIP_ROUND_UPWARDS:
         SCIPintervalSetRoundingModeUpwards();
         break;
      case SCIP_ROUND_NEAREST:
         SCIPintervalSetRoundingModeToNearest();
         break;
      default:
         break;
      }
   }
   nom = mpz_get_d(numerator(r->r).backend().data());
   denom = mpz_get_d(denominator(r->r).backend().data());
   //printf("computing %f/%f \n", nom, denom);

   //printf("compare with mpq value: %.*e \n", __DBL_DECIMAL_DIG__, mpq_get_d(r->r.backend().data()));

   realapprox = nom/denom;

   if( current != roundmode )
      SCIPintervalSetRoundingMode(current);

   // printf("%.*e , Roundmode %d \n",__DBL_DECIMAL_DIG__, realapprox, roundmode );

   return realapprox;
}

/** get the relaxation of a rational as a real, unfortunately you can't control the roundmode without using mpfr */
SCIP_Real RgetRealApprox(
   SCIP_Rational*        r                   /**< the rational */
   )
{
   assert(r != NULL);
   assert(r->r != NULL);

   if( r->isinf )
      return (r->r.sign() * SCIP_DEFAULT_INFINITY);

   return static_cast<SCIP_Real>(r->r);
}

void testNumericsRational(
   )
{
   if (std::numeric_limits<Rational>::is_specialized == false)
   {
   std::cout << "type " << typeid(Rational).name()  << " is not specialized for std::numeric_limits!" << std::endl;
   // ...
   }
   if(std::numeric_limits<Rational>::has_infinity)
   {
      std::cout << std::numeric_limits<Rational>::infinity() << std::endl;
   }
   else
   {
      std::cout << "type rational is not specialized for infinity" << std:: endl;
   }
   if( std::numeric_limits<Rational>::is_signed == true)
   {
      std::cout << "type rational is signed " << std::endl;
   }
   std::cout << "rounding style is " << std::numeric_limits<double>::round_style << std::endl;
}

}