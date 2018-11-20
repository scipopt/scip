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
SCIP_Rational* Rcreate(
   BMS_BLKMEM*           mem                 /**< block memory */
   )
{
   SCIP_Rational* retrat;

   BMSallocBlockMemory(mem, &retrat);

   retrat->isinf = FALSE;
   new (&retrat->r) Rational();

   return retrat;
}

/** allocate and create a rational from nominator and denominator */
SCIP_Rational* RcreateTemp(
   BMS_BUFMEM*           mem                 /**< block memory */
   )
{
   SCIP_Rational* retrat;

   BMSallocBufferMemory(mem, &retrat);

   retrat->isinf = FALSE;
   new (&retrat->r) Rational();

   return retrat;
}

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
   assert(*r != NULL);

   (*r)->r.~Rational();
   BMSfreeBlockMemory(mem, r);
}

/** delete a rational and free the allocated memory */
void RdeleteTemp(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational**       r                   /**< adress of the rational */
   )
{
   assert(*r != NULL);

   (*r)->r.~Rational();
   BMSfreeBufferMemory(mem, r);
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

/** set a rational to the value described by a string */
void RsetString(
   SCIP_Rational*        res,                /**< the result */
   const char*           desc                /**< the String describing the rational */
   )
{
   assert(res != NULL);
   assert(res->r != NULL);

   res->r = Rational(desc);

   if( 0 == strcmp(desc, "inf") )
   {
      res->r = 1;
      res->isinf = TRUE;
   }
   else if ( 0 == strcmp(desc, "-inf") )
   {
      res->r = -1;
      res->isinf = TRUE;
   }
}

/** set a rational to the value of another a real */
void RsetReal(
   SCIP_Rational*        r,
   SCIP_Real             real
   )
{
   assert(r != NULL);
   assert(r->r != NULL);

   r->isinf = false;
   r->r = real;
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
   if( 1.0 > quot )
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

    if( op1->isinf )
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

/** compute the minimum of two rationals */
void Rmax(
   SCIP_Rational*        ret,                /**< the result */
   SCIP_Rational*        r1,                 /**< the first rational */
   SCIP_Rational*        r2                  /**< the second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);
   assert(r1->r != NULL && r2->r != NULL);

   if( r1->isinf )
      ret->r = r1->r > 0 ? r1->r : r2->r;
   else if( r2->isinf )
      ret->r = r2->r > 0 ? r1->r : r1->r;
   else
      ret->r = max(r1->r, r2->r);
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

/** check if the rational is negative infinity */
SCIP_Bool RisIntegral(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return !(r->isinf) && (denominator(r->r) == 1);
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

/*
 * Dynamic Arrays
 */

/** calculate memory size for dynamically allocated arrays (copied from scip/set.c) */
static
int calcGrowSize(
   int                   initsize,           /**< initial size of array */
   SCIP_Real             growfac,            /**< growing factor of array */
   int                   num                 /**< minimum number of entries to store */
   )
{
   int size;

   assert(initsize >= 0);
   assert(growfac >= 1.0);
   assert(num >= 0);

   if( growfac == 1.0 )
      size = MAX(initsize, num);
   else
   {
      int oldsize;

      /* calculate the size with this loop, such that the resulting numbers are always the same (-> block memory) */
      initsize = MAX(initsize, 4);
      size = initsize;
      oldsize = size - 1;

      /* second condition checks against overflow */
      while( size < num && size > oldsize )
      {
         oldsize = size;
         size = (int)(growfac * size + initsize);
      }

      /* if an overflow happened, set the correct value */
      if( size <= oldsize )
         size = num;
   }

   assert(size >= initsize);
   assert(size >= num);

   return size;
}

/** creates a dynamic array of real values */
SCIP_RETCODE SCIPrationalarrayCreate(
   SCIP_RATIONALARRAY**  rationalarray,      /**< pointer to store the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(rationalarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, rationalarray) );
   (*rationalarray)->blkmem = blkmem;
   (*rationalarray)->vals = NULL;
   (*rationalarray)->valssize = 0;
   (*rationalarray)->firstidx = -1;
   (*rationalarray)->minusedidx = INT_MAX;
   (*rationalarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of real values */
SCIP_RETCODE SCIPrationalarrayCopy(
   SCIP_RATIONALARRAY**  rationalarray,      /**< pointer to store the copied real array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_RATIONALARRAY*   sourcerationalarray /**< dynamic real array to copy */
   )
{
   assert(rationalarray != NULL);
   assert(sourcerationalarray != NULL);

   SCIP_CALL( SCIPrationalarrayCreate(rationalarray, blkmem) );
   if( sourcerationalarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*rationalarray)->vals, sourcerationalarray->vals, \
                     sourcerationalarray->valssize) );
   }
   (*rationalarray)->valssize = sourcerationalarray->valssize;
   (*rationalarray)->firstidx = sourcerationalarray->firstidx;
   (*rationalarray)->minusedidx = sourcerationalarray->minusedidx;
   (*rationalarray)->maxusedidx = sourcerationalarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
SCIP_RETCODE SCIPrationalarrayFree(
   SCIP_RATIONALARRAY**  rationalarray   /**< pointer to the real array */
   )
{
   assert(rationalarray != NULL);
   assert(*rationalarray != NULL);

   RdeleteArray((*rationalarray)->blkmem, &(*rationalarray)->vals, (*rationalarray)->valssize);

   BMSfreeBlockMemory((*rationalarray)->blkmem, rationalarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPrationalarrayExtend(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(rationalarray != NULL);
   assert(rationalarray->minusedidx == INT_MAX || rationalarray->firstidx >= 0);
   assert(rationalarray->maxusedidx == INT_MIN || rationalarray->firstidx >= 0);
   assert(rationalarray->minusedidx == INT_MAX || rationalarray->minusedidx >= rationalarray->firstidx);
   assert(rationalarray->maxusedidx == INT_MIN || rationalarray->maxusedidx < rationalarray->firstidx + rationalarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   minidx = MIN(minidx, rationalarray->minusedidx);
   maxidx = MAX(maxidx, rationalarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   SCIPdebugMessage("extending rationalarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n",
      (void*)rationalarray, rationalarray->firstidx, rationalarray->valssize, rationalarray->minusedidx, rationalarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > rationalarray->valssize )
   {
      SCIP_Rational** newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = calcGrowSize(arraygrowinit, arraygrowfac, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(rationalarray->blkmem, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( rationalarray->firstidx != -1 )
      {
         for( i = 0; i < rationalarray->minusedidx - newfirstidx; ++i )
            newvals[i] = Rcreate(rationalarray->blkmem);

         /* check for possible overflow or negative value */
         assert(rationalarray->maxusedidx - rationalarray->minusedidx + 1 > 0);

         for( i = 0; i < rationalarray->maxusedidx - rationalarray->minusedidx + 1; ++i )
         {
            newvals[i + rationalarray->minusedidx - newfirstidx] =
               rationalarray->vals[i + rationalarray->minusedidx - rationalarray->firstidx];
         }

         for( i = rationalarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = Rcreate(rationalarray->blkmem);
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = Rcreate(rationalarray->blkmem);
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(rationalarray->blkmem, &rationalarray->vals, rationalarray->valssize);
      rationalarray->vals = newvals;
      rationalarray->valssize = newvalssize;
      rationalarray->firstidx = newfirstidx;
   }
   else if( rationalarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = rationalarray->valssize - nused;
      assert(nfree >= 0);
      rationalarray->firstidx = minidx - nfree/2;
      assert(rationalarray->firstidx <= minidx);
      assert(maxidx < rationalarray->firstidx + rationalarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < rationalarray->valssize; ++i )
         assert(RisZero(rationalarray->vals[i]));
#endif
   }
   else if( minidx < rationalarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = rationalarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + rationalarray->valssize);

      if( rationalarray->minusedidx <= rationalarray->maxusedidx )
      {
         int shift;

         assert(rationalarray->firstidx <= rationalarray->minusedidx);
         assert(rationalarray->maxusedidx < rationalarray->firstidx + rationalarray->valssize);

         /* shift used part of array to the right */
         shift = rationalarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = rationalarray->maxusedidx - rationalarray->firstidx; i >= rationalarray->minusedidx - rationalarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < rationalarray->valssize);
            Rset(rationalarray->vals[i + shift], rationalarray->vals[i]);
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            RsetReal(rationalarray->vals[rationalarray->minusedidx - rationalarray->firstidx + i], 0.0);
      }
      rationalarray->firstidx = newfirstidx;
   }
   else if( maxidx >= rationalarray->firstidx + rationalarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = rationalarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + rationalarray->valssize);

      if( rationalarray->minusedidx <= rationalarray->maxusedidx )
      {
         int shift;

         assert(rationalarray->firstidx <= rationalarray->minusedidx);
         assert(rationalarray->maxusedidx < rationalarray->firstidx + rationalarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - rationalarray->firstidx;
         assert(shift > 0);
         for( i = rationalarray->minusedidx - rationalarray->firstidx; i <= rationalarray->maxusedidx - rationalarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < rationalarray->valssize);
            Rset(rationalarray->vals[i - shift], rationalarray->vals[i]);
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            RsetReal(rationalarray->vals[rationalarray->maxusedidx - rationalarray->firstidx - i], 0.0);
      }
      rationalarray->firstidx = newfirstidx;
   }

   assert(minidx >= rationalarray->firstidx);
   assert(maxidx < rationalarray->firstidx + rationalarray->valssize);

   return SCIP_OKAY;
}

/* todo: exip this might be wrong */
/** clears a dynamic real array */
SCIP_RETCODE SCIPrationalarrayClear(
   SCIP_RATIONALARRAY*   rationalarray           /**< dynamic real array */
   )
{
   assert(rationalarray != NULL);

   SCIPdebugMessage("clearing rationalarray %p (firstidx=%d, size=%d, range=[%d,%d])\n",
      (void*)rationalarray, rationalarray->firstidx, rationalarray->valssize, rationalarray->minusedidx, rationalarray->maxusedidx);

   if( rationalarray->minusedidx <= rationalarray->maxusedidx )
   {
      assert(rationalarray->firstidx <= rationalarray->minusedidx);
      assert(rationalarray->maxusedidx < rationalarray->firstidx + rationalarray->valssize);
      assert(rationalarray->firstidx != -1);
      assert(rationalarray->valssize > 0);

      /* clear the used part of array */
      BMSclearMemoryArray(&rationalarray->vals[rationalarray->minusedidx - rationalarray->firstidx],
         rationalarray->maxusedidx - rationalarray->minusedidx + 1); /*lint !e866*/

      /* mark the array cleared */
      rationalarray->minusedidx = INT_MAX;
      rationalarray->maxusedidx = INT_MIN;
   }
   assert(rationalarray->minusedidx == INT_MAX);
   assert(rationalarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;

}

/** gets value of entry in dynamic array */
void SCIPrationalarrayGetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   idx,                /**< array index to get value for */
   SCIP_Rational*        result              /**< store the result */
   )
{
   assert(rationalarray != NULL);
   assert(idx >= 0);
   assert(idx >= rationalarray->minusedidx && idx <= rationalarray->maxusedidx);

   if( idx < rationalarray->minusedidx || idx > rationalarray->maxusedidx )
      RsetReal(result, 0.0);
   else
   {
      assert(rationalarray->vals != NULL);
      assert(idx - rationalarray->firstidx >= 0);
      assert(idx - rationalarray->firstidx < rationalarray->valssize);

      Rset(result, rationalarray->vals[idx - rationalarray->firstidx]);
   }
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPrationalarraySetVal(
   SCIP_RATIONALARRAY*   rationalarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   SCIP_Rational*        val                 /**< value to set array index to */
   )
{
   assert(rationalarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting rationalarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %g\n",
      (void*)rationalarray, rationalarray->firstidx, rationalarray->valssize, rationalarray->minusedidx, rationalarray->maxusedidx, idx, RgetRealApprox(val));

   if( !RisZero(val) )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPrationalarrayExtend(rationalarray, arraygrowinit, arraygrowfac, idx, idx) );
      assert(idx >= rationalarray->firstidx);
      assert(idx < rationalarray->firstidx + rationalarray->valssize);

      /* set the array value of the index */
      Rset(rationalarray->vals[idx - rationalarray->firstidx], val);

      /* update min/maxusedidx */
      rationalarray->minusedidx = MIN(rationalarray->minusedidx, idx);
      rationalarray->maxusedidx = MAX(rationalarray->maxusedidx, idx);
   }
   else if( idx >= rationalarray->firstidx && idx < rationalarray->firstidx + rationalarray->valssize )
   {
      /* set the array value of the index to zero */
      RsetReal(rationalarray->vals[idx - rationalarray->firstidx], 0.0);

      /* check, if we can tighten the min/maxusedidx */
      if( idx == rationalarray->minusedidx )
      {
         assert(rationalarray->maxusedidx >= 0);
         assert(rationalarray->maxusedidx < rationalarray->firstidx + rationalarray->valssize);
         do
         {
            Rdelete(rationalarray->blkmem, &rationalarray->vals[rationalarray->minusedidx]);
            rationalarray->minusedidx++;
         }
         while( rationalarray->minusedidx <= rationalarray->maxusedidx
            && RisZero(rationalarray->vals[rationalarray->minusedidx - rationalarray->firstidx]) );

         if( rationalarray->minusedidx > rationalarray->maxusedidx )
         {
            rationalarray->minusedidx = INT_MAX;
            rationalarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == rationalarray->maxusedidx )
      {
         assert(rationalarray->minusedidx >= 0);
         assert(rationalarray->minusedidx < rationalarray->maxusedidx);
         assert(rationalarray->maxusedidx < rationalarray->firstidx + rationalarray->valssize);
         do
         {
            Rdelete(rationalarray->blkmem, &rationalarray->vals[rationalarray->maxusedidx]);
            rationalarray->maxusedidx--;
            assert(rationalarray->minusedidx <= rationalarray->maxusedidx);
         }
         while( RisZero(rationalarray->vals[rationalarray->maxusedidx - rationalarray->firstidx]) );
      }
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPrationalarrayIncVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Rational*        incval              /**< value to increase array index */
   )
{
   assert(incval != NULL);

   if( RisZero(incval) )
      return SCIP_OKAY;
   else
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPrationalarrayExtend(rationalarray, arraygrowinit, arraygrowfac, idx, idx) );
      assert(idx >= rationalarray->firstidx);
      assert(idx < rationalarray->firstidx + rationalarray->valssize);

      /* set the array value of the index */
      Radd(rationalarray->vals[idx - rationalarray->firstidx], rationalarray->vals[idx - rationalarray->firstidx], incval);

      /* update min/maxusedidx */
      rationalarray->minusedidx = MIN(rationalarray->minusedidx, idx);
      rationalarray->maxusedidx = MAX(rationalarray->maxusedidx, idx);
   }
}


/** returns the minimal index of all stored non-zero elements */
int SCIPrationalarrayGetMinIdx(
   SCIP_RATIONALARRAY*   rationalarray           /**< dynamic real array */
   )
{
   assert(rationalarray != NULL);

   return rationalarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPrationalarrayGetMaxIdx(
   SCIP_RATIONALARRAY*   rationalarray           /**< dynamic real array */
   )
{
   assert(rationalarray != NULL);

   return rationalarray->maxusedidx;
}

}