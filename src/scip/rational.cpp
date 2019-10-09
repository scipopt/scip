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
#include "blockmemshell/memory.h"
#include "scip/rational.h"
#include "scip/multiprecision.hpp"
#include "scip/type_message.h"
#include "scip/type_retcode.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/intervalarith.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <boost/numeric/ublas/vector_sparse.hpp>

#ifdef WITH_GMP
#include <gmp.h>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/number.hpp>
#endif

extern "C"{

using std::vector;
using std::map;

struct SCIP_Rational{
   Rational val;
   unsigned int isinf:1;
   unsigned int fpexact:2;
};

struct SCIP_RationalArray
{
   map<int, Rational>    vals;
};

/** basis status for columns and rows */
enum SCIP_fpexact
{
   SCIP_FPEXACT_UNKNOWN = 0,
   SCIP_FPEXACT_TRUE = 1,
   SCIP_FPEXACT_FALSE = 2
};

/*
 * Creation methods
 */

static
long Rnumerator(
    SCIP_Rational*        rational
   )
{
   long result;
#ifdef SCIP_WITH_DEBUG_ADAPTOR
   result =  mpz_get_si(&(rational->val.backend().value().data())->_mp_num);
#else
   result = (boost::multiprecision::numerator(rational->val)).convert_to<long>();
#endif
   return result;
}

static
long Rdenominator(
    SCIP_Rational*        rational
   )
{
   long result;
#ifdef SCIP_WITH_DEBUG_ADAPTOR
   result mpz_get_si(&(rational->val.backend().value().data())->_mp_den);
#else
    result = (boost::multiprecision::denominator(rational->val)).convert_to<long>();
#endif
   return result;
}

/** allocate and create a rational from nominator and denominator */
SCIP_RETCODE RatCreate(
   SCIP_Rational**       rational            /**< pointer to the rational to create */
)
{
   SCIP_ALLOC( BMSallocMemory(rational) );

   (*rational)->isinf = FALSE;
   (*rational)->fpexact = SCIP_FPEXACT_TRUE;
   new (&(*rational)->val) Rational(0);

   return SCIP_OKAY;
}

/** allocate and create a rational from nominator and denominator */
SCIP_RETCODE RatCreateBlock(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   )
{
   SCIP_ALLOC( BMSallocBlockMemory(mem, rational) );

   (*rational)->isinf = FALSE;
   (*rational)->fpexact = SCIP_FPEXACT_TRUE;
   new (&(*rational)->val) Rational(0);

   return SCIP_OKAY;
}

/** allocate and create a rational from nominator and denominator */
SCIP_RETCODE RatCreateBuffer(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   )
{
   BMSallocBufferMemory(mem, rational);

   (*rational)->isinf = FALSE;
   (*rational)->fpexact = SCIP_FPEXACT_TRUE;
   new (&(*rational)->val) Rational(0);

   return SCIP_OKAY;
}

/** allocate and create a rational from a string in the format, e.g. "12/35" */
SCIP_RETCODE RatCreateString(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,            /**< pointer to the rational to create */
   char*                 desc                /**< the String describing the rational */
)
{
   SCIP_CALL( RatCreateBlock(mem, rational) );

   if( 0 == strcmp(desc, "inf") )
   {
      (*rational)->val  = 1;
      (*rational)->isinf = TRUE;
   }
   else if ( 0 == strcmp(desc, "-inf") )
   {
      (*rational)->val = -1;
      (*rational)->isinf = TRUE;
   }
   else
   {
      (*rational)->val  = Rational(desc);
      (*rational)->isinf = FALSE;
   }
   (*rational)->fpexact = SCIP_FPEXACT_UNKNOWN;
   return SCIP_OKAY;
}

/** create an array of rationals */
SCIP_RETCODE RatCreateBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      rational,           /**< pointer to the array to create */
   int                   size                /** the size of the array */
   )
{
   BMSallocBlockMemoryArray(mem, rational, size);

   for( int i = 0; i < size; ++i )
   {
      SCIP_CALL( RatCreateBlock(mem, &(*rational)[i]) );
      (*rational)[i]->fpexact = SCIP_FPEXACT_TRUE;
   }

   return SCIP_OKAY;
}

/** create an array of rationals */
SCIP_RETCODE RatCreateBufferArray(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational***      rational,           /**< pointer to the arrat to create */
   int                   size                /** the size of the array */
   )
{
   BMSallocBufferMemoryArray(mem, rational, size);

   for( int i = 0; i < size; ++i )
   {
      SCIP_CALL( RatCreateBuffer(mem, &(*rational)[i]) );
      (*rational)[i]->fpexact = SCIP_FPEXACT_TRUE;
   }

   return SCIP_OKAY;
}

/** copy an array of rationals */
SCIP_RETCODE RatCopyBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      result,             /**< address to copy to */
   SCIP_Rational**       src,                /**< src array */
   int                   len                 /**< size of src array */
   )
{
   int i;

   BMSduplicateBlockMemoryArray(mem, result, src, len);
 
   for( i = 0; i < len; ++i )
   {
      SCIP_CALL( RatCopy(mem, &(*result)[i], src[i]) );
   }

   return SCIP_OKAY;
}



/* create a copy of a rational */
SCIP_RETCODE RatCopy(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       result,             /**< pointer to the rational to create */
   SCIP_Rational*        src                 /**< rational to copy */
   )
{
   SCIP_CALL( RatCreateBlock(mem, result) );

   RatSet(*result, src);

   return SCIP_OKAY;
}


/** cretae a rational from an mpq_t */
#ifdef SCIP_WITH_GMP
SCIP_RETCODE RatCreateGMP(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,           /**< pointer to the rational to create */
   mpq_t                 numb                /**< the gmp rational */
   )
{

   BMSallocBlockMemory(mem, rational);
   new (&(*rational)->val) Rational(numb);
   (*rational)->isinf = FALSE;
   (*rational)->fpexact = SCIP_FPEXACT_UNKNOWN;

   return SCIP_OKAY;
}

/** get the underlying mpq_t* */
 mpq_t* RatGetGMP(
   SCIP_Rational*        rational            /**< the rational */
   )
{
   assert(rational != NULL);

   if( rational->isinf )
   {
      /** @todo exip: get proper inf value in here */
      RatSetReal(rational, 1e20 * rational->val.sign());
   }
#ifdef SCIP_WITH_DEBUG_ADAPTOR
   return &(rational->val.backend().value().data());
#else
   return &(rational->val.backend().data());
#endif
}

/** set value of a rational from gmp data */
void RatSetGMP(
   SCIP_Rational*        rational,           /**< the rational */
   const mpq_t           numb                /**< mpq_rational */
   )
{
   rational->val = numb;
   rational->isinf = FALSE;
   rational->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** init and set value of mpq array from rational array */
void RatSetGMPArray(
   mpq_t*                mpqaaray,           /** mpq array */
   SCIP_Rational**       ratarrray,          /** array of rationals */
   int                   len                 /** array length */
   )
{
   int i;
   for( int i = 0; i < len; i++ )
   {
      mpq_init(mpqaaray[i]);
      mpq_set(mpqaaray[i], *RatGetGMP(ratarrray[i]));
   }
}

/** set value of rational array from mpq array */
void RatSetArrayGMP(
   SCIP_Rational**       ratarray,           /** array of rationals */
   mpq_t*                mpqarray,           /** mpq array */
   int                   len                 /** array length */
   )
{
   int i;
   for( int i = 0; i < len; i++ )
   {
      RatSetGMP(ratarray[i], mpqarray[i]);
   }
}

/** clear array of mpqs */
void RatClearGMPArray(
   mpq_t*                mpqarray,           /** mpq array */
   int                   len                 /** array length */
   )
{
   int i;
   for( int i = 0; i < len; i++ )
   {
      mpq_clear(mpqarray[i]);
   }
}
#endif

/** free an array of rationals */
void RatFreeBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      ratbufarray,           /**< pointer to the array */
   int                   size                /**< size of the array */
   )
{
   assert(ratbufarray != NULL);

   for( int i = 0; i < size; ++i )
   {
      RatFreeBlock(mem, &((*ratbufarray)[i]));
   }

   BMSfreeBlockMemoryArrayNull(mem, ratbufarray, size);
}

/** free an array of rationals */
void RatFreeBufferArray(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational***      ratblockarray,           /**< pointer to the array */
   int                   size                /**< size of the array */
   )
{
   assert(ratblockarray != NULL);

   for( int i = size - 1; i >= 0; --i )
   {
      RatFreeBuffer(mem, &((*ratblockarray)[i]));
   }

   BMSfreeBufferMemoryArrayNull(mem, ratblockarray);
}

/** delete a rational and free the allocated memory */
void RatFree(
   SCIP_Rational**       rational            /**< adress of the rational */
   )
{
   assert(*rational != NULL);

   (*rational)->val.~Rational();
   BMSfreeMemory(rational);
}

/** delete a rational and free the allocated memory */
void RatFreeBlock(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational            /**< adress of the rational */
   )
{
   assert(*rational != NULL);

   (*rational)->val.~Rational();
   BMSfreeBlockMemory(mem, rational);
}

/** delete a rational and free the allocated memory */
void RatFreeBuffer(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational            /**< adress of the rational */
   )
{
   assert(*rational != NULL);

   (*rational)->val.~Rational();
   BMSfreeBufferMemory(mem, rational);
}

/** set a rational to the value of another rational */
void RatSet(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        src                 /**< the src */
   )
{
   assert(res != NULL);

   res->val = src->val;
   res->isinf = src->isinf;
   res->fpexact = src->fpexact;
}

/** set a rational to a nom/denom value */
void RatSetInt(
   SCIP_Rational*        res,                /**< the result */
   int                   nom,                /**< the nominator */
   int                   denom               /**< the denominator */
   )
{
   assert(res != NULL);
   assert(denom != 0);

   res->val = (nom/denom);
   res->isinf = FALSE;
   res->fpexact = SCIP_FPEXACT_UNKNOWN;

}

/** set a rational to the value described by a string */
void RatSetString(
   SCIP_Rational*        res,                /**< the result */
   char*                 desc                /**< the String describing the rational */
   )
{
   assert(res != NULL);

   if( 0 == strcmp(desc, "inf") )
   {
      res->val =  1;
      res->isinf = TRUE;
   }
   else if ( 0 == strcmp(desc, "-inf") )
   {
      res->val = -1;
      res->isinf = TRUE;
   }
   else
   {
      res->val = Rational(desc);
      res->isinf = FALSE;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** set a rational to the value of another a real */
void RatSetReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Real             real                /**< real to set from */
   )
{
   assert(res != NULL);

   res->isinf = FALSE;
   res->val = real;
   res->fpexact = SCIP_FPEXACT_TRUE;
}

/*
 * Computing methods
 */

/* transform rational into canonical form */
void RatCanonicalize(
   SCIP_Rational*        rational            /**< rational to put in canonical form */
   )
{
   mpq_canonicalize(rational->val.backend().data());
}

/** add two rationals and save the result in res*/
void RatAdd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   if( op1->isinf || op2->isinf )
   {
      RatSet(res, op1->isinf ? op1 : op2 );
      if( op1->val.sign() != op2->val.sign() && op1->isinf && op2->isinf )
      {
         SCIPerrorMessage("addition of pos and neg infinity not supported \n");
         SCIPABORT();
      }
   }
   else
   {
      res->isinf = FALSE;
      res->val = op1->val + op2->val;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** add a rational and a real and save the result in res*/
void RatAddReal(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   )
{
   assert(res != NULL && rat != NULL);
   if( rat->isinf )
      RatSet(res, rat);
   else
   {
      res->isinf = FALSE;
      res->val = rat->val + real;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/*** subtract two rationals and save the result in res*/
void RatDiff(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   if( op1->isinf || op2->isinf )
   {
      op1->isinf ? RatSet(res, op1) : RatNegate(res, op2);
      if( op1->val.sign() != op2->val.sign() && op1->isinf && op2->isinf )
      {
         SCIPerrorMessage("addition of pos and neg infinity not supported \n");
         SCIPABORT();
      }
   }
   else
   {
      res->isinf = FALSE;
      res->val = (op1->val) - (op2->val);
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** subtract a rational and a real and save the result in res*/
void RatDiffReal(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   )
{
   assert(res != NULL && rat != NULL);

   if( rat->isinf )
      RatSet(res, rat);
   else
   {
      res->isinf = FALSE;
      res->val = rat->val - real;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
void RatRelDiff(
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

   absval1 = abs(val1->val);
   absval2 = abs(val2->val);
   quot = max(absval1, absval2);
   if( 1.0 > quot )
      quot = 1.0;

   res->val = ((val1->val)-(val2->val))/quot;
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** multiply two rationals and save the result in res*/
void RatMult(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   if( op1->isinf || op2->isinf )
   {
      if( op1->val.is_zero() || op2->val.is_zero() )
      {
         res->val = 0;
         res->isinf = FALSE;
      }
      else
      {
         SCIPerrorMessage("multiplying with infinity might produce undesired behavior \n");
         res->val = op1->val.sign() * op2->val.sign();
         res->isinf = TRUE;
      }
   }
   else
   {
      res->val = op1->val * op2->val;
      res->isinf = FALSE;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** multiply a rational and a real and save the result in res*/
void RatMultReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   )
{
    assert(res != NULL && op1 != NULL);

    if( op1->isinf )
    {
       SCIPdebugMessage("multiplying with infinity might produce undesired behavior \n");
       if( op2 == 0.0 )
       {
          res->isinf = FALSE;
          res->val = 0;
       }
       else
       {
          op2 > 0 ? RatSet(res, op1) : RatNegate(res, op1);
       }
    }
    else
    {
       res->val = op1->val * op2;
       res->isinf = FALSE;
    }
    res->fpexact = SCIP_FPEXACT_UNKNOWN;
}


/** divide two rationals and save the result in res*/
void RatDiv(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);
   assert(!RatIsZero(op2));
   assert(!op1->isinf && !op2->isinf);

   res->val = op1->val / op2->val;
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** divide a rational and a real and save the result in res*/
void RatDivReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL);
   assert(!op1->isinf);
   assert(op2 != 0.0);

   RatMultReal(res, op1, 1.0 / op2 );
}

/* Computes res += op1 * op2 and saves the result in res */
void RatAddProd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);
   assert(!res->isinf);

   if( op1->isinf || op2->isinf )
   {
      if( op1->val.is_zero() || op2->val.is_zero() )
         return;
      else
      {
         SCIPerrorMessage("multiplying with infinity might produce undesired behavior \n");
         res->val = op1->val.sign() * op2->val.sign();
         res->isinf = TRUE;
         res->fpexact = SCIP_FPEXACT_FALSE;
      }
   }
   else
   {
      res->isinf = FALSE;
      res->val += op1->val * op2->val;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/* Computes res -= op1 * op2 and saves the result in res */
void RatDiffProd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);
   assert(!res->isinf);

   if( op1->isinf || op2->isinf )
   {
      if( op1->val.is_zero() || op2->val.is_zero() )
         return;
      else
      {
         SCIPerrorMessage("multiplying with infinity might produce undesired behavior \n");
         res->val = op1->val.sign() * op2->val.sign();
         res->isinf = TRUE;
         res->fpexact = SCIP_FPEXACT_FALSE;
      }
   }
   else
   {
      res->isinf = FALSE;
      res->val -= op1->val * op2->val;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** set res to -op */
void RatNegate(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   )
{
   assert(res != NULL && op != NULL);

   res->val = -op->val;
   res->isinf = op->isinf;
   res->fpexact = op->fpexact;
}


/** set res to Abs(op) */
void RatAbs(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   )
{
   assert(res != NULL && op != NULL);

   res->val = abs(op->val);
   res->isinf = op->isinf;
   res->fpexact = op->fpexact;
}


/** set res to 1/op */
void RatInvert(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   )
{
   assert(res != NULL && op != NULL);
   assert(!op->isinf);
   assert(!op->val.is_zero());

   res->val = 1 / op->val;
   res->isinf = FALSE;
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/*
 * Comparisoon methods
 */

/** compute the minimum of two rationals */
void RatMIN(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< the first rational */
   SCIP_Rational*        op2                 /**< the second rational */
   )
{
   SCIP_Bool positive;
   assert(op1 != NULL && op2 != NULL);

   if( op1->isinf )
   {
      if( op1->val > 0 )
         RatSet(res, op2);
      else
         RatSet(res, op1);
   }
   else if( op2->isinf )
   {
      if( op2->val > 0 )
         RatSet(res, op1);
      else
         RatSet(res, op2);
   }
   else
   {
      res->val = min(op1->val, op2->val);
      res->isinf = FALSE;
      res->fpexact = SCIP_FPEXACT_UNKNOWN;
   }
}

/** compute the minimum of two rationals */
void RatMAX(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< the first rational */
   SCIP_Rational*        op2                 /**< the second rational */
   )
{
   SCIP_Bool positive;
   assert(op1 != NULL && op2 != NULL);

   if( op1->isinf )
   {
      if( op1->val > 0 )
         RatSet(res, op1);
      else
         RatSet(res, op2);
   }
   else if( op2->isinf )
   {
      if( op2->val > 0 )
         RatSet(res, op2);
      else
         RatSet(res, op1);
   }
   else
   {
      res->val = max(op1->val, op2->val);
      res->isinf = FALSE;
      res->fpexact = SCIP_FPEXACT_UNKNOWN;
   }
}

/** check if two rationals are equal */
SCIP_Bool RatIsEqual(
   SCIP_Rational*        rat1,               /**< the first rational */
   SCIP_Rational*        rat2                /**< the second rational */
   )
{
   assert(rat1 != NULL && rat2 != NULL);

   if( rat1->val == rat2->val )
      return (rat1->isinf == rat2->isinf);

   return FALSE;
}

/** check if two rationals are equal */
SCIP_Bool RatIsAbsEqual(
   SCIP_Rational*        rat1,               /**< the first rational */
   SCIP_Rational*        rat2                /**< the second rational */
   )
{
   assert(rat1 != NULL && rat2 != NULL);

   if( abs(rat1->val) == abs(rat2->val) )
      return rat1->isinf == rat2->isinf;

   return FALSE;
}

/** check if a rational and a real are equal */
SCIP_Bool RatIsEqualReal(
   SCIP_Rational*        rat,                /**< the rational */
   SCIP_Real             real                /**< the real */
   )
{
   assert(rat != NULL);

   return !rat->isinf && rat->val == real;
}

/** check if real approx of rational and a real are equal */
SCIP_Bool RatIsApproxEqualReal(
   SCIP_Rational*        rat,               /**< the rational */
   SCIP_Real             real               /**< the real */
   )
{
   assert(rat != NULL);

   if( rat->isinf )
   {
      return FALSE;
   }
   else
      return RatApproxReal(rat) == real;
}

/** check if the first rational is greater than the second*/
SCIP_Bool RatIsGT(
   SCIP_Rational*        rat1,               /**< The first rational */
   SCIP_Rational*        rat2                /**< The second rational */
   )
{
   assert(rat1 != NULL && rat2 != NULL);

   if( rat1->isinf )
   {
      if( rat1->val < 0 )
         return FALSE;
      else if( rat2->isinf && rat2->val > 0 )
         return FALSE;
      else
         return TRUE;
   }
   else if( rat2->isinf )
   {
      if( rat2->val > 0 )
         return FALSE;
      else
         return TRUE;
   }
   else
   {
      return rat1->val > rat2->val;
   }
}

/** check if the first rational is greater than the second*/
SCIP_Bool RatIsAbsGT(
   SCIP_Rational*        rat1,                 /**< The first rational */
   SCIP_Rational*        rat2                  /**< The second rational */
   )
{
   assert(rat1 != NULL && rat2 != NULL);

   if( rat1->isinf && !rat2->isinf )
      return TRUE;
   else if( rat2->isinf )
      return FALSE;
   else
      return abs(rat1->val) > abs(rat2->val);
}

/** check if the first rational is smaller than the second*/
SCIP_Bool RatIsLT(
   SCIP_Rational*        rat1,               /**< The first rational */
   SCIP_Rational*        rat2                /**< The second rational */
   )
{
   assert(rat1 != NULL && rat2 != NULL);

   return RatIsGT(rat2, rat1);
}

/** check if the first rational is greater or equal than the second*/
EXTERN
SCIP_Bool RatIsGE(
   SCIP_Rational*        rat1,               /**< The first rational */
   SCIP_Rational*        rat2                /**< The second rational */
   )
{
   assert(rat1 != NULL && rat2 != NULL);

   if( rat1->isinf )
   {
      if( rat1->val > 0 )
         return TRUE;
      else if( rat2->isinf && (rat2->val) < 0 )
         return TRUE;
      else
         return FALSE;
   }
   else if( rat2->isinf )
   {
      if( rat2->val < 0 )
         return TRUE;
      else
         return FALSE;
   }
   else
   {
      return rat1->val >= rat2->val;
   }
}

/** check if the first rational is less or equal than the second*/
EXTERN
SCIP_Bool RatIsLE(
   SCIP_Rational*        rat1,               /**< The first rational */
   SCIP_Rational*        rat2                /**< The second rational */
   )
{
   assert(rat1 != NULL && rat2 != NULL);

   return RatIsGE(rat2, rat1);
}


/** check if the rational is zero */
SCIP_Bool RatIsZero(
   SCIP_Rational*        rational            /**< the rational to check */
   )
{
   assert(rational != NULL);

   return rational->val.is_zero();
}

/** check if the rational is positive */
SCIP_Bool RatIsPositive(
   SCIP_Rational*        rational            /**< the rational to check */
   )
{
   assert(rational != NULL);

   return rational->val.sign() > 0;
}

/** check if the rational is negative */
SCIP_Bool RatIsNegative(
   SCIP_Rational*        rational            /**< the rational to check */
   )
{
   assert(rational != NULL);

   return rational->val.sign() < 0;
}

/** check if the rational is positive infinity */
SCIP_Bool RatIsInfinity(
   SCIP_Rational*        rational             /**< the rational to check */
   )
{
   assert(rational != NULL);

   return rational->isinf && rational->val.sign() > 0;
}

/** check if the rational is negative infinity */
SCIP_Bool RatIsNegInfinity(
   SCIP_Rational*        rational             /**< the rational to check */
   )
{
   assert(rational != NULL);

   return rational->isinf && rational->val.sign() < 0;
}

/** check if the rational is negative infinity */
SCIP_Bool RatIsAbsInfinity(
   SCIP_Rational*        rational            /**< the rational to check */
   )
{
   assert(rational != NULL);
   assert(!rational->val.is_zero() || !rational->isinf);

   return rational->isinf;
}

/** check if the rational is negative infinity */
SCIP_Bool RatIsIntegral(
   SCIP_Rational*        rational             /**< the rational to check */
   )
{
   assert(rational != NULL);

#ifdef SCIP_WITH_DEBUG_ADAPTOR
   return !(rational->isinf) && (mpz_cmp_ui(&rational->val.backend().value().data()->_mp_den, 1) == 0);
#else
   return !(rational->isinf) && (mpz_cmp_ui(&rational->val.backend().data()->_mp_den, 1) == 0);
#endif
}

/** check if rational is exactly representable as double */
SCIP_Bool RatIsFpRepresentable(
   SCIP_Rational*        rational             /**< the rational to check */
   )
{
   assert(rational != NULL);
   if( rational->fpexact == SCIP_FPEXACT_TRUE )
   {
      assert(RatRoundReal(rational, SCIP_ROUND_DOWNWARDS) == RatRoundReal(rational, SCIP_ROUND_UPWARDS));
      return TRUE;
   }
   else if( rational->fpexact == SCIP_FPEXACT_FALSE )
   {
      return FALSE;
   }
   else
   {
      rational->fpexact = (RatRoundReal(rational, SCIP_ROUND_DOWNWARDS) 
         == RatRoundReal(rational, SCIP_ROUND_UPWARDS)) ? SCIP_FPEXACT_TRUE : SCIP_FPEXACT_FALSE;
   }

   return rational->fpexact == SCIP_FPEXACT_TRUE ? TRUE : FALSE;
}

/*
 * Printing/Conversion methods
 */

/** convert a Rational to a string for printing, returns the number of copied characters.
 * If return value is equal to strlen, it means the string was truncated.
 */
int RatToString(
   SCIP_Rational*        rational,           /**< the rational to print */
   char*                 str,                /**< the string to save the rational in */
   int                   strlen              /**< maximal length that can be copied to str */
   )
{
   int ret = 0;
   assert(rational != NULL);

   if( rational->isinf )
   {
      if( rational->val > 0 )
         ret = SCIPstrncpy(str, "inf", strlen);
      else
         ret = SCIPstrncpy(str, "-inf", strlen);
   }
   else
   {
      std::string s = rational->val.str();
      ret = SCIPstrncpy(str, s.c_str(), strlen);
   }
   if( ret == strlen )
   {
      SCIPerrorMessage("Rational string to long to fit in buffer");
      RatPrint(rational);
   }

   return ret;
}

/* allocates and returns a null-terminated string-representation of the rational
   warning! have to free this yourself, does not check infinity! only for debugging. */
const char* RatGetString(
   SCIP_Rational*        rational           /**< the rational to print */
   )
{
   assert(rational != NULL);
   return mpq_get_str(0, 10, rational->val.backend().data());
}

/** return the strlen of a rational number */
SCIP_Longint RatStrlen(
   SCIP_Rational*        rational           /** rational to consider */
   )
{
   assert(rational != NULL);
   if( rational->isinf )
   {
      if( rational->val > 0 )
         return 3;
      else
         return 4;
   }
   else
      return (SCIP_Longint) rational->val.str().length();
}

/** print a rational to command line (for debugging) */
void RatPrint(
   SCIP_Rational*        rational            /**< the rational to print */
   )
{
   assert(rational != NULL);
   if( rational->isinf )
      std::cout << rational->val << "inf" << std::endl;
   else
      std::cout << rational->val << std::endl;
}

/** print rational to file using message handler */
void RatMessage(
   SCIP_MESSAGEHDLR*     msg,                /**< message handler */
   FILE*                 file,               /**< file pointer */
   SCIP_Rational*        rational            /**< the rational to print */
   )
{
   char buf[SCIP_MAXSTRLEN];
   assert(rational != NULL);

   if( SCIP_MAXSTRLEN == RatToString(rational, buf, SCIP_MAXSTRLEN) )
      SCIPerrorMessage("WARNING: Rational does not fit in line \n.");

   SCIPmessageFPrintInfo(msg, file, "%s", buf);
}

/** get the relaxation of a rational as a real, unfortunately you can't control the roundmode without using mpfr */
SCIP_Real RatRoundReal(
   SCIP_Rational*        rational,           /**< the rational */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   )
{
   SCIP_Real realapprox;
   SCIP_Real nom, denom;
   SCIP_ROUNDMODE current;

   assert(rational != NULL);

   if( rational->isinf )
      return (rational->val.sign() * SCIP_DEFAULT_INFINITY);
   if( rational->fpexact == SCIP_FPEXACT_TRUE || roundmode == SCIP_ROUND_NEAREST )
      return RatApproxReal(rational);
   if( (roundmode == SCIP_ROUND_DOWNWARDS && RatIsPositive(rational)) || ((roundmode == SCIP_ROUND_UPWARDS) && RatIsNegative(rational)) )
   {
      realapprox = RatApproxReal(rational);
      SCIPdebugMessage("%.*e , Roundmode %d, R is positive: %d \n",__DBL_DECIMAL_DIG__, realapprox, roundmode, RatIsPositive(rational) );
      return realapprox;
   }


/* #ifdef SCIP_WITH_MPFR
   {
      mpfr_t valmpfr;
      mpq_t* val;

      val = RgetGMP(r);
      switch(roundmode)
      {
         case SCIP_ROUND_DOWNWARDS:
            mpfr_init_set_q(valmpfr, *val, MPFR_RNDD);
            break;
         case SCIP_ROUND_UPWARDS:
            mpfr_init_set_q(valmpfr, *val, MPFR_RNDU);
            break;
         case SCIP_ROUND_NEAREST:
            mpfr_init_set_q(valmpfr, *val, MPFR_RNDN);
            break;
         default:
            break;
      }
      realapprox = (SCIP_Real) mpfr_get_d(valmpfr, MPFR_RNDN);
      mpfr_clear(valmpfr);
   }
#else */

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

   nom = Rnumerator(rational);
   denom = Rdenominator(rational);

   SCIPdebugMessage("computing %f/%f \n", nom, denom);

   realapprox = nom/denom;
   //realapprox = RgetRealApprox(r);

   if( current != roundmode )
      SCIPintervalSetRoundingMode(current);
/* #endif
 */
   SCIPdebugMessage("%.*e , Roundmode %d \n",__DBL_DECIMAL_DIG__, realapprox, roundmode );

   return realapprox;
}

void RatRound(
   SCIP_Rational*        res,                /**< the resulting rounded integer */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   )
{
   mpz_t roundint;

   assert(src != NULL);
   assert(res != NULL);

   mpz_init(roundint);
   switch (roundmode)
   {
   case SCIP_ROUND_DOWNWARDS:
      mpz_fdiv_q(roundint, mpq_numref(*RatGetGMP(src)), mpq_denref(*RatGetGMP(src)));
      break;
   case SCIP_ROUND_UPWARDS:
      mpz_cdiv_q(roundint, mpq_numref(*RatGetGMP(src)), mpq_denref(*RatGetGMP(src)));
      break;
   case SCIP_ROUND_NEAREST:
   default:
      SCIPerrorMessage("roundmode not supported for integer-rounding \n");
      SCIPABORT();
      break;
   }
   mpq_set_z(*RatGetGMP(res), roundint);

   mpz_clear(roundint);
}

/** round rational to next integer in direction of roundmode */
SCIP_Bool RatRoundInteger(
   long int*             res,                /**< the resulting rounded long int */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   )
{
   SCIP_Bool success = FALSE;
   mpz_t roundint;

   assert(src != NULL);
   assert(res != NULL);

   mpz_init(roundint);
   switch (roundmode)
   {
   case SCIP_ROUND_DOWNWARDS:
      mpz_fdiv_q(roundint, mpq_numref(*RatGetGMP(src)), mpq_denref(*RatGetGMP(src)));
      break;
   case SCIP_ROUND_UPWARDS:
      mpz_cdiv_q(roundint, mpq_numref(*RatGetGMP(src)), mpq_denref(*RatGetGMP(src)));
      break;
   case SCIP_ROUND_NEAREST:
   default:
      SCIPerrorMessage("roundmode not supported for integer-rounding \n");
      SCIPABORT();
      break;
   }

   if( mpz_fits_slong_p(roundint) )
   {
      *res = mpz_get_si(roundint);
      success = TRUE;
   }
   mpz_clear(roundint);

   return success;
}

/** get the relaxation of a rational as a real, unfortunately you can't control the roundmode without using mpfr */
SCIP_Real RatApproxReal(
   SCIP_Rational*        rational            /**< the rational */
   )
{
   assert(rational != NULL);

   if( rational->isinf )
      return (rational->val.sign() * SCIP_DEFAULT_INFINITY);

   return mpq_get_d(rational->val.backend().data());
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

void testRuntimesRational(
   )
{
   SCIP_Rational* r;
   SCIP_Rational* rat2;

   clock_t startt, endt;
   int niterations = 1000000;
   int i;
   int nrep = 0;
   double runtime = 0;
   double runtime2 = 0;
   double addval;

   RatCreate(&r);
   RatCreate(&rat2);

   srand((unsigned int)time(NULL));

   printf("Testing time for performing tasks %d times\n", niterations);

   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
   }
   endt = clock();
   printf(" cpu time used for setting: %e \n", ((double) (endt - startt)) / CLOCKS_PER_SEC);

   runtime = 0;
   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      double val = ((float)rand())/RAND_MAX;
      RatSetReal(r, val);
      nrep += RatIsFpRepresentable(r) ? 1 : 0;
   }
   endt = clock();
   runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   printf(" cpu time used for checking fp-rep: %e \n", runtime);
   if( nrep != niterations )
   {
      printf(" error! fprep, %d iterations but only %d are representable \n", niterations, nrep);
   }

   runtime = 0;
   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
      addval += RatRoundReal(r, SCIP_ROUND_DOWNWARDS);
      addval += RatRoundReal(r, SCIP_ROUND_UPWARDS);
   }
   endt = clock();
   runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   printf(" cpu time used for rounding: %.17e, addval %e \n", runtime, addval);

   runtime = 0;
   addval = 0;
   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
      addval += RatApproxReal(r);
   }
   endt = clock();
   runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   printf(" cpu time used for apporx: %e, addval %e \n", runtime, addval);

   runtime = 0;
         startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
      RatSetReal(rat2, ((float)rand())/RAND_MAX);
      RatAdd(r, r, rat2);
   }
   endt = clock();
   runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   printf(" cpu time used for adding: %e \n", runtime);

   {
   SCIP_Real reals[niterations];
   runtime = 0;
   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      reals[i] = ((float)rand())/RAND_MAX;
      reals[i] += ((float)rand())/RAND_MAX;
   }
   endt = clock();
   runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   printf(" cpu time used for adding reals: %e \n", runtime);
   }

   runtime = 0;
   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
      RatSetReal(rat2, ((float)rand())/RAND_MAX);
      RatMult(r, r, rat2);
   }
   endt = clock();
   runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   printf(" cpu time used for multiplication: %e \n", runtime);

   RatFree(&r);
   RatFree(&rat2);
}

/*
 * Dynamic Arrays
 */

/** creates a dynamic array of real values */
SCIP_RETCODE SCIPrationalarrayCreate(
   SCIP_RATIONALARRAY**  rationalarray,      /**< pointer to store the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(rationalarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, rationalarray) );
   new (&(*rationalarray)->vals) map<int,Rational>();

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
   (*rationalarray)->vals = sourcerationalarray->vals;

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
SCIP_RETCODE SCIPrationalarrayFree(
   SCIP_RATIONALARRAY**  rationalarray,   /**< pointer to the real array */
   BMS_BLKMEM*           blkmem          /**< block memory */
   )
{
   assert(rationalarray != NULL);
   assert(*rationalarray != NULL);

   (*rationalarray)->vals.~map();
   BMSfreeBlockMemory(blkmem, rationalarray);

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
   auto search = rationalarray->vals.find(idx);
   if(  search != rationalarray->vals.end() )
      result->val = search->second;
   else
      result->val = 0;
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPrationalarraySetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic rational array */
   int                   idx,                /**< array index to set value for */
   SCIP_Rational*        val                       /**< value to set array index to */
   )
{
   assert(rationalarray != NULL);
   assert(idx >= 0);

   rationalarray->vals[idx] = val->val;

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPrationalarrayIncVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic real array */
   int                   idx,                /**< array index to increase value for */
    SCIP_Rational*  incval              /**< value to increase array index */
   )
{
   assert(incval != NULL);

   if( RatIsZero(incval) )
      return SCIP_OKAY;
   else
      rationalarray->vals[idx] += incval->val;

   return SCIP_OKAY;
}

}
