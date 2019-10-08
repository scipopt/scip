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
   Rational* r;
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
    SCIP_Rational*        r
   )
{
   long result;
#ifdef SCIP_WITH_DEBUG_ADAPTOR
   result =  mpz_get_si(&(r->r.backend().value().data())->_mp_num);
#else
   result = (boost::multiprecision::numerator(*r->r)).convert_to<long>();
#endif
   return result;
}

static
long Rdenominator(
    SCIP_Rational*        r
   )
{
   long result;
#ifdef SCIP_WITH_DEBUG_ADAPTOR
   result mpz_get_si(&(r->r.backend().value().data())->_mp_den);
#else
    result = (boost::multiprecision::denominator(*r->r)).convert_to<long>();
#endif
   return result;
}

/** allocate and create a rational from nominator and denominator */
SCIP_RETCODE RcreateNoMem(
   SCIP_Rational**       rational            /**< pointer to the rational to create */
)
{
   SCIP_ALLOC( BMSallocMemory(rational) );
   (*rational)->r = static_cast<Rational*>(BMSallocMemoryCPP(sizeof(Rational)));

   (*rational)->isinf = FALSE;
   (*rational)->fpexact = SCIP_FPEXACT_TRUE;
   new ((*rational)->r) Rational(0);

   return SCIP_OKAY;
}

/** allocate and create a rational from nominator and denominator */
SCIP_RETCODE Rcreate(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   )
{
   SCIP_ALLOC( BMSallocBlockMemory(mem, rational) );
   (*rational)->r = static_cast<Rational*>(BMSallocMemoryCPP(sizeof(Rational)));

   (*rational)->isinf = FALSE;
   (*rational)->fpexact = SCIP_FPEXACT_TRUE;
   new ((*rational)->r) Rational(0);

   return SCIP_OKAY;
}

/** allocate and create a rational from nominator and denominator */
SCIP_RETCODE RcreateTemp(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational            /**< pointer to the rational to create */
   )
{
   BMSallocBufferMemory(mem, rational);
   (*rational)->r = static_cast<Rational*>(BMSallocMemoryCPP(sizeof(Rational)));

   (*rational)->isinf = FALSE;
   (*rational)->fpexact = SCIP_FPEXACT_TRUE;
   new ((*rational)->r) Rational(0);

   return SCIP_OKAY;
}

/** allocate and create a rational from a string in the format, e.g. "12/35" */
SCIP_RETCODE RcreateString(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,            /**< pointer to the rational to create */
   char*                 desc                /**< the String describing the rational */
)
{
   BMSallocBlockMemory(mem, rational);
   (*rational)->r = static_cast<Rational*>(BMSallocMemoryCPP(sizeof(Rational)));
   if( 0 == strcmp(desc, "inf") )
   {
      new ((*rational)->r) Rational(1);
      (*rational)->isinf = TRUE;
   }
   else if ( 0 == strcmp(desc, "-inf") )
   {
      new ((*rational)->r) Rational(-1);
      (*rational)->isinf = TRUE;
   }
   else
   {
      new ((*rational)->r) Rational(desc);
      (*rational)->isinf = FALSE;
   }
   (*rational)->fpexact = SCIP_FPEXACT_UNKNOWN;
   return SCIP_OKAY;
}

/** create an array of rationals */
SCIP_RETCODE RcreateArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      rational,           /**< pointer to the array to create */
   int                   size                /** the size of the array */
   )
{
   BMSallocBlockMemoryArray(mem, rational, size);

   for( int i = 0; i < size; ++i )
   {
      SCIP_CALL( Rcreate(mem, &(*rational)[i]) );
      (*rational)[i]->fpexact = SCIP_FPEXACT_TRUE;
   }

   return SCIP_OKAY;
}

/** create an array of rationals */
SCIP_RETCODE RcreateArrayTemp(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational***      rational,           /**< pointer to the arrat to create */
   int                   size                /** the size of the array */
   )
{
   BMSallocBufferMemoryArray(mem, rational, size);

   for( int i = 0; i < size; ++i )
   {
      SCIP_CALL( RcreateTemp(mem, &(*rational)[i]) );
      (*rational)[i]->fpexact = SCIP_FPEXACT_TRUE;
   }

   return SCIP_OKAY;
}

/** copy an array of rationals */
SCIP_RETCODE RcopyArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      target,             /**< address to copy to */
   SCIP_Rational**       src,                /**< src array */
   int                   len                 /**< size of src array */
   )
{
   int i;

   BMSduplicateBlockMemoryArray(mem, target, src, len);

   for( i = 0; i < len; ++i )
   {
      SCIP_CALL( Rcopy(mem, &(*target)[i], src[i]) );
   }

   return SCIP_OKAY;
}



/* create a copy of a rational */
SCIP_RETCODE Rcopy(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,           /**< pointer to the rational to create */
   SCIP_Rational*        src                 /**< rational to copy */
   )
{
   Rcreate(mem, rational);

   Rset(*rational, src);
   return SCIP_OKAY;
}


/** cretae a rational from an mpq_t */
#ifdef SCIP_WITH_GMP
SCIP_Rational* RcreateGMP(
   BMS_BLKMEM*           mem,                /**< block memory */
   mpq_t                 numb                /**< the gmp rational */
   )
{
   SCIP_Rational* retrat;

   BMSallocBlockMemory(mem, &retrat);
   retrat->r = static_cast<Rational*>(BMSallocMemoryCPP(sizeof(Rational)));
   retrat->isinf = FALSE;
   new (retrat->r) Rational(numb);
   retrat->fpexact = SCIP_FPEXACT_UNKNOWN;

   assert(retrat != NULL);

   return retrat;
}

/** get the underlying mpq_t* */
 mpq_t* RgetGMP(
   SCIP_Rational*        r                   /**< the rational */
   )
{
   if( r->isinf )
   {
      /** @todo exip: get proper inf value in here */
      RsetReal(r, 1e20 * r->r->sign());
   }
#ifdef SCIP_WITH_DEBUG_ADAPTOR
   return &(r->r->backend().value().data());
#else
   return &(r->r->backend().data());
#endif
}

void RsetGMP(
   SCIP_Rational*        r,
   const mpq_t           numb
   )
{
   *r->r = numb;
   r->isinf = FALSE;
   r->fpexact = SCIP_FPEXACT_UNKNOWN;
}


void RsetGMPArray(
   mpq_t*                res,
   SCIP_Rational**       src,
   int                   len
   )
{
   int i;
   for( int i = 0; i < len; i++ )
   {
      mpq_init(res[i]);
      mpq_set(res[i], *RgetGMP(src[i]));
   }
}

void RsetArrayGMP(
   SCIP_Rational**       res,
   mpq_t*                src,
   int                   len
   )
{
   int i;
   for( int i = 0; i < len; i++ )
   {
      RsetGMP(res[i], src[i]);
   }
}

void RclearGMPArray(
   mpq_t*                ar,
   int                   len
   )
{
   int i;
   for( int i = 0; i < len; i++ )
   {
      mpq_clear(ar[i]);
   }
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
      Rdelete(mem, &((*array)[i]));
   }

   BMSfreeBlockMemoryArrayNull(mem, array, size);
}

/** free an array of rationals */
void RdeleteArrayTemp(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational***      array,              /**< pointer to the array */
   int                   size                /**< size of the array */
   )
{
   assert(array != NULL);

   for( int i = size - 1; i >= 0; --i )
   {
      RdeleteTemp(mem, &((*array)[i]));
   }

   BMSfreeBufferMemoryArrayNull(mem, array);
}

/** delete a rational and free the allocated memory */
void RdeleteNoMem(
   SCIP_Rational**       r                   /**< adress of the rational */
   )
{
   assert(*r != NULL);

   (*r)->r->~Rational();
   BMSfreeMemory(&((*r)->r));
   BMSfreeMemory(r);
}

/** delete a rational and free the allocated memory */
void Rdelete(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       r                   /**< adress of the rational */
   )
{
   assert(*r != NULL);

   (*r)->r->~Rational();
   BMSfreeMemory(&((*r)->r));
   BMSfreeBlockMemory(mem, r);
}

/** delete a rational and free the allocated memory */
void RdeleteTemp(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational**       r                   /**< adress of the rational */
   )
{
   assert(*r != NULL);

   (*r)->r->~Rational();
   BMSfreeMemory(&((*r)->r));
   BMSfreeBufferMemory(mem, r);
}

/** set a rational to the value of another rational */
void Rset(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        src                 /**< the src */
   )
{
   assert(res != NULL);
   assert(res->r != NULL);

   *(res->r) = *(src->r);
   res->isinf = src->isinf;
   res->fpexact = src->fpexact;
}

/** set a rational to a nom/denom value */
void RsetInt(
   SCIP_Rational*        res,                /**< the result */
   int                   nom,                /**< the nominator */
   int                   denom               /**< the denominator */
   )
{
   assert(res != NULL);
   assert(denom != 0);

   *(res->r) = (nom/denom);
   res->isinf = FALSE;
   res->fpexact = SCIP_FPEXACT_UNKNOWN;

}

/** set a rational to the value described by a string */
void RsetString(
   SCIP_Rational*        res,                /**< the result */
   char*                 desc                /**< the String describing the rational */
   )
{
   assert(res != NULL);

   if( 0 == strcmp(desc, "inf") )
   {
      *(res->r) =  1;
      res->isinf = TRUE;
   }
   else if ( 0 == strcmp(desc, "-inf") )
   {
      *(res->r) = -1;
      res->isinf = TRUE;
   }
   else
   {
      *(res->r) = Rational(desc);
      res->isinf = FALSE;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** set a rational to the value of another a real */
void RsetReal(
   SCIP_Rational*        r,
   SCIP_Real             real
   )
{
   assert(r != NULL);

   r->isinf = FALSE;
   *(r->r) = real;
   r->fpexact = SCIP_FPEXACT_TRUE;
}

/*
 * Computing methods
 */

/* transform rational into canonical form */
void Rcanonicalize(
   SCIP_Rational*        r                   /**< rational to put in canonical form */
   )
{
   mpq_canonicalize((*r->r).backend().data());
}

/** add two rationals and save the result in res*/
void Radd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   if( op1->isinf || op2->isinf )
   {
      Rset(res, op1->isinf ? op1 : op2 );
      if( op1->r->sign() != op2->r->sign() && op1->isinf && op2->isinf )
      {
         SCIPerrorMessage("addition of pos and neg infinity not supported \n");
         SCIPABORT();
      }
   }
   else
   {
      res->isinf = FALSE;
      *(res->r) = *(op1->r) + *(op2->r);
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** add a rational and a real and save the result in res*/
void RaddReal(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   )
{
   assert(res != NULL && rat != NULL);
   if( rat->isinf )
      Rset(res, rat);
   else
   {
      res->isinf = FALSE;
      *(res->r) = *(rat->r) + real;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/*** subtract two rationals and save the result in res*/
void Rdiff(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   if( op1->isinf || op2->isinf )
   {
      op1->isinf ? Rset(res, op1) : Rneg(res, op2);
      if( op1->r->sign() != op2->r->sign() && op1->isinf && op2->isinf )
      {
         SCIPerrorMessage("addition of pos and neg infinity not supported \n");
         SCIPABORT();
      }
   }
   else
   {
      res->isinf = FALSE;
      *(res->r) = (*(op1->r)) - (*(op2->r));
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** subtract a rational and a real and save the result in res*/
void RdiffReal(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   )
{
   assert(res != NULL && rat != NULL);

   if( rat->isinf )
      Rset(res, rat);
   else
   {
      res->isinf = FALSE;
      *(res->r) = *(rat->r) - real;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
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

   absval1 = abs(*(val1->r));
   absval2 = abs(*(val2->r));
   quot = max(absval1, absval2);
   if( 1.0 > quot )
      quot = 1.0;

   *res->r = ((val1->r)-(val2->r))/quot;
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** multiply two rationals and save the result in res*/
void Rmult(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);

   if( op1->isinf || op2->isinf )
   {
      if( op1->r->is_zero() || op2->r->is_zero() )
      {
         *(res->r) = 0;
         res->isinf = FALSE;
      }
      else
      {
         SCIPerrorMessage("multiplying with infinity might produce undesired behavior \n");
         *res->r = op1->r->sign() * op2->r->sign();
         res->isinf = TRUE;
      }
   }
   else
   {
      *(res->r) = *(op1->r) * *(op2->r);
      res->isinf = FALSE;
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
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
       SCIPdebugMessage("multiplying with infinity might produce undesired behavior \n");
       if( op2 == 0.0 )
       {
          res->isinf = FALSE;
          *(res)->r = 0;
       }
       else
       {
          op2 > 0 ? Rset(res, op1) : Rneg(res, op1);
       }
    }
    else
    {
       *(res->r) = *(op1->r) * op2;
       res->isinf = FALSE;
    }
    res->fpexact = SCIP_FPEXACT_UNKNOWN;
}


/** divide two rationals and save the result in res*/
void Rdiv(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);
   assert(!RisZero(op2));
   assert(!op1->isinf && !op2->isinf);

   *res->r = *op1->r / *op2->r;
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** divide a rational and a real and save the result in res*/
void RdivReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Real             op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL);
   assert(!op1->isinf);
   assert(op2 != 0.0);

   RmultReal(res, op1, 1.0 / op2 );
}

/* Computes res += op1 * op2 and saves the result in res */
void RaddProd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);
   assert(!res->isinf);

   if( op1->isinf || op2->isinf )
   {
      if( op1->r->is_zero() || op2->r->is_zero() )
         return;
      else
      {
         SCIPerrorMessage("multiplying with infinity might produce undesired behavior \n");
         *res->r = op1->r->sign() * op2->r->sign();
         res->isinf = TRUE;
         res->fpexact = SCIP_FPEXACT_FALSE;
      }
   }
   else
   {
      res->isinf = FALSE;
      *(res->r) += *(op1->r) * (*(op2->r));
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/* Computes res -= op1 * op2 and saves the result in res */
void RdiffProd(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op1,                /**< first operand */
   SCIP_Rational*        op2                 /**< second operand */
   )
{
   assert(res != NULL && op1 != NULL && op2 != NULL);
   assert(!res->isinf);

   if( op1->isinf || op2->isinf )
   {
      if( op1->r->is_zero() || op2->r->is_zero() )
         return;
      else
      {
         SCIPerrorMessage("multiplying with infinity might produce undesired behavior \n");
         *res->r = op1->r->sign() * op2->r->sign();
         res->isinf = TRUE;
         res->fpexact = SCIP_FPEXACT_FALSE;
      }
   }
   else
   {
      res->isinf = FALSE;
      *(res->r) -= *(op1->r) * (*(op2->r));
   }
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
}

/** set res to -op */
void Rneg(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   )
{
   assert(res != NULL && op != NULL);

   *(res->r) = -(*(op->r));
   res->isinf = op->isinf;
   res->fpexact = op->fpexact;
}


/** set res to Abs(op) */
void Rabs(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   )
{
   assert(res != NULL && op != NULL);

   *(res->r) = abs(*(op->r));
   res->isinf = op->isinf;
   res->fpexact = op->fpexact;
}


/** set res to 1/op */
void Rinv(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Rational*        op                  /**< operand */
   )
{
   assert(res != NULL && op != NULL);
   assert(!op->isinf);
   assert(!op->r->is_zero());

   *(res->r) = 1 / *(op->r);
   res->isinf = FALSE;
   res->fpexact = SCIP_FPEXACT_UNKNOWN;
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
   SCIP_Bool positive;
   assert(r1 != NULL && r2 != NULL);

   if( r1->isinf )
   {
      if( *(r1->r) > 0 )
         Rset(ret, r2);
      else
         Rset(ret, r1);
   }
   else if( r2->isinf )
   {
      if( *(r2->r) > 0 )
         Rset(ret, r1);
      else
         Rset(ret, r2);
   }
   else
   {
      *(ret->r) = min(*(r1->r), *(r2->r));
      ret->isinf = FALSE;
      ret->fpexact = SCIP_FPEXACT_UNKNOWN;
   }
}

/** compute the minimum of two rationals */
void Rmax(
   SCIP_Rational*        ret,                /**< the result */
   SCIP_Rational*        r1,                 /**< the first rational */
   SCIP_Rational*        r2                  /**< the second rational */
   )
{
   SCIP_Bool positive;
   assert(r1 != NULL && r2 != NULL);

   if( r1->isinf )
   {
      if( *(r1->r) > 0 )
         Rset(ret, r1);
      else
         Rset(ret, r2);
   }
   else if( r2->isinf )
   {
      if( *(r2->r) > 0 )
         Rset(ret, r2);
      else
         Rset(ret, r1);
   }
   else
   {
      *(ret->r) = max(*(r1->r), *(r2->r));
      ret->isinf = FALSE;
      ret->fpexact = SCIP_FPEXACT_UNKNOWN;
   }
}

/** check if two rationals are equal */
SCIP_Bool RisEqual(
   SCIP_Rational*        r1,                 /**< the first rational */
   SCIP_Rational*        r2                  /**< the second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);

   if( *(r1->r) == *(r2->r) )
      return (r1->isinf == r2->isinf);

   return FALSE;
}

/** check if two rationals are equal */
SCIP_Bool RisAbsEqual(
   SCIP_Rational*        r1,                 /**< the first rational */
   SCIP_Rational*        r2                  /**< the second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);

   if( abs(*(r1->r)) == abs(*(r2->r)) )
      return (r1->isinf == r2->isinf);

   return FALSE;
}

/** check if a rational and a real are equal */
SCIP_Bool RisEqualReal(
   SCIP_Rational*        r1,                 /**< the rational */
   SCIP_Real             r2                  /**< the real */
   )
{
   assert(r1 != NULL);

   return !r1->isinf && *(r1->r) == r2;
}

/** check if real approx of rational and a real are equal */
SCIP_Bool RApproxEqualReal(
   SCIP_Rational*        r1,                 /**< the rational */
   SCIP_Real             r2                  /**< the real */
   )
{
   assert(r1 != NULL);

   if( r1->isinf )
   {
      return FALSE;
   }
   else
      return RgetRealApprox(r1) == r2;
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
      if( *(r1->r) < 0 )
         return FALSE;
      else if( r2->isinf && *(r2->r) > 0 )
         return FALSE;
      else
         return TRUE;
   }
   else if( r2->isinf )
   {
      if( *(r2->r) > 0 )
         return FALSE;
      else
         return TRUE;
   }
   else
   {
      return *(r1->r) > *(r2->r);
   }
}

/** check if the first rational is greater than the second*/
SCIP_Bool RisAbsGT(
   SCIP_Rational*        r1,                 /**< The first rational */
   SCIP_Rational*        r2                  /**< The second rational */
   )
{
   assert(r1 != NULL && r2 != NULL);

   if( r1->isinf && !r2->isinf )
      return TRUE;
   else if( r2->isinf )
      return FALSE;
   else
      return abs(*(r1->r)) > abs(*(r2->r));
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
      if( *(r1->r) > 0 )
         return TRUE;
      else if( r2->isinf && *(r2->r) < 0 )
         return TRUE;
      else
         return FALSE;
   }
   else if( r2->isinf )
   {
      if( *(r2->r) < 0 )
         return TRUE;
      else
         return FALSE;
   }
   else
   {
      return *(r1->r) >= *(r2->r);
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

   return r->r->is_zero();
}

/** check if the rational is positive */
SCIP_Bool RisPositive(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return r->r->sign() > 0;
}

/** check if the rational is negative */
SCIP_Bool RisNegative(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return r->r->sign() < 0;
}

/** check if the rational is positive infinity */
SCIP_Bool RisInfinity(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return r->isinf && r->r->sign() > 0;
}

/** check if the rational is negative infinity */
SCIP_Bool RisNegInfinity(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

   return r->isinf && r->r->sign() < 0;
}

/** check if the rational is negative infinity */
SCIP_Bool RisAbsInfinity(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);
   assert(!r->r->is_zero() || !r->isinf);

   return r->isinf;
}

/** check if the rational is negative infinity */
SCIP_Bool RisIntegral(
   SCIP_Rational*        r                   /**< the rational to check */
   )
{
   assert(r != NULL);

#ifdef SCIP_WITH_DEBUG_ADAPTOR
   return !(r->isinf) && (mpz_cmp_ui(&r->r->backend().value().data()->_mp_den, 1) == 0);
#else
   return !(r->isinf) && (mpz_cmp_ui(&r->r->backend().data()->_mp_den, 1) == 0);
#endif
}

SCIP_Bool RisFpRepresentable(
   SCIP_Rational*        r
   )
{
   assert(r != NULL);
   if( r->fpexact == SCIP_FPEXACT_TRUE )
   {
      assert(RgetRealRelax(r, SCIP_ROUND_DOWNWARDS) == RgetRealRelax(r, SCIP_ROUND_UPWARDS));
      return TRUE;
   }
   else if( r->fpexact == SCIP_FPEXACT_FALSE )
   {
      return FALSE;
   }
   else
   {
      r->fpexact = (RgetRealRelax(r, SCIP_ROUND_DOWNWARDS) == RgetRealRelax(r, SCIP_ROUND_UPWARDS)) ? SCIP_FPEXACT_TRUE : SCIP_FPEXACT_FALSE;
   }

   return r->fpexact == SCIP_FPEXACT_TRUE ? TRUE : FALSE;
}

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
   )
{
   int ret = 0;
   assert(r != NULL);

   if( r->isinf )
   {
      if( *(r->r) > 0 )
         ret = SCIPstrncpy(str, "inf", strlen);
      else
         ret = SCIPstrncpy(str, "-inf", strlen);
   }
   else
   {
      std::string s = r->r->str();
      ret = SCIPstrncpy(str, s.c_str(), strlen);
   }

   return ret;
}

/* allocates and returns a null-terminated string-representation of the rational */
const char* RgetString(
   SCIP_Rational*        r
   )
{
   assert(r != NULL);
   return mpq_get_str(0, 10, r->r->backend().data());
}

/** return the strlen of a rational number */
SCIP_Longint Rstrlen(
   SCIP_Rational*        r                /** rational to consider */
   )
{
   assert(r != NULL);
   if( r->isinf )
   {
      if( *(r->r) > 0 )
         return 3;
      else
         return 4;
   }
   else
      return (SCIP_Longint) r->r->str().length();
}

/** print a rational to command line (for debugging) */
void Rprint(
   SCIP_Rational*        r,                  /**< the rational to print */
   FILE*                 file                /**< file to print to (or NULL for std output) */
   )
{
   assert(r != NULL);
   if( r->isinf )
   {
      if( file == NULL )
         std::cout << *(r->r) << "inf" << std::endl;
      else
         fprintf(file, "%d inf", r->r->sign());
   }
   else
   {
      if( file == NULL )
         std::cout << *(r->r) << std::endl;
      else
         fprintf(file, "%s", r->r->str().c_str());
   }
}

/** print rational to file using message handler */
void Rmessage(
   SCIP_MESSAGEHDLR*     msg,                /**< message handler */
   FILE*                 file,               /**< file pointer */
   SCIP_Rational*        r                   /**< the rational to print */
   )
{
   char buf[SCIP_MAXSTRLEN];
   assert(r != NULL);

   SCIPmessageFPrintInfo(msg, file, "%s", r->r->str().c_str());
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
      return (r->r->sign() * SCIP_DEFAULT_INFINITY);
   if( r->fpexact == SCIP_FPEXACT_TRUE || roundmode == SCIP_ROUND_NEAREST )
      return RgetRealApprox(r);
   if( (roundmode == SCIP_ROUND_DOWNWARDS && RisPositive(r)) || ((roundmode == SCIP_ROUND_UPWARDS) && RisNegative(r)) )
   {
      realapprox = RgetRealApprox(r);
      SCIPdebugMessage("%.*e , Roundmode %d, R is positive: %d \n",__DBL_DECIMAL_DIG__, realapprox, roundmode, RisPositive(r) );
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

   nom = Rnumerator(r);
   denom = Rdenominator(r);

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

void Rround(
   SCIP_Rational*        retval,             /**< the resulting rounded integer */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   )
{
   mpz_t roundint;

   assert(src != NULL);
   assert(retval != NULL);

   mpz_init(roundint);
   switch (roundmode)
   {
   case SCIP_ROUND_DOWNWARDS:
      mpz_fdiv_q(roundint, mpq_numref(*RgetGMP(src)), mpq_denref(*RgetGMP(src)));
      break;
   case SCIP_ROUND_UPWARDS:
      mpz_cdiv_q(roundint, mpq_numref(*RgetGMP(src)), mpq_denref(*RgetGMP(src)));
      break;
   case SCIP_ROUND_NEAREST:
   default:
      SCIPerrorMessage("roundmode not supported for integer-rounding \n");
      SCIPABORT();
      break;
   }
   mpq_set_z(*RgetGMP(retval), roundint);

   mpz_clear(roundint);
}

/** round rational to next integer in direction of roundmode */
SCIP_Bool RroundInteger(
   long int*             retval,             /**< the resulting rounded long int */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   )
{
   SCIP_Bool success = FALSE;
   mpz_t roundint;

   assert(src != NULL);
   assert(retval != NULL);

   mpz_init(roundint);
   switch (roundmode)
   {
   case SCIP_ROUND_DOWNWARDS:
      mpz_fdiv_q(roundint, mpq_numref(*RgetGMP(src)), mpq_denref(*RgetGMP(src)));
      break;
   case SCIP_ROUND_UPWARDS:
      mpz_cdiv_q(roundint, mpq_numref(*RgetGMP(src)), mpq_denref(*RgetGMP(src)));
      break;
   case SCIP_ROUND_NEAREST:
   default:
      SCIPerrorMessage("roundmode not supported for integer-rounding \n");
      SCIPABORT();
      break;
   }

   if( mpz_fits_slong_p(roundint) )
   {
      *retval = mpz_get_si(roundint);
      success = TRUE;
   }
   mpz_clear(roundint);

   return success;
}

/** get the relaxation of a rational as a real, unfortunately you can't control the roundmode without using mpfr */
SCIP_Real RgetRealApprox(
   SCIP_Rational*        r                   /**< the rational */
   )
{
   assert(r != NULL);

   if( r->isinf )
      return (r->r->sign() * SCIP_DEFAULT_INFINITY);

   return mpq_get_d(r->r->backend().data());
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
   SCIP_Rational* r2;

   clock_t startt, endt;
   int niterations = 1000000;
   int i;
   int nrep = 0;
   double runtime = 0;
   double runtime2 = 0;
   double addval;

   RcreateNoMem(&r);
   RcreateNoMem(&r2);

   srand((unsigned int)time(NULL));

   printf("Testing time for performing tasks %d times\n", niterations);

   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      RsetReal(r, ((float)rand())/RAND_MAX);
   }
   endt = clock();
   printf(" cpu time used for setting: %e \n", ((double) (endt - startt)) / CLOCKS_PER_SEC);

   runtime = 0;
   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      double val = ((float)rand())/RAND_MAX;
      RsetReal(r, val);
      nrep += RisFpRepresentable(r) ? 1 : 0;
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
      RsetReal(r, ((float)rand())/RAND_MAX);
      addval += RgetRealRelax(r, SCIP_ROUND_DOWNWARDS);
      addval += RgetRealRelax(r, SCIP_ROUND_UPWARDS);
   }
   endt = clock();
   runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   printf(" cpu time used for rounding: %.17e, addval %e \n", runtime, addval);

   runtime = 0;
   addval = 0;
   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      RsetReal(r, ((float)rand())/RAND_MAX);
      addval += RgetRealApprox(r);
   }
   endt = clock();
   runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   printf(" cpu time used for apporx: %e, addval %e \n", runtime, addval);

   runtime = 0;
         startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      RsetReal(r, ((float)rand())/RAND_MAX);
      RsetReal(r2, ((float)rand())/RAND_MAX);
      Radd(r, r, r2);
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
      RsetReal(r, ((float)rand())/RAND_MAX);
      RsetReal(r2, ((float)rand())/RAND_MAX);
      Rmult(r, r, r2);
   }
   endt = clock();
   runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   printf(" cpu time used for multiplication: %e \n", runtime);

   RdeleteNoMem(&r);
   RdeleteNoMem(&r2);
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
      *(result->r) = search->second;
   else
      *(result->r) = 0;
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

   rationalarray->vals[idx] = (*(val->r));

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

   if( RisZero(incval) )
      return SCIP_OKAY;
   else
      rationalarray->vals[idx] += *(incval->r);

   return SCIP_OKAY;
}

}
