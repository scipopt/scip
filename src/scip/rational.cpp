/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file rational.cpp
 * @brief rational wrapper file to use the multiprecision rationals class in SCIP
 * @ingroup INTERNALAPI
 * @author  Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "blockmemshell/memory.h"
#include "scip/rational.h"
#include "scip/struct_rational.h"
#include "scip/multiprecision.hpp"
#include "scip/type_message.h"
#include "scip/type_retcode.h"
#include "scip/type_rational.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/intervalarith.h"
#include "scip/set.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <numeric>
#include <string.h>
#include <ostream>

#ifdef SCIP_WITH_BOOST
#include <boost/format.hpp>
#ifdef SCIP_WITH_GMP
#include <gmp.h>
#include <boost/multiprecision/gmp.hpp>
#endif
#include <boost/multiprecision/number.hpp>
#else
#endif

extern "C" {

static const char posinf[4] = "inf";
static const char neginf[5] = "-inf";
static char stringbuf[SCIP_MAXSTRLEN];
static SCIP_Rational buffer;
static SCIP_Real infinity = SCIP_DEFAULT_INFINITY; /* values above this are considered to be infinite */

/*
 * Creation methods
 */
#ifdef SCIP_WITH_BOOST
static
SCIP_Longint Rnumerator(
   SCIP_Rational*        rational
   )
{
   long result;

   result = (boost::multiprecision::numerator(rational->val)).convert_to<long>();

   return result;
}

static
SCIP_Longint Rdenominator(
   SCIP_Rational*        rational
   )
{
   long result;

    result = (boost::multiprecision::denominator(rational->val)).convert_to<long>();

   return result;
}
#else
static
SCIP_Longint Rnumerator(
   SCIP_Rational*        rational
   )
{
   return rational->val.val;
}

static
SCIP_Longint Rdenominator(
   SCIP_Rational*        rational
   )
{
   return 1.0;
}
#endif

/** allocate and create a rational from nominator and denominator */
SCIP_RETCODE RatCreate(
   SCIP_Rational**       rational            /**< pointer to the rational to create */
)
{
   SCIP_ALLOC( BMSallocMemory(rational) );

   (*rational)->isinf = FALSE;
   (*rational)->isfprepresentable = SCIP_ISFPREPRESENTABLE_TRUE;
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
   (*rational)->isfprepresentable = SCIP_ISFPREPRESENTABLE_TRUE;
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
   (*rational)->isfprepresentable = SCIP_ISFPREPRESENTABLE_TRUE;
   new (&(*rational)->val) Rational(0);

   return SCIP_OKAY;
}

/** allocate and create a rational from a string in the format, e.g. "12/35" */
SCIP_RETCODE RatCreateString(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational**       rational,           /**< pointer to the rational to create */
   char*                 desc                /**< the String describing the rational */
)
{
   SCIP_CALL( RatCreateBlock(mem, rational) );

   if( 0 == strcmp(desc, "inf") )
   {
      (*rational)->val = 1;
      (*rational)->isinf = TRUE;
   }
   else if ( 0 == strcmp(desc, "-inf") )
   {
      (*rational)->val = -1;
      (*rational)->isinf = TRUE;
   }
   else
   {
      (*rational)->val = Rational(desc);
      (*rational)->isinf = FALSE;
   }
   (*rational)->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
   return SCIP_OKAY;
}

/** create an array of rationals */
SCIP_RETCODE RatCreateArray(
   SCIP_Rational***      rational,           /**< pointer to the array to create */
   int                   size                /** the size of the array */
   )
{
   BMSallocMemoryArray(rational, size);

   for( int i = 0; i < size; ++i )
   {
      SCIP_CALL( RatCreate(&(*rational)[i]) );
      (*rational)[i]->isfprepresentable = SCIP_ISFPREPRESENTABLE_TRUE;
   }

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
      (*rational)[i]->isfprepresentable = SCIP_ISFPREPRESENTABLE_TRUE;
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
      (*rational)[i]->isfprepresentable = SCIP_ISFPREPRESENTABLE_TRUE;
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



/** creates a copy of a rational */
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

#ifdef SCIP_WITH_BOOST
/** creates a rational from an mpq_t */
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
   (*rational)->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;

   return SCIP_OKAY;
}

/** gets the underlying mpq_t* */
mpq_t* RatGetGMP(
   SCIP_Rational*        rational            /**< the rational */
   )
{
   assert(rational != NULL);

   if( rational->isinf )
   {
      /** @todo exip: get proper inf value in here */
      RatSetReal(rational, 1e100 * rational->val.sign());
      rational->isinf = TRUE;
   }

   return &(rational->val.backend().data());
}

/** set value of a rational from gmp data */
void RatSetGMP(
   SCIP_Rational*        rational,           /**< the rational */
   const mpq_t           numb                /**< mpq_rational */
   )
{
   rational->val = numb;
   rational->isinf = FALSE;
   rational->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
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
#else
   /** get the underlying mpq_t* */
 mpq_t* RatGetGMP(
   SCIP_Rational*        rational            /**< the rational */
   )
{
   return 0;
}
#endif

/* transform rational into canonical form */
/** @todo exip: this does not work with cpp_rational currently */
void RatCanonicalize(
   SCIP_Rational*        rational            /**< rational to put in canonical form */
   )
{
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_BOOST)
   mpq_canonicalize(rational->val.backend().data());
#endif
}

/** free an array of rationals */
void RatFreeArray(
   SCIP_Rational***      ratarray,           /**< pointer to the array */
   int                   size                /**< size of the array */
   )
{
   assert(ratarray != NULL);

   for( int i = 0; i < size; ++i )
   {
      RatFree(&((*ratarray)[i]));
   }

   BMSfreeMemoryArrayNull(ratarray);
}

/** free an array of rationals */
void RatFreeBlockArray(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_Rational***      ratblockarray,      /**< pointer to the array */
   int                   size                /**< size of the array */
   )
{
   assert(ratblockarray != NULL);

   for( int i = 0; i < size; ++i )
   {
      RatFreeBlock(mem, &((*ratblockarray)[i]));
   }

   BMSfreeBlockMemoryArrayNull(mem, ratblockarray, size);
}

/** free an array of rationals */
void RatFreeBufferArray(
   BMS_BUFMEM*           mem,                /**< block memory */
   SCIP_Rational***      ratbufarray,        /**< pointer to the array */
   int                   size                /**< size of the array */
   )
{
   assert(ratbufarray != NULL);

   for( int i = size - 1; i >= 0; --i )
   {
      RatFreeBuffer(mem, &((*ratbufarray)[i]));
   }

   BMSfreeBufferMemoryArrayNull(mem, ratbufarray);
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
   res->isfprepresentable = src->isfprepresentable;
}

/** set a rational to a nom/denom value */
void RatSetInt(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Longint          nom,                /**< the nominator */
   SCIP_Longint          denom               /**< the denominator */
   )
{
   char buf[SCIP_MAXSTRLEN];

   assert(res != NULL);
   assert(denom != 0);

   if( denom < 0 )
   {
      nom *= -1;
      denom *= -1;
   }

   res->val = Rational(nom, denom);

   res->isinf = FALSE;
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
}

/** set a rational to the value described by a string */
void RatSetString(
   SCIP_Rational*        res,                /**< the result */
   const char*           desc                /**< the String describing the rational */
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
      std::string s(desc);
      /* case 1: string is given in nom/den format */
      if( s.find('.') == std::string::npos )
      {
         res->val = Rational(desc);
         res->isinf = FALSE;
      }
      /* case 2: string is given as base-10 decimal number */
      else
      {
         // std::cout << s << std::endl;
         if( s[0] == '.' )
            s.insert(0, "0");
         size_t pos = s.find('.');
         size_t exp = s.length() - 1 - pos;
         std::string den("1");
         for( int i = 0; i < exp; ++i )
            den.append("0");

         s.erase(pos, 1);
         assert(std::all_of(s.begin()+1, s.end(), ::isdigit));
         s.append("/");
         s.append(den);
         res->val = Rational(s);
         res->isinf = FALSE;
         // RatPrint(res);
      }
   }
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
}

/** set a rational to the value of another a real */
void RatSetReal(
   SCIP_Rational*        res,                /**< the result */
   SCIP_Real             real                /**< real to set from */
   )
{
   assert(res != NULL);

   if( REALABS(real) >= infinity )
   {
      res->isinf = TRUE;
      res->val = real > 0 ? 1 : -1;
      res->isfprepresentable = TRUE;
   }
   else
   {
      res->isinf = FALSE;
      res->val = real;
      res->isfprepresentable = SCIP_ISFPREPRESENTABLE_TRUE;
   }
}

/*
 * Computing methods
 */

/** add two rationals and save the result in res */
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
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
}

/** add a rational and a real and save the result in res */
void RatAddReal(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   )
{
   assert(res != NULL && rat != NULL);
   if( rat->isinf )
      RatSet(res, rat);
   else if( REALABS(real) >= infinity )
   {
      RatSetReal(res, real);
   }
   else
   {
      res->isinf = FALSE;
      res->val = rat->val + real;
   }
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
}

/*** subtract two rationals and save the result in res */
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
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
}

/** subtract a rational and a real and save the result in res */
void RatDiffReal(
   SCIP_Rational*        res,                /**< The result */
   SCIP_Rational*        rat,                /**< rational number */
   SCIP_Real             real                /**< real number */
   )
{
   assert(res != NULL && rat != NULL);

   RatAddReal(res, rat, -real);
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
   if( quot < 1.0 )
      quot = 1.0;

   res->val = ((val1->val)-(val2->val))/quot;
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
}

/** multiply two rationals and save the result in res */
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
         res->isfprepresentable = TRUE;
      }
      else
      {
         // SCIPerrorMessage("multiplying with infinity might produce undesired behavior \n");
         res->val = op1->val.sign() * op2->val.sign();
         res->isinf = TRUE;
      }
   }
   else
   {
      res->val = op1->val * op2->val;
      res->isinf = FALSE;
   }
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
}

/** multiplies a rational and a real and saves the result in res */
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
    else if( REALABS(op2) >= infinity )
    {
       RatSetReal(res, op2 * op1->val.sign());
    }
    else
    {
       res->val = op1->val * op2;
       res->isinf = FALSE;
    }
    res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
}


/** divides two rationals and saves the result in res */
/** @todo exip: should we allow infinity here? */
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
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
}

/** divides a rational by a real and saves the result in res */
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
         res->isfprepresentable = SCIP_ISFPREPRESENTABLE_FALSE;
      }
   }
   else
   {
      res->isinf = FALSE;
      res->val += op1->val * op2->val;
   }
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
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
         res->isfprepresentable = SCIP_ISFPREPRESENTABLE_FALSE;
      }
   }
   else
   {
      res->isinf = FALSE;
      res->val -= op1->val * op2->val;
   }
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
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
   res->isfprepresentable = op->isfprepresentable;
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
   res->isfprepresentable = op->isfprepresentable;
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
   res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
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
      res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
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
      res->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
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

   if( rat1->isinf && rat2->isinf )
      return (rat1->val.sign() == rat2->val.sign());

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
   if( rat1->isinf && rat2->isinf )
      return TRUE;

   return FALSE;
}

/** check if a rational and a real are equal */
SCIP_Bool RatIsEqualReal(
   SCIP_Rational*        rat,                /**< the rational */
   SCIP_Real             real                /**< the real */
   )
{
   assert(rat != NULL);

   if( REALABS(real) >= infinity && rat->isinf )
      return (real > 0 && RatIsPositive(rat));

   return !rat->isinf && rat->val == real;
}

/** check if real approx of rational and a real are equal */
SCIP_Bool RatIsApproxEqualReal(
   SCIP_SET*             set,                /**< SCIP set pointer */
   SCIP_Rational*        rat,                /**< the rational */
   SCIP_Real             real                /**< the real */
   )
{
   assert(rat != NULL);

   if( rat->isinf )
   {
      return RatIsPositive(rat) ? SCIPsetIsInfinity(set, real) : SCIPsetIsInfinity(set, -real);
   }
   else
      return SCIPsetIsEQ(set, real, RatApproxReal(rat));
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
   SCIP_Rational*        rat1,               /**< The first rational */
   SCIP_Rational*        rat2                /**< The second rational */
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
SCIP_EXPORT
SCIP_Bool RatIsGE(
   SCIP_Rational*        rat1,               /**< The first rational */
   SCIP_Rational*        rat2                /**< The second rational */
   )
{
   assert(rat1 != NULL && rat2 != NULL);

   if( RatIsEqual(rat1, rat2) )
      return TRUE;
   else
      return RatIsGT(rat1, rat2);
}

/** check if the first rational is less or equal than the second*/
SCIP_EXPORT
SCIP_Bool RatIsLE(
   SCIP_Rational*        rat1,               /**< The first rational */
   SCIP_Rational*        rat2                /**< The second rational */
   )
{
   assert(rat1 != NULL && rat2 != NULL);

   return RatIsGE(rat2, rat1);
}

/** check if the rational is greater than than the double */
SCIP_Bool RatIsGTReal(
   SCIP_Rational*        rat,                /**< The rational */
   SCIP_Real             real                /**< The real */
   )
{
   assert(rat != NULL);

   if( rat->isinf )
   {
      return (rat->val > 0) && (real < infinity);
   }
   else
   {
      return rat->val > real;
   }
}

/** check if the rational is greater or equal than than the double */
SCIP_Bool RatIsGEReal(
   SCIP_Rational*        rat,                /**< The rational */
   SCIP_Real             real                /**< The real */
   )
{
   assert(rat != NULL);

   if( rat->isinf )
   {
      return RatApproxReal(rat) == real;
   }
   else
   {
      return rat->val >= real;
   }
}

/** check if the rational is less than than the double */
SCIP_Bool RatIsLTReal(
   SCIP_Rational*        rat,                /**< The rational */
   SCIP_Real             real                /**< The real */
   )
{
   assert(rat != NULL);

   if( rat->isinf )
   {
      return (rat->val < 0) && (real > -infinity);
   }
   else
   {
      return rat->val < real;
   }
}

/** check if the rational is less or equal than than the double */
SCIP_Bool RatIsLEReal(
   SCIP_Rational*        rat,                /**< The rational */
   SCIP_Real             real                /**< The real */
   )
{
   assert(rat != NULL);

   if( rat->isinf )
   {
      return RatApproxReal(rat) == real;
   }
   else
   {
      return rat->val <= real;
   }
}



/** check if the rational is zero */
SCIP_Bool RatIsZero(
   SCIP_Rational*        rational            /**< the rational to check */
   )
{
   assert(rational != NULL);

   if( rational->val.is_zero() )
   {
      assert(!rational->isinf);
      return TRUE;
   }
   else
      return FALSE;
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
   SCIP_Rational*        rational            /**< the rational to check */
   )
{
   assert(rational != NULL);

   return rational->isinf && rational->val.sign() > 0;
}

/** check if the rational is negative infinity */
SCIP_Bool RatIsNegInfinity(
   SCIP_Rational*        rational            /**< the rational to check */
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
   SCIP_Rational*        rational            /**< the rational to check */
   )
{
   assert(rational != NULL);
   if( rational->isinf )
      return FALSE;
   else if( denominator(rational->val) == 1 )
      return TRUE;
   else if( numerator(rational->val) < denominator(rational->val) )
      return FALSE;
   else
   {
      RatCanonicalize(rational);
      return (denominator(rational->val) == 1);
   }
}

/** check if rational is exactly representable as double */
SCIP_Bool RatIsFpRepresentable(
   SCIP_Rational*        rational            /**< the rational to check */
   )
{
   assert(rational != NULL);
   if( rational->isfprepresentable == SCIP_ISFPREPRESENTABLE_TRUE )
   {
      assert(RatRoundReal(rational, SCIP_ROUND_DOWNWARDS) == RatRoundReal(rational, SCIP_ROUND_UPWARDS));
      return TRUE;
   }
   else if( rational->isfprepresentable == SCIP_ISFPREPRESENTABLE_FALSE )
   {
      return FALSE;
   }
   else
   {
      rational->isfprepresentable = (RatRoundReal(rational, SCIP_ROUND_DOWNWARDS)
         == RatRoundReal(rational, SCIP_ROUND_UPWARDS)) ? SCIP_ISFPREPRESENTABLE_TRUE : SCIP_ISFPREPRESENTABLE_FALSE;
   }

   return rational->isfprepresentable == SCIP_ISFPREPRESENTABLE_TRUE ? TRUE : FALSE;
}

/*
 * Printing/Conversion methods
 */

/** converts a rational to a string for printing, returns the number of copied characters.
 *
 *  @note If return value is equal to strlen, it means the string was truncated.
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
      if( rational->val.sign() > 0 )
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
      RatDebugMessage("Rational string to long to fit in buffer. Rational : %q \n", rational);
   }

   return ret;
}

/** returns the strlen of a rational number */
SCIP_Longint RatStrlen(
   SCIP_Rational*        rational            /** rational to consider */
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

/** prints rational to file using message handler */
void RatMessage(
   SCIP_MESSAGEHDLR*     msg,                /**< message handler */
   FILE*                 file,               /**< file pointer */
   SCIP_Rational*        rational            /**< the rational to print */
   )
{
   char buf[SCIP_MAXSTRLEN];
   assert(rational != NULL);

   if( SCIP_MAXSTRLEN == RatToString(rational, buf, SCIP_MAXSTRLEN) )
   {
      SCIPerrorMessage("WARNING: Rational does not fit in line \n.");
   }

   SCIPmessageFPrintInfo(msg, file, "%s", buf);
}

/** print a rational to command line (for debugging) */
void RatPrint(
   SCIP_Rational*        rational            /**< the rational to print */
   )
{
   assert(rational != NULL);
   if( rational->isinf )
      std::cout << rational->val.sign() << "inf" << "\n";
   else
      std::cout << rational->val << "\n";
}


std::ostream& operator<<(std::ostream& os, SCIP_Rational const & r) {
   if( r.isinf )
      os << r.val.sign() << "inf";
   else
      os << r.val;
}

/* convert va_arg format string into std:string */
static
std::string RatString(const char *format, va_list arguments)
{
   std::ostringstream stream;
   SCIP_Rational* rat;
   char* sval;
   SCIP_Real dval;
   int ival;
   char cval;

   while( *format != '\0' )
   {
      if(*format == '%' && *(format+1) != '%')
      {
         switch (*++format)
         {
         case 'q':
            rat = va_arg(arguments, SCIP_Rational*);
            stream << *rat;
            break;
         case 's':
            for (sval = va_arg(arguments, char *); *sval; sval++)
               stream << (*sval);
            break;
         case 'f':
            dval = va_arg(arguments, SCIP_Real);
            stream << boost::format("%f") % dval;
         case 'g':
            dval = va_arg(arguments, SCIP_Real);
            stream << boost::format("%g") % dval;
         case 'e':
            dval = va_arg(arguments, SCIP_Real);
            stream << boost::format("%e") % dval;
            break;
         case 'd':
         case 'i':
         case 'u':
            ival = va_arg(arguments, int);
            stream << ival;
            break;
         case 'c':
            cval = (char) va_arg(arguments, int);
            stream << cval;
            break;
         default:
            stream << (*format);
            break;
         }
      }
      else
      {
         stream << (*format);
      }
      ++format;
   }

   std::string ret = stream.str();
   return ret;
}

/* printf extension for rationals (not supporting all format options) */
void RatPrintf(const char *format, ...)
{
   SCIP_Rational* rat;
   char* sval;
   SCIP_Real dval;
   int ival;
   char cval;

   va_list arguments;
   va_start(arguments, format);
   while( *format != '\0' )
   {
      if(*format == '%' && *(format+1) != '%')
      {
         switch (*++format)
         {
         case 'q':
            rat = va_arg(arguments, SCIP_Rational*);
            RatPrint(rat);
            break;
         case 's':
            for (sval = va_arg(arguments, char *); *sval; sval++)
               putchar(*sval);
            break;
         case 'f':
            dval = va_arg(arguments, SCIP_Real);
            printf("%f", dval);
            break;
         case 'g':
            dval = va_arg(arguments, SCIP_Real);
            printf("%g", dval);
            break;
         case 'e':
            dval = va_arg(arguments, SCIP_Real);
            printf("%e", dval);
            break;
         case 'd':
         case 'i':
            ival = va_arg(arguments, int);
            printf("%d", ival);
            break;
         case 'u':
            ival = va_arg(arguments, int);
            printf("%u", ival);
            break;
         case 'c':
            cval = (char) va_arg(arguments, int);
            printf("%c", cval);
            break;
         default:
            putchar(*format);
            break;
         }
      }
      else
      {
         putchar(*format);
      }
      ++format;
   }

   va_end(arguments);
}


/** get the relaxation of a rational as a real, unfortunately you can't control the roundmode without using mpfr */
/** @todo exip: we might have to worry about incorrect results when huge coefficients occur */
SCIP_Real RatRoundReal(
   SCIP_Rational*        rational,           /**< the rational */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   )
{
   SCIP_Real realapprox;
   SCIP_Longint nom, denom;
   SCIP_ROUNDMODE current;

   assert(rational != NULL);

   if( rational->isinf )
      return (rational->val.sign() * infinity);
   if( rational->isfprepresentable == SCIP_ISFPREPRESENTABLE_TRUE || roundmode == SCIP_ROUND_NEAREST )
      return RatApproxReal(rational);

#ifdef FALSE
   {
      mpfr_t valmpfr;
      mpq_t* val;

      // RatCanonicalize(rational);

      val = RatGetGMP(rational);
      switch(roundmode)
      {
         case SCIP_ROUND_DOWNWARDS:
            mpfr_init_set_q(valmpfr, *val, MPFR_RNDD);
            realapprox = (SCIP_Real) mpfr_get_d(valmpfr, MPFR_RNDD);
            break;
         case SCIP_ROUND_UPWARDS:
            mpfr_init_set_q(valmpfr, *val, MPFR_RNDU);
            realapprox = (SCIP_Real) mpfr_get_d(valmpfr, MPFR_RNDU);
            break;
         case SCIP_ROUND_NEAREST:
            mpfr_init_set_q(valmpfr, *val, MPFR_RNDN);
            realapprox = (SCIP_Real) mpfr_get_d(valmpfr, MPFR_RNDN);
            break;
         default:
            break;
      }
      mpfr_clear(valmpfr);
   }
#else

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

   SCIPdebugMessage("computing %lld/%lld \n", nom, denom);

   realapprox = ((double) nom) / denom;
   //realapprox = RgetRealApprox(r);

   if( current != roundmode )
      SCIPintervalSetRoundingMode(current);
  #endif

   SCIPdebugMessage("%.*e , Roundmode %d \n",__DBL_DECIMAL_DIG__, realapprox, roundmode );

   return realapprox;
}

/** round a rational to the nearest integer and save it as a rational */
void RatRound(
   SCIP_Rational*        res,                /**< the resulting rounded integer */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   )
{
#ifdef SCIP_WITH_BOOST
   Integer roundint, rest;

   assert(src != NULL);
   assert(res != NULL);

   if( src->isinf )
      RatSet(res, src);
   else
   {
      roundint = 0;
      rest = 0;
      divide_qr(numerator(src->val), denominator(src->val), roundint, rest);
      if( rest != 0 )
      {
         switch (roundmode)
         {
         case SCIP_ROUND_DOWNWARDS:
            roundint = src->val.sign() > 0 ? roundint : roundint - 1;
            break;
         case SCIP_ROUND_UPWARDS:
            roundint = src->val.sign() > 0 ? roundint + 1 : roundint;
            break;
         case SCIP_ROUND_NEAREST:
            roundint = abs(rest) * 2 >= denominator(src->val) ? roundint + src->val.sign() : roundint;
            break;
         default:
            SCIPerrorMessage("roundmode not supported for integer-rounding \n");
            SCIPABORT();
            break;
         }
      }
      res->val = roundint;
   }
#endif
}

/** round rational to next integer in direction of roundmode, return FALSE 
 * if rational outside of long-range
 */
SCIP_Bool RatRoundInteger(
   SCIP_Longint*         res,                /**< the resulting rounded long int */
   SCIP_Rational*        src,                /**< the rational to round */
   SCIP_ROUNDMODE        roundmode           /**< the rounding direction */
   )
{
   SCIP_Bool success = FALSE;
#ifdef SCIP_WITH_BOOST
   Integer roundint, rest;

   assert(src != NULL);
   assert(res != NULL);
   assert(!src->isinf);

   divide_qr(numerator(src->val), denominator(src->val), roundint, rest);

   if( rest != 0 )
   {
      switch (roundmode)
      {
      case SCIP_ROUND_DOWNWARDS:
         roundint = src->val.sign() > 0 ? roundint : roundint - 1;
         break;
      case SCIP_ROUND_UPWARDS:
         roundint = src->val.sign() > 0 ? roundint + 1 : roundint;
         break;
      case SCIP_ROUND_NEAREST:
         roundint = abs(rest) * 2 >= denominator(src->val) ? roundint + src->val.sign() : roundint;
         break;
      default:
         SCIPerrorMessage("roundmode not supported for integer-rounding \n");
         SCIPABORT();
         break;
      }
   }
   *res = roundint.convert_to<SCIP_Longint>();
   if( *res == roundint )
      success = TRUE;
#endif
   return success;
}

/** get the relaxation of a rational as a real, unfortunately you can't control the roundmode without using mpfr */
SCIP_Real RatApproxReal(
   SCIP_Rational*        rational            /**< the rational */
   )
{
   SCIP_Real retval;
#ifdef SCIP_WITH_BOOST
   assert(rational != NULL);

   if( rational->isinf )
      return (rational->val.sign() * infinity);

   retval = rational->val.convert_to<SCIP_Real>();
#endif
   return retval;
}

/*
 * Vector arithmetic (to shorten code and provide benefits due to
 * usage of expression Templates)
 */
void RatScalarProduct(
   SCIP_Rational*        result,             /**< the resulting rational */
   SCIP_Rational**       array1,             /**< the first array */
   SCIP_Rational**       array2,             /**< the second array */
   int                   len                 /**< length of the arrays */
   )
{
   int i;
   result->isinf = FALSE;
   result->isfprepresentable = SCIP_ISFPREPRESENTABLE_UNKNOWN;
   Rational& rat = result->val;
   for( i = 0; i < len; ++i )
   {
      assert(array1[i] != NULL);
      assert(array2[i] != NULL);
      assert(!(array1[i]->isinf && array2[i]->isinf));
      rat += array1[i]->val * array2[i]->val;
   }
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
   new (&(*rationalarray)->vals) sparsevec();
   (*rationalarray)->firstidx = -1;

   return SCIP_OKAY;
}

/** creates a dynamic array of real values */
SCIP_RETCODE SCIPrationalarrayResize(
   SCIP_RATIONALARRAY*   rationalarray,      /**< pointer to store the real array */
   int                   newsize             /**< new size */
   )
{
   assert(rationalarray != NULL);

   rationalarray->vals.resize(newsize);

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of real values */
SCIP_RETCODE SCIPrationalarrayCopy(
   SCIP_RATIONALARRAY**  rationalarray,      /**< pointer to store the copied real array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_RATIONALARRAY*   sourcerationalarray /**< dynamic rational array to copy */
   )
{
   assert(rationalarray != NULL);
   assert(sourcerationalarray != NULL);

   SCIP_CALL( SCIPrationalarrayCreate(rationalarray, blkmem) );
   (*rationalarray)->vals = sourcerationalarray->vals;
   (*rationalarray)->firstidx = sourcerationalarray->firstidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
SCIP_RETCODE SCIPrationalarrayFree(
   SCIP_RATIONALARRAY**  rationalarray,      /**< pointer to the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(rationalarray != NULL);
   assert(*rationalarray != NULL);

   (*rationalarray)->vals.~sparsevec();
   BMSfreeBlockMemory(blkmem, rationalarray);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
void SCIPrationalarrayGetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic rational array */
   int                   idx,                /**< array index to get value for */
   SCIP_Rational*        result              /**< store the result */
   )
{
   assert(rationalarray != NULL);
   assert(idx >= 0);
   if( rationalarray->firstidx == -1 || idx < rationalarray->firstidx
      || idx >= rationalarray->vals.size() + rationalarray->firstidx )
      RatSetInt(result, 0, 1);
   else
      RatSet(result, &rationalarray->vals[idx - rationalarray->firstidx]);
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPrationalarraySetVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic rational array */
   int                   idx,                /**< array index to set value for */
   SCIP_Rational*        val                 /**< value to set array index to */
   )
{
   assert(rationalarray != NULL);
   assert(idx >= 0);

   if( rationalarray-> firstidx == -1 )
   {
      rationalarray->vals.push_back(*val);
      rationalarray->firstidx = idx;
   }
   if( idx < rationalarray->firstidx )
   {
      int ninserts = rationalarray->firstidx - idx;
      SCIP_Rational r;
      rationalarray->vals.insert(rationalarray->vals.begin(), ninserts, r);
      rationalarray->firstidx = idx;
      rationalarray->vals[0] = *val;
   }
   else if( idx >= rationalarray->vals.size() + rationalarray->firstidx )
   {
      int ninserts = idx - rationalarray->vals.size() - rationalarray->firstidx + 1;
      SCIP_Rational r;
      rationalarray->vals.insert(rationalarray->vals.end(), ninserts, r);
      rationalarray->vals[rationalarray->vals.size() - 1] = *val;
   }
   else
   {
      rationalarray->vals[idx - rationalarray->firstidx] = *val;
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPrationalarrayIncVal(
   SCIP_RATIONALARRAY*   rationalarray,      /**< dynamic rational array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Rational*        incval              /**< value to increase array index */
   )
{
   assert(incval != NULL);
   assert(!incval->isinf);

   if( RatIsZero(incval) )
      return SCIP_OKAY;
   else if( idx < rationalarray->firstidx || idx >= rationalarray->vals.size() + rationalarray->firstidx )
      SCIP_CALL( SCIPrationalarraySetVal(rationalarray, idx, incval) );
   else
   {
      rationalarray->vals[idx - rationalarray->firstidx].val += incval->val;
      rationalarray->vals[idx - rationalarray->firstidx].isfprepresentable = FALSE;
   }

   return SCIP_OKAY;
}

/** prints a rationalarray to std out */
SCIP_RETCODE SCIPrationalArrayPrint(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic rational array */
   )
{
   printf("Array with firstidx %d, length %d \n", rationalarray->firstidx, (int) rationalarray->vals.size());
   for( auto val : rationalarray->vals )
   {
      RatPrint(&val);
   }
   printf("\n");

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPrationalarrayGetMinIdx(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic rational array */
   )
{
   assert(rationalarray != NULL);

   return rationalarray->firstidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPrationalarrayGetMaxIdx(
   SCIP_RATIONALARRAY*   rationalarray       /**< dynamic rational array */
   )
{
   assert(rationalarray != NULL);

   return rationalarray->firstidx + rationalarray->vals.size() - 1;
}

}

/** set the infinity threshold to new value */
void RatSetInfinity(
   SCIP_Real             inf                 /**< new infinity value */
   )
{
   assert(inf > 0);
   infinity = inf;
}

/** return the infinity threshold for rationals */
SCIP_Real RatGetInfinity(
   void
   )
{
   return infinity;
}
