/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: intervalarith.c,v 1.63 2011/03/03 14:07:18 bzfviger Exp $"

/**@file   intervalarith.c
 * @brief  interval arithmetics for provable bounds
 * @author Tobias Achterberg
 * @author Stefan Vigerske
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/intervalarith.h"

#ifdef ROUNDING_FE
#define ROUNDING
/*
 * Linux rounding operations
 */

#include <fenv.h>

/** Linux rounding mode settings */
#define SCIP_ROUND_DOWNWARDS FE_DOWNWARD     /**< round always down */
#define SCIP_ROUND_UPWARDS   FE_UPWARD       /**< round always up */
#define SCIP_ROUND_NEAREST   FE_TONEAREST    /**< round always to nearest */
#define SCIP_ROUND_ZERO      FE_TOWARDZERO   /**< round always towards 0.0 */

/** returns whether rounding mode control is available */
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   )
{
   return TRUE;
}

/** sets rounding mode of floating point operations */
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE   roundmode           /**< rounding mode to activate */
   )
{
   if( fesetround(roundmode) != 0 )
   {
      SCIPerrorMessage("error setting rounding mode to %d\n", roundmode);
      SCIPABORT();
   }
}

/** gets current rounding mode of floating point operations */
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   )
{
   return (SCIP_ROUNDMODE)fegetround();
}

#endif



#ifdef ROUNDING_FP
#define ROUNDING
/*
 * OSF rounding operations
 */

#include <float.h>

/** OSF rounding mode settings */
#define SCIP_ROUND_DOWNWARDS FP_RND_RM       /**< round always down */
#define SCIP_ROUND_UPWARDS   FP_RND_RP       /**< round always up */
#define SCIP_ROUND_NEAREST   FP_RND_RN       /**< round always to nearest */
#define SCIP_ROUND_ZERO      FP_RND_RZ       /**< round always towards 0.0 */

/** returns whether rounding mode control is available */
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   )
{
   return TRUE;
}

/** sets rounding mode of floating point operations */
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE   roundmode           /**< rounding mode to activate */
   )
{
   if( write_rnd(roundmode) != 0 )
   {
      SCIPerrorMessage("error setting rounding mode to %d\n", roundmode);
      SCIPABORT();
   }
}

/** gets current rounding mode of floating point operations */
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   )
{
   return read_rnd();
}

#endif



#ifdef ROUNDING_MS
#define ROUNDING
/*
 * Microsoft compiler rounding operations
 */

#include <float.h>

/** Microsoft rounding mode settings */
#define SCIP_ROUND_DOWNWARDS RC_DOWN         /**< round always down */
#define SCIP_ROUND_UPWARDS   RC_UP           /**< round always up */
#define SCIP_ROUND_NEAREST   RC_NEAREST      /**< round always to nearest */
#define SCIP_ROUND_ZERO      RC_TRUNCATE     /**< round always towards zero */

/** returns whether rounding mode control is available */
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   )
{
   return TRUE;
}

/** sets rounding mode of floating point operations */
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE   roundmode           /**< rounding mode to activate */
   )
{
   if( (_controlfp(roundmode, _MCW_RC) & _MCW_RC) != roundmode )
   {
      SCIPerrorMessage("error setting rounding mode to %x\n", roundmode);
      SCIPABORT();
   }
}

/** gets current rounding mode of floating point operations */
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   )
{
	return _controlfp(0, 0) & _MCW_RC;
}
#endif



#ifndef ROUNDING
/*
 * rouding operations not available
 */
#define SCIP_ROUND_DOWNWARDS 0               /**< round always down */
#define SCIP_ROUND_UPWARDS   1               /**< round always up */
#define SCIP_ROUND_NEAREST   2               /**< round always to nearest */
#define SCIP_ROUND_ZERO      3               /**< round always towards zero */

/** returns whether rounding mode control is available */
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   )
{
   return FALSE;
}

/** sets rounding mode of floating point operations */
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE   roundmode           /**< rounding mode to activate */
   )
{  /*lint --e{715}*/
   SCIPwarningMessage("setting rounding mode not available - interval arithmetic is invalid!\n");
}

/** gets current rounding mode of floating point operations */
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   )
{
   return SCIP_ROUND_NEAREST;
}
#else
#undef ROUNDING
#endif


#if defined(__GNUC__)  /* gcc or icc compiler */

/** gets the negation of a double
 * Do this in a way that the compiler does not "optimize" it away, which usually does not considers rounding modes.
 * However, compiling with -frounding-math would allow to return -x here.
 */
static
double negate(
   /* we explicitely use double here, since I'm not sure the assembler code would work as it for other float's */
   double                x                   /**< number that should be negated */
   )
{
   /* The following line of code is taken from GAOL, http://sourceforge.net/projects/gaol. */
   asm volatile ("fldl %1; fchs; fstpl %0" : "=m" (x) : "m" (x));
   return x;
}

#elif defined(_MSC_VER)  /* cl or icl compiler */

/** gets the negation of a double
 * Do this in a way that the compiler does not "optimize" it away, which usually does not considers rounding modes.
 */
static
double negate(
   /* we explicitely use double here, since I'm not sure the assembler code would work as it for other float's */
   double                x                   /**< number that should be negated */
   )
{
   /* The following lines of code are taken from GAOL, http://sourceforge.net/projects/gaol. */
   __asm {
     fld x
     fchs
     fstp x
   }
   return x;
}

#else /* unknown compiler */

/** gets the negation of a double
 * Do this in a way that the compiler does not "optimize" it away, which usually does not considers rounding modes.
 */
static
SCIP_Real negate(
   SCIP_Real             x                   /**< number that should be negated */
   )
{
   SCIPwarningMessage("setting rounding mode not available - interval arithmetic is invalid!\n");
   return -x;
}

#endif

/* on systems where the function nextafter is not defined, we provide an implementation from Sun */
#ifdef NO_NEXTAFTER
/* The following implementation of the routine nextafter() comes with the following license:
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x

static
double nextafter(double x, double y)
{
   int   hx,hy,ix,iy;
   unsigned lx,ly;

   hx = __HI(x);     /* high word of x */
   lx = __LO(x);     /* low  word of x */
   hy = __HI(y);     /* high word of y */
   ly = __LO(y);     /* low  word of y */
   ix = hx&0x7fffffff;     /* |x| */
   iy = hy&0x7fffffff;     /* |y| */

   if(((ix>=0x7ff00000)&&((ix-0x7ff00000)|lx)!=0) ||   /* x is nan */
      ((iy>=0x7ff00000)&&((iy-0x7ff00000)|ly)!=0))     /* y is nan */
      return x+y;
   if(x==y) return x;      /* x=y, return x */
   if((ix|lx)==0) {        /* x == 0 */
       __HI(x) = hy&0x80000000;  /* return +-minsubnormal */
       __LO(x) = 1;
       y = x*x;
       if(y==x) return y; else return x;  /* raise underflow flag */
   }
   if(hx>=0) {          /* x > 0 */
       if(hx>hy||((hx==hy)&&(lx>ly))) {   /* x > y, x -= ulp */
      if(lx==0) hx -= 1;
      lx -= 1;
       } else {            /* x < y, x += ulp */
      lx += 1;
      if(lx==0) hx += 1;
       }
   } else {          /* x < 0 */
       if(hy>=0||hx>hy||((hx==hy)&&(lx>ly))){/* x < y, x -= ulp */
      if(lx==0) hx -= 1;
      lx -= 1;
       } else {            /* x > y, x += ulp */
      lx += 1;
      if(lx==0) hx += 1;
       }
   }
   hy = hx&0x7ff00000;
   if(hy>=0x7ff00000) return x+x;   /* overflow  */
   if(hy<0x00100000) {     /* underflow */
       y = x*x;
       if(y!=x) {    /* raise underflow flag */
      __HI(y) = hx; __LO(y) = lx;
      return y;
       }
   }
   __HI(x) = hx; __LO(x) = lx;
   return x;
}
#endif

/** sets rounding mode of floating point operations to downwards rounding */
void SCIPintervalSetRoundingModeDownwards(
   void
   )
{
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
}

/** sets rounding mode of floating point operations to upwards rounding */
void SCIPintervalSetRoundingModeUpwards(
   void
   )
{
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
}

/** sets rounding mode of floating point operations to nearest rounding */
void SCIPintervalSetRoundingModeToNearest(
   void
   )
{
   SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
}

/** sets rounding mode of floating point operations to towards zero rounding */
void SCIPintervalSetRoundingModeTowardsZero(
   void
   )
{
   SCIPintervalSetRoundingMode(SCIP_ROUND_ZERO);
}

/** negates a number in a way that the compiler does not optimize it away */
SCIP_Real SCIPintervalNegateReal(
   SCIP_Real             x                   /**< number to negate */
   )
{
   return negate((double)x);
}

/*
 * Interval arithmetic operations
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPintervalGetInf
#undef SCIPintervalGetSup
#undef SCIPintervalSet
#undef SCIPintervalSetBounds
#undef SCIPintervalSetEmpty
#undef SCIPintervalIsEmpty
#undef SCIPintervalSetEntire
#undef SCIPintervalIsEntire
#undef SCIPintervalIsPositiveInfinity
#undef SCIPintervalIsNegativeInfinity

/** returns infimum of interval */
SCIP_Real SCIPintervalGetInf(
   SCIP_INTERVAL         interval            /**< interval */
   )
{
   return interval.inf;
}

/** returns supremum of interval */
SCIP_Real SCIPintervalGetSup(
   SCIP_INTERVAL         interval            /**< interval */
   )
{
   return interval.sup;
}

/** stores given value as interval */
void SCIPintervalSet(
   SCIP_INTERVAL*        resultant,          /**< interval to store value into */
   SCIP_Real             value               /**< value to store */
   )
{
   assert(resultant != NULL);

   resultant->inf = value;
   resultant->sup = value;
}

/** stores given infimum and supremum as interval */
void SCIPintervalSetBounds(
   SCIP_INTERVAL*        resultant,          /**< interval to store value into */
   SCIP_Real             inf,                /**< value to store as infimum */
   SCIP_Real             sup                 /**< value to store as supremum */
   )
{
   assert(resultant != NULL);
   assert(inf <= sup);

   resultant->inf = inf;
   resultant->sup = sup;
}

/** sets interval to empty interval, which will be [infinity, -infinity] */
void SCIPintervalSetEmpty(
   SCIP_INTERVAL*        resultant           /**< resultant interval of operation */
   )
{
   assert(resultant != NULL);
   
   resultant->inf =  1.0;
   resultant->sup = -1.0;
}

/** indicates whether interval is empty, i.e., whether inf > sup */
SCIP_Bool SCIPintervalIsEmpty(
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   return operand.sup < operand.inf;
}

/** sets interval to entire [-infinity, +infinity] */
void SCIPintervalSetEntire(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant           /**< resultant interval of operation */
   )
{
   assert(resultant != NULL);

   resultant->inf = -infinity;
   resultant->sup =  infinity;
}

/** indicates whether interval is entire, i.e., whether inf <= -infinity and sup >= infinity */
SCIP_Bool SCIPintervalIsEntire(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   return operand.inf <= -infinity && operand.sup >= infinity;
}

/** indicates whether interval is positive infinity, i.e., [infinity, infinity] */
SCIP_Bool SCIPintervalIsPositiveInfinity(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   return operand.inf >=  infinity && operand.sup >= operand.inf;
}

/** indicates whether interval is negative infinity, i.e., [-infinity, -infinity] */
SCIP_Bool SCIPintervalIsNegativeInfinity(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   return operand.sup <= -infinity && operand.inf <= operand.sup;
}

/** indicates whether operand1 is contained in operand2 */
SCIP_Bool SCIPintervalIsSubsetEQ(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   /* the empty interval is contained everywhere */
   if( operand1.inf > operand1.sup )
      return TRUE;
   
   /* something not-empty is not contained in the empty interval */
   if( operand2.inf > operand2.sup )
      return FALSE;
   
   return (MAX(-infinity, operand1.inf) >= operand2.inf) &&
          (MIN( infinity, operand1.sup) <= operand2.sup);
}

/** indicates whether operand1 and operand2 are disjoint */
SCIP_Bool SCIPintervalAreDisjoint(
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   return (operand1.sup < operand2.inf) || (operand2.sup < operand1.inf);
}

/** intersection of two intervals */
void SCIPintervalIntersect(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(resultant != NULL);
   
   resultant->inf = MAX(operand1.inf, operand2.inf);
   resultant->sup = MIN(operand1.sup, operand2.sup);
}

/** interval enclosure of the union of two intervals */
void SCIPintervalUnify(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(resultant != NULL);
   
   if( operand1.inf > operand1.sup )
   {
      /* operand1 is empty */
      *resultant = operand2;
      return;
   }

   if( operand2.inf > operand2.sup )
   {
      /* operand2 is empty */
      *resultant = operand1;
      return;
   }

   resultant->inf = MIN(operand1.inf, operand2.inf);
   resultant->sup = MAX(operand1.sup, operand2.sup);
}

/** adds operand1 and operand2 and stores infimum of result in infimum of resultant */
void SCIPintervalAddInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_DOWNWARDS);
   assert(resultant != NULL);

   /* [a,...] + [-inf,...] = [-inf,...] for all a, in particular, [+inf,...] + [-inf,...] = [-inf,...] */
   if( operand1.inf <= -infinity || operand2.inf <= -infinity )
   {
      resultant->inf = -infinity;
   }
   /* [a,...] + [+inf,...] = [+inf,...] for all a > -inf */
   else if( operand1.inf >= infinity || operand2.inf >= infinity )
   {
      resultant->inf = infinity;
   }
   else
   {
      resultant->inf = operand1.inf + operand2.inf;
   }
}

/** adds operand1 and operand2 and stores supremum of result in supremum of resultant */
void SCIPintervalAddSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_UPWARDS);
   assert(resultant != NULL);

   /* [...,b] + [...,+inf] = [...,+inf] for all b, in particular, [...,-inf] + [...,+inf] = [...,+inf] */
   if( operand1.sup >= infinity || operand2.sup >= infinity )
   {
      resultant->sup = infinity;
   }
   /* [...,b] + [...,-inf] = [...,-inf] for all b < +inf */
   else if( operand1.sup <= -infinity || operand2.sup <= -infinity )
   {
      resultant->sup = -infinity;
   }
   else
   {
      resultant->sup = operand1.sup + operand2.sup;
   }
}

/** adds operand1 and operand2 and stores result in resultant */
void SCIPintervalAdd(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   roundmode = SCIPintervalGetRoundingMode();

   /* compute infimum of result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   SCIPintervalAddInf(infinity, resultant, operand1, operand2);

   /* compute supremum of result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   SCIPintervalAddSup(infinity, resultant, operand1, operand2);

   SCIPintervalSetRoundingMode(roundmode);
}

/** adds operand1 and scalar operand2 and stores result in resultant */
void SCIPintervalAddScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);

   roundmode = SCIPintervalGetRoundingMode();

   /* -inf + something >= -inf */
   if( operand1.inf <= -infinity || operand2 <= -infinity )
   {
      resultant->inf = -infinity;
   }
   else if( operand1.inf >= infinity || operand2 >= infinity )
   {
      /* inf + finite = inf, inf + inf = inf */
      resultant->inf = infinity;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
      resultant->inf = operand1.inf + operand2;
   }
   
   /* inf + something <= inf */
   if( operand1.sup >=  infinity || operand2 >= infinity )
   {
      resultant->sup =  infinity;
   }
   else if( operand1.sup <= -infinity || operand2 <= -infinity )
   {
      /* -inf + finite = -inf, -inf + (-inf) = -inf */
      resultant->sup = -infinity;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
      resultant->sup = operand1.sup + operand2;
   }
   
   SCIPintervalSetRoundingMode(roundmode);
}

/** adds vector operand1 and vector operand2 and stores result in vector resultant */
void SCIPintervalAddVectors(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< array of resultant intervals of operation */
   int                   length,             /**< length of arrays */
   SCIP_INTERVAL*        operand1,           /**< array of first operands of operation */
   SCIP_INTERVAL*        operand2            /**< array of second operands of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   int i;

   roundmode = SCIPintervalGetRoundingMode();

   /* compute infimums of resultant array */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   for( i = 0; i < length; ++i )
   {
      SCIPintervalAddInf(infinity, &resultant[i], operand1[i], operand2[i]);
   }
   /* compute supremums of result array */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   for( i = 0; i < length; ++i )
   {
      SCIPintervalAddSup(infinity, &resultant[i], operand1[i], operand2[i]);
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** substracts operand2 from operand1 and stores result in resultant */
void SCIPintervalSub(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   roundmode = SCIPintervalGetRoundingMode();

   if( operand1.inf <= -infinity || operand2.sup >=  infinity )
      resultant->inf = -infinity;
   /* [a,b] - [-inf,-inf] = [+inf,+inf] */
   else if( operand1.inf >= infinity || operand2.sup <= -infinity )
   {
      resultant->inf = infinity;
      resultant->sup = infinity;
      return;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
      resultant->inf = operand1.inf - operand2.sup;
   }
   
   if( operand1.sup >=  infinity || operand2.inf <= -infinity )
      resultant->sup =  infinity;
   /* [a,b] - [+inf,+inf] = [-inf,-inf] */
   else if( operand1.sup <= -infinity || operand2.inf >= infinity )
   {
      assert(resultant->inf == -infinity);  /* should be set above, since operand1.inf <= operand1.sup <= -infinity */  /*lint !e777*/
      resultant->sup = -infinity;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
      resultant->sup = operand1.sup - operand2.inf;
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** substracts scalar operand2 from operand1 and stores result in resultant */
void SCIPintervalSubScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIPintervalAddScalar(infinity, resultant, operand1, -operand2);
}

/** multiplies operand1 with operand2 and stores infimum of result in infimum of resultant */
void SCIPintervalMulInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   )
{
   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_DOWNWARDS); 

   if( operand1.inf >= infinity )
   {
      /* operand1 is infinity scalar */
      assert(operand1.sup >= infinity);
      SCIPintervalMulScalarInf(infinity, resultant, operand2, infinity);
   }
   else if( operand2.inf >= infinity )
   {
      /* operand2 is infinity scalar */
      assert(operand2.sup >=  infinity);
      SCIPintervalMulScalarInf(infinity, resultant, operand1, infinity);
   }
   else if( operand1.sup <= -infinity )
   {
      /* operand1 is -infinity scalar */
      assert(operand1.inf <= -infinity);
      SCIPintervalMulScalarInf(infinity, resultant, operand2, -infinity);
   }
   else if( operand2.sup <= -infinity )
   {
      /* operand2 is -infinity scalar */
      assert(operand2.inf <= -infinity);
      SCIPintervalMulScalarInf(infinity, resultant, operand1, -infinity);
   }
   else if( ( operand1.inf <= -infinity && operand2.sup > 0.0 ) 
      || ( operand1.sup > 0.0 && operand2.inf <= -infinity ) 
      || ( operand1.inf < 0.0 && operand2.sup >= infinity ) 
      || ( operand1.sup >= infinity && operand2.inf < 0.0 ) )
   {
      resultant->inf = -infinity;
   }
   else
   {
      SCIP_Real cand1;
      SCIP_Real cand2;
      SCIP_Real cand3;
      SCIP_Real cand4;
      
      cand1 = operand1.inf * operand2.inf;
      cand2 = operand1.inf * operand2.sup;
      cand3 = operand1.sup * operand2.inf;
      cand4 = operand1.sup * operand2.sup;
      resultant->inf = MIN(MIN(cand1, cand2), MIN(cand3, cand4));
   }
}

/** multiplies operand1 with operand2 and stores supremum of result in supremum of resultant */
void SCIPintervalMulSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   )
{
   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_UPWARDS); 

   if( operand1.inf >= infinity )
   {
      /* operand1 is infinity scalar */
      assert(operand1.sup >= infinity);
      SCIPintervalMulScalarSup(infinity, resultant, operand2, infinity);
   }
   else if( operand2.inf >= infinity )
   {
      /* operand2 is infinity scalar */
      assert(operand2.sup >=  infinity);
      SCIPintervalMulScalarSup(infinity, resultant, operand1, infinity);
   }
   else if( operand1.sup <= -infinity )
   {
      /* operand1 is -infinity scalar */
      assert(operand1.inf <= -infinity);
      SCIPintervalMulScalarSup(infinity, resultant, operand2, -infinity);
   }
   else if( operand2.sup <= -infinity )
   {
      /* operand2 is -infinity scalar */
      assert(operand2.inf <= -infinity);
      SCIPintervalMulScalarSup(infinity, resultant, operand1, -infinity);
   }
   else if( ( operand1.inf <= -infinity && operand2.inf < 0.0 ) 
      || ( operand1.inf < 0.0 && operand2.inf <= -infinity ) 
      || ( operand1.sup > 0.0 && operand2.sup >= infinity ) 
      || ( operand1.sup >= infinity && operand2.sup > 0.0 ) )
   {
      resultant->sup =  infinity;
   }
   else
   {
      SCIP_Real cand1;
      SCIP_Real cand2;
      SCIP_Real cand3;
      SCIP_Real cand4;
      
      cand1 = operand1.inf * operand2.inf;
      cand2 = operand1.inf * operand2.sup;
      cand3 = operand1.sup * operand2.inf;
      cand4 = operand1.sup * operand2.sup;
      resultant->sup = MAX(MAX(cand1, cand2), MAX(cand3, cand4));
   }
}

/** multiplies operand1 with operand2 and stores result in resultant */
void SCIPintervalMul(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   roundmode = SCIPintervalGetRoundingMode();

   /* compute infimum result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   SCIPintervalMulInf(infinity, resultant, operand1, operand2);

   /* compute supremum of result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   SCIPintervalMulSup(infinity, resultant, operand1, operand2);

   SCIPintervalSetRoundingMode(roundmode);
}

/** multiplies operand1 with scalar operand2 and stores infimum of result in infimum of resultant */
void SCIPintervalMulScalarInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation; can be +/- inf */
   )
{
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_DOWNWARDS);
   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);

   if( operand2 >= infinity )
   {
      /* result.inf defined by sign of operand1.inf */ 
      if( operand1.inf > 0 )
         resultant->inf = infinity;
      else if( operand1.inf < 0 )
         resultant->inf = -infinity;
      else
         resultant->inf = 0.0;
   }
   else if( operand2 <= -infinity )
   {
      /* result.inf defined by sign of operand1.sup */ 
      if( operand1.sup > 0 )
         resultant->inf = -infinity;
      else if( operand1.sup < 0 )
         resultant->inf = infinity;
      else
         resultant->inf = 0.0;
   }
   else if( operand2 == 0.0 )
   {
      resultant->inf = 0.0;
   }
   else if( operand2 > 0.0 )
   {
      if( operand1.inf <= -infinity )
         resultant->inf = -infinity;
      else if( operand1.inf >= infinity )
         resultant->inf =  infinity;
      else
         resultant->inf = operand1.inf * operand2;
   }
   else
   {
      if( operand1.sup >= infinity )
         resultant->inf = -infinity;
      else if( operand1.sup <= -infinity )
         resultant->inf =  infinity;
      else
         resultant->inf = operand1.sup * operand2;
   }
}

/** multiplies operand1 with scalar operand2 and stores supremum of result in supremum of resultant */
void SCIPintervalMulScalarSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation; can be +/- inf */
   )
{
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_UPWARDS);
   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);

   if( operand2 >= infinity )
   {
      /* result.sup defined by sign of operand1.sup */ 
      if( operand1.sup > 0 )
         resultant->sup = infinity;
      else if( operand1.sup < 0 )
         resultant->sup = -infinity;
      else
         resultant->sup = 0.0;
   }
   else if( operand2 <= -infinity )
   {
      /* result.sup defined by sign of operand1.inf */ 
      if( operand1.inf > 0 )
         resultant->sup = -infinity;
      else if( operand1.inf < 0 )
         resultant->sup = infinity;
      else
         resultant->sup = 0.0;
   }
   else if( operand2 == 0.0 )
   {
      resultant->sup = 0.0;
   }
   else if( operand2 > 0.0 )
   {
      if( operand1.sup >= infinity )
         resultant->sup = infinity;
      else if( operand1.sup <= -infinity )
         resultant->sup = -infinity;
      else
         resultant->sup = operand1.sup * operand2;
   }
   else
   {
      if( operand1.inf <= -infinity )
         resultant->sup = infinity;
      else if( operand1.inf >= infinity )
         resultant->sup = -infinity;
      else
         resultant->sup = operand1.inf * operand2;
   }
}

/** multiplies operand1 with scalar operand2 and stores result in resultant */
void SCIPintervalMulScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);

   roundmode = SCIPintervalGetRoundingMode();

   /* compute infimum result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   SCIPintervalMulScalarInf(infinity, resultant, operand1, operand2);

   /* compute supremum of result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   SCIPintervalMulScalarSup(infinity, resultant, operand1, operand2);

   SCIPintervalSetRoundingMode(roundmode);
}

/** divides operand1 by operand2 and stores result in resultant */
void SCIPintervalDiv(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_INTERVAL intmed;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   if( operand2.inf <= 0.0 && operand2.sup >= 0.0 )
   {  /* division by [0,0] or interval containing 0 gives [-inf, +inf] */
      resultant->inf = -infinity;
      resultant->sup =  infinity;
      return;
   }

   if( operand1.inf == 0.0 && operand1.sup == 0.0 )
   {  /* division of [0,0] by something nonzero */
      SCIPintervalSet(resultant, 0.0);
      return;
   }

   roundmode = SCIPintervalGetRoundingMode();

   /* division by nonzero: resultant = x * (1/y) */
   if( operand2.sup >=  infinity || operand2.sup <= -infinity )
   {
      intmed.inf = 0.0;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
      intmed.inf = 1.0 / operand2.sup;
   }
   if( operand2.inf <= -infinity || operand2.inf >= infinity )
   {
      intmed.sup = 0.0;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
      intmed.sup = 1.0 / operand2.inf;
   }
   SCIPintervalMul(infinity, resultant, operand1, intmed);
 
   SCIPintervalSetRoundingMode(roundmode);
}

/** divides operand1 by scalar operand2 and stores result in resultant */
void SCIPintervalDivScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);

   roundmode = SCIPintervalGetRoundingMode();
   
   if( operand2 >= infinity || operand2 <= -infinity )
   {
      /* division by +/-infinity is 0.0 */
      resultant->inf = 0.0;
      resultant->sup = 0.0;
   }
   else if( operand2 > 0.0 )
   {
      if( operand1.inf <= -infinity )
      {
         resultant->inf = -infinity;
      }
      else if( operand1.inf >= infinity )
      {
         /* infinity / + = infinity */
         resultant->inf = infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.inf / operand2;
      }
      if( operand1.sup >= infinity )
      {
         resultant->sup =  infinity;
      }
      else if( operand1.sup <= -infinity )
      {
         /* -infinity / + = -infinity */
         resultant->sup = -infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.sup / operand2;
      }
   }
   else if( operand2 < 0.0 )
   {
      if( operand1.sup >=  infinity )
      {
         resultant->inf = -infinity;
      }
      else if( operand1.sup <= -infinity )
      {
         /* -infinity / - = infinity */
         resultant->inf = infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.sup / operand2;
      }
      if( operand1.inf <= -infinity )
      {
         resultant->sup = infinity;
      }
      else if( operand1.inf >= infinity )
      {
         /* infinity / - = -infinity */
         resultant->sup = -infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.inf / operand2;
      }
   }
   else
   { /* division by 0.0 */
      if( operand1.inf >= 0 )
      {
         /* [+,+] / [0,0] = [+inf, +inf] */
         resultant->inf =  infinity;
         resultant->sup =  infinity;
      }
      else if( operand1.sup <= 0 )
      {
         /* [-,-] / [0,0] = [-inf, -inf] */
         resultant->inf = -infinity;
         resultant->sup = -infinity;
      }
      else
      {
         /* [-,+] / [0,0] = [-inf, +inf] */
         resultant->inf = -infinity;
         resultant->sup =  infinity;
      }
      return;
   }
  
   SCIPintervalSetRoundingMode(roundmode);
}

/** computes the scalar product of two vectors of intervals and stores result in resultant */
void SCIPintervalScalprod(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals; can have +/-inf entries */
   SCIP_INTERVAL*        operand2            /**< second vector as array of intervals; can have +/-inf entries */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_INTERVAL prod;
   int i;
   
   roundmode = SCIPintervalGetRoundingMode();

   resultant->inf = 0.0;
   resultant->sup = 0.0;

   /* compute infimum of resultant */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   for( i = 0; i < length && resultant->inf > -infinity; ++i )
   {
      SCIPintervalSetEntire(infinity, &prod);
      SCIPintervalMulInf(infinity, &prod, operand1[i], operand2[i]);
      SCIPintervalAddInf(infinity, resultant, *resultant, prod); 
   }
   assert(resultant->sup == 0.0);

   /* compute supremum of resultant */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   for( i = 0; i < length && resultant->sup < infinity ; ++i )
   {
      SCIPintervalSetEntire(infinity, &prod);
      SCIPintervalMulSup(infinity, &prod, operand1[i], operand2[i]);
      SCIPintervalAddSup(infinity, resultant, *resultant, prod); 
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** computes scalar product of a vector of intervals and a vector of scalars and stores infimum of result in infimum of 
 *  resultant 
 */
void SCIPintervalScalprodScalarsInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   )
{
   SCIP_INTERVAL prod;
   int i;

   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_DOWNWARDS);

   resultant->inf = 0.0;

   /* compute infimum of resultant */
   SCIPintervalSetEntire(infinity, &prod);
   for( i = 0; i < length && resultant->inf > -infinity; ++i )
   {
      SCIPintervalMulScalarInf(infinity, &prod, operand1[i], operand2[i]);
      assert(prod.sup >= infinity);
      SCIPintervalAddInf(infinity, resultant, *resultant, prod); 
   }
}

/** computes the scalar product of a vector of intervals and a vector of scalars and stores supremum of result in 
 *  supremum of resultant 
 */
void SCIPintervalScalprodScalarsSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   )
{
   SCIP_INTERVAL prod;
   int i;

   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_UPWARDS);

   resultant->sup = 0.0;

   /* compute supremum of resultant */
   SCIPintervalSetEntire(infinity, &prod);
   for( i = 0; i < length && resultant->sup < infinity; ++i )
   {
      SCIPintervalMulScalarSup(infinity, &prod, operand1[i], operand2[i]);
      assert(prod.inf <= -infinity);
      SCIPintervalAddSup(infinity, resultant, *resultant, prod); 
   }
}

/** computes the scalar product of a vector of intervals and a vector of scalars and stores result in resultant */
void SCIPintervalScalprodScalars(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   )
{
   SCIP_ROUNDMODE roundmode;

   roundmode = SCIPintervalGetRoundingMode();

   resultant->inf = 0.0;
   resultant->sup = 0.0;

   /* compute infimum of resultant */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   SCIPintervalScalprodScalarsInf(infinity, resultant, length, operand1, operand2);
   assert(resultant->sup == 0.0);

   /* compute supremum of resultant */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   SCIPintervalScalprodScalarsSup(infinity, resultant, length, operand1, operand2);

   SCIPintervalSetRoundingMode(roundmode);
}

/** squares operand and stores result in resultant */
void SCIPintervalSquare(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand.inf <= operand.sup);
  
   roundmode = SCIPintervalGetRoundingMode();

   if( operand.sup <= 0.0 )
   {  /* operand is left of 0.0 */
      if( operand.sup <= -infinity )
      {
         resultant->inf =  infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand.sup * operand.sup;
      }
      if( operand.inf <= -infinity )
      {
         resultant->sup = infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand.inf * operand.inf;
      }
   }
   else if( operand.inf >= 0.0 )
   {  /* operand is right of 0.0 */
      if( operand.inf >= infinity )
      {
         resultant->inf = infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand.inf * operand.inf;
      }
      if( operand.sup >= infinity )
      {
         resultant->sup = infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand.sup * operand.sup;
      }
   }
   else
   {  /* [-,+]^2 */
      resultant->inf = 0.0;
      if( operand.inf <= -infinity || operand.sup >= infinity )
      {
         resultant->sup = infinity;
      }
      else
      {
         SCIP_Real x;
         SCIP_Real y;

         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         x = operand.inf * operand.inf;
         y = operand.sup * operand.sup;
         resultant->sup = MAX(x, y);
      }
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** stores (positive part of) square root of operand in resultant
 * @caution we assume a correctly rounded sqrt(double) function when rounding is to nearest
 */
void SCIPintervalSquareRoot(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(operand.inf <= operand.sup);
   
   if( operand.sup < 0.0 )
   {
      SCIPintervalSetEmpty(resultant);
      return;
   }

   if( operand.inf == operand.sup )
   {
      if( operand.inf >= infinity )
      {
         resultant->inf = infinity;
         resultant->sup = infinity;
      }
      else
      {
         SCIP_Real tmp;

         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, noone should have changed rounding mode */
         tmp = sqrt(operand.inf);
         resultant->inf = nextafter(tmp, SCIP_REAL_MIN);
         resultant->sup = nextafter(tmp, SCIP_REAL_MAX);
      }

      return;
   }
  
   if( operand.inf <= 0.0 )
   {
      resultant->inf = 0.0;
   }
   else if( operand.inf >= infinity )
   {
      resultant->inf = infinity;
      resultant->sup = infinity;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, noone should have changed rounding mode */
      resultant->inf = nextafter(sqrt(operand.inf), SCIP_REAL_MIN);
   }
   
   if( operand.sup >= infinity )
   {
      resultant->sup = infinity;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, noone should have changed rounding mode */
      resultant->sup = nextafter(sqrt(operand.sup), SCIP_REAL_MAX);
   }
}

/** stores operand1 to the power of operand2 in resultant
 * 
 * uses SCIPintervalPowerScalar if operand2 is a scalar, otherwise computes exp(op2*log(op1)) */
void SCIPintervalPower(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{  /*lint --e{777}*/
   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   if( operand2.inf == operand2.sup )
   {  /* operand is number */
      SCIPintervalPowerScalar(infinity, resultant, operand1, operand2.inf);
      return;
   }
   
   /* resultant := log(op1) */
   SCIPintervalLog(infinity, resultant, operand1);
   if( SCIPintervalIsEmpty(*resultant) )
      return;
   
   /* resultant := op2 * resultant */
   SCIPintervalMul(infinity, resultant, operand2, *resultant);
   
   /* resultant := exp(resultant) */
   SCIPintervalExp(infinity, resultant, *resultant);
}

/** computes lower bound on power of a scalar operand1 to an integer operand2
 * both operands need to be finite numbers
 * need to have operand1 >= 0 and need to have operand2 >= 0 if operand1 == 0
 */
SCIP_Real SCIPintervalPowerScalarIntegerInf(
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_Real result;

   assert(operand1 >= 0.0);

   if( operand1 == 0.0 )
   {
      assert(operand2 >= 0);
      if( operand2 == 0 )
         return 1.0; /* 0^0 = 1 */
      else
         return 0.0; /* 0^positive = 0 */
   }

   /* 1^n = 1, x^0 = 1 */
   if( operand1 == 1.0 || operand2 == 0 )
      return 1.0;

   if( operand2 < 0 )
   {
      /* x^n = 1 / x^(-n) */
      result = SCIPintervalPowerScalarIntegerSup(operand1, -operand2);
      assert(result != 0.0);

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
      result = 1.0 / result;
      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {
      unsigned int n;
      SCIP_Real z;

      roundmode = SCIPintervalGetRoundingMode();

      result = 1.0;
      n = (unsigned int)operand2;
      z = operand1;

      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);

      /* use a binary exponentiation algorithm:
       * consider the binary respresentation of n: n = sum_i 2^i d_i with d_i \in {0,1}
       * then x^n = prod_{i:d_i=1} x^(2^i)
       * in the following, we loop over i=1,..., thereby storing x^(2^i) in z
       * whenever d_i is 1, we multiply result with x^(2^i) (the current z)
       * at the last (highest) i with d_i = 1 we stop, thus having x^n stored in result
       *
       * the binary representation of n and bit shifting is used for the loop
       */
      assert(n >= 1);
      do
      {
         if( n&1 ) /* n is odd (d_i=1), so multiply result with current z (=x^{2^i}) */
         {
            result = result * z;
            n >>= 1;
            if( n == 0 )
               break;
         }
         else
            n >>= 1;
         z = z * z;
      } while( TRUE );

      SCIPintervalSetRoundingMode(roundmode);
   }

   return result;
}

/** computes upper bound on power of a scalar operand1 to an integer operand2
 * both operands need to be finite numbers
 * need to have operand1 >= 0 and need to have operand2 >= 0 if operand1 == 0
 */
extern
SCIP_Real SCIPintervalPowerScalarIntegerSup(
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_Real result;

   assert(operand1 >= 0.0);

   if( operand1 == 0.0 )
   {
      assert(operand2 >= 0);
      if( operand2 == 0 )
         return 1.0; /* 0^0 = 1 */
      else
         return 0.0; /* 0^positive = 0 */
   }

   /* 1^n = 1, x^0 = 1 */
   if( operand1 == 1.0 || operand2 == 0 )
      return 1.0;

   if( operand2 < 0 )
   {
      /* x^n = 1 / x^(-n) */
      result = SCIPintervalPowerScalarIntegerInf(operand1, -operand2);
      assert(result != 0.0);

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
      result = 1.0 / result;
      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {
      unsigned int n;
      SCIP_Real z;

      roundmode = SCIPintervalGetRoundingMode();

      result = 1.0;
      n = (unsigned int)operand2;
      z = operand1;

      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);

      /* use a binary exponentiation algorithm... see comments in SCIPintervalPowerScalarIntegerInf */
      assert(n >= 1);
      do
      {
         if( n&1 )
         {
            result = result * z;
            n >>= 1;
            if( n == 0 )
               break;
         }
         else
            n >>= 1;
         z = z * z;
      } while( TRUE );

      SCIPintervalSetRoundingMode(roundmode);
   }

   return result;
}

/** computes bounds on power of a scalar operand1 to an integer operand2
 * both operands need to be finite numbers
 * need to have operand1 >= 0 and need to have operand2 >= 0 if operand1 == 0
 */
extern
void SCIPintervalPowerScalarInteger(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(operand1 >= 0.0);

   if( operand1 == 0.0 )
   {
      assert(operand2 >= 0);
      if( operand2 == 0 )
      {
         SCIPintervalSet(resultant, 1.0); /* 0^0 = 1 */
         return;
      }
      else
      {
         SCIPintervalSet(resultant, 0.0); /* 0^positive = 0 */
         return;
      }
   }

   /* 1^n = 1, x^0 = 1 */
   if( operand1 == 1.0 || operand2 == 0 )
   {
      SCIPintervalSet(resultant, 1.0);
      return;
   }

   if( operand2 < 0 )
   {
      /* x^n = 1 / x^(-n) */
      SCIPintervalPowerScalarInteger(resultant, operand1, -operand2);
      assert(resultant->inf > 0.0 || resultant->sup < 0.0);
      SCIPintervalReciprocal(SCIP_REAL_MAX, resultant, *resultant); /* value for infinity does not matter, since there should be no 0.0 in the interval, so just use something large enough */
   }
   else
   {
      unsigned int n;
      SCIP_Real z_inf;
      SCIP_Real z_sup;
      SCIP_Real result_sup;
      SCIP_Real result_inf;

      roundmode = SCIPintervalGetRoundingMode();

      result_inf = 1.0;
      result_sup = 1.0;
      z_inf = operand1;
      z_sup = operand1;
      n = (unsigned int)operand2;

      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);

      /* use a binary exponentiation algorithm... see comments in SCIPintervalPowerScalarIntegerInf
       * we compute lower and upper bounds within the same loop
       * to get correct lower bounds while rounding mode is upwards, we negate arguments */
      assert(n >= 1);
      do
      {
         if( n&1 )
         {
            result_inf = negate(negate(result_inf) * z_inf);
            result_sup = result_sup * z_sup;
            n >>= 1;
            if( n == 0 )
               break;
         }
         else
            n >>= 1;
         z_inf = negate(negate(z_inf) * z_inf);
         z_sup = z_sup * z_sup;
      } while( TRUE );

      SCIPintervalSetRoundingMode(roundmode);

      resultant->inf = result_inf;
      resultant->sup = result_sup;
      assert(resultant->inf <= resultant->sup);
   }
}

/** stores bounds on the power of a scalar operand1 to a scalar operand2 in resultant
 * both operands need to be finite numbers
 * need to have operand1 >= 0 or operand2 integer and need to have operand2 >= 0 if operand1 == 0
 * @caution we assume a correctly rounded pow(double) function when rounding is to nearest
 */
void SCIPintervalPowerScalarScalar(
   SCIP_INTERVAL*        resultant,          /**< resultant of operation */
   SCIP_Real             operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_Real result;

   assert(resultant != NULL);

   if( operand1 == 0.0 )
   {
      assert(operand2 >= 0);
      if( operand2 == 0 )
      {
         SCIPintervalSet(resultant, 1.0); /* 0^0 = 1 */
         return;
      }
      else
      {
         SCIPintervalSet(resultant, 0.0); /* 0^positive = 0 */
         return;
      }
   }

   /* 1^n = 1, x^0 = 1 */
   if( operand1 == 1.0 || operand2 == 0 )
   {
      SCIPintervalSet(resultant, 1.0);
      return;
   }

   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, noone should have changed rounding mode */
   result = pow(operand1, operand2);

   /* to get safe bounds, get the floating point numbers just below and above result */
   resultant->inf = nextafter(result, SCIP_REAL_MIN);
   resultant->sup = nextafter(result, SCIP_REAL_MAX);
}

/** stores operand1 to the power of the scalar operand2 in resultant
 * @caution we assume a correctly rounded pow(double) function when rounding is to nearest
 */
void SCIPintervalPowerScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{  /*lint --e{777}*/
   SCIP_Bool op2isint;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   
   if( operand2 == infinity )
   {
      /* 0^infinity =  0
       * +^infinity =  infinity
       * -^infinity = -infinity
       */
      if( operand1.inf < 0.0 )
         resultant->inf = -infinity;
      else
         resultant->inf = 0.0;
      if( operand1.sup > 0.0 )
         resultant->sup =  infinity;
      else
         resultant->sup = 0.0;
      return;
   }

   if( operand2 == 0.0 )
   { /* special case, since x^0 = 1 for x != 0, but 0^0 = 0 */
      if( operand1.inf == 0.0 && operand1.sup == 0.0 )
      {
         resultant->inf = 0.0;
         resultant->sup = 0.0;
      }
      else if( operand1.inf <= 0.0 || operand1.sup >= 0.0 )
      { /* 0.0 in x gives [0,1] */
         resultant->inf = 0.0;
         resultant->sup = 1.0;
      }
      else
      { /* 0.0 outside x gives [1,1] */
         resultant->inf = 1.0;
         resultant->sup = 1.0;
      }
      return;
   }

   if( operand2 == 0.0 )
   {
      /* x^0 = 1 */
      resultant->inf = 1.0;
      resultant->sup = 1.0;
      return;
   }

   if( operand2 == 1.0 )
   {
      /* x^1 = x */
      *resultant = operand1;
      return;
   }

   op2isint = (ceil(operand2) == operand2);
   
   if( !op2isint && operand1.inf < 0.0 )
   {  /* x^n with x negative not defined for n not integer*/
      operand1.inf = 0.0;
      if( operand1.sup < operand1.inf )
      {
         SCIPintervalSetEmpty(resultant);
         return;
      }
   }

   if( operand1.inf >= 0.0 )
   {  /* easy case: x^n with x>=0 */
      if( operand2 >= 0.0 )
      {
         if( operand1.inf >= infinity )
         {
            /* inf^+ = inf */
            resultant->inf = infinity;
         }
         else if( operand1.inf > 0.0 )
         {
            assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, noone should have changed rounding mode */
            resultant->inf = nextafter(pow(operand1.inf, operand2), SCIP_REAL_MIN);
         }
         else
            resultant->inf = 0.0;
         if( operand1.sup >= infinity )
         {
            resultant->sup = infinity;
         }
         else if( operand1.sup > 0.0 )
         {
            assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, noone should have changed rounding mode */
            resultant->sup = nextafter(pow(operand1.sup, operand2), SCIP_REAL_MAX);
         }
         else
            resultant->sup = 0.0;
      }
      else
      {
         if( operand1.sup >= infinity )
         {
            resultant->inf = 0.0;
         }
         else if( operand1.sup == 0.0 )
         {
            /* x^(negative even) = infinity for x->0 (from both sides),
             * but x^(negative odd) = -infinity for x->0 from left side */
            if( ceil(operand2/2) == operand2/2 )
               resultant->inf =  infinity;
            else
               resultant->inf = -infinity;
         }
         else
         {
            assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, noone should have changed rounding mode */
            resultant->inf = nextafter(pow(operand1.sup, operand2), SCIP_REAL_MIN);
         }
         if( operand1.inf == 0.0 )
         {
            /* 0^(negative) = infinity */
            resultant->sup = infinity;
         }
         else
         {
            assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, noone should have changed rounding mode */
            resultant->sup = nextafter(pow(operand1.inf, operand2), SCIP_REAL_MAX);
         }
      }
   }
   else if( operand1.sup <= 0.0 )
   {  /* more difficult case: x^n with x < 0; we now know, that n is integer */
      assert(op2isint);
      if( operand2 >= 0.0 && ceil(operand2/2) == operand2/2 )
      {
         /* x^n with n>=2 and even -> x^n is mon. decreasing for x < 0 */
         if( operand1.sup == -infinity )
         {
            /* (-inf)^n = inf */
            resultant->inf = infinity;
         }
         else
         {
            resultant->inf = SCIPintervalPowerScalarIntegerInf(-operand1.sup, (int)operand2);
         }
         if( operand1.inf <= -infinity )
         {
            resultant->sup = infinity;
         }
         else
         {
            resultant->sup = SCIPintervalPowerScalarIntegerSup(-operand1.inf, (int)operand2);
         }
      }
      else if( operand2 <= 0.0 && ceil(operand2/2) != operand2/2 )
      {
         /* x^n with n<=-1 and odd -> x^n = 1/x^(-n) is mon. decreasing for x<0 */
         if( operand1.sup == -infinity )
         {
            /* (-inf)^n = 1/(-inf)^(-n) = 1/(-inf) = 0 */
            resultant->inf = 0.0;
         }
         else if( operand1.sup == 0.0 )
         {
            /* x^n -> -infinity for x->0 from left */
            resultant->inf = -infinity;
         }
         else
         {
            resultant->inf = -SCIPintervalPowerScalarIntegerSup(-operand1.sup, (int)operand2);
         }
         if( operand1.inf <= -infinity )
         {
            /* (-inf)^n = 1/(-inf)^(-n) = 1/(-inf) = 0 */
            resultant->sup = 0.0;
         }
         else if( operand1.inf == 0.0 )
         {
            /* x^n -> infinity for x->0 from right */
            resultant->sup = infinity;
         }
         else
         {
            resultant->sup = -SCIPintervalPowerScalarIntegerInf(-operand1.inf, (int)operand2);
         }
      }
      else if( operand2 >= 0.0 )
      {
         /* x^n with n>0 and odd -> x^n is mon. increasing for x<0 */
         if( operand1.inf <= -infinity )
         {
            resultant->inf = -infinity;
         }
         else
         {
            resultant->inf = -SCIPintervalPowerScalarIntegerSup(-operand1.inf, (int)operand2);
         }
         if( operand1.sup <= -infinity )
         {
            resultant->sup = -infinity;
         }
         else
         {
            resultant->sup = -SCIPintervalPowerScalarIntegerInf(-operand1.sup, (int)operand2);
         }
      }
      else
      {
         /* x^n with n<0 and even -> x^n is mon. increasing for x<0 */
         if( operand1.inf <= -infinity )
         {
            resultant->inf = 0.0;
         }
         else if( operand1.inf == 0.0 )
         {
            /* x^n -> infinity for x->0 from both sides */
            resultant->inf = infinity;
         }
         else
         {
            resultant->inf = SCIPintervalPowerScalarIntegerSup(-operand1.inf, (int)operand2);
         }
         if( operand1.sup <= -infinity )
         {
            resultant->sup = 0.0;
         }
         else if( operand1.sup == 0.0 )
         {
            /* x^n -> infinity for x->0 from both sides */
            resultant->sup = infinity;
         }
         else
         {
            resultant->sup = SCIPintervalPowerScalarIntegerSup(-operand1.sup, (int)operand2);
         }
      }
      assert(resultant->inf <= resultant->sup);

   }
   else
   {  /* similar difficult case: x^n with x in [<0, >0], but n is integer */
      assert(op2isint); /* otherwise we had set operand1.inf == 0.0, which was handled in first case */
      if( operand2 >= 0.0 && operand2/2 == ceil(operand2/2) )
      {
         /* n even positive integer */
         resultant->inf = 0.0;
         if( operand1.inf == -infinity || operand1.sup == infinity )
         {
            resultant->sup = infinity;
         }
         else
         {
            resultant->sup = SCIPintervalPowerScalarIntegerSup(MAX(-operand1.inf, operand1.sup), (int)operand2);
         }
      }
      else if( operand2 <= 0.0 && ceil(operand2/2) != operand2/2 )
      {
         /* n even negative integer */
         resultant->sup = infinity;  /* since 0^n = infinity */
         if( operand1.inf == -infinity || operand1.sup == infinity )
         {
            resultant->inf = 0.0;
         }
         else
         {
            resultant->inf = SCIPintervalPowerScalarIntegerInf(MAX(-operand1.inf, operand1.sup), (int)operand2);
         }
      }
      else if( operand2 >= 0.0 )
      {
         /* n odd positive integer, so monotonically increasing function */
         if( operand1.inf == -infinity )
         {
            resultant->inf = -infinity;
         }
         else
         {
            resultant->inf = -SCIPintervalPowerScalarIntegerSup(-operand1.inf, (int)operand2);
         }
         if( operand1.sup == infinity )
         {
            resultant->sup = infinity;
         }
         else
         {
            resultant->sup = SCIPintervalPowerScalarIntegerSup(operand1.sup, (int)operand2);
         }
      }
      else
      {
         /* n odd negative integer:
          * x^n -> -infinity for x->0 from left
          * x^n ->  infinity for x->0 from right */
         resultant->inf = -infinity;
         resultant->sup =  infinity;
      }
   }
}

/** stores operand1 to the signed power of the scalar positive operand2 in resultant 
 * 
 * the signed power of x w.r.t. an exponent n >= 0 is given as sign(x) * abs(x)^n
 *
 * @caution we assume correctly rounded sqrt(double) and pow(double) functions when rounding is to nearest
 */
void SCIPintervalSignPowerScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   assert(resultant != NULL);

   assert(operand1.inf <= operand1.sup);
   assert(operand2     >= 0.0);

   if( operand2 == infinity )
   {
      /* 0^infinity =  0
       * +^infinity =  infinity
       *-+^infinity = -infinity
       */
      if( operand1.inf < 0.0 )
         resultant->inf = -infinity;
      else
         resultant->inf = 0.0;
      if( operand1.sup > 0.0 )
         resultant->sup =  infinity;
      else
         resultant->sup = 0.0;
      return;
   }

   if( operand2 == 0.0 )
   { /* special case, since x^0 = 1 for x != 0, but 0^0 = 0 */
      if( operand1.inf < 0.0 )
         resultant->inf = -1.0;
      else if( operand1.inf == 0.0 )
         resultant->inf =  0.0;
      else
         resultant->inf =  1.0;

      if( operand1.sup < 0.0 )
         resultant->sup = -1.0;
      else if( operand1.sup == 0.0 )
         resultant->sup =  0.0;
      else
         resultant->sup =  1.0;
      
      return;
   }

   if( operand2 == 1.0 )
   { /* easy case that should run fast */
      *resultant = operand1;
      return;
   }

   roundmode = SCIPintervalGetRoundingMode();

   if( operand2 == 2.0 )
   { /* common case where pow can easily be avoided */
      if( operand1.inf <= -infinity )
      {
         resultant->inf = -infinity;
      }
      else if( operand1.inf >= infinity )
      {
         resultant->inf =  infinity;
      }
      else if( operand1.inf > 0.0 )
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.inf * operand1.inf;
      }
      else
      {
         /* need upwards since we negate result of multiplication */
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->inf = negate(operand1.inf * operand1.inf);
      }

      if( operand1.sup >=  infinity )
      {
         resultant->sup =  infinity;
      }
      else if( operand1.sup <= -infinity )
      {
         resultant->sup = -infinity;
      }
      else if( operand1.sup > 0.0 )
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.sup * operand1.sup;
      }
      else
      {
         /* need downwards since we negate result of multiplication */
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->sup = negate(operand1.sup * operand1.sup);
      }
      assert(resultant->inf <= resultant->sup);
   }
   else if( operand2 == 0.5 )
   { /* another common case where pow can easily be avoided */
      if( operand1.inf <= -infinity )
      {
         resultant->inf = -infinity;
      }
      else if( operand1.inf >= infinity )
      {
         resultant->inf = infinity;
      }
      else if( operand1.inf >= 0.0 )
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->inf =  nextafter(sqrt( operand1.inf), SCIP_REAL_MIN);
      }
      else
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->inf = -nextafter(sqrt(-operand1.inf), SCIP_REAL_MAX);
      }

      if( operand1.sup >=  infinity )
      {
         resultant->sup =  infinity;
      }
      else if( operand1.sup <= -infinity )
      {
         resultant->sup = -infinity;
      }
      else if( operand1.sup > 0.0 )
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->sup =  nextafter(sqrt( operand1.sup), SCIP_REAL_MAX);
      }
      else
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->sup = -nextafter(sqrt(-operand1.sup), SCIP_REAL_MAX);
      }
      assert(resultant->inf <= resultant->sup);
   }
   else
   {
      if( operand1.inf <= -infinity )
      {
         resultant->inf = -infinity;
      }
      else if( operand1.inf >= infinity )
      {
         resultant->inf =  infinity;
      }
      else if( operand1.inf > 0.0 )
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->inf =  nextafter(pow( operand1.inf, operand2), SCIP_REAL_MIN);
      }
      else
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->inf = -nextafter(pow(-operand1.inf, operand2), SCIP_REAL_MAX);
      }

      if( operand1.sup >=  infinity )
      {
         resultant->sup =  infinity;
      }
      else if( operand1.sup <= -infinity )
      {
         resultant->sup = -infinity;
      }
      else if( operand1.sup > 0.0 )
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->sup =  nextafter(pow( operand1.sup, operand2), SCIP_REAL_MAX);
      }
      else
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->sup = -nextafter(pow(-operand1.sup, operand2), SCIP_REAL_MIN);
      }
   }
   
   SCIPintervalSetRoundingMode(roundmode);
}

/** computes the reciprocal of an interval
 */
void SCIPintervalReciprocal(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand.inf <= operand.sup);

   if( operand.inf == 0.0 && operand.sup == 0.0 )
   { /* 1/0 = [-inf,inf] */
      resultant->inf =  infinity;
      resultant->sup = -infinity;
      return;
   }

   roundmode = SCIPintervalGetRoundingMode();

   if( operand.inf >= 0.0 )
   {  /* 1/x with x >= 0 */
      if( operand.sup >= infinity )
      {
         resultant->inf = 0.0;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = 1.0 / operand.sup;
      }

      if( operand.inf >= infinity )
      {
         resultant->sup = 0.0;
      }
      else if( operand.inf == 0.0 )
      {
         resultant->sup = infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = 1.0 / operand.inf;
      }

      SCIPintervalSetRoundingMode(roundmode);
   }
   else if( operand.sup <= 0.0 )
   {  /* 1/x with x <= 0 */
      if( operand.sup <= -infinity )
      {
         resultant->inf = 0.0;
      }
      else if( operand.sup == 0.0 )
      {
         resultant->inf = -infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = 1.0 / operand.sup;
      }

      if( operand.inf <= -infinity )
      {
         resultant->sup = infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = 1.0 / operand.inf;
      }
      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {  /* 1/x with x in [-,+] is division by zero */
      resultant->inf = -infinity;
      resultant->sup =  infinity;
   }
}

/** stores exponential of operand in resultant
 * @caution we assume a correctly rounded exp(double) function when rounding is to nearest
 */
void SCIPintervalExp(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(operand.inf <= operand.sup);

   if( operand.sup <= -infinity )
   {
      resultant->inf = 0.0;
      resultant->sup = 0.0;
      return;
   }

   if( operand.inf >=  infinity )
   {
      resultant->inf = infinity;
      resultant->sup = infinity;
      return;
   }

   if( operand.inf == operand.sup )
   {
      SCIP_Real tmp;

      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
      tmp = exp(operand.inf);
      resultant->inf = nextafter(tmp, SCIP_REAL_MIN);
      resultant->sup = nextafter(tmp, SCIP_REAL_MAX);

      return;
   }

   if( operand.inf <= -infinity )
   {
      resultant->inf = 0.0;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
      resultant->inf = nextafter(exp(operand.inf), SCIP_REAL_MIN);
   }
  
   if( operand.sup >=  infinity )
   {
      resultant->sup = infinity;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
      resultant->sup = nextafter(exp(operand.sup), SCIP_REAL_MAX);
   }
}

/** stores natural logarithm of operand in resultant
 * @caution we assume a correctly rounded log(double) function when rounding is to nearest
 */
void SCIPintervalLog(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(operand.inf <= operand.sup);
  
   if( operand.sup < 0.0 )
   {
      SCIPintervalSetEmpty(resultant);
      return;
   }

   if( operand.inf == operand.sup )
   {
      if( operand.sup <= 0.0 )
      {
         resultant->inf = -infinity;
         resultant->sup = -infinity;
      }
      else
      {
         SCIP_Real tmp;

         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         tmp = log(operand.inf);
         resultant->inf = nextafter(tmp, SCIP_REAL_MIN);
         resultant->sup = nextafter(tmp, SCIP_REAL_MAX);
      }

      return;
   }
  
   if( operand.inf <= 0.0 )
   {
      resultant->inf = -infinity;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
      resultant->inf = nextafter(log(operand.inf), SCIP_REAL_MIN);
   }

   if( operand.sup >= infinity )
   {
      resultant->sup =  infinity;
   }
   else if( operand.sup == 0.0 )
   {
      resultant->sup = -infinity;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
      resultant->sup = nextafter(log(operand.sup), SCIP_REAL_MAX);
   }
}

/** stores minimum of operands in resultant */
void SCIPintervalMin(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   resultant->inf = MIN(operand1.inf, operand2.inf);
   resultant->sup = MIN(operand1.sup, operand2.sup);
}

/** stores maximum of operands in resultant */
void SCIPintervalMax(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   resultant->inf = MAX(operand1.inf, operand2.inf);
   resultant->sup = MAX(operand1.sup, operand2.sup);
}

/** stores absolute value of operand in resultant */
void SCIPintervalAbs(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(operand.inf <= operand.sup);
   
   if( operand.inf <= 0.0 && operand.sup >= 0.0)
   {
      resultant->inf = 0.0;
      resultant->sup = MAX(-operand.inf, operand.sup);
   }
   else if( operand.inf > 0.0 )
   {
      *resultant = operand;
   }
   else
   {
      resultant->inf = -operand.sup;
      resultant->sup = -operand.inf;
   }
}

/** stores sign of operand in resultant */
void SCIPintervalSign(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(operand.inf <= operand.sup);

   if( operand.sup < 0.0 )
   {
      resultant->inf = -1.0;
      resultant->sup = -1.0;
   }
   else if( operand.inf >= 0.0 )
   {
      resultant->inf =  1.0;
      resultant->sup =  1.0;
   }
   else
   {
      resultant->inf = -1.0;
      resultant->sup =  1.0;
   }
}

/** computes exact upper bound on \f$ a x^2 + b x \f$ for x in [xlb, xub], b an interval, and a scalar
 * 
 * Uses Algorithm 2.2 from Domes and Neumaier: Constraint propagation on quadratic constraints (2008) */
SCIP_Real SCIPintervalQuadUpperBound(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_Real             a,                  /**< coefficient of x^2 */
   SCIP_INTERVAL         b_,                 /**< coefficient of x */
   SCIP_INTERVAL         x                   /**< range of x */
   )
{
   SCIP_Real b;
   SCIP_Real u;
   
   assert(!SCIPintervalIsEmpty(x));
   assert(b_.inf <  infinity);
   assert(b_.sup > -infinity);
   assert( x.inf <  infinity);
   assert( x.sup > -infinity);

   /* handle b*x separately */
   if( a == 0.0 )
   {
      if( (b_.inf <= -infinity && x.inf <   0.0     ) ||
          (b_.inf <   0.0      && x.inf <= -infinity) ||
          (b_.sup >   0.0      && x.sup >=  infinity) ||
          (b_.sup >=  infinity && x.sup >   0.0     ) )
      {
         u = infinity;
      }
      else
      {
         SCIP_ROUNDMODE roundmode;
         SCIP_Real cand1, cand2, cand3, cand4;

         roundmode = SCIPintervalGetRoundingMode();
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         cand1 = b_.inf * x.inf;
         cand2 = b_.inf * x.sup;
         cand3 = b_.sup * x.inf;
         cand4 = b_.sup * x.sup;
         u = MAX(MAX(cand1, cand2), MAX(cand3, cand4));

         SCIPintervalSetRoundingMode(roundmode);
      }
      
      return u;
   }
   
   if( x.sup <= 0.0 )
   { /* change sign of x: enclose a*x^2 + [-bub, -blb]*(-x) for (-x) in [-xub, -xlb] */
      u = x.sup;
      x.sup = -x.inf;
      x.inf = -u;
      b = -b_.inf;
   }
   else
   {
      b = b_.sup;
   }

   if( x.inf >= 0.0 )
   {  /* upper bound for a*x^2 + b*x */
      SCIP_ROUNDMODE roundmode;
      SCIP_Real s,t;
      
      if( b >= infinity )
         return infinity;
    
      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
    
      u = MAX(x.inf * (a*x.inf + b), x.sup * (a*x.sup + b));
      s = b/2;
      t = s/negate(a);
      if( t > x.inf && negate(2*a)*x.sup > b && s*t > u )
         u = s*t;
      
      SCIPintervalSetRoundingMode(roundmode);
      return u;
   }
   else
   {
      SCIP_INTERVAL xlow = x;
      SCIP_Real cand1;
      SCIP_Real cand2;
      assert(x.inf < 0.0 && x.sup > 0);

      xlow.sup = 0;  /* so xlow is lower part of interval */ 
      x.inf = 0;     /* so x    is upper part of interval now */
      cand1 = SCIPintervalQuadUpperBound(infinity, a, b_, xlow);
      cand2 = SCIPintervalQuadUpperBound(infinity, a, b_, x);
      return MAX(cand1, cand2);
   }
}

/** stores range of quadratic term in resultant
 * 
 * given scalar a and intervals b and x, computes interval for \f$ a x^2 + b x \f$ */
void SCIPintervalQuad(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         xrng                /**< range of x */
   )
{
   SCIP_Real tmp;

   if( SCIPintervalIsEmpty(xrng) )
   {
      SCIPintervalSetEmpty(resultant);
      return;
   }
   if( sqrcoeff == 0.0 )
   {
      SCIPintervalMul(infinity, resultant, lincoeff, xrng);
      return;
   }

   resultant->sup =  SCIPintervalQuadUpperBound(infinity,  sqrcoeff, lincoeff, xrng);
  
   tmp = lincoeff.inf;
   lincoeff.inf = -lincoeff.sup;
   lincoeff.sup = -tmp;
   resultant->inf = -SCIPintervalQuadUpperBound(infinity, -sqrcoeff, lincoeff, xrng);
   
   assert(resultant->sup >= resultant->inf);
#if 0
   if( resultant->sup < resultant->inf )
   { /* in case upper bound is below -infinity */
      assert(resultant->inf <= -infinity);
      resultant->sup = resultant->inf;
   }
   else if( resultant->inf > resultant->sup )
   { /* in case lower bound is above  infinity */
      assert(resultant->sup >= infinity);
      resultant->inf = resultant->sup;
   }
#endif
}

/** computes interval with positive solutions of a quadratic equation with interval coefficients
 * 
 * Given intervals a, b, and c, this function computes an interval that contains all positive solutions of \f$ a x^2 + b x \geq c\f$.
 */
void SCIPintervalSolveUnivariateQuadExpressionPositive(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         rhs                 /**< right hand side of equation */
)
{
   assert(resultant != NULL);
  
   /* find x>=0 s.t. a.inf x^2 + b.inf x <= c.sup  -> -a.inf x^2 - b.inf x >= -c.sup */
   if( lincoeff.inf <= -infinity || rhs.sup >= infinity || sqrcoeff.inf <= -infinity )
   {
      resultant->inf = 0.0;
      resultant->sup = infinity;
   }
   else
   {
      SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar(infinity, resultant, -sqrcoeff.inf, -lincoeff.inf, -rhs.sup);
      SCIPdebugMessage("solve %g*x^2 + %g*x >= %g gives [%.20f, %.20f]\n", -sqrcoeff.inf, -lincoeff.inf, -rhs.sup, resultant->inf, resultant->sup);
   }
   
   /* find x>=0 s.t. a.sup x^2 + b.sup x >= c.inf */
   if( lincoeff.sup <  infinity && rhs.inf >  -infinity && sqrcoeff.sup <  infinity )
   {
      SCIP_INTERVAL res2;
      SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar(infinity, &res2, sqrcoeff.sup, lincoeff.sup, rhs.inf);
      SCIPdebugMessage("solve %g*x^2 + %g*x >= %g gives [%.20f, %.20f]\n", sqrcoeff.sup, lincoeff.sup, rhs.inf, res2.inf, res2.sup);
      SCIPdebugMessage("intersect [%.20f, %.20f] and [%.20f, %.20f]", resultant->inf, resultant->sup, res2.inf, res2.sup);
      /* intersect both results */
      SCIPintervalIntersect(resultant, *resultant, res2);
      SCIPdebugPrintf(" gives [%.20f, %.20f]\n", resultant->inf, resultant->sup);
   }
   /* else res2 = [0, infty] */
   
   if( resultant->inf >= infinity || resultant->sup <= -infinity )
   {
      SCIPintervalSetEmpty(resultant);
   }
}

/** computes positive solutions of a quadratic equation with scalar coefficients
 * 
 * Given scalar a, b, and c, this function computes an interval that contains all positive solutions of \f$ a x^2 + b x \geq c\f$.
 * Implements Algorithm 3.2 from Domes and Neumaier: Constraint propagation on quadratic constraints (2008).
 */
void SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             sqrcoeff,           /**< coefficient of x^2 */
   SCIP_Real             lincoeff,           /**< coefficient of x */
   SCIP_Real             rhs                 /**< right hand side of equation */
)
{
   SCIP_ROUNDMODE roundmode;
   SCIP_Real     b;
   SCIP_Real     delta;
   SCIP_Real     z;
   
   assert(resultant != NULL);
   assert(sqrcoeff <  infinity);
   assert(sqrcoeff > -infinity);
  
   resultant->inf = 0.0;
   resultant->sup = infinity;
   
   roundmode = SCIPintervalGetRoundingMode();

   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   b = lincoeff / 2.0;

   if( lincoeff >= 0.0 )
   { /* b >= 0.0 */
      /*SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS); */
      if( rhs > 0.0 )
      { /* b >= 0.0 and c > 0.0 */
         delta = b*b + sqrcoeff*rhs;
         if( delta < 0.0 || (sqrcoeff == 0.0 && lincoeff == 0.0) )
         {
            SCIPintervalSetEmpty(resultant);
         }
         else
         {
            SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
            z = nextafter(sqrt(delta), SCIP_REAL_MAX);
            SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
            z += b;
            resultant->inf = negate(negate(rhs)/z);
            if( sqrcoeff < 0.0 )
               resultant->sup = z / negate(sqrcoeff);
         }
      }
      else
      { /* b >= 0.0 and c <= 0.0 */
         if( sqrcoeff < 0.0 )
         {
            delta = b*b + sqrcoeff*rhs;
            SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
            z = nextafter(sqrt(delta), SCIP_REAL_MAX);
            SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
            z += b;
            resultant->sup = z / negate(sqrcoeff);
         }
      }
   }
   else
   { /* b < 0.0 */
      if( rhs > 0.0 )
      { /* b < 0.0 and c > 0.0 */
         if( sqrcoeff > 0.0 )
         {
            SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
            delta = b*b + sqrcoeff*rhs;
            SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
            z = nextafter(sqrt(delta), SCIP_REAL_MAX);
            SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
            z += negate(b);
            resultant->inf = z / sqrcoeff;
         }
         else
         {
            SCIPintervalSetEmpty(resultant);
         }
      }
      else
      { /* b < 0.0 and c <= 0.0 */
         /* SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS); */
         delta = b*b + sqrcoeff * rhs;
         if( delta >= 0.0 && sqrcoeff <= 0.0 )
         {
            SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
            z = nextafter(sqrt(delta), SCIP_REAL_MAX);
            SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
            z += negate(b);
            resultant->sup = negate(rhs/z);
         }
         /* @TODO actually we could generate a hole here
         if( delta >= 0.0 )
         {
            z = -b + sqrt(delta);
            resultant->sup = -(c/z);
            if( sqrcoeff > 0.0 )
               x2->inf = z/a;
         } */
      }
   }
   
   SCIPintervalSetRoundingMode(roundmode);
}

/** solves a quadratic equation with interval coefficients
 *
 * Given intervals a, b and c, this function computes an interval that contains all solutions of \f$ a x^2 + b x \in c\f$
 */
void SCIPintervalSolveUnivariateQuadExpression(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         rhs                 /**< right hand side of equation */
)
{
   SCIP_Real tmp;
   SCIP_INTERVAL xpos;
   SCIP_INTERVAL xneg;

   assert(resultant != NULL);

   if( sqrcoeff.inf == 0.0 && sqrcoeff.sup == 0.0 )
   { /* relatively easy case: x \in rhs / lincoeff */
      if( lincoeff.inf == 0.0 && lincoeff.sup == 0.0 )
      { /* equation became 0.0 \in rhs */
         if( rhs.inf <= 0.0 && rhs.sup >= 0.0 )
            SCIPintervalSetEntire(infinity, resultant);
         else
            SCIPintervalSetEmpty(resultant);
      }
      else
         SCIPintervalDiv(infinity, resultant, rhs, lincoeff);
      SCIPdebugMessage("  solving [%g,%g]*x in [%g,%g] gives [%g,%g]\n", SCIPintervalGetInf(lincoeff), SCIPintervalGetSup(lincoeff), SCIPintervalGetInf(rhs), SCIPintervalGetSup(rhs), SCIPintervalGetInf(*resultant), SCIPintervalGetSup(*resultant));
      return;
   }

   if( lincoeff.inf == 0.0 && lincoeff.sup == 0.0 )
   { /* easy case: x \in +/- sqrt(rhs/a) */
      SCIPintervalDiv(infinity, resultant, rhs, sqrcoeff);
      SCIPintervalSquareRoot(infinity, resultant, *resultant);
      resultant->inf = -resultant->sup;
      return;
   }

   /* find all x>=0 such that a*x^2+b*x = c */
   SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, &xpos, sqrcoeff, lincoeff, rhs);
   SCIPdebugMessage("  positive solutions of [%g,%g]*x^2 + [%g,%g]*x in [%g,%g] are [%g,%g]\n",
      sqrcoeff.inf, sqrcoeff.sup, lincoeff.inf, lincoeff.sup, rhs.inf, rhs.sup, xpos.inf, xpos.sup);

   tmp = lincoeff.inf;
   lincoeff.inf = -lincoeff.sup;
   lincoeff.sup = -tmp;

   /* find all x>=0 such that a*x^2-b*x = c */
   SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, &xneg, sqrcoeff, lincoeff, rhs);
   tmp = xneg.inf;
   xneg.inf = -xneg.sup;
   xneg.sup = -tmp;
   SCIPdebugMessage("  negative solutions of [%g,%g]*x^2 + [%g,%g]*x in [%g,%g] are [%g,%g]\n",
      sqrcoeff.inf, sqrcoeff.sup, lincoeff.inf, lincoeff.sup, rhs.inf, rhs.sup, xneg.inf, xneg.sup);

   SCIPintervalUnify(resultant, xpos, xneg);
   SCIPdebugMessage("  unify gives [%g,%g]\n", SCIPintervalGetInf(*resultant), SCIPintervalGetSup(*resultant));
}
