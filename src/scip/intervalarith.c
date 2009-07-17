/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: intervalarith.c,v 1.19 2009/07/17 15:31:23 bzfviger Exp $"

/**@file   intervalarith.c
 * @brief  interval arithmetics for provable bounds
 * @author Tobias Achterberg
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
enum RoundMode
{
   SCIP_ROUND_DOWNWARDS = FE_DOWNWARD,  /**< round always down */
   SCIP_ROUND_UPWARDS   = FE_UPWARD     /**< round always up */
};
typedef enum RoundMode ROUNDMODE;

/** sets rounding mode of floating point operations */
static
void setRoundingMode(
   ROUNDMODE        roundmode           /**< rounding mode to activate */
   )
{
   if( fesetround(roundmode) != 0 )
   {
      SCIPerrorMessage("error setting rounding mode to %d\n", roundmode);
      SCIPABORT();
   }
}

/** gets current rounding mode of floating point operations */
static
ROUNDMODE getRoundingMode(
   void
   )
{
   return (ROUNDMODE)fegetround();
}
#endif



#ifdef ROUNDING_FP
#define ROUNDING
/*
 * OSF rounding operations
 */

#include <float.h>

/** OSF rounding mode settings */
enum RoundMode
{
   SCIP_ROUND_DOWNWARDS = FP_RND_RM,    /**< round always down */
   SCIP_ROUND_UPWARDS   = FP_RND_RP     /**< round always up */
};
typedef enum RoundMode ROUNDMODE;

/** sets rounding mode of floating point operations */
static
void setRoundingMode(
   ROUNDMODE        roundmode           /**< rounding mode to activate */
   )
{
   if( write_rnd(roundmode) != 0 )
   {
      SCIPerrorMessage("error setting rounding mode to %d\n", roundmode);
      SCIPABORT();
   }
}

/** gets current rounding mode of floating point operations */
static
ROUNDMODE getRoundingMode(
   void
   )
{
   return read_rnd();
}
#endif



#ifndef ROUNDING
/*
 * rouding operations not available
 */
enum RoundMode
{
   SCIP_ROUND_DOWNWARDS = 0,            /**< round always down */
   SCIP_ROUND_UPWARDS   = 1             /**< round always up */
};
typedef enum RoundMode ROUNDMODE;

/** sets rounding mode of floating point operations */
static
void setRoundingMode(
   ROUNDMODE        roundmode           /**< rounding mode to activate */
   )
{  /*lint --e{715}*/
   SCIPwarningMessage("setting rounding mode not available - interval arithmetic is invalid!\n");
}

/** gets current rounding mode of floating point operations */
static
ROUNDMODE getRoundingMode(
   void
   )
{
   return SCIP_ROUND_DOWNWARDS;
}
#else
#undef ROUNDING
#endif




/*
 * Interval arithmetic operations
 */

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

/** adds operand1 and operand2 and stores result in resultant */
void SCIPintervalAdd(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   roundmode = getRoundingMode();

   if( operand1.inf <= -infinity || operand2.inf <= -infinity )
   {
      resultant->inf = -infinity;
   }
   else
   {
      setRoundingMode(SCIP_ROUND_DOWNWARDS);
      resultant->inf = operand1.inf + operand2.inf;
   }
   
   if( operand1.sup >=  infinity || operand2.sup >=  infinity )
   {
      resultant->sup =  infinity;
   }
   else
   {
      setRoundingMode(SCIP_ROUND_UPWARDS);
      resultant->sup = operand1.sup + operand2.sup;
   }
   
   setRoundingMode(roundmode);
}

/** substracts operand2 from operand1 and stores result in resultant */
void SCIPintervalSub(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   roundmode = getRoundingMode();

   if( operand1.inf <= -infinity || operand2.sup >=  infinity )
   {
      resultant->inf = -infinity;
   }
   else
   {
      setRoundingMode(SCIP_ROUND_DOWNWARDS);
      resultant->inf = operand1.inf - operand2.sup;
   }
   
   if( operand1.sup >=  infinity || operand2.inf <= -infinity )
   {
      resultant->sup =  infinity;
   }
   else
   {
      setRoundingMode(SCIP_ROUND_UPWARDS);
      resultant->sup = operand1.sup - operand2.inf;
   }

   setRoundingMode(roundmode);
}

/** multiplies operand1 with operand2 and stores result in resultant */
void SCIPintervalMul(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   ROUNDMODE roundmode;
   SCIP_Real cand1;
   SCIP_Real cand2;
   SCIP_Real cand3;
   SCIP_Real cand4;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   roundmode = getRoundingMode();

   if( (operand1.inf <= -infinity && operand2.sup > 0         ) ||
       (operand1.sup >   0        && operand2.inf <= -infinity) ||
       (operand1.inf <   0        && operand2.sup >=  infinity) ||
       (operand1.sup >=  infinity && operand2.inf <   0       ) )
   {
      resultant->inf = -infinity;
   }
   else
   {
      setRoundingMode(SCIP_ROUND_DOWNWARDS);
      cand1 = operand1.inf * operand2.inf;
      cand2 = operand1.inf * operand2.sup;
      cand3 = operand1.sup * operand2.inf;
      cand4 = operand1.sup * operand2.sup;
      resultant->inf = MIN(MIN(cand1, cand2), MIN(cand3, cand4));
   }
   
   if( (operand1.inf <= -infinity && operand2.inf <   0       ) ||
       (operand1.inf <   0        && operand2.inf <= -infinity) ||
       (operand1.sup >   0        && operand2.sup >=  infinity) ||
       (operand1.sup >=  infinity && operand2.sup >   0       ) )
   {
      resultant->sup =  infinity;
   }
   else
   {
      setRoundingMode(SCIP_ROUND_UPWARDS);
      cand1 = operand1.inf * operand2.inf;
      cand2 = operand1.inf * operand2.sup;
      cand3 = operand1.sup * operand2.inf;
      cand4 = operand1.sup * operand2.sup;
      resultant->sup = MAX(MAX(cand1, cand2), MAX(cand3, cand4));
   }
   
#if 0   
   if( operand1.inf >= 0.0 )
   {
      if( operand2.inf >= 0.0 )
      {
         /* [+,+] * [+,+] */
         setRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.inf * operand2.inf;
         setRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.sup * operand2.sup;
      }
      else if( operand2.sup >= 0.0 )
      {
         /* [+,+] * [-,+] */
         setRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.sup * operand2.inf;
         setRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.sup * operand2.sup;
      }
      else
      {
         /* [+,+] * [-,-] */
         setRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.sup * operand2.inf;
         setRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.inf * operand2.sup;
      }
   }
   else if( operand1.sup >= 0.0 )
   {
      if( operand2.inf >= 0.0 )
      {
         /* [-,+] * [+,+] */
         setRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.inf * operand2.sup;
         setRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.sup * operand2.sup;
      }
      else if( operand2.sup >= 0.0 )
      {
         SCIP_Real x;
         SCIP_Real y;

         /* [-,+] * [-,+] */
         setRoundingMode(SCIP_ROUND_DOWNWARDS);
         x = operand1.inf * operand2.sup;
         y = operand1.sup * operand2.inf;
         resultant->inf = MIN(x, y);
         setRoundingMode(SCIP_ROUND_UPWARDS);
         x = operand1.inf * operand2.inf;
         y = operand1.sup * operand2.sup;
         resultant->sup = MAX(x, y);
      }
      else
      {
         /* [-,+] * [-,-] */
         setRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.sup * operand2.inf;
         setRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.inf * operand2.inf;
      }
   }
   else
   {
      if( operand2.inf >= 0.0 )
      {
         /* [-,-] * [+,+] */
         setRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.inf * operand2.sup;
         setRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.sup * operand2.inf;
      }
      else if( operand2.sup >= 0.0 )
      {
         /* [-,-] * [-,+] */
         setRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.inf * operand2.sup;
         setRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.inf * operand2.inf;
      }
      else
      {
         /* [-,-] * [-,-] */
         setRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.sup * operand2.sup;
         setRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.inf * operand2.inf;
      }
   }
#endif
   
   setRoundingMode(roundmode);
}

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
