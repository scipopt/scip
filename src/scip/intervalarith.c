/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: intervalarith.c,v 1.4 2004/07/13 15:03:50 bzfpfend Exp $"

/**@file   intervalarith.c
 * @brief  interval arithmetics for provable bounds
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>

#include "message.h"
#include "intervalarith.h"



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
      errorMessage("error setting rounding mode to %d\n", roundmode);
      abort();
   }
}

/** gets current rounding mode of floating point operations */
static
ROUNDMODE getRoundingMode(
   void
   )
{
   return fegetround();
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
      errorMessage("error setting rounding mode to %d\n", roundmode);
      abort();
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
{
   warningMessage("setting rounding mode not available - interval arithmetic is invalid!\n");
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
   INTERVAL*        resultant,          /**< interval to store value into */
   Real             value               /**< value to store */
   )
{
   assert(resultant != NULL);

   resultant->inf = value;
   resultant->sup = value;
}

/** stores given infimum and supremum as interval */
void SCIPintervalSetBounds(
   INTERVAL*        resultant,          /**< interval to store value into */
   Real             inf,                /**< value to store as infimum */
   Real             sup                 /**< value to store as supremum */
   )
{
   assert(resultant != NULL);
   assert(inf <= sup);

   resultant->inf = inf;
   resultant->sup = sup;
}

/** adds operand1 and operand2 and stores result in resultant */
void SCIPintervalAdd(
   INTERVAL*        resultant,          /**< resultant interval of operation */
   INTERVAL         operand1,           /**< first operand of operation */
   INTERVAL         operand2            /**< second operand of operation */
   )
{
   ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   roundmode = getRoundingMode();

   setRoundingMode(SCIP_ROUND_DOWNWARDS);
   resultant->inf = operand1.inf + operand2.inf;
   setRoundingMode(SCIP_ROUND_UPWARDS);
   resultant->sup = operand1.sup + operand2.sup;

   setRoundingMode(roundmode);
}

/** substracts operand2 from operand1 and stores result in resultant */
void SCIPintervalSub(
   INTERVAL*        resultant,          /**< resultant interval of operation */
   INTERVAL         operand1,           /**< first operand of operation */
   INTERVAL         operand2            /**< second operand of operation */
   )
{
   ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   roundmode = getRoundingMode();

   setRoundingMode(SCIP_ROUND_DOWNWARDS);
   resultant->inf = operand1.inf - operand2.sup;
   setRoundingMode(SCIP_ROUND_UPWARDS);
   resultant->sup = operand1.sup - operand2.inf;

   setRoundingMode(roundmode);
}

/** multiplies operand1 with operand2 and stores result in resultant */
void SCIPintervalMul(
   INTERVAL*        resultant,          /**< resultant interval of operation */
   INTERVAL         operand1,           /**< first operand of operation */
   INTERVAL         operand2            /**< second operand of operation */
   )
{
   ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(operand1.inf <= operand1.sup);
   assert(operand2.inf <= operand2.sup);

   roundmode = getRoundingMode();

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
         Real x;
         Real y;

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

   setRoundingMode(roundmode);
}

/** returns infimum of interval */
Real SCIPintervalGetInf(
   INTERVAL         interval            /**< interval */
   )
{
   return interval.inf;
}

/** returns supremum of interval */
Real SCIPintervalGetSup(
   INTERVAL         interval            /**< interval */
   )
{
   return interval.sup;
}
