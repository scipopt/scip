/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   random.c
 * @brief  random number generator
 * @author Jakob Witzig
 */

#include "scip/random.h"


#define DEFAULT_XOR  362436000;
#define DEFAULT_MWC  521288629;
#define DEFAULT_CST  7654321;


SCIP_RETCODE SCIPrandomInit(
   SCIP_RANDGEN*         randgen,
   int                   initseed
   )
{
   assert(randgen != NULL);

   randgen->seed = (uint32_t)initseed;
   randgen->xor = DEFAULT_XOR;
   randgen->mwc = DEFAULT_MWC;
   randgen->cst = DEFAULT_CST;

   return SCIP_OKAY;
}

/** returns a random number between 0 and INT_MAX
 *
 *  implementation of KISS random number generator developed by George Marsaglia.
 *  KISS is combination of three different random number generators:
 *   - Lbinear congruential generator
 *   - Xorshift
 *   - Lag-1 Multiply-with-carry
 *
 *  KISS has a period of 2^123 and passes all statistical test part of BigCrush-Test of TestU01 [1].
 *
 *  [1] http://dl.acm.org/citation.cfm?doid=1268776.1268777
 */
static
int getRand(
   SCIP_RANDGEN*         randgen             /**< random number generator */
   )
{
   uint64_t t;

   /* linear congruential */
   randgen->seed = randgen->seed * (SCIP_Longint)1103515245 + 12345;

   /* Xorshift */
   randgen->xor ^= randgen->xor << 13;
   randgen->xor ^= randgen->xor >> 17;
   randgen->xor ^= randgen->xor << 5;

   /* Multiple-with-carry */
   t = 698769069ULL * randgen->mwc + randgen->cst;
   randgen->cst = t >> 32;
   randgen->mwc = (uint32_t) t;

   return (int)((randgen->seed + randgen->xor + randgen->mwc) % INT_MAX);
}

/** returns a random integer between minrandval and maxrandval */
int SCIPrandomGetInt(
   SCIP_RANDGEN*         randgen,            /**< random number generator */
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval          /**< maximal value to return */
   )
{
   SCIP_Real randnumber;

   randnumber = (SCIP_Real)getRand(randgen)/(INT_MAX+1.0);
   assert(randnumber >= 0.0);
   assert(randnumber < 1.0);

   /* we multiply minrandval and maxrandval separately by randnumber in order to avoid overflow if they are more than INT_MAX
    * apart
    */
   return (int) (minrandval*(1.0 - randnumber) + maxrandval*randnumber + randnumber);
}

/** returns a random real between minrandval and maxrandval */
SCIP_Real SCIPrandomGetReal(
   SCIP_RANDGEN*         randgen,            /**< random number generator */
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval          /**< maximal value to return */
   )
{
   SCIP_Real randnumber;

   randnumber = (SCIP_Real)getRand(randgen)/(SCIP_Real)INT_MAX;
   assert(randnumber >= 0.0);
   assert(randnumber <= 1.0);

   /* we multiply minrandval and maxrandval separately by randnumber in order to avoid overflow if they are more than
    * SCIP_REAL_MAX apart
    */
   return minrandval*(1.0 - randnumber) + maxrandval*randnumber;
}
