/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bandit_exp3ix.c
 * @ingroup OTHER_CFILES
 * @brief  methods for Exp.3-IX bandit selection
 * @author Antonia Chmiela
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/bandit.h"
#include "scip/bandit_exp3ix.h"
#include "scip/pub_bandit.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/scip_bandit.h"
#include "scip/scip_mem.h"
#include "scip/scip_randnumgen.h"

#define BANDIT_NAME "exp3ix"

/*
 * Data structures
 */

/** implementation specific data of Exp.3 bandit algorithm */
struct SCIP_BanditData
{
   SCIP_Real*            weights;            /**< exponential weight for each arm */
   SCIP_Real             weightsum;          /**< the sum of all weights */
   int                   iter;               /**< current iteration counter to compute parameters gamma_t and eta_t */
};

/*
 * Local methods
 */

/*
 * Callback methods of bandit algorithm
 */

/** callback to free bandit specific data structures */
SCIP_DECL_BANDITFREE(SCIPbanditFreeExp3IX)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   int nactions;
   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   BMSfreeBlockMemoryArray(blkmem, &banditdata->weights, nactions);

   BMSfreeBlockMemory(blkmem, &banditdata);

   SCIPbanditSetData(bandit, NULL);

   return SCIP_OKAY;
}

/** selection callback for bandit selector */
SCIP_DECL_BANDITSELECT(SCIPbanditSelectExp3IX)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   SCIP_RANDNUMGEN* rng;
   SCIP_Real* weights;
   SCIP_Real weightsum;
   int i;
   int nactions;
   SCIP_Real psum;
   SCIP_Real randnr;

   assert(bandit != NULL);
   assert(selection != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   rng = SCIPbanditGetRandnumgen(bandit);
   assert(rng != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   /* initialize some local variables to speed up probability computations */
   weightsum = banditdata->weightsum;
   weights = banditdata->weights;

   /* draw a random number between 0 and 1 */
   randnr = SCIPrandomGetReal(rng, 0.0, 1.0);

   /* loop over probability distribution until rand is reached
    * the loop terminates without looking at the last action,
    * which is then selected automatically if the target probability
    * is not reached earlier
    */
   psum = 0.0;
   for( i = 0; i < nactions - 1; ++i )
   {
      SCIP_Real prob;

      /* compute the probability for arm i */
      prob = weights[i] / weightsum;
      psum += prob;

      /* break and select element if target probability is reached */
      if( randnr <= psum )
         break;
   }

   /* select element i, which is the last action in case that the break statement hasn't been reached */
   *selection = i;

   return SCIP_OKAY;
}

/** compute gamma_t */
static
SCIP_Real SCIPcomputeGamma(
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   int                   t                   /**< current iteration */
   )
{
   return sqrt(log((SCIP_Real)nactions) / (4.0 * (SCIP_Real)t * (SCIP_Real)nactions));
}

/** update callback for bandit algorithm */
SCIP_DECL_BANDITUPDATE(SCIPbanditUpdateExp3IX)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   SCIP_Real etaparam;
   SCIP_Real lossestim;
   SCIP_Real prob;
   SCIP_Real weightsum;
   SCIP_Real newweightsum;
   SCIP_Real* weights;
   SCIP_Real gammaparam;
   int nactions;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   assert(selection >= 0);
   assert(selection < nactions);

   weights = banditdata->weights;
   weightsum = banditdata->weightsum;
   newweightsum = weightsum;
   gammaparam = SCIPcomputeGamma(nactions, banditdata->iter);
   etaparam = 2.0 * gammaparam;

   /* probability of selection */
   prob = weights[selection] / weightsum;

   /* estimated loss */
   lossestim = (1.0 - score) / (prob + gammaparam);
   assert(lossestim >= 0);

   /* update the observation for the current arm */
   newweightsum -= weights[selection];
   weights[selection] *= exp(-etaparam * lossestim);
   newweightsum += weights[selection];

   banditdata->weightsum = newweightsum;

   /* increase iteration counter */
   banditdata->iter += 1;

   return SCIP_OKAY;
}

/** reset callback for bandit algorithm */
SCIP_DECL_BANDITRESET(SCIPbanditResetExp3IX)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   SCIP_Real* weights;
   int nactions;
   int i;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);
   weights = banditdata->weights;

   assert(nactions > 0);

   /* initialize all weights with 1.0 */
   for( i = 0; i < nactions; ++i )
      weights[i] = 1.0;

   banditdata->weightsum = (SCIP_Real)nactions;

   /* set iteration counter to 1 */
   banditdata->iter = 1;

   return SCIP_OKAY;
}


/*
 * bandit algorithm specific interface methods
 */

/** direct bandit creation method for the core where no SCIP pointer is available */
SCIP_RETCODE SCIPbanditCreateExp3IX(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_BANDITVTABLE*    vtable,             /**< virtual function table for callback functions of Exp.3-IX */
   SCIP_BANDIT**         exp3ix,             /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial random seed */
   )
{
   SCIP_BANDITDATA* banditdata;

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &banditdata) );
   assert(banditdata != NULL);

   banditdata->iter = 1;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &banditdata->weights, nactions) );

   SCIP_CALL( SCIPbanditCreate(exp3ix, vtable, blkmem, bufmem, priorities, nactions, initseed, banditdata) );

   return SCIP_OKAY;
}

/** creates and resets an Exp.3-IX bandit algorithm using \p scip pointer */
SCIP_RETCODE SCIPcreateBanditExp3IX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         exp3ix,             /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial seed for random number generation */
   )
{
   SCIP_BANDITVTABLE* vtable;

   vtable = SCIPfindBanditvtable(scip, BANDIT_NAME);
   if( vtable == NULL )
   {
      SCIPerrorMessage("Could not find virtual function table for %s bandit algorithm\n", BANDIT_NAME);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbanditCreateExp3IX(SCIPblkmem(scip), SCIPbuffer(scip), vtable, exp3ix,
         priorities, nactions, SCIPinitializeRandomSeed(scip, initseed)) );

   return SCIP_OKAY;
}

/** returns probability to play an action */
SCIP_Real SCIPgetProbabilityExp3IX(
   SCIP_BANDIT*          exp3ix,             /**< bandit algorithm */
   int                   action              /**< index of the requested action */
   )
{
   SCIP_BANDITDATA* banditdata = SCIPbanditGetData(exp3ix);

   assert(banditdata->weightsum > 0.0);
   assert(SCIPbanditGetNActions(exp3ix) > 0);

   return banditdata->weights[action] / banditdata->weightsum;
}

/** include virtual function table for Exp.3-IX bandit algorithms */
SCIP_RETCODE SCIPincludeBanditvtableExp3IX(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BANDITVTABLE* vtable;

   SCIP_CALL( SCIPincludeBanditvtable(scip, &vtable, BANDIT_NAME,
         SCIPbanditFreeExp3IX, SCIPbanditSelectExp3IX, SCIPbanditUpdateExp3IX, SCIPbanditResetExp3IX) );
   assert(vtable != NULL);

   return SCIP_OKAY;
}
