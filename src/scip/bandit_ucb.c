/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bandit_ucb.c
 * @brief  methods for UCB bandit selection
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/bandit_ucb.h"

#define BANDIT_NAME "ucb"

/*
 * Data structures
 */

/** implementation specific data of UCB bandit algorithm */
struct SCIP_BanditData
{
   int                   nselections;        /**< counter for the number of selections */
   int*                  counter;            /**< array of counters how often every action has been chosen */
   int*                  startperm;          /**< indices for starting permutation */
   SCIP_Real*            meanscores;         /**< array of average scores for the actions */
   SCIP_Real             alpha;              /**< parameter to increase confidence width */
};


/*
 * Local methods
 */


/*
 * Callback methods of bandit algorithm
 */

/** callback to free bandit specific data structures */
static
SCIP_DECL_BANDITFREE(banditFreeUcb)
{  /*lint --e{715}*/

   SCIP_BANDITDATA* banditdata;
   int nactions;
   assert(scip != NULL);
   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   SCIPfreeBlockMemoryArray(scip, &banditdata->counter, nactions);
   SCIPfreeBlockMemoryArray(scip, &banditdata->startperm, nactions);
   SCIPfreeBlockMemoryArray(scip, &banditdata->meanscores, nactions);
   SCIPfreeBlockMemory(scip, &banditdata);

   SCIPbanditSetData(bandit, NULL);

   return SCIP_OKAY;
}

/** selection callback for bandit selector */
static
SCIP_DECL_BANDITSELECT(banditSelectUcb)
{  /*lint --e{715}*/

   SCIP_BANDITDATA* banditdata;
   int nactions;
   int* counter;
   SCIP_Real* meanscores;

   assert(scip != NULL);
   assert(bandit != NULL);
   assert(selection != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   counter = banditdata->counter;
   /* select the next uninitialized action from the start permutation */
   if( banditdata->nselections < nactions )
   {
      *selection = banditdata->startperm[banditdata->nselections];
      assert(counter[*selection] == 0);
   }
   else
   {
      /* select the action with the highest upper confidence bound */
      SCIP_Real widthfactor;
      SCIP_Real maxucb;
      int i;

      /* compute the confidence width factor that is common for all actions */
      widthfactor = banditdata->alpha * log(1 + banditdata->nselections);
      widthfactor = sqrt(widthfactor);
      maxucb = -1.0;

      /* loop over the actions and determine the maximum upper confidence bound.
       * The upper confidence bound of an action is the sum of its mean score
       * plus a confidence term that decreases with increasing number of observations of
       * this action.
       */
      for( i = 0; i < nactions; ++i )
      {
         SCIP_Real uppercb;
         SCIP_Real rootcount;
         assert(counter[i] > 0);

         /* compute the upper confidence bound for action i */
         uppercb = meanscores[i];
         rootcount = sqrt((SCIP_Real)counter[i]);
         uppercb += widthfactor / rootcount;
         assert(uppercb > 0);

         /* update maximum */
         if( uppercb > maxucb )
         {
            maxucb = uppercb;
            *selection = i;
         }
      }
   }

   assert(*selection > 0);
   assert(*selection <= nactions);
   banditdata->nselections++;

   return SCIP_OKAY;
}

/** update callback for bandit algorithm */
static
SCIP_DECL_BANDITUPDATE(banditUpdateUcb)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   int nactions;
   SCIP_Real delta;


   assert(scip != NULL);
   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   assert(selection >= 0);
   assert(selection < nactions);

   /* increase the mean by the incremental formula: A_n = A_n-1 + 1/n (a_n - A_n-1) */
   delta = score - banditdata->meanscores[selection];
   ++banditdata->counter[selection];
   banditdata->meanscores[selection] += delta / (SCIP_Real)banditdata->counter[selection];

   return SCIP_OKAY;
}

/** reset callback for bandit algorithm */
static
SCIP_DECL_BANDITRESET(banditResetUcb)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   int nactions;
   int i;

   assert(scip != NULL);
   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   nactions = SCIPbanditGetNActions(bandit);

   /* clear counters and scores */
   BMSclearMemoryArray(banditdata->counter, nactions);
   BMSclearMemoryArray(banditdata->meanscores, nactions);
   banditdata->nselections = 0;

   /* initialize start permutation as identity */
   for( i = 0; i < nactions; ++i )
      banditdata->startperm[i] = i;

   /* prepare the start permutation in decreasing order of priority */
   if( priorities != NULL )
   {
      SCIP_Real* prioritycopy;

      SCIP_CALL( SCIPduplicateBufferArray(scip, &prioritycopy, priorities, nactions) );

      SCIPsortDownRealInt(prioritycopy, banditdata->startperm, nactions);

      SCIPfreeBufferArray(scip, &prioritycopy);
   }
   else
   {
      /* use a random start permutation */
      SCIP_RANDNUMGEN* rng = SCIPbanditGetRandnumgen(bandit);
      assert(rng != NULL);

      SCIPrandomPermuteIntArray(rng, banditdata->startperm, 0, nactions);
   }

   return SCIP_OKAY;
}

/*
 * bandit algorithm specific interface methods
 */

/** returns the upper confidence bound of a selected action */
SCIP_Real SCIPgetConfidenceBoundUcb(
   SCIP_BANDIT*          ucb,                /**< UCB bandit algorithm */
   int                   action              /**< index of the queried action */
   )
{
   SCIP_Real uppercb;
   SCIP_BANDITDATA* banditdata;
   int nactions;

   assert(ucb != NULL);
   banditdata = SCIPbanditGetData(ucb);
   nactions = SCIPbanditGetNActions(ucb);
   assert(action < nactions);

   /* since only scores between 0 and 1 are allowed, 1.0 is a sure upper confidence bound */
   if( banditdata->nselections < nactions )
      return 1.0;

   /* the bandit algorithm must have picked every action once */
   assert(banditdata->counter[action] > 0);
   uppercb = banditdata->meanscores[action];
   uppercb += sqrt(banditdata->alpha * log(1.0 + banditdata->nselections) / (SCIP_Real)banditdata->counter[action]);

   return uppercb;
}

/** return start permutation of the UCB bandit algorithm */
int* SCIPgetStartPermutationUcb(
   SCIP_BANDIT*          ucb                 /**< UCB bandit algorithm */
   )
{
   SCIP_BANDITDATA* banditdata = SCIPbanditGetData(ucb);

   assert(banditdata != NULL);

   return banditdata->startperm;
}

/** create UCB bandit algorithm */
SCIP_RETCODE SCIPcreateBanditUcb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         ucb,                /**< pointer to store bandit algorithm */
   int                   nactions,           /**< the number of actions for this bandit algorithm */
   SCIP_Real             alpha               /**< parameter to increase confidence width */
   )
{
   SCIP_BANDITDATA* banditdata;
   SCIP_BANDITVTABLE* vtable;

   vtable = SCIPfindBanditvtable(scip, BANDIT_NAME);
   if( vtable == NULL )
   {
      SCIPerrorMessage("Could not find virtual function table for %s bandit algorithm\n", BANDIT_NAME);
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, &banditdata) );
   assert(banditdata != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &banditdata->counter, nactions) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &banditdata->startperm, nactions) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &banditdata->meanscores, nactions) );

   banditdata->alpha = alpha;

   SCIP_CALL( SCIPcreateBandit(scip, ucb, vtable, nactions, banditdata) );

   return SCIP_OKAY;
}

/** include virtual function table for UCB bandit algorithms */
SCIP_RETCODE SCIPincludeBanditvtableUcb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BANDITVTABLE* vtable;

   SCIP_CALL( SCIPincludeBanditvtable(scip, &vtable, BANDIT_NAME,
         banditFreeUcb, banditSelectUcb, banditUpdateUcb, banditResetUcb) );
   assert(vtable != NULL);

   return SCIP_OKAY;
}
