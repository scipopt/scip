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

/**@file   bandit_epsgreedy.c
 * @brief  implementation of epsilon greedy bandit algorithm
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "blockmemshell/memory.h"
#include "scip/bandit_epsgreedy.h"
#include "scip/scip.h"

#define BANDIT_NAME           "eps-greedy"
#define DEFAULT_WEIGHT 0.2
#define DEFAULT_INITSEED 999
/*
 * Data structures
 */

/** private data structure of epsilon greedy bandit algorithm */
struct SCIP_BanditData
{
   SCIP_Real             eps;                /**< epsilon parameter (between 0 and 1) to control epsilon greedy */
   SCIP_Real*            weights;            /**< weights for every action */
   int                   nselections;        /**< counter for the number of selection calls */
};

/*
 * Callback methods of bandit algorithm virtual function table
 */

/** callback to free bandit specific data structures */
SCIP_DECL_BANDITFREE(SCIPbanditFreeEpsgreedy)
{  /*lint --e{715}*/

   SCIP_BANDITDATA* banditdata;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   assert(banditdata->weights != NULL);

   BMSfreeBlockMemoryArray(blkmem, &banditdata->weights, SCIPbanditGetNActions(bandit));
   BMSfreeBlockMemory(blkmem, &banditdata);

   SCIPbanditSetData(bandit, NULL);

   return SCIP_OKAY;
}

/** selection callback for bandit algorithm */
SCIP_DECL_BANDITSELECT(SCIPbanditSelectEpsgreedy)
{  /*lint --e{715}*/

   SCIP_BANDITDATA* banditdata;
   SCIP_Real randnr;
   SCIP_Real curreps;
   SCIP_RANDNUMGEN* rng;
   int nactions;
   assert(bandit != NULL);
   assert(selection != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);
   rng = SCIPbanditGetRandnumgen(bandit);
   assert(rng != NULL);

   nactions = SCIPbanditGetNActions(bandit);

   /** roll the dice to check if the best element should be picked, or an element at random */
   randnr = SCIPrandomGetReal(rng, 0.0, 1.0);

   /* make epsilon decrease with an increasing number of selections */
   banditdata->nselections++;
   curreps = banditdata->eps * sqrt((SCIP_Real)nactions/(SCIP_Real)banditdata->nselections);

   /* select the best action seen so far */
   if( randnr >= curreps )
   {
      SCIP_Real* weights = banditdata->weights;
      int j;
      SCIP_Real maxreward;

      assert(weights != NULL);

      /** pick the element with the largest reward */
      maxreward = weights[0];
      *selection = 0;

      /* determine reward for every element */
      for( j = 1; j < nactions; ++j )
      {
         SCIP_Real reward = weights[j];

         if( maxreward < reward )
         {
            *selection = j;
            maxreward = reward;
         }
      }
   }
   else
   {
      /* play one of the actions at random */
      *selection = SCIPrandomGetInt(rng, 0, nactions - 1);
   }

   return SCIP_OKAY;
}

/** update callback for bandit algorithm */
SCIP_DECL_BANDITUPDATE(SCIPbanditUpdateEpsgreedy)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);

   /* use geometrically decreasing weights for older observations */
   banditdata->weights[selection] *= 0.9;
   banditdata->weights[selection] += 0.1 * score;

   return SCIP_OKAY;
}

/** reset callback for bandit algorithm */
SCIP_DECL_BANDITRESET(SCIPbanditResetEpsgreedy)
{  /*lint --e{715}*/
   SCIP_BANDITDATA* banditdata;
   SCIP_Real* weights;
   int w;
   int nactions;

   assert(bandit != NULL);

   banditdata = SCIPbanditGetData(bandit);
   assert(banditdata != NULL);

   weights = banditdata->weights;
   assert(weights != NULL);
   nactions = SCIPbanditGetNActions(bandit);
   assert(nactions > 0);

   /* reset weights to reflect priorities, if priorities are given */
   if( priorities != NULL )
   {
      SCIP_Real priosum;
      SCIP_Real normalization;

      /* compute priority sum for normalization */
      priosum = priorities[0];
      for( w = 1; w < nactions; ++w )
         priosum += priorities[w];

      /* use DEFAULT_WEIGHT as the standard epsilon greedy weight for uninitialized neighborhoods */
      normalization = DEFAULT_WEIGHT * nactions / priosum;

      /* reset weights */
      for( w = 0; w < nactions; ++w )
         weights[w] = priorities[w] * normalization;
   }
   else
   {
      /* reset weights */
      for( w = 0; w < nactions; ++w )
         weights[w] = DEFAULT_WEIGHT;
   }

   banditdata->nselections = 0;

   return SCIP_OKAY;
}

/*
 * interface methods of the Epsilon Greedy bandit algorithm
 */

/** internal method to create an epsilon greedy bandit algorithm */
SCIP_RETCODE SCIPbanditCreateEpsgreedy(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_BANDITVTABLE*    vtable,             /**< virtual function table with epsilon greedy callbacks */
   SCIP_BANDIT**         epsgreedy,          /**< pointer to store the epsilon greedy bandit algorithm */
   SCIP_Real             eps,                /**< probability for exploration between all actions */
   int                   nactions,           /**< the number of possible actions */
   unsigned int          initseed            /**< initial random seed */
   )
{
   SCIP_BANDITDATA* banditdata;

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &banditdata) );
   assert(banditdata != NULL);

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &banditdata->weights, nactions) );
   banditdata->eps = eps;
   banditdata->nselections = 0;

   SCIP_CALL( SCIPbanditCreate(epsgreedy, vtable, blkmem, nactions, initseed, banditdata) );

   return SCIP_OKAY;
}

/** create an epsilon greedy bandit algorithm */
SCIP_RETCODE SCIPcreateBanditEpsgreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         epsgreedy,          /**< pointer to store the epsilon greedy bandit algorithm */
   SCIP_Real             eps,                /**< probability for exploration between all actions */
   int                   nactions            /**< the number of possible actions */
   )
{
   SCIP_BANDITVTABLE* vtable;
   assert(scip != NULL);
   assert(epsgreedy != NULL);
   assert(eps <= 1.0);
   assert(eps >= 0.0);
   assert(nactions > 0);

   vtable = SCIPfindBanditvtable(scip, BANDIT_NAME);
   if( vtable == NULL )
   {
      SCIPerrorMessage("Could not find virtual function table for %s bandit algorithm\n", BANDIT_NAME);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbanditCreateEpsgreedy(SCIPblkmem(scip), vtable, epsgreedy, eps, nactions, SCIPinitializeRandomSeed(scip, DEFAULT_INITSEED)) );

   return SCIP_OKAY;
}

/** get weights array of epsilon greedy bandit algorithm */
SCIP_Real* SCIPgetWeightsEpsgreedy(
   SCIP_BANDIT*          epsgreedy           /**< epsilon greedy bandit algorithm */
   )
{
   SCIP_BANDITDATA* banditdata;
   assert(epsgreedy != NULL);
   banditdata = SCIPbanditGetData(epsgreedy);
   assert(banditdata != NULL);

   return banditdata->weights;
}

/** set epsilon parameter of epsilon greedy bandit algorithm */
void SCIPsetEpsilonEpsgreedy(
   SCIP_BANDIT*          epsgreedy,          /**< epsilon greedy bandit algorithm */
   SCIP_Real             eps                 /**< epsilon parameter (increase for more exploration) */
   )
{
   SCIP_BANDITDATA* banditdata;
   assert(epsgreedy != NULL);

   banditdata = SCIPbanditGetData(epsgreedy);

   banditdata->eps = eps;
}


/** creates the epsilon greedy bandit algorithm includes it in SCIP */
SCIP_RETCODE SCIPincludeBanditvtableEpsgreedy(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BANDITVTABLE* banditvtable;

   SCIP_CALL( SCIPincludeBanditvtable(scip, &banditvtable, BANDIT_NAME,
         SCIPbanditFreeEpsgreedy, SCIPbanditSelectEpsgreedy, SCIPbanditUpdateEpsgreedy, SCIPbanditResetEpsgreedy) );

   return SCIP_OKAY;
}
