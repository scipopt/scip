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

/**@file   iis_greedy.c
 * @brief  greedy deletion and addition filter heuristic to compute (I)ISs
 * @author Marc Pfetsch
 * @author Mark Turner
 *
 */

#include <assert.h>
#include <string.h>

#include "scip/scip_iis.h"
#include <scip/iis_greedy.h>

#define IIS_NAME                 "greedy"
#define IIS_DESC                 "greedy deletion or addition constraint deletion"
#define IIS_PRIORITY              8000
#define RANDSEED                  0x5EED

#define DEFAULT_BATCHSIZE        1
#define DEFAULT_MAXNNODESPERITER INT_MAX
#define DEFAULT_TIMELIMPERITER   1e+20
#define DEFAULT_ADDITIVE         TRUE
#define DEFAULT_CONSERVATIVE     TRUE
#define DEFAULT_DYNAMICREORDERING TRUE
#define DEFAULT_DELAFTERADD      TRUE

/*
 * Data structures
 */

/** IIS data */
struct SCIP_IISData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random generator for sorting constraints */
   SCIP_Real             timelimperiter;     /**< time limit per individual solve call */
   SCIP_Bool             minify;             /**< whether the computed IS should undergo a final deletion round to ensure an IIS */
   SCIP_Bool             additive;           /**< whether an additive approach instead of deletion based approach should be used */
   SCIP_Bool             conservative;       /**< should a node or time limit solve be counted as feasible when deleting constraints */
   SCIP_Bool             dynamicreordering;  /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */
   SCIP_Bool             delafteradd;        /**< should the deletion routine be performed after the addition routine (in the case of additive) */
   SCIP_Longint          maxnnodesperiter;   /**< maximum number of nodes per individual solve call */
   int                   batchsize;          /**< the number of constraints to delete or add per iteration */
};

/*
 * Local methods
 */

/** solve subproblem for deletionFilter */
static
SCIP_RETCODE deletionFilterConsSubproblem(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             silent,             /**< run silently? */
   SCIP_Bool*            valid,              /**< Whether the returned subscip is a valid (I)IS */
   SCIP_Real*            timelim,            /**< The global time limit on the IIS call */
   SCIP_Longint*         nodelim,            /**< The global node limit on the IIS call */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          maxnnodesperiter,   /**< maximum number of nodes per individual solve call */
   SCIP_Bool             conservative,       /**< whether we treat a subproblem to be feasible, if it reached its node limt */
   int                   ndelconss,          /**< the number of constraints that will be deleted */
   int*                  idxs,               /**< the indices of the constraints (in the conss array) that will be deleted */
   SCIP_CONS**           conss,              /**< the array of constraints (may be a superset of the current constraints) */
   SCIP_Bool*            stop,               /**< pointer to store whether we have to stop */
   SCIP_Bool*            alldeletionssolved  /**< pointer to store whether all the subscips solved */
   )
{
   SCIP_RETCODE retcode;
   SCIP_STATUS status;
   int j;

   assert( stop != NULL );

   *stop = FALSE;

   /* remove constraints */
   for (j = 0; j < ndelconss; ++j)
   {
      assert( SCIPconsIsInProb(conss[idxs[j]]) );
      SCIP_CALL( SCIPcaptureCons(scip, conss[idxs[j]]) );
      SCIP_CALL( SCIPdelCons(scip, conss[idxs[j]]) );
   }

   /* solve problem until first solution is found or infeasibility has been proven */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", MIN(*timelim, timelimperiter)) );
   if( *nodelim == -1 && maxnnodesperiter == -1 )
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", *nodelim) );
   else if( *nodelim == -1 || maxnnodesperiter == -1 )
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", MAX(*nodelim, maxnnodesperiter)) );
   else
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", MIN(*nodelim, maxnnodesperiter)) );
   retcode = SCIPsolve(scip);
   *timelim -= SCIPgetSolvingTime(scip);
   if( *nodelim != -1 )
   {
      *nodelim -= SCIPgetNTotalNodes(scip);
      if( *nodelim == -1 )
         *nodelim = -2;
   }
   if( retcode != SCIP_OKAY )
   {
      SCIP_CALL( SCIPfreeTransform(scip) );
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Error in sub-scip. Re-add deleted constraints. \n");
      for( j = 0; j < ndelconss; j++ )
      {
         SCIP_CALL( SCIPaddCons(scip, conss[idxs[j]]) );
         SCIP_CALL( SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      *alldeletionssolved = FALSE;
      return SCIP_OKAY;
   }
   status = SCIPgetStatus(scip);

   /* free transform */
   SCIP_CALL( SCIPfreeTransform(scip) );

   /* check status */
   switch ( status )
   {
   case SCIP_STATUS_USERINTERRUPT:    /* if an user interrupt occurred, just stop */
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "User interrupt. Stopping. \n");
      for( j = 0; j < ndelconss; j++ )
      {
         SCIP_CALL( SCIPaddCons(scip, conss[idxs[j]]) );
         SCIP_CALL( SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      *stop = TRUE;
      *alldeletionssolved = FALSE;
      break;

   case SCIP_STATUS_TIMELIMIT:        /* if we reached some limit */
   case SCIP_STATUS_NODELIMIT:
   case SCIP_STATUS_TOTALNODELIMIT:
   case SCIP_STATUS_STALLNODELIMIT:
   case SCIP_STATUS_MEMLIMIT:
   case SCIP_STATUS_GAPLIMIT:
   case SCIP_STATUS_RESTARTLIMIT:
      *alldeletionssolved = FALSE;
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Some limit reached. Removing batch if set as non-conservative. \n");
      if( !conservative )
         *valid = FALSE;
      for( j = 0; j < ndelconss; j++ )
      {
         if( conservative )
            SCIP_CALL( SCIPaddCons(scip, conss[idxs[j]]) );
         SCIP_CALL(SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      break;

   case SCIP_STATUS_INFEASIBLE:       /* if the problem is infeasible */
   case SCIP_STATUS_INFORUNBD:
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Subproblem infeasible. Remove constraint batch.\n");
      *valid = TRUE;
      for( j = 0; j < ndelconss; j++ )
      {
         SCIP_CALL(SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      break;

   case SCIP_STATUS_BESTSOLLIMIT:     /* we found a solution */
   case SCIP_STATUS_SOLLIMIT:
   case SCIP_STATUS_OPTIMAL:
   case SCIP_STATUS_UNBOUNDED:
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Found solution. Keeping constraint batch\n");
      for( j = 0; j < ndelconss; j++ )
      {
         SCIP_CALL( SCIPaddCons(scip, conss[idxs[j]]) );
         SCIP_CALL( SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      break;

   case SCIP_STATUS_UNKNOWN:
   default:
      *alldeletionssolved = FALSE;
      SCIPerrorMessage("unexpected return status %d. Exiting ...\n", status);
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** solve subproblem for additionFilter */
static
SCIP_RETCODE additionFilterConsSubproblem(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             silent,             /**< run silently? */
   SCIP_Real*            timelim,            /**< The global time limit on the IIS call */
   SCIP_Longint*         nodelim,            /**< The global node limit on the IIS call */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          maxnnodesperiter,   /**< maximum number of nodes per individual solve call */
   SCIP_Bool*            feasible,           /**< pointer to store whether the problem is now feasible */
   SCIP_Bool*            stop                /**< pointer to store whether we have to stop */
)
{
   SCIP_RETCODE retcode;
   SCIP_STATUS status;
   
   assert( stop != NULL );
   
   *stop = FALSE;
   *feasible = FALSE;
   
   /* solve problem until first solution is found or infeasibility has been proven */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", MIN(*timelim, timelimperiter)) );
   if( *nodelim == -1 && maxnnodesperiter == -1 )
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", *nodelim) );
   else if( *nodelim == -1 || maxnnodesperiter == -1 )
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", MAX(*nodelim, maxnnodesperiter)) );
   else
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", MIN(*nodelim, maxnnodesperiter)) );
   retcode = SCIPsolve(scip);
   if( retcode != SCIP_OKAY )
   {
      if ( !silent )
         SCIPinfoMessage(scip, NULL, "Error in sub-scip. Be safe and keep added constraints. \n");
      return SCIP_ERROR;
   }
   *timelim -= SCIPgetSolvingTime(scip);
   if( *nodelim != -1 )
   {
      *nodelim -= SCIPgetNTotalNodes(scip);
      if( *nodelim == -1 )
         *nodelim = -2;
   }
   status = SCIPgetStatus(scip);
   
   /* check status */
   switch ( status )
   {
      case SCIP_STATUS_TIMELIMIT:        /* if we reached the time limit, then stop */
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "Time limit exceeded. Added conss failed to induce infeasibility.\n");
         break;
      
      case SCIP_STATUS_USERINTERRUPT:    /* if an user interrupt occurred, just stop */
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "User interrupt. Stopping. \n");
         *stop = TRUE;
         break;
      
      case SCIP_STATUS_NODELIMIT:        /* if we reached some limit */
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "Some limit reached. Added batch failed to induce infeasibility. \n");
         break;
      
      case SCIP_STATUS_INFEASIBLE:       /* if the problem is infeasible */
      case SCIP_STATUS_INFORUNBD:
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "Subproblem infeasible. Final batch of constraints added.\n");
         break;
      
      case SCIP_STATUS_BESTSOLLIMIT:     /* we found a solution */
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_UNBOUNDED:
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "Found solution. Keep adding constraint batches\n");
         *feasible = TRUE;
         break;
      
      case SCIP_STATUS_UNKNOWN:
      default:
         SCIPerrorMessage("unexpected return status %d. Exiting ...\n", status);
         return SCIP_ERROR;
   }
   
   return SCIP_OKAY;
}


/** deletion filter to greedily remove constraints to obtain an (I)IS -- detailed function call */
static
SCIP_RETCODE deletionFilterBatchCons(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool*            valid,              /**< Whether the returned subscip is a valid (I)IS */
   SCIP_Real*            timelim,            /**< The global time limit on the IIS call */
   SCIP_Longint*         nodelim,            /**< The global node limit on the IIS call */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          maxnnodesperiter,   /**< maximum number of nodes per individual solve call */
   SCIP_Bool             conservative,       /**< should a node or time limit solve be counted as feasible when deleting constraints */
   int                   batchsize,          /**< the number of constraints to delete or add per iteration */
   SCIP_Bool*            alldeletionssolved  /**< pointer to store whether all the subscips solved */
   )
{
   SCIP_CONS** conss;
   int* order;
   int* idxs;
   int nconss;
   int nbatches;
   int i;
   
   /* Get constraint information */
   nconss = SCIPgetNOrigConss(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conss, nconss) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &conss, SCIPgetOrigConss(scip), nconss) );

   /* reset problem */
   SCIP_CALL( SCIPfreeTransform(scip) );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

   /* prepare random order */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nconss) );
   for (i = 0; i < nconss; ++i)
      order[i] = i;
   SCIPrandomPermuteIntArray(randnumgen, order, 0, nconss);
   
   /* Calculate the number of batches */
   nbatches = (int) SCIPceil(scip, (SCIP_Real) nconss / batchsize);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idxs, batchsize) );

   /* Loop through all batches of constraints in random order */
   for (i = 0; i < nbatches; ++i)
   {
      int minconsidx;
      int maxconsidx;
      int j;
      int ndelconss;
      SCIP_Bool stopiter = FALSE;
   
      minconsidx = i * batchsize;
      if( minconsidx >= nconss )
         break;
      maxconsidx = MIN((i + 1) * batchsize, nconss);
      ndelconss = maxconsidx - minconsidx;
      for( j = 0; j < ndelconss; j++ )
      {
         idxs[j] = order[i*batchsize+j];
      }

      /* treat subproblem */
      SCIP_CALL( deletionFilterConsSubproblem(scip, silent, valid, timelim, nodelim, timelimperiter,
                                              maxnnodesperiter, conservative, ndelconss, idxs, conss,
                                              &stopiter, alldeletionssolved) );
      if( *timelim <= 0 || *nodelim <= -2 )
         break;
      
      assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

      if ( stopiter )
         break;
   }
   
   SCIPfreeBlockMemoryArray(scip, &idxs, batchsize);
   SCIPfreeBlockMemoryArray(scip, &order, nconss);
   SCIPfreeBlockMemoryArray(scip, &conss, nconss);

   return SCIP_OKAY;
}

/** addition filter to greedily add constraints to obtain an (I)IS -- detailed function call */
static
SCIP_RETCODE additionFilterBatchCons(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool*            valid,              /**< Whether the returned subscip is a valid (I)IS */
   SCIP_Real*            timelim,            /**< The global time limit on the IIS call */
   SCIP_Longint*         nodelim,            /**< The global node limit on the IIS call */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          maxnnodesperiter,   /**< maximum number of nodes per individual solve call */
   SCIP_Bool             dynamicreordering,  /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */
   int                   batchsize           /**< the number of constraints to delete or add per iteration */
   )
{
   SCIP_CONS** conss;
   SCIP_SOL* sol;
   SCIP_SOL* copysol;
   SCIP_Bool* inIS;
   int* order;
   SCIP_Bool feasible;
   SCIP_Bool stopiter;
   SCIP_RETCODE retcode;
   SCIP_RESULT result;
   int nconss;
   int i;
   int j;
   int k;
   
   assert( *valid == TRUE );
   
   /* Get constraint information */
   nconss = SCIPgetNOrigConss(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conss, nconss) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &conss, SCIPgetOrigConss(scip), nconss) );
   
   /* Initialise information for whether a constraint is in the final infeasible system */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &inIS, nconss) );
   
   /* First capture and then delete all constraints */
   assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
   for( i = 0; i < nconss; ++i )
   {
      assert( SCIPconsIsInProb(conss[i]) );
      SCIP_CALL( SCIPcaptureCons(scip, conss[i]) );
      SCIP_CALL( SCIPdelCons(scip, conss[i]) );
      inIS[i] = FALSE;
   }
   *valid = FALSE;
   
   /* Prepare random order in which the constraints will be added back */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nconss) );
   for (i = 0; i < nconss; ++i)
      order[i] = i;
   SCIPrandomPermuteIntArray(randnumgen, order, 0, nconss);
   
   /* Continue to add constraints until an infeasible status is reached */
   i = 0;
   feasible = TRUE;
   stopiter = FALSE;
   result = SCIP_FEASIBLE;
   while( feasible && !stopiter )
   {
      /* Add the next batch of constraints */
      k = 0;
      for( j = i; j < nconss; ++j )
      {
         if( k >= batchsize )
            break;
         if( !inIS[order[j]] )
         {
            SCIP_CALL( SCIPaddCons(scip, conss[order[j]]) );
            SCIP_CALL( SCIPreleaseCons(scip, &conss[order[j]]) );
            inIS[order[j]] = TRUE;
            ++k;
         }
      }
      if( k == 0 )
         break;
      i = i + k;
      
      /* Solve the reduced problem */
      retcode = additionFilterConsSubproblem(scip, silent, timelim, nodelim, timelimperiter, maxnnodesperiter, &feasible, &stopiter);
      if( *timelim <= 0 || *nodelim <= -2 )
      {
         SCIP_CALL( SCIPfreeTransform(scip) );
         break;
      }
      
      if( retcode == SCIP_OKAY )
      {
         /* free transform and copy solution if there is one */
         sol = SCIPgetBestSol(scip);
         if( sol != NULL )
         {
            SCIP_CALL( SCIPcreateSolCopyOrig(scip, &copysol, sol) );
            SCIP_CALL( SCIPunlinkSol(scip, copysol) );
         }
         SCIP_CALL( SCIPfreeTransform(scip) );
         assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
   
         /* Add any other constraints that are also feasible for the current solution */
         if( feasible && !stopiter && copysol != NULL && dynamicreordering )
         {
            k = 0;
            for( j = i; j < nconss; ++j )
            {
               if( !inIS[order[j]] )
               {
                  SCIP_CALL( SCIPcheckCons(scip, conss[order[j]], copysol, FALSE, FALSE, FALSE, &result) );
                  if( result == SCIP_FEASIBLE )
                  {
                     SCIP_CALL( SCIPaddCons(scip, conss[order[j]]) );
                     SCIP_CALL( SCIPreleaseCons(scip, &conss[order[j]]) );
                     inIS[order[j]] = TRUE;
                     k++;
                  }
               }
            }
            if( k > 0 )
               SCIPinfoMessage(scip, NULL, "Added %d many constraints dynamically.\n", k);
         }
   
         if ( stopiter )
            break;
      }
      else
      {
         SCIP_CALL( SCIPfreeTransform(scip) );
         assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
         if ( stopiter )
            break;
      }
      
   }
   
   /* Release any cons not in the IS */
   for( i = 0; i < nconss; ++i )
   {
      if( !inIS[order[i]] )
         SCIP_CALL( SCIPreleaseCons(scip, &conss[order[i]]) );
         
   }
   
   SCIPfreeBlockMemoryArray(scip, &order, nconss);
   SCIPfreeBlockMemoryArray(scip, &inIS, nconss);
   SCIPfreeBlockMemoryArray(scip, &conss, nconss);
   *valid = TRUE;
   
   return SCIP_OKAY;
}

/*
 * Callback methods of IIS
 */


/** copy method for IIS plugin (called when SCIP copies plugins) */
static
SCIP_DECL_IISCOPY(iisCopyGreedy)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(iis != NULL);
   assert(strcmp(SCIPiisGetName(iis), IIS_NAME) == 0);
   
   /* call inclusion method of IIS */
   SCIP_CALL( SCIPincludeIISGreedy(scip) );
   
   return SCIP_OKAY;
}

/** destructor of IIS to free user data (called when SCIP is exiting) */
/**! [SnippetIISFreeGreedy] */
static
SCIP_DECL_IISFREE(iisFreeGreedy)
{  /*lint --e{715}*/
   SCIP_IISDATA* iisdata;
   
   iisdata = SCIPiisGetData(iis);
   SCIPfreeRandom(scip, &(iisdata->randnumgen));
   
   SCIPfreeBlockMemory(scip, &iisdata);
   
   SCIPiisSetData(iis, NULL);
   
   return SCIP_OKAY;
}
/**! [SnippetIISFreeGreedy] */

/** IIS generation method of IIS */
static
SCIP_DECL_IISGENERATE(iisGenerateGreedy)
{  /*lint --e{715}*/
   struct SCIP_IISData* iisdata;
   
   assert(iis != NULL);
   assert(result != NULL);
   
   *result = SCIP_SUCCESS;
   
   iisdata = SCIPiisGetData(iis);
   assert(iisdata != NULL);
   
   SCIP_CALL( SCIPgenerateIISGreedy(scip, valid, irreducible, timelim, nodelim, silent, iisdata->randnumgen,
                                    iisdata->timelimperiter, iisdata->additive, iisdata->conservative,
                                    iisdata->dynamicreordering, iisdata->delafteradd,
                                    iisdata->maxnnodesperiter, iisdata->batchsize, result) );
   
   return SCIP_OKAY;
}


/*
 * IIS specific interface methods
 */

/** creates the greedy IIS and includes it in SCIP */
SCIP_RETCODE SCIPincludeIISGreedy(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_IISDATA* iisdata;
   SCIP_IIS* iis;
   
   /* create greedy IIS data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &iisdata) );
   BMSclearMemory(iisdata);
   SCIP_CALL( SCIPcreateRandom(scip, &(iisdata)->randnumgen, RANDSEED, TRUE) );
   
   SCIP_CALL( SCIPincludeIISBasic(scip, &iis, IIS_NAME, IIS_DESC, IIS_PRIORITY, iisGenerateGreedy,
                                     iisdata) );
   
   assert(iis != NULL);
   
   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetIISCopy(scip, iis, iisCopyGreedy) );
   SCIP_CALL( SCIPsetIISFree(scip, iis, iisFreeGreedy) );
   
   /* add greedy IIS parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "iis/" IIS_NAME "/batchsize",
         "batch size of constraints that are removed or added in one iteration",
         &iisdata->batchsize, FALSE, DEFAULT_BATCHSIZE, 1, INT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IIS_NAME "/timelimperiter",
         "time limit of optimization process for each individual subproblem",
         &iisdata->timelimperiter, FALSE, DEFAULT_TIMELIMPERITER, 0.0, SCIP_INVALID/10.0, NULL, NULL) );
   
   SCIP_CALL( SCIPaddLongintParam(scip,
         "iis/" IIS_NAME "/maxnnodesperiter",
         "node limit of optimization process for each individual subproblem",
         &iisdata->maxnnodesperiter, FALSE, DEFAULT_MAXNNODESPERITER, 1, INT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IIS_NAME "/additive",
         "should an additive constraint approach be used instead of deletion",
         &iisdata->additive, FALSE, DEFAULT_ADDITIVE, NULL, NULL) );
   
   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IIS_NAME "/conservative",
         "should a hit limit (e.g. node  time) solve be counted as feasible when deleting constraints",
         &iisdata->conservative, TRUE, DEFAULT_CONSERVATIVE, NULL, NULL) );
   
   SCIP_CALL( SCIPaddBoolParam(scip,
      "iis/" IIS_NAME "/dynamicreordering",
      "should satisfied constraints outside the batch of an intermediate solve be added during the additive method",
      &iisdata->dynamicreordering, TRUE, DEFAULT_DYNAMICREORDERING, NULL, NULL) );
   
   SCIP_CALL( SCIPaddBoolParam(scip,
      "iis/" IIS_NAME "/delafteradd",
      "should the deletion routine be performed after the addition routine (in the case of additive)",
      &iisdata->delafteradd, TRUE, DEFAULT_DELAFTERADD, NULL, NULL) );
   
   return SCIP_OKAY;
}

/** perform a greedy addition or deletion algorithm to obtain an infeasible subsystem (IS).
 *
 *  This is the generation method for the greedy IIS rule.
 *  Depending on the parameter choices, constraints are either greedily added from an empty problem,
 *  or deleted from a complete problem. In the case of constraints being added, this is done until the problem
 *  becomes infeasible, after which one can then begin deleting constraints. In the case of deleting constraints,
 *  this is done until no more constraints (or batches of constraints) can be deleted without making
 *  the problem feasible.
 */
SCIP_RETCODE SCIPgenerateIISGreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            valid,              /**< Whether the returned subscip is a valid (I)IS */
   SCIP_Bool*            irreducible,        /**< Whether the returned subscip is a minimal IIS */
   SCIP_Real*            timelim,            /**< The global time limit on the IIS call */
   SCIP_Longint*         nodelim,            /**< The global node limit on the IIS call */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Bool             additive,           /**< whether an additive approach instead of deletion based approach should be used */
   SCIP_Bool             conservative,       /**< should a hit limit (e.g. node / time) solve be counted as feasible when deleting constraints */
   SCIP_Bool             dynamicreordering,  /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */
   SCIP_Bool             delafteradd,        /**< should the deletion routine be performed after the addition routine (in the case of additive) */
   SCIP_Longint          maxnnodesperiter,   /**< maximum number of nodes per individual solve call */
   int                   batchsize,          /**< the number of constraints to delete or add per iteration */
   SCIP_RESULT*          result              /**< pointer to store the result os the IIS run */
   )
{
   SCIP_Bool alldeletionssolved = TRUE;
   
   if( additive )
   {
      if( !silent )
         SCIPinfoMessage(scip, NULL, "Starting greedy addition algorithm\n");
      SCIP_CALL( additionFilterBatchCons(scip, valid, timelim, nodelim, silent, randnumgen, timelimperiter, maxnnodesperiter, dynamicreordering, batchsize) );
      *irreducible = FALSE;
      if( *timelim <= 0 || *nodelim <= -2 )
      {
         *result = SCIP_SUCCESS;
         return SCIP_OKAY;
      }
   }
   else
   {
      if( !silent )
         SCIPinfoMessage(scip, NULL, "Starting greedy deletion algorithm\n");
      SCIP_CALL( deletionFilterBatchCons(scip, valid, timelim, nodelim, silent, randnumgen, timelimperiter, maxnnodesperiter, conservative, batchsize, &alldeletionssolved) );
      if( alldeletionssolved && batchsize == 1 )
         *irreducible = TRUE;
      if( *timelim <= 0 || *nodelim <= -2 )
      {
         *result = SCIP_SUCCESS;
         return SCIP_OKAY;
      }
   }
   
   if( delafteradd && additive )
   {
      if( !silent )
         SCIPinfoMessage(scip, NULL, "Starting greedy deletion algorithm on reduced problem\n");
      SCIP_CALL( deletionFilterBatchCons(scip, valid, timelim, nodelim, silent, randnumgen, timelimperiter, maxnnodesperiter, conservative, batchsize, &alldeletionssolved) );
      if( alldeletionssolved && batchsize == 1 )
         *irreducible = TRUE;
      if( *timelim <= 0 || *nodelim <= -2 )
      {
         *result = SCIP_SUCCESS;
         return SCIP_OKAY;
      }
   }
   
   return SCIP_OKAY;
}
