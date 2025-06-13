/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   iisfinder_greedy.c
 * @brief  greedy deletion and addition filter heuristic to compute (I)ISs
 * @author Marc Pfetsch
 * @author Mark Turner
 * @author Paul Meinhold
 */

#include <assert.h>

#include "scip/iisfinder_greedy.h"

#define IISFINDER_NAME           "greedy"
#define IISFINDER_DESC           "greedy deletion or addition constraint deletion"
#define IISFINDER_PRIORITY        8000

#define DEFAULT_TIMELIMPERITER   1e+20 /**< time limit per individual solve call */
#define DEFAULT_NODELIMPERITER   -1L   /**< maximum number of nodes per individual solve call */

#define DEFAULT_ADDITIVE         TRUE  /**< whether an additive approach instead of deletion based approach should be used */
#define DEFAULT_CONSERVATIVE     TRUE  /**< should a solve that reached some limit be counted as feasible when deleting constraints */
#define DEFAULT_DYNAMICREORDERING TRUE /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */
#define DEFAULT_DELAFTERADD      TRUE  /**< should the deletion routine be performed after the addition routine (in the case of additive) */

#define DEFAULT_INITBATCHSIZE    16    /**< the initial batchsize for the first iteration */
#define DEFAULT_INITRELBATCHSIZE 0.03125 /**< the initial batchsize relative to the original problem for the first iteration */
#define DEFAULT_MAXBATCHSIZE     INT_MAX /**< the maximum batchsize per iteration */
#define DEFAULT_MAXRELBATCHSIZE  0.5   /**< the maximum batchsize relative to the original problem per iteration */
#define DEFAULT_BATCHINGFACTOR   2.0   /**< the factor with which the batchsize is multiplied in every update */
#define DEFAULT_BATCHINGOFFSET   0.0   /**< the offset which is added to the multiplied batchsize in every update */
#define DEFAULT_BATCHUPDATEINTERVAL 1  /**< the number of iterations to run with a constant batchsize before updating (1: always update) */


/*
 * Data structures
 */

/** IIS finder data */
struct SCIP_IISfinderData
{
   SCIP_Real             timelimperiter;     /**< time limit per individual solve call */
   SCIP_Longint          nodelimperiter;     /**< maximum number of nodes per individual solve call */

   SCIP_Bool             additive;           /**< whether an additive approach instead of deletion based approach should be used */
   SCIP_Bool             conservative;       /**< should a solve that reached some limit be counted as feasible when deleting constraints */
   SCIP_Bool             delafteradd;        /**< should the deletion routine be performed after the addition routine (in the case of additive) */
   SCIP_Bool             dynamicreordering;  /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */

   int                   initbatchsize;      /**< the initial batchsize for the first iteration */
   SCIP_Real             initrelbatchsize;   /**< the initial batchsize relative to the original problem for the first iteration */
   int                   maxbatchsize;       /**< the maximum batchsize per iteration */
   SCIP_Real             maxrelbatchsize;    /**< the maximum batchsize relative to the original problem per iteration */
   SCIP_Real             batchingfactor;     /**< the factor with which the batchsize is multiplied in every update */
   SCIP_Real             batchingoffset;     /**< the offset which is added to the multiplied batchsize in every update */
   int                   batchupdateinterval; /**< the number of iterations to run with a constant batchsize before updating (1: always update) */
};

/*
 * Local methods
 */

/* Set time and node limits on the subproblem */
static
SCIP_RETCODE setLimits(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip  */
   SCIP_Real             timelim,            /**< total time limit allowed of the whole call */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelim,            /**< node limit allowed for the whole call */
   SCIP_Longint          nodelimperiter      /**< maximum number of nodes per individual solve call */
   )
{
   SCIP_Real currtime;
   SCIP_Real mintimelim;
   SCIP_Longint globalnodelim;

   /* Set the time limit for the solve call. Take into account the global time limit, the current time used, and the time lim on each individual call */
   currtime = SCIPiisGetTime(iis);
   mintimelim = MIN(timelim - currtime, timelimperiter);
   mintimelim = MAX(mintimelim, 0);
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", mintimelim) );

   /* Set the node limit for the solve call. Take into account the global node limit, the current nodes processed, and the node lim on each individual call */
   if( nodelim == -1 )
      globalnodelim = -1;
   else
   {
      globalnodelim = nodelim - SCIPiisGetNNodes(iis);
      assert( globalnodelim >= 0 );
   }
   if( globalnodelim == -1 && nodelimperiter == -1 )
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", nodelim) );
   }
   else if( globalnodelim == -1 || nodelimperiter == -1 )
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", MAX(globalnodelim, nodelimperiter)) );
   }
   else
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", MIN(globalnodelim, nodelimperiter)) );
   }
   return SCIP_OKAY;
}

/* Revert bound changes made in the subproblem */
static
SCIP_RETCODE revertBndChgs(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_VAR**            vars,               /**< the array of variables whose bounds are changed */
   SCIP_Real*            bounds,             /**< the array of original bounds for the variables */
   int*                  idxs,               /**< the indices of the vars (in the vars array) that have been deleted */
   int                   ndelbounds,         /**< the number of bounds that will be deleted */
   SCIP_Bool             islb                /**< are the bounds that are being deleted LBs? */
   )
{
   int i;

   /* Iterate over the affected variables and restore the bounds back to their original values */
   for (i = 0; i < ndelbounds; ++i)
   {
      if( islb )
      {
         if( !SCIPisInfinity(scip, -bounds[i]) )
            SCIP_CALL( SCIPchgVarLb(scip, vars[idxs[i]], bounds[i]) );
      }
      else
      {
         if( !SCIPisInfinity(scip, bounds[i]) )
            SCIP_CALL( SCIPchgVarUb(scip, vars[idxs[i]], bounds[i]) );
      }
   }
   return SCIP_OKAY;
}

/* Revert deleted constraint changes made in the subproblem */
static
SCIP_RETCODE revertConssDeletions(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_CONS**           conss,              /**< the array of constraints where some have been deleted */
   int*                  idxs,               /**< the indices of the cons (in the conss array) that have been deleted */
   int                   ndelconss,          /**< the number of constraints that have been deleted */
   SCIP_Bool             addconss,           /**< Should the constraints be added back */
   SCIP_Bool             keepptrs            /**< Should the constraint pointers be kept for reuse */
   )
{
   int i;
   SCIP_CONS* copycons;

   assert( addconss || !keepptrs );

   for( i = 0; i < ndelconss; ++i )
   {
      if( addconss )
      {
         SCIP_CALL( SCIPaddCons(scip, conss[idxs[i]]) );
      }
      if( keepptrs )
      {
         copycons = conss[idxs[i]];
         assert(SCIPconsGetNUses(copycons) > 1);
         SCIP_CALL( SCIPreleaseCons(scip, &copycons) );
      }
      else
      {
         SCIP_CALL( SCIPreleaseCons(scip, &conss[idxs[i]]) );
      }
   }

   return SCIP_OKAY;
}

/* Update the batchsize accoring to the chosen rule */
static
SCIP_RETCODE updateBatchsize(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   initbatchsize,      /**< the initial batchsize */
   int                   maxbatchsize,       /**< the maximum batchsize per iteration */
   int                   iteration,          /**< the current iteration */
   SCIP_Bool             resettoinit,        /**< should the batchsize be reset to the initial batchsize? */
   SCIP_Real             batchingfactor,     /**< the factor with which the batchsize is multiplied in every update */
   SCIP_Real             batchingoffset,     /**< the offset which is added to the multiplied batchsize in every update */
   int                   batchupdateinterval, /**< the number of iterations to run with a constant batchsize before updating (1: always update) */
   int*                  batchsize           /**< the batchsize to be updated */
   )
{
   if( resettoinit )
      *batchsize = initbatchsize;
   else if( iteration % batchupdateinterval == 0 )
      *batchsize = (int)ceil(batchingfactor * (*batchsize) + batchingoffset);

   /* respect limits and maximum */
   *batchsize = MIN(*batchsize, maxbatchsize);
   *batchsize = MAX(*batchsize, 1);
   SCIPdebugMsg(scip, "Updated batchsize to %d\n", *batchsize);

   return SCIP_OKAY;
}

/** solve subproblem for deletionFilter */
static
SCIP_RETCODE deletionSubproblem(
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip */
   SCIP_CONS**           conss,              /**< The array of constraints (may be a superset of the current constraints) */
   SCIP_VAR**            vars,               /**< the array of vars */
   int*                  idxs,               /**< the indices of the constraints / vars that will be deleted / bounds removed */
   int                   ndels,              /**< the number of bounds that will be deleted */
   SCIP_Real             timelim,            /**< The global time limit on the IIS finder call */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS finder call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool             conservative,       /**< whether we treat a subproblem to be feasible, if it reaches some status that could result in infeasible, e.g. node limit */
   SCIP_Bool             delbounds,          /**< whether bounds should be deleted instead of constraints */
   SCIP_Bool             islb,               /**< are the bounds that are being deleted LBs? */
   SCIP_Bool*            deleted,            /**< have the deleted bounds or constraints stayed deleted */
   SCIP_Bool*            stop,               /**< pointer to store whether we have to stop */
   SCIP_Bool*            alldeletionssolved  /**< pointer to store whether all the subscips solved */
   )
{
   SCIP* scip;
   SCIP_Real* bounds = NULL;
   SCIP_RETCODE retcode;
   SCIP_STATUS status;
   SCIP_Bool chgmade = FALSE;
   int i;

   *deleted = FALSE;
   *stop = FALSE;
   scip = SCIPiisGetSubscip(iis);

   /* remove bounds or constraints */
   if( delbounds )
   {
      assert(vars != NULL);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &bounds, ndels) );
      for (i = 0; i < ndels; ++i)
      {
         if( islb )
         {
            bounds[i] = SCIPvarGetLbOriginal(vars[idxs[i]]);
            if( !SCIPisInfinity(scip, -bounds[i]) )
            {
               SCIP_CALL( SCIPchgVarLb(scip, vars[idxs[i]], -SCIPinfinity(scip)) );
               chgmade = TRUE;
            }
         }
         else
         {
            bounds[i] = SCIPvarGetUbOriginal(vars[idxs[i]]);
            if( !SCIPisInfinity(scip, bounds[i]) )
            {
               SCIP_CALL( SCIPchgVarUb(scip, vars[idxs[i]], SCIPinfinity(scip)) );
               chgmade = TRUE;
            }
         }
      }
   }
   else
   {
      assert(conss != NULL);

      if( ndels > 0 )
         chgmade = TRUE;
      for (i = 0; i < ndels; ++i)
      {
         assert( SCIPconsIsInProb(conss[idxs[i]]) );
         SCIP_CALL( SCIPcaptureCons(scip, conss[idxs[i]]) );
         SCIP_CALL( SCIPdelCons(scip, conss[idxs[i]]) );
      }
   }

   if( !chgmade )
   {
      if( delbounds )
         SCIPfreeBlockMemoryArray(scip, &bounds, ndels);
      return SCIP_OKAY;
   }

   /* solve problem until first solution is found or infeasibility has been proven */
   SCIP_CALL( setLimits(scip, iis, timelim, timelimperiter, nodelim, nodelimperiter) );
   retcode = SCIPsolve(scip);
   SCIPiisAddNNodes(iis, SCIPgetNTotalNodes(scip));

   if( retcode != SCIP_OKAY )
   {
      SCIP_CALL( SCIPfreeTransform(scip) );
      SCIPdebugMsg(scip, "Error in sub-scip with deleted constraints / bounds. Re-adding them.\n");
      if( delbounds )
      {
         SCIP_CALL( revertBndChgs(scip, vars, bounds, idxs, ndels, islb) );
         SCIPfreeBlockMemoryArray(scip, &bounds, ndels);
      }
      else
      {
         SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, TRUE, TRUE) );
      }
      *alldeletionssolved = FALSE;
      return SCIP_OKAY;
   }

   status = SCIPgetStatus(scip);
   SCIP_CALL( SCIPfreeTransform(scip) );

   /* check status and handle accordingly */
   switch ( status )
   {
      case SCIP_STATUS_USERINTERRUPT:    /* if a user interrupt occurred, just stop */
      case SCIP_STATUS_TERMINATE:
         SCIPdebugMsg(scip, "User interrupt. Stopping.\n");
         if( delbounds )
         {
            SCIP_CALL( revertBndChgs(scip, vars, bounds, idxs, ndels, islb) );
         }
         else
         {
            SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, TRUE, FALSE) );
         }
         *stop = TRUE;
         *alldeletionssolved = FALSE;
         break;

      case SCIP_STATUS_TIMELIMIT:        /* if we reached some status that may have ended up in an infeasible problem */
      case SCIP_STATUS_NODELIMIT:
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_DUALLIMIT:
         *alldeletionssolved = FALSE;
         SCIPdebugMsg(scip, "Some limit reached. Keeping bounds / constraints removed if non-conservative.\n");
         if( !conservative )
         {
            SCIPiisSetSubscipInfeasible(iis, FALSE);
            *deleted = TRUE;
         }
         if( conservative && delbounds )
         {
            SCIP_CALL( revertBndChgs(scip, vars, bounds, idxs, ndels, islb) );
         }
         if( !delbounds )
         {
            SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, conservative, conservative) );
         }
         break;

      case SCIP_STATUS_INFEASIBLE:       /* if the problem is infeasible */
         SCIPdebugMsg(scip, "Subproblem with bounds / constraints removed infeasible. Keep them removed.\n");
         SCIPiisSetSubscipInfeasible(iis, TRUE);
         if( !delbounds )
         {
            SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, FALSE, FALSE) );
         }
         *deleted = TRUE;
         break;

      case SCIP_STATUS_BESTSOLLIMIT:     /* we found a solution */
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_UNBOUNDED:
      case SCIP_STATUS_PRIMALLIMIT:
         SCIPdebugMsg(scip, "Found solution to subproblem with bounds / constraints removed. Add them back.\n");
         if( delbounds )
         {
            SCIP_CALL( revertBndChgs(scip, vars, bounds, idxs, ndels, islb) );
         }
         else
         {
            SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, TRUE, TRUE) );
         }
         break;

      case SCIP_STATUS_UNKNOWN:
      default:
         *alldeletionssolved = FALSE;
         SCIPerrorMessage("Unexpected return status %d in removed bounds subproblem. Exiting...\n", status);
         if( delbounds )
            SCIPfreeBlockMemoryArray(scip, &bounds, ndels);
         return SCIP_ERROR;
   }

   if( delbounds )
      SCIPfreeBlockMemoryArray(scip, &bounds, ndels);

   return SCIP_OKAY;
}

/** solve subproblem for additionFilter */
static
SCIP_RETCODE additionSubproblem(
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip  */
   SCIP_Real             timelim,            /**< The global time limit on the IIS finder call */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS finder call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool*            feasible,           /**< pointer to store whether the problem is feasible */
   SCIP_Bool*            stop                /**< pointer to store whether we have to stop */
   )
{
   SCIP* scip;
   SCIP_RETCODE retcode;
   SCIP_STATUS status;

   assert( stop != NULL );
   scip = SCIPiisGetSubscip(iis);

   *stop = FALSE;

   /* solve problem until first solution is found or infeasibility has been proven */
   SCIP_CALL( setLimits(scip, iis, timelim, timelimperiter, nodelim, nodelimperiter) );
   retcode = SCIPsolve(scip);

   if( retcode != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "Error in sub-scip with added constraints. Keep added constraints.\n");
      return SCIP_ERROR;
   }

   SCIPiisAddNNodes(iis, SCIPgetNTotalNodes(scip));
   status = SCIPgetStatus(scip);

   /* check status */
   switch ( status )
   {
      case SCIP_STATUS_USERINTERRUPT:    /* if an user interrupt occurred, just stop */
      case SCIP_STATUS_TERMINATE:
         SCIPdebugMsg(scip, "User interrupt. Stopping.\n");
         *stop = TRUE;
         break;

      case SCIP_STATUS_NODELIMIT:        /* if we reached some limit */
      case SCIP_STATUS_TIMELIMIT:
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_DUALLIMIT:
         SCIPdebugMsg(scip, "Some limit reached. Added constraint batch failed to induce infeasibility. Continue adding.\n");
         break;

      case SCIP_STATUS_INFEASIBLE:       /* if the problem is infeasible */
         SCIPdebugMsg(scip, "Subproblem with added constraints infeasible. Final batch of constraints added.\n");
         *feasible = FALSE;
         break;

      case SCIP_STATUS_BESTSOLLIMIT:     /* we found a solution */
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_UNBOUNDED:
      case SCIP_STATUS_PRIMALLIMIT:
         SCIPdebugMsg(scip, "Found solution of subproblem with added constraints. Keep adding constraint batches.\n");
         *feasible = TRUE;
         break;

      case SCIP_STATUS_UNKNOWN:
      default:
         SCIPerrorMessage("Unexpected return status %d. Exiting ...\n", status);
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/** Deletion filter to greedily remove constraints to obtain an (I)IS */
static
SCIP_RETCODE deletionFilterBatch(
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip  */
   SCIP_Real             timelim,            /**< The global time limit on the IIS finder call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS finder call */
   SCIP_Bool             removebounds,       /**< Whether the algorithm should remove bounds as well as constraints */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */

   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool             conservative,       /**< should a node or time limit solve be counted as feasible when deleting constraints */

   int                   initbatchsize,      /**< the initial batchsize for the first iteration */
   int                   maxbatchsize,       /**< the maximum batchsize per iteration */
   SCIP_Real             batchingfactor,     /**< the factor with which the batchsize is multiplied in every update */
   SCIP_Real             batchingoffset,     /**< the offset which is added to the multiplied batchsize in every update */
   int                   batchupdateinterval, /**< the number of iterations to run with a constant batchsize before updating (1: always update) */

   SCIP_Bool*            alldeletionssolved  /**< pointer to store whether all the subscips solved */
   )
{
   SCIP* scip;
   SCIP_CONS** origconss;
   SCIP_CONS** conss;
   SCIP_VAR** origvars;
   SCIP_VAR** vars;
   SCIP_RANDNUMGEN* randnumgen;
   int* order;
   int* idxs;
   SCIP_Bool stopiter;
   SCIP_Bool deleted;
   int nconss;
   int nvars;
   int batchindex;
   int batchsize;
   int iteration;
   int i;
   int k;

   /* get current subscip */
   scip = SCIPiisGetSubscip(iis);
   assert( scip != NULL );
   assert( SCIPiisIsSubscipInfeasible(iis) );

   /* get random generator */
   randnumgen = SCIPiisGetRandnumgen(iis);
   assert( randnumgen != NULL );

   /* get batch size */
   assert( initbatchsize >= 1 );
   assert( maxbatchsize >= 1 );
   initbatchsize = MIN(initbatchsize, maxbatchsize);

   /* allocate indices array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idxs, maxbatchsize) );

   /* get constraint information */
   nconss = SCIPgetNOrigConss(scip);
   origconss = SCIPgetOrigConss(scip);
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &conss, origconss, nconss) );

   /* reset problem */
   SCIP_CALL( SCIPfreeTransform(scip) );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

   /* prepare random order for constraints */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nconss) );
   for (i = 0; i < nconss; ++i)
      order[i] = i;
   SCIPrandomPermuteIntArray(randnumgen, order, 0, nconss);

   /* Loop through all batches of constraints in random order */
   i = 0;
   iteration = 0;
   deleted = FALSE;
   stopiter = FALSE;
   while( i < nconss )
   {
      /* update batchsize */
      SCIP_CALL( updateBatchsize(scip, initbatchsize, maxbatchsize, iteration, !deleted, batchingfactor, batchingoffset, batchupdateinterval, &batchsize) );

      k = 0;
      batchindex = i;
      while( i < nconss && k < batchsize )
      {
         assert( conss[order[i]] != NULL );
         if( SCIPconsGetNUses(conss[order[i]]) == 1 )
         {
            idxs[k] = order[i];
            k++;
         }
         i++;
      }

      /* treat subproblem */
      SCIP_CALL( deletionSubproblem(iis, conss, NULL, idxs, k, timelim, timelimperiter, nodelim, nodelimperiter,
            conservative, FALSE, FALSE, &deleted, &stopiter, alldeletionssolved) );
      if( !silent && deleted )
         SCIPiisfinderInfoMessage(iis, FALSE);

      if( timelim - SCIPiisGetTime(iis) <= 0 || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || stopiter )
         break;

      /* reset i to beginning of current batch if batch has not been deleted and k was large */
      if( !deleted && (k > initbatchsize) )
         i = batchindex;

      ++iteration;

      assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
   }

   SCIPfreeBlockMemoryArray(scip, &order, nconss);
   SCIPfreeBlockMemoryArray(scip, &conss, nconss);
   SCIPfreeBlockMemoryArray(scip, &idxs, maxbatchsize);

   if( timelim - SCIPiisGetTime(iis) <= 0 || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || stopiter )
      return SCIP_OKAY;

   /* Repeat the above procedure but for bounds instead of constraints */
   if( removebounds )
   {
      /* allocate indices array */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idxs, maxbatchsize) );

      nvars = SCIPgetNOrigVars(scip);
      origvars = SCIPgetOrigVars(scip);
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vars, origvars, nvars) );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nvars) );
      for (i = 0; i < nvars; ++i)
         order[i] = i;
      SCIPrandomPermuteIntArray(randnumgen, order, 0, nvars);

      i = 0;
      iteration = 0;
      deleted = FALSE;
      while( i < nvars )
      {
         /* update batchsize */
         SCIP_CALL( updateBatchsize(scip, initbatchsize, maxbatchsize, iteration, !deleted, batchingfactor, batchingoffset, batchupdateinterval, &batchsize) );

         k = 0;
         batchindex = i;
         /* Do not delete bounds of binary variables or bother with calculations of free variables */
         while( i < nvars && k < batchsize )
         {
            if( (SCIPvarGetType(vars[order[i]]) != SCIP_VARTYPE_BINARY) && (!SCIPisInfinity(scip, -SCIPvarGetLbOriginal(vars[order[i]])) || !SCIPisInfinity(scip, SCIPvarGetUbOriginal(vars[order[i]]))) )
            {
               idxs[k] = order[i];
               k++;
            }
            i++;
         }
         if( k == 0 )
            break;

         /* treat subproblem with LB deletions */
         SCIP_CALL( deletionSubproblem(iis, NULL, vars, idxs, k, timelim, timelimperiter, nodelim, nodelimperiter,
               conservative, TRUE, TRUE, &deleted, &stopiter, alldeletionssolved) );
         if( timelim - SCIPiisGetTime(iis) <= 0 || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || stopiter )
            break;
         if( !silent && deleted )
            SCIPiisfinderInfoMessage(iis, FALSE);

         /* treat subproblem with UB deletions */
         SCIP_CALL( deletionSubproblem(iis, NULL, vars, idxs, k, timelim, timelimperiter, nodelim, nodelimperiter,
               conservative, TRUE, FALSE, &deleted, &stopiter, alldeletionssolved) );

         if( timelim - SCIPiisGetTime(iis) <= 0 || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || stopiter )
            break;
         if( !silent && deleted )
            SCIPiisfinderInfoMessage(iis, FALSE);

         /* reset i to beginning of current batch if batch has not been deleted and k was large */
         if( !deleted && (k > initbatchsize) )
            i = batchindex;

         ++iteration;
      }

      SCIPfreeBlockMemoryArray(scip, &order, nvars);
      SCIPfreeBlockMemoryArray(scip, &vars, nvars);
      SCIPfreeBlockMemoryArray(scip, &idxs, maxbatchsize);
   }

   return SCIP_OKAY;
}

/** Addition filter to greedily add constraints to obtain an (I)IS */
static
SCIP_RETCODE additionFilterBatch(
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip  */
   SCIP_Real             timelim,            /**< The global time limit on the IIS finder call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS finder call */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */

   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool             dynamicreordering,  /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */

   int                   initbatchsize,      /**< the initial batchsize for the first iteration */
   int                   maxbatchsize,       /**< the maximum batchsize per iteration */
   SCIP_Real             batchingfactor,     /**< the factor with which the batchsize is multiplied in every update */
   SCIP_Real             batchingoffset,     /**< the offset which is added to the multiplied batchsize in every update */
   int                   batchupdateinterval /**< the number of iterations to run with a constant batchsize before updating (1: always update) */
   )
{
   SCIP* scip;
   SCIP_CONS** origconss;
   SCIP_CONS** conss;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_SOL* sol;
   SCIP_SOL* copysol;
   SCIP_Bool* inIS;
   int* order;
   SCIP_Bool feasible;
   SCIP_Bool stopiter;
   SCIP_RETCODE retcode;
   SCIP_RESULT result;
   int nconss;
   int batchsize;
   int iteration;
   int i;
   int j;
   int k;

   /* get current subscip */
   scip = SCIPiisGetSubscip(iis);
   assert( scip != NULL );
   assert( SCIPiisIsSubscipInfeasible(iis) );

   /* get random generator */
   randnumgen = SCIPiisGetRandnumgen(iis);
   assert( randnumgen != NULL );

   /* get batch size */
   assert( initbatchsize >= 1 );
   assert( maxbatchsize >= 1 );
   initbatchsize = MIN(initbatchsize, maxbatchsize);
   batchsize = initbatchsize;

   /* get constraint information */
   nconss = SCIPgetNOrigConss(scip);
   origconss = SCIPgetOrigConss(scip);
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &conss, origconss, nconss) );

   /* Initialise information for whether a constraint is in the final infeasible system */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &inIS, nconss) );

   /* First capture and then delete all constraints */
   assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
   for( i = 0; i < nconss; ++i )
   {
      assert( SCIPconsIsInProb(conss[i]) );
      if( SCIPconsGetNUses(conss[i]) > 1 )
      {
         inIS[i] = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPcaptureCons(scip, conss[i]) );
         SCIP_CALL( SCIPdelCons(scip, conss[i]) );
         inIS[i] = FALSE;
      }
   }
   SCIPiisSetSubscipInfeasible(iis, FALSE);

   /* Prepare random order in which the constraints will be added back */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nconss) );
   for (i = 0; i < nconss; ++i)
      order[i] = i;
   SCIPrandomPermuteIntArray(randnumgen, order, 0, nconss);

   /* Continue to add constraints until an infeasible status is reached */
   i = 0;
   iteration = 0;
   feasible = TRUE;
   stopiter = FALSE;
   while( feasible )
   {
      /* Add the next batch of constraints */
      k = 0;
      while( i < nconss && k < batchsize )
      {
         if( !inIS[order[i]] )
         {
            SCIP_CALL( SCIPaddCons(scip, conss[order[i]]) );
            SCIP_CALL( SCIPreleaseCons(scip, &conss[order[i]]) );
            inIS[order[i]] = TRUE;
            ++k;
         }
         i++;
      }
      if( k == 0 )
         break;

      /* Solve the reduced problem */
      retcode = additionSubproblem(iis, timelim, timelimperiter, nodelim, nodelimperiter, &feasible, &stopiter);
      if( !silent )
         SCIPiisfinderInfoMessage(iis, FALSE);
      if( stopiter || timelim - SCIPiisGetTime(iis) <= 0 || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) )
      {
         SCIP_CALL( SCIPfreeTransform(scip) );
         break;
      }

      if( retcode == SCIP_OKAY )
      {
         /* free transform and copy solution if there is one */
         copysol = NULL;
         sol = SCIPgetBestSol(scip);
         if( sol != NULL )
         {
            SCIP_CALL( SCIPcreateSolCopyOrig(scip, &copysol, sol) );
            SCIP_CALL( SCIPunlinkSol(scip, copysol) );
         }
         SCIP_CALL( SCIPfreeTransform(scip) );
         assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

         /* update batchsize if problem feasible */
         if( feasible )
         {
            SCIP_CALL( updateBatchsize(scip, initbatchsize, maxbatchsize, iteration, FALSE, batchingfactor, batchingoffset, batchupdateinterval, &batchsize) );
         }

         /* Add any other constraints that are also feasible for the current solution */
         if( feasible && (copysol != NULL) && dynamicreordering )
         {
            k = 0;
            for( j = i; j < nconss; ++j )
            {
               /* Don't dynamically add indicator constraints */
               if( !inIS[order[j]] && ( strcmp("indicator", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[order[j]]))) != 0 ) )
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
            {
               SCIPdebugMsg(scip, "Added %d constraints by reordering dynamically.\n", k);
               if( ! silent )
                  SCIPiisfinderInfoMessage(iis, FALSE);
            }
            SCIP_CALL( SCIPfreeSol(scip, &copysol) );
         }
      }
      else
      {
         SCIP_CALL( SCIPfreeTransform(scip) );
         assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
      }
      ++iteration;
   }

   /* Release any cons not in the IS */
   for( i = 0; i < nconss; ++i )
   {
      if( !inIS[order[i]] )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &conss[order[i]]) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &order, nconss);
   SCIPfreeBlockMemoryArray(scip, &inIS, nconss);
   SCIPfreeBlockMemoryArray(scip, &conss, nconss);
   if( feasible )
      SCIPiisSetSubscipInfeasible(iis, FALSE);
   else
      SCIPiisSetSubscipInfeasible(iis, TRUE);

   return SCIP_OKAY;
}

/** perform a greedy addition or deletion algorithm to obtain an infeasible subsystem (IS).
 *
 *  This is the generation method for the greedy IIS finder rule.
 *  Depending on the parameter choices, constraints are either greedily added from an empty problem,
 *  or deleted from a complete problem. In the case of constraints being added, this is done until the problem
 *  becomes infeasible, after which one can then begin deleting constraints. In the case of deleting constraints,
 *  this is done until no more constraints (or batches of constraints) can be deleted without making
 *  the problem feasible.
 */
static
SCIP_RETCODE execIISfinderGreedy(
   SCIP_IIS*             iis,                /**< IIS data structure */
   SCIP_IISFINDERDATA*   iisfinderdata,      /**< IIS finder data */
   SCIP_RESULT*          result              /**< pointer to store the result of the IIS finder run. SCIP_DIDNOTFIND if the algorithm failed, otherwise SCIP_SUCCESS. */
   )
{
   SCIP* scip = SCIPiisGetSubscip(iis);
   SCIP_Real timelim;
   SCIP_Longint nodelim;
   SCIP_Bool removebounds;
   SCIP_Bool silent;
   SCIP_Bool alldeletionssolved = TRUE;
   int nvars;
   int nconss;
   int maxbatchsize;
   int initbatchsize;

   assert( scip != NULL );
   assert( iisfinderdata != NULL );
   assert( result != NULL );

   SCIP_CALL( SCIPgetRealParam(scip, "iis/time", &timelim) );
   SCIP_CALL( SCIPgetLongintParam(scip, "iis/nodes", &nodelim) );
   SCIP_CALL( SCIPgetBoolParam(scip, "iis/removebounds", &removebounds) );
   SCIP_CALL( SCIPgetBoolParam(scip, "iis/silent", &silent) );

   nvars = SCIPgetNOrigVars(scip);
   nconss = SCIPgetNOrigConss(scip);
   maxbatchsize = MAX(nvars, nconss);
   initbatchsize = iisfinderdata->initrelbatchsize > 0.0
         ? (int)ceil(iisfinderdata->initrelbatchsize * maxbatchsize) : MIN(iisfinderdata->initbatchsize, maxbatchsize);
   maxbatchsize = (int)ceil(iisfinderdata->maxrelbatchsize * maxbatchsize);
   maxbatchsize = MIN(iisfinderdata->maxbatchsize, maxbatchsize);
   initbatchsize = MAX(initbatchsize, 1);
   maxbatchsize = MAX(maxbatchsize, 1);

   *result = SCIP_SUCCESS;

   if( iisfinderdata->additive )
   {
      if( !silent )
      {
         SCIPdebugMsg(scip, "----- STARTING GREEDY ADDITION ALGORITHM -----\n");
      }
      SCIP_CALL( additionFilterBatch(iis, timelim, nodelim, silent, iisfinderdata->timelimperiter,
            iisfinderdata->nodelimperiter, iisfinderdata->dynamicreordering, initbatchsize, maxbatchsize,
            iisfinderdata->batchingfactor, iisfinderdata->batchingoffset, iisfinderdata->batchupdateinterval) );
      SCIPiisSetSubscipIrreducible(iis, FALSE);
      if( timelim - SCIPiisGetTime(iis) <= 0 || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) )
         return SCIP_OKAY;
   }
   else
   {
      if( !silent )
      {
         SCIPdebugMsg(scip, "----- STARTING GREEDY DELETION ALGORITHM -----\n");
      }
      SCIP_CALL( deletionFilterBatch(iis, timelim, nodelim, removebounds, silent, iisfinderdata->timelimperiter,
            iisfinderdata->nodelimperiter, iisfinderdata->conservative, initbatchsize, maxbatchsize,
            iisfinderdata->batchingfactor, iisfinderdata->batchingoffset, iisfinderdata->batchupdateinterval,
            &alldeletionssolved) );
      if( timelim - SCIPiisGetTime(iis) <= 0 || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) )
         return SCIP_OKAY;
      if( alldeletionssolved && initbatchsize == 1 )
         SCIPiisSetSubscipIrreducible(iis, TRUE);
   }

   if( iisfinderdata->delafteradd && iisfinderdata->additive )
   {
      if( !silent )
      {
         SCIPdebugMsg(scip, "----- STARTING GREEDY DELETION ALGORITHM FOLLOWING COMPLETED ADDITION ALGORITHM -----\n");
      }
      SCIP_CALL( deletionFilterBatch(iis, timelim, nodelim, removebounds, silent, iisfinderdata->timelimperiter,
            iisfinderdata->nodelimperiter, iisfinderdata->conservative, initbatchsize, maxbatchsize,
            iisfinderdata->batchingfactor, iisfinderdata->batchingoffset, iisfinderdata->batchupdateinterval,
            &alldeletionssolved) );
      if( timelim - SCIPiisGetTime(iis) <= 0 || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) )
         return SCIP_OKAY;
      if( alldeletionssolved && initbatchsize == 1 )
         SCIPiisSetSubscipIrreducible(iis, TRUE);
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of IIS finder
 */


/** copy method for IIS finder plugin (called when SCIP copies plugins) */
static
SCIP_DECL_IISFINDERCOPY(iisfinderCopyGreedy)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(iisfinder != NULL);
   assert(strcmp(SCIPiisfinderGetName(iisfinder), IISFINDER_NAME) == 0);

   /* call inclusion method of IIS finder */
   SCIP_CALL( SCIPincludeIISfinderGreedy(scip) );

   return SCIP_OKAY;
}

/** destructor of IIS finder to free user data (called when SCIP is exiting) */
/**! [SnippetIISfinderFreeGreedy] */
static
SCIP_DECL_IISFINDERFREE(iisfinderFreeGreedy)
{  /*lint --e{715}*/
   SCIP_IISFINDERDATA* iisfinderdata;

   iisfinderdata = SCIPiisfinderGetData(iisfinder);

   SCIPfreeBlockMemory(scip, &iisfinderdata);

   SCIPiisfinderSetData(iisfinder, NULL);

   return SCIP_OKAY;
}
/**! [SnippetIISfinderFreeGreedy] */

/** IIS finder generation method of IIS */
static
SCIP_DECL_IISFINDEREXEC(iisfinderExecGreedy)
{  /*lint --e{715}*/
   SCIP_IISFINDERDATA* iisfinderdata;

   assert(iisfinder != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   iisfinderdata = SCIPiisfinderGetData(iisfinder);
   assert(iisfinderdata != NULL);

   SCIP_CALL( execIISfinderGreedy(iis, iisfinderdata, result) );

   return SCIP_OKAY;
}


/*
 * IIS finder specific interface methods
 */

/** creates the greedy IIS finder and includes it in SCIP */
SCIP_RETCODE SCIPincludeIISfinderGreedy(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_IISFINDERDATA* iisfinderdata;
   SCIP_IISFINDER* iisfinder;

   /* create greedy IIS finder data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &iisfinderdata) );
   BMSclearMemory(iisfinderdata);

   SCIP_CALL( SCIPincludeIISfinderBasic(scip, &iisfinder, IISFINDER_NAME, IISFINDER_DESC, IISFINDER_PRIORITY,
         iisfinderExecGreedy, iisfinderdata) );

   assert(iisfinder != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetIISfinderCopy(scip, iisfinder, iisfinderCopyGreedy) );
   SCIP_CALL( SCIPsetIISfinderFree(scip, iisfinder, iisfinderFreeGreedy) );

   /* add greedy IIS finder parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/timelimperiter",
         "time limit of optimization process for each individual subproblem",
         &iisfinderdata->timelimperiter, FALSE, DEFAULT_TIMELIMPERITER, 0.0, SCIP_INVALID/10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "iis/" IISFINDER_NAME "/nodelimperiter",
         "node limit of optimization process for each individual subproblem",
         &iisfinderdata->nodelimperiter, FALSE, DEFAULT_NODELIMPERITER, -1L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IISFINDER_NAME "/additive",
         "should an additive constraint approach be used instead of deletion",
         &iisfinderdata->additive, FALSE, DEFAULT_ADDITIVE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IISFINDER_NAME "/conservative",
         "should a hit limit (e.g. node  time) solve be counted as feasible when deleting constraints",
         &iisfinderdata->conservative, TRUE, DEFAULT_CONSERVATIVE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IISFINDER_NAME "/delafteradd",
         "should the deletion routine be performed after the addition routine (in the case of additive)",
         &iisfinderdata->delafteradd, TRUE, DEFAULT_DELAFTERADD, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IISFINDER_NAME "/dynamicreordering",
         "should satisfied constraints outside the batch of an intermediate solve be added during the additive method",
         &iisfinderdata->dynamicreordering, TRUE, DEFAULT_DYNAMICREORDERING, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "iis/" IISFINDER_NAME "/initbatchsize",
         "the initial batchsize for the first iteration, ignored if initrelbatchsize is positive",
         &iisfinderdata->initbatchsize, FALSE, DEFAULT_INITBATCHSIZE, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/initrelbatchsize",
         "the initial batchsize relative to the original problem for the first iteration (0.0: use initbatchsize)",
         &iisfinderdata->initrelbatchsize, FALSE, DEFAULT_INITRELBATCHSIZE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "iis/" IISFINDER_NAME "/maxbatchsize",
         "the maximum batchsize per iteration",
         &iisfinderdata->maxbatchsize, TRUE, DEFAULT_MAXBATCHSIZE, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/maxrelbatchsize",
         "the maximum batchsize relative to the original problem per iteration",
         &iisfinderdata->maxrelbatchsize, TRUE, DEFAULT_MAXRELBATCHSIZE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/batchingfactor",
         "the factor with which the batchsize is multiplied in every update",
         &iisfinderdata->batchingfactor, TRUE, DEFAULT_BATCHINGFACTOR, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/batchingoffset",
         "the offset which is added to the multiplied batchsize in every update",
         &iisfinderdata->batchingoffset, TRUE, DEFAULT_BATCHINGOFFSET, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "iis/" IISFINDER_NAME "/batchupdateinterval",
         "the number of iterations to run with a constant batchsize before updating (1: always update)",
         &iisfinderdata->batchupdateinterval, TRUE, DEFAULT_BATCHUPDATEINTERVAL, 1, INT_MAX, NULL, NULL) );


   return SCIP_OKAY;
}

/** perform the greedy deletion algorithm with singleton batches to obtain an irreducible infeasible subsystem (IIS) */
SCIP_RETCODE SCIPiisGreedyMinimize(
   SCIP_IIS*             iis                 /**< IIS data structure */
   )
{
   SCIP* scip = SCIPiisGetSubscip(iis);
   SCIP_Real timelim;
   SCIP_Longint nodelim;
   SCIP_Bool removebounds;
   SCIP_Bool silent;
   SCIP_Bool alldeletionssolved = TRUE;
   int nvars;
   int nconss;
   int maxbatchsize;

   assert( scip != NULL );

   if( !SCIPiisIsSubscipInfeasible(iis) )
   {
      SCIPerrorMessage("infeasible problem required\n");
      return SCIP_INVALIDDATA;
   }

   nvars = SCIPgetNOrigVars(scip);
   nconss = SCIPgetNOrigConss(scip);
   maxbatchsize = MAX(nvars, nconss);

   SCIP_CALL( SCIPgetRealParam(scip, "iis/time", &timelim) );
   SCIP_CALL( SCIPgetLongintParam(scip, "iis/nodes", &nodelim) );
   SCIP_CALL( SCIPgetBoolParam(scip, "iis/removebounds", &removebounds) );
   SCIP_CALL( SCIPgetBoolParam(scip, "iis/silent", &silent) );

   SCIP_CALL( deletionFilterBatch(iis, timelim, nodelim, removebounds, silent,
         DEFAULT_TIMELIMPERITER, DEFAULT_NODELIMPERITER, TRUE, 1, maxbatchsize,
         DEFAULT_BATCHINGFACTOR, DEFAULT_BATCHINGOFFSET, DEFAULT_BATCHUPDATEINTERVAL, &alldeletionssolved) );
   if( alldeletionssolved && SCIPiisGetTime(iis) < timelim && ( nodelim == -1 || SCIPiisGetNNodes(iis) < nodelim ) )
      SCIPiisSetSubscipIrreducible(iis, TRUE);

   return SCIP_OKAY;
}
