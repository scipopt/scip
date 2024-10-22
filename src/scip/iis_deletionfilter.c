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

/**@file   iis_deletionfilter.c
 * @brief  deletion filter heuristic to compute (I)ISs
 * @author Marc Pfetsch
 * @author Mark Turner
 *
 * An irreducibly infeasible subsystem (IIS) is a subset of the constraints that is infeasible and (set-wise) minimial
 * in this respect. The deletion filter heuristic greedily removes constraints and checks whether the remaining problem
 * is still infeasible. The method is based on
 *
 * O. Guieu and J. Chinneck, Analyzing infeasible mixed-integer and integer linear programs,@p
 * INFORMS J. Comput. 11, no. 1 (1999), pp. 63–77.
 *
 * We cannot guarantee that we are minimal at the end, so we might just obtain an infeasible subsystem (IS).
 */

#include <string.h>

#include <scip/scipdefplugins.h>
#include <scip/iis_deletionfilter.h>

/* default values */
#define DEFAULT_MINNNODES 50
#define DEFAULT_FACTORNODES 2.0
#define DEFAULT_BATCHSIZE 1

/** creates a sub-SCIP and sets parameters */
static
SCIP_RETCODE createSubscipCopy(
   SCIP*                 scip,               /**< main SCIP data structure */
   SCIP**                subscip             /**< pointer to store created sub-SCIP */
   )
{
   SCIP_Bool success;
   
   assert(scip != NULL);
   assert(subscip != NULL);
   
   /* create a new SCIP instance */
   SCIP_CALL( SCIPcreate(subscip) );
   
   /* create problem in sub-SCIP */
   SCIP_CALL( SCIPcopyOrig(scip, *subscip, NULL, NULL, "iis", TRUE, FALSE, FALSE, &success) );
   
   if( success == FALSE )
   {
      return SCIP_ERROR;
   }
   /* copy parameter settings */
   SCIP_CALL( SCIPcopyParamSettings(scip, *subscip) );
   
   /* avoid recursive calls */
   SCIP_CALL( SCIPsetSubscipsOff(*subscip, TRUE) );
   
   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(*subscip, SCIP_PARAMSETTING_FAST, TRUE) );
   
   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(*subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(*subscip, "timing/statistictiming", FALSE) );
#endif
   
   /* set parameter for solve to stop after finding a single solution (only need to prove feasibility) */
   SCIP_CALL( SCIPsetIntParam(*subscip, "limits/bestsol", 1) );
   
   return SCIP_OKAY;
}


/** solve subproblem for deletionFilter */
static
SCIP_RETCODE deletionFilterConsSubproblem(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             conservative,       /**< whether we treat a subproblem to be feasible, if it reached its node limt */
   int                   ndelconss,          /**< the number of constraints that will be deleted */
   int*                  idxs,               /**< the indices of the constraints (in the conss array) that will be deleted */
   SCIP_CONS**           conss,              /**< the array of constraints (may be a superset of the current constraints) */
   SCIP_Bool             silent,             /**< run silently? */
   SCIP_Bool*            stop                /**< pointer to store whether we have to stop */
   )
{
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
   SCIP_CALL( SCIPsolve(scip) );
   status = SCIPgetStatus(scip);

   /* free transform */
   SCIP_CALL( SCIPfreeTransform(scip) );

   /* check status */
   switch ( status )
   {
   case SCIP_STATUS_TIMELIMIT:        /* if we reached the time limit */
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Time limit exceeded. Not removing batch.\n");
      for( j = 0; j < ndelconss; j++ )
      {
         SCIP_CALL( SCIPaddCons(scip, conss[idxs[j]]) );
         SCIP_CALL( SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      break;

   case SCIP_STATUS_USERINTERRUPT:    /* if an user interrupt occurred, just stop */
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "User interrupt. Stopping. \n");
      for( j = 0; j < ndelconss; j++ )
      {
         SCIP_CALL( SCIPaddCons(scip, conss[idxs[j]]) );
         SCIP_CALL( SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      *stop = TRUE;
      break;

   case SCIP_STATUS_NODELIMIT:        /* if we reached the node limit */
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Node limit reached. Removing batch if set as non-conservative. \n");
      
      for( j = 0; j < ndelconss; j++ )
      {
         if( conservative )
            SCIP_CALL( SCIPaddCons(scip, conss[idxs[j]]) );
         SCIP_CALL(SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      break;

   case SCIP_STATUS_INFEASIBLE:       /* if the problem is infeasible */
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Subproblem infeasible. Remove constraint batch.\n");
      
      for( j = 0; j < ndelconss; j++ )
      {
         SCIP_CALL(SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      break;

   case SCIP_STATUS_BESTSOLLIMIT:     /* we found a solution */
   case SCIP_STATUS_OPTIMAL:
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Found solution. Keeping constraint batch\n");
      for( j = 0; j < ndelconss; j++ )
      {
         SCIP_CALL( SCIPaddCons(scip, conss[idxs[j]]) );
         SCIP_CALL( SCIPreleaseCons(scip, &conss[idxs[j]]) );
      }
      break;

   case SCIP_STATUS_UNKNOWN:
   case SCIP_STATUS_TOTALNODELIMIT:
   case SCIP_STATUS_STALLNODELIMIT:
   case SCIP_STATUS_MEMLIMIT:
   case SCIP_STATUS_GAPLIMIT:
   case SCIP_STATUS_SOLLIMIT:
   case SCIP_STATUS_RESTARTLIMIT:
   case SCIP_STATUS_UNBOUNDED:
   case SCIP_STATUS_INFORUNBD:
   default:
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
   SCIP_Bool*            feasible,           /**< pointer to store whether the problem is now feasible */
   SCIP_Bool*            stop                /**< pointer to store whether we have to stop */
)
{
   SCIP_STATUS status;
   
   assert( stop != NULL );
   
   *stop = FALSE;
   *feasible = FALSE;
   
   /* solve problem until first solution is found or infeasibility has been proven */
   SCIP_CALL( SCIPsolve(scip) );
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
      
      case SCIP_STATUS_NODELIMIT:        /* if we reached the node limit */
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "Node limit reached. Added batch failed to induce infeasibility. \n");
         break;
      
      case SCIP_STATUS_INFEASIBLE:       /* if the problem is infeasible */
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "Subproblem infeasible. Final batch of constraints added.\n");
         break;
      
      case SCIP_STATUS_BESTSOLLIMIT:     /* we found a solution */
      case SCIP_STATUS_OPTIMAL:
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "Found solution. Keep adding constraint batches\n");
         *feasible = TRUE;
         break;
      
      case SCIP_STATUS_UNKNOWN:
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
      case SCIP_STATUS_UNBOUNDED:
      case SCIP_STATUS_INFORUNBD:
      default:
         SCIPerrorMessage("unexpected return status %d. Exiting ...\n", status);
         return SCIP_ERROR;
   }
   
   return SCIP_OKAY;
}


/** deletion filter to greedily remove constraints to obtain an (I)IS -- detailed function call */
SCIP_RETCODE deletionFilterBatchCons(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             conservative,       /**< whether we treat a subproblem to be feasible, if it reached its node limt */
   SCIP_Bool             silent,             /**< run silently? */
   SCIP_Longint*         nnodes,             /**< pointer to store the total number of nodes needed (or NULL) */
   SCIP_Bool*            success             /**< pointer to store whether we have obtained an (I)IS */
   )
{
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_CLOCK* totalTimeClock = NULL;
   int nconss;
   int nbatches;
   SCIP_CONS** conss;
   int* order;
   int* idxs;
   int i;
   
   assert( success != NULL );

   *success = FALSE;
   if ( nnodes != NULL )
      *nnodes = 0;

   /* create and start clock */
   SCIP_CALL( SCIPcreateClock(scip, &totalTimeClock) );
   SCIP_CALL( SCIPstartClock(scip, totalTimeClock) );
   
   /* Get constraint information */
   nconss = SCIPgetNOrigConss(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conss, nconss) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &conss, SCIPgetOrigConss(scip), nconss) );

   /* reset problem */
   SCIP_CALL( SCIPfreeTransform(scip) );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

   /* prepare random order */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nconss) );
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 4678, FALSE) );
   for (i = 0; i < nconss; ++i)
      order[i] = i;
   SCIPrandomPermuteIntArray(randnumgen, order, 0, nconss);
   SCIPfreeRandom(scip, &randnumgen);
   
   /* Calculate the number of batches */
   nbatches = (int) SCIPceil(scip, (SCIP_Real) nconss / DEFAULT_BATCHSIZE);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idxs, DEFAULT_BATCHSIZE) );

   /* Loop through all batches of constraints in random order */
   for (i = 0; i < nbatches; ++i)
   {
      int minconsidx;
      int maxconsidx;
      int j;
      int ndelconss;
      SCIP_Bool stopiter = FALSE;
   
      minconsidx = i * DEFAULT_BATCHSIZE;
      if( minconsidx >= nconss )
         break;
      maxconsidx = MIN((i + 1) * DEFAULT_BATCHSIZE, nconss);
      ndelconss = maxconsidx - minconsidx;
      for( j = 0; j < ndelconss; j++ )
      {
         idxs[j] = order[i*DEFAULT_BATCHSIZE+j];
      }

      /* treat subproblem */
      SCIP_CALL( deletionFilterConsSubproblem(scip, conservative, ndelconss, idxs, conss,
            silent, &stopiter) );

      if ( nnodes != NULL )
         *nnodes += SCIPgetNTotalNodes(scip);
      assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

      if ( stopiter )
         break;
   }
   
   SCIPfreeBlockMemoryArray(scip, &idxs, DEFAULT_BATCHSIZE);
   SCIPfreeBlockMemoryArray(scip, &order, nconss);
   SCIPfreeBlockMemoryArray(scip, &conss, nconss);
   
   SCIP_CALL( SCIPfreeClock(scip, &totalTimeClock) );
   assert( totalTimeClock == NULL );

   return SCIP_OKAY;
}

/** addition filter to greedily add constraints to obtain an (I)IS -- detailed function call */
SCIP_RETCODE additionFilterBatchCons(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             silent,             /**< run silently? */
   SCIP_Longint*         nnodes,             /**< pointer to store the total number of nodes needed (or NULL) */
   SCIP_Bool*            success             /**< pointer to store whether we have obtained an (I)IS */
)
{
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_CLOCK* totalTimeClock = NULL;
   int nconss;
   SCIP_CONS** conss;
   SCIP_SOL* sol;
   SCIP_SOL* copysol;
   SCIP_Bool* inIS;
   int* order;
   SCIP_Bool feasible;
   SCIP_Bool stopiter;
   SCIP_RESULT result;
   int i;
   int j;
   int k;
   
   assert( success != NULL );
   
   *success = FALSE;
   if ( nnodes != NULL )
      *nnodes = 0;
   
   /* create and start clock */
   SCIP_CALL( SCIPcreateClock(scip, &totalTimeClock) );
   SCIP_CALL( SCIPstartClock(scip, totalTimeClock) );
   
   /* Get constraint information */
   nconss = SCIPgetNOrigConss(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conss, nconss) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &conss, SCIPgetOrigConss(scip), nconss) );
   
   /* Initialise information for whether a constraint is in the final infeasible system */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &inIS, nconss) );
   
   /* reset problem */
   SCIP_CALL( SCIPfreeTransform(scip) );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
   for( i = 0; i < nconss; ++i )
   {
      assert( SCIPconsIsInProb(conss[i]) );
      SCIP_CALL( SCIPcaptureCons(scip, conss[i]) );
      SCIP_CALL( SCIPdelCons(scip, conss[i]) );
      inIS[i] = FALSE;
   }
   
   /* prepare random order */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nconss) );
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 4678, FALSE) );
   for (i = 0; i < nconss; ++i)
      order[i] = i;
   SCIPrandomPermuteIntArray(randnumgen, order, 0, nconss);
   SCIPfreeRandom(scip, &randnumgen);
   i = 0;
   feasible = TRUE;
   stopiter = FALSE;
   result = SCIP_FEASIBLE;
   
   /* Continue to add constraints until an infeasible status is reached */
   while( feasible && !stopiter )
   {
      /* Add the next batch of constraints */
      k = 0;
      for( j = i; j < nconss; ++j )
      {
         if( k >= DEFAULT_BATCHSIZE )
            break;
         if( !inIS[order[j]] )
         {
            SCIP_CALL( SCIPaddCons(scip, conss[order[j]]) );
            SCIP_CALL( SCIPreleaseCons(scip, &conss[order[j]]) );
            inIS[order[j]] = TRUE;
            ++k;
         }
      }
      i = i + k;
      
      /* Solve the reduced problem */
      SCIP_CALL( additionFilterConsSubproblem(scip, silent, &feasible, &stopiter) );
   
      if ( nnodes != NULL )
         *nnodes += SCIPgetNTotalNodes(scip);
   
      /* free transform */
      // TODO: Figure out the correct way to handle the solution copying. Need to copy -> freeTransform -> checkCons -> deleteSol (not sure on last step)
      sol = SCIPgetBestSol(scip);
      if( sol != NULL )
      {
         SCIP_CALL( SCIPcreateSolCopyOrig(scip, &copysol, sol) );
         SCIP_CALL( SCIPunlinkSol(scip, copysol) );
      }
      SCIP_CALL( SCIPfreeTransform(scip) );
      assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
      
      /* Add any other constraints that are also feasible for the current solution */
      if( feasible && !stopiter && copysol != NULL )
      {
         for( j = i; j < nconss; ++j )
         {
            if( !inIS[order[j]] )
            {
               SCIP_CALL(SCIPcheckCons(scip, conss[order[j]], copysol, FALSE, FALSE, FALSE, &result) );
               if( result == SCIP_FEASIBLE )
               {
                  SCIP_CALL( SCIPaddCons(scip, conss[order[j]]) );
                  SCIP_CALL( SCIPreleaseCons(scip, &conss[order[j]]) );
                  inIS[order[j]] = TRUE;
               }
            }
         }
      }
   
      if ( stopiter )
         break;
      
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
   
   SCIP_CALL( SCIPfreeClock(scip, &totalTimeClock) );
   assert( totalTimeClock == NULL );
   
   return SCIP_OKAY;
}


/** run deletion filter to obtain an (I)IS */
SCIP_RETCODE SCIPrunDeletionFilter(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             silent,             /**< run silently? */
   int*                  sizeIS,             /**< pointer to store the size of the (I)IS */
   SCIP_Bool*            isIIS,              /**< pointer to store whether we found an IIS */
   SCIP_Bool*            success             /**< pointer to store whether we have obtained an (I)IS */
   )
{
   SCIP* subscip;
   
   SCIP_CALL( createSubscipCopy(scip, &subscip) );
   SCIP_CALL( deletionFilterBatchCons(subscip, TRUE, FALSE, NULL, success) );
   // SCIP_CALL( additionFilterBatchCons(subscip, FALSE, NULL, success) );
   SCIPinfoMessage(scip, NULL, "NCONSS: %d\n", SCIPgetNConss(subscip));

   return SCIP_OKAY;
}
