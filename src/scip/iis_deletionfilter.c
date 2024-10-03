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


/** solve subproblem for deletionFilter */
static
SCIP_RETCODE deletionFilterConsSubproblem(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             conservative,       /**< whether we treat a subproblem to be feasible, if it reached its node limt */
   int                   currentidx,         /**< current index of constraint to be deleted */
   int                   cnt,                /**< number of subproblems solved (for output) */
   SCIP_Bool             silent,             /**< run silently? */
   SCIP_Bool             finalrun,           /**< is this the final run? */
   SCIP_Bool*            IS,                 /**< output: array indicating which constraints take part in (I)IS */
   int*                  sizeIS,             /**< pointer to store the size of the (I)IS */
   SCIP_Bool*            isIIS,              /**< pointer to store whether we found an IIS */
   SCIP_Bool*            success,            /**< pointer to store whether we have obtained an (I)IS */
   SCIP_Bool*            feasible,           /**< pointer to store whether the subproblem was feasible */
   SCIP_Bool*            stop                /**< pointer to store whether we have to stop */
   )
{
   char consname[SCIP_MAXSTRLEN];
   SCIP_STATUS status;
   SCIP_CONS** conss;
   int nconss;
   int j;

   assert( feasible != NULL );
   assert( stop != NULL );
   assert( success != NULL );

   *success = TRUE;
   *feasible = FALSE;
   *stop = FALSE;

   /* transform problem and locally remove constraint */
   SCIP_CALL( SCIPtransformProb(scip) );

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   assert( nconss == SCIPgetNOrigConss(scip) );
   assert( 0 <= currentidx && currentidx < nconss );

   /* store constraint name */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "(%d: %s)", currentidx, SCIPconsGetName(conss[currentidx]));

   /* remove constraints */
   for (j = nconss - 1; j >= 0; --j)
   {
      if ( j == currentidx || ! IS[j] )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[j]) );
      }
   }

#ifdef SCIP_DEBUG
   if ( finalrun )
   {
      char name[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "delfilterres_%s.cip", SCIPgetProbName(scip));
      SCIP_CALL( SCIPwriteTransProblem(scip, name, "cip", FALSE) );
   }
   else
   {
      SCIP_CALL( SCIPwriteOrigProblem(scip, "deletionfilter.cip", "cip", FALSE) );
   }
#endif

   /* solve problem until first solution is found or infeasibility has been proved */
   SCIP_CALL( SCIPsolve(scip) );
   status = SCIPgetStatus(scip);

   /* free transform */
   SCIP_CALL( SCIPfreeTransform(scip) );

   /* check status */
   switch ( status )
   {
   case SCIP_STATUS_TIMELIMIT:        /* if we reached the time limit, then stop */
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "%d: Time limit exceeded (removed constraint <%s>).\n", cnt, consname);

      *success = FALSE;
      *isIIS = FALSE;  /* cannot be sure whether we have an IIS */
      *stop = TRUE;
      break;

   case SCIP_STATUS_USERINTERRUPT:    /* if an user interrupt occurred, just stop */
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "%d: User interrupt (removed constraint <%s>).\n", cnt, consname);

      *success = FALSE;
      *isIIS = FALSE;  /* cannot be sure whether we have an IIS */
      *stop = TRUE;
      break;

   case SCIP_STATUS_NODELIMIT:        /* if we reached the node limit */
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "%d: Node limit reached (removed constraint <%s>).\n", cnt, consname);

      if ( ! conservative && currentidx < nconss )
      {
         assert( IS[currentidx] );
         IS[currentidx] = FALSE;   /* no solution found -> remove constraint */
         --(*sizeIS);
      }
      *isIIS = FALSE;  /* cannot be sure whether we have an IIS */
      *success = FALSE;
      break;

   case SCIP_STATUS_INFEASIBLE:       /* if the problem is infeasible */
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "%d: Subproblem infeasible (removed constraint <%s>).\n", cnt, consname);

      if ( currentidx < nconss )
      {
         IS[currentidx] = FALSE;
         --(*sizeIS);
      }
      *success = TRUE;
      break;

   case SCIP_STATUS_BESTSOLLIMIT:     /* we found a solution */
   case SCIP_STATUS_OPTIMAL:
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "%d: Found solution (removed constraint <%s>).\n", cnt, consname);
      /* do nothing */
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
SCIP_RETCODE deletionFilterCons(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_Bool             conservative,       /**< whether we treat a subproblem to be feasible, if it reached its node limt */
   SCIP_Real             timelimit,          /**< total time limit */
   SCIP_Longint          nodelimit,          /**< maximal node limit for each run */
   SCIP_Bool             silent,             /**< run silently? */
   SCIP_Bool*            IS,                 /**< output: array indicating which constraints take part in (I)IS */
   int*                  sizeIS,             /**< pointer to store the size of the (I)IS */
   SCIP_Bool*            isIIS,              /**< pointer to store whether we found an IIS */
   SCIP_Longint*         nnodes,             /**< pointer to store the total number of nodes needed (or NULL) */
   SCIP_Bool*            success             /**< pointer to store whether we have obtained an (I)IS */
   )
{
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_CLOCK* totalTimeClock = NULL;
   SCIP_STATUS status;
   SCIP_Longint initnnodes = -1LL;
   SCIP_Longint oldnodelimit;
   SCIP_Real oldtimelimit;
   int oldverblevel;
   int olddispfreq;
   int oldbestsollimit;
   int nconss;
   int* order;
   int i;

   assert( IS != NULL );
   assert( sizeIS != NULL );
   assert( isIIS != NULL );
   assert( success != NULL );

   *success = FALSE;
   *isIIS = FALSE;
   if ( nnodes != NULL )
      *nnodes = 0;

   nconss = SCIPgetNOrigConss(scip);
   if ( ! silent )
   {
      SCIPinfoMessage(scip, NULL, "Number of constraints:\t\t%d\n", nconss);
      if ( nodelimit < SCIP_LONGINT_MAX )
         SCIPinfoMessage(scip, NULL, "Maximal node limit:\t\t%"SCIP_LONGINT_FORMAT"\n", nodelimit);
      if ( ! SCIPisInfinity(scip, timelimit) )
         SCIPinfoMessage(scip, NULL, "Time limit:\t\t\t%f\n", timelimit);
   }

   /* create and start clock */
   SCIP_CALL( SCIPcreateClock(scip, &totalTimeClock) );
   SCIP_CALL( SCIPstartClock(scip, totalTimeClock) );

   /* init (I)IS and problem */
   *sizeIS = nconss;
   for (i = 0; i < nconss; ++i)
      IS[i] = TRUE;

   /* get current parameters */
   SCIP_CALL( SCIPgetIntParam(scip, "display/verblevel", &oldverblevel) );
   SCIP_CALL( SCIPgetIntParam(scip, "display/freq", &olddispfreq) );
   SCIP_CALL( SCIPgetLongintParam(scip, "limits/nodes", &oldnodelimit) );
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &oldtimelimit) );
   SCIP_CALL( SCIPgetIntParam(scip, "limits/bestsol", &oldbestsollimit) );

   /* set output to console */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", 1000) );
#else
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", 1000) );
#endif

   /* set node limit */
   if ( nodelimit < SCIP_LONGINT_MAX )
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", nodelimit) );
   }

   /* set time limit */
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      SCIP_Real t;

      t = timelimit - SCIPgetClockTime(scip, totalTimeClock);
      if ( t <= 0.0 )
      {
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "Time limit exceeded in initial problem.\n");
         goto TERMINATE;
      }
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", t) );
   }

   /* stop solution process as soon as solution has been found */
   SCIP_CALL( SCIPsetIntParam(scip, "limits/bestsol", 1) );

   /* possibly ouput problem */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPwriteOrigProblem(scip, "deletionfilter.cip", "cip", FALSE) );
#endif

   /* Solve problem until infeasibility is proved (or first solution is found) and determine node number as estimate. */
   SCIP_CALL( SCIPsolve(scip) );
   status = SCIPgetStatus(scip);
   initnnodes = (SCIP_Longint) (DEFAULT_FACTORNODES * (SCIP_Real)SCIPgetNTotalNodes(scip));
   initnnodes = MAX(DEFAULT_MINNNODES, initnnodes);
   if ( nnodes != NULL )
      *nnodes += SCIPgetNTotalNodes(scip);

   /* check whether original problem is infeasible */
   if ( status != SCIP_STATUS_INFEASIBLE )
   {
      if ( ! silent )
      {
         if ( status == SCIP_STATUS_BESTSOLLIMIT || status == SCIP_STATUS_OPTIMAL )
         {
            SCIPinfoMessage(scip, NULL, "Original problem is feasible.\n");
            *sizeIS = 0;
            *isIIS = TRUE;
            *success = TRUE;
            for (i = 0; i < nconss; ++i)
               IS[i] = FALSE;
         }
         else
            SCIPinfoMessage(scip, NULL, "Status %d initial problem -> exit.\n", status);
      }
      goto TERMINATE;
   }

   /* set new node limit */
   if ( initnnodes < nodelimit )
   {
      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Node limit for subproblems:\t%"SCIP_LONGINT_FORMAT"\n\n", initnnodes);
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", initnnodes) );
   }

   /* ---------------------------------------------------------------------------------------------------*/
   /* reset problem */
   SCIP_CALL( SCIPfreeTransform(scip) );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

   /* Init status to be an IIS: If in one iteration, we obtained an indecisive solution status, we cannot be sure that
    * we found an IIS, so this will be reset. */
   *isIIS = TRUE;

   /* prepare random order */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nconss) );
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 4678, FALSE) );
   for (i = 0; i < nconss; ++i)
      order[i] = i;
   SCIPrandomPermuteIntArray(randnumgen, order, 0, nconss);
   SCIPfreeRandom(scip, &randnumgen);

   /* loop through all constraints in random order */
   for (i = 0; i < nconss; ++i)
   {
      SCIP_Bool stopiter = FALSE;
      SCIP_Bool feasible;
      int idx;

      idx = order[i];

      /* treat subproblem */
      SCIP_CALL( deletionFilterConsSubproblem(scip, conservative, idx, i,
            silent, FALSE, IS, sizeIS, isIIS, success, &feasible, &stopiter) );

      if ( nnodes != NULL )
         *nnodes += SCIPgetNTotalNodes(scip);
      assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

      if ( stopiter )
         break;
   }
   SCIPfreeBlockMemoryArray(scip, &order, nconss);

   /* one final check whether we have an IS if the loop above completed, we actually reduced the size (otherwise the status has been determined already) */
   if ( i >= nconss && *sizeIS < nconss && ! (*success) && ! (*isIIS) )
   {
      SCIP_Bool stopiter = FALSE;
      SCIP_Bool feasible = FALSE;

      if ( ! silent )
         SCIPinfoMessage(scip, NULL, "Solving final problem to determine whether we have an IS.\n");

      /* possibly reset node limit */
      if ( initnnodes < nodelimit )
      {
         SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", nodelimit) );
      }

      /* treat subproblem */
      *success = TRUE;
      SCIP_CALL( deletionFilterConsSubproblem(scip, conservative, nconss, nconss,
            silent, TRUE, IS, sizeIS, isIIS, success, &feasible, &stopiter) );
      if ( feasible )
         *success = FALSE;
      if ( nnodes != NULL )
         *nnodes += SCIPgetNTotalNodes(scip);
   }

   if ( ! silent )
   {
      SCIPinfoMessage(scip, NULL, "\n");
      if ( *success )
      {
         if ( *isIIS )
            SCIPinfoMessage(scip, NULL, "Size of found IIS: %d.\n", *sizeIS);
         else
            SCIPinfoMessage(scip, NULL, "Size of found  IS: %d.\n", *sizeIS);
      }
      else
         SCIPinfoMessage(scip, NULL, "Run unsuccessful.\n");

      SCIPinfoMessage(scip, NULL, "\nTotal used time:\t %f\n", SCIPgetClockTime(scip, totalTimeClock));
   }

 TERMINATE:
   /* reset parameters */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", oldverblevel) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", olddispfreq) );
   SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", oldnodelimit) );
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", oldtimelimit) );
   SCIP_CALL( SCIPsetIntParam(scip, "limits/bestsol", oldbestsollimit) );

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
   SCIP_Real timelimit;
   SCIP_Longint nodelimit;
   SCIP_Bool* IS;
   int nconss;

   /* determine node limit */
   nodelimit = SCIP_LONGINT_MAX;

   /* set time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      SCIP_Real t;

      t = timelimit - SCIPgetTotalTime(scip);
      if ( t <= 0.0 )
      {
         if ( ! silent )
            SCIPinfoMessage(scip, NULL, "Time limit already exceeded.\n");
         return SCIP_OKAY;
      }
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", t) );
   }

   nconss = SCIPgetNOrigConss(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &IS, nconss) );
   SCIP_CALL( deletionFilterCons(scip, TRUE, timelimit, nodelimit, FALSE, IS, sizeIS, isIIS, NULL, success) );
   SCIPfreeBlockMemoryArray(scip, &IS, nconss);

   return SCIP_OKAY;
}
