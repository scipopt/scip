/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_listscheduling.c
 * @brief  scheduling specific primal heuristic which is based on bidirectional serial generation scheme.
 * @author Jens Schulz
 *
 * @page HEUR List scheduling heuristic
 *
 * The heuristic performs a serial SGS (schedule generation scheme), see Kolisch and Hartmann 2006 \cite KolH06.
 * Therefore, the jobs are considered in a topological order (e.g., sorted by their earliest start) and are scheduled
 * according to that order as early as possible respecting the precedence and resource constraints.
 *
 * The serial generation scheme is extended to bidirectional SGS; see Li and Willis 1992 \cite LiW92. The first obtained
 * schedule is the so-called forward schedule.  Then, all jobs are sorted in non-increasing order of their completion
 * times in the forward schedule.  According to that ordering, a backward schedule is created by scheduling all jobs as
 * late as possible again with respect to precedence and resource constraints.  It gets clear from the way the algorithm
 * works, that if a feasible forward schedule has been found, a feasible backward schedule can be obtained, since no job
 * needs to be scheduled earlier as in the forward schedule.  Recreating a forward schedule by sorting the jobs
 * according to their start times in the backward schedule leads to a makespan not larger than the one in the first
 * forward schedule.
 *
 * @section REFERENCES References
 *
 * -# Rainer Kolisch and S&ouml;nke Hartmann. Experimental investigation of heuristics for resource-constrained
 *    project scheduling: An update. <em>European Journal of Operational Research</em>, 174(1):23&ndash;37, 2006.
 * -# K.Y. Li and R.J. Willis. An iterative scheduling technique for resource-constrained project
 *    scheduling. <em>European Journal of Operational Research</em>, 56(3):370&ndash;379, 1992.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_listscheduling.h"
#include "scip/pub_misc.h"

/**@name Properties of the heuristic
 *
 */

#define HEUR_NAME             "listscheduling"
#define HEUR_DESC             "scheduling specific primal heuristic which is based on bidirectional serial generation scheme"
#define HEUR_DISPCHAR         'x'
#define HEUR_PRIORITY         10000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         100000
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE | SCIP_HEURTIMING_BEFOREPRESOL
#define HEUR_USESSUBSCIP      FALSE      /**< does the heuristic use a secondary SCIP instance? */

/**@} */

/*
 * Data structures
 */


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_DIGRAPH*         precedencegraph;    /**< precedence graph of the jobs */
   int**                 resourcedemands;    /**< resource demands matrix (job i needs resourcedemands[i][j] units of resource j) */
   SCIP_VAR**            vars;               /**< array of start time variables */
   int*                  durations;          /**< array of duration for each job */
   int*                  capacities;         /**< array to store the capacities of all cum constraints */
   int                   njobs;              /**< number of jobs */
   int                   nresources;         /**< number of resources */
   SCIP_Bool             initialized;        /**< stores if initialization has already occurred */
};

/**@name Local methods
 *
 */

/** initializes heuristic data structures */
static
SCIP_RETCODE heurdataInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_DIGRAPH*         precedencegraph,    /**< precedence graph */
   SCIP_VAR**            vars,               /**< start time variables */
   int*                  durations,          /**< duration of the jobs independent of the resources */
   int**                 resourcedemands,    /**< resource demand matrix */
   int*                  capacities,         /**< capacities of the resources */
   int                   njobs,              /**< number if jobs */
   int                   nresources          /**< number of resources */
   )
{
   int j;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(!heurdata->initialized);

   heurdata->nresources = nresources;
   heurdata->njobs = njobs;

   if( njobs == 0 )
   {
      heurdata->resourcedemands = NULL;
      heurdata->capacities = NULL;
      heurdata->precedencegraph = NULL;
   }
   else
   {
      /* copy precedence graph */
      SCIP_CALL( SCIPdigraphCopy(&heurdata->precedencegraph, precedencegraph) );

      /* topological sort the precedence graph */
      SCIP_CALL( SCIPdigraphComputeUndirectedComponents(heurdata->precedencegraph, -1, NULL, NULL) );
      assert(SCIPdigraphGetNComponents(heurdata->precedencegraph) == 1);

      /* use the topological sorted for the variables */
      SCIP_CALL( SCIPdigraphTopoSortComponents(heurdata->precedencegraph) );
      SCIPdebug( SCIPdigraphPrintComponents(heurdata->precedencegraph, SCIPgetMessagehdlr(scip), NULL) );

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->capacities, capacities, nresources) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->vars, vars, njobs) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->durations, durations, njobs) );

      SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->resourcedemands, njobs) );
      for( j = 0; j < njobs; ++j )
      {
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->resourcedemands[j], resourcedemands[j], nresources) );
      }

      heurdata->initialized = TRUE;
   }

   return SCIP_OKAY;
}

/* frees heuristic data structures */
static
SCIP_RETCODE heurdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   int njobs;

   assert(scip != NULL);
   assert(heurdata != NULL);

   njobs = heurdata->njobs;

   if( njobs > 0 )
   {
      int j;

      for( j = 0; j < njobs; ++j )
      {
         SCIPfreeMemoryArray(scip, &heurdata->resourcedemands[j]);
      }

      SCIPfreeMemoryArray(scip, &heurdata->resourcedemands);
      SCIPfreeMemoryArray(scip, &heurdata->capacities);
      SCIPfreeMemoryArray(scip, &heurdata->vars);
      SCIPfreeMemoryArray(scip, &heurdata->durations);
      SCIPdigraphFree(&heurdata->precedencegraph);
   }

   heurdata->initialized = FALSE;
   heurdata->njobs = 0;

   return SCIP_OKAY;
}

/** constructs a solution with the given start values for the integer start variables */
static
SCIP_RETCODE constructSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to be constructed */
   SCIP_VAR**            vars,               /**< integer start variables */
   int*                  starttimes,         /**< start times for the integer start variables */
   int                   nvars               /**< number of integer start variables */
   )
{
   SCIP_VAR* var;
   SCIP_Real val;
   int v;

   /* set start time variables */
   for( v = 0; v < nvars; ++v )
   {
      /* get some values */
      var = vars[v];
      val = (SCIP_Real)starttimes[v];

      SCIP_CALL( SCIPsetSolVal(scip, sol, var, val) );
   }

   return SCIP_OKAY;
}

/** insert given job into the profiles */
static
void profilesInsertJob(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROFILE**        profiles,           /**< array of resource profiles */
   int                   nprofiles,          /**< number of profiles */
   int                   starttime,          /**< start time of the job */
   int                   duration,           /**< duration of the job */
   int*                  demands,            /**< profile depending demands */
   SCIP_Bool*            infeasible          /**< pointer to store if the insertion is infeasible */
   )
{
   int pos;
   int p;

   /* found a feasible start time, insert the job into all profiles */
   for( p = 0; p < nprofiles && !(*infeasible); ++p )
   {
      /* add job to resource profile */
      SCIPprofileInsertCore(profiles[p], starttime, starttime + duration, demands[p], &pos, infeasible);
   }
}


/** retruns for the given job (duration and demands) the earliest feasible start time w.r.t. all profiles */
static
int profilesFindEarliestFeasibleStart(
   SCIP_PROFILE**       profiles,           /**< array of resource profiles */
   int                   nprofiles,          /**< number of profiles */
   int                   est,                /**< earliest start time */
   int                   lst,                /**< latest start time */
   int                   duration,           /**< duration of the job */
   int*                  demands,            /**< profile depending demands */
   SCIP_Bool*            infeasible          /**< pointer to store if it is infeasible to do */
   )
{
   SCIP_Bool changed;
   int start;
   int r;

   assert(!(*infeasible));

   do
   {
      changed = FALSE;

      for( r = 0; r < nprofiles; ++r )
      {
         assert(est >= 0);

         /* get next possible time to start from the current earliest starting time */
         start = SCIPprofileGetEarliestFeasibleStart(profiles[r], est, lst, duration, demands[r], infeasible);

         /* stop if job cannot be inserted */
         if( *infeasible )
         {
            SCIPdebugMessage("Terminate after start: resource %d, est %d, duration %d, demand %d\n",
               r, est, duration, demands[r]);
            return -1;
         }

         /* check if the earliest start time changes */
         if( r > 0 && start > est )
            changed = TRUE;

         est = start;
      }
   }
   while( changed );

   SCIPdebugMessage("earliest feasible start time: %d\n", est);

   return est;
}

/** retruns for the given job (duration and demands) the earliest feasible start time w.r.t. all profiles */
static
int profilesFindLatestFeasibleStart(
   SCIP_PROFILE**       profiles,           /**< array of resource profiles */
   int                   nprofiles,          /**< number of profiles */
   int                   lst,                /**< latest start time */
   int                   duration,           /**< duration of the job */
   int*                  demands,            /**< profile depending demands */
   SCIP_Bool*            infeasible          /**< pointer to store if it is infeasible to do */
   )
{
   SCIP_Bool changed;
   int start;
   int r;

   do
   {
      changed = FALSE;

      for( r = 0; r < nprofiles; ++r )
      {
         /* get next latest possible time to start from */
         start = SCIPprofileGetLatestFeasibleStart(profiles[r], 0, lst, duration, demands[r], infeasible);

         if( *infeasible )
         {
            SCIPdebugMessage("Terminate after start: resource %d, lst %d, duration %d, demand %d\n",
               r, lst, duration, demands[r]);
            return -1;
         }

         assert(start <= lst);

         /* check if the earliest start time changes */
         if( r > 0 && start < lst )
            changed = TRUE;

         lst = start;
      }
   }
   while( changed );

   return lst;
}

/** collect earliest and latest start times for all variables in the order given in the variables array */
static
void collectEstLst(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  ests,               /**< array to store the earliest start times */
   int*                  lsts,               /**< array to store the latest start times */
   int                   nvars               /**< number of variables */
   )
{
   SCIP_VAR* var;
   int j;

   assert(ests != NULL);
   assert(lsts != NULL);

   /* initialize earliest and latest start times */
   for( j = 0; j < nvars; ++j )
   {
      var = vars[j];
      assert(var != NULL);

      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         ests[j] = (int)(SCIPgetSolVal(scip, NULL, var) + 0.5);
      }
      else
      {
         ests[j] = (int)(SCIPvarGetLbLocal(var) + 0.5);
      }
      lsts[j] = (int)(SCIPvarGetUbGlobal(var) + 0.5);
      assert(ests[j] <= lsts[j]);
   }
}

/** propagate the earliest start time of the given job via the precedence graph to all successors jobs*/
static
void propagateEst(
   SCIP_DIGRAPH*         precedencegraph,    /**< precedence graph */
   int*                  ests,               /**< array of earliest start time for each job */
   int*                  lsts,               /**< array of latest start times for each job */
   int                   pred,               /**< index of the job which earliest start time showed propagated */
   int                   duration,           /**< duration of the job */
   SCIP_Bool*            infeasible          /**< pointer to store if the propagate detected an infeasibility */
   )
{
   int* successors;
   int nsuccessors;
   int succ;
   int ect;
   int s;

   nsuccessors = SCIPdigraphGetNSuccessors(precedencegraph, pred);
   successors = SCIPdigraphGetSuccessors(precedencegraph, pred);

   /* compute earliest completion time */
   ect = ests[pred] + duration;

   for( s = 0; s < nsuccessors && !(*infeasible); ++s )
   {
      succ = successors[s];

      /* check if the new earliest start time is smaller than the latest start time of the job */
      if( ect > lsts[succ] )
         *infeasible = TRUE;
      else
         ests[succ] = MAX(ests[succ], ect);
   }
}

/** propagate the latest start time of the given job via the precedence graph w.r.t. all successors jobs */
static
void propagateLst(
   SCIP_DIGRAPH*         precedencegraph,    /**< precedence graph */
   int*                  lsts,               /**< array of latest start times for each job */
   int                   pred,               /**< index of the job which earliest start time showed propagated */
   int                   duration            /**< duration of the job */
   )
{
   int* successors;
   int nsuccessors;
   int succ;
   int s;

   nsuccessors = SCIPdigraphGetNSuccessors(precedencegraph, pred);
   successors = SCIPdigraphGetSuccessors(precedencegraph, pred);

   for( s = 0; s < nsuccessors; ++s )
   {
      succ = successors[s];
      lsts[pred] = MIN(lsts[pred], lsts[succ] - duration);
   }
}

/** perform forward scheduling, that is, assigned jobs (in the given ordering) to their earliest start time, propagate
 *  w.r.t. the precedence graph and resource profiles
 */
static
SCIP_RETCODE performForwardScheduling(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int*                  starttimes,         /**< array to store the start times for each job */
   int*                  lsts,               /**< array of latest start times for each job */
   int*                  perm,               /**< permutation defining the order of the jobs */
   int*                  makespan,           /**< pointer to store the makespan of the forward scheduling solution */
   SCIP_Bool*            infeasible          /**< pointer to store if an infeasibility was detected */
   )
{
   SCIP_PROFILE** profiles;
   int nresources;
   int njobs;
   int j;

   nresources = heurdata->nresources;
   njobs = heurdata->njobs;
   *makespan = 0;

   assert(*infeasible == FALSE);

   SCIPdebugMessage("perform forward scheduling\n");

   /* create resource profiles for checking the resource requirements */
   SCIP_CALL( SCIPallocBufferArray(scip, &profiles, nresources) );
   for( j = 0; j < nresources; ++j )
   {
      assert(heurdata->capacities[j] > 0);
      SCIP_CALL( SCIPprofileCreate(&profiles[j], heurdata->capacities[j]) );
   }

   for( j = 0; j < njobs && !(*infeasible); ++j )
   {
      int* demands;
      int duration;
      int idx;

      idx = perm[j];
      assert(idx >= 0 && idx < njobs);

      duration = heurdata->durations[idx];
      demands = heurdata->resourcedemands[idx];
      assert(demands != NULL);

      /* skip jobs which have a duration of zero */
      if( duration > 0 )
      {
         /* find earliest start time w.r.t to all resource profiles */
         starttimes[idx] = profilesFindEarliestFeasibleStart(profiles, nresources, starttimes[idx], lsts[idx], duration, demands, infeasible);

         if( *infeasible )
            break;

         /* adjust makespan */
         (*makespan) = MAX(*makespan, starttimes[idx] + duration);

         /* insert the job into the profiles */
         profilesInsertJob(scip, profiles, nresources, starttimes[idx], duration, demands, infeasible);
         if( *infeasible )
            break;

         /* propagate the new earliest start time of the job */
         propagateEst(heurdata->precedencegraph, starttimes, lsts, idx, duration, infeasible);
      }
      SCIPdebugMessage("job %d -> est %d\n", idx, starttimes[idx]);
   }

   /* free resource profiles */
   for( j = 0; j < nresources; ++j )
   {
      SCIPprofileFree(&profiles[j]);
   }
   SCIPfreeBufferArray(scip, &profiles);

   SCIPdebugMessage("forward scheduling: makespan %d, feasible %d\n", *makespan, !(*infeasible));

   return SCIP_OKAY;
}

/** perform backward scheduling, that is, schedule jobs (in the ordering of their latest completion time) to their and
 *  propagate w.r.t. the precedence graph and resource profiles
 */
static
SCIP_RETCODE performBackwardScheduling(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int*                  starttimes,         /**< array of latest start times for each job */
   int*                  perm,               /**< permutation defining the order of jobs */
   SCIP_Bool*            infeasible          /**< pointer to store if an infeasibility was detected */
   )
{
   SCIP_PROFILE** profiles;
   int* durations;
   int* demands;
   int duration;
   int nresources;
   int njobs;
   int idx;
   int j;

   nresources = heurdata->nresources;
   njobs = heurdata->njobs;
   durations = heurdata->durations;

   SCIPdebugMessage("perform forward scheduling\n");

   /* create resource profiles for checking the resource requirements */
   SCIP_CALL( SCIPallocBufferArray(scip, &profiles, nresources) );
   for( j = 0; j < nresources; ++j )
   {
      assert(heurdata->capacities[j] > 0);
      SCIP_CALL( SCIPprofileCreate(&profiles[j], heurdata->capacities[j]) );
   }

   for( j = 0; j < njobs; ++j )
   {
      idx = perm[j];
      assert(idx >= 0 && idx < njobs);

      duration = durations[idx];
      demands = heurdata->resourcedemands[idx];
      assert(demands != NULL);

      /* propagate the new latest start time */
      propagateLst(heurdata->precedencegraph, starttimes, idx, duration);

      /* check against the resource profiles if the duration is greater than zero */
      if( duration > 0 )
      {
         /* find earliest start time w.r.t to all resource profiles */
         starttimes[idx] = profilesFindLatestFeasibleStart(profiles, nresources, starttimes[idx], duration, demands, infeasible);
         if( *infeasible )
            break;

         /* insert the job into the profiles */
         profilesInsertJob(scip, profiles, nresources, starttimes[idx], duration, demands, infeasible);
         if( *infeasible )
            break;
      }

      SCIPdebugMessage("job %d -> est %d\n", idx, starttimes[idx]);

   }

   /* free resource profiles */
   for( j = 0; j < nresources; ++j )
   {
      SCIPprofileFree(&profiles[j]);
   }
   SCIPfreeBufferArray(scip, &profiles);

   return SCIP_OKAY;
}

static
SCIP_RETCODE getEstPermutation(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  starttimes,         /**< array of start times for each job */
   int*                  ests,               /**< earliest start times */
   int*                  durations,          /**< array of durations */
   int*                  perm,               /**< array to store the permutation w.r.t. earliest start time */
   int                   njobs               /**< number of jobs */
   )
{
   int* sortingkeys;
   int j;

   SCIP_CALL( SCIPallocBufferArray(scip, &sortingkeys, njobs) );

   for( j = 0; j < njobs; ++j )
   {
      perm[j] = j;
      sortingkeys[j] = starttimes[j];
      starttimes[j] = ests[j];
   }
   SCIPsortIntInt(sortingkeys, perm, njobs);

   SCIPfreeBufferArray(scip, &sortingkeys);

   return SCIP_OKAY;
}

static
SCIP_RETCODE getLctPermuataion(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  starttimes,         /**< array of start times for each job */
   int*                  durations,          /**< array of durations */
   int*                  perm,               /**< array to store the permutation w.r.t. latest completion time */
   int                   njobs               /**< number of jobs */
   )
{
   int* sortingkeys;
   int j;

   SCIP_CALL( SCIPallocBufferArray(scip, &sortingkeys, njobs) );

   for( j = 0; j < njobs; ++j )
   {
      perm[j] = j;
      sortingkeys[j] = starttimes[j] + durations[j];
   }
   SCIPsortDownIntInt(sortingkeys, perm, njobs);

   SCIPfreeBufferArray(scip, &sortingkeys);

   return SCIP_OKAY;
}


/** execution method of heuristic  */
static
SCIP_RETCODE executeHeuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< Heuristic data structure */
   SCIP_Bool*            result              /**< pointer to store whether solution is found or not */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   int* starttimes;
   int* ests;
   int* lsts;
   int* perm;
   SCIP_Bool infeasible;
   SCIP_Bool stored;
   int makespan;
   int njobs;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->initialized);

   vars = heurdata->vars;
   njobs = heurdata->njobs;
   infeasible = FALSE;

   /* create initialized permutation */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Real* solvals;
      int v;

      SCIP_CALL( SCIPallocBufferArray(scip, &perm, njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &solvals, njobs) );

      /* in case the LP relaxation was solved to optimality we use the LP solution as initialized permutation */
      for( v = 0; v < njobs; ++v )
      {
         solvals[v] = SCIPgetSolVal(scip, NULL, vars[v]);
         perm[v] = v;
      }
      SCIPsortRealInt(solvals, perm, njobs);

      SCIPfreeBufferArray(scip, &solvals);
   }
   else
   {
      int* component;

      /* in case the LP was not solved we use the topologically sorted variables w.r.t. precedences graph */
      SCIPdigraphGetComponent(heurdata->precedencegraph, 0, &component, NULL);
      SCIP_CALL( SCIPduplicateBufferArray(scip, &perm, component, njobs) );
   }

   /* collect earliest and latest start times for all variables in the order given in the variables array */
   SCIP_CALL( SCIPallocBufferArray(scip, &ests, njobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lsts, njobs) );
   collectEstLst(scip, vars, ests, lsts, njobs);

   /* initialize the start times with the earliest start times */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &starttimes, ests, njobs) );

   /* LIST schedule:
    *
    * STEP 1: perform forward scheduling, that is, shift all jobs to the left as much as possible in the given ordering
    */
   SCIP_CALL( performForwardScheduling(scip, heurdata, starttimes, lsts, perm, &makespan, &infeasible) );

   if( !infeasible )
   {
      SCIP_SOL* sol;

      /* get permutation w.r.t. latest completion time given by start times */
      SCIP_CALL( getLctPermuataion(scip, starttimes, heurdata->durations, perm, njobs) );

      /* backward scheduling w.r.t. latest completion time */
      performBackwardScheduling(scip, heurdata, starttimes, perm, &infeasible);

      if( !infeasible )
      {
         /* get permutation w.r.t. earliest start time given by the starttimes and reset the start time to the earliest start time */
         SCIP_CALL( getEstPermutation(scip, starttimes, ests, heurdata->durations, perm, njobs) );

         SCIP_CALL( performForwardScheduling(scip, heurdata, starttimes, lsts, perm, &makespan, &infeasible) );

         SCIP_CALL( SCIPcreateOrigSol(scip, &sol, heur) );

         SCIP_CALL( constructSolution(scip, sol, vars, starttimes, njobs) );

         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, TRUE, TRUE, TRUE, &stored) );

         if( stored )
            *result = SCIP_FOUNDSOL;
      }
   }

   SCIPfreeBufferArray(scip, &starttimes);
   SCIPfreeBufferArray(scip, &lsts);
   SCIPfreeBufferArray(scip, &ests);

   SCIPfreeBufferArray(scip, &perm);

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods
 *
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyListScheduling NULL


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeListScheduling)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( heurdataFree(scip, heurdata) );

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
#define heurInitListScheduling NULL

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitListScheduling NULL

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolListScheduling NULL

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolListScheduling NULL

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecListScheduling)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;                  /* Primal heuristic data */

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   (*result) = SCIP_DIDNOTRUN;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( !heurdata->initialized )
      return SCIP_OKAY;

   SCIPdebugMessage("execute heuristic <"HEUR_NAME">\n");

   if( heurdata->njobs == 0 )
      return SCIP_OKAY;

   (*result) = SCIP_DIDNOTFIND;

   SCIP_CALL( executeHeuristic(scip, heur, result) );

   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 */

/** creates the list scheduling primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurListScheduling(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create list scheduling primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdata->resourcedemands = NULL;
   heurdata->vars = NULL;
   heurdata->njobs = 0;
   heurdata->initialized = FALSE;

   heur = NULL;
   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY,
         HEUR_FREQ, HEUR_FREQOFS, HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurExecListScheduling, heurdata) );
   assert(heur != NULL);

   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeListScheduling) );

   return SCIP_OKAY;
}

/** initialize heuristic */
SCIP_RETCODE SCIPinitializeHeurListScheduling(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH*         precedencegraph,    /**< precedence graph */
   SCIP_VAR**            vars,               /**< start time variables */
   int*                  durations,          /**< duration of the jobs independent of the resources */
   int**                 resourcedemands,    /**< resource demand matrix */
   int*                  capacities,         /**< resource capacities */
   int                   njobs,              /**< number if jobs */
   int                   nresources          /**< number of resources */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   heur = SCIPfindHeur(scip, HEUR_NAME);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( heurdata->initialized )
   {
      /* free old heuristic data structure */
      SCIP_CALL( heurdataFree(scip, heurdata) );
   }

   /* initialize the heuristic data structure */
   SCIP_CALL( heurdataInit(scip, heurdata, precedencegraph, vars, durations, resourcedemands, capacities, njobs, nresources) );

   return SCIP_OKAY;
}

/**@} */
