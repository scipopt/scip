/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_optcumulative.h
 * @ingroup PRIMALHEURISTICS
 * @brief  heuristic for cumulative scheduling with optional activities
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_cumulative.h"
#include "heur_optcumulative.h"


#define HEUR_NAME             "optcumulative"
#define HEUR_DESC             "problem specific heuristic of cumulative scheduling problems with optional jobs"
#define HEUR_DISPCHAR         'q'
#define HEUR_PRIORITY         -1106000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      1000LL    /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MAXPROPROUNDS -1        /**< maximum number of propagation rounds during probing */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_VAR***           binvars;            /**< machnine job matrix (choice variables) */
   SCIP_VAR***           vars;               /**< machnine job matrix (start time variables) */
   int**                 durations;          /**< machnine job duration matrix */
   int**                 demands;            /**< machnine job demands matrix */
   int*                  machines;           /**< number of jobs for each machines */
   int*                  capacities;         /**< machine capacities */
   int                   nmachines;          /**< number of machines */
   int                   njobs;              /**< number of njobs */

   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   int                   maxproprounds;      /**< maximum number of propagation rounds during probing */
   SCIP_Bool             initialized;        /**< are the candidate list initialized? */
};

/*
 * Local methods
 */

/** reset heuristic data structure */
static
void heurdataReset(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< structure containing heurdata */
   )
{
   heurdata->vars = NULL;
   heurdata->binvars = NULL;
   heurdata->durations = NULL;
   heurdata->demands = NULL;
   heurdata->machines = NULL;
   heurdata->capacities = NULL;
   heurdata->nmachines = 0;
   heurdata->njobs = 0;

   heurdata->initialized = FALSE;
}

/** apply variable bound fixing during probing */
static
SCIP_RETCODE applyOptcumulativeFixings(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< structure containing heurdata */
   SCIP_Bool*            infeasible          /**< pointer to store whether problem is infeasible */
   )
{
   SCIP_VAR*** binvars;
   int* machines;
   int* possitions;
   int nmachines;
   int j;
   int m;

   binvars = heurdata->binvars;
   nmachines = heurdata->nmachines;
   machines = heurdata->machines;

   SCIP_CALL( SCIPallocBufferArray(scip, &possitions, nmachines) );
   BMSclearMemoryArray(possitions, nmachines);

   while( !(*infeasible) )
   {
      SCIP_VAR* var;
      SCIP_Real objval;
      int bestmachine;

      bestmachine = -1;
      objval = SCIPinfinity(scip);

      /* search over all machines and find the next cheapest job to assign */
      for( m = 0; m < nmachines; ++m )
      {
         int currentpos;

         currentpos = possitions[m];

         /* find next unfixed variable for the current machine */
         for( j = currentpos; j < machines[m]; ++j )
         {
            if( SCIPvarGetLbLocal(binvars[m][j]) + 0.5 < SCIPvarGetUbLocal(binvars[m][j]) )
               break;

            possitions[m]++;
         }

         currentpos = possitions[m];

         /* check if we have a variable left on that machine */
         if( currentpos < machines[m] )
         {
            assert(binvars[m][currentpos] != NULL);

            /* check if the objective coefficient is better than the best known one */
            if( SCIPvarGetObj(binvars[m][currentpos]) < objval )
            {
               objval = SCIPvarGetObj(binvars[m][currentpos]);
               bestmachine = m;
            }
         }
      }

      /* check if unsigned variable was left */
      if( bestmachine == -1 )
         break;

      assert(bestmachine < nmachines);
      assert(possitions[bestmachine] < machines[bestmachine]);

      var = binvars[bestmachine][possitions[bestmachine]];
      assert(var != NULL);
      assert(SCIPvarGetLbLocal(var) + 0.5 < SCIPvarGetUbLocal(var));

      possitions[bestmachine]++;

      SCIP_CALL( SCIPnewProbingNode(scip) );

      SCIP_CALL( SCIPfixVarProbing(scip, var, 1.0) );

      SCIPdebugMessage("variable <%s> objective coefficient <%g> fixed to 1.0 (%d pseudo cands)\n",
         SCIPvarGetName(var), SCIPvarGetObj(var), SCIPgetNPseudoBranchCands(scip));

      /* check if problem is already infeasible */
      SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, infeasible, NULL) );

      if( *infeasible )
      {
         /* backtrack */
         SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );

         /* after backtracking the variable might be already fixed to zero */
         if( SCIPvarGetUbLocal(var) > 0.5 )
         {
            SCIP_CALL( SCIPfixVarProbing(scip, var, 0.0) );
         }

         SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, infeasible, NULL) );
      }
   }

   SCIPfreeBufferArray(scip, &possitions);

   SCIPdebugMessage("probing ended with %sfeasible problem\n", (*infeasible) ? "in" : "");

   return SCIP_OKAY;
}

/** initialize the solution by assign the lower bound of the variable as solution value */
static
SCIP_RETCODE initializeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution to be initialize */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, vars[v], SCIPvarGetLbLocal(vars[v])) );
   }

   return SCIP_OKAY;
}

/** main procedure of the optcumulative heuristic */
static
SCIP_RETCODE applyOptcumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_Real lowerbound;
   SCIP_Real upperbound;
   SCIP_Real pseudoobj;
   SCIP_Bool infeasible;

   assert(heur != NULL);
   assert(heurdata != NULL);

   /* initialize default values */
   infeasible = FALSE;

   *result = SCIP_DIDNOTFIND;

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* apply the variable fixings */
   SCIP_CALL( applyOptcumulativeFixings(scip, heurdata, &infeasible) );

   lowerbound =  SCIPgetLowerbound(scip);
   upperbound =  SCIPgetUpperbound(scip);
   pseudoobj = SCIPgetPseudoObjval(scip);

   /* if a solution has been found --> fix all other variables by subscip if necessary */
   if( !infeasible && pseudoobj >= lowerbound && pseudoobj < upperbound )
   {
      SCIP_SOL* sol;
      SCIP_VAR** vars;
      SCIP_Real* lbs;
      SCIP_Real* ubs;
      int* durations;
      int* demands;
      SCIP_Bool unbounded;
      int njobs;
      int j;
      int m;

      /* create temporary solution */
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

      /* initialize the solution with the lower bound of all variables */
      SCIP_CALL( initializeSol(scip, sol) );

      njobs = heurdata->njobs;

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &durations, njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &demands, njobs) );

      SCIP_CALL( SCIPallocBufferArray(scip, &lbs, njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ubs, njobs) );

      for( m = 0; m < heurdata->nmachines && !infeasible; ++m )
      {
         int nvars;

         nvars = 0;

         for( j = 0; j < heurdata->machines[m]; ++j )
         {
            SCIP_VAR* var;

            var = heurdata->binvars[m][j];
            assert(var != NULL);

            /* check if job is assign to that machine */
            if( SCIPvarGetLbLocal(var) > 0.5 )
            {
               vars[nvars] = heurdata->vars[m][j];
               durations[nvars] = heurdata->durations[m][j];
               demands[nvars] = heurdata->demands[m][j];
               nvars++;

               SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
            }
         }

         SCIPdebugMessage("check machine %d (variables %d)\n", m, nvars);

         if( nvars == 0 )
            continue;

         /* solve the cumulative condition separately */
         SCIP_CALL( SCIPsolveCumulativeCondition(scip, nvars, vars, durations, demands,
               heurdata->capacities[m], 0, INT_MAX, lbs, ubs, heurdata->maxnodes, &infeasible, &unbounded) );
         assert(!unbounded);

         if( infeasible )
         {
            SCIPdebugMessage("infeasible :-(\n");
            break;
         }

         for( j = 0; j < nvars; ++j )
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, vars[j], lbs[j]) );
         }
      }

      SCIPfreeBufferArray(scip, &ubs);
      SCIPfreeBufferArray(scip, &lbs);

      SCIPfreeBufferArray(scip, &demands);
      SCIPfreeBufferArray(scip, &durations);
      SCIPfreeBufferArray(scip, &vars);

      /* try and free solution */
      if( !infeasible )
      {
         SCIP_Bool stored;

         SCIPdebugMessage("************ try solution\n");

         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, FALSE, TRUE, &stored) );

         if( stored )
            *result = SCIP_FOUNDSOL;
      }
      else
      {
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
      }
   }


   /* exit probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyOptcumulative)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of heuristic */
   SCIP_CALL( SCIPincludeHeurOptcumulative(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeOptcumulative)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int m;

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* release all variables */
   for( m = 0; m < heurdata->nmachines; ++m )
   {
      SCIPfreeMemoryArray(scip, &heurdata->vars[m]);
      SCIPfreeMemoryArray(scip, &heurdata->binvars[m]);
      SCIPfreeMemoryArray(scip, &heurdata->durations[m]);
      SCIPfreeMemoryArray(scip, &heurdata->demands[m]);
   }

   /* free varbounds array */
   SCIPfreeMemoryArrayNull(scip, &heurdata->vars);
   SCIPfreeMemoryArrayNull(scip, &heurdata->binvars);
   SCIPfreeMemoryArrayNull(scip, &heurdata->durations);
   SCIPfreeMemoryArrayNull(scip, &heurdata->demands);
   SCIPfreeMemoryArrayNull(scip, &heurdata->machines);
   SCIPfreeMemoryArrayNull(scip, &heurdata->capacities);

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
#define heurInitOptcumulative NULL

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitOptcumulative NULL

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolOptcumulative NULL

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolOptcumulative  NULL

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecOptcumulative)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNPseudoBranchCands(scip) == 0 )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( !heurdata->initialized )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMessage("apply optcumulative heuristic at node %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   *result = SCIP_DIDNOTFIND;

   /* try variable lower and upper bounds which respect to objective coefficients */
   SCIP_CALL( applyOptcumulative(scip, heur, heurdata, result) );

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the optcumulative primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurOptcumulative(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create optcumulative primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdataReset(scip, heurdata);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyOptcumulative,
         heurFreeOptcumulative, heurInitOptcumulative, heurExitOptcumulative,
         heurInitsolOptcumulative, heurExitsolOptcumulative, heurExecOptcumulative,
         heurdata) );

   /* add variable bounds primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxproprounds",
         "maximum number of propagation rounds during probing (-1 infinity)",
         &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX/4, NULL, NULL) );

   return SCIP_OKAY;
}

/** initialize the heuristics data structure */
SCIP_RETCODE SCIPinitHeurOptcumulative(
   SCIP*                 scip,               /**< original SCIP data structure */
   int                   nmachines,          /**< number of machines */
   int                   njobs,              /**< number of njobs */
   int*                  machines,           /**< number of jobs for each machines */
   SCIP_VAR***           binvars,            /**< machnine job matrix (choice variables) */
   SCIP_VAR***           vars,               /**< machnine job matrix (start time variables) */
   int**                 durations,          /**< machnine job duration matrix */
   int**                 demands,            /**< machnine job demands matrix */
   int*                  capacities          /**< machine capacities */
   )
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   int m;

   heur = SCIPfindHeur(scip, HEUR_NAME);

   if( heur == NULL )
   {
      SCIPerrorMessage("optcumulative heuristic not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* copy the problem data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->vars, nmachines) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->binvars, nmachines) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->durations, nmachines) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->demands, nmachines) );

   for( m = 0; m < nmachines; ++m )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->vars[m], vars[m], machines[m]) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->binvars[m], binvars[m], machines[m]) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->durations[m], durations[m], machines[m]) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->demands[m], demands[m], machines[m]) );

      /* sort variable w.r.t. their objective coefficient */
      SCIPsortPtrPtrIntInt((void**)heurdata->binvars[m], (void**)heurdata->vars[m],
         heurdata->durations[m], heurdata->demands[m], SCIPvarCompObj, machines[m]);
   }

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->machines, machines, nmachines) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->capacities, capacities, nmachines) );

   heurdata->nmachines = nmachines;

   heurdata->initialized = TRUE;

   return SCIP_OKAY;
}
