/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
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
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      1000LL    /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MAXPROPROUNDS -1        /**< maximum number of propagation rounds during probing */

/*
 * Data structures
 */

struct SCIP_Assignment
{
   SCIP_Bool**           vars;
   SCIP_Real**           solvals;
   SCIP_Bool*            feasibles;
   unsigned int*         keys;
   int*                  nones;
   int                   nassignments;
   int                   sassignments;
};
typedef struct SCIP_Assignment SCIP_ASSIGNMENT;


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

   SCIP_ASSIGNMENT**     machineassignments;
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
   heurdata->machineassignments = NULL;
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

   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);

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
#if 0
   int depth;
#endif

   assert(heur != NULL);
   assert(heurdata != NULL);

   /* initialize default values */
   infeasible = FALSE;

   *result = SCIP_DIDNOTFIND;
#if 0
   depth = SCIPgetDepth(scip);
#endif

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
      SCIP_ASSIGNMENT* machineassignment;
      int pos;

      SCIP_SOL* sol;
      SCIP_VAR** vars;
      SCIP_Real* lbs;
      SCIP_Real* ubs;
      int* durations;
      int* demands;
      SCIP_Bool unbounded;
      int njobs;
      int nvars;
      int j;
      int m;

      /* create temporary solution */
      SCIP_CALL( SCIPcreateOrigSol(scip, &sol, heur) );

      /* initialize the solution with the lower bound of all variables */
      SCIP_CALL( initializeSol(scip, sol) );

      njobs = heurdata->njobs;

      /* allocate memory for collecting the information for the single machines */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &durations, njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &demands, njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lbs, njobs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ubs, njobs) );

      nvars = -1;

      for( m = 0; m < heurdata->nmachines && !infeasible; ++m )
      {
         unsigned int key;
         int a;

         machineassignment = heurdata->machineassignments[m];

         pos = machineassignment->nassignments;

         /* realloc memory if not enough space left */
         if( machineassignment->nassignments == machineassignment->sassignments)
         {
            int oldsize;
            int newsize;

            oldsize = machineassignment->sassignments;
            newsize = SCIPcalcMemGrowSize(scip, pos + 1);

            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(machineassignment->vars), oldsize, newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(machineassignment->solvals), oldsize, newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(machineassignment->feasibles), oldsize, newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(machineassignment->keys), oldsize, newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(machineassignment->nones), oldsize, newsize) );

            machineassignment->sassignments = newsize;
         }
         assert(machineassignment->sassignments > pos);

         assert(njobs >= heurdata->machines[m]);
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &machineassignment->vars[pos], heurdata->machines[m]) ); /*lint !e866*/
         BMSclearMemoryArray(machineassignment->vars[pos], heurdata->machines[m]); /*lint !e866*/
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &machineassignment->solvals[pos], heurdata->machines[m]) ); /*lint !e866*/
         machineassignment->nassignments++;
         nvars = 0;
         key = 0;

         /* collect the jobs which are assign to that machine */
         for( j = 0; j < heurdata->machines[m]; ++j )
         {
            SCIP_VAR* binvar;

            binvar = heurdata->binvars[m][j];
            assert(binvar != NULL);

            /* check if job is assign to that machine */
            if( SCIPvarGetLbLocal(binvar) > 0.5 )
            {
               vars[nvars] = heurdata->vars[m][j];
               durations[nvars] = heurdata->durations[m][j];
               demands[nvars] = heurdata->demands[m][j];
               nvars++;

               machineassignment->vars[pos][j] = TRUE;
               key |= (1 << (j % 32)); /*lint !e701*/

               SCIP_CALL( SCIPsetSolVal(scip, sol, binvar, 1.0) );
            }
         }
         machineassignment->nones[pos] = nvars;
         machineassignment->keys[pos] = key;

         /* if none of the variables is assigned to that machine we skip it */
         if( nvars == 0 )
         {
            SCIPfreeBlockMemoryArray(scip, &machineassignment->vars[pos], heurdata->machines[m]); /*lint !e866*/
            SCIPfreeBlockMemoryArray(scip, &machineassignment->solvals[pos], heurdata->machines[m]); /*lint !e866*/
            machineassignment->nassignments--;
            continue;
         }

         /* check whether we already have try a subset of this variable combination */
         for( a = pos - 1; a >= 0; --a )
         {
            /* infeasible check */
            if( !machineassignment->feasibles[a]
               && nvars > machineassignment->nones[a] && ((~key & machineassignment->keys[a]) == 0) )
            {
               /* if we compare to an infeasible assignment, that assignment can be smaller or equal since a smaller
                * infeasible assignment induces a infeasibility for all assignments which include that assignment
                */

               /* do the expensive pairwise comparison */
               for( j = heurdata->machines[m] - 1; j >= 0; --j )
               {
                  /* at least the same variables in the old combination have to be assigned to 1 */
                  if( machineassignment->vars[pos][j] < machineassignment->vars[a][j] )
                     break;
               }
               /* we already tried this combination */
               if( j == -1 )
                  break;
            }
            /* feasible check */
            else if( machineassignment->feasibles[a] &&
               nvars < machineassignment->nones[a] && ((key & ~(machineassignment->keys[a])) == 0) )
            {
               /* if we compare to a feasible assignment, that assignment can be larger or equal since a larger feasible
                * assignment induces a feasibility for all assignments which is subset of that assignment
                */

               /* do the expensive pairwise comparison */
               for( j = heurdata->machines[m] - 1; j >= 0; --j )
               {
                  if( machineassignment->vars[pos][j] > machineassignment->vars[a][j] )
                     break;
               }
               /* we already tried this combination */
               if( j == -1 )
                  break;
            }
            else if( nvars == machineassignment->nones[a] && ((~key & machineassignment->keys[a]) == 0) )
            {
               /* do the expensive pairwise comparison */
               for( j = heurdata->machines[m] - 1; j >= 0; --j )
               {
                  if( machineassignment->vars[pos][j] != machineassignment->vars[a][j] )
                     break;
               }
               /* we already tried this combination */
               if( j == -1 )
                  break;
            }
         }

         if( a >= 0 )
         {
            SCIPdebugMessage("We already tried %s this combination, it was %s\n",
               machineassignment->nones[pos] > machineassignment->nones[a] ? "a subset of" : (machineassignment->nones[pos] > machineassignment->nones[a] ? "a superset of" : ""),
               machineassignment->feasibles[a] ? "feasible" : "infeasible");

            /* delete unnecessary data */
            SCIPfreeBlockMemoryArray(scip, &machineassignment->vars[pos], heurdata->machines[m]); /*lint !e866*/
            SCIPfreeBlockMemoryArray(scip, &machineassignment->solvals[pos], heurdata->machines[m]); /*lint !e866*/
            machineassignment->nassignments--;

            infeasible = !machineassignment->feasibles[a];

            if( infeasible )
               break;

            for( j = 0; j < heurdata->machines[m]; ++j )
            {
               if( machineassignment->vars[a][j] && SCIPvarGetLbLocal(heurdata->binvars[m][j]) > 0.5 )
               {
                  SCIP_CALL( SCIPsetSolVal(scip, sol, heurdata->vars[m][j], machineassignment->solvals[a][j]) );
               }
            }
         }
         else
         {
            SCIP_Real* objvals;
            SCIP_Real timelimit;
            SCIP_Real memorylimit;
            SCIP_Bool solved;
            SCIP_Bool error;
            int v;

            SCIPdebugMessage("check machine %d (variables %d)\n", m, nvars);

            SCIP_CALL( SCIPallocBufferArray(scip, &objvals, nvars) );

            for( v = 0; v < nvars; ++v )
            {
               SCIP_VAR* var;

               var = vars[v];
               assert(var != NULL);

               lbs[v] = SCIPvarGetLbLocal(var);
               ubs[v] = SCIPvarGetUbLocal(var);
               objvals[v] = SCIPvarGetObj(var);
            }

            /* check whether there is enough time and memory left */
            SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
            if( !SCIPisInfinity(scip, timelimit) )
               timelimit -= SCIPgetSolvingTime(scip);
            SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

            /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
            if( !SCIPisInfinity(scip, memorylimit) )
            {
               memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
               memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
            }

            /* solve the cumulative condition separately */
            SCIP_CALL( SCIPsolveCumulative(scip, nvars, lbs, ubs, objvals, durations, demands, heurdata->capacities[m], 0, INT_MAX,
                  timelimit, memorylimit, heurdata->maxnodes, &solved, &infeasible, &unbounded, &error) );
            assert(!unbounded);
            assert(!error);

            SCIPfreeBufferArray(scip, &objvals);

            machineassignment->feasibles[pos] = !infeasible;

            if( infeasible )
            {
               SCIPdebugMessage("infeasible :-(\n");
               break;
            }

            for( j = 0, v = 0; j < heurdata->machines[m]; ++j )
            {
               if( machineassignment->vars[pos][j] && SCIPvarGetLbLocal(heurdata->binvars[m][j]) > 0.5 )
               {
                  SCIP_CALL( SCIPsetSolVal(scip, sol, heurdata->vars[m][j], lbs[v]) );
                  machineassignment->solvals[pos][j] = lbs[v];
                  v++;
               }
            }
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

         SCIPdebugMessage("************ try solution <%g>\n", SCIPgetSolOrigObj(scip, sol));

         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, FALSE, FALSE, TRUE, &stored) );

         if( stored )
            *result = SCIP_FOUNDSOL;
      }
#if 0
      else
      {
         /* check that code */
         int v;

         SCIP_CALL( SCIPinitConflictAnalysis(scip) );

         for( v = 0; v < heurdata->machines[m]; ++v )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, heurdata->binvars[m][v]) );
            SCIP_CALL( SCIPaddConflictLb(scip, heurdata->vars[m][v], NULL) );
            SCIP_CALL( SCIPaddConflictUb(scip, heurdata->vars[m][v], NULL) );
         }

         /* analyze the conflict */
#if 0
         SCIP_CALL( SCIPanalyzeConflict(scip, depth, NULL) );
#endif
         SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
      }
#endif
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
   for( m = heurdata->nmachines - 1; m >= 0; --m )
   {
      int a;

      for( a = 0; a < heurdata->machineassignments[m]->nassignments; ++a )
      {
         SCIPfreeBlockMemoryArray(scip, &(heurdata->machineassignments[m]->vars[a]), heurdata->machines[m]); /*lint !e866*/
         SCIPfreeBlockMemoryArray(scip, &(heurdata->machineassignments[m]->solvals[a]), heurdata->machines[m]); /*lint !e866*/
      }

      SCIPfreeBlockMemoryArray(scip, &(heurdata->machineassignments[m]->nones), heurdata->machineassignments[m]->sassignments);
      SCIPfreeBlockMemoryArray(scip, &(heurdata->machineassignments[m]->keys), heurdata->machineassignments[m]->sassignments);
      SCIPfreeBlockMemoryArray(scip, &(heurdata->machineassignments[m]->feasibles), heurdata->machineassignments[m]->sassignments);
      SCIPfreeBlockMemoryArray(scip, &(heurdata->machineassignments[m]->solvals), heurdata->machineassignments[m]->sassignments);
      SCIPfreeBlockMemoryArray(scip, &(heurdata->machineassignments[m]->vars), heurdata->machineassignments[m]->sassignments);
      SCIPfreeBlockMemory(scip, &heurdata->machineassignments[m]); /*lint !e866*/

      SCIPfreeBlockMemoryArray(scip, &heurdata->vars[m], heurdata->machines[m]);
      SCIPfreeBlockMemoryArray(scip, &heurdata->binvars[m], heurdata->machines[m]);
      SCIPfreeBlockMemoryArray(scip, &heurdata->durations[m], heurdata->machines[m]);
      SCIPfreeBlockMemoryArray(scip, &heurdata->demands[m], heurdata->machines[m]);
   }

   /* free arrays */
   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->machineassignments, heurdata->nmachines);
   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->demands, heurdata->nmachines);
   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->durations, heurdata->nmachines);
   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->binvars, heurdata->nmachines);
   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->vars, heurdata->nmachines);

   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->capacities, heurdata->nmachines);
   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->machines, heurdata->nmachines);

   SCIPfreeBlockMemory(scip, &heurdata);
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

   SCIPdebugMessage("apply optcumulative heuristic at node %"SCIP_LONGINT_FORMAT"\n",
         SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

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
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
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
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->vars, nmachines) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->binvars, nmachines) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->durations, nmachines) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->demands, nmachines) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->machineassignments, nmachines) );

   for( m = 0; m < nmachines; ++m )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &heurdata->vars[m], vars[m], machines[m]) ); /*lint !e866*/
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &heurdata->binvars[m], binvars[m], machines[m]) ); /*lint !e866*/
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &heurdata->durations[m], durations[m], machines[m]) ); /*lint !e866*/
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &heurdata->demands[m], demands[m], machines[m]) ); /*lint !e866*/

      /* sort variable w.r.t. their objective coefficient */
      SCIPsortPtrPtrIntInt((void**)heurdata->binvars[m], (void**)heurdata->vars[m],
         heurdata->durations[m], heurdata->demands[m], SCIPvarCompObj, machines[m]);

      SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata->machineassignments[m]) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(heurdata->machineassignments[m]->vars), njobs) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(heurdata->machineassignments[m]->solvals), njobs) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(heurdata->machineassignments[m]->feasibles), njobs) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(heurdata->machineassignments[m]->keys), njobs) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(heurdata->machineassignments[m]->nones), njobs) );
      heurdata->machineassignments[m]->nassignments = 0;
      heurdata->machineassignments[m]->sassignments = njobs;
   }

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &heurdata->machines, machines, nmachines) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &heurdata->capacities, capacities, nmachines) );

   heurdata->nmachines = nmachines;
   heurdata->njobs = njobs;
   heurdata->initialized = TRUE;

   return SCIP_OKAY;
}
