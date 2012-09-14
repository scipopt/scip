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

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "heur_optcumulative.h"


#define HEUR_NAME             "optcumulative"
#define HEUR_DESC             "LNS heuristic uses the variable lower and upper bounds to determine the search neighborhood"
#define HEUR_DISPCHAR         'q'
#define HEUR_PRIORITY         -1106000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE           /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINFIXINGRATE 0.5       /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which optcumulative heuristic should at least improve the incumbent          */
#define DEFAULT_MINNODES      500LL     /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_NODESOFS      5000LL    /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_MAXPROPROUNDS -1        /* maximum number of propagation rounds during probing */

/* enable statistic output by defining macro STATISTIC_INFORMATION */
#ifdef STATISTIC_INFORMATION
#define STATISTIC(x)                x
#else
#define STATISTIC(x)             /**/
#endif

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_VAR***           vars;               /**< job machine choice matrix */
   int                   njobs;              /**< number of jobs */
   int*                  nmachines;          /**< number of machines */

   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          usednodes;          /**< nodes already used by optcumulative heuristic in earlier calls                         */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Real             minimprove;         /**< factor by which optcumulative heuristic should at least improve the incumbent          */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   int                   maxproprounds;      /**< maximum number of propagation rounds during probing */
   SCIP_Bool             initialized;        /**< are the candidate list initialized? */
   SCIP_Bool             applicable;         /**< is the heuristic applicable? */
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
   heurdata->njobs = 0;
   heurdata->nmachines = NULL;
   heurdata->initialized = FALSE;
   heurdata->applicable = FALSE;
}

/** apply variable bound fixing during probing */
static
SCIP_RETCODE applyOptcumulativeFixings(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< structure containing heurdata */
   SCIP_Bool*            infeasible          /**< pointer to store whether problem is infeasible */
   )
{
   SCIP_VAR*** vars;
   SCIP_VAR* var;
   int njobs;
   int nmachines;
   int j;
   int m;

   vars = heurdata->vars;
   njobs = heurdata->njobs;

   for( j = 0; j < njobs && !(*infeasible); ++j )
   {
      nmachines = heurdata->nmachines[j];

      for( m = 0; m < nmachines && !(*infeasible); ++m )
      {
         var = vars[j][m];

         /* skip variables which are already fixed */
         if( SCIPvarGetLbLocal(var) + 0.5 > SCIPvarGetUbLocal(var) )
            continue;

         SCIP_CALL( SCIPnewProbingNode(scip) );

         SCIP_CALL( SCIPfixVarProbing(scip, var, 1.0) );

         SCIPdebugMessage("fixing %d: variable <%s> objective coefficient <%g> (%d pseudo cands)\n",
            j, SCIPvarGetName(var), SCIPvarGetObj(var), SCIPgetNPseudoBranchCands(scip));

         /* check if problem is already infeasible */
         SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, infeasible, NULL) );

         if( !(*infeasible) )
         {
#if 0
            SCIP_Bool lperror;

            SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror) );

            fixed = TRUE;

            if( !lperror && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
#endif
               break;
         }
         else if( SCIPvarGetLbLocal(var) + 0.5 > SCIPvarGetUbLocal(var) )
            break;

         /* backtrack */
         SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );

         /* after backtracking should the variable be not be fixed */
         assert(SCIPvarGetLbLocal(var) + 0.5 < SCIPvarGetUbLocal(var));

         SCIP_CALL( SCIPfixVarProbing(scip, var, 0.0) );

         SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, infeasible, NULL) );
      }
   }

   SCIPdebugMessage("probing ended with %sfeasible problem\n", (*infeasible) ? "in" : "");

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_SOL*             newsol,             /**< working solution */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( subvars != NULL );
   assert( subsol != NULL );

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* subSCIP may have more variable than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert( nvars <= SCIPgetNOrigVars(subscip) );

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** main procedure of the optcumulative heuristic */
static
SCIP_RETCODE applyOptcumulative(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_Real             timelimit,          /**< timelimit for the subproblem                                   */
   SCIP_Real             memorylimit,        /**< memorylimit for the subproblem                                 */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem                    */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_SOL* newsol;
   SCIP_Bool infeasible;

   assert(heur != NULL);
   assert(heurdata != NULL);

   /* initialize default values */
   infeasible = FALSE;

   *result = SCIP_DIDNOTFIND;

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* create temporary solution */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

   /* apply the variable fixings */
   SCIP_CALL( applyOptcumulativeFixings(scip, heurdata, &infeasible) );

   /* if a solution has been found --> fix all other variables by subscip if necessary */
   if( !infeasible )
   {
#if 0
      SCIP_Bool success;

      SCIP_CALL( SCIPlinkCurrentSol(scip, newsol) );

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, FALSE, TRUE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

#else
      SCIP* subscip;
      SCIP_VAR** vars;
      SCIP_VAR** subvars;
      SCIP_HASHMAP* varmap;
      SCIP_Bool valid;
      int nvars;
      int i;

      valid = FALSE;
      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);

      /* create subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

      SCIP_CALL( SCIPcopy(scip, subscip, varmap, NULL, "_optcumulative", FALSE, FALSE, FALSE, &valid) );

      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

      for( i = 0; i < nvars; i++ )
         subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);

      /* free hash map */
      SCIPhashmapFree(&varmap);

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

      /* set limits for the subproblem */
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) );
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

      /* forbid call of heuristics and separators solving sub-CIPs */
      SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

      /* disable cutting plane separation */
      SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

      /* disable expensive presolving */
      SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

#ifdef SCIP_DEBUG
      /* for debugging optcumulative heuristic, enable MIP output */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
      SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#endif

      /* if there is already a solution, add an objective cutoff */
      if( SCIPgetNSols(scip) > 0 )
      {
         SCIP_Real upperbound;
         SCIP_Real minimprove;
         SCIP_Real cutoff;

         minimprove = heurdata->minimprove;
         cutoff = SCIPinfinity(scip);
         assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

         upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

         if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
         {
            cutoff = (1-minimprove)*SCIPgetUpperbound(scip) + minimprove*SCIPgetLowerbound(scip);
         }
         else
         {
            if( SCIPgetUpperbound ( scip ) >= 0 )
               cutoff = ( 1 - minimprove ) * SCIPgetUpperbound ( scip );
            else
               cutoff = ( 1 + minimprove ) * SCIPgetUpperbound ( scip );
         }
         cutoff = MIN(upperbound, cutoff);
         SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );
      }

      /* solve the subproblem */
      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
       */
#ifdef NDEBUG
      {
         SCIP_RETCODE retstat;
         retstat = SCIPpresolve(subscip);
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while presolving subMIP in optcumulative heuristic; subSCIP terminated with code <%d>\n", retstat);
         }
      }
#else
      SCIP_CALL( SCIPpresolve(subscip) );
#endif

      SCIPdebugMessage("optcumulative heuristic presolved subproblem: %d vars, %d cons\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip));

      /* after presolving, we should have at least reached a certain fixing rate over ALL variables (including continuous)
       * to ensure that not only the MIP but also the LP relaxation is easy enough
       */
      if( ( nvars - SCIPgetNVars(subscip) ) / (SCIP_Real)nvars >= heurdata->minfixingrate / 2.0 )
      {
         SCIP_SOL** subsols;
         SCIP_Bool success;
         int nsubsols;

         SCIPdebugMessage("solving subproblem: nstallnodes=%"SCIP_LONGINT_FORMAT", maxnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->maxnodes);

#ifdef NDEBUG
         {
            SCIP_RETCODE retstat;
            retstat = SCIPsolve(subscip);
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving subMIP in optcumulative heuristic; subSCIP terminated with code <%d>\n",retstat);
            }
         }
#else
         SCIP_CALL( SCIPsolve(subscip) );
#endif

         /* check, whether a solution was found; due to numerics, it might happen that not all solutions are feasible ->
          * try all solutions until one was accepted
          */
         nsubsols = SCIPgetNSols(subscip);
         subsols = SCIPgetSols(subscip);
         success = FALSE;

         for( i = 0; i < nsubsols && !success; ++i )
         {
            SCIP_CALL( createNewSol(scip, subscip, subvars, newsol, subsols[i], &success) );
         }
         if( success )
            *result = SCIP_FOUNDSOL;
      }

#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

      /* free subproblem */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
   }
#endif
   /* free solution */
   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

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
   int j;

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* release all variables */
   for( j = 0; j < heurdata->njobs; ++j )
   {
#if 0
      int m;

      for( m = 0; m < heurdata->nmachines; ++m )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &heurdata->vars[j][m]) );
      }
#endif
      SCIPfreeMemoryArrayNull(scip, &heurdata->vars[j]);
   }

   /* free varbounds array */
   SCIPfreeMemoryArrayNull(scip, &heurdata->vars);
   SCIPfreeMemoryArrayNull(scip, &heurdata->nmachines);

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
   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNPseudoBranchCands(scip) == 0 )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( !heurdata->initialized || !heurdata->applicable )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward variable bounds heuristic if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMessage("skipping "HEUR_NAME": nstallnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( timelimit < 10.0 || memorylimit <= 0.0 )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMessage("apply optcumulative heuristic at node %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   *result = SCIP_DIDNOTFIND;

   /* try variable lower and upper bounds which respect to objective coefficients */
   SCIP_CALL( applyOptcumulative(scip, heur, heurdata, timelimit, memorylimit, nstallnodes, result) );

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
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
         "minimum percentage of integer variables that have to be fixable ",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which "HEUR_NAME" heuristic should at least improve the incumbent  ",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxproprounds",
         "maximum number of propagation rounds during probing (-1 infinity)",
         &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX/4, NULL, NULL) );

   return SCIP_OKAY;
}

/** comparison method for sorting variables by non-decreasing index */
static
SCIP_DECL_SORTPTRCOMP(varComp)
{
   SCIP_Real diff;

   diff = SCIPvarGetObj((SCIP_VAR*)elem1) - SCIPvarGetObj((SCIP_VAR*)elem2);

   if( diff >= 0.5 )
      return 1;
   else if( diff <= -0.5 )
      return -1;
   else
      return 0;
}


/** initialize the heuristics data structure */
SCIP_RETCODE SCIPinitHeurOptcumulative(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_VAR***           vars,               /**< job machine matrix */
   int                   njobs,              /**< number of jobs */
   int                   nmachines           /**< number of machines */
   )
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   int j;
   int m;

   heur = SCIPfindHeur(scip, HEUR_NAME);

   if( heur == NULL )
   {
      SCIPerrorMessage("optcumulative heuristic not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   heurdata->njobs = njobs;

   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->vars, njobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->nmachines, njobs) );

   {
      SCIP_Real* sums;

      SCIP_CALL( SCIPallocBufferArray(scip, &sums, njobs) );

      for( j = 0; j < njobs; ++j )
      {
         heurdata->nmachines[j] = 0;
         sums[j] = 0.0;

         SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->vars[j], nmachines) );
         for( m = 0; m < nmachines; ++m )
         {
            if( vars[m][j] == NULL )
               continue;

            heurdata->vars[j][heurdata->nmachines[j]] = vars[m][j];

            SCIP_CALL( SCIPcaptureVar(scip, vars[m][j]) );

            heurdata->nmachines[j]++;
            sums[j] += SCIPvarGetObj(vars[m][j]);
         }

         if(  heurdata->nmachines[j] > 0 )
            sums[j] /= heurdata->nmachines[j];


         SCIPsortPtr( (void**)heurdata->vars[j], varComp, heurdata->nmachines[j]);

#ifdef SCIP_DEBUG
         printf("job %2d: %g <%s> ", j, SCIPvarGetObj(heurdata->vars[j][0]), SCIPvarGetName(heurdata->vars[j][0]));
         for(m = 1;  m < heurdata->nmachines[j]; ++m )
         {
            printf("%g <%s> ", SCIPvarGetObj(heurdata->vars[j][m]), SCIPvarGetName(heurdata->vars[j][m]));
            assert(SCIPvarGetObj(heurdata->vars[j][m-1]) <= SCIPvarGetObj(heurdata->vars[j][m]));
         }

         printf(" --> %g\n", sums[j]);
#endif
      }

      SCIPsortDownRealIntPtr(sums, heurdata->nmachines, (void**)heurdata->vars, njobs);

#ifdef SCIP_DEBUG
      for( j = 0; j < njobs; ++j )
      {
         printf("job %2d: %g ", j, SCIPvarGetObj(heurdata->vars[j][0]));
         for(m = 1;  m < heurdata->nmachines[j]; ++m )
         {
            printf("%g <%s> ", SCIPvarGetObj(heurdata->vars[j][m]), SCIPvarGetName(heurdata->vars[j][m]));
            assert(SCIPvarGetObj(heurdata->vars[j][m-1]) <= SCIPvarGetObj(heurdata->vars[j][m]));
         }

         printf(" --> %g\n", sums[j]);
      }
#endif

      SCIPfreeBufferArray(scip, &sums);

   }

   heurdata->initialized = TRUE;
   heurdata->applicable = TRUE;

   return SCIP_OKAY;
}
