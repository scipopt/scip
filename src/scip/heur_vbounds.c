/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_vbounds.c
 * @brief  LNS heuristic uses the variable lower and upper bounds to determine the search neighborhood
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Jens Schulz
 * @author Gerald Gamrath
 *
 * @todo allow smaller fixing rate for probing LP?
 * @todo allow smaller fixing rate after presolve if total number of variables is small (<= 1000)?
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_vbounds.h"

#ifdef SCIP_STATISTIC
#include "scip/clock.h"
#endif

#define HEUR_NAME             "vbounds"
#define HEUR_DESC             "LNS heuristic uses the variable lower and upper bounds to determine the search neighborhood"
#define HEUR_DISPCHAR         'V'
#define HEUR_PRIORITY         -1106000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE           /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINFIXINGRATE 0.25      /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which vbounds heuristic should at least improve the incumbent          */
#define DEFAULT_MINNODES      500LL     /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_NODESOFS      500LL     /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_MAXPROPROUNDS 2         /* maximum number of propagation rounds during probing */
#define DEFAULT_COPYCUTS      TRUE      /**< should all active cuts from the cutpool of the original scip be copied to
                                         *   constraints of the subscip
                                         */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_VAR**            vbvars;             /**< topological sorted variables with respect to the variable bounds */
   SCIP_BOUNDTYPE*       vbbounds;           /**< topological sorted variables with respect to the variable bounds */
   int                   nvbvars;            /**< number of variables in variable lower bound array */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          usednodes;          /**< nodes already used by vbounds heuristic in earlier calls                         */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Real             minimprove;         /**< factor by which vbounds heuristic should at least improve the incumbent          */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   int                   maxproprounds;      /**< maximum number of propagation rounds during probing */
   SCIP_Bool             initialized;        /**< are the candidate list initialized? */
   SCIP_Bool             applicable;         /**< is the heuristic applicable? */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in
                                              *   subproblem?
                                              */
};

/**@name Propagator defines
 *
 * @{
 *
 * The propagator works on indices representing a bound of a variable. This index will be called bound index in the
 * following. For a given active variable with problem index i (note that active variables have problem indices
 * between 0 and nactivevariable - 1), the bound index of its lower bound is 2*i, the bound index of its upper
 * bound is 2*i + 1. The other way around, a given bound index i corresponds to the variable with problem index
 * i/2 (rounded down), and to the lower bound, if i is even, to the upper bound if i is odd.
 * The following macros can be used to convert bound index into variable problem index and boundtype and vice versa.
 */
#define getLbIndex(idx) (2*(idx))
#define getUbIndex(idx) (2*(idx)+1)
#define getVarIndex(idx) ((idx)/2)
#define getBoundtype(idx) (((idx) % 2 == 0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER)
#define isIndexLowerbound(idx) ((idx) % 2 == 0)


/*
 * Hash map callback methods
 */

/*
 * Local methods
 */

/** reset heuristic data structure */
static
void heurdataReset(
   SCIP_HEURDATA*        heurdata            /**< structure containing heurdata */
   )
{
   heurdata->vbvars = NULL;
   heurdata->vbbounds = NULL;
   heurdata->nvbvars = 0;
   heurdata->initialized = FALSE;
   heurdata->applicable = FALSE;
}


/** performs depth-first-search in the implicitly given directed graph from the given start index */
static
SCIP_RETCODE dfs(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   startnode,          /**< node to start the depth-first-search */
   SCIP_Bool*            visited,            /**< array to store for each node, whether it was already visited */
   int*                  dfsstack,           /**< array of size number of nodes to store the stack;
                                              *   only needed for performance reasons */
   int*                  stacknextedge,      /**< array of size number of nodes to store the number of adjacent nodes
                                              *   already visited for each node on the stack; only needed for
                                              *   performance reasons */
   int*                  dfsnodes,           /**< array of nodes that can be reached starting at startnode, in reverse
                                              *   dfs order */
   int*                  ndfsnodes           /**< pointer to store number of nodes that can be reached starting at
                                              *   startnode */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* startvar;
   SCIP_VAR** vbvars;
   SCIP_Real* coefs;
   SCIP_Bool lower;
   int stacksize;
   int curridx;
   int idx;
   int nvbvars;
   int i;

   assert(startnode >= 0);
   assert(startnode < 2 * SCIPgetNVars(scip));
   assert(visited != NULL);
   assert(visited[startnode] == FALSE);
   assert(dfsstack != NULL);
   assert(dfsnodes != NULL);
   assert(ndfsnodes != NULL);

   vars = SCIPgetVars(scip);

   /* put start node on the stack */
   dfsstack[0] = startnode;
   stacknextedge[0] = 0;
   stacksize = 1;
   idx = -1;

   /* we run until no more bounds indices are on the stack, i.e. all changed bounds were propagated */
   while( stacksize > 0 )
   {
      /* get next node from stack */
      curridx = dfsstack[stacksize - 1];

      /* mark current node as visited */
      assert(visited[curridx] == (stacknextedge[stacksize - 1] > 0));
      visited[curridx] = TRUE;

      startvar = vars[getVarIndex(curridx)];
      lower = isIndexLowerbound(curridx);

      /* go over edges corresponding to varbounds */
      if( lower )
      {
         vbvars = SCIPvarGetVlbVars(startvar);
         coefs = SCIPvarGetVlbCoefs(startvar);
         nvbvars = SCIPvarGetNVlbs(startvar);
      }
      else
      {
         vbvars = SCIPvarGetVubVars(startvar);
         coefs = SCIPvarGetVubCoefs(startvar);
         nvbvars = SCIPvarGetNVubs(startvar);
      }

      /* iterate over all vbounds for the given bound */
      for( i = stacknextedge[stacksize - 1]; i < nvbvars; ++i )
      {
         if( !SCIPvarIsActive(vbvars[i]) )
            continue;

         idx = (SCIPisPositive(scip, coefs[i]) == lower) ? getLbIndex(SCIPvarGetProbindex(vbvars[i])) : getUbIndex(SCIPvarGetProbindex(vbvars[i]));
         assert(idx >= 0);

         /* break when the first unvisited node is reached */
         if( !visited[idx] )
            break;
      }

      /* we stopped because we found an unhandled node and not because we reached the end of the list */
      if( i < nvbvars )
      {
         assert(!visited[idx]);

         /* put the adjacent node onto the stack */
         dfsstack[stacksize] = idx;
         stacknextedge[stacksize] = 0;
         stacknextedge[stacksize - 1] = i + 1;
         stacksize++;
         assert(stacksize <= 2* SCIPgetNVars(scip));

         /* restart while loop, get next index from stack */
         continue;
      }

      /* the current node was completely handled, remove it from stack */
      stacksize--;

      if( (stacksize > 0 || nvbvars > 0) && SCIPvarGetType(startvar) != SCIP_VARTYPE_CONTINUOUS )
      {
         /* store node in the sorted nodes array */
         dfsnodes[(*ndfsnodes)] = curridx;
         (*ndfsnodes)++;
      }
      else
         visited[curridx] = FALSE;
   }

   return SCIP_OKAY;
}


/** sort the bounds of variables topologically */
static
SCIP_RETCODE topologicalSort(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  vbvars,             /**< array to store variable bounds in topological order */
   int*                  nvbvars             /**< array to store number of variable bounds in the graph */
   )
{
   int* dfsstack;
   int* stacknextedge;
   SCIP_Bool* inqueue;
   int nbounds;
   int i;

   assert(scip != NULL);

   nbounds = 2 * SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &dfsstack, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stacknextedge, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inqueue, nbounds) );
   BMSclearMemoryArray(inqueue, nbounds);

   /* while there are unvisited nodes, run dfs starting from one of these nodes; the dfs orders are stored in the
    * topoorder array, later dfs calls are just appended after the stacks of previous dfs calls, which gives us a
    * reverse topological order
    */
   for( i = 0; i < nbounds; ++i )
   {
      if( !inqueue[i] )
      {
         SCIP_CALL( dfs(scip, i, inqueue, dfsstack, stacknextedge, vbvars, nvbvars) );
      }
   }
   assert(*nvbvars <= nbounds);

   SCIPfreeBufferArray(scip, &inqueue);
   SCIPfreeBufferArray(scip, &stacknextedge);
   SCIPfreeBufferArray(scip, &dfsstack);

   return SCIP_OKAY;
}

/** initialize candidate lists */
static
SCIP_RETCODE initializeCandsLists(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< structure containing heurdata */
   )
{
   SCIP_VAR** vars;
   int* vbs;
   int nvars;
   int nvbs;
   int v;

   SCIPdebugMessage("initialize variable bound heuristic (%s)\n", SCIPgetProbName(scip));

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nvbs = 0;

   /* allocate memory for the arrays of the heurdata */
   SCIP_CALL( SCIPallocBufferArray(scip, &vbs, 2 * nvars) );

   /* create the topological sorted variable array with respect to the variable bounds */
   SCIP_CALL( topologicalSort(scip, vbs, &nvbs) );

   /* check if the candidate list contains enough candidates */
   if( nvbs >= heurdata->minfixingrate * nvars )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->vbvars, nvbs) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->vbbounds, nvbs) );

      /* capture variable candidate list */
      for( v = 0; v < nvbs; ++v )
      {
         heurdata->vbvars[v] = vars[getVarIndex(vbs[v])];
         heurdata->vbbounds[v] = getBoundtype(vbs[v]);

         SCIP_CALL( SCIPcaptureVar(scip, heurdata->vbvars[v]) );
      }

      heurdata->nvbvars = nvbs;
      heurdata->applicable = TRUE;
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &vbs);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->initialized = TRUE;

   SCIPstatisticMessage("vbvars %.3g (%s)\n",
      (nvbs * 100.0) / nvars, SCIPgetProbName(scip));

   return SCIP_OKAY;
}

/** apply variable bound fixing during probing */
static
SCIP_RETCODE applyVboundsFixings(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< structure containing heurdata */
   SCIP_VAR**            vars,               /**< variables to fix during probing */
   int                   nvbvars,            /**< number of variables in the variable bound graph */
   SCIP_Bool             forward,            /**< should fixings be done forward w.r.t. the vbound graph? */
   SCIP_Bool             tighten,            /**< should variables be fixed to cause other fixings? */
   SCIP_Bool             obj,                /**< should the objective be taken into account? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether problem is infeasible */
   SCIP_VAR**            lastvar,            /**< last fixed variable */
   SCIP_Bool*            fixedtolb           /**< was last fixed variable fixed to its lower bound? */
   )
{
   SCIP_VAR* var;
   SCIP_RETCODE retcode;
   SCIP_BOUNDTYPE bound;
   SCIP_Bool newnode = TRUE;
   int v;

   /* for each variable in topological order: start at best bound (MINIMIZE: neg coeff --> ub, pos coeff: lb) */
   for( v = 0; v < nvbvars && !(*infeasible); ++v )
   {
      var = forward ? vars[v] : vars[nvbvars - 1 - v];
      bound = forward ? heurdata->vbbounds[v] : heurdata->vbbounds[nvbvars - 1 - v];

      /*SCIPdebugMessage("topoorder[%d]: %s(%s) (%s)\n", v,
         bound == SCIP_BOUNDTYPE_UPPER ? "ub" : "lb", SCIPvarGetName(var),
         SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS ? "c" : "d");*/

      /* only check integer or binary variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      /* skip variables which are already fixed */
      if( SCIPvarGetLbLocal(var) + 0.5 > SCIPvarGetUbLocal(var) )
         continue;

      if( obj && ((SCIPvarGetObj(var) >= 0) == (bound == SCIP_BOUNDTYPE_LOWER)) )
         continue;

      if( newnode )
      {
         retcode = SCIPnewProbingNode(scip);
         if( retcode == SCIP_MAXDEPTHLEVEL )
            newnode = FALSE;
         else
         {
            SCIP_CALL( retcode );
         }
      }

      *lastvar = var;

      if( obj ? (tighten == (SCIPvarGetObj(var) >= 0)) : (tighten == (bound == SCIP_BOUNDTYPE_UPPER)) )
      {
         /* fix variable to lower bound */
         SCIP_CALL( SCIPfixVarProbing(scip, var, SCIPvarGetLbLocal(var)) );
         SCIPdebugMessage("fixing %d: variable <%s> to lower bound <%g> (%d pseudo cands)\n",
            v, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPgetNPseudoBranchCands(scip));
         *fixedtolb = TRUE;
      }
      else
      {
         assert((obj && (tighten == (SCIPvarGetObj(var) < 0)))
            || (!obj && (tighten == (bound == SCIP_BOUNDTYPE_LOWER))));
         /* fix variable to upper bound */
         SCIP_CALL( SCIPfixVarProbing(scip, var, SCIPvarGetUbLocal(var)) );
         SCIPdebugMessage("fixing %d: variable <%s> to upper bound <%g> (%d pseudo cands)\n",
            v, SCIPvarGetName(var), SCIPvarGetUbLocal(var), SCIPgetNPseudoBranchCands(scip));
         *fixedtolb = FALSE;
      }

      /* check if problem is already infeasible */
      SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, infeasible, NULL) );
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

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
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

/** main procedure of the vbounds heuristic */
static
SCIP_RETCODE applyVbounds(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            vbvars,             /**< variables to fix during probing */
   int                   nvbvars,            /**< number of variables to fix */
   SCIP_Bool             forward,            /**< should fixings be done forward w.r.t. the vbound graph? */
   SCIP_Bool             tighten,            /**< should variables be fixed to cause other fixings? */
   SCIP_Bool             obj,                /**< should the objective be taken into account? */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIPstatistic( SCIP_CLOCK* clock; )
   SCIP_VAR** vars;
   SCIP_VAR* lastfixedvar = NULL;
   SCIP_SOL* newsol;
   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */
   SCIP_LPSOLSTAT lpstatus;
   SCIP_Bool infeasible;
   SCIP_Bool lastfixedlower = TRUE;
   SCIP_Bool lperror;
   SCIP_Bool solvelp;
   SCIP_Bool foundsol = FALSE;
   int oldnpscands;
   int npscands;
   int nvars;
   SCIPstatistic( int nprevars = nvars; )

   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(nvbvars > 0);

   /* initialize default values */
   infeasible = FALSE;

   /* get variable data of original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   if( nvbvars < nvars * heurdata->minfixingrate )
      return SCIP_OKAY;

   if( *result == SCIP_DIDNOTRUN )
      *result = SCIP_DIDNOTFIND;

   oldnpscands = SCIPgetNPseudoBranchCands(scip);

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
      SCIPdebugMessage("skipping " HEUR_NAME ": nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPstatistic( SCIP_CALL( SCIPcreateClock(scip, &clock) ) );
   SCIPstatistic( SCIP_CALL( SCIPstartClock(scip, clock) ) );

   SCIPdebugMessage("apply variable bounds heuristic at node %lld on %d variable bounds\n",
      SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), nvbvars);

   /* check whether the LP should be solved at the current node in the tree to determine whether the heuristic
    * is allowed to solve an LP
    */
   solvelp = SCIPhasCurrentNodeLP(scip);

   if( !SCIPisLPConstructed(scip) && solvelp )
   {
      SCIP_Bool nodecutoff;

      SCIP_CALL( SCIPconstructLP(scip, &nodecutoff) );
      SCIP_CALL( SCIPflushLP(scip) );
      if( nodecutoff )
         return SCIP_OKAY;
   }


   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* create temporary solution */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

   /* apply the variable fixings */
   SCIP_CALL( applyVboundsFixings(scip, heurdata, vbvars, nvbvars, forward, tighten, obj, &infeasible,
         &lastfixedvar, &lastfixedlower) );

   /* try to repair probing */
   if( infeasible )
   {
      assert(lastfixedvar != NULL);

      SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );

      /* fix the last variable, which was fixed the reverse bound */
      SCIP_CALL( SCIPfixVarProbing(scip, lastfixedvar,
            lastfixedlower ? SCIPvarGetUbLocal(lastfixedvar) : SCIPvarGetLbLocal(lastfixedvar)) );

      /* propagate fixings */
      SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &infeasible, NULL) );

      SCIPdebugMessage("backtracking ended with %sfeasible problem\n", (infeasible ? "in" : ""));
   }

   /* check that we had enough fixings */
   npscands = SCIPgetNPseudoBranchCands(scip);

   SCIPdebugMessage("npscands=%d, oldnpscands=%d, heurdata->minfixingrate=%g\n", npscands, oldnpscands, heurdata->minfixingrate);

   /* check fixing rate */
   if( npscands > oldnpscands * (1 - heurdata->minfixingrate) )
   {
      SCIPdebugMessage("--> too few fixings\n");

      goto TERMINATE;
   }

   /*************************** Probing LP Solving ***************************/

   lpstatus = SCIP_LPSOLSTAT_ERROR;
   lperror = FALSE;
   /* solve lp only if the problem is still feasible */
   if( !infeasible && solvelp )
   {
      SCIPdebugMessage("starting solving vbound-lp at time %g\n", SCIPgetSolvingTime(scip));

      /* solve LP; errors in the LP solver should not kill the overall solving process, if the LP is just needed for a
       * heuristic.  hence in optimized mode, the return code is caught and a warning is printed, only in debug mode,
       * SCIP will stop.
       */
#ifdef NDEBUG
      {
         SCIP_Bool retstat;
         retstat = SCIPsolveProbingLP(scip, -1, &lperror, NULL);
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while solving LP in vbound heuristic; LP solve terminated with code <%d>\n",
               retstat);
         }
      }
#else
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) );
#endif
      SCIPdebugMessage("ending solving vbound-lp at time %g\n", SCIPgetSolvingTime(scip));

      lpstatus = SCIPgetLPSolstat(scip);

      SCIPdebugMessage(" -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
      SCIPdebugMessage(" -> error=%u, status=%d\n", lperror, lpstatus);
   }

   /* check if this is a feasible solution */
   if( lpstatus == SCIP_LPSOLSTAT_OPTIMAL && !lperror )
   {
      SCIP_Bool stored;
      SCIP_Bool success;

      /* copy the current LP solution to the working solution */
      SCIP_CALL( SCIPlinkLPSol(scip, newsol) );

      SCIP_CALL( SCIProundSol(scip, newsol, &success) );

      if( success )
      {
         SCIPdebugMessage("vbound heuristic found roundable primal solution: obj=%g\n",
            SCIPgetSolOrigObj(scip, newsol));

         /* check solution for feasibility, and add it to solution store if possible.
          * Neither integrality nor feasibility of LP rows have to be checked, because they
          * are guaranteed by the heuristic at this stage.
          */
#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &stored) );
#else
         SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, TRUE, FALSE, FALSE, &stored) );
#endif

         foundsol = TRUE;

         if( stored )
         {
            SCIPdebugMessage("found feasible solution:\n");
            *result = SCIP_FOUNDSOL;
         }
      }
   }
   else
   {
      SCIP_CALL( SCIPclearSol(scip, newsol) );
   }

   /*************************** END Probing LP Solving ***************************/

   /* if no solution has been found --> fix all other variables by subscip if necessary */
   if( !foundsol && lpstatus != SCIP_LPSOLSTAT_INFEASIBLE && lpstatus != SCIP_LPSOLSTAT_OBJLIMIT && !infeasible )
   {
      SCIP* subscip;
      SCIP_VAR** subvars;
      SCIP_HASHMAP* varmap;
      SCIP_Bool valid;
      int i;

      valid = FALSE;

      /* create subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

      SCIP_CALL( SCIPcopy(scip, subscip, varmap, NULL, "_vbounds", FALSE, FALSE, TRUE, &valid) );

      if( heurdata->copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmap, NULL, FALSE, NULL) );
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

      for( i = 0; i < nvars; i++ )
         subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);

      /* free hash map */
      SCIPhashmapFree(&varmap);

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

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

      /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
      if( timelimit <= 0.0 || memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
      {
         /* free subproblem */
         SCIPfreeBufferArray(scip, &subvars);
         SCIP_CALL( SCIPfree(&subscip) );

         goto TERMINATE;
      }

#ifndef SCIP_DEBUG
      /* disable statistic timing inside sub SCIP */
      SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

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
#if 0
      /* use best estimate node selection */
      if( SCIPfindNodesel(subscip, "uct") != NULL && !SCIPisParamFixed(subscip, "nodeselection/uct/stdpriority") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/uct/stdpriority", INT_MAX/4) );
      }
#endif
      /* use inference branching */
      if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
      }

      /* disable conflict analysis */
      if( !SCIPisParamFixed(subscip, "conflict/enable") )
      {
         SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", FALSE) );
      }

      /* employ a limit on the number of enforcement rounds in the quadratic constraint handlers; this fixes the issue that
       * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
       * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
       * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no decutions shall be
       * made for the original SCIP
       */
      if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
      }

#ifdef SCIP_DEBUG
      /* for debugging vbounds heuristic, enable MIP output */
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
       * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
#ifdef NDEBUG
      {
         SCIP_RETCODE retstat;
         retstat = SCIPpresolve(subscip);
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while presolving subMIP in vbounds heuristic; sub-SCIP terminated with code <%d>\n", retstat);
         }
      }
#else
      SCIP_CALL( SCIPpresolve(subscip) );
#endif

      SCIPdebugMessage("vbounds heuristic presolved subproblem: %d vars, %d cons\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip));

      /* after presolving, we should have at least reached a certain fixing rate over ALL variables (including continuous)
       * to ensure that not only the MIP but also the LP relaxation is easy enough
       */
      if( ( nvars - SCIPgetNVars(subscip) ) / (SCIP_Real)nvars >= heurdata->minfixingrate )
      {
         SCIP_SOL** subsols;
         SCIP_Bool success;
         int nsubsols;

         SCIPstatistic( nprevars = SCIPgetNVars(subscip) );

         SCIPdebugMessage("solving subproblem: nstallnodes=%" SCIP_LONGINT_FORMAT ", maxnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->maxnodes);

#ifdef NDEBUG
         {
            SCIP_RETCODE retstat;
            retstat = SCIPsolve(subscip);
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving subMIP in vbounds heuristic; sub-SCIP terminated with code <%d>\n",retstat);
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

 TERMINATE:

#ifdef SCIP_STATISTIC
   SCIP_CALL( SCIPstopClock(scip, clock) );
   SCIPstatisticMessage("vbound: forward=%d tighten=%d obj=%d nvars=%d presolnvars=%d ratio=%.2f infeas=%d found=%d time=%.2f\n",
      forward, tighten, obj, nvars, nprevars, (nvars - nprevars) / (SCIP_Real)nvars, infeasible,
      foundsol ? 1 : 0, SCIPclockGetTime(clock) );
#endif

   SCIPstatistic( SCIP_CALL( SCIPfreeClock(scip, &clock) ) );

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
SCIP_DECL_HEURCOPY(heurCopyVbounds)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of heuristic */
   SCIP_CALL( SCIPincludeHeurVbounds(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeVbounds)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolVbounds)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int v;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* release all variables */
   for( v = 0; v < heurdata->nvbvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &heurdata->vbvars[v]) );
   }

   /* free varbounds array */
   SCIPfreeMemoryArrayNull(scip, &heurdata->vbbounds);
   SCIPfreeMemoryArrayNull(scip, &heurdata->vbvars);

   /* reset heuristic data structure */
   heurdataReset(heurdata);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecVbounds)
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
   {
      SCIP_CALL( initializeCandsLists(scip, heurdata) );
   }

   if( !heurdata->applicable )
      return SCIP_OKAY;

   /* try variable bounds */
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, TRUE, TRUE, TRUE, result) );
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, TRUE, TRUE, FALSE, result) );
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, TRUE, FALSE, FALSE, result) );
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, TRUE, FALSE, TRUE, result) );
#if 0
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, FALSE, TRUE, TRUE, result) );
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, FALSE, TRUE, FALSE, result) );
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, FALSE, FALSE, TRUE, result) );
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, FALSE, FALSE, FALSE, result) );
#endif
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the vbounds primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create vbounds primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdataReset(heurdata);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecVbounds, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyVbounds) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeVbounds) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolVbounds) );

   /* add variable bounds primal heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "minimum percentage of integer variables that have to be fixable",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which " HEUR_NAME " heuristic should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxproprounds",
         "maximum number of propagation rounds during probing (-1 infinity)",
         &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX/4, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
