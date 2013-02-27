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

/**@file   heur_vbounds.c
 * @brief  LNS heuristic uses the variable lower and upper bounds to determine the search neighborhood
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Jens Schulz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_vbounds.h"


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
#define DEFAULT_MINFIXINGRATE 0.5       /* minimum percentage of integer variables that have to be fixed       */
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
   SCIP_VAR**            lbvars;             /**< topological sorted variables with respect to the variable lower bound */
   SCIP_VAR**            ubvars;             /**< topological sorted variables with respect to the variable upper bound */
   SCIP_VAR**            impvars;            /**< topological sorted variables with respect to the variable lower and upper
                                              *   bound and with a corresponding improving objective coefficient */
   int                   nlbvars;            /**< number of variables in variable lower bound array */
   int                   nubvars;            /**< number of variables in variable upper bound array */
   int                   nlbimpvars;         /**< number of variables in variable improving lower bound array */
   int                   nubimpvars;         /**< number of variables in variable improving upper bound array */

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

/*
 * Hash map callback methods
 */

/** hash key retrieval function for variables */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the indices of both variables are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqVar)
{  /*lint --e{715}*/
   if( key1 == key2 )
      return TRUE;
   return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVar)
{  /*lint --e{715}*/
   assert( SCIPvarGetIndex((SCIP_VAR*) key) >= 0 );
   return (unsigned int) SCIPvarGetIndex((SCIP_VAR*) key);
}


/*
 * Local methods
 */

/** reset heuristic data structure */
static
void heurdataReset(
   SCIP_HEURDATA*        heurdata            /**< structure containing heurdata */
   )
{
   heurdata->lbvars = NULL;
   heurdata->ubvars = NULL;
   heurdata->impvars = NULL;
   heurdata->nlbvars = 0;
   heurdata->nubvars = 0;
   heurdata->nlbimpvars = 0;
   heurdata->nubimpvars = 0;
   heurdata->initialized = FALSE;
   heurdata->applicable = FALSE;
}

/** gets the requested variables bounds */
static
void getVariableBounds(
   SCIP_VAR*             var,                /**< variable to get the variable bounds from */
   SCIP_VAR***           vbvars,             /**< pointer to store the variable bound array */
   int*                  nvbvars,            /**< pointer to store the number of variable bounds */
   SCIP_Bool             lowerbound          /**< variable lower bounds? (otherwise variable upper bound) */
   )
{
   if( lowerbound )
   {
      /* get variable lower bounds */
      (*vbvars) = SCIPvarGetVlbVars(var);
      (*nvbvars) = SCIPvarGetNVlbs(var);
   }
   else
   {
      /* get variable upper bounds */
      (*vbvars) = SCIPvarGetVubVars(var);
      (*nvbvars) = SCIPvarGetNVubs(var);
   }
}

/** perform depth-first-search from the given variable using the variable lower or upper bounds of the variable */
static
SCIP_RETCODE depthFirstSearch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to start the depth-first-search  */
   SCIP_HASHMAP*         varPosMap,          /**< mapping a variable to its position in the (used) variable array, or NULL */
   SCIP_VAR**            usedvars,           /**< array of variables which are involved in the propagation, or NULL */
   int*                  nusedvars,          /**< number of variables which are involved in the propagation, or NULL */
   SCIP_HASHTABLE*       connected,          /**< hash table storing if a node was already visited */
   SCIP_VAR**            sortedvars,         /**< array that will contain the topological sorted variables */
   int*                  nsortedvars,        /**< pointer to store the number of already collects variables in the sorted variables array */
   SCIP_Bool             lowerbound          /**< depth-first-search with respect to the variable lower bounds, otherwise variable upper bound */
   )
{
   SCIP_VAR** vbvars;
   SCIP_VAR* vbvar;
   SCIP_Real scalar;
   SCIP_Real constant;
   int nvbvars;
   int v;

   assert(scip != NULL);
   assert(var != NULL);
   assert(varPosMap == NULL || (varPosMap != NULL && usedvars != NULL && nusedvars != NULL));
   assert(sortedvars != NULL);
   assert(nsortedvars != NULL);
   assert(*nsortedvars >= 0);
   assert(SCIPvarGetProbindex(var) > -1);
   assert(SCIPhashtableExists(connected, var));

   /* mark variable as visited, remove variable from hash table */
   SCIP_CALL( SCIPhashtableRemove(connected, var) );

   /* get variable bounds */
   getVariableBounds(var, &vbvars, &nvbvars, lowerbound);

   SCIPdebugMessage("variable <%s> has %d variable %s bounds\n", SCIPvarGetName(var), nvbvars,
      lowerbound ? "lower" : "upper");

   for( v = 0; v < nvbvars; ++v )
   {
      vbvar = vbvars[v];
      assert(vbvar != NULL);

      scalar = 1.0;
      constant = 0.0;

      /* transform variable bound variable to an active variable if possible */
      SCIP_CALL( SCIPgetProbvarSum(scip, &vbvar, &scalar, &constant) );

      /* we could not resolve the variable bound variable to one active variable, therefore, ignore this variable bound */
      if( !SCIPvarIsActive(vbvar) )
         continue;

      /* insert variable bound variable into the hash table since they are involved in later propagation */
      if( varPosMap != NULL && !SCIPhashmapExists(varPosMap, vbvar) )
      {
         SCIPdebugMessage("insert variable <%s> with position %d into the hash map\n", SCIPvarGetName(vbvar), *nusedvars);
         SCIP_CALL( SCIPhashmapInsert(varPosMap, vbvar, (void*)(size_t)(*nusedvars)) );
         usedvars[*nusedvars] = vbvar;
         (*nusedvars)++;
      }

      /* check if the variable bound variable was already visited */
      if( SCIPhashtableExists(connected, vbvar) )
      {
         /* recursively call depth-first-search */
         SCIP_CALL( depthFirstSearch(scip, vbvar, varPosMap, usedvars, nusedvars, connected, sortedvars, nsortedvars, lowerbound) );
      }
   }

   /* store variable in the sorted variable array */
   sortedvars[(*nsortedvars)] = var;
   (*nsortedvars)++;

   /* insert variable bound variable into the hash table since they are involve in the later propagation */
   if( varPosMap != NULL && !SCIPhashmapExists(varPosMap, var) )
   {
      SCIPdebugMessage("insert variable <%s> with position %d into the hash map\n", SCIPvarGetName(var), *nusedvars);
      SCIP_CALL( SCIPhashmapInsert(varPosMap, var, (void*) (size_t)(*nusedvars)) );
      usedvars[*nusedvars] = var;
      (*nusedvars)++;
   }

   return SCIP_OKAY;
}

/** create a topological sorted variable array of the given variables and stores if (needed) the involved variables into
 *  the corresponding variable array and hash map
 *
 * @note: for all arrays and the hash map (if requested) you need to allocate enough memory before calling this method
 */
static
SCIP_RETCODE createTopoSortedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable which we want sort */
   int                   nvars,              /**< number of variables */
   SCIP_HASHMAP*         varPosMap,          /**< mapping a variable to its position in the (used) variable array, or NULL */
   SCIP_VAR**            usedvars,           /**< array of variables which are involved in the propagation, or NULL */
   int*                  nusedvars,          /**< number of variables which are involved in the propagation, or NULL */
   SCIP_VAR**            topovars,           /**< array where the topological sorted variables are stored */
   int*                  ntopovars,          /**< pointer to store the number of topological sorted variables */
   SCIP_Bool             lowerbound          /**< topological sorted with respect to the variable lower bounds, otherwise variable upper bound */
   )
{
   SCIP_VAR** sortedvars;
   SCIP_VAR** vbvars;
   SCIP_VAR* var;
   SCIP_HASHTABLE* connected;
   int nvbvars;
   int hashsize;
   int i;
   int v;

   assert(scip != NULL);
   assert(vars != NULL || nvars == 0);
   assert(varPosMap == NULL || (varPosMap != NULL && usedvars != NULL && nusedvars != NULL));
   assert(topovars != NULL);
   assert(ntopovars != NULL);

   SCIPdebugMessage("create topological sorted variable array with respect to variables %s bounds\n",
      lowerbound ? "lower" : "upper");

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedvars, nvars) );

   hashsize = SCIPcalcHashtableSize(5 * nvars);

   /* create hash table for variables which are (still) connected */
   SCIP_CALL( SCIPhashtableCreate(&connected, SCIPblkmem(scip), hashsize, SCIPvarGetHashkey, SCIPvarIsHashkeyEq, SCIPvarGetHashkeyVal, NULL) );

   /* detect isolated variables; mark all variables which have at least one entering or leaving arc as connected */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);

      if( !SCIPvarIsActive(var) )
         continue;

      /* get variable bounds */
      getVariableBounds(var, &vbvars, &nvbvars, lowerbound);

      if( nvbvars > 0 && !SCIPhashtableExists(connected, var) )
      {
         SCIP_CALL( SCIPhashtableInsert(connected, var) );
      }

      for( i = 0; i < nvbvars; ++i )
      {
         if( !SCIPvarIsActive(vbvars[i]) )
            continue;

         /* there is a leaving arc, hence, the variable/node is connected */
         assert(vbvars[i] != NULL);
         if( !SCIPhashtableExists(connected, vbvars[i]) )
         {
            SCIP_CALL( SCIPhashtableInsert(connected, vbvars[i]) );
         }
      }
   }

   /* loop over all "connected" variable and find for each connected component a "almost" topological sorted version */
   for( v = 0; v < nvars; ++v )
   {
      if( SCIPhashtableExists(connected, vars[v]) )
      {
         int nsortedvars;

         SCIPdebugMessage("start depth-first-search with variable <%s>\n", SCIPvarGetName(vars[v]));

         /* use depth first search to get a "almost" topological sorted variables for the connected component which
          * includes vars[v]
          */
         nsortedvars = 0;
         SCIP_CALL( depthFirstSearch(scip, vars[v], varPosMap, usedvars, nusedvars, connected, sortedvars, &nsortedvars, lowerbound) );

         SCIPdebugMessage("detected connected component of size <%d>\n", nsortedvars);

         /* copy variables */
         for( i = 0; i < nsortedvars; ++i )
         {
            topovars[(*ntopovars)] = sortedvars[i];
            (*ntopovars)++;
         }
      }
   }

   assert(*ntopovars <= nvars);
   SCIPdebugMessage("topological sorted array contains %d of %d variables (variable %s bound)\n",
      *ntopovars, nvars, lowerbound ? "lower" : "upper");

   /* free hash table */
   SCIPhashtableFree(&connected);

   /* free buffer memory */
   SCIPfreeBufferArray(scip, &sortedvars);

   return SCIP_OKAY;
}

/** initialize candidate lists */
static
SCIP_RETCODE initializeCandsLists(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< structure containing heurdata */
   )
{
   SCIP_VAR** allvars;
   SCIP_VAR** vars;
   SCIP_VAR** lbvars;
   SCIP_VAR** ubvars;
   SCIP_VAR** impvars;
   SCIP_HASHTABLE* collectedvars;
   int nallvars;
   int nvars;
   int nlbvars;
   int nubvars;
   int nlbimpvars;
   int nubimpvars;
   int v;

   SCIPdebugMessage("initialize variable bound heuristic (%s)\n", SCIPgetProbName(scip));

   allvars = SCIPgetVars(scip);
   nallvars = SCIPgetNVars(scip);
   nvars = 0;
   nlbvars = 0;
   nubvars = 0;
   nlbimpvars = 0;
   nubimpvars = 0;

   /* allocate memory for the arrays of the heurdata */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nallvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lbvars, nallvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubvars, nallvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &impvars, nallvars) );

   /* create the topological sorted variable array with respect to the variable lower bounds */
   SCIP_CALL( createTopoSortedVars(scip, allvars, nallvars, NULL, vars, &nvars, lbvars, &nlbvars, TRUE) );

   /* create the topological sorted variable array with respect to the variable upper bounds */
   SCIP_CALL( createTopoSortedVars(scip, allvars, nallvars, NULL, vars, &nvars, ubvars, &nubvars, FALSE) );

   /* create hash table for variables which are already collected */
   SCIP_CALL( SCIPhashtableCreate(&collectedvars, SCIPblkmem(scip), SCIPcalcHashtableSize(nallvars), hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );

   /* collect variables which improve the objective by fixing them to suggested bound  */
   for( v = 0; v < nlbvars; ++v )
   {
      if( SCIPisGE(scip, SCIPvarGetObj(lbvars[v]), 0.0) )
      {
         SCIP_CALL( SCIPhashtableInsert(collectedvars, lbvars[v]) );
         assert(nlbimpvars < nallvars);
         impvars[nlbimpvars] = lbvars[v];
         nlbimpvars++;
      }
   }

   for( v = 0; v < nubvars; ++v )
   {
      if( SCIPisLE(scip, SCIPvarGetObj(ubvars[v]), 0.0) && !SCIPhashtableExists(collectedvars, ubvars[v])  )
      {
         assert(nlbimpvars + nubimpvars < nallvars);
         impvars[nlbimpvars + nubimpvars] = ubvars[v];
         nubimpvars++;
      }
   }

   /* free hash table */
   SCIPhashtableFree(&collectedvars);

   /* check if the candidate lists contain enough candidates */
   if( nlbvars >= heurdata->minfixingrate * nallvars )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->lbvars, lbvars, nlbvars) );
      heurdata->nlbvars = nlbvars;
      heurdata->applicable = TRUE;

      /* capture variable candidate list */
      for( v = 0; v < nlbvars; ++v )
      {
         SCIP_CALL( SCIPcaptureVar(scip, heurdata->lbvars[v]) );
      }
   }
   if( nubvars >= heurdata->minfixingrate * nallvars )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->ubvars, ubvars, nubvars) );
      heurdata->nubvars = nubvars;
      heurdata->applicable = TRUE;

      /* capture variable candidate list */
      for( v = 0; v < nubvars; ++v )
      {
         SCIP_CALL( SCIPcaptureVar(scip, heurdata->ubvars[v]) );
      }
   }
   if( nlbvars > nlbimpvars && nubvars > nubimpvars && nlbimpvars + nubimpvars >= heurdata->minfixingrate * nallvars )
   {
      assert(nlbimpvars < INT_MAX - nubimpvars);
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &heurdata->impvars, impvars, nlbimpvars + nubimpvars) );
      heurdata->nlbimpvars = nlbimpvars;
      heurdata->nubimpvars = nubimpvars;
      heurdata->applicable = TRUE;

      /* capture variable candidate list */
      for( v = 0; v < nlbimpvars + nubimpvars; ++v )
      {
         SCIP_CALL( SCIPcaptureVar(scip, heurdata->impvars[v]) );
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &impvars);
   SCIPfreeBufferArray(scip, &ubvars);
   SCIPfreeBufferArray(scip, &lbvars);
   SCIPfreeBufferArray(scip, &vars);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->initialized = TRUE;

   SCIPstatisticMessage("lbvars %.3g, ubvars %.3g, impvars %.3g (%s)\n",
      (nlbvars * 100.0) / nallvars, (nubvars * 100.0) / nallvars,
      ((nlbimpvars + nubimpvars) * 100.0) / nallvars, SCIPgetProbName(scip));

   return SCIP_OKAY;
}

/** apply variable bound fixing during probing */
static
SCIP_RETCODE applyVboundsFixings(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< structure containing heurdata */
   SCIP_VAR**            vars,               /**< variables to fix during probing */
   int                   nlbvars,            /**< number of variables to use the lower bound */
   int                   nubvars,            /**< number of variables to use the upper bound */
   SCIP_SOL*             sol,                /**< working solution */
   SCIP_Bool*            infeasible,         /**< pointer to store whether problem is infeasible */
   SCIP_RESULT*          result              /**< pointer to store the result (solution found) */
   )
{
   SCIP_Bool success;
   int v;

   /* for each variable in topological order: start at best bound (MINIMIZE: neg coeff --> ub, pos coeff: lb) */
   for( v = 0; v < nlbvars + nubvars && !(*infeasible); ++v )
   {
      /* only check integer or binary variables */
      if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      /* skip variables which are already fixed */
      if( SCIPvarGetLbLocal(vars[v]) + 0.5 > SCIPvarGetUbLocal(vars[v]) )
         continue;

      if( v < nlbvars )
      {
         /* fix variable to lower bound */
         SCIP_CALL( SCIPfixVarProbing(scip, vars[v], SCIPvarGetLbLocal(vars[v])) );
         SCIPdebugMessage("fixing %d: variable <%s> to lower bound <%g> (%d pseudo cands)\n",
            v, SCIPvarGetName(vars[v]), SCIPvarGetLbLocal(vars[v]), SCIPgetNPseudoBranchCands(scip));
      }
      else
      {
         /* fix variable to lower bound */
         SCIP_CALL( SCIPfixVarProbing(scip, vars[v], SCIPvarGetUbLocal(vars[v])) );
         SCIPdebugMessage("fixing %d: variable <%s> to upper bound <%g> (%d pseudo cands)\n",
            v, SCIPvarGetName(vars[v]), SCIPvarGetUbLocal(vars[v]), SCIPgetNPseudoBranchCands(scip));
      }

      /* check if problem is already infeasible */
      SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, infeasible, NULL) );

      if( !(*infeasible) )
      {
         /* create solution from probing run and try to round it */
         SCIP_CALL( SCIPlinkCurrentSol(scip, sol) );
         SCIP_CALL( SCIProundSol(scip, sol, &success) );

         if( success )
         {
            SCIPdebugMessage("vbound heuristic found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, sol));

            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, TRUE, &success) );

            /* check, if solution was feasible and good enough */
            if( success )
            {
               SCIPdebugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;

               if( SCIPisStopped(scip) )
                  break;
            }
         }
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
   SCIP_VAR**            probvars,           /**< variables to fix during probing */
   int                   nlbvars,            /**< number of variables to use the lower bound */
   int                   nubvars,            /**< number of variables to use the upper bound */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_VAR** vars;
   SCIP_SOL* newsol;
   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */
   SCIP_Bool infeasible;
   int nvars;

   assert(heur != NULL);
   assert(heurdata != NULL);

   /* initialize default values */
   infeasible = FALSE;

   /* get variable data of original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   if( nlbvars + nubvars < nvars * heurdata->minfixingrate )
      return SCIP_OKAY;

   assert(nlbvars + nubvars > 0);

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

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMessage("apply variable bounds heuristic at node %lld on %d variable lower bound and %d variable upper bounds\n",
      SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), nlbvars, nubvars);

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* create temporary solution */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

   /* apply the variable fixings */
   SCIP_CALL( applyVboundsFixings(scip, heurdata, probvars, nlbvars, nubvars, newsol, &infeasible, result) );

   /* if a solution has been found --> fix all other variables by subscip if necessary */
   if( !infeasible )
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
         /** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
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
   for( v = 0; v < heurdata->nlbvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &heurdata->lbvars[v]) );
   }
   /* release all variables */
   for( v = 0; v < heurdata->nubvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &heurdata->ubvars[v]) );
   }
   /* release all variables */
   for( v = 0; v < heurdata->nlbimpvars + heurdata->nubimpvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &heurdata->impvars[v]) );
   }

   /* free varbounds array */
   SCIPfreeMemoryArrayNull(scip, &heurdata->impvars);
   SCIPfreeMemoryArrayNull(scip, &heurdata->lbvars);
   SCIPfreeMemoryArrayNull(scip, &heurdata->ubvars);

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

   *result = SCIP_DIDNOTFIND;

   /* try variable lower and upper bounds which respect to objective coefficients */
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->impvars, heurdata->nlbimpvars, heurdata->nubimpvars, result) );

   /* try variable lower bounds */
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->lbvars, heurdata->nlbvars, 0, result) );

   /* try variable upper bounds */
   SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->ubvars, 0, heurdata->nubvars, result) );

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
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
         "minimum percentage of integer variables that have to be fixable",
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
         "factor by which "HEUR_NAME" heuristic should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxproprounds",
         "maximum number of propagation rounds during probing (-1 infinity)",
         &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX/4, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
