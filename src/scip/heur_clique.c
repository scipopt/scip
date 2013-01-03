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

/**@file   heur_clique.c
 * @brief  LNS heuristic using a clique partition to restrict the search neighborhood
 * @brief  clique primal heuristic
 * @author Stefan Heinz
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/heur_clique.h"
#include "scip/cons_logicor.h"


#define HEUR_NAME             "clique"
#define HEUR_DESC             "LNS heuristic using a clique partition to restrict the search neighborhood"
#define HEUR_DISPCHAR         'Q'
#define HEUR_PRIORITY         -1000500
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE                       /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      5000LL                     /* maximum number of nodes to regard in the subproblem */
#define DEFAULT_MINFIXINGRATE 0.5                        /* minimum percentage of integer variables that have to be fixed */
#define DEFAULT_MINIMPROVE    0.01                       /* factor by which clique heuristic should at least improve the
                                                          * incumbent
                                                          */
#define DEFAULT_MINNODES      500LL                      /* minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS      500LL                      /* number of nodes added to the contingent of the total nodes */
#define DEFAULT_NODESQUOT     0.1                        /* subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_MAXPROPROUNDS 2                          /* maximum number of propagation rounds during probing */
#define DEFAULT_INITSEED      0                          /**< random seed value to initialize the random permutation
                                                          * value for variables
                                                          */
#define DEFAULT_MULTIPLIER    1.1                        /**< value to increase node number to determine the next run */
#define DEFAULT_COPYCUTS      TRUE                       /**< should all active cuts from the cutpool of the
                                                          *   original scip be copied to constraints of the subscip
                                                          */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Longint          usednodes;          /**< nodes already used by clique heuristic in earlier calls */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed */
   SCIP_Real             minimprove;         /**< factor by which clique heuristic should at least improve the incumbent */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   int                   maxproprounds;      /**< maximum number of propagation rounds during probing */
   SCIP_Longint          nnodefornextrun;    /**< node number for next run */
   SCIP_Real             multiplier;         /**< multiplier to determine next node number */
   int                   initseed;           /**< initial random seed value */
   unsigned int          seed;               /**< seed value for random number generator */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in
                                              *   subproblem?
                                              */
};

/*
 * Local methods
 */

/** comparison method for sorting variables by non-decreasing index */
static
SCIP_DECL_SORTPTRCOMP(varObjSort)
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   var1 = (SCIP_VAR*)elem1;
   var2 = (SCIP_VAR*)elem2;

   if( SCIPvarGetObj(var1) < SCIPvarGetObj(var2) )
      return -1;
   else if( SCIPvarGetObj(var1) > SCIPvarGetObj(var2) )
      return +1;
   else
      return 0;
}

/** sort the binary variable array w.r.t. the clique partition; thereby ensure the current order within the cliques are
 *  not changed
 */
static
SCIP_RETCODE stableSortBinvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            binvars,            /**< array of binary variables to sort */
   int                   nbinvars,           /**< number of binary variables */
   int*                  cliquepartition,    /**< clique partition to use */
   int                   ncliques            /**< number of cliques */
   )
{
   SCIP_VAR*** varpointers;
   SCIP_VAR** vars;
   int* cliquecount;
   int nextpos;
   int c;
   int v;
   int cliquenumber;

   assert(scip != NULL);
   assert(binvars != NULL);
   assert(cliquepartition != NULL);

   /* @note: we don't want to loose order from same clique numbers, so we need a stable sorting algorithm, or we first
    *       count all clique items and alloc temporary memory for a bucket sort */
   /* sort variables after clique-numbers */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquecount, ncliques) );
   BMSclearMemoryArray(cliquecount, ncliques);

   /* first we count for each clique the number of elements */
   for( v = nbinvars - 1; v >= 0; --v )
   {
      assert(0 <= cliquepartition[v] && cliquepartition[v] < ncliques);
      ++(cliquecount[cliquepartition[v]]);
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbinvars) );
#ifndef NDEBUG
   BMSclearMemoryArray(vars, nbinvars);
#endif
   SCIP_CALL( SCIPallocBufferArray(scip, &varpointers, ncliques) );

   nextpos = 0;
   /* now we initialize all start pointers for each clique, so they will be ordered */
   for( c = 0; c < ncliques; ++c )
   {
      /* to reach the goal that all variables of each clique will be standing next to each other we will initialize the
       * starting pointers for each clique by adding the number of each clique to the last clique starting pointer
       * e.g. clique1 has 4 elements and clique2 has 3 elements the the starting pointer for clique1 will be the pointer
       *      to vars[0], the starting pointer to clique2 will be the pointer to vars[4] and to clique3 it will be
       *      vars[7]
       *
       */
      varpointers[c] = (SCIP_VAR**) (vars + nextpos);
      assert(cliquecount[c] > 0);
      nextpos += cliquecount[c];
      assert(nextpos > 0);
   }
   assert(nextpos == nbinvars);

   /* now we copy all variable to the right order in our temporary variable array */
   for( v = 0; v < nbinvars; ++v )
   {
      *(varpointers[cliquepartition[v]]) = binvars[v];
      ++(varpointers[cliquepartition[v]]);
   }
#ifndef NDEBUG
   for( v = 0; v < nbinvars; ++v )
      assert(vars[v] != NULL);
#endif

   /* move all variables back to our variable array */
   BMScopyMemoryArray(binvars, vars, nbinvars);

   cliquenumber = 0;
   nextpos = cliquecount[0];

   c = 1;
   for( v = 0; v < nbinvars; ++v )
   {
      if( v == nextpos )
      {
         nextpos += cliquecount[c];
         ++c;
         ++cliquenumber;
      }
      cliquepartition[v] = cliquenumber;
   }
   assert(cliquepartition[v - 1] == ncliques - 1);

#ifndef NDEBUG
   for( v = 1; v < nbinvars; ++v )
      assert(SCIPvarGetObj(binvars[v - 1]) <= SCIPvarGetObj(binvars[v - 1]));
#endif

#if 0
#ifndef NDEBUG
   c = -1;
   for( v = 0; v < nbinvars; ++v )
   {
      if( cliquepartition[v] > c )
      {
         ++c;
#if 0
         if( c == 1 )
            break;
#endif
         assert(cliquepartition[v] == c);
         printf("Clique %d (%d elements): ",c, cliquecount[c]);
      }
#if 0
      printf("%s ", SCIPvarGetName(binvars[v]));
#endif
      if( v < nbinvars - 1 && cliquepartition[v + 1] > c )
      {
         printf("\n");
      }
   }
   printf("\n");
#endif
#endif

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &varpointers);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &cliquecount);

   return SCIP_OKAY;
}

/** apply clique fixing using probing */
static
SCIP_RETCODE applyCliqueFixings(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< structure containing heurdata */
   SCIP_VAR**            binvars,            /**< binary variables order w.r.t. to clique partition */
   int                   nbinvars,           /**< number of binary variables */
   int*                  cliquepartition,    /**< clique partition of all binary variables */
   int                   ncliques,           /**< number of cliques */
   SCIP_VAR**            onefixvars,         /**< array to store all variables which are stored to one */
   int*                  nonefixvars,        /**< pointer to store the number of variables fixed to one */
   SCIP_SOL*             sol,                /**< working solution */
   int*                  probingdepthofonefix,/**< pointer to store in which depth the last fixing to was applied  */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the propagation stopped with infeasibility */
   SCIP_RESULT*          result              /**< pointer to store the result (solution found) */
   )
{
#if 0
   SCIP_Bool success;
#endif
   SCIP_Bool alreadyone;
   SCIP_Bool allfixed;
   int bestpos;
   int v;
   int c;
#ifdef SCIP_DEBUG
   int nsolsround;
   int nsolstried;
#endif

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(binvars != NULL);
   assert(onefixvars != NULL);
   assert(nonefixvars != NULL);
   assert(sol != NULL);
   assert(probingdepthofonefix != NULL);
   assert(cutoff != NULL);
   assert(result != NULL);

   *cutoff = FALSE;
   *probingdepthofonefix = 0;

#ifdef SCIP_DEBUG
   nsolsround = 0;
   nsolstried = 0;
#endif
   v = 0;
   /* @todo maybe try to fix more than one variable to one in each probing node, to gain faster results */
   for( c = 0; c < ncliques; ++c )
   {
      alreadyone = FALSE;
      allfixed = TRUE;
      bestpos = nbinvars;

      /* find first unfixed variable in this clique */
      while( v < nbinvars && cliquepartition[v] == c )
      {
         if( SCIPvarGetLbLocal(binvars[v]) > 0.5 )
            alreadyone = TRUE;
         else if( allfixed && SCIPvarGetUbLocal(binvars[v]) > 0.5 )
         {
            bestpos = v;
            allfixed = FALSE;
         }

         ++v;
      }
      if( v == nbinvars && allfixed )
         break;

      /* if all clique variables are fixed, continue with the next clique */
      if( allfixed )
         continue;

      assert(bestpos < nbinvars);
      assert(c == cliquepartition[bestpos]);

      SCIP_CALL( SCIPnewProbingNode(scip) );

      v = bestpos;
      if( !alreadyone )
      {
         *probingdepthofonefix = SCIPgetProbingDepth(scip);
         onefixvars[(*nonefixvars)] = binvars[v];
         ++(*nonefixvars);

         /* fix best possible clique variable to 1 */
         SCIP_CALL( SCIPfixVarProbing(scip, binvars[v], 1.0) );
         SCIPdebugMessage("probing: fixing variable <%s> to 1\n", SCIPvarGetName(binvars[v]));

         ++v;
      }

      /* fix rest of unfixed clique variables to 0 */
      while( v < nbinvars && cliquepartition[v] == c )
      {
         if( SCIPvarGetUbLocal(binvars[v]) > 0.5 && SCIPvarGetLbLocal(binvars[v]) < 0.5 )
         {
            SCIP_CALL( SCIPfixVarProbing(scip, binvars[v], 0.0) );
            SCIPdebugMessage("probing: fixing variable <%s> to 0\n", SCIPvarGetName(binvars[v]));
         }
         ++v;
      }

      /* propagate fixings */
      SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, cutoff, NULL) );

      if( *cutoff )
         break;

      /* @todo need to be check if it's ok to always try to round and check the solution in each probing step */
#if 0

#ifdef SCIP_DEBUG
      ++nsolsround;
#endif
      /* create solution from probing run and try to round it */
      SCIP_CALL( SCIPlinkCurrentSol(scip, sol) );
      SCIP_CALL( SCIProundSol(scip, sol, &success) );

      if( success )
      {
         SCIPdebugMessage("clique heuristic found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, sol));

#ifdef SCIP_DEBUG
         ++nsolstried;
#endif
         /* try to add solution to SCIP */
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, TRUE, &success) );

         /* check, if solution was feasible and good enough */
         if( success )
         {
            SCIPdebugMessage(" -> solution was feasible and good enough\n");
            *result = SCIP_FOUNDSOL;
         }
      }
#endif

      if( SCIPisStopped(scip) )
         return SCIP_OKAY;

#if 0
      /* if the rest of all variables are in cliques with one variable stop */
      if( nbinvars - v == ncliques - c )
         break;
#endif
   }
   assert((*probingdepthofonefix > 0 && *nonefixvars > 0) || (*probingdepthofonefix == 0 && *nonefixvars == 0));
   assert(*cutoff || (nbinvars - v == ncliques - c) || (v == nbinvars && (c == ncliques || c == ncliques - 1)));

   SCIPdebugMessage("fixed %d of %d variables in probing\n", v, nbinvars);
   SCIPdebugMessage("applied %d of %d cliques in probing\n", c, ncliques);
   SCIPdebugMessage("probing was %sfeasible\n", (*cutoff) ? "in" : "");
#ifdef SCIP_DEBUG
   SCIPdebugMessage("clique heuristic rounded %d solutions and tried %d of them\n", nsolsround, nsolstried);
#endif
   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_SOL*             newsol,             /**< working solution */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   int nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);
   assert(success != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyClique)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurClique(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeClique)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitClique)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* reset heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* set the seed value to the initial random seed value */
   heurdata->seed = (unsigned int) heurdata->initseed;

   heurdata->usednodes = 0;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecClique)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   int nvars;
   SCIP_VAR** binvars;
   int nbinvars;
   int* cliquepartition;
   int ncliques;
#if 0
   SCIP_Longint tmpnnodes;
#endif
   SCIP_Bool cutoff;
   SCIP_Bool backtrackcutoff;
   SCIP_Bool lperror;

   int probingdepthofonefix;
   SCIP_VAR** onefixvars;
   int nonefixvars;
   SCIP_Bool enabledconflicts;
   SCIP_LPSOLSTAT lpstatus;
   SCIP_CONS* conflictcons;
   SCIP_Bool shortconflict;
   SCIP_Bool allfixsolfound;
   SCIP_Bool backtracked;
   char consname[SCIP_MAXSTRLEN];

   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */

   SCIP_SOL* sol;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   /* assert(SCIPhasCurrentNodeLP(scip)); */

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

#if 0
   if( heurdata->nnodefornextrun != SCIPgetNNodes(scip) )
      return SCIP_OKAY;
#endif
   /* get all binary variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );

   if( nbinvars < 2 )
   {
      heurdata->nnodefornextrun = INT_MAX;
      return SCIP_OKAY;
   }

   /* check for necessary information to apply this heuristic */
   if( SCIPgetNCliques(scip) == 0 && SCIPgetNImplications(scip) == 0 )
   {
      heurdata->nnodefornextrun = INT_MAX;
      return SCIP_OKAY;
   }

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

   *result = SCIP_DIDNOTFIND;

   onefixvars = NULL;
   sol = NULL;

   /* allocate memory */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &binvars, vars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquepartition, nbinvars) );

#if 1
   /* @todo change sorting after some attempts to random variable order */
   if( SCIPgetNNodes(scip) == 1 )
   {
      /* sort variables after increasing objective value */
      SCIPsortPtr((void**)binvars, varObjSort, nbinvars);
   }
   else
   {
      SCIPpermuteArray((void**)binvars, 0, nbinvars, &(heurdata->seed));
   }
#endif

   /* get clique partitions */
   SCIP_CALL( SCIPcalcCliquePartition(scip, binvars, nbinvars, cliquepartition, &ncliques) );
   /* @todo get negated clique partition and use this too, or maybe mix both */

   SCIPdebugMessage("found %d cliques\n", ncliques);

   /* disable conflict analysis, because we can it better than SCIP itself, cause we have more information */
   SCIP_CALL( SCIPgetBoolParam(scip, "conflict/enable", &enabledconflicts) );
   if( !SCIPisParamFixed(scip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "conflict/enable", FALSE) );
   }

   if( ncliques == nbinvars )
   {
      heurdata->nnodefornextrun = INT_MAX;
      goto TERMINATE;
   }

   /* sort the cliques together by respecting the current order (which is w.r.t. the objective coefficients */
   SCIP_CALL( stableSortBinvars(scip, binvars, nbinvars, cliquepartition, ncliques) );

   if( !SCIPisLPConstructed(scip) && SCIPhasCurrentNodeLP(scip) )
   {
      SCIP_Bool nodecutoff;

      SCIP_CALL( SCIPconstructLP(scip, &nodecutoff) );
      SCIP_CALL( SCIPflushLP(scip) );
      if( nodecutoff )
         goto TERMINATE;
   }

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* create a solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   /* allocate memory for all variables which will be fixed to one during probing */
   SCIP_CALL(SCIPallocBufferArray(scip, &onefixvars, ncliques) );
   nonefixvars = 0;

   /* apply fixings due to clique information */
   SCIP_CALL( applyCliqueFixings(scip, heurdata, binvars, nbinvars, cliquepartition, ncliques, onefixvars, &nonefixvars, sol, &probingdepthofonefix, &cutoff, result) );

   if( SCIPisStopped(scip) )
      goto TERMINATE;

   backtrackcutoff = FALSE;
   backtracked = FALSE;

   /* try to repair probing */
   if( cutoff && nonefixvars > 0)
   {
      assert(probingdepthofonefix > 0);

      SCIP_CALL( SCIPbacktrackProbing(scip, probingdepthofonefix - 1) );

      /* fix the last variable, which was fixed to 1 and led to the cutoff, to 0 */
      SCIP_CALL( SCIPfixVarProbing(scip, onefixvars[nonefixvars - 1], 0.0) );

      backtracked = TRUE;

      /* propagate fixings */
      SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &backtrackcutoff, NULL) );
   }

   /*************************** Probing LP Solving ***************************/

   lpstatus = SCIP_LPSOLSTAT_ERROR;
   lperror = FALSE;
   allfixsolfound = FALSE;
   /* solve lp only if the problem is still feasible */
   if( !backtrackcutoff && SCIPhasCurrentNodeLP(scip) )
   {
#if 1
      SCIPdebugMessage("starting solving clique-lp at time %g\n", SCIPgetSolvingTime(scip));

      /* solve LP; errors in the LP solver should not kill the overall solving process, if the LP is just needed for a
       * heuristic.  hence in optimized mode, the return code is caught and a warning is printed, only in debug mode,
       * SCIP will stop.
       */
#ifdef NDEBUG
      {
         SCIP_Bool retstat;
         retstat = SCIPsolveProbingLP(scip, -1, &lperror);
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while solving LP in clique heuristic; LP solve terminated with code <%d>\n",
               retstat);
         }
      }
#else
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror) );
#endif
      SCIPdebugMessage("ending solving clique-lp at time %g\n", SCIPgetSolvingTime(scip));

      lpstatus = SCIPgetLPSolstat(scip);

      SCIPdebugMessage(" -> new LP iterations: %"SCIP_LONGINT_FORMAT"\n", SCIPgetNLPIterations(scip));
      SCIPdebugMessage(" -> error=%u, status=%d\n", lperror, lpstatus);
   }

   /* check if this is a feasible solution */
   if( lpstatus == SCIP_LPSOLSTAT_OPTIMAL && !lperror )
   {
      SCIP_Bool stored;
      SCIP_Bool success;

      /* copy the current LP solution to the working solution */
      SCIP_CALL( SCIPlinkLPSol(scip, sol) );

      SCIP_CALL( SCIProundSol(scip, sol, &success) );

      if( success )
      {
         SCIPdebugMessage("clique heuristic found roundable primal solution: obj=%g\n",
            SCIPgetSolOrigObj(scip, sol));

         /* check solution for feasibility, and add it to solution store if possible.
          * Neither integrality nor feasibility of LP rows have to be checked, because they
          * are guaranteed by the heuristic at this stage.
          */
#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPtrySol(scip, sol, TRUE, TRUE, TRUE, TRUE, &stored) );
#else
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, TRUE, FALSE, FALSE, &stored) );
#endif
         if( stored )
         {
            SCIPdebugMessage("found feasible solution:\n");
            SCIPdebug( SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) ) );
            *result = SCIP_FOUNDSOL;
            allfixsolfound = TRUE;
         }
      }
   }
#endif

   /*************************** END Probing LP Solving ***************************/
   /*************************** Create Conflict ***************************/

   if( lpstatus == SCIP_LPSOLSTAT_INFEASIBLE || lpstatus == SCIP_LPSOLSTAT_OBJLIMIT || backtrackcutoff )
   {
      /* in case the last fixing in both direction led to infeasibility or to a reached objlimit than our conflict will
       * only include all variable before that last fixing
       */
      shortconflict = cutoff && (nonefixvars > 0);

      /* create own conflict */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "conf%d", SCIPgetNNodes(scip));

      /* get negated variables for our conflict */
      SCIP_CALL( SCIPgetNegatedVars(scip, nonefixvars, onefixvars, onefixvars) );

      /* create conflict constraint */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &conflictcons, consname, (shortconflict ? nonefixvars - 1 : nonefixvars), onefixvars,
            FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, SCIPgetCurrentNode(scip), conflictcons, NULL) );
      SCIPdebugPrintCons(scip, conflictcons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &conflictcons) );
   }

   /*************************** End Conflict ***************************/

   /*************************** Start Subscip Solving ***************************/

   /* if no solution has been found yet and the subproblem is still feasible --> fix all other variables by subscip if
    * necessary
    */
   if( !allfixsolfound && lpstatus != SCIP_LPSOLSTAT_INFEASIBLE && lpstatus != SCIP_LPSOLSTAT_OBJLIMIT && !backtrackcutoff )
   {
      SCIP* subscip;
      SCIP_VAR** subvars;
      SCIP_HASHMAP* varmap;
      SCIP_Bool valid;
      int i;

      valid = FALSE;

      /* create subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      /* allocate temporary memory for subscip variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

      SCIP_CALL( SCIPcopy(scip, subscip, varmap, NULL, "_clique", FALSE, FALSE, TRUE, &valid) );

      if( heurdata->copycuts )
      {
         /** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmap, NULL, FALSE, NULL) );
      }

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
      /* for debugging clique heuristic, enable MIP output */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
      SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#endif

      /* if there is already a solution, add an objective cutoff */
      if( SCIPgetNSols(scip) > 0 )
      {
         SCIP_Real upperbound;
         SCIP_Real minimprove;
         SCIP_Real cutoffbound;

         minimprove = heurdata->minimprove;
         cutoffbound = SCIPinfinity(scip);
         assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

         upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

         if( !SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
         {
            cutoffbound = (1-minimprove) * SCIPgetUpperbound(scip) + minimprove * SCIPgetLowerbound(scip);
         }
         else
         {
            if( SCIPgetUpperbound ( scip ) >= 0 )
               cutoffbound = (1 - minimprove) * SCIPgetUpperbound(scip);
            else
               cutoffbound = (1 + minimprove) * SCIPgetUpperbound(scip);
         }
         cutoffbound = MIN(upperbound, cutoffbound);
         SCIP_CALL( SCIPsetObjlimit(subscip, cutoffbound) );
         SCIPdebugMessage("setting objlimit for subscip to %g\n", cutoffbound);
      }

      SCIPdebugMessage("starting solving clique-submip at time %g\n", SCIPgetSolvingTime(scip));

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
            SCIPwarningMessage(scip, "Error while presolving subMIP in clique heuristic; sub-SCIP terminated with code <%d>\n", retstat);
         }
      }
#else
      SCIP_CALL( SCIPpresolve(subscip) );
#endif

      SCIPdebugMessage("clique heuristic presolved subproblem: %d vars, %d cons; fixing value = %g\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip), ((nvars - SCIPgetNVars(subscip)) / (SCIP_Real)nvars));

      /* after presolving, we should have at least reached a certain fixing rate over ALL variables (including continuous)
       * to ensure that not only the MIP but also the LP relaxation is easy enough
       */
      if( ((nvars - SCIPgetNVars(subscip)) / (SCIP_Real)nvars) >= (heurdata->minfixingrate / 2.0) )
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
               SCIPwarningMessage(scip, "Error while solving subMIP in clique heuristic; sub-SCIP terminated with code <%d>\n",retstat);
            }
         }
#else
         SCIP_CALL( SCIPsolve(subscip) );
#endif
         SCIPdebugMessage("ending solving clique-submip at time %g, status = %d\n", SCIPgetSolvingTime(scip), SCIPgetStatus(subscip));

         /* check, whether a solution was found; due to numerics, it might happen that not all solutions are feasible ->
          * try all solutions until one was accepted
          */
         nsubsols = SCIPgetNSols(subscip);
         subsols = SCIPgetSols(subscip);
         success = FALSE;

         for( i = 0; i < nsubsols && !success; ++i )
         {
            SCIP_CALL( createNewSol(scip, subscip, subvars, sol, subsols[i], &success) );
         }
         if( success )
            *result = SCIP_FOUNDSOL;

         /* if subscip was infeasible we can add a conflict too */
         if( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE )
         {
            /* in case the last fixing in both direction led to infeasibility or to a reached objlimit than our conflict will only include all variable before that last fixing */
            shortconflict = backtracked;

            /* create own conflict */
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "conf%d", SCIPgetNNodes(scip));

            /* get negated variables for our conflict */
            SCIP_CALL( SCIPgetNegatedVars(scip, nonefixvars, onefixvars, onefixvars) );

            /* create conflict constraint */
            SCIP_CALL( SCIPcreateConsLogicor(scip, &conflictcons, consname, (shortconflict ? nonefixvars - 1 : nonefixvars), onefixvars,
                  FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddConsNode(scip, SCIPgetCurrentNode(scip), conflictcons, NULL) );
            SCIPdebugPrintCons(scip, conflictcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &conflictcons) );
         }

      }

#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

      /* free subproblem */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
   }

   /*************************** End Subscip Solving ***************************/

 TERMINATE:

   /* reset the conflict analysis */
   if( !SCIPisParamFixed(scip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "conflict/enable", enabledconflicts) );
   }

   /* free conflict variables */
   if( onefixvars != NULL )
      SCIPfreeBufferArray(scip, &onefixvars);

   /* freeing solution */
   if( sol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &sol) );
   }

   /* end probing */
   if( SCIPinProbing(scip) )
   {
      SCIP_CALL( SCIPendProbing(scip) );
   }

   SCIPfreeBufferArray(scip, &cliquepartition);
   SCIPfreeBufferArray(scip, &binvars);

#if 0
   /* calculate next node number to run this heuristic */
   tmpnnodes = (SCIP_Longint) SCIPceil(scip, heurdata->nnodefornextrun * heurdata->multiplier);
   heurdata->nnodefornextrun = MIN(tmpnnodes, INT_MAX);
   SCIPdebugMessage("Next run will be at node %lld.\n", heurdata->nnodefornextrun);
#endif
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the clique primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurClique(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create clique primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecClique, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyClique) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeClique) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitClique) );

   /* add clique primal heuristic parameters */

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/multiplier",
         "value to increase nodenumber to determine the next run",
         &heurdata->multiplier, TRUE, DEFAULT_MULTIPLIER, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/initseed",
         "initial random seed value to permutate variables",
         &(heurdata->initseed), TRUE, DEFAULT_INITSEED, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
         "minimum percentage of integer variables that have to be fixable",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

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
