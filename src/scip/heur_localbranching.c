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

/**@file   heur_localbranching.c
 * @brief  Local branching heuristic according to Fischetti and Lodi
 * @author Timo Berthold
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_localbranching.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "localbranching"
#define HEUR_DESC             "local branching heuristic by Fischetti and Lodi"
#define HEUR_DISPCHAR         'L'
#define HEUR_PRIORITY         -1102000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NEIGHBORHOODSIZE  18    /* radius of the incumbents neighborhood to be searched                     */
#define DEFAULT_NODESOFS      1000      /* number of nodes added to the contingent of the total nodes               */
#define DEFAULT_MAXNODES      10000     /* maximum number of nodes to regard in the subproblem                      */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which localbranching should at least improve the incumbent     */
#define DEFAULT_MINNODES      1000      /* minimum number of nodes required to start the subproblem                 */
#define DEFAULT_NODESQUOT     0.05      /* contingent of sub problem nodes in relation to original nodes            */
#define DEFAULT_NWAITINGNODES 200       /* number of nodes without incumbent change that heuristic should wait      */
#define DEFAULT_USELPROWS    FALSE      /* should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructors of the constraints handlers are used    */
#define DEFAULT_COPYCUTS      TRUE      /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip
                                         */


#define EXECUTE               0
#define WAITFORNEWSOL         1


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait  */
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes           */
   int                   minnodes;           /**< minimum number of nodes required to start the subproblem             */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem                  */
   SCIP_Longint          usednodes;          /**< amount of nodes local branching used during all calls                */
   SCIP_Real             nodesquot;          /**< contingent of sub problem nodes in relation to original nodes        */
   SCIP_Real             minimprove;         /**< factor by which localbranching should at least improve the incumbent */
   int                   neighborhoodsize;   /**< radius of the incumbent's neighborhood to be searched                */
   int                   callstatus;         /**< current status of localbranching heuristic                           */
   SCIP_SOL*             lastsol;            /**< the last incumbent localbranching used as reference point            */
   int                   curneighborhoodsize;/**< current neighborhoodsize                                             */
   int                   curminnodes;        /**< current minimal number of nodes required to start the subproblem     */
   int                   emptyneighborhoodsize;/**< size of neighborhood that was proven to be empty                   */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?         */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?
                                              */
};


/*
 * Local methods
 */

/** copies the problem of scip to the problem of subscip - only necessary if uselprows is false */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< SCIP data structure of the original problem                      */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem                            */
   SCIP_VAR**            subvars             /**< variables of the subproblem                                      */
   )
{
   SCIP_ROW** rows;
   int nrows;
   int i;

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   for( i = 0; i < nrows; i++ )
   {
      SCIP_CONS* cons;
      SCIP_VAR** consvars;
      SCIP_COL** cols;
      SCIP_Real constant;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real* vals;
      int nnonz;
      int j;

      /* ignore rows that are only locally valid */
      if( SCIProwIsLocal(rows[i]) )
         continue;

      /* get the row's data */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIProwGetRhs(rows[i]) - constant;
      vals = SCIProwGetVals(rows[i]);
      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);

      assert(lhs <= rhs);

      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ )
         consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

      /* create new constraint and add it to subscip */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

      /* free memory */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** create the extra constraint of local branching and add it to subscip */
static
SCIP_RETCODE addLocalBranchingConstraint(
   SCIP*                 scip,               /**< SCIP data structure of the original problem   */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem         */
   SCIP_VAR**            subvars,            /**< variables of the subproblem                   */
   SCIP_HEURDATA*        heurdata            /**< heuristic's data structure                    */
   )
{
   SCIP_CONS* cons;                        /* local branching constraint to create */
   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;

   int nbinvars;
   int i;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real* consvals;
   char consname[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_localbranchcons", SCIPgetProbName(scip));

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, NULL, NULL, NULL) );
   bestsol = SCIPgetBestSol(scip);
   assert( bestsol != NULL );

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nbinvars) );

   /* set initial left and right hand sides of local branching constraint */
   lhs = (SCIP_Real)heurdata->emptyneighborhoodsize + 1.0;
   rhs = (SCIP_Real)heurdata->curneighborhoodsize;

   /* create the distance (to incumbent) function of the binary variables */
   for( i = 0; i < nbinvars; i++ )
   {
      SCIP_Real solval;

      solval = SCIPgetSolVal(scip, bestsol, vars[i]);
      assert( SCIPisFeasIntegral(scip,solval) );

      /* is variable i  part of the binary support of bestsol? */
      if( SCIPisFeasEQ(scip,solval,1.0) )
      {
         consvals[i] = -1.0;
         rhs -= 1.0;
         lhs -= 1.0;
      }
      else
         consvals[i] = 1.0;
      consvars[i] = subvars[i];
      assert( SCIPvarGetType(consvars[i]) == SCIP_VARTYPE_BINARY );
   }

   /* creates localbranching constraint and adds it to subscip */
   SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, consname, nbinvars, consvars, consvals,
         lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

   /* free local memory */
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< SCIP data structure  of the original problem      */
   SCIP*                 subscip,            /**< SCIP data structure  of the subproblem            */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< the Localbranching heuristic                      */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< pointer to store, whether new solution was found  */
   )
{
   SCIP_VAR** vars;
   int nvars;
   SCIP_SOL* newsol;
   SCIP_Real* subsolvals;

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( subvars != NULL );
   assert( subsol != NULL );

   /* copy the solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyLocalbranching)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurLocalbranching(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLocalbranching)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitLocalbranching)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* with a little abuse we initialize the heurdata as if localbranching would have finished its last step regularly */
   heurdata->callstatus = WAITFORNEWSOL;
   heurdata->lastsol = NULL;
   heurdata->usednodes = 0;
   heurdata->curneighborhoodsize = heurdata->neighborhoodsize;
   heurdata->curminnodes = heurdata->minnodes;
   heurdata->emptyneighborhoodsize = 0;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLocalbranching)
{  /*lint --e{715}*/
   SCIP_Longint maxnnodes;                   /* maximum number of subnodes                            */
   SCIP_Longint nsubnodes;                   /* nodelimit for subscip                                 */

   SCIP_HEURDATA* heurdata;
   SCIP* subscip;                            /* the subproblem created by localbranching              */
   SCIP_VAR** subvars;                       /* subproblem's variables                                */
   SCIP_SOL* bestsol;                        /* best solution so far                                  */

   SCIP_Real timelimit;                      /* timelimit for subscip (equals remaining time of scip) */
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem                   */
   SCIP_Real upperbound;
   SCIP_Real memorylimit;

   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_VAR** vars;

   int nvars;
   int i;

   SCIP_Bool success;

   SCIP_RETCODE retcode;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* there should be enough binary variables that a local branching constraint makes sense */
   if( SCIPgetNBinVars(scip) < 2*heurdata->neighborhoodsize )
      return SCIP_OKAY;

   *result = SCIP_DELAYED;

   /* only call heuristic, if an IP solution is at hand */
   if( SCIPgetNSols(scip) <= 0  )
         return SCIP_OKAY;

   /* only call heuristic, if the best solution comes from transformed problem */
   assert( SCIPgetBestSol(scip) != NULL );
   if( SCIPsolIsOriginal(SCIPgetBestSol(scip)) )
      return SCIP_OKAY;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip,SCIPgetBestSol(scip))  < heurdata->nwaitingnodes)
      return SCIP_OKAY;

   /* reset neighborhood and minnodes, if new solution was found */
   bestsol = SCIPgetBestSol(scip);
   if( heurdata->lastsol != bestsol )
   {
      heurdata->curneighborhoodsize = heurdata->neighborhoodsize;
      heurdata->curminnodes = heurdata->minnodes;
      heurdata->emptyneighborhoodsize = 0;
      heurdata->callstatus = EXECUTE;
      heurdata->lastsol = bestsol;
   }

   /* if no new solution was found and local branching also seems to fail, just keep on waiting */
   if( heurdata->callstatus == WAITFORNEWSOL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   maxnnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward local branching if it succeeded often */
   maxnnodes = (SCIP_Longint)(maxnnodes * (1.0 + 2.0*(SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur)+1.0)));
   maxnnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   maxnnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nsubnodes = maxnnodes - heurdata->usednodes;
   nsubnodes = MIN(nsubnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call sub problem solving */
   if( nsubnodes < heurdata->curminnodes )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("running localbranching heuristic ...\n");

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initializing the subproblem */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );
   success = FALSE;

   if( heurdata->uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

      /* get name of the original problem and add the string "_localbranchsub" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_localbranchsub", SCIPgetProbName(scip));

         /* create the subproblem */
      SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );
   }
   else
   {
      SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "localbranchsub", TRUE, FALSE, TRUE, &success) );

      if( heurdata->copycuts )
      {
         /** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, TRUE, NULL) );
      }
   }
   SCIPdebugMessage("Copying the plugins was %ssuccessful.\n", success ? "" : "not ");

   for (i = 0; i < nvars; ++i)
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* if the subproblem could not be created, free memory and return */
   if( !success )
   {
      *result = SCIP_DIDNOTRUN;
      goto TERMINATE;
   }

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifndef SCIP_DEBUG
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif

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
      goto TERMINATE;

   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nsubnodes) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", 1) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
   if( !SCIPisParamFixed(subscip, "conflict/useprop") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useinflp") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useinflp", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useboundlp", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/usesb") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/usepseudo") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) );
   }

   /* copy the original problem and add the local branching constraint */
   if( heurdata->uselprows )
   {
      SCIP_CALL( createSubproblem(scip, subscip, subvars) );
   }
   SCIP_CALL( addLocalBranchingConstraint(scip, subscip, subvars, heurdata) );

   /* add an objective cutoff */
   cutoff = SCIPinfinity(scip);
   assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

   upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
   if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
   {
      cutoff = (1-heurdata->minimprove)*SCIPgetUpperbound(scip) + heurdata->minimprove*SCIPgetLowerbound(scip);
   }
   else
   {
      if( SCIPgetUpperbound ( scip ) >= 0 )
         cutoff = ( 1 - heurdata->minimprove ) * SCIPgetUpperbound ( scip );
      else
         cutoff = ( 1 + heurdata->minimprove ) * SCIPgetUpperbound ( scip );
   }
   cutoff = MIN(upperbound, cutoff );
   SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );

   /* solve the subproblem */
   SCIPdebugMessage("solving local branching subproblem with neighborhoodsize %d and maxnodes %"SCIP_LONGINT_FORMAT"\n",
      heurdata->curneighborhoodsize, nsubnodes);
   retcode = SCIPsolve(subscip);

   /* Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in local branching heuristic; sub-SCIP terminated with code <%d>\n",retcode);
   }

   heurdata->usednodes += SCIPgetNNodes(subscip);
   SCIPdebugMessage("local branching used %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT" nodes\n",
      SCIPgetNNodes(subscip), nsubnodes);

   /* check, whether a solution was found */
   if( SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      }
      if( success )
      {
         SCIPdebugMessage("-> accepted solution of value %g\n", SCIPgetSolOrigObj(subscip, subsols[i]));
         *result = SCIP_FOUNDSOL;
      }
   }

   /* check the status of the sub-MIP */
   switch( SCIPgetStatus(subscip) )
   {
   case SCIP_STATUS_OPTIMAL:
   case SCIP_STATUS_BESTSOLLIMIT:
      heurdata->callstatus = WAITFORNEWSOL; /* new solution will immediately be installed at next call */
      SCIPdebugMessage(" -> found new solution\n");
      break;

   case SCIP_STATUS_NODELIMIT:
   case SCIP_STATUS_STALLNODELIMIT:
   case SCIP_STATUS_TOTALNODELIMIT:
      heurdata->callstatus = EXECUTE;
      heurdata->curneighborhoodsize = (heurdata->emptyneighborhoodsize + heurdata->curneighborhoodsize)/2;
      heurdata->curminnodes *= 2;
      SCIPdebugMessage(" -> node limit reached: reduced neighborhood to %d, increased minnodes to %d\n",
         heurdata->curneighborhoodsize, heurdata->curminnodes);
      if( heurdata->curneighborhoodsize <= heurdata->emptyneighborhoodsize )
      {
         heurdata->callstatus = WAITFORNEWSOL;
         SCIPdebugMessage(" -> new neighborhood was already proven to be empty: wait for new solution\n");
      }
      break;

   case SCIP_STATUS_INFEASIBLE:
   case SCIP_STATUS_INFORUNBD:
      heurdata->emptyneighborhoodsize = heurdata->curneighborhoodsize;
      heurdata->curneighborhoodsize += heurdata->curneighborhoodsize/2;
      heurdata->curneighborhoodsize = MAX(heurdata->curneighborhoodsize, heurdata->emptyneighborhoodsize + 2);
      heurdata->callstatus = EXECUTE;
      SCIPdebugMessage(" -> neighborhood is empty: increased neighborhood to %d\n", heurdata->curneighborhoodsize);
      break;

   case SCIP_STATUS_UNKNOWN:
   case SCIP_STATUS_USERINTERRUPT:
   case SCIP_STATUS_TIMELIMIT:
   case SCIP_STATUS_MEMLIMIT:
   case SCIP_STATUS_GAPLIMIT:
   case SCIP_STATUS_SOLLIMIT:
   case SCIP_STATUS_UNBOUNDED:
   default:
      heurdata->callstatus = WAITFORNEWSOL;
      SCIPdebugMessage(" -> unexpected sub-MIP status <%d>: waiting for new solution\n", SCIPgetStatus(subscip));
      break;
   }

 TERMINATE:
   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the localbranching primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurLocalbranching(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   SCIP_HEUR* heur;

   /* create Localbranching primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLocalbranching, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLocalbranching) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLocalbranching) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLocalbranching) );

   /* add localbranching primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/neighborhoodsize",
         "radius (using Manhattan metric) of the incumbent's neighborhood to be searched",
         &heurdata->neighborhoodsize, FALSE, DEFAULT_NEIGHBORHOODSIZE, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which localbranching should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
