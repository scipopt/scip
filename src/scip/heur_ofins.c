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

/**@file   heur_ofins.c
 * @brief  OFINS - Objective Function Induced Neighborhood Search. Primal heuristic for reoptimization
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "scip/heur_ofins.h"
#include "scip/scipdefplugins.h"       /* needed for the secondary SCIP instance */
#include "scip/pub_misc.h"

#define HEUR_NAME             "ofins"
#define HEUR_DESC             "primal heuristic for reoptimization, objective function induced neighborhood search"
#define HEUR_DISPCHAR         'A'
#define HEUR_PRIORITY         60000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

/* default values for OFINS-specific plugins */
#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem */
#define DEFAULT_MAXCHGRATE    0.50      /* maximum percentage of changed objective coefficients */
#define DEFAULT_COPYCUTS      TRUE      /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip
                                         */
#define DEFAULT_MAXCHANGE     0.04     /* maximal rate of change per coefficient to get fixed */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which OFINS should at least improve the incumbent          */
#define DEFAULT_ADDALLSOLS   FALSE      /* should all subproblem solutions be added to the original SCIP? */
#define DEFAULT_MINNODES      50LL      /* minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS      500LL     /* number of nodes added to the contingent of the total nodes */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_LPLIMFAC      2.0       /* factor by which the limit on the number of LP depends on the node limit */

/* event handler properties */
#define EVENTHDLR_NAME         "Ofins"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Real             maxchangerate;           /**< maximal rate of changed coefficients in the objective function */
   SCIP_Longint          maxnodes;                /**< maximum number of nodes to regard in the subproblem */
   SCIP_Bool             copycuts;                /**< should all active cuts from cutpool be copied to constraints in subproblem? */
   SCIP_Bool             addallsols;              /**< should all subproblem solutions be added to the original SCIP? */
   SCIP_Longint          minnodes;                /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;                /**< number of nodes added to the contingent of the total nodes */
   SCIP_Real             maxchange;               /**< maximal rate of change per coefficient to get fixed */
   SCIP_Real             minimprove;              /**< factor by which OFINS should at least improve the incumbent */
   SCIP_Real             nodesquot;               /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             nodelimit;               /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;                /**< factor by which the limit on the number of LP depends on the node limit */
};

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecOfins)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP */
   if( SCIPgetNLPs(scip) > heurdata->lplimfac * heurdata->nodelimit )
   {
      SCIPdebugMessage("interrupt after %" SCIP_LONGINT_FORMAT " LPs\n",SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

/** creates a subproblem by fixing a number of variables */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_Bool*            chgcoeffs           /**< array indicating which coefficients have changed */
   )
{
   SCIP_VAR** vars;
   SCIP_SOL* sol;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);

   /* get binary variables; all continuous variables remain unfixed */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNBinVars(scip);

   /* get optimal solution of the last iteration */
   sol = SCIPgetReoptLastOptSol(scip);
   assert(sol != NULL);

   /* change bounds of variables of the subproblem */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real solval;

      if( !chgcoeffs[i] )
      {
         assert(subvars[i] != NULL);

         solval = SCIPgetSolVal(scip, sol, vars[i]);

         /* perform the bound change */
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], solval) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], solval) );
      }
   }

   /* set an objective limit */
   SCIPdebugMessage("set objective limit of %g to sub-SCIP\n", SCIPgetUpperbound(scip));
   SCIP_CALL( SCIPsetObjlimit(subscip, SCIPgetUpperbound(scip)) );

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< RENS heuristic structure                            */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;                         /* the original problem's number of variables      */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
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

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** main procedure of the OFINS heuristic, creates and solves a sub-SCIP */
static
SCIP_RETCODE applyOfins(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_HEURDATA*        heurdata,           /**< euristic's private data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem */
   SCIP_Bool*            chgcoeffs           /**< array of changed coefficients */
   )
{
   SCIP* subscip;                                 /* the subproblem created by OFINS */
   SCIP_HASHMAP* varmapfw;                        /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_VAR** vars;                               /* source problem's variables */
   SCIP_VAR** subvars;                            /* subproblem's variables */
   SCIP_EVENTHDLR* eventhdlr;                     /* event handler for LP events */

   SCIP_Real timelimit;                           /* time limit for OFINS subproblem */
   SCIP_Real memorylimit;                         /* memory limit for OFINS subproblem */

   int nvars;                                     /* number of source problem's variables */
   int i;

   SCIP_SOL** subsols;
   int nsubsols;

   SCIP_Bool valid;
   SCIP_Bool success;
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);
   assert(chgcoeffs != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("+---+ Start OFINS heuristic +---+\n");

   /* check whether there is enough time and memory left */
   timelimit = 0.0;
   memorylimit = 0.0;
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
      SCIPdebugMessage("-> not enough memory left\n");
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* get variable data */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   eventhdlr = NULL;

   valid = FALSE;

   /* copy complete SCIP instance */
   SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "ofins", TRUE, FALSE, TRUE, &valid) );

   if( heurdata->copycuts )
   {
      /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
      SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, TRUE, NULL) );
   }

   SCIPdebugMessage("Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");

   /* create event handler for LP events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecOfins, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   for( i = 0; i < nvars; i++ )
   {
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);
     assert(subvars[i] != NULL);
   }
   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
   SCIP_CALL( createSubproblem(scip, subscip, subvars, chgcoeffs) );
   SCIPdebugMessage("OFINS subproblem: %d vars, %d cons\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip));

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* set limits for the subproblem */
   heurdata->nodelimit = heurdata->maxnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
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
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", FALSE) );
   }

#ifdef SCIP_DEBUG
   /* for debugging RENS, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#endif

   /* presolve the subproblem */
   retcode = SCIPpresolve(subscip);

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while presolving subproblem in %s heuristic; sub-SCIP terminated with code <%d>\n", HEUR_NAME, retcode);

      /* free */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
      return SCIP_OKAY;
   }

   SCIPdebugMessage("%s presolved subproblem: %d vars, %d cons\n", HEUR_NAME, SCIPgetNVars(subscip), SCIPgetNConss(subscip));

   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

   /* solve the subproblem */
   SCIPdebugMessage("solving subproblem: nstallnodes=%" SCIP_LONGINT_FORMAT ", maxnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->maxnodes);
   retcode = SCIPsolve(subscip);

   SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in RENS heuristic; sub-SCIP terminated with code <%d>\n", retcode);
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);
   success = FALSE;
   for( i = 0; i < nsubsols && (!success || heurdata->addallsols); i++ )
   {
      SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      if( success )
         *result = SCIP_FOUNDSOL;
   }

   SCIPstatisticPrintf("%s statistic: fixed %6.3f integer variables, needed %6.1f seconds, %" SCIP_LONGINT_FORMAT " nodes, solution %10.4f found at node %" SCIP_LONGINT_FORMAT "\n",
      HEUR_NAME, 0.0, SCIPgetSolvingTime(subscip), SCIPgetNNodes(subscip), success ? SCIPgetPrimalbound(scip) : SCIPinfinity(scip),
      nsubsols > 0 ? SCIPsolGetNodenum(SCIPgetBestSol(subscip)) : -1 );

   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/* put your local methods here, and declare them static */


/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeOfins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecOfins)
{/*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   SCIP_Bool* chgcoeffs;
   SCIP_Longint nstallnodes;
   int nchgcoefs;
   int nvars;
   int v;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* only call heuristic, if reoptimization is enabled */
   if( !SCIPisReoptEnabled(scip) )
      return SCIP_OKAY;

   /* only call the heuristic, if we are in run >= 2 */
   if( SCIPgetNReoptRuns(scip) <= 1 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward OFINS if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-SCIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMessage("skipping OFINS: nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* get variable data and check which coefficient has changed  */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNBinVars(scip);
   nchgcoefs = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &chgcoeffs, nvars) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_VAR* origvar;
      SCIP_Real constant;
      SCIP_Real scalar;
      SCIP_Real newcoef;
      SCIP_Real oldcoef;
      SCIP_Real newcoefabs;
      SCIP_Real oldcoefabs;
      SCIP_Real frac;

      /* we only want to count variables that are unfixed after the presolving */
      assert(SCIPvarGetStatus(vars[v]) != SCIP_VARSTATUS_ORIGINAL);
      assert(SCIPvarIsActive(vars[v]));

      origvar = vars[v];
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      newcoef = SCIPvarGetObj(origvar);
      SCIP_CALL( SCIPgetReoptOldObjCoef(scip, origvar, SCIPgetNReoptRuns(scip)-1, &oldcoef) );
      newcoefabs = fabs(newcoef);
      oldcoefabs = fabs(oldcoef);

      frac = SCIP_INVALID;

      /* if both coefficients are zero nothing has changed */
      if( SCIPisZero(scip, newcoef) && SCIPisZero(scip, oldcoef) )
      {
         frac = 0;
      }
      /* if exactly one coefficient is zero, the other need to be close to zero */
      else if( SCIPisZero(scip, newcoef) || SCIPisZero(scip, oldcoef) )
      {
         assert(SCIPisZero(scip, newcoef) != SCIPisZero(scip, oldcoef));
         if( !SCIPisZero(scip, newcoef) )
            frac = MIN(1, newcoefabs);
         else
            frac = MIN(1, oldcoefabs);
      }
      /* if both coefficients have the same sign we calculate the quotient
       * MIN(newcoefabs, oldcoefabs)/MAX(newcoefabs, oldcoefabs)
       */
      else if( SCIPisPositive(scip, newcoef) == SCIPisPositive(scip, oldcoef) )
      {
         frac = MIN(newcoefabs, oldcoefabs)/MAX(newcoefabs, oldcoefabs);
      }
      /* if both coefficients have a different sign, we set frac = 1 */
      else
      {
         assert((SCIPisPositive(scip, newcoef) && SCIPisNegative(scip, oldcoef))
             || (SCIPisNegative(scip, newcoef) && SCIPisPositive(scip, oldcoef)));

         frac = 1;
      }

      if( frac > heurdata->maxchange )
      {
         chgcoeffs[v] = TRUE;
         nchgcoefs++;
      }
      else
         chgcoeffs[v] = FALSE;
   }

   SCIPdebugMessage("%d (rate %.4f) changed coefficients\n", nchgcoefs, nchgcoefs/((SCIP_Real)nvars));

   /* we only want to run the heuristic, if there at least 3 changed coefficients.
    * if the number of changed coefficients is 2 the trivialnegation heuristic will construct an
    * optimal solution without solving a MIP.
    */
   if( nchgcoefs < 3 )
      goto TERMINATE;

   /* run the heuristic, if not too many coefficients have changed */
   if( nchgcoefs/((SCIP_Real)SCIPgetNBinVars(scip)) > heurdata->maxchangerate )
      goto TERMINATE;

   SCIP_CALL( applyOfins(scip, heur, heurdata, result, nstallnodes, chgcoeffs) );

  TERMINATE:
   SCIPfreeBufferArray(scip, &chgcoeffs);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the ofins primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurOfins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create ofins primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   assert(heurdata != NULL);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecOfins, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeOfins) );

   /* add ofins primal heuristic parameters */

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxchangerate",
         "maximal rate of changed coefficients",
         &heurdata->maxchangerate, FALSE, DEFAULT_MAXCHGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxchange",
         "maximal rate of change per coefficient to get fixed",
         &heurdata->maxchange, FALSE, DEFAULT_MAXCHANGE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/addallsols",
         "should all subproblem solutions be added to the original SCIP?",
         &heurdata->addallsols, TRUE, DEFAULT_ADDALLSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which RENS should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
