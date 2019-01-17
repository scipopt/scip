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

/**@file   heur_trustregion.c
 * @brief  Large neighbourhood search heuristic for Benders' decomposition based on trust region methods
 * @author Stephen J. Maher
 *
 * The Trust Region heuristic draws upon trust region methods for solving optimisation problems, especially in the
 * context of Benders' decomposition. This heuristic has been developed to improve the heuristic performance of the
 * Benders' decomposition algorithm within SCIP.
 *
 * The Trust Region heuristic copies the original SCIP instance and adds a constraint to penalise changes from the
 * incumbent solution. Consider a problem that includes a set of binary variables \f$\mathcal{B}\f$. Given a feasible
 * solution \f$\hat{x}\f$ to the original problem, we define the set \f$\mathcal{B}^{+}\f$ as the index set for the
 * binary variables that are 1 in the input solution and \f$\mathcal{B}^{-}\f$ as the index set for binary variables
 * that are 0. The trust region constraint, which is added to the sub-SCIP, is given by
 *
 * \f[
 *    \sum_{i \in \mathcal{B}^{+}}(1 - x_{i}) + \sum_{i \in \mathcal{B}^{-}}x_{i} \le \theta
 * \f]
 *
 * The variable \f$\theta\f$ measure the distance, in terms of the binary variables, of candidate solutions to the input
 * solution.
 *
 * In addition, an upper bounding constraint is explicitly added to enforce a minimum improvement from the heuristic,
 * given by \f$f(x) \le f(\hat{x}) - \epsilon\f$. The parameter \f$\epsilon \ge 0\f$ denotes the minimum improvement
 * that must be achieved by the heuristic.
 *
 * The objective function is then modified to \f$f(x) + M\theta\f$, where \f$M\f$ is a parameter for penalising the
 * distance of solutions from the input solution \f$\hat{x}\f$.
 *
 * If a new incumbent solution is found by this heuristic, then the Trust Region heuristic is immediately
 * re-executed with this new incumbent solution.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/heuristics.h"
#include "scip/heur_trustregion.h"
#include "scip/pub_event.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"
#include <string.h>

#define HEUR_NAME             "trustregion"
#define HEUR_DESC             "LNS heuristic for Benders' decomposition based on trust region methods"
#define HEUR_DISPCHAR         '?'
#define HEUR_PRIORITY         -1102000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MINBINVARS    10        /* the minimum number of binary variables necessary to run the heuristic */
#define DEFAULT_NODESOFS      1000      /* number of nodes added to the contingent of the total nodes */
#define DEFAULT_MAXNODES      10000     /* maximum number of nodes to regard in the subproblem */
#define DEFAULT_MINNODES      100       /* minimum number of nodes required to start the subproblem */
#define DEFAULT_NODESQUOT     0.05      /* contingent of sub problem nodes in relation to original nodes */
#define DEFAULT_LPLIMFAC      1.5       /* factor by which the limit on the number of LP depends on the node limit */
#define DEFAULT_NWAITINGNODES 1         /* number of nodes without incumbent change that heuristic should wait */
#define DEFAULT_USELPROWS     FALSE     /* should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE      /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip */
#define DEFAULT_BESTSOLLIMIT   3         /* limit on number of improving incumbent solutions in sub-CIP */
#define DEFAULT_USEUCT        FALSE     /* should uct node selection be used at the beginning of the search? */

#define DEFAULT_VIOLPENALTY   100       /* the penalty for violating the trust region */
#define DEFAULT_OBJMINIMPROVE 1e-2      /* the minimum improvement in the objective function value */
#define DEFAULT_NWAITCALLS    3         /* the number of calls to wait before resolving with the same solution */

/* event handler properties */
#define EVENTHDLR_NAME         "Trustregion"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"


#define EXECUTE               0
#define WAITFORNEWSOL         1


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   int                   minnodes;           /**< minimum number of nodes required to start the subproblem */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Longint          usednodes;          /**< amount of nodes trust region used during all calls */
   SCIP_Real             nodesquot;          /**< contingent of sub problem nodes in relation to original nodes */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   int                   minbinvars;         /**< minimum number of binary variables necessary to run the heuristic */
   int                   callstatus;         /**< current status of trustregion heuristic */
   SCIP_SOL*             lastsol;            /**< the last incumbent trustregion used as reference point */
   int                   curminnodes;        /**< current minimal number of nodes required to start the subproblem */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows? */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem? */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP */
   SCIP_Bool             useuct;             /**< should uct node selection be used at the beginning of the search? */
   SCIP_Real             violpenalty;        /**< the penalty for violating the trust region */
   SCIP_Real             objminimprove;      /**< the minimum improvement in the objective function value */
};


/*
 * Local methods
 */

/** create the extra constraint of trust region and add it to subscip */
static
SCIP_RETCODE addTrustRegionConstraints(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem */
   SCIP_VAR**            subvars,            /**< variables of the subproblem */
   SCIP_HEURDATA*        heurdata            /**< heuristic's data structure */
   )
{
   SCIP_VAR* violvar;
   SCIP_CONS* cons;                        /* trust region constraint to create */
   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;

   int nvars;
   int nbinvars;
   int i;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real* consvals;
   char name[SCIP_MAXSTRLEN];

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   bestsol = SCIPgetBestSol(scip);
   assert( bestsol != NULL );

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars + 1) );

   /* set initial left and right hand sides of trust region constraint */
   lhs = 0.0;
   rhs = 0.0;

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

   /* adding the violation variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "violationvar", i);
   SCIP_CALL( SCIPcreateVarBasic(subscip, &violvar, name, 0.0, SCIPinfinity(subscip), heurdata->violpenalty, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(subscip, violvar) );
   consvars[nbinvars] = violvar;
   consvals[nbinvars] = -1.0;

   /* creates trustregion constraint and adds it to subscip */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_trustregioncons", SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, name, nbinvars + 1, consvars, consvals,
         lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

   /* create the upper bounding constraint */
   lhs = -SCIPinfinity(subscip);
   rhs = SCIPgetSolTransObj(scip, bestsol) - heurdata->objminimprove;

   /* if the objective function is integer, then the floor of the RHS is taken */
   if( SCIPisObjIntegral(scip) )
      rhs = SCIPfeasFloor(scip, rhs);

   /* adding the coefficients to the upper bounding constriant */
   for( i = 0; i < nvars; i++ )
   {
      consvals[i] = SCIPvarGetObj(subvars[i]);
      consvars[i] = subvars[i];
   }

   /* creates trustregion constraint and adds it to subscip */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_upperboundcons", SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, name, nvars, consvars, consvals,
         lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

   /* releasing the violation variable */
   SCIP_CALL( SCIPreleaseVar(subscip, &violvar) );

   /* free local memory */
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< SCIP data structure  of the original problem */
   SCIP*                 subscip,            /**< SCIP data structure  of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEUR*            heur,               /**< the Trustregion heuristic */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_Bool*            success             /**< pointer to store, whether new solution was found */
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

   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecTrustregion)
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
      SCIPdebugMsg(scip, "interrupt after  %" SCIP_LONGINT_FORMAT " LPs\n",SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTrustregion)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTrustregion(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTrustregion)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitTrustregion)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* with a little abuse we initialize the heurdata as if trustregion would have finished its last step regularly */
   heurdata->callstatus = WAITFORNEWSOL;
   heurdata->lastsol = NULL;
   heurdata->usednodes = 0;
   heurdata->curminnodes = heurdata->minnodes;

   return SCIP_OKAY;
}

/** todo setup And Solve Subscip */
static
SCIP_RETCODE setupAndSolveSubscipTrustregion(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< the subproblem created by trustregion */
   SCIP_HEUR*            heur,               /**< trustregion heuristic */
   SCIP_Longint          nsubnodes,          /**< nodelimit for subscip */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_VAR** subvars;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_HEURDATA* heurdata;
   SCIP_HASHMAP* varmapfw;
   SCIP_VAR** vars;

   int nvars;
   int i;

   SCIP_Bool success;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );
   success = FALSE;

   /* create a problem copy as sub SCIP */
   SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "trustregion", NULL, NULL, 0, heurdata->uselprows,
         heurdata->copycuts, &success, NULL) );

   SCIPdebugMsg(scip, "Copying SCIP was %ssuccessful.\n", success ? "" : "not ");

   /* if the subproblem could not be created, free memory and return */
   if( !success )
   {
      *result = SCIP_DIDNOTRUN;
      goto TERMINATE;
   }

   /* create event handler for LP events */
   eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecTrustregion, NULL) );
   if( eventhdlr == NULL )
   {
      /* free hash map */
      SCIPhashmapFree(&varmapfw);

      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   for (i = 0; i < nvars; ++i)
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* set limits for the subproblem */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );
   heurdata->nodelimit = nsubnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nsubnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", MAX(10, nsubnodes/10)) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->bestsollimit) );

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

   /* activate uct node selection at the top of the tree */
   if( heurdata->useuct && SCIPfindNodesel(subscip, "uct") != NULL && !SCIPisParamFixed(subscip, "nodeselection/uct/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/uct/stdpriority", INT_MAX/2) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* enable conflict analysis, disable analysis of boundexceeding LPs, and restrict conflict pool */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", TRUE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      SCIP_CALL( SCIPsetCharParam(subscip, "conflict/useboundlp", 'o') );
   }
   if( !SCIPisParamFixed(subscip, "conflict/maxstoresize") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "conflict/maxstoresize", 100) );
   }

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handler; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no deductions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 500) );
   }

   SCIP_CALL( addTrustRegionConstraints(scip, subscip, subvars, heurdata) );

   /* catch LP events of sub-SCIP */
   if( !heurdata->uselprows )
   {
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );
   }

   /* solve the subproblem */
   SCIPdebugMsg(scip, "solving trust region subproblem with maxnodes %" SCIP_LONGINT_FORMAT "\n", nsubnodes);

   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/trysol/priority", 100000) );

   /* Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* drop LP events of sub-SCIP */
   if( !heurdata->uselprows )
   {
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   heurdata->usednodes += SCIPgetNNodes(subscip);
   SCIPdebugMsg(scip, "trust region used %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT " nodes\n",
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
      for( i = 0; i < nsubsols; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );

         if( success )
            break;
      }
      if( success )
      {
         SCIPdebugMsg(scip, "-> accepted solution of value %g\n", SCIPgetSolOrigObj(subscip, subsols[i]));
         *result = SCIP_FOUNDSOL;
      }
   }

   /* check the status of the sub-MIP */
   switch( SCIPgetStatus(subscip) )
   {
   case SCIP_STATUS_OPTIMAL:
   case SCIP_STATUS_BESTSOLLIMIT:
      heurdata->callstatus = WAITFORNEWSOL; /* new solution will immediately be installed at next call */
      SCIPdebugMsg(scip, " -> found new solution\n");
      break;

   case SCIP_STATUS_NODELIMIT:
   case SCIP_STATUS_STALLNODELIMIT:
   case SCIP_STATUS_TOTALNODELIMIT:
      heurdata->callstatus = EXECUTE;
      heurdata->curminnodes *= 2;
      break;

   case SCIP_STATUS_INFEASIBLE:
   case SCIP_STATUS_INFORUNBD:
      heurdata->callstatus = WAITFORNEWSOL;
      break;

   case SCIP_STATUS_UNKNOWN:
   case SCIP_STATUS_USERINTERRUPT:
   case SCIP_STATUS_TERMINATE:
   case SCIP_STATUS_TIMELIMIT:
   case SCIP_STATUS_MEMLIMIT:
   case SCIP_STATUS_GAPLIMIT:
   case SCIP_STATUS_SOLLIMIT:
   case SCIP_STATUS_RESTARTLIMIT:
   case SCIP_STATUS_UNBOUNDED:
   default:
      heurdata->callstatus = WAITFORNEWSOL;
      SCIPdebugMsg(scip, " -> unexpected sub-MIP status <%d>: waiting for new solution\n", SCIPgetStatus(subscip));
      break;
   }

 TERMINATE:
   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTrustregion)
{  /*lint --e{715}*/
   SCIP_Longint maxnnodes;
   SCIP_Longint nsubnodes;

   SCIP_HEURDATA* heurdata;
   SCIP* subscip;

   SCIP_SOL* bestsol;

   SCIP_Bool success;
   SCIP_RETCODE retcode;

   int prevnwaitingnodes;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* there should be enough binary variables that a trust region constraint makes sense */
   if( SCIPgetNBinVars(scip) < heurdata->minbinvars )
      return SCIP_OKAY;

   *result = SCIP_DELAYED;

   /* only call heuristic, if an IP solution is at hand */
   if( SCIPgetNSols(scip) <= 0  )
      return SCIP_OKAY;

   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);

   /* only call heuristic, if the best solution comes from transformed problem */
   if( SCIPsolIsOriginal(bestsol) )
      return SCIP_OKAY;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip, bestsol)  < heurdata->nwaitingnodes)
      return SCIP_OKAY;

   /* only call heuristic, if the best solution does not come from trivial heuristic */
   if( SCIPsolGetHeur(bestsol) != NULL && strcmp(SCIPheurGetName(SCIPsolGetHeur(bestsol)), "trivial") == 0 )
      return SCIP_OKAY;

   /* reset minnodes if new solution was found */
   if( heurdata->lastsol != bestsol )
   {
      heurdata->curminnodes = heurdata->minnodes;
      heurdata->callstatus = EXECUTE;
      heurdata->lastsol = bestsol;
   }

   /* if no new solution was found and trust region also seems to fail, just keep on waiting */
   if( heurdata->callstatus == WAITFORNEWSOL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   maxnnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward trust region if it found solutions often.
    * In this case, the trust region heuristic is designed for Benders' decomposition and solutions found may not be
    * added by this heuristic but by trysol. So we don't reward finding best solutions, but finding any solution. */
   maxnnodes = (SCIP_Longint)(maxnnodes * (1.0 + 2.0*(SCIPheurGetNSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur)+1.0)));
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

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   /* abort if no time is left or not enough memory to create a copy of SCIP */
   if( !success )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMsg(scip, "running trust region heuristic ...\n");

   SCIP_CALL( SCIPcreate(&subscip) );

   retcode = setupAndSolveSubscipTrustregion(scip, subscip, heur, nsubnodes, result);

   SCIP_CALL( SCIPfree(&subscip) );

   /* if a solution is found, then we execute the trust region heuristic again */
   if( bestsol != SCIPgetBestSol(scip) )
   {
      prevnwaitingnodes = heurdata->nwaitingnodes;
      heurdata->nwaitingnodes = 0;
      SCIP_CALL( heurExecTrustregion(scip, heur, heurtiming, nodeinfeasible, result) );
      heurdata->nwaitingnodes = prevnwaitingnodes;
   }

   return retcode;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the trustregion primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTrustregion(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Trustregion primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTrustregion, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTrustregion) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTrustregion) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitTrustregion) );

   /* add trustregion primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minbinvars",
         "the number of binary variables necessary to run the heuristic",
         &heurdata->minbinvars, FALSE, DEFAULT_MINBINVARS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/bestsollimit",
         "limit on number of improving incumbent solutions in sub-CIP",
         &heurdata->bestsollimit, FALSE, DEFAULT_BESTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useuct",
         "should uct node selection be used at the beginning of the search?",
         &heurdata->useuct, TRUE, DEFAULT_USEUCT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/violpenalty",
         "the penalty for each change in the binary variables from the candidate solution",
         &heurdata->violpenalty, FALSE, DEFAULT_VIOLPENALTY, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/objminimprove",
         "the minimum improvement in the objective function value",
         &heurdata->objminimprove, FALSE, DEFAULT_OBJMINIMPROVE, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
