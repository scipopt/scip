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

/**@file   heur_rins.c
 * @brief  LNS heuristic that combines the incumbent with the LP optimum
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/heur_rins.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "rins"
#define HEUR_DESC             "relaxation induced neighborhood search by Danna, Rothberg, and Le Pape"
#define HEUR_DISPCHAR         'N'
#define HEUR_PRIORITY         -1101000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          5
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NODESOFS      500       /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_MAXNODES      5000      /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINNODES      500       /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which RINS should at least improve the incumbent          */
#define DEFAULT_MINFIXINGRATE 0.0       /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_NWAITINGNODES 200       /* number of nodes without incumbent change that heuristic should wait */
#define DEFAULT_USELPROWS    FALSE      /* should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE      /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip
                                         */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   int                   minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed       */
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   SCIP_Real             minimprove;         /**< factor by which RINS should at least improve the incumbent          */
   SCIP_Longint          usednodes;          /**< nodes already used by RINS in earlier calls                         */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?        */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?
                                              */
};

/*
 * Local methods
 */

/** creates a subproblem for subscip by fixing a number of variables */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure                                  */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                        */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                               */
   SCIP_Real             minfixingrate,      /**< percentage of integer variables that have to be fixed         */
   SCIP_Bool             uselprows,          /**< should subproblem be created out of the rows in the LP rows?   */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_SOL* bestsol;                        /* incumbent solution of the original problem */
   SCIP_VAR** vars;                          /* original scip variables                    */
   SCIP_ROW** rows;                          /* original scip rows                         */
   SCIP_Real fixingrate;

   int nrows;
   int nvars;
   int nbinvars;
   int nintvars;
   int i;
   int fixingcounter;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);

   fixingcounter = 0;

   /* change bounds of variables of the subproblem */
   for( i = 0; i <  nbinvars + nintvars; i++ )
   {
      SCIP_Real lpsolval;
      SCIP_Real solval;

      /* get the current LP solution and the incumbent solution for each variable */
      lpsolval = SCIPvarGetLPSol(vars[i]);
      solval = SCIPgetSolVal(scip, bestsol, vars[i]);

      /* iff both solutions are equal, variable is fixed to that value in the subproblem, otherwise it is just copied */
      if( SCIPisFeasEQ(scip, lpsolval, solval) )
      {
         /* perform the bound change */
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], solval) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], solval) );
         fixingcounter++;
      }
   }

   fixingrate = 0.0;

   /* abort, if all variables were fixed */
   if( fixingcounter == nbinvars + nintvars )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   else
      fixingrate = (SCIP_Real)fixingcounter / (SCIP_Real)(MAX(nbinvars + nintvars, 1));

   /* abort, if the amount of fixed variables is insufficient */
   if( fixingrate < minfixingrate )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   if( uselprows )
   {
      /* get the rows and their number */
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

      /* copy all rows to linear constraints */
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

         assert( lhs <= rhs );

         /* allocate memory array to be filled with the corresponding subproblem variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );
         for( j = 0; j < nnonz; j++ )
            consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

         /* create a new linear constraint and add it to the subproblem */
         SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

         /* free temporary memory */
         SCIPfreeBufferArray(scip, &consvars);
      }
   }

   *success = TRUE;
   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< RINS heuristic structure                            */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
)
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( subvars != NULL );
   assert( subsol != NULL );

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

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRins)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRins(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRins)
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
SCIP_DECL_HEURINIT(heurInitRins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* initialize data */
   heurdata->usednodes = 0;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRins)
{  /*lint --e{715}*/
   SCIP_Longint nstallnodes;

   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP* subscip;                            /* the subproblem created by RINS      */
   SCIP_VAR** vars;                          /* original problem's variables        */
   SCIP_VAR** subvars;                       /* subproblem's variables              */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */

   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem */
   SCIP_Real upperbound;

   int nvars;
   int nbinvars;
   int nintvars;
   int i;

   SCIP_Bool success;
   SCIP_RETCODE retcode;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );
   assert( SCIPhasCurrentNodeLP(scip) );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution and a feasible solution are at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL || SCIPgetNSols(scip) <= 0  )
      return SCIP_OKAY;

   /* only call heuristic, if the best solution comes from transformed problem */
   assert( SCIPgetBestSol(scip) != NULL );
   if( SCIPsolIsOriginal(SCIPgetBestSol(scip)) )
      return SCIP_OKAY;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip,SCIPgetBestSol(scip))  < heurdata->nwaitingnodes)
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward RINS if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* check whether discrete variables are available */
   if( nbinvars == 0 && nintvars == 0 )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* initializing the subproblem */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

   if( heurdata->uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

      /* get name of the original problem and add the string "_rinssub" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_rinssub", SCIPgetProbName(scip));

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

      /* create the subproblem */
      SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );
   }
   else
   {
      SCIP_Bool valid;
      valid = FALSE;

      SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "rins", TRUE, FALSE, TRUE, &valid) );

      if( heurdata->copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, TRUE, NULL) );
      }

      SCIPdebugMessage("Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");
   }

   for( i = 0; i < nvars; i++ )
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   success = FALSE;

   /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
   SCIP_CALL( createSubproblem(scip, subscip, subvars, heurdata->minfixingrate, heurdata->uselprows, &success) );

   if( !success )
   {
      *result = SCIP_DIDNOTRUN;
      goto TERMINATE;
   }

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
      goto TERMINATE;

   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nstallnodes) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", 3) );
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
   retcode = SCIPsolve(subscip);

   /* Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in RINS heuristic; sub-SCIP terminated with code <%d>\n",retcode);
   }

   heurdata->usednodes += SCIPgetNNodes(subscip);

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
         *result = SCIP_FOUNDSOL;
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

/** creates the RINS primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Rins primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRins, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRins) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRins) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRins) );

   /* add RINS primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which "HEUR_NAME" should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
         "minimum percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
