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

/**@file   sepa_rapidlearning.c
 * @ingroup SEPARATORS
 * @brief  rapidlearning separator
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#ifndef NDEBUG
#include <string.h>
#endif

#include "scip/sepa_rapidlearning.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_var.h"

#define SEPA_NAME              "rapidlearning"
#define SEPA_DESC               "rapid learning heuristic and separator"
#define SEPA_PRIORITY          -1200000
#define SEPA_FREQ                    -1 
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_APPLYCONFLICTS     TRUE /**< should the found conflicts be applied in the original SCIP?                 */
#define DEFAULT_APPLYBDCHGS        TRUE /**< should the found global bound deductions be applied in the original SCIP?   
					 *   apply only if conflicts and incumbent solution will be copied too
					 */
#define DEFAULT_APPLYINFERVALS     TRUE /**< should the inference values be used as initialization in the original SCIP? */
#define DEFAULT_APPLYPRIMALSOL     TRUE /**< should the incumbent solution be copied to the original SCIP?               */
#define DEFAULT_APPLYSOLVED        TRUE /**< should a solved status ba copied to the original SCIP?                      */

#define DEFAULT_MAXNVARS          10000 /**< maximum problem size (variables) for which rapid learning will be called */
#define DEFAULT_MAXNCONSS         10000 /**< maximum problem size (constraints) for which rapid learning will be called */

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_Bool             applyconflicts;     /**< should the found conflicts be applied in the original SCIP?                 */
   SCIP_Bool             applybdchgs;        /**< should the found global bound deductions be applied in the original SCIP?   */
   SCIP_Bool             applyinfervals;     /**< should the inference values be used as initialization in the original SCIP? */
   SCIP_Bool             applyprimalsol;     /**< should the incumbent solution be copied to the original SCIP?               */
   SCIP_Bool             applysolved;        /**< should a solved status ba copied to the original SCIP?                      */
   int                   maxnvars;           /**< maximum problem size (variables) for which rapid learning will be called   */
   int                   maxnconss;          /**< maximum problem size (constraints) for which rapid learning will be called */
};

/** creates a subproblem for subscip by fixing a number of variables */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure                                  */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                        */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                               */
   SCIP_HASHMAP*         varmapfw,           /**< mapping of SCIP variables to subSCIP variables                */
   SCIP_HASHMAP*         varmapbw,           /**< mapping of subSCIP variables to SCIP variables                */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_CONSHDLR** conshdlrs;
   SCIP_VAR** vars;                          /* original SCIP variables */
   int nvars;
   int i; 
 
   char consname[SCIP_MAXSTRLEN];
   
   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nvars == SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip));

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_rapidsub", SCIPgetProbName(scip));

   /* create the subproblem */
   SCIP_CALL( SCIPcreateProb(subscip, consname, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create the variables of the subproblem */
   for( i = 0; i < nvars; i++ )
   {     
      assert(SCIPvarGetLbLocal(vars[i]) == SCIPvarGetLbGlobal(vars[i]));
      assert(SCIPvarGetUbLocal(vars[i]) == SCIPvarGetUbGlobal(vars[i]));
      SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]), SCIPvarGetLbLocal(vars[i]),
            SCIPvarGetUbLocal(vars[i]), SCIPvarGetObj(vars[i]), 
            SCIPvarGetType(vars[i]) == SCIP_VARTYPE_IMPLINT ? SCIP_VARTYPE_INTEGER : SCIPvarGetType(vars[i]),
            SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, subvars[i]) );
      SCIP_CALL( SCIPhashmapInsert(varmapfw, vars[i], subvars[i]) );
   }

   conshdlrs = SCIPgetConshdlrs(scip);
 
   /* copy problem: loop through all constraint handlers */  
   for( i = 0; i < SCIPgetNConshdlrs(scip); ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONS** conss;
      SCIP_CONS* cons;
      SCIP_CONS* conscopy;
      SCIP_Bool succeed;
      int nconss;
      int c;

      conshdlr = conshdlrs[i];

      SCIPdebugMessage("rapid learning separator attempting to copy %d %s constraints\n", SCIPconshdlrGetNConss(conshdlr), SCIPconshdlrGetName(conshdlr));

      conss = SCIPconshdlrGetConss(conshdlr);
      nconss = SCIPconshdlrGetNConss(conshdlr);

      /* copy problem: loop through all constraints of one type */  
      for( c = 0; c < nconss; ++c )
      {
         cons = conss[c];
         assert(cons != NULL);

         /* copy each constraint */
         SCIP_CALL( SCIPcopyCons(subscip, &conscopy, NULL, conshdlr, scip, cons, varmapfw,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), TRUE, SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
               FALSE, &succeed) );

         /* add the copied constraint to subSCIP, print a warning if conshdlr does not support copying */
         if( succeed )
         {
            SCIP_CALL( SCIPaddCons(subscip, conscopy) );
            SCIP_CALL( SCIPreleaseCons(subscip, &conscopy) );
         }
         else
         {
            SCIPdebugMessage("failed to copy constraint %s\n", SCIPconsGetName(cons));
         }
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
   SCIP_HEUR*            heur,               /**< trysol heuristic structure                          */
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
   assert( heur != NULL );
   assert( subsol != NULL );
   assert( success != NULL );

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert( nvars == SCIPgetNOrigVars(subscip) );  
 
   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );
       
   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* check feasible of new solution and pass it to trysol heuristic */
   SCIP_CALL( SCIPcheckSol(scip, newsol, TRUE, TRUE, TRUE, success) );
   if( *success )
   {
      SCIPdebugMessage("Solution checking successful.\n");
      SCIP_CALL( SCIPheurPassSolTrySol(scip, heur, newsol) );
   }

   /* free allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeRapidlearning)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   SCIPfreeMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#define sepaInitRapidlearning NULL

/** deinitialization method of separator (called before transformed problem is freed) */
#define sepaExitRapidlearning NULL

/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolRapidlearning NULL

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#define sepaExitsolRapidlearning NULL

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpRapidlearning)
{/*lint --e{715}*/
   SCIP* subscip;                            /* the subproblem created by rapid learning       */
   SCIP_SEPADATA* sepadata;                  /* separator's private data                       */

   SCIP_VAR** vars;                          /* original problem's variables                   */
   SCIP_VAR** subvars;                       /* subproblem's variables                         */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to subSCIP variables */    
   SCIP_HASHMAP* varmapbw;                   /* mapping of subSCIP variables to SCIP variables */

   SCIP_CONSHDLR** conshdlrs;                /* array of constraint handler's that might that might obtain conflicts */
   int* oldnconss;                           /* number of constraints without rapid learning conflicts               */

   SCIP_Longint nodelimit;                   /* node limit for the subproblem                  */
   SCIP_Real timelimit;                      /* time limit for the subproblem                  */
   SCIP_Real memorylimit;                    /* memory limit for the subproblem                */

   int nconshdlrs;                           /* size of conshdlr and oldnconss array                      */
   int nfixedvars;                           /* number of variables that could be fixed by rapid learning */
   int nvars;                                /* number of variables                                       */           
   int restartnum;                           /* maximal number of conflicts that should be created        */
   int restarts;                             /* maximal number of restarts that should be performed       */
   int i;                                    /* counter                                                   */

   SCIP_Bool success;                        /* was problem creation / copying constraint successful? */

#ifdef NDEBUG
   SCIP_RETCODE retstat;                     /* used for catching subSCIP errors in debug mode */
#endif

   int nconflicts;                          /* statistic: number of conflicts applied         */
   int nbdchgs;                             /* statistic: number of bound changes applied     */
   int n1startinfers;                       /* statistic: number of one side infer values     */
   int n2startinfers;                       /* statistic: number of both side infer values    */

   SCIP_Bool soladded;                      /* statistic: was a new incumbent found?          */
   SCIP_Bool dualboundchg;                  /* statistic: was a new dual bound found?         */

   assert(sepa != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   
   /* only run when still not fixed binary variables exists */
   if( SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;

   /* only run for binary programs */
   if( SCIPgetNBinVars(scip) != SCIPgetNVars(scip) )
      return SCIP_OKAY;

   /* if the separator should be exclusive to the root node, this prevents multiple calls due to restarts */
   if(  SCIPsepaGetFreq(sepa) == 0 && SCIPsepaGetNCalls(sepa) > 0)
      return SCIP_OKAY;

   /* call separator at most once per node */
   if( SCIPsepaGetNCallsAtNode(sepa) > 0 )
      return SCIP_OKAY;

   /* get separator's data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* do not call rapid learning, if the probelm is too big */
   if( SCIPgetNVars(scip) > sepadata->maxnvars || SCIPgetNConss(scip) > sepadata->maxnconss )
      return SCIP_OKAY; 

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

   *result = SCIP_DIDNOTFIND;
   
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initializing the subproblem */  
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) ); 
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* mimic an FD solver: DFS, no LP solving, 1-FUIP instead of all-FUIP */
   SCIP_CALL( SCIPsetIntParam(subscip, "lp/solvefreq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "conflict/fuiplevels", 1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/dfs/stdpriority", INT_MAX/4) ); 
   SCIP_CALL( SCIPsetBoolParam(subscip, "constraints/disableenfops", TRUE) );
   SCIP_CALL( SCIPsetIntParam(subscip, "propagating/pseudoobj/freq", -1) );

   /* use inference branching */
   SCIP_CALL( SCIPsetBoolParam(subscip, "branching/inference/useweightedsum", FALSE) );

   /* only create short conflicts */
   SCIP_CALL( SCIPsetRealParam(subscip, "conflict/maxvarsfac", 0.05) );
  
   /* set limits for the subproblem */
   nodelimit = SCIPgetNLPIterations(scip);
   nodelimit = MAX(500, nodelimit);
   nodelimit = MIN(5000, nodelimit);
   restarts = 0;
   restartnum = 1000;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit/5) ); 
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/restarts", restarts) );
   SCIP_CALL( SCIPsetIntParam(subscip, "conflict/restartnum", restartnum) );

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/crossover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/undercover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rins/freq", -1) ); 
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rens/freq", -1) ); 
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/localbranching/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/mutation/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/dins/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/rapidlearning/freq", -1) );

   /* disable cut separation in sub problem */
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxroundsroot", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxcuts", 0) ); 
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxcutsroot", 0) );
   
   /* disable expensive presolving */
   SCIP_CALL( SCIPsetIntParam(subscip, "presolving/probing/maxrounds", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "constraints/linear/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "constraints/setppc/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "constraints/logicor/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetRealParam(subscip, "constraints/linear/maxaggrnormscale", 0.0) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifndef SCIP_DEBUG
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&varmapbw, SCIPblkmem(scip), nvars) );

   /* copy the problem */
   success = FALSE;
   SCIP_CALL( createSubproblem(scip, subscip, subvars, varmapfw, varmapbw, &success) );

   if( !success )
   {
     *result = SCIP_DIDNOTRUN;
     goto TERMINATE;
   } 

   /* add an objective cutoff */
   SCIP_CALL( SCIPsetObjlimit(subscip, SCIPgetUpperbound(scip)) );

   /* store reversing mapping of variables */
   SCIP_CALL( SCIPtransformProb(subscip) );
   for( i = 0; i < nvars; ++i)
   {  
      SCIP_CALL( SCIPhashmapInsert(varmapbw, SCIPvarGetTransVar(subvars[i]), vars[i]) );
   }

   /** allocate memory for constraints storage. Each constraint that will be created from now on will be a conflict.
    *  Therefore, we need to remember oldnconss to get the conflicts from the FD search. 
    */
   nconshdlrs = 4;
   SCIP_CALL( SCIPallocBufferArray(scip, &conshdlrs, nconshdlrs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oldnconss, nconshdlrs) );

   /* store number of constraints before rapid learning search */
   conshdlrs[0] = SCIPfindConshdlr(subscip, "bounddisjunction");
   conshdlrs[1] = SCIPfindConshdlr(subscip, "setppc");
   conshdlrs[2] = SCIPfindConshdlr(subscip, "linear");
   conshdlrs[3] = SCIPfindConshdlr(subscip, "logicor");

   /* redundant constraints might be eliminated in presolving */
   SCIP_CALL( SCIPpresolve(subscip));

   for( i = 0; i < nconshdlrs; ++i)
   {
      if( conshdlrs[i] != NULL )
         oldnconss[i] = SCIPconshdlrGetNConss(conshdlrs[i]);
   }

   nfixedvars = SCIPgetNFixedVars(scip);
   
   /* solve the subproblem.
    * Errors in the sub problem solver should not kill the overall solving process.
    * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
    */
#ifdef NDEBUG
   retstat = SCIPsolve(subscip);

   if( retstat != SCIP_OKAY )
   { 
      SCIPwarningMessage("Error while solving subMIP in rapid learning separator; subSCIP terminated with code <%d>\n", retstat);
   }
#else
   SCIP_CALL( SCIPsolve(subscip) );
#endif

 
   /* abort solving, if limit of applied conflicts is reached */
   if( SCIPgetNConflictConssApplied(subscip) >= restartnum && restarts == 0 )
   {
      SCIPdebugMessage("finish after %lld successful conflict calls.\n", SCIPgetNConflictConssApplied(subscip)); 
   }
   /* if the first 20% of the solution process were successful, proceed */
   else if( (sepadata->applyprimalsol && SCIPgetNSols(subscip) > 0 && SCIPisFeasLT(scip, SCIPgetUpperbound(subscip), SCIPgetUpperbound(scip) ) )
      || (sepadata->applybdchgs && SCIPgetNFixedVars(subscip) > nfixedvars)
      || (sepadata->applyconflicts && SCIPgetNConflictConssApplied(subscip) > 0) ) 
   {
      SCIPdebugMessage("proceed solving after the first 20%% of the solution process, since:\n");

      if( SCIPgetNSols(subscip) > 0 && SCIPisFeasLE(scip, SCIPgetUpperbound(subscip), SCIPgetUpperbound(scip) ) )
      {
         SCIPdebugMessage("   - there was a better solution (%f < %f)\n",SCIPgetUpperbound(subscip), SCIPgetUpperbound(scip));
      }
      if( SCIPgetNFixedVars(subscip) > nfixedvars )
      {
         SCIPdebugMessage("   - there were %d variables fixed\n", SCIPgetNFixedVars(scip)-nfixedvars );
      }
      if( SCIPgetNConflictConssFound(subscip) > 0 )
      {
         SCIPdebugMessage("   -  there were %lld conflict constraints created\n", SCIPgetNConflictConssApplied(subscip));
      }

      /* set node limit to 100% */
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) ); 

#ifdef NDEBUG
      retstat = SCIPsolve(subscip);

      if( retstat != SCIP_OKAY )
      { 
         SCIPwarningMessage("Error while solving subMIP in rapid learning separator; subSCIP terminated with code <%d>\n", retstat);
      }
#else
  
      SCIP_CALL( SCIPsolve(subscip) );
#endif
   }
   else
   {
     SCIPdebugMessage("do not proceed solving after the first 20%% of the solution process.\n");
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

  /* check, whether a solution was found */
   if( sepadata->applyprimalsol && SCIPgetNSols(subscip) > 0 && SCIPfindHeur(scip, "trysol") != NULL )
   {
      SCIP_HEUR* heurtrysol;
      SCIP_SOL** subsols;
      int nsubsols;

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until was declared to be feasible 
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      soladded = FALSE;
      heurtrysol = SCIPfindHeur(scip, "trysol");

      /* sequentially add solutions to trysol heuristic */
      for( i = 0; i < nsubsols && !soladded; ++i )
      {
	 SCIPdebugMessage("Try to create new solution by copying subscip solution.\n");
         SCIP_CALL( createNewSol(scip, subscip, subvars, heurtrysol, subsols[i], &soladded) );
      }
   }

   /* if the sub problem was solved completely, we update the dual bound */
   dualboundchg = FALSE;
   if( sepadata->applysolved && (SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL || SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE) )
   {
      /* we need to multiply the dualbound with the scaling facting and add the offset, 
       * because thiinformation was   disregarded in the subscip */
      SCIPdebugMessage("Update old dualbound %g to new dualbound %g.\n", SCIPgetDualbound(scip), SCIPgetTransObjscale(scip) * SCIPgetDualbound(subscip) + SCIPgetTransObjoffset(scip));
      SCIP_CALL( SCIPupdateLocalDualbound(scip, SCIPgetDualbound(subscip) * SCIPgetTransObjscale(scip) + SCIPgetTransObjoffset(scip)) );
      dualboundchg = TRUE;
   }

   /* check, whether conflicts were created */
   nconflicts = 0;
   if( sepadata->applyconflicts && SCIPgetNConflictConssApplied(subscip) > 0 )
   {
      /* loop over all constraint handlers that might contain conflict constraints */
      for( i = 0; i < nconshdlrs; ++i)
      {
         /* copy constraints that have been created in FD run */
         if( conshdlrs[i] != NULL && SCIPconshdlrGetNConss(conshdlrs[i]) > oldnconss[i] )
         {
	    SCIP_CONS** conss;
            int c;
            int nconss;
            
	    nconss = SCIPconshdlrGetNConss(conshdlrs[i]);
	    conss = SCIPconshdlrGetConss(conshdlrs[i]);

            /* loop over all constraints that have been added in subSCIP run, these are the conflicts */            
            for( c = oldnconss[i]; c < nconss; ++c)
            {
               SCIP_CONS* cons;
               SCIP_CONS* conscopy;
               
               cons = conss[c];
               assert(cons != NULL);        
               
               SCIP_CALL( SCIPcopyCons(scip, &conscopy, NULL, conshdlrs[i], subscip, cons, varmapbw,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), TRUE, SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
                     SCIPconsIsRemovable(cons), FALSE, &success) );

               if( success )
               {
                  nconflicts++;
                  SCIP_CALL( SCIPaddCons(scip, conscopy) );
                  SCIP_CALL( SCIPreleaseCons(scip, &conscopy) );
               }
               else
               {
                  SCIPdebugMessage("failed to copy constraint %s\n", SCIPconsGetName(cons));
               }
            }
         }
      }   
   }

   /* check, whether tighter global bounds were detected */
   nbdchgs = 0;
   if( sepadata->applybdchgs )
      for( i = 0; i < nvars; ++i )
      {
	 SCIP_Bool infeasible;
	 SCIP_Bool tightened;
	 
	 assert(SCIPisLE(scip, SCIPvarGetLbGlobal(vars[i]), SCIPvarGetLbGlobal(subvars[i]))); 
	 assert(SCIPisLE(scip, SCIPvarGetLbGlobal(subvars[i]), SCIPvarGetUbGlobal(subvars[i])));
	 assert(SCIPisLE(scip, SCIPvarGetUbGlobal(subvars[i]), SCIPvarGetUbGlobal(vars[i])));  
	 
	 /* update the bounds of the original SCIP, if a better bound was proven in the subSCIP */
	 SCIP_CALL( SCIPtightenVarUb(scip, vars[i], SCIPvarGetUbGlobal(subvars[i]), FALSE, &infeasible, &tightened) );
	 if( tightened ) 
	   nbdchgs++;
	 
	 SCIP_CALL( SCIPtightenVarLb(scip, vars[i], SCIPvarGetLbGlobal(subvars[i]), FALSE, &infeasible, &tightened) );
	 if( tightened )
	   nbdchgs++;   
      }

   n1startinfers = 0;
   n2startinfers = 0;

   /* install start values for inference branching */
   if( sepadata->applyinfervals )
   {
      for( i = 0; i < nvars; ++i )
      {
         SCIP_Longint downval;
         SCIP_Longint upval;
         
         /* copy downwards inference value to original SCIP */
         downval = 0;
         if( SCIPvarGetNBranchings(subvars[i], SCIP_BRANCHDIR_DOWNWARDS) > 0 )
            downval = (SCIP_Longint) (SCIPvarGetNInferences(subvars[i], SCIP_BRANCHDIR_DOWNWARDS) / SCIPvarGetNBranchings(subvars[i], SCIP_BRANCHDIR_DOWNWARDS));
         
         /* copy upwards inference value to original SCIP */
         upval = 0;
         if( SCIPvarGetNBranchings(subvars[i], SCIP_BRANCHDIR_UPWARDS) > 0 )
            upval = (SCIP_Longint) (SCIPvarGetNInferences(subvars[i], SCIP_BRANCHDIR_UPWARDS) / SCIPvarGetNBranchings(subvars[i], SCIP_BRANCHDIR_UPWARDS));

         /* memorize statistics */
         if( downval != 0 || upval != 0 )
            n1startinfers++;

         if( downval != 0 && upval != 0 )
            n2startinfers++;

         SCIP_CALL( SCIPsetVarNInferencesInitial(scip, vars[i], downval, upval) );
      }   
   }
   
   SCIPdebugMessage("Rapidlearning added %d conflicts, changed %d bounds, %s primal solution, %s dual bound improvement.\n", nconflicts, nbdchgs, soladded ? "found" : "no", 
      dualboundchg ? "found" : "no");

   SCIPdebugMessage("Infervalues initialized on one side: %5.2f %% of variables, %5.2f %% on both sides\n", n1startinfers/(SCIP_Real)nvars, n2startinfers/(SCIP_Real)nvars);

   /* change result pointer */
   if( nconflicts > 0 || dualboundchg )
      *result = SCIP_CONSADDED;
   else if( nbdchgs > 0 )
      *result = SCIP_REDUCEDDOM;
  
   /* free local data */
   SCIPfreeBufferArray(scip, &oldnconss);
   SCIPfreeBufferArray(scip, &conshdlrs);

 TERMINATE:
   
   SCIPhashmapFree(&varmapbw);
   SCIPhashmapFree(&varmapfw);

   /* free subproblem */
   SCIP_CALL( SCIPfreeTransform(subscip) );
   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
   }
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );
  
   return SCIP_OKAY;
}

/** arbitrary primal solution separation method of separator */
#define sepaExecsolRapidlearning NULL

/*
 * separator specific interface methods
 */

/** creates the rapidlearning separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaRapidlearning(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create rapidlearning separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_DELAY,
         sepaFreeRapidlearning, sepaInitRapidlearning, sepaExitRapidlearning, 
         sepaInitsolRapidlearning, sepaExitsolRapidlearning,
         sepaExeclpRapidlearning, sepaExecsolRapidlearning,
         sepadata) );

   /* add rapidlearning separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/rapidlearning/applyconflicts",
         "should the found conflicts be applied in the original SCIP?",
         &sepadata->applyconflicts, TRUE, DEFAULT_APPLYCONFLICTS, NULL, NULL) );
  
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/rapidlearning/applybdchgs",
         "should the found global bound deductions be applied in the original SCIP?",
         &sepadata->applybdchgs, TRUE, DEFAULT_APPLYBDCHGS, NULL, NULL) );
  
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/rapidlearning/applyinfervals",
         "should the inference values be used as initialization in the original SCIP?",
         &sepadata->applyinfervals, TRUE, DEFAULT_APPLYINFERVALS, NULL, NULL) );
  
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/rapidlearning/applyprimalsol",
         "should the incumbent solution be copied to the original SCIP?",
         &sepadata->applyprimalsol, TRUE, DEFAULT_APPLYPRIMALSOL, NULL, NULL) );
  
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/rapidlearning/applysolved",
         "should a solved status ba copied to the original SCIP?",
         &sepadata->applysolved, TRUE, DEFAULT_APPLYSOLVED, NULL, NULL) );
 
   SCIP_CALL( SCIPaddIntParam(scip, "separating/rapidlearning/maxnvars",
         "maximum problem size (variables) for which rapid learning will be called",
         &sepadata->maxnvars, TRUE, DEFAULT_MAXNVARS, 0, INT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddIntParam(scip, "separating/rapidlearning/maxnconss",
         "maximum problem size (constraints) for which rapid learning will be called",
         &sepadata->maxnconss, TRUE, DEFAULT_MAXNCONSS, 0, INT_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}
