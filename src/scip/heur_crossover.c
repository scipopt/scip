/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_crossover.c,v 1.1 2005/11/02 14:28:06 bzfberth Exp $"

/**@file   heur_crossover.c
 * @brief  crossover primal heuristic
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/heur_crossover.h"

#define HEUR_NAME             "crossover"
#define HEUR_DESC             "LNS heuristic that fixes all variables that are identic in a couple of solutions"
#define HEUR_DISPCHAR         'C'
#define HEUR_PRIORITY         -1000000
#define HEUR_FREQ             30
#define HEUR_FREQOFS          10
#define HEUR_MAXDEPTH         -1
#define HEUR_PSEUDONODES      FALSE     /* call heuristic at nodes where only a pseudo solution exist?   */
#define HEUR_DURINGPLUNGING   TRUE      /* call heuristic during plunging?                               */
#define HEUR_DURINGLPLOOP     FALSE     /* call heuristic during the LP price-and-cut loop? */
#define HEUR_AFTERNODE        TRUE      /* call heuristic after or before the current node was solved?   */

#define DEFAULT_NODESOFS      500       /* number of nodes added to the contingent of the total nodes    */
#define DEFAULT_MAXNODES      5000      /* maximum number of nodes to regard in the subproblem           */
#define DEFAULT_MINNODES      500       /* minimum number of nodes to regard in the subproblem           */
#define DEFAULT_MINFIXINGRATE 0.666     /* minimum percentage of integer variables that have to be fixed */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_NUSEDSOLS     3         /* number of solutions that will be taken into account           */
#define DEFAULT_NWAITINGNODES 200       /* number of nodes without incumbent change that heuristic should wait */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             lastsol;            /**< worst solution taken into account during the last run         */
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes    */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem           */
   int                   minnodes;           /**< minimum number of nodes to regard in the subproblem           */
   int                   nusedsols;          /**< number of solutions that will be taken into account           */
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed */
   SCIP_Longint          usednodes;          /**< nodes already used by crossover in earlier calls      */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
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
   int                   nusedsols,          /**< number of solutions to be taken into account                  */
   SCIP_Real             minfixingrate,      /**< percentage of integer variables that have to be fixed         */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_VAR** vars;                          /* original scip variables                */
   SCIP_ROW** rows;                          /* original scip rows                     */
   SCIP_SOL** sols;                          /* pool of solutions                      */
   SCIP_Real fixingrate;                     /* percentage of variables that are fixed */
   int nsols;
   int nrows;
   int nvars; 
   int nbinvars;
   int nintvars;
   int i;
   int fixingcounter;
   char consname[SCIP_MAXSTRLEN];  

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   sols = SCIPgetSols(scip);
   nsols = SCIPgetNSols(scip);

   assert( sols != NULL );
   assert( nusedsols > 1);

   /* get name of the original problem and add the string "_crossoversub" */
   sprintf(consname, "%s_crossoversub", SCIPgetProbName(scip));

   /* create the subproblem */
   SCIP_CALL( SCIPcreateProb(subscip, consname, NULL, NULL, NULL, NULL, NULL, NULL) );
   fixingcounter = 0;

   /* create the binary and general integer variables of the subproblem */
   for( i = 0; i < nbinvars + nintvars; i++ )
   {
      int j;

      /* get the current LP solution and the incumbent solution for each variable */
      SCIP_Real solval;
      SCIP_Bool fixable = TRUE;
      int nused;

      solval = 0.0;
      for( j = 0, nused = 0; nused < nusedsols && j < nsols; j++ )
      {
         SCIP_Real varsolval;

         if( SCIPsolGetOrigin(sols[j]) == SCIP_SOLORIGIN_ORIGINAL )
            continue;
         varsolval = SCIPgetSolVal(scip, sols[j], vars[i]);
         if( nused == 0 )
            solval = varsolval;
         else if( REALABS(solval - varsolval) > 0.5 )
         {
            fixable = FALSE;
            break;
         }
         nused++;
      }

      /* check if we actually had enough solutions in the transformed problem space. Abort, if not */
      if( fixable && nused < nusedsols )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      /* iff both solutions are equal, variable is fixed to that value in the subproblem, otherwise it is just copied */
      if( fixable )
      {
         SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]), solval,
               solval, SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]),
               SCIPvarIsInitial(vars[i]), SCIPvarIsRemoveable(vars[i]), NULL, NULL, NULL, NULL) );
         fixingcounter++;
      }
      else 
      {
         SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]), SCIPvarGetLbGlobal(vars[i]),
               SCIPvarGetUbGlobal(vars[i]), SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]),
               SCIPvarIsInitial(vars[i]), SCIPvarIsRemoveable(vars[i]), NULL, NULL, NULL, NULL) );
      }
      SCIP_CALL( SCIPaddVar(subscip, subvars[i]) );
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
     
   /* create the continuous variables of the subproblem */  
   for( i = nbinvars + nintvars; i < nvars; i++ )
   {
      SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]), SCIPvarGetLbGlobal(vars[i]),
            SCIPvarGetUbGlobal(vars[i]), SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]),
            SCIPvarIsInitial(vars[i]), SCIPvarIsRemoveable(vars[i]), NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, subvars[i]) );
   }

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
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      
      /* free temporary memory */
      SCIPfreeBufferArray(scip, &consvars);
   }

   *success = TRUE;
   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_HEUR*            heur,               /**< crossover heuristic structure               */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
)
{
   SCIP_VAR** subvars;                       /* the subproblem's variables                      */
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
   SCIP_SOL*  subsol;                        /* incumbent of the subproblem                     */
        
   assert( scip != NULL );
   assert( subscip != NULL );

   subsol = SCIPgetBestSol(subscip);
   assert( subsol != NULL );

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert( nvars == SCIPgetNOrigVars(subscip) );  
 
   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   subvars = SCIPgetOrigVars(subscip); 
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );
       
   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeCrossover)
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
SCIP_DECL_HEURINIT(heurInitCrossover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->lastsol = NULL;

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitCrossover NULL

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolCrossover NULL

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolCrossover NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecCrossover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;                  /* heuristic's data                                    */
   SCIP* subscip;                            /* the subproblem created by crossover         */
   SCIP_VAR** vars;                          /* original problem's variables                        */
   SCIP_VAR** subvars;                       /* subproblem's variables                              */
   SCIP_SOL** sols;
   SCIP_Real timelimit;                      /* timelimit for the subproblem                        */
   int nvars;                                /* number of original problem's variables              */
   int nusedsols;                            /* number of solutions that will be taken into account */
   SCIP_Bool success;
   int i;   
   SCIP_Longint maxnnodes;                  
   SCIP_Longint nsubnodes;                   /* node limit for the subproblem                       */
     
   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );
   assert( SCIPhasCurrentNodeLP(scip) );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );
   
   nusedsols = heurdata->nusedsols;
   assert( nusedsols > 1);

   *result = SCIP_DELAYED;

   /* only call heuristic, if enough solutions are at hand */
   if( SCIPgetNSols(scip) < nusedsols  )
      return SCIP_OKAY;

   /* only recall heuristic, if at least one new good solution was found in the meantime */
   sols = SCIPgetSols(scip);
   if( sols[nusedsols-1] == heurdata->lastsol )
      return SCIP_OKAY;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip,SCIPgetBestSol(scip))  < heurdata->nwaitingnodes)
      return SCIP_OKAY;
   
   *result = SCIP_DIDNOTRUN;
   
   /* calculate the maximal number of branching nodes until heuristic is aborted */
   maxnnodes = heurdata->nodesquot * SCIPgetNNodes(scip);

   /* reward crossover if it succeeded often */
   maxnnodes *= 1.0 + 2.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0);
   maxnnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   maxnnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nsubnodes = maxnnodes - heurdata->usednodes;
   nsubnodes = MIN(nsubnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nsubnodes < heurdata->minnodes )
       return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initializing the subproblem */  
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) ); 
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );
 
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
  
   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nsubnodes) ); 
 
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit - SCIPgetTotalTime(scip) + 10.0) );

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rins/freq", -1) ); 
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/localbranching/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/crossover/freq", -1) );

   /* disable cut separation in sub problem */
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxroundsroot", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxcuts", 0) ); 
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxcutsroot", 0) );
   
   /* use pseudo cost branching without strong branching */
   SCIP_CALL( SCIPsetIntParam(subscip, "branching/pscost/priority", INT_MAX) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetIntParam(subscip, "presolving/probing/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "constraints/linear/maxpresolpairrounds", 0) );
   SCIP_CALL( SCIPsetRealParam(subscip, "constraints/linear/maxaggrnormscale", 0.0) );

   /* disable conflict analysis */
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) ); 
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/uselp", FALSE) ); 
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) ); 
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) );

   success = FALSE;
   /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
   createSubproblem(scip, subscip, subvars, nusedsols, heurdata->minfixingrate, &success);
   heurdata->lastsol = sols[nusedsols-1];
  
   /* if creation of subscip was aborted (e.g. due to number of fixings), free subscip and abort */
   if( !success )
   {
      int nbinvars;
      int nintvars;
      *result = SCIP_DIDNOTRUN;
      SCIP_CALL( SCIPfreeTransform(subscip) );
      SCIP_CALL( SCIPgetVarsData(subscip, NULL, NULL, &nbinvars, &nintvars, NULL, NULL) );
      for( i = 0; i < nbinvars + nintvars; i++ )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
      }
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
      return SCIP_OKAY;
   }
   
   /* add an objective cutoff */
   SCIP_CALL( SCIPsetObjlimit(subscip, SCIPgetSolTransObj(scip, sols[0]) - SCIPsumepsilon(scip)) );

   /* solve the subproblem */
   SCIP_CALL( SCIPsolve(subscip) );
   
   heurdata->usednodes += SCIPgetNNodes(subscip);

   /* check, whether a solution was found */
   if( SCIPgetNSols(subscip) > 0 )
   {
      success = FALSE;
      SCIP_CALL( createNewSol(scip, subscip, heur, &success) );
      if( success )
      {
         *result = SCIP_FOUNDSOL;
      }
   }
   
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

/*
 * primal heuristic specific interface methods
 */

/** creates the crossover primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurCrossover(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */ 
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING, HEUR_DURINGLPLOOP, HEUR_AFTERNODE,
         heurFreeCrossover, heurInitCrossover, heurExitCrossover,
         heurInitsolCrossover, heurExitsolCrossover, heurExecCrossover,
         heurdata) );
  
   /* add crossover primal heuristic parameters */ 
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/crossover/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/crossover/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/crossover/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/crossover/nusedsols",
         "number of solutions to be takten into account",
         &heurdata->nusedsols, DEFAULT_NUSEDSOLS, 2, INT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/crossover/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, DEFAULT_NWAITINGNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/crossover/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );
   
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/crossover/minfixingrate",
         "minimum percentage of integer variables that have to be fixed ",
         &heurdata->minfixingrate, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );
   
   return SCIP_OKAY;
}
