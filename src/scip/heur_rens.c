/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_rens.c,v 1.4 2006/06/20 20:24:01 bzfpfend Exp $"

/**@file   heur_rens.c
 * @brief  RENS primal heuristic
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/heur_rens.h"

#define HEUR_NAME             "rens"
#define HEUR_DESC             "LNS exploring fractional neighborhood of relaxation's optimum"
#define HEUR_DISPCHAR         'E'
#define HEUR_PRIORITY         -1100000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE

#define DEFAULT_BINARYBOUNDS  FALSE     /* should general integers get binary bounds [floor(.),ceil(.)] ?      */ 
#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINFIXINGRATE 0.5       /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which RENS should at least improve the incumbent          */
#define DEFAULT_MINNODES      500LL     /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_NODESOFS      500LL     /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem       */



/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;          /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;          /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          nodesofs;          /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          usednodes;         /**< nodes already used by RENS in earlier calls                         */

   SCIP_Real             minfixingrate;     /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Real             minimprove;        /**< factor by which RENS should at least improve the incumbent          */
   SCIP_Real             nodesquot;         /**< subproblem nodes in relation to nodes of the original problem       */
   SCIP_Real             nsuccesses;        /**< number of RENS-calls, where a real improvement was achieved         */

   SCIP_Bool             binarybounds;      /**< should general integers get binary bounds [floor(.),ceil(.)] ?      */
};


/*
 * Local methods
 */

/** creates a subproblem for subscip by fixing a number of variables */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                         */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                                */
   SCIP_Real             minfixingrate,      /**< percentage of integer variables that have to be fixed          */
   SCIP_Bool             binarybounds,       /**< should general integers get binary bounds [floor(.),ceil(.)] ? */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully  */
   )
{
   SCIP_VAR** vars;                          /* original scip variables                    */
   SCIP_ROW** rows;                          /* original scip rows                         */
   SCIP_Real fixingrate;

   int nrows;
   int nvars;   
   int nbinvars;
   int nintvars;
   int i; 
   int fixingcounter;

   char consname[SCIP_MAXSTRLEN];
   
   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* get name of the original problem and add the string "_renssub" */
   sprintf(consname, "%s_renssub", SCIPgetProbName(scip));

   /* create the subproblem */
   SCIP_CALL( SCIPcreateProb(subscip, consname, NULL, NULL, NULL, NULL, NULL, NULL) );
   fixingcounter = 0;

   /* create the variables of the subproblem */
   for( i = 0; i < nbinvars + nintvars; i++ )
   {
      SCIP_Real lpsolval;

      /* get the current LP solution and the incumbent solution for each variable */
      lpsolval = SCIPvarGetLPSol(vars[i]);
      
      /* iff both solutions are equal, variable is fixed to that value in the subproblem, otherwise it is just copied */
      if( SCIPisFeasIntegral(scip, lpsolval) )
      {
         SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]), lpsolval,
               lpsolval, SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]),
               SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL, NULL) );
         fixingcounter++;
      }
      else 
      {
         /* optionally, general integers can get binary bounds, too */
         if( i >= nbinvars && binarybounds )
         {
            SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]), SCIPfeasFloor(scip,lpsolval),
                  SCIPfeasCeil(scip,lpsolval), SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]),
                  SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL, NULL) );
         }
         else
         {
            SCIP_CALL( SCIPcreateVar(subscip, &subvars[i], SCIPvarGetName(vars[i]), SCIPvarGetLbGlobal(vars[i]),
                  SCIPvarGetUbGlobal(vars[i]), SCIPvarGetObj(vars[i]), SCIPvarGetType(vars[i]),
                  SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL, NULL) );
         }
      }
      SCIP_CALL( SCIPaddVar(subscip, subvars[i]) );
   }
      
   fixingrate = 0.0;

   /* abort, if all variables were fixed (which should not happen) */
   if( fixingcounter == nbinvars + nintvars )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   else
      fixingrate = fixingcounter / (SCIP_Real)(MAX(nbinvars + nintvars, 1));

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
            SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL, NULL) );
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
      SCIP_CALL( SCIPallocBufferArray(subscip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ ) 
         consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      
      /* free temporary memory */
      SCIPfreeBufferArray(subscip, &consvars);
   }

   *success = TRUE;
   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_HEUR*            heur,               /**< RENS heuristic structure                            */
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
SCIP_DECL_HEURFREE(heurFreeRens)
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
SCIP_DECL_HEURINIT(heurInitRens)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->nsuccesses = 0;

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitRens NULL

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolRens NULL

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolRens NULL

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRens)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP* subscip;                            /* the subproblem created by RENS      */
   SCIP_VAR** vars;                          /* original problem's variables        */
   SCIP_VAR** subvars;                       /* subproblem's variables              */
  
   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem */
   
   SCIP_Bool success;
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */

   int nvars;                     
   int i;   

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );
   assert( SCIPhasCurrentNodeLP(scip) );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;
   
   *result = SCIP_DIDNOTRUN;
 
   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));
   
   /* reward RENS if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (heurdata->nsuccesses+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
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

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initializing the subproblem */  
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) ); 
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );
 
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
 
   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) ); 
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) ); 
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/crossover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rens/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rins/freq", -1) ); 
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/localbranching/freq", -1) );

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
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useinflp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useboundlp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) );

   success = FALSE;

   /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
   SCIP_CALL( createSubproblem(scip, subscip, subvars, heurdata->minfixingrate, heurdata->binarybounds, &success) );

   /* if the subproblem could not be created, free memory and return */
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

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      cutoff = SCIPinfinity(scip);
      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );   
      
      if( !SCIPisInfinity(scip,SCIPgetLowerbound(scip)) )
      {
         SCIP_Real upperbound;
         cutoff = (1-heurdata->minimprove)*SCIPgetUpperbound(scip) - heurdata->minimprove*SCIPgetLowerbound(scip);
         upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
         cutoff = MIN(upperbound, cutoff );
      }
      else
         cutoff = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );
   }

   /* solve the subproblem */
   SCIP_CALL( SCIPpresolve(subscip) );

   /* after presolving, we should have at least reached a certain fixing rate over ALL variables (including continuous)
    * to ensure that not only the MIP but also the LP relaxation is easy enough
    */
   if( ( nvars - SCIPgetNVars(subscip) ) / (SCIP_Real)nvars >= heurdata->minfixingrate / 2.0 )
   {
      SCIP_CALL( SCIPsolve(subscip) );
      
      /* check, whether a solution was found */
      if( SCIPgetNSols(subscip) > 0 )
      {
         success = FALSE;
         SCIP_CALL( createNewSol(scip, subscip, heur, &success) );
         if( success )
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

/** creates the rens primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRens(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING,
         heurFreeRens, heurInitRens, heurExitRens, 
         heurInitsolRens, heurExitsolRens, heurExecRens,
         heurdata) );

   /* add rens primal heuristic parameters */
 
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/rens/minfixingrate",
         "minimum percentage of integer variables that have to be fixed ",
         &heurdata->minfixingrate, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );
   
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/rens/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
 
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/rens/binarybounds",
         "should general integers get binary bounds [floor(.),ceil(.)] ?",
         &heurdata->binarybounds, DEFAULT_BINARYBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/rens/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/rens/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
 
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/rens/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/rens/minimprove",
         "factor by which RENS should at least improve the incumbent  ",
         &heurdata->minimprove, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );
   
   return SCIP_OKAY;
}
