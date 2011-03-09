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

/**@file   heur_rens.c
 * @ingroup PRIMALHEURISTICS
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
#include "scip/pub_misc.h"

#define HEUR_NAME             "rens"
#define HEUR_DESC             "LNS exploring fractional neighborhood of relaxation's optimum"
#define HEUR_DISPCHAR         'E'
#define HEUR_PRIORITY         -1100000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE

#define DEFAULT_BINARYBOUNDS  TRUE      /* should general integers get binary bounds [floor(.),ceil(.)] ?      */ 
#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINFIXINGRATE 0.5       /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which RENS should at least improve the incumbent          */
#define DEFAULT_MINNODES      500LL     /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_NODESOFS      500LL     /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_USELPROWS     TRUE      /* should subproblem be created out of the rows in the LP rows, 
                                         * otherwise, the copy constructor of the constraints handlers are used*/



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
   SCIP_Bool             binarybounds;      /**< should general integers get binary bounds [floor(.),ceil(.)] ?      */
   SCIP_Bool             uselprows;         /**< should subproblem be created out of the rows in the LP rows?        */
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
   SCIP_HASHMAP*         varmapfw,           /**< mapping of SCIP variables to subSCIP variables                */
   SCIP_Real             minfixingrate,      /**< percentage of integer variables that have to be fixed          */
   SCIP_Bool             binarybounds,       /**< should general integers get binary bounds [floor(.),ceil(.)] ? */
   SCIP_Bool             uselprows,          /**< should subproblem be created out of the rows in the LP rows?   */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully  */
   )
{
   SCIP_VAR** vars;                          /* original scip variables                    */
   SCIP_Real fixingrate;

   int nvars;   
   int nbinvars;
   int nintvars;
   int i; 
   int fixingcounter;

   char consname[SCIP_MAXSTRLEN];
   
   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* get name of the original problem and add the string "_renssub" */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_renssub", SCIPgetProbName(scip));

   /* create the subproblem */
   SCIP_CALL( SCIPcreateProb(subscip, consname, NULL, NULL, NULL, NULL, NULL, NULL) );
   fixingcounter = 0;

   /* create the variables of the subproblem */
   for( i = 0; i < nbinvars + nintvars; i++ )
   {
      SCIP_Real lpsolval;

      /* get the current LP solution for each variable */
      lpsolval = SCIPvarGetLPSol(vars[i]);
      
      /* iff variable is integer, it is fixed in the subproblem, otherwise it is just copied */
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

      /* insert variable into mapping between SCIP and the subSCIP */
      SCIP_CALL( SCIPhashmapInsert(varmapfw, vars[i], subvars[i]) );
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

      /* insert variable into mapping between SCIP and the subSCIP */
      SCIP_CALL( SCIPhashmapInsert(varmapfw, vars[i], subvars[i]) );
   }

   if( uselprows )
   {
      SCIP_ROW** rows;                          /* original scip rows                         */
      int nrows;
      
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
   }
   else
   {
      /* use the copy constructor of each constraint handler to create subSCIP */
      
      /**@todo it might be also of interest to copy the cuts, like in case of restart (see cons_linear.c) */

      /* copy problem: loop through all constraint handlers */  
      for( i = 0; i < SCIPgetNConshdlrs(scip); ++i )
      {
         SCIP_CONSHDLR* conshdlr;
         SCIP_CONS* cons;
         SCIP_CONS* conscopy;
         SCIP_Bool succeed;
         int c;
         
         conshdlr = SCIPgetConshdlrs(scip)[i];
         
         SCIPdebugMessage("rens heuristic attempting to copy %d %s constraints\n", SCIPconshdlrGetNConss(conshdlr), SCIPconshdlrGetName(conshdlr));
         
         /* copy problem: loop through all constraints of one type */  
         for( c = 0; c < SCIPconshdlrGetNConss(conshdlr); ++c )
         {
            cons = SCIPconshdlrGetConss(conshdlr)[c];
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
   SCIP_HEUR*            heur,               /**< RENS heuristic structure                            */
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
   assert( nvars == SCIPgetNOrigVars(subscip) );  
 
   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );
       
   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** main procedure of the RENS heuristic, creates and solves a subMIP */
SCIP_RETCODE SCIPapplyRens(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< result data structure                                          */
   SCIP_Real             timelimit,          /**< timelimit for the subproblem                                   */        
   SCIP_Real             memorylimit,        /**< memorylimit for the subproblem                                 */
   SCIP_Real             minfixingrate,      /**< minimum percentage of integer variables that have to be fixed  */
   SCIP_Real             minimprove,         /**< factor by which RENS should at least improve the incumbent     */
   SCIP_Longint          maxnodes,           /**< maximum number of  nodes for the subproblem                    */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem                    */
   SCIP_Bool             binarybounds,       /**< should general integers get binary bounds [floor(.),ceil(.)]?  */
   SCIP_Bool             uselprows           /**< should subproblem be created out of the rows in the LP rows?   */
   )
{
   SCIP* subscip;                            /* the subproblem created by RENS      */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to subSCIP variables */    
   SCIP_VAR** vars;                          /* original problem's variables        */
   SCIP_VAR** subvars;                       /* subproblem's variables              */
  
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem */
   
   SCIP_Bool success;

   int nvars;                     
   int i;   

#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initializing the subproblem */  
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) ); 
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );
 
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
 
   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) ); 
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", maxnodes) ); 
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rens/freq", -1) );
#if 1
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/undercover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/crossover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/oneopt/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rins/freq", -1) ); 
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/localbranching/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/mutation/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/dins/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/rapidlearning/freq", -1) ); 

   /* use best estimate node selection */
   SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) ); 

   /* disable cut separation in sub problem */
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxroundsroot", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxcuts", 0) ); 
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/maxcutsroot", 0) );

   /* use inference branching */
   SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetIntParam(subscip, "presolving/probing/maxrounds", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "constraints/linear/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "constraints/setppc/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "constraints/logicor/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetRealParam(subscip, "constraints/linear/maxaggrnormscale", 0.0) );

   /* disable conflict analysis */
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useinflp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useboundlp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) );
#endif

#ifdef SCIP_DEBUG
   /* for debugging RENS, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#endif

   success = FALSE;

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );

   /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
   SCIP_CALL( createSubproblem(scip, subscip, subvars, varmapfw, minfixingrate, binarybounds, uselprows, &success) );
   SCIPdebugMessage("RENS subproblem: %d vars, %d cons, success=%u\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip), success);
   
   /* free hash map */
   SCIPhashmapFree(&varmapfw);
   
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
      SCIP_Real upperbound;
      cutoff = SCIPinfinity(scip);
      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );   

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
      {
         cutoff = (1-minimprove)*SCIPgetUpperbound(scip) + minimprove*SCIPgetLowerbound(scip);
      }
      else
      {
         if ( SCIPgetUpperbound ( scip ) >= 0 )
            cutoff = ( 1 - minimprove ) * SCIPgetUpperbound ( scip );
         else
            cutoff = ( 1 + minimprove ) * SCIPgetUpperbound ( scip );
      }
      cutoff = MIN(upperbound, cutoff);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );
   }

   /* solve the subproblem */
   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
    */
#ifdef NDEBUG
   retstat = SCIPpresolve(subscip);
   if( retstat != SCIP_OKAY )
   { 
      SCIPwarningMessage("Error while presolving subMIP in RENS heuristic; subSCIP terminated with code <%d>\n",retstat);
   }
#else
   SCIP_CALL( SCIPpresolve(subscip) );
#endif



   SCIPdebugMessage("RENS presolved subproblem: %d vars, %d cons, success=%u\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip), success);

#if 0
   {
   char fname[SCIP_MAXSTRLEN];
   (void) SCIPsnprintf(fname, SCIP_MAXSTRLEN, "test/%s.lp",SCIPgetProbName(scip));
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", 1) );
   SCIP_CALL( SCIPsolve(subscip) );
   SCIP_CALL( SCIPwriteMIP(subscip,fname,TRUE,TRUE) );
   SCIPmessagePrintInfo("wrote RENS-subMIP to file <%s>\n",fname);
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", maxnodes ) );
   }
#endif

   /* after presolving, we should have at least reached a certain fixing rate over ALL variables (including continuous)
    * to ensure that not only the MIP but also the LP relaxation is easy enough
    */
   if( ( nvars - SCIPgetNVars(subscip) ) / (SCIP_Real)nvars >= minfixingrate / 2.0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;

      SCIPdebugMessage("solving subproblem: nstallnodes=%"SCIP_LONGINT_FORMAT", maxnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, maxnodes);

#ifdef NDEBUG
      retstat = SCIPsolve(subscip);
      if( retstat != SCIP_OKAY )
      { 
         SCIPwarningMessage("Error while solving subMIP in RENS heuristic; subSCIP terminated with code <%d>\n",retstat);
      }
#else
      SCIP_CALL( SCIPsolve(subscip) );
#endif


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
   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */

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
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMessage("skipping RENS: nstallnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

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

   SCIP_CALL( SCIPapplyRens(scip, heur, result, timelimit, memorylimit, heurdata->minfixingrate, heurdata-> minimprove,
         heurdata->maxnodes, nstallnodes, heurdata->binarybounds, heurdata->uselprows) );   

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
         "minimum percentage of integer variables that have to be fixable ",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );
   
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/rens/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
 
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/rens/binarybounds",
         "should general integers get binary bounds [floor(.),ceil(.)] ?",
         &heurdata->binarybounds, TRUE, DEFAULT_BINARYBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/rens/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/rens/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
 
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/rens/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/rens/minimprove",
         "factor by which RENS should at least improve the incumbent  ",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );
   
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/rens/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );
   
   return SCIP_OKAY;
}
