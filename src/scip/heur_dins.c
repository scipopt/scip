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

/**@file   heur_dins.c
 * @brief  DINS primal heuristic (according to Ghosh)
 * @author Robert Waniek
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/heur_dins.h"

#define HEUR_NAME             "dins"
#define HEUR_DESC             "distance induced neighborhood search by Ghosh"
#define HEUR_DISPCHAR         'D'
#define HEUR_PRIORITY         -1105000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NODESOFS      5000LL    /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINNODES      500LL     /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which DINS should at least improve the incumbent          */
#define DEFAULT_NODESQUOT     0.05      /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_NWAITINGNODES 0LL       /* number of nodes without incumbent change that heuristic should wait */
#define DEFAULT_NEIGHBORHOODSIZE  18    /* radius of the incumbents neighborhood to be searched                */
#define DEFAULT_SOLNUM        5         /* number of pool-solutions to be checked for flag array update        */
#define DEFAULT_USELPROWS    FALSE      /* should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE      /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip        */


/*
 * Data structures
 */

/** DINS primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   SCIP_Real             minimprove;         /**< factor by which DINS should at least improve the incumbent          */
   SCIP_Longint          usednodes;          /**< nodes already used by DINS in earlier calls                         */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   int                   neighborhoodsize;   /**< radius of the incumbent's neighborhood to be searched               */
   SCIP_Bool*            delta;              /**< stores whether a variable kept its value from root LP all the time  */
   int                   deltalength;        /**< if there are no binary variables, we need no flag array             */
   SCIP_Longint          lastnsolsfound;     /**< solutions found until the last call of DINS                         */
   int                   solnum;             /**< number of pool-solutions to be checked for flag array update        */
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
   SCIP*                 scip,               /**< SCIP data structure of the original problem                    */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem                          */
   SCIP_VAR**            vars,               /**< variables of the original problem                              */
   SCIP_VAR**            subvars,            /**< variables of the subproblem                                    */
   int                   nbinvars,           /**< number of binary variables of problem and subproblem           */
   int                   nintvars,           /**< number of general integer variables of problem and subproblem  */
   SCIP_Bool             uselprows           /**< should subproblem be created out of the rows in the LP rows?   */
   )
{
   SCIP_SOL* bestsol;

   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(vars != NULL);
   assert(subvars != NULL);

   /* get the best MIP-solution known so far */
   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);

   /* create the rebounded general integer variables of the subproblem */
   for( i = nbinvars; i < nbinvars + nintvars; i++ )
   {
      SCIP_Real mipsol;
      SCIP_Real lpsol;

      SCIP_Real lbglobal;
      SCIP_Real ubglobal;

      /* get the bounds for each variable */
      lbglobal = SCIPvarGetLbGlobal(vars[i]);
      ubglobal = SCIPvarGetUbGlobal(vars[i]);

      assert(SCIPvarGetType(vars[i]) == SCIP_VARTYPE_INTEGER);
      /* get the current LP solution for each variable */
      lpsol = SCIPvarGetLPSol(vars[i]);
      /* get the current MIP solution for each variable */
      mipsol = SCIPgetSolVal(scip, bestsol, vars[i]);

      /* if the solution values differ by 0.5 or more, the variable is rebounded, otherwise it is just copied */
      if( REALABS(lpsol-mipsol) >= 0.5 )
      {
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real range;

         lb = lbglobal;
         ub = ubglobal;

         /* create a equally sized range around lpsol for general integers: bounds are lpsol +- (mipsol-lpsol) */
         range = 2*lpsol-mipsol;

         if( mipsol >= lpsol )
         {
            range = SCIPfeasCeil(scip, range);
            lb = MAX(lb, range);

            /* when the bound new upper bound is equal to the current MIP solution, we set both bounds to the integral bound (without eps) */
            if( SCIPisFeasEQ(scip, mipsol, lb) )
               ub = lb;
            else
               ub = mipsol;
         }
         else
         {
            range = SCIPfeasFloor(scip, range);
            ub = MIN(ub, range);

            /* when the bound new upper bound is equal to the current MIP solution, we set both bounds to the integral bound (without eps) */
            if( SCIPisFeasEQ(scip, mipsol, ub) )
               lb = ub;
            else
               lb = mipsol;
         }

         /* the global domain of variables might have been reduced since incumbent was found: adjust lb and ub accordingly */
         lb = MAX(lb, lbglobal);
         ub = MIN(ub, ubglobal);

         /* perform the bound change */
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], lb) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], ub) );
      }
      else
      {
         /* the global domain of variables might have been reduced since incumbent was found: adjust it accordingly */
         mipsol = MAX(mipsol, lbglobal);
         mipsol = MIN(mipsol, ubglobal);

         /* hard fixing for general integer variables with abs(mipsol-lpsol) < 0.5 */
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], mipsol) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], mipsol) );
      }
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

         assert(lhs <= rhs);

         /* allocate memory array to be filled with the corresponding subproblem variables */
         SCIP_CALL( SCIPallocBufferArray(subscip, &consvars, nnonz) );
         for( j = 0; j < nnonz; j++ )
            consvars[j] = subvars [ SCIPvarGetProbindex(SCIPcolGetVar(cols[j])) ];

         /* create a new linear constraint and add it to the subproblem */
         SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

         /* free temporary memory */
         SCIPfreeBufferArray(subscip, &consvars);
      }
   }

   return SCIP_OKAY;
}

/** create the extra constraint of local branching and add it to subscip */
static
SCIP_RETCODE addLocalBranchingConstraint(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem       */
   SCIP_VAR**            subvars,            /**< variables of the subproblem                 */
   SCIP_HEURDATA*        heurdata,           /**< heuristic's data structure                  */
   SCIP_Bool*            fixed               /**< TRUE --> include variable in LB constraint  */
   )
{
   SCIP_CONS* cons;                     /* local branching constraint to create          */
   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;

   SCIP_Real* consvals;
   SCIP_Real solval;
   SCIP_Real lhs;
   SCIP_Real rhs;

   char consname[SCIP_MAXSTRLEN];

   int nbinvars;
   int i;

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_dinsLBcons", SCIPgetProbName(scip));

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, NULL, NULL, NULL) );
   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nbinvars) );

   /* set initial left and right hand sides of local branching constraint */
   lhs = 0.0;
   rhs = (SCIP_Real) heurdata->neighborhoodsize;

   /* create the distance function of the binary variables (to incumbent solution) */
   for( i = 0; i < nbinvars; i++ )
   {
      consvars[i] = subvars[i];
      assert(SCIPvarGetType(consvars[i]) == SCIP_VARTYPE_BINARY);
      if( fixed[i] )
      {
         consvals[i]=0.0;
         continue;
      }

      solval = SCIPgetSolVal(scip, bestsol, vars[i]);
      assert(SCIPisFeasIntegral(scip, solval));

      /* is variable i part of the binary support of the current solution? */
      if( SCIPisFeasEQ(scip, solval, 1.0) )
      {
         consvals[i] = -1.0;
         rhs -= 1.0;
         lhs -= 1.0;
      }
      else
         consvals[i] = 1.0;
   }

   /* creates local branching constraint and adds it to subscip */
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
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< DINS heuristic structure                            */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(heur != NULL);
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
   if( *success )
   {
      SCIPdebugMessage("DINS successfully found new solution\n");
   }

   SCIPfreeBufferArray(scip, &subsolvals);
   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyDins)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurDins(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeDins)
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


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolDins)
{
   SCIP_HEURDATA* heurdata;
   int i;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->lastnsolsfound = 0;

   /* create flag array */
   heurdata->deltalength = SCIPgetNBinVars(scip);

   /* no binvars => no flag array needed */
   if( heurdata->deltalength > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(heurdata->delta), heurdata->deltalength) );
      for( i = 0; i < heurdata->deltalength; i++ )
         heurdata->delta[i] = TRUE;
   }
   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolDins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free flag array if exist */
   if( heurdata->deltalength > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &(heurdata->delta), heurdata->deltalength);
   }
   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecDins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP* subscip;                            /* the subproblem created by DINS                               */
   SCIP_VAR** subvars;                       /* subproblem's variables                                       */
   SCIP_VAR** vars;                          /* variables of the original problem                            */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_SOL* bestsol;                        /* best solution known so far                                   */
   SCIP_SOL** sols;                          /* list of known solutions                                      */

   SCIP_Bool* fixed;                         /* fixing flag array                                            */
   SCIP_Bool* delta;                         /* flag array if variable value changed during solution process */


   SCIP_Longint maxnnodes;                   /* maximum number of subnodes                                   */
   SCIP_Longint nsubnodes;                   /* nodelimit for subscip                                        */
   SCIP_Longint nsolsfound;

   SCIP_Real timelimit;                      /* timelimit for subscip (equals remaining time of scip)        */
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem                          */
   SCIP_Real upperbound;
   SCIP_Real memorylimit;                    /* memory limit for solution process of subscip                 */
   SCIP_Real lpsolval;
   SCIP_Real rootlpsolval;
   SCIP_Real mipsolval;
   SCIP_Real solval;

   int ufcount;                              /* counts the number of true fixing flag entries                */
   int nvars;                                /* number of variables in original SCIP                         */
   int nbinvars;                             /* number of binary variables in original SCIP                  */
   int nintvars;                             /* number of general integer variables in original SCIP         */
   int nsols;                                /* number of known solutions                                    */
   int nsubsols;
   int checklength;
   int i;
   int j;

   SCIP_Bool success;                        /* used to store whether new solution was found or not          */
   SCIP_Bool infeasible;                     /* stores whether the hard fixing of a variables was feasible or not */

   SCIP_RETCODE retcode;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DELAYED;

   /* only call heuristic, if a CIP solution is at hand */
   if( SCIPgetNSols(scip) <= 0 )
      return SCIP_OKAY;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   delta = heurdata->delta;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip, SCIPgetBestSol(scip)) < heurdata->nwaitingnodes )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* determine the node limit for the current process */
   maxnnodes = (SCIP_Longint) (heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward DINS if it succeeded often */
   maxnnodes = (SCIP_Longint) (maxnnodes * (1.0 + 2.0 * (SCIPheurGetNBestSolsFound(heur)+1.0) / (SCIPheurGetNCalls(heur) + 1.0)));

   /* count the setup costs for the sub-MIP as 100 nodes */
   maxnnodes -= 100 * SCIPheurGetNCalls(heur);
   maxnnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nsubnodes = maxnnodes - heurdata->usednodes;
   nsubnodes = MIN(nsubnodes , heurdata->maxnodes);

   /* check whether we have enough nodes left to call sub problem solving */
   if( nsubnodes < heurdata->minnodes )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
     return SCIP_OKAY;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   assert(nbinvars <= nvars);

   /* do not run heuristic if only continuous variables are present */
   if( nbinvars == 0 && nintvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

   success = FALSE;
   if( heurdata->uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

      /* get name of the original problem and add the string "_dinssub" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_dinssub", SCIPgetProbName(scip));

      /* create the subproblem */
      SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );
   }
   else
   {
      SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "dins", TRUE, FALSE, TRUE, &success) );

      if( heurdata->copycuts )
      {
         /** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, TRUE, NULL) );
      }

      SCIPdebugMessage("Copying the SCIP instance was %ssuccessful.\n", success ? "" : "not ");
   }

   for( i = 0; i < nvars; i++ )
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* create variables and rebound them if their bounds differ by more than 0.5 */
   SCIP_CALL( createSubproblem(scip, subscip, vars, subvars, nbinvars, nintvars, heurdata->uselprows) );
   SCIPdebugMessage("DINS subproblem: %d vars (%d binvars & %d intvars), %d cons\n",
      SCIPgetNVars(subscip), SCIPgetNBinVars(subscip) , SCIPgetNIntVars(subscip) , SCIPgetNConss(subscip));

   *result = SCIP_DIDNOTFIND;

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

   /* get the best MIP-solution known so far */
   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);

   /* get solution pool and number of solutions in pool */
   sols = SCIPgetSols(scip);
   nsols = SCIPgetNSols(scip);
   nsolsfound = SCIPgetNSolsFound(scip);
   checklength = MIN(nsols, heurdata->solnum);
   assert(sols != NULL);
   assert(nsols > 0);

   /* create fixing flag array */
   SCIP_CALL( SCIPallocBufferArray(scip, &fixed, nbinvars) );

   /* if new binary variables have been created, e.g., due to column generation, reallocate the delta array */
   if( heurdata->deltalength < nbinvars )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, nbinvars);
      assert(newsize >= nbinvars);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &heurdata->delta, heurdata->deltalength, newsize) );
      delta = heurdata->delta;

      /* initialize new part of delta array */
      for( i = heurdata->deltalength; i < newsize; i++ )
         delta[i] = TRUE;

      heurdata->deltalength = newsize;
   }

   /* fixing for binary variables */
   /* hard fixing for some with mipsol(s)=lpsolval=rootlpsolval and preparation for soft fixing for the remaining */
   ufcount = 0;
   for( i = 0; i < nbinvars; i++ )
   {
      /* soft fixing if the variable somewhen changed its value or the relaxations differ by adding a local branching constraint */
      fixed[i] = FALSE;

      /* get the current LP solution for each variable */
      lpsolval = SCIPvarGetLPSol(vars[i]);
      /* get the current MIP solution for each variable */
      mipsolval = SCIPgetSolVal(scip, bestsol, vars[i]);
      /* get the root LP solution for each variable */
      rootlpsolval = SCIPvarGetRootSol(vars[i]);

      if( SCIPisFeasEQ(scip, lpsolval, mipsolval) && SCIPisFeasEQ(scip, mipsolval, rootlpsolval) )
      {
         /* update delta */
         if( nsols > 1 && heurdata->lastnsolsfound != nsolsfound && delta[i] ) /* no need to update delta[i] if already FALSE */
         {
            /* no need to update delta[i] if already FALSE or sols[i] already checked on previous run or worse than DINS-solution of last run */
            for( j = 0; delta[i] && j < checklength && SCIPgetSolHeur(scip, sols[j]) != heur ; j++ )
            {
               solval = SCIPgetSolVal(scip, sols[j], vars[i]);
               delta[i] = delta[i] && SCIPisFeasEQ(scip, mipsolval, solval);
            }
         }

         /* hard fixing if rootlpsolval=nodelpsolval=mipsolval(s) and delta (is TRUE) */
         if( delta[i] && SCIPisFeasEQ(scip, mipsolval, lpsolval) && SCIPisFeasEQ(scip, mipsolval, rootlpsolval)
            && SCIPisFeasEQ(scip, rootlpsolval, lpsolval) )
         {
            SCIP_CALL( SCIPfixVar(subscip, subvars[i], mipsolval, &infeasible, &success) );
            fixed[i] = !infeasible;
            if( !success )
            {
               SCIPdebugMessage("variable %d was already fixed\n", i);
            }
            if( infeasible )
            {
               SCIPdebugMessage("fixing of variable %d to value %f was infeasible\n", i, mipsolval);
            }
         }
      }
      if( !fixed[i] )
         ufcount++;
   }

   /* store the number of found solutions for next run */
   heurdata->lastnsolsfound = nsolsfound;

   /* perform prepared softfixing for all unfixed vars if the number of unfixed vars is larger than the neighborhoodsize (otherwise it will be useless) */
   if( ufcount > heurdata->neighborhoodsize )
   {
      SCIP_CALL( addLocalBranchingConstraint(scip, subscip, subvars, heurdata, fixed) );
   }

   /* free fixing flag array */
   SCIPfreeBufferArray(scip, &fixed);

   /* add an objective cutoff */
   cutoff = SCIPinfinity(scip);
   assert(!SCIPisInfinity(scip, SCIPgetUpperbound(scip)));

   if( !SCIPisInfinity(scip, -1.0*SCIPgetLowerbound(scip)) )
   {
      cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip) + heurdata->minimprove * SCIPgetLowerbound(scip);
      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      cutoff = MIN(upperbound, cutoff);
   }
   else
   {
      if( SCIPgetUpperbound(scip) >= 0 )
         cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip);
      else
         cutoff = (1 + heurdata->minimprove) * SCIPgetUpperbound(scip);
      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      cutoff = MIN(upperbound, cutoff);
   }
   SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );

   /* solve the subproblem */
   SCIPdebugMessage("solving DINS sub-MIP with neighborhoodsize %d and maxnodes %"SCIP_LONGINT_FORMAT"\n", heurdata->neighborhoodsize, nsubnodes);
   retcode = SCIPsolve(subscip);

   /* Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in DINS heuristic; sub-SCIP terminated with code <%d>\n", retcode);
   }

   heurdata->usednodes += SCIPgetNNodes(subscip);
   nsubsols = SCIPgetNSols(subscip);
   SCIPdebugMessage("DINS used %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT" nodes and found %d solutions\n", SCIPgetNNodes(subscip), nsubnodes, nsubsols);

   /* check, whether a  (new) solution was found */
   if( nsubsols > 0 )
   {
      SCIP_SOL** subsols;

      /* check, whether a solution was found; due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted */
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

/** creates the DINS primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurDins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Dins primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecDins, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyDins) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeDins) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolDins) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolDins) );

   /* add DINS primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, FALSE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/solnum",
         "number of pool-solutions to be checked for flag array update (for hard fixing of binary variables)",
         &heurdata->solnum, FALSE, DEFAULT_SOLNUM, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/neighborhoodsize",
         "radius (using Manhattan metric) of the incumbent's neighborhood to be searched",
         &heurdata->neighborhoodsize, FALSE, DEFAULT_NEIGHBORHOODSIZE, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which "HEUR_NAME" should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
