/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   benderscut_opt.c
 * @brief  Generates a standard Benders' decomposition optimality cut
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/benderscut_opt.h"
#include "scip/pub_benders.h"
#include "scip/pub_benderscut.h"
#include "scip/misc_benders.h"

#include "scip/cons_linear.h"


#define BENDERSCUT_NAME             "optimality"
#define BENDERSCUT_DESC             "Standard Benders' decomposition optimality cut"
#define BENDERSCUT_PRIORITY         0


/*
 * Data structures
 */

/* TODO: fill in the necessary compression data */

/** Benders' decomposition cuts data */
#if 0
struct SCIP_BenderscutData
{
};
#endif


/*
 * Local methods
 */

/* computing as standard Benders' optimality cut from the dual solutions of the LP */
static
SCIP_RETCODE computeStandardOptimalityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the subproblem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_CONS*            cut                 /**< the cut that is generated from the pricing problem */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** fixedvars;
   SCIP_CONS** conss;
   int nvars;
   int nfixedvars;
   int nconss;
   SCIP_Real dualsol;
   SCIP_Real addval;    /* the value that must be added to the lhs */
   SCIP_Real lhs;       /* the left hand side of the cut */
   int i;

#ifndef NDEBUG
   SCIP_Real verifyobj = 0;
   SCIP_Real checkobj = 0;
#endif

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(cut != NULL);

   nvars = SCIPgetNVars(subproblem);
   vars = SCIPgetVars(subproblem);
   nfixedvars = SCIPgetNFixedVars(subproblem);
   fixedvars = SCIPgetFixedVars(subproblem);

   nconss = SCIPgetNConss(subproblem);
   conss = SCIPgetConss(subproblem);


   /* looping over all constraints and setting the coefficients of the cut */
   for( i = 0; i < nconss; i++ )
   {
      dualsol = BDconsGetDualsol(subproblem, conss[i]);

      assert( !SCIPisInfinity(subproblem, dualsol) && !SCIPisInfinity(subproblem, -dualsol) );

      if( SCIPisZero(subproblem, dualsol) )
         continue;

      lhs = SCIPgetLhsLinear(masterprob, cut);


      if( SCIPisPositive(subproblem, dualsol) )
         addval = dualsol*BDconsGetLhs(subproblem, conss[i]);
      else if( SCIPisNegative(subproblem, dualsol) )
         addval = dualsol*BDconsGetRhs(subproblem, conss[i]);

      lhs += addval;

      /* Update the lhs of the cut */
      SCIP_CALL( SCIPchgLhsLinear(masterprob, cut, lhs) );
   }

   /* looping over all variables to update the coefficients in the computed cut. */
   for( i = 0; i < nvars + nfixedvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_VAR* mastervar;
      SCIP_Real redcost;


      if( i < nvars )
         var = vars[i];
      else
         var = fixedvars[i - nvars];

      /* retreiving the master problem variable for the given subproblem variable. */
      mastervar = SCIPgetBendersMasterVar(masterprob, benders, var);

      var = SCIPvarGetProbvar(var);

      redcost = SCIPgetVarRedcost(subproblem, var);

#ifndef NDEBUG
      checkobj += SCIPvarGetUnchangedObj(var)*SCIPvarGetSol(var, TRUE);
#endif

      /* checking whether the subproblem variable has a corresponding master variable. */
      if( mastervar != NULL )
      {
         SCIP_Real coef;

         //coef = -1.0*(SCIPvarGetObj(var) + redcost);
         coef = 1.0*(SCIPvarGetObj(var) + redcost);

#ifndef NDEBUG
         verifyobj -= SCIPvarGetObj(var)*SCIPvarGetSol(var, TRUE);
#endif

         SCIP_CALL( SCIPaddCoefLinear(masterprob, cut, mastervar, coef) );
      }
      else
      {
         if( !SCIPisZero(subproblem, redcost) )
         {
             addval = 0;

             /* get current lhs of the subproblem cut */
             lhs = SCIPgetLhsLinear(masterprob, cut);

             if( SCIPisPositive(subproblem, redcost) )
                addval = redcost*SCIPvarGetLbLocal(var);
             else if( SCIPisNegative(subproblem, redcost) )
                addval = redcost*SCIPvarGetUbLocal(var);

             lhs += addval;

             /* Update lhs */
             SCIP_CALL( SCIPchgLhsLinear(masterprob, cut, lhs) );
         }
      }
   }

#ifndef NDEBUG
   lhs = SCIPgetLhsLinear(masterprob, cut);
   verifyobj += lhs;

   /* need to generate a solution to verify the cut. Not sure how to do this yet. */
   //SCIP_CALL( GCGcreateSolFromGcgCol(subproblem, &sol, probnr, gcgcol) );
   //verifyobj -= SCIPgetActivityLinear(masterprob, cut, sol);
#endif

   //assert(SCIPisFeasEQ(origprob, SCIPgetLPObjval(origprob), verifyobj ));

   assert(cut != NULL);

   return SCIP_OKAY;
}


/* adds the auxiliary variable to the generated cut. If this is the first optimality cut for the subproblem, then the
 * auxiliary variable is first created and added to the master problem. */
static
SCIP_RETCODE addAuxiliaryVariableToCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_CONS*            cut,                /**< the cut that is generated from the pricing problem */
   int                   probnumber,         /**< the number of the pricing problem */
   SCIP_Bool*            optimal             /**< flag to indicate whether the current subproblem is optimal for the master */
   )
{
   SCIP_VAR* auxiliaryvar;
   SCIP_SOL* bestsol;               /* the current solution to the master problem */
   SCIP_Real auxiliaryvarval;

   assert(masterprob != NULL);
   assert(benders != NULL);
   assert(cut != NULL);

   (*optimal) = FALSE;

   auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, probnumber);

   /* retrieving the best solution to check the value of the auxiliary variable against the subproblem solution */
   bestsol = SCIPgetBestSol(masterprob);

   auxiliaryvarval = SCIPgetSolVal(masterprob, bestsol, auxiliaryvar);

   /* if the value of the auxiliary variable in the master problem is greater or equal to the subproblem objective,
    * then a cut is not added by the subproblem.
    * TODO: Need to use a epsilon tolerance for this check. */
   if( SCIPisGE(masterprob, auxiliaryvarval, SCIPbendersGetSubprobObjval(benders, probnumber)) )
      (*optimal) = TRUE;
   else
   {
      /* adding the auxiliary variable to the generated cut */
      SCIP_CALL( SCIPaddCoefLinear(masterprob, cut, auxiliaryvar, 1.0) );
   }

   return SCIP_OKAY;
}


/* generates and applies Benders' cuts */
static
SCIP_RETCODE generateAndApplyBendersCuts(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< the benders' decomposition cut method */
   int                   probnumber,         /**< the number of the pricing problem */
   SCIP_RESULT*          result              /**< the result from solving the subproblems */
   )
{
   SCIP_CONS* cut;                  /* the cut that will be generated from the solution to the pricing problem */
   char cutname[SCIP_MAXSTRLEN];    /* the name of the generated cut */
   SCIP_Bool optimal;               /* flag to indicate whether the current subproblem is optimal for the master */

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(result != NULL);
   assert(SCIPgetStatus(subproblem) == SCIP_STATUS_OPTIMAL);

   /* setting the name of the generated cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "optimalitycut_%d_%d", probnumber,
      SCIPbenderscutGetNFound(benderscut, probnumber) );

   /* creating the constraint for the cut */
   SCIP_CALL( SCIPcreateConsBasicLinear(masterprob, &cut, cutname, 0, NULL, NULL, 0.0, SCIPinfinity(masterprob)) );


   /* computing the coefficients of the optimality cut */
   SCIP_CALL( computeStandardOptimalityCut(masterprob, subproblem, benders, cut) );

   /* adding the auxiliary variable to the optimality cut */
   SCIP_CALL( addAuxiliaryVariableToCut(masterprob, benders, cut, probnumber, &optimal) );

   /* if the current subproblem is optimal for the master, then we do not add a constraint. */
   if( optimal )
   {
      SCIPinfoMessage(masterprob, NULL, "No cut added for subproblem %d\n", probnumber);
      SCIP_CALL( SCIPreleaseCons(masterprob, &cut) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPprintCons(masterprob, cut, NULL) );
   SCIPinfoMessage(masterprob, NULL, "\n");

   /* adding the constraint to the master problem */
   SCIP_CALL( SCIPaddCons(masterprob, cut) );

   SCIP_CALL( SCIPreleaseCons(masterprob, &cut) );

   (*result) = SCIP_CONSADDED;


   return SCIP_OKAY;
}
/*
 * Callback methods of Benders' decomposition cuts
 */

/* TODO: Implement all necessary Benders' decomposition cuts methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for Benders' decomposition cuts plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BENDERSCUTCOPY(benderscutCopyOpt)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of opt Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutCopyOpt NULL
#endif

/** destructor of Benders' decomposition cuts to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_BENDERSCUTFREE(benderscutFreeOpt)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of opt Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutFreeOpt NULL
#endif


/** initialization method of Benders' decomposition cuts (called after problem was transformed) */
#if 0
static
SCIP_DECL_BENDERSCUTINIT(benderscutInitOpt)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of opt Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutInitOpt NULL
#endif


/** deinitialization method of Benders' decomposition cuts (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BENDERSCUTEXIT(benderscutExitOpt)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of opt Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutExitOpt NULL
#endif


/** solving process initialization method of Benders' decomposition cuts (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BENDERSCUTINITSOL(benderscutInitsolOpt)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of opt Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutInitsolOpt NULL
#endif


/** solving process deinitialization method of Benders' decomposition cuts (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BENDERSCUTEXITSOL(benderscutExitsolOpt)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of opt Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutExitsolOpt NULL
#endif


/** execution method of Benders' decomposition cuts */
static
SCIP_DECL_BENDERSCUTEXEC(benderscutExecOpt)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(benders != NULL);
   assert(benderscut != NULL);
   assert(result != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* only generate optimality cuts if the subproblem is optimal */
   if( SCIPgetStatus(SCIPbendersSubproblem(benders, probnumber)) == SCIP_STATUS_OPTIMAL )
   {
      /* generating a cut for a given subproblem */
      SCIP_CALL( generateAndApplyBendersCuts(scip, SCIPbendersSubproblem(benders, probnumber), benders, benderscut,
            probnumber, result) );
   }

   return SCIP_OKAY;
}


/*
 * Benders' decomposition cuts specific interface methods
 */

/** creates the opt Benders' decomposition cuts and includes it in SCIP */
SCIP_RETCODE SCIPincludeBenderscutOpt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_BENDERSCUTDATA* benderscutdata;
   SCIP_BENDERSCUT* benderscut;

   assert(benders != NULL);

   /* create opt Benders' decomposition cuts data */
   benderscutdata = NULL;

   benderscut = NULL;

   /* include Benders' decomposition cuts */
#if 0
   /* use SCIPincludeBenderscut() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenderscut(scip, benders, BENDERSCUT_NAME, BENDERSCUT_DESC, BENDERSCUT_PRIORITY,
         benderscutCopyOpt, benderscutFreeXyz, benderscutInitXyz, benderscutExitXyz, benderscutInitsolXyz,
         benderscutExitsolOpt, benderscutExecXyz, benderscutdata) );
#else
   /* use SCIPincludeBenderscutBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenderscutBasic(scip, benders, &benderscut, BENDERSCUT_NAME, BENDERSCUT_DESC, BENDERSCUT_PRIORITY,
         benderscutExecOpt, benderscutdata) );

   assert(benderscut != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBenderscutCopy(scip, benderscut, benderscutCopyOpt) );
   SCIP_CALL( SCIPsetBenderscutFree(scip, benderscut, benderscutFreeOpt) );
   SCIP_CALL( SCIPsetBenderscutInit(scip, benderscut, benderscutInitOpt) );
   SCIP_CALL( SCIPsetBenderscutExit(scip, benderscut, benderscutExitOpt) );
   SCIP_CALL( SCIPsetBenderscutInitsol(scip, benderscut, benderscutInitsolOpt) );
   SCIP_CALL( SCIPsetBenderscutExitsol(scip, benderscut, benderscutExitsolOpt) );
#endif

   /* add opt Benders' decomposition cuts parameters */
   /* TODO: (optional) add Benders' decomposition cuts specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
