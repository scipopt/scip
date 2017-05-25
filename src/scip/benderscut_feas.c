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

/**@file   benderscut_feas.c
 * @brief  Standard feasibility cuts for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/benderscut_feas.h"
#include "scip/pub_benders.h"
#include "scip/pub_benderscut.h"
#include "scip/misc_benders.h"

#include "scip/cons_linear.h"


#define BENDERSCUT_NAME             "feas"
#define BENDERSCUT_DESC             "Standard feasibility cuts for Benders' decomposition"
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

/* computing as standard Benders' feasibility cut from the dual solutions of the LP */
/* NOTE: The cut must be created before being passed to this function */
static
SCIP_RETCODE computeStandardFeasibilityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,        /**< the SCIP instance of the pricing problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
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
   SCIP_Real lhs;       /* the left hand side of the cut */
   SCIP_Real addval;    /* the value that must be added to the lhs */
   int i;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_Real activity;
   SCIP_Real* farkascoefs;    // the coefficients of the farkas proof
   SCIP_Real farkasact = 0;   // the activities of the farkas proof
   SCIP_Real farkaslhs = 0;   // the lhs of the farkas proof
   int nconsvars;
   int j;

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

   SCIP_CALL( SCIPallocBufferArray(subproblem, &farkascoefs, nvars + nfixedvars) );
   for( i = 0; i < nvars + nfixedvars; i++ )
      farkascoefs[i] = 0;

   /* looping over all constraints and setting the coefficients of the cut */
   for( i = 0; i < nconss; i++ )
   {
      dualsol = BDconsGetDualfarkas(subproblem, conss[i]);

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

      farkaslhs += addval;

      nconsvars = BDconsGetNVars(subproblem, conss[i]);
      SCIP_CALL( SCIPallocBufferArray(subproblem, &consvars, nconsvars) );
      SCIP_CALL( SCIPallocBufferArray(subproblem, &consvals, nconsvars) );
      SCIP_CALL( BDconsGetVars(subproblem, conss[i], consvars, nconsvars) );
      SCIP_CALL( BDconsGetVals(subproblem, conss[i], consvals, nconsvars) );

      /* loop over all variables with non-zero coefficient */
      for( j = 0; j < nconsvars; j++ )
      {
         SCIP_VAR* mastervar;
         SCIP_VAR* consvar;
         SCIP_Real consval;

         consvar = consvars[j];
         consval = consvals[j];

         /* retreiving the master problem variable for the given subproblem variable. */
         mastervar = SCIPgetBendersMasterVar(masterprob, benders, consvars[j]);

         /* TODO: Do we need the problem variable? */
         consvar = SCIPvarGetProbvar(consvars[j]);

         //assert(!BDoriginalVarIsLinking(consvar));

         /* update the coefficient in the farkas activity */
         farkascoefs[SCIPvarGetProbindex(consvar)] += dualsol * consval;

         /* if the variable is a master variable, then it will be on the rhs of the constraint.
          * In computing the contribution of the fixed variables, we don't need to solution value because this is
          * given by the upper bound of the variable. */
         if( mastervar != NULL )
            farkaslhs -= dualsol * consval * SCIPvarGetUbLocal(consvar);
      }

      SCIPfreeBufferArray(subproblem, &consvars);
      SCIPfreeBufferArray(subproblem, &consvals);
   }

   /* looping over all variables to update the coefficients in the computed cut. */
   for( i = 0; i < nvars + nfixedvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_VAR* mastervar;

      if( i < nvars )
         var = vars[i];
      else
         var = fixedvars[i - nvars];

      /* retreiving the master problem variable for the given subproblem variable. */
      mastervar = SCIPgetBendersMasterVar(masterprob, benders, var);

      var = SCIPvarGetProbvar(var);

      //dualsol = farkascoefs[SCIPvarGetProbindex(var)];
      dualsol = SCIPgetVarFarkasCoef(subproblem, var);

      if( SCIPisZero(subproblem, dualsol) )
         continue;

      /* checking whether the original variable is a linking variable.
       * If this is the case, then the corresponding master variable is added to the generated cut.
       * If the pricing variable is not a linking variable, then the farkas dual value is added to the lhs */
      if( mastervar != NULL )
      {
         SCIP_CALL( SCIPaddCoefLinear(masterprob, cut, mastervar, dualsol) );
      }
      else
      {
         assert(SCIPisNegative(subproblem, dualsol));

         addval = 0;

         /* get current lhs of the subproblem cut */
         lhs = SCIPgetLhsLinear(masterprob, cut);

         if( SCIPisPositive(subproblem, dualsol) )
            addval = dualsol*SCIPvarGetUbLocal(var);
         else if( SCIPisNegative(subproblem, dualsol) )
            addval = dualsol*SCIPvarGetLbLocal(var);

         lhs += addval;


#ifndef NDEBUG
         farkasact += addval;
#endif

         /* Update lhs */
         SCIP_CALL( SCIPchgLhsLinear(masterprob, cut, lhs) );
      }
   }

#ifndef NDEBUG
   lhs = SCIPgetLhsLinear(masterprob, cut);
   activity = SCIPgetActivityLinear(masterprob, cut, sol);
   printf("LHS: %g Activity: %g\n", lhs, activity);
   assert(activity < lhs);
#endif


   assert(cut != NULL);


#ifndef NDEBUG
   /* TODO: Not sure about how to generate the solution for the first assert. Need to check */
   //assert(SCIPgetActivityLinear(masterprob, cut, pricingsol) < SCIPgetLhsLinear(masterprob, cut));
   printf("FARKAS - act: %g lhs: %g\n", farkasact, farkaslhs);
   assert(farkasact < farkaslhs);
   SCIPfreeBufferArray(subproblem, &farkascoefs);
#endif


   return SCIP_OKAY;
}


/* generates and applies Benders' cuts */
static
SCIP_RETCODE generateAndApplyBendersCuts(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< the benders' decomposition cut method */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the number of the pricing problem */
   SCIP_RESULT*          result              /**< the result from solving the subproblems */
   )
{
   SCIP_CONS* cut;                  /* the cut that will be generated from the solution to the pricing problem */
   char cutname[SCIP_MAXSTRLEN];    /* the name of the generated cut */

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(result != NULL);
   assert(SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE);

   /* setting the name of the generated cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "feasibilitycut_%d_%d", probnumber,
      SCIPbenderscutGetNFound(benderscut, probnumber) );

   /* creating the constraint for the cut */
   SCIP_CALL( SCIPcreateConsBasicLinear(masterprob, &cut, cutname, 0, NULL, NULL, 0.0, SCIPinfinity(masterprob)) );

   if( SCIPgetNLPIterations(subproblem) == 0 )
      SCIPinfoMessage(masterprob, NULL, "No iterations in pricing problem %d\n", probnumber);

   /* computing the coefficients of the feasibility cut */
   SCIP_CALL( computeStandardFeasibilityCut(masterprob, subproblem, benders, sol, cut) );

   //SCIP_CALL( SCIPprintCons(masterprob, cut, NULL) );
   //SCIPinfoMessage(masterprob, NULL, "\n");

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
SCIP_DECL_BENDERSCUTCOPY(benderscutCopyFeas)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of feas Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutCopyFeas NULL
#endif

/** destructor of Benders' decomposition cuts to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_BENDERSCUTFREE(benderscutFreeFeas)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of feas Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutFreeFeas NULL
#endif


/** initialization method of Benders' decomposition cuts (called after problem was transformed) */
#if 0
static
SCIP_DECL_BENDERSCUTINIT(benderscutInitFeas)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of feas Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutInitFeas NULL
#endif


/** deinitialization method of Benders' decomposition cuts (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BENDERSCUTEXIT(benderscutExitFeas)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of feas Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutExitFeas NULL
#endif


/** solving process initialization method of Benders' decomposition cuts (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BENDERSCUTINITSOL(benderscutInitsolFeas)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of feas Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutInitsolFeas NULL
#endif


/** solving process deinitialization method of Benders' decomposition cuts (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BENDERSCUTEXITSOL(benderscutExitsolFeas)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of feas Benders' decomposition cuts not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutExitsolFeas NULL
#endif


/** execution method of Benders' decomposition cuts */
static
SCIP_DECL_BENDERSCUTEXEC(benderscutExecFeas)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(benders != NULL);
   assert(benderscut != NULL);
   assert(result != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* only generate feasibility cuts if the subproblem is infeasible */
   if( SCIPgetStatus(SCIPbendersSubproblem(benders, probnumber)) == SCIP_STATUS_INFEASIBLE )
   {
      /* generating a cut for a given subproblem */
      SCIP_CALL( generateAndApplyBendersCuts(scip, SCIPbendersSubproblem(benders, probnumber), benders, benderscut,
            sol, probnumber, result) );
   }

   return SCIP_OKAY;
}


/*
 * Benders' decomposition cuts specific interface methods
 */

/** creates the Standard Feasibility Benders' decomposition cuts and includes it in SCIP */
SCIP_RETCODE SCIPincludeBenderscutFeas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_BENDERSCUTDATA* benderscutdata;
   SCIP_BENDERSCUT* benderscut;

   assert(benders != NULL);

   /* create feas Benders' decomposition cuts data */
   benderscutdata = NULL;

   benderscut = NULL;

   /* include Benders' decomposition cuts */
#if 0
   /* use SCIPincludeBenderscut() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenderscut(scip, benders, BENDERSCUT_NAME, BENDERSCUT_DESC, BENDERSCUT_PRIORITY,
         benderscutCopyFeas, benderscutFreeFeas, benderscutInitFeas, benderscutExitFeas, benderscutInitsolFeas,
         benderscutExitsolFeas, benderscutExecFeas, benderscutdata) );
#else
   /* use SCIPincludeBenderscutBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenderscutBasic(scip, benders, &benderscut, BENDERSCUT_NAME, BENDERSCUT_DESC, BENDERSCUT_PRIORITY,
         benderscutExecFeas, benderscutdata) );

   assert(benderscut != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBenderscutCopy(scip, benderscut, benderscutCopyFeas) );
   SCIP_CALL( SCIPsetBenderscutFree(scip, benderscut, benderscutFreeFeas) );
   SCIP_CALL( SCIPsetBenderscutInit(scip, benderscut, benderscutInitFeas) );
   SCIP_CALL( SCIPsetBenderscutExit(scip, benderscut, benderscutExitFeas) );
   SCIP_CALL( SCIPsetBenderscutInitsol(scip, benderscut, benderscutInitsolFeas) );
   SCIP_CALL( SCIPsetBenderscutExitsol(scip, benderscut, benderscutExitsolFeas) );
#endif

   /* add feas Benders' decomposition cuts parameters */
   /* TODO: (optional) add Benders' decomposition cuts specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
