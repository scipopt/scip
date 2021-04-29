/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   benderscut_opt.c
 * @ingroup OTHER_CFILES
 * @brief  Generates a standard Benders' decomposition optimality cut
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/exprinterpret.h"
#include "nlpi/pub_expr.h"
#include "scip/benderscut_opt.h"
#include "scip/cons_linear.h"
#include "scip/pub_benderscut.h"
#include "scip/pub_benders.h"
#include "scip/pub_lp.h"
#include "scip/pub_nlp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_linear.h"
#include "scip/pub_var.h"
#include "scip/scip.h"
#include <string.h>

#define BENDERSCUT_NAME             "optimality"
#define BENDERSCUT_DESC             "Standard Benders' decomposition optimality cut"
#define BENDERSCUT_PRIORITY      5000
#define BENDERSCUT_LPCUT            TRUE

#define SCIP_DEFAULT_ADDCUTS             FALSE  /** Should cuts be generated, instead of constraints */

/*
 * Data structures
 */

/** Benders' decomposition cuts data */
struct SCIP_BenderscutData
{
   SCIP_Bool             addcuts;            /**< should cuts be generated instead of constraints */
};


/*
 * Local methods
 */

/** in the case of numerical troubles, the LP is resolved with solution polishing activated */
static
SCIP_RETCODE polishSolution(
   SCIP*                 subproblem,         /**< the SCIP data structure */
   SCIP_Bool*            success             /**< TRUE is the resolving of the LP was successful */
   )
{
   int oldpolishing;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   assert(subproblem != NULL);
   assert(SCIPinProbing(subproblem));

   (*success) = FALSE;

   /* setting the solution polishing parameter */
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/solutionpolishing", &oldpolishing) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/solutionpolishing", 2) );

   /* resolving the probing LP */
   SCIP_CALL( SCIPsolveProbingLP(subproblem, -1, &lperror, &cutoff) );

   if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_OPTIMAL )
      (*success) = TRUE;

   /* resetting the solution polishing parameter */
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/solutionpolishing", oldpolishing) );

   return SCIP_OKAY;
}

/** verifying the activity of the cut when master variables are within epsilon of their upper or lower bounds
 *
 *  When setting up the Benders' decomposition subproblem, master variables taking values that are within epsilon
 *  greater than their upper bound or less than their lower bound are set to their upper and lower bounds respectively.
 *  As such, there can be a difference between the subproblem dual solution objective and the optimality cut activity,
 *  when computed using the master problem solution directly. This check is to verify whether this difference is an
 *  actual error or due to the violation of the upper and lower bounds when setting up the Benders' decomposition
 *  subproblem.
 */
static
SCIP_RETCODE checkSetupTolerances(
   SCIP*                 masterprob,         /**< the SCIP data structure */
   SCIP_SOL*             sol,                /**< the master problem solution */
   SCIP_VAR**            vars,               /**< pointer to array of variables in the generated cut with non-zero coefficient */
   SCIP_Real*            vals,               /**< pointer to array of coefficients of the variables in the generated cut */
   SCIP_Real             lhs,                /**< the left hand side of the cut */
   SCIP_Real             checkobj,           /**< the objective of the subproblem computed from the dual solution */
   int                   nvars,              /**< the number of variables in the cut */
   SCIP_Bool*            valid               /**< returns true is the cut is valid */
   )
{
   SCIP_Real verifyobj;
   int i;

   assert(masterprob != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   /* initialising the verify objective with the left hand side of the optimality cut */
   verifyobj = lhs;

   /* computing the activity of the cut from the master solution and the constraint values */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real solval;

      solval = SCIPgetSolVal(masterprob, sol, vars[i]);

      /* checking whether the solution value is less than or greater than the variable bounds */
      if( !SCIPisLT(masterprob, solval, SCIPvarGetUbLocal(vars[i])) )
         solval = SCIPvarGetUbLocal(vars[i]);
      else if( !SCIPisGT(masterprob, solval, SCIPvarGetLbLocal(vars[i])) )
         solval = SCIPvarGetLbLocal(vars[i]);

      verifyobj -= solval*vals[i];
   }

   (*valid) = SCIPisFeasEQ(masterprob, checkobj, verifyobj);

   return SCIP_OKAY;
}

/** when solving NLP subproblems, numerical issues are addressed by tightening the feasibility tolerance */
static
SCIP_RETCODE resolveNLPWithTighterFeastol(
   SCIP*                 subproblem,         /**< the SCIP data structure */
   SCIP_Real             multiplier,         /**< the amount by which to decrease the tolerance */
   SCIP_Bool*            success             /**< TRUE is the resolving of the LP was successful */
   )
{
   SCIP_NLPSOLSTAT nlpsolstat;
#ifdef SCIP_DEBUG
   SCIP_NLPTERMSTAT nlptermstat;
#endif
#ifdef SCIP_MOREDEBUG
   SCIP_SOL* nlpsol;
#endif
   SCIP_Real feastol;
   SCIP_Real objtol;

   assert(subproblem != NULL);
   assert(SCIPinProbing(subproblem));

   (*success) = FALSE;

#ifdef SCIP_MOREDEBUG
   SCIP_CALL( SCIPsetNLPIntPar(subproblem, SCIP_NLPPAR_VERBLEVEL, 1) );
#endif

   SCIP_CALL( SCIPsetNLPIntPar(subproblem, SCIP_NLPPAR_ITLIM, INT_MAX) );

   /* getting the feasibility tolerance currently used for the NLP */
   SCIP_CALL( SCIPgetNLPRealPar(subproblem, SCIP_NLPPAR_FEASTOL, &feastol) );
   SCIP_CALL( SCIPgetNLPRealPar(subproblem, SCIP_NLPPAR_RELOBJTOL, &objtol) );

   /* setting the feasibility tolerance to 0.01x the current tolerance */
   SCIP_CALL( SCIPsetNLPRealPar(subproblem, SCIP_NLPPAR_FEASTOL, feastol*multiplier) );
   SCIP_CALL( SCIPsetNLPRealPar(subproblem, SCIP_NLPPAR_RELOBJTOL, objtol*multiplier) );

   SCIP_CALL( SCIPsolveNLP(subproblem) );

   nlpsolstat = SCIPgetNLPSolstat(subproblem);
#ifdef SCIP_DEBUG
   nlptermstat = SCIPgetNLPTermstat(subproblem);
   SCIPdebugMsg(subproblem, "NLP solstat %d termstat %d\n", nlpsolstat, nlptermstat);
#endif

   if( nlpsolstat == SCIP_NLPSOLSTAT_LOCOPT || nlpsolstat == SCIP_NLPSOLSTAT_GLOBOPT
      || nlpsolstat == SCIP_NLPSOLSTAT_FEASIBLE )
   {
#ifdef SCIP_MOREDEBUG
      SCIP_CALL( SCIPcreateNLPSol(subproblem, &nlpsol, NULL) );
      SCIP_CALL( SCIPprintSol(subproblem, nlpsol, NULL, FALSE) );
      SCIP_CALL( SCIPfreeSol(subproblem, &nlpsol) );
#endif

      (*success) = TRUE;
   }

   /* resetting the feasibility tolerance to 0.01x the current tolerance */
   SCIP_CALL( SCIPsetNLPRealPar(subproblem, SCIP_NLPPAR_FEASTOL, feastol) );

   return SCIP_OKAY;
}

/** adds a variable and value to the constraint/row arrays */
static
SCIP_RETCODE addVariableToArray(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_VAR***           vars,               /**< pointer to the array of variables in the generated cut with non-zero coefficient */
   SCIP_Real**           vals,               /**< pointer to the array of coefficients of the variables in the generated cut */
   SCIP_VAR*             addvar,             /**< the variable that will be added to the array */
   SCIP_Real             addval,             /**< the value that will be added to the array */
   int*                  nvars,              /**< the number of variables in the variable array */
   int*                  varssize            /**< the length of the variable size */
   )
{
   assert(masterprob != NULL);
   assert(vars != NULL);
   assert(*vars != NULL);
   assert(vals != NULL);
   assert(*vals != NULL);
   assert(addvar != NULL);
   assert(nvars != NULL);
   assert(varssize != NULL);

   if( *nvars >= *varssize )
   {
      *varssize = SCIPcalcMemGrowSize(masterprob, *varssize + 1);
      SCIP_CALL( SCIPreallocBufferArray(masterprob, vars, *varssize) );
      SCIP_CALL( SCIPreallocBufferArray(masterprob, vals, *varssize) );
   }
   assert(*nvars < *varssize);

   (*vars)[*nvars] = addvar;
   (*vals)[*nvars] = addval;
   (*nvars)++;

   return SCIP_OKAY;
}

/** returns the variable solution either from the NLP or from the primal vals array */
static
SCIP_Real getNlpVarSol(
   SCIP_VAR*             var,                /**< the variable for which the solution is requested */
   SCIP_Real*            primalvals,         /**< the primal solutions for the NLP, can be NULL */
   SCIP_HASHMAP*         var2idx             /**< mapping from variable of the subproblem to the index in the dual arrays, can be NULL */
   )
{
   SCIP_Real varsol;
   int idx;

   assert(var != NULL);
   assert((primalvals == NULL && var2idx == NULL) || (primalvals != NULL && var2idx != NULL));

   if( var2idx != NULL && primalvals != NULL )
   {
      assert(SCIPhashmapExists(var2idx, (void*)var) );
      idx = SCIPhashmapGetImageInt(var2idx, (void*)var);
      varsol = primalvals[idx];
   }
   else
      varsol = SCIPvarGetNLPSol(var);

   return varsol;
}

/** computes a standard Benders' optimality cut from the dual solutions of the LP */
static
SCIP_RETCODE computeStandardLPOptimalityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the subproblem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_VAR***           vars,               /**< pointer to array of variables in the generated cut with non-zero coefficient */
   SCIP_Real**           vals,               /**< pointer to array of coefficients of the variables in the generated cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   SCIP_Real*            rhs,                /**< the right hand side of the cut */
   int*                  nvars,              /**< the number of variables in the cut */
   int*                  varssize,           /**< the number of variables in the array */
   SCIP_Real*            checkobj,           /**< stores the objective function computed from the dual solution */
   SCIP_Bool*            success             /**< was the cut generation successful? */
   )
{
   SCIP_VAR** subvars;
   SCIP_VAR** fixedvars;
   int nsubvars;
   int nfixedvars;
   SCIP_Real dualsol;
   SCIP_Real addval;
   int nrows;
   int i;

   (*checkobj) = 0;

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(vars != NULL);
   assert(*vars != NULL);
   assert(vals != NULL);
   assert(*vals != NULL);

   (*success) = FALSE;

   /* looping over all LP rows and setting the coefficients of the cut */
   nrows = SCIPgetNLPRows(subproblem);
   for( i = 0; i < nrows; i++ )
   {
      SCIP_ROW* lprow;

      lprow = SCIPgetLPRows(subproblem)[i];
      assert(lprow != NULL);

      dualsol = SCIProwGetDualsol(lprow);
      assert( !SCIPisInfinity(subproblem, dualsol) && !SCIPisInfinity(subproblem, -dualsol) );

      if( SCIPisZero(subproblem, dualsol) )
         continue;

      if( dualsol > 0.0 )
         addval = dualsol*SCIProwGetLhs(lprow);
      else
         addval = dualsol*SCIProwGetRhs(lprow);

      (*lhs) += addval;

      /* if the bound becomes infinite, then the cut generation terminates. */
      if( SCIPisInfinity(masterprob, (*lhs)) || SCIPisInfinity(masterprob, -(*lhs))
         || SCIPisInfinity(masterprob, addval) || SCIPisInfinity(masterprob, -addval))
      {
         (*success) = FALSE;
         SCIPdebugMsg(masterprob, "Infinite bound when generating optimality cut. lhs = %g addval = %g.\n", (*lhs), addval);
         return SCIP_OKAY;
      }
   }

   nsubvars = SCIPgetNVars(subproblem);
   subvars = SCIPgetVars(subproblem);
   nfixedvars = SCIPgetNFixedVars(subproblem);
   fixedvars = SCIPgetFixedVars(subproblem);

   /* looping over all variables to update the coefficients in the computed cut. */
   for( i = 0; i < nsubvars + nfixedvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_VAR* mastervar;
      SCIP_Real redcost;

      if( i < nsubvars )
         var = subvars[i];
      else
         var = fixedvars[i - nsubvars];

      /* retrieving the master problem variable for the given subproblem variable. */
      SCIP_CALL( SCIPgetBendersMasterVar(masterprob, benders, var, &mastervar) );

      redcost = SCIPgetVarRedcost(subproblem, var);

      (*checkobj) += SCIPvarGetUnchangedObj(var)*SCIPvarGetSol(var, TRUE);

      /* checking whether the subproblem variable has a corresponding master variable. */
      if( mastervar != NULL )
      {
         SCIP_Real coef;

         coef = -1.0*(SCIPvarGetObj(var) + redcost);

         if( !SCIPisZero(masterprob, coef) )
         {
            /* adding the variable to the storage */
            SCIP_CALL( addVariableToArray(masterprob, vars, vals, mastervar, coef, nvars, varssize) );
         }
      }
      else
      {
         if( !SCIPisZero(subproblem, redcost) )
         {
            addval = 0;

            if( SCIPisPositive(subproblem, redcost) )
               addval = redcost*SCIPvarGetLbLocal(var);
            else if( SCIPisNegative(subproblem, redcost) )
               addval = redcost*SCIPvarGetUbLocal(var);

            (*lhs) += addval;

            /* if the bound becomes infinite, then the cut generation terminates. */
            if( SCIPisInfinity(masterprob, (*lhs)) || SCIPisInfinity(masterprob, -(*lhs))
               || SCIPisInfinity(masterprob, addval) || SCIPisInfinity(masterprob, -addval))
            {
               (*success) = FALSE;
               SCIPdebugMsg(masterprob, "Infinite bound when generating optimality cut.\n");
               return SCIP_OKAY;
            }
         }
      }
   }

   assert(SCIPisInfinity(masterprob, (*rhs)));
   /* the rhs should be infinite. If it changes, then there is an error */
   if( !SCIPisInfinity(masterprob, (*rhs)) )
   {
      (*success) = FALSE;
      SCIPdebugMsg(masterprob, "RHS is not infinite. rhs = %g.\n", (*rhs));
      return SCIP_OKAY;
   }

   (*success) = TRUE;

   return SCIP_OKAY;
}

/** computes a standard Benders' optimality cut from the dual solutions of the NLP */
static
SCIP_RETCODE computeStandardNLPOptimalityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the subproblem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_VAR***           vars,               /**< pointer to array of variables in the generated cut with non-zero coefficient */
   SCIP_Real**           vals,               /**< pointer to array of coefficients of the variables in the generated cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   SCIP_Real*            rhs,                /**< the right hand side of the cut */
   int*                  nvars,              /**< the number of variables in the cut */
   int*                  varssize,           /**< the number of variables in the array */
   SCIP_Real             objective,          /**< the objective function of the subproblem */
   SCIP_Real*            primalvals,         /**< the primal solutions for the NLP, can be NULL */
   SCIP_Real*            consdualvals,       /**< dual variables for the constraints, can be NULL */
   SCIP_Real*            varlbdualvals,      /**< the dual variables for the variable lower bounds, can be NULL */
   SCIP_Real*            varubdualvals,      /**< the dual variables for the variable upper bounds, can be NULL */
   SCIP_HASHMAP*         row2idx,            /**< mapping between the row in the subproblem to the index in the dual array, can be NULL */
   SCIP_HASHMAP*         var2idx,            /**< mapping from variable of the subproblem to the index in the dual arrays, can be NULL */
   SCIP_Real*            checkobj,           /**< stores the objective function computed from the dual solution */
   SCIP_Bool*            success             /**< was the cut generation successful? */
   )
{
   SCIP_EXPRINT* exprinterpreter;
   SCIP_VAR** subvars;
   SCIP_VAR** fixedvars;
   int nsubvars;
   int nfixedvars;
   SCIP_Real dirderiv;
   SCIP_Real dualsol;
   int nrows;
   int idx;
   int i;

   (*checkobj) = 0;

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(SCIPisNLPConstructed(subproblem));
   assert(SCIPgetNLPSolstat(subproblem) <= SCIP_NLPSOLSTAT_FEASIBLE || consdualvals != NULL);
   assert(SCIPhasNLPSolution(subproblem) || consdualvals != NULL);

   (*success) = FALSE;

   if( !(primalvals == NULL && consdualvals == NULL && varlbdualvals == NULL && varubdualvals == NULL && row2idx == NULL && var2idx == NULL)
      && !(primalvals != NULL && consdualvals != NULL && varlbdualvals != NULL && varubdualvals != NULL && row2idx != NULL && var2idx != NULL) ) /*lint !e845*/
   {
      SCIPerrorMessage("The optimality cut must generated from either a SCIP instance or all of the dual solutions and indices must be supplied");
      (*success) = FALSE;

      return SCIP_ERROR;
   }

   nsubvars = SCIPgetNNLPVars(subproblem);
   subvars = SCIPgetNLPVars(subproblem);
   nfixedvars = SCIPgetNFixedVars(subproblem);
   fixedvars = SCIPgetFixedVars(subproblem);

   /* our optimality cut implementation assumes that SCIP did not modify the objective function and sense,
    * that is, that the objective function value of the NLP corresponds to the value of the auxiliary variable
    * if that wouldn't be the case, then the scaling and offset may have to be considered when adding the
    * auxiliary variable to the cut (cons/row)?
    */
   assert(SCIPgetTransObjoffset(subproblem) == 0.0);
   assert(SCIPgetTransObjscale(subproblem) == 1.0);
   assert(SCIPgetObjsense(subproblem) == SCIP_OBJSENSE_MINIMIZE);

   (*lhs) = objective;
   assert(!SCIPisInfinity(subproblem, REALABS(*lhs)));

   (*rhs) = SCIPinfinity(masterprob);

   dirderiv = 0.0;

   SCIP_CALL( SCIPexprintCreate(SCIPblkmem(subproblem), &exprinterpreter) );

   /* looping over all NLP rows and setting the corresponding coefficients of the cut */
   nrows = SCIPgetNNLPNlRows(subproblem);
   for( i = 0; i < nrows; i++ )
   {
      SCIP_NLROW* nlrow;

      nlrow = SCIPgetNLPNlRows(subproblem)[i];
      assert(nlrow != NULL);

      if( row2idx != NULL && consdualvals != NULL )
      {
         assert(SCIPhashmapExists(row2idx, (void*)nlrow) );
         idx = SCIPhashmapGetImageInt(row2idx, (void*)nlrow);
         dualsol = consdualvals[idx];
      }
      else
         dualsol = SCIPnlrowGetDualsol(nlrow);
      assert( !SCIPisInfinity(subproblem, dualsol) && !SCIPisInfinity(subproblem, -dualsol) );

      if( SCIPisZero(subproblem, dualsol) )
         continue;

      SCIP_CALL( SCIPaddNlRowGradientBenderscutOpt(masterprob, subproblem, benders, nlrow, exprinterpreter,
            -dualsol, primalvals, var2idx, &dirderiv, vars, vals, nvars, varssize) );
   }

   SCIP_CALL( SCIPexprintFree(&exprinterpreter) );

   /* looping over all variable bounds and updating the corresponding coefficients of the cut; compute checkobj */
   for( i = 0; i < nsubvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_VAR* mastervar;
      SCIP_Real coef;

      var = subvars[i];

      (*checkobj) += SCIPvarGetObj(var) * getNlpVarSol(var, primalvals, var2idx);

      /* retrieving the master problem variable for the given subproblem variable. */
      SCIP_CALL( SCIPgetBendersMasterVar(masterprob, benders, var, &mastervar) );

      if( var2idx != NULL && varubdualvals != NULL && varlbdualvals != NULL )
      {
         assert(SCIPhashmapExists(var2idx, (void*)var) );
         idx = SCIPhashmapGetImageInt(var2idx, (void*)var);
         dualsol = varubdualvals[idx] - varlbdualvals[idx];
      }
      else
         dualsol = SCIPgetNLPVarsUbDualsol(subproblem)[i] - SCIPgetNLPVarsLbDualsol(subproblem)[i];

      /* checking whether the subproblem variable has a corresponding master variable. */
      if( mastervar == NULL || dualsol == 0.0 )
         continue;

      coef = -dualsol;

      /* adding the variable to the storage */
      SCIP_CALL( addVariableToArray(masterprob, vars, vals, mastervar, coef, nvars, varssize) );

      dirderiv += coef * getNlpVarSol(var, primalvals, var2idx);
   }

   for( i = 0; i < nfixedvars; i++ )
      *checkobj += SCIPvarGetUnchangedObj(fixedvars[i]) * getNlpVarSol(fixedvars[i], primalvals, var2idx);

   *lhs += dirderiv;

   /* if the side became infinite or dirderiv was infinite, then the cut generation terminates. */
   if( SCIPisInfinity(masterprob, *lhs) || SCIPisInfinity(masterprob, -*lhs)
      || SCIPisInfinity(masterprob, dirderiv) || SCIPisInfinity(masterprob, -dirderiv))
   {
      (*success) = FALSE;
      SCIPdebugMsg(masterprob, "Infinite bound when generating optimality cut. lhs = %g dirderiv = %g.\n", lhs, dirderiv);
      return SCIP_OKAY;
   }

   (*success) = TRUE;

   return SCIP_OKAY;
}


/** Adds the auxiliary variable to the generated cut. If this is the first optimality cut for the subproblem, then the
 *  auxiliary variable is first created and added to the master problem.
 */
static
SCIP_RETCODE addAuxiliaryVariableToCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_VAR**            vars,               /**< the variables in the generated cut with non-zero coefficient */
   SCIP_Real*            vals,               /**< the coefficients of the variables in the generated cut */
   int*                  nvars,              /**< the number of variables in the cut */
   int                   probnumber          /**< the number of the pricing problem */
   )
{
   SCIP_VAR* auxiliaryvar;

   assert(masterprob != NULL);
   assert(benders != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, probnumber);

   vars[(*nvars)] = auxiliaryvar;
   vals[(*nvars)] = 1.0;
   (*nvars)++;

   return SCIP_OKAY;
}


/*
 * Callback methods of Benders' decomposition cuts
 */

/** destructor of Benders' decomposition cuts to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BENDERSCUTFREE(benderscutFreeOpt)
{  /*lint --e{715}*/
   SCIP_BENDERSCUTDATA* benderscutdata;

   assert( benderscut != NULL );
   assert( strcmp(SCIPbenderscutGetName(benderscut), BENDERSCUT_NAME) == 0 );

   /* free Benders' cut data */
   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert( benderscutdata != NULL );

   SCIPfreeBlockMemory(scip, &benderscutdata);

   SCIPbenderscutSetData(benderscut, NULL);

   return SCIP_OKAY;
}


/** execution method of Benders' decomposition cuts */
static
SCIP_DECL_BENDERSCUTEXEC(benderscutExecOpt)
{  /*lint --e{715}*/
   SCIP* subproblem;
   SCIP_BENDERSCUTDATA* benderscutdata;
   SCIP_Bool nlprelaxation;
   SCIP_Bool addcut;
   char cutname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(benders != NULL);
   assert(benderscut != NULL);
   assert(result != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* retrieving the Benders' cut data */
   benderscutdata = SCIPbenderscutGetData(benderscut);

   /* if the cuts are generated prior to the solving stage, then rows can not be generated. So constraints must be
    * added to the master problem.
    */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE )
      addcut = FALSE;
   else
      addcut = benderscutdata->addcuts;

   /* setting the name of the generated cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "optimalitycut_%d_%" SCIP_LONGINT_FORMAT, probnumber,
      SCIPbenderscutGetNFound(benderscut) );

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   if( subproblem == NULL )
   {
      SCIPdebugMsg(scip, "The subproblem %d is set to NULL. The <%s> Benders' decomposition cut can not be executed.\n",
         probnumber, BENDERSCUT_NAME);

      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* setting a flag to indicate whether the NLP relaxation should be used to generate cuts */
   nlprelaxation = SCIPisNLPConstructed(subproblem) && SCIPgetNNlpis(subproblem);

   /* only generate optimality cuts if the subproblem LP or NLP is optimal,
    * since we use the dual solution of the LP/NLP to construct the optimality cut
    */
   if( SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING &&
      ((!nlprelaxation && SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_OPTIMAL) ||
       (nlprelaxation && SCIPgetNLPSolstat(subproblem) <= SCIP_NLPSOLSTAT_FEASIBLE)) )
   {
      /* generating a cut for a given subproblem */
      SCIP_CALL( SCIPgenerateAndApplyBendersOptCut(scip, subproblem, benders, benderscut, sol, probnumber, cutname,
            SCIPbendersGetSubproblemObjval(benders, probnumber), NULL, NULL, NULL, NULL, NULL, NULL, type, addcut,
            FALSE, result) );

      /* if it was not possible to generate a cut, this could be due to numerical issues. So the solution to the LP is
       * resolved and the generation of the cut is reattempted. For NLPs, we do not have such a polishing yet.
       */
      if( (*result) == SCIP_DIDNOTFIND )
      {
         SCIP_Bool success;

         SCIPdebugMsg(scip, "Numerical trouble generating optimality cut for subproblem %d.", probnumber);

         if( !nlprelaxation )
         {
            SCIPdebugMsg(scip, "Attempting to polish the LP solution to find an alternative dual extreme point.\n");

            SCIP_CALL( polishSolution(subproblem, &success) );

            /* only attempt to generate a cut if the solution polishing was successful */
            if( success )
            {
               SCIP_CALL( SCIPgenerateAndApplyBendersOptCut(scip, subproblem, benders, benderscut, sol, probnumber, cutname,
                     SCIPbendersGetSubproblemObjval(benders, probnumber), NULL, NULL, NULL, NULL, NULL, NULL, type, addcut,
                     FALSE, result) );
            }
         }
         else
         {
            SCIP_Real multiplier = 0.01;

            SCIPdebugMsg(scip, "Attempting to resolve the NLP with a tighter feasibility tolerance to find an "
               "alternative dual extreme point.\n");

            while( multiplier > 1e-06 && (*result) == SCIP_DIDNOTFIND )
            {
               SCIP_CALL( resolveNLPWithTighterFeastol(subproblem, multiplier, &success) );

               if( success )
               {
                  SCIP_CALL( SCIPgenerateAndApplyBendersOptCut(scip, subproblem, benders, benderscut, sol, probnumber, cutname,
                        SCIPbendersGetSubproblemObjval(benders, probnumber), NULL, NULL, NULL, NULL, NULL, NULL, type, addcut,
                        FALSE, result) );
               }

               multiplier *= 0.1;
            }
         }
      }
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
   char paramname[SCIP_MAXSTRLEN];

   assert(benders != NULL);

   /* create opt Benders' decomposition cuts data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &benderscutdata) );

   benderscut = NULL;

   /* include Benders' decomposition cuts */
   SCIP_CALL( SCIPincludeBenderscutBasic(scip, benders, &benderscut, BENDERSCUT_NAME, BENDERSCUT_DESC,
         BENDERSCUT_PRIORITY, BENDERSCUT_LPCUT, benderscutExecOpt, benderscutdata) );

   assert(benderscut != NULL);

   /* setting the non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBenderscutFree(scip, benderscut, benderscutFreeOpt) );

   /* add opt Benders' decomposition cuts parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/benderscut/%s/addcuts",
      SCIPbendersGetName(benders), BENDERSCUT_NAME);
   SCIP_CALL( SCIPaddBoolParam(scip, paramname,
         "should cuts be generated and added to the cutpool instead of global constraints directly added to the problem.",
         &benderscutdata->addcuts, FALSE, SCIP_DEFAULT_ADDCUTS, NULL, NULL) );

   return SCIP_OKAY;
}

/** Generates a classical Benders' optimality cut using the dual solutions from the subproblem or the input arrays. If
 *  the dual solutions are input as arrays, then a mapping between the array indices and the rows/variables is required.
 *  This method can also be used to generate a feasibility cut, if a problem to minimise the infeasibilities has been solved
 *  to generate the dual solutions
 */
SCIP_RETCODE SCIPgenerateAndApplyBendersOptCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< the benders' decomposition cut method */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the number of the pricing problem */
   char*                 cutname,            /**< the name for the cut to be generated */
   SCIP_Real             objective,          /**< the objective function of the subproblem */
   SCIP_Real*            primalvals,         /**< the primal solutions for the NLP, can be NULL */
   SCIP_Real*            consdualvals,       /**< dual variables for the constraints, can be NULL */
   SCIP_Real*            varlbdualvals,      /**< the dual variables for the variable lower bounds, can be NULL */
   SCIP_Real*            varubdualvals,      /**< the dual variables for the variable upper bounds, can be NULL */
   SCIP_HASHMAP*         row2idx,            /**< mapping between the row in the subproblem to the index in the dual array, can be NULL */
   SCIP_HASHMAP*         var2idx,            /**< mapping from variable of the subproblem to the index in the dual arrays, can be NULL */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_Bool             addcut,             /**< should the Benders' cut be added as a cut or constraint */
   SCIP_Bool             feasibilitycut,     /**< is this called for the generation of a feasibility cut */
   SCIP_RESULT*          result              /**< the result from solving the subproblems */
   )
{
   SCIP_CONSHDLR* consbenders;
   SCIP_CONS* cons;
   SCIP_ROW* row;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nvars;
   int varssize;
   int nmastervars;
   SCIP_Bool optimal;
   SCIP_Bool success;

   SCIP_Real checkobj;
   SCIP_Real verifyobj;

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(benderscut != NULL);
   assert(result != NULL);
   assert((primalvals == NULL && consdualvals == NULL && varlbdualvals == NULL && varubdualvals == NULL
         && row2idx == NULL && var2idx == NULL)
      || (primalvals != NULL && consdualvals != NULL && varlbdualvals != NULL && varubdualvals != NULL
         && row2idx != NULL && var2idx != NULL));

   row = NULL;
   cons = NULL;

   /* retrieving the Benders' decomposition constraint handler */
   consbenders = SCIPfindConshdlr(masterprob, "benders");

   /* checking the optimality of the original problem with a comparison between the auxiliary variable and the
    * objective value of the subproblem */
   if( feasibilitycut )
      optimal = FALSE;
   else
   {
      SCIP_CALL( SCIPcheckBendersSubproblemOptimality(masterprob, benders, sol, probnumber, &optimal) );
   }

   if( optimal )
   {
      (*result) = SCIP_FEASIBLE;
      SCIPdebugMsg(masterprob, "No cut added for subproblem %d\n", probnumber);
      return SCIP_OKAY;
   }

   /* allocating memory for the variable and values arrays */
   nmastervars = SCIPgetNVars(masterprob) + SCIPgetNFixedVars(masterprob);
   SCIP_CALL( SCIPallocClearBufferArray(masterprob, &vars, nmastervars) );
   SCIP_CALL( SCIPallocClearBufferArray(masterprob, &vals, nmastervars) );
   lhs = 0.0;
   rhs = SCIPinfinity(masterprob);
   nvars = 0;
   varssize = nmastervars;

   if( SCIPisNLPConstructed(subproblem) && SCIPgetNNlpis(subproblem) )
   {
      /* computing the coefficients of the optimality cut */
      SCIP_CALL( computeStandardNLPOptimalityCut(masterprob, subproblem, benders, &vars, &vals, &lhs, &rhs, &nvars,
            &varssize, objective, primalvals, consdualvals, varlbdualvals, varubdualvals, row2idx,
            var2idx, &checkobj, &success) );
   }
   else
   {
      /* computing the coefficients of the optimality cut */
      SCIP_CALL( computeStandardLPOptimalityCut(masterprob, subproblem, benders, &vars, &vals, &lhs, &rhs, &nvars,
            &varssize, &checkobj, &success) );
   }

   /* if success is FALSE, then there was an error in generating the optimality cut. No cut will be added to the master
    * problem. Otherwise, the constraint is added to the master problem.
    */
   if( !success )
   {
      (*result) = SCIP_DIDNOTFIND;
      SCIPdebugMsg(masterprob, "Error in generating Benders' optimality cut for problem %d.\n", probnumber);
   }
   else
   {
      /* creating an empty row or constraint for the Benders' cut */
      if( addcut )
      {
         SCIP_CALL( SCIPcreateEmptyRowConshdlr(masterprob, &row, consbenders, cutname, lhs, rhs, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddVarsToRow(masterprob, row, nvars, vars, vals) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsBasicLinear(masterprob, &cons, cutname, nvars, vars, vals, lhs, rhs) );
         SCIP_CALL( SCIPsetConsDynamic(masterprob, cons, TRUE) );
         SCIP_CALL( SCIPsetConsRemovable(masterprob, cons, TRUE) );
      }

      /* computing the objective function from the cut activity to verify the accuracy of the constraint */
      verifyobj = 0.0;
      if( addcut )
      {
         verifyobj += SCIProwGetLhs(row) - SCIPgetRowSolActivity(masterprob, row, sol);
      }
      else
      {
         verifyobj += SCIPgetLhsLinear(masterprob, cons) - SCIPgetActivityLinear(masterprob, cons, sol);
      }

      /* it is possible that numerics will cause the generated cut to be invalid. This cut should not be added to the
       * master problem, since its addition could cut off feasible solutions. The success flag is set of false, indicating
       * that the Benders' cut could not find a valid cut.
       */
      if( !feasibilitycut && !SCIPisFeasEQ(masterprob, checkobj, verifyobj) )
      {
         SCIP_Bool valid;

         /* the difference in the checkobj and verifyobj could be due to the setup tolerances. This is checked, and if
          * so, then the generated cut is still valid
          */
         SCIP_CALL( checkSetupTolerances(masterprob, sol, vars, vals, lhs, checkobj, nvars, &valid) );

         if( !valid )
         {
            success = FALSE;
            SCIPdebugMsg(masterprob, "The objective function and cut activity are not equal (%g != %g).\n", checkobj,
               verifyobj);

#ifdef SCIP_DEBUG
            /* we only need to abort if cut strengthen is not used. If cut strengthen has been used in this round and the
             * cut could not be generated, then another subproblem solving round will be executed
             */
            if( !SCIPbendersInStrengthenRound(benders) )
            {
#ifdef SCIP_MOREDEBUG
               int i;

               for( i = 0; i < nvars; i++ )
                  printf("<%s> %g %g\n", SCIPvarGetName(vars[i]), vals[i], SCIPgetSolVal(masterprob, sol, vars[i]));
#endif
               SCIPABORT();
            }
#endif
         }
      }

      if( success )
      {
         /* adding the auxiliary variable to the optimality cut */
         if( !feasibilitycut )
         {
            SCIP_CALL( addAuxiliaryVariableToCut(masterprob, benders, vars, vals, &nvars, probnumber) );
         }

         /* adding the constraint to the master problem */
         if( addcut )
         {
            SCIP_Bool infeasible;

            /* adding the auxiliary variable coefficient to the row */
            if( !feasibilitycut )
            {
               SCIP_CALL( SCIPaddVarToRow(masterprob, row, vars[nvars - 1], vals[nvars - 1]) );
            }

            if( type == SCIP_BENDERSENFOTYPE_LP || type == SCIP_BENDERSENFOTYPE_RELAX )
            {
               SCIP_CALL( SCIPaddRow(masterprob, row, FALSE, &infeasible) );
               assert(!infeasible);
            }
            else
            {
               assert(type == SCIP_BENDERSENFOTYPE_CHECK || type == SCIP_BENDERSENFOTYPE_PSEUDO);
               SCIP_CALL( SCIPaddPoolCut(masterprob, row) );
            }

            (*result) = SCIP_SEPARATED;
         }
         else
         {
            /* adding the auxiliary variable coefficient to the constraint */
            if( !feasibilitycut )
            {
               SCIP_CALL( SCIPaddCoefLinear(masterprob, cons, vars[nvars - 1], vals[nvars - 1]) );
            }

            SCIPdebugPrintCons(masterprob, cons, NULL);

            SCIP_CALL( SCIPaddCons(masterprob, cons) );

            (*result) = SCIP_CONSADDED;
         }

         /* storing the data that is used to create the cut */
         SCIP_CALL( SCIPstoreBendersCut(masterprob, benders, vars, vals, lhs, rhs, nvars) );
      }
      else
      {
         (*result) = SCIP_DIDNOTFIND;
         SCIPdebugMsg(masterprob, "Error in generating Benders' optimality cut for problem %d.\n", probnumber);
      }

      /* releasing the row or constraint */
      if( addcut )
      {
         /* release the row */
         SCIP_CALL( SCIPreleaseRow(masterprob, &row) );
      }
      else
      {
         /* release the constraint */
         SCIP_CALL( SCIPreleaseCons(masterprob, &cons) );
      }
   }

   SCIPfreeBufferArray(masterprob, &vals);
   SCIPfreeBufferArray(masterprob, &vars);

   return SCIP_OKAY;
}


/** adds the gradient of a nonlinear row in the current NLP solution of a subproblem to a linear row or constraint in the master problem
 *
 * Only computes gradient w.r.t. master problem variables.
 * Computes also the directional derivative, that is, mult times gradient times solution.
 */
SCIP_RETCODE SCIPaddNlRowGradientBenderscutOpt(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the subproblem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_EXPRINT*         exprint,            /**< expressions interpreter */
   SCIP_Real             mult,               /**< multiplier */
   SCIP_Real*            primalvals,         /**< the primal solutions for the NLP, can be NULL */
   SCIP_HASHMAP*         var2idx,            /**< mapping from variable of the subproblem to the index in the dual arrays, can be NULL */
   SCIP_Real*            dirderiv,           /**< storage to add directional derivative */
   SCIP_VAR***           vars,               /**< pointer to array of variables in the generated cut with non-zero coefficient */
   SCIP_Real**           vals,               /**< pointer to array of coefficients of the variables in the generated cut */
   int*                  nvars,              /**< the number of variables in the cut */
   int*                  varssize            /**< the number of variables in the array */
   )
{
   SCIP_EXPRTREE* tree;
   SCIP_VAR* var;
   SCIP_VAR* mastervar;
   SCIP_Real coef;
   int i;

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(nlrow != NULL);
   assert(exprint != NULL);
   assert((primalvals == NULL && var2idx == NULL) || (primalvals != NULL && var2idx != NULL));
   assert(mult != 0.0);
   assert(dirderiv != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   /* linear part */
   for( i = 0; i < SCIPnlrowGetNLinearVars(nlrow); i++ )
   {
      var = SCIPnlrowGetLinearVars(nlrow)[i];
      assert(var != NULL);

      /* retrieving the master problem variable for the given subproblem variable. */
      SCIP_CALL( SCIPgetBendersMasterVar(masterprob, benders, var, &mastervar) );
      if( mastervar == NULL )
         continue;

      coef = mult * SCIPnlrowGetLinearCoefs(nlrow)[i];

      /* adding the variable to the storage */
      SCIP_CALL( addVariableToArray(masterprob, vars, vals, mastervar, coef, nvars, varssize) );

      *dirderiv += coef * getNlpVarSol(var, primalvals, var2idx);
   }

   /* quadratic part */
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); i++ )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_VAR* mastervar1;
      SCIP_VAR* mastervar2;
      SCIP_Real coef1;
      SCIP_Real coef2;

      assert(SCIPnlrowGetQuadElems(nlrow)[i].idx1 < SCIPnlrowGetNQuadVars(nlrow));
      assert(SCIPnlrowGetQuadElems(nlrow)[i].idx2 < SCIPnlrowGetNQuadVars(nlrow));

      var1  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx1];
      var2  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx2];

      /* retrieving the master problem variables for the given subproblem variables. */
      SCIP_CALL( SCIPgetBendersMasterVar(masterprob, benders, var1, &mastervar1) );
      SCIP_CALL( SCIPgetBendersMasterVar(masterprob, benders, var2, &mastervar2) );

      coef1 = mult * SCIPnlrowGetQuadElems(nlrow)[i].coef * getNlpVarSol(var2, primalvals, var2idx);
      coef2 = mult * SCIPnlrowGetQuadElems(nlrow)[i].coef * getNlpVarSol(var1, primalvals, var2idx);

      /* adding the variable to the storage */
      if( mastervar1 != NULL )
      {
         SCIP_CALL( addVariableToArray(masterprob, vars, vals, mastervar1, coef1, nvars, varssize) );
      }
      if( mastervar2 != NULL )
      {
         SCIP_CALL( addVariableToArray(masterprob, vars, vals, mastervar2, coef2, nvars, varssize) );
      }

      if( mastervar1 != NULL )
         *dirderiv += coef1 * getNlpVarSol(var1, primalvals, var2idx);

      if( mastervar2 != NULL )
         *dirderiv += coef2 * getNlpVarSol(var2, primalvals, var2idx);
   }

   /* tree part */
   tree = SCIPnlrowGetExprtree(nlrow);
   if( tree != NULL )
   {
      SCIP_Real* treegrad;
      SCIP_Real* x;
      SCIP_Real val;

      SCIP_CALL( SCIPallocBufferArray(subproblem, &x, SCIPexprtreeGetNVars(tree)) );
      SCIP_CALL( SCIPallocBufferArray(subproblem, &treegrad, SCIPexprtreeGetNVars(tree)) );

      /* compile expression tree, if not done before */
      if( SCIPexprtreeGetInterpreterData(tree) == NULL )
      {
         SCIP_CALL( SCIPexprintCompile(exprint, tree) );
      }

      /* sets the solution value */
      for( i = 0; i < SCIPexprtreeGetNVars(tree); ++i )
         x[i] = getNlpVarSol(SCIPexprtreeGetVars(tree)[i], primalvals, var2idx);

      SCIP_CALL( SCIPexprintGrad(exprint, tree, x, TRUE, &val, treegrad) );

      /* update corresponding gradient entry */
      for( i = 0; i < SCIPexprtreeGetNVars(tree); ++i )
      {
         var = SCIPexprtreeGetVars(tree)[i];
         assert(var != NULL);

         /* retrieving the master problem variable for the given subproblem variable. */
         SCIP_CALL( SCIPgetBendersMasterVar(masterprob, benders, var, &mastervar) );
         if( mastervar == NULL )
            continue;

         coef = mult * treegrad[i];

         /* adding the variable to the storage */
         SCIP_CALL( addVariableToArray(masterprob, vars, vals, mastervar, coef, nvars, varssize) );

         *dirderiv += coef * getNlpVarSol(var, primalvals, var2idx);
      }

      SCIPfreeBufferArray(subproblem, &treegrad);
      SCIPfreeBufferArray(subproblem, &x);
   }

   return SCIP_OKAY;
}
