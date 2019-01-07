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

/**@file   cons_expr_nlhdlr_perspective.c
 * @brief  perspective nonlinear handler
 * @author Ksenia Bestuzheva
 */


#include <string.h>

#include "scip/cons_varbound.h"
#include "scip/cons_expr_nlhdlr_perspective.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_sum.h"
#include "scip/scip_sol.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "perspective"
#define NLHDLR_DESC         "perspective handler for expressions"
#define NLHDLR_PRIORITY     75

/*
 * Data structures
 */

/** data structure to store information of a semicontinuous variable */
struct SCIP_SCVarData
{
   SCIP_VAR*             bvar;               /**< the binary variable on which the variable domain depends */
   SCIP_Real             val0;               /**< var value when bvar = 0 */
};
typedef struct SCIP_SCVarData SCIP_SCVARDATA;

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_EXPRCURV         curvature;          /**< curvature of the expression */
   SCIP_CONSEXPR_EXPR**  onoffterms;         /**< on/off terms for which we apply perspective cuts */
   SCIP_Real*            onoffcoefs;         /**< coefficients of onoffterms */
   SCIP_VAR**            bvars;              /**< binary vars associated with onoffterms */
   SCIP_CONSEXPR_EXPR**  convterms;          /**< convex terms for which we apply gradient cuts */
   SCIP_Real*            convcoefs;          /**< coefficients of convterms */
   int                   nonoffterms;        /**< number of on/off expressions */
   int                   nconvterms;         /**< number of convterms */
   int                   convtermssize;      /**< size of the convterms array */
   SCIP_CONSEXPR_EXPR**  varexprs;           /**< variable expressions */
   int                   nvarexprs;          /**< total number of variable expressions */
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_HASHMAP*         scvars;             /**< maps semicontinuous variables to their on/off bounds */
};

/*
 * Local methods
 */

/** checks if a variable is semicontinuous and, if needed, updates the hashmap
 *
 * A variable is semicontinuous if its bounds depend on the binary variable bvar and bvar == 0 => var = v_off for some
 * real constant v_off. If the bvar is not specified, find the first binary variable that var depends on.
 */
static
SCIP_RETCODE varIsSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< the variable to check */
   SCIP_VAR**            bvar,               /**< the binary variable */
   SCIP_HASHMAP*         scvars,             /**< semicontinuous variable information */
   SCIP_Bool*            result              /**< buffer to store whether var is semicontinuous */
   )
{
   SCIP_SCVARDATA* scv;
   SCIP_CONSHDLR* varboundconshdlr;
   SCIP_CONS** vbconss;
   SCIP_SCVARDATA* olddata;
   int nvbconss, c;
   SCIP_Real lb0, ub0, lb1, ub1;

   assert(scip != NULL);
   assert(var != NULL);
   assert(scvars != NULL);
   assert(result != NULL);

   *result = FALSE;

   olddata = (SCIP_SCVARDATA*) SCIPhashmapGetImage(scvars, (void*)var);
   if( olddata != NULL )
   {
      if( *bvar == NULL )
         *bvar = olddata->bvar;
      *result = (olddata->bvar == *bvar);
      return SCIP_OKAY;
   }

   lb0 = -SCIPinfinity(scip);
   ub0 = SCIPinfinity(scip);
   lb1 = -SCIPinfinity(scip);
   ub1 = SCIPinfinity(scip);

   varboundconshdlr = SCIPfindConshdlr(scip, "varbound");
   /* TODO @bzfkensi in some unittests, we don't have the varbound conshdlr; so is it correct to return FALSE here? */
   if( varboundconshdlr == NULL )
      return SCIP_OKAY;

   nvbconss = SCIPconshdlrGetNConss(varboundconshdlr);
   vbconss = SCIPconshdlrGetConss(varboundconshdlr);

   for( c = 0; c < nvbconss; ++c )
   {
      /* look for varbounds containing var and, if specified, bvar */
      if( SCIPgetVarVarbound(scip, vbconss[c]) != var || (*bvar != NULL && SCIPgetVbdvarVarbound(scip, vbconss[c]) != *bvar) )
         continue;

      if( *bvar == NULL )
         *bvar = SCIPgetVbdvarVarbound(scip, vbconss[c]);

      /* try to update the bound information */
      if( SCIPgetLhsVarbound(scip, vbconss[c]) > lb0 )
         lb0 = SCIPgetLhsVarbound(scip, vbconss[c]);
      if( SCIPgetRhsVarbound(scip, vbconss[c]) < ub0 )
         ub0 = SCIPgetRhsVarbound(scip, vbconss[c]);
      if( SCIPgetLhsVarbound(scip, vbconss[c]) - SCIPgetVbdcoefVarbound(scip, vbconss[c]) > lb1 )
         lb1 = SCIPgetLhsVarbound(scip, vbconss[c]) - SCIPgetVbdcoefVarbound(scip, vbconss[c]);
      if( SCIPgetRhsVarbound(scip, vbconss[c]) - SCIPgetVbdcoefVarbound(scip, vbconss[c]) > ub1 )
         ub1 = SCIPgetRhsVarbound(scip, vbconss[c]) - SCIPgetVbdcoefVarbound(scip, vbconss[c]);
   }

   /* if any bound is dominated by the global bound, replace it */
   if( SCIPvarGetLbGlobal(var) >= lb0 )
      lb0 = SCIPvarGetLbGlobal(var);
   if( SCIPvarGetUbGlobal(var) <= ub0 )
      ub0 = SCIPvarGetUbGlobal(var);
   if( SCIPvarGetLbGlobal(var) >= lb1 )
      lb1 = SCIPvarGetLbGlobal(var);
   if( SCIPvarGetUbGlobal(var) <= ub1 )
      ub1 = SCIPvarGetUbGlobal(var);

   SCIPdebugMsg(scip, "onoff bounds: lb0 = %f, ub0 = %f, lb1 = %f, ub1 = %f\n", lb0, ub0, lb1, ub1);

   /* the 'off' domain of the variable should reduce to a single point and be different from the 'on' domain */
   if( lb0 == ub0 && (lb0 != lb1 || ub0 != ub1) )  /*lint !e777*/
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &scv) );
      scv->bvar = *bvar;
      scv->val0 = lb0;
      if( lb0 != 0.0 )
      {
         SCIPdebugMsg(scip, "var %s has a non-zero off value %f\n", SCIPvarGetName(var), lb0);
      }
      SCIP_CALL( SCIPhashmapInsert(scvars, var, scv) );
      *result = TRUE;
   }

   return SCIP_OKAY;
}

/** adds an expression to the hashmap linking binary vars to on/off expressions */
static
SCIP_RETCODE addOnoffTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_HASHMAP*         onoffterms,         /**< hashmap linking binary vars to on/off terms */
   int*                  nonoffterms,        /**< number of on/off terms */
   SCIP_Real             coef,               /**< coef of the added term */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expr to add */
   SCIP_VAR*             bvar                /**< the binary variable */
)
{
   SCIP_EXPRCURV ecurv;
   SCIP_EXPRCURV pcurv;
   SCIP_CONSEXPR_EXPR* persp;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(onoffterms != NULL);
   assert(expr != NULL);
   assert(bvar != NULL);

   ecurv = SCIPexprcurvMultiply(coef, SCIPgetConsExprExprCurvature(expr));
   if( SCIPhashmapExists(onoffterms, (void*)bvar) )
   {
      persp = (SCIP_CONSEXPR_EXPR*) SCIPhashmapGetImage(onoffterms, (void*)bvar);
      pcurv = SCIPgetConsExprExprCurvature(persp);
      if( pcurv == SCIP_EXPRCURV_LINEAR && (ecurv == SCIP_EXPRCURV_CONVEX || ecurv == SCIP_EXPRCURV_CONCAVE) )
         SCIPsetConsExprExprCurvature(persp, ecurv);
      else if( ecurv == SCIP_EXPRCURV_UNKNOWN || (ecurv != SCIP_EXPRCURV_LINEAR && ecurv != pcurv) )
         SCIPsetConsExprExprCurvature(persp, SCIP_EXPRCURV_UNKNOWN);
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, persp, expr, coef) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &persp, 1, &expr, &coef, 0.0) );
      SCIPsetConsExprExprCurvature(persp, ecurv);
      SCIP_CALL( SCIPhashmapInsert(onoffterms, (void*)bvar, (void*)persp) );
   }
   ++*nonoffterms;

   return SCIP_OKAY;
}

/** adds an expression to the array of convex expressions */
static
SCIP_RETCODE addConvTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nonlinear handler expression data */
   SCIP_Real             coef,               /**< coefficient of expr in the original expression */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression to be added */
)
{
   int newsize;

   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(expr != NULL);

   if( nlhdlrexprdata->nconvterms + 1 > nlhdlrexprdata->convtermssize )
   {
      newsize = SCIPcalcMemGrowSize(scip, nlhdlrexprdata->nconvterms + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->convterms,  nlhdlrexprdata->convtermssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->convcoefs, nlhdlrexprdata->convtermssize, newsize) );
      nlhdlrexprdata->convtermssize = newsize;
   }
   assert(nlhdlrexprdata->nconvterms + 1 <= nlhdlrexprdata->convtermssize);

   nlhdlrexprdata->convcoefs[nlhdlrexprdata->nconvterms] = coef;
   nlhdlrexprdata->convterms[nlhdlrexprdata->nconvterms] = expr;
   nlhdlrexprdata->nconvterms++;

   return SCIP_OKAY;
}

/** constructs gradient linearization of a given expression and adds it to rowprep */
static
SCIP_RETCODE addGradientLinearisation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_ROWPREP*         rowprep,            /**< a rowprep where the linearization is stored */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be linearized */
   SCIP_Real             coef,               /**< coefficient of expr in the original expression */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_Bool*            success             /**< indicates whether the linearization could be computed */
)
{
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_Real constant;
   int i, v, nvars;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(rowprep != NULL);
   assert(expr != NULL);
   assert(success != NULL);

   /* compute gradient */
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );

   /* gradient evaluation error -> skip */
   if( SCIPgetConsExprExprDerivative(expr) == SCIP_INVALID ) /*lint !e777*/
   {
      *success = FALSE;
      SCIPdebugMsg(scip, "gradient evaluation error for %p\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* get g(x*) */
   constant = SCIPgetConsExprExprValue(expr);

   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", constant, (void*)expr);
      return SCIP_OKAY;
   }

   /* compute gradient cut */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip)) );
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, varexprs, &nvars) );
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real derivative;
      SCIP_Real val;

      assert(varexprs[i] != NULL);
      assert(SCIPisConsExprExprVar(varexprs[i]));

      /* get the variable of the variable expression */
      var = SCIPgetConsExprExprVarVar(varexprs[i]);
      assert(var != NULL);

      /* get solution value */
      val = SCIPgetSolVal(scip, sol, var);

      /* avoid overhead of SCIPgetConsExprExprPartialDiff by accessing the derivative directly */
      derivative = SCIPgetConsExprExprDerivative(varexprs[i]);
      assert(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, var) == derivative); /*lint !e777*/

      /* evaluation error or too large values -> skip */
      if( SCIPisInfinity(scip, REALABS(derivative * val)) )
      {
         *success = FALSE;
         SCIPdebugMsg(scip, "evaluation error / too large values (%g %g) for %s in %p\n", derivative, val,
                      SCIPvarGetName(var), (void*)expr);
         for( v = 0; v < nvars; ++v )
         {
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
         }
         SCIPfreeBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip));
         return SCIP_OKAY;
      }

      /* - grad(g(x*))_i x*_i */
      constant -= derivative * val;

      /* grad(g(x*))_i x_i */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef*derivative) );
   }

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
   }
   SCIPfreeBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip));

   /* add constant */
   SCIPaddRowprepConstant(rowprep, coef*constant);

   return SCIP_OKAY;
}

/** constructs perspective linearization of a given expression and adds it to rowprep */
static
SCIP_RETCODE addPerspectiveLinearisation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_HASHMAP*         scvars,             /**< hashmap linking semicontinuous vars to their data */
   SCIP_ROWPREP*         rowprep,            /**< a rowprep where the linearization is stored */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be linearized */
   SCIP_Real             coef,               /**< coefficient of expr */
   SCIP_VAR*             bvar,               /**< binary variable */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_Bool*            success             /**< indicates whether the linearization could be computed */
   )
{
   SCIP_SOL* sol0;
   SCIP_Real* vals0;
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_VAR** vars;
   SCIP_Real scalar_prod, fval, fval0;
   int nvars, v;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(scvars != NULL);
   assert(rowprep != NULL);
   assert(expr != NULL);
   assert(bvar != NULL);
   assert(success != NULL);

   /* add the cut: auxvar >= (x - x0) \nabla f(sol) + (f(sol) - f(x0) - (sol - x0) \nabla f(sol)) z + f(x0),
    * where x is semicontinuous, z is binary and x0 is the value of x when z = 0
    */

   SCIP_CALL( SCIPcreateSol(scip, &sol0, NULL) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip)) );
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, varexprs, &nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vals0, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, nvars) );
#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "bvar = \n");
   SCIPprintVar(scip, bvar, NULL);
   SCIPdebugMsg(scip, NULL, "pexpr = \n");
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
#endif

   /* get x0 */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_SCVARDATA* vardata;
      vars[v] = SCIPgetConsExprExprVarVar(varexprs[v]);
      vardata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(scvars, (void*)vars[v]);
      vals0[v] = vardata->val0;
   }

   /* set x to x0 in sol0 */
   SCIP_CALL( SCIPsetSolVals(scip, sol0, nvars, vars, vals0) );

   /* get f(x0) */
   SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, expr, sol0, 0) );
   fval0 = SCIPgetConsExprExprValue(expr);
   SCIP_CALL( SCIPfreeSol(scip, &sol0) );

   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(fval0)) )
   {
      *success = FALSE;
      for(v = 0; v < nvars; ++v)
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
      }
      SCIPfreeBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip));
      SCIPfreeBlockMemoryArray(scip, &vals0, nvars);
      SCIPfreeBlockMemoryArray(scip, &vars, nvars);
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", fval0, (void*)expr);
      return SCIP_OKAY;
   }

   /* get f(sol) */
   SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, expr, sol, 0) );
   fval = SCIPgetConsExprExprValue(expr);

   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(fval)) )
   {
      *success = FALSE;
      for(v = 0; v < nvars; ++v)
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
      }
      SCIPfreeBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip));
      SCIPfreeBlockMemoryArray(scip, &vals0, nvars);
      SCIPfreeBlockMemoryArray(scip, &vars, nvars);
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", fval, (void*)expr);
      return SCIP_OKAY;
   }

   /* add (f(sol) - f(x0))z + f(x0) */
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, bvar, coef*(fval - fval0)) );
   SCIPaddRowprepConstant(rowprep, coef*fval0);

   /* compute gradient */
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );

   /* gradient evaluation error -> skip */
   if( SCIPgetConsExprExprDerivative(expr) == SCIP_INVALID ) /*lint !e777*/
   {
      *success = FALSE;
      for(v = 0; v < nvars; ++v)
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
      }
      SCIPfreeBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip));
      SCIPfreeBlockMemoryArray(scip, &vals0, nvars);
      SCIPfreeBlockMemoryArray(scip, &vars, nvars);
      SCIPdebugMsg(scip, "gradient evaluation error for %p\n", (void*)expr);
      return SCIP_OKAY;
   }

   scalar_prod = 0.0;
   for(v = 0; v < nvars; ++v)
   {
      SCIP_VAR* var;
      var = SCIPgetConsExprExprVarVar(varexprs[v]);

      /* add xi f'xi(sol) */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef*SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, var)) );
      /* add -x0i f'xi(sol) */
      SCIPaddRowprepConstant(rowprep, -coef*vals0[v]*SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, var));

      /* compute -(soli - x0i) f'xi(sol) */
      scalar_prod -= (SCIPgetSolVal(scip, sol, var) - vals0[v])*SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, var);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
   }
   SCIPfreeBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip));
   SCIPfreeBlockMemoryArray(scip, &vals0, nvars);
   SCIPfreeBlockMemoryArray(scip, &vars, nvars);

   /* add -(sol - x0) \nabla f(sol)) z */
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, bvar, coef*scalar_prod) );

   return SCIP_OKAY;
}

/** frees nlhdlrexprdata structure */
static
SCIP_RETCODE freeNlhdlrExprData(
   SCIP*                         scip,            /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata   /**< nlhdlr expression data */
   )
{
   int c;

   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->convterms), nlhdlrexprdata->convtermssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->convcoefs), nlhdlrexprdata->convtermssize);

   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->onoffterms), nlhdlrexprdata->nonoffterms);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->onoffcoefs), nlhdlrexprdata->nonoffterms);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->bvars), nlhdlrexprdata->nonoffterms);

   if( nlhdlrexprdata->varexprs != NULL )
   {
      for( c = 0; c < nlhdlrexprdata->nvarexprs; ++c )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(nlhdlrexprdata->varexprs[c])) );
      }
      SCIPfreeBlockMemoryArray(scip, &nlhdlrexprdata->varexprs, nlhdlrexprdata->nvarexprs);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrPerspective)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrPerspective(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}


/** callback to free data of handler */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataPerspective)
{ /*lint --e{715}*/
   SCIP_HASHMAPENTRY* entry;
   SCIP_SCVARDATA* data;
   int c;

   if( (*nlhdlrdata)->scvars != NULL )
   {
      for( c = 0; c < SCIPhashmapGetNEntries((*nlhdlrdata)->scvars); ++c )
      {
         entry = SCIPhashmapGetEntry((*nlhdlrdata)->scvars, c);
         if( entry != NULL )
         {
            data = (SCIP_SCVARDATA*) SCIPhashmapEntryGetImage(entry);
            SCIPfreeBlockMemory(scip, &data);
         }
      }
      SCIPhashmapFree(&((*nlhdlrdata)->scvars));
   }

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}


/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataPerspective)
{  /*lint --e{715}*/
   SCIP_CALL( freeNlhdlrExprData(scip, *nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}



/** callback to be called in initialization */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINIT(nlhdlrInitPerspective)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSHDLR* linconshdlr;
   SCIP_CONSHDLR* varboundconshdlr;
   int c;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* look for important constraint handlers */
   linconshdlr = SCIPfindConshdlr(scip, "linear");
   varboundconshdlr = SCIPfindConshdlr(scip, "varbound");

   for( c = 0; c < SCIPgetNConss(scip); ++c )
   {
      SCIP_CONS* cons = SCIPgetConss(scip)[c];
      assert(cons != NULL);

      if( SCIPconsGetHdlr(cons) == linconshdlr )
      {
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");

         /* TODO extract information with interface functions from cons_linear.h */
      }
      else if( SCIPconsGetHdlr(cons) == varboundconshdlr )
      {
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");

         /* TODO extract information with interface functions from cons_varbound.h */
      }
      else
      {
         /* TODO check whether other constraint handlers are important */
      }
   }

   /* TODO store data in nlhdlrdata */

   return SCIP_OKAY;
}
#else
#define nlhdlrInitPerspective NULL
#endif


/** callback to be called in deinitialization */
static
SCIP_DECL_CONSEXPR_NLHDLREXIT(nlhdlrExitPerspective)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** callback to detect structure in expression tree
 *
 * We are looking for expressions of the form: \sum\limits_{i=1}^p g_i(x_i) + g_0(x_0), where:
 *  each vector x_i has a single fixed value x^{off}_i when a binary var b_i is 0;
 *  g_i, i=1,..,p are nonlinear and either all convex or all concave;
 *  g_0 is either linear or has the same curvature as g_i, i=1,..,p;
 *  p != 0.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectPerspective)
{ /*lint --e{715}*/
   SCIP_EXPRCURV child_curvature;
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_VAR* var;
   SCIP_VAR* bvar;
   int nvars, v, c, nchildren;
   SCIP_CONSEXPR_EXPR** children;
   SCIP_CONSEXPR_EXPR* pexpr;
   SCIP_Real* coefs;
   SCIP_Bool expr_is_onoff;
   SCIP_Bool var_is_sc;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_HASHMAP* onoffterms;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrdata != NULL);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "Called perspective detect, expr = %p: \n", expr);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
#endif

   *success = FALSE;

   if( SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_UNKNOWN || SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_LINEAR )
   {
      SCIPdebugMsg(scip, "curvature of expr %p is %s\n", expr, SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_LINEAR ? "linear" : "unknown");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   (*nlhdlrexprdata)->curvature = SCIPgetConsExprExprCurvature(expr);
   SCIPdebugMsg(scip, "expr %p is %s\n", expr, (*nlhdlrexprdata)->curvature == SCIP_EXPRCURV_CONVEX ? "convex" : "concave");

   SCIP_CALL( SCIPgetConsExprExprNVars(scip, conshdlr, expr, &nvars) );
   if( nlhdlrdata->scvars == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&(nlhdlrdata->scvars), SCIPblkmem(scip), nvars) );
   }

   /* prepare the list of terms */
   if( strcmp("sum", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr))) == 0 )
   {
      children = SCIPgetConsExprExprChildren(expr);
      nchildren = SCIPgetConsExprExprNChildren(expr);
      coefs = SCIPgetConsExprExprSumCoefs(expr);
   }
   else
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &children, 1) );
      *children = expr;
      nchildren = 1;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &coefs, 1) );
      *coefs = 1.0;
   }
   SCIP_CALL( SCIPhashmapCreate(&onoffterms, SCIPblkmem(scip), nchildren) );

   /* this loop collects terms that satisfy the conditions for g_i(x_i) and creates an entry in onoffterms for each bvar
    * that has any on/off terms depending on it; the image associated with bvar is the sum of all such on/off terms.
    * all other terms are stored in convterms
    */
   for( c = 0; c < nchildren; ++c )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip)) );
      SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, children[c], varexprs, &nvars) );
      bvar = NULL;
      expr_is_onoff = TRUE;

      /* a constant is not an on/off expression */
      if( nvars == 0 )
         expr_is_onoff = FALSE;

      /* all variables of an on/off term should be semicontinuous and depend on the same binary var */
      for( v = 0; v < nvars; ++v )
      {
         var = SCIPgetConsExprExprVarVar(varexprs[v]);
         SCIP_CALL( varIsSemicontinuous(scip, var, &bvar, nlhdlrdata->scvars, &var_is_sc) );
         if( !var_is_sc )
         {
            expr_is_onoff = FALSE;
            break;
         }
      }

      if( !expr_is_onoff || SCIPgetConsExprExprCurvature(children[c]) != (*nlhdlrexprdata)->curvature )
      {
         /* if the term is not on/off, add it to convterms */
         SCIP_CALL( addConvTerm(scip, *nlhdlrexprdata, coefs[c], children[c]) );
      }
      else
      {
#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "Adding on/off term: ");
         SCIPprintConsExprExpr(scip, conshdlr, children[c], NULL);
#endif
         /* if the term satisfies the requirements for g_i(x_i), add it to onoffterms */
         SCIP_CALL( addOnoffTerm(scip, conshdlr, onoffterms, &(*nlhdlrexprdata)->nonoffterms, coefs[c], children[c], bvar) );
      }

      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
      }
      SCIPfreeBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip));
   }

   if( strcmp("sum", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr))) != 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &children, 1);
      SCIPfreeBlockMemoryArrayNull(scip, &coefs, 1);
   }

   /* check curvature of the terms */
   if( !SCIPhashmapIsEmpty(onoffterms) )
   {
      *success = TRUE;
      for( c = 0; c < (*nlhdlrexprdata)->nconvterms; ++c )
      {
         child_curvature = SCIPexprcurvMultiply((*nlhdlrexprdata)->convcoefs[c], SCIPgetConsExprExprCurvature((*nlhdlrexprdata)->convterms[c]));
         if( child_curvature != (*nlhdlrexprdata)->curvature && child_curvature != SCIP_EXPRCURV_LINEAR )
         {
            SCIPdebugMsg(scip, "Non-convex cont term, curv %d, expr curv %d: detect will return false\n", child_curvature, (*nlhdlrexprdata)->curvature);
            *success = FALSE;
            break;
         }
      }
      SCIPdebugMsg(scip, "Found on/off terms: ");
      for( c = 0; c < SCIPhashmapGetNEntries(onoffterms); ++c )
      {
         SCIP_HASHMAPENTRY* entry = SCIPhashmapGetEntry(onoffterms, c);
         if( entry != NULL )
         {
            pexpr = (SCIP_CONSEXPR_EXPR*) SCIPhashmapEntryGetImage(entry);
#ifdef SCIP_DEBUG
            SCIPprintVar(scip, (SCIP_VAR*)SCIPhashmapEntryGetOrigin(entry), NULL);
            SCIPdebugMsg(scip, ": ");
            SCIPprintConsExprExpr(scip, conshdlr, pexpr, NULL);
#endif
            child_curvature = SCIPgetConsExprExprCurvature(pexpr);
            if( child_curvature != (*nlhdlrexprdata)->curvature )
            {
               SCIPdebugMsg(scip, "Non-convex onoff term, curv %d, expr curv %d: detect will return false\n", child_curvature, (*nlhdlrexprdata)->curvature);
               *success = FALSE;
            }
         }
      }
   }

   if( *success )
   {
      int nterm;

      /* depending on curvature, set enforcemethods */
      if( (*nlhdlrexprdata)->curvature == SCIP_EXPRCURV_CONVEX )
      {
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
         *enforcedbelow = TRUE;
         SCIPdebugMsg(scip, "detected expr to be convex -> can enforce expr <= auxvar\n");
      }
      else if( (*nlhdlrexprdata)->curvature == SCIP_EXPRCURV_CONCAVE )
      {
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
         *enforcedabove = TRUE;
         SCIPdebugMsg(scip, "detected expr to be concave -> can enforce expr >= auxvar\n");
      }
      /* save varexprs to nlhdlrexprdata */
      SCIP_CALL( SCIPgetConsExprExprNVars(scip, conshdlr, expr, &(*nlhdlrexprdata)->nvarexprs) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->varexprs, (*nlhdlrexprdata)->nvarexprs) );
      SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, (*nlhdlrexprdata)->varexprs, &((*nlhdlrexprdata)->nvarexprs)) );

      /* move the on/off terms to an array */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->onoffterms, (*nlhdlrexprdata)->nonoffterms) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->onoffcoefs, (*nlhdlrexprdata)->nonoffterms) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->bvars, (*nlhdlrexprdata)->nonoffterms) );
      nterm = 0;
      for( c = 0; c < SCIPhashmapGetNEntries(onoffterms); ++c )
      {
         SCIP_HASHMAPENTRY* entry = SCIPhashmapGetEntry(onoffterms, c);
         if( entry != NULL )
         {
            bvar = (SCIP_VAR*)SCIPhashmapEntryGetOrigin(entry);
            pexpr = (SCIP_CONSEXPR_EXPR*) SCIPhashmapEntryGetImage(entry);
            children = SCIPgetConsExprExprChildren(pexpr);
            nchildren = SCIPgetConsExprExprNChildren(pexpr);
            coefs = SCIPgetConsExprExprSumCoefs(pexpr);

            for( int nexpr = 0; nexpr < nchildren; ++nexpr )
            {
               assert(nterm < (*nlhdlrexprdata)->nonoffterms);
               (*nlhdlrexprdata)->onoffterms[nterm] = children[nexpr];
               (*nlhdlrexprdata)->onoffcoefs[nterm] = coefs[nexpr];
               (*nlhdlrexprdata)->bvars[nterm] = bvar;
               nterm++;
            }
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &pexpr) );
         }
      }
   }
   else
   {
      for( c = 0; c < SCIPhashmapGetNEntries(onoffterms); ++c )
      {
         SCIP_HASHMAPENTRY* entry = SCIPhashmapGetEntry(onoffterms, c);
         if( entry != NULL )
         {
            pexpr = (SCIP_CONSEXPR_EXPR*) SCIPhashmapEntryGetImage(entry);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &pexpr) );
         }
      }
      SCIP_CALL( freeNlhdlrExprData(scip, *nlhdlrexprdata) );
      SCIPfreeBlockMemory(scip, nlhdlrexprdata);
      *nlhdlrexprdata = NULL;
   }

   SCIPhashmapFree(&onoffterms);

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalauxPerspective)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(auxvalue != NULL);

   *auxvalue = SCIPgetConsExprExprValue(expr);

   return SCIP_OKAY;
}


/** callback to detect structure in expression tree */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(nlhdlrInitSepaPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#endif


/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(nlhdlrExitSepaPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#endif


/** nonlinear handler separation callback
 *
 * Applies perspective linearization to on/off terms and gradient linearization to everything else.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRSEPA(nlhdlrSepaPerspective)
{ /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;
   SCIP_ROW* row;
   SCIP_VAR* auxvar;
   int i;
   SCIP_CONSEXPR_EXPR* pexpr;
   SCIP_Bool success;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_Real pcoef;
   SCIP_VAR* bvar;

   *result = SCIP_DIDNOTFIND;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "sepa method of perspective nonlinear handler called for expr %p: ", expr);
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   assert(scip != NULL);
   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(result != NULL);
   assert(ncuts != NULL);
   assert(nlhdlrdata != NULL);

   *ncuts = 0;

   /* if estimating on non-convex side, then do nothing */
   if( ( overestimate && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONVEX) ||
       (!overestimate && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONCAVE) )
   {
      SCIPdebugMsg(scip, "Estimating on non-convex side, do nothing\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );
   auxvar = SCIPgetConsExprExprAuxVar(expr);
   assert(auxvar != NULL);
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, auxvar, -1.0) );

   if( SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) )
      SCIPaddRowprepConstant(rowprep, SCIPgetConsExprExprSumConstant(expr));

   success = TRUE; /* think positive */

   /* handle convex terms */
   for( i = 0; i < nlhdlrexprdata->nconvterms && success; ++i )
   {
      SCIP_CALL( addGradientLinearisation(scip, conshdlr, rowprep, nlhdlrexprdata->convterms[i], nlhdlrexprdata->convcoefs[i], sol, &success) );
   }

   /* handle on/off terms */
   for( i = 0; i < nlhdlrexprdata->nonoffterms && success; ++i )
   {
      pexpr = nlhdlrexprdata->onoffterms[i];
      pcoef = nlhdlrexprdata->onoffcoefs[i];
      bvar = nlhdlrexprdata->bvars[i];
      SCIP_CALL( addPerspectiveLinearisation(scip, conshdlr, nlhdlrdata->scvars, rowprep, pexpr, pcoef, bvar, sol, &success) );
   }

   /* merge coefficients that belong to same variable */
   SCIPmergeRowprepTerms(scip, rowprep);

   rowprep->local = FALSE;
   if( success )
   {
      SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, mincutviolation, NULL, &success) );
   }

   /* if cut looks good (numerics ok and cutting off solution), then turn into row and add to sepastore */
   if( success )
   {
      SCIP_Bool infeasible;
      SCIP_CALL( SCIPgetRowprepRowConshdlr(scip, &row, rowprep, conshdlr) );
#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "Separating sol point\n");
      for( int v = 0; v < nlhdlrexprdata->nvarexprs; ++v )
      {
         SCIP_VAR* var = SCIPgetConsExprExprVarVar(nlhdlrexprdata->varexprs[v]);
         SCIPwriteVarName(scip, NULL, var, TRUE);
         SCIPinfoMessage(scip, NULL, ": %f\n",  SCIPgetSolVal(scip, sol, var));
      }
      SCIPinfoMessage(scip, NULL, "by perspective cut ");
      SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif

      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

      if( infeasible )
      {
         *result = SCIP_CUTOFF;
      }
      else
      {
         *result = SCIP_SEPARATED;
         ++*ncuts;
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}


/** nonlinear handler under/overestimation callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimatePerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#endif


/** nonlinear handler interval evaluation callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#endif


/** nonlinear handler callback for reverse propagation */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(nlhdlrReversepropPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#endif


/** nonlinear handler callback for branching scores */
static
SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE(nlhdlrBranchscorePerspective)
{ /*lint --e{715}*/
   SCIP_Real violation;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX || SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONCAVE);
   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
   assert(auxvalue == SCIPgetConsExprExprValue(expr)); /* given auxvalue should have been computed by nlhdlrEvalAuxConvex */  /*lint !e777*/
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->varexprs != NULL);
   assert(nlhdlrexprdata->nvarexprs > 0);
   assert(success != NULL);

   /* we separate only convex functions here, so there should be little use for branching
    * if violations are small or there are numerical issues, then we will not have generated a cut, though
    * in that case, we will still branch, that is, register branchscores for all depending var exprs
    */

   /* compute violation */
   if( auxvalue == SCIP_INVALID ) /*lint !e777*/
      violation = SCIPinfinity(scip); /* evaluation error -> we should branch */
   else if( SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX  )
      violation = auxvalue - SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr));
   else
      violation = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr)) - auxvalue;

   /* if violation is not on the side that we need to enforce, then no need for branching */
   if( violation <= 0.0 )
      return SCIP_OKAY;

   /* TODO try to figure out which variables appear linear and skip them here */
   for( i = 0; i < nlhdlrexprdata->nvarexprs; ++i )
   {
      assert(nlhdlrexprdata->varexprs[i] != NULL);
      assert(SCIPisConsExprExprVar(nlhdlrexprdata->varexprs[i]));

      SCIPaddConsExprExprBranchScore(scip, nlhdlrexprdata->varexprs[i], brscoretag, violation);
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/** nonlinear handler callback for reformulation */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE(nlhdlrReformulatePerspective)
{ /*lint --e{715}*/

   /* TODO detect structure */

   /* TODO create expression and store it in refexpr */

   /* set refexpr to expr and capture it if no reformulation is possible */
   *refexpr = expr;
   SCIPcaptureConsExprExpr(*refexpr);

   return SCIP_OKAY;
}
#endif

/*
 * nonlinear handler specific interface methods
 */

/** includes Perspective nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrPerspective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSEXPR_NLHDLR* nlhdlr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   /* create nonlinear handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );
   BMSclearMemory(nlhdlrdata);


   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectPerspective, nlhdlrEvalauxPerspective, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrPerspective);
   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrFreehdlrdataPerspective);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrFreeExprDataPerspective);
   SCIPsetConsExprNlhdlrInitExit(scip, nlhdlr, nlhdlrInitPerspective, nlhdlrExitPerspective);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, nlhdlrSepaPerspective, NULL, NULL);
   SCIPsetConsExprNlhdlrBranchscore(scip, nlhdlr, nlhdlrBranchscorePerspective);

   return SCIP_OKAY;
}
