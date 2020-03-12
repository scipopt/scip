/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr.c
 * @brief  tests basic nonlinear handler methods
 *
 * This test implements a nonlinear handler for bivariate quadratic expressions that are either convex or concave.
 * We do some simplifying assumptions, e.g., assume that the arguments are actual variables and not non-variable expressions.
 * Also separation is only implemented for the convex side.
 * These assumptions are ok for the example that is executed by the test.
 *
 * The test constructs a problem with two convex quadratic constraints.
 * Convexity is not exploited when using the separation methods of the expressions alone, as the quadratic terms
 * are considered separately. Thus, if no other nonlinear handlers take care of this, then only with this nonlinear handler
 * we can solve this problem by separation only (Kelley' cutting plane).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#define SCIP_PRIVATE_ROWPREP
#include "scip/cons_quadratic.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_var.h"

/*
 * NONLINEAR HANDLER (i.e. SCIP code)
 */

struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_Bool             initialized;        /**< whether handler has been initialized and not yet de-initialized */
};

/** compact storage for variables and coefficients in a bivariate quadratic term that is either convex or concave */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_VAR*             varx;               /**< first variable */
   SCIP_VAR*             vary;               /**< second variable */
   SCIP_Real             xcoef;              /**< coefficient of first variable linear term */
   SCIP_Real             ycoef;              /**< coefficient of second variable linear term */
   SCIP_Real             xycoef;             /**< coefficient of bilinear term */
   SCIP_Real             xxcoef;             /**< coefficient of first variable square term */
   SCIP_Real             yycoef;             /**< coefficient of second variable square term */
   SCIP_Real             constant;           /**< constant term */
   SCIP_Bool             convex;             /**< whether convex or concave */

   SCIP_CONSEXPR_EXPR*   exprx;               /**< expression corresponding to first variable */
   SCIP_CONSEXPR_EXPR*   expry;               /**< expression corresponding to second variable */
};


static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(freeHdlrData)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrdata != NULL);
   assert(*nlhdlrdata != NULL);
   assert(!(*nlhdlrdata)->initialized);

   SCIPfreeMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(freeExprData)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);

   SCIPfreeMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRINIT(initHdlr)
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(nlhdlr != NULL);

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   assert(!nlhdlrdata->initialized);

   nlhdlrdata->initialized = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLREXIT(exitHdlr)
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(nlhdlr != NULL);

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   assert(nlhdlrdata->initialized);

   nlhdlrdata->initialized = FALSE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(detectHdlr)
{
   SCIP_CONSEXPR_EXPRHDLR* powhdlr;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_CONSEXPR_NLHDLREXPRDATA exprdata;
   int c;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* only look at sum expressions */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr) )
      return SCIP_OKAY;

   powhdlr = SCIPfindConsExprExprHdlr(conshdlr, "pow");
   assert(powhdlr != NULL);

   BMSclearMemory(&exprdata);
   exprdata.constant = SCIPgetConsExprExprSumConstant(expr);

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      child = SCIPgetConsExprExprChildren(expr)[c];
      assert(child != NULL);

      if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrVar(conshdlr) )
      {
         SCIP_VAR* var;

         var = SCIPgetConsExprExprVarVar(child);
         assert(var != NULL);

         if( var == exprdata.varx )
         {
            exprdata.xcoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( var == exprdata.vary )
         {
            exprdata.ycoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.varx == NULL )
         {
            exprdata.varx = var;
            exprdata.exprx = child;
            assert(exprdata.xcoef == 0.0);
            exprdata.xcoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.vary == NULL )
         {
            exprdata.vary = var;
            exprdata.expry = child;
            assert(exprdata.ycoef == 0.0);
            exprdata.ycoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else
         {
            /* more than two variables -> give up */
            return SCIP_OKAY;
         }
      }
      else if( SCIPgetConsExprExprHdlr(child) == powhdlr )
      {
         SCIP_VAR* var;

         if( SCIPgetConsExprExprPowExponent(child) != 2.0 )
            return SCIP_OKAY;  /* only allow for exponent 2 (square terms) */

         assert(SCIPgetConsExprExprNChildren(child) == 1);
         child = SCIPgetConsExprExprChildren(child)[0];

         if( SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrVar(conshdlr) )
            return SCIP_OKAY;  /* only allow for variable as base of power */

         var = SCIPgetConsExprExprVarVar(child);
         if( var == exprdata.varx )
         {
            exprdata.xxcoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( var == exprdata.vary )
         {
            exprdata.yycoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.varx == NULL )
         {
            assert(exprdata.xxcoef == 0.0);
            exprdata.varx = var;
            exprdata.exprx = child;
            exprdata.xxcoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.vary == NULL )
         {
            assert(exprdata.yycoef == 0.0);
            exprdata.vary = var;
            exprdata.expry = child;
            exprdata.yycoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else
         {
            /* more than two variables -> give up */
            return SCIP_OKAY;
         }
      }
      else if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrProduct(conshdlr) )
      {
         SCIP_VAR* var1;
         SCIP_VAR* var2;

         if( SCIPgetConsExprExprNChildren(child) != 2 )
            return SCIP_OKAY; /* only allow for bilinear term */

         if( SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(child)[0]) != SCIPgetConsExprExprHdlrVar(conshdlr) )
            return SCIP_OKAY; /* only allow for variables in bilinear term */
         if( SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(child)[1]) != SCIPgetConsExprExprHdlrVar(conshdlr) )
            return SCIP_OKAY; /* only allow for variables in bilinear term */

         var1 = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[0]);
         var2 = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[1]);
         assert(var1 != var2);

         if( (var1 == exprdata.varx && var2 == exprdata.vary) || (var1 == exprdata.vary && var2 == exprdata.varx)  )
         {
            exprdata.xycoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( ((var1 == exprdata.varx) || (var2 == exprdata.varx)) && exprdata.vary == NULL )
         {
            assert(exprdata.xycoef == 0.0);
            exprdata.vary = (var1 == exprdata.varx) ? var2 : var1;
            exprdata.expry = (var1 == exprdata.varx) ? SCIPgetConsExprExprChildren(child)[1] : SCIPgetConsExprExprChildren(child)[0];
            exprdata.xycoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.varx == NULL )
         {
            assert(exprdata.xycoef == 0.0);
            assert(exprdata.vary == NULL);
            exprdata.varx = var1;
            exprdata.exprx = SCIPgetConsExprExprChildren(child)[0];
            exprdata.vary = var2;
            exprdata.expry = SCIPgetConsExprExprChildren(child)[0];
            exprdata.xycoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else
         {
            /* more than two variables */
            return SCIP_OKAY;
         }
      }
      else
      {
         /* unknown expression type */
         return SCIP_OKAY;
      }
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, " -> %gx%+gy%+gxy%+gx^2%+gy^2%+g (x=<%s>, y=<%s>)\n", exprdata.xcoef, exprdata.ycoef, exprdata.xycoef, exprdata.xxcoef, exprdata.yycoef, exprdata.constant, SCIPvarGetName(exprdata.varx), SCIPvarGetName(exprdata.vary));
#endif

   /* separable function is not of interest (for this unittest) */
   if( exprdata.xycoef == 0.0 )
      return SCIP_OKAY;

   /* check if convex or concave */
   if( exprdata.xxcoef >= 0.0 && exprdata.yycoef >= 0.0 && SCIPisGE(scip, 4.0 * exprdata.xxcoef * exprdata.yycoef, exprdata.xycoef*exprdata.xycoef) )
      exprdata.convex = TRUE;
   else if( exprdata.xxcoef <= 0.0 && exprdata.yycoef <= 0.0 && SCIPisGE(scip, 4.0 * exprdata.xxcoef * exprdata.yycoef, exprdata.xycoef*exprdata.xycoef) )
      exprdata.convex = FALSE;
   else
      return SCIP_OKAY; /* indefinite */

   /* communicate that we will enforce one side by separation */
   *enforcemethods |= exprdata.convex ? SCIP_CONSEXPR_EXPRENFO_SEPABELOW : SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
   *enforcedbelow |= exprdata.convex;
   *enforcedabove |= !exprdata.convex;
   *success = TRUE;
   SCIP_CALL( SCIPduplicateMemory(scip, nlhdlrexprdata, &exprdata) );

   /* communicate that we will also do inteval and reverseprop */
   *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_INTEVAL;
   *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_REVERSEPROP;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(evalauxHdlr)
{
   SCIP_Real xval;
   SCIP_Real yval;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->exprx) == nlhdlrexprdata->varx);
   assert(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->expry) == nlhdlrexprdata->vary);
   assert(auxvalue != NULL);

   xval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->varx);
   yval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vary);

   *auxvalue = nlhdlrexprdata->constant;
   *auxvalue += nlhdlrexprdata->xcoef * xval;
   *auxvalue += nlhdlrexprdata->xxcoef * xval * xval;
   *auxvalue += nlhdlrexprdata->ycoef * yval;
   *auxvalue += nlhdlrexprdata->yycoef * yval * yval;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRENFO(enfoHdlr)
{
   SCIP_VAR* auxvar;
   SCIP_ROWPREP* rowprep = NULL;
   SCIP_Bool infeasible;
   SCIP_Real xval;
   SCIP_Real yval;
   SCIP_Real xcoef;
   SCIP_Real ycoef;
   SCIP_Real side;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->exprx) == nlhdlrexprdata->varx);
   assert(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->expry) == nlhdlrexprdata->vary);

   *result = SCIP_DIDNOTFIND;

   /* get auxiliary variable */
   auxvar = SCIPgetConsExprExprAuxVar(expr);
   assert(auxvar != NULL);

   xval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->varx);
   yval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vary);

   /* adjust the reference point */
   xval = MIN(MAX(xval, SCIPvarGetLbLocal(nlhdlrexprdata->varx)), SCIPvarGetUbLocal(nlhdlrexprdata->varx));
   yval = MIN(MAX(yval, SCIPvarGetLbLocal(nlhdlrexprdata->vary)), SCIPvarGetUbLocal(nlhdlrexprdata->vary));

   if( nlhdlrexprdata->convex == !overestimate )
   {
      /* convex side -> linearize
       * constant + xcoef * x + ycoef * y + xxcoef*xval*xval + yycoef*yval*yval + xycoef*xval*yval
       *   + (2*xxcoef*xval + xycoef*yval) * (x - xval) + (2*yycoef*yval + xycoef*xval) * (y - yval)
       * = (xcoef + 2*xxcoef*xval + xycoef*yval) * x + (ycoef * 2*yycoef*yval + xycoef*xval) * y
       *   - xxcoef*xval*xval - yycoef*yval*yval - xycoef*xval*yval + constant
       */
      xcoef = nlhdlrexprdata->xcoef + 2.0 * nlhdlrexprdata->xxcoef * xval + nlhdlrexprdata->xycoef * yval;
      ycoef = nlhdlrexprdata->ycoef + 2.0 * nlhdlrexprdata->yycoef * yval + nlhdlrexprdata->xycoef * xval;
      side = nlhdlrexprdata->xxcoef * xval * xval + nlhdlrexprdata->yycoef * yval * yval + nlhdlrexprdata->xycoef * xval * yval - nlhdlrexprdata->constant;

      SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, FALSE) );
      SCIPaddRowprepSide(rowprep, side);
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, nlhdlrexprdata->varx, xcoef) );
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, nlhdlrexprdata->vary, ycoef) );
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, auxvar, -1.0) );

      /* SCIP_CALL( SCIPprintRow(scip, cut, NULL) ); */
   }
   else
   {
      /* we do not implement separation for this side (and did not advertise to do so in detect) */
   }

   if( rowprep == NULL )
      return SCIP_OKAY;

   /* check whether its violation and numerical properties are ok (and maybe improve) */
   /* TODO if not allowweakcuts, maybe use SCIPcleanupRowprep2 */
   SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, SCIPgetLPFeastol(scip), NULL, &success) );

   if( success )
   {
      SCIP_ROW* cut;

      SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "testhdlrcut_cvx");
      SCIP_CALL( SCIPgetRowprepRowConshdlr(scip, &cut, rowprep, conshdlr) );

      assert(-SCIPgetRowSolFeasibility(scip, cut, sol) >= SCIPgetLPFeastol(scip));

      /* add cut */
      SCIP_CALL( SCIPaddRow(scip, cut, FALSE, &infeasible) );
      *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;

#ifdef SCIP_DEBUG
      if( *result == SCIP_CUTOFF )
      {
         SCIPdebugMsg(scip, "add cut makes node infeasible!\n");
      }
      else
      {
         SCIPdebugMsg(scip, "add cut with violation %e\n", -SCIPgetRowSolFeasibility(scip, cut, sol));
      }
      SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

      /* release cut */
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }

   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(intevalHdlr)
{
   assert(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->exprx) == nlhdlrexprdata->varx);
   assert(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->expry) == nlhdlrexprdata->vary);

   SCIPintervalQuadBivar(SCIP_INTERVAL_INFINITY, interval, nlhdlrexprdata->xxcoef, nlhdlrexprdata->yycoef, nlhdlrexprdata->xycoef, nlhdlrexprdata->xcoef, nlhdlrexprdata->ycoef, SCIPgetConsExprExprActivity(scip, nlhdlrexprdata->exprx), SCIPgetConsExprExprActivity(scip, nlhdlrexprdata->expry));
   SCIPintervalAddScalar(SCIP_INTERVAL_INFINITY, interval, *interval, nlhdlrexprdata->constant);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(reversepropHdlr)
{
   SCIP_INTERVAL childbounds;
   SCIP_INTERVAL rhs;

   assert(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->exprx) == nlhdlrexprdata->varx);
   assert(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->expry) == nlhdlrexprdata->vary);

   rhs = SCIPgetConsExprExprActivity(scip, expr);
   SCIPintervalSubScalar(SCIP_INTERVAL_INFINITY, &rhs, rhs, nlhdlrexprdata->constant);

   /* solve conv({x in xbnds : xxcoef*x^2 + yycoef*y^2 + xycoef*x*y + xcoef*x + ycoef*y \in rhs, y \in ybnds}) */
   SCIPintervalSolveBivariateQuadExpressionAllScalar(SCIP_INTERVAL_INFINITY, &childbounds, nlhdlrexprdata->xxcoef, nlhdlrexprdata->yycoef, nlhdlrexprdata->xycoef, nlhdlrexprdata->xcoef, nlhdlrexprdata->ycoef, rhs, SCIPgetConsExprExprActivity(scip, nlhdlrexprdata->exprx), SCIPgetConsExprExprActivity(scip, nlhdlrexprdata->expry));
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, nlhdlrexprdata->exprx, childbounds, force, reversepropqueue, infeasible, nreductions) );

   if( !*infeasible )
   {
      /* solve conv({y in ybnds : yycoef*y^2 + xxcoef*x^2 + xycoef*y*x + ycoef*y + xcoef*x \in rhs, x \in xbnds}) */
      SCIPintervalSolveBivariateQuadExpressionAllScalar(SCIP_INTERVAL_INFINITY, &childbounds, nlhdlrexprdata->yycoef, nlhdlrexprdata->xxcoef, nlhdlrexprdata->xycoef, nlhdlrexprdata->ycoef, nlhdlrexprdata->xcoef, rhs, SCIPgetConsExprExprActivity(scip, nlhdlrexprdata->expry), SCIPgetConsExprExprActivity(scip, nlhdlrexprdata->exprx));
      SCIP_CALL( SCIPtightenConsExprExprInterval(scip, nlhdlrexprdata->expry, childbounds, force, reversepropqueue, infeasible, nreductions) );
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(copyHdlr)
{
   SCIP_CONSEXPR_NLHDLR* targetnlhdlr;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourceconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), "testhdlr") == 0);

   SCIP_CALL( SCIPallocClearMemory(targetscip, &nlhdlrdata) );

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(targetscip, targetconsexprhdlr, &targetnlhdlr,
      SCIPgetConsExprNlhdlrName(sourcenlhdlr), SCIPgetConsExprNlhdlrDesc(sourcenlhdlr), SCIPgetConsExprNlhdlrPriority(sourcenlhdlr), detectHdlr, evalauxHdlr, nlhdlrdata) );
   SCIPsetConsExprNlhdlrFreeHdlrData(targetscip, targetnlhdlr, freeHdlrData);
   SCIPsetConsExprNlhdlrFreeExprData(targetscip, targetnlhdlr, freeExprData);
   SCIPsetConsExprNlhdlrCopyHdlr(targetscip, targetnlhdlr, copyHdlr);
   SCIPsetConsExprNlhdlrInitExit(targetscip, targetnlhdlr, initHdlr, exitHdlr);
   SCIPsetConsExprNlhdlrSepa(targetscip, targetnlhdlr, NULL, enfoHdlr, NULL, NULL);
   SCIPsetConsExprNlhdlrProp(targetscip, targetnlhdlr, intevalHdlr, reversepropHdlr);

   return SCIP_OKAY;
}

/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;


/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   const char* input1 = "[expr] <test>: (<x>+<y>-0)^2 + (<y>-0)^2 <= 1.5;";
   const char* input2 = "[expr] <test>: (<x>+<y>-1)^2 + (<y>-1)^2 <= 1.0;";

   /* create scip with all plugins */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 5.0, -1.5, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 5.0, -2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input1,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );

   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input2,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

Test(conshdlr, nlhdlr, .init = setup, .fini = teardown,
   .description = "test basic functionality of nonlinear handler of the cons_expr constraint handler."
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPR_NLHDLR* nlhdlr;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPallocClearMemory(scip, &nlhdlrdata) );

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, conshdlr, &nlhdlr, "testhdlr", "tests nonlinear handler functionality", 0, detectHdlr, evalauxHdlr, nlhdlrdata) );

   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, freeHdlrData);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, freeExprData);
   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, copyHdlr);
   SCIPsetConsExprNlhdlrInitExit(scip, nlhdlr, initHdlr, exitHdlr);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, enfoHdlr, NULL, NULL);
   SCIPsetConsExprNlhdlrProp(scip, nlhdlr, intevalHdlr, reversepropHdlr);

   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
   /* SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", 1e-6) ); */
/*
   SCIP_CALL( SCIPpresolve(scip) );
   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );
   SCIP_CALL( SCIPprintTransProblem(scip, NULL, NULL, FALSE) );
*/
   SCIP_CALL( SCIPsolve(scip) );
   /* SCIP_CALL( SCIPprintBestSol(scip, NULL, TRUE) ); */

   cr_assert(SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL, "not solved to optimality");
   cr_assert(SCIPisFeasEQ(scip, SCIPgetPrimalbound(scip), -1.93649230212515), "optimal value not correct, expected -1.93649230212515, but got %.20g", SCIPgetPrimalbound(scip));
   cr_assert(SCIPgetNNodes(scip) <= 1, "convex NLP should be solved without branching, but took %" SCIP_LONGINT_FORMAT " nodes", SCIPgetNNodes(scip));
}
