/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
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

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/*
 * NONLINEAR HANDLER
 */

struct SCIP_NlhdlrData
{
   SCIP_Bool             initialized;        /**< whether handler has been initialized and not yet de-initialized */
};

/** compact storage for variables and coefficients in a bivariate quadratic term that is either convex or concave */
struct SCIP_NlhdlrExprData
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

   SCIP_EXPR*            exprx;               /**< expression corresponding to first variable */
   SCIP_EXPR*            expry;               /**< expression corresponding to second variable */
};


static
SCIP_DECL_NLHDLRFREEHDLRDATA(freeHdlrData)
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
SCIP_DECL_NLHDLRFREEEXPRDATA(freeExprData)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);

   SCIPfreeMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_NLHDLRINIT(initHdlr)
{
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(nlhdlr != NULL);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   assert(!nlhdlrdata->initialized);

   nlhdlrdata->initialized = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_NLHDLREXIT(exitHdlr)
{
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(nlhdlr != NULL);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   assert(nlhdlrdata->initialized);

   nlhdlrdata->initialized = FALSE;

   return SCIP_OKAY;
}

static
SCIP_DECL_NLHDLRDETECT(detectHdlr)
{
   SCIP_EXPR* child;
   SCIP_NLHDLREXPRDATA exprdata;
   int c;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(enforcing != NULL);
   assert(participating != NULL);

   /* only look at sum expressions */
   if( !SCIPisExprSum(scip, expr) )
      return SCIP_OKAY;

   BMSclearMemory(&exprdata);
   exprdata.constant = SCIPgetConstantExprSum(expr);

   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      child = SCIPexprGetChildren(expr)[c];
      assert(child != NULL);

      if( SCIPisExprVar(scip, child) )
      {
         SCIP_VAR* var;

         var = SCIPgetVarExprVar(child);
         assert(var != NULL);

         if( var == exprdata.varx )
         {
            exprdata.xcoef += SCIPgetCoefsExprSum(expr)[c];
         }
         else if( var == exprdata.vary )
         {
            exprdata.ycoef += SCIPgetCoefsExprSum(expr)[c];
         }
         else if( exprdata.varx == NULL )
         {
            exprdata.varx = var;
            exprdata.exprx = child;
            assert(exprdata.xcoef == 0.0);
            exprdata.xcoef = SCIPgetCoefsExprSum(expr)[c];
         }
         else if( exprdata.vary == NULL )
         {
            exprdata.vary = var;
            exprdata.expry = child;
            assert(exprdata.ycoef == 0.0);
            exprdata.ycoef = SCIPgetCoefsExprSum(expr)[c];
         }
         else
         {
            /* more than two variables -> give up */
            return SCIP_OKAY;
         }
      }
      else if( SCIPisExprPower(scip, child) )
      {
         SCIP_VAR* var;

         if( SCIPgetExponentExprPow(child) != 2.0 )
            return SCIP_OKAY;  /* only allow for exponent 2 (square terms) */

         assert(SCIPexprGetNChildren(child) == 1);
         child = SCIPexprGetChildren(child)[0];

         if( !SCIPisExprVar(scip, child) )
            return SCIP_OKAY;  /* only allow for variable as base of power */

         var = SCIPgetVarExprVar(child);
         if( var == exprdata.varx )
         {
            exprdata.xxcoef += SCIPgetCoefsExprSum(expr)[c];
         }
         else if( var == exprdata.vary )
         {
            exprdata.yycoef += SCIPgetCoefsExprSum(expr)[c];
         }
         else if( exprdata.varx == NULL )
         {
            assert(exprdata.xxcoef == 0.0);
            exprdata.varx = var;
            exprdata.exprx = child;
            exprdata.xxcoef = SCIPgetCoefsExprSum(expr)[c];
         }
         else if( exprdata.vary == NULL )
         {
            assert(exprdata.yycoef == 0.0);
            exprdata.vary = var;
            exprdata.expry = child;
            exprdata.yycoef = SCIPgetCoefsExprSum(expr)[c];
         }
         else
         {
            /* more than two variables -> give up */
            return SCIP_OKAY;
         }
      }
      else if( SCIPisExprProduct(scip, child) )
      {
         SCIP_VAR* var1;
         SCIP_VAR* var2;

         if( SCIPexprGetNChildren(child) != 2 )
            return SCIP_OKAY; /* only allow for bilinear term */

         if( !SCIPisExprVar(scip, SCIPexprGetChildren(child)[0]) )
            return SCIP_OKAY; /* only allow for variables in bilinear term */
         if( !SCIPisExprVar(scip, SCIPexprGetChildren(child)[1]) )
            return SCIP_OKAY; /* only allow for variables in bilinear term */

         var1 = SCIPgetVarExprVar(SCIPexprGetChildren(child)[0]);
         var2 = SCIPgetVarExprVar(SCIPexprGetChildren(child)[1]);
         assert(var1 != var2);

         if( (var1 == exprdata.varx && var2 == exprdata.vary) || (var1 == exprdata.vary && var2 == exprdata.varx)  )
         {
            exprdata.xycoef += SCIPgetCoefsExprSum(expr)[c];
         }
         else if( ((var1 == exprdata.varx) || (var2 == exprdata.varx)) && exprdata.vary == NULL )
         {
            assert(exprdata.xycoef == 0.0);
            exprdata.vary = (var1 == exprdata.varx) ? var2 : var1;
            exprdata.expry = (var1 == exprdata.varx) ? SCIPexprGetChildren(child)[1] : SCIPexprGetChildren(child)[0];
            exprdata.xycoef = SCIPgetCoefsExprSum(expr)[c];
         }
         else if( exprdata.varx == NULL )
         {
            assert(exprdata.xycoef == 0.0);
            assert(exprdata.vary == NULL);
            exprdata.varx = var1;
            exprdata.exprx = SCIPexprGetChildren(child)[0];
            exprdata.vary = var2;
            exprdata.expry = SCIPexprGetChildren(child)[0];
            exprdata.xycoef = SCIPgetCoefsExprSum(expr)[c];
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
   SCIP_CALL( SCIPprintExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, " -> %gx%+gy%+gxy%+gx^2%+gy^2%+g (x=<%s>, y=<%s>)\n", exprdata.xcoef, exprdata.ycoef, exprdata.xycoef, exprdata.xxcoef, exprdata.yycoef, exprdata.constant, SCIPvarGetName(exprdata.varx), exprdata.vary != NULL ? SCIPvarGetName(exprdata.vary) : "na");
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

   /* communicate that we will participate on one side by separation */
   *participating = exprdata.convex ? SCIP_NLHDLR_METHOD_SEPABELOW : SCIP_NLHDLR_METHOD_SEPAABOVE;
   /* communicate that we will also do inteval and reverseprop */
   *participating |= SCIP_NLHDLR_METHOD_ACTIVITY;
   /* everything where we participate, we do thoroughly */
   *enforcing |= *participating;

   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, exprdata.exprx, FALSE, TRUE, FALSE, FALSE) );
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, exprdata.expry, FALSE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPduplicateMemory(scip, nlhdlrexprdata, &exprdata) );

   return SCIP_OKAY;
}

static
SCIP_DECL_NLHDLREVALAUX(evalauxHdlr)
{
   SCIP_Real xval;
   SCIP_Real yval;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->exprx) == nlhdlrexprdata->varx);
   assert(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->expry) == nlhdlrexprdata->vary);
   assert(auxvalue != NULL);

   xval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->varx);
   yval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vary);

   *auxvalue = nlhdlrexprdata->constant;
   *auxvalue += nlhdlrexprdata->xcoef * xval;
   *auxvalue += nlhdlrexprdata->xxcoef * xval * xval;
   *auxvalue += nlhdlrexprdata->ycoef * yval;
   *auxvalue += nlhdlrexprdata->yycoef * yval * yval;
   *auxvalue += nlhdlrexprdata->xycoef * xval * yval;

   return SCIP_OKAY;
}

static
SCIP_DECL_NLHDLRENFO(enfoHdlr)
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
   assert(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->exprx) == nlhdlrexprdata->varx);
   assert(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->expry) == nlhdlrexprdata->vary);

   *result = SCIP_DIDNOTFIND;

   /* get auxiliary variable */
   auxvar = SCIPgetExprAuxVarNonlinear(expr);
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
      SCIProwprepAddSide(rowprep, side);
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
   SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIPgetLPFeastol(scip), NULL, &success) );

   if( success )
   {
      SCIP_ROW* cut;

      SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "testhdlrcut_cvx");
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
SCIP_DECL_NLHDLRINTEVAL(intevalHdlr)
{
   assert(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->exprx) == nlhdlrexprdata->varx);
   assert(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->expry) == nlhdlrexprdata->vary);

   SCIPintervalQuadBivar(SCIP_INTERVAL_INFINITY, interval, nlhdlrexprdata->xxcoef, nlhdlrexprdata->yycoef, nlhdlrexprdata->xycoef, nlhdlrexprdata->xcoef, nlhdlrexprdata->ycoef, SCIPexprGetActivity(nlhdlrexprdata->exprx), SCIPexprGetActivity(nlhdlrexprdata->expry));
   SCIPintervalAddScalar(SCIP_INTERVAL_INFINITY, interval, *interval, nlhdlrexprdata->constant);

   return SCIP_OKAY;
}

static
SCIP_DECL_NLHDLRREVERSEPROP(reversepropHdlr)
{
   SCIP_INTERVAL childbounds;
   SCIP_INTERVAL rhs;

   assert(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->exprx) == nlhdlrexprdata->varx);
   assert(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->expry) == nlhdlrexprdata->vary);

   rhs = bounds;
   SCIPintervalSubScalar(SCIP_INTERVAL_INFINITY, &rhs, rhs, nlhdlrexprdata->constant);

   /* solve conv({x in xbnds : xxcoef*x^2 + yycoef*y^2 + xycoef*x*y + xcoef*x + ycoef*y \in rhs, y \in ybnds}) */
   SCIPintervalSolveBivariateQuadExpressionAllScalar(SCIP_INTERVAL_INFINITY, &childbounds, nlhdlrexprdata->xxcoef, nlhdlrexprdata->yycoef, nlhdlrexprdata->xycoef, nlhdlrexprdata->xcoef, nlhdlrexprdata->ycoef, rhs, SCIPexprGetActivity(nlhdlrexprdata->exprx), SCIPexprGetActivity(nlhdlrexprdata->expry));
   SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, nlhdlrexprdata->exprx, childbounds, infeasible, nreductions) );

   if( !*infeasible )
   {
      /* solve conv({y in ybnds : yycoef*y^2 + xxcoef*x^2 + xycoef*y*x + ycoef*y + xcoef*x \in rhs, x \in xbnds}) */
      SCIPintervalSolveBivariateQuadExpressionAllScalar(SCIP_INTERVAL_INFINITY, &childbounds, nlhdlrexprdata->yycoef, nlhdlrexprdata->xxcoef, nlhdlrexprdata->xycoef, nlhdlrexprdata->ycoef, nlhdlrexprdata->xcoef, rhs, SCIPexprGetActivity(nlhdlrexprdata->expry), SCIPexprGetActivity(nlhdlrexprdata->exprx));
      SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, nlhdlrexprdata->expry, childbounds, infeasible, nreductions) );
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_NLHDLRCOPYHDLR(copyHdlr)
{
   SCIP_NLHDLR* targetnlhdlr;
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), "testhdlr") == 0);

   SCIP_CALL( SCIPallocClearMemory(targetscip, &nlhdlrdata) );

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(targetscip, &targetnlhdlr,
      SCIPnlhdlrGetName(sourcenlhdlr), SCIPnlhdlrGetDesc(sourcenlhdlr),
      SCIPnlhdlrGetDetectPriority(sourcenlhdlr), SCIPnlhdlrGetEnfoPriority(sourcenlhdlr),
      detectHdlr, evalauxHdlr, nlhdlrdata) );
   SCIPnlhdlrSetFreeHdlrData(targetnlhdlr, freeHdlrData);
   SCIPnlhdlrSetFreeExprData(targetnlhdlr, freeExprData);
   SCIPnlhdlrSetCopyHdlr(targetnlhdlr, copyHdlr);
   SCIPnlhdlrSetInitExit(targetnlhdlr, initHdlr, exitHdlr);
   SCIPnlhdlrSetSepa(targetnlhdlr, NULL, enfoHdlr, NULL, NULL);
   SCIPnlhdlrSetProp(targetnlhdlr, intevalHdlr, reversepropHdlr);

   return SCIP_OKAY;
}

/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;


/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   const char* input1 = "[nonlinear] <test>: (<x>+<y>-0)^2 + (<y>-0)^2 <= 1.5;";
   const char* input2 = "[nonlinear] <test>: (<x>+<y>-1)^2 + (<y>-1)^2 <= 1.0;";

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
   SCIP_CALL( SCIPparseCons(scip, &cons, input1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, input2, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
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
   .description = "test basic functionality of nonlinear handler of the nonlinear constraint handler."
   )
{
   SCIP_NLHDLR* nlhdlr;
   SCIP_NLHDLRDATA* nlhdlrdata;

   SCIP_CALL( SCIPallocClearMemory(scip, &nlhdlrdata) );

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, "testhdlr",
         "tests nonlinear handler functionality", 10000, 10000, detectHdlr, evalauxHdlr, nlhdlrdata) );

   SCIPnlhdlrSetFreeHdlrData(nlhdlr, freeHdlrData);
   SCIPnlhdlrSetFreeExprData(nlhdlr, freeExprData);
   SCIPnlhdlrSetCopyHdlr(nlhdlr, copyHdlr);
   SCIPnlhdlrSetInitExit(nlhdlr, initHdlr, exitHdlr);
   SCIPnlhdlrSetSepa(nlhdlr, NULL, enfoHdlr, NULL, NULL);
   SCIPnlhdlrSetProp(nlhdlr, intevalHdlr, reversepropHdlr);

#ifndef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
#endif
   /* SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", 1e-6) ); */
/*
   SCIP_CALL( SCIPpresolve(scip) );
   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );
   SCIP_CALL( SCIPprintTransProblem(scip, NULL, NULL, FALSE) );
*/
   SCIP_CALL( SCIPsolve(scip) );
   /* SCIP_CALL( SCIPprintBestSol(scip, NULL, TRUE) ); */

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );
#endif

   cr_assert(SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL, "not solved to optimality");
   cr_assert(SCIPisFeasEQ(scip, SCIPgetPrimalbound(scip), -1.93649230212515), "optimal value not correct, expected -1.93649230212515, but got %.20g", SCIPgetPrimalbound(scip));
   cr_assert(SCIPgetNNodes(scip) <= 1, "convex NLP should be solved without branching, but took %" SCIP_LONGINT_FORMAT " nodes", SCIPgetNNodes(scip));
}
