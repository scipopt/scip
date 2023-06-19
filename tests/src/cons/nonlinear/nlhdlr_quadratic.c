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

/**@file   nlhdlr_quadratic.c
 * @brief  tests quadratic nonlinear handler methods
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "include/scip_test.h"
/* XXX: need the consdata struct because we don't have getNlhdlrs or findNlhdlrs; I don't add those function because I'm unsure
 * we actually need them
 */
#include "scip/nlhdlr_quadratic.c"
#include "scip/nlhdlr.c"
#include "scip/cons_nonlinear.c"
#include "scip/nlhdlr_default.h"
#include "scip/expr_abs.h"
#include "scip/expr_exp.h"
#include "scip/expr_log.h"
#include "scip/expr_pow.h"
#include "scip/expr_product.h"
#include "scip/expr_sum.h"
#include "scip/expr_trig.h"
#include "scip/expr_var.h"
#include "scip/expr_value.h"
#include "scip/nlpi_ipopt.h"
#include "scip/struct_expr.h"

/*
 * TEST
 */

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* w;
static SCIP_VAR* z;
static SCIP_VAR* s;
static SCIP_VAR* t;

static SCIP_CONSHDLR* conshdlr;
static SCIP_NLHDLR* nlhdlr = NULL;

static RAYS* myrays = NULL;

#define EXPECTFEQ(a,b) cr_expect_float_eq(a, b, 1e-6, "%s = %g != %g (dif %g)", #a, a, b, ABS(a-b))

/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_nonlinear */
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );

   /* include some expr handlers */
   SCIP_CALL( SCIPincludeExprhdlrAbs(scip) );
   SCIP_CALL( SCIPincludeExprhdlrExp(scip) );
   SCIP_CALL( SCIPincludeExprhdlrLog(scip) );
   SCIP_CALL( SCIPincludeExprhdlrVar(scip) );
   SCIP_CALL( SCIPincludeExprhdlrValue(scip) );
   SCIP_CALL( SCIPincludeExprhdlrSum(scip) );
   SCIP_CALL( SCIPincludeExprhdlrPow(scip) );
   SCIP_CALL( SCIPincludeExprhdlrProduct(scip) );
   SCIP_CALL( SCIPincludeExprhdlrSin(scip) );
   SCIP_CALL( SCIPincludeExprhdlrCos(scip) );

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert_not_null(conshdlr);

   /* include quadratic and default nlhdlr */
   SCIP_CALL( SCIPincludeNlhdlrQuadratic(scip) );
   SCIP_CALL( SCIPincludeNlhdlrDefault(scip) );

   /* get quadratic nlhdl */
   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, "quadratic");
   cr_assert_not_null(nlhdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to SOLVING stage */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &s, "s", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &t, "t", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPaddVar(scip, s) );
   SCIP_CALL( SCIPaddVar(scip, t) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* free rays */
   if( myrays != NULL )
      freeRays(scip, &myrays);

   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &s) );
   SCIP_CALL( SCIPreleaseVar(scip, &t) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

static
void checkNQuad(
   SCIP_EXPR* expr,
   int enlin,
   int enquad,
   int enbilin
   )
{
   int nlinexprs;
   int nquadexprs;
   int nbilinexprs;

   SCIPexprGetQuadraticData(expr, NULL, &nlinexprs, NULL, NULL, &nquadexprs, &nbilinexprs, NULL, NULL);

   cr_expect_eq(nlinexprs, enlin, "Expecting %d linear expr, got %d\n", enlin, nlinexprs);
   cr_expect_eq(nquadexprs, enquad, "Expecting %d quadratic terms, got %d\n", enquad, nquadexprs);
   cr_expect_eq(nbilinexprs, enbilin, "Expecting %d bilinear terms, got %d\n", enbilin, nbilinexprs);
}

static
void checkQuadTerm(
   SCIP_EXPR* expr,
   int nterm,
   SCIP_EXPR* eexpr,
   SCIP_VAR* evar,
   SCIP_Real elincoef,
   SCIP_Real esqrcoef
   )
{
   SCIP_EXPR* qexpr;
   SCIP_Real lincoef ;
   SCIP_Real sqrcoef;

   SCIPexprGetQuadraticQuadTerm(expr, nterm, &qexpr, &lincoef, &sqrcoef, NULL, NULL, NULL);

   cr_expect_eq(elincoef, lincoef, "Expecting lincoef %g in quad term, got %g\n", elincoef, lincoef);
   cr_expect_eq(esqrcoef, sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", esqrcoef, sqrcoef);

   if( evar != NULL )
   {
      SCIP_VAR* var;
      cr_expect(SCIPisExprVar(scip, qexpr));

      var = SCIPgetVarExprVar(qexpr);
      cr_expect_eq(evar, var, "Expecting var %s in quad term, got %s\n", SCIPvarGetName(evar), SCIPvarGetName(var));
   }

   if( eexpr != NULL )
   {
      cr_expect_eq(qexpr, eexpr);
   }
}

/* detects x^2 + x as quadratic expression */
Test(nlhdlrquadratic, detectandfree1, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_EXPR* expr;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_NLHDLR_METHOD participatingexpected;
   SCIP_NLHDLR_METHOD enforcingexpected;
   SCIP_Bool infeasible;
   SCIP_Bool success;

   /* skip when no ipopt */
   if( ! SCIPisIpoptAvailableIpopt() )
      return;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[nonlinear] <test>: <x>^2 + <x> <= 1", TRUE, TRUE,
            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   success = FALSE;
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_assert(!infeasible);

   /* get expr and work with it */
   expr = SCIPgetExprNonlinear(cons);

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectQuadratic(scip, conshdlr, nlhdlr, expr, NULL, &enforcing, &participating, &nlhdlrexprdata) );
   cr_assert_not_null(nlhdlrexprdata);
   cr_assert_eq(nlhdlrexprdata->qexpr, expr);

   /* x^2 + x <= something is convex and so convex nlhdlr should take care of it; x^2 + x >= something is nonconvex
    * and so intersection cuts participates
    */
   participatingexpected = /*SCIP_NLHDLR_METHOD_SEPAABOVE |*/ SCIP_NLHDLR_METHOD_ACTIVITY;
   cr_expect_eq(participating, participatingexpected, "participating expecting %d got %d\n", participatingexpected, participating);

   enforcingexpected = SCIP_NLHDLR_METHOD_ACTIVITY;
   cr_expect_eq(enforcing, enforcingexpected, "enforcing expecting %d got %d\n", enforcingexpected, enforcing);

   nlhdlrFreeexprdataQuadratic(scip, nlhdlr, expr, &nlhdlrexprdata);

   checkNQuad(expr, 0, 1, 0);
   checkQuadTerm(expr, 0, NULL, x, 1.0, 1.0);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects x^2 + 2*x cos(y x^2) + cos(y x^2)^2 <= 1 as convex quadratic expression:
 * simplify yields x^2 + 2 * x cos(x^2 y) + cos(x^2 y)^2 <= 1 --> should detect x^2 + 2 x * w + w^2
 */
Test(nlhdlrquadratic, detectandfree2, .init = setup, .fini = teardown)
{
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_NLHDLR_METHOD participatingexpected;
   SCIP_NLHDLR_METHOD enforcingexpected;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_EXPR* cosexpr;
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Bool infeasible;

   /* skip when no ipopt */
   if( ! SCIPisIpoptAvailableIpopt() )
      return;

   /* create expression, simplify it and find common subexpressions*/
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[nonlinear] <test>: <x>^2 + 2 * <x> * cos(<y> * <x>^2) + cos(<y> * <x>^2)^2 <= 1", TRUE, TRUE,
            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   success = FALSE;
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_assert(!infeasible);

   /* get expr and work with it */
   expr = SCIPgetExprNonlinear(cons);

   /* get cosine expression */
   cr_assert_eq(SCIPexprGetNChildren(expr), 3);
   cosexpr = SCIPexprGetChildren(expr)[1]; /*  x * cos(x^2 y) */
   cosexpr = SCIPexprGetChildren(cosexpr)[1]; /* cos(x^2 y) */
   cr_assert_str_eq(SCIPexprhdlrGetName(SCIPexprGetHdlr(cosexpr)), "cos", "expecting cos got %s\n",
         SCIPexprhdlrGetName(SCIPexprGetHdlr(cosexpr)));

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectQuadratic(scip, conshdlr, nlhdlr, expr, NULL, &enforcing, &participating, &nlhdlrexprdata) );
   cr_assert_not_null(nlhdlrexprdata);
   participatingexpected = /*SCIP_NLHDLR_METHOD_SEPAABOVE | */SCIP_NLHDLR_METHOD_ACTIVITY;
   enforcingexpected = SCIP_NLHDLR_METHOD_ACTIVITY;
   cr_expect_eq(participating, participatingexpected, "part expecting %d got %d\n", participatingexpected, participating);
   cr_expect_eq(enforcing, enforcingexpected, "enfo expecting %d got %d\n", enforcingexpected, enforcing);

   nlhdlrFreeexprdataQuadratic(scip, nlhdlr, expr, &nlhdlrexprdata);

   checkNQuad(expr, 0, 2, 1);

   checkQuadTerm(expr, 0, NULL, x, 0.0, 1.0);
   checkQuadTerm(expr, 1, cosexpr, NULL, 0.0, 1.0);

#if DEFAULT_USEINTERCUTS
//   cr_expect(SCIPgetExprNAuxvarUsesNonlinear(quad.expr) > 0, "cos expr should have auxiliary variable!\n");
#endif


   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;
   SCIP_Real coef;
   SCIP_VAR* var;
   SCIPexprGetQuadraticBilinTerm(expr, 0, &expr1, &expr2, &coef, NULL, NULL);

   cr_assert_not_null(expr1);
   cr_assert_not_null(expr2);

   cr_expect(SCIPisExprVar(scip, expr1));
   var = SCIPgetVarExprVar(expr1);

   cr_expect_eq(var, x, "Expecting %s as first factor, got %s\n", SCIPvarGetName(x), SCIPvarGetName(var));
   cr_expect_eq(expr2, cosexpr);
   cr_expect_eq(2.0, coef, "Expecting bilinear coef of %g, got %g\n", 2.0, coef);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* properly detect quadratic expression in exp(abs(log(x^2 + 2 * x*y + y^2))) <= 1 */
Test(nlhdlrquadratic, detectandfree3, .init = setup, .fini = teardown)
{
   SCIP_EXPR* expr;
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Bool infeasible;

   /* create expression and simplify it */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[nonlinear] <test>: exp(abs(log(<x>^2 + 2 * <x> * <y> + <y> + 2 * <y>^2))) <= 1", TRUE, TRUE,
            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* adds locks which are needed for detectNlhdlrs */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_assert_not(infeasible);

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1) );

   /* get expr and work with it */
   expr = SCIPgetExprNonlinear(cons);

   /* expr is exponential expr */
   cr_assert_eq(SCIPexprGetNChildren(expr), 1);
   cr_assert_str_eq(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "exp", "expecting exp got %s\n",
         SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));

   /* expr is abs expr */
   expr = SCIPexprGetChildren(expr)[0];
   cr_assert_eq(SCIPexprGetNChildren(expr), 1);
   cr_assert_str_eq(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "abs", "expecting abs got %s\n",
         SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));

   /* expr is log expr */
   expr = SCIPexprGetChildren(expr)[0];
   cr_assert_eq(SCIPexprGetNChildren(expr), 1);
   cr_assert_str_eq(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "log", "expecting log got %s\n",
         SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));

   /* expr is sum expr */
   expr = SCIPexprGetChildren(expr)[0];
   cr_assert_eq(SCIPexprGetNChildren(expr), 4);
   cr_assert_str_eq(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "sum", "expecting sum got %s\n",
         SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));

#if 0
   /* it should be identified that child should not have aux vars because of locks */
   for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];
      cr_expect_null(child->auxvar);
   }
#endif

   checkNQuad(expr, 0, 2, 1);

   /* quadratic terms */
   checkQuadTerm(expr, 0, NULL, x, 0.0, 1.0);
   checkQuadTerm(expr, 1, NULL, y, 1.0, 2.0);


   ///* bilinear term */
   //SCIP_QUADEXPR_BILINTERM bilin;
   //cr_expect_eq(1, expr->enfos[0]->nlhdlrexprdata->quaddata->nbilinexprterms);
   //bilin = expr->enfos[0]->nlhdlrexprdata->quaddata->bilinexprterms[0];
   //cr_assert_not_null(bilin.expr1);
   //cr_assert_not_null(bilin.expr2);
   //cr_expect_eq(2.0, bilin.coef, "Expecting bilincoef %g in quad term, got %g\n", 2.0, bilin.coef);
   //cr_expect_eq(SCIPgetExprAuxVarNonlinear(bilin.expr1), y);
   //cr_expect_eq(SCIPgetExprAuxVarNonlinear(bilin.expr2), x);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableNonlinear(scip, conshdlr, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* x^2 + y^2 + w*z should not be propagated by this nlhandler */
Test(nlhdlrquadratic, notpropagablequadratic1, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x>^2 + <y>^2 + <w>*<z>", NULL, NULL, NULL) );
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   cr_expect_not(changed);
   cr_expect_not(infeasible);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = simplified;

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectQuadratic(scip, conshdlr, nlhdlr, expr, FALSE, &enforcing, &participating, &nlhdlrexprdata) );
#if DEFAULT_USEINTERCUTS
   /* should have detected separation only */
   cr_expect_not_null(nlhdlrexprdata);
   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_SEPABOTH, "got %d\n", participating);
   cr_expect_eq(enforcing, SCIP_NLHDLR_METHOD_NONE);
#else
   /* shouldn't have detected anything -> provides nothing */
   cr_expect_null(nlhdlrexprdata);
   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_NONE, "got %d\n", participating);
   cr_expect_eq(enforcing, SCIP_NLHDLR_METHOD_NONE);
#endif
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* log^2 x + sin^2 y + cos^2 z should not be handled by this nlhandler when intersection cuts are not available */
Test(nlhdlrquadratic, notpropagable2, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"log(<x>)^2 + sin(<y>)^2 + cos(<z>)^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   cr_expect(changed);
   cr_expect_not(infeasible);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = simplified;

   /* detect */
   SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/quadratic/useintersectioncuts", FALSE) );
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectQuadratic(scip, conshdlr, nlhdlr, expr, FALSE, &enforcing, &participating, &nlhdlrexprdata) );

   /* shouldn't have detected anything -> provides nothing */
   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_NONE);
   cr_assert_eq(enforcing, SCIP_NLHDLR_METHOD_NONE);
   cr_expect_null(nlhdlrexprdata);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* x^2 + y^2 + z^2 * x, should only provide propagation:
 * Note: we use this expression because variables are automatically detected to be
 * common subexpressions. Since we cannot call detect common subexpression to a given expression
 * as easily as calling simplify, we content with this work around
 * The alternative would be to create a constraint and canonilize it, then get the expression
 * and call the detection method of the quadratic to this expression. This is the cleanest way
 * and probably the way it should be done (TODO)
 */
Test(nlhdlrquadratic, onlyPropagation, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;
   SCIP_CONS* cons;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x>^2 + <y>^2 + <z>^2 * <x>", NULL, NULL, NULL) );
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   cr_expect(changed);
   cr_expect_not(infeasible);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = simplified;
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "cons", expr, -100.0, 100.0) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   expr = SCIPgetExprNonlinear(cons);

   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectQuadratic(scip, conshdlr, nlhdlr, expr, FALSE, &enforcing, &participating, &nlhdlrexprdata) );

   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_ACTIVITY, "got %d\n", participating);
   cr_expect_eq(enforcing, SCIP_NLHDLR_METHOD_ACTIVITY, "got %d\n", enforcing);
   cr_expect_not_null(nlhdlrexprdata);

   checkNQuad(expr, 0, 3, 1);

   SCIP_CALL( nlhdlrFreeexprdataQuadratic(scip, nlhdlr, expr, &nlhdlrexprdata) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

#if SCIP_DISABLED_CODE
/* TODO: it might be that test needs to be updated to new propagation algorithms */
/* test propagation of yz - xz = z(y - x) */
Test(nlhdlrquadratic, factorize, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"-<z>*<x> + <y>*<z>", NULL, NULL, NULL) );
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   cr_expect(changed);
   cr_expect_not(infeasible);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = simplified;
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* detect */
   participating = SCIP_NLHDLR_METHOD_NONE;
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectQuadratic(scip, conshdlr, nlhdlr, expr, FALSE, &enforcing, &participating, &nlhdlrexprdata) );

   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_ACTIVITY /*ALL*/, "got %d\n", participating);
   cr_expect_eq(enforcing, SCIP_NLHDLR_METHOD_ACTIVITY, "got %d\n", enforcing);
   cr_expect_not_null(nlhdlrexprdata);

   cr_expect(nlhdlrexprdata->quaddata->nquadexprs == 3);

   /* no auxiliary variables should have been created */
   cr_expect(SCIPgetExprNAuxvarUsesNonlinear(nlhdlrexprdata->quaddata->quadexprterms[0].expr) == 0);
   cr_expect(SCIPgetExprNAuxvarUsesNonlinear(nlhdlrexprdata->quaddata->quadexprterms[1].expr) == 0);
   cr_expect(SCIPgetExprNAuxvarUsesNonlinear(nlhdlrexprdata->quaddata->quadexprterms[2].expr) == 0);

   /* check internal structure */
   cr_expect_eq(expr->children[0]->children[0], nlhdlrexprdata->quaddata->quadexprterms[0].expr); /* x should be first */
   cr_expect_eq(expr->children[0]->children[1], nlhdlrexprdata->quaddata->quadexprterms[1].expr); /* then z */
   cr_expect_eq(expr->children[1]->children[0], nlhdlrexprdata->quaddata->quadexprterms[2].expr); /* finally y */
   cr_expect_eq(nlhdlrexprdata->quaddata->bilinexprterms[0].expr1, nlhdlrexprdata->quaddata->bilinexprterms[1].expr1); /* z should be the first on both */

   /* interval evaluate */
   SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &interval);
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, expr, 0, NULL, NULL) );
   SCIP_CALL( nlhdlrIntevalQuadratic(scip, nlhdlr, expr, nlhdlrexprdata, &interval, NULL, NULL) );

   cr_expect_float_eq(interval.inf, matinf, 1e-7); cr_expect_leq(interval.inf, matinf);
   cr_expect_float_eq(interval.sup, matsup, 1e-7); cr_expect_geq(interval.sup, matsup);

   /* test reverse propagation */
   {
      SCIP_Bool infeasible = FALSE;
      int nreductions = 0;
      exprinterval.inf = 35;
      exprinterval.sup = 35;
      SCIP_CALL( SCIPdismantleExpr(scip, NULL, expr) );
      SCIP_CALL( nlhdlrReversepropQuadratic(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, exprinterval, &infeasible, &nreductions) );
      SCIP_CALL( SCIPdismantleExpr(scip, NULL, expr) );
      cr_expect_eq(nreductions, 2);
      cr_expect_not(infeasible);
      cr_expect_float_eq(SCIPvarGetLbLocal(z), -0.0741996, 1e-7);
      cr_expect_float_eq(SCIPvarGetUbLocal(x), -0.928007, 1e-6);
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
#endif

/* We test propagation with the following quadratic expression
 * x^2 - 3.1*x*y + 12.2*w*z + w^2 + 1.3*x*z - 4.8754*z^2 - 0.5*y^2 - 17.1*x + 22.02*y + 5*z - w
 * In the nlhdlr data we store the variables according to the order x < z < y < w because of repetition
 * FORWARD (or INTEVAL)
 * to compute an approximation of the max and min of this quadratic function,
 * it replaces it with the sum of univariate interval quadratics:
 * qx(x) + qy(y) + qw(w) + qz(z)
 * where
 * qx(x) = x^2 - 3.1*x*y - 17.1*x  + 1.3*x*z  -> x^2 + x*[-3.1*y - 17.1 + 1.3 * z]
 * qz(z) = -4.8754*z^2 + 5*z + 12.2*z*w       -> -4.8754*z^2 + z * [5 + 12.2*w]
 * qy(y) = -0.5*y^2 + 22.02*y                 -> -0.5*y^2 + 22.02*y
 * qw(w) = w^2 - w                            -> w^2 - 1.0 * w
 * Then computes the maximum and minimum of each q and adds them.
 *
 * BACKWARDS see mathematica code below. Now that we use the intervals of qx, qy, qz, qw
 * as computed by the forward propagation and not just doing the interval evaluation!
 *
 * Mathematicas code to generate answer: not very smart
 *
 *
Ix=Interval[{-101/100,101/100}];Iy=Interval[{7/100,9/100}];Iz=Interval[{-9/10,7/10}];Iw=Interval[{149/100,151/100}];
q[x_,y_,z_,w_]=x^2-3.1*x*y+12.2*z*w+w^2+16.3*z*x-4.8754*z^2-0.5*y^2-17.1*x+22.02*y+5*z-w;
qx[x_,y_,z_]=x^2+x*(-3.1*y-17.1+16.3*z);
qz[z_,w_]=-4.8754*z^2+z*(5+12.2*w);
qy[y_]=-0.5*y^2+22.02*y;
qw[w_]=w^2-1*w;
SetPrecision[q[Ix,Iy,Iz,Iw],15]
QXU=MaxValue[{qx[x,y,z],Element[{x},Ix]&&Element[{y},Iy]&&Element[{z},Iz]},{x,y,z}];
QYU=MaxValue[{qy[y],Element[{y},Iy]},{y}];
QWU=MaxValue[{qw[w],Element[{w},Iw]},{w}];
QZU=MaxValue[{qz[z,w],Element[{z},Iz]&&Element[{w},Iw]},{z,w}];
QXL=MinValue[{qx[x,y,z],Element[{x},Ix]&&Element[{y},Iy]&&Element[{z},Iz]},{x,y,z}];
QYL=MinValue[{qy[y],Element[{y},Iy]},{y}];
QWL=MinValue[{qw[w],Element[{w},Iw]},{w}];
QZL=MinValue[{qz[z,w],Element[{z},Iz]&&Element[{w},Iw]},{z,w}];
SetPrecision[QXL+QYL+QWL+QZL,20] (*computes lower bound*)
SetPrecision[QXU+QYU+QWU+QZU,20] (*computes upper bound*)
QX=Interval[{QXL,QXU}];QY=Interval[{QYL,QYU}];QW=Interval[{QWL,QWU}];QZ=Interval[{QZL,QZU}];
(*reverse propagation*)
Iq=Interval[{35,35}];
(* Propagating qx: first standard then bilinear terms *)
rhsx=Iq-(QY+QW+QZ); rx=Reduce[Exists[{y,z},Element[{y},Iy]&&Element[{z},Iz]&&Element[{qx[x,y,z]},rhsx]],Reals];
newIx=Interval[{MinValue[x,(rx)&&Element[{x},Ix],x],MaxValue[x,(rx)&&Element[{x},Ix],x]}]; (*improves x*)
Print["new Ix : ",newIx]

bilinxrhs = Interval[{MinValue[y/x-x,Element[{x},newIx]&&Element[{y},rhsx],{x,y}],MaxValue[y/x-x,Element[{x},newIx]&&Element[{y},rhsx],{x,y}]}];
Print["bilinxrhs ", bilinxrhs, " compare to ", rhsx/newIx - newIx]
consxforyz=Element[{-3.1*y-17.1+16.3*z},bilinxrhs];
newIy=Interval[{MinValue[y,consxforyz&&Element[{y},Iy]&&Element[{z},Iz],{y,z}],MaxValue[y,consxforyz&&Element[{y},Iy]&&Element[{z},Iz],{y,z}]}]; (* doesn't improve y*)
Print["new Iy : ",newIy]
newIz=Interval[{MinValue[z,consxforyz&&Element[{y},Iy]&&Element[{z},Iz],{y,z}],MaxValue[z,consxforyz&&Element[{y},Iy]&&Element[{z},Iz],{y,z}]}]; (*improves z*)
Print["new Iz : ",newIz]
(* Propagating qz *)
rhsz=Iq-(QX+QY+QW); rz=Reduce[Exists[{w},Element[{w},Iw]&&Element[{qz[z,w]},rhsz]],Reals];
newIz=Interval[{MinValue[z,(rz)&&Element[{z},newIz],z],MaxValue[z,(rz)&&Element[{z},newIz],z]}];
Print["new Iz : ",newIz]

(* 0 is in newIz so we don't do the following *)
(*conszforw=Element[{5+12.2*w},rhsz/newIz - newIz];
newIw=Interval[{MinValue[w,conszforw&&Element[{w},Iw],{w}],MaxValue[w,conszforw&&Element[{w},Iw],{w}]}]; (* doesn't improve y*)
Print["new Iw : ",newIw]*)
newIw=Iw;
(* Propagating qy *)
rhsy=Iq-(QX+QW+QZ); ry=Reduce[Element[{qy[y]},rhsy],Reals];
newIy=Interval[{MinValue[y,(ry)&&Element[{y},newIy],y],MaxValue[y,(ry)&&Element[{y},newIy],y]}];
Print["new Iy : ",newIy]
(* Propagating qw *)
rhsw=Iq-(QX+QZ+QY); rw=Reduce[Element[{qw[w]},rhsw],Reals];
newIw=Interval[{MinValue[w,(rw)&&Element[{w},newIw],w],MaxValue[w,(rw)&&Element[{w},newIw],w]}];
Print["new Iw : ",newIw]
 */
Test(nlhdlrquadratic, propagation_inteval, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_INTERVAL interval;
   SCIP_CONS* cons;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;

   /* set bounds: is important to do it here so that interval evaluations work */
   SCIPchgVarLb(scip, x, -1.01); SCIPchgVarUb(scip, x, 1.01);
   SCIPchgVarLb(scip, y,  0.07); SCIPchgVarUb(scip, y, 0.09);
   SCIPchgVarLb(scip, z,  -0.9); SCIPchgVarUb(scip, z,  0.7);
   SCIPchgVarLb(scip, w,  1.49); SCIPchgVarUb(scip, w, 1.51);

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x>^2 - 3.1*<x>*<y> + 12.2*<z>*<w> + <w>^2 + 16.3*<z>*<x> - 4.8754*<z>^2 - 0.5*<y>^2 - 17.1*<x> + 22.02*<y> + 5*<z> - <w>", NULL, NULL, NULL) );
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   cr_expect(changed);
   cr_expect_not(infeasible);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = simplified;
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* the reverse propagation test requires that var-exprs get their activity updated immediately when the var bounds are updated
    * this happens when the var boundchange event is processed, which only happens for variables that are used in a constraint
    */
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "cons", expr, -SCIPinfinity(scip), SCIPinfinity(scip)) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   expr = SCIPgetExprNonlinear(cons);
   SCIP_CALL( SCIPaddCons(scip, cons) ); /* to register events */

   /* detect */
   SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/quadratic/useintersectioncuts", FALSE) );

   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectQuadratic(scip, conshdlr, nlhdlr, expr, FALSE, &enforcing, &participating, &nlhdlrexprdata) );

   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_ACTIVITY, "got %d\n", participating);
   cr_expect_eq(enforcing, SCIP_NLHDLR_METHOD_ACTIVITY, "got %d\n", enforcing);
   cr_expect_not_null(nlhdlrexprdata);

   ///* check internal sorting of factors */
   //for( int i = 0; i < nlhdlrexprdata->quaddata->nbilinexprterms; ++ i)
   //{
   //   /* x always first */
   //   cr_expect(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->quaddata->bilinexprterms[i].expr2) != x);
   //   /* w never first */
   //   cr_expect(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->quaddata->bilinexprterms[i].expr1) != w);

   //   /* z can only be second to x */
   //   if( SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->quaddata->bilinexprterms[i].expr2) == z )
   //      cr_expect(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->quaddata->bilinexprterms[i].expr1) == x);

   //   /* y can only be first to w */
   //   if( SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->quaddata->bilinexprterms[i].expr1) == y )
   //      cr_expect(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->quaddata->bilinexprterms[i].expr1) == w);
   //}

   /* interval evaluate */
   SCIP_CALL( SCIPevalExprActivity(scip, expr) );
   SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &interval);
   SCIP_CALL( nlhdlrIntevalQuadratic(scip, nlhdlr, expr, nlhdlrexprdata, &interval, NULL, NULL) );

   EXPECTFEQ(interval.inf, -54.1092139000245);
   EXPECTFEQ(interval.sup, 50.1438939955093);

   /* test reverse propagation */
   {
      int nreductions = 0;
      SCIP_INTERVAL bounds;
      infeasible = FALSE;
      SCIPintervalSet(&bounds, 35);
      //SCIP_CALL( SCIPdismantleExpr(scip, NULL, expr) );
      SCIP_CALL( nlhdlrReversepropQuadratic(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, bounds, &infeasible, &nreductions) );
      //SCIP_CALL( SCIPdismantleExpr(scip, NULL, expr) );
      cr_expect_eq(nreductions, 3, "expecting 3 got %d\n", nreductions); /* three because the z improved twice */
      cr_expect_not(infeasible);

      EXPECTFEQ(SCIPvarGetLbLocal(z), -0.0485777477946283);
      EXPECTFEQ(SCIPvarGetUbLocal(z), 0.0198745061962769);
      EXPECTFEQ(SCIPvarGetUbLocal(x), -0.559537393062365);
   }


   SCIP_CALL( nlhdlrFreeexprdataQuadratic(scip, nlhdlr, expr, &nlhdlrexprdata) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* code to generate an answer:
 * MinValue[y/x + 1/2*x, Element[{x}, Ix] && Element[{y}, Ir], {x, y}]
 * MaxValue[y/x + 1/2*x, Element[{x}, Ix] && Element[{y}, Ir], {x, y}]
 * for some intervals Ix and Ir (one can also change the coef 1/2 above)
 */
Test(nlhdlrquadratic, bilin_rhs_range, .init = setup, .fini = teardown)
{
   SCIP_INTERVAL rhs;
   SCIP_INTERVAL exprdom;
   SCIP_INTERVAL range;

   rhs.inf = -8.0;
   rhs.sup = -4.0;

   exprdom.inf = 1.0;
   exprdom.sup = 10.0;

   computeRangeForBilinearProp(exprdom, 0.25, rhs, &range);
   EXPECTFEQ(range.inf, -8.25);
   EXPECTFEQ(range.sup, -2.0);

   exprdom.inf = -10.0;
   exprdom.sup = -1.0;

   computeRangeForBilinearProp(exprdom, 0.25, rhs, &range);
   EXPECTFEQ(range.inf, 2.0);
   EXPECTFEQ(range.sup, 8.25);

   rhs.inf = -8.0;
   rhs.sup = 4.0;

   computeRangeForBilinearProp(exprdom, 0.25, rhs, &range);
   EXPECTFEQ(range.inf, -3.75);
   EXPECTFEQ(range.sup, 8.25);

   computeRangeForBilinearProp(exprdom, 0.0, rhs, &range);
   EXPECTFEQ(range.inf, -4.0);
   EXPECTFEQ(range.sup, 8.0);

   rhs.inf = 0.0;
   rhs.sup = 4.0;

   computeRangeForBilinearProp(exprdom, -0.5, rhs, &range);
   EXPECTFEQ(range.inf, -5.4);
   EXPECTFEQ(range.sup, -0.5);

   exprdom.inf = -1.0;
   exprdom.sup = 10.0;

   computeRangeForBilinearProp(exprdom, 0.25, rhs, &range);
   cr_expect(SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, range));
}

/* test propagation of x*y + z + z^2; this is interesting to see how reverse propagation handles the term x*y.
 *
 * x in [1,5], y in [1, +inf], and z in [4, 5].
 *
 * Backward propagation:
 *
Ix = Interval[{1, 5}]; Iy = Interval[{1, 5}]; Iz = Interval[{4, 5}];
cons = x*y + z - z^2 == 10;
xmax = MaxValue[{x, Element[{x}, Ix] && y >= 1 && Element[{z}, Iz] && cons}, {x, y, z}]
xmin = MinValue[{x, Element[{x}, Ix] && y >= 1 && Element[{z}, Iz] && cons}, {x, y, z}]
ymax = MaxValue[{y, Element[{x}, Ix] && y >= 1 && Element[{z}, Iz] && cons}, {x, y, z}]
ymin = MinValue[{y, Element[{x}, Ix] && y >= 1 && Element[{z}, Iz] && cons}, {x, y, z}]
zmax = MaxValue[{z, Element[{x}, Ix] && y >= 1 && Element[{z}, Iz] && cons}, {x, y, z}]
zmin = MinValue[{z, Element[{x}, Ix] && y >= 1 && Element[{z}, Iz] && cons}, {x, y, z}]
 *
 */
/* FIXME: nlhdlr_quadratic now treats x*y as an argument, so detect and propagation calls for the x*y expression need to be added here */
#if SCIP_DISABLED_CODE
Test(nlhdlrquadratic, propagation_freq1vars, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;
   SCIP_INTERVAL interval;

   /* set bounds: is important to do it here so that interval evaluations work */
   SCIPchgVarLb(scip, x, 1.0); SCIPchgVarUb(scip, x, 5.0);
   SCIPchgVarLb(scip, y, 1.0);
   SCIPchgVarLb(scip, z, 4.0); SCIPchgVarUb(scip, z, 5.0);

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<y>*<x> + <z> - <z>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPprintExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPsimplifyExpr(scip, conshdlr, expr, &simplified, &changed, &infeasible) );
   cr_expect(changed);
   cr_expect_not(infeasible);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = simplified;
   SCIP_CALL( SCIPprintExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");


   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectQuadratic(scip, conshdlr, nlhdlr, expr, FALSE, &enforcing, &participating, &nlhdlrexprdata) );

   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_ACTIVITY, "got %d\n", participating);
   cr_expect_eq(enforcing, SCIP_NLHDLR_METHOD_ACTIVITY, "got %d\n", enforcing);
   cr_expect_not_null(nlhdlrexprdata);

   /* interval evaluate */
   SCIP_CALL( SCIPevalExprActivity(scip, conshdlr, expr, &interval, FALSE, FALSE) );
   //SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &interval);
   //SCIP_CALL( nlhdlrIntevalQuadratic(scip, nlhdlr, expr, nlhdlrexprdata, &interval, NULL, FALSE, NULL) );

   //cr_expect_float_eq(interval.inf, matinf, 1e-7, "got %f, expected %f\n", interval.inf, matinf); cr_expect_leq(interval.inf, matinf);
   //cr_expect_float_eq(interval.sup, matsup, 1e-7, "got %f, expected %f\n", interval.sup, matsup); cr_expect_geq(interval.sup, matsup);

   /* test reverse propagation */
   {
      int nreductions = 0;
      infeasible = FALSE;
      SCIPintervalSet(&interval, 10);
      SCIP_CALL( SCIPdismantleExpr(scip, NULL, expr) );
      SCIP_CALL( nlhdlrReversepropQuadratic(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, interval, &infeasible, &nreductions) );
      SCIP_CALL( SCIPdismantleExpr(scip, NULL, expr) );
      cr_expect_eq(nreductions, 2);
      cr_expect_not(infeasible);
   }

   /* check result */
   {
      int i;
      SCIP_Real expectedlowerbounds[] = {1.0, 4.4, 4.0};
      SCIP_Real expectedupperbounds[] = {5.0, 30.0, 5.0};
      SCIP_VAR* vars[] = {x, y, z};

      for( i = 0; i < (int)(sizeof(vars)/sizeof(SCIP_VAR*)); ++i )
      {
         cr_expect_float_eq(SCIPvarGetLbLocal(vars[i]), expectedlowerbounds[i], 1e-7, "var %s expecting %g, got %g\n",
               SCIPvarGetName(vars[i]), expectedlowerbounds[i], SCIPvarGetLbLocal(vars[i]));
         cr_expect_float_eq(SCIPvarGetUbLocal(vars[i]), expectedupperbounds[i], 1e-7, "var %s expecting %g, got %g\n",
               SCIPvarGetName(vars[i]), expectedupperbounds[i], SCIPvarGetUbLocal(vars[i]));
      }
   }

   /* register enforcer info in expr and free */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(expr->enfos), 1) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &(expr->enfos[0])) );
   expr->enfos[0]->nlhdlr = nlhdlr;
   expr->enfos[0]->nlhdlrexprdata = nlhdlrexprdata;
   expr->nenfos = 1;
   expr->enfos[0]->issepainit = FALSE;

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
#endif

/* Intersection cuts tests
 * =======================
 *
 * A way to generate simple LPs with a given set of rays is as follows.
 * Collect the rays as columns in a matrix R.
 * The LP that represents the cone with rays R is:
 * -R^{-1} x <= 0
 *
 * The cone with rays R and appex vertex is then given by:
 * -R^{-1} (x - vertex) <= 0
 *
 * Then, the cut in the slack variables is given by 1/interpoint of the ray with the S-free set >= 1.
 * Let us call this cut alpha^t s >= 1.
 * Thus, in the original variables the cut is
 * alpha^t (R^{-1} (x - vertex)) >= 1
 *
 */

TestSuite(interCuts, .init = setup, .fini = teardown);

static
void simplifyAndDetect(
   SCIP_CONS**           cons,
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata, /**< nlhdlr expression data */
   const char*           str                 /**< string to parse for constraint */
   )
{
   SCIP_EXPR* expr;
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_Bool enforcing;
   SCIP_Bool participating;
   SCIP_Bool infeasible;
   SCIP_Bool success;

   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, cons, str, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success)
         );
   cr_assert(success);

   infeasible = TRUE;
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_assert(!infeasible);
   expr = SCIPgetExprNonlinear(*cons);

   /* SCIP_CALL( SCIPprintExpr(scip, conshdlr, expr, NULL) ); */
   /* SCIPinfoMessage(scip, NULL, "\n"); */

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   nlhdlrdata->useintersectioncuts = TRUE;
   SCIP_CALL( nlhdlrDetectQuadratic(scip, conshdlr, nlhdlr, expr, *cons, &enforcing, &participating, nlhdlrexprdata) );

   cr_expect_not_null(nlhdlrexprdata);
}

static
void registerAndFree(
   SCIP_CONS*           cons,
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata  /**< nlhdlr expression data */
   )
{
   SCIP_EXPR* expr;
   SCIP_EXPR_OWNERDATA* ownerdata;

   expr = SCIPgetExprNonlinear(cons);
   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ownerdata->enfos), 1) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &(ownerdata->enfos[0])) );
   ownerdata->enfos[0]->nlhdlr = nlhdlr;
   ownerdata->enfos[0]->nlhdlrexprdata = nlhdlrexprdata;
   ownerdata->nenfos = 1;
   ownerdata->enfos[0]->issepainit = FALSE;

   ///* if there is an nlhdlr, then there must also be an auxvar */
   //SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, expr, NULL) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

static
void testRays(
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata, /**< nlhdlr expression data */
   SCIP_VAR* auxvar,
   SCIP_VAR** vars,
   SCIP_VAR** consvars,
   SCIP_Real* expectedrayscoefs,
   int* expectedraysidx,
   int* expectedlppos,
   int* expectedbegin,
   int expectednrays,
   int expectednnonz,
   int expectedncols
   )
{
   SCIP_EXPR* quaddata;
   SCIP_EXPR** linexprs;
   SCIP_Bool success;
   SCIP_COL** cols;
   int nquadexprs;
   int nlinexprs;
   int ncols;

   quaddata = nlhdlrexprdata->qexpr;

   ncols = SCIPgetNLPCols(scip);
   cols = SCIPgetLPCols(scip);

   cr_expect(ncols == expectedncols);
   SCIPexprGetQuadraticData(quaddata, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, NULL, NULL, NULL);

   /* first check that vars in lp and constraint are sorted as vars and consvars, respectively */
   for( int i = 0; i < ncols; ++i )
   {
      cr_assert_eq(SCIPcolGetVar(cols[i]), vars[i], "expected %s got %s\n", SCIPvarGetName(vars[i]), SCIPvarGetName(SCIPcolGetVar(cols[i])));
   }
   for( int i = 0; i < nquadexprs; ++i )
   {
      SCIP_EXPR* expr;
      SCIPexprGetQuadraticQuadTerm(quaddata, i, &expr, NULL, NULL, NULL, NULL, NULL);

      printf("%d got %s exp %s\n", i, SCIPvarGetName(SCIPgetExprAuxVarNonlinear(expr)), SCIPvarGetName(consvars[i]));
      cr_assert_eq(SCIPgetExprAuxVarNonlinear(expr), consvars[i]);
   }
   for( int i = 0; i < nlinexprs; ++i )
   {
      printf("%d got %s exp %s\n", i, SCIPvarGetName(SCIPgetExprAuxVarNonlinear(linexprs[i])), SCIPvarGetName(consvars[i + nquadexprs]));
      cr_assert_eq(SCIPgetExprAuxVarNonlinear(linexprs[i]), consvars[i + nquadexprs]);
   }

   /*
    * create rays and test them
    */
   SCIP_CALL( createAndStoreSparseRays(scip, nlhdlrexprdata, auxvar, &myrays, &success) );
   cr_expect(success);
   cr_expect_eq(myrays->nrays, expectednrays, "e %d g %d\n", expectednrays, myrays->nrays);
   cr_expect_eq(myrays->raysbegin[myrays->nrays], expectednnonz, "e %d g %d\n", expectednnonz, myrays->raysbegin[myrays->nrays]);

   for( int i = 0; i < expectednnonz; ++i )
   {
      cr_expect_float_eq(myrays->rays[i], expectedrayscoefs[i], 1e-9, "%d-th entry: expected %g, got %g\n", i,
            expectedrayscoefs[i], myrays->rays[i]);
   }
   cr_expect_arr_eq(myrays->raysidx, expectedraysidx, expectednnonz * sizeof(int));
   cr_expect_arr_eq(myrays->lpposray, expectedlppos, expectednrays * sizeof(int));
   cr_expect_arr_eq(myrays->raysbegin, expectedbegin, (expectednrays + 1) * sizeof(int));
}

/* builds lp so that every var is nonbasic at lower */
static
void buildAndSolveSimpleProbingLP(void)
{
   SCIP_Bool cutoff;
   SCIP_Bool lperror;

   SCIP_CALL( SCIPchgVarObjProbing(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarObjProbing(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarObjProbing(scip, z, 1.0) );
   SCIP_CALL( SCIPchgVarObjProbing(scip, w, 1.0) );

   SCIP_CALL( SCIPchgVarLbProbing(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarLbProbing(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarLbProbing(scip, z, 0.0) );
   SCIP_CALL( SCIPchgVarLbProbing(scip, w, 0.0) );

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
   //cr_expect_not(cutoff); /* cutoff == TRUE since the LP is solved to optimality with a obtimal value greater than cutoffbound */
   cr_expect_not(lperror);

   /* all variables should be nonbasic */
   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(x)), SCIP_BASESTAT_LOWER);
   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(y)), SCIP_BASESTAT_LOWER);
   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(z)), SCIP_BASESTAT_LOWER);
   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(w)), SCIP_BASESTAT_LOWER);
}

/* builds lp so that x, y, and z are basic and the rays are
 * what rays should be
 * x      1     0     0
 * y  =   0 s1  1 s2  1 s3
 * z      0     2     0
 */
static
void buildAndSolveSimpleProbingLP2(void)
{
   SCIP_ROW* row1;
   SCIP_ROW* row2;
   SCIP_ROW* row3;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   /* - x <= -1 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row1, "row1", -SCIPinfinity(scip), -1.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, x,  -1.0) );

   SCIP_CALL( SCIPaddRowProbing(scip, row1) );
   SCIP_CALL( SCIPreleaseRow(scip, &row1) );

   /* -z/2 <= -1/2 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row2, "row2", -SCIPinfinity(scip), -0.5, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row2, z,  -0.5) );

   SCIP_CALL( SCIPaddRowProbing(scip, row2) );
   SCIP_CALL( SCIPreleaseRow(scip, &row2) );

   /* -y + z/2 <= 3/2 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row3, "row3", -SCIPinfinity(scip), 1.5, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row3, y,  -1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row3, z,  0.5) );

   SCIP_CALL( SCIPaddRowProbing(scip, row3) );
   SCIP_CALL( SCIPreleaseRow(scip, &row3) );

   /* set objective coefficient */
   SCIP_CALL( SCIPchgVarObjProbing(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarObjProbing(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarObjProbing(scip, z, 1.0) );



   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
   cr_assert_not(lperror);
   /* interestingly this fails because the cutoffbound of scip is .0001 or something TODO: figure out why */
   //cr_assert_not(cutoff);
   SCIP_CALL( SCIPprintSol(scip, NULL, NULL, TRUE) );

   /* all variables should be nonbasic */
   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(x)), SCIP_BASESTAT_BASIC, "got %d\n", SCIPcolGetBasisStatus(SCIPvarGetCol(x)));
   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(y)), SCIP_BASESTAT_BASIC);
   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(z)), SCIP_BASESTAT_BASIC);

}

/* function to test cuts: coefficients should be normalize so to have norm 1 */
static
void testCut(
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata, /**< nlhdlr expression data */
   SCIP_CONS*            cons,
   SCIP_Bool             overestimate,
   SCIP_Real*            expectedcoefs,
   SCIP_VAR**            expectedvars,
   SCIP_Real             expectedlhs,
   int                   expectedncoefs
   )
{
   SCIP_EXPR* expr;
   SCIP_ROWPREP* rowprep;
   SCIP_COL** cols;
   SCIP_ROW* cut;
   SCIP_Real* coefs;
   SCIP_Real enorm;
   SCIP_Real cutnorm;
   SCIP_Real side;
   int nnonz;
   SCIP_Bool success;

   /* compute expected cut's norm: this norm does not have to be equal to cutnorm! */
   enorm = 0.0;
   for( int i = 0; i < expectedncoefs; i++ )
      enorm += SQR( expectedcoefs[i] );
   enorm = SQRT( enorm );
   cr_assert(enorm > 0);

   expr = SCIPgetExprNonlinear(cons);

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_LEFT, TRUE) );

   SCIP_CALL( generateIntercut(scip, expr, SCIPnlhdlrGetData(nlhdlr), nlhdlrexprdata, cons, NULL, rowprep, overestimate, &success ) );

   cr_expect(success);

   /* create cut from rowprep and get some info */
   SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, rowprep, cons) );

   cutnorm = SCIProwGetNorm(cut);
   nnonz = SCIProwGetNNonz(cut);
   cols = SCIProwGetCols(cut);
   coefs = SCIProwGetVals(cut);

   SCIPprintRow(scip, cut, NULL);

   /* ~~check same number of coefficients~~
    * currently not, as with HiGHS a coef of 2e-9 is obtained, which messes up nonzero count but still passes the coef checks below
    */
   /* cr_expect_eq(nnonz, expectedncoefs, "expected %d coefs, but got %d\n", expectedncoefs, nnonz); */

   /* check norm not zero */
   cr_assert(cutnorm > 0);

   /* check side */
   side = SCIProwGetLhs(cut);
   cr_expect(!SCIPisInfinity(scip, -side));
   cr_expect_float_eq(side / cutnorm, expectedlhs / enorm, 1e-6, "expecting lhs %g, got %g\n", expectedlhs / enorm, side
         / cutnorm);

    /* check coefficients */
   for( int j = 0; j < nnonz; j++ )
   {
      SCIP_Real expcoeffound = 0.0;
       for( int i = 0; i < expectedncoefs; i++ )
          if( SCIPcolGetVar(cols[j]) == expectedvars[i] )
            expcoeffound = expectedcoefs[i];
      cr_expect_float_eq(coefs[j] / cutnorm, expcoeffound / enorm, 1e-6, "expecting cut coef %g, got %g\n",
         expcoeffound / enorm, coefs[j] / cutnorm);
   }
   /* and another round to check that there isn't an expected coef missing */
   for( int i = 0; i < expectedncoefs; i++ )
   {
      SCIP_Real coeffound = 0.0;
      for( int j = 0; j < nnonz; j++ )
         if( SCIPcolGetVar(cols[j]) == expectedvars[i] )
            coeffound = coefs[j];
      cr_expect_float_eq(coeffound / cutnorm, expectedcoefs[i] / enorm, 1e-6, "expecting cut coef %g, got %g\n",
         expectedcoefs[i] / enorm, coeffound / cutnorm);
   }

   /* free cut */
   SCIPreleaseRow(scip, &cut);
   SCIPfreeRowprep(scip, &rowprep);
}

/* test that stored rays are
 * x    1     0     0     0
 * y  = 0 x + 1 y + 0 w + 0 z
 * z    0     0     0     1
 */
Test(interCuts, testRays1)
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;

   /* simplify and detect quadratic structure in: 6xy + 2x^2 -2z^2 + 2 <= 0 */
   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: 6.0*<x>*<y> + 2.0*<x>^2 - 2.0*<z>^2 + 2 <= 0");

   /*
    * build LP so that every var is non-basic
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );
   buildAndSolveSimpleProbingLP();


   /* rays should be
    * x    1     0     0     0
    * y  = 0 x + 1 y + 0 w + 0 z
    * z    0     0     0     1
    */
   {
      SCIP_Real expectedrayscoefs[3] = {1.0, 1.0, 1.0};
      int expectedraysidx[3] = {0, 1, 2};
      int expectedlppos[3] = {0, 1, 3};
      int expectedbegin[4] = {0, 1, 2, 3};
      int expectednrays = 3;
      int expectednnonz = 3;
      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[3] = {x, y, z};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

/* test that stored rays are
 * x    1     0     0     0
 * y  = 0 x + 1 y + 0 w + 0 z
 */
Test(interCuts, testRays2)
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;

   /* simplify and detect quadratic structure in: 6xy + 2x^2 + 2 <= 0 */
   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: 6.0*<x>*<y> + 2.0*<x>^2 + 2 <= 0");

   /*
    * build LP so that every var is non-basic
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );
   buildAndSolveSimpleProbingLP();


   /* rays should be
    * x    1     0     0     0
    * y  = 0 x + 1 y + 0 w + 0 z
    */
   {
      int expectednnonz = 2;
      SCIP_Real expectedrayscoefs[2] = {1.0, 1.0};
      int expectedraysidx[2] = {0, 1};
      int expectednrays = 2;
      int expectedlppos[2] = {0, 1};
      int expectedbegin[3] = {0, 1, 2};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[2] = {x, y};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}
/* test that stored rays are
 * y    0     1     0     0
 * z  = 0 x + 0 y + 0 w + 1 z
 */
Test(interCuts, testRays3)
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;

   /* simplify and detect quadratic structure in: 6zy + 2z^2 + 2 <= 0 */
   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: 6.0*<z>*<y> + 2.0*<z>^2 + 2 <= 0");

   /*
    * build LP so that every var is non-basic
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );
   buildAndSolveSimpleProbingLP();


   /* rays should be
    * y    0     1     0     0
    * z  = 0 x + 0 y + 0 w + 1 z
    */
   {
      int expectednnonz = 2;
      SCIP_Real expectedrayscoefs[2] = {1.0, 1.0};
      int expectedraysidx[2] = {0, 1};
      int expectednrays = 2;
      int expectedlppos[2] = {1, 3};
      int expectedbegin[3] = {0, 1, 2};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[2] = {y, z};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}
/* test that stored rays are
 * z  = 0 x + 0 y + 0 w + 1 z
 */
Test(interCuts, testRays4)
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;

   /* simplify and detect quadratic structure in: 6z + 2z^2 + 2 <= 0 */
   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: 6.0*<z> + 2.0*<z>^2 + 2 <= 0");

   /*
    * build LP so that every var is non-basic
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );
   buildAndSolveSimpleProbingLP();


   /* rays should be
    * z  = 0 x + 0 y + 0 w + 1 z
    */
   {
      int expectednnonz = 1;
      SCIP_Real expectedrayscoefs[1] = {1.0};
      int expectedraysidx[2] = {0};
      int expectednrays = 1;
      int expectedlppos[1] = {3};
      int expectedbegin[2] = {0, 1};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[1] = {z};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

/* test with tableau:
 *
 * x    -3/2     0    -47/12   -1     0
 * y  +    1 w+  0 z+  0    s+  2 t = 0
 *
 * where w and s are nonbasic at lower and z and t are nonbasic at upper
 * Thus, rays should be
 *
 * x    3/2    0    47/12  -1
 * y  =  -1 w  0fz  0    s  2ft
 * w      1    0    0       0
 * z      0   -1fz  0       0
 * s      0    0    1       0
 * t      0    0    0      -1ft
 *
 * z and t are going to be linear variables
 *
 * Note SCIP will add slack variables nonetheless, so this will contribute ray entries
 * the slack variables coefficients are +1 in the tableau.
 * I don't know if one can really tell at which bound they are going to be fixed; this depends
 * on the LP solver we use. Apparently, for CPLEX, the row is going to be active at its lhs,
 * which means that the base status is going to be at lower (see lpi.h). However, this means that
 * the slack variable is active at its upper bound! This is different for XPRESS for example.
 *
 */
Test(interCuts, testRays5)
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;

   /* simplify and detect quadratic structure in: x*y + x*w - z + w*s + t + 2 <= 0 */
   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: <x>*<y> + <x>*<w> - <z> + <w>*<s> + <t> + 2 <= 0");

   /* build LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   /* create and solve LP */
   {
      SCIP_ROW* row1;
      SCIP_ROW* row2;
      SCIP_Bool lperror;

      /* x  +  -3/2 w  +  0z    -47/12s   -1t = 0 */
      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row1, "row1", 0.0, 0.0, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, x,  1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, w,  -3.0/2) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, s, -47.0/12) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, t, -1.0) );

      SCIP_CALL( SCIPaddRowProbing(scip, row1) );
      SCIP_CALL( SCIPreleaseRow(scip, &row1) );

      /* y  +    1 w+  0 z+  0    s+  2 t = 0 */
      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row2, "row2", 0.0, 0.0, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarToRow(scip, row2, y,  1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, row2, w,  1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, row2, t,  2.0) );

      SCIP_CALL( SCIPaddRowProbing(scip, row2) );
      SCIP_CALL( SCIPreleaseRow(scip, &row2) );

      /* lower bound and objective of w,s so that they are nonbasic at lower */
      SCIP_CALL( SCIPchgVarLbProbing(scip, w, 0.0) );
      SCIP_CALL( SCIPchgVarLbProbing(scip, s, 0.0) );
      SCIP_CALL( SCIPchgVarObjProbing(scip, w, 1.0) );
      SCIP_CALL( SCIPchgVarObjProbing(scip, s, 1.0) );

      /* upper bound and objective of t,z so that they are nonbasic at upper */
      SCIP_CALL( SCIPchgVarUbProbing(scip, z, 0.0) );
      SCIP_CALL( SCIPchgVarUbProbing(scip, t, 0.0) );
      SCIP_CALL( SCIPchgVarObjProbing(scip, z, -1.0) );
      SCIP_CALL( SCIPchgVarObjProbing(scip, t, -1.0) );


      /* give bounds to x and y to prevent basestat zero (needed for soplex/xpress apparently) */
      SCIP_CALL( SCIPchgVarLbProbing(scip, x, -1.0e10) ); SCIP_CALL( SCIPchgVarUbProbing(scip, x, 1.0e10) );
      //SCIP_CALL( SCIPchgVarObjProbing(scip, x, -1.0e-5) );
      SCIP_CALL( SCIPchgVarLbProbing(scip, y, -1.0e10) ); SCIP_CALL( SCIPchgVarUbProbing(scip, y, 1.0e10) );
      //SCIP_CALL( SCIPchgVarObjProbing(scip, y,  1.0e-5) );

      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
      //cr_assert_not(cutoff); /* cutoff == TRUE since the LP is solved to optimality with a obtimal value greater than cutoffbound */
      cr_assert_not(lperror);
      SCIP_CALL( SCIPprintSol(scip, NULL, NULL, FALSE) );

      /* all variables should be nonbasic */
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(x)), SCIP_BASESTAT_BASIC, "got %d\n", SCIPcolGetBasisStatus(SCIPvarGetCol(x)));
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(y)), SCIP_BASESTAT_BASIC);
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(w)), SCIP_BASESTAT_LOWER);
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(s)), SCIP_BASESTAT_LOWER);
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(z)), SCIP_BASESTAT_UPPER);
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(t)), SCIP_BASESTAT_UPPER);
   }

   /* what rays should be
    * x    3/2    0    47/12  -1
    * y  =  -1 w  0fz  0    s  2ft
    * w      1    0    0       0
    * z      0   -1fz  0       0
    * s      0    0    1       0
    * t      0    0    0      -1ft
    *
    * the quadraitc: x*y + x*w - z + w*s + t + 2
    */
   {
      SCIP_ROW** rows;
      int nrows;
      int expectednnonz = 11;
      SCIP_Real factor1;
      SCIP_Real factor2;

      /* get basis status of slack variables to know what rays to expect */
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

      factor1 = SCIProwGetBasisStatus(rows[0]) == SCIP_BASESTAT_LOWER ? 1.0 : -1.0;
      factor2 = SCIProwGetBasisStatus(rows[1]) == SCIP_BASESTAT_LOWER ? 1.0 : -1.0;

      SCIP_Real expectedrayscoefs[11] = {
         3.0/2, -1.0, 1.0,
         -1.0,
         47.0/12, 1.0,
         -1.0, 2.0, -1.0,
         /* slacks entries */
         factor1 * 1.0,
         factor2 * 1.0
      };
      int expectedraysidx[11] = {
         0, 1, 2,
         4,
         0, 3,
         0, 1, 5,
         /* slacks entries */
         0,
         1
      };
      int expectednrays = 6;
      int expectedlppos[6] = {2, 3, 4, 5, -1, -2};
      int expectedbegin[7] = {0, 3, 4, 6, 9, 10, 11};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[6] = {x, y, w, s, z, t};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

/* test where all vars in quadratic basic. quadratic on x y and z, and linear on x
 * y*z + z^2 + x <= 0
 *
 * Rays:
 * x     1       1.3        0
 * y = 0.5 s1 + -0.5 s2 + 1.1 s3
 * z    -1      -0.7     -1.2
 *
 * It also tests case 4 coefficients of equations of rays and roots.
 * 1st ray should intersect in case 4b,
 * the other 2 intersect in case 4a
 */
Test(interCuts, testRays6)
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_Real vb[2];
   SCIP_Real vzlp[2];
   SCIP_Real wcoefs[2];
   SCIP_Real wzlp;
   SCIP_Real kappa;
   SCIP_Real coefs4a[5];
   SCIP_Real coefs4b[5];
   SCIP_Real coefscond[3];
   SCIP_Real root;
   SCIP_Bool success;
   SCIP_Real* eigenvalues;

   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: <y>*<z> + <z>^2 + <x> + 2.0 <= 2.0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /* build LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   /* create and solve LP */
   {
      SCIP_ROW* row1;
      SCIP_ROW* row2;
      SCIP_ROW* row3;
      SCIP_Bool lperror;

      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row1, "row1", -SCIPinfinity(scip), -31.0/18, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, x, -137.0/72) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, y, -13.0/6) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, z, -143.0/72) );

      SCIP_CALL( SCIPaddRowProbing(scip, row1) );
      SCIP_CALL( SCIPreleaseRow(scip, &row1) );

      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row2, "row2", -SCIPinfinity(scip), 5.0/9, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarToRow(scip, row2, x,  25.0/36) );
      SCIP_CALL( SCIPaddVarToRow(scip, row2, y,   5.0/3) );
      SCIP_CALL( SCIPaddVarToRow(scip, row2, z,  55.0/36) );

      SCIP_CALL( SCIPaddRowProbing(scip, row2) );
      SCIP_CALL( SCIPreleaseRow(scip, &row2) );

      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row3, "row3", -SCIPinfinity(scip), 35.0/18, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarToRow(scip, row3, x,  85.0/72) );
      SCIP_CALL( SCIPaddVarToRow(scip, row3, y,   5.0/6) );
      SCIP_CALL( SCIPaddVarToRow(scip, row3, z, 115.0/72) );

      SCIP_CALL( SCIPaddRowProbing(scip, row3) );
      SCIP_CALL( SCIPreleaseRow(scip, &row3) );

      /* set objective so that lp sol is finite */
      SCIP_CALL( SCIPchgVarObjProbing(scip, x, 100.0) );
      SCIP_CALL( SCIPchgVarObjProbing(scip, y, 118.0) );
      SCIP_CALL( SCIPchgVarObjProbing(scip, z, 101.0) );



      /* give bounds to x and y to prevent basestat zero */
      //SCIP_CALL( SCIPchgVarLbProbing(scip, x, -1.0e10) ); SCIP_CALL( SCIPchgVarUbProbing(scip, x, 1.0e10) );
      //SCIP_CALL( SCIPchgVarObjProbing(scip, x, -1.0e-5) );
      //SCIP_CALL( SCIPchgVarLbProbing(scip, y, -1.0e10) ); SCIP_CALL( SCIPchgVarUbProbing(scip, y, 1.0e10) );
      //SCIP_CALL( SCIPchgVarObjProbing(scip, y,  1.0e-5) );

      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
      cr_assert_not(lperror);

      /* interestingly this fails because the cutoffbound of scip is .5 something TODO: figure out why */
      /* cr_expect_not(cutoff); */
      SCIP_CALL( SCIPprintSol(scip, NULL, NULL, TRUE) );

      /* check nonbasic statuc */
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(x)), SCIP_BASESTAT_BASIC, "got %d\n", SCIPcolGetBasisStatus(SCIPvarGetCol(x)));
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(y)), SCIP_BASESTAT_BASIC);
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(z)), SCIP_BASESTAT_BASIC);
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(w)), SCIP_BASESTAT_ZERO);
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(s)), SCIP_BASESTAT_ZERO);
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(t)), SCIP_BASESTAT_ZERO);

      cr_expect_float_eq(SCIPvarGetLPSol(x), 1.0, 1e-12);
      cr_expect_float_eq(SCIPvarGetLPSol(y), -1.0, 1e-12);
      cr_expect_float_eq(SCIPvarGetLPSol(z), 1.0, 1e-12);
   }

   /* what rays should be
    * y =  1/2 s1 + -0.5 s2 + 1.1 s3
    * z     -1      -0.7     -1.2
    * x      1       1.3        0
    *
    * y*z + z^2 + x <= 0
    */
   {
      int expectednnonz = 8;
      SCIP_Real expectedrayscoefs[8] = {
         0.5, -1.0, 1.0,
         -0.5, -0.7, 1.3,
         1.1, -1.2
      };
      int expectedraysidx[8] = {
         0, 1, 2,
         0, 1, 2,
         0, 1
      };
      int expectednrays = 3;
      int expectedlppos[3] = {-1, -2, -3};
      int expectedbegin[4] = {0, 3, 6, 8};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[3] = {y, z, x} ;

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }

   /* check coefficients of equations and roots */
   SCIP_CALL( intercutsComputeCommonQuantities(scip, nlhdlrexprdata, NULL, 1.0, NULL, vb, vzlp, wcoefs, &wzlp, &kappa) );

   /* w should be w(x,y,z) = x */
   cr_expect_float_eq(wcoefs[0], 0.0, 1e-9);
   cr_expect_float_eq(wcoefs[1], 0.0, 1e-9);
   cr_expect_float_eq(wzlp, 1.0, 1e-12);
   cr_expect_float_eq(kappa, 0.0, 1e-9);

   /* rays case 4a */
   {
      SCIP_Real expectedcoefs4a[15] ={
         0.3977475644174330, -0.4571067811865475, 0.3535533905932738, 0.01843405788072983, 1.163423134802327,
         0.4302996026178613, 0.1050252531694167, 0.3535533905932738, 0.08811293459634177, 1.163423134802327,
         0.4508846494072806, -0.7985281374238570, 0.3535533905932738, -0.3861570698336348, 1.163423134802327};
      SCIP_Real expectedroots4a[3] = { 2.335548281119806, 1.661274085559792, 1.662222003744240};
      SCIP_Real expectedcoefs4b[15] = {
         0.05223665235168156, -0.1616116523516816, 0.125, -0.4785533905932738, 1.353553390593274,
         0.002757575950825020, 0.03713203435596426, 0.125, -0.5474873734152916, 1.353553390593274,
         0.1594117965644036, -0.2823223304703363, 0.125, -0.4492640687119285, 1.353553390593274};
      SCIP_Real expectedroots4b[3] = { 1.0 + SQRT( 2.0 ), 5.0 / 3, 5.0 * (1 + SQRT( 2.0 )) / 6};

      for( int nray = 0; nray < myrays->nrays; ++nray )
      {
         SCIP_CALL( computeRestrictionToRay(scip, nlhdlrexprdata, 1.0, TRUE, &myrays->rays[myrays->raysbegin[nray]],
                  &myrays->raysidx[myrays->raysbegin[nray]], myrays->raysbegin[nray + 1] - myrays->raysbegin[nray], vb, vzlp,
                  wcoefs, wzlp, kappa, coefs4a, coefs4b, coefscond, &success) );
         cr_expect(success);

         /* check coefficients */
         for( int i = 0; i < 5; ++i )
         {
            cr_expect_float_eq(coefs4a[i], expectedcoefs4a[5*nray + i], 1e-12, "case 4a: coefs4a %i for ray %d: got %g, exp %g dif %g\n", i,
                  nray, coefs4a[i], expectedcoefs4a[5*nray + i], ABS(coefs4a[i]- expectedcoefs4a[5*nray + i]));

            cr_expect_float_eq(coefs4b[i], expectedcoefs4b[5*nray + i], 1e-12, "case 4b: coefs4b %i for ray %d: got %g, exp %g dif %g\n", i,
                  nray, coefs4b[i], expectedcoefs4b[5*nray + i], ABS(coefs4b[i]- expectedcoefs4b[5*nray + i]));
         }

         /* check roots */
         printf("computing root with A, B, C, D, E = %g, %g, %g, %g, %g\n",coefs4a[0], coefs4a[1], coefs4a[2], coefs4a[3], coefs4a[4]);
         root = computeRoot(scip, coefs4a);

         /* high tolerance because we do a bin search to ensure that phi(root) <= 0 and small, we don't solve exactly */
         cr_expect_float_eq(root, expectedroots4a[nray], 1e-5, "case 4a: root for ray %d: got %g, exp %g\n", nray, root,
               expectedroots4a[nray]);
         root = computeRoot(scip, coefs4b);
         cr_expect_float_eq(root, expectedroots4b[nray], 1e-5, "case 4b: root for ray %d: got %g, exp %g dif %g\n", nray, root,
               expectedroots4b[nray], root - expectedroots4b[nray]);
      }
   }


   /* check some individual rays */
   {
      /* test ray {y, z, x} = {0, 0, 1}: there is finite intersection in case 4a, but infinte in 4b. overall, there is
       * no intersection
       */
      SCIP_Real testraycoef[1] = {1.0};
      int testrayidx[1] = {2};
      int testraynnonz = 1;
      SCIP_Real expectedroot4a = 2.0 + 4.0 * SQRT( 2.0 ) + 2.0 * SQRT( 10.0 + 6.0 * SQRT( 2.0 ) );
      SCIP_Real expectedroot4b = SCIPinfinity(scip);

      printf("testing ray with finite/infinte intersection\n");

      SCIP_CALL( computeRestrictionToRay(scip, nlhdlrexprdata, 1.0, TRUE, testraycoef, testrayidx, testraynnonz,
               vb, vzlp, wcoefs, wzlp, kappa, coefs4a, coefs4b, coefscond, &success) );
      cr_expect(success);
      root = computeRoot(scip, coefs4a);

      cr_expect_float_eq(root, expectedroot4a, 1e-12, "case 4a: root for custom ray: got %g, exp %g\n", root,
            expectedroot4a);

      root = computeRoot(scip, coefs4b);
      cr_expect_float_eq(root, expectedroot4b, 1e-12, "case 4b: root for custom ray: got %g, exp %g\n", root,
            expectedroot4b);
   }
   {
      /* test ray {y, z, x} = {0, 0.136472436, 1} */
      SCIP_Real testraycoef[2] = {34118109.0 / 250000000, 1.0};

      int testrayidx[2] = {1,2};
      int testraynnonz = 2;
      SCIP_Real expectedroot4a = 2.0249091949668 * 1e9;
      SCIP_Real expectedroot4b = SCIPinfinity(scip);

      printf("testing ray with large root\n");

      /* 4a */
      SCIP_CALL( computeRestrictionToRay(scip, nlhdlrexprdata, 1.0, TRUE, testraycoef, testrayidx, testraynnonz,
               vb, vzlp, wcoefs, wzlp, kappa, coefs4a, coefs4b, coefscond, &success) );
      cr_expect(success);
      root = computeRoot(scip, coefs4a);

      printf("computing root with A, B, C, D, E = %.15f, %.15f, %.15f, %.15f, %.15f\n",coefs4a[0], coefs4a[1], coefs4a[2], coefs4a[3], coefs4a[4]);
      /* TODO: the error that we have is around 300, implementing a binary search should help improve; however, 300 is
       * not going to make much of a difference for the cut coefficient */
      cr_expect_float_eq(root, expectedroot4a, 1000, "case 4a: root for custom ray: got %.15f, exp %.15f dif %g\n", root,
            expectedroot4a, ABS(root - expectedroot4a));

      /* 4b */
      root = computeRoot(scip, coefs4b);

      cr_expect_float_eq(root, expectedroot4b, 1e-12, "case 4b: root for custom ray: got %g, exp %g\n", root,
            expectedroot4b);
   }

   /* check cut */
   {
      int       expectednvars    = 3;
      SCIP_Real expectedcoefs[3] = {-1.63242263171138, -2.91415905378622, -5.07771472287715};
      SCIP_VAR* expectedvars[3]  = {x, y, z};
      SCIP_Real expectedlhs      = 1.00399143674532;

      testCut(nlhdlrexprdata, cons, FALSE, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }


   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

/* test when aux var is present and all nonbasic */
Test(interCuts, testRaysAuxvar1)
{
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_EXPR* expr;
   SCIP_VAR* auxvar;

   /* simplify and detect quadratic structure in: x - 6z + 2z^2 + 2 <= 0 */
   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: <x> - 6.0*<z> + 2.0*<z>^2 + 2 <= 0");

   /* create aux variable */
   expr = SCIPgetExprNonlinear(cons);
   cr_assert_not_null(expr);
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, expr, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( initSepa(scip, conshdlr, &cons, 1, &infeasible) );

   /*
    * build LP so that every var is non-basic
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   /* set bounds of auxvar */
   auxvar = SCIPgetExprAuxVarNonlinear(expr);
   cr_assert_not_null(auxvar);

   SCIP_CALL( SCIPchgVarObjProbing(scip, auxvar, 1.0) );
   SCIP_CALL( SCIPchgVarLbProbing(scip, auxvar, 0.0) );

   buildAndSolveSimpleProbingLP();

   /* rays should be
    * z       = 0 x + 0 y + 0 w + 1 z + 0 auxvar
    * x       = 1 x + 0 y + 0 w + 0 z + 0 auxvar
    * auxvar  = 0 x + 0 y + 0 w + 0 z + 1 auxvar
    */
   {
      int expectednnonz = 3;
      SCIP_Real expectedrayscoefs[3] = {1.0, 1.0, 1.0};
      int expectedraysidx[3] = {1, 0, 2};
      int expectednrays = 3;
      int expectedlppos[3] = {0, 3, 6};
      int expectedbegin[4] = {0, 1, 2, 3};

      int expectedncols = 7;
      SCIP_VAR* vars[7] = {x, y, w, z, s, t, auxvar};
      SCIP_VAR* consvars[3] = {z, x, auxvar};

      testRays(nlhdlrexprdata, auxvar, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

#if SCIP_DISABLED_CODE
/* Gurobi finds a different basis here since the optimal solution of this problem is degenerate.
 * We disable it for now since it fails for Gurobi. TODO: fix that or run only for specific lp solvers */
/* test when aux var is present and auxvar basic */
Test(interCuts, testRaysAuxvar2)
{
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_EXPR* expr;
   SCIP_VAR* auxvar;
   SCIP_Real* eigenvalues;

   /* simplify and detect quadratic structure in: x - 6z + 2z^2 + 2 <= 0*/
   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: <x> - 6.0*<z> - 2.0*<z>^2 + 2.0 <= 0.0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /* for this example simplify changes the cons to -x +6z +2z^2 -2
    * we are testing the constraint induced by the auxiliary variable:
    * -x +6z +2z^2 -2 = auxvar
    *  The solutin of the LP is 0, so the violated constraint is
    * -x +6z +2z^2 -2 >= auxvar
    *  thus we need to overestimate.
    */

   /* create aux variable */
   expr = SCIPgetExprNonlinear(cons);
   cr_assert_not_null(expr);
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, expr, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( initSepa(scip, conshdlr, &cons, 1, &infeasible) );

   //SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, expr, NULL) );
   auxvar = SCIPgetExprAuxVarNonlinear(expr);
   cr_assert_not_null(auxvar);

   /* build LP:
    * I want
    *      x = 2 z - 1.1 slack1 - 1.2 slack_2
    * auxvar = 1 z + 0.1 slack1 - 1.3 slack_2
    * so LP must be the one below
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   SCIP_ROW* row1;
   SCIP_ROW* row2;
   SCIP_Bool lperror;

   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row1, "row1", -SCIPinfinity(scip), 0.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, z,      -28.0 / 31) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, x,       26.0 / 31) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, auxvar, -24.0 / 31) );

   SCIP_CALL( SCIPaddRowProbing(scip, row1) );
   SCIP_CALL( SCIPreleaseRow(scip, &row1) );

   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row2, "row2", -SCIPinfinity(scip), 0.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row2, z,      -26.0 / 31) );
   SCIP_CALL( SCIPaddVarToRow(scip, row2, x,        2.0 / 31) );
   SCIP_CALL( SCIPaddVarToRow(scip, row2, auxvar,  22.0 / 31) );

   SCIP_CALL( SCIPaddRowProbing(scip, row2) );
   SCIP_CALL( SCIPreleaseRow(scip, &row2) );

   /* we want z to be nonbasic at lower */
   SCIP_CALL( SCIPchgVarLbProbing(scip, z, 0.0) );

   /* set objective to force x and aux to be basic */
   SCIP_CALL( SCIPchgVarObjProbing(scip, z,       1.0) );
   SCIP_CALL( SCIPchgVarObjProbing(scip, x,      -1.0/3) );
   SCIP_CALL( SCIPchgVarObjProbing(scip, auxvar,  1.0/4) );


   /* SCIP_CALL( SCIPwriteLP(scip, "probing.lp") ); */

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );

   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(z)), SCIP_BASESTAT_LOWER);
   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(x)), SCIP_BASESTAT_BASIC, "got %d\n", SCIPcolGetBasisStatus(SCIPvarGetCol(x)));
   cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(auxvar)), SCIP_BASESTAT_BASIC);

   /* rays should be
    *      x = 2 z - 1.1 slack1 - 1.2 slack_2
    * auxvar = 1 z + 0.1 slack1 - 1.3 slack_2
    * so
    * z       = 1 z +   0 slack1 +   0 slack2
    * x       = 2 z - 1.1 slack1 - 1.2 slack2
    * auxvar  = 1 z + 0.1 slack1 - 1.3 slack2
    */
   {
      int expectednnonz = 7;
      SCIP_Real expectedrayscoefs[7] = {1.0, 2.0, 1.0, -1.1, 0.1, -1.2, -1.3};
      int expectedraysidx[7] = {0, 1, 2, 1, 2, 1, 2};
      int expectednrays = 3;
      int expectedlppos[3] = {3, -1,  -2};
      int expectedbegin[4] = {0, 3, 5, 7};

      int expectedncols = 7;
      SCIP_VAR* vars[7] = {x, y, w, z, s, t, auxvar};
      SCIP_VAR* consvars[3] = {z, x, auxvar};

      testRays(nlhdlrexprdata, auxvar, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }

   {
      int       expectednvars    = 3;
      SCIP_Real expectedcoefs[3] = {-0.140028008402801, 0.980196058819607, -0.140028008402801};
      SCIP_VAR* expectedvars[3]  = {x, z, auxvar};
      SCIP_Real expectedlhs      = 0.280056016805602;

      nlhdlrexprdata->cons = NULL; /* force to use aux var */
      testCut(nlhdlrexprdata, cons, TRUE /* overestimate */, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}
#endif

Test(interCuts, cut1, .description = "test cut for Case 2")
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_Real* eigenvalues;

   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: 6.0*<x>*<y> + 2.0*<x>^2 - 2.0*<y>^2 + 2 <= 0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /*
    * build LP so that every var is non-basic
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   buildAndSolveSimpleProbingLP();

   /* check cut */
   {
      int       expectednvars    = 2;
      SCIP_Real expectedcoefs[2] = {0.6335517491618165545, 1.1838022718621540656};
      SCIP_VAR* expectedvars[2]  = {x, y};
      SCIP_Real expectedlhs      = 1.0;

      testCut(nlhdlrexprdata, cons, FALSE, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }

   SCIP_CALL( SCIPendProbing(scip) );


   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}
Test(interCuts, cut2, .description = "test cut for Case 1")
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_Real* eigenvalues;

   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: (<x> - <y>)^2 - <z>*<x> + <x> -<z>^2 <= -1.0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /*
    * build LP so that every var is non-basic
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   buildAndSolveSimpleProbingLP();


   /* check cut without strengthening */
   {
      int       expectednvars    = 2;
      SCIP_Real expectedcoefs[2] = {0.136963531501410, 0.99057609048405};
      SCIP_VAR* expectedvars[2]  = {y, z};
      SCIP_Real expectedlhs      = 0.98106166069077;

      SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/quadratic/usestrengthening", FALSE) );
      testCut(nlhdlrexprdata, cons, FALSE, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }

   SCIP_CALL( SCIPendProbing(scip) );


   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}
Test(interCuts, cut3, .description = "test cut for Case 3")
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_Real* eigenvalues;

   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: -(<x> - <y>)^2 + (2 + <x> - <w>)^2 + (<x> + <y> + <z>)^2 - (<w> + <z>)^2 <= 0.25");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /*
    * build LP so that every var is non-basic
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   buildAndSolveSimpleProbingLP();


   /* check cut */
   {
      int       expectednvars    = 4;
      SCIP_Real expectedcoefs[4] = {0.0902290773616422, 0.389427393655165, 0.874367751710887, 0.275110983854626};
      SCIP_VAR* expectedvars[4]  = {x, y, w, z};
      SCIP_Real expectedlhs      = 0.788874349522027;

      SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/quadratic/usestrengthening", FALSE) );
      testCut(nlhdlrexprdata, cons, FALSE, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }

   SCIP_CALL( SCIPendProbing(scip) );


   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

/*
 * Tests for Strengthening
 */
Test(interCuts, strength1, .description = "test strengthening case 1")
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_Real* eigenvalues;

   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: (<x> - <y>)^2 -<z>*<x> + <x> + 1 - <z>^2 <= 0.0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /* build LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   buildAndSolveSimpleProbingLP2();

   /* what rays should be
    * x      1     0     0
    * y  =   0 s1  1 s2  1 s3
    * z      0     2     0
    */
   {
      int expectednnonz = 4;
      SCIP_Real expectedrayscoefs[4] = {
         1.0,
         1.0, 2.0,
         1.0
      };
      int expectedraysidx[4] = {
         0,
         1, 2,
         1
      };
      int expectednrays = 3;
      int expectedlppos[3] = {-1, -2, -3};
      int expectedbegin[4] = {0, 1, 3, 4};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[3] = {x, y, z};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }


   /* check cut */
   {
      int       expectednvars    = 3;
      SCIP_Real expectedcoefs[3] = {-0.4241094145918607, 0.5633142925617439, 0.7090897067721483};
      SCIP_VAR* expectedvars[3]  = {x, y, z};
      SCIP_Real expectedlhs      = 0.43067995082696664;

      SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/quadratic/usestrengthening", TRUE) );
      testCut(nlhdlrexprdata, cons, FALSE, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }

   SCIP_CALL( SCIPendProbing(scip) );


   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

Test(interCuts, strength2, .description = "test strengthening case 2")
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_Real* eigenvalues;

   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: (<x> - <y>)^2 - 2*<z>*<x> - <x> - <z>^2 <= -3.0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /* build LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   buildAndSolveSimpleProbingLP2();

   /* what rays should be
    * x      1     0     0
    * y  =   0 s1  1 s2  1 s3
    * z      0     2     0
    */
   {
      int expectednnonz = 4;
      SCIP_Real expectedrayscoefs[4] = {
         1.0,
         1.0, 2.0,
         1.0
      };
      int expectedraysidx[4] = {
         0,
         1, 2,
         1
      };
      int expectednrays = 3;
      int expectedlppos[3] = {-1, -2, -3};
      int expectedbegin[4] = {0, 1, 3, 4};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[3] = {x, y, z};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }


   /* check cut */
   {
      int       expectednvars    = 3;
      SCIP_Real expectedcoefs[3] = {-0.02813776917115283967, 0.6065759643564169028, 0.7945274478651783946};
      SCIP_VAR* expectedvars[3]  = {x, y, z};
      SCIP_Real expectedlhs      = 0.6638190109054805748;

      SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/quadratic/usestrengthening", TRUE) );
      /* we need to overestimate because simplify will multiply cons by -1 */
      testCut(nlhdlrexprdata, cons, TRUE, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }

   SCIP_CALL( SCIPendProbing(scip) );


   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

Test(interCuts, strength3, .description = "test strengthening case 3")
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_Real* eigenvalues;

   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: (<x> - <y>)^2 -<z>*<x> + <x> - <z>^2 <= 0.0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /* build LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   buildAndSolveSimpleProbingLP2();

   /* what rays should be
    * x      1     0     0
    * y  =   0 s1  1 s2  1 s3
    * z      0     2     0
    */
   {
      int expectednnonz = 4;
      SCIP_Real expectedrayscoefs[4] = {
         1.0,
         1.0, 2.0,
         1.0
      };
      int expectedraysidx[4] = {
         0,
         1, 2,
         1
      };
      int expectednrays = 3;
      int expectedlppos[3] = {-1, -2, -3};
      int expectedbegin[4] = {0, 1, 3, 4};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[3] = {x, y, z};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }


   /* check cut */
   {
      int       expectednvars    = 3;
      SCIP_Real expectedcoefs[3] = {-0.4180242003191831481, 0.6139777515027494797, 0.6695424472034132719};
      SCIP_VAR* expectedvars[3]  = {x, y, z};
      SCIP_Real expectedlhs      = 0.1747724718238416967;

      SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/quadratic/usestrengthening", TRUE) );
      testCut(nlhdlrexprdata, cons, FALSE, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }

   SCIP_CALL( SCIPendProbing(scip) );


   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

Test(interCuts, strength4, .description = "test strengthening case 4")
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_Real* eigenvalues;

   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: <x>^2 -<z>*<x> + <y>*<x> + <x>- <z> + 2.0 <= 0.0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /* build LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   buildAndSolveSimpleProbingLP2();

   /* what rays should be
    * x      1     0     0
    * y  =   0 s1  1 s2  1 s3
    * z      0     2     0
    */
   {
      int expectednnonz = 4;
      SCIP_Real expectedrayscoefs[4] = {
         1.0,
         1.0, 2.0,
         1.0
      };
      int expectedraysidx[4] = {
         0,
         1, 2,
         1
      };
      int expectednrays = 3;
      int expectedlppos[3] = {-1, -2, -3};
      int expectedbegin[4] = {0, 1, 3, 4};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[3] = {x, y, z};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }


   /* check cut */
   {
      int       expectednvars    = 3;
      SCIP_Real expectedcoefs[3] = {-0.21313592330512526, 0.10270409145863757, 0.9716094625900509};
      SCIP_VAR* expectedvars[3]  = {x, y, z};
      SCIP_Real expectedlhs      = 1.3375271971557974;

      SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/quadratic/usestrengthening", TRUE) );
      testCut(nlhdlrexprdata, cons, FALSE, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }

   SCIP_CALL( SCIPendProbing(scip) );


   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

/* The first ray and third ray are in the recession cone of C
 * The second ray intersects in 4a
 * After strengthening the third ray, the new ray is in the recession cone of C due to 4b
 */
Test(interCuts, strength4ab, .description = "more complicated test strengthening case 4")
{
   SCIP_Bool cutoff;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_Real* eigenvalues;

   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: <x>^2 - 3*<y>^2 + 3*<x>*<y> + <z> + 1.0 <= 0.0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);
   SCIP_CALL( SCIPstartProbing(scip) );

   /* build LP
    * 7 x/18 - 8 y/9 - z/2 <= -1/9,
    * x/2 - z/2 <= 0,
    * -2 x/9 + 2 y/9 <= -2/9
    */
   {
      SCIP_ROW* row1;
      SCIP_ROW* row2;
      SCIP_ROW* row3;
      SCIP_Bool lperror;

      /* 7 x/18 - 8 y/9 - z/2 <= -1/9 */
      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row1, "row1", -SCIPinfinity(scip), -1.0 / 9, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, x,   7.0 / 18) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, y,  -8.0 / 9) );
      SCIP_CALL( SCIPaddVarToRow(scip, row1, z,  -1.0 / 2) );

      SCIP_CALL( SCIPaddRowProbing(scip, row1) );
      SCIP_CALL( SCIPreleaseRow(scip, &row1) );

      /* x/2 -z/2 <= 0 */
      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row2, "row2", -SCIPinfinity(scip), 0.0, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarToRow(scip, row2, x,   1.0 / 2) );
      SCIP_CALL( SCIPaddVarToRow(scip, row2, z,  -1.0 / 2) );

      SCIP_CALL( SCIPaddRowProbing(scip, row2) );
      SCIP_CALL( SCIPreleaseRow(scip, &row2) );

      /* 2 x/9 + 2 y/9 <= -2/9 */
      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row3, "row3", -SCIPinfinity(scip), -2.0 / 9, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarToRow(scip, row3, x,  -2.0 / 9) );
      SCIP_CALL( SCIPaddVarToRow(scip, row3, y,   2.0 / 9) );

      SCIP_CALL( SCIPaddRowProbing(scip, row3) );
      SCIP_CALL( SCIPreleaseRow(scip, &row3) );

      /* set objective coefficient */
      SCIP_CALL( SCIPchgVarObjProbing(scip, x, 0.0) );
      SCIP_CALL( SCIPchgVarObjProbing(scip, y, 0.0) );
      SCIP_CALL( SCIPchgVarObjProbing(scip, z, 1.0) );


      /* SCIP_CALL( SCIPwriteLP(scip, "probing.lp") ); */

      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
      cr_assert_not(lperror);
      /* interestingly this fails because the cutoffbound of scip is .0001 or something TODO: figure out why */
      //cr_assert_not(cutoff);
      SCIP_CALL( SCIPprintSol(scip, NULL, NULL, TRUE) );

      /* all variables should be nonbasic */
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(x)), SCIP_BASESTAT_BASIC, "got %d\n", SCIPcolGetBasisStatus(SCIPvarGetCol(x)));
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(y)), SCIP_BASESTAT_BASIC);
      cr_expect_eq(SCIPcolGetBasisStatus(SCIPvarGetCol(z)), SCIP_BASESTAT_BASIC);

   }

   /* what rays should be
    * x      1     -1        4
    * y  =   1 s1  -1 s2  -1/2 s3
    * z      1      1        4
    */
   {
      int expectednnonz = 9;
      SCIP_Real expectedrayscoefs[9] = {
         1.0, 1.0, 1.0,
         -1.0, -1.0, 1.0,
         4.0, -0.5, 4.0
      };
      int expectedraysidx[9] = {
         0, 1, 2,
         0, 1, 2,
         0, 1, 2
      };
      int expectednrays = 3;
      int expectedlppos[3] = {-1, -2, -3};
      int expectedbegin[4] = {0, 3, 6, 9};

      int expectedncols = 6;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_VAR* consvars[3] = {x, y, z};

      testRays(nlhdlrexprdata, NULL, vars, consvars, expectedrayscoefs, expectedraysidx, expectedlppos, expectedbegin,
            expectednrays, expectednnonz, expectedncols);
   }


   /* check cut */
   {
      int       expectednvars    = 2;
      SCIP_Real expectedcoefs[2] = {-0.191886095440581, -0.981417203016418};
      SCIP_VAR* expectedvars[2]  = {x, y};
      SCIP_Real expectedlhs      = 0.551754189836395;

      SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/quadratic/usestrengthening", TRUE) );
      testCut(nlhdlrexprdata, cons, FALSE, expectedcoefs, expectedvars, expectedlhs, expectednvars);
   }

   SCIP_CALL( SCIPendProbing(scip) );


   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}

/* test that stored rays are
 * x    1     0     0     0
 * y  = 0 x + 1 y + 0 w + 0 z
 * z    0     0     0     1
 */
Test(interCuts, testBoundRays1)
{
   SCIP_Bool cutoff;
   SCIP_Bool lperror;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_SOL* vertex;
   SCIP_Real* eigenvalues;

   /* simplify and detect quadratic structure in: 6xy + 2x^2 -2z^2 + 2 <= 0 */
   simplifyAndDetect(&cons, &nlhdlrexprdata, "[nonlinear] <test>: 6.0*<x>*<y> + 2.0*<x>^2 - 2.0*<z>^2 + 2 <= 0");

   /* skip if eigenvalue decomposition does not exist */
   SCIPexprGetQuadraticData(nlhdlrexprdata->qexpr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      registerAndFree(cons, nlhdlrexprdata);
      return;
   }

   /*
    * build LP
    */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* the nonlinear constraint was not added, so it shouldn't add weird constraints to LP */
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   /* add bounds to vars */
   SCIPchgVarLbProbing(scip, x, -1.0); SCIPchgVarUbProbing(scip, x, 10.0);
   SCIPchgVarLbProbing(scip, y, -2.0); SCIPchgVarUbProbing(scip, y, 9.0);
   SCIPchgVarLbProbing(scip, z, 1.0);

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
   //cr_expect_not(cutoff); /* cutoff == TRUE since the LP is solved to optimality with a obtimal value greater than cutoffbound */
   cr_expect_not(lperror);

   /* choose solution to separate */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 4.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 6.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );

   /* the vertex we separate is {-1.0, 9.0, -, 1.0, -, -}
    * and bound rays are
    * x    1      0     0
    * y  = 0 x + -1 y + 0 z
    * z    0      0     1
    */
   {
      SCIP_Real expectedvertexcoefs[6] = {-1.0, 9.0, 0.0, 1.0, 0.0, 0.0};
      SCIP_Real expectedrayscoefs[8] = {1.0, -1.0, 1.0};
      int expectednrays = 3;
      int expectednnonz = 3;
      SCIP_VAR* vars[6] = {x, y, w, z, s, t};
      SCIP_Bool success;

      /*
      * find and test nearest vertex and rays
      */
      SCIP_CALL( SCIPcreateSol(scip, &vertex, NULL) );
      SCIP_CALL( findVertexAndGetRays(scip, nlhdlrexprdata, sol, vertex, NULL, &myrays, &success) );

      for( int i = 0; i < 6; ++i )
      {
         cr_expect_float_eq(SCIPgetSolVal(scip, vertex, vars[i]), expectedvertexcoefs[i], 1e-9, "%d-th entry: expected %g, got %g\n", i,
            expectedvertexcoefs[i], SCIPgetSolVal(scip, vertex, vars[i]));
      }

      cr_expect_eq(myrays->nrays, expectednrays, "e %d g %d\n", expectednrays, myrays->nrays);
      cr_expect_eq(myrays->raysbegin[myrays->nrays], expectednnonz, "e %d g %d\n", expectednnonz, myrays->raysbegin[myrays->nrays]);

      for( int i = 0; i < expectednnonz; ++i )
      {
         cr_expect_float_eq(myrays->rays[i], expectedrayscoefs[i], 1e-9, "%d-th entry: expected %g, got %g\n", i,
            expectedrayscoefs[i], myrays->rays[i]);
      }
      /* cr_expect_arr_eq(myrays->raysidx, expectedraysidx, expectednnonz * sizeof(int));
      cr_expect_arr_eq(myrays->lpposray, expectedlppos, expectednrays * sizeof(int));
      cr_expect_arr_eq(myrays->raysbegin, expectedbegin, (expectednrays + 1) * sizeof(int)); */

      SCIP_CALL( SCIPfreeSol(scip, &vertex) );
   }

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   /* register enforcer info in expr and free */
   registerAndFree(cons, nlhdlrexprdata);
}
