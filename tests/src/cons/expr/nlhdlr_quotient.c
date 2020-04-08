/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_quotient.c
 * @brief  tests quotient nonlinear handler methods
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr.c"
#include "scip/cons_expr_nlhdlr_quotient.c"


/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* w;
static SCIP_CONSEXPR_NLHDLR* nlhdlr;
static SCIP_CONSHDLR* conshdlr;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers and nonlinear handlers; get quadratic handler and conshdlr */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   nlhdlr = SCIPfindConsExprNlhdlr(conshdlr, "quotient");
   cr_assert_not_null(nlhdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 1.5, 5.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -4.0, 0.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 1.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -4.0, -1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPaddVar(scip, w) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   BMScheckEmptyMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* test suite */
TestSuite(nlhdlrquotient, .init = setup, .fini = teardown);

/** checks whether the values in nlhdlrexprdata are as expected */
static
void checkData(
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< the nlhdlr expression data */
   SCIP_VAR*             numvar,             /**< expected numerator variable */
   SCIP_Real             numcoef,            /**< expected numerator coefficient */
   SCIP_Real             numconst,           /**< expected numerator constant */
   SCIP_VAR*             denomvar,           /**< expected denominator variable */
   SCIP_Real             denomcoef,          /**< expected denominator coefficient */
   SCIP_Real             denomconst,         /**< expected denominator constant */
   SCIP_Real             constant            /**< expected constant */
   )
{
   cr_expect_not_null(nlhdlrexprdata->numvar);
   cr_expect_not_null(nlhdlrexprdata->denomvar);
   cr_expect(SCIPisEQ(scip, numcoef, nlhdlrexprdata->numcoef));
   cr_expect(SCIPisEQ(scip, numconst, nlhdlrexprdata->numconst));
   cr_expect(SCIPisEQ(scip, denomcoef, nlhdlrexprdata->denomcoef));
   cr_expect(SCIPisEQ(scip, denomconst, nlhdlrexprdata->denomconst));
   cr_expect(SCIPisEQ(scip, constant, nlhdlrexprdata->constant));
}

/* checks whether the values of the cut are as expected */
static
void checkCut(
   SCIP_ROW*             cut,                /**< the cut to check */
   SCIP_VAR**            vars,               /**< array of expected variables */
   SCIP_Real*            vals,               /**< array of expected coefficients */
   SCIP_Real             lhs,                /**< expected lhs value */
   SCIP_Real             rhs,                /**< expected rhs value */
   int                   nvars               /**< expected number of variables */
   )
{
   SCIP_Bool* foundvar;
   int i;
   int j;


   cr_assert_not_null(cut);
   cr_assert_not_null(vars);
   cr_assert_not_null(vals);

   (void) SCIPprintRow(scip, cut, NULL);
   cr_expect_eq(SCIProwGetNNonz(cut), nvars, "expected %d variable(s) in cut, but got %d\n",
      nvars, SCIProwGetNNonz(cut));
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(cut), lhs), "expected lhs = %f, but got %f\n", lhs, SCIProwGetLhs(cut));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(cut), rhs), "expected rhs = %f, but got %f\n", rhs, SCIProwGetRhs(cut));

   SCIP_CALL( SCIPallocClearBufferArray(scip, &foundvar, nvars) );

   for( i = 0; i < nvars; ++i )
   {
      for( j = 0; j < nvars; ++j )
      {
         if( SCIPcolGetVar(SCIProwGetCols(cut)[j]) == vars[i] )
         {
            foundvar[i] = TRUE;

            cr_expect(SCIPisEQ(scip, SCIProwGetVals(cut)[j], vals[i]), "expected coef of %s to be %f, but got %f\n",
               SCIPvarGetName(vars[i]), vals[i], SCIProwGetVals(cut)[j]);

            break;
         }
      }

      cr_expect(foundvar[i], "did not find variable %s in the cut", SCIPvarGetName(vars[i]));
   }

   SCIPfreeBufferArray(scip, &foundvar);
}

/* detects x / y */
Test(nlhdlrquotient, detectandfree1, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: <x> / <y> <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 1.0, 0.0, y, 1.0, 0.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects (4x + 1) / (-3x - 3) */
Test(nlhdlrquotient, detectandfree2, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: (4*<x> + 1) / (-3*<x> - 3) <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 4.0, 1.0, x, -3.0, -3.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects log((4x + 3) / (x + 1)) */
Test(nlhdlrquotient, detectandfree3, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: log((4*<x> + 3) / (<x> + 1)) <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons)->children[0];

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 4.0, 3.0, x, 1.0, 1.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects that (4x + 2y + 3) / (x + 1)) is invalid */
Test(nlhdlrquotient, detectandfree4, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: (4*<x> + 2*<y> + 3) / (<x> + 1) <= 10",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_null(nlhdlrexprdata);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects log(x) / log(y)) */
Test(nlhdlrquotient, detectandfree5, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR* auxvar1;
   SCIP_VAR* auxvar2;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: log(<x>) / abs(<y>) <= 10",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   auxvar1 = SCIPgetConsExprExprAuxVar(expr->children[0]);
   auxvar2 = SCIPgetConsExprExprAuxVar(expr->children[1]->children[0]);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, auxvar1, 1.0, 0.0, auxvar2, 1.0, 0.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects (4x + 1) / (-3x + 3) + 2 after simplification */
Test(nlhdlrquotient, detectandfree6, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: ((4*<x> + 1) / (-3*<x> - 3) + 2) <= 3",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* adding a constraint also adds locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_expect_not(infeasible);

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 4.0, 1.0, x, -3.0, -3.0, 2.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* tests interval evaluation for ((+/-)4x + 1) / (-3x + 3) - 2 */
Test(nlhdlrquotient, inteval, .description = "tests interval evaluation of simple quotient expression")
{
   SCIP_INTERVAL varbnds;
   SCIP_INTERVAL result;

   /* test interval including 0 in denominator*/

   varbnds.inf = 0.0;
   varbnds.sup = 2.0;

   result = intEvalQuotient(scip, varbnds, 4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, result));

   /* test positive denominator part for monotone increasing expression */

   varbnds.inf = 2.0;
   varbnds.sup = 9.0;

   result = intEvalQuotient(scip, varbnds, 4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, -5.0));
   cr_expect(SCIPisEQ(scip, result.sup, -37.0 / 24.0 - 2.0), "expected %f, but got %f\n",
      -37.0 / 24.0 - 2.0, result.sup);

   /* test negative denominator part for monotone increasing expression */

   varbnds.inf = -1.0;
   varbnds.sup = 0.0;

   result = intEvalQuotient(scip, varbnds, 4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, -2.5));
   cr_expect(SCIPisEQ(scip, result.sup, 1.0 / 3.0 - 2.0));


   /* test positive denominator part for monotone decreasing expression */

   varbnds.inf = 2.0;
   varbnds.sup = 9.0;

   result = intEvalQuotient(scip, varbnds, -4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, 35.0 / 24.0 - 2.0));
   cr_expect(SCIPisEQ(scip, result.sup, 7.0 / 3.0 - 2.0));

   /* test negative denominator part for monotone decreasing expression */

   varbnds.inf = -1.0;
   varbnds.sup = 0.0;

   result = intEvalQuotient(scip, varbnds, -4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, 1.0 / 3.0 - 2.0));
   cr_expect(SCIPisEQ(scip, result.sup, 5.0 / 6.0 - 2.0));
}

/* tests reverse propagation for univariate quotients */
Test(nlhdlrquotient, reverseprop, .description = "tests reverse propagation simple univariate quotient expressions")
{
   SCIP_INTERVAL bnds;
   SCIP_INTERVAL result;

   /* x / (x + 1) in [-3,-1] => x in [-0.75,0.5]*/
   SCIPintervalSetBounds(&bnds, -3.0, -1.0);
   result = reversepropQuotient(bnds, 1.0, 0.0, 1.0, 1.0, 0.0);
   cr_expect(SCIPisEQ(scip, result.inf, -0.75));
   cr_expect(SCIPisEQ(scip, result.sup, -0.5));

   /* x / (x + 1) in [-2,1] => x in [-2/3,+inf]*/
   SCIPintervalSetBounds(&bnds, -2.0, 1.0);
   result = reversepropQuotient(bnds, 1.0, 0.0, 1.0, 1.0, 0.0);
   cr_expect(SCIPisEQ(scip, result.inf, -2.0/3.0));
   cr_expect(SCIPisInfinity(scip, result.sup));

   /* (-5x + 2) / (3*x + 3) + 6 in [3,5] => x in [-inf,+inf]*/
   SCIPintervalSetBounds(&bnds, 3.0, 5.0);
   result = reversepropQuotient(bnds, -5.0, 2.0, 3.0, 3.0, 6.0);
   cr_expect(SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, result));

   /* (-5x + 2) / (3*x + 3) + 6 in [-2,-1] => x in [-23/16,-26/19]*/
   SCIPintervalSetBounds(&bnds, -2.0, -1.0);
   result = reversepropQuotient(bnds, -5.0, 2.0, 3.0, 3.0, 6.0);
   cr_expect(SCIPisEQ(scip, result.inf, -23.0/16.0));
   cr_expect(SCIPisEQ(scip, result.sup, -26.0/19.0));
}

/* separates x = 2 for (4x + 1) / (-3x + 3) + 2 */
Test(nlhdlrquotient, separation1, .description = "separates simple univariate quotient expression")
{
   SCIP_ROWPREP* cutprep;
   SCIP_ROW* cut;
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_SOL* sol;
   SCIP_VAR* auxvar;
   SCIP_VAR* cutvar;
   SCIP_Real cutcoef;
   SCIP_Real tmpinf;
   SCIP_Real tmpsup;
   SCIP_Bool infeasible;
   SCIP_Bool success;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: ((4*<x> + 1) / (-3*<x> + 3) - 2) <= 3",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_expect_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* create auxvar for expression */
   SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, expr, &auxvar) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );

   /**
    *  test overestimation
    */

   SCIP_CALL( SCIPcreateRowprep(scip, &cutprep, SCIP_SIDETYPE_LEFT, TRUE) );
   SCIP_CALL( sepaUnivariate(scip, sol, x, 4.0, 1.0, -3.0, 3.0, -2.0, TRUE, cutprep, &success) );

   cr_assert(success);
   cr_assert_not_null(cutprep);

   SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, cutprep, cons) );

   cutvar = x;
   cutcoef = 5.0 / 3.0;

   checkCut(cut, &cutvar, &cutcoef, 5.0 + cutcoef * 2.0, SCIPinfinity(scip), 1);

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   SCIPfreeRowprep(scip, &cutprep);

   /*
    *  test underestimation
    */

   SCIP_CALL( SCIPcreateRowprep(scip, &cutprep, SCIP_SIDETYPE_RIGHT, TRUE) );
   SCIP_CALL( sepaUnivariate(scip, sol, x, 4.0, 1.0, -3.0, 3.0, -2.0, FALSE, cutprep, &success) );

   cr_assert(success);
   cr_assert_not_null(cutprep);

   SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, cutprep, cons) );

   tmpinf = -7.0 / 1.5 - 2.0;
   tmpsup = -21.0 / 12.0 - 2.0;
   cutvar = x;
   cutcoef = (tmpsup - tmpinf) / 3.5;
   checkCut(cut, &cutvar, &cutcoef, -SCIPinfinity(scip), -tmpsup + cutcoef * 5.0, 1);

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   SCIPfreeRowprep(scip, &cutprep);

   /* free the rest */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* separates y = -1 for (4y + 1) / (-3y + 3) + 2 */
Test(nlhdlrquotient, separation2, .description = "separates simple univariate quotient expression")
{
   SCIP_ROWPREP* cutprep;
   SCIP_ROW* cut;
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_SOL* sol;
   SCIP_VAR* auxvar;
   SCIP_VAR* cutvar;
   SCIP_Real cutcoef;
   SCIP_Real tmpinf;
   SCIP_Real tmpsup;
   SCIP_Bool infeasible;
   SCIP_Bool success;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: ((4*<y> + 1) / (-3*<y> + 3) - 2) <= 3",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_expect_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* create auxvar for expression */
   SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, expr, &auxvar) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -1.0) );

   /*
    *  test overestimation
    */

   SCIP_CALL( SCIPcreateRowprep(scip, &cutprep, SCIP_SIDETYPE_LEFT, TRUE) );
   SCIP_CALL( sepaUnivariate(scip, sol, y, 4.0, 1.0, -3.0, 3.0, -2.0, TRUE, cutprep, &success) );

   cr_assert(success);
   cr_assert_not_null(cutprep);

   SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, cutprep, cons) );

   tmpinf = -3.0;
   tmpsup = 1.0 / 3.0 - 2.0;
   cutvar = y;
   cutcoef = (tmpsup - tmpinf) / 4.0;

   checkCut(cut, &cutvar, &cutcoef, -tmpsup, SCIPinfinity(scip), 1);

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   SCIPfreeRowprep(scip, &cutprep);

   /*
    *  test underestimation
    */

   SCIP_CALL( SCIPcreateRowprep(scip, &cutprep, SCIP_SIDETYPE_RIGHT, TRUE) );
   SCIP_CALL( sepaUnivariate(scip, sol, y, 4.0, 1.0, -3.0, 3.0, -2.0, FALSE, cutprep, &success) );

   cr_assert(success);
   cr_assert_not_null(cutprep);

   SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, cutprep, cons) );

   cutvar = y;
   cutcoef = 5.0 / 12.0;
   checkCut(cut, &cutvar, &cutcoef, -SCIPinfinity(scip), 2.5 - cutcoef, 1);

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   SCIPfreeRowprep(scip, &cutprep);

   /* free the rest */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* separates (x,y) = (3,2) for x/y for x in [1,4] and z in [1.5,5] */
Test(nlhdlrquotient, separation3, .description = "separates simple bivariate quotient expression")
{
   SCIP_Real vals[3];
   SCIP_Bool success;

   /*
    * test overestimation
    */

   SCIP_CALL( estimateBivariate(scip, 1.0, 4.0, 1.5, 5.0, -SCIPinfinity(scip), SCIPinfinity(scip), 3.0, 2.0, 0.0, TRUE,
      &vals[0], &vals[1], &vals[2], &success) );
   cr_expect(success);
   cr_expect(SCIPisEQ(scip, vals[0], 2.0 / 3.0), "got %g expected %g", vals[0], 2.0 / 3.0);
   cr_expect(SCIPisEQ(scip, vals[1], -1.0 / 7.5), "got %g expected %g", vals[1], -1.0 / 7.5);
   cr_expect(SCIPisEQ(scip, vals[2], 1.0 / 5.0), "got %g expected %g", vals[2], 1.0 / 5.0);

   /*
    * test underestimation
    */

   SCIP_CALL( estimateBivariate(scip, 1.0, 4.0, 1.5, 5.0, -SCIPinfinity(scip), SCIPinfinity(scip), 3.0, 2.0, 0.0, FALSE,
      &vals[0], &vals[1], &vals[2], &success) );
   cr_expect(success);
   cr_expect(SCIPisEQ(scip, vals[0], 5.0 / 9.0), "got %g expected %g", vals[0], 5.0 / 9.0);
   cr_expect(SCIPisEQ(scip, vals[1], -25.0 / 36.0), "got %g expected %g", vals[1], -25.0 / 36.0);
   cr_expect(SCIPisEQ(scip, vals[2], 10.0 / 9.0), "got %g expected %g", vals[2], 10.0 / 9.0);
}

/* separates (x,y) = (-3,2) for x/y for x in [-4,-1] and y in [1.5,5] */
Test(nlhdlrquotient, separation4, .description = "separates simple bivariate quotient expression")
{
   SCIP_Real vals[3];
   SCIP_Bool success;

   /*
    * test overestimation
    */

   SCIP_CALL( estimateBivariate(scip, -4.0, -1.0, 1.5, 5.0, -SCIPinfinity(scip), SCIPinfinity(scip), -3.0, 2.0, 0.0, TRUE,
      &vals[0], &vals[1], &vals[2], &success) );
   cr_expect(success);
   cr_expect(SCIPisEQ(scip, vals[0], 5.0 / 9.0), "got %g expected %g", vals[0], 5.0 / 9.0);
   cr_expect(SCIPisEQ(scip, vals[1], 25.0 / 36.0), "got %g expected %g", vals[1], 25.0 / 36.0);
   cr_expect(SCIPisEQ(scip, vals[2], -10.0 / 9.0), "got %g expected %g", vals[2], -10.0 / 9.0);

   /*
    * test underestimation
    */

   SCIP_CALL( estimateBivariate(scip, -4.0, -1.0, 1.5, 5.0, -SCIPinfinity(scip), SCIPinfinity(scip), -3.0, 2.0, 0.0, FALSE,
      &vals[0], &vals[1], &vals[2], &success) );
   cr_expect(success);
   cr_expect(SCIPisEQ(scip, vals[0], 2.0 / 3.0), "got %g expected %g", vals[0], 2.0 / 3.0);
   cr_expect(SCIPisEQ(scip, vals[1], 2.0 / 15.0), "got %g expected %g", vals[1], 2.0 / 15.0);
   cr_expect(SCIPisEQ(scip, vals[2], -1.0 / 5.0), "got %g expected %g", vals[2], -1.0 / 5.0);
}

/* TODO add a test for the case that 0 is in the domain of the numerator */