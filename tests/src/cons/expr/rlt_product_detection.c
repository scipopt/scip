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

/**@file   rlt.c
 * @brief  tests hidden product detection for rlt cuts
 * @author Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_DEBUG

#include "scip/sepa_rlt.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SEPA* sepa;
static SCIP_VAR* x1o;
static SCIP_VAR* x2o;
static SCIP_VAR* x3o;
static SCIP_VAR* x4o;
static SCIP_VAR* y12o;
static SCIP_VAR* b1o;
static SCIP_VAR* b2o;
static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* x3;
static SCIP_VAR* x4;
static SCIP_VAR* y12;
static SCIP_VAR* b1;
static SCIP_VAR* b2;
static SCIP_VAR* auxvar;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   SCIP_CALL( SCIPincludeConshdlrLinear(scip) );

   /* include rlt separator */
   SCIP_CALL( SCIPincludeSepaRlt(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* get rlt separator */
   sepa = SCIPfindSepa(scip, "rlt");
   assert(sepa != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1o, "x1", -1.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2o, "x2", -6.0, -3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3o, "x3", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4o, "x4", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &b1o, "b1", 0, 1, 1, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &b2o, "b2", 0, 1, 1, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y12o, "y12", 2.0, 4.0, -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x1o) );
   SCIP_CALL( SCIPaddVar(scip, x2o) );
   SCIP_CALL( SCIPaddVar(scip, x3o) );
   SCIP_CALL( SCIPaddVar(scip, x4o) );
   SCIP_CALL( SCIPaddVar(scip, b1o) );
   SCIP_CALL( SCIPaddVar(scip, b2o) );
   SCIP_CALL( SCIPaddVar(scip, y12o) );

   /* get SCIP into SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* get transformed vars and release vars */
   SCIP_CALL( SCIPgetTransformedVar(scip, x1o, &x1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x2o, &x2) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x3o, &x3) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x4o, &x4) );
   SCIP_CALL( SCIPgetTransformedVar(scip, b1o, &b1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, b2o, &b2) );
   SCIP_CALL( SCIPgetTransformedVar(scip, y12o, &y12) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4o) );
   SCIP_CALL( SCIPreleaseVar(scip, &b1o) );
   SCIP_CALL( SCIPreleaseVar(scip, &b2o) );
   SCIP_CALL( SCIPreleaseVar(scip, &y12o) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

Test(rlt_selection, implrels, .init = setup, .fini = teardown, .description = "test extracting products from two implied relations")
{
   SCIP_ROW* row1;
   SCIP_ROW* row2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   SCIP_SEPADATA* sepadata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR* vars[3];
   SCIP_Real coefs1[3], coefs2[3];
   int c;
   SCIP_Bool cutoff;
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_CONSEXPR_EXPR** linchildren;
   SCIP_Real* lincoefs;

   vars[0] = x1; /* w or y */
   vars[1] = x2; /* y or w */
   vars[2] = b1; /* x */

   coefs1[0] = 2.0; /* a1 or c1 */
   coefs1[1] = 1.5; /* c1 or a1 */
   coefs1[2] = 3.0; /* b1 */

   coefs2[0] = 4.0; /* a2 or c2 */
   coefs2[1] = 2.0; /* c2 or a2 */
   coefs2[2] = -1.0; /* b2 */

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   sepadata->isinitialround = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear constraints */

   /* 0 <= 2x1 + 1.5x2 + 3b1 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons1, "c1", 3, vars, coefs1, 0.0, 1.0) );

   /* 4x1 + 2x2 - b1 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons2, "c2", 3, vars, coefs2, -SCIPinfinity(scip), 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons1) );
   SCIP_CALL( SCIPaddCons(scip, cons2) );

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* create separator data - this will also detect the products */
   SCIP_CALL( createSepaData(scip, sepadata) );

   /* the product should be: xy >= (1/(a1c2 - c1a2))*(a1a2w + (a2(b1 - d1) + a1d2)x + a1c2y - a1d2)
    * or, substituting the given values: xy <= (1/(2*2 - 1.5*4))*(2*4w + (4(3 - 1) + 2*1)x + 2*2y - 2*1)
    * xy <= -4w - 5x - 2y + 1
    * the second product:
    * xy >= (1/(1.5*4 - 2*2))*(1.5*2w + (2(3 - 1) + 1.5*1)x + 1.5*4y - 1.5*1)
    * xy >= 1.5w + 2.75x + 3y - 0.75 */

   /* check the numbers */
   cr_expect_eq(sepadata->nbilinvars, 3, "\nExpected 3 bilinear vars, got %d", sepadata->nbilinvars);
   cr_expect_eq(sepadata->nbilinterms, 2, "\nExpected 2 bilinear terms, got %d", sepadata->nbilinterms);
   cr_expect_eq(sepadata->nlinexprs[0], 4, "\nExpected 4 linear expressions for product 0, got %d", sepadata->nlinexprs[0]);
   cr_expect_eq(sepadata->nlinexprs[1], 4, "\nExpected 4 linear expressions for product 1, got %d", sepadata->nlinexprs[1]);

   /* check the product expressions */
   cr_assert(sepadata->bilinterms[0] != NULL);
   cr_assert(sepadata->bilinterms[1] != NULL);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[0]) == 2);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[1]) == 2);

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[0]);
   cr_expect_eq(var, b1, "Var 0 of product 0 should be b1, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[1]);
   cr_expect_eq(var, x2, "Var 1 of product 0 should be x2, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[0]);
   cr_expect_eq(var, b1, "Var 0 of product 1 should be b1, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[1]);
   cr_expect_eq(var, x1, "Var 1 of product 1 should be x1, got %s", SCIPvarGetName(var));

   /* check the (sorted) linear expressions and sides */
   cr_assert(sepadata->linexprs[0][0] != NULL);
   cr_assert(sepadata->linexprs[1][0] != NULL);

   /* first product: b1x2 */
   /* first linear expression (from first implied relation and x1 upper bound) */
   cr_expect_eq(SCIPgetConsExprExprNChildren(sepadata->linexprs[0][0]), 2, "Linexpr 0 for product b1x2 should have 2 children, got %d",
   SCIPgetConsExprExprNChildren(sepadata->linexprs[0][0]));

   linchildren = SCIPgetConsExprExprChildren(sepadata->linexprs[0][0]);
   lincoefs = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[0][0]);

   var = SCIPgetConsExprExprVarVar(linchildren[0]);
   cr_expect_eq(lincoefs[0], -8.0, "Coefficient of x = b1 should be -8.0, got %f", lincoefs[0]);
   cr_expect_eq(SCIPgetConsExprExprVarVar(linchildren[0]), b1, "Var 0 of linexpr 1 for product b1x2 should be b1, got %s",
   SCIPvarGetName(SCIPgetConsExprExprVarVar(linchildren[0])));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[0][0])[1]);
   cr_expect_eq(var, x1, "Var 1 of linexpr 1 for product b1x2 should be x1, got %s", SCIPvarGetName(var));

   cr_assert(sepadata->linoverestimate[0][0]);
   cr_assert(!sepadata->linunderestimate[0][0]);

   /* second linear expression (from two implied relations) */
   cr_expect_eq(SCIPgetConsExprExprNChildren(sepadata->linexprs[0][1]), 3, "Linexpr 1 for product b1x2 should have 3 children, got %d",
      SCIPgetConsExprExprNChildren(sepadata->linexprs[0][1]));

   linchildren = SCIPgetConsExprExprChildren(sepadata->linexprs[0][1]);
   lincoefs = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[0][1]);

   var = SCIPgetConsExprExprVarVar(linchildren[0]);
   cr_expect_eq(lincoefs[0], -5.0, "Coefficient of x = b1 should be -5.0, got %f", lincoefs[0]);
   cr_expect_eq(SCIPgetConsExprExprVarVar(linchildren[0]), b1, "Var 0 of linexpr 1 for product b1x2 should be b1, got %s",
      SCIPvarGetName(SCIPgetConsExprExprVarVar(linchildren[0])));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[0][1])[1]);
   cr_expect_eq(lincoefs[1], -4.0, "Coefficient of w = x1 should be -4.0, got %f", lincoefs[1]);
   cr_expect_eq(var, x1, "Var 1 of linexpr 1 for product b1x2 should be x1, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[0][1])[2]);
   cr_expect_eq(lincoefs[2], -2.0, "Coefficient of y = x2 should be -2.0, got %f", lincoefs[2]);
   cr_expect_eq(var, x2, "Var 2 of linexpr 1 for product b1x2 should be x2, got %s", SCIPvarGetName(var));

   cr_expect_eq(SCIPgetConsExprExprSumConstant(sepadata->linexprs[0][1]), 1.0, "Linexpr 1 for product b1x2 should have constant 1.0, got %f",
      SCIPgetConsExprExprSumConstant(sepadata->linexprs[0][1]));

   cr_assert(sepadata->linoverestimate[0][1]);
   cr_assert(!sepadata->linunderestimate[0][1]);

   /* third linear expression (from first implied relation and x1 lower bound) */
   cr_expect_eq(SCIPgetConsExprExprNChildren(sepadata->linexprs[0][2]), 3, "Linexpr 2 for product b1x2 should have 3 children, got %d",
      SCIPgetConsExprExprNChildren(sepadata->linexprs[0][2]));

   linchildren = SCIPgetConsExprExprChildren(sepadata->linexprs[0][2]);
   lincoefs = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[0][2]);

   cr_expect_eq(lincoefs[2], 1.0, "Coefficient of y = x2 should be 1.0, got %f", lincoefs[2]);
   cr_expect_eq(SCIPgetConsExprExprSumConstant(sepadata->linexprs[0][2]), 0.0, "Linexpr 2 for product b1x2 should have constant 0.0, got %f",
      SCIPgetConsExprExprSumConstant(sepadata->linexprs[0][2]));

   cr_assert(!sepadata->linunderestimate[0][2]);
   cr_assert(sepadata->linoverestimate[0][2]);

   /* fourth linear expression (from implied relation and bound) */
   cr_expect_eq(SCIPgetConsExprExprNChildren(sepadata->linexprs[0][3]), 3, "Linexpr 1 for product b1x2 should have 3 children, got %d",
   SCIPgetConsExprExprNChildren(sepadata->linexprs[0][3]));

   linchildren = SCIPgetConsExprExprChildren(sepadata->linexprs[0][3]);
   lincoefs = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[0][3]);

   var = SCIPgetConsExprExprVarVar(linchildren[0]);
   cr_expect_eq(lincoefs[0], -9.5, "Coefficient of x = b1 should be -9.5, got %f", lincoefs[0]);
   cr_expect_eq(SCIPgetConsExprExprVarVar(linchildren[0]), b1, "Var 0 of linexpr 1 for product b1x2 should be b1, got %s",
   SCIPvarGetName(SCIPgetConsExprExprVarVar(linchildren[0])));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[0][3])[1]);
   cr_expect_eq(lincoefs[1], 2.0, "Coefficient of w = x1 should be 2.0, got %f", lincoefs[1]);
   cr_expect_eq(var, x1, "Var 1 of linexpr 1 for product b1x2 should be x1, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[0][3])[2]);
   cr_expect_eq(lincoefs[2], 1.0, "Coefficient of y = x2 should be 1.0, got %f", lincoefs[2]);
   cr_expect_eq(var, x2, "Var 2 of linexpr 1 for product b1x2 should be x2, got %s", SCIPvarGetName(var));

   cr_expect_eq(SCIPgetConsExprExprSumConstant(sepadata->linexprs[0][3]), -0.5, "Linexpr 3 for product b1x2 should have constant -0.5, got %f",
   SCIPgetConsExprExprSumConstant(sepadata->linexprs[0][3]));

   cr_assert(!sepadata->linoverestimate[0][3]);
   cr_assert(sepadata->linunderestimate[0][3]);

   /* second product: b1x1 */
   /* second linear expression */
   cr_expect_eq(SCIPgetConsExprExprNChildren(sepadata->linexprs[1][3]), 3, "Linexpr 3 for product b1x1 should have 3 children, got %d",
      SCIPgetConsExprExprNChildren(sepadata->linexprs[1][3]));

   lincoefs = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[1][3]);

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[1][3])[0]);
   cr_expect_eq(lincoefs[0], 2.75, "Coefficient of x = b1 should be 2.75, got %f", coef);
   cr_expect_eq(var, b1, "Var 0 of linexpr 3 for product b1x1 should be b1, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[1][3])[1]);
   cr_expect_eq(lincoefs[1], 3.0, "Coefficient of y = x1 should be 3.0, got %f", coef);
   cr_expect_eq(var, x1, "Var 1 of linexpr 3 for product b1x1 should be x1, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[1][3])[2]);
   cr_expect_eq(lincoefs[2], 1.5, "Coefficient of w = x2 should be 1.5, got %f", coef);
   cr_expect_eq(var, x2, "Var 2 of linexpr 3 for product b1x1 should be x2, got %s", SCIPvarGetName(var));

   cr_expect_eq(SCIPgetConsExprExprSumConstant(sepadata->linexprs[1][3]), -0.75, "Linexpr 3 for product b1x1 should have constant -0.75, got %f",
      SCIPgetConsExprExprSumConstant(sepadata->linexprs[1][3]));

   SCIP_CALL( freeSepaData(scip, sepadata) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(rlt_selection, implrelbnd, .init = setup, .fini = teardown, .description = "test extracting products from an implied relation and an implied bound")
{
   SCIP_ROW* row1;
   SCIP_ROW* row2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   SCIP_SEPADATA* sepadata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR* vars[3];
   SCIP_Real coefs1[3];
   int c, nbdchgs;
   SCIP_Bool cutoff, infeasible;
   SCIP_VAR* var;
   SCIP_Real coef;

   vars[0] = x1; /* w or y */
   vars[1] = x2; /* y or w */
   vars[2] = b1; /* x */

   coefs1[0] = 2.0; /* a1 or c1 */
   coefs1[1] = 2.0; /* c1 or a1 */
   coefs1[2] = 3.0; /* b1 */

   /* d1 = 1, d2 = 3 */
   /* a1 = 2, c1 = 2, b1 = 3, */
   /* a2 = 1, c2 = 0, b2 = -M or */

   /* a1 = 2, c1 = 2, b1 = 3, */
   /* a2 = 0, c2 = 1, b2 = -M */

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   sepadata->isinitialround = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear relations */

   /* 0 <= 2x1 + 2x2 + 3b1 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons1, "c1", 3, vars, coefs1, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons1) );

   /* b1 == 0  =>  x1 <= 3 */
   /* x1 - Mb1 <= 3 */
   SCIP_CALL( SCIPaddVarImplication(scip, b1, FALSE, x1, SCIP_BOUNDTYPE_UPPER, 3.0, &infeasible, &nbdchgs) );

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* create separator data - this will also detect the products */
   SCIP_CALL( createSepaData(scip, sepadata) ); /*a1 = 2; den = a1c2 - c1a2 = 2*0 - 2.0 < 0*/

   /* the product should be: xy >= (1/(a1c2 - c1a2))*(a1a2w + (a2(b1 - d1) + a1d2)x + a1c2y - a1d2)
    * or, substituting the given values: xy <= (-1/2)*(2w + (3 - 1 + 2*3)x - 2*3)
    * xy <= -w - 4x + 3
    * the second product: xy >= 3x + 1y - 3 (x = b1, y = x1)            y(x-1) >= 3(x-1)    y <= 3
    * /

   /* check the numbers */
   cr_expect_eq(sepadata->nbilinvars, 3, "\nExpected 3 bilinear vars, got %d", sepadata->nbilinvars);
   cr_expect_eq(sepadata->nbilinterms, 2, "\nExpected 2 bilinear terms, got %d", sepadata->nbilinterms);
   cr_expect_eq(sepadata->nlinexprs[0], 3, "\nExpected 3 linear expression for product 0, got %d", sepadata->nlinexprs[0]);
   cr_expect_eq(sepadata->nlinexprs[1], 3, "\nExpected 3 linear expression for product 1, got %d", sepadata->nlinexprs[1]);

   /* check the product expressions */
   cr_assert(sepadata->bilinterms[0] != NULL);
   cr_assert(sepadata->bilinterms[1] != NULL);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[0]) == 2);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[1]) == 2);

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[0]);
   cr_expect_eq(var, b1, "Var 0 of product 0 should be b1, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[1]);
   cr_expect_eq(var, x2, "Var 1 of product 0 should be x2, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[0]);
   cr_expect_eq(var, b1, "Var 0 of product 1 should be b1, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[1]);
   cr_expect_eq(var, x1, "Var 1 of product 1 should be x1, got %s", SCIPvarGetName(var));

   /* check the linear expressions obtained from the implied relation and the implied bound */
   cr_assert(sepadata->linexprs[0][1] != NULL);
   cr_assert(sepadata->linexprs[1][0] != NULL);

   /* first linear expression */
   cr_expect_eq(SCIPgetConsExprExprNChildren(sepadata->linexprs[0][1]), 2, "Linexpr for product b1x2 should have 2 children, got %d",
      SCIPgetConsExprExprNChildren(sepadata->linexprs[0][1]));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[0][1])[0]);
   coef = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[0][1])[0];
   cr_expect_eq(var, b1, "Var 0 of linexpr for product b1x2 should be b1, got %s", SCIPvarGetName(var));
   cr_expect_eq(coef, -4, "Coefficient of x = b1 should be -4, got %f", coef);

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[0][1])[1]);
   coef = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[0][1])[1];
   cr_expect_eq(var, x1, "Var 1 of linexpr for product b1x2 should be x1, got %s", SCIPvarGetName(var));
   cr_expect_eq(coef, -1.0, "Coefficient of w = x1 should be -1.0, got %f", coef);

   cr_expect_eq(SCIPgetConsExprExprSumConstant(sepadata->linexprs[0][1]), 3.0, "Linexpr for product b1x2 should have constant 3.0, got %f",
      SCIPgetConsExprExprSumConstant(sepadata->linexprs[0][1]));

   /* second linear expression */
   cr_expect_eq(SCIPgetConsExprExprNChildren(sepadata->linexprs[1][0]), 2, "Linexpr for product b1x1 should have 2 children, got %d",
      SCIPgetConsExprExprNChildren(sepadata->linexprs[1][0]));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[1][0])[0]);
   coef = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[1][0])[0];
   cr_expect_eq(var, b1, "Var 0 of linexpr for product b1x1 should be b1, got %s", SCIPvarGetName(var));
   cr_expect_eq(coef, 3.0, "Coefficient of x = b1 should be 3.0, got %f", coef);

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->linexprs[1][0])[1]);
   coef = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[1][0])[1];
   cr_expect_eq(var, x1, "Var 1 of linexpr for product b1x1 should be x1, got %s", SCIPvarGetName(var));
   cr_expect_eq(coef, 1.0, "Coefficient of y = x1 should be 1.0, got %f", coef);

   cr_expect_eq(SCIPgetConsExprExprSumConstant(sepadata->linexprs[1][0]), -3.0, "Linexpr for product b1x1 should have constant -3.0, got %f",
      SCIPgetConsExprExprSumConstant(sepadata->linexprs[1][0]));

   cr_assert(sepadata->linoverestimate[0][0]);
   cr_assert(sepadata->linunderestimate[1][0]);

   SCIP_CALL( freeSepaData(scip, sepadata) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(rlt_selection, implrelclique, .init = setup, .fini = teardown, .description = "test extracting products from an implied relation and a clique")
{
   SCIP_ROW* row1;
   SCIP_ROW* row2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   SCIP_SEPADATA* sepadata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR* vars[3];
   SCIP_Real coefs1[3];
   int c, nbdchgs;
   SCIP_Bool cutoff, infeasible;
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_VAR* clique_vars[2];
   SCIP_Bool clique_vals[2];

   vars[0] = x1; /* w or y */
   vars[1] = b1; /* y or w */
   vars[2] = b2; /* x */

   coefs1[0] = 2.0; /* a1 or c1 */
   coefs1[1] = 2.0; /* c1 or a1 */
   coefs1[2] = 3.0; /* b1 */

   /* d1 = 1, d2 = 3 */
   /* a1 = 2, c1 = 2, b1 = 3, */
   /* a2 = 1, c2 = 0, b2 = -M or */

   /* a1 = 2, c1 = 2, b1 = 3, */
   /* a2 = 0, c2 = 1, b2 = -M */

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   sepadata->isinitialround = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear relations */

   /* 0 <= 2x1 + 2b1 + 3b2 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons1, "c1", 3, vars, coefs1, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons1) );

   /* add a clique b1 + b2 <= 1 */
   clique_vars[0] = b1;
   clique_vars[1] = b2;
   clique_vals[0] = 1;
   clique_vals[1] = 0;
   SCIP_CALL( SCIPaddClique(scip, clique_vars, clique_vals, 2, FALSE, &infeasible, &nbdchgs) );

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* create separator data - this will also detect the products */
   SCIP_CALL( createSepaData(scip, sepadata) ); /*a1 = 2; den = a1c2 - c1a2 = 2*0 - 2.0 < 0*/

   /* the product should be: xy >= (1/(a1c2 - c1a2))*(a1a2w + (a2(b1 - d1) + a1d2)x + a1c2y - a1d2)
    * or, substituting the given values: xy <= (-1/2)*(2w + (3 - 1 + 2*3)x - 2*3)
    * xy <= -w - 4x + 3
    * the second product: xy >= 3x + 1y - 3 (x = b1, y = x1)
    * /

   /* check the numbers */
//   cr_expect_eq(sepadata->nbilinvars, 3, "\nExpected 3 bilinear vars, got %d", sepadata->nbilinvars);
//   cr_expect_eq(sepadata->nbilinterms, 2, "\nExpected 2 bilinear terms, got %d", sepadata->nbilinterms);
//   cr_expect_eq(sepadata->nlinexprs[0], 3, "\nExpected 3 linear expression for product 0, got %d", sepadata->nlinexprs[0]);
//   cr_expect_eq(sepadata->nlinexprs[1], 3, "\nExpected 3 linear expression for product 1, got %d", sepadata->nlinexprs[1]);

//   /* check the product expressions */
//   cr_assert(sepadata->bilinterms[0] != NULL);
//   cr_assert(sepadata->bilinterms[1] != NULL);
//   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[0]) == 2);
//   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[1]) == 2);

//   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[0]);
//   cr_expect_eq(var, b1, "Var 0 of product 0 should be b1, got %s", SCIPvarGetName(var));
//   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[1]);
//   cr_expect_eq(var, x2, "Var 1 of product 0 should be x2, got %s", SCIPvarGetName(var));

//   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[0]);
//   cr_expect_eq(var, b1, "Var 0 of product 1 should be b1, got %s", SCIPvarGetName(var));
//   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[1]);
//   cr_expect_eq(var, x1, "Var 1 of product 1 should be x1, got %s", SCIPvarGetName(var));


   SCIP_CALL( freeSepaData(scip, sepadata) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIPfreeBuffer(scip, &sepadata);
}