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

static
void checkLinearisation(
   SCIP_SEPADATA*        sepadata,
   int                   prodidx,
   int                   linidx,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   int                   nvars,
   SCIP_Real             cst,
   SCIP_Bool             underestimate,
   SCIP_Bool             overestimate
   )
{
   SCIP_CONSEXPR_EXPR* linexpr;
   SCIP_CONSEXPR_EXPR* linchildren[3];
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_Bool changed;
   SCIP_Bool infeasible;
   int i;
   int datanvars;
   SCIP_CONSEXPR_EXPR** datachildren;
   SCIP_Real* datacoefs;
   SCIP_Real datacst;
   int compare;

   datacoefs = SCIPgetConsExprExprSumCoefs(sepadata->linexprs[prodidx][linidx]);
   datanvars = SCIPgetConsExprExprNChildren(sepadata->linexprs[prodidx][linidx]);
   datacst = SCIPgetConsExprExprSumConstant(sepadata->linexprs[prodidx][linidx]);

   cr_expect(datanvars == nvars, "linearisation expression [%d][%d] should have %d vars, got %d", prodidx, linidx, nvars, datanvars);

   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &linchildren[i], vars[i]) );
      if( SCIPisEQ(scip, vals[i], datacoefs[i]) )
      {
         vals[i] = datacoefs[i];
      }
   }

   if( SCIPisEQ(scip, cst, datacst) )
   {
      cst = datacst;
   }

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &linexpr, nvars, linchildren, vals, cst) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, linexpr, &simplified, &changed, &infeasible) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &linexpr) );
   linexpr = simplified;

   compare = SCIPcompareConsExprExprs(linexpr, sepadata->linexprs[prodidx][linidx]);

   if( compare != 0 )
   {
      SCIPinfoMessage(scip, NULL, "\n\nprinting the two non-matching exprs (compare = %d):\n", compare);
      SCIPprintConsExprExpr(scip, conshdlr, linexpr, NULL);
      SCIPinfoMessage(scip, NULL, "\n");
      SCIPprintConsExprExpr(scip, conshdlr, sepadata->linexprs[prodidx][linidx], NULL);
      SCIPinfoMessage(scip, NULL, "\n");
      SCIPinfoMessage(scip, NULL, "\n");
   }

   cr_expect(compare == 0, "linearisation expression [%d][%d] doesn't match expected", prodidx, linidx);

   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &linchildren[i]) );
   }

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &linexpr) );

   cr_expect(sepadata->linunderestimate[prodidx][linidx] == underestimate,
         "linearisation expression [%d][%d] should %sunderestimate", prodidx, linidx, underestimate ? "" : "NOT ");
   cr_expect(sepadata->linoverestimate[prodidx][linidx] == overestimate,
         "linearisation expression [%d][%d] should %soverestimate", prodidx, linidx, overestimate ? "" : "NOT ");
}

Test(rlt_product_detection, implrels, .init = setup, .fini = teardown, .description = "test extracting products from two implied relations")
{
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   SCIP_SEPADATA* sepadata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR* vars[3];
   SCIP_Real coefs1[3];
   SCIP_Real coefs2[3];
   SCIP_Bool cutoff;
   SCIP_VAR* var;
   SCIP_Real* coefs_expected;
   SCIP_VAR** vars_expected;

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
   sepadata->detecthidden = TRUE;
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
   /* first linear expression (from first implied relation (rhs) and x1 upper bound, binary variable = b1) */
   vars_expected = (SCIP_VAR*[2]) {b1, x1};
   coefs_expected = (SCIP_Real[2]) {-8.0, -4.0/3.0};
   checkLinearisation(sepadata, 0, 0, vars_expected, coefs_expected, 2, 20.0/3.0, FALSE, TRUE);

   /* second linear expression (from two implied relations, binary variable = x1) */
   vars_expected = (SCIP_VAR*[3]) {b1, x1, x2};
   coefs_expected = (SCIP_Real[3]) {-5.0, -4.0, -2.0};
   checkLinearisation(sepadata, 0, 1, vars_expected, coefs_expected, 3, 1.0, FALSE, TRUE);

   /* third linear expression (from first implied relation (lhs) and x1 lower bound, binary variable = x1) */
   vars_expected = (SCIP_VAR*[3]) {b1, x1, x2};
   coefs_expected = (SCIP_Real[3]) {4.0/3.0, 4.0/3.0, 1.0};
   checkLinearisation(sepadata, 0, 2, vars_expected, coefs_expected, 3, 0.0, FALSE, TRUE);

   /* fourth linear expression (from implied relation and bound) */
   vars_expected = (SCIP_VAR*[3]) {b1, x1, x2};
   coefs_expected = (SCIP_Real[3]) {-9.5, 2.0, 1.0};
   checkLinearisation(sepadata, 0, 3, vars_expected, coefs_expected, 3, -0.5, TRUE, FALSE);

   /* second product: b1x1 */
   /* second linear expression */
   vars_expected = (SCIP_VAR*[3]) {b1, x1, x2};
   coefs_expected = (SCIP_Real[3]) {2.75, 3.0, 1.5};
   checkLinearisation(sepadata, 1, 3, vars_expected, coefs_expected, 3, -0.75, TRUE, FALSE);

   SCIP_CALL( freeSepaData(scip, sepadata) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(rlt_product_detection, implrelbnd, .init = setup, .fini = teardown, .description = "test extracting products from an implied relation and an implied bound")
{
   SCIP_CONS* cons;
   SCIP_SEPADATA* sepadata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR* vars[3];
   SCIP_Real coefs[3];
   int nbdchgs;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_VAR* var;
   SCIP_Real* coefs_expected;
   SCIP_VAR** vars_expected;

   vars[0] = x1; /* w or y */
   vars[1] = x2; /* y or w */
   vars[2] = b1; /* x */

   coefs[0] = 2.0; /* a1 or c1 */
   coefs[1] = 2.0; /* c1 or a1 */
   coefs[2] = 3.0; /* b1 */

   /* d1 = 1, d2 = 3 */
   /* a1 = 2, c1 = 2, b1 = 3, */
   /* a2 = 1, c2 = 0, b2 = -M or */

   /* a1 = 2, c1 = 2, b1 = 3, */
   /* a2 = 0, c2 = 1, b2 = -M */

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   sepadata->isinitialround = TRUE;
   sepadata->detecthidden = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear relations */

   /* 0 <= 2x1 + 2x2 + 3b1 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 3, vars, coefs, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

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

   cr_expect_eq(sepadata->nlinexprs[0], 3, "\nExpected 3 linear expressions for product b1x2, got %d", sepadata->nlinexprs[0]);
   cr_expect_eq(sepadata->nlinexprs[1], 2, "\nExpected 2 linear expressions for product b1x1, got %d", sepadata->nlinexprs[1]);

   /* check the linear expressions obtained from the implied relation and the implied bound */
   cr_assert(sepadata->linexprs[0][1] != NULL);
   cr_assert(sepadata->linexprs[1][0] != NULL);

   /* first linear expression */
   vars_expected = (SCIP_VAR*[2]) {b1, x1};
   coefs_expected = (SCIP_Real[2]) {-4.0, -1.0};
   checkLinearisation(sepadata, 0, 1, vars_expected, coefs_expected, 2, 3.0, FALSE, TRUE);

   SCIP_CALL( freeSepaData(scip, sepadata) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(rlt_product_detection, implrelclique, .init = setup, .fini = teardown, .description = "test extracting products from an implied relation and a clique")
{
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   SCIP_SEPADATA* sepadata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR* vars[3];
   SCIP_Real coefs1[3];
   int nbdchgs;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_VAR* var;
   SCIP_VAR* clique_vars[2];
   SCIP_Bool clique_vals[2];
   SCIP_Real* coefs_expected;
   SCIP_VAR** vars_expected;

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
   sepadata->detecthidden = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear relations */

   /* 0 <= 2x1 + 2b1 + 3b2 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons1, "c1", 3, vars, coefs1, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons1) );
   scip->stat->nnz = 3; /* set this manually so that the clique is not rejected */

   /* add a clique b1 + (1-b2) <= 1 */
   clique_vars[0] = b1;
   clique_vars[1] = b2;
   clique_vals[0] = 1;
   clique_vals[1] = 0;
   SCIP_CALL( SCIPaddClique(scip, clique_vars, clique_vals, 2, FALSE, &infeasible, &nbdchgs) );

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) ); /* TODO remove this? */

   /* create separator data - this will also detect the products */
   SCIP_CALL( createSepaData(scip, sepadata) ); /*a1 = 2; den = a1c2 - c1a2 = 2*0 - 2.0 < 0*/

   /* the product should be: xy >= (1/(a1c2 - c1a2))*(a1a2w + (a2(b1 - d1) + a1d2)x + a1c2y - a1d2)
    * or, substituting the given values: xy <= (-1/2)*(2w + (3 - 1 + 2*3)x - 2*3)
    * xy <= -w - 4x + 3
    * the second product: xy >= 3x + 1y - 3 (x = b1, y = x1)
    * /

   /* check the numbers */
   cr_expect_eq(sepadata->nbilinvars, 3, "\nExpected 3 bilinear vars, got %d", sepadata->nbilinvars);
   cr_expect_eq(sepadata->nbilinterms, 3, "\nExpected 3 bilinear terms, got %d", sepadata->nbilinterms);
   cr_expect_eq(sepadata->nlinexprs[0], 3, "\nExpected 3 linear expressions for product 0, got %d", sepadata->nlinexprs[0]);
   cr_expect_eq(sepadata->nlinexprs[1], 4, "\nExpected 4 linear expressions for product 1, got %d", sepadata->nlinexprs[1]);
   cr_expect_eq(sepadata->nlinexprs[2], 3, "\nExpected 3 linear expressions for product 2, got %d", sepadata->nlinexprs[2]);

   /* check the product expressions */
   cr_assert(sepadata->bilinterms[0] != NULL);
   cr_assert(sepadata->bilinterms[1] != NULL);
   cr_assert(sepadata->bilinterms[2] != NULL);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[0]) == 2);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[1]) == 2);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[2]) == 2);

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[0]);
   cr_expect_eq(var, b1, "Var 0 of product 0 should be b1, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[1]);
   cr_expect_eq(var, x1, "Var 1 of product 0 should be x1, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[0]);
   cr_expect_eq(var, b1, "Var 0 of product 1 should be b1, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[1]);
   cr_expect_eq(var, b2, "Var 1 of product 1 should be b2, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[2])[0]);
   cr_expect_eq(var, b2, "Var 0 of product 2 should be b2, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[2])[1]);
   cr_expect_eq(var, x1, "Var 1 of product 2 should be x1, got %s", SCIPvarGetName(var));

   /* t_b1t_x1 <= 1.5t_b2 - 1.5t_b1 + t_x1 (from constraint and clique, (x,w,y) = (b1,b2,x1)) */
   vars_expected = (SCIP_VAR*[3]) {b1, b2, x1};
   coefs_expected = (SCIP_Real[3]) {-1.5, 1.5, 1.0};
   checkLinearisation(sepadata, 0, 2, vars_expected, coefs_expected, 3, 0.0, FALSE, TRUE);

   /* t_b2t_x1 <= -t_b1 - t_b2 (from constraint and clique, (x,w,y) = (b2,b1,x1)) */
   vars_expected = (SCIP_VAR*[2]) {b1, b2};
   coefs_expected = (SCIP_Real[3]) {-1.0, -1.0};
   checkLinearisation(sepadata, 2, 1, vars_expected, coefs_expected, 2, 0.0, FALSE, TRUE);

   SCIP_CALL( freeSepaData(scip, sepadata) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(rlt_product_detection, implbnd, .init = setup, .fini = teardown, .description = "test extracting products from an implied bound and an unconditional relation")
{
   SCIP_CONS* cons;
   SCIP_SEPADATA* sepadata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR* vars[2];
   SCIP_Real coefs[2];
   int nbdchgs;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_VAR* var;
   SCIP_Real* coefs_expected;
   SCIP_VAR** vars_expected;

   vars[0] = x1;
   vars[1] = x2;

   coefs[0] = 2.0;
   coefs[1] = 2.0;

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   sepadata->isinitialround = TRUE;
   sepadata->detecthidden = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear relations */

   /* 0 <= 2x1 + 2x2 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 2, vars, coefs, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* b1 == 0  =>  x1 <= 1 */
   /* x1 <= -2b1 + 1 */
   SCIPaddVarVub(scip, x1, b1, -2.0, 1.0, &infeasible, &nbdchgs);

   /* x1 <= b2 + 0.0 */
   SCIPaddVarVub(scip, x1, b2, 1.0, 0.0, &infeasible, &nbdchgs);

   SCIPinfoMessage(scip, NULL, "\nvar x1 has %d vubs:", SCIPvarGetNVubs(x1));
   for( int i = 0; i < SCIPvarGetNVubs(x1); ++i )
   {
      SCIPinfoMessage(scip, NULL, "\nx1 <= %f%s + %f", SCIPvarGetVubCoefs(x1)[i], SCIPvarGetName(SCIPvarGetVubVars(x1)[i]), SCIPvarGetVubConstants(x1)[i]);
   }

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* create separator data - this will also detect the products */
   SCIP_CALL( createSepaData(scip, sepadata) );

   /* check the numbers */
   cr_expect_eq(sepadata->nbilinvars, 3, "\nExpected 3 bilinear vars, got %d", sepadata->nbilinvars);
   cr_expect_eq(sepadata->nbilinterms, 3, "\nExpected 3 bilinear terms, got %d", sepadata->nbilinterms);
   cr_expect_eq(sepadata->nlinexprs[0], 2, "\nExpected 3 linear expressions for product 0, got %d", sepadata->nlinexprs[0]);
   cr_expect_eq(sepadata->nlinexprs[1], 1, "\nExpected 3 linear expressions for product 1, got %d", sepadata->nlinexprs[1]);
   cr_expect_eq(sepadata->nlinexprs[2], 1, "\nExpected 3 linear expressions for product 2, got %d", sepadata->nlinexprs[2]);

   /* check the product expressions */
   cr_assert(sepadata->bilinterms[0] != NULL);
   cr_assert(sepadata->bilinterms[1] != NULL);
   cr_assert(sepadata->bilinterms[2] != NULL);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[0]) == 2);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[1]) == 2);
   cr_assert(SCIPgetConsExprExprNChildren(sepadata->bilinterms[2]) == 2);

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[0]);
   cr_expect_eq(var, b1, "Var 0 of product 0 should be b1, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[0])[1]);
   cr_expect_eq(var, b2, "Var 1 of product 0 should be b2, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[0]);
   cr_expect_eq(var, b1, "Var 0 of product 1 should be b1, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[1])[1]);
   cr_expect_eq(var, x2, "Var 1 of product 1 should be x2, got %s", SCIPvarGetName(var));

   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[2])[0]);
   cr_expect_eq(var, b2, "Var 0 of product 1 should be b2, got %s", SCIPvarGetName(var));
   var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(sepadata->bilinterms[2])[1]);
   cr_expect_eq(var, x2, "Var 1 of product 1 should be x2, got %s", SCIPvarGetName(var));

   /* check the linear expressions obtained from the implied relation and the implied bound */
   cr_assert(sepadata->linexprs[0][1] != NULL);
   cr_assert(sepadata->linexprs[1][0] != NULL);
   cr_assert(sepadata->linexprs[2][0] != NULL);

   /* first linear expression (from two implied bounds on x1) */
   /* should be t_b1t_b2 <= -1t_x1 + -1t_b1 + 1t_b2 + 0 */
   vars_expected = (SCIP_VAR*[3]) {b1, b2, x1};
   coefs_expected = (SCIP_Real[3]) {-1.0, 1.0, -1.0};
   checkLinearisation(sepadata, 0, 0, vars_expected, coefs_expected, 3, 0.0, FALSE, TRUE);

   /* second linear expression (from two implied bounds on x1) */
   /* should be t_b2t_b1 <= -0.5t_x1 + 0.5t_b2 */
   vars_expected = (SCIP_VAR*[2]) {b2, x1};
   coefs_expected = (SCIP_Real[2]) {0.5, -0.5};
   checkLinearisation(sepadata, 0, 1, vars_expected, coefs_expected, 2, 0.0, FALSE, TRUE);

   /* third linear expression (from implied bound on x1 (with b1) and unconditional) */
   /* should be t_b1t_x2 >= 1t_x1 + 1.5t_b1 + 1t_x2 + -0.5 */
   vars_expected = (SCIP_VAR*[3]) {b1, x1, x2};
   coefs_expected = (SCIP_Real[3]) {1.5, 1.0, 1.0};
   checkLinearisation(sepadata, 1, 0, vars_expected, coefs_expected, 3, -0.5, TRUE, FALSE);

   /* fourth linear expression (from implied bound on x1 (with b2) and unconditional) */
   /* should be t_b2t_x2 <= -1t_x1 + 0.5t_b2 */
   vars_expected = (SCIP_VAR*[2]) {b2, x1};
   coefs_expected = (SCIP_Real[2]) {0.5, -1.0};
   checkLinearisation(sepadata, 2, 0, vars_expected, coefs_expected, 2, 0.0, FALSE, TRUE);

   SCIP_CALL( freeSepaData(scip, sepadata) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(rlt_product_detection, rowlist, .init = setup, .fini = teardown, .description = "test creating linked lists of rows")
{
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   SCIP_CONS* cons3;
   SCIP_CONS* cons4;
   SCIP_CONS* cons5;
   SCIP_VAR* vars[3];
   SCIP_Real coefs1[3];
   SCIP_Bool cutoff;
   SCIP_VAR* var;
   SCIP_Real coef;
   int nrows;
   int i;
   int nvars_in_2rels;
   int r;
   SCIP_ROW** prob_rows;
   SCIP_ROW* row1;
   int* row_list;
   SCIP_HASHTABLE* hashtable2;
   SCIP_HASHTABLE* hashtable3;
   SCIP_VAR** vars_in_2rels;
   SCIP_VAR*** related_vars;
   int* nrelated_vars;
   HASHDATA* foundhashdata;

   /* create linear relations */

   vars[0] = x1;
   vars[1] = x2;
   vars[2] = b1;

   /* 0 <= 2x1 + 2x2 <= 1 */
   coefs1[0] = 2.0; coefs1[1] = 2.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons1, "c1", 2, vars, coefs1, 0.0, 1.0) );

   /* 0 <= 2x1 + 1.5x2 + 3b1 <= 1 */
   coefs1[0] = 2.0; coefs1[1] = 1.5; coefs1[2] = 3.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons2, "c2", 3, vars, coefs1, 0.0, 1.0) );

   /* 4x1 + 2x2 - b1 <= 1 */
   coefs1[0] = 4.0; coefs1[1] = 2.0; coefs1[2] = -1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons3, "c3", 3, vars, coefs1, -SCIPinfinity(scip), 1.0) );

   /* x1 + 3x2 + b1 <= 2 */
   coefs1[0] = 1.0; coefs1[1] = 3.0; coefs1[2] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons4, "c4", 3, vars, coefs1, -SCIPinfinity(scip), 2.0) );

   /* x1 - 2x2 <= 1 */
   coefs1[0] = 1.0; coefs1[1] = -2.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons5, "c5", 2, vars, coefs1, -SCIPinfinity(scip), 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons1) );
   SCIP_CALL( SCIPaddCons(scip, cons2) );
   SCIP_CALL( SCIPaddCons(scip, cons3) );
   SCIP_CALL( SCIPaddCons(scip, cons4) );
   SCIP_CALL( SCIPaddCons(scip, cons5) );

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* create the rows */
   nrows = 5;
   SCIP_CALL( getInitialRows(scip, &prob_rows, &nrows) );
   SCIPsortPtr((void**)prob_rows, SCIProwComp, nrows);

   /* create tables of implied and unconditional relations */
   SCIP_CALL( SCIPhashtableCreate(&hashtable3, SCIPblkmem(scip), nrows, SCIPhashGetKeyStandard,
      hashdataKeyEqConss, hashdataKeyValConss, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&hashtable2, SCIPblkmem(scip), nrows, SCIPhashGetKeyStandard,
      hashdataKeyEqConss, hashdataKeyValConss, (void*) scip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &row_list, nrows) );

   /* allocate the array of variables that appear in 2-var relations */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars_in_2rels, SCIPgetNVars(scip)) );
   /* allocate the array of arrays of variables that appear in 2-var relations with each variable */
   SCIP_CALL( SCIPallocBufferArray(scip, &related_vars, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nrelated_vars, SCIPgetNVars(scip)) );

   /* create separator data - this will also detect the products */
   SCIP_CALL( createRelationTables(scip, prob_rows, nrows, hashtable2, hashtable3, vars_in_2rels, related_vars,
      nrelated_vars, &nvars_in_2rels, row_list) );

   for( i = 0; i < SCIPhashtableGetNEntries(hashtable3); ++i )
   {
      foundhashdata = (HASHDATA*)SCIPhashtableGetEntry(hashtable3, i);
      if( foundhashdata == NULL )
         continue;

      SCIPdebugMsg(scip, "(%s, %s, %s): ", SCIPvarGetName(foundhashdata->vars[0]),
         SCIPvarGetName(foundhashdata->vars[1]), SCIPvarGetName(foundhashdata->vars[2]));

      SCIPinfoMessage(scip, NULL, "\nrow array:");
      r = foundhashdata->firstrow;
      while( r != -1 )
      {
         assert(r < nrows && r >= 0);
         row1 = prob_rows[r];
         SCIPinfoMessage(scip, NULL, "%d; ", SCIProwGetIndex(row1));
         r = row_list[r];
      }

      SCIPfreeBufferArray(scip, &(foundhashdata->vars));
      SCIPfreeBuffer(scip, &foundhashdata);
   }

   for( i = 0; i < SCIPhashtableGetNEntries(hashtable2); ++i )
   {
      foundhashdata = (HASHDATA*)SCIPhashtableGetEntry(hashtable2, i);
      if( foundhashdata == NULL )
         continue;

      SCIPdebugMsg(scip, "(%s, %s): ", SCIPvarGetName(foundhashdata->vars[0]),
      SCIPvarGetName(foundhashdata->vars[1]));

      SCIPinfoMessage(scip, NULL, "\nrow array:");
      r = foundhashdata->firstrow;
      while( r != -1 )
      {
         assert(r < nrows && r >= 0);
         row1 = prob_rows[r];
         SCIPinfoMessage(scip, NULL, "%d; ", SCIProwGetIndex(row1));
         r = row_list[r];
      }

      SCIPfreeBufferArray(scip, &(foundhashdata->vars));
      SCIPfreeBuffer(scip, &foundhashdata);
   }

   for( i = 0; i < nvars_in_2rels; ++i )
   {
      SCIPfreeBufferArray(scip, &related_vars[i]);
   }

   SCIPfreeBufferArray(scip, &nrelated_vars);
   SCIPfreeBufferArray(scip, &related_vars);
   SCIPfreeBufferArray(scip, &vars_in_2rels);

   SCIPfreeBufferArray(scip, &row_list);

   SCIPhashtableFree(&hashtable2);
   SCIPhashtableFree(&hashtable3);

   SCIPfreeBufferArray(scip, &prob_rows);

   SCIP_CALL( SCIPreleaseCons(scip, &cons5) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons4) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons3) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
}
