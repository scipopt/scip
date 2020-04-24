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

/**@file   readers.c
 * @brief  tests readers
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "scip/cons_expr.c"
#include <libgen.h>
#include "include/scip_test.h"

Test(readers, pip)
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR** children;
   SCIP_CONSEXPR_EXPR** grandchildren;
   char filename[SCIP_MAXSTRLEN];
   SCIP_CONSEXPR_EXPRHDLR* sumhdlr;

   /* get file to read: test.mps that lives in the same directory as this file */
   (void)SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s", __FILE__);
   dirname(filename);
   strcat(filename, "/test.pip");
   printf("Reading %s\n", filename);

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPreadProb(scip, filename, NULL));

   /* check that vars are what we expect */
   cr_expect_eq(SCIPgetNVars(scip), 5);
   vars = SCIPgetVars(scip);

   cr_expect_str_eq(SCIPvarGetName(vars[0]), "x");
   cr_expect_str_eq(SCIPvarGetName(vars[1]), "y");
   cr_expect_str_eq(SCIPvarGetName(vars[2]), "z");
   cr_expect_str_eq(SCIPvarGetName(vars[3]), "objconst");
   cr_expect_str_eq(SCIPvarGetName(vars[4]), "nonlinobjvar");

   /* check that cons are what we expect */
   cr_expect_eq(SCIPgetNConss(scip), 3);
   conss = SCIPgetConss(scip);

   cr_expect_str_eq(SCIPconsGetName(conss[0]), "nonlinobj");
   cr_expect_str_eq(SCIPconsGetName(conss[1]), "e1");
   cr_expect_str_eq(SCIPconsGetName(conss[2]), "e2");

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   sumhdlr = SCIPgetConsExprExprHdlrSum(conshdlr);
   cr_assert_not_null(sumhdlr);

   /* check objective coefficients */
   cr_expect_eq(SCIPvarGetObj(vars[0]), 0.0);
   cr_expect_eq(SCIPvarGetObj(vars[1]), 0.0);
   cr_expect_eq(SCIPvarGetObj(vars[2]), 0.0);
   cr_expect_eq(SCIPvarGetObj(vars[3]), 11.0);
   cr_expect_eq(SCIPvarGetObj(vars[4]), 1.0);

   /*
    * check whether the first constraint is nonlinear and of the form + x - y^3 z^2 - nonlinobj <= 0
    */

   cr_expect_eq(SCIPconsGetHdlr(conss[0]), SCIPfindConshdlr(scip, "expr"));
   expr = SCIPgetExprConsExpr(scip, conss[0]);
   cr_expect_not_null(expr);
   children = SCIPgetConsExprExprChildren(expr);
   cr_expect_not_null(children);
   cr_expect_eq(SCIPgetConsExprExprNChildren(expr), 3);

   /* check sides */
   cr_expect(SCIPisInfinity(scip, -SCIPgetLhsConsExpr(scip, conss[0])));
   cr_expect_eq(SCIPgetRhsConsExpr(scip, conss[0]), 0.0);

   /* check constant and coefficients of the sum expression */
   cr_expect_eq(SCIPgetConsExprExprSumConstant(expr), 0.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[0], 1.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[1], -1.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[2], -1.0);

   /* check whether first and third child is a variable */
   cr_expect(SCIPisConsExprExprVar(children[0]));
   cr_expect(SCIPisConsExprExprVar(children[2]));

   /* check whether second child is a product; both grandchildren need to be power expressions */
   cr_expect(SCIPgetConsExprExprHdlr(children[1]) == SCIPfindConsExprExprHdlr(conshdlr, "prod"));
   grandchildren = SCIPgetConsExprExprChildren(children[1]);
   cr_expect_eq(SCIPgetConsExprExprNChildren(children[1]), 2);
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[0]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[1]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[0]), 3.0);
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[1]), 2.0);

   /*
    * check whether the second constraint is linear and of the form x + 2y + 3z <= 1
    */

   cr_expect_eq(SCIPconsGetHdlr(conss[1]), SCIPfindConshdlr(scip, "linear"));
   cr_expect(SCIPisInfinity(scip, -SCIPgetLhsLinear(scip, conss[1])));
   cr_expect_eq(SCIPgetRhsLinear(scip, conss[1]), 1.0);
   cr_expect_eq(SCIPgetNVarsLinear(scip, conss[1]), 3);
   cr_expect_eq(SCIPgetValsLinear(scip, conss[1])[0], 1.0);
   cr_expect_eq(SCIPgetValsLinear(scip, conss[1])[1], 2.0);
   cr_expect_eq(SCIPgetValsLinear(scip, conss[1])[2], 3.0);

   /*
    * check whether the third constraint is nonlinear and of the form x^2 y^3 z^4 + x + y + 2 = 10
    */

   cr_expect_eq(SCIPconsGetHdlr(conss[2]), SCIPfindConshdlr(scip, "expr"));
   expr = SCIPgetExprConsExpr(scip, conss[2]);
   cr_expect_not_null(expr);
   children = SCIPgetConsExprExprChildren(expr);
   cr_expect_not_null(children);
   cr_expect_eq(SCIPgetConsExprExprNChildren(expr), 3);

   /* check sides */
   cr_expect_eq(SCIPgetLhsConsExpr(scip, conss[2]), 10.0);
   cr_expect_eq(SCIPgetRhsConsExpr(scip, conss[2]), 10.0);

   /* check constant and coefficients of the sum expression */
   cr_expect_eq(SCIPgetConsExprExprSumConstant(expr), 7.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[0], 4.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[1], 5.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[2], 6.0);

   /* check whether first child is a product; all three grandchildren need to be power expressions */
   cr_expect(SCIPgetConsExprExprHdlr(children[0]) == SCIPfindConsExprExprHdlr(conshdlr, "prod"));
   grandchildren = SCIPgetConsExprExprChildren(children[0]);
   cr_expect_eq(SCIPgetConsExprExprNChildren(children[0]), 3);
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[0]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[1]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect(SCIPgetConsExprExprHdlr(grandchildren[2]) == SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[0]), 2.0);
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[1]), 3.0);
   cr_expect_eq(SCIPgetConsExprExprPowExponent(grandchildren[2]), 4.0);
}

Test(readers, mps1)
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPR_EXPR* expr;
   char filename[SCIP_MAXSTRLEN];
   SCIP_CONSEXPR_EXPRHDLR* sumhdlr;

   /* get file to read: test.mps that lives in the same directory as this file */
   (void)SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s", __FILE__);
   dirname(filename);
   strcat(filename, "/test.mps");
   printf("Reading %s\n", filename);

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPreadProb(scip, filename, NULL));

   /* check that vars are what we expect */
   cr_expect_eq(SCIPgetNVars(scip), 4);
   vars = SCIPgetVars(scip);

   cr_expect_str_eq(SCIPvarGetName(vars[0]), "x1");
   cr_expect_str_eq(SCIPvarGetName(vars[1]), "x2");
   cr_expect_str_eq(SCIPvarGetName(vars[2]), "x3");
   cr_expect_str_eq(SCIPvarGetName(vars[3]), "qmatrixvar");

   /* check that cons are what we expect */
   cr_expect_eq(SCIPgetNConss(scip), 4);
   conss = SCIPgetConss(scip);

   cr_expect_str_eq(SCIPconsGetName(conss[0]), "c0");
   cr_expect_str_eq(SCIPconsGetName(conss[1]), "c1");
   cr_expect_str_eq(SCIPconsGetName(conss[2]), "c2");
   cr_expect_str_eq(SCIPconsGetName(conss[3]), "qmatrix");

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   sumhdlr = SCIPgetConsExprExprHdlrSum(conshdlr);
   cr_assert_not_null(sumhdlr);

   /* check some things from the constraints which should be
    * c0: x1 + 9 x1^2 + 0.5 x1 * x2 + 0.5 x2 * x1 = 0.1
    * c1: x2 + 2 x3^2 <= 0.2
    * c2: 0.5 x1 + x1^2 + x2^2 - x3^2 >= 0.3
    * qmatrix: 2x1*x2 + 0.1x2*x3 -0.1 x3*x1 <= qmatrixvar
    */
   SCIP_Real lhs[4] = {0.1, -SCIPinfinity(scip), 0.3, -SCIPinfinity(scip)};
   SCIP_Real rhs[4] = {0.1, 0.2, SCIPinfinity(scip),  0.0};
   int expectednnonz[4] = {4, 2, 4, 4};
   /* cons expr creates first the quadratic terms and then the linear terms; note that for whatever reason the quadratic
    * part of the objective is divided by two */
   SCIP_Real expectedcoeffs[4][4] = {{9.0, 0.5, 0.5, 1.0}, {2.0, 1.0, SCIPinfinity(scip), SCIPinfinity(scip)},
      {1.0,1.0,-1.0, 0.5}, {1.0,0.05,-0.05,-1.0} };

   for( int i = 0; i < 4; ++i )
   {
      cr_assert_eq(conshdlr, SCIPconsGetHdlr(conss[i]));
      cr_expect_eq(SCIPgetLhsConsExpr(scip, conss[i]), lhs[i], "lhs cons %d: expected %g, got %g\n", i, lhs[i], SCIPgetLhsConsExpr(scip, conss[i]));
      cr_expect_eq(SCIPgetRhsConsExpr(scip, conss[i]), rhs[i], "rhs cons %d: expected %g, got %g\n", i, rhs[i], SCIPgetRhsConsExpr(scip, conss[i]));

      expr = SCIPgetExprConsExpr(scip, conss[i]);

      cr_assert_not_null(expr);
      cr_assert_eq(expr->exprhdlr, sumhdlr);
      cr_expect_eq(SCIPgetConsExprExprNChildren(expr), expectednnonz[i]);
      cr_expect_eq(SCIPgetConsExprExprSumConstant(expr), 0.0);
      for( int j = 0; j < expectednnonz[i]; ++j )
         cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[j], expectedcoeffs[i][j], "i,j = %d,%d: expected %g, got %g\n",
               i, j, expectedcoeffs[i][j], SCIPgetConsExprExprSumCoefs(expr)[j]);
   }

   /* check constraint c1 */
   expr = SCIPgetExprConsExpr(scip, conss[1]);
   cr_expect_eq(expr->children[0]->exprhdlr, SCIPgetConsExprExprHdlrPower(conshdlr));
   cr_expect_eq(expr->children[0]->children[0]->exprhdlr, SCIPgetConsExprExprHdlrVar(conshdlr));
   cr_expect_eq(SCIPgetConsExprExprVarVar(expr->children[0]->children[0]), vars[2]);

   cr_expect_eq(expr->children[1]->exprhdlr, SCIPgetConsExprExprHdlrVar(conshdlr));
   cr_expect_eq(SCIPgetConsExprExprVarVar(expr->children[1]), vars[1]);
}

Test(readers, zimpl)
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPR_EXPR* expr;
   char filename[SCIP_MAXSTRLEN];
   SCIP_CONSEXPR_EXPRHDLR* sumhdlr;

   /* get file to read: test.zpl that lives in the same directory as this file */
   (void)SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s", __FILE__);
   dirname(filename);
   strcat(filename, "/test.zpl");
   printf("Reading %s\n", filename);

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   if( SCIPfindReader(scip, "zplreader") == NULL )
      return;

   SCIP_CALL( SCIPreadProb(scip, filename, NULL));

   /* check that vars are what we expect; zimpl will create 2 auxiliary variables, hence we expect 5 */
   cr_expect_eq(SCIPgetNVars(scip), 5, "\nexpected 5 variables, got %d", SCIPgetNVars(scip));
   vars = SCIPgetVars(scip);

   cr_expect_str_eq(SCIPvarGetName(vars[0]), "X1");
   cr_expect_str_eq(SCIPvarGetName(vars[1]), "X2");
   cr_expect_str_eq(SCIPvarGetName(vars[2]), "X3");
   cr_expect_str_eq(SCIPvarGetName(vars[3]), "@@polyfun_1_t_0");
   cr_expect_str_eq(SCIPvarGetName(vars[4]), "@@polyfun_1_r_1");

   /* check that cons are what we expect:  */
   cr_expect_eq(SCIPgetNConss(scip), 5);
   conss = SCIPgetConss(scip);

   cr_expect_str_eq(SCIPconsGetName(conss[0]), "quad_1");
   cr_expect_str_eq(SCIPconsGetName(conss[1]), "poly_1");
   cr_expect_str_eq(SCIPconsGetName(conss[2]), "polyfun_1_a_0");
   cr_expect_str_eq(SCIPconsGetName(conss[3]), "polyfun_1_b_1");
   cr_expect_str_eq(SCIPconsGetName(conss[4]), "polyfun_1");

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   sumhdlr = SCIPgetConsExprExprHdlrSum(conshdlr);
   cr_assert_not_null(sumhdlr);

   /* check the constraints which should be
    * quad: X2^2 + 4*X1^2 == 0;
    * poly: X2^2 + 2 * X1 * X2^3 - X1^2 >= -0.2;
    * <polyfun_1_a_0>:  X2^4 - @@polyfun_1_t_0 == 0.0
    * <polyfun_1_b_1>: -@@polyfun_1_r_1 + 0.434294 * ln(@@polyfun_1_t_0) == 0.0
    * <polyfun_1>:     2 * @@polyfun_1_r_1^2 + 2 * X2^3 * X1 + X1^4 <= 1.0
    */
   SCIP_Real infty = SCIPinfinity(scip);
   SCIP_Real lhs[5] = {0.0, -0.2,  0.0, 0.0, -infty};
   SCIP_Real rhs[5] = {0.0, infty, 0.0, 0.0, 1.0};
   int expectednnonz[5] = {2, 3, 2, 2, 3};
   SCIP_Real expectedcoeffs[5][3] = {{1.0, 4.0, infty}, {1.0, 2.0, -1.0}, {1.0, -1.0, infty},
                                      {-1.0, 1.0 / log(10.0), infty}, {2.0, 2.0, 1.0}};
   enum exprhdlrtype {POW, PRODUCT, VAR, LOG, NONE};
   enum exprhdlrtype exprhdlrtypes[5][3] = {{POW, POW, NONE}, {POW, PRODUCT, POW}, {POW, VAR, NONE}, {VAR, LOG, NONE},
                                            {POW, PRODUCT, POW}};
   SCIP_VAR* expectedvars[5][3] = {{vars[1], vars[0], NULL}, {vars[1], vars[0], vars[0]}, {vars[1], vars[3], NULL},
                                   {vars[4], vars[3], NULL}, {vars[4], vars[0], vars[0]}};
   SCIP_Real expectedexps[5][3] = {{2.0, 2.0, infty}, {2.0, infty, 2.0}, {4.0, infty, infty}, {infty, infty, infty},
                                   {2.0, infty, 4.0}};

   for( int i = 0; i < 5; ++i )
   {
      cr_assert_eq(conshdlr, SCIPconsGetHdlr(conss[i]));
      cr_assert(SCIPisEQ(scip, SCIPgetLhsConsExpr(scip, conss[i]), lhs[i]));
      cr_assert(SCIPisEQ(scip, SCIPgetRhsConsExpr(scip, conss[i]), rhs[i]));

      expr = SCIPgetExprConsExpr(scip, conss[i]);

      cr_assert_not_null(expr);
      cr_assert_eq(expr->exprhdlr, sumhdlr);
      cr_expect_eq(SCIPgetConsExprExprNChildren(expr), expectednnonz[i]);
      cr_expect_eq(SCIPgetConsExprExprSumConstant(expr), 0.0);
      for( int j = 0; j < expectednnonz[i]; ++j )
      {
         SCIP_CONSEXPR_EXPR* childexpr = NULL;
         SCIP_VAR* childvar;

         cr_expect_eq(SCIPgetConsExprExprSumCoefs(expr)[j], expectedcoeffs[i][j], "i,j = %d,%d: expected %g, got %g\n",
                      i, j, expectedcoeffs[i][j], SCIPgetConsExprExprSumCoefs(expr)[j]);
         switch( exprhdlrtypes[i][j] )
         {
         case POW:
            cr_expect_eq(expr->children[j]->exprhdlr, SCIPgetConsExprExprHdlrPower(conshdlr));
            cr_expect_eq(SCIPgetConsExprExprPowExponent(expr->children[j]), expectedexps[i][j]);
            childexpr = expr->children[j]->children[0];
            break;
         case PRODUCT:
            cr_expect_eq(expr->children[j]->exprhdlr, SCIPgetConsExprExprHdlrProduct(conshdlr));
            cr_expect_eq(SCIPgetConsExprExprNChildren(expr->children[j]), 2);

            /* there is only one product expression, check the power child here */
            /* for some reason the order of product children is non-deterministic, so check which child is power */
            childexpr = expr->children[j]->children[0]->exprhdlr == SCIPgetConsExprExprHdlrPower(conshdlr) ?
                        expr->children[j]->children[0] : expr->children[j]->children[1];

            cr_expect_eq(childexpr->exprhdlr, SCIPgetConsExprExprHdlrPower(conshdlr));
            cr_expect_eq(SCIPgetConsExprExprPowExponent(childexpr), 3);
            childvar = SCIPgetConsExprExprVarVar(childexpr->children[0]);
            cr_expect_eq(childvar, vars[1]);

            /* save the other expr to childexpr */
            childexpr = expr->children[j]->children[0]->exprhdlr == SCIPgetConsExprExprHdlrPower(conshdlr) ?
                        expr->children[j]->children[1] : expr->children[j]->children[0];
               break;
         case VAR:
            cr_expect_eq(expr->children[j]->exprhdlr, SCIPgetConsExprExprHdlrVar(conshdlr));
            childexpr = expr->children[j];
               break;
         case LOG:
            cr_expect_eq(expr->children[j]->exprhdlr, SCIPgetConsExprExprHdlrLog(conshdlr));
            childexpr = expr->children[j]->children[0];
               break;
         default:
            cr_assert_fail("\nshouldn't have reached this, i, j = %d, %d", i, j);
         }

         cr_assert(childexpr != NULL);
         cr_expect_eq(childexpr->exprhdlr, SCIPgetConsExprExprHdlrVar(conshdlr));
         cr_expect_eq(SCIPgetConsExprExprVarVar(childexpr), expectedvars[i][j]);
      }
   }
}
