/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   readers.c
 * @brief  tests readers that create nonlinear constraints
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

void test_readers_pip(void)
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_EXPR* expr;
   SCIP_EXPR** children;
   SCIP_EXPR** grandchildren;
   char filename[SCIP_MAXSTRLEN];

   /* get file to read: test.mps that lives in the same directory as this file */
   TESTsetTestfilename(filename, __FILE__, "test.pip");
   printf("Reading %s\n", filename);

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPreadProb(scip, filename, NULL));

   /* check that vars are what we expect */
   SOFT_ASSERT_EQUAL(SCIPgetNVars(scip), 5);
   vars = SCIPgetVars(scip);

   SOFT_ASSERT_EQUAL_STRING("x", SCIPvarGetName(vars[0]));
   SOFT_ASSERT_EQUAL_STRING("y", SCIPvarGetName(vars[1]));
   SOFT_ASSERT_EQUAL_STRING("z", SCIPvarGetName(vars[2]));
   SOFT_ASSERT_EQUAL_STRING("objconst", SCIPvarGetName(vars[3]));
   SOFT_ASSERT_EQUAL_STRING("nonlinobjvar", SCIPvarGetName(vars[4]));

   /* check that cons are what we expect */
   SOFT_ASSERT_EQUAL(SCIPgetNConss(scip), 3);
   conss = SCIPgetConss(scip);

   SOFT_ASSERT_EQUAL_STRING("nonlinobj", SCIPconsGetName(conss[0]));
   SOFT_ASSERT_EQUAL_STRING("e1", SCIPconsGetName(conss[1]));
   SOFT_ASSERT_EQUAL_STRING("e2", SCIPconsGetName(conss[2]));

   /* check objective coefficients */
   SOFT_ASSERT_EQUAL(SCIPvarGetObj(vars[0]), 0.0);
   SOFT_ASSERT_EQUAL(SCIPvarGetObj(vars[1]), 0.0);
   SOFT_ASSERT_EQUAL(SCIPvarGetObj(vars[2]), 0.0);
   SOFT_ASSERT_EQUAL(SCIPvarGetObj(vars[3]), 11.0);
   SOFT_ASSERT_EQUAL(SCIPvarGetObj(vars[4]), 1.0);

   /*
    * check whether the first constraint is nonlinear and of the form + x - y^3 z^2 - nonlinobj <= 0
    */

   SOFT_ASSERT_EQUAL(SCIPconsGetHdlr(conss[0]), SCIPfindConshdlr(scip, "nonlinear"));
   expr = SCIPgetExprNonlinear(conss[0]);
   SOFT_ASSERT_NOT_NULL(expr);
   children = SCIPexprGetChildren(expr);
   SOFT_ASSERT_NOT_NULL(children);
   SOFT_ASSERT_EQUAL(SCIPexprGetNChildren(expr), 3);

   /* check sides */
   SOFT_ASSERT(SCIPisInfinity(scip, -SCIPgetLhsNonlinear(conss[0])));
   SOFT_ASSERT_EQUAL(SCIPgetRhsNonlinear(conss[0]), 0.0);

   /* check constant and coefficients of the sum expression */
   SOFT_ASSERT(SCIPisExprSum(scip, expr));
   SOFT_ASSERT_EQUAL(SCIPgetConstantExprSum(expr), 0.0);
   SOFT_ASSERT_EQUAL(SCIPgetCoefsExprSum(expr)[0], 1.0);
   SOFT_ASSERT_EQUAL(SCIPgetCoefsExprSum(expr)[1], -1.0);
   SOFT_ASSERT_EQUAL(SCIPgetCoefsExprSum(expr)[2], -1.0);

   /* check whether first and third child is a variable */
   SOFT_ASSERT(SCIPisExprVar(scip, children[0]));
   SOFT_ASSERT(SCIPisExprVar(scip, children[2]));

   /* check whether second child is a product; both grandchildren need to be power expressions */
   SOFT_ASSERT(SCIPisExprProduct(scip, children[1]));
   grandchildren = SCIPexprGetChildren(children[1]);
   SOFT_ASSERT_EQUAL(SCIPexprGetNChildren(children[1]), 2);
   SOFT_ASSERT(SCIPisExprPower(scip, grandchildren[0]));
   SOFT_ASSERT(SCIPisExprPower(scip, grandchildren[1]));
   SOFT_ASSERT_EQUAL(SCIPgetExponentExprPow(grandchildren[0]), 3.0);
   SOFT_ASSERT_EQUAL(SCIPgetExponentExprPow(grandchildren[1]), 2.0);

   /*
    * check whether the second constraint is linear and of the form x + 2y + 3z <= 1
    */

   SOFT_ASSERT_EQUAL(SCIPconsGetHdlr(conss[1]), SCIPfindConshdlr(scip, "linear"));
   SOFT_ASSERT(SCIPisInfinity(scip, -SCIPgetLhsLinear(scip, conss[1])));
   SOFT_ASSERT_EQUAL(SCIPgetRhsLinear(scip, conss[1]), 1.0);
   SOFT_ASSERT_EQUAL(SCIPgetNVarsLinear(scip, conss[1]), 3);
   SOFT_ASSERT_EQUAL(SCIPgetValsLinear(scip, conss[1])[0], 1.0);
   SOFT_ASSERT_EQUAL(SCIPgetValsLinear(scip, conss[1])[1], 2.0);
   SOFT_ASSERT_EQUAL(SCIPgetValsLinear(scip, conss[1])[2], 3.0);

   /*
    * check whether the third constraint is nonlinear and of the form x^2 y^3 z^4 + x + y + 2 = 10
    */

   SOFT_ASSERT_EQUAL(SCIPconsGetHdlr(conss[2]), SCIPfindConshdlr(scip, "nonlinear"));
   expr = SCIPgetExprNonlinear(conss[2]);
   SOFT_ASSERT_NOT_NULL(expr);
   children = SCIPexprGetChildren(expr);
   SOFT_ASSERT_NOT_NULL(children);
   SOFT_ASSERT_EQUAL(SCIPexprGetNChildren(expr), 3);

   /* check sides */
   SOFT_ASSERT_EQUAL(SCIPgetLhsNonlinear(conss[2]), 10.0);
   SOFT_ASSERT_EQUAL(SCIPgetRhsNonlinear(conss[2]), 10.0);

   /* check constant and coefficients of the sum expression */
   SOFT_ASSERT(SCIPisExprSum(scip, expr));
   SOFT_ASSERT_EQUAL(SCIPgetConstantExprSum(expr), 7.0);
   SOFT_ASSERT_EQUAL(SCIPgetCoefsExprSum(expr)[0], 4.0);
   SOFT_ASSERT_EQUAL(SCIPgetCoefsExprSum(expr)[1], 5.0);
   SOFT_ASSERT_EQUAL(SCIPgetCoefsExprSum(expr)[2], 6.0);

   /* check whether first child is a product; all three grandchildren need to be power expressions */
   SOFT_ASSERT(SCIPisExprProduct(scip, children[0]));
   grandchildren = SCIPexprGetChildren(children[0]);
   SOFT_ASSERT_EQUAL(SCIPexprGetNChildren(children[0]), 3);
   SOFT_ASSERT(SCIPisExprPower(scip, grandchildren[0]));
   SOFT_ASSERT(SCIPisExprPower(scip, grandchildren[1]));
   SOFT_ASSERT(SCIPisExprPower(scip, grandchildren[2]));
   SOFT_ASSERT_EQUAL(SCIPgetExponentExprPow(grandchildren[0]), 2.0);
   SOFT_ASSERT_EQUAL(SCIPgetExponentExprPow(grandchildren[1]), 3.0);
   SOFT_ASSERT_EQUAL(SCIPgetExponentExprPow(grandchildren[2]), 4.0);
}

void test_readers_mps1(void)
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_EXPR* expr;
   char filename[SCIP_MAXSTRLEN];

   /* get file to read: test.mps that lives in the same directory as this file */
   TESTsetTestfilename(filename, __FILE__, "test.mps");
   printf("Reading %s\n", filename);

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPreadProb(scip, filename, NULL));

   /* check that vars are what we expect */
   SOFT_ASSERT_EQUAL(SCIPgetNVars(scip), 4);
   vars = SCIPgetVars(scip);

   SOFT_ASSERT_EQUAL_STRING("x1", SCIPvarGetName(vars[0]));
   SOFT_ASSERT_EQUAL_STRING("x2", SCIPvarGetName(vars[1]));
   SOFT_ASSERT_EQUAL_STRING("x3", SCIPvarGetName(vars[2]));
   SOFT_ASSERT_EQUAL_STRING("qmatrixvar", SCIPvarGetName(vars[3]));

   /* check that cons are what we expect */
   SOFT_ASSERT_EQUAL(SCIPgetNConss(scip), 4);
   conss = SCIPgetConss(scip);

   SOFT_ASSERT_EQUAL_STRING("c0", SCIPconsGetName(conss[0]));
   SOFT_ASSERT_EQUAL_STRING("c1", SCIPconsGetName(conss[1]));
   SOFT_ASSERT_EQUAL_STRING("c2", SCIPconsGetName(conss[2]));
   SOFT_ASSERT_EQUAL_STRING("qmatrix", SCIPconsGetName(conss[3]));

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
      SOFT_ASSERT_EQUAL(SCIPfindConshdlr(scip, "nonlinear"), SCIPconsGetHdlr(conss[i]));
      SOFT_ASSERT_EQUAL(SCIPgetLhsNonlinear(conss[i]), lhs[i], "lhs cons %d: expected %g, got %g\n", i, lhs[i], SCIPgetLhsNonlinear(conss[i]));
      SOFT_ASSERT_EQUAL(SCIPgetRhsNonlinear(conss[i]), rhs[i], "rhs cons %d: expected %g, got %g\n", i, rhs[i], SCIPgetRhsNonlinear(conss[i]));

      expr = SCIPgetExprNonlinear(conss[i]);
      TEST_ASSERT_NOT_NULL(expr);
      SOFT_ASSERT(SCIPisExprSum(scip, expr));
      SOFT_ASSERT_EQUAL(SCIPexprGetNChildren(expr), expectednnonz[i]);
      SOFT_ASSERT_EQUAL(SCIPgetConstantExprSum(expr), 0.0);
      for( int j = 0; j < expectednnonz[i]; ++j )
         SOFT_ASSERT_EQUAL(SCIPgetCoefsExprSum(expr)[j], expectedcoeffs[i][j], "i,j = %d,%d: expected %g, got %g\n",
               i, j, expectedcoeffs[i][j], SCIPgetCoefsExprSum(expr)[j]);
   }

   /* check constraint c1 */
   expr = SCIPgetExprNonlinear(conss[1]);
   SOFT_ASSERT(SCIPisExprPower(scip, SCIPexprGetChildren(expr)[0]));;
   SOFT_ASSERT(SCIPisExprVar(scip, SCIPexprGetChildren(SCIPexprGetChildren(expr)[0])[0]));
   SOFT_ASSERT_EQUAL(SCIPgetVarExprVar(SCIPexprGetChildren(SCIPexprGetChildren(expr)[0])[0]), vars[2]);

   SOFT_ASSERT(SCIPisExprVar(scip, SCIPexprGetChildren(expr)[1]));
   SOFT_ASSERT_EQUAL(SCIPgetVarExprVar(SCIPexprGetChildren(expr)[1]), vars[1]);
}

void test_readers_zimpl(void)
{
   SCIP* scip;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_EXPR* expr;
   char filename[SCIP_MAXSTRLEN];

   /* get file to read: test.zpl that lives in the same directory as this file */
   TESTsetTestfilename(filename, __FILE__, "test.zpl");
   printf("Reading %s\n", filename);

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   if( SCIPfindReader(scip, "zplreader") == NULL )
      return;

   SCIP_CALL( SCIPreadProb(scip, filename, NULL));

   /* check that vars are what we expect; zimpl will create 2 auxiliary variables, hence we expect 5 */
   SOFT_ASSERT_EQUAL(SCIPgetNVars(scip), 5, "\nexpected 5 variables, got %d", SCIPgetNVars(scip));
   vars = SCIPgetVars(scip);

   SOFT_ASSERT_EQUAL_STRING("X1", SCIPvarGetName(vars[0]));
   SOFT_ASSERT_EQUAL_STRING("X2", SCIPvarGetName(vars[1]));
   SOFT_ASSERT_EQUAL_STRING("X3", SCIPvarGetName(vars[2]));
   SOFT_ASSERT_EQUAL_STRING("@@polyfun_1_t_0", SCIPvarGetName(vars[3]));
   SOFT_ASSERT_EQUAL_STRING("@@polyfun_1_r_1", SCIPvarGetName(vars[4]));

   /* check that cons are what we expect:  */
   SOFT_ASSERT_EQUAL(SCIPgetNConss(scip), 5);
   conss = SCIPgetConss(scip);

   SOFT_ASSERT_EQUAL_STRING("quad_1", SCIPconsGetName(conss[0]));
   SOFT_ASSERT_EQUAL_STRING("poly_1", SCIPconsGetName(conss[1]));
   SOFT_ASSERT_EQUAL_STRING("polyfun_1_a_0", SCIPconsGetName(conss[2]));
   SOFT_ASSERT_EQUAL_STRING("polyfun_1_b_1", SCIPconsGetName(conss[3]));
   SOFT_ASSERT_EQUAL_STRING("polyfun_1", SCIPconsGetName(conss[4]));

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
      TEST_ASSERT_EQUAL(SCIPfindConshdlr(scip, "nonlinear"), SCIPconsGetHdlr(conss[i]));
      TEST_ASSERT(SCIPisEQ(scip, SCIPgetLhsNonlinear(conss[i]), lhs[i]));
      TEST_ASSERT(SCIPisEQ(scip, SCIPgetRhsNonlinear(conss[i]), rhs[i]));

      expr = SCIPgetExprNonlinear(conss[i]);
      TEST_ASSERT_NOT_NULL(expr);

      SOFT_ASSERT(SCIPisExprSum(scip, expr));
      SOFT_ASSERT_EQUAL(SCIPexprGetNChildren(expr), expectednnonz[i]);
      SOFT_ASSERT_EQUAL(SCIPgetConstantExprSum(expr), 0.0);
      for( int j = 0; j < expectednnonz[i]; ++j )
      {
         SCIP_EXPR* childj;
         SCIP_EXPR* childexpr = NULL;
         SCIP_VAR* childvar;

         SOFT_ASSERT_EQUAL(SCIPgetCoefsExprSum(expr)[j], expectedcoeffs[i][j], "i,j = %d,%d: expected %g, got %g\n",
                      i, j, expectedcoeffs[i][j], SCIPgetCoefsExprSum(expr)[j]);
         childj = SCIPexprGetChildren(expr)[j];
         switch( exprhdlrtypes[i][j] )
         {
            case POW:
               SOFT_ASSERT(SCIPisExprPower(scip, childj));
               SOFT_ASSERT_EQUAL(SCIPgetExponentExprPow(childj), expectedexps[i][j]);
               childexpr = SCIPexprGetChildren(childj)[0];
               break;
            case PRODUCT:
               SOFT_ASSERT(SCIPisExprProduct(scip, childj));
               SOFT_ASSERT_EQUAL(SCIPexprGetNChildren(childj), 2);

               /* there is only one product expression, check the power child here */
               /* for some reason the order of product children is non-deterministic, so check which child is power */
               childexpr = SCIPisExprPower(scip, SCIPexprGetChildren(childj)[0]) ?
                  SCIPexprGetChildren(childj)[0] : SCIPexprGetChildren(childj)[1];

               SOFT_ASSERT(SCIPisExprPower(scip, childexpr));
               SOFT_ASSERT_EQUAL(SCIPgetExponentExprPow(childexpr), 3);
               childvar = SCIPgetVarExprVar(SCIPexprGetChildren(childexpr)[0]);
               SOFT_ASSERT_EQUAL(childvar, vars[1]);

               /* save the other expr to childexpr */
               childexpr = SCIPisExprPower(scip, SCIPexprGetChildren(childj)[0]) ?
                  SCIPexprGetChildren(childj)[1] : SCIPexprGetChildren(childj)[0];
               break;
            case VAR:
               SOFT_ASSERT(SCIPisExprVar(scip, childj));
               childexpr = childj;
               break;
            case LOG:
               SOFT_ASSERT(SCIPisExprLog(scip, childj));
               childexpr = SCIPexprGetChildren(childj)[0];
               break;
            default:
               TEST_FAIL_MESSAGE(scip_test_msg("\nshouldn't have reached this, i, j = %d, %d", i, j));
         }

         TEST_ASSERT(childexpr != NULL);
         SOFT_ASSERT(SCIPisExprVar(scip, childexpr));
         SOFT_ASSERT_EQUAL(SCIPgetVarExprVar(childexpr), expectedvars[i][j]);
      }
   }
}

void setUp(void) { }

void tearDown(void) { }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_readers_pip);
   RUN_TEST(test_readers_mps1);
   RUN_TEST(test_readers_zimpl);
   return UNITY_END();
}
