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

/**@file   free.c
 * @brief  unit test for checking freeing expressions
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/expr_var.h"
#include "scip/expr_sum.h"
#include "scip/expr_product.h"
#include "scip/expr_value.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
#define BIG 100

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* create expression x + (-2)*x/y*(-5) */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

/***** TEST SUITE: all tests of the form Test(free, xyz) belong to the same suite and share the setup and teardown *****/
TestSuite(free, .init = setup, .fini = teardown);

/***** ACTUAL TESTS *****/

Test(free, simple_5xy)
{
   SCIP_EXPR* expr_x;
   SCIP_EXPR* expr_y;
   SCIP_EXPR* expr_5;
   SCIP_EXPR* expr_xy5;

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_y, y, NULL, NULL) );

   /* create expression for constant 5 */
   SCIP_CALL( SCIPcreateExprValue(scip, &expr_5, 5.0, NULL, NULL) );

   /* create expression for product of 5, x, and y */
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr_xy5, 1, &expr_x, 2.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_xy5, expr_y) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_xy5, expr_5) );

   /* release leaf expressions (this should not free them yet, as they are captured by prod_xy5) */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_5) );

   /* release product expression (this should free the product and its children) */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_xy5) );
   /* in the teardown we check, after scip is freed, that there are no leaks */
}

Test(free, long_expr)
{
   int i;
   SCIP_EXPR* exprs[BIG];
   SCIP_Real  coefs[BIG];
   SCIP_EXPR* expr_x;
   SCIP_EXPR* sumexpr;

   /* create expressions for variables x */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );

   for( i = 0; i < BIG; i++ )
   {
      exprs[i] = expr_x;
      coefs[i] = i;
   }

   /* create expression for sum */
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, BIG, exprs, coefs, -1.0, NULL, NULL) );

   cr_expect_eq(SCIPexprGetNUses(expr_x), BIG + 1);

   /* release leaf expressions (this should not free them yet, as they are captured by sumexpr) */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );

   cr_expect_eq(SCIPexprGetNUses(exprs[0]), BIG);

   /* release sum expression (this should free the sum and its children) */
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   /* in the teardown we check, after scip is freed, that there are no leaks */
}

Test(free, deep_expr)
{
   int i;
   SCIP_EXPR* sumexprs[BIG];
   SCIP_EXPR* expr_x;
   SCIP_EXPR* sumexpr;

   /* create expressions for variables x */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );

   SCIP_CALL( SCIPcreateExprSum(scip, &sumexprs[0], 1, &expr_x, NULL, 0.0, NULL, NULL) );
   for( i = 1; i < BIG; i++ )
   {
      /* create expressions for sum */
      SCIP_CALL( SCIPcreateExprSum(scip, &sumexprs[i], 1, &sumexprs[i-1], NULL, 1.0 * i, NULL, NULL) );
      cr_expect_eq(SCIPexprGetNUses(sumexprs[i]), 1);
      cr_expect_eq(SCIPexprGetNUses(sumexprs[i-1]), 2);
   }
   printf("finish big loop\n");

   /* create expressions for sum */
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 1, &sumexprs[BIG-1], NULL, 1.0 * BIG, NULL, NULL) );

   /* check nuses */
   cr_expect_eq(SCIPexprGetNUses(sumexpr), 1);
   cr_expect_eq(SCIPexprGetNUses(expr_x), 2);

   /* release leaf expressions (this should not free them yet, as they are captured by sumexpr) */
   for( i = 0; i < BIG; i++ )
   {
      cr_expect_eq(SCIPexprGetNUses(sumexprs[i]), 2);
      SCIP_CALL( SCIPreleaseExpr(scip, &sumexprs[i]) );
   }
   cr_expect_eq(SCIPexprGetNUses(expr_x), 2);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );

   /* release sum expression (this should free all expressions) */
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   /* in the teardown we check, after scip is freed, that there are no leaks */
}

Test(free, long_and_deep_expr)
{
   SCIP_EXPR* expr_x;
   SCIP_EXPR* expr_y;
   SCIP_EXPR* prodexprs[BIG];
   SCIP_EXPR* xysum[3];
   SCIP_EXPR* exprs[BIG];
   SCIP_EXPR* sumexpr;
   SCIP_EXPR* crazyexpr;
   SCIP_Real  coefs[BIG];
   int i;

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_y, y, NULL, NULL) );

   /* create long sum */
   for( i = 0; i < BIG; i++ )
   {
      exprs[i] = expr_x;
      coefs[i] = i;
   }
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, BIG, exprs, coefs, -1.0, NULL, NULL) );

   /* create deep product of ys and sums */
   xysum[0] = expr_x; xysum[1] = expr_y; xysum[2] = sumexpr;
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexprs[0], 3, (SCIP_EXPR**)&xysum, 1.0, NULL, NULL) );
   for( i = 1; i < BIG; i++ )
   {
      if( BIG % 2 == 1 )
      {
         xysum[0] = sumexpr; xysum[1] = prodexprs[i-1]; xysum[2] = expr_y;
         SCIP_CALL( SCIPcreateExprProduct(scip, &prodexprs[i], 3, (SCIP_EXPR**)&xysum, 1.0 * i, NULL, NULL) );
      }
      else
      {
         exprs[0] = prodexprs[i-1]; exprs[1] = prodexprs[i-1]; exprs[2] = expr_y;
         SCIP_CALL( SCIPcreateExprSum(scip, &prodexprs[i], BIG, exprs, coefs, -1.0, NULL, NULL) );
      }
   }

   /* create expressions for crazy expr */
   xysum[0] = sumexpr; xysum[1] = prodexprs[BIG-1]; xysum[2] = sumexpr;
   SCIP_CALL( SCIPcreateExprSum(scip, &crazyexpr, 3, (SCIP_EXPR**)&xysum, NULL, 1.0 * BIG, NULL, NULL) );

   /* release leaf expressions (this should not free them yet, as they are captured by crazyexpr) */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   for( i = 0; i < BIG; i++ )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &prodexprs[i]) );
   }

   /* release crazy expression (this should free the product and its children) */
   SCIP_CALL( SCIPreleaseExpr(scip, &crazyexpr) );
}
