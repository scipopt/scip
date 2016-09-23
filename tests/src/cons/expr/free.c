/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   free.c
 * @brief  unit test for checking freeing expressions
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_value.h"

/* TESTS */
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static const int BIG = 100;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* currently expr constraints cannot be created */
   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

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

/***** TEST SUITE: all tests of the form Test(walk, xxx) belong to the same suite and share the setup and teardown *****/
TestSuite(free, .init = setup, .fini = teardown);

/***** ACTUAL TESTS *****/

Test(free, simple_5xy)
{
   SCIP_CONSEXPR_EXPR* expr_x;
   SCIP_CONSEXPR_EXPR* expr_y;
   SCIP_CONSEXPR_EXPR* expr_5;
   SCIP_CONSEXPR_EXPR* expr_xy5;

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );

   /* create expression for constant 5 */
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr_5, 5.0) );

   /* create expression for product of 5, x, and y */
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_xy5, 1, &expr_x, 2.0) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_xy5, expr_y) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_xy5, expr_5) );

   /* release leaf expressions (this should not free them yet, as they are captured by prod_xy5) */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_5) );

   /* release product expression (this should free the product and its children) */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_xy5) );
   /* in the teardown we check, after scip is freed, that there are no leaks */
}

Test(free, long_expr)
{
   int i;
   SCIP_CONSEXPR_EXPR* exprs[BIG];
   SCIP_Real           coefs[BIG];
   SCIP_CONSEXPR_EXPR* expr_x;
   SCIP_CONSEXPR_EXPR* sumexpr;

   /* create expressions for variables x */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );

   for( i = 0; i < BIG; i++ )
   {
      exprs[i] = expr_x;
      coefs[i] = i;
   }

   /* create expression for sum */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, BIG, exprs, coefs, -1.0) );

   cr_expect_eq(SCIPgetConsExprExprNUses(expr_x), BIG + 1);

   /* release leaf expressions (this should not free them yet, as they are captured by sumexpr) */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );

   cr_expect_eq(SCIPgetConsExprExprNUses(exprs[0]), BIG);

   /* release sum expression (this should free the sum and its children) */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   /* in the teardown we check, after scip is freed, that there are no leaks */
}

Test(free, deep_expr)
{
   int i;
   SCIP_CONSEXPR_EXPR* sumexprs[BIG];
   SCIP_CONSEXPR_EXPR* expr_x;
   SCIP_CONSEXPR_EXPR* sumexpr;

   /* create expressions for variables x */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexprs[0], 1, &expr_x, NULL, 0.0) );
   for( i = 1; i < BIG; i++ )
   {
      /* create expressions for sum */
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexprs[i], 1, &sumexprs[i-1], NULL, 1.0 * i) );
      cr_expect_eq(SCIPgetConsExprExprNUses(sumexprs[i]), 1);
      cr_expect_eq(SCIPgetConsExprExprNUses(sumexprs[i-1]), 2);
   }
   printf("finish big loop\n");

   /* create expressions for sum */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, 1, &sumexprs[BIG-1], NULL, 1.0 * BIG) );

   /* check nuses */
   cr_expect_eq(SCIPgetConsExprExprNUses(sumexpr), 1);
   cr_expect_eq(SCIPgetConsExprExprNUses(expr_x), 2);

   /* release leaf expressions (this should not free them yet, as they are captured by sumexpr) */
   for( i = 0; i < BIG; i++ )
   {
      cr_expect_eq(SCIPgetConsExprExprNUses(sumexprs[i]), 2);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexprs[i]) );
   }
   cr_expect_eq(SCIPgetConsExprExprNUses(expr_x), 2);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );

   /* release sum expression (this should free all expressions) */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   /* in the teardown we check, after scip is freed, that there are no leaks */
}

Test(free, long_and_deep_expr)
{
   SCIP_CONSEXPR_EXPR* expr_x;
   SCIP_CONSEXPR_EXPR* expr_y;
   SCIP_CONSEXPR_EXPR* prodexprs[BIG];
   SCIP_CONSEXPR_EXPR* xysum[3];
   SCIP_CONSEXPR_EXPR* exprs[BIG];
   SCIP_CONSEXPR_EXPR* sumexpr;
   SCIP_CONSEXPR_EXPR* crazyexpr;
   SCIP_Real           coefs[BIG];
   int i;

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );

   /* create long sum */
   for( i = 0; i < BIG; i++ )
   {
      exprs[i] = expr_x;
      coefs[i] = i;
   }
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, BIG, exprs, coefs, -1.0) );

   /* create deep product of ys and sums */
   xysum[0] = expr_x; xysum[1] = expr_y; xysum[2] = sumexpr;
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &prodexprs[0], 3, (SCIP_CONSEXPR_EXPR**)&xysum, 1.0) );
   for( i = 1; i < BIG; i++ )
   {
      if( BIG % 2 == 1 )
      {
         xysum[0] = sumexpr; xysum[1] = prodexprs[i-1]; xysum[2] = expr_y;
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &prodexprs[i], 3, (SCIP_CONSEXPR_EXPR**)&xysum, 1.0 * i) );
      }
      else
      {
         exprs[0] = prodexprs[i-1]; exprs[1] = prodexprs[i-1]; exprs[2] = expr_y;
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &prodexprs[i], BIG, exprs, coefs, -1.0) );
      }
   }

   /* create expressions for crazy expr */
   xysum[0] = sumexpr; xysum[1] = prodexprs[BIG-1]; xysum[2] = sumexpr;
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &crazyexpr, 3, (SCIP_CONSEXPR_EXPR**)&xysum,
            NULL, 1.0 * BIG) );

   /* release leaf expressions (this should not free them yet, as they are captured by crazyexpr) */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   for( i = 0; i < BIG; i++ )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexprs[i]) );
   }

   /* release crazy expression (this should free the product and its children) */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &crazyexpr) );
}
