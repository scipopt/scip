/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   hash.c
 * @brief  tests hashes of expressions
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_abs.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_CONSEXPR_EXPR* xexpr;
static SCIP_CONSEXPR_EXPR* yexpr;

/* helper function for testing */
static
void checkHashkey(
   SCIP_CONSEXPR_EXPR*   expr1,              /**< first expression to be tested */
   SCIP_CONSEXPR_EXPR*   expr2               /**< second expression to be tested (might be NULL) */
   )
{
   unsigned int hashkey1;
   unsigned int hashkey2;
   cr_assert_not_null(expr1);

   SCIPinfoMessage(scip, NULL, "hash key of expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr1, NULL) );
   SCIP_CALL( SCIPgetConsExprExprHashkey(scip, expr1, &hashkey1));
   SCIPinfoMessage(scip, NULL, " = %u\n", hashkey1);
   cr_assert_neq(hashkey1, 0);

   if( expr2 != NULL )
   {
      SCIPinfoMessage(scip, NULL, "hash key of expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr2, NULL) );
      SCIP_CALL( SCIPgetConsExprExprHashkey(scip, expr2, &hashkey2));
      SCIPinfoMessage(scip, NULL, " = %u\n", hashkey2);
      cr_assert_eq(hashkey1, hashkey2);
      cr_assert_neq(hashkey2, 0);
   }
}

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

TestSuite(hash, .init = setup, .fini = teardown);

Test(hash, hashEqualExpressions)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* expr2;

   /* sum expressions */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr, 0, NULL, NULL, 2.5) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, xexpr, 14.3) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, yexpr, -2.3) );

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr2, 0, NULL, NULL, 2.5) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr2, yexpr, -2.3) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr2, xexpr, 14.3) );

   checkHashkey(expr, expr2);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* product expressions */
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr, 0, NULL, 2.5) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, xexpr) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, yexpr) );

   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr2, 0, NULL, 2.5) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, yexpr) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, xexpr) );

   checkHashkey(expr, expr2);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "2*<x> + 3.3*<y>", NULL, &expr)) );
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "3.3*<y> + 2*<x>", NULL, &expr2)) );

   checkHashkey(expr, expr2);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(hash, hashSingleExpr)
{
   SCIP_CONSEXPR_EXPR* expr;
   int i;

   /* value expressions */
   for( i = -2; i < 2; ++i )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr, i) );
      checkHashkey(expr, NULL);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* variable expressions */
   checkHashkey(xexpr, NULL);
   checkHashkey(yexpr, NULL);

   /* absolute expression */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "abs(<x>)", NULL, &expr)) );
   checkHashkey(expr, NULL);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* logarithmic expression */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "exp(<x>)", NULL, &expr)) );
   checkHashkey(expr, NULL);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* exponential expression */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "log(<x>)", NULL, &expr)) );
   checkHashkey(expr, NULL);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* complicated expression */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "abs(exp(<x>*<y>^2/<x>^4) - log(2*<x>)*(3+5*<x>-2*<y>)*(<x>+<y>)^(-3.5)) + 2", NULL, &expr)) );
   checkHashkey(expr, NULL);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
