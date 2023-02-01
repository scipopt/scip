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

/**@file   hash.c
 * @brief  tests hashes of expressions
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/expr_sum.h"
#include "scip/expr_product.h"
#include "scip/expr_var.h"
#include "scip/expr_value.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_EXPR* xexpr;
static SCIP_EXPR* yexpr;

/* helper function for testing */
static
void checkHashkey(
   SCIP_EXPR*   expr1,              /**< first expression to be tested */
   SCIP_EXPR*   expr2               /**< second expression to be tested (might be NULL) */
   )
{
   unsigned int hashkey1;
   unsigned int hashkey2;
   cr_assert_not_null(expr1);

   SCIPinfoMessage(scip, NULL, "hash key of expression: ");
   SCIP_CALL( SCIPprintExpr(scip, expr1, NULL) );
   SCIP_CALL( SCIPhashExpr(scip, expr1, &hashkey1));
   SCIPinfoMessage(scip, NULL, " = %u\n", hashkey1);
   cr_assert_neq(hashkey1, 0);

   if( expr2 != NULL )
   {
      SCIPinfoMessage(scip, NULL, "hash key of expression: ");
      SCIP_CALL( SCIPprintExpr(scip, expr2, NULL) );
      SCIP_CALL( SCIPhashExpr(scip, expr2, &hashkey2));
      SCIPinfoMessage(scip, NULL, " = %u\n", hashkey2);
      cr_assert_eq(hashkey1, hashkey2);
      cr_assert_neq(hashkey2, 0);
   }
}

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &yexpr, y, NULL, NULL) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

TestSuite(hash, .init = setup, .fini = teardown);

Test(hash, hashEqualExpressions)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* expr2;

   /* sum expressions */
   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 0, NULL, NULL, 2.5, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr, xexpr, 14.3) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr, yexpr, -2.3) );

   SCIP_CALL( SCIPcreateExprSum(scip, &expr2, 0, NULL, NULL, 2.5, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr2, yexpr, -2.3) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr2, xexpr, 14.3) );

   checkHashkey(expr, expr2);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* product expressions */
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr, 0, NULL, 2.5, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr, xexpr) );
   SCIP_CALL( SCIPappendExprChild(scip, expr, yexpr) );

   SCIP_CALL( SCIPcreateExprProduct(scip, &expr2, 0, NULL, 2.5, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr, yexpr) );
   SCIP_CALL( SCIPappendExprChild(scip, expr, xexpr) );

   checkHashkey(expr, expr2);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   SCIP_CALL( (SCIPparseExpr(scip, &expr, "2*<x> + 3.3*<y>", NULL, NULL, NULL)) );
   SCIP_CALL( (SCIPparseExpr(scip, &expr2, "3.3*<y> + 2*<x>", NULL, NULL, NULL)) );

   checkHashkey(expr, expr2);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(hash, hashSingleExpr)
{
   SCIP_EXPR* expr;
   int i;

   /* value expressions */
   for( i = -2; i < 2; ++i )
   {
      SCIP_CALL( SCIPcreateExprValue(scip, &expr, i, NULL, NULL) );
      checkHashkey(expr, NULL);
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }

   /* variable expressions */
   checkHashkey(xexpr, NULL);
   checkHashkey(yexpr, NULL);

   /* absolute expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "abs(<x>)", NULL, NULL, NULL) );
   checkHashkey(expr, NULL);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* logarithmic expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "exp(<x>)", NULL, NULL, NULL) );
   checkHashkey(expr, NULL);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* exponential expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "log(<x>)", NULL, NULL, NULL) );
   checkHashkey(expr, NULL);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* complicated expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "abs(exp(<x>*<y>^2/<x>^4) - log(2*<x>)*(3+5*<x>-2*<y>)*(<x>+<y>)^(-3.5)) + 2", NULL, NULL, NULL) );
   checkHashkey(expr, NULL);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
