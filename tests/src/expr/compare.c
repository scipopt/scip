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

/**@file   compare.c
 * @brief  tests expression comparison methods
 * @note   expressions should be simplified!
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/expr_sum.h"
#include "scip/expr_product.h"
#include "scip/expr_var.h"
#include "scip/expr_value.h"
#include "scip/expr_pow.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_EXPR* expr_x; /* x */
static SCIP_EXPR* expr_y; /* y */
static SCIP_EXPR* expr_z; /* z */
static SCIP_EXPR* expr_posvalue; /* 1.3 */
static SCIP_EXPR* expr_negvalue; /* -12.58 */
static SCIP_EXPR* expr_sum; /* x + y */
static SCIP_EXPR* expr_prod; /* x*y*(x+y) */
static SCIP_EXPR* expr_halfx; /* 0.5 * x */
static SCIP_EXPR* expr_sqrtx; /* \sqrt x */
static SCIP_EXPR* expr_half_sqrx; /* 0.5 * x^2 */
static SCIP_EXPR* expr_sqrx; /* x^2 */
static SCIP_EXPR* expr_signsqrx; /* signpower(x,2) */
static SCIP_EXPR* expr_sum_fracpow; /* (x+y)^3.2 */
static SCIP_EXPR* expr_subprod_fracpow; /* x*(y*(x+y))^2.8 */

static
void setup(void)
{
   SCIP_EXPR* tmp;
   SCIP_EXPR* tmp2;
   SCIP_Real aux;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 2.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -1.0, 2.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create expressions for negative and positive values */
   SCIP_CALL( SCIPcreateExprValue(scip, &expr_posvalue, 1.3, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprValue(scip, &expr_negvalue, -12.58, NULL, NULL) );

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_y, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_z, z, NULL, NULL) );

   /* create expression sum: x+y */
   SCIP_CALL( SCIPcreateExprSum(scip, &expr_sum, 1, &expr_x, NULL, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr_sum, expr_y, 1.0) );

   /* create expression sum_fracpow: (x+y)^3.2 */
   SCIP_CALL( SCIPcreateExprPow(scip, &expr_sum_fracpow, expr_sum, 3.2, NULL, NULL) );

   /* create expression halfx: 0.5*x */
   /* SCIP_CALL( SCIPcreateExprProduct(scip, conshdlr, &expr_halfx, 1, &expr_x, NULL, 0.5) ); */ /*this is not simplified*/
   aux = 0.5;
   SCIP_CALL( SCIPcreateExprSum(scip, &expr_halfx, 1, &expr_x, &aux, 0.0, NULL, NULL) );

   /* create expression sqrtx: sqrt(x) */
   SCIP_CALL( SCIPcreateExprPow(scip, &expr_sqrtx, expr_x, 0.5, NULL, NULL) );

   /* create expression sqrx: x^2 */
   SCIP_CALL( SCIPcreateExprPow(scip, &expr_sqrx, expr_x, 2.0, NULL, NULL) );

   /* create expression signsqrx: signpower(x,2) */
   SCIP_CALL( SCIPcreateExprSignpower(scip, &expr_signsqrx, expr_x, 2.0, NULL, NULL) );

   /* create expression half_srqx: 0.5*x^2 */
   aux = 0.5;
   SCIP_CALL( SCIPcreateExprSum(scip, &expr_half_sqrx, 1, &expr_sqrx, &aux, 0.0, NULL, NULL) );

   /* create expression prod: x*y*(x+y) */
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr_prod, 1, &expr_x, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_prod, expr_y) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_prod, expr_sum) );

   /* create expression subprod_fracpower: x*(y*(x+y))^2.8 */
   aux = 2.8;
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr_subprod_fracpow, 1, &expr_x, 1.0, NULL, NULL) );  /* x */
   SCIP_CALL( SCIPcreateExprProduct(scip, &tmp, 1, &expr_y, 1.0, NULL, NULL) );  /* y */
   SCIP_CALL( SCIPappendExprChild(scip, tmp, expr_sum) );  /* y -> y * (x+y) */
   SCIP_CALL( SCIPcreateExprPow(scip, &tmp2, tmp, aux, NULL, NULL) );  /* (y*(x+y))^aux */
   SCIP_CALL( SCIPappendExprChild(scip, expr_subprod_fracpow, tmp2) );   /* x -> x * (y*(x+y))^aux  */

   SCIP_CALL( SCIPreleaseExpr(scip, &tmp) );
   SCIP_CALL( SCIPreleaseExpr(scip, &tmp2) );
}

static
void teardown(void)
{
   /* release */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_z) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_negvalue) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_posvalue) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sum) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_prod) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_subprod_fracpow) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sum_fracpow) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_halfx) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sqrx) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_signsqrx) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_half_sqrx) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sqrtx) );

   /* release scip */
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   /* check for leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

TestSuite(compare, .init = setup, .fini = teardown);

Test(compare, linear_expressions)
{
   SCIP_EXPR* lin_expr1;
   SCIP_EXPR* lin_expr2;

   cr_expect( SCIPcompareExpr(scip, expr_x, expr_x) == 0 );
   cr_expect( SCIPcompareExpr(scip, expr_y, expr_y) == 0 );
   cr_expect( SCIPcompareExpr(scip, expr_x, expr_y) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_y, expr_x) == 1 );

   SCIP_CALL( SCIPcreateExprSum(scip, &lin_expr1, 1, &expr_x, NULL, 2.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &lin_expr2, 1, &expr_x, NULL, 3.0, NULL, NULL) );
   cr_expect( SCIPcompareExpr(scip, lin_expr1, lin_expr2) == -1 );
   cr_expect( SCIPcompareExpr(scip, lin_expr2, lin_expr1) == 1 );
   cr_expect( SCIPcompareExpr(scip, lin_expr1, lin_expr1) == 0 );

   SCIP_CALL( SCIPappendExprSumExpr(scip, lin_expr1, expr_y, 3.0) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, lin_expr2, expr_y, 3.0) );
   cr_expect( SCIPcompareExpr(scip, lin_expr1, lin_expr2) == -1 );
   cr_expect( SCIPcompareExpr(scip, lin_expr2, lin_expr1) == 1 );
   cr_expect( SCIPcompareExpr(scip, lin_expr1, lin_expr1) == 0 );

   SCIP_CALL( SCIPreleaseExpr(scip, &lin_expr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &lin_expr2) );
}

Test(compare, values)
{
   cr_expect( SCIPcompareExpr(scip, expr_negvalue, expr_negvalue) == 0 );
   cr_expect( SCIPcompareExpr(scip, expr_posvalue, expr_negvalue) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_negvalue, expr_posvalue) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_posvalue, expr_x) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_posvalue, expr_y) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_posvalue, expr_sum) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_posvalue, expr_prod) == -1 );
}

Test(compare, exponents)
{
   cr_expect( SCIPcompareExpr(scip, expr_halfx, expr_x) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_x, expr_halfx) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_sqrtx, expr_x) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_x, expr_sqrtx) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_sqrtx, expr_sqrx) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_sqrx, expr_sqrtx) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_x, expr_sqrx) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_sqrx, expr_x) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_prod, expr_sum_fracpow) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_sum_fracpow, expr_prod) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_sum_fracpow, expr_sum_fracpow) == 0 );

   /* compare signed and unsigned power (same exponents) */
   cr_expect( SCIPcompareExpr(scip, expr_sqrx, expr_signsqrx) == -1 );
}

Test(compare, sums_and_products)
{
   cr_expect( SCIPcompareExpr(scip, expr_sum, expr_prod) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_prod, expr_sum) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_prod, expr_prod) == 0 );
   cr_expect( SCIPcompareExpr(scip, expr_sum, expr_sum) == 0 );
   cr_expect( SCIPcompareExpr(scip, expr_x, expr_sum) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_y, expr_sum) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_x, expr_prod) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_y, expr_prod) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_sum, expr_x) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_sum, expr_y) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_prod, expr_x) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_prod, expr_y) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_prod, expr_subprod_fracpow) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_subprod_fracpow, expr_prod) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_subprod_fracpow, expr_subprod_fracpow) == 0 );
   cr_expect( SCIPcompareExpr(scip, expr_subprod_fracpow, expr_z) == -1 );
   cr_expect( SCIPcompareExpr(scip, expr_z, expr_subprod_fracpow) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_sqrx, expr_half_sqrx) == 1 );
   cr_expect( SCIPcompareExpr(scip, expr_half_sqrx, expr_sqrx) == -1 );
}
