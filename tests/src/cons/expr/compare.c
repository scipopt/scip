/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   compare.c
 * @brief  tests comparison methods
 * @note   expressions should be simplified!
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_pow.h"


#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_CONSEXPR_EXPR* expr_x; /* x */
static SCIP_CONSEXPR_EXPR* expr_y; /* y */
static SCIP_CONSEXPR_EXPR* expr_z; /* z */
static SCIP_CONSEXPR_EXPR* expr_posvalue; /* 1.3 */
static SCIP_CONSEXPR_EXPR* expr_negvalue; /* -12.58 */
static SCIP_CONSEXPR_EXPR* expr_sum; /* x + y */
static SCIP_CONSEXPR_EXPR* expr_prod; /* x */
static SCIP_CONSEXPR_EXPR* expr_halfx; /* 0.5 * x */
static SCIP_CONSEXPR_EXPR* expr_sqrtx; /* \sqrt x */
static SCIP_CONSEXPR_EXPR* expr_half_sqrx; /* 0.5 * x^2 */
static SCIP_CONSEXPR_EXPR* expr_sqrx; /* x^2 */
static SCIP_CONSEXPR_EXPR* expr_sum_fracpow; /* (x+y)^3.2 */
static SCIP_CONSEXPR_EXPR* expr_subprod_fracpow; /* x*(y*(x+y))^2.8 */

static
void setup(void)
{
   SCIP_CONSEXPR_EXPR* tmp;
   SCIP_CONSEXPR_EXPR* tmp2;
   SCIP_Real aux;

   conshdlr = NULL;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 2.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -1.0, 2.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create expressions for negative and positive values */
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr_posvalue, 1.3) );
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr_negvalue, -12.58) );

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_z, z) );

   /* create expression sum: x+y */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_sum, 1, &expr_x, NULL, 0.0) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr_sum, expr_y, 1.0) );

   /* create expression sum_fracpow: (x+y)^3.2 */
   //SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_sum_fracpow, 1, &expr_sum, &aux, 1.0) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &expr_sum_fracpow, expr_sum, 3.2) );

   /* create expression halfx: 0.5*x */
   /* SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_halfx, 1, &expr_x, NULL, 0.5) ); */ /*this is not simplified*/
   aux = 0.5;
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_halfx, 1, &expr_x, &aux, 0.0) );

   /* create expression sqrtx: sqrt(x) */
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &expr_sqrtx, expr_x, 0.5) );

   /* create expression sqrx: x^2 */
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &expr_sqrx, expr_x, 2.0) );

   /* create expression half_srqx: 0.5*x^2 */
   aux = 0.5;
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_half_sqrx, 1, &expr_sqrx, &aux, 0.0) );

   /* create expression prod: x*y*(x+y) */
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_prod, 1, &expr_x, 1.0) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_prod, expr_y) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_prod, expr_sum) );

   /* create expression subprod_fracpower: x*(y*(x+y))^2.8 */
   aux = 2.8;
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_subprod_fracpow, 1, &expr_x, 1.0) );  /* x */
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &tmp, 1, &expr_y, 1.0) );  /* y */
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, tmp, expr_sum) );  /* y -> y * (x+y) */
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &tmp2, tmp, aux) );  /* (y*(x+y))^aux */
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_subprod_fracpow, tmp2) );   /* x -> x * (y*(x+y))^aux  */

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &tmp) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &tmp2) );
}

static
void teardown(void)
{
   /* release */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_z) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_negvalue) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_posvalue) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sum) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_prod) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_subprod_fracpow) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sum_fracpow) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_halfx) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sqrx) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_half_sqrx) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sqrtx) );

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
   SCIP_CONSEXPR_EXPR* lin_expr1;
   SCIP_CONSEXPR_EXPR* lin_expr2;

   cr_expect( SCIPcompareConsExprExprs(expr_x, expr_x) == 0 );
   cr_expect( SCIPcompareConsExprExprs(expr_y, expr_y) == 0 );
   cr_expect( SCIPcompareConsExprExprs(expr_x, expr_y) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_y, expr_x) == 1 );

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &lin_expr1, 1, &expr_x, NULL, 2.0) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &lin_expr2, 1, &expr_x, NULL, 3.0) );
   cr_expect( SCIPcompareConsExprExprs(lin_expr1, lin_expr2) == -1 );
   cr_expect( SCIPcompareConsExprExprs(lin_expr2, lin_expr1) == 1 );
   cr_expect( SCIPcompareConsExprExprs(lin_expr1, lin_expr1) == 0 );

   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, lin_expr1, expr_y, 3.0) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, lin_expr2, expr_y, 3.0) );
   cr_expect( SCIPcompareConsExprExprs(lin_expr1, lin_expr2) == -1 );
   cr_expect( SCIPcompareConsExprExprs(lin_expr2, lin_expr1) == 1 );
   cr_expect( SCIPcompareConsExprExprs(lin_expr1, lin_expr1) == 0 );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &lin_expr1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &lin_expr2) );
}

Test(compare, values)
{
   cr_expect( SCIPcompareConsExprExprs(expr_negvalue, expr_negvalue) == 0 );
   cr_expect( SCIPcompareConsExprExprs(expr_posvalue, expr_negvalue) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_negvalue, expr_posvalue) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_posvalue, expr_x) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_posvalue, expr_y) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_posvalue, expr_sum) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_posvalue, expr_prod) == -1 );
}

Test(compare, exponents)
{
   cr_expect( SCIPcompareConsExprExprs(expr_halfx, expr_x) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_x, expr_halfx) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_sqrtx, expr_x) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_x, expr_sqrtx) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_sqrtx, expr_sqrx) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_sqrx, expr_sqrtx) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_x, expr_sqrx) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_sqrx, expr_x) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_prod, expr_sum_fracpow) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_sum_fracpow, expr_prod) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_sum_fracpow, expr_sum_fracpow) == 0 );
}

Test(compare, sums_and_products)
{
   cr_expect( SCIPcompareConsExprExprs(expr_sum, expr_prod) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_prod, expr_sum) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_prod, expr_prod) == 0 );
   cr_expect( SCIPcompareConsExprExprs(expr_sum, expr_sum) == 0 );
   cr_expect( SCIPcompareConsExprExprs(expr_x, expr_sum) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_y, expr_sum) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_x, expr_prod) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_y, expr_prod) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_sum, expr_x) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_sum, expr_y) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_prod, expr_x) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_prod, expr_y) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_prod, expr_subprod_fracpow) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_subprod_fracpow, expr_prod) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_subprod_fracpow, expr_subprod_fracpow) == 0 );
   cr_expect( SCIPcompareConsExprExprs(expr_subprod_fracpow, expr_z) == -1 );
   cr_expect( SCIPcompareConsExprExprs(expr_z, expr_subprod_fracpow) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_sqrx, expr_half_sqrx) == 1 );
   cr_expect( SCIPcompareConsExprExprs(expr_half_sqrx, expr_sqrx) == -1 );
}
