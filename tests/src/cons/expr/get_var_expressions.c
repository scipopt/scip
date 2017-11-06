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

/**@file   get_var_expressions.c
 * @brief  tests for methods to get expressions of variables
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* w;
static SCIP_CONSEXPR_EXPR* expr;

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

   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -3.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -2.0, 2.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -3.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -3.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3", NULL, &expr) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

TestSuite(getvars, .init = setup, .fini = teardown);

Test(getvars, expression_not_containing_all_vars)
{
   SCIP_CONSEXPR_EXPR* varexprs[3];
   int nvarexprs;
   int i;
   
   /* note that this captures the variable expressions */
   SCIP_CALL( getVarExprs(scip, expr, varexprs, &nvarexprs) );
   cr_assert_eq(nvarexprs, 3);

   for( i = 0; i < nvarexprs; ++i )
   {
      cr_assert_not_null(varexprs[i]);
      cr_assert_str_eq(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(varexprs[i])), "var");

      /* release variable expression */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[i]) );
   }
}

Test(getvars, expression_containing_all_vars)
{
   SCIP_CONSEXPR_EXPR* varexprs[4];
   SCIP_CONSEXPR_EXPR* wexpr;
   SCIP_CONSEXPR_EXPR* sumexpr;
   int nvarexprs;
   int i;
   
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, 0, NULL, NULL, 0) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &wexpr, w) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, sumexpr, wexpr, 1.0) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, sumexpr, expr, 1.0) );

   /* note that this captures the variable expressions */
   SCIP_CALL( getVarExprs(scip, sumexpr, varexprs, &nvarexprs) );
   cr_assert_eq(nvarexprs, 4, "Expecting 4, got %d\n", nvarexprs);

   for( i = 0; i < nvarexprs; ++i )
   {
      cr_assert_not_null(varexprs[i]);
      cr_assert_str_eq(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(varexprs[i])), "var");

      /* release variable expression */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[i]) );
   }

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &wexpr) );
}
