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

/**@file   getvarexprs.c
 * @brief  tests for methods to get expressions of variables
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr.h"
#include "scip/scipdefplugins.h"
#include "scip/expr_sum.h"
#include "scip/expr_var.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_EXPR* expr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* w;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

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

   SCIP_CALL( SCIPparseExpr(scip, &expr, "1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3", NULL, NULL, NULL) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
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
   SCIP_EXPR* varexprs[3];
   int nvarexprs;
   int i;

   /* note that this captures the variable expressions */
   SCIP_CALL( SCIPgetExprVarExprs(scip, expr, varexprs, &nvarexprs) );
   cr_assert_eq(nvarexprs, 3);

   for( i = 0; i < nvarexprs; ++i )
   {
      cr_assert_not_null(varexprs[i]);
      cr_assert(SCIPisExprVar(scip, varexprs[i]));

      /* release variable expression */
      SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[i]) );
   }
}

Test(getvars, expression_containing_all_vars)
{
   SCIP_EXPR* varexprs[4];
   SCIP_EXPR* wexpr;
   SCIP_EXPR* sumexpr;
   int nvarexprs;
   int i;

   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 0, NULL, NULL, 0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &wexpr, w, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, sumexpr, wexpr, 1.0) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, sumexpr, expr, 1.0) );

   /* note that this captures the variable expressions */
   SCIP_CALL( SCIPgetExprVarExprs(scip, sumexpr, varexprs, &nvarexprs) );
   cr_assert_eq(nvarexprs, 4, "Expecting 4, got %d\n", nvarexprs);

   for( i = 0; i < nvarexprs; ++i )
   {
      cr_assert_not_null(varexprs[i]);
      cr_assert(SCIPisExprVar(scip, varexprs[i]));

      /* release variable expression */
      SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[i]) );
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &wexpr) );
}
