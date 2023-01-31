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

/**@file   parse.c
 * @brief  tests expression parsing
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/struct_expr.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
SCIP_RETCODE check_nuses(
   SCIP_EXPR*            rootexpr            /**< expression from which to check number of uses */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;

   assert(scip != NULL);
   assert(rootexpr != NULL);

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, FALSE) );

   for( expr = rootexpr; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      if( expr->nuses != 1 )
      {
         printf("following expression is captured too many times (%d, expected 1)\n", expr->nuses);
         SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
         assert(expr->nuses == 1);
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z(  ", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

static
void teardown(void)
{

   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(parse, .init = setup, .fini = teardown);

Test(parse, simple)
{
   SCIP_EXPR* expr_xy5;
   const char* input = "<x>[C] / <y>[I] *(-5)";

   /* create expression for product of -5, x, and y */
   cr_expect_eq(SCIPparseExpr(scip, &expr_xy5, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   /* print expression */
   SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
   SCIP_CALL( SCIPprintExpr(scip, expr_xy5, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* check that the expression is captured correctly */
   SCIP_CALL( check_nuses(expr_xy5) );

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_xy5) );
}

Test(parse, simple2)
{
   SCIP_EXPR* crazyexpr;
   const char* input = "-<x>[C] * <y>[I] ^(-1) + (<x>[C]+<y>[C])^2";

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &crazyexpr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   /* print expression */
   SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
   SCIP_CALL( SCIPprintExpr(scip, crazyexpr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &crazyexpr) );
}

Test(parse, signpower)
{
   SCIP_EXPR* crazyexpr;
   SCIP_EXPR* child;
   const char* input = "signpower(<x>^2   , 2.5)";

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &crazyexpr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   cr_assert(SCIPisExprSignpower(scip, crazyexpr));
   cr_assert_eq(SCIPgetExponentExprPow(crazyexpr), 2.5);

   child = SCIPexprGetChildren(crazyexpr)[0];
   cr_assert(SCIPisExprPower(scip, child));
   cr_assert_eq(SCIPgetExponentExprPow(child), 2.0);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &crazyexpr) );
}

/* disable undefined-behavior sanitizer because of intentional division by zero */
#if defined(__GNUC__) && __GNUC__ * 100 + __GNUC_MINOR__ * 10 >= 490 && !defined(__INTEL_COMPILER)
__attribute__((no_sanitize_undefined))
#endif
Test(parse, eval)
{
   SCIP_EXPR* crazyexpr;
   SCIP_SOL* crazysol;
   const char* input = "<x>*<y>^2/<x>^4 - 2*<x>*(3+5*<x>-2*<y>)*(<x>+<y>)^(-3.5)";
   const SCIP_Real vals[3][2] = { {1.0, 2.0}, {0.0, 0.0}, {42.0, 17.0} };
   int p;

#define CRAZYEVAL(x, y) ((x)*pow(y,2)/pow(x,4) - 2*(x)*(3+5*(x)-2*(y)) * pow((x)+(y), -3.5))

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &crazyexpr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   /* print expression */
   SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
   SCIP_CALL( SCIPprintExpr(scip, crazyexpr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* test expression by evaluating it */
   SCIP_CALL( SCIPcreateSol(scip, &crazysol, NULL) );

   for( p = 0; p < 3; ++p )
   {
      SCIP_Real expvalue = CRAZYEVAL(vals[p][0], vals[p][1]);

      SCIP_CALL( SCIPsetSolVal(scip, crazysol, x, vals[p][0]) );
      SCIP_CALL( SCIPsetSolVal(scip, crazysol, y, vals[p][1]) );

      SCIP_CALL( SCIPevalExpr(scip, crazyexpr, crazysol, 0) );
      SCIPinfoMessage(scip, NULL, "value for x=%g y=%g is %g, expected: %g\n", vals[p][0], vals[p][1], SCIPexprGetEvalValue(crazyexpr), expvalue);
      if( SCIPexprGetEvalValue(crazyexpr) == SCIP_INVALID )
         cr_expect(!SCIPisFinite(expvalue));
      else
         cr_expect(SCIPisEQ(scip, SCIPexprGetEvalValue(crazyexpr), expvalue));
   }

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &crazyexpr) );

   /* release solution */
   SCIP_CALL( SCIPfreeSol(scip, &crazysol) );
}

Test(parse, unusual_var_name)
{
   SCIP_EXPR* expr;
   const char* input = "(<x> - <y>) /   <z(  >^2";

   /* parse */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   /* print expression */
   SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* check that the expression is capture correctly */
   SCIP_CALL( check_nuses(expr) );

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(parse, invalid_expressions)
{
   SCIP_EXPR* e;

   SCIPmessageSetErrorPrinting(NULL, NULL);

   /* there is no variable with name "xx" */
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"<xx>", NULL, NULL, NULL), SCIP_READERROR);
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"5/<donothave> ", NULL, NULL, NULL), SCIP_READERROR);
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"<x> +-*5 ", NULL, NULL, NULL), SCIP_READERROR);
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"<x> / (<y>-5 ", NULL, NULL, NULL), SCIP_READERROR);
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"donothave(<x>) ", NULL, NULL, NULL), SCIP_READERROR);
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"donothave(<x> ", NULL, NULL, NULL), SCIP_READERROR);
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"val(1) ", NULL, NULL, NULL), SCIP_READERROR);

#ifdef FAILING_TESTS
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"<x>+2<y> ", NULL, NULL, NULL), SCIP_READERROR);
#endif

   SCIPmessageSetErrorPrintingDefault();
}

Test(parse, misc)
{
   SCIP_EXPR* e;

   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"-5+3*<x>", NULL, NULL, NULL), SCIP_OKAY);
   cr_expect(SCIPisExprSum(scip, e));
   cr_expect_eq(SCIPgetConstantExprSum(e), -5.0);
   cr_expect_eq(SCIPgetCoefsExprSum(e)[0], 3.0);
   cr_expect_eq(SCIPexprGetNChildren(e), 1);
   cr_expect(SCIPisExprVar(scip, SCIPexprGetChildren(e)[0]));
   cr_expect_eq(SCIPgetVarExprVar(SCIPexprGetChildren(e)[0]), x);
   SCIP_CALL( SCIPreleaseExpr(scip, &e) );

   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"<x>", NULL, NULL, NULL), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseExpr(scip, &e) );
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"<x> +5*<y>", NULL, NULL, NULL), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseExpr(scip, &e) );
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"<x> + <x>*5*<y>", NULL, NULL, NULL), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseExpr(scip, &e) );
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"<x> +5", NULL, NULL, NULL), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseExpr(scip, &e) );
   cr_expect_eq(SCIPparseExpr(scip, &e, (char*)"<x> +5 + <y>", NULL, NULL, NULL), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseExpr(scip, &e) );
}
