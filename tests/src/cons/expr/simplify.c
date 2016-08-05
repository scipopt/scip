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

/**@file   cmain.c
 * @brief  unit test for checking cons_expr
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "scip/scip.h"
#include "scip/cons_expr.h"

#include "include/scip_test.h"


/* specify parameters */
ParameterizedTestParameters(simplify /* test suite */, simplify_test /* test name */)
{
   /* needs! to be static */
   static const char* expressions[] =
   {
      "1+2*2+3",
      "1*2*2*3",
      "2*<x>",
      "<y> + <x> + 1",
      "<y> + <x> + 1 +2*(<y> + <x> + 1)",
      "<x>*<y> + 1 +0.5*(0+2*<x>*<y>)",
      "<x>*<y> + 0.5*(0+2*<x>*<y>)",
      "<x>*<y> + 0.5*(0+2*<x>*<y>) -<x>*0.5*<y>*2",
      "0+0",
      "-<x>+2*<y>-<y>-<y>",
      "-<x>+2*<y>+2*(0+0.5*<x>-<y>)",
      "<x>*<x>",
      "(2*<x>)*(2*<x>)",
      "<x>*<x>^2",
      "(<x>^0.5)^2",
      "(<y>^2)^2",
      "1*2*(<x>+<y>)*<x>*4*0*5",
      "(<x>^0.25)^2*(<x>^0.25)^2",
      "(<x>)^0.25*(<x>)^0.25*(<x>)^0.25*(<x>)^0.25",
      "(<x>^0.2)^1.25*(<x>^0.2)^1.25*(<x>^0.2)^1.25*(<x>^0.2)^1.25",
      "<x>^0.5 * (<x>^0.8)^(-0.625)",
      "(<x>^0.5*<y>^0.5)^0.5*(<x>^0.5*<y>^0.5)^0.5 * <y>",
      "2+3*<x>^2-<y>*<x>*4*(<x>+<y>)",
      "2*<x>*<y>",
      "<x>^10 + <x>^9.5 + <x>^9 + <x>^8.5 + <x>^8 + <x>^7.5",
      "<x>/<x>",
      "(2*<x>)^2",
      "(2*<x>)*(2*<x>) - 4 * <x>^2"
   };

   /* type of the parameter; the parameter; number of parameters */
   return cr_make_param_array(const char*, expressions, sizeof(expressions)/sizeof(const char *));
}

struct test_data
{
   SCIP* scip;
   SCIP_SOL* sol1;
   SCIP_SOL* sol2;
};


static struct test_data tdata;

static
void setup(void)
{
   SCIP_VAR* x;
   SCIP_VAR* y;

   SCIP_CALL( SCIPcreate(&tdata.scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(tdata.scip) );


   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(tdata.scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(tdata.scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(tdata.scip, &y, "y", 0.0, 2.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(tdata.scip, x) );
   SCIP_CALL( SCIPaddVar(tdata.scip, y) );

   /* create solutions: this should be on test; but to show how to define "test data" */
   SCIP_CALL( SCIPcreateSol(tdata.scip, &tdata.sol1, NULL) );
   SCIP_CALL( SCIPcreateSol(tdata.scip, &tdata.sol2, NULL) );

   SCIP_CALL( SCIPsetSolVal(tdata.scip, tdata.sol1, x, 1.2) );
   SCIP_CALL( SCIPsetSolVal(tdata.scip, tdata.sol1, y, -2.4) );
   SCIP_CALL( SCIPsetSolVal(tdata.scip, tdata.sol2, x, 0.7463) );
   SCIP_CALL( SCIPsetSolVal(tdata.scip, tdata.sol2, y, 12.037) );

   /* release variables */
   SCIP_CALL( SCIPreleaseVar(tdata.scip, &y) );
   SCIP_CALL( SCIPreleaseVar(tdata.scip, &x) );
}

static
void teardown(void)
{
   /* release everything */
   SCIP_CALL( SCIPfreeSol(tdata.scip, &tdata.sol2) );
   SCIP_CALL( SCIPfreeSol(tdata.scip, &tdata.sol1) );
   SCIP_CALL( SCIPfree(&tdata.scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}


/* actual test; we get one parameter as argument
 * NOTE: only one parameter as argument, so build a struct if you need several parameters
 */
ParameterizedTest(const char** expression, simplify, simplify_test, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real values[2];
   SCIP* scip = tdata.scip;
   SCIP_SOL* sol1 = tdata.sol1;
   SCIP_SOL* sol2 = tdata.sol2;

   //printf("parametrized test: processing expression %s\n", *expression);

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(tdata.scip, "expr");
   cr_assert_not_null(conshdlr, "No expression constraint handler!!!");

   /* parse expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, *expression, NULL, &expr) );

   /* evaluate */
   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol1, 0) );
   values[0] = SCIPgetConsExprExprValue(expr);
   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol2, 0) );
   values[1] = SCIPgetConsExprExprValue(expr);

   /* simplify */
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, &expr) );

   /* test its evaluate the same; expect because we want to release the expression even if this fails */
   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol1, 0) );
   cr_expect_float_eq(values[0], SCIPgetConsExprExprValue(expr), SCIPfeastol(scip));

   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol2, 0) );
   cr_assert_float_eq(values[1], SCIPgetConsExprExprValue(expr), SCIPfeastol(scip));

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
