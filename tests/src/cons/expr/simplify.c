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
#include "scip/nodesel_bfs.h"

#include "include/scip_test.h"

struct expr_type
{
   const char* expr;
   const char* type;
};

/* specify parameters */
ParameterizedTestParameters(simplify /* test suite */, simplify_test /* test name */)
{
   /* needs! to be static */
   static const struct expr_type expressions[] =
   {
      {"1+2*2+3", "val"},
      {"1*2*2*3", "val"},
      {"2*<x>", "sum"},
      {"<y> + <x> + 1", "sum"},
      {"<y> + <x> + 1 +2*(<y> + <x> + 1)", "sum"},
      {"<x>*<y> + 1 +0.5*(0+2*<x>*<y>)", "sum"},
      {"<x>*<y> + 0.5*(0+2*<x>*<y>)", "sum"},
      {"<x>*<y> + 0.5*(0+2*<x>*<y>) -<x>*0.5*<y>*2", "prod"},
      {"0+0", "val"},
      {"-<x>+2*<y>-<y>-<y>", "sum"},
      {"-<x>+2*<y>+2*(0+0.5*<x>-<y>)", "val"},
      {"<x>*<x>", "prod"},
      {"(2*<x>)*(2*<x>)", "sum"},
      {"<x>*<x>^2", "prod"},
      {"(<x>^0.5)^2", "var"},
      {"(<y>^2)^2", "prod"},
      {"1*2*(<x>+<y>)*<x>*4*0*5", "val"},
      {"(<x>^0.25)^2*(<x>^0.25)^2", "var"},
      {"(<x>)^0.25*(<x>)^0.25*(<x>)^0.25*(<x>)^0.25", "var"},
      {"(<x>^0.2)^1.25*(<x>^0.2)^1.25*(<x>^0.2)^1.25*(<x>^0.2)^1.25", "var"},
      {"<x>^0.5 * (<x>^0.8)^(-0.625)", "prod"}, // this should probably be simplified further
      {"(<x>^0.5*<y>^0.5)^0.5*(<x>^0.5*<y>^0.5)^0.5 * <y>", "prod"},
      {"2+3*<x>^2-<y>*<x>*4*(<x>+<y>)", "sum"},
      {"2*<x>*<y>", "sum"},
      {"<x>^10 + <x>^9.5 + <x>^9 + <x>^8.5 + <x>^8 + <x>^7.5", "sum"},
      {"<x>/<x>", "val"},
      {"(2*<x>)^2", "sum"},
      {"(<x> + <y>)^2", "sum"},
      {"(<x> + <y> + 2)^2", "sum"},
      {"(<x> + <y>)^2 - <x>^2 - 2*<x>*<y>", "prod"},
      {"(<x> + <y>)^2 - <x>^2 - 2*<x>*<y> - <y>^2", "val"},
      {"(<x> + <y> + <fixvar>)^2 - <x>^2 - 2*<x>*<y> - <y>^2 - <fixvar>^2 -2*<x>*<fixvar> - 2*<fixvar>*<y>", "val"},
      {"(1 + (<x>*<y>^2)^0.5)^2 - 1 - 2*(<x>*<y>^2)^0.5 -<y>^2*<x>", "val"},
      {"((<x>-<y>)^2 + (<x>+<y>)^2)^2/4 - <x>^4 - 2*(<x>*<y>)^2 - <y>^4", "val"},
      {"(2*<x>)*(2*<x>) - 4 * <x>^2", "val"},
      {"abs(-3.0)", "val"},
      {"log(exp(1.0))", "val"},
      {"exp(-3.0)", "val"},
      {"log(abs(-3.0))", "val"},
      {"log(3.0)", "val"}
      //{"<fixvar>", "val"}
      //{"<fixvar>^2", "val"}
   };

   /* type of the parameter; the parameter; number of parameters */
   return cr_make_param_array(const struct expr_type, expressions, sizeof(expressions)/sizeof(struct expr_type));
}

static SCIP_SOL* sol1;
static SCIP_SOL* sol2;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* fixvar;
static SCIP_VAR* agg1;
static SCIP_VAR* agg2;
static SCIP_VAR* multiagg1;
static SCIP_VAR* multiagg2;

/* helper method */
static
void parseSimplifyCheck(SCIP* scip, const char* input, const char* type)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real values[2];

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr, "No expression constraint handler!!!");

   /* parse expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, input, NULL, &expr) );

   /* evaluate */
   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol1, 0) );
   values[0] = SCIPgetConsExprExprValue(expr);
   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol2, 0) );
   values[1] = SCIPgetConsExprExprValue(expr);

   /* simplify */
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, expr, &simplified) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   expr = simplified;

   /* check type of simplified expression */
   cr_expect_str_eq(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), type);

   /* test it evaluates to the same; expect because we want to release the expression even if this fails */
   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol1, 0) );
   cr_expect_float_eq(values[0], SCIPgetConsExprExprValue(expr), SCIPfeastol(scip));

   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol2, 0) );
   cr_expect_float_eq(values[1], SCIPgetConsExprExprValue(expr), SCIPfeastol(scip));

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExec)
{  /*lint --e{715}*/
   SCIP_Bool infeasible;
   SCIP_Bool success;

   /* TEST multi-aggregated variable */
   {
      SCIP_Real scalars[2];
      SCIP_VAR* vars[2];
      SCIP_Real constant;

      /* multi aggregate: multiagg1 = 2.1 * <x>  - 1.5 * <y> - 0.3
       * -> multiagg1 = 2.1 * [1.2, 0.7463] - 1.5 * [-2.4, 12.037] - 0.3 = [5.82, -16.78827] */
      scalars[0] = 2.1; scalars[1] = -1.5; constant = -0.3;
      vars[0] = x; vars[1] = y;
      SCIP_CALL( SCIPmultiaggregateVar(scip, multiagg1, 2, vars, scalars, constant, &infeasible, &success) );
      cr_assert(!infeasible && success);

      /* set sol value; this wouldn't work if multiagg1 was created here, since sol is original, whereas multiagg1 would be trans */
      SCIP_CALL( SCIPsetSolVal(scip, sol1, multiagg1, 5.82) );
      SCIP_CALL( SCIPsetSolVal(scip, sol2, multiagg1, -16.78827) );

      /* simplify */
      parseSimplifyCheck(scip, "<t_multiagg1>", "sum");
   }

   /* TEST fixed variable */
   {
      printf("fixing var status %d\n", SCIPvarGetStatus(fixvar));
      SCIP_CALL( SCIPfixVar(scip, fixvar, 1.5, &infeasible, &success) );
      printf("%s is active? %d \n", SCIPvarGetName(fixvar), SCIPvarIsActive(fixvar));
      cr_assert(!infeasible && success);

      /* simplify */
      parseSimplifyCheck(scip, "<t_fixvar>", "val");
   }

   /* TEST multi-aggregated variable to fixed variabe */
   {
      SCIP_Real scalars[1];
      SCIP_Real constant;

      /* multi aggregate: multiagg2 = 2.1 * <fixvar> - 0.3 -> multiagg2 = 2.1 * 1.5 - 0.3 = 2.85 */
      scalars[0] = 2.1; constant = -0.3;
      printf("status pre multiagg %d\n", SCIPvarGetStatus(multiagg2));
      SCIP_CALL( SCIPmultiaggregateVar(scip, multiagg2, 1, &fixvar, scalars, constant, &infeasible, &success) );
      printf("status post multiagg %d\n", SCIPvarGetStatus(multiagg2));
      cr_assert(!infeasible && success);

      /* set sol value */
      SCIP_CALL( SCIPsetSolVal(scip, sol1, multiagg2, 2.85) );
      SCIP_CALL( SCIPsetSolVal(scip, sol2, multiagg2, 2.85) );

      /* simplify */
      parseSimplifyCheck(scip, "<t_multiagg2>", "val");
   }

   /* TEST aggregated variable */
   {
      SCIP_Bool redundant;
      /* aggregate: agg1 = 5.7/7.5 + 3.3/7.5 x = 0.76 + 0.44 x -> agg1 = 0.76 + 0.44 *[1.2, 0.7463] = [1.288, 1.088372] */
      SCIP_CALL( SCIPaggregateVars(scip, agg1, x, 7.5, -3.3, 5.7, &infeasible, &redundant, &success) );
      cr_assert(!infeasible && success);
      cr_assert(!SCIPvarIsActive(SCIPvarGetTransVar(agg1)));

      /* set sol value */
      SCIP_CALL( SCIPsetSolVal(scip, sol1, agg1, 1.288) );
      SCIP_CALL( SCIPsetSolVal(scip, sol2, agg1, 1.088372) );

      /* simplify */
      parseSimplifyCheck(scip, "<t_agg1>", "sum");
   }

   /* TEST aggregated variable to single var*/
   {
      SCIP_Bool redundant;
      /* aggregate: agg2 = x -> agg2 = [1.288, 1.088372] */
      SCIP_CALL( SCIPaggregateVars(scip, agg2, x, 1.0, -1.0, 0.0, &infeasible, &redundant, &success) );
      cr_assert(!infeasible && success);
      cr_assert(!SCIPvarIsActive(SCIPvarGetTransVar(agg2)));

      /* set sol value */
      SCIP_CALL( SCIPsetSolVal(scip, sol1, agg2, 1.2) );
      SCIP_CALL( SCIPsetSolVal(scip, sol2, agg2, 0.7463) );

      /* simplify */
      parseSimplifyCheck(scip, "<t_agg2>", "var");
   }
   return SCIP_OKAY;
}

static SCIP* scip;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* include our presol which is needed for some tests */
   SCIP_CALL( SCIPincludePresolBasic(scip, NULL, "presol", "presol", 1, 1, SCIP_PRESOLTIMING_FAST,
         presolExec,
         NULL) );

   /* needed for presolve */
   SCIP_CALL( SCIPincludeNodeselBfs(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 2.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &fixvar, "fixvar", 1.5, 1.5, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &agg1, "agg1", -20.0, 23.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &agg2, "agg2", -20.0, 23.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &multiagg1, "multiagg1", -10.3, 10.0, 1.12, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &multiagg2, "multiagg2", -10.3, 10.0, 1.12, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, fixvar) );
   SCIP_CALL( SCIPaddVar(scip, agg1) );
   SCIP_CALL( SCIPaddVar(scip, agg2) );
   SCIP_CALL( SCIPaddVar(scip, multiagg1) );
   SCIP_CALL( SCIPaddVar(scip, multiagg2) );

   /* create solutions */
   SCIP_CALL( SCIPcreateSol(scip, &sol1, NULL) );
   SCIP_CALL( SCIPcreateSol(scip, &sol2, NULL) );

   SCIP_CALL( SCIPsetSolVal(scip, sol1, x, 1.2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol1, y, -2.4) );
   SCIP_CALL( SCIPsetSolVal(scip, sol1, fixvar, 1.5) );

   SCIP_CALL( SCIPsetSolVal(scip, sol2, x, 0.7463) );
   SCIP_CALL( SCIPsetSolVal(scip, sol2, y, 12.037) );
   SCIP_CALL( SCIPsetSolVal(scip, sol2, fixvar, 1.5) );
}

static
void teardown(void)
{
   /* release everything */
   SCIP_CALL( SCIPreleaseVar(scip, &multiagg2) );
   SCIP_CALL( SCIPreleaseVar(scip, &multiagg1) );
   SCIP_CALL( SCIPreleaseVar(scip, &agg2) );
   SCIP_CALL( SCIPreleaseVar(scip, &agg1) );
   SCIP_CALL( SCIPreleaseVar(scip, &fixvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIP_CALL( SCIPfreeSol(scip, &sol2) );
   SCIP_CALL( SCIPfreeSol(scip, &sol1) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

/***** TEST SUITE: all tests of the form Test(simplify, xxx) belong to the same suite and share the setup and teardown *****/
TestSuite(simplify, .init = setup, .fini = teardown);

/* actual test; we get one parameter as argument */
ParameterizedTest(const struct expr_type* expression, simplify, simplify_test)
{
   //printf("received %s and %s\n", expression->expr, expression->type);
   parseSimplifyCheck(scip, expression->expr, expression->type);
}

/* non-parametrized test, which calls presolve */
Test(simplify, remove_fix_variables)
{
   SCIP_CALL( SCIPpresolve(scip) );
   SCIP_CALL( SCIPfreeTransform(scip) ); /* why do I need to call this one before freeing the sols? */
}

///* parameterized test don't work with --single :/ */
//Test(simplify, debug)
//{
//   parseSimplifyCheck(scip, "(<x> + <y> + <fixvar>)^2 - <x>^2 - 2*<x>*<y> - <y>^2 - <fixvar>^2 -2*<x>*<fixvar> - 2*<fixvar>*<y>", "val");
//}
