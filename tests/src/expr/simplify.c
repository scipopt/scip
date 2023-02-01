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

/**@file   simplify.c
 * @brief  unit test for expression simplification
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

struct expr_type
{
   const char expr[1000];
   const char type[1000];
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
      {"<x>*<x>", "pow"},
      {"(2*<x>)*(2*<x>)", "sum"},
      {"<x>*<x>^2", "pow"},
      {"(<x>^0.5)^2", "pow"},  /* not simplified because of implicit x >= 0 constraint */
      {"(<x>^2)^0.5", "abs"},
      {"(<x>^3)^0.5", "pow"},  /* should be |x|^1.5 */
      {"(<x>^4)^0.5", "pow"},
      {"(<x>^2)^1.5", "pow"},  /* should be |x|^3 */
      {"(<y>^2)^2", "pow"},
      {"1*2*(<x>+<y>)*<x>*4*0*5", "val"},
      {"(<x>^0.25)^2*(<x>^0.25)^2", "var"},
      {"(<x>)^0.25*(<x>)^0.25*(<x>)^0.25*(<x>)^0.25", "var"},
      {"(<x>^0.2)^1.25*(<x>^0.2)^1.25*(<x>^0.2)^1.25*(<x>^0.2)^1.25", "var"},
      {"<x>^0.5 * (<x>^0.8)^(-0.625)", "val"},
      {"<x>^0.5 * (<x>^2)^(-0.5)", "pow"},
      {"(<x>^0.5*<y>^0.5)^0.5*(<x>^0.5*<y>^0.5)^0.5 * <y>", "prod"},
      {"2+3*<x>^2-<y>*<x>*4*(<x>+<y>)", "sum"},
      {"2*<x>*<y>", "sum"},
      {"<x>^10 + <x>^9.5 + <x>^9 + <x>^8.5 + <x>^8 + <x>^7.5", "sum"},
      {"<x>/<x>", "val"},
      {"(2*<x>)^2", "sum"},
      {"(<x> + <y>)^2", "sum"},
      {"(<x> + <y> + 2)^2", "sum"},
      {"(<x> + <y>)*(<x> + 1)^2 - <x>^3 - <x>^2*<y> - 2*<x>^2 - 2*<y>*<x> -<y> -<x>", "val"},
      {"(<x> + <y>)^2 - <x>^2 - 2*<x>*<y>", "pow"},
      {"-<x>^2 + (<x> + <y>)^2 - 2*<x>*<y> - <y>^2", "val"}, // order is important to test the internal algorithms
      {"<x>^2 * (1.0 / (<x>^2 * <y>)) * <y>", "val"}, // order is important to test the internal algorithms
      {"(<x> + <y> + <fixvar>)^2 - <x>^2 - 2*<x>*<y> - <y>^2 - <fixvar>^2 -2*<x>*<fixvar> - 2*<fixvar>*<y>", "val"},
      {"(1 + (<x>*<y>^2)^0.5)^2 - 1 - 2*(<x>*<y>^2)^0.5 -<y>^2*<x>", "val"},
      {"((<x>-<y>)^2 + (<x>+<y>)^2)^2/4 - <x>^4 - 2*(<x>*<y>)^2 - <y>^4", "val"},
      {"(2*<x>)*(2*<x>) - 4 * <x>^2", "val"},
      {"abs(-3.0)", "val"},
      {"log(exp(1.0))", "val"},
      {"exp(-3.0)", "val"},
      {"log(abs(-3.0))", "val"},
      {"log(3.0)", "val"},
      {"((<x>^0.5)^0.5 + <y> + 1)*(<x> + (<x>^0.5)^0.5 + 2)", "sum"},
      {"((<x>+<y>)^0.5 - <y>)*((<x> + <y>)^0.5 + 1)", "sum"},
      {"(<x>*<y>^0.5 + 1)*(<y>^1.5/<x> - 1)", "sum"},
      {"(2*<x>+1)*(<x>+1)", "sum"},
      {"(2*<x>+1)*(<x> + 1)", "sum"},
      {"((<x>+<y>)^0.5 - <y>)*((<x> + <y>)^0.5 + 1)", "sum"},
      {"((-2*<x>)^0.5+1)*((-2*<x>)^0.5 + 1)", "sum"},
      {"((2*<x>)^0.5+1)*((2*<x>)^0.5 + 1)", "sum"},
      {"((2*<x>)^0.5+1)*((1*<x>)^0.5 + 1)", "sum"},
      {"<x>*(<x>+1)", "sum"},
      {"2*<x>*(<x>+1)", "sum"},
      {"(<x>^0.5)^0.5*((<x>^0.5)^0.5+2) - 2*(<x>^0.5)^0.5 - <x>^0.5", "val"},
      {"(25.0 * <x>^2)^0.5", "sum"},
      {"exp(<x>)*exp(<y>)", "exp"},
      {"exp(<x>)*exp(-<x>)", "val"},
      {"exp(<x>^2)*<x>*exp(-<x>^2)", "var"},
      {"<x>*exp(<x>^2)*exp(-<x>^2)", "var"},
      {"<x>*exp(<x>^2)*exp(-<x>^2)*<x>", "pow"},
      {"2+exp(<x>*<y>)*exp(-<y>*<x>)", "val"},
      {"exp(<x>)^2", "exp"},
      {"2+2*<x>*exp(<x>*<y>)-2", "prod"},
      {"2+2*<x>*cos(<x>*<y>)-2", "sum"},
      {"2+(1+1)*<x>*exp(<x>*<y>)*2-2", "prod"},
      {"10.0*exp(<x>)", "exp"},
      {"-10.0*exp(<x>)", "sum"}
      /*TODO,
      {"<x>*abs(<x>)", "pow"}
      {"<x>*abs(<x>)^0.875", "pow"}*/
      /*{"<fixvar>", "val"}*/
      /*{"<fixvar>^2", "val"}*/
   };

   /* alloc memory */
  struct expr_type* gexpr = (struct expr_type*)cr_malloc(sizeof(expressions));

   for( unsigned int i = 0; i < sizeof(expressions)/sizeof(struct expr_type); ++i )
   {
      strcpy((char *)gexpr[i].type, (char *)expressions[i].type);
      strcpy((char *)gexpr[i].expr, (char *)expressions[i].expr);
   }

   /* type of the parameter; the parameter; number of parameters */
   return cr_make_param_array(const struct expr_type, gexpr, sizeof(expressions)/sizeof(struct expr_type));
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
void parseSimplifyCheck(SCIP* scip, const char* input, const char* type, SCIP_EXPR** simplifiedexpr)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_EXPR* simplified_again;
   SCIP_Real values[2];
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   /* parse expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, input, NULL, NULL, NULL) );

   /* evaluate */
   SCIP_CALL( SCIPevalExpr(scip, expr, sol1, 0) );
   values[0] = SCIPexprGetEvalValue(expr);
   SCIP_CALL( SCIPevalExpr(scip, expr, sol2, 0) );
   values[1] = SCIPexprGetEvalValue(expr);

#if 0
   SCIP_CALL( SCIPshowExpr(scip, expr) );
   fprintf(stderr, " simplifying!\n");
#else
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, " simplifying!\n");
#endif
   /* simplify */
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   cr_assert_not(infeasible);

#if 0
   fprintf(stderr,"done simplifying!\n");
   printf("original\n");
   SCIP_CALL( SCIPprintExpr(scip, expr, 0) );
   SCIPinfoMessage(scip, 0, "\n");
   printf("simplified\n");
   SCIP_CALL( SCIPprintExpr(scip, simplified, 0) );
   SCIPinfoMessage(scip, 0, "\n");
   //SCIP_CALL( SCIPshowExpr(scip, simplified) );
   printf("\n");
#endif

   /* check type of simplified expression */
   cr_expect_str_eq(SCIPexprhdlrGetName(SCIPexprGetHdlr(simplified)), type);

   /* test it evaluates to the same; expect because we want to release the expression even if this fails */
   SCIP_CALL( SCIPevalExpr(scip, simplified, sol1, 0) );
   cr_expect_float_eq(values[0], SCIPexprGetEvalValue(simplified), SCIPfeastol(scip), "expecting %f got %f",
         values[0], SCIPexprGetEvalValue(simplified));

   SCIP_CALL( SCIPevalExpr(scip, simplified, sol2, 0) );
   cr_expect_float_eq(values[1], SCIPexprGetEvalValue(simplified), SCIPfeastol(scip), "expecting %f got %f",
         values[1], SCIPexprGetEvalValue(simplified));

   /* test that the same expression is obtained when simplifying again */
   /*printf("~~~~~~~~~~~~~ simplify again ~~~~~~~~~~~~~~~~~~~~ \n");*/
   SCIP_CALL( SCIPsimplifyExpr(scip, simplified, &simplified_again, &changed, &infeasible, NULL, NULL) );
   cr_expect_eq(SCIPcompareExpr(scip, simplified, simplified_again), 0);
   cr_expect_not(changed);
   cr_assert_not(infeasible);


   /* release expressions */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &simplified_again) );

   if( simplifiedexpr != NULL )
      *simplifiedexpr = simplified;
   else
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );
   }

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
      parseSimplifyCheck(scip, "<t_multiagg1>", "sum", NULL);
   }

   /* TEST fixed variable */
   {
      printf("fixing var status %d\n", SCIPvarGetStatus(fixvar));
      SCIP_CALL( SCIPfixVar(scip, fixvar, 1.5, &infeasible, &success) );
      printf("%s is active? %d \n", SCIPvarGetName(fixvar), SCIPvarIsActive(fixvar));
      cr_assert(!infeasible && success);

      /* simplify */
      parseSimplifyCheck(scip, "<t_fixvar>", "val", NULL);
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
      parseSimplifyCheck(scip, "<t_multiagg2>", "val", NULL);
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
      parseSimplifyCheck(scip, "<t_agg1>", "sum", NULL);
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
      parseSimplifyCheck(scip, "<t_agg2>", "var", NULL);
   }
   return SCIP_OKAY;
}

static SCIP* scip;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include our presol which is needed for some tests */
   SCIP_CALL( SCIPincludePresolBasic(scip, NULL, "presol", "presol", 1, 1, SCIP_PRESOLTIMING_FAST, presolExec, NULL) );

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
   /*fprintf(stderr,"received %s and %s\n", expression->expr, expression->type);*/
   parseSimplifyCheck(scip, expression->expr, expression->type, NULL);
}

/* to debug parameterized test, since it doesn't work with --single :/ */
Test(simplify, debug)
{
   parseSimplifyCheck(scip, "-<x>+2*<y>+2*(0+0.5*<x>-<y>)", "val", NULL);
}


/* non-parametrized test, which calls presolve: tests aggregation */
Test(simplify, remove_fix_variables)
{
   SCIP_CALL( SCIPpresolve(scip) );
   SCIP_CALL( SCIPfreeTransform(scip) ); /* why do I need to call this one before freeing the sols? */
}

/* further simplification tests */
Test(simplify, more_simplification_tests)
{
   SCIP_EXPR* simplified = NULL;

   parseSimplifyCheck(scip, "<x>*exp(<x>)*exp(<x>)", "prod", &simplified);
   cr_assert_not_null(simplified);

   cr_expect_eq(SCIPexprGetNChildren(simplified), 2, "got %d", SCIPexprGetNChildren(simplified));
   SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );

   parseSimplifyCheck(scip, "<x>*exp(<x>)*exp(<y>)", "prod", &simplified);
   cr_assert_not_null(simplified);

   cr_expect_eq(SCIPexprGetNChildren(simplified), 2, "got %d", SCIPexprGetNChildren(simplified));
   SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );

   parseSimplifyCheck(scip, "<x>*exp(-<x>)*exp(<x>+<y>)", "prod", &simplified);
   cr_assert_not_null(simplified);

   cr_expect_eq(SCIPexprGetNChildren(simplified), 2, "got %d", SCIPexprGetNChildren(simplified));
   SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );

   parseSimplifyCheck(scip, "(exp(<x>/2.0))^2.0", "exp", &simplified);
   cr_assert_not_null(simplified);
   cr_expect(SCIPisExprVar(scip, SCIPexprGetChildren(simplified)[0]));
   SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );
}
