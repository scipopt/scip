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

/**@file   monotonicity.c
 * @brief  tests monotonicity expression handler callbacks
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */  /*lint !750 */

#include <string.h>
#include <math.h>

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_EXPR* expr;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   expr = NULL;
}

static
void teardown(void)
{
   if( expr != NULL )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }

   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/** helper function to create an expression */
static
SCIP_RETCODE createExpr(
   const char*           input,              /**< input creating an expression */
   const char*           exprhdlrname        /**< target expression handler name */
   )
{
   SCIP_EXPR* origexpr;
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   /* release previous expression */
   if( expr != NULL )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }

   /* create and print expression */
   cr_expect_eq(SCIPparseExpr(scip, &origexpr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   /* simplify expression */
   SCIP_CALL( SCIPsimplifyExpr(scip, origexpr, &expr, &changed, &infeasible, NULL, NULL) );
   SCIP_CALL( SCIPreleaseExpr(scip, &origexpr) );

   /* check name of the corresponding expression handler */
   exprhdlr = SCIPexprGetHdlr(expr);
   cr_assert(exprhdlr != NULL);
   cr_expect(strcmp(SCIPexprhdlrGetName(exprhdlr), exprhdlrname) == 0, "expect expression handler %s, got %s\n",
      exprhdlrname, SCIPexprhdlrGetName(exprhdlr));

   /* print simplified expression */
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   return SCIP_OKAY;
}

/** helper function to set bounds */
static
SCIP_RETCODE chgBounds(
   SCIP_VAR*            var,                /**< variable */
   SCIP_Real            lb,                 /**< new lower bound */
   SCIP_Real            ub                  /**< new upper bound */
   )
{
   cr_assert(lb <= ub);
   SCIP_CALL( SCIPchgVarLbGlobal(scip, var, lb) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, var, ub) );

   return SCIP_OKAY;
}

/**< helper function for creating an expression and checking its monotonicity */
static
SCIP_RETCODE testMonotonicity(
   int                  i,                  /**< i-th children of the root expression */
   SCIP_MONOTONE        expectedres         /**< expected result */
   )
{
   SCIP_MONOTONE monotonicity;

   cr_assert(i < SCIPexprGetNChildren(expr));

   /* evaluate monotonicity */
   SCIP_CALL( SCIPcallExprMonotonicity(scip, expr, i, &monotonicity) );

   /* check monotonicity */
   cr_expect(monotonicity == expectedres, "got %d, expected %d", monotonicity, expectedres);

   return SCIP_OKAY;
}

TestSuite(monotonicity, .init = setup, .fini = teardown);

/* check for abs expression */
Test(monotonicity, abs)
{
   SCIP_CALL( createExpr("abs(<x>[C])", "abs") );

   SCIP_CALL( chgBounds(x, -3.0, 2.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_UNKNOWN) );

   SCIP_CALL( chgBounds(x, 1.0, 2.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );

   SCIP_CALL( chgBounds(x, -2.0, -1.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
}

/* check for cosine expression */
Test(monotonicity, cos)
{
   SCIP_CALL( createExpr("cos(<x>[C])", "cos") );

   SCIP_CALL( chgBounds(x, 0.0, M_PI) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );

   SCIP_CALL( chgBounds(x, -2.0*SCIPepsilon(scip), M_PI) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_UNKNOWN) );

   SCIP_CALL( chgBounds(x, -5.0*M_PI, -4.0*M_PI) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );

   SCIP_CALL( chgBounds(x, -6.0*M_PI, -5.0*M_PI) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
}

/* check for entropy expression */
Test(monotonicity, entropy)
{
   SCIP_CALL( createExpr("entropy(<x>[C])", "entropy") );

   SCIP_CALL( chgBounds(x, exp(-1.0) - 0.1, exp(-1.0) + 0.1) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_UNKNOWN) );

   SCIP_CALL( chgBounds(x, -10.0, exp(-1.0)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );

   SCIP_CALL( chgBounds(x, exp(-1.0), 10.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
}

/* check for exp expression */
Test(monotonicity, exp)
{
   SCIP_CALL( createExpr("exp(<x>[C])", "exp") );

   SCIP_CALL( chgBounds(x, -SCIPinfinity(scip), SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
}

/* check for log expression */
Test(monotonicity, log)
{
   SCIP_CALL( createExpr("log(<x>[C])", "log") );

   SCIP_CALL( chgBounds(x, -SCIPinfinity(scip), SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
}

/* check for pow expressions */
Test(monotonicity, pow)
{
   SCIP_CALL( createExpr("<x>[C]^2", "pow") );

   SCIP_CALL( chgBounds(x, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( chgBounds(x, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
   SCIP_CALL( chgBounds(x, -1.0, 1.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_UNKNOWN) );

   SCIP_CALL( createExpr("<x>[C]^3", "pow") );

   SCIP_CALL( chgBounds(x, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( chgBounds(x, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( chgBounds(x, -1.0, 1.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );

   SCIP_CALL( createExpr("<x>[C]^(0.5)", "pow") );

   SCIP_CALL( chgBounds(x, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( chgBounds(x, -1.0, 1.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );

   SCIP_CALL( createExpr("<x>[C]^(-1)", "pow") );

   SCIP_CALL( chgBounds(x, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
   SCIP_CALL( chgBounds(x, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
   SCIP_CALL( chgBounds(x, -1.0, 1.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_UNKNOWN) );

   SCIP_CALL( createExpr("<x>[C]^(-2)", "pow") );

   SCIP_CALL( chgBounds(x, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
   SCIP_CALL( chgBounds(x, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( chgBounds(x, -1.0, 1.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_UNKNOWN) );

   SCIP_CALL( createExpr("<x>[C]^(-0.5)", "pow") );

   SCIP_CALL( chgBounds(x, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
   SCIP_CALL( chgBounds(x, -1.0, 1.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );

   SCIP_CALL( createExpr("signpower(<x>,2)", "signpower") );

   SCIP_CALL( chgBounds(x, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( chgBounds(x, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( chgBounds(x, -1.0, 1.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
}

/* check for product expression with two factors */
Test(monotonicity, prod_two)
{
   SCIP_CALL( createExpr("<x>[C] * <y>[C]", "prod") );

   /*
    * all positive
    */
   SCIP_CALL( chgBounds(x, 1.0, SCIPinfinity(scip)) );
   SCIP_CALL( chgBounds(y, 3.0, SCIPinfinity(scip)) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( testMonotonicity(1, SCIP_MONOTONE_INC) );

   /*
    * all negative
    */
   SCIP_CALL( chgBounds(x, -2.0, -1.0) );
   SCIP_CALL( chgBounds(y, -3.0, -2.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
   SCIP_CALL( testMonotonicity(1, SCIP_MONOTONE_DEC) );

   /*
    * mixed
    */
   SCIP_CALL( chgBounds(x, 1.0, 3.0) );
   SCIP_CALL( chgBounds(y, -1.0, 2.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_UNKNOWN) );
   SCIP_CALL( testMonotonicity(1, SCIP_MONOTONE_INC) );
}

/* check for product expression with three factors */
Test(monotonicity, prod_three)
{
   SCIP_CALL( createExpr("<x>[C] * <y>[C] * <z>[C]", "prod") );

   /*
    * no nonpositive
    */
   SCIP_CALL( chgBounds(x, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( chgBounds(y, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( chgBounds(z, 0.0, SCIPinfinity(scip)) );

   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( testMonotonicity(1, SCIP_MONOTONE_INC) );
   SCIP_CALL( testMonotonicity(2, SCIP_MONOTONE_INC) );

   /*
    * two negatives
    */
   SCIP_CALL( chgBounds(x, -SCIPinfinity(scip), -1.0) );
   SCIP_CALL( chgBounds(y, -SCIPinfinity(scip), -1.0) );
   SCIP_CALL( chgBounds(z, -SCIPinfinity(scip), -1.0) );

   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( testMonotonicity(1, SCIP_MONOTONE_INC) );
   SCIP_CALL( testMonotonicity(2, SCIP_MONOTONE_INC) );

   /*
    * mixed
    */
   SCIP_CALL( chgBounds(x, 0.0, 1.0) );
   SCIP_CALL( chgBounds(y, -3.0, -2.0) );
   SCIP_CALL( chgBounds(z, -4.0, -2.0) );

   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( testMonotonicity(1, SCIP_MONOTONE_DEC) );
   SCIP_CALL( testMonotonicity(2, SCIP_MONOTONE_DEC) );

   SCIP_CALL( chgBounds(x, 0.0, 1.0) );
   SCIP_CALL( chgBounds(y, 2.0, 4.0) );
   SCIP_CALL( chgBounds(z, -1.0, 1.0) );

   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_UNKNOWN) );
   SCIP_CALL( testMonotonicity(1, SCIP_MONOTONE_UNKNOWN) );
   SCIP_CALL( testMonotonicity(2, SCIP_MONOTONE_INC) );
}

/* check for sin expression */
Test(monotonicity, sin)
{
   SCIP_CALL( createExpr("sin(<x>[C])", "sin") );

   SCIP_CALL( chgBounds(x, -M_PI/2.0, M_PI/2.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );

   SCIP_CALL( chgBounds(x, -M_PI/2.0 - 2.0*SCIPepsilon(scip), M_PI/2.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_UNKNOWN) );

   SCIP_CALL( chgBounds(x, M_PI/2.0, 3.0*M_PI/2.0) );
   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_DEC) );
}

/* check for sum expression */
Test(monotonicity, sum)
{
   SCIP_CALL( createExpr("(<x>[C])^2 + 2.0 * <y>[C] - 3.0 * <z>[C]", "sum") );

   SCIP_CALL( testMonotonicity(0, SCIP_MONOTONE_INC) );
   SCIP_CALL( testMonotonicity(1, SCIP_MONOTONE_INC) );
   SCIP_CALL( testMonotonicity(2, SCIP_MONOTONE_DEC) );
}

/* check for value expression */
Test(monotonicity, value)
{
   SCIP_CALL( createExpr("-1.3", "val") );
   SCIP_CALL( testMonotonicity(-1, SCIP_MONOTONE_CONST) );
}

/* check for var expression */
Test(monotonicity, var)
{
   SCIP_CALL( createExpr("<x>[C]", "var") );
   SCIP_CALL( testMonotonicity(-1, SCIP_MONOTONE_INC) );
}
