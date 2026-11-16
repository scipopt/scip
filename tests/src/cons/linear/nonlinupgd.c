/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   nonlinupgd.c
 * @brief  tests linear constraint upgrade of linear nonlinear constraints
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/* upgrades a linear nonlinear constraint to a linear constraint */
void test_nonlinupgd_linear(void)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_CONS* lincons = NULL;
   SCIP_CONS* cons;
   SCIP_Bool changed;
   SCIP_Bool infeasible;
   int nupgdconss = 0;
   int i;

   const char* input = "1.0 * <x> + 2.0 * <y> - 3.0 * <z> + 0.5";

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL) );
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "test", simplified, -2.0, 2.0) );

   SCIP_CALL( upgradeConsNonlinear(scip, cons, 3, &nupgdconss, &lincons, 1) );

   TEST_ASSERT_EQUAL(nupgdconss, 1);
   TEST_ASSERT_NOT_NULL(lincons);
   SOFT_ASSERT(SCIPgetNVarsLinear(scip, lincons) == 3);
   SOFT_ASSERT(SCIPgetLhsLinear(scip, lincons) == -2.5);
   SOFT_ASSERT(SCIPgetRhsLinear(scip, lincons) == 1.5);

   /* check coefficients; soft asserts so that a wrong coefficient is reported
    * for every variable at once instead of only for the first one
    */
   for( i = 0; i < SCIPgetNVarsLinear(scip, lincons); ++i )
   {
      SCIP_VAR* var = SCIPgetVarsLinear(scip, lincons)[i];
      SCIP_Real coef = SCIPgetValsLinear(scip, lincons)[i];

      if( var == x )
         SOFT_ASSERT(coef == 1.0);
      else if( var == y )
         SOFT_ASSERT(coef == 2.0);
      else if( var == z )
         SOFT_ASSERT(coef == -3.0);
   }

   /* release constraints and expressions */
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* tries to upgrade a quadratic nonlinear constraint to a linear constraint, which should fail */
void test_nonlinupgd_quadratic(void)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_CONS* lincons = NULL;
   SCIP_CONS* cons;
   SCIP_Bool changed;
   SCIP_Bool infeasible;
   int nupgdconss = 0;

   const char* input = "<x>^2 + <y>";

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL) );
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "test", simplified, -2.0, 2.0) );

   SCIP_CALL( upgradeConsNonlinear(scip, cons, 2, &nupgdconss, &lincons, 1) );

   TEST_ASSERT_EQUAL(nupgdconss, 0);
   TEST_ASSERT_NULL(lincons);

   /* release constraints and expressions */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_nonlinupgd_linear);
   RUN_TEST(test_nonlinupgd_quadratic);
   return UNITY_END();
}
