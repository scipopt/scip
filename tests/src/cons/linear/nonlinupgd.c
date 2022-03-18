/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
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

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(nonlinupgd, .init = setup, .fini = teardown);

/* upgrades a linear nonlinear constraint to a linear constraint */
Test(nonlinupgd, linear)
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

   cr_assert_eq(nupgdconss, 1);
   cr_assert_not_null(lincons);
   cr_expect(SCIPgetNVarsLinear(scip, lincons) == 3);
   cr_expect(SCIPgetLhsLinear(scip, lincons) == -2.5);
   cr_expect(SCIPgetRhsLinear(scip, lincons) == 1.5);

   /* check coefficients */
   for( i = 0; i < SCIPgetNVarsLinear(scip, lincons); ++i )
   {
      SCIP_VAR* var = SCIPgetVarsLinear(scip, lincons)[i];
      SCIP_Real coef = SCIPgetValsLinear(scip, lincons)[i];

      if( var == x )
         cr_expect(coef == 1.0);
      else if( var == y )
         cr_expect(coef == 2.0);
      else if( var == z )
         cr_expect(coef == -3.0);
   }

   /* release constraints and expressions */
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* tries to upgrade a quadratic nonlinear constraint to a linear constraint, which should fail */
Test(nonlinupgd, quadratic)
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

   cr_assert_eq(nupgdconss, 0);
   cr_assert_null(lincons);

   /* release constraints and expressions */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
