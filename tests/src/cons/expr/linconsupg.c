/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   linconsupg.c
 * @brief  tests linear constraint upgrade of linear expression constraints
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_linear.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include expression and linear constraint handlers; the expression constraint handler needs to be added first */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) );

   /* get expression constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

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

TestSuite(linconsupg, .init = setup, .fini = teardown);

/* upgrades a linear expression constraint to a linear constraint */
Test(linconsupg, linear)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_CONS* lincons;
   SCIP_CONS* cons;
   SCIP_Bool changed;
   SCIP_Bool infeasible;
   int i;

   const char* input = "1.0 * <x> + 2.0 * <y> - 3.0 * <z> + 0.5";

   /* create expression constraint */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr, &simplified, &changed, &infeasible) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "test", simplified, -2.0, 2.0) );

   /* get an equivalent linear constraint */
   SCIP_CALL( SCIPgetLinearConsExpr(scip, cons, &lincons) );

   cr_assert(lincons != NULL);
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
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplified) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tries to upgrade a quadratic expression constraint to a linear constraint, which should fail */
Test(linconsupg, quadratic)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_CONS* lincons;
   SCIP_CONS* cons;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   const char* input = "<x>^2 + <y>";

   /* create expression constraint */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr, &simplified, &changed, &infeasible) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "test", simplified, -2.0, 2.0) );

   /* get an equivalent linear constraint, which should not be possible  */
   SCIP_CALL( SCIPgetLinearConsExpr(scip, cons, &lincons) );
   cr_assert(lincons == NULL);

   /* release constraints and expressions */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplified) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
