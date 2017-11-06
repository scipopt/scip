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

#include "scip/scip.h"
#include "scip/cons_expr.c"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;

static
void setup(void)
{

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_VAR* x;
   SCIP_VAR* y;

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );

}

static
void teardown(void)
{
   /* free scip and check for memory leaks */
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

TestSuite(commonSubexpr, .init = setup, .fini = teardown);

Test(commonSubexpr, singleConstraints)
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* children[2];

   SCIPinfoMessage(scip, NULL, "test single constraint\n");

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x>^2 * <y>^2", NULL, &children[0]) );
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x>^2 * <y>^2", NULL, &children[1]) );
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr, 2, children, 1.0) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, 1.0) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[0]) );

   SCIP_CALL( replaceCommonSubexpressions(scip, &cons, 1) );
   cr_assert(SCIPgetExprConsExpr(scip, cons)->children[0] == SCIPgetExprConsExpr(scip, cons)->children[1]);

   /* this should not change anything */
   SCIP_CALL( replaceCommonSubexpressions(scip, &cons, 1) );
   assert(SCIPgetExprConsExpr(scip, cons)->children[0] == SCIPgetExprConsExpr(scip, cons)->children[1]);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(commonSubexpr, multipleConstraints)
{
   SCIP_CONS* conss[4];
   SCIP_CONSEXPR_EXPR* expr;

   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "<x> * <y>", NULL, &expr)) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &conss[0], "cons", expr, -1.0, 1.0) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "exp(<x> * <y>)", NULL, &expr)) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &conss[1], "cons", expr, -1.0, 1.0) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "abs(exp(<x> * <y>))", NULL, &expr)) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &conss[2], "cons", expr, -1.0, 1.0) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "log(abs(exp(<x> * <y>)))", NULL, &expr)) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &conss[3], "cons", expr, -1.0, 1.0) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   SCIP_CALL( replaceCommonSubexpressions(scip, conss, 4) );

   cr_assert(SCIPgetExprConsExpr(scip, conss[0]) == SCIPgetExprConsExpr(scip, conss[1])->children[0]);
   cr_assert(SCIPgetExprConsExpr(scip, conss[0]) == SCIPgetExprConsExpr(scip, conss[2])->children[0]->children[0]);
   cr_assert(SCIPgetExprConsExpr(scip, conss[0]) == SCIPgetExprConsExpr(scip, conss[3])->children[0]->children[0]->children[0]);

   cr_assert(SCIPgetExprConsExpr(scip, conss[1]) == SCIPgetExprConsExpr(scip, conss[2])->children[0]);
   cr_assert(SCIPgetExprConsExpr(scip, conss[1]) == SCIPgetExprConsExpr(scip, conss[3])->children[0]->children[0]);

   cr_assert(SCIPgetExprConsExpr(scip, conss[2]) == SCIPgetExprConsExpr(scip, conss[3])->children[0]);

   SCIP_CALL( SCIPreleaseCons(scip, &conss[3]) );
   SCIP_CALL( SCIPreleaseCons(scip, &conss[2]) );
   SCIP_CALL( SCIPreleaseCons(scip, &conss[1]) );
   SCIP_CALL( SCIPreleaseCons(scip, &conss[0]) );
}
