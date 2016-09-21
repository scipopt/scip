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

/**@file   walk.c
 * @brief  unit test for checking walks in cons_expr
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_sumprod.h"
#include "scip/cons_expr_value.h"


/*
 * SCIP CODE
 */

static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(walk_printnode)
{
   assert(expr != NULL);

   printf("stage %d curchild %d type %s\n", stage, SCIPgetConsExprExprWalkCurrentChild(expr), SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)));

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

typedef struct
{
   SCIP_CONSEXPR_EXPR* e[10];
   int next;
} EXPRCOLLECT;

static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(walk_collect)
{
   EXPRCOLLECT* collect;
   assert(expr != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   collect = (EXPRCOLLECT*)data;
   collect->e[collect->next++] = expr;

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(walk_count)
{
   int* nnodes;
   assert(expr != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   nnodes = (int *)data;
   ++*nnodes;

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(walk_count_all)
{
   int nnodes;
   int *ntotalnodes;
   assert(expr != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_LEAVEEXPR);

   /* count number of nodes in sub-expression */
   nnodes = 0;
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, walk_count, NULL, NULL, NULL, &nnodes) );

   /* add to total nodes */
   ntotalnodes = (int *)data;
   *ntotalnodes += nnodes;

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

/* TESTS */
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_CONSEXPR_EXPR* expr_x;
static SCIP_CONSEXPR_EXPR* expr_y;
static SCIP_CONSEXPR_EXPR* expr_5;
static SCIP_CONSEXPR_EXPR* expr_xy5;
static SCIP_CONSEXPR_EXPR* expr_sum;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* currently expr constraints cannot be created */
   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* create expression x + (-2)*x/y*(-5) */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );

   /* create expression for constant -5 */
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr_5, -5.0) );

   /* create expression for product of -5, x, and y, and constant factor -2 */
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_xy5, 1, &expr_x, -2.0) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_xy5, expr_y) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_xy5, expr_5) );

   /* create expression for sum of x and product (expr_xy5) */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_sum, 1, &expr_x, NULL, 0) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr_sum, expr_xy5, 1.0) );
}

static
void teardown(void)
{
   /* release expressions, vars and free SCIP */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_5) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_xy5) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sum) );

   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

/***** TEST SUITE: all tests of the form Test(walk, xxx) belong to the same suite and share the setup and teardown *****/
TestSuite(walk, .init = setup, .fini = teardown);

/***** ACTUAL TESTS *****/

Test(walk, collection)
{
   EXPRCOLLECT collect;

   /* collect expression during walk (in initnode stage) and check that they come in the expected order */
   collect.next = 0;
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, walk_collect, NULL, NULL, NULL, &collect) );

   cr_assert(collect.next == 6, "Should have collected 6 expressions!");
   cr_assert(collect.e[0] == expr_sum);
   if( SCIPgetConsExprExprChildren(expr_sum)[0] == expr_x )
   {
      /* if expr_sum holds the children in the same order as they have been given to SCIPcreateConsExprExprSum */
      cr_assert(collect.e[1] == expr_x);
      cr_assert(collect.e[2] == expr_xy5);
      cr_assert(collect.e[3] == SCIPgetConsExprExprChildren(expr_xy5)[0]);
      cr_assert(collect.e[4] == SCIPgetConsExprExprChildren(expr_xy5)[1]);
      cr_assert(collect.e[5] == SCIPgetConsExprExprChildren(expr_xy5)[2]);
   }
   else
   {
      /* if expr_sum swapped the children */
      cr_assert(collect.e[1] == expr_xy5);
      cr_assert(collect.e[2] == SCIPgetConsExprExprChildren(expr_xy5)[0]);
      cr_assert(collect.e[3] == SCIPgetConsExprExprChildren(expr_xy5)[1]);
      cr_assert(collect.e[4] == SCIPgetConsExprExprChildren(expr_xy5)[2]);
      cr_assert(collect.e[5] == expr_x);
   }
}

Test(walk, walk_in_walk)
{
   int nnodes;

   /* returns sum_{expr in expr_sum} 1 */
   nnodes = 0;
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, walk_count, NULL, NULL, NULL, &nnodes) );
   cr_assert(nnodes == 6);

   /* returns sum_{expr in expr_sum} nchild(expr) by recursively calling walk_count */
   nnodes = 0;
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, NULL, NULL, NULL, walk_count_all, &nnodes) );
   cr_assert(nnodes == 14);
}

#include "walk.sol"

Test(walk, print)
{
   cr_redirect_stdout();

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, walk_printnode, walk_printnode, walk_printnode, walk_printnode, NULL) );

   fflush(stdout);

   cr_assert_stdout_eq_str(sol);
}
