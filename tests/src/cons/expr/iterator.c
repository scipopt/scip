/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   iterator.c
 * @brief  unit test for expression iterators
 * @author Benjamin Mueller
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr_iterator.h"
#include "scip/cons_expr.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_CONSEXPR_ITERATOR* bfsiter;
static SCIP_CONSEXPR_ITERATOR* dfsiter;
static SCIP_CONSEXPR_EXPR* expr;

static
void setup(void)
{
   static SCIP_VAR* xo;
   static SCIP_VAR* yo;
   static SCIP_VAR* zo;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &xo, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yo, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zo, "z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, xo) );
   SCIP_CALL( SCIPaddVar(scip, yo) );
   SCIP_CALL( SCIPaddVar(scip, zo) );

   /* goto presolving */
   TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE);

   /* get transformed vars and release vars */
   SCIP_CALL( SCIPgetTransformedVar(scip, xo, &x) );
   SCIP_CALL( SCIPgetTransformedVar(scip, yo, &y) );
   SCIP_CALL( SCIPgetTransformedVar(scip, zo, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &xo) );
   SCIP_CALL( SCIPreleaseVar(scip, &yo) );
   SCIP_CALL( SCIPreleaseVar(scip, &zo) );

   /* create iterator */
   SCIP_CALL( SCIPexpriteratorCreate(&bfsiter, SCIPblkmem(scip), SCIP_CONSEXPRITERATOR_BFS) );
   SCIP_CALL( SCIPexpriteratorCreate(&dfsiter, SCIPblkmem(scip), SCIP_CONSEXPRITERATOR_DFS) );

   /* NULL expression in order to free it in teardown() */
   expr = NULL;
}

static
void teardown(void)
{
   /* free expression if it has been created */
   if( expr != NULL )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   SCIPexpriteratorFree(&dfsiter);
   SCIPexpriteratorFree(&bfsiter);

   /* free scip and check for memory leaks */
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

TestSuite(iterator, .init = setup, .fini = teardown);

/* test BFS iterator on a tree containing single expression */
Test(iterator, bfs_single)
{
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x>", NULL, &expr) );

   cr_expect(SCIPexpriteratorInit(bfsiter, expr) == expr);
   cr_expect(SCIPexpriteratorIsEnd(bfsiter));
   cr_expect(SCIPexpriteratorGetNext(bfsiter) == NULL);

   /* reinitialize again */
   cr_expect(SCIPexpriteratorInit(bfsiter, expr) == expr);
}

/* test BFS iterator on a tree expression */
Test(iterator, bfs_tree)
{
   SCIP_CONSEXPR_EXPR* exprs[6];
   SCIP_CONSEXPR_EXPR* tmp;
   int i = 0;

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "sin((<t_x> + <t_y> + <t_z>)^2)", NULL, &expr) );

   exprs[0] = expr; /* sin */
   exprs[1] = SCIPgetConsExprExprChildren(exprs[0])[0]; /* pow */
   exprs[2] = SCIPgetConsExprExprChildren(exprs[1])[0]; /* sum */
   exprs[3] = SCIPgetConsExprExprChildren(exprs[2])[0]; /* x */
   exprs[4] = SCIPgetConsExprExprChildren(exprs[2])[1]; /* y */
   exprs[5] = SCIPgetConsExprExprChildren(exprs[2])[2]; /* z */

   for( tmp = SCIPexpriteratorInit(bfsiter, expr); !SCIPexpriteratorIsEnd(bfsiter); tmp = SCIPexpriteratorGetNext(bfsiter) )
   {
      cr_expect(tmp == exprs[i]);
      ++i;
   }
}

/* test DFS iterator on a tree containing a single expression */
Test(iterator, dfs_single)
{

}

/* test DFS iterator on a tree expression */
Test(iterator, dfs_tree)
{

}
