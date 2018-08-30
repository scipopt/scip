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
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_exp.h"
#include "scip/cons_expr_sin.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_CONSEXPR_EXPR* expr;
static SCIP_CONSEXPR_ITERATOR* it;

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

   /* create an iterator */
   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );

   /* get transformed vars and release vars */
   SCIP_CALL( SCIPgetTransformedVar(scip, xo, &x) );
   SCIP_CALL( SCIPgetTransformedVar(scip, yo, &y) );
   SCIP_CALL( SCIPgetTransformedVar(scip, zo, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &xo) );
   SCIP_CALL( SCIPreleaseVar(scip, &yo) );
   SCIP_CALL( SCIPreleaseVar(scip, &zo) );

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

   SCIPexpriteratorFree(&it);

   /* free scip and check for memory leaks */
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

TestSuite(iterator, .init = setup, .fini = teardown);

/* test BFS iterator on a tree containing single expression */
Test(iterator, bfs_single)
{

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x>", NULL, &expr) );

   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_BFS, TRUE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);

   /* reinitialize again */
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_BFS, TRUE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr);
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

   /* loop over the whole tree; please enjoy the beauty of this code */
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_BFS, TRUE) );
   for( tmp = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); tmp = SCIPexpriteratorGetNext(it) )
   {
      cr_expect(tmp == exprs[i]);
      ++i;
   }
}

/* test BFS iterator on an expression with common sub-expressions */
Test(iterator, bfs_general)
{
   SCIP_CONSEXPR_EXPR* expr_prod;
   SCIP_CONSEXPR_EXPR* expr_sin;
   SCIP_CONSEXPR_EXPR* expr_exp;
   SCIP_CONSEXPR_EXPR* expr_sum;
   SCIP_CONSEXPR_EXPR* expr_x;
   SCIP_CONSEXPR_EXPR* expr_y;
   SCIP_Real coef = 1.0;

   /* create and store expressions for exp(x+y) * sin(x+y)
    *
    * expected BFS order if allowing multiple visits: product, exp, sin, sum, sum, x, y, x, y
    * expected BFS order if disallowing multiple visits: product, exp, sin, sum, x, y
    */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_sum, 1, &expr_x, &coef, 0.0) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr_sum, expr_y, 1.0) );
   SCIP_CALL( SCIPcreateConsExprExprExp(scip, conshdlr, &expr_exp, expr_sum) );
   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr_sin, expr_sum) );
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_prod, 1, &expr_exp, 1.0) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_prod, expr_sin) );

   SCIP_CALL( SCIPexpriteratorInit(it, expr_prod, SCIP_CONSEXPRITERATOR_BFS, TRUE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_exp);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sin);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sum);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sum);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_x);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_y);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_x);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_y);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));

   SCIP_CALL( SCIPexpriteratorInit(it, expr_prod, SCIP_CONSEXPRITERATOR_BFS, FALSE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_exp);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sin);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sum);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_x);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_y);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_prod) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sin) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_exp) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sum) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
}

/* test RTOPOLOGICAL iterator on a tree containing a single expression */
Test(iterator, rtopological_single)
{
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x>", NULL, &expr) );

   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_RTOPOLOGIC, TRUE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);

   /* reinitialize again */
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_RTOPOLOGIC, TRUE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr);
}

/* test RTOPOLOGICAL iterator on a tree expression */
Test(iterator, rtopological_tree)
{
   SCIP_CONSEXPR_EXPR* exprs[6];
   SCIP_CONSEXPR_EXPR* tmp;
   int targetidx[6] = {3, 4, 5, 2, 1, 0};
   int i = 0;

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "sin((<t_x> + <t_y> + <t_z>)^2)", NULL, &expr) );

   exprs[0] = expr; /* sin */
   exprs[1] = SCIPgetConsExprExprChildren(exprs[0])[0]; /* pow */
   exprs[2] = SCIPgetConsExprExprChildren(exprs[1])[0]; /* sum */
   exprs[3] = SCIPgetConsExprExprChildren(exprs[2])[0]; /* x */
   exprs[4] = SCIPgetConsExprExprChildren(exprs[2])[1]; /* y */
   exprs[5] = SCIPgetConsExprExprChildren(exprs[2])[2]; /* z */

   /* loop over the whole tree; please enjoy the beauty of this code */
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_RTOPOLOGIC, TRUE) );
   for( tmp = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); tmp = SCIPexpriteratorGetNext(it) )
   {
      cr_expect(tmp == exprs[targetidx[i]]);
      ++i;
   }
}

/* test RTOPOLOGICAL iterator on an expression with common sub-expressions */
Test(iterator, rtopological_general)
{
   SCIP_CONSEXPR_EXPR* expr_prod;
   SCIP_CONSEXPR_EXPR* expr_sin;
   SCIP_CONSEXPR_EXPR* expr_exp;
   SCIP_CONSEXPR_EXPR* expr_sum;
   SCIP_CONSEXPR_EXPR* expr_x;
   SCIP_CONSEXPR_EXPR* expr_y;
   SCIP_Real coef = 1.0;

   /* create and store expressions for exp(x+y) * sin(x+y)
    *
    * expected RTOPOLOGICAL order if allowing multiple visits: x, y, sum, exp, x, y, sum, sin, product
    * expected RTOPOLOGICAL order if disallowing multiple visits: x, y, sum, exp, sin, product
    */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_sum, 1, &expr_x, &coef, 0.0) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr_sum, expr_y, 1.0) );
   SCIP_CALL( SCIPcreateConsExprExprExp(scip, conshdlr, &expr_exp, expr_sum) );
   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr_sin, expr_sum) );
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_prod, 1, &expr_exp, 1.0) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_prod, expr_sin) );

   SCIP_CALL( SCIPexpriteratorInit(it, expr_prod, SCIP_CONSEXPRITERATOR_RTOPOLOGIC, TRUE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr_x);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_y);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sum);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_exp);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_x);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_y);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sum);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sin);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_prod);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));

   SCIP_CALL( SCIPexpriteratorInit(it, expr_prod, SCIP_CONSEXPRITERATOR_RTOPOLOGIC, FALSE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr_x);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_y);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sum);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_exp);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sin);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_prod);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_prod) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sin) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_exp) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sum) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
}


/* test DFS iterators on a tree containing a single expression */
Test(iterator, dfs_single)
{
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<t_x>", NULL, &expr) );

   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);

   /* reinitialize again */
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr);


   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);

   /* reinitialize again */
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr);

}

/* test DFS iterator on a tree expression */
Test(iterator, dfs_tree)
{
   SCIP_CONSEXPR_EXPR* exprs[6];
   SCIP_CONSEXPR_EXPR* tmp;
   int targetidx[6] = { 0, 1, 2, 3, 4, 5 };
   int i = 0;

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "sin((<t_x> + <t_y> + <t_z>)^2)", NULL, &expr) );

   exprs[0] = expr; /* sin */
   exprs[1] = SCIPgetConsExprExprChildren(exprs[0])[0]; /* pow */
   exprs[2] = SCIPgetConsExprExprChildren(exprs[1])[0]; /* sum */
   exprs[3] = SCIPgetConsExprExprChildren(exprs[2])[0]; /* x */
   exprs[4] = SCIPgetConsExprExprChildren(exprs[2])[1]; /* y */
   exprs[5] = SCIPgetConsExprExprChildren(exprs[2])[2]; /* z */

   /* loop over the whole tree; please enjoy the beauty of this code */
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   for( tmp = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); tmp = SCIPexpriteratorGetNext(it) )
   {
      cr_expect(tmp == exprs[targetidx[i]]);
      ++i;
   }

   /* loop over the whole tree without revisits; same beauty as before */
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   for( tmp = SCIPexpriteratorGetCurrent(it), i = 0; !SCIPexpriteratorIsEnd(it); tmp = SCIPexpriteratorGetNext(it) )
   {
      cr_expect(tmp == exprs[targetidx[i]]);
      ++i;
   }
}

/* test DFS iterators on an expression with common sub-expressions */
Test(iterator, dfs_general)
{
   SCIP_CONSEXPR_EXPR* expr_prod;
   SCIP_CONSEXPR_EXPR* expr_sin;
   SCIP_CONSEXPR_EXPR* expr_exp;
   SCIP_CONSEXPR_EXPR* expr_sum;
   SCIP_CONSEXPR_EXPR* expr_x;
   SCIP_CONSEXPR_EXPR* expr_y;
   SCIP_Real coef = 1.0;

   /* create and store expressions for exp(x+y) * sin(x+y)
    *
    * expected DFS order if allowing multiple visits: product, exp, sum, x, y, sin, sum, x, y
    * expected DFS order if disallowing multiple visits: product, exp, sum, x, y, sin
    */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_sum, 1, &expr_x, &coef, 0.0) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr_sum, expr_y, 1.0) );
   SCIP_CALL( SCIPcreateConsExprExprExp(scip, conshdlr, &expr_exp, expr_sum) );
   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr_sin, expr_sum) );
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_prod, 1, &expr_exp, 1.0) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_prod, expr_sin) );

   SCIP_CALL( SCIPexpriteratorInit(it, expr_prod, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_exp);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sum);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_x);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_y);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sin);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sum);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_x);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_y);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));

   SCIP_CALL( SCIPexpriteratorInit(it, expr_prod, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   cr_expect(SCIPexpriteratorGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_exp);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sum);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_x);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_y);
   cr_expect(SCIPexpriteratorGetNext(it) == expr_sin);
   cr_expect(SCIPexpriteratorGetNext(it) == NULL);
   cr_expect(SCIPexpriteratorIsEnd(it));

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_prod) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sin) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_exp) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sum) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
}
