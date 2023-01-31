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

/**@file   iterator.c
 * @brief  unit test for expression iterators
 * @author Benjamin Mueller
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_EXPR* expr;
static SCIP_EXPRITER* it;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "t_x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "t_y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "t_z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create an iterator */
   SCIP_CALL( SCIPcreateExpriter(scip, &it) );

   /* NULL expression in order to free it in teardown() */
   expr = NULL;
}

static
void teardown(void)
{
   /* free expression if it has been created */
   if( expr != NULL )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }

   SCIPfreeExpriter(&it);

   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );

   /* free scip and check for memory leaks */
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

TestSuite(iterator, .init = setup, .fini = teardown);

/* test BFS iterator on a tree containing single expression */
Test(iterator, bfs_single)
{
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x>", NULL, NULL, NULL) );

   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_BFS, TRUE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));
   cr_expect(SCIPexpriterGetNext(it) == NULL);

   /* reinitialize again */
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_BFS, TRUE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr);
}

/* test BFS iterator on a tree expression */
Test(iterator, bfs_tree)
{
   SCIP_EXPR* exprs[6];
   SCIP_EXPR* tmp;
   int i = 0;

   SCIP_CALL( SCIPparseExpr(scip, &expr, "sin((<t_x> + <t_y> + <t_z>)^2)", NULL, NULL, NULL) );

   exprs[0] = expr; /* sin */
   exprs[1] = SCIPexprGetChildren(exprs[0])[0]; /* pow */
   exprs[2] = SCIPexprGetChildren(exprs[1])[0]; /* sum */
   exprs[3] = SCIPexprGetChildren(exprs[2])[0]; /* x */
   exprs[4] = SCIPexprGetChildren(exprs[2])[1]; /* y */
   exprs[5] = SCIPexprGetChildren(exprs[2])[2]; /* z */

   /* loop over the whole tree; please enjoy the beauty of this code */
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_BFS, TRUE) );
   for( tmp = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); tmp = SCIPexpriterGetNext(it) )
   {
      cr_expect(tmp == exprs[i]);
      ++i;
   }
}

/* test BFS iterator on an expression with common sub-expressions */
Test(iterator, bfs_general)
{
   SCIP_EXPR* expr_prod;
   SCIP_EXPR* expr_sin;
   SCIP_EXPR* expr_exp;
   SCIP_EXPR* expr_sum;
   SCIP_EXPR* expr_x;
   SCIP_EXPR* expr_y;
   SCIP_Real coef = 1.0;

   /* create and store expressions for exp(x+y) * sin(x+y)
    *
    * expected BFS order if allowing multiple visits: product, exp, sin, sum, sum, x, y, x, y
    * expected BFS order if disallowing multiple visits: product, exp, sin, sum, x, y
    */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_y, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &expr_sum, 1, &expr_x, &coef, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr_sum, expr_y, 1.0) );
   SCIP_CALL( SCIPcreateExprExp(scip, &expr_exp, expr_sum, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSin(scip, &expr_sin, expr_sum, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr_prod, 1, &expr_exp, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_prod, expr_sin) );

   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_BFS, TRUE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_sin);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_x);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetNext(it) == expr_x);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));

   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_BFS, FALSE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_sin);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_x);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_prod) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sin) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_exp) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sum) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );
}

/* test RTOPOLOGICAL iterator on a tree containing a single expression */
Test(iterator, rtopological_single)
{
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x>", NULL, NULL, NULL) );

   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_RTOPOLOGIC, TRUE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));
   cr_expect(SCIPexpriterGetNext(it) == NULL);

   /* reinitialize again */
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_RTOPOLOGIC, TRUE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr);
}

/* test RTOPOLOGICAL iterator on a tree expression */
Test(iterator, rtopological_tree)
{
   SCIP_EXPR* exprs[6];
   SCIP_EXPR* tmp;
   int targetidx[6] = {3, 4, 5, 2, 1, 0};
   int i = 0;

   SCIP_CALL( SCIPparseExpr(scip, &expr, "sin((<t_x> + <t_y> + <t_z>)^2)", NULL, NULL, NULL) );

   exprs[0] = expr; /* sin */
   exprs[1] = SCIPexprGetChildren(exprs[0])[0]; /* pow */
   exprs[2] = SCIPexprGetChildren(exprs[1])[0]; /* sum */
   exprs[3] = SCIPexprGetChildren(exprs[2])[0]; /* x */
   exprs[4] = SCIPexprGetChildren(exprs[2])[1]; /* y */
   exprs[5] = SCIPexprGetChildren(exprs[2])[2]; /* z */

   /* loop over the whole tree; please enjoy the beauty of this code */
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_RTOPOLOGIC, TRUE) );
   for( tmp = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); tmp = SCIPexpriterGetNext(it) )
   {
      cr_assert(i < 6);
      cr_expect(tmp == exprs[targetidx[i]]);
      ++i;
   }

   /* loop over the whole tree using DFS with stop at leave; should be same beauty as RTOPOLOGICAL */
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);
   for( tmp = SCIPexpriterGetCurrent(it), i = 0; !SCIPexpriterIsEnd(it); tmp = SCIPexpriterGetNext(it) )
   {
      cr_assert(i < 6);
      cr_expect(tmp == exprs[targetidx[i]]);
      ++i;
   }

}

/* test RTOPOLOGICAL iterator on an expression with common sub-expressions */
Test(iterator, rtopological_general)
{
   SCIP_EXPR* expr_prod;
   SCIP_EXPR* expr_sin;
   SCIP_EXPR* expr_exp;
   SCIP_EXPR* expr_sum;
   SCIP_EXPR* expr_x;
   SCIP_EXPR* expr_y;
   SCIP_Real coef = 1.0;

   /* create and store expressions for exp(x+y) * sin(x+y)
    *
    * expected RTOPOLOGICAL order if allowing multiple visits: x, y, sum, exp, x, y, sum, sin, product
    * expected RTOPOLOGICAL order if disallowing multiple visits: x, y, sum, exp, sin, product
    */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_y, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &expr_sum, 1, &expr_x, &coef, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr_sum, expr_y, 1.0) );
   SCIP_CALL( SCIPcreateExprExp(scip, &expr_exp, expr_sum, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSin(scip, &expr_sin, expr_sum, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr_prod, 1, &expr_exp, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_prod, expr_sin) );

   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_RTOPOLOGIC, TRUE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr_x);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_x);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_sin);
   cr_expect(SCIPexpriterGetNext(it) == expr_prod);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));

   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_RTOPOLOGIC, FALSE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr_x);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_sin);
   cr_expect(SCIPexpriterGetNext(it) == expr_prod);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_prod) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sin) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_exp) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sum) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );
}


/* test DFS iterators on a tree containing a single expression */
Test(iterator, dfs_single)
{
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x>", NULL, NULL, NULL) );

   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));
   cr_expect(SCIPexpriterGetNext(it) == NULL);


   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));
   cr_expect(SCIPexpriterGetNext(it) == NULL);

   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_ENTEREXPR);
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ALLSTAGES);
   cr_expect(SCIPexpriterGetCurrent(it) == expr);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_ENTEREXPR);
   cr_expect(SCIPexpriterGetNext(it) == expr);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_LEAVEEXPR);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));
}

/* test DFS iterator on a tree expression */
Test(iterator, dfs_tree)
{
   SCIP_EXPR* exprs[6];
   SCIP_EXPR* tmp;
   int targetidx[6] = { 0, 1, 2, 3, 4, 5 };
   int targetidx2[12] = { 0, 1, 2, 3, 3, 4, 4, 5, 5, 2, 1, 0 };
   int i = 0;

   SCIP_CALL( SCIPparseExpr(scip, &expr, "sin((<t_x> + <t_y> + <t_z>)^2)", NULL, NULL, NULL) );

   exprs[0] = expr; /* sin */
   exprs[1] = SCIPexprGetChildren(exprs[0])[0]; /* pow */
   exprs[2] = SCIPexprGetChildren(exprs[1])[0]; /* sum */
   exprs[3] = SCIPexprGetChildren(exprs[2])[0]; /* x */
   exprs[4] = SCIPexprGetChildren(exprs[2])[1]; /* y */
   exprs[5] = SCIPexprGetChildren(exprs[2])[2]; /* z */

   /* loop over the whole tree; please enjoy the beauty of this code */
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   for( tmp = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); tmp = SCIPexpriterGetNext(it) )
   {
      cr_assert(i < 6);
      cr_expect(tmp == exprs[targetidx[i]]);
      ++i;
   }

   /* loop over the whole tree without revisits; same beauty as before */
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );
   for( tmp = SCIPexpriterGetCurrent(it), i = 0; !SCIPexpriterIsEnd(it); tmp = SCIPexpriterGetNext(it) )
   {
      cr_assert(i < 6);
      cr_expect(tmp == exprs[targetidx[i]]);
      ++i;
   }

   /* loop over the whole tree without revisits but stop on enter & leave; more beauty */
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_LEAVEEXPR);
   for( tmp = SCIPexpriterGetCurrent(it), i = 0; !SCIPexpriterIsEnd(it); tmp = SCIPexpriterGetNext(it), ++i )
   {
      cr_assert(i < 12);
      cr_expect(tmp == exprs[targetidx2[i]]);

      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_ENTEREXPR :
         {
            SCIP_EXPRITER_USERDATA userdata = { .ptrval = tmp };
            SCIPexpriterSetCurrentUserData(it, userdata);
            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR :
            cr_expect(SCIPexpriterGetCurrentUserData(it).ptrval == tmp);
            break;

         default:
            cr_assert(0, "unexpected stage");
            break;
      }
   }
}

/* test DFS iterators on an expression with common sub-expressions */
Test(iterator, dfs_general)
{
   SCIP_EXPR* expr_prod;
   SCIP_EXPR* expr_sin;
   SCIP_EXPR* expr_exp;
   SCIP_EXPR* expr_sum;
   SCIP_EXPR* expr_x;
   SCIP_EXPR* expr_y;
   SCIP_Real coef = 1.0;

   /* create and store expressions for exp(x+y) * sin(x+y)
    *
    * expected DFS order if allowing multiple visits: product, exp, sum, x, y, sin, sum, x, y
    * expected DFS order if disallowing multiple visits: product, exp, sum, x, y, sin
    */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_y, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &expr_sum, 1, &expr_x, &coef, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr_sum, expr_y, 1.0) );
   SCIP_CALL( SCIPcreateExprExp(scip, &expr_exp, expr_sum, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSin(scip, &expr_sin, expr_sum, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr_prod, 1, &expr_exp, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_prod, expr_sin) );

   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_DFS, TRUE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriterGetParentDFS(it) == NULL);
   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetParentDFS(it) == expr_prod);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetParentDFS(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_x);
   cr_expect(SCIPexpriterGetParentDFS(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetParentDFS(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_sin);
   cr_expect(SCIPexpriterGetParentDFS(it) == expr_prod);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetParentDFS(it) == expr_sin);
   cr_expect(SCIPexpriterGetNext(it) == expr_x);
   cr_expect(SCIPexpriterGetParentDFS(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetParentDFS(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));

   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_DFS, FALSE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_x);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetNext(it) == expr_sin);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));

   /* same again, but stop when visiting a child or visited a child only */
   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_VISITEDCHILD);

   cr_expect(SCIPexpriterGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITINGCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 0);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_exp);

   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITINGCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 0);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_sum);

   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITINGCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 0);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_x);

   /* next will not stop at x, because it doesn't have a child */
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITEDCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 0);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_x);

   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITINGCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 1);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_y);

   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITEDCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 1);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_y);

   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITEDCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 0);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_sum);

   cr_expect(SCIPexpriterGetNext(it) == expr_prod);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITEDCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 0);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_exp);

   cr_expect(SCIPexpriterGetNext(it) == expr_prod);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITINGCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 1);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_sin);

   /* next will not stop at sin, since all its children have been visited already */
   cr_expect(SCIPexpriterGetNext(it) == expr_prod);
   cr_expect(SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITEDCHILD);
   cr_expect(SCIPexpriterGetChildIdxDFS(it) == 1);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_sin);

   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));


   /* now try out skip in enterexpr stage */
   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_DFS, FALSE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterSkipDFS(it) == expr_sin); /* if skip all children of sum (x,y), then we should be at sin next */
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));

   /* now try out skip in enterexpr stage, allow revisits */
   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_DFS, TRUE) );
   cr_expect(SCIPexpriterGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterSkipDFS(it) == expr_sin); /* if skip all children of sum (x,y), then we should be at sin next */
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetNext(it) == expr_x);
   cr_expect(SCIPexpriterGetNext(it) == expr_y);
   cr_expect(SCIPexpriterGetNext(it) == NULL);
   cr_expect(SCIPexpriterIsEnd(it));

   /* now try out skip in visitingchild stage */
   SCIP_CALL( SCIPexpriterInit(it, expr_prod, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD);
   cr_expect(SCIPexpriterGetCurrent(it) == expr_prod);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_exp);
   cr_expect(SCIPexpriterGetNext(it) == expr_sum);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_x);
   cr_expect(SCIPexpriterSkipDFS(it) == expr_sum);  /* we skip over x and so should be sum now, looking at the next child */
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_y);
   cr_expect(SCIPexpriterGetNext(it) == expr_prod);
   cr_expect(SCIPexpriterGetChildExprDFS(it) == expr_sin);
   cr_expect(SCIPexpriterGetNext(it) == NULL);  /* sum (as child of sin) already visited, so we are done */
   cr_expect(SCIPexpriterIsEnd(it));


   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_prod) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sin) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_exp) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sum) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );
}

Test(iterator, walk_in_walk)
{
   SCIP_EXPR* expr_x;
   SCIP_EXPR* expr_y;
   SCIP_EXPR* expr_5;
   SCIP_EXPR* expr_xy5;
   SCIP_EXPR* expr_sum;
   SCIP_EXPRITER* it2;
   int nnodes;

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_y, y, NULL, NULL) );

   /* create expression for constant -5 */
   SCIP_CALL( SCIPcreateExprValue(scip, &expr_5, -5.0, NULL, NULL) );

   /* create expression for product of -5, x, and y, and constant factor -2 */
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr_xy5, 1, &expr_x, -2.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_xy5, expr_y) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_xy5, expr_5) );

   /* create expression for sum of x and product (expr_xy5) */
   SCIP_CALL( SCIPcreateExprSum(scip, &expr_sum, 1, &expr_x, NULL, 0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr_sum, expr_xy5, 1.0) );

   /* returns sum_{expr in expr_sum} 1 */
   nnodes = 0;
   SCIP_CALL( SCIPexpriterInit(it, expr_sum, SCIP_EXPRITER_DFS, TRUE) );
   while( !SCIPexpriterIsEnd(it) )
   {
      ++nnodes;
      SCIPexpriterGetNext(it);
   }
   cr_assert(nnodes == 6);

   /* returns sum_{expr in expr_sum} nchild(expr) by repeatedly counting */
   nnodes = 0;
   SCIPcreateExpriter(scip, &it2);
   SCIP_CALL( SCIPexpriterInit(it, expr_sum, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);
   while( !SCIPexpriterIsEnd(it) )
   {
      SCIP_CALL( SCIPexpriterInit(it2, SCIPexpriterGetCurrent(it), SCIP_EXPRITER_DFS, TRUE) );
      while( !SCIPexpriterIsEnd(it2) )
      {
         ++nnodes;
         SCIPexpriterGetNext(it2);
      }
      SCIPexpriterGetNext(it);
   }
   cr_assert(nnodes == 14);
   SCIPfreeExpriter(&it2);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_5) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_xy5) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sum) );
}

#include "iterator_walk.sol"

Test(iterator, print)
{
   SCIP_EXPR* expr_x;
   SCIP_EXPR* expr_y;
   SCIP_EXPR* expr_5;
   SCIP_EXPR* expr_xy5;
   SCIP_EXPR* expr_sum;

   /* create expressions for variables x and y */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_x, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &expr_y, y, NULL, NULL) );

   /* create expression for constant -5 */
   SCIP_CALL( SCIPcreateExprValue(scip, &expr_5, -5.0, NULL, NULL) );

   /* create expression for product of -5, x, and y, and constant factor -2 */
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr_xy5, 1, &expr_x, -2.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_xy5, expr_y) );
   SCIP_CALL( SCIPappendExprChild(scip, expr_xy5, expr_5) );

   /* create expression for sum of x and product (expr_xy5) */
   SCIP_CALL( SCIPcreateExprSum(scip, &expr_sum, 1, &expr_x, NULL, 0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr_sum, expr_xy5, 1.0) );

   cr_redirect_stdout();

   SCIP_CALL( SCIPexpriterInit(it, expr_sum, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ALLSTAGES);

   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      printf("stage %d curchild %d type %s\n",
         SCIPexpriterGetStageDFS(it),
         SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITINGCHILD || SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_VISITEDCHILD ? SCIPexpriterGetChildIdxDFS(it) : -1,
         SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));

   fflush(stdout);

   cr_assert_stdout_eq_str(sol);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr_x) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_y) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_5) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_xy5) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr_sum) );
}
