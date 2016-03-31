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

#include <string.h>

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_sumprod.h"

/** macro to check the return of tests
 *
 *  @note assumes the existence of SCIP_RETCODE retcode
 */
#define CHECK_TEST(x)                            \
   do                                            \
   {                                             \
      retcode = (x);                             \
      if( retcode != SCIP_OKAY )                 \
      {                                          \
         printf("Unit test " #x " failed\n");    \
         SCIPprintError(retcode);                \
         return -1;                              \
      }                                          \
   }                                             \
   while( FALSE )

/* BIG = 75000 produces a seg fault because of stack overflow */
#define BIG 10000

/*
 * Local methods for tests
 */


/*
 * TESTS
 */

static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(walk_printnode)
{
   assert(expr != NULL);

   printf("stage %d node %p curchild %d %s\n", stage, (void*)expr, SCIPgetConsExprExprWalkCurrentChild(expr), SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)));

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

/** test freeing method */
static
SCIP_RETCODE testWalk(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* currently expr constraints cannot be created */
   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create expression x*y*(-5) and walk it
    * TODO: do this on a bit more interesting expression (at least one more depth)
    */
   {
      SCIP_CONSEXPR_EXPR* expr_x;
      SCIP_CONSEXPR_EXPR* expr_y;
      SCIP_CONSEXPR_EXPR* expr_5;
      SCIP_CONSEXPR_EXPR* expr_xy5;
      EXPRCOLLECT collect;

      /* create expressions for variables x and y */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );

      /* create expression for constant 5 */
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr_5, -5.0) );

      /* create expression for product of 5, x, and y (TODO should have something to add children to an existing product expr) */
      {
         SCIP_CONSEXPR_EXPR* xy5[3] = {expr_x, expr_y, expr_5};
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_xy5, 3, xy5, NULL, 2.0) );
      }

      /* release leaf expressions (this should not free them yet, as they are captured by prod_xy5) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_5) );

      /* print expression */
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr_xy5, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* print walk */
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_xy5, walk_printnode, walk_printnode, walk_printnode, walk_printnode, NULL) );

      /* collect expression during walk (in initnode stage) and check that they come in the expected order */
      collect.next = 0;
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_xy5, walk_collect, NULL, NULL, NULL, &collect) );
      assert(collect.next == 4);
      assert(collect.e[0] == expr_xy5);
      assert(collect.e[1] == SCIPgetConsExprExprChildren(expr_xy5)[0]);
      assert(collect.e[2] == SCIPgetConsExprExprChildren(expr_xy5)[1]);
      assert(collect.e[3] == SCIPgetConsExprExprChildren(expr_xy5)[2]);


      /* release product expression (this should free the product and its children) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_xy5) );
   }

   return SCIP_OKAY;
}

/** test freeing method */
static
SCIP_RETCODE testFree(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;

   conshdlr = NULL;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* currently expr constraints cannot be created */
   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create expression 5*x*y and free it */
   {
      SCIP_CONSEXPR_EXPR* expr_x;
      SCIP_CONSEXPR_EXPR* expr_y;
      SCIP_CONSEXPR_EXPR* expr_5;
      SCIP_CONSEXPR_EXPR* expr_xy5;

      /* create expressions for variables x and y */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );

      /* create expression for constant 5 */
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr_5, 5.0) );

      /* create expression for product of 5, x, and y (TODO should have something to add children to an existing product expr) */
      {
         SCIP_CONSEXPR_EXPR* xy5[3] = {expr_x, expr_y, expr_5};
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_xy5, 3, xy5, NULL, 2.0) );
      }

      /* release leaf expressions (this should not free them yet, as they are captured by prod_xy5) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_5) );

      /* release product expression (this should free the product and its children) */
      printf("freeing simple expression\n");
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_xy5) );
   }

   /* create expression sum x_i where x_i = x for all i and free it */
   {
      int i;
      SCIP_CONSEXPR_EXPR* exprs[BIG];
      SCIP_Real           coefs[BIG];
      SCIP_CONSEXPR_EXPR* expr_x;
      SCIP_CONSEXPR_EXPR* sumexpr;

      /* create expressions for variables x */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );

      for( i = 0; i < BIG; i++ )
      {
         exprs[i] = expr_x;
         coefs[i] = i;
      }
      printf("finish big loop\n");

      /* create expression for sum */
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, BIG, exprs, coefs, -1.0) );

      /* release leaf expressions (this should not free them yet, as they are captured by sumexpr) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );

      /* release sum expression (this should free the sum and its children) */
      printf("freeing large expression\n");
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   }

   /* create expression sum x_i where x_i = x for all i but now it is deep and free it */
   {
      int i;
      SCIP_CONSEXPR_EXPR* sumexprs[BIG];
      SCIP_CONSEXPR_EXPR* expr_x;
      SCIP_CONSEXPR_EXPR* sumexpr;

      /* create expressions for variables x */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );

      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexprs[0], 1, &expr_x, NULL, 0.0) );
      for( i = 1; i < BIG; i++ )
      {
         /* create expressions for sum */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexprs[i], 1, &sumexprs[i-1], NULL, 1.0 * i) );
      }
      printf("finish big loop\n");

      /* create expressions for sum */
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, 1, &sumexprs[BIG-1], NULL, 1.0 * BIG) );

      /* release leaf expressions (this should not free them yet, as they are captured by sumexpr) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
      for( i = 0; i < BIG; i++ )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexprs[i]) );
      }

      /* release sum expression (this should free the sum and its children) */
      printf("freeing deep expression\n");
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   }

   /* create expression 5*x*y and free it */
   {
      SCIP_CONSEXPR_EXPR* expr_x;
      SCIP_CONSEXPR_EXPR* expr_y;
      SCIP_CONSEXPR_EXPR* prodexprs[BIG];
      SCIP_CONSEXPR_EXPR* xysum[3];
      SCIP_CONSEXPR_EXPR* exprs[BIG];
      SCIP_CONSEXPR_EXPR* sumexpr;
      SCIP_CONSEXPR_EXPR* crazyexpr;
      SCIP_Real           coefs[BIG];
      int i;

      /* create expressions for variables x and y */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );

      /* create long sum */
      for( i = 0; i < BIG; i++ )
      {
         exprs[i] = expr_x;
         coefs[i] = i;
      }
      printf("finish big loop\n");

      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, BIG, exprs, coefs, -1.0) );

      /* create deep product of y * x * long sum */
      xysum[0] = expr_x; xysum[1] = expr_y; xysum[2] = sumexpr;
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &prodexprs[0], 3, (SCIP_CONSEXPR_EXPR**)&xysum, NULL, 1.0) );

      for( i = 1; i < BIG; i++ )
      {
         /* create expressions for product */
         if( BIG % 2 == 1 )
         {
            xysum[0] = sumexpr; xysum[1] = prodexprs[i-1]; xysum[2] = expr_y;
            SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &prodexprs[i], 3, (SCIP_CONSEXPR_EXPR**)&xysum, NULL, 1.0 * i) );
         }
         else
         {
            exprs[0] = prodexprs[i-1]; exprs[1] = prodexprs[i-1]; exprs[2] = expr_y;
            SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &prodexprs[i], BIG, exprs, coefs, -1.0) );
         }

      }

      /* create expressions for crazy expr */
      xysum[0] = sumexpr; xysum[1] = prodexprs[BIG-1]; xysum[2] = sumexpr;
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &crazyexpr, 3, (SCIP_CONSEXPR_EXPR**)&xysum, NULL, 1.0 * BIG) );

      /* release leaf expressions (this should not free them yet, as they are captured by sumexpr) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
      for( i = 0; i < BIG; i++ )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexprs[i]) );
      }

      /* release product expression (this should free the product and its children) */
      printf("freeing complicated expression\n");
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &crazyexpr) );
   }

   /* give up scip */
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   /* test something? */
   BMScheckEmptyMemory();

   return SCIP_OKAY;
}



/** main function */
int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-consexpr ===========\n");
   printf("=opt=  unittest-consexpr 0\n\n");

   CHECK_TEST( testFree() );

   CHECK_TEST( testWalk() );

   /* for automatic testing output the following */
   printf("SCIP Status        : all tests passed\n");
   printf("Ignore the following:\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return 0;
}
