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

/* we assume that asserts are always executed */
#ifdef NDEBUG
#undef NDEBUG
#endif

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
#define BIG 100

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


   /* create expression x + (-2)*x/y*(-5) and walk it */
   {
      SCIP_CONSEXPR_EXPR* expr_x;
      SCIP_CONSEXPR_EXPR* expr_y;
      SCIP_CONSEXPR_EXPR* expr_5;
      SCIP_CONSEXPR_EXPR* expr_xy5;
      SCIP_CONSEXPR_EXPR* expr_sum;
      EXPRCOLLECT collect;

      /* create expressions for variables x and y */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );

      /* create expression for constant -5 */
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr_5, -5.0) );

      /* create expression for product of -5, x, and y, and constant factor -2 (TODO should have something to add children to an existing product expr) */
      {
         SCIP_Real exponents[3] = {1.0, -1.0, 1.0};
         SCIP_CONSEXPR_EXPR* xy5[3] = {expr_x, expr_y, expr_5};
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_xy5, 3, xy5, exponents, -2.0) );
      }

      /* create expression for sum of x and product (expr_xy5) */
      {
         SCIP_CONSEXPR_EXPR* terms[2] = {expr_x, expr_xy5};
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_sum, 2, terms, NULL, 0) );
      }

      /* release leaf expressions (this should not free them yet, as they are captured by expr_xy5) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_5) );

      /* print expression */
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr_sum, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* print walk */
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, walk_printnode, walk_printnode, walk_printnode, walk_printnode, NULL) );

      /* collect expression during walk (in initnode stage) and check that they come in the expected order */
      collect.next = 0;
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, walk_collect, NULL, NULL, NULL, &collect) );
      assert(collect.next == 6);
      assert(collect.e[0] == expr_sum);
      if( SCIPgetConsExprExprChildren(expr_sum)[0] == expr_x )
      {
         /* if expr_sum holds the children in the same order as they have been given to SCIPcreateConsExprExprSum */
         assert(collect.e[1] == expr_x);
         assert(collect.e[2] == expr_xy5);
         assert(collect.e[3] == SCIPgetConsExprExprChildren(expr_xy5)[0]);
         assert(collect.e[4] == SCIPgetConsExprExprChildren(expr_xy5)[1]);
         assert(collect.e[5] == SCIPgetConsExprExprChildren(expr_xy5)[2]);
      }
      else
      {
         /* if expr_sum swapped the children */
         assert(collect.e[1] == expr_xy5);
         assert(collect.e[2] == SCIPgetConsExprExprChildren(expr_xy5)[0]);
         assert(collect.e[3] == SCIPgetConsExprExprChildren(expr_xy5)[1]);
         assert(collect.e[4] == SCIPgetConsExprExprChildren(expr_xy5)[2]);
         assert(collect.e[5] == expr_x);
      }

      /* release product expression (this should free the product and its children) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sum) );
   }

   /* test walk in walk; counts subexpressions of x + (-2)*x/y*(-5) */
   {
      SCIP_CONSEXPR_EXPR* expr_x;
      SCIP_CONSEXPR_EXPR* expr_y;
      SCIP_CONSEXPR_EXPR* expr_5;
      SCIP_CONSEXPR_EXPR* expr_xy5;
      SCIP_CONSEXPR_EXPR* expr_sum;
      int nnodes;

      /* create expressions for variables x and y */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_x, x) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr_y, y) );

      /* create expression for constant -5 */
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr_5, -5.0) );

      /* create expression for product of -5, x, and y, and constant factor -2 */
      {
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_xy5, 1, &expr_x, NULL, -2.0) );
         SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_xy5, expr_y, -1.0) );
         SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr_xy5, expr_5, 1.0) );
      }

      /* create expression for sum of x and product (expr_xy5) */
      {
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_sum, 1, &expr_x, NULL, 0) );
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr_sum, expr_xy5, 1.0) );
      }

      /* release leaf expressions (this should not free them yet, as they are captured by expr_xy5) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_5) );
      assert(expr_y != NULL);
      assert(expr_x != NULL);
      assert(expr_5 != NULL);

      /* returns sum_{expr in expr_sum} 1 */
      nnodes = 0;
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, walk_count, NULL, NULL, NULL, &nnodes) );
      assert(nnodes == 6);

      /* returns sum_{expr in expr_sum} nchild(expr) by recursively calling walk_count */
      nnodes = 0;
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, NULL, NULL, NULL, walk_count_all, &nnodes) );
      assert(nnodes == 14);

      /* release product expression (this should free the product and its children) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_sum) );
   }

   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** test parsing method */
static
SCIP_RETCODE testParse(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;

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
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z(  ", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create expression x/y*(5) from string */
   {
      SCIP_CONSEXPR_EXPR* expr_xy5;
      const char* input = "<x>[C] / <y>[I] *(-5)";

      /* create expression for product of -5, x, and y */
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr_xy5) == SCIP_OKAY);

      /* print expression */
      SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr_xy5, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* release expression */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_xy5) );
   }

   /* create more crazy expressions from string */
   {
      SCIP_CONSEXPR_EXPR* crazyexpr;
      const char* input = "-<x>[C] * <y>[I] ^(-1) + (<x>[C]+<y>[C])^2";

      /* create expression */
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &crazyexpr) == SCIP_OKAY);

      /* print expression */
      SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
      SCIP_CALL( SCIPprintConsExprExpr(scip, crazyexpr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* release expression (this should free the expression and its children) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &crazyexpr) );
   }

   /* create even more crazy expression from string */
   {
      SCIP_CONSEXPR_EXPR* crazyexpr;
      SCIP_SOL* crazysol;
      const char* input = "<x>*<y>^2/<x>^4 - 2*<x>*(3+5*<x>-2*<y>)*(<x>+<y>)^(-3.5)";
      const SCIP_Real vals[3][2] = { {1.0, 2.0}, {0.0, 0.0}, {42.0, 17.0} };
      int p;

#define CRAZYEVAL(x, y) ((x)*pow(y,2)/pow(x,4) - 2*(x)*(3+5*(x)-2*(y)) * pow((x)+(y), -3.5))

      /* create expression */
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &crazyexpr) == SCIP_OKAY);

      /* print expression */
      SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
      SCIP_CALL( SCIPprintConsExprExpr(scip, crazyexpr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* test expression by evaluating it */
      SCIP_CALL( SCIPcreateSol(scip, &crazysol, NULL) );

      for( p = 0; p < 3; ++p )
      {
         SCIP_Real expvalue = CRAZYEVAL(vals[p][0], vals[p][1]);

         SCIP_CALL( SCIPsetSolVal(scip, crazysol, x, vals[p][0]) );
         SCIP_CALL( SCIPsetSolVal(scip, crazysol, y, vals[p][1]) );

         SCIP_CALL( SCIPevalConsExprExpr(scip, crazyexpr, crazysol, 0) );
         SCIPinfoMessage(scip, NULL, "value for x=%g y=%g is %g, expected: %g\n", vals[p][0], vals[p][1], SCIPgetConsExprExprValue(crazyexpr), expvalue);
         if( SCIPgetConsExprExprValue(crazyexpr) == SCIP_INVALID )
            assert(!SCIPisFinite(expvalue));
         else
            assert(SCIPisEQ(scip, SCIPgetConsExprExprValue(crazyexpr), expvalue));
      }

      /* release expression (this should free the expression and its children) */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &crazyexpr) );

      /* release solution */
      SCIP_CALL( SCIPfreeSol(scip, &crazysol) );
   }

   /* create expression from string with unusual variable name */
   {
      SCIP_CONSEXPR_EXPR* expr;
      const char* input = "(<x> - <y>) /   <z(  >^2";

      /* parse */
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) == SCIP_OKAY);

      /* print expression */
      SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* release expression */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* create constraint holding x/y*(5) <= 1 from string */
   {
      SCIP_CONS* consexpr_xy5;
      SCIP_Bool success;
      const char* input = "[expr] <test1>: <x>[C] / <y>[I] *(5) >= 1;";

      /* parse constraint */
      success = FALSE;
      SCIP_CALL( SCIPparseCons(scip, &consexpr_xy5, input,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      assert(success);

      /* print constraint */
      SCIPinfoMessage(scip, NULL, "printing constraint %s after parsing from string:", input);
      SCIP_CALL( SCIPprintCons(scip, consexpr_xy5, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &consexpr_xy5) );
   }

   /* create constraint holding 1 <= x/y*(5) - z <= 2 from string */
   {
      SCIP_CONS* cons;
      SCIP_Bool success;
      const char* input = "[expr] <test2>: 1 <= <x>[C] / <y>[I] *(5) - <x> <= 2;";

      /* parse constraint */
      success = FALSE;
      SCIP_CALL( SCIPparseCons(scip, &cons, input,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      assert(success);

      /* print constraint */
      SCIPinfoMessage(scip, NULL, "printing constraint %s after parsing from string:", input);
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* try to create expressions from invalid strings */
   {
      SCIP_CONSEXPR_EXPR* e;

      SCIPmessageSetErrorPrinting(NULL, NULL);

      /* there is no variable with name "xx" */
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<xx>", NULL, &e) == SCIP_READERROR);

      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"5/<donothave> ", NULL, &e) == SCIP_READERROR);

      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> +-*5 ", NULL, &e) == SCIP_READERROR);

      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> / (<y>-5 ", NULL, &e) == SCIP_READERROR);

      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"donothave(<x>) ", NULL, &e) == SCIP_READERROR);

      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"donothave(<x> ", NULL, &e) == SCIP_READERROR);

      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"val(1) ", NULL, &e) == SCIP_READERROR);

      SCIPmessageSetErrorPrintingDefault();
   }

   /* these expressions shouldn't fail */
   {
      SCIP_CONSEXPR_EXPR* e;

      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"-5+3*<x>", NULL, &e) == SCIP_OKAY);
      assert(SCIPgetConsExprExprHdlr(e) == SCIPgetConsExprExprHdlrSum(conshdlr));
      assert(SCIPgetConsExprExprSumConstant(e) == -5.0);
      assert(SCIPgetConsExprExprSumCoefs(e)[0] == 3.0);
      assert(SCIPgetConsExprExprNChildren(e) == 1);
      assert(SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(e)[0]) == SCIPgetConsExprExprHdlrVar(conshdlr));
      assert(SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(e)[0]) == x);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );

      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x>", NULL, &e) == SCIP_OKAY);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> +5*<y>", NULL, &e) == SCIP_OKAY);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> + <x>*5*<y>", NULL, &e) == SCIP_OKAY);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> +5", NULL, &e) == SCIP_OKAY);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> +5 + <y>", NULL, &e) == SCIP_OKAY);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
   }

   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

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

/* creates expression for f(x,y) = 0.5 * ( (x^2*y^(-1)*5^(-4))^2 * (2*x + 1)^(-1) ) */
static
SCIP_RETCODE createExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< cons_expr constraint handler */
   SCIP_VAR*             x,                  /**< problem variable */
   SCIP_VAR*             y,                  /**< problem variable */
   SCIP_CONSEXPR_EXPR**  xexpr,              /**< pointer to store variable expression */
   SCIP_CONSEXPR_EXPR**  yexpr,              /**< pointer to store variable expression */
   SCIP_CONSEXPR_EXPR**  const_expr,         /**< pointer to store constant expression */
   SCIP_CONSEXPR_EXPR**  prodexpr,           /**< pointer to store product expression */
   SCIP_CONSEXPR_EXPR**  sumexpr,            /**< pointer to store sum expression */
   SCIP_CONSEXPR_EXPR**  mainexpr            /**< pointer to store full expression */
   )
{
   SCIP_CONSEXPR_EXPR* exprs[] = {NULL, NULL, NULL};
   SCIP_Real exponents[3] = {2.0, -1.0, -4.0};
   SCIP_Real coef = 2.0;

   /* create variable and constant expressions */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, const_expr, 5.0) );

   /* create sum and product expression */
   exprs[0] = *xexpr;
   exprs[1] = *yexpr;
   exprs[2] = *const_expr;
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, prodexpr, 3, exprs, exponents, 1.0) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, sumexpr, 1, exprs, &coef, 1.0) );

   /* create main expression */
   exprs[0] = *prodexpr;
   exprs[1] = *sumexpr;
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, mainexpr, 2, exprs, exponents, 0.5) );

   return SCIP_OKAY;
}

/* helper function to check evaluation of expression created with createExpr() */
 static
SCIP_RETCODE checkExprEval(
   SCIP*                scip,                /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*  xexpr,               /**< variable expression */
   SCIP_CONSEXPR_EXPR*  yexpr,               /**< variable expression */
   SCIP_CONSEXPR_EXPR*  const_expr,          /**< constant expression */
   SCIP_CONSEXPR_EXPR*  prodexpr,            /**< product expression */
   SCIP_CONSEXPR_EXPR*  sumexpr,             /**< sum expression */
   SCIP_CONSEXPR_EXPR*  mainexpr,            /**< full expression */
   SCIP_Real            x,                   /**< x value used for evaluation */
   SCIP_Real            y,                   /**< y value used for evaluation */
   unsigned int         tag                  /**< tag used for evaluation */
   )
{
   SCIP_Real prodval;
   SCIP_Real sumval;

   prodval = pow(x,2)*pow(y,-1)*pow(5,-4);
   sumval = 2*x + 1;

   /* check values */
   if( !SCIPisEQ(scip, SCIPgetConsExprExprValue(mainexpr), 0.5 * pow(prodval,2) * pow(sumval,-1) )
      || !SCIPisEQ(scip, SCIPgetConsExprExprValue(sumexpr), sumval)
      || !SCIPisEQ(scip, SCIPgetConsExprExprValue(prodexpr), prodval)
      || !SCIPisEQ(scip, SCIPgetConsExprExprValue(xexpr), x)
      || !SCIPisEQ(scip, SCIPgetConsExprExprValue(yexpr), y)
      || !SCIPisEQ(scip, SCIPgetConsExprExprValue(const_expr), 5.0) )
      return SCIP_ERROR;

   /* check tags */
   if( SCIPgetConsExprExprEvalTag(mainexpr) != tag
      || SCIPgetConsExprExprEvalTag(sumexpr) != tag
      || SCIPgetConsExprExprEvalTag(prodexpr) != tag
      || SCIPgetConsExprExprEvalTag(xexpr) != tag
      || SCIPgetConsExprExprEvalTag(yexpr) != tag
      || SCIPgetConsExprExprEvalTag(const_expr) != tag )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

/* Test expression evaluation method */
static
SCIP_RETCODE testExpreval(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_SOL* sol;
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_CONSEXPR_EXPR* yexpr;
   SCIP_CONSEXPR_EXPR* const_expr;
   SCIP_CONSEXPR_EXPR* prodexpr;
   SCIP_CONSEXPR_EXPR* sumexpr;
   SCIP_CONSEXPR_EXPR* mainexpr;
   int i;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 10.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create all expressions */
   SCIP_CALL( createExpr(scip, conshdlr, x, y, &xexpr, &yexpr, &const_expr, &prodexpr, &sumexpr, &mainexpr) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 4.0) );

   /* evaluate main expression and check values for sub-expressions */
   printf("evaluate and check expressions\n");
   SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 1) );
   SCIP_CALL( checkExprEval(scip, xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, 2.0, 4.0, 1) );

   /* modify solution and evaluate expression with the same tag again; values should not change */
   printf("evaluate and check expressions with a modified solution but the same tag\n");
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -5.0) );
   SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 1) );
   SCIP_CALL( checkExprEval(scip, xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, 2.0, 4.0, 1) );

   /* evaluate expression with a different tag; values should have changed */
   printf("evaluate expression with a new tag\n");
   SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 2) );
   SCIP_CALL( checkExprEval(scip, xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, -2.0, -5.0, 2) );

   /* evaluate solution with zero tag */
   printf("evaluate expression with a zero tag\n");
   for( i = 1; i < 100; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, i*i) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, -5.0/i) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 0) );
      SCIP_CALL( checkExprEval(scip, xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, i*i, -5.0 / i, 0) );
   }

   /* mainexpr is not defined for x = -1 or y = 0; the result should be SCIP_INVALID */
   printf("evaluate expression for an undefined point\n");
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 0.0) );
   SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 0) );
   if( SCIPgetConsExprExprValue(mainexpr) != SCIP_INVALID
      || SCIPgetConsExprExprValue(prodexpr) != SCIP_INVALID )
      return SCIP_ERROR;

   /* set values for variable expression explicitly */
   printf("evaluate expression after setting value for variable expressions\n");
   for( i = 1; i < 100; ++i )
   {
      SCIPsetConsExprExprEvalValue(xexpr, i*i, i);
      SCIPsetConsExprExprEvalValue(yexpr, 1.0 / i, i);
      SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, NULL, i) );
      SCIP_CALL( checkExprEval(scip, xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, i*i, 1.0 / i, i) );
   }

   /* free solution */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   /* release all expressions */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &const_expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &mainexpr) );

   /* release SCIP */
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** helper function to check interval evaluation of an expression */
static
SCIP_RETCODE checkExprIntEval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to evaluate */
   SCIP_Real             targetinf,          /**< target infimum */
   SCIP_Real             targetsup,          /**< target supremum */
   SCIP_Bool             empty,              /**< should the result be empty? */
   unsigned int          targettag           /**< target tag*/
   )
{
   SCIP_INTERVAL interval;

   assert(expr != NULL);

   interval = SCIPgetConsExprExprInterval(expr);

   /* check evaluation tag */
   if( targettag != SCIPgetConsExprExprEvalIntervalTag(expr) )
      return SCIP_ERROR;

   /* check if interval is and should be empty */
   if( empty )
      return SCIPintervalIsEmpty(SCIPinfinity(scip), interval) ? SCIP_OKAY : SCIP_ERROR;

   /* check interval */
   if(  SCIPintervalIsEmpty(SCIPinfinity(scip), interval)
      || !SCIPisEQ(scip, SCIPintervalGetInf(interval) - targetinf, 0.0)
      || !SCIPisEQ(scip, SCIPintervalGetSup(interval) - targetsup, 0.0) )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

/* Test expression interval evaluation method */
static
SCIP_RETCODE testExprinteval(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_CONSEXPR_EXPR* yexpr;
   SCIP_CONSEXPR_EXPR* const_expr;
   SCIP_CONSEXPR_EXPR* prodexpr;
   SCIP_CONSEXPR_EXPR* sumexpr;
   SCIP_CONSEXPR_EXPR* mainexpr;
   SCIP_CONSEXPR_EXPR* sqrtexpr;
   SCIP_INTERVAL interval;
   SCIP_Real xlb, xub, ylb, yub;
   SCIP_Real exponent;
   int i;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -5.0, 5.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create different kind of expressions */
   exponent = 0.5;
   SCIP_CALL( createExpr(scip, conshdlr, x, y, &xexpr, &yexpr, &const_expr, &prodexpr, &sumexpr, &mainexpr) );
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &sqrtexpr, 1, &xexpr, &exponent, 1.0) );

   /*
    * check interval evaluation method for constant expressions
    */
   printf("check interval-evaluation of constant expressions\n");
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, const_expr, 0) );
   SCIP_CALL( checkExprIntEval(scip, const_expr, 5.0, 5.0, FALSE, 0) );

   /*
    * check interval evaluation method for variable expressions
    */
   printf("check interval-evaluation of variable expressions\n");
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, xexpr, 0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, yexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, xexpr, 0.0, 10.0, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, yexpr, -5.0, 5.0, FALSE, 0) );

   /*
    * check interval evaluation method for sum expression
    */
   printf("check interval-evaluation of sum expression\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 7.5) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sumexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, sumexpr, 5.0, 16.0, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, xexpr, 2.0, 7.5, FALSE, 0) );

   /*
    * check interval evaluation method for product expression: (x^2 / (y*5^(-4)))
    */
   printf("check interval-evaluation of product expression\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.5) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 2.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, prodexpr, pow(5,-4.0) / 8.0 , pow(5,-4.0), FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, xexpr, 0.5, 1.0, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, yexpr, 1.0, 2.0, FALSE, 0) );

   /*
    * check interval-evaluation for a complicated expression: 0.5 * ( (x^2*y^(-1)*5^(-4))^2 * (2x + 1)^(-1) )
    */
   printf("check interval evaluation of a complicated expression\n");
   for( xub = 0.0; xub <= 10.0; xub += 1.0 )
      for( xlb = 0.0; xlb <= xub; xlb += 1.0 )
         for( yub = 1.0; yub <= 10.0; yub += 1.0 )
            for( ylb = 1.0; ylb <= yub; ylb += 1.0 )
            {
               SCIP_Real inf, sup, a, b;

               SCIP_CALL( SCIPchgVarLb(scip, x, xlb) );
               SCIP_CALL( SCIPchgVarUb(scip, x, xub) );
               SCIP_CALL( SCIPchgVarLb(scip, y, ylb) );
               SCIP_CALL( SCIPchgVarUb(scip, y, yub) );

               /* compute range of mainexpr */
               inf = (pow(5.0,-8) / 2.0) * MIN(pow(xlb,4), pow(xub,4)) * pow(yub, -2);
               sup = (pow(5.0,-8) / 2.0) * MAX(pow(xlb,4), pow(xub,4)) * pow(ylb, -2);
               a = MIN(1.0 / (1.0 + 2 * xlb), 1.0 / (1.0 + 2 * xub));
               b = MAX(1.0 / (1.0 + 2 * xlb), 1.0 / (1.0 + 2 * xub));
               inf *= b < 0 ? b : a;
               sup *= b < 0 ? a : b;

               SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, 0) );

               /* check all expressions */
               SCIP_CALL( checkExprIntEval(scip, mainexpr, MIN(inf,sup), MAX(inf,sup), FALSE, 0) );
               SCIP_CALL( checkExprIntEval(scip, sumexpr, 2*xlb + 1.0, 2*xub + 1.0, FALSE, 0) );
               inf = MIN(xlb*xlb, xub*xub) * (1.0/yub) * pow(5,-4);
               sup = MAX(xlb*xlb, xub*xub) * (1.0/ylb) * pow(5,-4);
               SCIP_CALL( checkExprIntEval(scip, prodexpr, inf, sup, FALSE, 0) );
               SCIP_CALL( checkExprIntEval(scip, xexpr, xlb, xub, FALSE, 0) );
               SCIP_CALL( checkExprIntEval(scip, yexpr, ylb, yub, FALSE, 0) );
               SCIP_CALL( checkExprIntEval(scip, const_expr, 5.0, 5.0, FALSE, 0) );
            }

   /*
    * check if interval evaluation for 1/(1+2*x)^3 or 1/y leads to infinite bounds
    */
   printf("check interval-evaluation for expressions containing functions like 1/f(x)\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, -0.5) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -0.5) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0) );

   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0) );

   SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 0.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0) );

   /* (1/y)^2 should lead to [0,inf] */
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, SCIPinfinity(scip), FALSE, 0) );

   /* (1/y)^2 should lead to [0,inf] but because of 1/(1+2*x)^3 we should get [-inf,inf] */
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 10.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0) );

   /*
    * check if interval evaluation aborts for some cases like (-1)^2
    */
   printf("check interval-evaluation for undefined expressions like (-1)^2\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, sqrtexpr, 0, 0, TRUE, 0) );

   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -0.5) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, sqrtexpr, 0, 0, TRUE, 0) );

   /* the result of sqrt([-1,2]) should be [0,sqrt(2)] and not an empty interval */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, sqrtexpr, 0, sqrt(2), FALSE, 0) );

   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, 0) );
   SCIP_CALL( checkExprIntEval(scip, sqrtexpr, 0.0, 1.0, FALSE, 0) );

   /*
    * check interval evaluation tags
    */
   printf("check interval tag behavior\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* do the expression store the box tag correctly? */
   for( i = 0; i < 10; ++i )
   {
      int tag = i % 3;
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, tag) );
      SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, 0.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, sumexpr, 1.0, 1.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, xexpr, 0.0, 0.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, yexpr, 1.0, 1.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, const_expr, 5.0, 5.0, FALSE, tag) );
   }

   /* set another tag for some subexpression; result should not change */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, 1) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, sumexpr, 1.0, 1.0, FALSE, 1) );

   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, 2) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sumexpr, 2) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, 1) ); /* this should not trigger a reevaluation */
   SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, 2) );
   SCIP_CALL( checkExprIntEval(scip, sumexpr, 1.0, 1.0, FALSE, 2) );

   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, 3) ); /* this should trigger a reevaluation */
   SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, 0.0, FALSE, 3) );
   SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, 3) );
   SCIP_CALL( checkExprIntEval(scip, sumexpr, 1.0, 1.0, FALSE, 3) );

   /* manipulate evaluation interval */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, 1) );
   SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, xexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, yexpr, 1.0, 1.0, FALSE, 1) );

   /* set bounds of x to [1,1] in xexpr */
   SCIPintervalSetBounds(&interval, 1.0, 1.0);
   SCIPsetConsExprExprEvalInterval(xexpr, &interval, 2);
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, 2) ); /* should use [1,1] for xexpr */
   SCIP_CALL( checkExprIntEval(scip, prodexpr, pow(5.0,-4), pow(5.0,-4), FALSE, 2) );
   SCIP_CALL( checkExprIntEval(scip, xexpr, 1.0, 1.0, FALSE, 2) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, 3) ); /* should use [0,0] for xexpr */
   SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, 3) );
   SCIP_CALL( checkExprIntEval(scip, xexpr, 0.0, 0.0, FALSE, 3) );

   /* release all expressions */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &const_expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &mainexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sqrtexpr) );

   /* release SCIP */
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

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

   CHECK_TEST( testWalk() );

   CHECK_TEST( testFree() );

   CHECK_TEST( testExpreval() );

   CHECK_TEST( testExprinteval() );

   CHECK_TEST( testParse() );

   /* for automatic testing output the following */
   printf("SCIP Status        : all tests passed\n");
   printf("Ignore the following:\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return 0;
}
