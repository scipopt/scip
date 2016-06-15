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
#include "scip/cons_expr_exp.h"
#include "scip/cons_expr_log.h"
#include "scip/cons_expr_abs.h"
#include "scip/struct_cons_expr.h"

#include "scip/nodesel_bfs.h" /* to be able to transform a problem */


/* declaration as in cons_expr.c */
SCIP_RETCODE replaceCommonSubexpressions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< total number of constraints */
   );

/* declaration as in cons_expr.c */
SCIP_RETCODE forwardPropCons(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONS*              cons,             /**< constraint to propagate */
   SCIP_Bool               intersect,        /**< should the new expr. bounds be intersected with the previous ones? */
   SCIP_Bool*              cutoff            /**< buffer to store whether the current node can be cutoff */
   );

/* declaration as in cons_expr.c */
SCIP_RETCODE reversePropConss(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONS**             conss,            /**< constraints to propagate */
   int                     nconss,           /**< total number of constraints to propagate */
   SCIP_Bool*              cutoff,           /**< buffer to store whether the current node can be cutoff */
   int*                    ntightenings      /**< buffer to store the number of (variable) tightenings */
   );

/* declaration as in cons_expr.c */
SCIP_RETCODE getVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_CONSEXPR_EXPR**    varexprs,         /**< array to store all variable expressions; array needs to be large enough */
   int*                    nvarexprs         /**< buffer to store the total number of variable expressions */
   );

/** macro to check the return of tests
 *
 *  @note assumes the existence of SCIP_RETCODE retcode
 */
#define CHECK_TEST(x)                            \
   do                                            \
   {                                             \
      printf("-----Start test " #x "\n");        \
      retcode = (x);                             \
      if( retcode != SCIP_OKAY )                 \
      {                                          \
         printf("Unit test " #x " failed\n");    \
         SCIPprintError(retcode);                \
         return -1;                              \
      }                                          \
      printf("-----Finish test " #x "\n\n");     \
   }                                             \
   while( FALSE )

/** macro to check interval of an expression */
#define CHECK_EXPRINTERVAL(scip,expr,a,b) (SCIPisEQ((scip), SCIPgetConsExprExprInterval((expr)).inf, (a)) && SCIPisEQ((scip), SCIPgetConsExprExprInterval((expr)).sup, (b)))

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
SCIP_DECL_CONSEXPREXPRWALK_VISIT(check_nuses)
{
   assert(expr != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   if( expr->nuses > 1 )
   {
      printf("following expression is captured too many times\n");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
      assert(expr->nuses == 1);
   }

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

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
         SCIP_CONSEXPR_EXPR* xy5[3];

         xy5[0] = expr_x;
         xy5[1] = expr_y;
         xy5[2] = expr_5;

         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr_xy5, 3, xy5, exponents, -2.0) );
      }

      /* create expression for sum of x and product (expr_xy5) */
      {
         SCIP_CONSEXPR_EXPR* terms[2];

         terms[0] = expr_x;
         terms[1] = expr_xy5;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr_sum, 2, terms, NULL, 0) );
      }

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

      /* release expressions */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_5) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_xy5) );
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

      /* release all but expr_sum */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_x) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_y) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_5) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_xy5) );

      /* returns sum_{expr in expr_sum} 1 */
      nnodes = 0;
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, walk_count, NULL, NULL, NULL, &nnodes) );
      assert(nnodes == 6);

      /* returns sum_{expr in expr_sum} nchild(expr) by recursively calling walk_count */
      nnodes = 0;
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_sum, NULL, NULL, NULL, walk_count_all, &nnodes) );
      assert(nnodes == 14);

      /* release expr_sum (this should free its children) */
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

      /* check that the expression is capture correctly */
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_xy5, check_nuses, NULL, NULL, NULL, NULL) );
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

      /* check that the expression is capture correctly */
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, check_nuses, NULL, NULL, NULL, NULL) );
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
         SCIP_CONSEXPR_EXPR* xy5[3];

         xy5[0] = expr_x;
         xy5[1] = expr_y;
         xy5[2] = expr_5;

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
   SCIP_CONSEXPR_PRINTDOTDATA* dotdata;
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

   /* evaluate main expression, print it, and check values for sub-expressions */
   printf("evaluate and check expressions\n");
   SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 1) );
   SCIP_CALL( SCIPprintConsExprExprDotInit(scip, &dotdata, NULL, SCIP_CONSEXPR_PRINTDOT_EXPRSTRING | SCIP_CONSEXPR_PRINTDOT_EVALTAG) );
   SCIP_CALL( SCIPprintConsExprExprDot(scip, dotdata, mainexpr) );
   SCIP_CALL( SCIPprintConsExprExprDotFinal(scip, &dotdata) );
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
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, const_expr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, const_expr, 5.0, 5.0, FALSE, 0) );

   /*
    * check interval evaluation method for variable expressions
    */
   printf("check interval-evaluation of variable expressions\n");
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, xexpr, FALSE, 0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, yexpr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, xexpr, 0.0, 10.0, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, yexpr, -5.0, 5.0, FALSE, 0) );

   /*
    * check interval evaluation method for sum expression
    */
   printf("check interval-evaluation of sum expression\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 7.5) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sumexpr, FALSE, 0) );
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
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 0) );
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

               SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0) );

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
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0) );

   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0) );

   SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 0.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0) );

   /* (1/y)^2 should lead to [0,inf] */
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, SCIPinfinity(scip), FALSE, 0) );

   /* (1/y)^2 should lead to [0,inf] but because of 1/(1+2*x)^3 we should get [-inf,inf] */
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 10.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0) );

   /*
    * check if interval evaluation aborts for some cases like (-1)^2
    */
   printf("check interval-evaluation for undefined expressions like (-1)^2\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, sqrtexpr, 0, 0, TRUE, 0) );

   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -0.5) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, sqrtexpr, 0, 0, TRUE, 0) );

   /* the result of sqrt([-1,2]) should be [0,sqrt(2)] and not an empty interval */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, FALSE, 0) );
   SCIP_CALL( checkExprIntEval(scip, sqrtexpr, 0, sqrt(2), FALSE, 0) );

   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, FALSE, 0) );
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
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, 0.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, sumexpr, 1.0, 1.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, xexpr, 0.0, 0.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, yexpr, 1.0, 1.0, FALSE, tag) );
      SCIP_CALL( checkExprIntEval(scip, const_expr, 5.0, 5.0, FALSE, tag) );
   }

   /* set another tag for some subexpression; result should not change */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, sumexpr, 1.0, 1.0, FALSE, 1) );

   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 2) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sumexpr, FALSE, 2) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 1) ); /* this should not trigger a reevaluation */
   SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, 2) );
   SCIP_CALL( checkExprIntEval(scip, sumexpr, 1.0, 1.0, FALSE, 2) );

   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 3) ); /* this should trigger a reevaluation */
   SCIP_CALL( checkExprIntEval(scip, mainexpr, 0.0, 0.0, FALSE, 3) );
   SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, 3) );
   SCIP_CALL( checkExprIntEval(scip, sumexpr, 1.0, 1.0, FALSE, 3) );

   /* manipulate evaluation interval */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, prodexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, xexpr, 0.0, 0.0, FALSE, 1) );
   SCIP_CALL( checkExprIntEval(scip, yexpr, 1.0, 1.0, FALSE, 1) );

   /* set bounds of x to [1,1] in xexpr */
   SCIPintervalSetBounds(&interval, 1.0, 1.0);
   SCIPsetConsExprExprEvalInterval(xexpr, &interval, 2);
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 2) ); /* should use [1,1] for xexpr */
   SCIP_CALL( checkExprIntEval(scip, prodexpr, pow(5.0,-4), pow(5.0,-4), FALSE, 2) );
   SCIP_CALL( checkExprIntEval(scip, xexpr, 1.0, 1.0, FALSE, 2) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 3) ); /* should use [0,0] for xexpr */
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


/** test duplicating expressions */
static
SCIP_RETCODE testDuplicate(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create expression 1.1*x*y/z + 3.2*x^2*y^(-5)*z + 0.5*z^3 from string */
   {
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_CONSEXPR_EXPR* duplicate;
      const char* input = "1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3";

      /* create expression from input string */
      assert(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) == SCIP_OKAY);

      /* print expression */
      SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* duplicate expression */
      SCIP_CALL( SCIPduplicateConsExprExpr(scip, expr, &duplicate) );

      /* print duplicate */
      SCIPinfoMessage(scip, NULL, "printing duplicate: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, duplicate, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate both expressions on different points */

      /* release expressions */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &duplicate) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** test copy of cons expression */
static
SCIP_RETCODE testCopy(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create constraint 1.1*x*y/z + 3.2*x^2*y^(-5)*z + 0.5*z^3 from string */
   {
      SCIP* subscip;
      SCIP_CONS* consexpr;
      SCIP_CONS* copyconsexpr;
      SCIP_Bool success;
      SCIP_Bool valid;
      const char* input = "[expr] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 == 2;";

      /* parse constraint and add it to SCIP */
      success = FALSE;
      SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      assert(success);

      SCIP_CALL( SCIPaddCons(scip, consexpr) );

      /* print constraint */
      SCIPinfoMessage(scip, NULL, "printing constraint %s after parsing from string:", input);
      SCIP_CALL( SCIPprintCons(scip, consexpr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* copy constraint: it seems scip needs to be transformed before copy :-/ */
      /* transform problem; we need a nodeselector for this */
      SCIP_CALL( SCIPincludeNodeselBfs(scip) );
      SCIP_CALL( SCIPtransformProb(scip) );

      SCIP_CALL( SCIPcreate(&subscip) );
      SCIP_CALL( SCIPcopy(scip, subscip, NULL, NULL, "copytest_", TRUE, FALSE, FALSE, &valid) );
      /*SCIP_CALL( SCIPcopyConss(scip, subscip, NULL, NULL, TRUE, FALSE, &valid) );*/
      /*SCIP_CALL( SCIPgetConsCopy(scip, subscip, consexpr, &copyconsexpr, conshdlr, NULL, NULL, "copycons",
               TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, &valid) );*/
      assert(valid);
      assert(SCIPgetNConss(subscip) == 1);

      /* get and print copied constraint; note that the subscip needs to be used to print copyconsexpr */
      copyconsexpr = SCIPgetConss(subscip)[0];
      SCIPinfoMessage(subscip, NULL, "printing copy: ");
      SCIP_CALL( SCIPprintCons(subscip, copyconsexpr, NULL) );
      SCIPinfoMessage(subscip, NULL, "\n");


      /* release transformed problem */
      SCIP_CALL( SCIPfreeTransform(scip) );

      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );

      /* release subscip */
      SCIP_CALL( SCIPfree(&subscip) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** test transforming of cons expression */
static
SCIP_RETCODE testTransform(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create constraint 1.1*x*y/z + 3.2*x^2*y^(-5)*z + 0.5*z^3 from string */
   {
      SCIP_CONS* consexpr;
      SCIP_CONS* transconsexpr;
      SCIP_Bool success;
      const char* input = "[expr] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 == 2;";

      /* parse constraint */
      success = FALSE;
      SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      assert(success);

      /* print constraint */
      SCIPinfoMessage(scip, NULL, "printing constraint %s after parsing from string:", input);
      SCIP_CALL( SCIPprintCons(scip, consexpr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* transform problem; we need a nodeselector for this */
      SCIP_CALL( SCIPincludeNodeselBfs(scip) );
      SCIP_CALL( SCIPtransformProb(scip) );

      /* transform constraint */
      SCIP_CALL( SCIPtransformCons(scip, consexpr, &transconsexpr) );

      /* print transformed constraint */
      SCIPinfoMessage(scip, NULL, "printing transform: ");
      SCIP_CALL( SCIPprintCons(scip, transconsexpr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* free transproblem */
      SCIP_CALL( SCIPreleaseCons(scip, &transconsexpr) );
      SCIP_CALL( SCIPfreeTransform(scip) );

      /* release constraints */
      SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** test CONSCHECK callback of cons expression */
static
SCIP_RETCODE testCheck(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 5.0, 2.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 5.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create constraint 1.1*x*y/z + 3.2*x^2*y^(-5)*z + 0.5*z^3 from string */
   {
      SCIP_SOL* sol;
      SCIP_CONS* consexpr;
      SCIP_Bool success;
      const char* input = "[expr] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

      /* parse constraint */
      success = FALSE;
      SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      assert(success);
      SCIP_CALL( SCIPaddCons(scip, consexpr) );

      /* create solution */
      SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

      /* create an infeasible solution */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
      SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, FALSE, FALSE, FALSE, &success) );
      if( success )
      {
         assert(FALSE);
         return SCIP_ERROR;
      }

      /* create a feasible solution */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
      SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, FALSE, FALSE, FALSE, &success) );
      if( !success )
      {
         assert(FALSE);
         return SCIP_ERROR;
      }

      /* create an undefined solution */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
      SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, FALSE, FALSE, FALSE, &success) );
      if( success )
      {
         assert(FALSE);
         return SCIP_ERROR;
      }

      /* release solution */
      SCIP_CALL( SCIPfreeSol(scip, &sol) );

      /* release constraints */
      SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** test exponential expressions */
static
SCIP_RETCODE testExp(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_SOL* sol;
   SCIP_VAR* x;
   SCIP_VAR* y;
   int i;

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
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   /* easy exponential expression */
   {
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_INTERVAL interval;
      const char* input = "exp(<x>[C]) + exp(<x>[C])";

      SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr)) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = -10; i <= 10; ++i )
      {
         /* evaluate expression */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
         SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
         assert(SCIPisRelEQ(scip, SCIPgetConsExprExprValue(expr), exp(i) + exp(i)));

         /* propagate expression */
         SCIP_CALL( SCIPchgVarLb(scip, x, i) );
         SCIP_CALL( SCIPchgVarUb(scip, x, i + 1.0 / (ABS(i) + 1)) );
         SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0) );
         interval = SCIPgetConsExprExprInterval(expr);
         assert(SCIPisRelEQ(scip, SCIPintervalGetInf(interval), 2*exp(i)));
         assert(SCIPisRelEQ(scip, SCIPintervalGetSup(interval), 2*exp(i + 1.0 / (ABS(i) + 1))));
      }

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* complicated exponential expression */
   {
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_INTERVAL interval;
      const char* input = "exp(exp(<x>[C])) * exp(<y>[C])^2";

      SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr)) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = 1; i <= 10; ++i )
      {
         /* evaluate expression */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) 1.0 / i) );
         SCIP_CALL( SCIPsetSolVal(scip, sol, y, (SCIP_Real) i) );
         SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
         assert(SCIPisRelEQ(scip, SCIPgetConsExprExprValue(expr), exp(exp(1.0 / i)) * exp(2*i)));

         /* propagate expression */
         SCIP_CALL( SCIPchgVarLb(scip, x, -1.0 / i) );
         SCIP_CALL( SCIPchgVarUb(scip, x,  1.0 / i) );
         SCIP_CALL( SCIPchgVarLb(scip, y, i) );
         SCIP_CALL( SCIPchgVarUb(scip, y, i + 1.0 / i) );
         SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0) );
         interval = SCIPgetConsExprExprInterval(expr);
         assert(SCIPisRelEQ(scip, SCIPintervalGetInf(interval), exp(exp(-1.0 / i)) * exp(2*i)));
         assert(SCIPisRelEQ(scip, SCIPintervalGetSup(interval), exp(exp(1.0 / i)) * exp(2*i + 2.0 / i)));
      }

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* free allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** test logarithmic expressions */
static
SCIP_RETCODE testLog(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_SOL* sol;
   SCIP_VAR* x;
   SCIP_VAR* y;
   int i;

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
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   /* easy logarithmic expression */
   {
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_INTERVAL interval;
      const char* input = "log(<x>[C]) + log(<x>[C])";

      SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr)) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = -10; i <= 10; ++i )
      {
         SCIP_Real xlb, xub;

         /* evaluate expression */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
         SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );

         if( i <= 0 )
            assert(SCIPgetConsExprExprValue(expr) == SCIP_INVALID);
         else
            assert(SCIPisRelEQ(scip, SCIPgetConsExprExprValue(expr), log(i) + log(i)));

         /* propagate expression */
         xlb = i;
         xub = i + 1.0 / (ABS(i) + 1);
         SCIP_CALL( SCIPchgVarLb(scip, x, xlb) );
         SCIP_CALL( SCIPchgVarUb(scip, x, xub) );
         SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0) );
         interval = SCIPgetConsExprExprInterval(expr);

         /* interval is empty if both bounds are non-positive */
         if( xub <= 0 )
            assert(SCIPintervalIsEmpty(SCIPinfinity(scip), interval));
         else
         {
            assert(SCIPisRelEQ(scip, SCIPintervalGetSup(interval), 2*log(xub)));

            if( xlb <= 0 )
               assert(SCIPisRelEQ(scip, SCIPintervalGetInf(interval), -SCIPinfinity(scip)));
            else
               assert(SCIPisRelEQ(scip, SCIPintervalGetInf(interval), 2*log(xlb)));
         }
      }

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* complicated logarithmic expression */
   {
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_INTERVAL interval;
      const char* input = "log(log(exp(<x>[C]) * exp(<y>[C])))";

      SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr)) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = 1; i <= 10; ++i )
      {
         /* evaluate expression */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
         SCIP_CALL( SCIPsetSolVal(scip, sol, y, (SCIP_Real) i + 1) );
         SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
         assert(SCIPisRelEQ(scip, SCIPgetConsExprExprValue(expr), log(2*i + 1) ));

         /* propagate expression */
         SCIP_CALL( SCIPchgVarLb(scip, x,  1.0 / i) );
         SCIP_CALL( SCIPchgVarUb(scip, x,  2.0 / i) );
         SCIP_CALL( SCIPchgVarLb(scip, y,  3.0 / i) );
         SCIP_CALL( SCIPchgVarUb(scip, y,  4.0 / i) );
         SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0) );
         interval = SCIPgetConsExprExprInterval(expr);
         assert(SCIPisRelEQ(scip, SCIPintervalGetInf(interval), log(1.0 / i + 3.0 / i)));
         assert(SCIPisRelEQ(scip, SCIPintervalGetSup(interval), log(2.0 / i + 4.0 / i)));
      }

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* free allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** test absolute expressions */
static
SCIP_RETCODE testAbs(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_SOL* sol;
   SCIP_VAR* x;
   SCIP_VAR* y;
   int i;

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
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   /* easy absolute expression */
   {
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_INTERVAL interval;
      const char* input = "abs(<x>[C]) + abs(<x>[C])";

      SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr)) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = -10; i <= 10; ++i )
      {
         /* evaluate expression */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
         SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );

         assert(SCIPisRelEQ(scip, SCIPgetConsExprExprValue(expr),  2 * REALABS(i)));

         /* propagate expression */
         SCIP_CALL( SCIPchgVarLb(scip, x, -REALABS(i)) );
         SCIP_CALL( SCIPchgVarUb(scip, x, i*i) );
         SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0) );
         interval = SCIPgetConsExprExprInterval(expr);

         assert(SCIPisRelEQ(scip, SCIPintervalGetInf(interval), 0));
         assert(SCIPisRelEQ(scip, SCIPintervalGetSup(interval), 2 * i*i));
      }

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* complicated abolute expression */
   {
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_INTERVAL interval;
      const char* input = "abs(abs(abs(<x>[C]) + abs(<y>[C])) * abs(<x>[C])^3 * abs(<y>[C]))";

      SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr)) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = -10; i <= 10; ++i )
      {
         /* evaluate expression */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
         SCIP_CALL( SCIPsetSolVal(scip, sol, y, (SCIP_Real) i) );
         SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
         assert(SCIPisRelEQ(scip, SCIPgetConsExprExprValue(expr), 2 * pow(REALABS(i),5)));

         /* propagate expression */
         SCIP_CALL( SCIPchgVarLb(scip, x, -REALABS(i)) );
         SCIP_CALL( SCIPchgVarUb(scip, x, REALABS(i)) );
         SCIP_CALL( SCIPchgVarLb(scip, y, -REALABS(i)) );
         SCIP_CALL( SCIPchgVarUb(scip, y, REALABS(i)) );
         SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0) );
         interval = SCIPgetConsExprExprInterval(expr);
         assert(SCIPisRelEQ(scip, SCIPintervalGetInf(interval), 0));
         assert(SCIPisRelEQ(scip, SCIPintervalGetSup(interval), 2 * pow(REALABS(i),5)));
      }

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* free allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** helper function to test correctness of the computed hashkeys */
static
SCIP_RETCODE checkHashkey(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr1,              /**< first expression to be tested */
   SCIP_CONSEXPR_EXPR*   expr2               /**< second expression to be tested (might be NULL) */
   )
{
   assert(expr1 != NULL);

   SCIPinfoMessage(scip, NULL, "hash key of expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr1, NULL) );
   SCIPinfoMessage(scip, NULL, " = %u\n", SCIPgetConsExprExprHashkey(scip, expr1));

   if( expr2 != NULL )
   {
      SCIPinfoMessage(scip, NULL, "hash key of expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr2, NULL) );
      SCIPinfoMessage(scip, NULL, " = %u\n", SCIPgetConsExprExprHashkey(scip, expr2));
      assert(SCIPgetConsExprExprHashkey(scip, expr1) == SCIPgetConsExprExprHashkey(scip, expr2));
   }

   return SCIP_OKAY;
}

/** test absolute expressions */
static
SCIP_RETCODE testHash(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_CONSEXPR_EXPR* yexpr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* expr2;
   int i;

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
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );


   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );

   /* value expressions */
   for( i = -2; i < 2; ++i )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr, i) );
      SCIP_CALL( checkHashkey(scip, expr, NULL) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* variable expressions */
   SCIP_CALL( checkHashkey(scip, xexpr, NULL) );
   SCIP_CALL( checkHashkey(scip, yexpr, NULL) );

   /* sum expressions */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr, 0, NULL, NULL, 2.5) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, xexpr, 14.3) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, yexpr, -2.3) );

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr2, 0, NULL, NULL, 2.5) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr2, yexpr, -2.3) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr2, xexpr, 14.3) );

   SCIP_CALL( checkHashkey(scip, expr, expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* product expressions */
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr, 0, NULL, NULL, 2.5) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, xexpr, 14.3) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, yexpr, -2.3) );

   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr2, 0, NULL, NULL, 2.5) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, yexpr, -2.3) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, xexpr, 14.3) );

   SCIP_CALL( checkHashkey(scip, expr, expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* absolute expression */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "abs(<x>)", NULL, &expr)) );
   SCIP_CALL( checkHashkey(scip, expr, NULL) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* logarithmic expression */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "exp(<x>)", NULL, &expr)) );
   SCIP_CALL( checkHashkey(scip, expr, NULL) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* exponential expression */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "log(<x>)", NULL, &expr)) );
   SCIP_CALL( checkHashkey(scip, expr, NULL) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* complicated expression */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "abs(exp(<x>*<y>^2/<x>^4) - log(2*<x>)*(3+5*<x>-2*<y>)*(<x>+<y>)^(-3.5)) + 2", NULL, &expr)) );
   SCIP_CALL( checkHashkey(scip, expr, NULL) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* equal expressions */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "<x> * <y>", NULL, &expr)) );
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "<y> * <x>", NULL, &expr2)) );
   SCIP_CALL( checkHashkey(scip, expr, expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "2*<x> + 3.3*<y>", NULL, &expr)) );
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "3.3*<y> + 2*<x>", NULL, &expr2)) );
   SCIP_CALL( checkHashkey(scip, expr, expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* free allocated memory */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** test absolute expressions */
static
SCIP_RETCODE testCommonSubexpr(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS* conss[4];
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_CONSEXPR_EXPR* expr;

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
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* single constraint */
   {
      SCIP_CONSEXPR_EXPR* children[2];
      SCIP_Real exponents[2] = {2, 2};

      SCIPinfoMessage(scip, NULL, "test single constraint\n");

      SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "<x> * <y>", NULL, &children[0])) );
      SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, "<x> * <y>", NULL, &children[1])) );
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr, 2, children, exponents, 1) );

      SCIP_CALL( SCIPcreateConsExprBasic(scip, &conss[0], "cons", expr, -1.0, 1.0) );

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[1]) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[0]) );

      SCIP_CALL( replaceCommonSubexpressions(scip, conss, 1) );
      assert(SCIPgetExprConsExpr(scip, conss[0])->children[0] == SCIPgetExprConsExpr(scip, conss[0])->children[1]);

      /* this should not change anything */
      SCIP_CALL( replaceCommonSubexpressions(scip, conss, 1) );
      assert(SCIPgetExprConsExpr(scip, conss[0])->children[0] == SCIPgetExprConsExpr(scip, conss[0])->children[1]);

      SCIP_CALL( SCIPreleaseCons(scip, &conss[0]) );
   }

   /* multiple constraint */
   {
      SCIPinfoMessage(scip, NULL, "test multiple constraints\n");

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

      assert(SCIPgetExprConsExpr(scip, conss[0]) == SCIPgetExprConsExpr(scip, conss[1])->children[0]);
      assert(SCIPgetExprConsExpr(scip, conss[0]) == SCIPgetExprConsExpr(scip, conss[2])->children[0]->children[0]);
      assert(SCIPgetExprConsExpr(scip, conss[0]) == SCIPgetExprConsExpr(scip, conss[3])->children[0]->children[0]->children[0]);

      assert(SCIPgetExprConsExpr(scip, conss[1]) == SCIPgetExprConsExpr(scip, conss[2])->children[0]);
      assert(SCIPgetExprConsExpr(scip, conss[1]) == SCIPgetExprConsExpr(scip, conss[3])->children[0]->children[0]);

      assert(SCIPgetExprConsExpr(scip, conss[2]) == SCIPgetExprConsExpr(scip, conss[3])->children[0]);

      SCIP_CALL( SCIPreleaseCons(scip, &conss[3]) );
      SCIP_CALL( SCIPreleaseCons(scip, &conss[2]) );
      SCIP_CALL( SCIPreleaseCons(scip, &conss[1]) );
      SCIP_CALL( SCIPreleaseCons(scip, &conss[0]) );
   }

   /* free allocated memory */
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}


/** test forward and reverse propagation */
static
SCIP_RETCODE testPropagation(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONS* cons;
   SCIP_Bool cutoff;
   int ntightenings;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* currently expr constraints cannot be created */
   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -2.0, 2.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -3.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -3.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* sum */
   {
      SCIP_CONSEXPR_EXPR* exprs[2];
      SCIP_Real coeffs[2] = {2.0, -1.0};

      SCIPinfoMessage(scip, NULL, "test sum expression\n");

      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &exprs[0], x) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &exprs[1], y) );
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr, 2, exprs, coeffs, 0.5) );

      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, 0.5, 1.5) );

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      SCIP_CALL( reversePropConss(scip, &cons, 1, &cutoff, &ntightenings) );
      assert(!cutoff);

      assert(SCIPisEQ(scip, expr->interval.inf, 0.5));
      assert(SCIPisEQ(scip, expr->interval.sup, 1.5));
      assert(SCIPisEQ(scip, expr->children[0]->interval.inf, -1.5));
      assert(SCIPisEQ(scip, expr->children[0]->interval.sup, 1.0));
      assert(SCIPisEQ(scip, expr->children[1]->interval.inf, -3.0));
      assert(SCIPisEQ(scip, expr->children[1]->interval.sup, 1.0));

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[1]) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[0]) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* prod */
   {
      SCIP_CONSEXPR_EXPR* exprs[2];
      SCIP_Real coeffs[2] = {2.0, -1.0};

      SCIPinfoMessage(scip, NULL, "test prod expression\n");

      SCIP_CALL( SCIPchgVarLb(scip, x, 1.0) );
      SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );
      SCIP_CALL( SCIPchgVarLb(scip, y, 2.0) );
      SCIP_CALL( SCIPchgVarUb(scip, y, 4.0) );

      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &exprs[0], x) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &exprs[1], y) );
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr, 2, exprs, coeffs, 0.5) );

      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, 0.0, 1.0) );

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      SCIP_CALL( reversePropConss(scip, &cons, 1, &cutoff, &ntightenings) );
      assert(!cutoff);

      assert(SCIPisEQ(scip, expr->interval.inf, 1.0 / 8.0));
      assert(SCIPisEQ(scip, expr->interval.sup, 1.0));
      assert(SCIPisEQ(scip, expr->children[0]->interval.inf, 1.0));
      assert(SCIPisEQ(scip, expr->children[0]->interval.sup, SQRT(8)));
      assert(SCIPisEQ(scip, expr->children[1]->interval.inf, 2.0));
      assert(SCIPisEQ(scip, expr->children[1]->interval.sup, 4.0));

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[1]) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[0]) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* abs */
   {
      SCIPinfoMessage(scip, NULL, "test abs expression\n");

      SCIP_CALL( SCIPchgVarLb(scip, x, -3.0) );
      SCIP_CALL( SCIPchgVarUb(scip, x, 4.0) );

      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "abs(<x>)", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, 1.0, 2.5) );

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      SCIP_CALL( reversePropConss(scip, &cons, 1, &cutoff, &ntightenings) );
      assert(!cutoff);

      assert(SCIPisEQ(scip, expr->children[0]->interval.inf, -2.5));
      assert(SCIPisEQ(scip, expr->children[0]->interval.sup, 2.5));

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* exp */
   {
      SCIPinfoMessage(scip, NULL, "test exp expression\n");

      SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
      SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );

      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "exp(<x>)", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, 2.0) );

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      SCIP_CALL( reversePropConss(scip, &cons, 1, &cutoff, &ntightenings) );
      assert(!cutoff);

      assert(SCIPisEQ(scip, expr->interval.inf, exp(-1)));
      assert(SCIPisEQ(scip, expr->interval.sup, 2.0));
      assert(SCIPisEQ(scip, expr->children[0]->interval.inf, -1.0));
      assert(SCIPisEQ(scip, expr->children[0]->interval.sup, log(2.0)));

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* log */
   {
      SCIPinfoMessage(scip, NULL, "test log expression\n");

      SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
      SCIP_CALL( SCIPchgVarUb(scip, x, 7.0) );

      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "log(<x>)", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, 1.0) );

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      SCIP_CALL( reversePropConss(scip, &cons, 1, &cutoff, &ntightenings) );
      assert(!cutoff);

      assert(SCIPisEQ(scip, expr->interval.inf, -1.0));
      assert(SCIPisEQ(scip, expr->interval.sup, 1.0));
      assert(SCIPisEQ(scip, expr->children[0]->interval.inf, exp(-1.0)));
      assert(SCIPisEQ(scip, expr->children[0]->interval.sup, exp(1.0)));

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* more complicated expression */
   {
      SCIP_CONSEXPR_EXPR* xexpr;
      SCIP_CONSEXPR_EXPR* yexpr;
      SCIP_CONSEXPR_EXPR* zexpr;
      SCIP_CONSEXPR_EXPR* rootexpr;
      SCIP_CONSEXPR_EXPR* prodexpr;
      SCIP_CONSEXPR_EXPR* sumexpr;
      SCIP_CONSEXPR_EXPR* logexpr;
      SCIP_CONSEXPR_EXPR* exprs[2];
      SCIP_Real coeffs[2];

      /* set variable bounds */
      SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
      SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
      SCIP_CALL( SCIPchgVarLb(scip, y, 2.0) );
      SCIP_CALL( SCIPchgVarUb(scip, y, 3.0) );
      SCIP_CALL( SCIPchgVarLb(scip, z, 1.0) );
      SCIP_CALL( SCIPchgVarUb(scip, z, 2.0) );

      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &zexpr, z) );

      /*
       * create constraint 0 <= (-x^2 + log(y)) / z <= 2
       */

      /* log(y) */
      SCIP_CALL( SCIPcreateConsExprExprLog(scip, conshdlr, &logexpr, yexpr) );

      /* x^2 */
      coeffs[0] = 2.0;
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &prodexpr, 1, &xexpr, coeffs, 1.0) );

      /* log(y) - x^2 */
      coeffs[0] = 1.0;
      exprs[0] = logexpr;
      coeffs[1] = -1.0;
      exprs[1] = prodexpr;
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, 2, exprs, coeffs, 0.0) );

      /* (-x^2 + log(y)) / z */
      coeffs[0] = 1.0;
      exprs[0] = sumexpr;
      coeffs[1] = -1.0;
      exprs[1] = zexpr;
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &rootexpr, 2, exprs, coeffs, 1.0) );

      SCIPinfoMessage(scip, NULL, "test more complicated expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, rootexpr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* create constraint */
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", rootexpr, 0.0, 2.0) );

      /* apply forward propagation */
      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      assert(CHECK_EXPRINTERVAL(scip, xexpr, -1, 1));
      assert(CHECK_EXPRINTERVAL(scip, yexpr, 2, 3));
      assert(CHECK_EXPRINTERVAL(scip, zexpr, 1, 2));
      assert(CHECK_EXPRINTERVAL(scip, logexpr, log(2), log(3)));
      assert(CHECK_EXPRINTERVAL(scip, prodexpr, 0, 1));
      assert(CHECK_EXPRINTERVAL(scip, sumexpr, log(2) - 1, log(3)));
      assert(CHECK_EXPRINTERVAL(scip, rootexpr, 0, log(3)));

      /* apply reverse propagation */
      SCIP_CALL( reversePropConss(scip, &cons, 1, &cutoff, &ntightenings) );
      assert(!cutoff);

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &rootexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &logexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &zexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   }

   /* test expressions with unbounded sub-expressions */
   {
      /* set variable bounds */
      SCIP_CALL( SCIPchgVarLb(scip, x, -SCIPinfinity(scip)) );
      SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
      SCIP_CALL( SCIPchgVarLb(scip, y, -SCIPinfinity(scip)) );
      SCIP_CALL( SCIPchgVarUb(scip, y, 0.0) );
      SCIP_CALL( SCIPchgVarLb(scip, z, 1.0) );
      SCIP_CALL( SCIPchgVarUb(scip, z, 2.0) );

      /*
       *  -5 <= x * y
       */
      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "2 * <x> * <y>", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -5.0, SCIPinfinity(scip)) );

      SCIPinfoMessage(scip, NULL, "test expressions with unbounded sub-expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      assert(CHECK_EXPRINTERVAL(scip, expr, -5.0, SCIPinfinity(scip)));

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

      /*
       *  -inf <= 2 + x - y <= inf
       */
      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "2 + <x> - <y>", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -SCIPinfinity(scip), SCIPinfinity(scip)) );

      SCIPinfoMessage(scip, NULL, "test expressions with unbounded sub-expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      assert(CHECK_EXPRINTERVAL(scip, expr, -SCIPinfinity(scip), SCIPinfinity(scip)));

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

      /*
       *  -inf <= x * y * z <= inf
       */
      SCIP_CALL( SCIPchgVarUb(scip, x, 0.0) );

      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x> * <y> * <z>", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -SCIPinfinity(scip), SCIPinfinity(scip)) );

      SCIPinfoMessage(scip, NULL, "test expressions with unbounded sub-expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      assert(CHECK_EXPRINTERVAL(scip, expr, 0.0, SCIPinfinity(scip)));

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* test if forward propagation uses stored expression bounds */
   {
      SCIPinfoMessage(scip, NULL, "test if forward propagation uses stored bounds\n");

      /* set variable bounds */
      SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
      SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
      SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) );
      SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x> + <y>", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, 0.5) );

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(!cutoff);
      assert(CHECK_EXPRINTERVAL(scip, expr, 0.0, 0.5));

      /* change variable expressions */
      expr->children[0]->interval.inf = -1.0;
      expr->children[0]->interval.sup = 0.2;
      expr->children[1]->interval.inf = -1.0;
      expr->children[1]->interval.sup = 0.2;

      /* new interval should be [0,1] intersected with [-2, 0.4] */
      SCIP_CALL( forwardPropCons(scip, cons, TRUE, &cutoff) );
      assert(!cutoff);
      assert(CHECK_EXPRINTERVAL(scip, expr, 0.0, 0.4));

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* infeasible after forward propagation */
   {
      SCIPinfoMessage(scip, NULL, "test expressions leading to an empty interval after forward propagation\n");

      /* set variable bounds */
      SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
      SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
      SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) );
      SCIP_CALL( SCIPchgVarUb(scip, y, 3.0) );

      /*
       * -7.0 <= 2 * x * y <= -6.1
       */
      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "2 * <x> * <y>", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -7.0, -6.1) );

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(cutoff);
      assert(SCIPintervalIsEmpty(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)));

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

      /*
       * -1.0 <= 1 + x / (1 + y) <= -0.1
       */
      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "1 + <x> / (1 + <y>)", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, -0.1) );

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(SCIPintervalIsEmpty(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)));
      assert(cutoff);

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

      /*
       * 0.0 <= 1 + expr(-5 * x + y^2) <= 0.9
       */
      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "1 + exp(-5 * <x> + <y>^2)", NULL, &expr) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -10.0, -0.1) );

      SCIP_CALL( forwardPropCons(scip, cons, FALSE, &cutoff) );
      assert(SCIPintervalIsEmpty(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)));
      assert(cutoff);

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   }

   /* infeasible after reverse propagation */
   {
      SCIP_CONSEXPR_EXPR* exprs[2];
      SCIP_CONSEXPR_EXPR* expr1;
      SCIP_CONSEXPR_EXPR* expr2;
      SCIP_CONSEXPR_EXPR* xexpr;
      SCIP_CONSEXPR_EXPR* yexpr;
      SCIP_CONS* cons1;
      SCIP_CONS* cons2;

      SCIPinfoMessage(scip, NULL, "test expressions leading to an empty interval after reverse propagation\n");

      /* set variable bounds */
      SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
      SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
      SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) );
      SCIP_CALL( SCIPchgVarUb(scip, y, 2.0) );

      /*
       * -1.0 <= x * y <= 1.0 and 3.5 <= x + y <= 5.0
       */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
      exprs[0] = xexpr;
      exprs[1] = yexpr;

      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr1, 2, exprs, NULL, 1.0) );
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr2, 2, exprs, NULL, 0.0) );

      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons1, "cons1", expr1, -1.0, 1.0) );
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons2, "cons2", expr2, 3.5, 5.0) );

      /* apply forward propagation for both constraints */
      SCIP_CALL( forwardPropCons(scip, cons1, FALSE, &cutoff) );
      assert(!cutoff);
      assert(CHECK_EXPRINTERVAL(scip, expr1, 0.0, 1.0));
      SCIP_CALL( forwardPropCons(scip, cons2, FALSE, &cutoff) );
      assert(!cutoff);
      assert(CHECK_EXPRINTERVAL(scip, expr2, 3.5, 4.0));

      /* reverse propagation of cons2 should lead to new bounds on x and y */
      SCIP_CALL( reversePropConss(scip, &cons2, 1, &cutoff, &ntightenings) );
      assert(!cutoff);
      assert(CHECK_EXPRINTERVAL(scip, expr2->children[0], 1.5, 2.0));
      assert(CHECK_EXPRINTERVAL(scip, expr2->children[1], 1.5, 2.0));

      /* reverse propagation of cons1 should lead to an empty interval for x */
      SCIP_CALL( reversePropConss(scip, &cons1, 1, &cutoff, &ntightenings) );
      assert(cutoff);
      assert(SCIPintervalIsEmpty(SCIPinfinity(scip), expr1->children[0]->interval));

      SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr1) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   }

   /* free allocated memory */
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** test collecting of variable expression */
static
SCIP_RETCODE testGetVarExprs(void)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* w;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* wexpr;
   SCIP_CONSEXPR_EXPR* sumexpr;
   SCIP_CONSEXPR_EXPR** varexprs;
   int nvarexprs;
   int i;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* currently expr constraints cannot be created */
   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_getvars") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -3.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -2.0, 2.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -3.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -3.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, SCIPgetNVars(scip)) );

   /*
    * test expression not containing all variables
    */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3", NULL, &expr) );

   SCIP_CALL( getVarExprs(scip, expr, varexprs, &nvarexprs) );
   assert(nvarexprs == 3);

   for( i = 0; i < nvarexprs; ++i )
   {
      assert(varexprs[i] != NULL);
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(varexprs[i])), "var") == 0);
   }

   /*
    * test expression containing all variables
    */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, 0, NULL, NULL, 0) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &wexpr, w) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, sumexpr, wexpr, 1.0) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, sumexpr, expr, 1.0) );

   SCIP_CALL( getVarExprs(scip, sumexpr, varexprs, &nvarexprs) );
   assert(nvarexprs == 4);
   assert(strcmp(SCIPvarGetName(SCIPgetConsExprExprVarVar(varexprs[0])), "w") == 0);

   for( i = 0; i < nvarexprs; ++i )
   {
      assert(varexprs[i] != NULL);
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(varexprs[i])), "var") == 0);
   }

   SCIPfreeBufferArray(scip, &varexprs);

   /* free allocated memory */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
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

   CHECK_TEST( testDuplicate() );

   CHECK_TEST( testTransform() );

   CHECK_TEST( testCopy() );

   CHECK_TEST( testCheck() );

   CHECK_TEST( testExp() );

   CHECK_TEST( testLog() );

   CHECK_TEST( testAbs() );

   CHECK_TEST( testHash() );

   CHECK_TEST( testCommonSubexpr() );

   CHECK_TEST( testPropagation() );

   CHECK_TEST( testGetVarExprs() );

   /* for automatic testing output the following */
   printf("SCIP Status        : all tests passed\n");
   printf("Ignore the following:\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return 0;
}
