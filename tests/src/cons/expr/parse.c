/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   parse.c
 * @brief  tests expression parsing
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_sumprod.h"
#include "scip/cons_expr_var.h"
#include "scip/struct_cons_expr.h"

static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(check_nuses)
{
   int expectednuses;

   assert(expr != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   /* root is captured by the walker!!! */
   expectednuses = 1;
   if( SCIPgetConsExprExprWalkParent(expr) == NULL )
      expectednuses = 2;

   if( expr->nuses != expectednuses )
   {
      printf("following expression is captured too many times (%d, expected %d)\n",
            expr->nuses, expectednuses);
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
      assert(expr->nuses == expectednuses);
   }

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

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

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z(  ", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

static
void teardown(void)
{

   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(parse, .init = setup, .fini = teardown);

Test(parse, simple)
{
   SCIP_CONSEXPR_EXPR* expr_xy5;
   const char* input = "<x>[C] / <y>[I] *(-5)";

   /* create expression for product of -5, x, and y */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr_xy5), SCIP_OKAY);

   /* print expression */
   SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr_xy5, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* check that the expression is capture correctly */
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr_xy5, check_nuses, NULL, NULL, NULL, NULL) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr_xy5) );
}

Test(parse, simple2)
{
   SCIP_CONSEXPR_EXPR* crazyexpr;
   const char* input = "-<x>[C] * <y>[I] ^(-1) + (<x>[C]+<y>[C])^2";

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &crazyexpr), SCIP_OKAY);

   /* print expression */
   SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
   SCIP_CALL( SCIPprintConsExprExpr(scip, crazyexpr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &crazyexpr) );
}

Test(parse, eval)
{
   SCIP_CONSEXPR_EXPR* crazyexpr;
   SCIP_SOL* crazysol;
   const char* input = "<x>*<y>^2/<x>^4 - 2*<x>*(3+5*<x>-2*<y>)*(<x>+<y>)^(-3.5)";
   const SCIP_Real vals[3][2] = { {1.0, 2.0}, {0.0, 0.0}, {42.0, 17.0} };
   int p;

#define CRAZYEVAL(x, y) ((x)*pow(y,2)/pow(x,4) - 2*(x)*(3+5*(x)-2*(y)) * pow((x)+(y), -3.5))

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &crazyexpr), SCIP_OKAY);

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
         cr_expect(!SCIPisFinite(expvalue));
      else
         cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprValue(crazyexpr), expvalue));
   }

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &crazyexpr) );

   /* release solution */
   SCIP_CALL( SCIPfreeSol(scip, &crazysol) );
}

Test(parse, unusual_var_name)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "(<x> - <y>) /   <z(  >^2";

   /* parse */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   /* print expression */
   SCIPinfoMessage(scip, NULL, "printing expression %s after parsing from string: ", input);
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* check that the expression is capture correctly */
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, check_nuses, NULL, NULL, NULL, NULL) );
   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(parse, constraint_with_spaces)
{
      SCIP_CONS* consexpr_xy5;
      SCIP_Bool success;
      const char* input = "[expr] <test1>: <x>[C] / <y>[I] *(5) >= 1;";

      /* parse constraint */
      success = FALSE;
      SCIP_CALL( SCIPparseCons(scip, &consexpr_xy5, input,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      cr_expect(success);

      /* print constraint */
      SCIPinfoMessage(scip, NULL, "printing constraint %s after parsing from string:", input);
      SCIP_CALL( SCIPprintCons(scip, consexpr_xy5, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &consexpr_xy5) );
}

Test(parse, constraint_with_sides)
{
      SCIP_CONS* cons;
      SCIP_Bool success;
      const char* input = "[expr] <test2>: 1 <= <x>[C] / <y>[I] *(5) - <x> <= 2;";

      /* parse constraint */
      success = FALSE;
      SCIP_CALL( SCIPparseCons(scip, &cons, input,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      cr_expect(success);

      /* print constraint */
      SCIPinfoMessage(scip, NULL, "printing constraint %s after parsing from string:", input);
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(parse, invalid_expressions)
{
   SCIP_CONSEXPR_EXPR* e;

   SCIPmessageSetErrorPrinting(NULL, NULL);

   /* there is no variable with name "xx" */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<xx>", NULL, &e), SCIP_READERROR);
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"5/<donothave> ", NULL, &e), SCIP_READERROR);
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> +-*5 ", NULL, &e), SCIP_READERROR);
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> / (<y>-5 ", NULL, &e), SCIP_READERROR);
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"donothave(<x>) ", NULL, &e), SCIP_READERROR);
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"donothave(<x> ", NULL, &e), SCIP_READERROR);
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"val(1) ", NULL, &e), SCIP_READERROR);

#ifdef FAILING_TESTS
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x>+2<y> ", NULL, &e), SCIP_READERROR);
#endif

   SCIPmessageSetErrorPrintingDefault();
}

Test(parse, misc)
{
   SCIP_CONSEXPR_EXPR* e;

   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"-5+3*<x>", NULL, &e), SCIP_OKAY);
   cr_expect_eq(SCIPgetConsExprExprHdlr(e), SCIPgetConsExprExprHdlrSum(conshdlr));
   cr_expect_eq(SCIPgetConsExprExprSumConstant(e), -5.0);
   cr_expect_eq(SCIPgetConsExprExprSumCoefs(e)[0], 3.0);
   cr_expect_eq(SCIPgetConsExprExprNChildren(e), 1);
   cr_expect_eq(SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(e)[0]), SCIPgetConsExprExprHdlrVar(conshdlr));
   cr_expect_eq(SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(e)[0]), x);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );

   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x>", NULL, &e), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> +5*<y>", NULL, &e), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> + <x>*5*<y>", NULL, &e), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> +5", NULL, &e), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x> +5 + <y>", NULL, &e), SCIP_OKAY);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &e) );
}
