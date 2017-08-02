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

/**@file   expreval.c
 * @brief  tests expression's pointwise and interval evaluation
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_pow.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;

/* creates scip, problem, includes expression constraint handler, creates  and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* currently expr constraints cannot be created */
   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1.0e-8) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* free allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

TestSuite(evalexpr, .init = setup, .fini = teardown);

/** TESTS **/
Test(evalexpr, absolute, .description = "Tests expression evaluation for absolute expressions.")
{
   int i;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_INTERVAL interval;
   const char* inputs[2] = {"abs(<x>[C]) + abs(<x>[C])",
      "abs(abs(abs(<x>[C]) + abs(<y>[C])) * abs(<x>[C])^3 * abs(<y>[C]))"};

   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)inputs[0], NULL, &expr)) );
   SCIPinfoMessage(scip, NULL, "testing expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* evaluate expression for different points */
   for( i = -10; i <= 10; ++i )
   {
      /* evaluate expression at i */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );

      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(expr),  2.0 * ABS(i)));

      /* evaluate expression at interval [-|i|, i^2]*/
      SCIP_CALL( SCIPchgVarLb(scip, x, -ABS(i)) );
      SCIP_CALL( SCIPchgVarUb(scip, x, i*i) );
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0, SCIPepsilon(scip)) );
      interval = SCIPgetConsExprExprInterval(expr);

      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), 0.0));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), 2.0 * i*i));
   }
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* test a more complicated expression |x + y| |x|^3 |y|*/
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)inputs[1], NULL, &expr)) );
   SCIPinfoMessage(scip, NULL, "testing expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* evaluate expression for different points */
   for( i = -10; i <= 10; ++i )
   {
      /* evaluate expression at i */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, (SCIP_Real) i) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(expr), 2.0 * pow(ABS(i),5)));

      /* evaluate expression at interval [-|i|, |i|]^2 */
      SCIP_CALL( SCIPchgVarLb(scip, x, -ABS(i)) );
      SCIP_CALL( SCIPchgVarUb(scip, x, ABS(i)) );
      SCIP_CALL( SCIPchgVarLb(scip, y, -ABS(i)) );
      SCIP_CALL( SCIPchgVarUb(scip, y, ABS(i)) );
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0, SCIPepsilon(scip)) );
      interval = SCIPgetConsExprExprInterval(expr);
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), 0.0));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), 2.0 * pow(ABS(i),5)));
   }
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(evalexpr, exponential, .description = "Tests expression evaluation for exponential expressions.")
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_INTERVAL interval;
   const char* inputs[2] = {"exp(<x>[C]) + exp(<x>[C])", "exp(exp(<x>[C])) * exp(<y>[C])^2"};
   int i;

   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)inputs[0], NULL, &expr)) );
   SCIPinfoMessage(scip, NULL, "testing expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* evaluate expression for different points */
   for( i = -10; i <= 10; ++i )
   {
      /* evaluate expression at i */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(expr), exp(i) + exp(i)));

      /* evaluate expression at interval [i, i + 1/(|i| + 1)] */
      SCIP_CALL( SCIPchgVarLb(scip, x, i) );
      SCIP_CALL( SCIPchgVarUb(scip, x, i + 1.0 / (ABS(i) + 1)) );
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0, SCIPepsilon(scip)) );
      interval = SCIPgetConsExprExprInterval(expr);
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), 2*exp(i)));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), 2*exp(i + 1.0 / (ABS(i) + 1))));
   }
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* complicated exponential expression e^(2y + e^x) */
   SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)inputs[1], NULL, &expr)) );
   SCIPinfoMessage(scip, NULL, "testing expression: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* evaluate expression for different points */
   for( i = 1; i <= 10; ++i )
   {
      /* evaluate expression at (1/i, i) */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) 1.0 / i) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, (SCIP_Real) i) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(expr), exp(exp(1.0 / i)) * exp(2*i)));

      /* evaluate expression at interval [-1/i, 1/i] x [i, i + 1/i] */
      SCIP_CALL( SCIPchgVarLb(scip, x, -1.0 / i) );
      SCIP_CALL( SCIPchgVarUb(scip, x,  1.0 / i) );
      SCIP_CALL( SCIPchgVarLb(scip, y, i) );
      SCIP_CALL( SCIPchgVarUb(scip, y, i + 1.0 / i) );
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0, SCIPepsilon(scip)) );
      interval = SCIPgetConsExprExprInterval(expr);
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), exp(exp(-1.0 / i)) * exp(2*i)));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), exp(exp(1.0 / i)) * exp(2*i + 2.0 / i)));
   }
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(evalexpr, logarithm, .description = "Tests expression evaluation for logarithmic expressions.")
{
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_INTERVAL interval;
      const char* inputs[2] = {"log(<x>[C]) + log(<x>[C])", "log(log(exp(<x>[C]) * exp(<y>[C])))"};
      int i;

      SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)inputs[0], NULL, &expr) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = -10; i <= 10; ++i )
      {
         SCIP_Real xlb, xub;

         /* evaluate expression at i */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
         SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );

         if( i <= 0 )
            cr_expect_eq(SCIPgetConsExprExprValue(expr), SCIP_INVALID);
         else
            cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(expr), log(i) + log(i)));

         /* evaluate expression at interval [i, i + 1/(|i|+1)] */
         xlb = i;
         xub = i + 1.0 / (ABS(i) + 1);
         SCIP_CALL( SCIPchgVarLb(scip, x, xlb) );
         SCIP_CALL( SCIPchgVarUb(scip, x, xub) );
         SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0, SCIPepsilon(scip)) );
         interval = SCIPgetConsExprExprInterval(expr);

         /* interval is empty if both bounds are non-positive */
         if( xub <= 0 )
            cr_expect(SCIPintervalIsEmpty(SCIPinfinity(scip), interval));
         else
         {
            cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), 2*log(xub)));

            if( xlb <= 0 )
               cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), -SCIPinfinity(scip)));
            else
               cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), 2*log(xlb)));
         }
      }
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

      /* complicated logarithmic expression: log( log( exp(x) * exp(y) ) ) = log( x + y ) */
      SCIP_CALL( (SCIPparseConsExprExpr(scip, conshdlr, (char*)inputs[1], NULL, &expr)) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = 1; i <= 10; ++i )
      {
         /* evaluate expression at (i, i + 1) */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
         SCIP_CALL( SCIPsetSolVal(scip, sol, y, (SCIP_Real) i + 1) );
         SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
         cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(expr), log(2*i + 1) ));

         /* evaluate expression at intervar [1/i, 2/i] x [3/i, 4/i] */
         SCIP_CALL( SCIPchgVarLb(scip, x,  1.0 / i) );
         SCIP_CALL( SCIPchgVarUb(scip, x,  2.0 / i) );
         SCIP_CALL( SCIPchgVarLb(scip, y,  3.0 / i) );
         SCIP_CALL( SCIPchgVarUb(scip, y,  4.0 / i) );
         SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0, SCIPepsilon(scip)) );
         interval = SCIPgetConsExprExprInterval(expr);
         cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), log(1.0 / i + 3.0 / i)));
         cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), log(2.0 / i + 4.0 / i)));
      }

      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(evalexpr, power, .description = "Tests expression evaluation for power expressions.")
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_INTERVAL interval;

   SCIP_CALL( SCIPchgVarLb(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 10.0) );

   /* create expressions x^2.5 */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &expr, xexpr, 2.5) );

   /* evaluate expression at 2.0 */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );
   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
   cr_assert(SCIPisEQ(scip, SCIPgetConsExprExprValue(expr), pow(2.0, 2.5)));

   /* evaluate expression for an undefined point: -1.0 */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -1.0) );
   SCIP_CALL( SCIPevalConsExprExpr(scip, expr, sol, 0) );
   cr_assert(SCIPgetConsExprExprValue(expr) == SCIP_INVALID);

   /* evaluate expression at interval [1.0, 3.0] */
   SCIP_CALL( SCIPchgVarLb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0, 0.0) );
   interval = SCIPgetConsExprExprInterval(expr);
   cr_assert(SCIPisEQ(scip, interval.inf, 1.0));
   cr_assert(SCIPisEQ(scip, interval.sup, pow(3.0, 2.5)));

   /* evaluate expression at an undefined interval [-2.0, -1.0] => resulting interval should be empty */
   SCIP_CALL( SCIPchgVarLb(scip, x, -2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, expr, FALSE, 0, 0.0) );
   interval = SCIPgetConsExprExprInterval(expr);
   cr_assert(SCIPintervalIsEmpty(SCIPinfinity(scip), interval));

   /* free expressions */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
}

/* creates expression for f(x,y) = 0.5 * ( (x^2*y^(-1)*5^(-4))^2 * (2*x + 1)^(-1) ) */
static
SCIP_RETCODE createExpr(
   SCIP_CONSEXPR_EXPR**  xexpr,              /**< pointer to store variable expression */
   SCIP_CONSEXPR_EXPR**  yexpr,              /**< pointer to store variable expression */
   SCIP_CONSEXPR_EXPR**  const_expr,         /**< pointer to store constant expression */
   SCIP_CONSEXPR_EXPR**  prodexpr,           /**< pointer to store product expression */
   SCIP_CONSEXPR_EXPR**  sumexpr,            /**< pointer to store sum expression */
   SCIP_CONSEXPR_EXPR**  mainexpr            /**< pointer to store full expression */
   )
{
   SCIP_CONSEXPR_EXPR* exprs[] = {NULL, NULL, NULL};
   SCIP_Real coef = 2.0;

   /* create variable and constant expressions */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, const_expr, 5.0) );

   /* create sum expression */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, sumexpr, 1, xexpr, &coef, 1.0) );  /* 2*x+1 */

   /* create product expression */
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &exprs[0], *xexpr,  2.0) );  /* x^2 */
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &exprs[1], *yexpr, -1.0) );  /* y^(-1) */
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &exprs[2], *const_expr, -4.0) );  /* 5^(-4) */
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, prodexpr, 3, exprs, 1.0) );  /* x^2*y^(-1)*5^(-4) */

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[0]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[2]) );

   /* create main expression */
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &exprs[0], *prodexpr,  2.0) );  /* (x^2*y^(-1)*5^(-4))^2 */
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &exprs[1], *sumexpr, -1.0) );  /* (2*x + 1)^(-1) */
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, mainexpr, 2, exprs, 0.5) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[0]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprs[1]) );

   return SCIP_OKAY;
}

/* helper function to check evaluation of expression created with createExpr() */
static
void checkExprEval(
   SCIP_CONSEXPR_EXPR*  xexpr,               /**< variable expression */
   SCIP_CONSEXPR_EXPR*  yexpr,               /**< variable expression */
   SCIP_CONSEXPR_EXPR*  const_expr,          /**< constant expression */
   SCIP_CONSEXPR_EXPR*  prodexpr,            /**< product expression */
   SCIP_CONSEXPR_EXPR*  sumexpr,             /**< sum expression */
   SCIP_CONSEXPR_EXPR*  mainexpr,            /**< full expression */
   SCIP_Real            xval,                /**< x value used for evaluation */
   SCIP_Real            yval,                /**< y value used for evaluation */
   unsigned int         tag                  /**< tag used for evaluation */
   )
{
   SCIP_Real prodval;
   SCIP_Real sumval;

   prodval = pow(xval,2)*pow(yval,-1)*pow(5,-4);
   sumval = 2*xval + 1;

   /* check values */
   cr_expect_float_eq(SCIPgetConsExprExprValue(mainexpr), 0.5 * pow(prodval,2) * pow(sumval,-1), SCIPepsilon(scip));
   cr_expect_float_eq(SCIPgetConsExprExprValue(sumexpr), sumval, SCIPepsilon(scip));
   cr_expect_float_eq(SCIPgetConsExprExprValue(prodexpr), prodval, SCIPepsilon(scip));
   cr_expect_float_eq(SCIPgetConsExprExprValue(xexpr), xval, SCIPepsilon(scip));
   cr_expect_float_eq(SCIPgetConsExprExprValue(yexpr), yval, SCIPepsilon(scip));
   cr_expect_float_eq(SCIPgetConsExprExprValue(const_expr), 5.0, SCIPepsilon(scip));

   /* check tags */
   cr_expect_eq(SCIPgetConsExprExprEvalTag(mainexpr), tag);
   cr_expect_eq(SCIPgetConsExprExprEvalTag(sumexpr), tag);
   cr_expect_eq(SCIPgetConsExprExprEvalTag(prodexpr), tag);
   cr_expect_eq(SCIPgetConsExprExprEvalTag(xexpr), tag);
   cr_expect_eq(SCIPgetConsExprExprEvalTag(yexpr), tag);
   cr_expect_eq(SCIPgetConsExprExprEvalTag(const_expr), tag);
}

Test(evalexpr, complicated, .description = "Tests expression evaluation for a large complicated expression.")
{
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_CONSEXPR_EXPR* yexpr;
   SCIP_CONSEXPR_EXPR* const_expr;
   SCIP_CONSEXPR_EXPR* prodexpr;
   SCIP_CONSEXPR_EXPR* sumexpr;
   SCIP_CONSEXPR_EXPR* mainexpr;
   SCIP_CONSEXPR_PRINTDOTDATA* dotdata;
   int i;

   /* create all expressions */
   SCIP_CALL( createExpr(&xexpr, &yexpr, &const_expr, &prodexpr, &sumexpr, &mainexpr) );

   /* initialize solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 4.0) );

   /* evaluate main expression, print it, and check values for sub-expressions */
   printf("evaluate and check expressions\n");
   SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 1) );
   SCIP_CALL( SCIPprintConsExprExprDotInit(scip, &dotdata, NULL, SCIP_CONSEXPR_PRINTDOT_EXPRSTRING | SCIP_CONSEXPR_PRINTDOT_EVALTAG) );
   SCIP_CALL( SCIPprintConsExprExprDot(scip, dotdata, mainexpr) );
   SCIP_CALL( SCIPprintConsExprExprDotFinal(scip, &dotdata) );
   checkExprEval(xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, 2.0, 4.0, 1);

   /* modify solution and evaluate expression with the same tag again; values should not change */
   printf("evaluate and check expressions with a modified solution but the same tag\n");
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -5.0) );
   SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 1) );
   checkExprEval(xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, 2.0, 4.0, 1);

   /* evaluate expression with a different tag; values should have changed */
   printf("evaluate expression with a new tag\n");
   SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 2) );
   checkExprEval(xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, -2.0, -5.0, 2);

   /* evaluate solution with zero tag */
   printf("evaluate expression with a zero tag\n");
   for( i = 1; i < 100; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, i*i) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, -5.0/i) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 0) );
      checkExprEval(xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, i*i, -5.0 / i, 0);
   }

   /* mainexpr is not defined for x = -1 or y = 0; the result should be SCIP_INVALID */
   printf("evaluate expression for an undefined point\n");
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 0.0) );
   SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprValue(mainexpr), SCIP_INVALID);

   /* set values for variable expression explicitly */
   printf("evaluate expression after setting value for variable expressions\n");
   for( i = 1; i < 3; ++i )
   {
      SCIPsetConsExprExprEvalValue(xexpr, i*i, i);
      SCIPsetConsExprExprEvalValue(yexpr, 1.0 / i, i);
      SCIP_CALL( SCIPevalConsExprExpr(scip, mainexpr, NULL, i) );
      checkExprEval(xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, i*i, 1.0 / i, i);
   }

   /* release all expressions */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &const_expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &mainexpr) );
}

/* helper function to check interval evaluation of an expression */
static
void checkExprIntEval(
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
   cr_expect_eq(targettag, SCIPgetConsExprExprEvalIntervalTag(expr));

   /* check if interval is and should be empty */
   cr_expect_eq(empty, SCIPintervalIsEmpty(SCIPinfinity(scip), interval));

   /* check interval */
   if( !empty )
   {
      cr_expect_float_eq(SCIPintervalGetInf(interval), targetinf, SCIPepsilon(scip));
      cr_expect_float_eq(SCIPintervalGetSup(interval), targetsup, SCIPepsilon(scip));
   }
}

/* Test expression interval evaluation method */
Test(evalexprInterval, complicated_interval, .description = "Tests expression interval evaluation for a large complicated expression.")
{
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

   /* currently expr constraints cannot be created */
   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -5.0, 5.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create different kind of expressions */
   exponent = 0.5;
   SCIP_CALL( createExpr(&xexpr, &yexpr, &const_expr, &prodexpr, &sumexpr, &mainexpr) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &sqrtexpr, xexpr, exponent) );

   /*
    * check interval evaluation method for constant expressions
    */
   printf("check interval-evaluation of constant expressions\n");
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, const_expr, FALSE, 0, 0.0) );
   checkExprIntEval(const_expr, 5.0, 5.0, FALSE, 0);

   /*
    * check interval evaluation method for variable expressions
    */
   printf("check interval-evaluation of variable expressions\n");
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, xexpr, FALSE, 0, 0.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, yexpr, FALSE, 0, 0.0) );
   checkExprIntEval(xexpr, 0.0, 10.0, FALSE, 0);
   checkExprIntEval(yexpr, -5.0, 5.0, FALSE, 0);

   /*
    * check interval evaluation method for sum expression
    */
   printf("check interval-evaluation of sum expression\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 7.5) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sumexpr, FALSE, 0, 0.0) );
   checkExprIntEval(sumexpr, 5.0, 16.0, FALSE, 0);
   checkExprIntEval(xexpr, 2.0, 7.5, FALSE, 0);

   /*
    * check interval evaluation method for product expression: (x^2 / (y*5^(-4)))
    */
   printf("check interval-evaluation of product expression\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.5) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 2.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 0, 0.0) );
   checkExprIntEval(prodexpr, pow(5,-4.0) / 8.0 , pow(5,-4.0), FALSE, 0);
   checkExprIntEval(xexpr, 0.5, 1.0, FALSE, 0);
   checkExprIntEval(yexpr, 1.0, 2.0, FALSE, 0);

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

               SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0, 0.0) );

               /* check all expressions */
               checkExprIntEval(mainexpr, MIN(inf,sup), MAX(inf,sup), FALSE, 0);
               checkExprIntEval(sumexpr, 2*xlb + 1.0, 2*xub + 1.0, FALSE, 0);
               inf = MIN(xlb*xlb, xub*xub) * (1.0/yub) * pow(5,-4);
               sup = MAX(xlb*xlb, xub*xub) * (1.0/ylb) * pow(5,-4);
               checkExprIntEval(prodexpr, inf, sup, FALSE, 0);
               checkExprIntEval(xexpr, xlb, xub, FALSE, 0);
               checkExprIntEval(yexpr, ylb, yub, FALSE, 0);
               checkExprIntEval(const_expr, 5.0, 5.0, FALSE, 0);
            }

   /*
    * check if interval evaluation for 1/(1+2*x)^3 or 1/y leads to infinite bounds
    */
   printf("check interval-evaluation for expressions containing functions like 1/f(x)\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, -0.5) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -0.5) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0, 0.0) );
   checkExprIntEval(mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0);

   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0, 0.0) );
   checkExprIntEval(mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0);

   SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 0.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0, 0.0) );
   checkExprIntEval(mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0);

   /* (1/y)^2 should lead to [0,inf] */
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0, 0.0) );
   checkExprIntEval(mainexpr, 0.0, SCIPinfinity(scip), FALSE, 0);

   /* (1/y)^2 should lead to [0,inf] but because of 1/(1+2*x)^3 we should get [-inf,inf] */
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 10.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 0, 0.0) );
   checkExprIntEval(mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, 0);

   /*
    * check if interval evaluation aborts for some cases like (-1)^2
    */
   printf("check interval-evaluation for undefined expressions like (-1)^2\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, FALSE, 0, 0.0) );
   checkExprIntEval(sqrtexpr, 0, 0, TRUE, 0);

   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -0.5) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, FALSE, 0, 0.0) );
   checkExprIntEval(sqrtexpr, 0, 0, TRUE, 0);

   /* the result of sqrt([-1,2]) should be [0,sqrt(2)] and not an empty interval */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, FALSE, 0, 0.0) );
   checkExprIntEval(sqrtexpr, 0, sqrt(2), FALSE, 0);

   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sqrtexpr, FALSE, 0, 0.0) );
   checkExprIntEval(sqrtexpr, 0.0, 1.0, FALSE, 0);

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
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, tag, 0.0) );
      checkExprIntEval(mainexpr, 0.0, 0.0, FALSE, tag);
      checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, tag);
      checkExprIntEval(sumexpr, 1.0, 1.0, FALSE, tag);
      checkExprIntEval(xexpr, 0.0, 0.0, FALSE, tag);
      checkExprIntEval(yexpr, 1.0, 1.0, FALSE, tag);
      checkExprIntEval(const_expr, 5.0, 5.0, FALSE, tag);
   }

   /* set another tag for some subexpression; result should not change */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 1, 0.0) );
   checkExprIntEval(mainexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval(sumexpr, 1.0, 1.0, FALSE, 1);

   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 2, 0.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, sumexpr, FALSE, 2, 0.0) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 1, 0.0) ); /* this should not trigger a reevaluation */
   checkExprIntEval( mainexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval( prodexpr, 0.0, 0.0, FALSE, 2);
   checkExprIntEval( sumexpr, 1.0, 1.0, FALSE, 2);

   SCIP_CALL( SCIPevalConsExprExprInterval(scip, mainexpr, FALSE, 3, 0.0) ); /* this should trigger a reevaluation */
   checkExprIntEval(mainexpr, 0.0, 0.0, FALSE, 3);
   checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, 3);
   checkExprIntEval(sumexpr, 1.0, 1.0, FALSE, 3);

   /* manipulate evaluation interval */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 1, 0.0) );
   checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval(xexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval(yexpr, 1.0, 1.0, FALSE, 1);

   /* set bounds of x to [1,1] in xexpr */
   SCIPintervalSetBounds(&interval, 1.0, 1.0);
   SCIPsetConsExprExprEvalInterval(xexpr, &interval, 2);
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 2, 0.0) ); /* should use [1,1] for xexpr */
   checkExprIntEval(prodexpr, pow(5.0,-4), pow(5.0,-4), FALSE, 2);
   checkExprIntEval(xexpr, 1.0, 1.0, FALSE, 2);
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, prodexpr, FALSE, 3, 0.0) ); /* should use [0,0] for xexpr */
   checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, 3);
   checkExprIntEval(xexpr, 0.0, 0.0, FALSE, 3);

   /* release all expressions */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &const_expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &mainexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sqrtexpr) );
}
