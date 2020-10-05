/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_erf.c
 * @brief  handler for gaussian error function expressions
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_erf.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr.h"

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "erf"
#define EXPRHDLR_DESC         "gaussian error function"
#define EXPRHDLR_PRECEDENCE   79000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(131071.0)

/*
 * Data structures
 */

/*
 * Local methods
 */

/** evaluates the gaussian error function at a given point */
static
SCIP_Real errorf(
   SCIP_Real             x                   /**< point to evaluate */
   )
{
   SCIP_Real a1 = +0.254829592;
   SCIP_Real a2 = -0.284496736;
   SCIP_Real a3 = +1.421413741;
   SCIP_Real a4 = -1.453152027;
   SCIP_Real a5 = +1.061405429;
   SCIP_Real p  = +0.3275911;
   int sign  = (x >= 0.0) ? 1 : -1;
   SCIP_Real t = 1.0 / (1.0 + p * REALABS(x));
   SCIP_Real y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x*x);

   return sign*y;
}

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrErf)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrErf(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

/** simplifies an erf expression */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyErf)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   /* check for value expression */
   if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrValue(conshdlr) )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, errorf(SCIPgetConsExprExprValueValue(child))) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

/** expression compare callback */
static
SCIP_DECL_CONSEXPR_EXPRCOMPARE(compareErf)
{  /*lint --e{715}*/
   assert(expr1 != NULL);
   assert(expr2 != NULL);

   SCIPerrorMessage("method of erf expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return 0;
}

/** expression parse callback */
static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseErf)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseConsExprExpr(scip, consexprhdlr, string, endstring, &childexpr) );
   assert(childexpr != NULL);

   /* create gaussian error function expression */
   SCIP_CALL( SCIPcreateConsExprExprErf(scip, consexprhdlr, expr, childexpr) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the gaussian error function expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalErf)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = errorf(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalErf)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
   {
      SCIP_Real childinf = SCIPintervalGetInf(childinterval);
      SCIP_Real childsup = SCIPintervalGetSup(childinterval);
      SCIP_Real inf = childinf <= -SCIP_INTERVAL_INFINITY ? -1.0 : errorf(childinf);
      SCIP_Real sup = childsup >= +SCIP_INTERVAL_INFINITY ? +1.0 : errorf(childsup);
      assert(inf <= sup);
      SCIPintervalSetBounds(interval, inf, sup);
   }

   return SCIP_OKAY;
}

/** expression under/overestimation callback */
static
SCIP_DECL_CONSEXPR_EXPRESTIMATE(estimateErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_EXPRREVERSEPROP(reversepropErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** erf hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= childrenhashes[0];

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_CONSEXPR_EXPRCURVATURE(curvatureErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   /* expression is
    *  - convex if child is convex and child <= 0
    *  - concave if child is concave and child >= 0
    */
   if( exprcurvature == SCIP_EXPRCURV_CONVEX )
   {
      *success = TRUE;
      *childcurv = SCIP_EXPRCURV_CONVEX;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(monotonicityErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);

   *result = SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEGRALITY(integralityErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   *isintegral = FALSE;

   return SCIP_OKAY;
}

/** creates the handler for erf expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrErf(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   /* include expression handler */
   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
         EXPRHDLR_PRECEDENCE, evalErf, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrErf, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, NULL, NULL, estimateErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCurvature(scip, consexprhdlr, exprhdlr, curvatureErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrMonotonicity(scip, consexprhdlr, exprhdlr, monotonicityErf) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntegrality(scip, consexprhdlr, exprhdlr, integralityErf) );

   return SCIP_OKAY;
}

/** creates a erf expression */
SCIP_RETCODE SCIPcreateConsExprExprErf(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< childr expression */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(consexprhdlr != NULL);
   assert(expr != NULL);

   exprhdlr = SCIPfindConsExprExprHdlr(consexprhdlr, EXPRHDLR_NAME);

   if( exprhdlr != NULL )
   {
      SCIPerrorMessage("could not find %s expression handler -> abort\n", EXPRHDLR_NAME);
      SCIPABORT();
      return SCIP_ERROR;
   }

   /* create expression data */
   exprdata = NULL;
   /* create expression */
   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, 1, &child) );

   return SCIP_OKAY;
}
