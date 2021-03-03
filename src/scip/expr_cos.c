/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   expr_cos.c
 * @brief  handler for cosine expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */  /*lint !750 */

#include <string.h>
#include <math.h>
#include "scip/expr_cos.h"
#include "scip/expr_sin.h"
#include "scip/expr_value.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2         1.57079632679489661923
#endif
#ifndef M_PI_4
#define M_PI_4         0.785398163397448309616
#endif

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "cos"
#define EXPRHDLR_DESC         "cosine expression"
#define EXPRHDLR_PRECEDENCE   92000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(82463.0)

/*
 * Local methods
 */

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrCos)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprHdlrCos(scip) );

   return SCIP_OKAY;
}

/** simplifies a cos expression
 *  Evaluates the sine value function when its child is a value expression
 *  TODO: add further simplifications
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyCos)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /* check for value expression */
   if( SCIPisExprValue(scip, child) )
   {
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, cos(SCIPgetValueExprValue(child)), ownercreate,
               ownercreatedata) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_EXPRPARSE(parseCos)
{  /*lint --e{715}*/
   SCIP_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseExpr(scip, &childexpr, string, endstring, ownercreate, ownercreatedata) );
   assert(childexpr != NULL);

   /* create cosine expression */
   SCIP_CALL( SCIPcreateExprCos(scip, expr, childexpr, ownercreate, ownercreatedata) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the cosine expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_EXPREVAL(evalCos)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = cos(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffCos)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(expr != NULL);
   assert(childidx == 0);
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(child)), "val") != 0);

   *val = -sin(SCIPexprGetEvalValue(child));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalCos)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
      SCIPintervalCos(SCIP_INTERVAL_INFINITY, interval, childinterval);

   return SCIP_OKAY;
}

/** separation initialization callback */
static
SCIP_DECL_EXPRINITESTIMATES(initEstimateCos)
{
   SCIP_Real childlb;
   SCIP_Real childub;

   childlb = bounds[0].inf;
   childub = bounds[0].sup;

   /* no need for cut if child is fixed */
   if( SCIPisRelEQ(scip, childlb, childub) )
      return SCIP_OKAY;

   /* compute cuts */
   SCIP_CALL( SCIPcomputeInitialCutsTrig(scip, expr, childlb, childub, ! overestimate, coefs, constant, nreturned) );

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_EXPRESTIMATE(estimateCos)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);
   assert(success != NULL);


   *success = SCIPcomputeEstimatorsTrig(scip, expr, coefs, constant, refpoint[0], localbounds[0].inf,
         localbounds[0].sup, ! overestimate);
   *islocal = TRUE;  /* TODO there are cases where cuts would be globally valid */

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropCos)
{  /*lint --e{715}*/
   SCIP_INTERVAL newbounds;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   /* bounds should have been intersected with activity, which is within [-1,1] */
   assert(SCIPintervalGetInf(bounds) >= -1.0);
   assert(SCIPintervalGetSup(bounds) <= 1.0);

   /* get the child interval */
   newbounds = childrenbounds[0];

   /* shift child interval to match sine */
   SCIPintervalAddScalar(SCIP_INTERVAL_INFINITY, &newbounds, newbounds, M_PI_2);  /* TODO use bounds on Pi/2 instead of approximation of Pi/2 */

   /* compute the new child interval */
   SCIP_CALL( SCIPcomputeRevPropIntervalSin(scip, bounds, newbounds, &newbounds) );

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newbounds) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* shift the new interval back */
   SCIPintervalAddScalar(SCIP_INTERVAL_INFINITY, &childrenbounds[0], newbounds, -M_PI_2);  /* TODO use bounds on Pi/2 instead of approximation of Pi/2 */

   return SCIP_OKAY;
}

/** cos hash callback */
static
SCIP_DECL_EXPRHASH(hashCos)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= childrenhashes[0];

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureCos)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_INTERVAL childinterval;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(exprcurvature != SCIP_EXPRCURV_UNKNOWN);
   assert(childcurv != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   SCIP_CALL( SCIPevalExprActivity(scip, child) );
   childinterval = SCIPexprGetActivity(child);

   /* TODO rewrite SCIPcomputeCurvatureSin so it provides the reverse operation */
   *success = TRUE;
   if( SCIPcomputeCurvatureSin(SCIP_EXPRCURV_CONCAVE, childinterval.inf + M_PI_2, childinterval.sup + M_PI_2) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_CONCAVE;
   else if( SCIPcomputeCurvatureSin(SCIP_EXPRCURV_CONVEX, childinterval.inf + M_PI_2, childinterval.sup + M_PI_2) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_CONVEX;
   else if( SCIPcomputeCurvatureSin(SCIP_EXPRCURV_LINEAR, childinterval.inf + M_PI_2, childinterval.sup + M_PI_2) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_LINEAR;
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityCos)
{  /*lint --e{715}*/
   SCIP_INTERVAL interval;
   SCIP_Real inf;
   SCIP_Real sup;
   int k;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx == 0);

   assert(SCIPexprGetChildren(expr)[0] != NULL);
   SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(expr)[0]) );
   interval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   *result = SCIP_MONOTONE_UNKNOWN;
   inf = SCIPintervalGetInf(interval);
   sup = SCIPintervalGetSup(interval);

   /* expression is not monotone because the interval is too large */
   if( sup - inf > M_PI )
      return SCIP_OKAY;

   /* compute k s.t. PI * k <= interval.inf <= PI * (k+1) */
   k = (int)floor(inf/M_PI);
   assert(M_PI * k <= inf);
   assert(M_PI * (k+1) >= inf);

   /* check whether [inf,sup] are contained in an interval for which the cosine function is monotone */
   if( sup <= M_PI * (k+1) )
      *result = ((k % 2 + 2) % 2) == 0 ? SCIP_MONOTONE_DEC : SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** creates the handler for cos expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeExprHdlrCos(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   /* include expression handler */
   SCIP_CALL( SCIPincludeExprHdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE, evalCos, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrCos, NULL);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyCos);
   SCIPexprhdlrSetParse(exprhdlr, parseCos);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalCos);
   SCIPexprhdlrSetEstimate(exprhdlr, initEstimateCos, estimateCos);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropCos);
   SCIPexprhdlrSetHash(exprhdlr, hashCos);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffCos, NULL, NULL);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureCos);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityCos);

   return SCIP_OKAY;
}

/** creates a cos expression */
SCIP_RETCODE SCIPcreateExprCos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindExprHdlr(scip, EXPRHDLR_NAME) != NULL);

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPfindExprHdlr(scip, EXPRHDLR_NAME), NULL, 1, &child, ownercreate,
            ownercreatedata) );

   return SCIP_OKAY;
}
