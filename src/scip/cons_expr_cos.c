/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_cos.c
 * @brief  handler for cosine expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */  /*lint !750 */

#include <string.h>
#include <math.h>
#include "scip/cons_expr_cos.h"
#include "scip/cons_expr_sin.h"
#include "scip/cons_expr_value.h"

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
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrCos)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrCos(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

/** simplifies a cos expression
 *  Evaluates the sine value function when its child is a value expression
 *  TODO: add further simplifications
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyCos)
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
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr,
            COS(SCIPgetConsExprExprValueValue(child))) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseCos)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseConsExprExpr(scip, consexprhdlr, string, endstring, &childexpr) );
   assert(childexpr != NULL);

   /* create cosine expression */
   SCIP_CALL( SCIPcreateConsExprExprCos(scip, consexprhdlr, expr, childexpr) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the cosine expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalCos)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = COS(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffCos)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(expr != NULL);
   assert(childidx == 0);
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);

   *val = -SIN(SCIPgetConsExprExprValue(child));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalCos)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
      SCIPintervalCos(SCIP_INTERVAL_INFINITY, interval, childinterval);

   return SCIP_OKAY;
}

/** separation initialization callback */
static
SCIP_DECL_CONSEXPR_EXPRINITSEPA(initSepaCos)
{
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Bool success;
   SCIP_VAR* childvar;

   SCIP_ROWPREP* cuts[5];   /* 0: secant, 1: left tangent, 2: right tangent, 3: left mid tangent, 4: right mid tangent */
   int i;

   *infeasible = FALSE;

   childvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[0]);
   assert(childvar != NULL);

   childlb = SCIPvarGetLbLocal(childvar);
   childub = SCIPvarGetUbLocal(childvar);

   /* no need for cut if child is fixed */
   if( SCIPisRelEQ(scip, childlb, childub) )
      return SCIP_OKAY;

   /* compute underestimating cuts */
   if( underestimate )
   {
      SCIP_CALL( SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &cuts[0], &cuts[1], &cuts[2], &cuts[3], &cuts[4],
            childlb, childub, TRUE) );

      for( i = 0; i < 5; ++i )
      {
         /* only the cuts which could be created are added */
         if( !*infeasible && cuts[i] != NULL )
         {
            SCIP_CALL( SCIPcleanupRowprep(scip, cuts[i], NULL, SCIP_CONSEXPR_CUTMAXRANGE, 0.0, NULL, &success) );

            if( success && cuts[i]->nvars == 2 )
            {
               /* make a SCIP_ROW and add to LP */
               SCIP_ROW* row;

               SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, cuts[i], cons) );
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }

            SCIPfreeRowprep(scip, &cuts[i]);
         }
      }
   }

   /* compute overestimating cuts */
   if( overestimate && *infeasible )
   {
      SCIP_CALL( SCIPcomputeInitialCutsTrig(scip, conshdlr, expr, &cuts[0], &cuts[1], &cuts[2], &cuts[3], &cuts[4],
            childlb, childub, FALSE) );

      for( i = 0; i < 5; ++i )
      {
         /* only the cuts which could be created are added */
         if( !*infeasible && cuts[i] != NULL )
         {
            SCIP_CALL( SCIPcleanupRowprep(scip, cuts[i], NULL, SCIP_CONSEXPR_CUTMAXRANGE, 0.0, NULL, &success) );

            if( success && cuts[i]->nvars == 2 )
            {
               /* make a SCIP_ROW and add to LP */
               SCIP_ROW* row;

               SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, cuts[i], cons) );
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }

            SCIPfreeRowprep(scip, &cuts[i]);
         }
      }
   }

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_CONSEXPR_EXPRESTIMATE(estimateCos)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* childvar;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);
   assert(success != NULL);

   /* get expression data */
   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childvar = SCIPgetConsExprExprAuxVar(child);
   assert(childvar != NULL);

   *success = SCIPcomputeEstimatorsTrig(scip, conshdlr, expr, coefs, constant, SCIPgetSolVal(scip, sol, childvar),
                                        SCIPvarGetLbLocal(childvar), SCIPvarGetUbLocal(childvar), !overestimate);
   *islocal = TRUE;  /* TODO there are cases where cuts would be globally valid */

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_EXPRREVERSEPROP(reversepropCos)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_INTERVAL newbounds;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(nreductions != NULL);
   assert(SCIPintervalGetInf(SCIPgetConsExprExprActivity(scip, expr)) >= -1.0);
   assert(SCIPintervalGetSup(SCIPgetConsExprExprActivity(scip, expr)) <= 1.0);

   *nreductions = 0;

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   /* get the child interval and shift it to match sine */
   newbounds = SCIPgetConsExprExprActivity(scip, child);
   newbounds.inf += M_PI_2;
   newbounds.sup += M_PI_2;

   /* compute the new child interval */
   SCIP_CALL( SCIPcomputeRevPropIntervalSin(scip, SCIPgetConsExprExprActivity(scip, expr), newbounds, &newbounds) );

   /* shift the new interval back */
   newbounds.inf -= M_PI_2;
   newbounds.sup -= M_PI_2;

   /* try to tighten the bounds of the child node */
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, child, newbounds, force, reversepropqueue, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** cos hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashCos)
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
SCIP_DECL_CONSEXPR_EXPRCURVATURE(curvatureCos)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_INTERVAL childinterval;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(exprcurvature != SCIP_EXPRCURV_UNKNOWN);
   assert(childcurv != NULL);
   assert(success != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childinterval = SCIPgetConsExprExprActivity(scip, child);

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
SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(monotonicityCos)
{  /*lint --e{715}*/
   SCIP_INTERVAL interval;
   SCIP_Real inf;
   SCIP_Real sup;
   int k;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx == 0);

   assert(SCIPgetConsExprExprChildren(expr)[0] != NULL);
   interval = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(expr)[0]);

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
SCIP_RETCODE SCIPincludeConsExprExprHdlrCos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   /* include expression handler */
   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
         EXPRHDLR_PRECEDENCE, evalCos, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrCos, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, initSepaCos, NULL, estimateCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCurvature(scip, consexprhdlr, exprhdlr, curvatureCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrMonotonicity(scip, consexprhdlr, exprhdlr, monotonicityCos) );

   return SCIP_OKAY;
}

/** creates a cos expression */
SCIP_RETCODE SCIPcreateConsExprExprCos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< single child */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindConsExprExprHdlr(consexprhdlr, EXPRHDLR_NAME) != NULL);

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPfindConsExprExprHdlr(consexprhdlr, EXPRHDLR_NAME), NULL, 1,
         &child) );

   return SCIP_OKAY;
}
