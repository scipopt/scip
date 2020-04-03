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

/**@file   cons_expr_entropy.c
 * @brief  handler for -x*log(x) expressions
 * @author Benjamin Mueller
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_entropy.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr.h"

#include <string.h>

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "entropy"
#define EXPRHDLR_DESC         "expression handler for -x*log(x)"
#define EXPRHDLR_PRECEDENCE   0
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(7477.0)

/*
 * Data structures
 */

/*
 * Local methods
 */

/** helper function for reverseProp() which returns an x* in [xmin,xmax] s.t. the distance -x*log(x) and a given target
 *  value is minimized; the function assumes that -x*log(x) is monotone on [xmin,xmax];
 */
static
SCIP_Real reversePropBinarySearch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             xmin,               /**< smallest possible x */
   SCIP_Real             xmax,               /**< largest possible x */
   SCIP_Bool             increasing,         /**< -x*log(x) is increasing or decreasing on [xmin,xmax] */
   SCIP_Real             targetval           /**< target value */
   )
{
   SCIP_Real xminval = (xmin == 0.0) ? 0.0 : -xmin * log(xmin);
   SCIP_Real xmaxval = (xmax == 0.0) ? 0.0 : -xmax * log(xmax);
   int i;

   assert(xmin <= xmax);
   assert(increasing ? xminval <= xmaxval : xminval >= xmaxval);

   /* function can not achieve -x*log(x) -> return xmin or xmax */
   if( SCIPisGE(scip, xminval, targetval) && SCIPisGE(scip, xmaxval, targetval) )
      return increasing ? xmin : xmax;
   else if( SCIPisLE(scip, xminval, targetval) && SCIPisLE(scip, xmaxval, targetval) )
      return increasing ? xmax : xmin;

   /* binary search */
   for( i = 0; i < 1000; ++i )
   {
      SCIP_Real x = (xmin + xmax) / 2.0;
      SCIP_Real xval = (x == 0.0) ? 0.0 : -x * log(x);

      /* found the corresponding point -> skip */
      if( SCIPisEQ(scip, xval, targetval) )
         return x;
      else if( SCIPisLT(scip, xval, targetval) )
      {
         if( increasing )
            xmin = x;
         else
            xmax = x;
      }
      else
      {
         if( increasing )
            xmax = x;
         else
            xmin = x;
      }
   }

   return SCIP_INVALID;
}

/** helper function for reverse propagation; needed for proper unittest */
static
SCIP_RETCODE reverseProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTERVAL         exprinterval,       /**< bounds on the expression */
   SCIP_INTERVAL         childinterval,      /**< bounds on the interval of the child */
   SCIP_INTERVAL*        interval            /**< resulting interval */
   )
{
   SCIP_INTERVAL childentropy;
   SCIP_INTERVAL intersection;
   SCIP_INTERVAL tmp;
   SCIP_Real childinf;
   SCIP_Real childsup;
   SCIP_Real extremum;
   SCIP_Real boundinf;
   SCIP_Real boundsup;

   assert(scip != NULL);
   assert(interval != NULL);

   /* check whether domain is empty, i.e., bounds on -x*log(x) > 1/e */
   if( SCIPisGT(scip, SCIPintervalGetInf(exprinterval), exp(-1.0))
      || SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
   {
      SCIPintervalSetEmpty(interval);
      return SCIP_OKAY;
   }

   /* compute the intersection between entropy([childinf,childsup]) and [expr.inf, expr.sup] */
   SCIPintervalEntropy(SCIP_INTERVAL_INFINITY, &childentropy, childinterval);
   SCIPintervalIntersect(&intersection, childentropy, exprinterval);

   /* intersection empty -> infeasible */
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, intersection) )
   {
      SCIPintervalSetEmpty(interval);
      return SCIP_OKAY;
   }

   /* intersection = childentropy -> nothing can be learned */
   if( SCIPintervalIsSubsetEQ(SCIP_INTERVAL_INFINITY, childentropy, intersection) )
   {
      *interval = childinterval;
      return SCIP_OKAY;
   }

   childinf = MAX(0.0, SCIPintervalGetInf(childinterval)); /*lint !e666*/
   childsup = SCIPintervalGetSup(childinterval);
   extremum = exp(-1.0);
   boundinf = SCIP_INVALID;
   boundsup = SCIP_INVALID;

   /*
    * check whether lower bound of child can be improved
    */
   SCIPintervalSet(&tmp, childinf);
   SCIPintervalEntropy(SCIP_INTERVAL_INFINITY, &tmp, tmp);

   /* entropy(childinf) < intersection.inf -> consider [childinf, MIN(childsup, extremum)] */
   if( SCIPintervalGetInf(intersection) > -SCIP_INTERVAL_INFINITY && SCIPintervalGetSup(tmp) - SCIPintervalGetInf(intersection) < -SCIPepsilon(scip) )
   {
      boundinf = reversePropBinarySearch(scip, childinf, MIN(extremum, childsup), TRUE,
         SCIPintervalGetInf(intersection));
   }
   /* entropy(childinf) > intersection.sup -> consider [MAX(childinf,extremum), childsup] */
   else if( SCIPintervalGetSup(intersection) < SCIP_INTERVAL_INFINITY && SCIPintervalGetInf(tmp) - SCIPintervalGetSup(intersection) > SCIPepsilon(scip) )
   {
      boundinf = reversePropBinarySearch(scip, MAX(childinf, extremum), childsup, FALSE,
         SCIPintervalGetSup(intersection));
   }
   /* using a strict greater-than here because we expect a tightening because we saw an at-least-epsilon-potential above */
   assert(boundinf == SCIP_INVALID || boundinf > childinf); /*lint !e777*/


   /*
    * check whether upper bound of child can be improved
    */
   if( childsup < SCIP_INTERVAL_INFINITY )
   {
      SCIPintervalSet(&tmp, childsup);
      SCIPintervalEntropy(SCIP_INTERVAL_INFINITY, &tmp, tmp);
   }
   else
      SCIPintervalSetBounds(&tmp, -SCIP_INTERVAL_INFINITY, -SCIP_INTERVAL_INFINITY);  /* entropy(inf) = -inf */

   /* entropy(childsup) < intersection.inf -> consider [MAX(childinf,extremum), childsup] */
   if( SCIPintervalGetInf(intersection) > -SCIP_INTERVAL_INFINITY && SCIPintervalGetSup(tmp) - SCIPintervalGetInf(intersection) < -SCIPepsilon(scip) )
   {
      boundsup = reversePropBinarySearch(scip, MAX(childinf, extremum), childsup, FALSE,
         SCIPintervalGetInf(intersection));
   }
   /* entropy(childsup) > intersection.sup -> consider [childinf, MIN(childsup,extremum)] */
   else if( SCIPintervalGetSup(intersection) < SCIP_INTERVAL_INFINITY && SCIPintervalGetInf(tmp) - SCIPintervalGetSup(intersection) > SCIPepsilon(scip) )
   {
      boundsup = reversePropBinarySearch(scip, childinf, MIN(childsup, extremum), TRUE,
         SCIPintervalGetSup(intersection));
   }
   /* using a strict smaller-than here because we expect a tightening because we saw an at-least-epsilon-potential above */
   assert(boundsup == SCIP_INVALID || boundsup < childsup); /*lint !e777*/

   if( boundinf != SCIP_INVALID ) /*lint !e777*/
   {
      childinf = MAX(childinf, boundinf);
   }
   if( boundsup != SCIP_INVALID ) /*lint !e777*/
   {
      childsup = boundsup;
   }
   assert(childinf <= childsup); /* infeasible case has been handled already */

   /* set the resulting bounds */
   SCIPintervalSetBounds(interval, childinf, childsup);

   return SCIP_OKAY;
}

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrEntropy)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrEntropy(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

/** simplifies an entropy expression */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyEntropy)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   /* check for value expression */
   if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrValue(conshdlr) )
   {
      SCIP_Real childvalue = SCIPgetConsExprExprValueValue(child);

      /* TODO how to handle a negative value? */
      assert(childvalue >= 0.0);

      if( childvalue == 0.0 || childvalue == 1.0 )
      {
         SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, 0.0) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, -childvalue * log(childvalue)) );
      }
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

/** expression data copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataEntropy)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPgetConsExprExprData(sourceexpr) == NULL);

   *targetexprdata = NULL;
   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataEntropy)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPsetConsExprExprData(expr, NULL);
   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseEntropy)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseConsExprExpr(scip, consexprhdlr, string, endstring, &childexpr) );
   assert(childexpr != NULL);

   /* create entropy expression */
   SCIP_CALL( SCIPcreateConsExprExprEntropy(scip, consexprhdlr, expr, childexpr) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the entropy expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}


/** expression (point-) evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalEntropy)
{  /*lint --e{715}*/
   SCIP_Real childvalue;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   childvalue = SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]);

   if( childvalue < 0.0 )
   {
      SCIPdebugMsg(scip, "invalid evaluation of entropy expression\n");
      *val = SCIP_INVALID;
   }
   else if( childvalue == 0.0 || childvalue == 1.0 )
   {
      /* -x*log(x) = 0 iff x in {0,1} */
      *val = 0.0;
   }
   else
   {
      *val = -childvalue * log(childvalue);
   }

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffEntropy)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_Real childvalue;

   assert(expr != NULL);
   assert(childidx == 0);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);

   childvalue = SCIPgetConsExprExprValue(child);

   /* derivative is not defined for x = 0 */
   if( childvalue <= 0.0 )
      *val = SCIP_INVALID;
   else
      *val = -1.0 - log(childvalue);

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalEntropy)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
      SCIPintervalEntropy(SCIP_INTERVAL_INFINITY, interval, childinterval);

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_CONSEXPR_EXPRESTIMATE(estimateEntropy)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* childvar;
   SCIP_Real refpoint;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);
   assert(success != NULL);

   *success = FALSE;

   /* get expression data */
   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childvar = SCIPgetConsExprExprAuxVar(child);
   assert(childvar != NULL);

   refpoint = SCIPgetSolVal(scip, sol, childvar);

   /* reference point is outside the domain of f(x) = -x*log(x) */
   if( refpoint < 0.0 )
      return SCIP_OKAY;

   /* use secant for underestimate (locally valid) */
   if( !overestimate )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real vallb;
      SCIP_Real valub;

      lb = SCIPvarGetLbLocal(childvar);
      ub = SCIPvarGetUbLocal(childvar);

      if( lb < 0.0 || SCIPisInfinity(scip, ub) || SCIPisEQ(scip, lb, ub) )
         return SCIP_OKAY;

      assert(lb >= 0.0 && ub >= 0.0);
      assert(ub - lb != 0.0);

      vallb = (lb == 0.0) ? 0.0 : -lb * log(lb);
      valub = (ub == 0.0) ? 0.0 : -ub * log(ub);

      *coefs = (valub - vallb) / (ub - lb);
      *constant = valub - *coefs * ub;
      assert(SCIPisEQ(scip, *constant, vallb - *coefs * lb));

      *islocal = TRUE;
   }
   /* use gradient cut for underestimate (globally valid) */
   else
   {
      /* no gradient cut possible if reference point is too close at 0 */
      if( SCIPisZero(scip, refpoint) )
         return SCIP_OKAY;

      /* -x*(1+log(x*)) + x* <= -x*log(x) */
      *coefs = -(1.0 + log(refpoint));
      *constant = refpoint;

      *islocal = FALSE;
      *branchcand = FALSE;
   }

   /* give up if the constant or coefficient is too large */
   if( SCIPisInfinity(scip, REALABS(*constant)) || SCIPisInfinity(scip, REALABS(*coefs)) )
      return SCIP_OKAY;

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_EXPRREVERSEPROP(reversepropEntropy)
{  /*lint --e{715}*/
   SCIP_INTERVAL newinterval;
   SCIP_INTERVAL exprinterval;
   SCIP_INTERVAL childinterval;

   SCIP_CONSEXPR_EXPR* child;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(nreductions != NULL);

   *nreductions = 0;

   child = SCIPgetConsExprExprChildren(expr)[0];
   childinterval = SCIPgetConsExprExprActivity(scip, child);
   exprinterval = SCIPgetConsExprExprActivity(scip, expr);

   /* compute resulting intervals */
   SCIP_CALL( reverseProp(scip, exprinterval, childinterval, &newinterval) );

   /* try to tighten the bounds of the child node */
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, conshdlr, child, newinterval, force, reversepropqueue, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** entropy hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashEntropy)
{  /*lint --e{715}*/
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
SCIP_DECL_CONSEXPR_EXPRCURVATURE(curvatureEntropy)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(childcurv != NULL);
   assert(success != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   /* to be concave, the child needs to be concave, too; we cannot be convex or linear */
   if( exprcurvature == SCIP_EXPRCURV_CONCAVE )
   {
      *childcurv = SCIP_EXPRCURV_CONCAVE;
      *success = TRUE;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(monotonicityEntropy)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_INTERVAL childbounds;
   SCIP_Real childinf;
   SCIP_Real childsup;
   SCIP_Real brpoint = exp(-1.0);

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx == 0);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   childbounds = SCIPgetConsExprExprActivity(scip, child);
   childinf = SCIPintervalGetInf(childbounds);
   childsup = SCIPintervalGetSup(childbounds);

   if( childsup <= brpoint )
      *result = SCIP_MONOTONE_INC;
   else if( childinf >= brpoint )
      *result = SCIP_MONOTONE_DEC;
   else
      *result = SCIP_MONOTONE_UNKNOWN;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEGRALITY(integralityEntropy)
{  /*lint --e{715}*/

   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   /* TODO it is possible to check for the special case that the child is integral and its bounds are [0,1]; in
    * this case the entropy expression can only achieve 0 and is thus integral
    */
   *isintegral = FALSE;

   return SCIP_OKAY;
}

/** creates the handler for x*log(x) expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrEntropy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLRDATA* exprhdlrdata;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   /* create expression handler data */
   exprhdlrdata = NULL;

   /* include expression handler */
   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
         EXPRHDLR_PRECEDENCE, evalEntropy, exprhdlrdata) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrEntropy, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataEntropy, freedataEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, NULL, NULL, estimateEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCurvature(scip, consexprhdlr, exprhdlr, curvatureEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrMonotonicity(scip, consexprhdlr, exprhdlr, monotonicityEntropy) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntegrality(scip, consexprhdlr, exprhdlr, integralityEntropy) );

   return SCIP_OKAY;
}

/** creates an x*log(x) expression */
SCIP_RETCODE SCIPcreateConsExprExprEntropy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< child expression */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(consexprhdlr != NULL);
   assert(expr != NULL);
   assert(child != NULL);

   exprhdlr = SCIPfindConsExprExprHdlr(consexprhdlr, EXPRHDLR_NAME);
   assert(exprhdlr != NULL);

   /* create expression data */
   exprdata = NULL;

   /* create expression */
   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, 1, &child) );

   return SCIP_OKAY;
}
