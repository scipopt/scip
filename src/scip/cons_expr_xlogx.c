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

/**@file   cons_expr_xlogx.c
 * @brief  handler for x*log(x) expressions
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_xlogx.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr.h"

#include <string.h>

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "xlogx"
#define EXPRHDLR_DESC         "expression handler for x*log(x)"
#define EXPRHDLR_PRECEDENCE   0
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(7477.0)

/*
 * Data structures
 */

/*
 * Local methods
 */

/** helper function to separate a given point; needed for proper unittest */
static
SCIP_RETCODE separatePointXlogx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< xlogx expression */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_ROW**            cut                 /**< pointer to store the row */
   )
{
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* auxvar;
   SCIP_VAR* childvar;
   SCIP_Real refpoint;
   SCIP_Real activity;
   SCIP_Real violation;
   SCIP_Real coef;
   SCIP_Real constant;
   SCIP_Bool overestimate;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(cut != NULL);

   *cut = NULL;

   /* get linearization variable */
   auxvar = SCIPgetConsExprExprLinearizationVar(expr);
   assert(auxvar != NULL);

   /* get expression data */
   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childvar = SCIPgetConsExprExprLinearizationVar(child);
   assert(childvar != NULL);

   refpoint = SCIPgetSolVal(scip, sol, childvar);

   /* reference point is outside the domain of f(x) = x*log(x) */
   if( refpoint < 0.0 )
      return SCIP_OKAY;

   activity = (refpoint == 0.0) ? 0.0 : refpoint * log(refpoint);
   violation = activity - SCIPgetSolVal(scip, sol, auxvar);
   overestimate = SCIPisLT(scip, violation, 0.0);

   /* use secant for overestimate (locally valid) */
   if( overestimate )
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

      vallb = (lb == 0.0) ? 0.0 : lb * log(lb);
      valub = (ub == 0.0) ? 0.0 : ub * log(ub);

      coef = (valub - vallb) / (ub - lb);
      constant = valub - coef * ub;
      assert(SCIPisEQ(scip, constant, vallb - coef * lb));
   }
   /* use gradient cut for underestimate (globally valid) */
   else
   {
      /* no gradient cut possible if reference point is too close at 0 */
      if( SCIPisZero(scip, refpoint) )
         return SCIP_OKAY;

      /* x*(1+log(x*)) - x* <= x*log(x) */
      coef = (1.0 + log(refpoint));
      constant = -refpoint;
   }

   /* create cut */
   SCIP_CALL( SCIPcreateRowCons(scip, cut, conshdlr, "xlogx_cut", 0, NULL, NULL,
         overestimate ? -constant : -SCIPinfinity(scip),
         overestimate ? SCIPinfinity(scip) : -constant,
         overestimate, FALSE, FALSE) );

   SCIP_CALL( SCIPaddVarToRow(scip, *cut, auxvar, -1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, *cut, childvar, coef) );

   return SCIP_OKAY;
}

/** helper function for reverseProp() which returns an x* in [xmin,xmax] s.t. the distance x*log(x) and a given target
 *  value is minimized; the function assumes that x*log(x) is monotone on [xmin,xmax];
 */
static
SCIP_Real reversePropBinarySearch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             xmin,               /**< smallest possible x */
   SCIP_Real             xmax,               /**< largest possible x */
   SCIP_Bool             increasing,         /**< x*log(x) is increasing or decreasing on [xmin,xmax] */
   SCIP_Real             targetval           /**< target value */
   )
{
   SCIP_Real xminval = (xmin == 0.0) ? 0.0 : xmin * log(xmin);
   SCIP_Real xmaxval = (xmax == 0.0) ? 0.0 : xmax * log(xmax);
   int i;

   assert(increasing ? xminval <= xmaxval : xminval >= xmaxval);

   /* function can not achieve x*log(x) -> return xmin or xmax */
   if( SCIPisGE(scip, xminval, targetval) && SCIPisGE(scip, xmaxval, targetval) )
      return increasing ? xmin : xmax;
   else if( SCIPisLE(scip, xminval, targetval) && SCIPisLE(scip, xmaxval, targetval) )
      return increasing ? xmax : xmin;

   /* binary search */
   for( i = 0; i < 1000; ++i )
   {
      SCIP_Real x = (xmin + xmax) / 2.0;
      SCIP_Real xval = (x == 0.0) ? 0.0 : x * log(x);

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
   SCIP_Real childsup;
   SCIP_Real childinf;
   SCIP_Real exprsup;
   SCIP_Real exprinf;
   SCIP_Real bound;

   assert(scip != NULL);
   assert(interval != NULL);

   /* check whether domain is empty, i.e., bounds on x*log(x) < -1/e */
   if( SCIPisLT(scip, SCIPintervalGetSup(exprinterval), -exp(-1.0))
      || SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval) )
   {
      SCIPintervalSetEmpty(interval);
      return SCIP_OKAY;
   }

   childsup = SCIPintervalGetSup(childinterval);
   childinf = SCIPintervalGetInf(childinterval);
   exprsup = SCIPintervalGetSup(exprinterval);
   exprinf = SCIPintervalGetInf(exprinterval);

   /*
    * consider bounds implied by upper bound on the expression
    */

   bound = reversePropBinarySearch(scip, exp(-1.0), childsup, TRUE, exprsup);
   assert(bound <= childsup);
   childsup = MIN(bound, childsup);

   if( SCIPisLT(scip, exprsup, 0.0) )
   {
      bound = reversePropBinarySearch(scip, childinf, exp(-1.0), FALSE, exprsup);
      assert(bound >= childinf);
      childinf = MAX(childinf, bound);
   }

   /*
    * consider bounds implied by lower bound on the expression
    */

   /* lower bound on expression can only imply a better lower bound on the child's interval */
   if( SCIPisGT(scip, exprinf, 0.0) && childinf >= exp(-1.0) )
   {
      bound = reversePropBinarySearch(scip, exp(-1.0), childsup, TRUE, exprinf);
      assert(bound >= childinf);
      childinf = MAX(childinf, bound);
   }

   /* set the resulting bounds */
   SCIPintervalSetBounds(interval, childinf, childsup);

   return SCIP_OKAY;
}

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrXlogx)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrXlogx(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

/** simplifies an xlogx expression */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyXlogx)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_CONSHDLR* conshdlr;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

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
         SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, childvalue*log(childvalue)) );
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
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataXlogx)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPgetConsExprExprData(sourceexpr) == NULL);

   *targetexprdata = NULL;
   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataXlogx)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPsetConsExprExprData(expr, NULL);
   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_CONSEXPR_EXPRPRINT(printXlogx)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print function with opening parenthesis */
         SCIPinfoMessage(scip, file, "xlogx(");
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      {
         assert(SCIPgetConsExprExprWalkCurrentChild(expr) == 0);
         break;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
      {
         /* print closing parenthesis */
         SCIPinfoMessage(scip, file, ")");
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      default: ;
   }

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalXlogx)
{  /*lint --e{715}*/
   SCIP_Real childvalue;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   childvalue = SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]);

   if( childvalue < 0.0 )
   {
      SCIPdebugMsg(scip, "invalid evaluation of xlogx expression\n");
      *val = SCIP_INVALID;
   }
   else if( childvalue == 0.0 || childvalue == 1.0 )
   {
      /* x*log(x) = 0 iff x in {0,1} */
      *val = 0.0;
   }
   else
   {
      *val = childvalue * log(childvalue);
   }

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffXlogx)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_Real childvalue;

   assert(expr != NULL);
   assert(idx == 0);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPgetConsExprExprChildren(expr)[idx];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);

   childvalue = SCIPgetConsExprExprValue(child);

   /* derivative is not defined for x = 0 */
   if( childvalue <= 0.0 )
      *val = SCIP_INVALID;
   else
      *val = 1.0 + childvalue;

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalXlogx)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

   SCIPintervalLog(SCIPinfinity(scip), interval, childinterval);
   SCIPintervalMul(SCIPinfinity(scip), interval, *interval, childinterval);

   return SCIP_OKAY;
}

/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaXlogx)
{  /*lint --e{715}*/
   SCIP_ROW* cut;
   SCIP_Bool infeasible;

   cut = NULL;
   *ncuts = 0;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( separatePointXlogx(scip, conshdlr, expr, sol, &cut) );

   /* failed to compute a cut */
   if( cut == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPmassageConsExprExprCut(scip, &cut, sol, minviolation) );

   /* cut violation or numerics were too bad */
   if( cut == NULL )
      return SCIP_OKAY;

   /* add cut */
   SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
   *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
   *ncuts += 1;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "add cut with violation %e\n", violation);
   SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropXlogx)
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
   childinterval = SCIPgetConsExprExprInterval(child);
   exprinterval = SCIPgetConsExprExprInterval(expr);

   /* compute resulting intervals */
   SCIP_CALL( reverseProp(scip, exprinterval, childinterval, &newinterval) );

   /* try to tighten the bounds of the child node */
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, child, newinterval, force, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** xlogx hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashXlogx)
{  /*lint --e{715}*/
   unsigned int childhash;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   *hashkey = EXPRHDLR_HASHKEY;

   assert(SCIPhashmapExists(expr2key, (void*)SCIPgetConsExprExprChildren(expr)[0]));
   childhash = (unsigned int)(size_t)SCIPhashmapGetImage(expr2key, SCIPgetConsExprExprChildren(expr)[0]);

   *hashkey ^= childhash;

   return SCIP_OKAY;
}

/** creates the handler for x*log(x) expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrXlogx(
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
         EXPRHDLR_PRECEDENCE, evalXlogx, exprhdlrdata) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrXlogx, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataXlogx, freedataXlogx) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyXlogx) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printXlogx) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalXlogx) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaXlogx) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropXlogx) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashXlogx) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffXlogx) );

   return SCIP_OKAY;
}

/** creates an x*log(x) expression */
SCIP_RETCODE SCIPcreateConsExprExprXlogx(
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

   exprhdlr = SCIPgetConsExprExprHdlrXlogx(consexprhdlr);
   assert(exprhdlr != NULL);

   /* create expression data */
   exprdata = NULL;

   /* create expression */
   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, 1, &child) );

   return SCIP_OKAY;
}
