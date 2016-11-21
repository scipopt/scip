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

/**@file   cons_expr_pow.c
 * @brief  power expression handler
 * @author Benjamin Mueller
 *
 * @todo initsepaPow
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_value.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_sum.h"

#define POW_PRECEDENCE  55000
#define POW_HASHKEY     SCIPcalcFibHash(21163.0)

/*
 * Data structures
 */

struct SCIP_ConsExpr_ExprData
{
   SCIP_Real  exponent;     /**< exponent */
};

/*
 * Local methods
 */

static
SCIP_RETCODE createData(
   SCIP*                    scip,            /**< SCIP data structure */
   SCIP_CONSEXPR_EXPRDATA** exprdata,        /**< pointer where to store expression data */
   SCIP_Real                exponent         /**< exponent of the power expression */
   )
{
   assert(exprdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, exprdata) );

   (*exprdata)->exponent = exponent;

   return SCIP_OKAY;
}


/*
 * Callback methods of expression handler
 */

/** the order of two pow is base_1^expo_1 < base_2^expo_2 if and only if
 * base_1 < base2 or, base_1 = base_2 and expo_1 < expo_2
 */
static
SCIP_DECL_CONSEXPR_EXPRCMP(comparePow)
{  /*lint --e{715}*/
   SCIP_Real expo1;
   SCIP_Real expo2;
   int compareresult;

   compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[0],
              SCIPgetConsExprExprChildren(expr2)[0]);
   if( compareresult != 0 )
      return compareresult;

   expo1 = SCIPgetConsExprExprPowExponent(expr1);
   expo2 = SCIPgetConsExprExprPowExponent(expr2);

   return  expo1 == expo2 ? 0 : expo1 < expo2 ? -1 : 1;
}

/** simplifies a pow expression.
 * Evaluates the power function when its child is a value expression
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyPow)
{
   SCIP_CONSEXPR_EXPR* base;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   base = SCIPgetConsExprExprChildren(expr)[0];
   assert(base != NULL);

   exponent = SCIPgetConsExprExprPowExponent(expr);
   /* when exponent is inteer, round exponent so that is actually an integer
    * TODO: should this go in the createConsExprExprPow? */
   if( SCIPisIntegral(scip, exponent) )
      exponent = SCIPround(scip, exponent);

   SCIPdebugPrintf("[simplifyPow] simplifying power with expo %g\n", exponent);

   /* enforces POW1 */
   if( exponent == 0.0 )
   {
      SCIPdebugPrintf("[simplifyPow] POW1\n");
      /* TODO: more checks? */
      if( SCIPgetConsExprExprHdlr(base) == SCIPgetConsExprExprHdlrValue(conshdlr) &&
            SCIPgetConsExprExprValueValue(base) == 0.0 )
      {
         assert(0);
      }
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, 1.0) );
      return SCIP_OKAY;
   }

   /* enforces POW2 */
   if( exponent == 1.0 )
   {
      SCIPdebugPrintf("[simplifyPow] POW2\n");
      *simplifiedexpr = base;
      SCIPcaptureConsExprExpr(*simplifiedexpr);
      return SCIP_OKAY;
   }

   /* enforces POW3 */
   if( SCIPgetConsExprExprHdlr(base) == SCIPgetConsExprExprHdlrValue(conshdlr) )
   {
      SCIP_Real baseval;

      SCIPdebugPrintf("[simplifyPow] POW3\n");
      baseval = SCIPgetConsExprExprValueValue(base);

      /* TODO check if those are all important asserts */
      assert(baseval >= 0.0 || fmod(exponent, 1.0) == 0.0);
      assert(baseval != 0.0 || exponent != 0.0);

      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, pow(baseval, exponent)) );
      return SCIP_OKAY;
   }

   if( SCIPisIntegral(scip, exponent) )
   {
      SCIP_CONSEXPR_EXPR* aux;
      SCIP_CONSEXPR_EXPR* simplifiedaux;

      /* enforces POW5
       * given (pow n (prod 1.0 expr_1 ... expr_k)) we distribute the exponent:
       * -> (prod 1.0 (pow n expr_1) ... (pow n expr_k))
       * notes: - since base is simplified and its coefficient is 1.0 (SP8)
       *        - n is an integer (excluding 1 and 0; see POW1-2 above)
       */
      if( SCIPgetConsExprExprHdlr(base) == SCIPgetConsExprExprHdlrProduct(conshdlr) )
      {
         int i;
         SCIP_CONSEXPR_EXPR* auxproduct;

         /* create empty product */
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &auxproduct, 0, NULL, 1.0) );

         for( i = 0; i < SCIPgetConsExprExprNChildren(base); ++i )
         {
            /* create (pow n expr_i) and simplify */
            SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &aux,
                     SCIPgetConsExprExprChildren(base)[i], exponent) );
            SCIP_CALL( simplifyPow(scip, aux, &simplifiedaux) );
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );

            /* append (pow n expr_i) to product */
            SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, auxproduct, simplifiedaux) );
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplifiedaux) );
         }

         /* simplify (prod 1.0 (pow n expr_1) ... (pow n expr_k))
          * TODO: ideally we would call simplifyProduct directly, since we know its children are simplified! */
         SCIP_CALL( SCIPsimplifyConsExprExpr(scip, auxproduct, simplifiedexpr) ); /*FIXME*/
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &auxproduct) );
         return SCIP_OKAY;
      }

      /* enforces POW6
       * given (pow n (sum 0.0 coef expr)) we can move `pow` inside `sum`:
       * (pow n (sum 0.0 coef expr) ) -> (sum 0.0 coef^n (pow n expr))
       * notes: - since base is simplified and its constant is 0, then coef != 1.0 (SS7)
       *        - n is an integer (excluding 1 and 0; see POW1-2 above)
       */
      if( SCIPgetConsExprExprHdlr(base) == SCIPgetConsExprExprHdlrSum(conshdlr)
            && SCIPgetConsExprExprNChildren(base) == 1
            && SCIPgetConsExprExprSumConstant(base) == 0.0 )
      {
         SCIP_Real newcoef;

         SCIPdebugPrintf("[simplifyPow] seing a sum with one term, exponent %g\n", exponent);
         /* assert SS7 holds */
         assert(SCIPgetConsExprExprSumCoefs(base)[0] != 1.0);

         /* create (pow n expr) and simplify it
          * note: we call simplifyPow directly, since we know that `expr` is simplified */
         newcoef = pow(SCIPgetConsExprExprSumCoefs(base)[0], exponent);
         SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &aux, SCIPgetConsExprExprChildren(base)[0], exponent) );
         SCIP_CALL( simplifyPow(scip, aux, &simplifiedaux) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );

         /* create (sum (pow n expr)) and simplify it
          * TODO: ideally we would call simplifySum directly, since we know its child is simplified! */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &aux, 1, &simplifiedaux, &newcoef, 0.0) );
         SCIP_CALL( SCIPsimplifyConsExprExpr(scip, aux, simplifiedexpr) ); /*FIXME*/
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplifiedaux) );
         return SCIP_OKAY;
      }

      /* enforces POW7
       * (const + sum alpha_i expr_i)^2 = sum alpha_i^2 expr_i^2
       * + sum_{j < i} 2*alpha_i alpha_j expr_i expr_j
       * + sum const alpha_i expr_i
       * TODO: put some limits on the number of children of the sum being expanded
       */
      if( SCIPgetConsExprExprHdlr(base) == SCIPgetConsExprExprHdlrSum(conshdlr) && exponent == 2 )
      {
         int i;
         int nchildren;
         int nexpandedchildren;
         SCIP_CONSEXPR_EXPR* expansion;
         SCIP_CONSEXPR_EXPR** expandedchildren;
         SCIP_Real* coefs;
         SCIP_Real constant;

         SCIPdebugPrintf("[simplifyPow] expanding sum^%g\n", exponent);

         nchildren = SCIPgetConsExprExprNChildren(base);
         nexpandedchildren = nchildren * (nchildren + 1) / 2 + nchildren;
         SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nexpandedchildren) );
         SCIP_CALL( SCIPallocBufferArray(scip, &expandedchildren, nexpandedchildren) );

         for( i = 0; i < nchildren; ++i )
         {
            int j;
            SCIP_CONSEXPR_EXPR* expansionchild;
            SCIP_CONSEXPR_EXPR* prodchildren[2];
            prodchildren[0] = SCIPgetConsExprExprChildren(base)[i];

            /* create and simplify expr_i * expr_j */
            for( j = 0; j < i; ++j )
            {
               prodchildren[1] = SCIPgetConsExprExprChildren(base)[j];
               coefs[i*(i+1)/2 + j] = 2 * SCIPgetConsExprExprSumCoefs(base)[i] * SCIPgetConsExprExprSumCoefs(base)[j];

               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expansionchild, 2, prodchildren, 1.0) );
               SCIP_CALL( SCIPsimplifyConsExprExpr(scip, expansionchild, &expandedchildren[i*(i+1)/2 + j]) ); /* FIXME: call simplifyProduct */
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expansionchild) );
            }
            /* create and simplify expr_i * expr_i */
            prodchildren[1] = SCIPgetConsExprExprChildren(base)[i];
            coefs[i*(i+1)/2 + i] = SCIPgetConsExprExprSumCoefs(base)[i] * SCIPgetConsExprExprSumCoefs(base)[i];

            SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expansionchild, 2, prodchildren, 1.0) );
            SCIP_CALL( SCIPsimplifyConsExprExpr(scip, expansionchild, &expandedchildren[i*(i+1)/2 + i]) ); /* FIXME: call simplifyProduct */
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expansionchild) );
         }
         /* create const * alpha_i expr_i */
         for( i = 0; i < nchildren; ++i )
         {
            coefs[i + nexpandedchildren - nchildren] = 2 * SCIPgetConsExprExprSumConstant(base) * SCIPgetConsExprExprSumCoefs(base)[i];
            expandedchildren[i + nexpandedchildren - nchildren] = SCIPgetConsExprExprChildren(base)[i];
         }

         constant = SCIPgetConsExprExprSumConstant(base);
         constant *= constant;
         /* create sum of all the above and simplify it with simplifySum since all of its children are simplified! */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expansion, nexpandedchildren,
                  expandedchildren, coefs, constant) );
         SCIP_CALL( SCIPsimplifyConsExprExpr(scip, expansion, simplifiedexpr) ); /* FIXME: call simplifySum */

         /* release eveything */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expansion) );
         /* release the *created* expanded children */
         for( i = 0; i < nexpandedchildren - nchildren; ++i )
         {
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expandedchildren[i]) );
         }
         SCIPfreeBufferArray(scip, &coefs);
         SCIPfreeBufferArray(scip, &expandedchildren);

         return SCIP_OKAY;
      }

      /* enforces POW8
       * given (pow n (pow expo expr)) we distribute the exponent:
       * -> (pow n*expo expr)
       * notes: n is an integer (excluding 1 and 0; see POW1-2 above)
       */
      /* FIXME: use SCIPgetConsExprExprHdlrPow */
      if( strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(base)), "pow") == 0 )
      {
         SCIP_Real newexponent;

         newexponent = SCIPgetConsExprExprPowExponent(base) * exponent;
         SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &aux, SCIPgetConsExprExprChildren(base)[0], newexponent) );
         SCIP_CALL( simplifyPow(scip, aux, simplifiedexpr) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );

         return SCIP_OKAY;
      }
   }

   SCIPdebugPrintf("[simplifyPow] power is simplified\n");
   *simplifiedexpr = expr;

   /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
   SCIPcaptureConsExprExpr(*simplifiedexpr);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrPow)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrPow(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataPow)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* sourceexprdata;

   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   sourceexprdata = SCIPgetConsExprExprData(sourceexpr);
   assert(sourceexprdata != NULL);

   *targetexprdata = NULL;

   SCIP_CALL( createData(targetscip, targetexprdata, sourceexprdata->exponent) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataPow)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIPfreeBlockMemory(scip, &exprdata);
   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

/** @todo: use precedence for better printing */
static
SCIP_DECL_CONSEXPR_EXPRPRINT(printPow)
{
   assert(expr != NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print function with opening parenthesis */
         SCIPinfoMessage(scip, file, "(");
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
         SCIPinfoMessage(scip, file, ")^%g", SCIPgetConsExprExprPowExponent(expr));
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      default: ;
   }

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalPow)
{  /*lint --e{715}*/
   SCIP_Real exponent;
   SCIP_Real base;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   exponent = SCIPgetConsExprExprPowExponent(expr);
   base = SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]);

   *val = pow(base, exponent);

   /* if there is a domain, pole, or range error, pow() should return some kind of NaN, infinity, or HUGE_VAL
    * we could also work with floating point exceptions or errno, but I am not sure this would be thread-safe
    */
   if( !SCIPisFinite(*val) || *val == HUGE_VAL || *val == -HUGE_VAL )
      *val = SCIP_INVALID;

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffPow)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_Real exponent;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);
   assert(idx >= 0 && idx < SCIPgetConsExprExprNChildren(expr));
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID);

   child = SCIPgetConsExprExprChildren(expr)[idx];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);

   exponent = SCIPgetConsExprExprPowExponent(expr);
   assert(exponent != 1.0 && exponent != 0.0);

   /* x^exponent is not differentiable for x = 0 and exponent in ]0,1[ */
   if( exponent > 0.0 && exponent < 1.0 && SCIPgetConsExprExprValue(child) == 0.0 )
      *val = SCIP_INVALID;
   else
      *val = exponent * pow(SCIPgetConsExprExprValue(child), exponent - 1.0);

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalPow)
{
   SCIP_INTERVAL childinterval;
   SCIP_Real exponent;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

   exponent = SCIPgetConsExprExprPowExponent(expr);

   SCIPintervalPowerScalar(SCIPinfinity(scip), interval, childinterval, exponent);

   return SCIP_OKAY;
}

/** helper function to separate a given point; needed for proper unittest */
static
SCIP_RETCODE separatePointPow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_ROW**            cut                 /**< pointer to store the row */
   )
{
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* auxvar;
   SCIP_VAR* childvar;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real exponent;
   SCIP_Bool overestimate;
   SCIP_Real violation;
   SCIP_Real refpoint;
   SCIP_Real lincoef;
   SCIP_Real linconstant;
   SCIP_Bool islocal;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "pow") == 0);
   assert(cut != NULL);

   exponent = SCIPgetConsExprExprPowExponent(expr);
   assert(exponent != 1.0 && exponent != 0.0); /* this should have been simplified */

   /* get expression data */
   auxvar = SCIPgetConsExprExprLinearizationVar(expr);
   assert(auxvar != NULL);
   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childvar = SCIPgetConsExprExprLinearizationVar(child);
   assert(childvar != NULL);

   *cut = NULL;
   refpoint = SCIPgetSolVal(scip, sol, childvar);

   /* we can not generate a cut at +/- infinity */
   if( SCIPisInfinity(scip, REALABS(refpoint)) )
      return SCIP_OKAY;

   /* compute the violation; this determines whether we need to over- or underestimate */
   violation = pow(refpoint, exponent) - SCIPgetSolVal(scip, sol, auxvar);

   /* check if there is a violation */
   if( SCIPisEQ(scip, violation, 0.0) )
      return SCIP_OKAY;

   overestimate = SCIPisLT(scip, violation, 0.0);
   success = TRUE;
   lincoef = 0.0;
   linconstant = 0.0;

   /* adjust the reference points */
   childlb = SCIPvarGetLbLocal(childvar);
   childub = SCIPvarGetUbLocal(childvar);
   refpoint = SCIPisLT(scip, refpoint, childlb) ? childlb : refpoint;
   refpoint = SCIPisGT(scip, refpoint, childub) ? childub : refpoint;
   assert(SCIPisLE(scip, refpoint, childub) && SCIPisGE(scip, refpoint, childlb));

   /* quadratic case */
   if( exponent == 2.0 )
   {
      if( overestimate )
      {
         SCIPaddSquareSecant(scip, 1.0, childlb, childub, refpoint, &lincoef, &linconstant, &success);
         islocal = TRUE; /* secants are only valid locally */
      }
      else
      {
         SCIPaddSquareLinearization(scip, 1.0, refpoint, SCIPvarIsIntegral(childvar), &lincoef, &linconstant, &success);
         islocal = FALSE; /* linearizations are globally valid */
      }

      /* @todo  allow lhs/rhs of +/- infinity? */
      if( success && !SCIPisInfinity(scip, REALABS(linconstant)) )
      {
         SCIP_CALL( SCIPcreateRowCons(scip, cut, conshdlr, islocal ? "square_secant" : "square_linearization", 0, NULL, NULL, -SCIPinfinity(scip),
               SCIPinfinity(scip), islocal, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *cut, childvar, lincoef) );
         SCIP_CALL( SCIPaddVarToRow(scip, *cut, auxvar, -1.0) );

         if( overestimate )
         {
            SCIP_CALL( SCIPchgRowLhs(scip, *cut, -linconstant) );
         }
         else
         {
            SCIP_CALL( SCIPchgRowRhs(scip, *cut, -linconstant) );
         }
      }
   }
   else
   {
      /* @todo can not be handled so far */
   }

   return SCIP_OKAY;
}

/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaPow)
{
   SCIP_ROW* cut;
   SCIP_Bool infeasible;

   cut = NULL;
   *ncuts = 0;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( separatePointPow(scip, conshdlr, expr, sol, &cut) );

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
   SCIPdebugMsg(scip, "add cut with violation %e\n", minviolation);
   SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** expression reverse propagaton callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropPow)
{
   SCIP_INTERVAL interval;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(nreductions != NULL);

   *nreductions = 0;

   /* not possible to learn bounds if expression interval is unbounded in both directions */
   if( SCIPintervalIsEntire(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)) )
      return SCIP_OKAY;

   exponent = SCIPgetConsExprExprPowExponent(expr);

   /* f = pow(c0, alpha) -> c0 = pow(f, 1/alpha) */
   SCIPintervalPowerScalarInverse(SCIPinfinity(scip), &interval,
      SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]), exponent, SCIPgetConsExprExprInterval(expr));

   /* try to tighten the bounds of the child node */
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, SCIPgetConsExprExprChildren(expr)[0], interval, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashPow)
{
   unsigned int childhash;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   *hashkey = POW_HASHKEY;

   assert(SCIPhashmapExists(expr2key, (void*)SCIPgetConsExprExprChildren(expr)[0]));
   childhash = (unsigned int)(size_t)SCIPhashmapGetImage(expr2key, SCIPgetConsExprExprChildren(expr)[0]);

   *hashkey ^= childhash;

   return SCIP_OKAY;
}

/** creates the handler for power expression and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrPow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "pow", "power expression",
         POW_PRECEDENCE, evalPow, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrPow, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataPow, freedataPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, comparePow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffPow) );

   return SCIP_OKAY;
}

/** creates a power expression */
SCIP_RETCODE SCIPcreateConsExprExprPow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child,              /**< single child */
   SCIP_Real             exponent            /**< exponent of the power expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindConsExprExprHdlr(consexprhdlr, "pow") != NULL);

   SCIP_CALL( createData(scip, &exprdata, exponent) );
   assert(exprdata != NULL);

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPfindConsExprExprHdlr(consexprhdlr, "pow"), exprdata, 1, &child) );

   return SCIP_OKAY;
}

/** gets the exponent of a power expression */
SCIP_Real SCIPgetConsExprExprPowExponent(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   return exprdata->exponent;
}
