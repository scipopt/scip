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

/** simplifies a pow expression.
 * Evaluates the power function when its child is a value expression
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyPow)
{
   SCIP_CONSEXPR_EXPR* child;
   SCIP_CONSHDLR* conshdlr;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   exponent = SCIPgetConsExprExprPowExponent(expr);

   /* check for value expression */
   if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrValue(conshdlr) )
   {
      SCIP_Real base;

      base = SCIPgetConsExprExprValueValue(child);

      /* TODO check if those are all important asserts */
      assert(base >= 0.0 || fmod(exponent, 2.0) == 1.0);
      assert(base != 0.0 || exponent != 0.0);

      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, pow(base, exponent)) );
   }
   else if( exponent == 1.0 )
   {
      /* @todo is it necessary to copy or capture the children? */
      *simplifiedexpr = child;
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }

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
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPgetConsExprExprData(sourceexpr) == NULL);

   *targetexprdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataPow)
{
   assert(expr != NULL);

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printPow)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print function with opening parenthesis */
         SCIPinfoMessage(scip, file, "pow(");
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

static
SCIP_DECL_CONSEXPR_EXPRPARSE(parsePow)
{
   assert(expr != NULL);

   /* @todo the parsing needs to be implemented in cons_expr.c if we want to use x^y instead of pow(x,y) */

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPREVAL(evalPow)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   if( SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) <= 0.0 )
   {
      SCIPdebugMessage("invalid evaluation of power expression\n");
      *val = SCIP_INVALID;
   }
   else
   {
      SCIP_Real exponent;
      SCIP_Real base;

      exponent = SCIPgetConsExprExprPowExponent(expr);
      base = SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]);

      /* TODO check if those are all important asserts */
      assert(base >= 0.0 || fmod(exponent, 2.0) == 1.0);
      assert(base != 0.0 || exponent != 0.0);

      *val = pow(base, exponent);
   }

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalPow)
{
   SCIP_INTERVAL childinterval;
   SCIP_Real exponent;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
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

   /* compute the violation; this determines whether we need to over- or underestimate */
   violation = pow(refpoint, exponent) - SCIPgetSolVal(scip, sol, auxvar);

   /* check if there is a violation */
   if( SCIPisEQ(scip, violation, 0.0) )
      return SCIP_OKAY;

   overestimate = SCIPisLT(scip, violation, 0.0);
   success = FALSE;

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
   SCIPdebugMessage("add cut with violation %e\n", violation);
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
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parsePow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropPow) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashPow) );

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
