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

/**@file   expr_value.c
 * @brief  constant value expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/expr_value.h"

#define EXPRHDLR_NAME            "val"
#define EXPRHDLR_DESC            "constant value"
#define EXPRHDLR_PRECEDENCE      10000
#define EXPRHDLR_HASHKEY         SCIPcalcFibHash(36787.0)

/** the order of two values is the real order */
static
SCIP_DECL_EXPRCOMPARE(compareValue)
{  /*lint --e{715}*/
   SCIP_Real val1;
   SCIP_Real val2;

   val1 = SCIPgetValueExprValue(expr1);
   val2 = SCIPgetValueExprValue(expr2);

   return val1 < val2 ? -1 : val1 == val2 ? 0 : 1; /*lint !e777*/
}

static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrValue)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprHdlrValue(scip) );

   return SCIP_OKAY;
}

static
SCIP_DECL_EXPRCOPYDATA(copydataValue)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   *targetexprdata = SCIPexprGetData(sourceexpr);

   return SCIP_OKAY;
}

static
SCIP_DECL_EXPRFREEDATA(freedataValue)
{  /*lint --e{715}*/
   assert(expr != NULL);

   /* nothing much to do, as currently the data is the pointer */
   SCIPexprSetData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_EXPRPRINT(printValue)
{  /*lint --e{715}*/
   assert(expr != NULL);

   if( stage == SCIP_EXPRITER_ENTEREXPR )
   {
      SCIP_Real v = SCIPgetValueExprValue(expr);
      if( v < 0.0 && EXPRHDLR_PRECEDENCE <= parentprecedence )
      {
         SCIPinfoMessage(scip, file, "(%g)", v);
      }
      else
      {
         SCIPinfoMessage(scip, file, "%g", v);
      }
   }

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalValue)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   memcpy(val, &exprdata, sizeof(SCIP_Real));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffValue)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);

   /* this should never happen because variable expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalValue)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   SCIP_Real val;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   memcpy(&val, &exprdata, sizeof(SCIP_Real));

   SCIPintervalSet(interval, val);

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_EXPRHASH(hashValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);
   assert(hashkey != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(SCIPgetValueExprValue(expr));

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);

   *result = SCIP_MONOTONE_CONST;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityValue)
{  /*lint --e{715}*/

   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   *isintegral = EPSISINT(SCIPgetValueExprValue(expr), 0.0); /*lint !e835 !e666*/

   return SCIP_OKAY;
}

/** creates the handler for constant value expression and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprHdlrValue(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeExprHdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE, evalValue, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrValue, NULL);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataValue, freedataValue);
   SCIPexprhdlrSetCompare(exprhdlr, compareValue);
   SCIPexprhdlrSetPrint(exprhdlr, printValue);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalValue);
   SCIPexprhdlrSetHash(exprhdlr, hashValue);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffValue, NULL, NULL);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureValue);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityValue);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityValue);

   return SCIP_OKAY;
}

/** creates constant value expression */
SCIP_RETCODE SCIPcreateExprValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_Real             value,              /**< value to be stored */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata  /**< data to pass to ownerdatacreate */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(SCIPisFinite(value));

   assert(sizeof(SCIP_Real) <= sizeof(SCIP_EXPRDATA*));
   memcpy(&exprdata, &value, sizeof(SCIP_Real));

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPgetExprHdlrValue(scip), exprdata, 0, NULL, ownerdatacreate, ownerdatacreatedata) );

   return SCIP_OKAY;
}

/* from pub_expr.h */

/** gets the value of a constant value expression */
SCIP_Real SCIPgetValueExprValue(
   SCIP_EXPR*   expr                /**< expression */
   )
{
   SCIP_EXPRDATA* exprdata;
   SCIP_Real v;

   assert(sizeof(SCIP_Real) <= sizeof(SCIP_EXPRDATA*));

   exprdata = SCIPexprGetData(expr);
   memcpy(&v, &exprdata, sizeof(SCIP_Real));

   return v;
}
