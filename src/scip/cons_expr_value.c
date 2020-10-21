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

/**@file   cons_expr_value.c
 * @brief  constant value expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_value.h"

#define EXPRHDLR_NAME            "val"
#define EXPRHDLR_DESC            "constant value"
#define EXPRHDLR_PRECEDENCE      10000
#define EXPRHDLR_HASHKEY         SCIPcalcFibHash(36787.0)

/** the order of two values is the real order */
static
SCIP_DECL_CONSEXPR_EXPRCOMPARE(compareValue)
{  /*lint --e{715}*/
   SCIP_Real val1;
   SCIP_Real val2;

   val1 = SCIPgetConsExprExprValueValue(expr1);
   val2 = SCIPgetConsExprExprValueValue(expr2);

   return val1 < val2 ? -1 : val1 == val2 ? 0 : 1; /*lint !e777*/
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrValue)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrValue(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataValue)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   *targetexprdata = SCIPgetConsExprExprData(sourceexpr);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataValue)
{  /*lint --e{715}*/
   assert(expr != NULL);

   /* nothing much to do, as currently the data is the pointer */
   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printValue)
{  /*lint --e{715}*/
   assert(expr != NULL);

   if( stage == SCIP_CONSEXPRITERATOR_ENTEREXPR )
   {
      SCIP_Real v = SCIPgetConsExprExprValueValue(expr);
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
SCIP_DECL_CONSEXPR_EXPREVAL(evalValue)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   memcpy(val, &exprdata, sizeof(SCIP_Real));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffValue)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);

   /* this should never happen because variable expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalValue)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_Real val;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   memcpy(&val, &exprdata, sizeof(SCIP_Real));

   SCIPintervalSet(interval, val);

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 0);
   assert(hashkey != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(SCIPgetConsExprExprValueValue(expr));

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_CONSEXPR_EXPRCURVATURE(curvatureValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(success != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 0);

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(monotonicityValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 0);

   *result = SCIP_MONOTONE_CONST;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEGRALITY(integralityValue)
{  /*lint --e{715}*/

   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   *isintegral = EPSISINT(SCIPgetConsExprExprValueValue(expr), 0.0); /*lint !e835 !e666*/

   return SCIP_OKAY;
}

/** creates the handler for constant value expression and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE, evalValue, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrValue, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataValue, freedataValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrDiff(scip, consexprhdlr, exprhdlr, bwdiffValue, NULL, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCurvature(scip, consexprhdlr, exprhdlr, curvatureValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrMonotonicity(scip, consexprhdlr, exprhdlr, monotonicityValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntegrality(scip, consexprhdlr, exprhdlr, integralityValue) );

   return SCIP_OKAY;
}

/** creates constant value expression */
SCIP_RETCODE SCIPcreateConsExprExprValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_Real             value               /**< value to be stored */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(consexprhdlr != NULL);
   assert(expr != NULL);
   assert(SCIPisFinite(value));

   assert(sizeof(SCIP_Real) <= sizeof(SCIP_CONSEXPR_EXPRDATA*));
   memcpy(&exprdata, &value, sizeof(SCIP_Real));

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPgetConsExprExprHdlrValue(consexprhdlr), exprdata, 0, NULL) );

   return SCIP_OKAY;
}


/** gets the value of a constant value expression */
SCIP_Real SCIPgetConsExprExprValueValue(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_Real v;

   assert(sizeof(SCIP_Real) <= sizeof(SCIP_CONSEXPR_EXPRDATA*));

   exprdata = SCIPgetConsExprExprData(expr);
   memcpy(&v, &exprdata, sizeof(SCIP_Real));

   return v;
}
