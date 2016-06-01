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

/**@file   cons_expr_constant.c
 * @brief  constant value expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_value.h"

#define VALUE_PRECEDENCE     10000
#define VALUE_HASHKEY        SCIPcalcFibHash(36787)

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrValue)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrValue(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataValue)
{
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   *targetexprdata = SCIPgetConsExprExprData(sourceexpr);
   assert(*targetexprdata != NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printValue)
{
   assert(expr != NULL);

   if( stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR )
   {
      SCIP_Real v = SCIPgetConsExprExprValueValue(expr);
      if( v < 0.0 && VALUE_PRECEDENCE <= SCIPgetConsExprExprWalkParentPrecedence(expr)  )
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
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   memcpy(val, &exprdata, sizeof(SCIP_Real));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalValue)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_Real val;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   memcpy(&val, &exprdata, sizeof(SCIP_Real));

   SCIPintervalSet(interval, val);

   return SCIP_OKAY;
}

/** expression reverse propagaton callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropValue)
{
   SCIPerrorMessage("Unexpected call of reverse propagation callback of value expression handler.");
   SCIPABORT();

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashValue)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 0);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   *hashkey = VALUE_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(SCIPgetConsExprExprValueValue(expr));

   return SCIP_OKAY;
}

/** creates the handler for constant value expression and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "val", "constant value", VALUE_PRECEDENCE, evalValue, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrValue, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataValue, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashValue) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropValue) );

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
