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

/**@file   cons_expr_exp.c
 * @brief  exponential expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_exp.h"

#define EXP_PRECEDENCE  60000

/*
 * Data structures
 */


/*
 * Local methods
 */

/*
 * Callback methods of expression handler
 */


static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrExp)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrExp(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataExp)
{
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPgetConsExprExprData(sourceexpr) == NULL);

   *targetexprdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataExp)
{
   assert(expr != NULL);

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printExp)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print function with opening parenthesis */
         SCIPinfoMessage(scip, file, "exp(");
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

      default: ;
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseExp)
{
   SCIP_CONSEXPR_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseConsExprExpr(scip, consexprhdlr, string, endstring, &childexpr) );
   assert(childexpr != NULL);

   /* create exponential expression */
   SCIP_CALL( SCIPcreateConsExprExprExp(scip, consexprhdlr, expr, childexpr) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the exponential expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPREVAL(evalExp)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID);

   *val = exp(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]));

   return SCIP_OKAY;
}


/** expression hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashExp)
{
   int childhash;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(expr->nchildren == 1);

   expr->hashkey = SCIPcalcFibHash(SCIPgetConsExprExprHdlrPrecedence(expr->exprhdlr));

   childhash = SCIPgetConsExprExprHashkey(SCIPgetConsExprExprChildren(expr)[0]);
   assert(childhash >= 0);

   expr->hashkey ^= childhash;
   assert(expr->hashkey >= 0);

   return SCIP_OKAY;
}


/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalExp)
{
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

   SCIPintervalExp(SCIPinfinity(scip), interval, childinterval);

   return SCIP_OKAY;
}

/** creates the handler for exponential expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrExp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "exp", "exponential function",
         EXP_PRECEDENCE, evalExp, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrExp, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataExp, freedataExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashExp) );

   return SCIP_OKAY;
}

/** creates an exponential expression */
SCIP_RETCODE SCIPcreateConsExprExprExp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< single child */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindConsExprExprHdlr(consexprhdlr, "exp") != NULL);

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPfindConsExprExprHdlr(consexprhdlr, "exp"), NULL, 1, &child) );

   return SCIP_OKAY;
}
