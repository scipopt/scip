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

/**@file   cons_expr_var.c
 * @brief  variable expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_var.h"

/** the order of two variable is given by their indices
 * @note: this is affected by permutations in the problem! */
static
SCIP_DECL_CONSEXPR_EXPRCMP(compareVar)
{
   int index1;
   int index2;

   index1 = SCIPvarGetIndex(SCIPgetConsExprExprVarVar(expr1));
   index2 = SCIPvarGetIndex(SCIPgetConsExprExprVarVar(expr2));

   return index1 < index2 ? -1 : index1 == index2 ? 0 : 1;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrVar)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrVar(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataVar)
{
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   if( mapvar == NULL )
   {
      /* identical mapping: just copy data pointer */
      assert(targetscip == sourcescip);

      *targetexprdata = SCIPgetConsExprExprData(sourceexpr);
      assert(*targetexprdata != NULL);

      SCIP_CALL( SCIPcaptureVar(targetscip, (SCIP_VAR*)*targetexprdata) );
   }
   else
   {
      /* call mapvar callback (captures targetvar) */
      SCIP_CALL( (*mapvar)(targetscip, (SCIP_VAR**)targetexprdata, sourcescip, (SCIP_VAR*)SCIPgetConsExprExprData(sourceexpr), mapvardata) );
      assert(*targetexprdata != NULL);
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataVar)
{
   SCIP_VAR* var;

   assert(expr != NULL);

   var = (SCIP_VAR*)SCIPgetConsExprExprData(expr);
   assert(var != NULL);

   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printVar)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);

   if( stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR )
   {
      SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName((SCIP_VAR*)SCIPgetConsExprExprData(expr)));
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPREVAL(evalVar)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);

   *val = SCIPgetSolVal(scip, sol, (SCIP_VAR*)SCIPgetConsExprExprData(expr));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalVar)
{
   SCIP_VAR* var;

   assert(expr != NULL);

   var = (SCIP_VAR*) SCIPgetConsExprExprData(expr);
   assert(var != NULL);

   SCIPintervalSetBounds(interval, MIN(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)),
      MAX(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   return SCIP_OKAY;
}

/** creates the handler for variable expression and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "var", "variable", 0, evalVar, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrVar, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataVar, freedataVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalVar) );

   return SCIP_OKAY;
}

/** creates a variable expression */
SCIP_RETCODE SCIPcreateConsExprExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_VAR*             var                 /**< variable to be stored */
   )
{
   assert(consexprhdlr != NULL);
   assert(expr != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcaptureVar(scip, var) );

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPgetConsExprExprHdlrVar(consexprhdlr), (SCIP_CONSEXPR_EXPRDATA*)var, 0, NULL) );

   return SCIP_OKAY;
}

/** gets the variable of a variable expression */
SCIP_VAR* SCIPgetConsExprExprVarVar(
   SCIP_CONSEXPR_EXPR*   expr                /**< variable expression */
   )
{
   assert(expr != NULL);

   return (SCIP_VAR*)SCIPgetConsExprExprData(expr);
}
