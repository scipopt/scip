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
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_var.h"

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrVar)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrVar(scip, consexprhdlr) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataVar)
{
   assert(targetscip == sourcescip);  /* TODO if this is a copy from one SCIP to another, we need to get the variable mapping */
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   *targetexprdata = SCIPgetConsExprExprData(sourceexpr);
   assert(*targetexprdata != NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printVar)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);

   SCIPinfoMessage(scip, file, "%s", SCIPvarGetName((SCIP_VAR*)SCIPgetConsExprExprData(expr)));

   return SCIP_OKAY;
}


/** creates the handler for variable expression and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "var", "variable", NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrVar, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataVar, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printVar) );

   return SCIP_OKAY;
}

/** creates the data of a variable expression */
SCIP_RETCODE SCIPcreateConsExprExprVar(
   SCIP*                    scip,            /**< SCIP data structure */
   SCIP_CONSHDLR*           consexprhdlr,    /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprhdlr,        /**< variable expression handler */
   SCIP_CONSEXPR_EXPRDATA** exprdata,        /**< pointer where to store data of expression */
   SCIP_VAR*                var              /**< variable to be stored */
   )
{
   assert(exprdata != NULL);

   *exprdata = (void*)var;

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
