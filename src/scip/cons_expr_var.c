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
 * @brief  variable operand handler
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_var.h"

static
SCIP_DECL_CONSEXPR_OPERANDCOPYHDLR(copyhdlrVar)
{
   SCIP_CALL( SCIPincludeExprOperandVar(scip, consexprhdlr) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_OPERANDCOPYDATA(copydataVar)
{
   assert(targetscip == sourcescip);  /* if this is a copy from one SCIP to another, we need to get the variable mapping */
   assert(targetoperanddata != NULL);

   *targetoperanddata = sourceoperanddata;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_OPERANDPRINT(printVar)
{
   assert(operanddata != NULL);

   SCIPinfoMessage(scip, file, "%s", SCIPvarGetName((SCIP_VAR*)(operanddata)));

   return SCIP_OKAY;
}


/** creates the handler for variable operands and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeExprOperandVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CALL( SCIPincludeOperandHdlrConshdlrExpr(
      scip, consexprhdlr, "var", "variable",
      copyhdlrVar, NULL, copydataVar, NULL, printVar,
      NULL) );

   return SCIP_OKAY;
}

/** creates the data of a variable operand */
SCIP_RETCODE SCIPcreateExprOperandVar(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR*  operandhdlr,        /**< variable operand handler */
   SCIP_CONSEXPR_OPERANDDATA** operanddata,        /**< pointer where to store data of operand */
   SCIP_VAR*                   var                 /**< variable to be stored */
   )
{
   assert(operanddata != NULL);

   *operanddata = (void*)var;

   return SCIP_OKAY;
}
