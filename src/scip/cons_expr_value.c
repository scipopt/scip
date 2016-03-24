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
 * @brief  constant value operand handler
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_value.h"

static
SCIP_DECL_CONSEXPR_OPERANDCOPYHDLR(copyhdlrValue)
{
   SCIP_CALL( SCIPincludeOperandHdlrValue(scip, consexprhdlr) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_OPERANDCOPYDATA(copydataValue)
{
   assert(targetoperanddata != NULL);
   assert(sourceexpr != NULL);

   *targetoperanddata = SCIPgetConsExprExprOperatorData(sourceexpr);
   assert(*targetoperanddata != NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_OPERANDPRINT(printValue)
{
   assert(expr != NULL);

   SCIPinfoMessage(scip, file, "%g", SCIPgetOperandValueValue(SCIPgetConsExprExprOperatorData(expr)));

   return SCIP_OKAY;
}


/** creates the handler for constant value operands and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeOperandHdlrValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr;

   SCIP_CALL( SCIPincludeOperandHdlrBasic(scip, consexprhdlr, &ophdlr, "val", "constant value", NULL) );
   assert(ophdlr != NULL);

   SCIP_CALL( SCIPsetOperandHdlrCopyFreeHdlr(scip, consexprhdlr, ophdlr, copyhdlrValue, NULL) );
   SCIP_CALL( SCIPsetOperandHdlrCopyFreeData(scip, consexprhdlr, ophdlr, copydataValue, NULL) );
   SCIP_CALL( SCIPsetOperandHdlrPrint(scip, consexprhdlr, ophdlr, printValue) );

   return SCIP_OKAY;
}

/** creates the data of a constant value operand */
SCIP_RETCODE SCIPcreateOperandValue(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR*  operandhdlr,        /**< constant value operand handler */
   SCIP_CONSEXPR_OPERANDDATA** operanddata,        /**< pointer where to store data of operand */
   SCIP_Real                   value               /**< value to be stored */
   )
{
   assert(operanddata != NULL);
   assert(sizeof(SCIP_Real) <= sizeof(SCIP_CONSEXPR_OPERANDDATA*));

   memcpy(operanddata, &value, sizeof(SCIP_Real));

   return SCIP_OKAY;
}

/** gets the value of a constant value operand */
SCIP_Real SCIPgetOperandValueValue(
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< operand data */
   )
{
   SCIP_Real v;

   assert(sizeof(SCIP_Real) <= sizeof(SCIP_CONSEXPR_OPERANDDATA*));

   memcpy(&v, operanddata, sizeof(SCIP_Real));

   return v;
}
