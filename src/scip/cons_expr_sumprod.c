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

/**@file   cons_expr_sumprod.c
 * @brief  sum and product operand handlers
 * @author Stefan Vigerske
 *
 * Implementation of the sum operator, representing a summation of a constant
 * and the arguments, each multiplied by a coefficients, i.e., sum_i a_i*x_i + constant.
 * Implementation of the product operator, representing a signomial term,
 * i.e., coef * prod_i x_i^e_i.
 * As both operands store similar data, we implement them in the same C file.
 * The data (a_i and constant, or e_i and coef) is currently stored as a SCIP_Real
 * array of length nchildren + 1, storing the constant/coef in the first position,
 * and the a_i/e_i afterwards.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_sumprod.h"

static
SCIP_RETCODE createData(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_OPERANDDATA** operanddata,        /**< pointer where to store data of operand */
   int                         ncoefficients,      /**< number of coefficients (i.e., number of children) */
   SCIP_Real*                  coefficients,       /**< array with coefficients for all operands (or NULL if all 1.0) */
   SCIP_Real                   constant            /**< constant term of sum */
   )
{
   SCIP_Real* opdata;

   assert(operanddata != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &opdata, ncoefficients + 1) );
   opdata[0] = constant;

   if( coefficients != NULL )
   {
      memcpy(opdata+1, coefficients, ncoefficients * sizeof(SCIP_Real));
   }
   else
   {
      int i;
      for( i = 1; i <= ncoefficients; ++i )
         opdata[i] = 1.0;
   }

   *operanddata = (SCIP_CONSEXPR_OPERANDDATA*) opdata;

   return SCIP_OKAY;
}


static
SCIP_DECL_CONSEXPR_OPERANDCOPYHDLR(copyhdlrSum)
{
   SCIP_CALL( SCIPincludeOperandHdlrSum(scip, consexprhdlr) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_OPERANDCOPYHDLR(copyhdlrProduct)
{
   SCIP_CALL( SCIPincludeOperandHdlrProduct(scip, consexprhdlr) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_OPERANDCOPYDATA(copydataSumProduct)
{
   SCIP_Real* opdata;

   assert(targetoperanddata != NULL);

   SCIP_CALL( SCIPduplicateBlockMemoryArray(targetscip, &opdata, (SCIP_Real*)sourceoperanddata, nchildren + 1) );

   *targetoperanddata = (SCIP_CONSEXPR_OPERANDDATA*)opdata;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_OPERANDFREEDATA(freedataSumProduct)
{
   SCIP_Real* opdata;

   assert(operanddata != NULL);
   opdata = (SCIP_Real*)operanddata;

   SCIPfreeBlockMemoryArray(scip, &opdata, nchildren + 1);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_OPERANDPRINT(printSumProduct)
{
   int i;

   assert(operanddata != NULL);

   SCIPinfoMessage(scip, file, "%s[%g", SCIPgetOperandHdlrName(operandhdlr), *(SCIP_Real*)operanddata);

   for( i = 0; i < nchildren; ++i )
   {
      SCIPinfoMessage(scip, file, "%c%g", i ? ',' : ';', ((SCIP_Real*)operanddata)[i+1]);
   }

   SCIPinfoMessage(scip, file, "]");

   return SCIP_OKAY;
}


/** creates the handler for sum operands and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeOperandHdlrSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr;

   SCIP_CALL( SCIPincludeOperandHdlrBasic(scip, consexprhdlr, &ophdlr, "sum", "summation with coefficients and a constant", NULL) );
   assert(ophdlr != NULL);

   SCIP_CALL( SCIPsetOperandHdlrCopyFreeHdlr(scip, consexprhdlr, ophdlr, copyhdlrSum, NULL) );
   SCIP_CALL( SCIPsetOperandHdlrCopyFreeData(scip, consexprhdlr, ophdlr, copydataSumProduct, freedataSumProduct) );
   SCIP_CALL( SCIPsetOperandHdlrPrint(scip, consexprhdlr, ophdlr, printSumProduct) );

   return SCIP_OKAY;
}

/** creates the data of a summation operand */
SCIP_RETCODE SCIPcreateOperandSum(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR*  operandhdlr,        /**< variable operand handler */
   SCIP_CONSEXPR_OPERANDDATA** operanddata,        /**< pointer where to store data of operand */
   int                         ncoefficients,      /**< number of coefficients (i.e., number of children) */
   SCIP_Real*                  coefficients,       /**< array with coefficients for all operands (or NULL if all 1.0) */
   SCIP_Real                   constant            /**< constant term of sum */
   )
{
   SCIP_CALL( createData(scip, operanddata, ncoefficients, coefficients, constant) );

   return SCIP_OKAY;
}

/** gets the coefficients of a summation operand */
EXTERN
SCIP_Real* SCIPgetOperandSumCoefs(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< sum operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< data of operand */
   )
{
   assert(operanddata != NULL);

   return ((SCIP_Real*)operanddata) + 1;
}

/** gets the constant of a summation operand */
EXTERN
SCIP_Real SCIPgetOperandSumConstant(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< sum operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< data of operand */
   )
{
   assert(operanddata != NULL);

   return *(SCIP_Real*)operanddata;
}


/** creates the handler for product operands and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeOperandHdlrProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr;

   SCIP_CALL( SCIPincludeOperandHdlrBasic(scip, consexprhdlr, &ophdlr, "prod", "product of children with exponents (actually a signomial)", NULL) );
   assert(ophdlr != NULL);

   SCIP_CALL( SCIPsetOperandHdlrCopyFreeHdlr(scip, consexprhdlr, ophdlr, copyhdlrProduct, NULL) );
   SCIP_CALL( SCIPsetOperandHdlrCopyFreeData(scip, consexprhdlr, ophdlr, copydataSumProduct, freedataSumProduct) );
   SCIP_CALL( SCIPsetOperandHdlrPrint(scip, consexprhdlr, ophdlr, printSumProduct) );

   return SCIP_OKAY;
}

/** creates the data of a product operand */
SCIP_RETCODE SCIPcreateOperandProduct(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR*  operandhdlr,        /**< variable operand handler */
   SCIP_CONSEXPR_OPERANDDATA** operanddata,        /**< pointer where to store data of operand */
   int                         nexponents,         /**< number of exponents (i.e., number of children) */
   SCIP_Real*                  exponents,          /**< array with exponents for all operands (or NULL if all 1.0) */
   SCIP_Real                   constant            /**< constant coefficient of product */
   )
{
   SCIP_CALL( createData(scip, operanddata, nexponents, exponents, constant) );

   return SCIP_OKAY;
}

/** gets the exponents of a product operand */
SCIP_Real* SCIPgetOperandProductExponents(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< product operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< data of operand */
   )
{
   assert(operanddata != NULL);

   return ((SCIP_Real*)operanddata) + 1;
}

/** gets the constant coefficient of a product operand */
SCIP_Real SCIPgetOperandProductCoef(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< product operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< data of operand */
   )
{
   assert(operanddata != NULL);

   return *(SCIP_Real*)operanddata;
}
