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
 * @brief  sum and product expression handlers
 * @author Stefan Vigerske
 *
 * Implementation of the sum expression, representing a summation of a constant
 * and the arguments, each multiplied by a coefficients, i.e., sum_i a_i*x_i + constant.
 * Implementation of the product expression, representing a signomial term,
 * i.e., coef * prod_i x_i^e_i.
 * As both expressions store similar data, we implement them in the same C file.
 * The data (a_i and constant, or e_i and coef) is currently stored as a SCIP_Real
 * array of length nchildren + 1, storing the constant/coef in the first position,
 * and the a_i/e_i afterwards.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_sumprod.h"

static
SCIP_RETCODE createData(
   SCIP*                    scip,            /**< SCIP data structure */
   SCIP_CONSEXPR_EXPRDATA** exprdata,        /**< pointer where to store expression data */
   int                      ncoefficients,   /**< number of coefficients (i.e., number of children) */
   SCIP_Real*               coefficients,    /**< array with coefficients for all children (or NULL if all 1.0) */
   SCIP_Real                constant         /**< constant term of sum */
   )
{
   SCIP_Real* edata;

   assert(exprdata != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &edata, ncoefficients + 1) );
   edata[0] = constant;

   if( coefficients != NULL )
   {
      memcpy(edata+1, coefficients, ncoefficients * sizeof(SCIP_Real));
   }
   else
   {
      int i;
      for( i = 1; i <= ncoefficients; ++i )
         edata[i] = 1.0;
   }

   *exprdata = (SCIP_CONSEXPR_EXPRDATA*) edata;

   return SCIP_OKAY;
}


static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrSum)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrSum(scip, consexprhdlr) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrProduct)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrProduct(scip, consexprhdlr) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataSumProduct)
{
   SCIP_Real* sourceexprdata;
   SCIP_Real* exprdata;

   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   sourceexprdata = (SCIP_Real*)SCIPgetConsExprExprData(sourceexpr);
   assert(sourceexprdata != NULL);

   SCIP_CALL( SCIPduplicateBlockMemoryArray(targetscip, &exprdata, sourceexprdata, SCIPgetConsExprExprNChildren(sourceexpr) + 1) );

   *targetexprdata = (SCIP_CONSEXPR_EXPRDATA*)exprdata;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataSumProduct)
{
   SCIP_Real* exprdata;

   assert(expr != NULL);

   exprdata = (SCIP_Real*)SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIPfreeBlockMemoryArray(scip, &exprdata, SCIPgetConsExprExprNChildren(expr) + 1);

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printSum)
{
   SCIP_Real* exprdata;

   assert(expr != NULL);

   exprdata = (SCIP_Real*)SCIPgetConsExprExprData(expr);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print constant, if nonzero */
         if( exprdata[0] != 0.0 )
         {
            SCIPinfoMessage(scip, file, "%g", exprdata[SCIPgetConsExprExprWalkCurrentChild(expr)]);
         }
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      {
         /* print coefficient (with plus or minus sign) and opening parenthesis */
         SCIPinfoMessage(scip, file, "%+g * (", exprdata[SCIPgetConsExprExprWalkCurrentChild(expr)]);
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
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
SCIP_DECL_CONSEXPR_EXPRPRINT(printProduct)
{
   SCIP_Real* exprdata;

   assert(expr != NULL);

   exprdata = (SCIP_Real*)SCIPgetConsExprExprData(expr);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print constant, if nonzero */
         if( exprdata[0] != 0.0 )
         {
            SCIPinfoMessage(scip, file, "%g", exprdata[SCIPgetConsExprExprWalkCurrentChild(expr)]);
         }
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      {
         /* print coefficient (with plus or minus sign) and opening parenthesis */
         SCIPinfoMessage(scip, file, " * (%g) * (", exprdata[SCIPgetConsExprExprWalkCurrentChild(expr)]);
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      {
         /* print closing parenthesis */
         SCIPinfoMessage(scip, file, ")");
         break;
      }

      default: ;
   }

   return SCIP_OKAY;
}

/** creates the handler for sum expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "sum", "summation with coefficients and a constant", 100000, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrSum, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataSumProduct, freedataSumProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printSum) );

   return SCIP_OKAY;
}

/** creates a sum expression */
SCIP_RETCODE SCIPcreateConsExprExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR**  children,           /**< children */
   SCIP_Real*            coefficients,       /**< array with coefficients for all children (or NULL if all 1.0) */
   SCIP_Real             constant            /**< constant term of sum */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   SCIP_CALL( createData(scip, &exprdata, nchildren, coefficients, constant) );

   SCIP_CALL( SCIPcreateConsExprExpr(scip, consexprhdlr, expr, SCIPgetConsExprExprHdlrSum(consexprhdlr), exprdata, nchildren, children) );

   return SCIP_OKAY;
}

/** gets the coefficients of a summation expression */
SCIP_Real* SCIPgetConsExprExprSumCoefs(
   SCIP_CONSEXPR_EXPR*   expr                /**< sum expression */
   )
{
   assert(expr != NULL);

   return ((SCIP_Real*)SCIPgetConsExprExprData(expr)) + 1;
}

/** gets the constant of a summation expression */
SCIP_Real SCIPgetConsExprExprSumConstant(
   SCIP_CONSEXPR_EXPR*   expr                /**< sum expression */
   )
{
   assert(expr != NULL);

   return *(SCIP_Real*)SCIPgetConsExprExprData(expr);
}


/** creates the handler for product expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "prod", "product of children with exponents (actually a signomial)", 50000, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrProduct, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataSumProduct, freedataSumProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printProduct) );

   return SCIP_OKAY;
}

/** creates a product expression */
SCIP_RETCODE SCIPcreateConsExprExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR**  children,           /**< children */
   SCIP_Real*            exponents,          /**< array with exponents for all children (or NULL if all 1.0) */
   SCIP_Real             constant            /**< constant coefficient of product */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   SCIP_CALL( createData(scip, &exprdata, nchildren, exponents, constant) );

   SCIP_CALL( SCIPcreateConsExprExpr(scip, consexprhdlr, expr, SCIPgetConsExprExprHdlrProduct(consexprhdlr), exprdata, nchildren, children) );

   return SCIP_OKAY;
}

/** gets the exponents of a product expression */
SCIP_Real* SCIPgetConsExprExprProductExponents(
   SCIP_CONSEXPR_EXPR*   expr                /**< product expression */
   )
{
   assert(expr != NULL);

   return ((SCIP_Real*)SCIPgetConsExprExprData(expr)) + 1;
}

/** gets the constant coefficient of a product expression */
SCIP_Real SCIPgetConsExprExprProductCoef(
   SCIP_CONSEXPR_EXPR*   expr                /**< product expression */
   )
{
   assert(expr != NULL);

   return *(SCIP_Real*)SCIPgetConsExprExprData(expr);
}
