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

/**@file   type_cons_expr.h
 * @brief  (public) types of expression constraint
 * @author Stefan Vigerske
 *
 * These are in particular types that define the expressions in cons_expr
 * and that need to be accessed by the linear estimation plugins of cons_expr.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_CONS_EXPR_H__
#define __SCIP_TYPE_CONS_EXPR_H__

#ifdef __cplusplus
extern "C" {
#endif


/** expression operator handler copy method
 *
 * the method includes the expression operator handler into a expression constraint handler
 *
 * This method is usually called when doing a copy of an expression constraint handler.
 *
 *  input:
 *  - scip              : target SCIP main data structure
 *  - consexprhdlr      : target expression constraint handler
 *  - sourcescip        : source SCIP main data structure
 *  - sourceconsexprhdlr : expression constraint handler in source SCIP
 *  - sourceoperandhdlr : expression operand handler in source SCIP
 */
#define SCIP_DECL_CONSEXPR_OPERANDCOPYHDLR(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* consexprhdlr, \
   SCIP* sourcescip, \
   SCIP_CONSHDLR* sourceconsexprhdlr, \
   SCIP_CONSEXPR_OPERANDHDLR* sourceoperandhdlr)

/** expression operator handler free method
 *
 * the method frees the data of an expression operand handler
 *
 *  input:
 *  - scip          : SCIP main data structure
 *  - consexprhdlr  : expression constraint handler
 *  - operandhdlr   : expression operand handler
 *  - operandhdlrdata : operand handler data to be freed
 */
#define SCIP_DECL_CONSEXPR_OPERANDFREEHDLR(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* consexprhdlr, \
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr, \
   SCIP_CONSEXPR_OPERANDHDLRDATA* operandhdlrdata)

/** expression operator data copy method
 *
 * the method copies the data of an expression operand
 *
 * This method is called when creating copies of an expression within
 * the same or between different SCIP instances.
 *
 *  input:
 *  - targetscip         : target SCIP main data structure
 *  - targetconsexprhdlr : expression constraint handler in target SCIP
 *  - targetoperandhdlr  : operand handler in target SCIP
 *  - targetoperanddata  : pointer to store the copied operand data
 *  - sourcescip         : source SCIP main data structure
 *  - sourceconsexprhdlr : expression constraint handler in source SCIP
 *  - sourceoperandhdlr  : expression operand handler in source SCIP
 *  - sourceoperanddata  : operand data (in source SCIP) to be copied
 *  - nchildren          : number of children in the corresponding expression
 */
#define SCIP_DECL_CONSEXPR_OPERANDCOPYDATA(x) SCIP_RETCODE x (\
   SCIP* targetscip, \
   SCIP_CONSHDLR* targetconsexprhdlr, \
   SCIP_CONSEXPR_OPERANDHDLR* targetoperandhdlr, \
   SCIP_CONSEXPR_OPERANDDATA** targetoperanddata, \
   SCIP* sourcescip, \
   SCIP_CONSHDLR* sourceconsexprhdlr, \
   SCIP_CONSEXPR_OPERANDHDLR* sourceoperandhdlr, \
   SCIP_CONSEXPR_OPERANDDATA* sourceoperanddata, \
   int nchildren)

/** expression operator data free method
 *
 * the method frees the data of an expression operand
 *
 *  input:
 *  - scip          : SCIP main data structure
 *  - consexprhdlr  : expression constraint handler
 *  - operandhdlr   : expression operand handler
 *  - operanddata   : operand data to be freed
 *  - nchildren     : number of children in the corresponding expression
 */
#define SCIP_DECL_CONSEXPR_OPERANDFREEDATA(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* consexprhdlr, \
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr, \
   SCIP_CONSEXPR_OPERANDDATA* operanddata, \
   int nchildren)

/** expression operator print method
 *
 * the method prints the data of an expression operand
 *
 * input:
 *  - scip          : SCIP main data structure
 *  - consexprhdlr  : expression constraint handler
 *  - operandhdlr   : expression operand handler
 *  - operanddata   : operand data to print
 *  - nchildren     : number of children in corresponding expression
 *  - file          : the file to print to
 */
#define SCIP_DECL_CONSEXPR_OPERANDPRINT(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* consexprhdlr, \
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr, \
   SCIP_CONSEXPR_OPERANDDATA* operanddata, \
   int nchildren, \
   FILE* file)

/** variability of expression operands */
enum SCIP_ConsExpr_Variability
{
   SCIP_CONSEXPR_INVARIATE    = 0,   /**< no children */
   SCIP_CONSEXPR_UNIVARIATE   = 1,   /**< exactly one child */
   SCIP_CONSEXPR_BIVARIATE    = 2,   /**< exactly two children */
   SCIP_CONSEXPR_MULTIVARIATE = 3    /**< arbitrary number of children */
};

typedef enum  SCIP_ConsExpr_Variability  SCIP_CONSEXPR_VARIABILITY; /**< variability of expression operands */
typedef union SCIP_ConsExpr_Children     SCIP_CONSEXPR_CHILDREN;    /**< storage type for children of an expression */
typedef struct SCIP_ConsExpr_OperandData SCIP_CONSEXPR_OPERANDDATA; /**< expression operand data */
typedef struct SCIP_ConsExpr_OperandHdlr SCIP_CONSEXPR_OPERANDHDLR; /**< expression operand handler */
typedef struct SCIP_ConsExpr_OperandHdlrData SCIP_CONSEXPR_OPERANDHDLRDATA; /**< expression operand handler data */

typedef struct SCIP_ConsExpr_Expr        SCIP_CONSEXPR_EXPR;        /**< expression */


#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_CONS_EXPR_H__ */
