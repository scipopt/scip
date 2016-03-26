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


/** expression handler copy callback
 *
 * the method includes the expression handler into a expression constraint handler
 *
 * This method is usually called when doing a copy of an expression constraint handler.
 *
 *  input:
 *  - scip              : target SCIP main data structure
 *  - consexprhdlr      : target expression constraint handler
 *  - sourcescip        : source SCIP main data structure
 *  - sourceconsexprhdlr : expression constraint handler in source SCIP
 *  - sourceexprhdlr    : expression handler in source SCIP
 */
#define SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* consexprhdlr, \
   SCIP* sourcescip, \
   SCIP_CONSHDLR* sourceconsexprhdlr, \
   SCIP_CONSEXPR_EXPRHDLR* sourceexprhdlr)

/** expression handler free callback
 *
 * the callback frees the data of an expression handler
 *
 *  input:
 *  - scip          : SCIP main data structure
 *  - consexprhdlr  : expression constraint handler
 *  - exprhdlr      : expression handler
 *  - exprhdlrdata  : expression handler data to be freed
 */
#define SCIP_DECL_CONSEXPR_EXPRFREEHDLR(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* consexprhdlr, \
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr, \
   SCIP_CONSEXPR_EXPRHDLRDATA** exprhdlrdata)

/** expression data copy callback
 *
 * the method copies the data of an expression
 *
 * This method is called when creating copies of an expression within
 * the same or between different SCIP instances. It is given the
 * source expression which data shall be copied. It expects
 * that *targetexprdata will be set. This data will then be used
 * to create a new expression.
 *
 *  input:
 *  - targetscip         : target SCIP main data structure
 *  - targetconsexprhdlr : expression constraint handler in target SCIP
 *  - targetexprhdlr     : expression handler in target SCIP
 *  - targetexprdata     : pointer to store the copied expression data
 *  - sourcescip         : source SCIP main data structure
 *  - sourceconsexprhdlr : expression constraint handler in source SCIP
 *  - sourceexpr         : expression in source SCIP which data is to be copied
 */
#define SCIP_DECL_CONSEXPR_EXPRCOPYDATA(x) SCIP_RETCODE x (\
   SCIP* targetscip, \
   SCIP_CONSHDLR* targetconsexprhdlr, \
   SCIP_CONSEXPR_EXPRHDLR* targetexprhdlr, \
   SCIP_CONSEXPR_EXPRDATA** targetexprdata, \
   SCIP* sourcescip, \
   SCIP_CONSHDLR* sourceconsexprhdlr, \
   SCIP_CONSEXPR_EXPR* sourceexpr)

/** expression data free callback
 *
 * The method frees the data of an expression.
 * It assumes that expr->exprdata will be set to NULL.
 *
 *  input:
 *  - scip          : SCIP main data structure
 *  - expr          : the expression which data to be freed
 */
#define SCIP_DECL_CONSEXPR_EXPRFREEDATA(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr)

/** expression print callback
 *
 * the method prints the data of an expression
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression which data is to be printed
 *  - file : the file to print to
 */
#define SCIP_DECL_CONSEXPR_EXPRPRINT(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   FILE* file)


typedef struct SCIP_ConsExpr_ExprData     SCIP_CONSEXPR_EXPRDATA;     /**< expression data */
typedef struct SCIP_ConsExpr_ExprHdlr     SCIP_CONSEXPR_EXPRHDLR;     /**< expression handler */
typedef struct SCIP_ConsExpr_ExprHdlrData SCIP_CONSEXPR_EXPRHDLRDATA; /**< expression handler data */
typedef struct SCIP_ConsExpr_Expr         SCIP_CONSEXPR_EXPR;         /**< expression */


#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_CONS_EXPR_H__ */
