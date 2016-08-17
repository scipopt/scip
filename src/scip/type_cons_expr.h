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
 * @author Benjamin Mueller
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

/** expression simplify callback
 *
 * the method receives the expression to be simplified and a pointer to store the simplified expression
 *
 * input:
 *  - scip           : SCIP main data structure
 *  - expr           : expression to simplify
 * output:
 *  - simplifiedexpr : the simplified expression
 */
#define SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(x) SCIP_RETCODE x (\
   SCIP*                 scip,               \
   SCIP_CONSEXPR_EXPR*   expr,               \
   SCIP_CONSEXPR_EXPR**  simplifiedexpr)

/** expression compare callback
 *
 * the method receives two expressions, expr1 and expr2. Must return
 * -1 if expr1 < expr2
 * 0  if expr1 = expr2
 * 1  if expr1 > expr2
 *
 * input:
 *  - expr1 : first expression to compare
 *  - expr2 : second expression to compare
 */
#define SCIP_DECL_CONSEXPR_EXPRCMP(x) int x (\
   SCIP_CONSEXPR_EXPR* expr1, \
   SCIP_CONSEXPR_EXPR* expr2)

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
   SCIP_CONSHDLR* sourceconsexprhdlr, \
   SCIP_CONSEXPR_EXPRHDLR* sourceexprhdlr, \
   SCIP_Bool* valid)

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

/** variable mapping callback for expression data callback
 *
 * The method maps a variable (in a source SCIP instance) to a variable
 * (in a target SCIP instance) and captures the target variable.
 *
 *  input:
 *  - targetscip         : target SCIP main data structure
 *  - targetvar          : pointer to store the mapped variable
 *  - sourcescip         : source SCIP main data structure
 *  - sourcevar          : variable to be mapped
 *  - mapvardata         : data of callback
 */
#define SCIP_DECL_CONSEXPR_EXPRCOPYDATA_MAPVAR(x) SCIP_RETCODE x (\
   SCIP* targetscip, \
   SCIP_VAR** targetvar, \
   SCIP* sourcescip, \
   SCIP_VAR* sourcevar, \
   void* mapvardata \
   )

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
 *  - targetexprhdlr     : expression handler in target SCIP
 *  - targetexprdata     : pointer to store the copied expression data
 *  - sourcescip         : source SCIP main data structure
 *  - sourceexpr         : expression in source SCIP which data is to be copied,
 *  - mapvar             : variable mapping callback for use by variable expression handler
 *  - mapvardata         : data of variable mapping callback
 */
#define SCIP_DECL_CONSEXPR_EXPRCOPYDATA(x) SCIP_RETCODE x (\
   SCIP* targetscip, \
   SCIP_CONSEXPR_EXPRHDLR* targetexprhdlr, \
   SCIP_CONSEXPR_EXPRDATA** targetexprdata, \
   SCIP* sourcescip, \
   SCIP_CONSEXPR_EXPR* sourceexpr, \
   SCIP_DECL_CONSEXPR_EXPRCOPYDATA_MAPVAR(mapvar), \
   void* mapvardata)

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
 * the method prints an expression
 * it is called during an expression walk at different stages of the walk
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression which data is to be printed
 *  - stage: stage of expression print walk
 *  - file : the file to print to
 */
#define SCIP_DECL_CONSEXPR_EXPRPRINT(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPREXPRWALK_STAGE stage, \
   FILE* file)

/** expression parse callback
 *
 * the method parses an expression
 * it is called when parsing a constraint and an operator with the expr handler name is found
 *
 * input:
 *  - scip         : SCIP main data structure
 *  - consexprhdlr : expression constraint handler
 *  - string       : string containing expression to be parse
 *
 *  output:
 *  - endstring    : pointer to store the position of string after parsing
 *  - expr         : pointer to store the parsed expression
 *  - success      : pointer to store whether the parsing was successful or not
 */
#define SCIP_DECL_CONSEXPR_EXPRPARSE(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* consexprhdlr, \
   const char* string, \
   const char** endstring, \
   SCIP_CONSEXPR_EXPR** expr, \
   SCIP_Bool* success)

/** expression (point-) evaluation callback
 *
 * The method evaluates an expression by taking the values of its children into account.
 * We might extend this later to store (optionally) also information for
 * gradient and Hessian computations.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - val : buffer where to store value
 *  - expr : expression to be evaluated
 *  - sol : solution that is evaluated (can be NULL)
 */
#define SCIP_DECL_CONSEXPR_EXPREVAL(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_Real* val, \
   SCIP_SOL* sol)

/** expression callback for reverse propagation
 *
 * The method propagates each child of an expression by taking the intervals of all other children into account. The
 * tighter interval is stored inside the interval variable of the corresponding child expression.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to be evaluated
 *  - infeasible: buffer to store whether an expression's bounds were propagated to an empty interval
 *  - nreductions : buffer to store the number of interval reductions of all children
 */
#define SCIP_DECL_CONSEXPR_REVERSEPROP(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_Bool* infeasible, \
   int* nreductions)

/** expression (interval-) evaluation callback
 *
 * The method evaluates an expression by taking the intervals of its children into account.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - interval : buffer where to store interval
 *  - expr : expression to be evaluated
 *  - varboundrelax : a suggested amount by which to relax variable bounds
 */
#define SCIP_DECL_CONSEXPR_EXPRINTEVAL(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_INTERVAL* interval, \
   SCIP_Real varboundrelax)

/** separation initialization method of an expression handler (called during CONSINITLP)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : expression constraint handler
 *  - expr            : expression
 *
 *  output:
 *  - infeasible      : pointer to store whether an infeasibility was detected while building the LP
 */
#define SCIP_DECL_CONSEXPR_EXPRINITSEPA(x) SCIP_RETCODE x (\
      SCIP* scip, \
      SCIP_CONSHDLR* conshdlr, \
      SCIP_CONSEXPR_EXPR* expr, \
      SCIP_Bool* infeasible)

/** separation deinitialization method of an expression handler (called during CONSEXITSOL)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : expression constraint handler
 *  - expr            : expression
 */
#define SCIP_DECL_CONSEXPR_EXPREXITSEPA(x) SCIP_RETCODE x (\
      SCIP* scip, \
      SCIP_CONSEXPR_EXPR* expr)

/** expression separation callback
 *
 * The method tries to separate a given point by using linearization variables stored at each expression.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression
 *  - sol  : solution to be separated (NULL for the LP solution)
 *  - minefficacy : minimal efficacy of a cut if it should be added to the LP
 *  - result : pointer to store the result
 *  - ncuts : pointer to store the number of added cuts
 */
#define SCIP_DECL_CONSEXPR_EXPRSEPA(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_SOL* sol, \
   SCIP_Real minefficacy, \
   SCIP_RESULT* result, \
   int* ncuts)

/** expression hash callback
 *
 * The method hashes an expression by taking the hash keys of its children into account.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to be hashed
 *  - expr2key : hash map containing keys for sub-expressions
 *  - hashkey: pointer to store the hash key
 */
#define SCIP_DECL_CONSEXPR_EXPRHASH(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_HASHMAP* expr2key, \
   unsigned int* hashkey)

/** stages of expression walker in which the walker callbacks are called */
typedef enum
{
   SCIP_CONSEXPREXPRWALK_ENTEREXPR,          /**< an expression is visited the first time (before any of its children are visited) */
   SCIP_CONSEXPREXPRWALK_VISITINGCHILD,      /**< a child of an expression is to be visited */
   SCIP_CONSEXPREXPRWALK_VISITEDCHILD,       /**< a child of an expression has been visited */
   SCIP_CONSEXPREXPRWALK_LEAVEEXPR           /**< an expression is to be left (all of its children have been processed) */
} SCIP_CONSEXPREXPRWALK_STAGE;

/** feedback from expression walker callback to expression walker to direct the walk
 *
 * The return code SCIP_CONSEXPREXPRWALK_SKIP is only allowed in the stages SCIP_CONSEXPREXPRWALK_ENTERNODE and SCIP_CONSEXPREXPRWALK_VISITINGCHILD.
 */
typedef enum
{
   SCIP_CONSEXPREXPRWALK_CONTINUE,           /**< continue the walk */
   SCIP_CONSEXPREXPRWALK_SKIP,               /**< skip this node (if in ENTEREXPR stage) or the next child (if in VISITINGCHILD stage) */
   SCIP_CONSEXPREXPRWALK_ABORT               /**< abort the walk */
} SCIP_CONSEXPREXPRWALK_RESULT;

/** user data storage type for expression walker */
typedef union
{
   SCIP_Real             realval;            /**< a floating-point value */
   int                   intval;             /**< an integer value */
   void*                 ptrval;             /**< a pointer */
} SCIP_CONSEXPREXPRWALK_IO;


/** expression graph walk callback
 *
 * input:
 *  - scip   : SCIP main data structure
 *  - expr   : expression node that is visited
 *  - stage  : the current stage of the expression walker
 *  - data   : pointer to user data
 *  - result : buffer to store how to proceed in the walk
 */
#define SCIP_DECL_CONSEXPREXPRWALK_VISIT(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPREXPRWALK_STAGE stage, \
   void* data, \
   SCIP_CONSEXPREXPRWALK_RESULT* result)

/** @name bitflags that customize what is printed by dot-format printer
 * @{
 */
#define SCIP_CONSEXPR_PRINTDOT_EXPRSTRING   0x1u /**< print the math. function that the expression represents (e.g., "c0+c1") */
#define SCIP_CONSEXPR_PRINTDOT_EXPRHDLR     0x2u /**< print expression handler name */
#define SCIP_CONSEXPR_PRINTDOT_NUSES        0x4u /**< print number of uses (reference counting) */
#define SCIP_CONSEXPR_PRINTDOT_NLOCKS       0x8u /**< print number of locks */
#define SCIP_CONSEXPR_PRINTDOT_EVALVALUE   0x10u /**< print evaluation value */
#define SCIP_CONSEXPR_PRINTDOT_EVALTAG     0x30u /**< print evaluation value and tag */
#define SCIP_CONSEXPR_PRINTDOT_INTERVAL    0x40u /**< print interval value */
#define SCIP_CONSEXPR_PRINTDOT_INTERVALTAG 0xC0u /**< print interval value and tag */

/** print everything */
#define SCIP_CONSEXPR_PRINTDOT_ALL SCIP_CONSEXPR_PRINTDOT_EXPRSTRING | SCIP_CONSEXPR_PRINTDOT_EXPRHDLR | SCIP_CONSEXPR_PRINTDOT_NUSES | SCIP_CONSEXPR_PRINTDOT_NLOCKS | SCIP_CONSEXPR_PRINTDOT_EVALTAG | SCIP_CONSEXPR_PRINTDOT_INTERVALTAG

/** type for printdot bitflags
 * @todo find a better name
 */
typedef unsigned int SCIP_CONSEXPR_PRINTDOT_WHAT;
/** @} */

typedef struct SCIP_ConsExpr_ExprData     SCIP_CONSEXPR_EXPRDATA;     /**< expression data */
typedef struct SCIP_ConsExpr_ExprHdlr     SCIP_CONSEXPR_EXPRHDLR;     /**< expression handler */
typedef struct SCIP_ConsExpr_ExprHdlrData SCIP_CONSEXPR_EXPRHDLRDATA; /**< expression handler data */
typedef struct SCIP_ConsExpr_Expr         SCIP_CONSEXPR_EXPR;         /**< expression */

typedef struct SCIP_ConsExpr_PrintDotData SCIP_CONSEXPR_PRINTDOTDATA; /**< printing a dot file data */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_CONS_EXPR_H__ */
