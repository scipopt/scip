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

/** expression handler copy callback
 *
 * the method includes the expression handler into a expression constraint handler
 *
 * This method is usually called when doing a copy of an expression constraint handler.
 *
 *  input:
 *  - scip              : target SCIP main data structure
 *  - consexprhdlr      : target expression constraint handler
 *  - sourceconsexprhdlr : expression constraint handler in source SCIP
 *  - sourceexprhdlr    : expression handler in source SCIP
 *  - valid             : to store indication whether the expression handler was copied
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
 *  - expr : expression to be evaluated
 *  - val : buffer where to store value
 *  - sol : solution that is evaluated (can be NULL)
 */
#define SCIP_DECL_CONSEXPR_EXPREVAL(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_Real* val, \
   SCIP_SOL* sol)

/** derivative evaluation callback
 *
 * The method computes the derivative of an expression using backward automatic differentiation.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to be evaluated
 *  - childidx : index of the child
 *  - val : buffer to store the partial derivative w.r.t. the i-th children
 */
#define SCIP_DECL_CONSEXPR_EXPRBWDIFF(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   int childidx, \
   SCIP_Real* val)

/** expression callback for reverse propagation
 *
 * The method propagates each child of an expression by taking the intervals of all other children into account. The
 * tighter interval is stored inside the interval variable of the corresponding child expression.
 * SCIPtightenConsExprExprInterval() shall be used to tighten a childs interval.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to be evaluated
 *  - reversepropqueue : expression queue in reverse propagation, to be passed on to SCIPtightenConsExprExprInterval
 *  - infeasible: buffer to store whether an expression's bounds were propagated to an empty interval
 *  - nreductions : buffer to store the number of interval reductions of all children
 *  - force : force tightening even if it is below the bound strengthening tolerance
 */
#define SCIP_DECL_CONSEXPR_REVERSEPROP(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_QUEUE* reversepropqueue, \
   SCIP_Bool* infeasible, \
   int* nreductions, \
   SCIP_Bool force)

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
 *  - overestimate    : whether the expression needs to be overestimated
 *  - underestimate   : whether the expression needs to be underestimated
 *
 *  output:
 *  - infeasible      : pointer to store whether an infeasibility was detected while building the LP
 */
#define SCIP_DECL_CONSEXPR_EXPRINITSEPA(x) SCIP_RETCODE x (\
      SCIP* scip, \
      SCIP_CONSHDLR* conshdlr, \
      SCIP_CONSEXPR_EXPR* expr, \
      SCIP_Bool overestimate, \
      SCIP_Bool underestimate, \
      SCIP_Bool* infeasible)

/** separation deinitialization method of an expression handler (called during CONSEXITSOL)
 *
 *  input:
 *  - scip            : SCIP main data structure
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
 *  - overestimate : whether the expression needs to be over- or underestimated
 *  - minviolation : minimal violation of a cut if it should be added to the LP
 *  - result : pointer to store the result
 *  - ncuts : pointer to store the number of added cuts
 */
#define SCIP_DECL_CONSEXPR_EXPRSEPA(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_SOL* sol, \
   SCIP_Bool overestimate, \
   SCIP_Real minviolation, \
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

/** expression branching score callback
 *
 * The method adds branching scores to its children if it finds that the value of the
 * linearization variables does not coincide with the value of the expression in the given solution.
 * It shall use the function SCIPaddConsExprExprBranchScore() to add a branching score to its children.
 * It shall return TRUE in success if no branching is necessary or branching scores have been added.
 * If returning FALSE in success, then other scoring methods will be applied, e.g., a fallback that
 * adds a score to every child.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to be hashed
 *  - sol  : solution (NULL for the LP solution)
 *  - brscoretag : value to be passed on to SCIPaddConsExprExprBranchScore()
 *  - success: buffer to store whether the branching score callback was successful
 */
#define SCIP_DECL_CONSEXPR_EXPRBRANCHSCORE(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_SOL* sol, \
   unsigned int brscoretag, \
   SCIP_Bool* success)

/** expression curvature detection callback
 *
 * The method computes the curvature of an given expression. It assumes that the interval evaluation of the expression
 * has been called before and the expression has been simplified.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - conshdlr : expression constraint handler
 *  - expr : expression to check the curvature for
 *  - curvature : buffer to store the curvature of the expression
 */
#define SCIP_DECL_CONSEXPR_EXPRCURVATURE(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_EXPRCURV* curvature)

/** expression monotonicity detection callback
 *
 * The method computes the monotonicity of an expression with respect to a given child.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to check the monotonicity for
 *  - childidx : index of the considered child expression
 *  - result : buffer to store the monotonicity
 */
#define SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   int childidx, \
   SCIP_MONOTONE* result)

/** expression integrality detection callback
 *
 * The method computes whether an expression evaluates always to an integral value.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to check the integrality for
 *  - isintegral : buffer to store whether expr is integral
 */
#define SCIP_DECL_CONSEXPR_EXPRINTEGRALITY(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_Bool* isintegral)

/* maybe should make this a parameter (was cutmaxrange in other conshdlr)
 * maybe should derive this from the current feastol (e.g., 10/feastol)
 */
#define SCIP_CONSEXPR_CUTMAXRANGE 1.0e7

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
   int                   intvals[2];         /**< two integer values */
   void*                 ptrval;             /**< a pointer */
} SCIP_CONSEXPREXPRWALK_IO;

/** mode for expression iterator */
typedef enum
{
   SCIP_CONSEXPRITERATOR_RTOPOLOGIC,         /**< reverse topological order */
   SCIP_CONSEXPRITERATOR_BFS                 /**< breadth-first search */
} SCIP_CONSEXPRITERATOR_TYPE;

/** monotonicity of an expression */
typedef enum
{
   SCIP_MONOTONE_UNKNOWN      = 0,          /**< unknown */
   SCIP_MONOTONE_INC          = 1,          /**< increasing */
   SCIP_MONOTONE_DEC          = 2,          /**< decreasing */
   SCIP_MONOTONE_CONST        = SCIP_MONOTONE_INC | SCIP_MONOTONE_DEC /**< constant */

} SCIP_MONOTONE;

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

#define SCIP_CONSEXPR_EXPRENFO_NONE           0x0u /**< no enforcement */
#define SCIP_CONSEXPR_EXPRENFO_SEPABELOW      0x1u /**< separation for expr <= auxvar, thus might estimate expr from below */
#define SCIP_CONSEXPR_EXPRENFO_SEPAABOVE      0x2u /**< separation for expr >= auxvar, thus might estimate expr from above */
#define SCIP_CONSEXPR_EXPRENFO_SEPABOTH       (SCIP_CONSEXPR_EXPRENFO_SEPABELOW | SCIP_CONSEXPR_EXPRENFO_SEPAABOVE)  /**< separation for expr == auxvar */
#define SCIP_CONSEXPR_EXPRENFO_INTEVAL        0x4u /**< interval evaluation */
#define SCIP_CONSEXPR_EXPRENFO_REVERSEPROP    0x8u /**< reverse propagation */
#define SCIP_CONSEXPR_EXPRENFO_BRANCHSCORE    0x10u /**< setting branching scores */
#define SCIP_CONSEXPR_EXPRENFO_ALL            (SCIP_CONSEXPR_EXPRENFO_SEPABOTH | SCIP_CONSEXPR_EXPRENFO_INTEVAL | SCIP_CONSEXPR_EXPRENFO_REVERSEPROP | SCIP_CONSEXPR_EXPRENFO_BRANCHSCORE) /**< all enforcement methods */

/** type for exprenfo bitflags */
typedef unsigned int SCIP_CONSEXPR_EXPRENFO_METHOD;

typedef struct SCIP_ConsExpr_ExprData     SCIP_CONSEXPR_EXPRDATA;     /**< expression data */
typedef struct SCIP_ConsExpr_ExprHdlr     SCIP_CONSEXPR_EXPRHDLR;     /**< expression handler */
typedef struct SCIP_ConsExpr_ExprHdlrData SCIP_CONSEXPR_EXPRHDLRDATA; /**< expression handler data */
typedef struct SCIP_ConsExpr_Expr         SCIP_CONSEXPR_EXPR;         /**< expression */
typedef struct SCIP_ConsExpr_ExprEnfo     SCIP_CONSEXPR_EXPRENFO;     /**< expression enforcement data */

typedef struct SCIP_ConsExpr_PrintDotData SCIP_CONSEXPR_PRINTDOTDATA; /**< printing a dot file data */

/** nonlinear handler copy callback
 *
 * the method includes the nonlinear handler into a expression constraint handler
 *
 * This method is usually called when doing a copy of an expression constraint handler.
 *
 *  input:
 *  - targetscip          : target SCIP main data structure
 *  - targetconsexprhdlr  : target expression constraint handler
 *  - sourceconsexprhdlr  : expression constraint handler in source SCIP
 *  - sourcenlhdlr        : nonlinear handler in source SCIP
 */
#define SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(x) SCIP_RETCODE x (\
   SCIP* targetscip, \
   SCIP_CONSHDLR* targetconsexprhdlr, \
   SCIP_CONSHDLR* sourceconsexprhdlr, \
   SCIP_CONSEXPR_NLHDLR* sourcenlhdlr)

/** callback to free data of handler
 *
 * - scip SCIP data structure
 * - nlhdlr nonlinear handler
 * - nlhdlrdata nonlinear handler data to be freed
 */
#define SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(x) SCIP_RETCODE x (\
   SCIP* scip, \
	SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_NLHDLRDATA** nlhdlrdata)

/** callback to free expression specific data
 *
 * - scip SCIP data structure
 * - nlhdlr nonlinear handler
 * - nlhdlrexprdata nonlinear handler expression data to be freed
 */
#define SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata)

/** callback to be called in initialization
 *
 * - scip SCIP data structure
 * - nlhdlr nonlinear handler
 */
#define SCIP_DECL_CONSEXPR_NLHDLRINIT(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr)

/** callback to be called in deinitialization
 *
 * - scip SCIP data structure
 * - nlhdlr nonlinear handler
 */
#define SCIP_DECL_CONSEXPR_NLHDLREXIT(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr)

/** callback to detect structure in expression tree
 *
 * The nonlinear handler shall analyze the current expression and decide whether it wants to contribute
 * in enforcing the relation between this expression (expr) and its auxiliary variable (auxvar).
 * We distinguish the relations expr <= auxvar (denoted as "below") and expr >= auxvar (denoted as "above").
 * Parameters enforcedbelow and enforcedabove indicate on input whether nonlinear handlers for these
 * relations already exist, or none is necessary.
 * Parameter enforcemethods indicates on input which enforcement methods are already provided by some
 * nonlinear handler.
 *
 * If the detect callback decides to become active at an expression, it shall
 * - set enforcedbelow to TRUE if it will enforce expr <= auxvar
 * - set enforcedabove to TRUE if it will enforce expr >= auxvar
 * - signal the enforcement methods it aims to provide by setting the corresponding bit in enforcemethods
 * - set success to TRUE
 *
 * A nonlinear handler can also return TRUE in success if it will not enforce any relation between expr and auxvar.
 * This can be useful for nonlinear handlers that do not implement a complete enforcement, e.g.,
 * a handler that only contributes cutting planes in some situations.
 * Note, that all (non-NULL) enforcement callbacks of the nonlinear handler are potentially called,
 * not only those that are signaled via enforcemethods.
 *
 * A nonlinear handler can still enforce if both enforcedbelow and enforcedabove are TRUE on input.
 * For example, another nonlinear handler may implement propagation and branching, while this handler could
 * provide separation. In this case, the detect callback should update the enforcemethods argument and
 * set success to TRUE.
 *
 * If a nonlinear handler decides to become active in an expression (success == TRUE), then it shall
 * create auxiliary variables for those subexpressions where they will be required.
 *
 * - scip SCIP data structure
 * - conshdlr expr-constraint handler
 * - nlhdlr nonlinear handler
 * - expr expression to analyze
 * - enforcemethods enforcement methods that are provided by some nonlinear handler (to be updated by detect callback)
 * - enforcedbelow indicates whether an enforcement method for expr <= auxvar exists (to be updated by detect callback) or is not necessary
 * - enforcedabove indicates whether an enforcement method for expr >= auxvar exists (to be updated by detect callback) or is not necessary
 * - success buffer to store whether the nonlinear handler should be called for this expression
 * - nlhdlrexprdata nlhdlr's expr data to be stored in expr, can only be set to non-NULL if success is set to TRUE
 */
#define SCIP_DECL_CONSEXPR_NLHDLRDETECT(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_EXPRENFO_METHOD* enforcemethods, \
   SCIP_Bool* enforcedbelow, \
   SCIP_Bool* enforcedabove, \
   SCIP_Bool* success, \
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata)

/** auxiliary evaluation callback of nonlinear handler
 *
 * Evaluates the expression w.r.t. the auxiliary variables that were introduced by the nonlinear handler (if any)
 * The method is used to determine the violation of the relation that the nonlinear
 * handler attempts to enforce. During enforcement, this violation value is used to
 * decide whether separation or branching score callbacks should be called.
 *
 * It can be assumed that the expression itself has been evaluated in the given sol.
 */
#define SCIP_DECL_CONSEXPR_NLHDLREVALAUX(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_Real* auxvalue, \
   SCIP_SOL* sol)


/** separation initialization method of a nonlinear handler (called during CONSINITLP)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : expression constraint handler
 *  - nlhdlr          : nonlinear handler
 *  - nlhdlrexprdata  : exprdata of nonlinear handler
 *  - expr            : expression
 *  - overestimate    : whether the expression needs to be overestimated
 *  - underestimate   : whether the expression needs to be underestimated
 *
 *  output:
 *  - infeasible      : pointer to store whether an infeasibility was detected while building the LP
 */
#define SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(x) SCIP_RETCODE x (\
      SCIP* scip, \
      SCIP_CONSHDLR* conshdlr, \
      SCIP_CONSEXPR_NLHDLR* nlhdlr, \
      SCIP_CONSEXPR_EXPR* expr, \
      SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
      SCIP_Bool overestimate, \
      SCIP_Bool underestimate, \
      SCIP_Bool* infeasible)

/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nlhdlr          : nonlinear handler
 *  - nlhdlrexprdata  : exprdata of nonlinear handler
 *  - expr            : expression
 */
#define SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(x) SCIP_RETCODE x (\
      SCIP* scip, \
      SCIP_CONSEXPR_NLHDLR* nlhdlr, \
      SCIP_CONSEXPR_EXPR* expr, \
      SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata)

/** nonlinear handler separation callback
 *
 * The method tries to separate a given point by means of the nonlinear handler.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - conshdlr : cons expr handler
 *  - nlhdlr : nonlinear handler
 *  - expr : expression
 *  - nlhdlrexprdata : expression specific data of the nonlinear handler
 *  - sol : solution to be separated (NULL for the LP solution)
 *  - auxvalue : current value of expression w.r.t. auxiliary variables as obtained from EVALAUX
 *  - overestimate : whether the expression needs to be over- or underestimated
 *  - minviolation :  minimal violation of a cut if it should be added to the LP
 *  - separated : whether another nonlinear handler already added a cut for this expression
 *  - result : pointer to store the result
 *  - ncuts : pointer to store the number of added cuts
 */
#define SCIP_DECL_CONSEXPR_NLHDLRSEPA(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_SOL* sol, \
   SCIP_Real auxvalue, \
   SCIP_Bool overestimate, \
   SCIP_Real minviolation, \
   SCIP_Bool separated, \
   SCIP_RESULT* result, \
   int* ncuts)

/** nonlinear handler interval evaluation callback
 *
 * The methods computes an interval that contains the image (range) of the expression.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - nlhdlr : nonlinear handler
 *  - expr : expression
 *  - nlhdlrexprdata : expression specific data of the nonlinear handler
 *  - interval : buffer where to store interval (on input: current interval for expr, on output: computed interval for expr)
 *  - varboundrelax : a suggested amount by which to relax variable bounds
 */
#define SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_INTERVAL* interval, \
   SCIP_Real varboundrelax)

/** nonlinear handler callback for reverse propagation
 *
 * The method propagates bounds over the arguments of an expression.
 * The arguments of an expression are other expressions and the tighter intervals should be stored inside the interval variable
 * of the corresponding argument (expression) by using SCIPtightenConsExprExprInterval().
 *
 * input:
 *  - scip : SCIP main data structure
 *  - nlhdlr : nonlinear handler
 *  - expr : expression
 *  - nlhdlrexprdata : expression specific data of the nonlinear handler
 *  - reversepropqueue : expression queue in reverse propagation, to be passed on to SCIPtightenConsExprExprInterval
 *  - infeasible: buffer to store whether an expression's bounds were propagated to an empty interval
 *  - nreductions : buffer to store the number of interval reductions of all children
 *  - force : force tightening even if it is below the bound strengthening tolerance
 */
#define SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_QUEUE* reversepropqueue, \
   SCIP_Bool* infeasible, \
   int* nreductions, \
   SCIP_Bool force)

/** nonlinear handler callback for branching scores
 *
 * The method adds branching scores to successors if it finds that this is how to enforce
 * the relation between the auxiliary variable and the value of the expression in the given solution.
 * It shall use the function SCIPaddConsExprExprBranchScore() to add a branching score to its successors.
 * It shall return TRUE in success if no branching is necessary or branching scores have been added.
 * If returning FALSE in success, then other scoring methods will be applied.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - nlhdlr : nonlinear handler
 *  - expr : expression to be hashed
 *  - nlhdlrexprdata : expression specific data of the nonlinear handler
 *  - sol  : solution (NULL for the LP solution)
 *  - auxvalue : current value of expression w.r.t. auxiliary variables as obtained from EVALAUX
 *  - brscoretag : value to be passed on to SCIPaddConsExprExprBranchScore()
 *  - success: buffer to store whether the branching score callback was successful
 */
#define SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_SOL* sol, \
   SCIP_Real auxvalue, \
   unsigned int brscoretag, \
   SCIP_Bool* success)

typedef struct SCIP_ConsExpr_Nlhdlr      SCIP_CONSEXPR_NLHDLR;        /**< nonlinear handler */
typedef struct SCIP_ConsExpr_NlhdlrData  SCIP_CONSEXPR_NLHDLRDATA;    /**< nonlinear handler data */
typedef struct SCIP_ConsExpr_NlhdlrExprData SCIP_CONSEXPR_NLHDLREXPRDATA;  /**< nonlinear handler data for a specific expression */
typedef struct SCIP_ConsExpr_Iterator    SCIP_CONSEXPR_ITERATOR;      /**< expression tree iterator */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_CONS_EXPR_H__ */
