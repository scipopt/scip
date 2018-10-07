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

/* maybe should make this a parameter (was cutmaxrange in other conshdlr)
 * maybe should derive this from the current feastol (e.g., 10/feastol)
 */
#define SCIP_CONSEXPR_CUTMAXRANGE 1.0e7

typedef struct SCIP_ConsExpr_ExprData  SCIP_CONSEXPR_EXPRDATA;     /**< expression data */
typedef struct SCIP_ConsExpr_Expr      SCIP_CONSEXPR_EXPR;         /**< expression */

/** monotonicity of an expression */
typedef enum
{
   SCIP_MONOTONE_UNKNOWN      = 0,          /**< unknown */
   SCIP_MONOTONE_INC          = 1,          /**< increasing */
   SCIP_MONOTONE_DEC          = 2,          /**< decreasing */
   SCIP_MONOTONE_CONST        = SCIP_MONOTONE_INC | SCIP_MONOTONE_DEC /**< constant */

} SCIP_MONOTONE;

/** bilinear term data */
struct SCIP_ConsExpr_BilinTerm
{
   SCIP_VAR*             x;                  /**< first variable */
   SCIP_VAR*             y;                  /**< second variable */
   SCIP_VAR*             auxvar;             /**< auxiliary variable for the product of x and y */
   int                   nlockspos;          /**< number of positive expression locks */
   int                   nlocksneg;          /**< number of negative expression locks */
};
typedef struct SCIP_ConsExpr_BilinTerm SCIP_CONSEXPR_BILINTERM;    /**< bilinear term data */

/** callback that returns bounds for a given variable as used in interval evaluation
 *
 * Implements a relaxation scheme for variable bounds and translates between different infinity values.
 *
 *  input:
 *  - scip           : SCIP main data structure
 *  - var            : variable for which to obtain bounds
 *  - intevalvardata : data that belongs to this callback
 *  output:
 *  - returns an interval that contains the current variable bounds, but might be (slightly) larger
 */
#define SCIP_DECL_CONSEXPR_INTEVALVAR(x) SCIP_INTERVAL x (\
   SCIP* scip, \
   SCIP_VAR* var, \
   void* intevalvardata \
   )

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
#define SCIP_DECL_CONSEXPR_MAPVAR(x) SCIP_RETCODE x (\
   SCIP* targetscip, \
   SCIP_VAR** targetvar, \
   SCIP* sourcescip, \
   SCIP_VAR* sourcevar, \
   void* mapvardata \
   )

/**@name Expression Handler */
/**@{ */

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
 * This callback must be implemented for expressions that have data.
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
   SCIP_DECL_CONSEXPR_MAPVAR(mapvar), \
   void* mapvardata)

/** expression data free callback
 *
 * The method frees the data of an expression.
 * It assumes that expr->exprdata will be set to NULL.
 *
 * This callback must be implemented for expressions that have data.
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
 * it is called while iterating over the expression graph at different stages
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression which data is to be printed
 *  - stage: stage of expression graph iteration
 *  - currentchild: index of current child if in stage visitingchild or visitedchild
 *  - parentprecedence: precedence of parent
 *  - file : the file to print to
 */
#define SCIP_DECL_CONSEXPR_EXPRPRINT(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPRITERATOR_STAGE stage, \
   int currentchild, \
   unsigned int parentprecedence, \
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
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr, \
   const char* string, \
   const char** endstring, \
   SCIP_CONSEXPR_EXPR** expr, \
   SCIP_Bool* success)

/** expression curvature detection callback
 *
 * The method returns whether an expression can have a desired curvature under conditions on the
 * curvature of the children.
 * That is, the method shall return TRUE in success and requirements on the curvature for each child
 * which will suffice for this expression to be convex (or concave, or linear, as specified by caller)
 * w.r.t. the current activities of all children.
 * It can return "unknown" for a child's curvature if its curvature does not matter (though that's
 * rarely the case).
 *
 * The method assumes that the activity evaluation of the expression has been called before
 * and the expression has been simplified.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - conshdlr : expression constraint handler
 *  - expr : expression to check the curvature for
 *  - exprcurvature : desired curvature of this expression
 *  - success: buffer to store whether the desired curvature be obtained
 *  - childcurv: array to store required curvature for each child
 */
#define SCIP_DECL_CONSEXPR_EXPRCURVATURE(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_EXPRCURV exprcurvature, \
   SCIP_Bool* success, \
   SCIP_EXPRCURV* childcurv )

/** expression monotonicity detection callback
 *
 * The method computes the monotonicity of an expression with respect to a given child.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - conshdlr: cons_expr constraint handler
 *  - expr : expression to check the monotonicity for
 *  - childidx : index of the considered child expression
 *  - result : buffer to store the monotonicity
 */
#define SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
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

/** expression hash callback
 *
 * The method hashes an expression by taking the hashes of its children into account.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to be hashed
 *  - hashkey : pointer to store the hash value
 *  - childrenhashes : array with hash values of children
 */
#define SCIP_DECL_CONSEXPR_EXPRHASH(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   unsigned int* hashkey, \
   unsigned int* childrenhashes)

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
#define SCIP_DECL_CONSEXPR_EXPRCOMPARE(x) int x (\
   SCIP_CONSEXPR_EXPR* expr1, \
   SCIP_CONSEXPR_EXPR* expr2)

/** derivative evaluation callback
 *
 * The method computes the partial derivative of expr w.r.t its child at childidx.
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

/** derivative evaluation callback for forward mode differentiation
 *
 * The method evaluates the directional derivative of an expression by taking the directional derivative of its children into account and their values.
 * Equivalently, it computes the total derivative, w.r.t its children, of expr.
 *
 * The directional derivative of a child can be accessed via SCIPgetConsExprExprDot
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to be evaluated
 *  - dot : buffer to store derivative value
 *  - direction : direction of the derivative (useful only for var expressions)
 */
#define SCIP_DECL_CONSEXPR_EXPRFWDIFF(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_Real* dot, \
   SCIP_SOL* direction)

/** derivative evaluation callback for hessian directions (backward over forward)
 *
 * The method computes the total derivative, w.r.t its children, of the partial derivative of expr w.r.t childidx
 * Equivalently, it computes the partial derivative w.r.t childidx of the total derivative
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression to be evaluated
 *  - childidx : index of the child
 *  - bardot : buffer to store derivative value
 *  - direction : direction of the derivative (useful only for var expressions)
 */
#define SCIP_DECL_CONSEXPR_EXPRBWFWDIFF(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   int childidx, \
   SCIP_Real* bardot, \
   SCIP_SOL* direction)

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

/** expression (interval-) evaluation callback
 *
 * The method evaluates an expression by taking the intervals of its children into account.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - interval : buffer where to store interval
 *  - expr : expression to be evaluated
 *  - intevalvar : callback to be called when interval evaluating a variable
 *  - intevalvardata : data to be passed to intevalvar callback
 */
#define SCIP_DECL_CONSEXPR_EXPRINTEVAL(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_INTERVAL* interval, \
   SCIP_DECL_CONSEXPR_INTEVALVAR((*intevalvar)), \
   void* intevalvardata)

/** expression under/overestimation callback
 *
 * The method tries to compute a linear under- or overestimator that is as tight as possible
 * at a given point by using auxiliary variables stored in all children.
 * If the value of the estimator in the solution is smaller (larger) than targetvalue
 * when underestimating (overestimating), then no estimator needs to be computed.
 * Note, that targetvalue can be infinite if any estimator will be accepted.
 * If successful, it shall store the coefficient of the i-th child in entry coefs[i] and
 * the constant part in \par constant.
 * The callback shall set branchcand[i] to FALSE if branching in the i-th child would not
 * improve the estimator. That is, branchcand[i] will be initialized to TRUE for all children.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - expr : expression
 *  - sol  : solution at which to estimate (NULL for the LP solution)
 *  - overestimate : whether the expression needs to be over- or underestimated
 *  - targetvalue : a value that the estimator shall exceed, can be +/-infinity
 *  - coefs : array to store coefficients of estimator
 *  - constant : buffer to store constant part of estimator
 *  - islocal : buffer to store whether estimator is valid locally only
 *  - success : buffer to indicate whether an estimator could be computed
 *  - branchcand: array to indicate which children (not) to consider for branching
 */
#define SCIP_DECL_CONSEXPR_EXPRESTIMATE(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_SOL* sol, \
   SCIP_Bool overestimate, \
   SCIP_Real targetvalue, \
   SCIP_Real* coefs, \
   SCIP_Real* constant, \
   SCIP_Bool* islocal, \
   SCIP_Bool* success, \
   SCIP_Bool* branchcand)

/** expression simplify callback
 *
 * the method receives the expression to be simplified and a pointer to store the simplified expression
 *
 * input:
 *  - scip           : SCIP main data structure
 *  - consexprhdlr   : expression constraint handler
 *  - expr           : expression to simplify
 * output:
 *  - simplifiedexpr : the simplified expression
 */
#define SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(x) SCIP_RETCODE x (\
   SCIP*                 scip,               \
   SCIP_CONSHDLR*        conshdlr,           \
   SCIP_CONSEXPR_EXPR*   expr,               \
   SCIP_CONSEXPR_EXPR**  simplifiedexpr)

/** expression callback for reverse propagation
 *
 * The method propagates given bounds over the children of an expression.
 * The tighter interval should be passed to the corresponding child expression by using
 * SCIPtightenConsExprExprInterval().
 *
 * input:
 *  - scip : SCIP main data structure
 *  - conshdlr: expr constraint handler
 *  - expr : expression
 *  - bounds : the bounds on the expression that should be propagated
 *  - infeasible: buffer to store whether an expression's bounds were propagated to an empty interval
 *  - nreductions : buffer to store the number of interval reductions of all children
 */
#define SCIP_DECL_CONSEXPR_EXPRREVERSEPROP(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_INTERVAL bounds, \
   SCIP_Bool* infeasible, \
   int* nreductions )

/** separation initialization method of an expression handler (called during CONSINITLP)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : expression constraint handler
 *  - cons            : expression constraint
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
      SCIP_CONS* cons, \
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

typedef struct SCIP_ConsExpr_ExprHdlr     SCIP_CONSEXPR_EXPRHDLR;     /**< expression handler */
typedef struct SCIP_ConsExpr_ExprHdlrData SCIP_CONSEXPR_EXPRHDLRDATA; /**< expression handler data */

/** @} */  /* expression handler */


/** @name expression iterator
 * @{
 */

/** maximal number of iterators that can be active on an expression graph concurrently
 *
 * How often an expression graph iteration can be started within an active iteration, plus one.
 */
#define SCIP_CONSEXPRITERATOR_MAXNACTIVE 5

/** stages of expression DFS iteration */
#define SCIP_CONSEXPRITERATOR_ENTEREXPR     1u /**< an expression is visited the first time (before any of its children are visited) */
#define SCIP_CONSEXPRITERATOR_VISITINGCHILD 2u /**< a child of an expression is to be visited */
#define SCIP_CONSEXPRITERATOR_VISITEDCHILD  4u /**< a child of an expression has been visited */
#define SCIP_CONSEXPRITERATOR_LEAVEEXPR     8u /**< an expression is to be left (all of its children have been processed) */
#define SCIP_CONSEXPRITERATOR_ALLSTAGES     (SCIP_CONSEXPRITERATOR_ENTEREXPR | SCIP_CONSEXPRITERATOR_VISITINGCHILD | SCIP_CONSEXPRITERATOR_VISITEDCHILD | SCIP_CONSEXPRITERATOR_LEAVEEXPR)

/** type to represent stage of DFS iterator */
typedef unsigned int SCIP_CONSEXPRITERATOR_STAGE;

/** user data storage type for expression iteration */
typedef union
{
   SCIP_Real             realval;            /**< a floating-point value */
   int                   intval;             /**< an integer value */
   int                   intvals[2];         /**< two integer values */
   unsigned int          uintval;            /**< an unsigned integer value */
   void*                 ptrval;             /**< a pointer */
} SCIP_CONSEXPRITERATOR_USERDATA;

/** mode for expression iterator */
typedef enum
{
   SCIP_CONSEXPRITERATOR_RTOPOLOGIC,         /**< reverse topological order */
   SCIP_CONSEXPRITERATOR_BFS,                /**< breadth-first search */
   SCIP_CONSEXPRITERATOR_DFS                 /**< depth-first search */
} SCIP_CONSEXPRITERATOR_TYPE;

typedef struct SCIP_ConsExpr_Expr_IterData SCIP_CONSEXPR_EXPR_ITERDATA; /**< expression tree iterator data for a specific expression */
typedef struct SCIP_ConsExpr_Iterator      SCIP_CONSEXPR_ITERATOR;      /**< expression tree iterator */

/** @} */

/** @name expression printing
 * @{
 */

#define SCIP_CONSEXPR_PRINTDOT_EXPRSTRING   0x1u /**< print the math. function that the expression represents (e.g., "c0+c1") */
#define SCIP_CONSEXPR_PRINTDOT_EXPRHDLR     0x2u /**< print expression handler name */
#define SCIP_CONSEXPR_PRINTDOT_NUSES        0x4u /**< print number of uses (reference counting) */
#define SCIP_CONSEXPR_PRINTDOT_NLOCKS       0x8u /**< print number of locks */
#define SCIP_CONSEXPR_PRINTDOT_EVALVALUE   0x10u /**< print evaluation value */
#define SCIP_CONSEXPR_PRINTDOT_EVALTAG     0x30u /**< print evaluation value and tag */
#define SCIP_CONSEXPR_PRINTDOT_ACTIVITY    0x40u /**< print activity value */
#define SCIP_CONSEXPR_PRINTDOT_ACTIVITYTAG 0xC0u /**< print activity value and corresponding tag */

/** print everything */
#define SCIP_CONSEXPR_PRINTDOT_ALL SCIP_CONSEXPR_PRINTDOT_EXPRSTRING | SCIP_CONSEXPR_PRINTDOT_EXPRHDLR | SCIP_CONSEXPR_PRINTDOT_NUSES | SCIP_CONSEXPR_PRINTDOT_NLOCKS | SCIP_CONSEXPR_PRINTDOT_EVALTAG | SCIP_CONSEXPR_PRINTDOT_ACTIVITYTAG


typedef unsigned int                      SCIP_CONSEXPR_PRINTDOT_WHAT; /**< type for printdot bitflags */
typedef struct SCIP_ConsExpr_PrintDotData SCIP_CONSEXPR_PRINTDOTDATA;  /**< printing a dot file data */

/** @} */

/** @name expression enforcement */
#define SCIP_CONSEXPR_EXPRENFO_NONE           0x0u /**< no enforcement */
#define SCIP_CONSEXPR_EXPRENFO_SEPABELOW      0x1u /**< separation for expr <= auxvar, thus might estimate expr from below */
#define SCIP_CONSEXPR_EXPRENFO_SEPAABOVE      0x2u /**< separation for expr >= auxvar, thus might estimate expr from above */
#define SCIP_CONSEXPR_EXPRENFO_SEPABOTH       (SCIP_CONSEXPR_EXPRENFO_SEPABELOW | SCIP_CONSEXPR_EXPRENFO_SEPAABOVE)  /**< separation for expr == auxvar */
#define SCIP_CONSEXPR_EXPRENFO_ACTIVITY       0x4u /**< activity computation (interval evaluation) and propagation (reverse propagation) */
#define SCIP_CONSEXPR_EXPRENFO_ALL            (SCIP_CONSEXPR_EXPRENFO_SEPABOTH | SCIP_CONSEXPR_EXPRENFO_ACTIVITY) /**< all enforcement methods */

typedef unsigned int                  SCIP_CONSEXPR_EXPRENFO_METHOD; /**< exprenfo bitflags */
typedef struct SCIP_ConsExpr_ExprEnfo SCIP_CONSEXPR_EXPRENFO;        /**< expression enforcement data */

/** @} */

/** @name Nonlinear Handler
 * @{
 */

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
 * - expr expression
 * - nlhdlrexprdata nonlinear handler expression data to be freed
 */
#define SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
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
 * in enforcing the relation between this expression (expr) and its auxiliary variable (auxvar) via
 * linear under- or overestimation, cut generation, and/or activity computation and propagation.
 *
 * We distinguish the following enforcement methods:
 * - SCIP_CONSEXPR_EXPRENFO_SEPABELOW: linear underestimation or cut generation for the relation expr <= auxvar (denoted as "below")
 * - SCIP_CONSEXPR_EXPRENFO_SEPAABOVE: linear overestimation or cut generation for the relation expr >= auxvar (denoted as "above")
 * - SCIP_CONSEXPR_EXPRENFO_ACTIVITY: domain propagation (i.e., constant under/overestimation) for the relation expr == auxvar.
 *
 * On input, parameter 'enforcing' indicates for any of these methods, whether
 * - it is not necessary to have such a method, e.g., because no auxvar will exist for expr, or no one uses or set activities of this expression,
 *   or due to analysis of the expression, expr >= auxvar is not necessary to be satisfied,
 * - or there already exists a nonlinear handler that will provide this method in an "enforcement" sense, that is,
 *   it believes that no one else could provide this method in a stronger sense. (This is mainly used by the default nlhdlr to check whether
 *   it should still reach out to the exprhdlr or whether is dominated by some nonlinear handler.)
 *
 * The DETECT callback shall augment the 'enforcing' bitmask by setting the enforcement methods it wants to provide in an "enforcement" sense.
 *
 * Additionally, the 'participating' bitmask shall be set if the nonlinear handler wants to be called on this expression at all.
 * Here, it shall set all methods that it wants to provide, which are those set in enforcing, but additionally those where it wants
 * to participate but leave enforcement to another nlhdlr.
 * This can be useful for nonlinear handlers that do not implement a complete enforcement, e.g., a handler that only contributes
 * cutting planes in some situations only.
 *
 * A nonlinear handler will be called only for those callbacks that it mentioned in participating, which is
 * - ENFO and/or ESTIMATE will be called with overestimate==FALSE if SCIP_CONSEXPR_EXPRENFO_SEPABELOW has been set
 * - ENFO and/or ESTIMATE will be called with overestimate==TRUE if SCIP_CONSEXPR_EXPRENFO_SEPAABOVE has been set
 * - INTEVAL and/or REVERSEPROP will be called if SCIP_CONSEXPR_EXPRENFO_ACTIVITY has been set
 * If SCIP_CONSEXPR_EXPRENFO_SEPABELOW or SCIP_CONSEXPR_EXPRENFO_SEPAABOVE has been set, then at least one of the
 * callbacks ENFO and ESTIMATE need to be implemented. Also EVALAUX will be called in this case.
 * If SCIP_CONSEXPR_EXPRENFO_ACTIVITY has been set, then at least one of INTEVAL and REVERSEPROP needs to be implemented.
 * If the nlhdlr chooses not to participate, then it must not return nlhdlrexprdata and can leave participating at its
 * initial value (SCIP_CONSEXPR_EXPRENFO_NONE).
 *
 * Additionally, a nonlinear handler that decides to participate in any of the enforcement methods must call
 * @ref SCIPregisterConsExprExprUsage() for every subexpression that it will use and indicate whether
 * - it will use an auxiliary variables,
 * - it will use activity for some subexpressions when computing estimators or cuts, and
 * - it will use activity for some subexpressions when for INTEVAL or REVERSEPROP.
 *
 * @note Auxiliary variables do not exist in subexpressions during detect and are not created by a call to @ref SCIPregisterConsExprExprUsage().
 *   They will be available when the INITSEPA callback is called.
 *
 * - scip SCIP data structure
 * - conshdlr expr-constraint handler
 * - nlhdlr nonlinear handler
 * - expr expression to analyze
 * - cons the constraint that expression defines, or NULL when the expr does not define any constraint, that is, when it is not the root of an expression of a constraint
 * - enforcing enforcement methods that are provided by some nonlinear handler (to be updated by detect callback)
 * - participating enforcement methods that this nonlinear handler should be called for (to be set by detect callback)
 * - nlhdlrexprdata nlhdlr's expr data to be stored in expr, can only be set to non-NULL if success is set to TRUE
 */
#define SCIP_DECL_CONSEXPR_NLHDLRDETECT(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONS* cons, \
   SCIP_CONSEXPR_EXPRENFO_METHOD* enforcing, \
   SCIP_CONSEXPR_EXPRENFO_METHOD* participating, \
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
 *  - intevalvar : callback to be called when interval evaluating a variable
 *  - intevalvardata : data to be passed to intevalvar callback
 */
#define SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_INTERVAL* interval, \
   SCIP_DECL_CONSEXPR_INTEVALVAR((*intevalvar)), \
   void* intevalvardata)

/** nonlinear handler callback for reverse propagation
 *
 * The method propagates the given bounds over the arguments of an expression.
 * The arguments of an expression are other expressions and the tighter intervals should be passed
 * to the corresponding argument (expression) via SCIPtightenConsExprExprInterval().
 *
 * input:
 *  - scip : SCIP main data structure
 *  - conshdlr: expr constraint handler
 *  - nlhdlr : nonlinear handler
 *  - expr : expression
 *  - nlhdlrexprdata : expression specific data of the nonlinear handler
 *  - bounds : the bounds on the expression that should be propagated
 *  - infeasible: buffer to store whether an expression's bounds were propagated to an empty interval
 *  - nreductions : buffer to store the number of interval reductions of all children
 */
#define SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_INTERVAL bounds, \
   SCIP_Bool* infeasible, \
   int* nreductions )

/** separation initialization method of a nonlinear handler (called during CONSINITLP)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : expression constraint handler
 *  - cons            : expression constraint
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
      SCIP_CONS* cons, \
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

/** nonlinear handler separation and enforcement callback
 *
 * The method tries to separate the given solution from the set defined by either
 *   expr - auxvar <= 0 (if !overestimate)
 * or
 *   expr - auxvar >= 0 (if  overestimate),
 * where auxvar = SCIPgetConsExprExprAuxVar(expr).
 *
 * It can do so by
 * - separation, i.e., finding an affine hyperplane (a cut) that separates
 *   the given point,
 * - bound tightening, i.e., changing bounds on a variable so that the given point is
 *   outside the updated domain,
 * - adding branching scores to potentially split the current problem into 2 subproblems
 * If parameter inenforcement is FALSE, then only the first option (separation) is allowed.
 *
 * If the NLHDLR always separates by computing a linear under- or overestimator of expr,
 * then it could be advantageous to implement the NLHDLRESTIMATE callback instead.
 *
 * Note, that the NLHDLR may also choose to separate for a relaxation of the mentioned sets,
 * e.g., expr <= upperbound(auxvar)  or  expr >= lowerbound(auxvar).
 * This is especially useful in situations where expr is the root expression of a constraint
 * and it is sufficient to satisfy lhs <= expr <= rhs. The cons_expr core ensures that
 * lhs <= lowerbound(auxvar) and upperbound(auxvar) <= rhs.
 *
 * cons_expr core may call this callback first with allowweakcuts=FALSE and repeat later with
 * allowweakcuts=TRUE, if it didn't succeed to enforce a solution without using weak cuts.
 * If in enforcement and the NLHDLR cannot enforce by separation or bound tightening, it should register
 * branching scores for those expressions where branching may help to compute tighter cuts in children.
 *
 * The nlhdlr must set result to SCIP_SEPARATED if it added a cut, to SCIP_REDUCEDDOM if it added
 * a bound change, and to SCIP_BRANCHED if it added branching scores.
 * Otherwise, it may set result to SCIP_DIDNOTRUN or SCIP_DIDNOTFIND.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - conshdlr : cons expr handler
 *  - cons : expression constraint
 *  - nlhdlr : nonlinear handler
 *  - expr : expression
 *  - nlhdlrexprdata : expression specific data of the nonlinear handler
 *  - sol : solution to be separated (NULL for the LP solution)
 *  - auxvalue : current value of expression w.r.t. auxiliary variables as obtained from EVALAUX
 *  - overestimate : whether the expression needs to be over- or underestimated
 *  - allowweakcuts : whether we should only look for "strong" cuts, or anything that separates is fine
 *  - separated : whether another nonlinear handler already added a cut for this expression
 *  - inenforcement: whether we are in enforcement, or only in separation
 *  - result : pointer to store the result
 */
#define SCIP_DECL_CONSEXPR_NLHDLRENFO(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONS* cons, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_SOL* sol, \
   SCIP_Real auxvalue, \
   SCIP_Bool overestimate, \
   SCIP_Bool allowweakcuts, \
   SCIP_Bool separated, \
   SCIP_Bool addbranchscores, \
   SCIP_RESULT* result)

/** nonlinear handler under/overestimation callback
 *
 * The method tries to compute a linear under- or overestimator that is as tight as possible
 * at a given point.
 * If the value of the estimator in the solution is smaller (larger) than targetvalue
 * when underestimating (overestimating), then no estimator needs to be computed.
 * Note, that targetvalue can be infinite if any estimator will be accepted.
 * If successful, it shall store the estimator in a given rowprep data structure and set the
 * rowprep->local flag accordingly.
 * It is assumed that the sidetype of the rowprep is not changed by the callback.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - conshdlr : constraint handler
 *  - nlhdlr : nonlinear handler
 *  - expr : expression
 *  - nlhdlrexprdata : expression data of nonlinear handler
 *  - sol  : solution at which to estimate (NULL for the LP solution)
 *  - auxvalue : current value of expression w.r.t. auxiliary variables as obtained from EVALAUX
 *  - overestimate : whether the expression needs to be over- or underestimated
 *  - targetvalue : a value the estimator shall exceed, can be +/-infinity
 *  - rowprep : a rowprep where to store the estimator
 *  - rowpreps: an array where to store the estimators
 *  - success : buffer to indicate whether an estimator could be computed
 *  - addbranchscores: indicates whether to register branching scores
 *  - addedbranchscores: buffer to store whether the branching score callback was successful
 */
#define SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(x) SCIP_RETCODE x (\
   SCIP* scip, \
   SCIP_CONSHDLR* conshdlr, \
   SCIP_CONSEXPR_NLHDLR* nlhdlr, \
   SCIP_CONSEXPR_EXPR* expr, \
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_SOL* sol, \
   SCIP_Real auxvalue, \
   SCIP_Bool overestimate, \
   SCIP_Real targetvalue, \
   SCIP_PTRARRAY* rowpreps, \
   SCIP_Bool* success, \
   SCIP_Bool addbranchscores, \
   SCIP_Bool* addedbranchscores)

typedef struct SCIP_ConsExpr_Nlhdlr         SCIP_CONSEXPR_NLHDLR;          /**< nonlinear handler */
typedef struct SCIP_ConsExpr_NlhdlrData     SCIP_CONSEXPR_NLHDLRDATA;      /**< nonlinear handler data */
typedef struct SCIP_ConsExpr_NlhdlrExprData SCIP_CONSEXPR_NLHDLREXPRDATA;  /**< nonlinear handler data for a specific expression */

/** @} */

/** evaluation callback for (vertex-polyhedral) functions used as input for facet computation of its envelopes
 *
 * input:
 * - args the point to be evaluated
 * - nargs the number of arguments of the function (length of array args)
 * - funcdata user-data of function evaluation callback
 *
 * return:
 * - value of function in point x or SCIP_INVALID if could not be evaluated
 */
#define SCIP_DECL_VERTEXPOLYFUN(f) SCIP_Real f (SCIP_Real* args, int nargs, void* funcdata)

/** maximum dimension of vertex-polyhedral function for which we can try to compute a facet of its convex or concave envelope */
#define SCIP_MAXVERTEXPOLYDIM 14


typedef struct SCIP_ConsExpr_QuadExpr      SCIP_CONSEXPR_QUADEXPR;      /**< representation of expression as quadratic */


/** storage for a linear row in preparation
 *
 * Uses to assemble data that could eventually make a SCIP_ROW.
 * @note Only one-sided rows are allowed here.
 */
struct SCIP_RowPrep
{
   SCIP_VAR**            vars;               /**< variables */
   SCIP_Real*            coefs;              /**< coefficients of variables */
   int                   nvars;              /**< number of variables (= number of coefficients) */
   int                   varssize;           /**< length of variables array (= lengths of coefficients array) */
   SCIP_Real             side;               /**< side */
   SCIP_SIDETYPE         sidetype;           /**< type of side */
   SCIP_Bool             local;              /**< whether the row is only locally valid (i.e., for the current node) */
   char                  name[SCIP_MAXSTRLEN]; /**< row name */

   SCIP_Bool             recordmodifications;/**< whether to remember variables which coefficients were modified during cleanup */
   SCIP_VAR**            modifiedvars;       /**< variables which coefficient were modified by cleanup */
   int                   nmodifiedvars;      /**< number of variables which coefficient was modified */
   int                   modifiedvarssize;   /**< length of modifiedvars array */
   SCIP_Bool             modifiedside;       /**< whether the side was modified (relaxed) by cleanup */
};
typedef struct SCIP_RowPrep SCIP_ROWPREP;

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_CONS_EXPR_H__ */
