/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_nlhdlr.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions related to nonlinear handlers of nonlinear constraints
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 *
 *  This file defines the interface for nonlinear handlers.
 *
 *  - \ref NLHDLRS "List of available nonlinear handlers"
 */

#ifndef SCIP_TYPE_NLHDLR_H_
#define SCIP_TYPE_NLHDLR_H_

#include "scip/type_expr.h"
#include "scip/type_cons.h"
#include "scip/type_misc.h"

#define SCIP_NLHDLR_METHOD_NONE       0x0u   /**< no enforcement */
#define SCIP_NLHDLR_METHOD_SEPABELOW  0x1u   /**< separation for expr <= auxvar, thus might estimate expr from below */
#define SCIP_NLHDLR_METHOD_SEPAABOVE  0x2u   /**< separation for expr >= auxvar, thus might estimate expr from above */
#define SCIP_NLHDLR_METHOD_SEPABOTH   (SCIP_NLHDLR_METHOD_SEPABELOW | SCIP_NLHDLR_METHOD_SEPAABOVE)  /**< separation for expr == auxvar */
#define SCIP_NLHDLR_METHOD_ACTIVITY   0x4u   /**< activity computation (interval evaluation) and propagation (reverse propagation) */
#define SCIP_NLHDLR_METHOD_ALL        (SCIP_NLHDLR_METHOD_SEPABOTH | SCIP_NLHDLR_METHOD_ACTIVITY) /**< all enforcement methods */

typedef unsigned int SCIP_NLHDLR_METHOD; /**< nlhdlr methods bitflags */

/** nonlinear handler copy callback
 *
 * The method includes the nonlinear handler into a nonlinear constraint handler.
 *
 * This method is usually called when doing a copy of a nonlinear constraint handler.
 *
 * - targetscip      : target SCIP main data structure
 * - targetconshdlr  : target nonlinear constraint handler
 * - sourceconshdlr  : nonlinear constraint handler in source SCIP
 * - sourcenlhdlr    : nonlinear handler in source SCIP
 */
#define SCIP_DECL_NLHDLRCOPYHDLR(x) SCIP_RETCODE x (\
   SCIP*          targetscip,     \
   SCIP_CONSHDLR* targetconshdlr, \
   SCIP_CONSHDLR* sourceconshdlr, \
   SCIP_NLHDLR*   sourcenlhdlr)

/** callback to free data of handler
 *
 * - scip       : SCIP data structure
 * - nlhdlr     : nonlinear handler
 * - nlhdlrdata : nonlinear handler data to be freed
 */
#define SCIP_DECL_NLHDLRFREEHDLRDATA(x) SCIP_RETCODE x (\
   SCIP*             scip,   \
   SCIP_NLHDLR*      nlhdlr, \
   SCIP_NLHDLRDATA** nlhdlrdata)

/** callback to free expression specific data
 *
 * - scip           : SCIP data structure
 * - nlhdlr         : nonlinear handler
 * - expr           : expression
 * - nlhdlrexprdata : nonlinear handler expression data to be freed
 */
#define SCIP_DECL_NLHDLRFREEEXPRDATA(x) SCIP_RETCODE x (\
   SCIP*                 scip,   \
   SCIP_NLHDLR*          nlhdlr, \
   SCIP_EXPR*            expr,   \
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata)

/** callback to be called in initialization
 *
 * - scip   : SCIP data structure
 * - nlhdlr : nonlinear handler
 */
#define SCIP_DECL_NLHDLRINIT(x) SCIP_RETCODE x (\
   SCIP*        scip, \
   SCIP_NLHDLR* nlhdlr)

/** callback to be called in deinitialization
 *
 * - scip   : SCIP data structure
 * - nlhdlr : nonlinear handler
 */
#define SCIP_DECL_NLHDLREXIT(x) SCIP_RETCODE x (\
   SCIP*        scip, \
   SCIP_NLHDLR* nlhdlr)

/** callback to detect structure in expression tree
 *
 * The nonlinear handler shall analyze the current expression and decide whether it wants to contribute
 * in enforcing the relation between this expression (expr) and its auxiliary variable (auxvar) via
 * linear under- or overestimation, cut generation, and/or activity computation and propagation.
 *
 * We distinguish the following enforcement methods:
 * - SCIP_NLHDLR_METHOD_SEPABELOW: linear underestimation or cut generation for the relation expr <= auxvar (denoted as "below")
 * - SCIP_NLHDLR_METHOD_SEPAABOVE: linear overestimation or cut generation for the relation expr >= auxvar (denoted as "above")
 * - SCIP_NLHDLR_METHOD_ACTIVITY: domain propagation (i.e., constant under/overestimation) for the relation expr == auxvar.
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
 * - ENFO and/or ESTIMATE will be called with overestimate==FALSE if SCIP_NLHDLR_METHOD_SEPABELOW has been set
 * - ENFO and/or ESTIMATE will be called with overestimate==TRUE if SCIP_NLHDLR_METHOD_SEPAABOVE has been set
 * - INTEVAL and/or REVERSEPROP will be called if SCIP_NLHDLR_METHOD_ACTIVITY has been set
 * If SCIP_NLHDLR_METHOD_SEPABELOW or SCIP_NLHDLR_METHOD_SEPAABOVE has been set, then at least one of the
 * callbacks ENFO and ESTIMATE need to be implemented. Also EVALAUX will be called in this case.
 * If SCIP_NLHDLR_METHOD_ACTIVITY has been set, then at least one of INTEVAL and REVERSEPROP needs to be implemented.
 * If the nlhdlr chooses not to participate, then it must not return nlhdlrexprdata and can leave participating at its
 * initial value (SCIP_NLHDLR_METHOD_NONE).
 *
 * Additionally, a nonlinear handler that decides to participate in any of the enforcement methods must call
 * @ref SCIPregisterExprUsageNonlinear() for every subexpression that it will use and indicate whether
 * - it will use an auxiliary variable,
 * - it will use activity for some subexpressions when computing estimators or cuts, and
 * - it will use activity for some subexpressions when for INTEVAL or REVERSEPROP.
 *
 * @note Auxiliary variables do not exist in subexpressions during detect and are not created by a call to @ref SCIPregisterExprUsageNonlinear().
 *   They will be available when the INITSEPA callback is called.
 *
 * - scip           : SCIP data structure
 * - conshdlr       : nonlinear constraint handler
 * - nlhdlr         : nonlinear handler
 * - expr           : expression to analyze
 * - cons           : the constraint that expression defines, or NULL when the expr does not define any constraint, that is, when it is not the root of an expression of a constraint
 * - enforcing      : enforcement methods that are provided by some nonlinear handler (to be updated by detect callback)
 * - participating  : enforcement methods that this nonlinear handler should be called for (to be set by detect callback)
 * - nlhdlrexprdata : nlhdlr's expr data to be stored in expr, can only be set to non-NULL if success is set to TRUE
 */
#define SCIP_DECL_NLHDLRDETECT(x) SCIP_RETCODE x (\
   SCIP*                 scip,          \
   SCIP_CONSHDLR*        conshdlr,      \
   SCIP_NLHDLR*          nlhdlr,        \
   SCIP_EXPR*            expr,          \
   SCIP_CONS*            cons,          \
   SCIP_NLHDLR_METHOD*   enforcing,     \
   SCIP_NLHDLR_METHOD*   participating, \
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata)

/** auxiliary evaluation callback of nonlinear handler
 *
 * Evaluates the expression w.r.t. the auxiliary variables that were introduced by the nonlinear handler (if any).
 * The method is used to determine the violation of the relation that the nonlinear
 * handler attempts to enforce. During enforcement, this violation value is used to
 * decide whether separation or branching score callbacks should be called.
 *
 * It can be assumed that the expression itself has been evaluated in the given sol.
 *
 * - scip           : SCIP data structure
 * - nlhdlr         : nonlinear handler
 * - expr           : expression to evaluate
 * - nlhdlrexprdata : expression specific data of the nonlinear handler
 * - auxvalue       : buffer to store value of expression w.r.t. auxiliary variables
 * - sol            : point to evaluate
 */
#define SCIP_DECL_NLHDLREVALAUX(x) SCIP_RETCODE x (\
   SCIP*                scip,           \
   SCIP_NLHDLR*         nlhdlr,         \
   SCIP_EXPR*           expr,           \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_Real*           auxvalue,       \
   SCIP_SOL*            sol)

/** nonlinear handler interval evaluation callback
 *
 * The method computes an interval that contains the image (range) of the expression.
 *
 * - scip           : SCIP main data structure
 * - nlhdlr         : nonlinear handler
 * - expr           : expression
 * - nlhdlrexprdata : expression specific data of the nonlinear handler
 * - interval       : buffer where to store interval (on input: current interval for expr, on output: computed interval for expr)
 * - intevalvar     : callback to be called when interval evaluating a variable
 * - intevalvardata : data to be passed to intevalvar callback
 */
#define SCIP_DECL_NLHDLRINTEVAL(x) SCIP_RETCODE x (\
   SCIP*                scip,                \
   SCIP_NLHDLR*         nlhdlr,              \
   SCIP_EXPR*           expr,                \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata,      \
   SCIP_INTERVAL*       interval,            \
   SCIP_DECL_EXPR_INTEVALVAR((*intevalvar)), \
   void*                intevalvardata)

/** nonlinear handler callback for reverse propagation
 *
 * The method propagates the given bounds over the arguments of an expression.
 * The arguments of an expression are other expressions and the tighter intervals should be passed
 * to the corresponding argument (expression) via SCIPtightenExprIntervalNonlinear().
 *
 * - scip           : SCIP main data structure
 * - conshdlr       : nonlinear constraint handler
 * - nlhdlr         : nonlinear handler
 * - expr           : expression
 * - nlhdlrexprdata : expression specific data of the nonlinear handler
 * - bounds         : the bounds on the expression that should be propagated
 * - infeasible     : buffer to store whether an expression's bounds were propagated to an empty interval
 * - nreductions    : buffer to store the number of interval reductions of all children
 */
#define SCIP_DECL_NLHDLRREVERSEPROP(x) SCIP_RETCODE x (\
   SCIP*                scip,           \
   SCIP_CONSHDLR*       conshdlr,       \
   SCIP_NLHDLR*         nlhdlr,         \
   SCIP_EXPR*           expr,           \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_INTERVAL        bounds,         \
   SCIP_Bool*           infeasible,     \
   int*                 nreductions)

/** separation initialization method of a nonlinear handler (called during CONSINITLP)
 *
 * - scip           : SCIP main data structure
 * - conshdlr       : nonlinear constraint handler
 * - cons           : nonlinear constraint
 * - nlhdlr         : nonlinear handler
 * - nlhdlrexprdata : exprdata of nonlinear handler
 * - expr           : expression
 * - overestimate   : whether the expression needs to be overestimated
 * - underestimate  : whether the expression needs to be underestimated
 * - infeasible     : pointer to store whether an infeasibility was detected while building the LP
 */
#define SCIP_DECL_NLHDLRINITSEPA(x) SCIP_RETCODE x (\
   SCIP*                scip,           \
   SCIP_CONSHDLR*       conshdlr,       \
   SCIP_CONS*           cons,           \
   SCIP_NLHDLR*         nlhdlr,         \
   SCIP_EXPR*           expr,           \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_Bool            overestimate,   \
   SCIP_Bool            underestimate,  \
   SCIP_Bool*           infeasible)

/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL)
 *
 * - scip            : SCIP main data structure
 * - nlhdlr          : nonlinear handler
 * - nlhdlrexprdata  : exprdata of nonlinear handler
 * - expr            : expression
 */
#define SCIP_DECL_NLHDLREXITSEPA(x) SCIP_RETCODE x (\
   SCIP*                scip,   \
   SCIP_NLHDLR*         nlhdlr, \
   SCIP_EXPR*           expr,   \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata)

/** nonlinear handler separation and enforcement callback
 *
 * The method tries to separate the given solution from the set defined by either
 *   expr - auxvar <= 0 (if !overestimate)
 * or
 *   expr - auxvar >= 0 (if  overestimate),
 * where auxvar = SCIPgetExprAuxVarNonlinear(expr).
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
 * - scip           : SCIP main data structure
 * - conshdlr       : cons nonlinear handler
 * - cons           : nonlinear constraint
 * - nlhdlr         : nonlinear handler
 * - expr           : expression
 * - nlhdlrexprdata : expression specific data of the nonlinear handler
 * - sol            : solution to be separated (NULL for the LP solution)
 * - auxvalue       : current value of expression w.r.t. auxiliary variables as obtained from EVALAUX
 * - overestimate   : whether the expression needs to be over- or underestimated
 * - allowweakcuts  : whether we should only look for "strong" cuts, or anything that separates is fine
 * - separated      : whether another nonlinear handler already added a cut for this expression
 * - inenforcement  : whether we are in enforcement, or only in separation
 * - result         : pointer to store the result
 */
#define SCIP_DECL_NLHDLRENFO(x) SCIP_RETCODE x (\
   SCIP*                scip,            \
   SCIP_CONSHDLR*       conshdlr,        \
   SCIP_CONS*           cons,            \
   SCIP_NLHDLR*         nlhdlr,          \
   SCIP_EXPR*           expr,            \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata,  \
   SCIP_SOL*            sol,             \
   SCIP_Real            auxvalue,        \
   SCIP_Bool            overestimate,    \
   SCIP_Bool            allowweakcuts,   \
   SCIP_Bool            separated,       \
   SCIP_Bool            addbranchscores, \
   SCIP_RESULT*         result)

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
 *  - scip              : SCIP main data structure
 *  - conshdlr          : constraint handler
 *  - nlhdlr            : nonlinear handler
 *  - expr              : expression
 *  - nlhdlrexprdata    : expression data of nonlinear handler
 *  - sol               : solution at which to estimate (NULL for the LP solution)
 *  - auxvalue          : current value of expression w.r.t. auxiliary variables as obtained from EVALAUX
 *  - overestimate      : whether the expression needs to be over- or underestimated
 *  - targetvalue       : a value the estimator shall exceed, can be +/-infinity
 *  - rowprep           : a rowprep where to store the estimator
 *  - rowpreps          : an array where to store the estimators
 *  - success           : buffer to indicate whether an estimator could be computed
 *  - addbranchscores   : indicates whether to register branching scores
 *  - addedbranchscores : buffer to store whether the branching score callback was successful
 */
#define SCIP_DECL_NLHDLRESTIMATE(x) SCIP_RETCODE x (\
   SCIP*                scip,            \
   SCIP_CONSHDLR*       conshdlr,        \
   SCIP_NLHDLR*         nlhdlr,          \
   SCIP_EXPR*           expr,            \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata,  \
   SCIP_SOL*            sol,             \
   SCIP_Real            auxvalue,        \
   SCIP_Bool            overestimate,    \
   SCIP_Real            targetvalue,     \
   SCIP_PTRARRAY*       rowpreps,        \
   SCIP_Bool*           success,         \
   SCIP_Bool            addbranchscores, \
   SCIP_Bool*           addedbranchscores)

typedef struct SCIP_Nlhdlr         SCIP_NLHDLR;          /**< nonlinear handler */
typedef struct SCIP_NlhdlrData     SCIP_NLHDLRDATA;      /**< nonlinear handler data */
typedef struct SCIP_NlhdlrExprData SCIP_NLHDLREXPRDATA;  /**< nonlinear handler data for a specific expression */

#endif /* SCIP_TYPE_NLHDLR_H_ */
