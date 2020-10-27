/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_expr.h
 * @brief  public functions to work with algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_PUB_EXPR_H_
#define SCIP_PUB_EXPR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "scip/def.h"
#include "scip/type_expr.h"

// TODO everything with SCIP* move to scip_expr.h

/**@name Expression Handler Methods */
/**@{ */

/** set the expression handler callbacks to copy and free an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetCopyFreeHdlr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOPYHDLR((*copyhdlr)),      /**< handler copy callback (can be NULL) */
   SCIP_DECL_EXPRFREEHDLR((*freehdlr))       /**< handler free callback (can be NULL) */
);

/** set the expression handler callbacks to copy and free expression data */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetCopyFreeData(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOPYDATA((*copydata)),      /**< expression data copy callback (can be NULL for expressions without data) */
   SCIP_DECL_EXPRFREEDATA((*freedata))       /**< expression data free callback (can be NULL if data does not need to be freed) */
);

/** set the print callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetPrint(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRPRINT((*print))             /**< print callback (can be NULL) */
);

/** set the parse callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetParse(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRPARSE((*parse))             /**< parse callback (can be NULL) */
);

/** set the curvature detection callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetCurvature(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCURVATURE((*curvature))     /**< curvature detection callback (can be NULL) */
);

/** set the monotonicity detection callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetMonotonicity(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRMONOTONICITY((*monotonicity)) /**< monotonicity detection callback (can be NULL) */
);

/** set the integrality detection callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetIntegrality(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINTEGRALITY((*integrality)) /**< integrality detection callback (can be NULL) */
);

/** set the hash callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetHash(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRHASH((*hash))               /**< hash callback (can be NULL) */
);

/** set the compare callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetCompare(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOMPARE((*compare))         /**< compare callback (can be NULL) */
);

/** set derivative evaluation callbacks of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetDiff(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRBWDIFF((*bwdiff)),          /**< backward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRFWDIFF((*fwdiff)),          /**< forward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRBWFWDIFF((*bwfwdiff))       /**< backward-forward derivative evaluation callback (can be NULL) */
);

/** set the interval evaluation callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetIntEval(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINTEVAL((*inteval))         /**< interval evaluation callback (can be NULL) */
);

/** set the simplify callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetSimplify(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRSIMPLIFY((*simplify))       /**< simplify callback (can be NULL) */
);

/** set the reverse propagation callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetReverseProp(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
);

/** set the separation and estimation callbacks of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrSetSepa(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINITSEPA((*initsepa)),      /**< separation initialization callback (can be NULL) */
   SCIP_DECL_EXPREXITSEPA((*exitsepa)),      /**< separation deinitialization callback (can be NULL) */
   SCIP_DECL_EXPRESTIMATE((*estimate))       /**< estimator callback (can be NULL) */
);

/** gives the name of an expression handler */
SCIP_EXPORT
const char* SCIPexprhdlrGetName(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
);

/** gives the description of an expression handler (can be NULL) */
SCIP_EXPORT
const char* SCIPexprhdlrGetDescription(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
);

/** gives the precedence of an expression handler */
SCIP_EXPORT
unsigned int SCIPexprhdlrGetPrecedence(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
);

/** gives the data of an expression handler */
SCIP_EXPORT
SCIP_EXPRHDLRDATA* SCIPexprhdlrGetData(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
);

/** returns whether expression handler implements the print callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasPrint(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the backward differentiation callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasBwdiff(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the interval evaluation callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasIntEval(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the estimator callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasEstimate(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the simplification callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasSimplify(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the curvature callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasCurvature(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the reverse propagation callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasReverseProp(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the initialization callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasInitSepa(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the deinitialization callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasExitSepa(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** calls the print callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRPRINT(SCIPexprhdlrPrint);

/** calls the parse callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRPARSE(SCIPexprhdlrParse);

/** calls the expression hash callback */
SCIP_EXPORT
SCIP_DECL_EXPRHASH(SCIPexprhdlrHash);

/** calls the expression compare callback */
SCIP_EXPORT
SCIP_DECL_EXPRCOMPARE(SCIPexprhdlrCompare);

/** calls the backward-differentiation callback of an expression handler
 *
 * further, allows to different w.r.t. given expression and children values
 */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrBwdiff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   int                   childidx,           /**< index of child w.r.t. which to compute derivative */
   SCIP_Real*            derivative,         /**< buffer to store value of derivative */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_Real             exprval             /**< value for expression, used only if childrenvals is not NULL */
);

/** calls the backward-forward differentiation callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRBWFWDIFF(SCIPexprhdlrBwfwdiff);

/** calls the forward differentiation callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRFWDIFF(SCIPexprhdlrFwdiff);

/** calls the evaluation callback of an expression handler
 *
 * further, allows to evaluate w.r.t. given children values
 */
SCIP_EXPORT
SCIP_RETCODE SCIPexprhdlrEval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            val,                /**< buffer store value of expression */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_SOL*             sol                 /**< solution that is evaluated (used by the var-expression) */
);

/** calls the expression interval evaluation callback */
SCIP_EXPORT
SCIP_DECL_EXPRINTEVAL(SCIPexprhdlrIntEval);

/** calls estimator method of expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRESTIMATE(SCIPexprhdlrEstimate);

/** calls the simplification method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRSIMPLIFY(SCIPexprhdlrSimplify);

/** calls the curvature check method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRCURVATURE(SCIPexprhdlrCurvature);

/** calls the expression callback for reverse propagation */
SCIP_EXPORT
SCIP_DECL_EXPRREVERSEPROP(SCIPexprhdlrReverseProp);

/** calls the separation initialization method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRINITSEPA(SCIPexprhdlrInitSepa);

/** calls the separation deinitialization method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPREXITSEPA(SCIPexprhdlrExitSepa);

/** increments the branching score count of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrIncrementNBranchScore(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** @} */



/**@name Expression Methods */
/**@{ */

/** gets the number of times the expression is currently captured */
SCIP_EXPORT
int SCIPexprGetNUses(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the number of children of an expression */
SCIP_EXPORT
int SCIPexprGetNChildren(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the children of an expression (can be NULL if no children) */
SCIP_EXPORT
SCIP_EXPR** SCIPexprGetChildren(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gets the expression handler of an expression
 *
 * This identifies the type of the expression (sum, variable, ...).
 */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPexprGetHdlr(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gets the expression data of an expression */
SCIP_EXPORT
SCIP_EXPRDATA* SCIPexprGetData(
   SCIP_EXPR*            expr                /**< expression */
   );

/** sets the expression data of an expression
 *
 * The pointer to possible old data is overwritten and the
 * freedata-callback is not called before.
 * This function is intended to be used by expression handler.
 */
SCIP_EXPORT
void SCIPexprSetData(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRDATA*        exprdata            /**< expression data to be set (can be NULL) */
   );

/** gives the value from the last evaluation of an expression (or SCIP_INVALID if there was an eval error) */
SCIP_EXPORT
SCIP_Real SCIPexprGetEvalValue(
   SCIP_EXPR*            expr                /**< expression */
   );

/* TODO make private or remove */
/** sets the evaluation value */
SCIP_EXPORT
void SCIPexprSetEvalValue(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             value,              /**< value to set */
   unsigned int          tag                 /**< tag of solution that was evaluated, or 0 */
   );

/** gives the evaluation tag from the last evaluation, or 0 */
SCIP_EXPORT
unsigned int SCIPexprGetEvalTag(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the derivative stored in an expression (or SCIP_INVALID if there was an evaluation error) */
SCIP_EXPORT
SCIP_Real SCIPexprGetDerivative(
   SCIP_EXPR*     expr              /**< expression */
   );

/** gives the value of directional derivative from the last evaluation of a directional derivative of expression (or SCIP_INVALID if there was an error) */
SCIP_EXPORT
SCIP_Real SCIPexprGetDot(
   SCIP_EXPR*     expr              /**< expression */
   );

/** returns the difftag stored in an expression
 *
 * can be used to check whether partial derivative value is valid
 */
SCIP_EXPORT
unsigned int SCIPexprGetDiffTag(
   SCIP_EXPR*     expr              /**< expression */
   );

/** returns the activity of the expression
 *
 * TODO this is conshdlr specific:
 * The caller needs to make sure that the activity is valid.
 * For expression and nonlinear handlers, this is made sure when the following callbacks are called:
 * - interval evaluation (intervals for children only)
 * - reverse propagation
 * - estimate and enforce (for exprs where activity usage was signaled during nlhdlr detect)
 */
SCIP_EXPORT
SCIP_INTERVAL SCIPexprGetActivity(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the tag associated with the activity of the expression
 *
 * TODO this is conshdlr specific:
 * Can be compared with SCIPgetConsExprCurBoundsTag() and SCIPgetConsExprLastBoundRelaxTag()
 * to check whether the activity currently stored in this expression is current and valid, respectively.
 */
SCIP_EXPORT
unsigned int SCIPexprGetActivityTag(
   SCIP_EXPR*            expr                /**< expression */
   );

/** compare expressions
 * @return -1, 0 or 1 if expr1 <, =, > expr2, respectively
 * @note: The given expressions are assumed to be simplified.
 */
SCIP_EXPORT
int SCIPexprCompare(
   SCIP_EXPR*            expr1,              /**< first expression */
   SCIP_EXPR*            expr2               /**< second expression */
   );

/** returns the curvature of an expression
 *
 *  @note Call SCIPcomputeExprCurvature before calling this function.
 */
SCIP_EXPORT
SCIP_EXPRCURV SCIPexprGetCurvature(
   SCIP_EXPR*            expr                /**< expression */
   );

/** sets the curvature of an expression */
SCIP_EXPORT
void SCIPexprSetCurvature(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRCURV         curvature           /**< curvature of the expression */
   );

/** returns whether an expression is integral */
SCIP_EXPORT
SCIP_Bool SCIPexprIsIntegral(
   SCIP_EXPR*            expr                /**< expression */
   );

/** @} */

/**@name Quadratic expression functions */
/**@{ */

/** gives the coefficients and expressions that define a quadratic expression
 *
 * It can return the constant part, the number, arguments, and coefficients of the purely linear part
 * and the number of quadratic terms and bilinear terms.
 * Note that for arguments that appear in the quadratic part, a linear coefficient is
 * stored with the quadratic term.
 * Use SCIPexprGetQuadraticQuadTerm() and SCIPexprGetQuadraticBilinTerm()
 * to access the data for a quadratic or bilinear term.
 *
 * It can also return the eigenvalues and the eigenvectors of the matrix Q when the quadratic is written
 * as x^T Q x + b^T x + c^T y + d, where c^T y defines the purely linear part.
 * Note, however, that to have access to them one needs to call SCIPcomputeExprQuadraticCurvature()
 * with storeeigeninfo equals to TRUE. If the eigen infortmation was not stored or it failed to be computed,
 * eigenvalues and eigenvectors will be set to NULL.
 *
 *
 * This function returns pointers to internal data in linexprs and lincoefs.
 * The user must not change this data.
 */
SCIP_EXPORT
void SCIPexprGetQuadraticData(
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_Real*            constant,           /**< buffer to store constant term, or NULL */
   int*                  nlinexprs,          /**< buffer to store number of expressions that appear linearly, or NULL */
   SCIP_EXPR***          linexprs,           /**< buffer to store pointer to array of expressions that appear linearly, or NULL */
   SCIP_Real**           lincoefs,           /**< buffer to store pointer to array of coefficients of expressions that appear linearly, or NULL */
   int*                  nquadexprs,         /**< buffer to store number of expressions in quadratic terms, or NULL */
   int*                  nbilinexprs,        /**< buffer to store number of bilinear expressions terms, or NULL */
   SCIP_Real**           eigenvalues,        /**< buffer to store pointer to array of eigenvalues of Q, or NULL */
   SCIP_Real**           eigenvectors        /**< buffer to store pointer to array of eigenvectors of Q, or NULL */
   );

/** gives the data of a quadratic expression term
 *
 * For a term a*expr^2 + b*expr + sum_i (c_i * expr * otherexpr_i), returns
 * expr, a, b, the number of summands, and indices of bilinear terms in the quadratic expressions bilinexprterms.
 *
 * This function returns pointers to internal data in adjbilin.
 * The user must not change this data.
 */
SCIP_EXPORT
void SCIPexprGetQuadraticQuadTerm(
   SCIP_EXPR*            expr,               /**< quadratic expression */
   int                   termidx,            /**< index of quadratic term */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to argument expression (the 'x') of this term, or NULL */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of variable, or NULL */
   SCIP_Real*            sqrcoef,            /**< buffer to store square coefficient of variable, or NULL */
   int*                  nadjbilin,          /**< buffer to store number of bilinear terms this variable is involved in, or NULL */
   int**                 adjbilin,           /**< buffer to store pointer to indices of associated bilinear terms, or NULL */
   SCIP_EXPR**           sqrexpr             /**< buffer to store pointer to square expression (the 'x^2') of this term or NULL if no square expression, or NULL */
   );

/** gives the data of a bilinear expression term
 *
 * For a term a*expr1*expr2, returns
 * expr1, expr2, a, and the position of the quadratic expression term that uses expr2 in the quadratic expressions quadexprterms.
 */
SCIP_EXPORT
void SCIPexprGetQuadraticBilinTerm(
   SCIP_EXPR*            expr,               /**< quadratic expression */
   int                   termidx,            /**< index of bilinear term */
   SCIP_EXPR**           expr1,              /**< buffer to store first factor, or NULL */
   SCIP_EXPR**           expr2,              /**< buffer to store second factor, or NULL */
   SCIP_Real*            coef,               /**< buffer to coefficient, or NULL */
   int*                  pos2,               /**< buffer to position of expr2 in quadexprterms array of quadratic expression, or NULL */
   SCIP_EXPR**           prodexpr            /**< buffer to store pointer to expression that is product if first and second factor, or NULL */
   );

/** returns whether all expressions that are used in a quadratic expression are variable expression
 *
 * @return TRUE iff all linexprs and quadexprterms[.].expr are variable expressions
 */
SCIP_EXPORT
SCIP_Bool SCIPexprAreQuadraticExprsVariables(
   SCIP_EXPR*            expr,               /**< quadratic expression */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* SCIP_PUB_EXPR_H_ */
