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

/**@file   cons_expr.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for expression constraints (in particular, nonlinear constraints)
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_H__
#define __SCIP_CONS_EXPR_H__


#include "scip/scip.h"
#include "scip/type_cons_expr.h"


#ifdef __cplusplus
extern "C" {
#endif

/**@name Expression Handler Methods */
/**@{ */

/** creates the handler for an expression handler and includes it into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprExprHdlrBasic(
   SCIP*                       scip,         /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR**    exprhdlr,     /**< buffer where to store expression handler */
   const char*                 name,         /**< name of expression handler (must not be NULL) */
   const char*                 desc,         /**< description of expression handler (can be NULL) */
   unsigned int                precedence,   /**< precedence of expression operation (used for printing) */
   SCIP_DECL_CONSEXPR_EXPREVAL((*eval)),     /**< point evaluation callback (can not be NULL) */
   SCIP_CONSEXPR_EXPRHDLRDATA* data          /**< data of expression handler (can be NULL) */
   );

/** set the expression handler callbacks to copy and free an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeHdlr(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOPYHDLR((*copyhdlr)), /**< handler copy callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRFREEHDLR((*freehdlr))  /**< handler free callback (can be NULL) */
);

/** set the expression handler callbacks to copy and free expression data */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOPYDATA((*copydata)), /**< expression data copy callback (can be NULL for expressions without data) */
   SCIP_DECL_CONSEXPR_EXPRFREEDATA((*freedata))  /**< expression data free callback (can be NULL if data does not need to be freed) */
);

/** set the print callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrPrint(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRPRINT((*print))    /**< print callback (can be NULL) */
);

/** set the parse callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrParse(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRPARSE((*parse))    /**< parse callback (can be NULL) */
);

/** set the curvature detection callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrCurvature(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCURVATURE((*curvature)) /**< curvature detection callback (can be NULL) */
);

/** set the monotonicity detection callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrMonotonicity(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRMONOTONICITY((*monotonicity)) /**< monotonicity detection callback (can be NULL) */
);

/** set the integrality detection callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrIntegrality(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRINTEGRALITY((*integrality)) /**< integrality detection callback (can be NULL) */
);

/** set the hash callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrHash(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRHASH((*hash))      /**< hash callback (can be NULL) */
);

/** set the compare callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrCompare(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOMPARE((*compare))/**< compare callback (can be NULL) */
);

/** set the derivative evaluation callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrBwdiff(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRBWDIFF((*bwdiff))  /**< derivative evaluation callback (can be NULL) */
);

/** set the interval evaluation callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrIntEval(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRINTEVAL((*inteval))/**< interval evaluation callback (can be NULL) */
);

/** set the simplify callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrSimplify(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRSIMPLIFY((*simplify))  /**< simplify callback (can be NULL) */
);

/** set the reverse propagation callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrReverseProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRREVERSEPROP((*reverseprop))/**< reverse propagation callback (can be NULL) */
);

/** set the separation and estimation callbacks of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRINITSEPA((*initsepa)), /**< separation initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPREXITSEPA((*exitsepa)), /**< separation deinitialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRESTIMATE((*estimate))  /**< estimator callback (can be NULL) */
);

/** gives expression handlers */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR** SCIPgetConsExprExprHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
);

/** gives number of expression handlers */
SCIP_EXPORT
int SCIPgetConsExprExprNHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
);

/** returns an expression handler of a given name (or NULL if not found) */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPfindConsExprExprHdlr(
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   const char*                name           /**< name of expression handler */
   );

/** returns expression handler for variable expressions */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrVar(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for constant value expressions */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrValue(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for sum expressions */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrSum(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for product expressions */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrProduct(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for power expressions */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrPower(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for signed power expressions */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrSignPower(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for exponential expressions */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrExponential(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for logarithm expressions */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrLog(
        SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
);

/** gives the name of an expression handler */
SCIP_EXPORT
const char* SCIPgetConsExprExprHdlrName(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
);

/** gives the description of an expression handler (can be NULL) */
SCIP_EXPORT
const char* SCIPgetConsExprExprHdlrDescription(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
);

/** gives the precedence of an expression handler */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprHdlrPrecedence(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
);

/** gives the data of an expression handler */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLRDATA* SCIPgetConsExprExprHdlrData(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr      /**< expression handler */
);

/** returns whether expression handler implements the print callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrPrint(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the backward differentiation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrBwdiff(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the interval evaluation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrIntEval(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the estimator callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrEstimate(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the simplification callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrSimplify(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the curvature callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrCurvature(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the reverse propagation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrReverseProp(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the initialization callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrInitSepa(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the deinitialization callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrExitSepa(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the branching score callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrBranchingScore(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** calls the print callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRPRINT(SCIPprintConsExprExprHdlr);

/** calls the parse callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRPARSE(SCIPparseConsExprExprHdlr);

/** calls the expression hash callback */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRHASH(SCIPhashConsExprExprHdlr);

/** calls the expression compare callback */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRCOMPARE(SCIPcompareConsExprExprHdlr);

/** calls the backward-differentiation callback of an expression handler
 *
 * further, allows to different w.r.t. given expression and children values
 */
SCIP_EXPORT
SCIP_RETCODE SCIPbwdiffConsExprExprHdlr(
   SCIP*                      scip,         /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*        expr,         /**< expression */
   int                        childidx,     /**< index of child w.r.t. which to compute derivative */
   SCIP_Real*                 derivative,   /**< buffer to store value of derivative */
   SCIP_Real*                 childrenvals, /**< values for children, or NULL if values stored in children should be used */
   SCIP_Real                  exprval       /**< value for expression, used only if childrenvals is not NULL */
);

/** calls the evaluation callback of an expression handler
 *
 * further, allows to evaluate w.r.t. given children values
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalConsExprExprHdlr(
   SCIP*                      scip,         /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*        expr,         /**< expression */
   SCIP_Real*                 val,          /**< buffer store value of expression */
   SCIP_Real*                 childrenvals, /**< values for children, or NULL if values stored in children should be used */
   SCIP_SOL*                  sol           /**< solution that is evaluated (used by the var-expression) */
);

/** calls the expression interval evaluation callback */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRINTEVAL(SCIPintevalConsExprExprHdlr);

/** calls estimator method of expression handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRESTIMATE(SCIPestimateConsExprExprHdlr);

/** calls the simplification method of an expression handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(SCIPsimplifyConsExprExprHdlr);

/** calls the curvature check method of an expression handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRCURVATURE(SCIPcurvatureConsExprExprHdlr);

/** calls the expression callback for reverse propagation */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRREVERSEPROP(SCIPreversepropConsExprExprHdlr);

/** calls the separation initialization method of an expression handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRINITSEPA(SCIPinitsepaConsExprExprHdlr);

/** calls the separation deinitialization method of an expression handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPREXITSEPA(SCIPexitsepaConsExprExprHdlr);

/** increments the branching score count of an expression handler */
SCIP_EXPORT
void SCIPincrementConsExprExprHdlrNBranchScore(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr
   );

/** returns whether we are ok to branch on auxiliary variables
 *
 * Currently returns whether depth of node in B&B tree is at least value of constraints/expr/branching/aux parameter.
 */
SCIP_EXPORT
SCIP_Bool SCIPgetConsExprBranchAux(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
);

/** @} */

/**@name Expression Methods */
/**@{ */

/** creates and captures an expression with given expression data and children */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_CONSEXPR_EXPRDATA* exprdata,         /**< expression data (expression assumes ownership) */
   int                     nchildren,        /**< number of children */
   SCIP_CONSEXPR_EXPR**    children          /**< children (can be NULL if nchildren is 0) */
   );

/** creates and captures an expression with given expression data and up to two children */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExpr2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_CONSEXPR_EXPRDATA* exprdata,         /**< expression data */
   SCIP_CONSEXPR_EXPR*     child1,           /**< first child (can be NULL) */
   SCIP_CONSEXPR_EXPR*     child2            /**< second child (can be NULL) */
   );

/** creates and captures an expression from a node in an (old-style) expression graph */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExpr3(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_EXPRGRAPH*         exprgraph,        /**< expression graph */
   SCIP_EXPRGRAPHNODE*     node              /**< expression graph node */
   );

/** creates and captures an expression representing a quadratic function */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExprQuadratic(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   int                     nlinvars,         /**< number of linear terms */
   SCIP_VAR**              linvars,          /**< array with variables in linear part */
   SCIP_Real*              lincoefs,         /**< array with coefficients of variables in linear part */
   int                     nquadterms,       /**< number of quadratic terms */
   SCIP_VAR**              quadvars1,        /**< array with first variables in quadratic terms */
   SCIP_VAR**              quadvars2,        /**< array with second variables in quadratic terms */
   SCIP_Real*              quadcoefs         /**< array with coefficients of quadratic terms */
   );

/** creates and captures an expression representing a monomial */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExprMonomial(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   int                     nfactors,         /**< number of factors in monomial */
   SCIP_VAR**              vars,             /**< variables in the monomial */
   SCIP_Real*              exponents         /**< exponent in each factor, or NULL if all 1.0 */
   );

/** appends child to the children list of expr */
SCIP_EXPORT
SCIP_RETCODE SCIPappendConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_EXPR*   child               /**< expression to be appended */
   );

/** remove all children of expr
 *
 * only use if you really know what you are doing
 */
SCIP_EXPORT
SCIP_RETCODE SCIPremoveConsExprExprChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** overwrites/replaces a child of an expressions
 *
 * @note the old child is released and the newchild is captured, unless they are the same (=same pointer)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPreplaceConsExprExprChild(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression which is going to replace a child */
   int                     childidx,         /**< index of child being replaced */
   SCIP_CONSEXPR_EXPR*     newchild          /**< the new child */
   );

/** duplicates the given expression */
SCIP_EXPORT
SCIP_RETCODE SCIPduplicateConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< original expression */
   SCIP_CONSEXPR_EXPR**  copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_Bool             copychildren        /**< whether children (and all successors) should be copied, too */
   );

/** gets the number of times the expression is currently captured */
SCIP_EXPORT
int SCIPgetConsExprExprNUses(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   );

/** captures an expression (increments usage count) */
SCIP_EXPORT
void SCIPcaptureConsExprExpr(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression to be captured */
   );

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_EXPORT
SCIP_RETCODE SCIPreleaseConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to expression to be released */
   );

/** gives the number of children of an expression */
SCIP_EXPORT
int SCIPgetConsExprExprNChildren(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** gives the children of an expression (can be NULL if no children) */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR** SCIPgetConsExprExprChildren(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** gets the handler of an expression
 *
 * This identifies the type of the expression (sum, variable, ...).
 */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlr(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** gets the expression data of an expression */
SCIP_EXPORT
SCIP_CONSEXPR_EXPRDATA* SCIPgetConsExprExprData(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** returns whether an expression is a variable expression */
SCIP_EXPORT
SCIP_Bool SCIPisConsExprExprVar(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** returns whether an expression is a value expression */
SCIP_EXPORT
SCIP_Bool SCIPisConsExprExprValue(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** returns the variable used for linearizing a given expression (return value might be NULL)
 *
 * @note for variable expression it returns the corresponding variable
 */
SCIP_EXPORT
SCIP_VAR* SCIPgetConsExprExprAuxVar(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** sets the expression data of an expression
 *
 * The pointer to possible old data is overwritten and the
 * freedata-callback is not called before.
 * This function is intended to be used by expression handler.
 */
SCIP_EXPORT
void SCIPsetConsExprExprData(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_CONSEXPR_EXPRDATA* exprdata          /**< expression data to be set (can be NULL) */
   );

/** print an expression as info-message */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be printed */
   FILE*                   file              /**< file to print to, or NULL for stdout */
   );

/** initializes printing of expressions in dot format to a give FILE* pointer */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExprDotInit(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_PRINTDOTDATA** dotdata,     /**< buffer to store dot printing data */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_CONSEXPR_PRINTDOT_WHAT whattoprint   /**< info on what to print for each expression */
   );

/** initializes printing of expressions in dot format to a file with given filename */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExprDotInit2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_PRINTDOTDATA** dotdata,     /**< buffer to store dot printing data */
   const char*             filename,         /**< name of file to print to */
   SCIP_CONSEXPR_PRINTDOT_WHAT whattoprint   /**< info on what to print for each expression */
   );

SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExprDot(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_PRINTDOTDATA* dotdata,      /**< data as initialized by \ref SCIPprintConsExprExprDotInit() */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to be printed */
   );

/** finishes printing of expressions in dot format */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExprDotFinal(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_PRINTDOTDATA** dotdata      /**< buffer where dot printing data has been stored */
   );

/** shows a single expression by use of dot and gv
 *
 * This function is meant for debugging purposes.
 * It's signature is kept as simple as possible to make it
 * easily callable from gdb, for example.
 *
 * It prints the expression into a temporary file in dot format, then calls dot to create a postscript file, then calls ghostview (gv) to show the file.
 * SCIP will hold until ghostscript is closed.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPshowConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to be printed */
   );

/** prints structure of an expression a la Maple's dismantle */
SCIP_EXPORT
SCIP_RETCODE SCIPdismantleConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to dismantle */
   );

/** Creates an expression from a string.
 * We specify the grammar that defines the syntax of an expression. Loosely speaking, a Base will be any "block",
 * a Factor is a Base to a power, a Term is a product of Factors and an Expression is a sum of terms
 * The actual definition:
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * Term       -> Factor { ("*" | "/" ) Factor }
 * Factor     -> Base [ "^" "number" | "^(" "number" ")" ]
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 * where [a|b] means a or b or none, (a|b) means a or b, {a} means 0 or more a.
 *
 * Note that Op and OpExpression are undefined. Op corresponds to the name of an expression handler and
 * OpExpression to whatever string the expression handler accepts (through its parse method).
 *
 * See also @ref parseExpr.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPparseConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   const char*           exprstr,            /**< string with the expr to parse */
   const char**          finalpos,           /**< buffer to store the position of exprstr where we finished reading, or NULL if not of interest */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to store the expr parsed */
   );

/** evaluate an expression in a point
 *
 * Iterates over expressions to also evaluate children, if necessary.
 * Value can be received via SCIPgetConsExprExprEvalValue().
 * If an evaluation error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 *
 * If a nonzero \p soltag is passed, then only (sub)expressions are
 * reevaluated that have a different solution tag. If a soltag of 0
 * is passed, then subexpressions are always reevaluated.
 * The tag is stored together with the value and can be received via
 * SCIPgetConsExprExprEvalTag().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be evaluated */
   SCIP_SOL*               sol,              /**< solution to be evaluated */
   unsigned int            soltag            /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/** gives the value from the last evaluation of an expression (or SCIP_INVALID if there was an eval error) */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprValue(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   );

/** sets the evaluation value */
SCIP_EXPORT
void SCIPsetConsExprExprEvalValue(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_Real               value,            /**< value to set */
   unsigned int            tag               /**< tag of solution that was evaluated, or 0 */
   );

/** gives the evaluation tag from the last evaluation, or 0 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprEvalTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   );

/** @name Differentiation methods
 * Automatic differentiation Backward mode:
 * Given a function, say, f(s(x,y),t(x,y)) there is a common mnemonic technique to compute its partial derivatives,
 * using a tree diagram. Suppose we want to compute the partial derivative of f w.r.t x. Write the function as a tree:
 * f
 * |-----|
 * s     t
 * |--|  |--|
 * x  y  x  y
 * The weight of an edge between two nodes represents the partial derivative of the parent w.r.t the children, eg,
 * f
 * |   is d_s f [where d is actually \f$ \partial \f$]
 * s
 * The weight of a path is the product of the weight of the edges in the path.
 * The partial derivative of f w.r.t. x is then the sum of the weights of all paths connecting f with x:
 * df/dx = d_s f * d_x s + d_t f * d_x t
 *
 * We follow this method in order to compute the gradient of an expression (root) at a given point (point).
 * Note that an expression is a DAG representation of a function, but there is a 1-1 correspondence between paths
 * in the DAG and path in a tree diagram of a function.
 * Initially, we set root->derivative to 1.0.
 * Then, traversing the tree in Depth First (see SCIPexpriteratorInit), for every expr that *has* children,
 * we store in its i-th child
 * child[i]->derivative = the derivative of expr w.r.t that child evaluated at point * expr->derivative
 * Example:
 * f->derivative = 1.0
 * s->derivative = d_s f * f->derivative = d_s f
 * x->derivative = d_x s * s->derivative = d_x s * d_s f
 * However, when the child is a variable expressions, we actually need to initialize child->derivative to 0.0
 * and afterwards add, instead of overwrite the computed value.
 * The complete example would then be:
 * f->derivative = 1.0, x->derivative = 0.0, y->derivative = 0.0
 * s->derivative = d_s f * f->derivative = d_s f
 * x->derivative += d_x s * s->derivative = d_x s * d_s f
 * y->derivative += d_t s * s->derivative = d_t s * d_s f
 * t->derivative = d_t f * f->derivative = d_t f
 * x->derivative += d_x t * t->derivative = d_x t * d_t f
 * y->derivative += d_t t * t->derivative = d_t t * d_t f
 *
 * At the end we have: x->derivative == d_x s * d_s f + d_x t * d_t f, y->derivative == d_t s * d_s f + d_t t * d_t f
 *
 * @{
 */

/** computes the gradient for a given point
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * Value can be received via SCIPgetConsExprExprPartialDiff().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeConsExprExprGradient(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be evaluated */
   SCIP_SOL*               sol,              /**< solution to be evaluated (NULL for the current LP solution) */
   unsigned int            soltag            /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/** returns the partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error)
 *
 * @note expression must belong to a constraint
 */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprPartialDiff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   );

/** returns the derivative stored in an expression (or SCIP_INVALID if there was an evaluation error) */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprDerivative(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   );

/** returns the difftag stored in an expression
 *
 * can be used to check whether partial derivative value is valid
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprDiffTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   );

/**@} */  /* end of differentiation methods */

/** returns the activity of the expression
 *
 * The caller needs to make sure that the activity is valid.
 * For expression and nonlinear handlers, this is made sure when the following callbacks are called:
 * - interval evaluation (intervals for children only)
 * - reverse propagation
 * - estimate and enforce (for exprs where activity usage was signaled during nlhdlr detect)
 */
SCIP_EXPORT
SCIP_INTERVAL SCIPgetConsExprExprActivity(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   );

/** returns the tag associated with the activity of the expression
 *
 * Can be compared with SCIPgetConsExprCurBoundsTag() and SCIPgetConsExprLastBoundRelaxTag()
 * to check whether the activity currently stored in this expression is current and valid, respectively.
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprActivityTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   );

/** possibly reevaluates and then returns the activity of the expression
 *
 * Reevaluate activity if currently stored is not valid (some bound was relaxed since last evaluation).
 * If validsufficient is set to FALSE, then it will also reevaluate activity if a bound tightening was happening
 * since last evaluation.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalConsExprExprActivity(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler, or NULL */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_INTERVAL*          activity,         /**< interval where to store expression */
   SCIP_Bool               validsufficient   /**< whether any valid activity is sufficient */
   );

/** returns bounds on the expression
 *
 * This gives an intersection of bounds from
 * - activity calculation (\ref SCIPgetConsExprExprActivity), if valid,
 * - auxiliary variable, if present,
 * - stored by \ref SCIPtightenConsExprExprInterval during domain propagation
 *
 * @note The returned interval can be empty!
 */
SCIP_EXPORT
SCIP_INTERVAL SCIPgetConsExprExprBounds(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   );

/** informs the expression about new bounds that can be used for reverse-propagation and to tighten bounds of
 * corresponding (auxiliary) variable (if any)
 *
 * @attention this function should only be called during domain propagation in cons_expr
 */
SCIP_EXPORT
SCIP_RETCODE SCIPtightenConsExprExprInterval(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be tightened */
   SCIP_INTERVAL           newbounds,        /**< new bounds for the expression */
   SCIP_Bool*              cutoff,           /**< buffer to store whether a cutoff was detected */
   int*                    ntightenings      /**< buffer to add the total number of tightenings, or NULL */
   );

/** mark constraints that include this expression to be propagated again
 *
 * This can be used by, e.g., nlhdlrs, to trigger a new propagation of constraints without
 * a change of variable bounds, e.g., because new information on the expression is available
 * that could potentially lead to tighter expression activity values.
 *
 * Note, that this call marks also constraints for propagation which only share some variable
 * with this expression.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPmarkConsExprExprPropagate(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to propagate again */
   );

/** increments the curboundstag and resets lastboundrelax in constraint handler data
 *
 * @note This method is not intended for normal use.
 *   These tags are maintained by the event handler for variable bound change events.
 *   This method is used by some unittests.
 */
SCIP_EXPORT
void SCIPincrementConsExprCurBoundsTag(
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_Bool               boundrelax        /**< indicates whether a bound was relaxed, i.e., lastboundrelax should be set too */
   );

/** adds violation-branching score to an expression
 *
 * Adds a score to the expression-specific violation-branching score, thereby marking it as branching candidate.
 * The expression must either be a variable expression or have an aux-variable.
 * In the latter case, branching on auxiliary variables must have been enabled.
 * In case of doubt, use SCIPaddConsExprExprsViolScore(). Roughly, the difference between these functions is that the current
 * function adds the violscore to the expression directly, while SCIPaddConsExprExprsViolScore() will split the
 * violation score among all the given expressions according to constraints/expr/branching/violsplit. See
 * SCIPaddConsExprExprsViolScore() for more details.
 */
SCIP_EXPORT
void SCIPaddConsExprExprViolScore(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression where to add branching score */
   SCIP_Real               violscore         /**< violation score to add to expression */
   );

/** adds violation-branching score to a set of expressions, distributing the score among all the expressions.
 *
 * Each expression must either be a variable expression or have an aux-variable.
 * If branching on aux-variables is disabled, then the violation branching score will be distributed among all among the
 * variables present in exprs
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddConsExprExprsViolScore(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_CONSEXPR_EXPR**    exprs,            /**< expressions where to add branching score */
   int                     nexprs,           /**< number of expressions */
   SCIP_Real               violscore,        /**< violation score to add to expression */
   SCIP_SOL*               sol,              /**< current solution */
   SCIP_Bool*              success           /**< buffer to store whether at least one violscore was added */
   );

/** gives violation-branching score stored in expression, or 0.0 if no valid score has been stored */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprViolScore(
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   );

/** returns the hash value of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprHash(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   unsigned int*           hashval           /**< pointer to store the hash value */
   );

/** compare expressions
 * @return -1, 0 or 1 if expr1 <, =, > expr2, respectively
 * @note: The given expressions are assumed to be simplified.
 */
SCIP_EXPORT
int SCIPcompareConsExprExprs(
   SCIP_CONSEXPR_EXPR*   expr1,              /**< first expression */
   SCIP_CONSEXPR_EXPR*   expr2               /**< second expression */
   );

/** simplifies an expression
 *
 * The given expression will be released and overwritten with the simplified expression.
 * To keep the expression, duplicate it via SCIPduplicateConsExprExpr before calling this method.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsimplifyConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSEXPR_EXPR*     rootexpr,         /**< expression to be simplified */
   SCIP_CONSEXPR_EXPR**    simplified,       /**< buffer to store simplified expression */
   SCIP_Bool*              changed,          /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*              infeasible        /**< buffer to store whether infeasibility has been detected */
   );

/** sets the curvature of an expression */
SCIP_EXPORT
void SCIPsetConsExprExprCurvature(
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_EXPRCURV         curvature           /**< curvature of the expression */
   );

/** returns the curvature of an expression
 *
 *  @note Call SCIPcomputeConsExprExprCurvature before calling this function.
 */
SCIP_EXPORT
SCIP_EXPRCURV SCIPgetConsExprExprCurvature(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** computes the curvature of a given expression and all its subexpressions
 *
 *  @note this function also evaluates all subexpressions w.r.t. current variable bounds
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeConsExprExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** computes the monotonicity of an expression w.r.t. to a given child */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprMonotonicity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   int                   childidx,           /**< index of child */
   SCIP_MONOTONE*        monotonicity        /**< buffer to store monotonicity */
   );

/** returns the number of positive rounding locks of an expression */
SCIP_EXPORT
int SCIPgetConsExprExprNLocksPos(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** returns the number of negative rounding locks of an expression */
SCIP_EXPORT
int SCIPgetConsExprExprNLocksNeg(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** computes integrality information of a given expression and all its subexpressions; the integrality information can
 * be accessed via SCIPisConsExprExprIntegral()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeConsExprExprIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** returns whether an expression is integral */
SCIP_EXPORT
SCIP_Bool SCIPisConsExprExprIntegral(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** number of nonlinear handlers whose activity computation and propagation methods depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprNPropUsesActivity(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** number of nonlinear handlers whose separation methods (estimate or enforcement) depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprNSepaUsesActivity(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** number of nonlinear handlers whose separation methods (estimate or enforcement) use auxiliary variable of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprNAuxvarUses(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** method to be called by a nlhdlr during NLHDLRDETECT to notify expression that it will be used
 *
 * - if useauxvar is enabled, then ensures that an auxiliary variable will be created in INITLP
 * - if useactivityforprop or useactivityforsepa{below,above} is enabled, then ensured that activity will be updated for expr
 * - if useactivityforprop is enabled, then increments the count returned by \ref SCIPgetConsExprExprNPropUsesActivity
 * - if useactivityforsepa{below,above} is enabled, then increments the count returned by \ref SCIPgetConsExprExprNSepaUsesActivity
 *   and also increments this count for all variables in the expression.
 *
 * The distinction into useactivityforprop and useactivityforsepa{below,above} is to recognize variables which domain influences
 * under/overestimators. Domain propagation routines (like OBBT) may invest more work for these variables.
 * The distinction into useactivityforsepabelow and useactivityforsepaabove is to recognize whether a nlhdlr that called this method
 * will use activity of expr in enfomethod sepabelow or enfomethod sepaabove.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPregisterConsExprExprUsage(
   SCIP*                 scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,         /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,             /**< expression */
   SCIP_Bool             useauxvar,        /**< whether an auxiliary variable will be used for estimate or cut generation */
   SCIP_Bool             useactivityforprop, /**< whether activity of expr will be used by domain propagation or activity calculation (inteval) */
   SCIP_Bool             useactivityforsepabelow, /**< whether activity of expr will be used by underestimation */
   SCIP_Bool             useactivityforsepaabove  /**< whether activity of expr will be used by overestimation */
   );

/** returns the total number of variables in an expression
 *
 * The function counts variables in common sub-expressions only once.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprNVars(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   int*                    nvars             /**< buffer to store the total number of variables */
   );

/** returns all variable expressions contained in a given expression; the array to store all variable expressions needs
 * to be at least of size the number of unique variables in the expression which is given by SCIPgetConsExprExprNVars()
 * and can be bounded by SCIPgetNVars().
 *
 * @note function captures variable expressions
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_CONSEXPR_EXPR**    varexprs,         /**< array to store all variable expressions */
   int*                    nvarexprs         /**< buffer to store the total number of variable expressions */
   );

/** computes absolute violation for auxvar relation in an expression w.r.t. original variables
 *
 * Assume the expression is f(x), where x are original (i.e., not auxiliary) variables.
 * Assume that f(x) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(x) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(x) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(x).
 *
 * If necessary, f is evaluated in the given solution. If that fails (domain error),
 * then viol is set to SCIPinfinity and both violover and violunder are set to TRUE.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprAbsOrigViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(x) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(x) is violated, or NULL */
   );

/** computes absolute violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(w) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(w) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(w).
 *
 * If the given value of f(w) is SCIP_INVALID, then viol is set to SCIPinfinity and
 * both violover and violunder are set to TRUE.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprAbsAuxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   );

/** computes relative violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * Taking the absolute violation from SCIPgetConsExprExprAbsAuxViolation, this function returns
 * the absolute violation divided by max(1,|f(w)|).
 *
 * If the given value of f(w) is SCIP_INVALID, then viol is set to SCIPinfinity and
 * both violover and violunder are set to TRUE.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprRelAuxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   );

/** @} */



/**@name Expression Constraint Handler Methods */
/**@{ */

/** gets the index an expression iterator can use to store iterator specific data in an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPactivateConsExprExprHdlrIterator(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   int*                       iterindex       /**< buffer to store iteration index */
   );

/** returns the index that an expression iterator used to store iterator specific data in an expression */
SCIP_EXPORT
void SCIPdeactivateConsExprExprHdlrIterator(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   int                        iterindex       /**< iteration index that is not used anymore */
   );

/** get a new tag that can be used to mark an expression as visited */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprHdlrNewVisitedTag(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** gets tag indicating current local variable bounds */
SCIP_EXPORT
unsigned int SCIPgetConsExprCurBoundsTag(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** gets the curboundstag at the last time where variable bounds were relaxed */
SCIP_EXPORT
unsigned int SCIPgetConsExprLastBoundRelaxTag(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** returns the hashmap that is internally used to map variables to their corresponding variable expressions */
SCIP_EXPORT
SCIP_HASHMAP* SCIPgetConsExprVarHashmap(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** notifies conshdlr that a variable expression is to be freed
 *
 * the conshdlr will then update its var2expr hashmap
 *
 * @note To be called only by var-exprhdlr.
 * @note Temporary method that will be replaced by ownerdata-free
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnotifyConsExprExprVarFreed(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*        varexpr         /**< variable expression to be freed */
   );

/** collects all bilinear terms for a given set of constraints
 *
 * @note This method should only be used for unit tests that depend on SCIPgetConsExprBilinTerms(),
 *       SCIPgetConsExprBilinTerm() or SCIPgetConsExprBilinTermIdx().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcollectConsExprBilinTerms(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_CONS**                conss,          /**< expression constraints */
   int                        nconss          /**< total number of expression constraints */
   );

/** returns the total number of bilinear terms that are contained in all expression constraints
 *
 *  @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 */
SCIP_EXPORT
int SCIPgetConsExprNBilinTerms(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** returns all bilinear terms that are contained in all expression constraints
 *
 * @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @note The value of the auxiliary variable of a bilinear term might be NULL, which indicates that the term does not have an auxiliary variable.
 */
SCIP_EXPORT
SCIP_CONSEXPR_BILINTERM* SCIPgetConsExprBilinTerms(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** returns the index of the bilinear term representing the product of the two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns -1 if the variables do not appear bilinearly.
 */
SCIP_EXPORT
int SCIPgetConsExprBilinTermIdx(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_VAR*                  x,              /**< first variable */
   SCIP_VAR*                  y               /**< second variable */
   );

/** returns the bilinear term that representing the product of two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns NULL if the variables do not appear bilinearly.
 */
SCIP_EXPORT
SCIP_CONSEXPR_BILINTERM* SCIPgetConsExprBilinTerm(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_VAR*                  x,              /**< first variable */
   SCIP_VAR*                  y               /**< second variable */
   );

/** evaluates an auxiliary expression for a bilinear term */
SCIP_EXPORT
SCIP_Real SCIPevalConsExprBilinAuxExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< first variable of the bilinear term */
   SCIP_VAR*             y,                  /**< second variable of the bilinear term */
   SCIP_CONSEXPR_AUXEXPR* auxexpr,           /**< auxiliary expression */
   SCIP_SOL*             sol                 /**< solution at which to evaluate (can be NULL) */
   );

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertBilinearTermExisting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   int                   nlockspos,          /**< number of positive expression locks */
   int                   nlocksneg           /**< number of negative expression locks */
   );

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertBilinearTermImplicit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   SCIP_Real             coefx,              /**< coefficient of x in the auxiliary expression */
   SCIP_Real             coefy,              /**< coefficient of y in the auxiliary expression */
   SCIP_Real             coefaux,            /**< coefficient of auxvar in the auxiliary expression */
   SCIP_Real             cst,                /**< constant of the auxiliary expression */
   SCIP_Bool             overestimate        /**< whether the auxiliary expression overestimates the bilinear product */
   );

/** returns the number of enforcements for an expression */
SCIP_EXPORT
int SCIPgetConsExprExprNEnfos(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

/** returns the data for one of the enforcements of an expression */
SCIP_EXPORT
void SCIPgetConsExprExprEnfoData(
   SCIP_CONSEXPR_EXPR*   expr,                         /**< expression */
   int                   idx,                          /**< position of enforcement in enfos array */
   SCIP_CONSEXPR_NLHDLR** nlhdlr,                      /**< buffer to store nlhldr */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata,      /**< buffer to store nlhdlr data for expression, or NULL */
   SCIP_CONSEXPR_EXPRENFO_METHOD* nlhdlrparticipation, /**< buffer to store methods where nonlinear handler participates, or NULL */
   SCIP_Bool*            sepabelowusesactivity,        /**< buffer to store whether sepabelow uses activity of some expression, or NULL */
   SCIP_Bool*            sepaaboveusesactivity,        /**< buffer to store whether sepaabove uses activity of some expression, or NULL */
   SCIP_Real*            auxvalue                      /**< buffer to store current auxvalue, or NULL */
   );

/** sets the auxiliary value of expression for one of the enforcements of an expression */
SCIP_EXPORT
void SCIPsetConsExprExprEnfoAuxValue(
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   int                   idx,                /**< position of enforcement in enfos array */
   SCIP_Real             auxvalue            /**< the new value of auxval */
   );

/** upgrading method for expression constraints into more specific constraints
 *
 * the method might upgrade an expression constraint into a set of upgrade constraints
 * the caller provided an array upgdconss to store upgrade constraints
 * the length of upgdconss is given by upgdconsssize
 * if an upgrade is not possible, set *nupgdconss to zero
 * if more than upgdconsssize many constraints shall replace cons, the function
 * should return the required number as negated value in *nupgdconss
 * i.e., if cons should be replaced by 3 constraints, the function should set
 * *nupgdconss to -3 and return with SCIP_OKAY
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cons            : the nonlinear constraint to upgrade
 *  - nvarexprs       : total number of variable expressions in the expression constraint
 *  - nupgdconss      : pointer to store number of constraints that replace this constraint
 *  - upgdconss       : array to store constraints that replace this constraint
 *  - upgdconsssize   : length of the provided upgdconss array
 */
#define SCIP_DECL_EXPRCONSUPGD(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONS* cons, int nvarexprs, \
      int* nupgdconss, SCIP_CONS** upgdconss, int upgdconsssize)


/** creates the handler for expr constraints and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrExpr(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes an expression constraint upgrade method into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_EXPRCONSUPGD((*exprconsupgd)),  /**< method to call for upgrading expression constraint, or NULL */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method by active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   );

/** creates and captures a expr constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   );

/** creates and captures a expr constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** creates and captures a quadratic expression constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   );

/** @} */

/**@name Expression Constraint Methods */
/**@{ */

/** returns the expression of the given expression constraint */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR* SCIPgetExprConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the root curvature of the given expression constraint
 *
 * @note The curvature information are computed during CONSINITSOL.
 */
SCIP_EXPORT
SCIP_EXPRCURV SCIPgetCurvatureConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns representation of the expression of the given expression constraint as quadratic form, if possible
 *
 * Only sets *quaddata to non-NULL if the whole expression is quadratic (in the non-extended formulation) and non-linear.
 * That is, the expr in each SCIP_CONSEXPR_QUADEXPRTERM will be a variable expressions and
 * \ref SCIPgetConsExprExprVarVar() can be used to retrieve the variable.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetQuadExprConsExpr(
   SCIP*                    scip,               /**< SCIP data structure */
   SCIP_CONS*               cons,               /**< constraint data */
   SCIP_CONSEXPR_QUADEXPR** quaddata            /**< buffer to store pointer to quaddata, if quadratic; stores NULL, otherwise */
   );

/** gets the expr constraint as a nonlinear row representation. */
SCIP_EXPORT
SCIP_RETCODE SCIPgetNlRowConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   );

/** gets the left hand side of an expression constraint */
SCIP_EXPORT
SCIP_Real SCIPgetLhsConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the right hand side of an expression constraint */
SCIP_EXPORT
SCIP_Real SCIPgetRhsConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** adds coef * var to expression constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddLinearTermConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             coef,               /**< coefficient */
   SCIP_VAR*             var                 /**< variable */
   );

/** gets absolute violation of expression constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * If this value is at most SCIPfeastol(scip), the constraint would be considered feasible.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetAbsViolationConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   );

/** gets scaled violation of expression constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * The scaling that is applied to the absolute violation of the constraint
 * depends on the setting of parameter constraints/expr/violscale.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetRelViolationConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   );

/** gives the unique index of an expression constraint
 *
 * Each expression constraint gets an index assigned when it is created.
 * This index never changes and is unique among all expression constraints
 * within the same SCIP instance.
 * Thus, it can be used to sort a set of expression constraints.
 */
SCIP_EXPORT
int SCIPgetConsExprIndex(
   SCIP_CONS*            cons                /**< constraint data */
   );

/** compares two expression constraints by its index
 *
 * Usable as compare operator in array sort functions.
 */
SCIP_EXPORT
int SCIPcompareConsExprIndex(
   void*                 cons1,
   void*                 cons2
   );

/** returns an equivalent linear constraint if possible */
SCIP_EXPORT
SCIP_RETCODE SCIPgetLinearConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_CONS**           lincons             /**< buffer to store linear constraint data */
   );

/** returns a variable that appears linearly that may be decreased without making any other constraint infeasible */
SCIP_EXPORT
SCIP_RETCODE SCIPgetLinvarMayDecreaseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   );

/** returns a variable that appears linearly that may be increased without making any other constraint infeasible */
SCIP_EXPORT
SCIP_RETCODE SCIPgetLinvarMayIncreaseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   );

/** detects nonlinear handlers that can handle the expressions and creates needed auxiliary variables
 *
 *  @note this method is only used for testing purposes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdetectConsExprNlhdlrs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check for auxiliary variables */
   int                   nconss              /**< total number of constraints */
   );

/** add the cut and maybe report branchscores */
SCIP_EXPORT
SCIP_RETCODE SCIPprocessConsExprRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_NLHDLR* nlhdlr,             /**< nonlinear handler which provided the estimator */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_ROWPREP*         rowprep,            /**< cut to be added */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             auxvalue,           /**< current value of expression w.r.t. auxiliary variables as obtained from EVALAUX */
   SCIP_Bool             allowweakcuts,      /**< whether we should only look for "strong" cuts, or anything that separates is fine */
   SCIP_Bool             branchscoresuccess, /**< buffer to store whether the branching score callback of the estimator was successful */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement, or only in separation */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result */
   );

/** checks whether an expression is quadratic and returns the corresponding coefficients
 *
 * An expression is quadratic if it is either a square (of some expression), a product (of two expressions),
 * or a sum of terms where at least one is a square or a product.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_QUADEXPR** quaddata         /**< buffer to store pointer to quadratic representation of expression, if it is quadratic, otherwise stores NULL */
   );

/** @} */

/**@name Nonlinear Handler Methods */
/**@{ */

/** creates the nonlinearity handler and includes it into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprNlhdlrBasic(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_NLHDLR**      nlhdlr,             /**< buffer where to store nonlinear handler */
   const char*                 name,               /**< name of nonlinear handler (must not be NULL) */
   const char*                 desc,               /**< description of nonlinear handler (can be NULL) */
   int                         detectpriority,     /**< detection priority of nonlinear handler */
   int                         enfopriority,       /**< enforcement priority of nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRDETECT((*detect)),  /**< structure detection callback of nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLREVALAUX((*evalaux)), /**< auxiliary evaluation callback of nonlinear handler */
   SCIP_CONSEXPR_NLHDLRDATA*   data                /**< data of nonlinear handler (can be NULL) */
   );

/** set the copy handler callback of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrCopyHdlr(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR((*copy)) /**< copy callback (can be NULL) */
);

/** set the nonlinear handler callback to free the nonlinear handler data */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrFreeHdlrData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA((*freehdlrdata)) /**< handler free callback (can be NULL) */
);

/** set the expression handler callback to free expression specific data of nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrFreeExprData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA((*freeexprdata)) /**< nonlinear handler expression data free callback (can be NULL if data does not need to be freed) */
);

/** set the initialization and deinitialization callback of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrInitExit(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRINIT((*init)),   /**< initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLREXIT((*exit))    /**< deinitialization callback (can be NULL) */
);

/** set the propagation callbacks of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRINTEVAL((*inteval)), /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
);

/** set the enforcement callbacks of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRINITSEPA((*initsepa)), /**< separation initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRENFO((*enfo)),         /**< enforcement callback (can be NULL if estimate is not NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRESTIMATE((*estimate)), /**< estimation callback (can be NULL if sepa is not NULL) */
   SCIP_DECL_CONSEXPR_NLHDLREXITSEPA((*exitsepa))  /**< separation deinitialization callback (can be NULL) */
);

/** gives name of nonlinear handler */
SCIP_EXPORT
const char* SCIPgetConsExprNlhdlrName(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr         /**< nonlinear handler */
);

/** gives description of nonlinear handler, can be NULL */
SCIP_EXPORT
const char* SCIPgetConsExprNlhdlrDesc(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr         /**< nonlinear handler */
);

/** gives detection priority of nonlinear handler */
SCIP_EXPORT
int SCIPgetConsExprNlhdlrDetectPriority(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr         /**< nonlinear handler */
);

/** gives enforcement priority of nonlinear handler */
SCIP_EXPORT
int SCIPgetConsExprNlhdlrEnfoPriority(
        SCIP_CONSEXPR_NLHDLR*      nlhdlr    /**< nonlinear handler */
);

/** returns a nonlinear handler of a given name (or NULL if not found) */
SCIP_EXPORT
SCIP_CONSEXPR_NLHDLR* SCIPfindConsExprNlhdlr(
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   const char*                name           /**< name of nonlinear handler */
   );

/** gives handler data of nonlinear handler */
SCIP_EXPORT
SCIP_CONSEXPR_NLHDLRDATA* SCIPgetConsExprNlhdlrData(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr         /**< nonlinear handler */
);

/** gives nonlinear handler expression data
 *
 * @return NULL if expr has not been detected by nlhdlr or nlhdlr did not store data
 */
SCIP_EXPORT
SCIP_CONSEXPR_NLHDLREXPRDATA* SCIPgetConsExprNlhdlrExprData(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_CONSEXPR_EXPR*        expr           /**< expression */
);

/** returns whether nonlinear handler implements the interval evaluation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrInteval(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the reverse propagation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrReverseProp(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the separation initialization callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrInitSepa(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the separation deinitialization callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrExitSepa(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the enforcement callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrEnfo(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the estimator callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrEstimate(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** call the detect callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRDETECT(SCIPdetectConsExprNlhdlr);

/** call the auxiliary evaluation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(SCIPevalauxConsExprNlhdlr);

/** calls the interval evaluation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(SCIPintevalConsExprNlhdlr);

/** calls the reverse propagation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(SCIPreversepropConsExprNlhdlr);

/** calls the separation initialization callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(SCIPinitsepaConsExprNlhdlr);

/** calls the separation deinitialization callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(SCIPexitsepaConsExprNlhdlr);

/** calls the enforcement callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRENFO(SCIPenfoConsExprNlhdlr);

/** calls the estimator callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(SCIPestimateConsExprNlhdlr);

/** @} */

/**@name Quadratic expression functions */
/**@{ */

/** gives the coefficients and expressions that define a quadratic expression
 *
 * It can return the constant part, the number, arguments, and coefficients of the purely linear part
 * and the number of quadratic terms and bilinear terms.
 * Note that for arguments that appear in the quadratic part, a linear coefficient is
 * stored with the quadratic term.
 * Use SCIPgetConsExprQuadraticQuadTermData() and SCIPgetConsExprQuadraticBilinTermData()
 * to access the data for a quadratic or bilinear term.
 *
 * This function returns pointers to internal data in linexprs and lincoefs.
 * The user must not change this data.
 */
SCIP_EXPORT
void SCIPgetConsExprQuadraticData(
   SCIP_CONSEXPR_QUADEXPR*       quaddata,         /**< quadratic coefficients data */
   SCIP_Real*                    constant,         /**< buffer to store constant term, or NULL */
   int*                          nlinexprs,        /**< buffer to store number of expressions that appear linearly, or NULL */
   SCIP_CONSEXPR_EXPR***         linexprs,         /**< buffer to store pointer to array of expressions that appear linearly, or NULL */
   SCIP_Real**                   lincoefs,         /**< buffer to store pointer to array of coefficients of expressions that appear linearly, or NULL */
   int*                          nquadexprs,       /**< buffer to store number of expressions in quadratic terms, or NULL */
   int*                          nbilinexprs       /**< buffer to store number of bilinear expressions terms, or NULL */
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
void SCIPgetConsExprQuadraticQuadTermData(
   SCIP_CONSEXPR_QUADEXPR*       quaddata,         /**< quadratic coefficients data */
   int                           termidx,          /**< index of quadratic term */
   SCIP_CONSEXPR_EXPR**          expr,             /**< buffer to store pointer to argument expression (the 'x') of this term, or NULL */
   SCIP_Real*                    lincoef,          /**< buffer to store linear coefficient of variable, or NULL */
   SCIP_Real*                    sqrcoef,          /**< buffer to store square coefficient of variable, or NULL */
   int*                          nadjbilin,        /**< buffer to store number of bilinear terms this variable is involved in, or NULL */
   int**                         adjbilin,         /**< buffer to store pointer to indices of associated bilinear terms, or NULL */
   SCIP_CONSEXPR_EXPR**          sqrexpr           /**< buffer to store pointer to square expression (the 'x^2') of this term or NULL if no square expression, or NULL */
   );

/** gives the data of a bilinear expression term
 *
 * For a term a*expr1*expr2, returns
 * expr1, expr2, a, and the position of the quadratic expression term that uses expr2 in the quadratic expressions quadexprterms.
 */
SCIP_EXPORT
void SCIPgetConsExprQuadraticBilinTermData(
   SCIP_CONSEXPR_QUADEXPR*       quaddata,         /**< quadratic coefficients data */
   int                           termidx,          /**< index of bilinear term */
   SCIP_CONSEXPR_EXPR**          expr1,            /**< buffer to store first factor, or NULL */
   SCIP_CONSEXPR_EXPR**          expr2,            /**< buffer to store second factor, or NULL */
   SCIP_Real*                    coef,             /**< buffer to coefficient, or NULL */
   int*                          pos2,             /**< buffer to position of expr2 in quadexprterms array of quadratic expression, or NULL */
   SCIP_CONSEXPR_EXPR**          prodexpr          /**< buffer to store pointer to expression that is product if first and second factor, or NULL */
   );

/** returns whether all expressions that are used in a quadratic expression are variable expression
 *
 * @return TRUE iff all linexprs and quadexprterms[.].expr in quaddata are variable expressions
 */
SCIP_EXPORT
SCIP_Bool SCIPareConsExprQuadraticExprsVariables(
   SCIP_CONSEXPR_QUADEXPR*       quaddata          /**< quadratic coefficients data */
   );

/** evaluates quadratic term in a solution w.r.t. auxiliary variables
 *
 * \note This assumes that for every expr used in the quadratic data, an auxiliary variable is available.
 */
SCIP_EXPORT
SCIP_Real SCIPevalConsExprQuadraticAux(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_QUADEXPR* quaddata,         /**< quadratic form */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL for LP solution */
   );

/** prints quadratic expression */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression conshdlr */
   SCIP_CONSEXPR_QUADEXPR* quaddata          /**< quadratic form */
   );

/** Checks the curvature of the quadratic function, x^T Q x + b^T x stored in quaddata
 *
 * For this, it builds the matrix Q and computes its eigenvalues using LAPACK; if Q is
 * - semidefinite positive -> provided is set to sepaunder
 * - semidefinite negative -> provided is set to sepaover
 * - otherwise -> provided is set to none
 *
 * If assumevarfixed is given and some entries of x correspond to variables present in
 * this hashmap, then the corresponding rows and columns are ignored in the matrix Q.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprQuadraticCurvature(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_QUADEXPR* quaddata,         /**< quadratic coefficients data */
   SCIP_EXPRCURV*          curv,             /**< pointer to store the curvature of quadratics */
   SCIP_HASHMAP*           assumevarfixed    /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   );

/** @} */


/** computes a facet of the convex or concave envelope of a vertex polyhedral function
 *
 * If \f$ f(x) \f$ is vertex-polyhedral, then \f$ g \f$ is a convex underestimator if and only if
 * \f$ g(v^i) \leq f(v^i), \forall i \f$, where \f$ \{ v^i \}_{i = 1}^{2^n} \subseteq \mathbb R^n \f$ are the vertices
 * of the domain of \f$ x \f$, \f$ [\ell,u] \f$. Hence, we can compute a linear underestimator by solving the following
 * LP (we don't necessarily get a facet of the convex envelope, see Technical detail below):
 *
 * \f{align*}{
 *              \max \, & \alpha^T x^* + \beta \\
 *     s.t. \; & \alpha^T v^i + \beta \le f(v^i), \, \forall i = 1, \ldots, 2^n
 * \f}
 *
 * In principle, one would need to update the LP whenever the domain changes. However, \f$ [\ell,u] = T([0, 1]^n) \f$,
 * where \f$ T \f$ is an affine linear invertible transformation given by \f$ T(y)_i = (u_i - \ell_i) y_i + \ell_i \f$.
 * Working with the change of variables \f$ x = T(y) \f$ allows us to keep the constraints of the LP, even if the domain
 * changes. Indeed, after the change of variables, the problem is: find an affine underestimator \f$ g \f$ such that \f$
 * g(T(y)) \le f(T(y)) \f$, for all \f$ y \in [0, 1]^n \f$. Now \f$ f(T(y)) \f$ is componentwise affine, but still
 * satisfies that \f$ g \f$ is a valid underestimator if and only if \f$ g(T(u)) \leq f(T(u)), \forall u \in \{0, 1\}^n
 * \f$. So we now look for \f$ \bar g(y) := g(T(y)) = g(((u_i - \ell_i) y_i + \ell_i)_i) = \bar \alpha^T y + \bar \beta
 * \f$, where \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i \f$ and \f$ \bar \beta = \sum_i \alpha_i \ell_i + \beta \f$. So
 * we find \f$ \bar g \f$ by solving the LP:
 *
 * \f{align*}{
 *              \max \, & \bar \alpha^T T^{-1}(x^*) + \bar \beta \\
 *     s.t. \; & \bar \alpha^T u + \bar \beta \le f(T(u)), \, \forall u \in \{0, 1\}^n
 * \f}
 *
 * and recover \f$ g \f$ by solving \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i, \bar \beta = \sum_i \alpha_i \ell_i +
 * \beta \f$. Notice that \f$ f(T(u^i)) = f(v^i) \f$ so the right hand side doesn't change after the change of variables.
 *
 * Furthermore, the LP has more constraints than variables, so we solve its dual:
 * \f{align*}{
 *              \min \, & \sum_i \lambda_i f(v^i) \\
 *     s.t. \; & \sum_i \lambda_i u^i = T^{-1}(x^*) \\
 *             & \sum_i \lambda_i = 1 \\
 *             & \forall i, \, \lambda_i \geq 0
 * \f}
 *
 * In case we look for an overestimate, we do exactly the same, but have to maximize in the dual LP instead
 * of minimize.
 *
 * #### Technical and implementation details
 * -# \f$ U \f$ has exponentially many variables, so we only apply this separator for \f$ n \leq 14 \f$.
 * -# If the bounds are not finite, there is no underestimator. Also, \f$ T^{-1}(x^*) \f$ must be in the domain,
 * otherwise the dual is infeasible.
 * -# After a facet is computed, we check whether it is a valid facet (i.e. we check \f$ \alpha^T v + \beta \le f(v) \f$
 *  for every vertex \f$ v \f$). If we find a violation of at most ADJUSTFACETFACTOR * SCIPlpfeastol, then we weaken \f$
 *  \beta \f$ by this amount, otherwise, we discard the cut.
 * -# If a variable is fixed within tolerances, we replace it with its value and compute the facet of the remaining
 * expression. Note that since we are checking the cut for validity, this will never produce wrong result.
 * -# If \f$ x^* \f$ is in the boundary of the domain, then the LP has infinitely many solutions, some of which might
 * have very bad numerical properties. For this reason, we perturb \f$ x^* \f$ to be in the interior of the region.
 * Furthermore, for some interior points, there might also be infinitely many solutions (e.g. for \f$ x y \f$ in \f$
 * [0,1]^2 \f$ any point \f$ (x^*, y^*) \f$ such that \f$ y^* = 1 - x^* \f$ has infinitely many solutions). For this
 * reason, we perturb any given \f$ x^* \f$. The idea is to try to get a facet of the convex/concave envelope. This only
 * happens when the solution has \f$ n + 1 \f$ non zero \f$ \lambda \f$'s (i.e. the primal has a unique solution).
 * -# We need to compute \f$ f(v^i) \f$ for every vertex of \f$ [\ell,u] \f$. A vertex is encoded by a number between 0
 * and \f$ 2^n - 1 \f$, via its binary representation (0 bit is lower bound, 1 bit is upper bound), so we can compute
 * all these values by iterating between 0 and \f$ 2^n - 1 \f$.
 * -# To check that the computed cut is valid we do the following: we use a gray code to loop over the vertices
 * of the box domain w.r.t. unfixed variables in order to evaluate the underestimator. To ensure the validity of the
 * underestimator, we check whether \f$ \alpha v^i + \beta \le f(v^i) \f$ for every vertex \f$ v^i \f$ and adjust
 * \f$ \beta \f$ if the maximal violation is small.
 *
 * @todo the solution is a facet if all variables of the primal have positive reduced costs (i.e. the solution is
 * unique). In the dual, this means that there are \f$ n + 1 \f$ variables with positive value. Can we use this or some
 * other information to handle any of both cases (point in the boundary or point in the intersection of polytopes
 * defining different pieces of the convex envelope)? In the case where the point is in the boundary, can we use that
 * information to maybe solve another to find a facet? How do the polytopes defining the pieces where the convex
 * envelope is linear looks like, i.e, given a point in the interior of a facet of the domain, does the midpoint of the
 * segment joining \f$ x^* \f$ with the center of the domain, always belongs to the interior of one of those polytopes?
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeFacetVertexPolyhedral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_DECL_VERTEXPOLYFUN((*function)),     /**< pointer to vertex polyhedral function */
   void*                 fundata,            /**< data for function evaluation (can be NULL) */
   SCIP_Real*            xstar,              /**< point to be separated */
   SCIP_Real*            box,                /**< box where to compute facet: should be lb_1, ub_1, lb_2, ub_2... */
   int                   nallvars,           /**< half of the length of box */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an array of length at least nallvars */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
);

/** given three points, constructs coefficient of equation for hyperplane generated by these three points
 * Three points a, b, and c are given.
 * Computes coefficients alpha, beta, gamma, and delta, such that a, b, and c, satisfy
 * alpha * x1 + beta * x2 + gamma * x3 = delta and gamma >= 0.0.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeHyperplaneThreePoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             a1,                 /**< first coordinate of a */
   SCIP_Real             a2,                 /**< second coordinate of a */
   SCIP_Real             a3,                 /**< third coordinate of a */
   SCIP_Real             b1,                 /**< first coordinate of b */
   SCIP_Real             b2,                 /**< second coordinate of b */
   SCIP_Real             b3,                 /**< third coordinate of b */
   SCIP_Real             c1,                 /**< first coordinate of c */
   SCIP_Real             c2,                 /**< second coordinate of c */
   SCIP_Real             c3,                 /**< third coordinate of c */
   SCIP_Real*            alpha,              /**< coefficient of first coordinate */
   SCIP_Real*            beta,               /**< coefficient of second coordinate */
   SCIP_Real*            gamma_,             /**< coefficient of third coordinate */
   SCIP_Real*            delta               /**< constant right-hand side */
   );

#ifdef __cplusplus
}
#endif

#endif
