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
   SCIP_DECL_CONSEXPR_EXPRSEPA((*sepa)),     /**< separation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRESTIMATE((*estimate))  /**< estimator callback (can be NULL) */
);

/** set the branching score callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrBranchscore(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRBRANCHSCORE((*brscore)) /**< branching score callback (can be NULL) */
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

/** returns whether expression handler implements the separation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrSepa(
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

/** calls separator method of expression handler to separate a given solution */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRSEPA(SCIPsepaConsExprExprHdlr);

/** calls the expression branching score callback */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_EXPRBRANCHSCORE(SCIPbranchscoreConsExprExprHdlr);

/** increments the branching score count of an expression handler */
SCIP_EXPORT
void SCIPincrementConsExprExprHdlrNBranchScore(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr
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

/** appends child to the children list of expr */
SCIP_EXPORT
SCIP_RETCODE SCIPappendConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_EXPR*   child               /**< expression to be appended */
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

/** returns the partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error) */
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

/**@} */  /* end of differentiation methods */

/** returns the activity of the expression
 *
 * The caller needs to make sure that the activity is valid.
 * For expression and nonlinear handlers, this is made sure when the following callbacks are called:
 * - interval evaluation (intervals for children only)
 * - reverse propagation
 * - monotonicity computation
 * - convexity detection
 * - structure detection
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

/** tightens the bounds of an expression and stores the result in the expression interval; variables in variable
 *  expression will be tightened immediately if SCIP is in a stage above SCIP_STAGE_TRANSFORMED
 *
 *  If a reversepropqueue is given, then the expression will be added to the queue if its bounds could be tightened without detecting infeasibility.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPtightenConsExprExprInterval(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be tightened */
   SCIP_INTERVAL           newbounds,        /**< new bounds for the expression */
   SCIP_Bool               force,            /**< force tightening even if below bound strengthening tolerance */
   SCIP_QUEUE*             reversepropqueue, /**< reverse propagation queue, or NULL if not in reverse propagation */
   SCIP_Bool*              cutoff,           /**< buffer to store whether a node's bounds were propagated to an empty interval */
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

/** adds branching score to an expression
 *
 * Adds a score to the expression-specific branching score.
 * The branchscoretag argument is used to identify whether the score in the expression needs to be reset before adding a new score.
 * In an expression with children, the scores are distributed to its children.
 * In an expression that is a variable, the score may be used to identify a variable for branching.
 */
SCIP_EXPORT
void SCIPaddConsExprExprBranchScore(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression where to add branching score */
   unsigned int            branchscoretag,   /**< tag to identify current branching scores */
   SCIP_Real               branchscore       /**< branching score to add to expression */
   );

/** returns the hash value of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprHash(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   unsigned int*           hashval           /**< pointer to store the hash value */
   );

/** creates and gives the auxiliary variable for a given expression
 *
 * @note if auxiliary variable already present for that expression, then only returns this variable
 * @note for a variable expression it returns the corresponding variable
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExprAuxVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR**            auxvar              /**< buffer to store pointer to auxiliary variable, or NULL */
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

/** reformulate an expression; this functions works similar as SCIPsimplifyConsExprExpr() but instead of calling the
 *  simplify callback of an expression handler it iterates through all nonlinear handlers and uses the reformulation
 *  callback
 */
SCIP_EXPORT
SCIP_RETCODE SCIPreformulateConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSEXPR_EXPR*     rootexpr,         /**< expression to be simplified */
   SCIP_CONSEXPR_EXPR**    refrootexpr,      /**< buffer to store reformulated expression */
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

/** returns the monotonicity of an expression w.r.t. to a given child */
SCIP_EXPORT
SCIP_MONOTONE SCIPgetConsExprExprMonotonicity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   int                   childidx            /**< index of child */
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
 * to be at least of size the number of unique variables in the expression which is given by SCIpgetConsExprExprNVars()
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

/** collects all bilinear terms for a given set of constraints
 *
 * @note This method should only be used for unit tests that depend on SCIPgetConsExprBilinTerms()
 *       or SCIPgetConsExprBilinTermAuxar().
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
 * @note The value of auxvars[i] might be NULL, which indicates that xs[i] * ys[i] does not have an auxiliary variable.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprBilinTerms(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_VAR**                 xs,             /**< array to store first variables (of size >= SCIPgetConsExprNBilinTerms()) */
   SCIP_VAR**                 ys,             /**< array to store second variables (of size >= SCIPgetConsExprNBilinTerms()) */
   SCIP_VAR**                 auxvars         /**< array to store auxiliary variables (of size >= SCIPgetConsExprNBilinTerms()) */
   );

/** returns the auxiliary variable of a bilinear term, if it exists
 *
 * @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprBilinTermAuxar(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_VAR*                  x,              /**< first variable */
   SCIP_VAR*                  y,              /**< second variable */
   SCIP_VAR**                 auxvar,         /**< pointer to store auxiliary variable (might be NULL) */
   SCIP_Bool*                 found           /**< pointer to store whether the bilinear term xy exists */
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
 *  - nupgdconss      : pointer to store number of constraints that replace this constraint
 *  - upgdconss       : array to store constraints that replace this constraint
 *  - upgdconsssize   : length of the provided upgdconss array
 */
#define SCIP_DECL_EXPRCONSUPGD(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONS* cons, \
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
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
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

/** @} */

/**@name Expression Constraint Methods */
/**@{ */

/** returns the expression of the given expression constraint */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR* SCIPgetExprConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
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
   int                   nconss,             /**< total number of constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected while creating the auxiliary vars */
   );

/** @} */

/**@name Nonlinear Handler Methods */
/**@{ */

/** creates the nonlinearity handler and includes it into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprNlhdlrBasic(
   SCIP*                       scip,         /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_NLHDLR**      nlhdlr,       /**< buffer where to store nonlinear handler */
   const char*                 name,         /**< name of nonlinear handler (must not be NULL) */
   const char*                 desc,         /**< description of nonlinear handler (can be NULL) */
   int                         priority,     /**< priority of nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRDETECT((*detect)), /**< structure detection callback of nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLREVALAUX((*evalaux)), /**< auxiliary evaluation callback of nonlinear handler */
   SCIP_CONSEXPR_NLHDLRDATA*   data          /**< data of nonlinear handler (can be NULL) */
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

/** set the reformulate callback of a nonlinear handler */
void SCIPsetConsExprNlhdlrReformulate(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE((*reformulate)) /**< reformulation callback */
   );

/** set the propagation callbacks of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRINTEVAL((*inteval)), /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
);

/** set the separation callbacks of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRINITSEPA((*initsepa)), /**< separation initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRSEPA((*sepa)),         /**< separation callback (can be NULL if estimate is not NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRESTIMATE((*estimate)), /**< estimation callback (can be NULL if sepa is not NULL) */
   SCIP_DECL_CONSEXPR_NLHDLREXITSEPA((*exitsepa))  /**< separation deinitialization callback (can be NULL) */
);

/** set the branching score callback of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrBranchscore(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE((*branchscore)) /**< branching score callback */
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

/** gives priority of nonlinear handler */
SCIP_EXPORT
int SCIPgetConsExprNlhdlrPriority(
   SCIP_CONSEXPR_NLHDLR*      nlhdlr         /**< nonlinear handler */
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

/** returns whether nonlinear handler implements the reformulation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrReformulate(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
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

/** returns whether nonlinear handler implements the separation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrSepa(
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

/** calls the reformulation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE(SCIPreformulateConsExprNlhdlr);

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

/** calls the separation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRSEPA(SCIPsepaConsExprNlhdlr);

/** calls the estimator callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(SCIPestimateConsExprNlhdlr);

/** calls the nonlinear handler branching score callback */
SCIP_EXPORT
SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE(SCIPbranchscoreConsExprNlHdlr);

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

#ifdef __cplusplus
}
#endif

#endif
