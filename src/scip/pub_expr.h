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

/**@name Expression Handler Methods */
/**@{ */

// TODO everything with SCIP* move to scip_expr.h
// TODO some cons_nonlinear specific move to cons_nonlinear.h


/** creates the handler for an expression handler and includes it into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprExprHdlrBasic(
   SCIP*                       scip,         /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,     /**< expression constraint handler */
   SCIP_EXPRHDLR**    exprhdlr,     /**< buffer where to store expression handler */
   const char*                 name,         /**< name of expression handler (must not be NULL) */
   const char*                 desc,         /**< description of expression handler (can be NULL) */
   unsigned int                precedence,   /**< precedence of expression operation (used for printing) */
   SCIP_DECL_EXPREVAL((*eval)),     /**< point evaluation callback (can not be NULL) */
   SCIP_EXPRHDLRDATA* data          /**< data of expression handler (can be NULL) */
   );

/** set the expression handler callbacks to copy and free an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeHdlr(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_EXPRCOPYHDLR((*copyhdlr)), /**< handler copy callback (can be NULL) */
   SCIP_DECL_EXPRFREEHDLR((*freehdlr))  /**< handler free callback (can be NULL) */
);

/** set the expression handler callbacks to copy and free expression data */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_EXPRCOPYDATA((*copydata)), /**< expression data copy callback (can be NULL for expressions without data) */
   SCIP_DECL_EXPRFREEDATA((*freedata))  /**< expression data free callback (can be NULL if data does not need to be freed) */
);

/** set the print callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrPrint(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRPRINT((*print))    /**< print callback (can be NULL) */
);

/** set the parse callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrParse(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRPARSE((*parse))    /**< parse callback (can be NULL) */
);

/** set the curvature detection callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrCurvature(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRCURVATURE((*curvature)) /**< curvature detection callback (can be NULL) */
);

/** set the monotonicity detection callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrMonotonicity(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRMONOTONICITY((*monotonicity)) /**< monotonicity detection callback (can be NULL) */
);

/** set the integrality detection callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrIntegrality(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRINTEGRALITY((*integrality)) /**< integrality detection callback (can be NULL) */
);

/** set the hash callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrHash(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRHASH((*hash))      /**< hash callback (can be NULL) */
);

/** set the compare callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrCompare(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRCOMPARE((*compare))/**< compare callback (can be NULL) */
);

/** set derivative evaluation callbacks of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrDiff(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRBWDIFF((*bwdiff)),  /**< backward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRFWDIFF((*fwdiff)),  /**< forward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRBWFWDIFF((*bwfwdiff))/**< backward-forward derivative evaluation callback (can be NULL) */
);

/** set the interval evaluation callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrIntEval(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRINTEVAL((*inteval))/**< interval evaluation callback (can be NULL) */
);

/** set the simplify callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrSimplify(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRSIMPLIFY((*simplify))  /**< simplify callback (can be NULL) */
);

/** set the reverse propagation callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrReverseProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRREVERSEPROP((*reverseprop))/**< reverse propagation callback (can be NULL) */
);

/** set the separation and estimation callbacks of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetConsExprExprHdlrSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_EXPRINITSEPA((*initsepa)), /**< separation initialization callback (can be NULL) */
   SCIP_DECL_EXPREXITSEPA((*exitsepa)), /**< separation deinitialization callback (can be NULL) */
   SCIP_DECL_EXPRESTIMATE((*estimate))  /**< estimator callback (can be NULL) */
);

/** gives expression handlers */
SCIP_EXPORT
SCIP_EXPRHDLR** SCIPgetConsExprExprHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
);

/** gives number of expression handlers */
SCIP_EXPORT
int SCIPgetConsExprExprNHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
);

/** returns an expression handler of a given name (or NULL if not found) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPfindConsExprExprHdlr(
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   const char*                name           /**< name of expression handler */
   );

/** returns expression handler for variable expressions */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetConsExprExprHdlrVar(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for constant value expressions */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetConsExprExprHdlrValue(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** gives the value of directional derivative from the last evaluation of a directional derivative of expression (or SCIP_INVALID if there was an error) */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprDot(
   SCIP_EXPR*     expr              /**< expression */
   );

/** returns expression handler for sum expressions */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetConsExprExprHdlrSum(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for product expressions */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetConsExprExprHdlrProduct(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for power expressions */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetConsExprExprHdlrPower(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for signed power expressions */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetConsExprExprHdlrSignPower(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for exponential expressions */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetConsExprExprHdlrExponential(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   );

/** returns expression handler for logarithm expressions */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetConsExprExprHdlrLog(
        SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
);

/** gives the name of an expression handler */
SCIP_EXPORT
const char* SCIPgetConsExprExprHdlrName(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
);

/** gives the description of an expression handler (can be NULL) */
SCIP_EXPORT
const char* SCIPgetConsExprExprHdlrDescription(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
);

/** gives the precedence of an expression handler */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprHdlrPrecedence(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
);

/** gives the data of an expression handler */
SCIP_EXPORT
SCIP_EXPRHDLRDATA* SCIPgetConsExprExprHdlrData(
   SCIP_EXPRHDLR*    exprhdlr      /**< expression handler */
);

/** returns whether expression handler implements the print callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrPrint(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the backward differentiation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrBwdiff(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the interval evaluation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrIntEval(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the estimator callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrEstimate(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the simplification callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrSimplify(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the curvature callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrCurvature(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the reverse propagation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrReverseProp(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the initialization callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrInitSepa(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the deinitialization callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrExitSepa(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** returns whether expression handler implements the branching score callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprExprHdlrBranchingScore(
   SCIP_EXPRHDLR*    exprhdlr       /**< expression handler */
   );

/** calls the print callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRPRINT(SCIPprintConsExprExprHdlr);

/** calls the parse callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRPARSE(SCIPparseConsExprExprHdlr);

/** calls the expression hash callback */
SCIP_EXPORT
SCIP_DECL_EXPRHASH(SCIPhashConsExprExprHdlr);

/** calls the expression compare callback */
SCIP_EXPORT
SCIP_DECL_EXPRCOMPARE(SCIPcompareConsExprExprHdlr);

/** calls the backward-differentiation callback of an expression handler
 *
 * further, allows to different w.r.t. given expression and children values
 */
SCIP_EXPORT
SCIP_RETCODE SCIPbwdiffConsExprExprHdlr(
   SCIP*                      scip,         /**< SCIP data structure */
   SCIP_EXPR*        expr,         /**< expression */
   int                        childidx,     /**< index of child w.r.t. which to compute derivative */
   SCIP_Real*                 derivative,   /**< buffer to store value of derivative */
   SCIP_Real*                 childrenvals, /**< values for children, or NULL if values stored in children should be used */
   SCIP_Real                  exprval       /**< value for expression, used only if childrenvals is not NULL */
);

/** calls the backward-forward differentiation callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPbwfwdiffConsExprExprHdlr(
   SCIP*                      scip,         /**< SCIP data structure */
   SCIP_EXPR*        expr,         /**< expression */
   int                        childidx,     /**< index of child w.r.t. which to compute derivative */
   SCIP_Real*                 derivative    /**< buffer to store value of the backward-forward derivative */
);

/** calls the forward differentiation callback of an expression handler */
SCIP_EXPORT
SCIP_RETCODE SCIPfwdiffConsExprExprHdlr(
   SCIP*                      scip,         /**< SCIP data structure */
   SCIP_EXPR*        expr,         /**< expression */
   SCIP_Real*                 derivative    /**< buffer to store value of the forward derivative */
);


/** calls the evaluation callback of an expression handler
 *
 * further, allows to evaluate w.r.t. given children values
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalConsExprExprHdlr(
   SCIP*                      scip,         /**< SCIP data structure */
   SCIP_EXPR*        expr,         /**< expression */
   SCIP_Real*                 val,          /**< buffer store value of expression */
   SCIP_Real*                 childrenvals, /**< values for children, or NULL if values stored in children should be used */
   SCIP_SOL*                  sol           /**< solution that is evaluated (used by the var-expression) */
);

/** calls the expression interval evaluation callback */
SCIP_EXPORT
SCIP_DECL_EXPRINTEVAL(SCIPintevalConsExprExprHdlr);

/** calls estimator method of expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRESTIMATE(SCIPestimateConsExprExprHdlr);

/** calls the simplification method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRSIMPLIFY(SCIPsimplifyConsExprExprHdlr);

/** calls the curvature check method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRCURVATURE(SCIPcurvatureConsExprExprHdlr);

/** calls the expression callback for reverse propagation */
SCIP_EXPORT
SCIP_DECL_EXPRREVERSEPROP(SCIPreversepropConsExprExprHdlr);

/** calls the separation initialization method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRINITSEPA(SCIPinitsepaConsExprExprHdlr);

/** calls the separation deinitialization method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPREXITSEPA(SCIPexitsepaConsExprExprHdlr);

/** increments the branching score count of an expression handler */
SCIP_EXPORT
void SCIPincrementConsExprExprHdlrNBranchScore(
   SCIP_EXPRHDLR*    exprhdlr
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
   SCIP_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_EXPRDATA* exprdata,         /**< expression data (expression assumes ownership) */
   int                     nchildren,        /**< number of children */
   SCIP_EXPR**    children          /**< children (can be NULL if nchildren is 0) */
   );

/** creates and captures an expression with given expression data and up to two children */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExpr2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_EXPRDATA* exprdata,         /**< expression data */
   SCIP_EXPR*     child1,           /**< first child (can be NULL) */
   SCIP_EXPR*     child2            /**< second child (can be NULL) */
   );

/** creates and captures an expression from a node in an (old-style) expression graph */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExpr3(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_EXPRGRAPH*         exprgraph,        /**< expression graph */
   SCIP_EXPRGRAPHNODE*     node              /**< expression graph node */
   );

/** creates and captures an expression representing a quadratic function */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExprQuadratic(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_EXPR**    expr,             /**< pointer where to store expression */
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
   SCIP_EXPR**    expr,             /**< pointer where to store expression */
   int                     nfactors,         /**< number of factors in monomial */
   SCIP_VAR**              vars,             /**< variables in the monomial */
   SCIP_Real*              exponents         /**< exponent in each factor, or NULL if all 1.0 */
   );

/** appends child to the children list of expr */
SCIP_EXPORT
SCIP_RETCODE SCIPappendConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_EXPR*   child               /**< expression to be appended */
   );

/** remove all children of expr
 *
 * only use if you really know what you are doing
 */
SCIP_EXPORT
SCIP_RETCODE SCIPremoveConsExprExprChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr                /**< expression */
   );

/** overwrites/replaces a child of an expressions
 *
 * @note the old child is released and the newchild is captured, unless they are the same (=same pointer)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPreplaceConsExprExprChild(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPR*     expr,             /**< expression which is going to replace a child */
   int                     childidx,         /**< index of child being replaced */
   SCIP_EXPR*     newchild          /**< the new child */
   );

/** duplicates the given expression */
SCIP_EXPORT
SCIP_RETCODE SCIPduplicateConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< original expression */
   SCIP_EXPR**  copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_Bool             copychildren        /**< whether children (and all successors) should be copied, too */
   );

/** gets the number of times the expression is currently captured */
SCIP_EXPORT
int SCIPgetConsExprExprNUses(
   SCIP_EXPR*   expr               /**< expression */
   );

/** captures an expression (increments usage count) */
SCIP_EXPORT
void SCIPcaptureConsExprExpr(
   SCIP_EXPR*   expr                /**< expression to be captured */
   );

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_EXPORT
SCIP_RETCODE SCIPreleaseConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**  expr                /**< pointer to expression to be released */
   );

/** gives the number of children of an expression */
SCIP_EXPORT
int SCIPgetConsExprExprNChildren(
   SCIP_EXPR*   expr                /**< expression */
   );

/** gives the children of an expression (can be NULL if no children) */
SCIP_EXPORT
SCIP_EXPR** SCIPgetConsExprExprChildren(
   SCIP_EXPR*   expr                /**< expression */
   );

/** gets the handler of an expression
 *
 * This identifies the type of the expression (sum, variable, ...).
 */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetConsExprExprHdlr(
   SCIP_EXPR*   expr                /**< expression */
   );

/** gets the expression data of an expression */
SCIP_EXPORT
SCIP_EXPRDATA* SCIPgetConsExprExprData(
   SCIP_EXPR*   expr                /**< expression */
   );

/** returns whether an expression is a variable expression */
SCIP_EXPORT
SCIP_Bool SCIPisConsExprExprVar(
   SCIP_EXPR*   expr                /**< expression */
   );

/** returns whether an expression is a value expression */
SCIP_EXPORT
SCIP_Bool SCIPisConsExprExprValue(
   SCIP_EXPR*   expr                /**< expression */
   );

/** returns the variable used for linearizing a given expression (return value might be NULL)
 *
 * @note for variable expression it returns the corresponding variable
 */
SCIP_EXPORT
SCIP_VAR* SCIPgetConsExprExprAuxVar(
   SCIP_EXPR*   expr                /**< expression */
   );

/** sets the expression data of an expression
 *
 * The pointer to possible old data is overwritten and the
 * freedata-callback is not called before.
 * This function is intended to be used by expression handler.
 */
SCIP_EXPORT
void SCIPsetConsExprExprData(
   SCIP_EXPR*     expr,             /**< expression */
   SCIP_EXPRDATA* exprdata          /**< expression data to be set (can be NULL) */
   );

/** print an expression as info-message */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_EXPR*     expr,             /**< expression to be printed */
   FILE*                   file              /**< file to print to, or NULL for stdout */
   );

/** initializes printing of expressions in dot format to a give FILE* pointer */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExprDotInit(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_EXPRPRINTDOTDATA** dotdata,     /**< buffer to store dot printing data */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_EXPRPRINTDOT_WHAT whattoprint   /**< info on what to print for each expression */
   );

/** initializes printing of expressions in dot format to a file with given filename */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExprDotInit2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_EXPRPRINTDOTDATA** dotdata,     /**< buffer to store dot printing data */
   const char*             filename,         /**< name of file to print to */
   SCIP_EXPRPRINTDOT_WHAT whattoprint   /**< info on what to print for each expression */
   );

SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExprDot(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDOTDATA* dotdata,      /**< data as initialized by \ref SCIPprintConsExprExprDotInit() */
   SCIP_EXPR*     expr              /**< expression to be printed */
   );

/** finishes printing of expressions in dot format */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprExprDotFinal(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDOTDATA** dotdata      /**< buffer where dot printing data has been stored */
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
   SCIP_EXPR*     expr              /**< expression to be printed */
   );

/** prints structure of an expression a la Maple's dismantle */
SCIP_EXPORT
SCIP_RETCODE SCIPdismantleConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_EXPR*     expr              /**< expression to dismantle */
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
   SCIP_EXPR**  expr                /**< pointer to store the expr parsed */
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
   SCIP_EXPR*     expr,             /**< expression to be evaluated */
   SCIP_SOL*               sol,              /**< solution to be evaluated */
   unsigned int            soltag            /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/** gives the value from the last evaluation of an expression (or SCIP_INVALID if there was an eval error) */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprValue(
   SCIP_EXPR*     expr              /**< expression */
   );

/** sets the evaluation value */
SCIP_EXPORT
void SCIPsetConsExprExprEvalValue(
   SCIP_EXPR*     expr,             /**< expression */
   SCIP_Real               value,            /**< value to set */
   unsigned int            tag               /**< tag of solution that was evaluated, or 0 */
   );

/** gives the evaluation tag from the last evaluation, or 0 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprEvalTag(
   SCIP_EXPR*     expr              /**< expression */
   );

/** @name Differentiation methods
 * Gradients (Automatic differentiation Backward mode)
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
 * child[i]->derivative = the derivative of expr w.r.t child evaluated at point * expr->derivative
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
 * Note that, to compute this, we only need to know, for each expression, its partial derivatives w.r.t a given child
 * at a point. This is what the callback SCIP_DECL_EXPRBWDIFF should return.
 * Indeed, from child[i]->derivative = the derivative of expr w.r.t child evaluated at point * expr->derivative,
 * note that at the moment of processing a child, we already know expr->derivative, so the only
 * missing piece of information is 'the derivative of expr w.r.t child evaluated at point'.
 *
 * An equivalent way of interpreting the procedure is that expr->derivative stores the derivative of the root w.r.t expr.
 * This way, x->derivative and y->derivative will contain the partial derivatives of root w.r.t to the variable,
 * that is, the gradient. Note, however, that this analogy is only correct for leave expressions, since
 * the derivative value of an intermediate expression gets overwritten.
 *
 *
 * Hessian (Automatic differentiation Backward on Forward mode)
 * Computing the Hessian is more complicated since it is the derivative of the gradient, which is a function with more than one output.
 * We compute the Hessian by computing 'directions' of the Hessian, that is H*u for different 'u'
 * This is easy in general, since it is the gradient of the *scalar* function `grad f^T u`, that is, the directional derivative of f
 * in the direction u, D_u f.
 * This is easily computed via the so called forward mode.
 * Just as expr->derivative stores the partial derivative of the root w.r.t expr,
 * expr->dot stores the directional derivative of expr in the direction 'u'.
 * Then, by the chain rule, expr->dot = sum_(c : children) d_c expr * c->dot.
 * Starting with x_i->dot = u_i, we can compute expr->dot for every expression at the same time we evaluate expr.
 * Computing expr->dot is the purpose of the callback SCIP_DECL_EXPRFWDIFF.
 * Obviously, when this callback is called, the dots of all children are known
 * (just like evaluation, where the value of all children are known).
 *
 * Once we have this information, we compute the gradient of this function, following the same idea as before.
 * We define expr->bardot to be the directional derivative in direction u of the partial derivative of the root w.r.t expr `grad f^T u` w.r.t expr,
 * that is D_u (d_expr f) = D_u (expr->derivative).
 *
 * This way, x_i->bardot = D_u (d_(x_i) f) = e_i^T H_f u. Hence vars->bardot contain H_f u.
 * By the chain rule, product rule, and definition we have
 *
 * expr->bardot = D_u (d_expr f) =
 * D_u ( d_parent f * d_expr parent ) =
 * D_u( parent->derivative * d_expr parent ) =
 * d_expr parent * D_u (parent->derivative) + parent->derivative * D_u (d_expr parent) =
 * parent->bardot * d_expr parent + parent->derivative * D_u (d_expr parent)
 *
 * Note that we have computed parent->bardot and parent->derivative at this point,
 * while (d_expr parent) is the return of SCIP_DECL_EXPRBWDIFF.
 * Hence the only information we need to compute is D_u (d_expr parent).
 * This is the purpose of the callback SCIP_DECL_EXPRBWFWDIFF.
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
   SCIP_EXPR*     expr,             /**< expression to be evaluated */
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
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   );

/** returns the var's coordinate of Hu partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error)
 *
 * @note expression must belong to a constraint
 */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprPartialDiffGradientDir(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< root expression of constraint used in the last SCIPcomputeConsExprHessianDir() call */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   );

/** computes the hessian * v at a given point
 *
 * Evaluates children, if necessary.
 * Value can be received via SCIPgetConsExprExprPartialDiffGradientDir()
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeConsExprHessianDir(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONS*              cons,             /**< constraint for which we will compute directional derivative */
   SCIP_SOL*               sol,              /**< solution to be evaluated (NULL for the current LP solution) */
   SCIP_SOL*               direction,        /**< direction */
   unsigned int            soltag            /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/** returns the derivative stored in an expression (or SCIP_INVALID if there was an evaluation error) */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprDerivative(
   SCIP_EXPR*     expr              /**< expression */
   );

/** returns the difftag stored in an expression
 *
 * can be used to check whether partial derivative value is valid
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprDiffTag(
   SCIP_EXPR*     expr              /**< expression */
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
   SCIP_EXPR*     expr              /**< expression */
   );

/** returns the tag associated with the activity of the expression
 *
 * Can be compared with SCIPgetConsExprCurBoundsTag() and SCIPgetConsExprLastBoundRelaxTag()
 * to check whether the activity currently stored in this expression is current and valid, respectively.
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprActivityTag(
   SCIP_EXPR*     expr              /**< expression */
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
   SCIP_EXPR*     expr,             /**< expression */
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
   SCIP_EXPR*     expr              /**< expression */
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
   SCIP_EXPR*     expr,             /**< expression to be tightened */
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
   SCIP_EXPR*     expr              /**< expression to propagate again */
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
   SCIP_EXPR*     expr,             /**< expression where to add branching score */
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
   SCIP_EXPR**    exprs,            /**< expressions where to add branching score */
   int                     nexprs,           /**< number of expressions */
   SCIP_Real               violscore,        /**< violation score to add to expression */
   SCIP_SOL*               sol,              /**< current solution */
   SCIP_Bool*              success           /**< buffer to store whether at least one violscore was added */
   );

/** gives violation-branching score stored in expression, or 0.0 if no valid score has been stored */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprViolScore(
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_EXPR*     expr              /**< expression */
   );

/** returns the hash value of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprHash(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPR*     expr,             /**< expression */
   unsigned int*           hashval           /**< pointer to store the hash value */
   );

/** compare expressions
 * @return -1, 0 or 1 if expr1 <, =, > expr2, respectively
 * @note: The given expressions are assumed to be simplified.
 */
SCIP_EXPORT
int SCIPcompareConsExprExprs(
   SCIP_EXPR*   expr1,              /**< first expression */
   SCIP_EXPR*   expr2               /**< second expression */
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
   SCIP_EXPR*     rootexpr,         /**< expression to be simplified */
   SCIP_EXPR**    simplified,       /**< buffer to store simplified expression */
   SCIP_Bool*              changed,          /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*              infeasible        /**< buffer to store whether infeasibility has been detected */
   );

/** sets the curvature of an expression */
SCIP_EXPORT
void SCIPsetConsExprExprCurvature(
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_EXPRCURV         curvature           /**< curvature of the expression */
   );

/** returns the curvature of an expression
 *
 *  @note Call SCIPcomputeConsExprExprCurvature before calling this function.
 */
SCIP_EXPORT
SCIP_EXPRCURV SCIPgetConsExprExprCurvature(
   SCIP_EXPR*   expr                /**< expression */
   );

/** computes the curvature of a given expression and all its subexpressions
 *
 *  @note this function also evaluates all subexpressions w.r.t. current variable bounds
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeConsExprExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr                /**< expression */
   );

/** computes the monotonicity of an expression w.r.t. to a given child */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprExprMonotonicity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< expression */
   int                   childidx,           /**< index of child */
   SCIP_MONOTONE*        monotonicity        /**< buffer to store monotonicity */
   );

/** returns the number of positive rounding locks of an expression */
SCIP_EXPORT
int SCIPgetConsExprExprNLocksPos(
   SCIP_EXPR*   expr                /**< expression */
   );

/** returns the number of negative rounding locks of an expression */
SCIP_EXPORT
int SCIPgetConsExprExprNLocksNeg(
   SCIP_EXPR*   expr                /**< expression */
   );

/** computes integrality information of a given expression and all its subexpressions; the integrality information can
 * be accessed via SCIPisConsExprExprIntegral()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeConsExprExprIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   expr                /**< expression */
   );

/** returns whether an expression is integral */
SCIP_EXPORT
SCIP_Bool SCIPisConsExprExprIntegral(
   SCIP_EXPR*   expr                /**< expression */
   );

/** number of nonlinear handlers whose activity computation and propagation methods depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprNPropUsesActivity(
   SCIP_EXPR*   expr                /**< expression */
   );

/** number of nonlinear handlers whose separation methods (estimate or enforcement) depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprNSepaUsesActivity(
   SCIP_EXPR*   expr                /**< expression */
   );

/** number of nonlinear handlers whose separation methods (estimate or enforcement) use auxiliary variable of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
SCIP_EXPORT
unsigned int SCIPgetConsExprExprNAuxvarUses(
   SCIP_EXPR*   expr                /**< expression */
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
   SCIP_EXPR*   expr,             /**< expression */
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
   SCIP_EXPR*     expr,             /**< expression */
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
   SCIP_EXPR*     expr,             /**< expression */
   SCIP_EXPR**    varexprs,         /**< array to store all variable expressions */
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
   SCIP_EXPR*   expr,               /**< expression */
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
   SCIP_EXPR*   expr,               /**< expression */
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
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   );

/** @} */

/**@name Expression iterator */
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
 * It can also return the eigenvalues and the eigenvectors of the matrix Q when the quadratic is written
 * as x^T Q x + b^T x + c^T y + d, where c^T y defines the purely linear part.
 * Note, however, that to have access to them one needs to call SCIPgetConsExprQuadraticCurvature()
 * with storeeigeninfo equals to TRUE. If the eigen infortmation was not stored or it failed to be computed,
 * eigenvalues and eigenvectors will be set to NULL.
 *
 *
 * This function returns pointers to internal data in linexprs and lincoefs.
 * The user must not change this data.
 */
SCIP_EXPORT
void SCIPgetConsExprQuadraticData(
   SCIP_EXPR*           expr,
//   SCIP_QUADEXPR*       quaddata,         /**< quadratic coefficients data */
   SCIP_Real*                    constant,         /**< buffer to store constant term, or NULL */
   int*                          nlinexprs,        /**< buffer to store number of expressions that appear linearly, or NULL */
   SCIP_EXPR***         linexprs,         /**< buffer to store pointer to array of expressions that appear linearly, or NULL */
   SCIP_Real**                   lincoefs,         /**< buffer to store pointer to array of coefficients of expressions that appear linearly, or NULL */
   int*                          nquadexprs,       /**< buffer to store number of expressions in quadratic terms, or NULL */
   int*                          nbilinexprs,      /**< buffer to store number of bilinear expressions terms, or NULL */
   SCIP_Real**                   eigenvalues,      /**< buffer to store pointer to array of eigenvalues of Q, or NULL */
   SCIP_Real**                   eigenvectors      /**< buffer to store pointer to array of eigenvectors of Q, or NULL */
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
   SCIP_QUADEXPR*       quaddata,         /**< quadratic coefficients data */
   int                           termidx,          /**< index of quadratic term */
   SCIP_EXPR**          expr,             /**< buffer to store pointer to argument expression (the 'x') of this term, or NULL */
   SCIP_Real*                    lincoef,          /**< buffer to store linear coefficient of variable, or NULL */
   SCIP_Real*                    sqrcoef,          /**< buffer to store square coefficient of variable, or NULL */
   int*                          nadjbilin,        /**< buffer to store number of bilinear terms this variable is involved in, or NULL */
   int**                         adjbilin,         /**< buffer to store pointer to indices of associated bilinear terms, or NULL */
   SCIP_EXPR**          sqrexpr           /**< buffer to store pointer to square expression (the 'x^2') of this term or NULL if no square expression, or NULL */
   );

/** gives the data of a bilinear expression term
 *
 * For a term a*expr1*expr2, returns
 * expr1, expr2, a, and the position of the quadratic expression term that uses expr2 in the quadratic expressions quadexprterms.
 */
SCIP_EXPORT
void SCIPgetConsExprQuadraticBilinTermData(
   SCIP_QUADEXPR*       quaddata,         /**< quadratic coefficients data */
   int                           termidx,          /**< index of bilinear term */
   SCIP_EXPR**          expr1,            /**< buffer to store first factor, or NULL */
   SCIP_EXPR**          expr2,            /**< buffer to store second factor, or NULL */
   SCIP_Real*                    coef,             /**< buffer to coefficient, or NULL */
   int*                          pos2,             /**< buffer to position of expr2 in quadexprterms array of quadratic expression, or NULL */
   SCIP_EXPR**          prodexpr          /**< buffer to store pointer to expression that is product if first and second factor, or NULL */
   );

/** returns whether all expressions that are used in a quadratic expression are variable expression
 *
 * @return TRUE iff all linexprs and quadexprterms[.].expr in quaddata are variable expressions
 */
SCIP_EXPORT
SCIP_Bool SCIPareConsExprQuadraticExprsVariables(
   SCIP_QUADEXPR*       quaddata          /**< quadratic coefficients data */
   );

/** evaluates quadratic term in a solution w.r.t. auxiliary variables
 *
 * \note This assumes that for every expr used in the quadratic data, an auxiliary variable is available.
 */
SCIP_EXPORT
SCIP_Real SCIPevalConsExprQuadraticAux(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_QUADEXPR* quaddata,         /**< quadratic form */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL for LP solution */
   );

/** prints quadratic expression */
SCIP_EXPORT
SCIP_RETCODE SCIPprintConsExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression conshdlr */
   SCIP_QUADEXPR* quaddata          /**< quadratic form */
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
   SCIP_QUADEXPR* quaddata,         /**< quadratic coefficients data */
   SCIP_EXPRCURV*          curv,             /**< pointer to store the curvature of quadratics */
   SCIP_HASHMAP*           assumevarfixed,   /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   SCIP_Bool               storeeigeninfo    /**< whether the eigenvalues and eigenvectors should be stored */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* SCIP_PUB_EXPR_H_ */
