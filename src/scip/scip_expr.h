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

/**@file   scip_expr.h
 * @brief  public functions to work with algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_SCIP_EXPR_H_
#define SCIP_SCIP_EXPR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "scip/type_scip.h"
#include "scip/type_expr.h"

/**@name Expression Handler Methods */
/**@{ */

/** creates the handler for an expression handler and includes it into SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprHdlr(
   SCIP*                 scip,         /**< SCIP data structure */
   SCIP_EXPRHDLR**       exprhdlr,     /**< buffer where to store expression handler */
   const char*           name,         /**< name of expression handler (must not be NULL) */
   const char*           desc,         /**< description of expression handler (can be NULL) */
   unsigned int          precedence,   /**< precedence of expression operation (used for printing) */
   SCIP_DECL_EXPREVAL((*eval)),        /**< point evaluation callback (must not be NULL) */
   SCIP_EXPRHDLRDATA*    data          /**< data of expression handler (can be NULL) */
   );

/** gives expression handlers */
SCIP_EXPORT
SCIP_EXPRHDLR** SCIPgetExprHdlrs(
   SCIP*                      scip           /**< SCIP data structure */
);

/** gives number of expression handlers */
SCIP_EXPORT
int SCIPgetNExprHdlrs(
   SCIP*                      scip           /**< SCIP data structure */
);

/** returns an expression handler of a given name (or NULL if not found) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPfindExprHdlr(
   SCIP*                      scip,          /**< SCIP data structure */
   const char*                name           /**< name of expression handler */
   );

/** returns expression handler for variable expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprHdlrVar(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for constant value expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprHdlrValue(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for sum expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprHdlrSum(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for product expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprHdlrProduct(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for power expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprHdlrPower(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for signed power expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprHdlrSignPower(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for exponential expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprHdlrExponential(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for logarithm expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprHdlrLogarithm(
   SCIP*                      scip           /**< SCIP data structure */
);

/** calls the print callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRPRINT(SCIPcallExprhdlrPrint);

/** calls the parse callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRPARSE(SCIPcallExprhdlrParse);

/** calls the expression hash callback */
SCIP_EXPORT
SCIP_DECL_EXPRHASH(SCIPcallExprhdlrHash);

/** calls the expression compare callback */
SCIP_EXPORT
SCIP_DECL_EXPRCOMPARE(SCIPcallExprhdlrCompare);

/** calls the backward-differentiation callback of an expression handler
 *
 * further, allows to different w.r.t. given expression and children values
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcallExprhdlrBwdiff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   int                   childidx,           /**< index of child w.r.t. which to compute derivative */
   SCIP_Real*            derivative,         /**< buffer to store value of derivative */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_Real             exprval             /**< value for expression, used only if childrenvals is not NULL */
);

/** calls the backward-forward differentiation callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRBWFWDIFF(SCIPcallExprhdlrBwfwdiff);

/** calls the forward differentiation callback of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRFWDIFF(SCIPcallExprhdlrFwdiff);

/** calls the evaluation callback of an expression handler
 *
 * further, allows to evaluate w.r.t. given children values
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcallExprhdlrEval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            val,                /**< buffer store value of expression */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_SOL*             sol                 /**< solution that is evaluated (used by the var-expression) */
);

/** calls the expression interval evaluation callback */
SCIP_EXPORT
SCIP_DECL_EXPRINTEVAL(SCIPcallExprhdlrIntEval);

/** calls estimator method of expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRESTIMATE(SCIPcallExprhdlrEstimate);

/** calls the intitial estimators method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRINITESTIMATES(SCIPcallExprhdlrInitEstimates);

/** calls the simplification method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRSIMPLIFY(SCIPcallExprhdlrSimplify);

/** calls the curvature check method of an expression handler */
SCIP_EXPORT
SCIP_DECL_EXPRCURVATURE(SCIPcallExprhdlrCurvature);

/** calls the expression callback for reverse propagation */
SCIP_EXPORT
SCIP_DECL_EXPRREVERSEPROP(SCIPcallExprhdlrReverseProp);

/** @} */



/**@name Expression Methods */
/**@{ */

/** creates and captures an expression with given expression data and children */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data (expression assumes ownership) */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children            /**< children (can be NULL if nchildren is 0) */
   );

/** creates and captures an expression with given expression data and up to two children */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExpr2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data */
   SCIP_EXPR*            child1,             /**< first child (can be NULL) */
   SCIP_EXPR*            child2              /**< second child (can be NULL) */
   );

/** creates and captures an expression representing a quadratic function */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs           /**< array with coefficients of quadratic terms */
   );

/** creates and captures an expression representing a monomial */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprMonomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nfactors,           /**< number of factors in monomial */
   SCIP_VAR**            vars,               /**< variables in the monomial */
   SCIP_Real*            exponents           /**< exponent in each factor, or NULL if all 1.0 */
   );

/** appends child to the children list of expr */
SCIP_EXPORT
SCIP_RETCODE SCIPappendExprChild(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPR*            child               /**< expression to be appended */
   );

/** overwrites/replaces a child of an expressions
 *
 * @note the old child is released and the newchild is captured, unless they are the same (=same pointer)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPreplaceExprChild(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPR*              expr,             /**< expression which is going to replace a child */
   int                     childidx,         /**< index of child being replaced */
   SCIP_EXPR*              newchild          /**< the new child */
   );

/** remove all children of expr
 *
 * @attention only use if you really know what you are doing
 */
SCIP_EXPORT
SCIP_RETCODE SCIPremoveExprChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** duplicates the given expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcopyExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_Bool             copychildren        /**< whether children (and all successors) should be copied, too */
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
SCIP_RETCODE SCIPparseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer to store the expr parsed */
   const char*           exprstr,            /**< string with the expr to parse */
   const char**          finalpos            /**< buffer to store the position of exprstr where we finished reading, or NULL if not of interest */
   );

/** captures an expression (increments usage count) */
SCIP_EXPORT
void SCIPcaptureExpr(
   SCIP_EXPR*            expr                /**< expression to be captured */
   );

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_EXPORT
SCIP_RETCODE SCIPreleaseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr                /**< pointer to expression to be released */
   );

/** returns whether an expression is a variable expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a value expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** print an expression as info-message */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be printed */
   FILE*                 file                /**< file to print to, or NULL for stdout */
   );

/** initializes printing of expressions in dot format to a give FILE* pointer */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExprDotInit(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDATA**    printdata,        /**< buffer to store dot printing data */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_EXPRPRINT_WHAT     whattoprint       /**< info on what to print for each expression */
   );

/** initializes printing of expressions in dot format to a file with given filename */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExprDotInit2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDATA**    printdata,        /**< buffer to store dot printing data */
   const char*             filename,         /**< name of file to print to */
   SCIP_EXPRPRINT_WHAT     whattoprint       /**< info on what to print for each expression */
   );

SCIP_EXPORT
SCIP_RETCODE SCIPprintExprDot(
   SCIP*                  scip,              /**< SCIP data structure */
   SCIP_EXPRPRINTDATA*    printdata,         /**< data as initialized by \ref SCIPprintExprDotInit() */
   SCIP_EXPR*             expr               /**< expression to be printed */
   );

/** finishes printing of expressions in dot format */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExprDotFinal(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDATA**    printdata         /**< buffer where dot printing data has been stored */
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
SCIP_RETCODE SCIPshowExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression to be printed */
   );

/** prints structure of an expression a la Maple's dismantle */
SCIP_EXPORT
SCIP_RETCODE SCIPdismantleExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   SCIP_EXPR*            expr                /**< expression to dismantle */
   );

/** evaluate an expression in a point
 *
 * Iterates over expressions to also evaluate children, if necessary.
 * Value can be received via SCIPexprGetEvalValue().
 * If an evaluation error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 *
 * If a nonzero \p soltag is passed, then only (sub)expressions are
 * reevaluated that have a different solution tag. If a soltag of 0
 * is passed, then subexpressions are always reevaluated.
 * The tag is stored together with the value and can be received via
 * SCIPexprGetEvalTag().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
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
 * Value can be received via SCIPgetExprPartialDiffNonlinear().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeExprGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated (NULL for the current LP solution) */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/**@} */  /* end of differentiation methods */

/** possibly reevaluates and then returns the activity of the expression
 *
 * Reevaluate activity if currently stored is no longer uptodate (some bound was changed since last evaluation).
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalExprActivity(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPR*              expr,             /**< expression */
   SCIP_INTERVAL*          activity          /**< interval where to store expression */
   );

/** compute the hash value of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPhashExpr(
   SCIP*                 scip,             /**< SCIP data structure */
   SCIP_EXPR*            expr,             /**< expression */
   unsigned int*         hashval           /**< pointer to store the hash value */
   );

/** simplifies an expression
 *
 * The given expression will be released and overwritten with the simplified expression.
 * To keep the expression, duplicate it via SCIPcopyExpr before calling this method.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsimplifyExpr(
   SCIP*                 scip,             /**< SCIP data structure */
   SCIP_EXPR*            rootexpr,         /**< expression to be simplified */
   SCIP_EXPR**           simplified,       /**< buffer to store simplified expression */
   SCIP_Bool*            changed,          /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*            infeasible        /**< buffer to store whether infeasibility has been detected */
   );

/** computes the curvature of a given expression and all its subexpressions
 *
 *  @note this function also evaluates all subexpressions w.r.t. current variable bounds
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** computes the monotonicity of an expression w.r.t. to a given child */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeExprMonotonicity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   int                   childidx,           /**< index of child */
   SCIP_MONOTONE*        monotonicity        /**< buffer to store monotonicity */
   );

/** computes integrality information of a given expression and all its subexpressions
 *
 * the integrality information can be accessed via SCIPexprIsIntegral()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeExprIntegrality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the total number of variables in an expression
 *
 * The function counts variables in common sub-expressions only once.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetExprNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   int*                  nvars               /**< buffer to store the total number of variables */
   );

/** returns all variable expressions contained in a given expression
 *
 * the array to store all variable expressions needs
 * to be at least of size the number of unique variables in the expression which is given by SCIPgetExprNVars()
 * and can be bounded by SCIPgetNVars().
 *
 * @note function captures variable expressions
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetExprVarExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPR**           varexprs,           /**< array to store all variable expressions */
   int*                  nvarexprs           /**< buffer to store the total number of variable expressions */
   );

/** @} */


/**@name Expression iterator */
/**@{ */

/** gets the index an expression iterator can use to store iterator specific data in an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPactivateExprIterator(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  iterindex           /**< buffer to store iteration index */
   );

/** returns the index that an expression iterator used to store iterator specific data in an expression */
SCIP_EXPORT
void SCIPdeactivateExprIterator(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   iterindex           /**< iteration index that is not used anymore */
   );

/** get a new tag that can be used to mark an expression as visited */
SCIP_EXPORT
unsigned int SCIPgetExprIteratorNewVisitedTag(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @} */



/**@name Quadratic expression functions */
/**@{ */

/** checks whether an expression is quadratic
 *
 * An expression is quadratic if it is either a square (of some expression), a product (of two expressions),
 * or a sum of terms where at least one is a square or a product.
 *
 * Use \ref SCIPexprGetQuadraticData to get data about the representation as quadratic.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcheckExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool*            isquadratic         /**< buffer to store result */
   );

/** evaluates quadratic term in a solution
 *
 * \note This requires that every expr used in the quadratic data is a variable expression.
 */
SCIP_EXPORT
SCIP_Real SCIPevalExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL for LP solution */
   );

/** prints quadratic expression */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< quadratic expression */
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
SCIP_RETCODE SCIPcomputeExprQuadraticCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_EXPRCURV*        curv,               /**< pointer to store the curvature of quadratics */
   SCIP_HASHMAP*         assumevarfixed,     /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   SCIP_Bool             storeeigeninfo      /**< whether the eigenvalues and eigenvectors should be stored */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* SCIP_SCIP_EXPR_H_ */
