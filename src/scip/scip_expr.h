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
   SCIP_EXPRHDLR**       exprhdlr,     /**< buffer where to store created expression handler */
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

/** duplicates the given expression (including children) */
SCIP_EXPORT
SCIP_RETCODE SCIPcopyExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_DECL_EXPR_MAPVAR((*mapvar)),         /**< variable mapping function, or NULL for identity mapping */
   void*                 mapvardata,         /**< data of variable mapping function */
   SCIP_DECL_EXPR_MAPEXPR((*mapexpr)),       /**< expression mapping function, or NULL for creating new expressions */
   void*                 mapexprdata,        /**< data of expression mapping function */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call on expression copy to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata, /**< data to pass to ownerdatacreate */
   SCIP_DECL_EXPR_OWNERDATAFREE((*ownerdatafree)),     /**< function to call when freeing expression, e.g., to free ownerdata */
   );

/** duplicates the given expression without its children */
SCIP_EXPORT
SCIP_RETCODE SCIPcopyExprShallow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr            /**< buffer to store (shallow) duplicate of expr */
   );

/** creates an expression from a string
 *
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
 * See also @ref parseExpr in expr.c.
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

/** returns whether an expression is a sum expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a product expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a power expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprPower(
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

/** main part of printing an expression in dot format */
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
 * Then, traversing the tree in Depth First (see SCIPexpriterInit), for every expr that *has* children,
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

/** evaluates gradient of an expression for a given point
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffNonlinear().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalExprGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   SCIP_SOL*             sol,                /**< solution to be evaluated (NULL for the current LP solution) */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/** evaluates Hessian-vector product of an expression for a given point and direction
 *
 * Evaluates children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffGradientDirNonlinear().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalExprHessianDir(
   SCIP*                 scip,             /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   SCIP_SOL*             sol,              /**< solution to be evaluated (NULL for the current LP solution) */
   unsigned int          soltag,           /**< tag that uniquely identifies the solution (with its values), or 0. */
   SCIP_SOL*             direction         /**< direction */
   );

/**@} */  /* end of differentiation methods */

/** possibly reevaluates and then returns the activity of the expression
 *
 * Reevaluate activity if currently stored is no longer uptodate (some bound was changed since last evaluation).
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalExprActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** compute the hash value of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPhashExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   unsigned int*         hashval             /**< pointer to store the hash value */
   );

/** simplifies an expression
 *
 * The given expression will be released and overwritten with the simplified expression.
 * To keep the expression, duplicate it via SCIPcopyExpr before calling this method.
 *
 * This is largely inspired in Joel Cohen's
 * Computer algebra and symbolic computation: Mathematical methods
 * In particular Chapter 3
 * The other fountain of inspiration is the current simplifying methods in expr.c.
 *
 * Note: The things to keep in mind when adding simplification rules are the following.
 * I will be using the product expressions as an example.
 * There are mainly 3 parts of the simplification process. You need to decide
 * at which stage the simplification rule makes sense.
 * 1. Simplify each factor (simplifyFactor): At this stage we got the children of the product expression.
 * At this point, each child is simplified when viewed as a stand-alone
 * expression, but not necessarily when viewed as child of a product
 * expression. Rules like SP2, SP7, etc are enforced at this point.
 * 2. Multiply the factors (mergeProductExprlist): At this point rules like SP4, SP5 and SP14 are enforced.
 * 3. Build the actual simplified product expression (buildSimplifiedProduct):
 * At this point rules like SP10, SP11, etc are enforced.
 *
 * **During step 1. and 2. do not forget to set the flag changed to TRUE when something actually changes**
 *
 * Definition of simplified expressions
 * ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * An expression is simplified if it
 * - is a value expression
 * - is a var expression
 * - is a product expression such that
 *    SP1:  every child is simplified
 *    SP2:  no child is a product
 *    SP4:  no two children are the same expression (those should be multiplied)
 *    SP5:  the children are sorted [commutative rule]
 *    SP7:  no child is a value
 *    SP8:  its coefficient is 1.0 (otherwise should be written as sum)
 *    SP10: it has at least two children
 *    ? at most one child is an abs
 *    SP11: no two children are expr*log(expr)
 *    (TODO: we could handle more complicated stuff like x*y*log(x) -> - y * entropy(x), but I am not sure this should
 *    happen at the simplifcation level, or (x*y) * log(x*y), which currently simplifies to x * y * log(x*y))
 *    SP12: if it has two children, then neither of them is a sum (expand sums)
 *    SP13: no child is a sum with a single term
 *    SP14: at most one child is an exp
 * - is a (signed)power expression such that
 *   TODO: Some of these criteria are too restrictive for signed powers; for example, the exponent does not need to be
 *   an integer for signedpower to distribute over a product (POW5, POW6, POW8). Others can also be improved
 *    POW1: exponent is not 0
 *    POW2: exponent is not 1
 *    POW3: its child is not a value
 *    POW4: its child is simplified
 *    POW5: if exponent is integer, its child is not a product
 *    POW6: if exponent is integer, its child is not a sum with a single term ((2*x)^2 -> 4*x^2)
 *    POW7: if exponent is 2, its child is not a sum (expand sums)
 *    POW8: its child is not a power unless (x^n)^m with n*m being integer and n or m fractional and n not being even integer
 *    POW9: its child is not a sum with a single term with a positive coefficient: (25*x)^0.5 -> 5 x^0.5
 *    POW10: its child is not a binary variable: b^e and e > 0 --> b, b^e and e < 0 --> fix b to 1
 *    POW11: its child is not an exponential: exp(expr)^e --> exp(e * expr)
 * - is a signedpower expression such that
 *   TODO: Some of these criteria are too restrictive for signed powers; for example, the exponent does not need to be
 *   an integer for signedpower to distribute over a product (SPOW5, SPOW6, SPOW8). Others can also be improved
 *    SPOW1: exponent is not 0
 *    SPOW2: exponent is not 1
 *    SPOW3: its child is not a value
 *    SPOW4: its child is simplified
 *    SPOW5: (TODO) do we want to distribute signpowers over products like we do powers?
 *    SPOW6: exponent is not an odd integer: (signpow odd expr) -> (pow odd expr)
 *    SPOW8: if exponent is integer, its child is not a power
 *    SPOW9: its child is not a sum with a single term: (25*x)^0.5 -> 5 x^0.5
 *    SPOW10: its child is not a binary variable: b^e and e > 0 --> b, b^e and e < 0 --> fix b to 1
 *    SPOW11: its child is not an exponential: exp(expr)^e --> exp(e * expr)
 *    SPOW?: TODO: what happens when child is another signed power?
 *    SPOW?: if child >= 0 -> transform to normal power; if child < 0 -> transform to - normal power
 * - is a sum expression such that
 *    SS1: every child is simplified
 *    SS2: no child is a sum
 *    SS3: no child is a value (values should go in the constant of the sum)
 *    SS4: no two children are the same expression (those should be summed up)
 *    SS5: the children are sorted [commutative rule]
 *    SS6: it has at least one child
 *    SS7: if it consists of a single child, then either constant is != 0.0 or coef != 1
 *    SS8: no child has coefficient 0
 *    SS9: if a child c is a product that has an exponential expression as one of its factors, then the coefficient of c is +/-1.0
 *    SS10: if a child c is an exponential, then the coefficient of c is +/-1.0 (TODO)
 *    x if it consists of a single child, then its constant != 0.0 (otherwise, should be written as a product)
 * - it is a function with simplified arguments, but not all of them can be values
 * ? a logarithm doesn't have a product as a child
 * ? the exponent of an exponential is always 1
 *
 * ORDERING RULES (see SCIPexprCompare())
 * ^^^^^^^^^^^^^^
 * These rules define a total order on *simplified* expressions.
 * There are two groups of rules, when comparing equal type expressions and different type expressions
 * Equal type expressions:
 * OR1: u,v value expressions: u < v <=> val(u) < val(v)
 * OR2: u,v var expressions: u < v <=> SCIPvarGetIndex(var(u)) < SCIPvarGetIndex(var(v))
 * OR3: u,v are both sum or product expression: < is a lexicographical order on the terms
 * OR4: u,v are both pow: u < v <=> base(u) < base(v) or, base(u) == base(v) and expo(u) < expo(v)
 * OR5: u,v are u = FUN(u_1, ..., u_n), v = FUN(v_1, ..., v_m): u < v <=> For the first k such that u_k != v_k, u_k < v_k,
 *      or if such a k doesn't exist, then n < m.
 *
 * Different type expressions:
 * OR6: u value, v other: u < v always
 * OR7: u sum, v var or func: u < v <=> u < 0+v
 *      In other words, u = sum_{i = 1}^n alpha_i u_i, then u < v <=> u_n < v or if u_n = v and alpha_n < 1
 * OR8: u product, v pow, sum, var or func: u < v <=> u < 1*v
 *      In other words, u = Pi_{i = 1}^n u_i,  then u < v <=> u_n < v
 *      @note: since this applies only to simplified expressions, the form of the product is correct. Simplified products
 *             do *not* have constant coefficients
 * OR9: u pow, v sum, var or func: u < v <=> u < v^1
 * OR10: u var, v func: u < v always
 * OR11: u func, v other type of func: u < v <=> name(type(u)) < name(type(v))
 * OR12: none of the rules apply: u < v <=> ! v < u
 * Examples:
 * OR12: x < x^2 ?:  x is var and x^2 product, so none applies.
 *       Hence, we try to answer x^2 < x ?: x^2 < x <=> x < x or if x = x and 2 < 1 <=> 2 < 1 <=> False, so x < x^2 is True
 *       x < x^-1 --OR12--> ~(x^-1 < x) --OR9--> ~(x^-1 < x^1) --OR4--> ~(x < x or -1 < 1) --> ~True --> False
 *       x*y < x --OR8--> x*y < 1*x --OR3--> y < x --OR2--> False
 *       x*y < y --OR8--> x*y < 1*y --OR3--> y < x --OR2--> False
 *
 * Algorithm
 * ^^^^^^^^^
 * The recursive version of the algorithm is
 *
 * EXPR simplify(expr)
 *    for c in 1..expr->nchildren
 *       expr->children[c] = simplify(expr->children[c])
 *    end
 *    return expr->exprhdlr->simplify(expr)
 * end
 *
 * Important: Whatever is returned by a simplify callback **has** to be simplified.
 * Also, all children of the given expression **are** already simplified
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsimplifyExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            rootexpr,           /**< expression to be simplified */
   SCIP_EXPR**           simplified,         /**< buffer to store simplified expression */
   SCIP_Bool*            changed,            /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*            infeasible          /**< buffer to store whether infeasibility has been detected */
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

/** get the monotonicity of an expression w.r.t. to a given child */
SCIP_EXPORT
SCIP_RETCODE SCIPgetExprMonotonicity(
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


/**@name Expression Iterator Methods */
/**@{ */

/** creates an expression iterator */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExpriter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRITER**       iterator            /**< buffer to store expression iterator */
   );

/** frees an expression iterator */
SCIP_EXPORT
void SCIPfreeExpriter(
   SCIP_EXPRITER**       iterator            /**< pointer to the expression iterator */
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
