/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlpi/pub_expr.h
 * @brief  public methods for expressions, expression trees, expression graphs, and related stuff
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __NLPI_PUB_EXPR_H__
#define __NLPI_PUB_EXPR_H__

#include "scip/def.h"
#include "scip/pub_message.h"
#include "scip/intervalarith.h"
#include "blockmemshell/memory.h"
#include "nlpi/type_expr.h"
#include "nlpi/type_exprinterpret.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@name Expression curvature methods */
/**@{ */

/* gives curvature for a sum of two functions with given curvature */
extern
SCIP_EXPRCURV SCIPexprcurvAdd(
   SCIP_EXPRCURV         curv1,              /**< curvature of first summand */
   SCIP_EXPRCURV         curv2               /**< curvature of second summand */
   );

/** gives the curvature for the negation of a function with given curvature */
extern
SCIP_EXPRCURV SCIPexprcurvNegate(
   SCIP_EXPRCURV         curvature           /**< curvature of function */
   );

/* gives curvature for a functions with given curvature multiplied by a constant factor */
extern
SCIP_EXPRCURV SCIPexprcurvMultiply(
   SCIP_Real             factor,             /**< constant factor */
   SCIP_EXPRCURV         curvature           /**< curvature of other factor */
   );

/* gives curvature for base^exponent for given bounds and curvature of base-function and constant exponent */
extern
SCIP_EXPRCURV SCIPexprcurvPower(
   SCIP_INTERVAL         basebounds,         /**< bounds on base function */
   SCIP_EXPRCURV         basecurv,           /**< curvature of base function */
   SCIP_Real             exponent            /**< exponent */
   );

/* gives curvature for a monomial with given curvatures and bounds for each factor */
extern
SCIP_EXPRCURV SCIPexprcurvMonomial(
   int                   nfactors,           /**< number of factors in monomial */
   SCIP_Real*            exponents,          /**< exponents in monomial, or NULL if all 1.0 */
   int*                  factoridxs,         /**< indices of factors, or NULL if identity mapping */
   SCIP_EXPRCURV*        factorcurv,         /**< curvature of each factor */
   SCIP_INTERVAL*        factorbounds        /**< bounds of each factor */
   );

/** gives name as string for a curvature */
extern
const char* SCIPexprcurvGetName(
   SCIP_EXPRCURV         curv                /**< curvature */
   );

/**@} */

/**@name Expression operand methods */
/**@{ */

/** gives the name of an operand */
extern
const char* SCIPexpropGetName(
   SCIP_EXPROP           op                  /**< expression operand */
   );

/** gives the number of children of a simple operand
 * @return -1 for invalid operands and -2 for complex operands (those where the number of children depends on the expression)
 */
extern
int SCIPexpropGetNChildren(
   SCIP_EXPROP           op                  /**< expression operand */
   );

/**@} */

/**@name Expression methods */
/**@{ */

/** creates a simple expression */
extern
SCIP_RETCODE SCIPexprCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   SCIP_EXPROP           op,                 /**< operand of expression */
   ...                                       /**< arguments of operand */
   );

/** copies an expression including its children */
extern
SCIP_RETCODE SCIPexprCopyDeep(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           targetexpr,         /**< buffer to store pointer to copied expression */
   SCIP_EXPR*            sourceexpr          /**< expression to copy */
   );

/** frees an expression including its children */
extern
void SCIPexprFreeDeep(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr                /**< pointer to expression to free */
   );

/** gives operator of expression */
extern
SCIP_EXPROP SCIPexprGetOperator(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives number of children of an expression */
extern
int SCIPexprGetNChildren(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives pointer to array with children of an expression */
extern
SCIP_EXPR** SCIPexprGetChildren(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives index belonging to a SCIP_EXPR_VARIDX or SCIP_EXPR_PARAM operand */
extern
int SCIPexprGetOpIndex(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives real belonging to a SCIP_EXPR_CONST operand */ 
extern
SCIP_Real SCIPexprGetOpReal(
   SCIP_EXPR* expr                           /**< expression */
   );

/** gives void* belonging to a complex operand */
extern
void* SCIPexprGetOpData(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives exponent belonging to a SCIP_EXPR_REALPOWER expression */
extern
SCIP_Real SCIPexprGetRealPowerExponent(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives exponent belonging to a SCIP_EXPR_INTPOWER expression */
extern
int SCIPexprGetIntPowerExponent(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives exponent belonging to a SCIP_EXPR_SIGNPOWER expression */
extern
SCIP_Real SCIPexprGetSignPowerExponent(
   SCIP_EXPR*            expr                /**< expression */
   );

/** creates a SCIP_EXPR_LINEAR expression that is (affine) linear in its children: constant + sum_i coef_i child_i */
extern
SCIP_RETCODE SCIPexprCreateLinear(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children of expression */
   SCIP_Real*            coefs,              /**< coefficients of children */
   SCIP_Real             constant            /**< constant part */
   );

/** gives linear coefficients belonging to a SCIP_EXPR_LINEAR expression */
extern
SCIP_Real* SCIPexprGetLinearCoefs(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives constant belonging to a SCIP_EXPR_LINEAR expression  */
extern
SCIP_Real SCIPexprGetLinearConstant(
   SCIP_EXPR*            expr                /**< expression */
   );

/** creates a SCIP_EXPR_QUADRATIC expression: constant + sum_i coef_i child_i + sum_i coef_i child1_i child2_i */
SCIP_RETCODE SCIPexprCreateQuadratic(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children of expression */
   SCIP_Real             constant,           /**< constant */
   SCIP_Real*            lincoefs,           /**< linear coefficients of children, or NULL if all 0.0 */
   int                   nquadelems,         /**< number of quadratic elements */
   SCIP_QUADELEM*        quadelems           /**< quadratic elements specifying coefficients and child indices */
   );

/** gives constant belonging to a SCIP_EXPR_QUADRATIC expression */
extern
SCIP_Real SCIPexprGetQuadConstant(
   SCIP_EXPR*            expr                /**< quadratic expression */
   );

/** gives linear coefficients belonging to a SCIP_EXPR_QUADRATIC expression
 * can be NULL if all coefficients are 0.0 */
extern
SCIP_Real* SCIPexprGetQuadLinearCoefs(
   SCIP_EXPR*            expr                /**< quadratic expression */
   );

/** gives quadratic elements belonging to a SCIP_EXPR_QUADRATIC expression */
extern
SCIP_QUADELEM* SCIPexprGetQuadElements(
   SCIP_EXPR*            expr                /**< quadratic expression */
   );

/** gives number of quadratic elements belonging to a SCIP_EXPR_QUADRATIC expression */
extern
int SCIPexprGetNQuadElements(
   SCIP_EXPR*            expr                /**< quadratic expression */
   );

/** ensures that quadratic elements of a quadratic expression are sorted */
extern
void SCIPexprSortQuadElems(
   SCIP_EXPR*            expr                /**< quadratic expression */
   );

/** creates a SCIP_EXPR_POLYNOMIAL expression from an array of monomials: constant + sum_i monomial_i */
extern
SCIP_RETCODE SCIPexprCreatePolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children of expression */
   int                   nmonomials,         /**< number of monomials */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< monomials */
   SCIP_Real             constant,           /**< constant part */
   SCIP_Bool             copymonomials       /**< should monomials by copied or ownership be assumed? */
   );

/** gives the monomials belonging to a SCIP_EXPR_POLYNOMIAL expression */
extern
SCIP_EXPRDATA_MONOMIAL** SCIPexprGetMonomials(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the number of monomials belonging to a SCIP_EXPR_POLYNOMIAL expression */
extern
int SCIPexprGetNMonomials(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the constant belonging to a SCIP_EXPR_POLYNOMIAL expression */
extern
SCIP_Real SCIPexprGetPolynomialConstant(
   SCIP_EXPR*            expr                /**< expression */
   );

/** adds an array of monomials to a SCIP_EXPR_POLYNOMIAL expression */
extern
SCIP_RETCODE SCIPexprAddMonomials(
   BMS_BLKMEM*           blkmem,             /**< block memory of expression */
   SCIP_EXPR*            expr,               /**< expression */
   int                   nmonomials,         /**< number of monomials to add */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< the monomials to add */
   SCIP_Bool             copymonomials       /**< should monomials by copied or ownership be assumed? */
   );

/** changes the constant in a SCIP_EXPR_POLYNOMIAL expression */
extern
void SCIPexprChgPolynomialConstant(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             constant            /**< new value for constant */
   );

/** multiplies each summand of a polynomial by a given constant */
extern
void SCIPexprMultiplyPolynomialByConstant(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_Real             factor              /**< constant factor */
   );

/** multiplies each summand of a polynomial by a given monomial */
extern
SCIP_RETCODE SCIPexprMultiplyPolynomialByMonomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_EXPRDATA_MONOMIAL* factor,           /**< monomial factor */
   int*                  childmap            /**< map children in factor to children in expr, or NULL for 1:1 */
   );

/** multiplies this polynomial by a polynomial
 * factor needs to be different from expr */
extern
SCIP_RETCODE SCIPexprMultiplyPolynomialByPolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_EXPR*            factor,             /**< polynomial factor */
   int*                  childmap            /**< map children in factor to children in expr, or NULL for 1:1 */
   );

/** takes a power of the polynomial
 * exponent need to be an integer
 * polynomial need to be a monomial, if exponent is negative
 */
extern
SCIP_RETCODE SCIPexprPolynomialPower(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   int                   exponent            /**< exponent of power operation */
   );

/** merges monomials in a polynomial expression that differ only in coefficient into a single monomial
 * eliminates monomials with coefficient between -eps and eps
 */
extern
void SCIPexprMergeMonomials(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_Real             eps,                /**< threshold under which numbers are treat as zero */
   SCIP_Bool             mergefactors        /**< whether to merge factors in monomials too */
   );

/** creates a monomial */
extern
SCIP_RETCODE SCIPexprCreateMonomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_MONOMIAL** monomial,        /**< buffer where to store pointer to new monomial */
   SCIP_Real             coef,               /**< coefficient of monomial */
   int                   nfactors,           /**< number of factors in monomial */
   int*                  childidxs,          /**< indices of children corresponding to factors, or NULL if identity */
   SCIP_Real*            exponents           /**< exponent in each factor, or NULL if all 1.0 */
   );

/** frees a monomial */
extern
void SCIPexprFreeMonomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_MONOMIAL** monomial         /**< pointer to monomial that should be freed */
   );

/** gets coefficient of a monomial */
extern
SCIP_Real SCIPexprGetMonomialCoef(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   );

/** gets number of factors of a monomial */
extern
int SCIPexprGetMonomialNFactors(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   );

/** gets indices of children corresponding to factors of a monomial */
extern
int* SCIPexprGetMonomialChildIndices(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   );

/** gets exponents in factors of a monomial */
extern
SCIP_Real* SCIPexprGetMonomialExponents(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   );

/** ensures that factors in a monomial are sorted */
extern
void SCIPexprSortMonomialFactors(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   );

/** finds a factor corresponding to a given child index in a monomial
 * note that if the factors have not been merged, the position of some factor corresponding to a given child is given
 * returns TRUE if a factor is found, FALSE if not
 */
extern
SCIP_Bool SCIPexprFindMonomialFactor(
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   int                   childidx,           /**< index of the child which factor to search for */
   int*                  pos                 /**< buffer to store position of factor */
   );

/** checks if two monomials are equal */
extern
SCIP_Bool SCIPexprAreMonomialsEqual(
   SCIP_EXPRDATA_MONOMIAL*  monomial1,       /**< first monomial */
   SCIP_EXPRDATA_MONOMIAL*  monomial2,       /**< second monomial */
   SCIP_Real             eps                 /**< threshold under which numbers are treated as 0.0 */
   );

/** adds factors to a monomial */
extern
SCIP_RETCODE SCIPexprAddMonomialFactors(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   int                   nfactors,           /**< number of factors to add */
   int*                  childidxs,          /**< indices of children corresponding to factors */
   SCIP_Real*            exponents           /**< exponent in each factor */
   );

/** changes coefficient of monomial */
extern
void SCIPexprChgMonomialCoef(
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   SCIP_Real             newcoef             /**< new coefficient */
   );

/** multiplies a monomial with a monomial */
extern
SCIP_RETCODE SCIPexprMultiplyMonomialByMonomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   SCIP_EXPRDATA_MONOMIAL* factor,           /**< factor monomial */
   int*                  childmap            /**< map to apply to children in factor, or NULL for 1:1 */
   );

/** replaces the monomial by a power of the monomial
 * allows only integers as exponent
 */
extern
void SCIPexprMonomialPower(
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   int                   exponent            /**< integer exponent of power operation */
   );

/** merges factors that correspond to the same child by adding exponents
 * eliminates factors with exponent between -eps and eps
 */
extern
void SCIPexprMergeMonomialFactors(
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   SCIP_Real             eps                 /**< threshold under which numbers are treated as 0.0 */
   );

/** ensures that monomials of a polynomial are sorted */
extern
void SCIPexprSortMonomials(
   SCIP_EXPR*            expr                /**< polynomial expression */
   );

/** indicates whether the expression contains a SCIP_EXPR_PARAM */
extern
SCIP_Bool SCIPexprHasParam(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gets maximal degree of expression, or SCIP_EXPR_DEGREEINFINITY if not a polynomial */
extern
SCIP_RETCODE SCIPexprGetMaxDegree(
   SCIP_EXPR*            expr,               /**< expression */
   int*                  maxdegree           /**< buffer to store maximal degree */
   );

/** counts usage of variables in expression */
extern
void SCIPexprGetVarsUsage(
   SCIP_EXPR*            expr,               /**< expression to update */
   int*                  varsusage           /**< array with counters of variable usage */
   );

/** compares whether two expressions are the same
 * inconclusive, i.e., may give FALSE even if expressions are equivalent (x*y != y*x)
 */
extern
SCIP_Bool SCIPexprAreEqual(
   SCIP_EXPR*            expr1,              /**< first expression */
   SCIP_EXPR*            expr2,              /**< second expression */
   SCIP_Real             eps                 /**< threshold under which numbers are assumed to be zero */
   );

/** aims at simplifying an expression and splitting of a linear expression
 * if linear variables are split off, expression interpreter data, if stored in the tree, is freed
 */
extern
SCIP_RETCODE SCIPexprSimplify(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             eps,                /**< threshold, under which positive values are treat as 0 */
   int                   maxexpansionexponent,/**< maximal exponent for which we still expand non-monomial polynomials */
   int                   nvars,              /**< number of variables in expression */
   int*                  nlinvars,           /**< buffer to store number of linear variables in linear part, or NULL if linear part should not be separated */
   int*                  linidxs,            /**< array to store indices of variables in expression tree which belong to linear part, or NULL */
   SCIP_Real*            lincoefs            /**< array to store coefficients of linear part, or NULL */
   );

/** evaluates an expression w.r.t. a point */
extern
SCIP_RETCODE SCIPexprEval(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            varvals,            /**< values for variables, can be NULL if the expression is constant */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression is not parameterized */
   SCIP_Real*            val                 /**< buffer to store value */
   );

/** evaluates an expression w.r.t. an interval */
extern
SCIP_RETCODE SCIPexprEvalInt(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values for variables, can be NULL if the expression is constant */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression is not parameterized */
   SCIP_INTERVAL*        val                 /**< buffer to store value */
   );

/** tries to determine the curvature type of an expression w.r.t. given variable domains */
extern
SCIP_RETCODE SCIPexprCheckCurvature(
   SCIP_EXPR*            expr,               /**< expression to check */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        varbounds,          /**< domains of variables */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression is not parameterized */
   SCIP_EXPRCURV*        curv,               /**< buffer to store curvature of expression */
   SCIP_INTERVAL*        bounds              /**< buffer to store bounds on expression */
   );

/** substitutes variables (SCIP_EXPR_VARIDX) by expressions
 * Note that only the children of the given expr are checked!
 * A variable with index i is replaced by a copy of substexprs[i], if that latter is not NULL
 * if substexprs[i] == NULL, then the variable expression i is not touched
 */
extern
SCIP_RETCODE SCIPexprSubstituteVars(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR*            expr,               /**< expression, which of the children may be replaced */
   SCIP_EXPR**           substexprs          /**< array of substitute expressions; single entries can be NULL */
   );

/** updates variable indices in expression tree */
extern
void SCIPexprReindexVars(
   SCIP_EXPR*            expr,               /**< expression to update */
   int*                  newindices          /**< new indices of variables */
   );

/** updates parameter indices in expression tree */
extern
void SCIPexprReindexParams(
   SCIP_EXPR*            expr,               /**< expression to update */
   int*                  newindices          /**< new indices of variables */
   );

/** prints an expression */
extern
void SCIPexprPrint(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file for printing, or NULL for stdout */
   const char**          varnames,           /**< names of variables, or NULL for default names */
   const char**          paramnames,         /**< names of parameters, or NULL for default names */
   SCIP_Real*            paramvals           /**< values of parameters, or NULL for not printing */
   );


/**@} */

/**@name Expression tree methods */
/**@{ */

/** creates an expression tree */
extern
SCIP_RETCODE SCIPexprtreeCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRTREE**       tree,               /**< buffer to store address of created expression tree */
   SCIP_EXPR*            root,               /**< pointer to root expression, not copied deep !, can be NULL */
   int                   nvars,              /**< number of variables in variable mapping */
   int                   nparams,            /**< number of parameters in expression */
   SCIP_Real*            params              /**< values for parameters, or NULL (if NULL but nparams > 0, then params is initialized with zeros) */
   );

/** copies an expression tree */
extern
SCIP_RETCODE SCIPexprtreeCopy(
   BMS_BLKMEM*           blkmem,             /**< block memory that should be used in new expression tree */
   SCIP_EXPRTREE**       targettree,         /**< buffer to store address of copied expression tree */
   SCIP_EXPRTREE*        sourcetree          /**< expression tree to copy */
   );

/** frees an expression tree */
extern
SCIP_RETCODE SCIPexprtreeFree(
   SCIP_EXPRTREE**       tree                /**< pointer to expression tree that is freed */
   );

/** returns root expression of an expression tree */
extern
SCIP_EXPR* SCIPexprtreeGetRoot(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** returns number of variables in expression tree */
extern
int SCIPexprtreeGetNVars(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** returns number of parameters in expression tree */
extern
int SCIPexprtreeGetNParams(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** returns values of parameters or NULL if none */
extern
SCIP_Real* SCIPexprtreeGetParamVals(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** sets value of a single parameter in expression tree */
extern
void SCIPexprtreeSetParamVal(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   paramidx,           /**< index of parameter */
   SCIP_Real             paramval            /**< new value of parameter */
   );

/** sets number and values of all parameters in expression tree */
extern
SCIP_RETCODE SCIPexprtreeSetParams(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nparams,            /**< number of parameters */
   SCIP_Real*            paramvals           /**< values of parameters, can be NULL if nparams == 0 */
   );

/** gets data of expression tree interpreter, or NULL if not set */
extern
SCIP_EXPRINTDATA* SCIPexprtreeGetInterpreterData(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** sets data of expression tree interpreter */
extern
void SCIPexprtreeSetInterpreterData(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_EXPRINTDATA*     interpreterdata     /**< expression interpreter data */
   );

/** frees data of expression tree interpreter, if any */
extern
SCIP_RETCODE SCIPexprtreeFreeInterpreterData(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** indicates whether there are parameterized constants (SCIP_EXPR_PARAM) in expression tree */
extern
SCIP_Bool SCIPexprtreeHasParam(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** Gives maximal degree of expression in expression tree.
 * If constant expression, gives 0,
 * if linear expression, gives 1,
 * if polynomial expression, gives its maximal degree,
 * otherwise (nonpolynomial nonconstant expressions) gives at least SCIP_EXPR_DEGREEINFINITY.
 */
extern
SCIP_RETCODE SCIPexprtreeGetMaxDegree(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int*                  maxdegree           /**< buffer to store maximal degree */
   );

/** gives the number of usages for each variable in the expression tree */
extern
void SCIPexprtreeGetVarsUsage(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int*                  varsusage           /**< array where to store for each variable how often it is used in the tree */
   );

/** aims at simplifying an expression and splitting of a linear expression
 * if linear variables are split off, expression interpreter data, if stored in the tree, is freed
 */
extern
SCIP_RETCODE SCIPexprtreeSimplify(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Real             eps,                /**< threshold, under which positive values are treat as 0 */
   int                   maxexpansionexponent,/**< maximal exponent for which we still expand non-monomial polynomials */
   int*                  nlinvars,           /**< buffer to store number of linear variables in linear part, or NULL if linear part should not be separated */
   int*                  linidxs,            /**< array to store indices of variables in expression tree which belong to linear part, or NULL */
   SCIP_Real*            lincoefs            /**< array to store coefficients of linear part, or NULL */
   );

/** adds an expression to the root expression of the tree
 * the root is replaced with an SCIP_EXPR_PLUS expression which has the previous root and the given expression as children
 */
extern
SCIP_RETCODE SCIPexprtreeAddExpr(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_EXPR*            expr,               /**< expression to add to tree */
   SCIP_Bool             copyexpr            /**< whether expression should be copied */
   );

/** evaluates an expression tree w.r.t. a point */
extern
SCIP_RETCODE SCIPexprtreeEval(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values for variables */
   SCIP_Real*            val                 /**< buffer to store expression tree value */
   );

/** evaluates an expression tree w.r.t. an interval */
extern
SCIP_RETCODE SCIPexprtreeEvalInt(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< intervals for variables */
   SCIP_INTERVAL*        val                 /**< buffer to store expression tree value */
   );

/** tries to determine the curvature type of an expression tree w.r.t. given variable domains */
extern
SCIP_RETCODE SCIPexprtreeCheckCurvature(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varbounds,          /**< domains of variables */
   SCIP_EXPRCURV*        curv,               /**< buffer to store curvature of expression */
   SCIP_INTERVAL*        bounds              /**< buffer to store bounds on expression, or NULL if not needed */
   );

/** substitutes variables (SCIP_EXPR_VARIDX) in an expression tree by expressions
 * A variable with index i is replaced by a copy of substexprs[i], if that latter is not NULL
 * if substexprs[i] == NULL, then the variable expression i is not touched
 */
extern
SCIP_RETCODE SCIPexprtreeSubstituteVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_EXPR**           substexprs          /**< array of substitute expressions; single entries can be NULL */
   );

/** prints an expression tree */
extern
void SCIPexprtreePrint(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file for printing, or NULL for stdout */
   const char**          varnames,           /**< names of variables, or NULL for default names */
   const char**          paramnames          /**< names of parameters, or NULL for default names */
   );

/**@} */

/**@name Quadratic element methods */
/**@{ */

/** sorts an array of quadratic elements
 * The elements are sorted such that the first index is increasing and
 * such that among elements with the same first index, the second index is increasing.
 * For elements with same first and second index, the order is not defined.
 */
extern
void SCIPquadelemSort(
   SCIP_QUADELEM*        quadelems,          /**< array of quadratic elements */
   int                   nquadelems          /**< number of quadratic elements */
   );

/** Finds an index pair in a sorted array of quadratic elements.
 * If (idx1,idx2) is found in quadelems, then returns TRUE and stores position of quadratic element in *pos.
 * If (idx1,idx2) is not found in quadelems, then returns FALSE and stores position where a quadratic element with these indices would be inserted in *pos.
 * Assumes that idx1 <= idx2.
 */
extern
SCIP_Bool SCIPquadelemSortedFind(
   SCIP_QUADELEM*        quadelems,          /**< array of quadratic elements */
   int                   idx1,               /**< index of first  variable in element to search for */
   int                   idx2,               /**< index of second variable in element to search for */
   int                   nquadelems,         /**< number of quadratic elements in array */
   int*                  pos                 /**< buffer to store position of found quadratic element, or position where it would be inserted */
   );

/** Adds quadratic elements with same index and removes elements with coefficient 0.0.
 * Assumes that elements have been sorted before.
 */
extern
void SCIPquadelemSqueeze(
   SCIP_QUADELEM*        quadelems,          /**< array of quadratic elements */
   int                   nquadelems,         /**< number of quadratic elements */
   int*                  nquadelemsnew       /**< pointer to store new (reduced) number of quadratic elements */
   );

/**@} */

/**@name Expression graph node methods */
/**@{ */

/** creates an expression graph node */
extern
SCIP_RETCODE SCIPexprgraphCreateNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   SCIP_EXPROP           op,                 /**< operator type of expression */
   ...
   );

/** creates an expression graph node for a linear expression */
extern
SCIP_RETCODE SCIPexprgraphCreateNodeLinear(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   int                   ncoefs,             /**< number of coefficients */
   SCIP_Real*            coefs,              /**< coefficients of linear expression */
   SCIP_Real             constant            /**< constant of linear expression */
   );

/** creates an expression graph node for a quadratic expression */
extern
SCIP_RETCODE SCIPexprgraphCreateNodeQuadratic(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   int                   nchildren,          /**< number of children */
   SCIP_Real*            lincoefs,           /**< linear coefficients for children, or NULL */
   int                   nquadelems,         /**< number of quadratic elements */
   SCIP_QUADELEM*        quadelems,          /**< quadratic elements, or NULL if nquadelems == 0 */
   SCIP_Real             constant            /**< constant */
   );

/** creates an expression graph node for a polynomial expression */
extern
SCIP_RETCODE SCIPexprgraphCreateNodePolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   int                   nmonomials,         /**< number of monomials */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< monomials */
   SCIP_Real             constant,           /**< constant of polynomial */
   SCIP_Bool             copymonomials       /**< whether to copy monomials or to assume ownership */
   );

/** adds monomials to an expression graph node that is a polynomial expression */
extern
SCIP_RETCODE SCIPexprgraphNodePolynomialAddMonomials(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE*   node,               /**< store expression graph node with polynomial operator */
   int                   nmonomials,         /**< number of monomials */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< monomials */
   SCIP_Bool             copymonomials       /**< whether to copy monomials or to assume ownership */
   );

/** given a node of an expression graph, splitup a linear part which variables are not used somewhere else in the same expression
 * E.g., if the expression is 1 + x + y + y^2, one gets 1 + x and the node remains at y + y^2.
 * If the node is a linear expression, it may be freed.
 * If it is not linear, the node may change, i.e., the remaining nonlinear part may be stored in a new node.
 * It is assumed that the user had captured the node.
 * It is assumed that the expression graph has been simplified before.
 */
extern
SCIP_RETCODE SCIPexprgraphNodeSplitOffLinear(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE**  node,               /**< expression graph node where to splitup linear part */
   int                   linvarssize,        /**< length of linvars and lincoefs arrays */
   int*                  nlinvars,           /**< buffer to store length of linear term that have been splitup */
   void**                linvars,            /**< buffer to store variables of linear part */
   SCIP_Real*            lincoefs,           /**< buffer to store coefficients of linear part */
   SCIP_Real*            constant            /**< buffer to store constant part */
   );

/** moves parents from a one node to another node
 * in other words, replaces the child srcnode by targetnode in all parents of srcnode
 * srcnode may be freed, if not captured
 * it is assumes that targetnode represents the same expression as srcnode
 */
extern
SCIP_RETCODE SCIPexprgraphMoveNodeParents(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE**  srcnode,            /**< node which parents to move */
   SCIP_EXPRGRAPHNODE*   targetnode          /**< node where to move parents to */
   );

/** captures node, i.e., increases number of uses */
extern
void SCIPexprgraphCaptureNode(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node to capture */
   );

/** releases node, i.e., decreases number of uses
 * node is freed if no parents and no other uses
 * children are recursively released if they have no other parents
 * nodes that are removed are also freed
 * if node correspond to a variable, then the variable is removed from the expression graph
 * similar for constants
 */
extern
SCIP_RETCODE SCIPexprgraphReleaseNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE**  node                /**< expression graph node to release */
   );

/** frees a node of an expression graph */
extern
void SCIPexprgraphFreeNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node                /**< pointer to expression graph node that should be freed */
   );

/** enables a node and recursively all its children in an expression graph */
extern
void SCIPexprgraphEnableNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node to enable */
   );

/** disables a node and recursively all children which have no enabled parents in an expression graph */
extern
void SCIPexprgraphDisableNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node to enable */
   );

/** returns whether a node is currently enabled */
extern
SCIP_Bool SCIPexprgraphIsNodeEnabled(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node to enable */
   );

/** gets number of children of a node in an expression graph */
extern
int SCIPexprgraphGetNodeNChildren(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gets children of a node in an expression graph */
extern
SCIP_EXPRGRAPHNODE** SCIPexprgraphGetNodeChildren(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gets number of parents of a node in an expression graph */
extern
int SCIPexprgraphGetNodeNParents(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gets parents of a node in an expression graph */
extern
SCIP_EXPRGRAPHNODE** SCIPexprgraphGetNodeParents(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** returns whether the node has siblings in the expression graph */
extern
SCIP_Bool SCIPexprgraphHasNodeSibling(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gets depth of node in expression graph */
extern
int SCIPexprgraphGetNodeDepth(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gets position of node in expression graph at its depth level */
extern
int SCIPexprgraphGetNodePosition(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gets operator of a node in an expression graph */
extern
SCIP_EXPROP SCIPexprgraphGetNodeOperator(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives index belonging to a SCIP_EXPR_VARIDX or SCIP_EXPR_PARAM operand */
extern
int SCIPexprgraphGetNodeOperatorIndex(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives real belonging to a SCIP_EXPR_CONST operand */
extern
SCIP_Real SCIPexprgraphGetNodeOperatorReal(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives variable belonging to a SCIP_EXPR_VARIDX expression */
extern
void* SCIPexprgraphGetNodeVar(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives exponent belonging to a SCIP_EXPR_REALPOWER expression */
extern
SCIP_Real SCIPexprgraphGetNodeRealPowerExponent(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives exponent belonging to a SCIP_EXPR_INTPOWER expression */
extern
int SCIPexprgraphGetNodeIntPowerExponent(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives exponent belonging to a SCIP_EXPR_SIGNPOWER expression */
extern
SCIP_Real SCIPexprgraphGetNodeSignPowerExponent(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives linear coefficients belonging to a SCIP_EXPR_LINEAR expression */
extern
SCIP_Real* SCIPexprgraphGetNodeLinearCoefs(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives constant belonging to a SCIP_EXPR_LINEAR expression  */
extern
SCIP_Real SCIPexprgraphGetNodeLinearConstant(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives constant belonging to a SCIP_EXPR_QUADRATIC expression */
extern
SCIP_Real SCIPexprgraphGetNodeQuadraticConstant(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives linear coefficients belonging to a SCIP_EXPR_QUADRATIC expression, or NULL if all coefficients are 0.0 */
extern
SCIP_Real* SCIPexprgraphGetNodeQuadraticLinearCoefs(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives quadratic elements belonging to a SCIP_EXPR_QUADRATIC expression */
extern
SCIP_QUADELEM* SCIPexprgraphGetNodeQuadraticQuadElements(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives number of quadratic elements belonging to a SCIP_EXPR_QUADRATIC expression */
extern
int SCIPexprgraphGetNodeQuadraticNQuadElements(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives the monomials belonging to a SCIP_EXPR_POLYNOMIAL expression */
extern
SCIP_EXPRDATA_MONOMIAL** SCIPexprgraphGetNodePolynomialMonomials(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives the number of monomials belonging to a SCIP_EXPR_POLYNOMIAL expression */
extern
int SCIPexprgraphGetNodePolynomialNMonomials(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gives the constant belonging to a SCIP_EXPR_POLYNOMIAL expression */
extern
SCIP_Real SCIPexprgraphGetNodePolynomialConstant(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gets bounds of a node in an expression graph */
extern
SCIP_INTERVAL SCIPexprgraphGetNodeBounds(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gets value of expression associated to node from last evaluation call */
extern
SCIP_Real SCIPexprgraphGetNodeVal(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** gets curvature of expression associated to node from last curvature check call */
extern
SCIP_EXPRCURV SCIPexprgraphGetNodeCurvature(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** returns whether all children of an expression graph node are variable nodes
 * gives TRUE for nodes without children
 */
extern
SCIP_Bool SCIPexprgraphAreAllNodeChildrenVars(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** returns whether the node has an ancestor which has a nonlinear expression operand */
extern
SCIP_Bool SCIPexprgraphHasNodeNonlinearAncestor(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** prints an expression graph node */
extern
void SCIPexprgraphPrintNode(
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file to print to, or NULL for stdout */
   );

/** tightens the bounds in a node of the graph
 * preparation for reverse propagation
 * sets bound status to SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTRECENT if tightening is strong enough and not cutoff
 */
extern
void SCIPexprgraphTightenNodeBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< node in expression graph with no parents */
   SCIP_INTERVAL         nodebounds,         /**< new bounds for node */
   SCIP_Real             minstrength,        /**< minimal required relative bound strengthening in a node to trigger a propagation into children nodes */
   SCIP_Bool*            cutoff              /**< buffer to store whether a node's bounds were propagated to an empty interval */
   );

/** ensures that bounds and curvature information in a node is uptodate
 * assumes that bounds and curvature in children are uptodate
 */
extern
SCIP_RETCODE SCIPexprgraphUpdateNodeBoundsCurvature(
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Real             minstrength,        /**< minimal required relative bound strengthening to trigger a bound recalculation in parent nodes */
   SCIP_Bool             clearreverseprop    /**< whether to reset bound tightenings from reverse propagation */
   );

/**@} */

/**@name Expression graph methods */
/**@{ */

/** creates an empty expression graph */
extern
SCIP_RETCODE SCIPexprgraphCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPH**      exprgraph,          /**< buffer to store pointer to expression graph */
   int                   varssizeinit,       /**< minimal initial size for variables array, or -1 to choose automatically */
   int                   depthinit,          /**< minimal initial depth of expression graph, or -1 to choose automatically */
   SCIP_DECL_EXPRGRAPHVARADDED((*exprgraphvaradded)), /** callback method to invoke when a variable has been added to the expression graph, or NULL if not needed */
   SCIP_DECL_EXPRGRAPHVARREMOVE((*exprgraphvarremove)), /** callback method to invoke when a variable will be removed from the expression graph, or NULL if not needed */
   SCIP_DECL_EXPRGRAPHVARCHGIDX((*exprgraphvarchgidx)), /** callback method to invoke when a variable changes its index in the expression graph, or NULL if not needed */
   void*                 userdata            /**< user data to pass to callback functions */
   );

/** frees an expression graph */
extern
SCIP_RETCODE SCIPexprgraphFree(
   SCIP_EXPRGRAPH**      exprgraph           /**< pointer to expression graph that should be freed */
   );

/** adds an expression graph node to an expression graph
 * expression graph assumes ownership of node
 * children are notified about new parent
 * depth will be chosen to be the maximum of mindepth and the depth of all children plus one
 */
extern
SCIP_RETCODE SCIPexprgraphAddNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node to add */
   int                   mindepth,           /**< minimal depth in expression graph where to add node, e.g., 0 or smaller to choose automatically */
   int                   nchildren,          /**< number of children */
   SCIP_EXPRGRAPHNODE**  children            /**< children nodes, or NULL if no children */
   );

/** adds variables to an expression graph, if not existing yet
 * also already existing nodes are enabled
 */
extern
SCIP_RETCODE SCIPexprgraphAddVars(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   int                   nvars,              /**< number of variables to add */
   void**                vars,               /**< variables to add */
   SCIP_EXPRGRAPHNODE**  varnodes            /**< array to store nodes corresponding to variables, or NULL if not of interest */
   );

/** adds a constant to an expression graph, if not existing yet
 * also already existing nodes are enabled
 */
SCIP_RETCODE SCIPexprgraphAddConst(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             constant,           /**< constant to add */
   SCIP_EXPRGRAPHNODE**  constnode           /**< buffer to store pointer to expression graph node corresponding to constant */
   );

/** adds sum of expression trees into expression graph
 * node will also be captured
 */
extern
SCIP_RETCODE SCIPexprgraphAddExprtreeSum(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   int                   nexprtrees,         /**< number of expression trees to add */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees that should be added */
   SCIP_Real*            coefs,              /**< coefficients of expression trees, or NULL if all 1.0 */
   SCIP_EXPRGRAPHNODE**  rootnode,           /**< buffer to store expression graph node corresponding to root of expression tree */
   SCIP_Bool*            rootnodeisnew       /**< buffer to indicate whether the node in *rootnode has been newly created for this expression tree (otherwise, expression tree was already in graph) */
   );

/** replaces variable in expression graph by a linear sum of variables
 * variables will be added if not in the graph yet
 */
extern
SCIP_RETCODE SCIPexprgraphReplaceVarByLinearSum(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   void*                 var,                /**< variable to replace */
   int                   ncoefs,             /**< number of coefficients in linear term */
   SCIP_Real*            coefs,              /**< coefficients in linear term, or NULL if ncoefs == 0 */
   void**                vars,               /**< variables in linear term */
   SCIP_Real             constant            /**< constant offset */
   );

/** finds expression graph node corresponding to a variable */
extern
SCIP_Bool SCIPexprgraphFindVarNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   void*                 var,                /**< variable to search for */
   SCIP_EXPRGRAPHNODE**  varnode             /**< buffer to store node corresponding to variable, if found, or NULL if not found */
   );

/** finds expression graph node corresponding to a constant */
extern
SCIP_Bool SCIPexprgraphFindConstNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             constant,           /**< constant to search for */
   SCIP_EXPRGRAPHNODE**  constnode           /**< buffer to store node corresponding to constant, if found, or NULL if not found */
   );

/** get current maximal depth of expression graph */
extern
int SCIPexprgraphGetDepth(
   SCIP_EXPRGRAPH*       exprgraph           /**< expression graph */
   );

/** gets array with number of nodes at each depth of expression graph */
extern
int* SCIPexprgraphGetNNodes(
   SCIP_EXPRGRAPH*       exprgraph           /**< expression graph */
   );

/** gets nodes of expression graph, one array per depth */
extern
SCIP_EXPRGRAPHNODE*** SCIPexprgraphGetNodes(
   SCIP_EXPRGRAPH*       exprgraph           /**< expression graph */
   );

/** gets number of variables in expression graph */
extern
int SCIPexprgraphGetNVars(
   SCIP_EXPRGRAPH*       exprgraph           /**< pointer to expression graph that should be freed */
   );

/** gets array of variables in expression graph */
extern
void** SCIPexprgraphGetVars(
   SCIP_EXPRGRAPH*       exprgraph           /**< pointer to expression graph that should be freed */
   );

/** gets array of expression graph nodes corresponding to variables */
extern
SCIP_EXPRGRAPHNODE** SCIPexprgraphGetVarNodes(
   SCIP_EXPRGRAPH*       exprgraph           /**< pointer to expression graph that should be freed */
   );

/** prints an expression graph in dot format */
extern
SCIP_RETCODE SCIPexprgraphPrintDot(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   const char**          varnames            /**< variable names, or NULL for generic names */
   );

/** evaluates nodes of expression graph for given values of variables */
extern
SCIP_RETCODE SCIPexprgraphEval(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real*            varvals             /**< values for variables */
   );

/** sets value for a single variable given as expression graph node */
extern
void SCIPexprgraphSetVarNodeValue(
   SCIP_EXPRGRAPHNODE*   varnode,            /**< expression graph node corresponding to variable */
   SCIP_Real             value               /**< new value for variable */
   );

/** sets bounds for variables */
extern
void SCIPexprgraphSetVarsBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_INTERVAL*        varbounds           /**< new bounds for variables */
   );

/** sets bounds for a single variable */
extern
void SCIPexprgraphSetVarBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   void*                 var,                /**< variable */
   SCIP_INTERVAL         varbounds           /**< new bounds of variable */
   );

/** sets bounds for a single variable given as expression graph node */
extern
void SCIPexprgraphSetVarNodeBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   varnode,            /**< expression graph node corresponding to variable */
   SCIP_INTERVAL         varbounds           /**< new bounds of variable */
   );

/** sets lower bound for a single variable given as expression graph node */
extern
void SCIPexprgraphSetVarNodeLb(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   varnode,            /**< expression graph node corresponding to variable */
   SCIP_Real             lb                  /**< new lower bound for variable */
   );

/** sets upper bound for a single variable given as expression graph node */
extern
void SCIPexprgraphSetVarNodeUb(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   varnode,            /**< expression graph node corresponding to variable */
   SCIP_Real             ub                  /**< new upper bound for variable */
   );

/** gets bounds that are stored for all variables */
extern
SCIP_INTERVAL* SCIPexprgraphGetVarsBounds(
   SCIP_EXPRGRAPH*       exprgraph           /**< expression graph */
   );

/** propagates bound changes in variables forward through the expression graph */
extern
SCIP_RETCODE SCIPexprgraphPropagateVarBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Bool             clearreverseprop,   /**< whether to reset bound tightenings from reverse propagation */
   SCIP_Bool*            domainerror         /**< buffer to store whether a node with empty bounds has been found, propagation is interrupted in this case */
   );

/** propagates bound changes in nodes backward through the graph
 * new bounds are not stored in varbounds, but only in nodes corresponding to variables
 * NOTE: it is assumed that SCIPexprgraphPropagateVarBounds was called before if variable bounds were relaxed
 */
extern
void SCIPexprgraphPropagateNodeBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Real             minstrength,        /**< minimal required relative bound strengthening in a node to trigger a propagation into children nodes */
   SCIP_Bool*            cutoff              /**< buffer to store whether a node's bounds were propagated to an empty interval */
   );

/** updates curvature information in expression graph nodes w.r.t. currently stored variable bounds
 * implies update of bounds in expression graph
 */
extern
SCIP_RETCODE SCIPexprgraphCheckCurvature(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Bool             clearreverseprop    /**< whether to reset bound tightenings from reverse propagation */
   );

/** aims at simplifying an expression graph
 * a domain error can occur when variables were fixed to values for which a parent expression is not defined (e.g., 0^(-1) or log(-1))
 */
extern
SCIP_RETCODE SCIPexprgraphSimplify(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Real             eps,                /**< threshold, under which positive values are treat as 0 */
   int                   maxexpansionexponent,/**< maximal exponent for which we still expand non-monomial polynomials */
   SCIP_Bool*            havechange,         /**< buffer to indicate whether the graph has been modified */
   SCIP_Bool*            domainerror         /**< buffer to indicate whether a domain error has been encountered, i.e., some expressions turned into NaN */
   );

/** creates an expression tree from a given node in an expression graph */
extern
SCIP_RETCODE SCIPexprgraphGetTree(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   rootnode,           /**< expression graph node that should represent root of expression tree */
   SCIP_EXPRTREE**       exprtree            /**< buffer to store pointer to created expression tree */
   );

/** creates a sum of expression trees with pairwise disjoint variables from a given node in an expression graph
 * Giving SCIPexprgraphGetNodeNChildren() for exprtreesize is always sufficient.
 */
extern
SCIP_RETCODE SCIPexprgraphGetSeparableTrees(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node which represents expression to get */
   int                   exprtreessize,      /**< length of exprtrees and exprtreecoefs arrays, need to be at least one */
   int*                  nexprtrees,         /**< buffer to store number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< array where to store expression trees */
   SCIP_Real*            exprtreecoefs       /**< array where to store coefficients of expression trees */
   );

/** gives the number of summands which the expression of an expression graph node consists of */
extern
int SCIPexprgraphGetSumTreesNSummands(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   );

/** creates a sum of expression trees, possibly sharing variables, from a given node in an expression graph */
extern
SCIP_RETCODE SCIPexprgraphGetSumTrees(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node which represents expression to get */
   int                   exprtreessize,      /**< length of exprtrees and exptreecoefs arrays, should be at least SCIPexprgraphGetSumTreesNSummands() */
   int*                  nexprtrees,         /**< buffer to store number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< array where to store expression trees */
   SCIP_Real*            exprtreecoefs       /**< array where to store coefficients of expression trees */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif /* __NLPI_PUB_EXPR_H__ */
