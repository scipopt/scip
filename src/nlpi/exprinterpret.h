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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   exprinterpret.h
 * @brief  methods to interpret (evaluate) an expression "fast"
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 * Realized similar to LPI: one implementation of an interpreter is linked in.
 */

/* @todo product Gradient times vector
   @todo product Hessian times vector
   @todo product Hessian of Lagrangian times vector
   @todo sparse Hessian of expression
   @todo sparse Hessian of Lagrangian (sets of expressions and quadratic parts)?
*/

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPRINTERPRET_H__
#define __SCIP_EXPRINTERPRET_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_expr.h"
#include "nlpi/type_exprinterpret.h"
#include "scip/intervalarith.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup EXPRINTS
 * @{
 */

/** gets name and version of expression interpreter */
SCIP_EXPORT
const char* SCIPexprintGetName(void);

/** gets descriptive text of expression interpreter */
SCIP_EXPORT
const char* SCIPexprintGetDesc(void);

/** gets capabilities of expression interpreter (using bitflags) */
SCIP_EXPORT
SCIP_EXPRINTCAPABILITY SCIPexprintGetCapability(void);

/** creates an expression interpreter object */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRINT**        exprint             /**< buffer to store pointer to expression interpreter */
   );

/** frees an expression interpreter object */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintFree(
   SCIP_EXPRINT**        exprint             /**< expression interpreter that should be freed */
   );

/** compiles an expression and stores compiled data in expression */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintCompile(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the capability to evaluate an expression by the expression interpreter
 *
 * In cases of user-given expressions, higher order derivatives may not be available for the user-expression,
 * even if the expression interpreter could handle these. This method allows to recognize that, e.g., the
 * Hessian for an expression is not available because it contains a user expression that does not provide
 * Hessians.
 */
SCIP_EXPORT
SCIP_EXPRINTCAPABILITY SCIPexprintGetExprCapability(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** frees interpreter data */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintFreeData(
   SCIP_EXPRINTDATA**    interpreterdata     /**< interpreter data that should freed */
   );

/** evaluates an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintEval(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Real*            val                 /**< buffer to store value of expression */
   );

/** evaluates an expression on intervals */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintEvalInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables */
   SCIP_INTERVAL*        val                 /**< buffer to store interval value of expression */
   );

/** computes value and gradient of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintGrad(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to a point evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store expression value */
   SCIP_Real*            gradient            /**< buffer to store expression gradient, need to have length at least SCIPgetExprNVars(expr) */
   );

/** computes interval value and interval gradient of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintGradInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable interval values changed since last call to an interval evaluation routine? */
   SCIP_INTERVAL*        val,                /**< buffer to store expression interval value */
   SCIP_INTERVAL*        gradient            /**< buffer to store expression interval gradient, need to have length at least SCIPgetExprNVars(expr) */
   );

/** gives sparsity pattern of lower-triangular part of hessian
 *
 * Since the AD code might need to do a forward sweep, you should pass variable values in here.
 *
 * Result will have (*colidxs)[i] <= (*rowidixs)[i] for i=0..*nnz.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintHessianSparsity(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            varvals,            /**< values of variables */
   int**                 rowidxs,            /**< buffer to return array with row indices of Hessian elements */
   int**                 colidxs,            /**< buffer to return array with column indices of Hessian elements */
   int*                  nnz                 /**< buffer to return length of arrays */
   );

/** computes value and hessian of an expression
 *
 * Returned arrays rowidxs and colidxs and number of elements nnz are the same as given by SCIPexprintHessianSparsity().
 * Returned array hessianvals will contain the corresponding Hessian elements.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintHessian(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store function value */
   int**                 rowidxs,            /**< buffer to return array with row indices of Hessian elements */
   int**                 colidxs,            /**< buffer to return array with column indices of Hessian elements */
   SCIP_Real**           hessianvals,        /**< buffer to return array with Hessian elements */
   int*                  nnz                 /**< buffer to return length of arrays */
   );

/** gives sparsity pattern of hessian
 *
 * NOTE: this function might be replaced later by something nicer.
 * Since the AD code might need to do a forward sweep, you should pass variable values in here.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintHessianSparsityDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Bool*            sparsity            /**< buffer to store sparsity pattern of Hessian, sparsity[i+n*j] indicates whether entry (i,j) is nonzero in the hessian */
   );

/** computes value and dense hessian of an expression
 *
 *  The full hessian is computed (lower left and upper right triangle).
 */
SCIP_EXPORT
SCIP_RETCODE SCIPexprintHessianDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store function value */
   SCIP_Real*            hessian             /**< buffer to store hessian values, need to have size at least n*n */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPRINTERPRET_H__ */
