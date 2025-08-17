/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_exactlinear.h
 * @ingroup CONSHDLRS
 * @brief  Constraint handler for linear constraints in their most general form, \f$lhs <= a^T x <= rhs\f$.
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXACTLINEAR_H__
#define __SCIP_CONS_EXACTLINEAR_H__

#include "scip/def.h"
#include "scip/intervalarith.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_lpexact.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"
#include "scip/type_certificate.h"
#include "scip/type_rational.h"

#ifdef __cplusplus
extern "C" {
#endif



/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrExactLinear(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Exact Linear Constraints
 *
 * This constraint handler handles linear constraints in their most general form
 * \f[
 *   lhs \leq \sum_{i=1}^n a_i x_i \leq rhs
 * \f]
 * in a numerically exact way, where \f$a_i \in Q, i = 1,\dots,n\f$, \f$lhs\in Q \cup \{-\infty\}\f$, \f$rhs\in Q \cup \{\infty\}\f$,
 * and \f$x_i, i = 1,\dots,n\f$ which can be binary, integer, or continuous decision variables.
 *
 * @{
 */

/** creates and captures a linear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_RATIONAL**       vals,               /**< array with coefficients of constraint entries */
   SCIP_RATIONAL*        lhs,                /**< left hand side of constraint */
   SCIP_RATIONAL*        rhs,                /**< right hand side of constraint */
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

/** creates and captures a linear constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLinear(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsLinear() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_RATIONAL**       vals,               /**< array with coefficients of constraint entries */
   SCIP_RATIONAL*        lhs,                /**< left hand side of constraint */
   SCIP_RATIONAL*        rhs                 /**< right hand side of constraint */
   );

/** creates a linear constraint from an exact linear constraint by rounding values to floating-point and captures it */
SCIP_EXPORT
SCIP_RETCODE SCIPcopyConsExactLinear(
   SCIP*                 scip,               /**< target SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to store the created target constraint */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in source variable array */
   SCIP_VAR**            sourcevars,         /**< source variables of the linear constraints */
   SCIP_INTERVAL*        sourcecoefs,        /**< coefficient array of the linear constraint, or NULL if all coefficients are one */
   SCIP_Real             lhs,                /**< left hand side of the linear constraint */
   SCIP_Real             rhs,                /**< right hand side of the linear constraint */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            valid               /**< pointer to store if the copying was valid */
   );

/** adds coefficient to linear constraint (if it is not zero) */
SCIP_EXPORT
SCIP_RETCODE SCIPaddCoefExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_RATIONAL*        val                 /**< coefficient of constraint entry */
   );

/** changes coefficient of variable in linear constraint; deletes the variable if coefficient is zero; adds variable if
 *  not yet contained in the constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint and variable.
 *
 *  @note This method requires linear time to search for occurences of the variable in the constraint data.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgCoefExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_RATIONAL*        val                 /**< new coefficient of constraint entry */
   );

/** deletes variable from linear constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint and variable.
 *
 *  @note This method requires linear time to search for occurences of the variable in the constraint data.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdelCoefExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   );

/** gets left hand side of linear constraint */
SCIP_EXPORT
SCIP_RATIONAL* SCIPgetLhsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets right hand side of linear constraint */
SCIP_EXPORT
SCIP_RATIONAL* SCIPgetRhsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** changes left hand side of linear constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPchgLhsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_RATIONAL*        lhs                 /**< new left hand side */
   );

/** changes right hand side of linear constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPchgRhsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_RATIONAL*        rhs                 /**< new right hand side */
   );

/** gets the number of variables in the linear constraint */
SCIP_EXPORT
int SCIPgetNVarsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the array of variables in the linear constraint; the user must not modify this array! */
SCIP_EXPORT
SCIP_VAR** SCIPgetVarsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the array of coefficient values in the linear constraint; the user must not modify this array! */
SCIP_EXPORT
SCIP_INTERVAL* SCIPgetValsRealExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the array of coefficient values in the linear constraint; the user must not modify this array! */
SCIP_EXPORT
SCIP_RATIONAL** SCIPgetValsExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the activity of the linear constraint in the given solution
 *
 *  @note if the activity comprises positive and negative infinity contributions, the result is currently undefined
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetActivityExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol,                /**< solution, or NULL to use current node's solution */
   SCIP_RATIONAL*        ret                 /**< pointer to store result */
   );

/** gets the feasibility of the linear constraint in the given solution */
SCIP_EXPORT
SCIP_RETCODE SCIPgetFeasibilityExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol,                /**< solution, or NULL to use current node's solution */
   SCIP_RATIONAL*        ret                 /**< pointer to store the result */
   );

/** gets the dual solution of the linear constraint in the current LP
 *
 *  @note this method currently returns the value from the floating-point LP
 */
SCIP_EXPORT
void SCIPgetFpDualsolExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_RATIONAL*        ret                 /**< pointer to store the result */
   );

/** gets the dual Farkas value of the linear constraint in the current infeasible LP
 *
 *  @note this method currently returns the value from the floating-point LP
 */
SCIP_EXPORT
void SCIPgetFpDualfarkasExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_RATIONAL*        ret                 /**< pointer to store the result */
   );

/** returns the linear relaxation of the given linear constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_EXPORT
SCIP_ROW* SCIPgetRowExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the exact linear relaxation of the given linear constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_EXPORT
SCIP_ROWEXACT* SCIPgetRowExactExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** prints the certificate for a given original exact linear constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPcertifyConsOrigExactLinear(
   SCIP*                 scip,
   SCIP_CONSHDLR*        conshdlr,
   SCIP_CONS*            cons
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
