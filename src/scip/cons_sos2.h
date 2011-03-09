/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_sos2.h
 * @brief  constraint handler for SOS type 2 constraints
 * @author Marc Pfetsch
 *
 * A specially ordered set of type 2 (SOS2) is a sequence of variables such that at most two
 * variables are nonzero and if two variables are nonzero they must be adjacent in the specified
 * sequence. Note that it is in principle allowed that a variables appears twice, but it then can be
 * fixed to 0 if it is at least two apart in the sequence.
 *
 * This constraint is useful when considering a piecewise affine approximation of a univariate
 * (nonlinear) function \f$: [a,b] \rightarrow R\f$: Let \f$x_1 < \ldots < x_n\f$ be points in
 * \f$[a,b]\f$ and introduce variables \f$\lambda_1, \ldots, \lambda_n\f$. To evaluate \f$f(x')\f$
 * at some point \f$x' \in [a,b]\f$ one can use the following constraints:
 * \f[
 *     \lambda_1 + \cdots + \lambda_n = 1,\quad x' = x_1 \lambda_1 + \cdots + x_n \lambda_n.
 * \f]
 * The value of \f$f(x')\f$ can the be approximated as
 * \f[
 *     f(x_1) \lambda_1 + \cdots + f(x_n) \lambda_n.
 * \f]
 * To get a valid piecewise affine approximation, \f$\lambda_1, \ldots, \lambda_n\f$ have to obey an
 * SOS constraint of type 2.
 *
 * This implementation of this constraint handler is based on classical ideas, see e.g.@n
 *  "Special Facilities in General Mathematical Programming System for
 *  Non-Convex Problems Using Ordered Sets of Variables"@n
 *  E. Beale and J. Tomlin, Proc. 5th IFORS Conference, 447-454 (1970)
 *
 * The order of the variables is determined as follows:
 *
 * - If the constraint is created with SCIPcreateConsSOS2() and weights are given, the weights
 *   determine the order (decreasing weights). Additional variables can be added with
 *   SCIPaddVarSOS2(), which adds a variable with given weight.
 *
 * - If an empty constraint is created and then variables are added with SCIPaddVarSOS2(), weights
 *   are needed and stored.
 *
 * - All other calls ignore the weights, i.e., if an nonempty constraint is created or variables are
 *   added with SCIPappendVarSOS2().
 *
 * @todo Allow to adapt the order of the constraints, e.g. by priorities. This for instance
 *   determines the branching order.
 * @todo Separate the following cuts for each pair of variables x, y of at least distance 2 in the
 *   SOS2 constraint: \f$ \min \{l_x, l_y\} \leq x + y \leq \max \{u_x, u_y\}\f$, where \f$l_x, u_x,
 *   l_y, u_y\f$ are the lower and upper bounds of x and y, respectively.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_SOS2_H__
#define __SCIP_CONS_SOS2_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for SOS2 constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrSOS2(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures an SOS2 constraint */
extern
SCIP_RETCODE SCIPcreateConsSOS2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            weights,            /**< weights determining the variable order, or NULL if natural order should be used */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** adds variable to SOS2 constraint, the position is determined by the given weight */
extern
SCIP_RETCODE SCIPaddVarSOS2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_Real             weight              /**< weight determining position of variable */
   );

/** appends variable to SOS2 constraint */
extern
SCIP_RETCODE SCIPappendVarSOS2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   );

/** gets number of variables in SOS2 constraint */
extern
int SCIPgetNVarsSOS2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets array of variables in SOS2 constraint */
SCIP_VAR** SCIPgetVarsSOS2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets array of weights in SOS2 constraint (or NULL if not existent) */
SCIP_Real* SCIPgetWeightsSOS2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

#ifdef __cplusplus
}
#endif

#endif
