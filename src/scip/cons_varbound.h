/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_varbound.h
 * @brief  constraint handler for variable bound constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_VARBOUND_H__
#define __SCIP_CONS_VARBOUND_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for variable bound constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrVarbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a variable bound constraint: lhs <= x + c*y <= rhs */
extern
SCIP_RETCODE SCIPcreateConsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             vbdcoef,            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs,                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs,                /**< right hand side of variable bound inequality */
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
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** gets left hand side of variable bound constraint lhs <= x + c*y <= rhs */
extern
SCIP_Real SCIPgetLhsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets right hand side of variable bound constraint lhs <= x + c*y <= rhs */
extern
SCIP_Real SCIPgetRhsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets bounded variable x of variable bound constraint lhs <= x + c*y <= rhs */
extern
SCIP_VAR* SCIPgetVarVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets bounding variable y of variable bound constraint lhs <= x + c*y <= rhs */
extern
SCIP_VAR* SCIPgetVbdvarVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets bound coefficient c of variable bound constraint lhs <= x + c*y <= rhs */
extern
SCIP_Real SCIPgetVbdcoefVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the dual solution of the variable bound constraint in the current LP */
extern
SCIP_Real SCIPgetDualsolVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the dual Farkas value of the variable bound constraint in the current infeasible LP */
extern
SCIP_Real SCIPgetDualfarkasVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the linear relaxation of the given variable bound constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
extern
SCIP_ROW* SCIPgetRowVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

#ifdef __cplusplus
}
#endif

#endif
