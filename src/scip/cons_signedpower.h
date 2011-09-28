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

/**@file   cons_signedpower.h
 * @brief  constraint handler for signedpower constraints
 * @author Stefan Vigerske
 * 
 * This constraint handler handles constraints of the form
 * \f[
 *   lhs \leq sign(x+offset) |x+offset|^n + c z \leq rhs
 * \f]
 * for n > 1.0 a rational number, c and offset arbitrary, and x and z variables.
 * x can have -offset in the interior of its domain.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_SIGNEDPOWER_H__
#define __SCIP_CONS_SIGNEDPOWER_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for signedpower constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrSignedpower(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a signedpower constraint */
extern
SCIP_RETCODE SCIPcreateConsSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             x,                  /**< nonlinear variable x in constraint */
   SCIP_VAR*             z,                  /**< linear variable z in constraint */
   SCIP_Real             exponent,           /**< exponent n of |x+offset|^n term in constraint */
   SCIP_Real             xoffset,            /**< offset in |x+offset|^n term in constraint */
   SCIP_Real             zcoef,              /**< coefficient of z in constraint */
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
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** gets the signedpower constraint as a nonlinear row representation */
extern
SCIP_RETCODE SCIPgetNlRowSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< a buffer where to store pointer to nonlinear row */
);

/** gets nonlinear variable x in signedpower constraint */
extern
SCIP_VAR* SCIPgetNonlinearVarSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< signedpower constraint */
);

/** gets linear variable z in signedpower constraint */
extern
SCIP_VAR* SCIPgetLinearVarSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< signedpower constraint */
);

/** gets exponent in power term in signedpower constraint */
extern
SCIP_Real SCIPgetExponentSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< signedpower constraint */
);

/** gets offset in power term in signedpower constraint */
extern
SCIP_Real SCIPgetOffsetSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< signedpower constraint */
);

/** gets coefficient of linear variable in signedpower constraint */
extern
SCIP_Real SCIPgetCoefLinearSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< signedpower constraint */
);

/** gets left hand side in signedpower constraint */
extern
SCIP_Real SCIPgetLhsSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< signedpower constraint */
);

/** gets right hand side in signedpower constraint */
extern
SCIP_Real SCIPgetRhsSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< signedpower constraint */
);

/** gets the absolute violation of a signedpower constraint by a solution */
extern
SCIP_Real SCIPgetViolationSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< signedpower constraint */
   SCIP_SOL*             sol                 /**< LP solution */
   );

#ifdef __cplusplus
}
#endif

#endif
