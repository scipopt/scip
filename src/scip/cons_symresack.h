/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_symresack.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for symresack constraints
 * @author Christopher Hojny
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_SYMRESACK_H__
#define __SCIP_CONS_SYMRESACK_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Symresack Constraints
 *
 * @{
 *
 * Given a permutation that acts on the order of the variables of a (mixed) 0/1-program
 * such that the permutation is a symmetry of the program, this constraint handler can
 * be used to handle the symmetries corresponding to the permutation. The symmetries
 * are handled by enforcing that a binary solution is lexicographically not smaller than
 * its permutation.
 *
 * Moreover, the constraint handler checks whether each cycle of the permutation is
 * contained in a set packing or partitioning constraint. In this case, the symresack
 * is upgraded to a ppsymresack and strong symmetry handling inequalities are added during
 * the initialization of the constraint handler.
 *
 * Furthermore, the constraint handler checks whether the permutation is a composition of
 * 2-cycles. In this case, the symresack is a so-called orbisack and will be treated by
 * specialized methods of the orbisack constraint handler.
 */

/** creates the handler for symresack constraints and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrSymresack(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** creates and captures a symresack constraint
 *
 *  In a presolving step, we check whether the permutation acts only on binary points. Otherwise, we eliminate
 *  the non-binary variables from the permutation.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in problem */
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
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   SCIP_Bool*            success             /**< pointer to store whether permutation was created */
   );

/** creates and captures a symresack constraint
 *  in its most basic variant, i.e., with all constraint flags set to their default values, which can be set
 *  afterwards using SCIPsetConsFLAGNAME() in scip.h
 *
 *  In a presolving step, we remove all fixed points and cycles that act on non-binary variables of the permutation
 *
 *  @see SCIPcreateConsSymresack() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in problem */
   SCIP_Bool*            success             /**< pointer to store whether permutation was created */
   );

#ifdef __cplusplus
}
#endif

#endif
