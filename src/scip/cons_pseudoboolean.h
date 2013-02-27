/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_pseudoboolean.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for pseudoboolean constraints
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 * The constraint handler deals with pseudo boolean constraints. These are constraints of the form
 * \f[
 * \mbox{lhs} \leq \sum_{k=0}^m c_k \cdot x_k  +  \sum_{i=0}^n c_i \cdot \prod_{j \in I_i} x_j \leq \mbox{rhs}
 * \f]
 * where all \f$x\f$ are binary and all \f$c\f$ are integer.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_PSEUDOBOOLEAN_H__
#define __SCIP_CONS_PSEUDOBOOLEAN_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for pseudoboolean constraints and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrPseudoboolean(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a pseudoboolean constraint
 *
 *  @note linear and nonlinear terms can be added using SCIPaddCoefPseudoboolean() and SCIPaddTermPseudoboolean(),
 *        respectively
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            linvars,            /**< variables of the linear part, or NULL */
   int                   nlinvars,           /**< number of variables of the linear part */
   SCIP_Real*            linvals,            /**< coefficients of linear part, or NULL */
   SCIP_VAR***           terms,              /**< nonlinear terms of variables, or NULL */
   int                   nterms,             /**< number of terms of variables of nonlinear term */
   int*                  ntermvars,          /**< number of variables in nonlinear terms, or NULL */
   SCIP_Real*            termvals,           /**< coefficients of nonlinear parts, or NULL */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures a pseudoboolean constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values, which can be set
 *  afterwards using SCIPsetConsFLAGNAME() in scip.h
 *
 *  @see SCIPcreateConsPseudoboolean() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            linvars,            /**< variables of the linear part, or NULL */
   int                   nlinvars,           /**< number of variables of the linear part */
   SCIP_Real*            linvals,            /**< coefficients of linear part, or NULL */
   SCIP_VAR***           terms,              /**< nonlinear terms of variables, or NULL */
   int                   nterms,             /**< number of terms of variables of nonlinear term */
   int*                  ntermvars,          /**< number of variables in nonlinear terms, or NULL */
   SCIP_Real*            termvals,           /**< coefficients of nonlinear parts, or NULL */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** adds linear term pseudo boolean constraint (if it is not zero) */
EXTERN
SCIP_RETCODE SCIPaddCoefPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   );

/** adds nonlinear term to pseudo boolean constraint (if it is not zero) */
EXTERN
SCIP_RETCODE SCIPaddTermPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR**            vars,               /**< variables of the nonlinear term */
   int                   nvars,              /**< number of variables of the nonlinear term */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   );

/** gets left hand side of pseudo boolean constraint */
EXTERN
SCIP_Real SCIPgetLhsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets right hand side of pseudoboolean constraint */
EXTERN
SCIP_Real SCIPgetRhsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets indicator variable of pseudoboolean constraint, or NULL if there is no */
EXTERN
SCIP_VAR* SCIPgetIndVarPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets linear variables of pseudoboolean constraint */
EXTERN
SCIP_VAR** SCIPgetLinearVarsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets linear coefficients of pseudoboolean constraint */
EXTERN
SCIP_Real* SCIPgetLinearValsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets number of linear variables of pseudoboolean constraint */
EXTERN
int SCIPgetNLinearVarsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets array with non-linear term variables and corresponding number of variable in each non-linear term iff
 *  'termvarssize' (size of array for all non-linear term variables) is big enough, or will return the needed size for
 *  the array for all non-linear variables in termvarssize otherwise, of a given pseudoboolean constraint
 */
EXTERN
SCIP_RETCODE SCIPgetTermVarsDataPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR***           termvars,           /**< array for all arrays of term variables */
   int*                  ntermvars,          /**< array for number of variable in each term */
   int*                  termvarssize        /**< pointer to store size of array which is needed for all term variables */
   );

/** gets non-linear coefficients of pseudoboolean constraint */
EXTERN
SCIP_Real* SCIPgetTermValsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets number of non-linear coefficients of pseudoboolean constraint */
EXTERN
int SCIPgetNTermValsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** changes left hand side of pseudoboolean constraint */
EXTERN
SCIP_RETCODE SCIPchgLhsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of pseudoboolean constraint */
EXTERN
SCIP_RETCODE SCIPchgRhsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right hand side */
   );

#ifdef __cplusplus
}
#endif

#endif
