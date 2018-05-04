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

/**@file   cons_benders.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for benders decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_BENDERS_H__
#define __SCIP_CONS_BENDERS_H__


#include "scip/scip.h"
#include "scip/benders.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for benders constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             twophase            /**< should the two phase method be used? */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Benders Constraints
 *
 * @{
 */

/** the methods for the enforcement of a SCIP solution */
EXTERN
SCIP_RETCODE SCIPconsBendersEnforceSolution(
   SCIP*                 scip,               /**< the SCIP instance */
   SCIP_SOL*             sol,                /**< the primal solution to enforce, or NULL for the current LP/pseudo sol */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_RESULT*          result,             /**< the result of the enforcement */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should integrality be considered when checking the subproblems */
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
