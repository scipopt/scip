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

/**@file   pub_misc_linear.h
 * @ingroup INTERNALAPI
 * @brief  internal miscellaneous methods for linear constraints
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MISC_LINEAR_H__
#define __SCIP_MISC_LINEAR_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns the right-hand side of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_EXPORT
SCIP_Real SCIPconsGetRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which right-hand side is queried */
   SCIP_Bool*            success             /**< pointer to store whether a valid right-hand side was returned */
   );

/** returns the left-hand side of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_EXPORT
SCIP_Real SCIPconsGetLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get left hand side for */
   SCIP_Bool*            success             /**< pointer to store whether a valid left-hand side was returned */
   );

/** returns the value array of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which the coefficients are wanted */
   SCIP_Real*            vals,               /**< array to store the coefficients of the constraint */
   int                   varssize,           /**< available slots in vals array needed to check if the array is large enough */
   SCIP_Bool*            success             /**< pointer to store whether the coefficients are successfully copied */
   );

/** returns the dual farkas solution of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the dual farkas solution
 */
SCIP_EXPORT
void SCIPconsGetDualfarkas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get left hand side for */
   SCIP_Real*            dualfarkas,         /**< pointer to store the dual farkas solution */
   SCIP_Bool*            success             /**< pointer to store whether the dual farkas solution is successfully returned */
   );

/** returns the dual solution of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the dual solution
 */
void SCIPconsGetDualsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get left hand side for */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution */
   SCIP_Bool*            success             /**< pointer to store whether the dual solution is successfully returned */
   );

/** returns the row of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *  or NULL of no row is awailable
 */
SCIP_EXPORT
SCIP_ROW* SCIPconsGetRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   );

/** adds the given variable to the input constraint.
 *  If the constraint is setppc or logicor the value is ignored. If the constraint is knapsack, then the value is
 *  converted to an int. A warning is passed if the SCIP_Real is not an integer.
 *  TODO: Allow val to be a pointer.
 */
SCIP_RETCODE SCIPconsAddCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which row is queried */
   SCIP_VAR*             var,                /**< variable of the constraint entry */
   SCIP_Real             val                 /**< the coefficient of the constraint entry */
   );

#ifdef __cplusplus
}
#endif

#endif
