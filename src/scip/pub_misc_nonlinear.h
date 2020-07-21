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

/**@file   pub_misc_nonlinear.h
 * @ingroup INTERNALAPI
 * @brief  internal miscellaneous methods for nonlinear constraints
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MISC_NONLINEAR_H__
#define __SCIP_MISC_NONLINEAR_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_cons.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns the right-hand side of an arbitrary SCIP constraint that can be represented as a single nonlinear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_Real SCIPconsNonlinearGetRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which right-hand side is queried */
   SCIP_Bool*            success             /**< pointer to store whether a valid right-hand side was returned */
   );

/** returns the left-hand side of an arbitrary SCIP constraint that can be represented as a single nonlinear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_Real SCIPconsNonlinearGetLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get left-hand side for */
   SCIP_Bool*            success             /**< pointer to store whether a valid left-hand side was returned */
   );

/** adds the given variable to the input constraint. */
SCIP_RETCODE SCIPconsNonlinearAddLinearCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which row is queried */
   SCIP_VAR*             var,                /**< variable of the constraint entry */
   SCIP_Real             val                 /**< the coefficient of the constraint entry */
   );

#ifdef __cplusplus
}
#endif

#endif
