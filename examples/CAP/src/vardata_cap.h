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

/**@file   vardata_cap.h
 * @brief  Variable data containing the ids of constraints in which the variable appears
 * @author Stephen J. Maher
 *
 * This file implements the handling of the variable data which is attached to each file. See SCIP_VarData and \ref CAP_PRICER.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VARDATA_CAP__
#define __SCIP_VARDATA_CAP__

#include "scip/scip.h"

enum varprob
{
   MASTER = 0,
   SUBPROB = 1
};
typedef enum varprob VARPROB;

/** create variable data */
extern
SCIP_RETCODE SCIPvardataCreateCAP(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   VARPROB               prob,               /**< the problem where this variable exists */
   int                   nvars               /**< the number of corresponding variables,
                                                  this is the number of subproblems or 1 if the variable is a subproblem var */
   );

/** creates variable */
SCIP_RETCODE SCIPcreateVarCAP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< the lower bound */
   SCIP_Real             ub,                 /**< the upper bound */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          type                /**< the variable type */
   );

/** adds a variable to the variable mapping array */
extern
void SCIPvardataAddVarMapping(
   SCIP_VARDATA*         vardata,            /**< variable data */
   SCIP_VAR*             mappedvar,          /**< the variable that is mapped to this variable */
   int                   probnum             /**< the problem number that the mapped var belongs to */
   );

/** returns the mapped variable for a given problem */
extern
SCIP_VAR* SCIPvardataGetMappedVar(
   SCIP_VARDATA*         vardata,            /**< variable data */
   int                   probnum             /**< the problem number for the required mapped var, -1 if the master prob */
   );

/** returns the mapped variable for a given problem */
extern
VARPROB SCIPvardataGetProb(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

#endif
