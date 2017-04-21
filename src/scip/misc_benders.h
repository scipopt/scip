/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    benders_misc.h
 * @ingroup PUBLICMETHODS
 * @brief   various SCIP helper methods for Benders' decomposition
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SCIP_BENDERS_MISC_H__
#define SCIP_BENDERS_MISC_H__

#include "scip/scip.h"
#include "scip/cons_setppc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** constraint types */
typedef enum  {
   linear, knapsack, varbound, setpacking, setcovering, setpartitioning,
   logicor, sos1, sos2, unknown, nconsTypeItems
} consType;


/** returns TRUE if variable is relevant, FALSE otherwise */
extern
SCIP_Bool BDisVarRelevant(
   SCIP_VAR*             var                 /**< variable to test */
   );

/** returns the type of an arbitrary SCIP constraint */
extern
consType BDconsGetType(
   SCIP_CONS*            cons                /**< constraint to get type for */
   );

/** returns the rhs of an arbitrary SCIP constraint */
extern
SCIP_Real BDconsGetRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   );

/** returns the lhs of an arbitrary SCIP constraint */
extern
SCIP_Real BDconsGetLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   );

/** returns the dual farkas sol of an arbitrary SCIP constraint */
extern
SCIP_Real BDconsGetDualfarkas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   );

/** returns the dual sol of an arbitrary SCIP constraint */
extern
SCIP_Real BDconsGetDualsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   );

/** returns the number of variables in an arbitrary SCIP constraint */
extern
int BDconsGetNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get number of variables */
   );

/** returns the variable array of an arbitrary SCIP constraint */
extern
SCIP_RETCODE BDconsGetVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get variables from */
   SCIP_VAR**            vars,               /**< array where variables are stored */
   int                   nvars               /**< size of storage array */
   );

/** returns the value array of an arbitrary SCIP constraint */
extern
SCIP_RETCODE BDconsGetVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get values from */
   SCIP_Real*            vals,               /**< array where values are stored */
   int                   nvals               /**< size of storage array */
   );

/** returns true if the constraint should be a master constraint and false otherwise */
SCIP_Bool BDgetConsIsSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SETPPCTYPE*      setppctype          /**< returns the type of the constraints */
   );

/** returns TRUE or FALSE, depending whether we are in the root node or not */
extern
SCIP_Bool BDisRootNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
SCIP_RETCODE BDincludeDialogsGraph(
   SCIP* scip
   );

#ifdef __cplusplus
}
#endif

#endif /* SCIP_BENDERS_MISC_H_ */
