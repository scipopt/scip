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

/**@file   struct_history.h
 * @brief  datastructures for branching and inference history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_HISTORY_H__
#define __SCIP_STRUCT_HISTORY_H__


#include "scip/def.h"
#include "scip/type_history.h"

#ifdef __cplusplus
extern "C" {
#endif

/** branching and inference history information for single variable */
struct SCIP_History
{
   SCIP_Real             pscostcount[2];     /**< nr of (partial) summands in down/upwards pseudo costs (may be fractional) */
   SCIP_Real             pscostsum[2];       /**< sum of (partial) pseudo cost values for down/upwards branching */
   SCIP_Real             vsids[2];           /**< degree of how often the variable was reason for a conflict */
   SCIP_Real             conflengthsum[2];   /**< overall length of all active conflicts for which the variable gave reason */
   SCIP_Real             inferencesum[2];    /**< degree of how often branching on the variable lead to inference of another bound */
   SCIP_Real             cutoffsum[2];       /**< degree of how often branching on the variable lead to an infeasible sub problem */
   SCIP_Longint          nactiveconflicts[2];/**< number of active conflicts for which the variable gave reason */
   SCIP_Longint          nbranchings[2];     /**< nr of times, the variable changed its bounds due to branching */
   SCIP_Longint          branchdepthsum[2];  /**< sum of depth levels, at which the branching bound changes took place */
};

#ifdef __cplusplus
}
#endif

#endif
