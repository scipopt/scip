/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_history.h,v 1.9 2005/02/14 13:35:52 bzfpfend Exp $"

/**@file   struct_history.h
 * @brief  datastructures for branching and inference history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_HISTORY_H__
#define __STRUCT_HISTORY_H__


#include "scip/def.h"
#include "scip/type_history.h"


/** branching and inference history information for single variable */
struct History
{
   Real             pscostcount[2];     /**< nr of (partial) summands in down/upwards pseudo costs (may be fractional) */
   Real             pscostsum[2];       /**< sum of (partial) pseudo cost values for down/upwards branching */
   Longint          nbranchings[2];     /**< nr of times, the variable changed its bounds due to branching */
   Longint          ninferences[2];     /**< nr of times, branching on the variable lead to inference of another bound */
   Longint          ncutoffs[2];        /**< nr of times, branching on the variable lead to an infeasible sub problem */
   Longint          branchdepthsum[2];  /**< sum of depth levels, at which the branching bound changes took place */
};


#endif
