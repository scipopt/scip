/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_history.h,v 1.3 2004/04/06 15:21:07 bzfpfend Exp $"

/**@file   struct_history.h
 * @brief  datastructures for branching and inference history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_HISTORY_H__
#define __STRUCT_HISTORY_H__


#include "def.h"
#include "type_history.h"


/** branching and inference history information for single variable */
struct History
{
   Real             pscostcount[2];     /**< nr of (partial) summands in down/upwards pseudo costs (may be fractional) */
   Real             pscostsum[2];       /**< sum of (partial) pseudo cost values for down/upwards branching */
   Longint          nbranchings;        /**< nr of times, the variable changed its bounds due to branching */
   Longint          ninferences;        /**< nr of times, branching on the variable lead to inference of another bound */
};


#endif
