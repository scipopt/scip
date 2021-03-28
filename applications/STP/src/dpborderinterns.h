/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dptermsinterns.h
 * @brief  Dynamic programming internals for Steiner tree (sub-) problems with small number of terminals
 * @author Daniel Rehfeldt
 *
 * Internal methods and data structures for DP.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef APPLICATIONS_STP_SRC_DPBORDERINTERNS_H_
#define APPLICATIONS_STP_SRC_DPBORDERINTERNS_H_


#include "scip/scip.h"
#include "graph.h"
#include "stpvector.h"

/** DP border structure */
struct dynamic_programming_border
{
   int                   nnodes;             /**< number of nodes of underlying graph */
};




#endif /* APPLICATIONS_STP_SRC_DPBORDERINTERNS_H_ */
