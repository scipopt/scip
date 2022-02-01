/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dpborder.h
 * @brief  Dynamic programming solver for Steiner tree (sub-) problems with small border
 * @author Daniel Rehfeldt
 *
 * This file implements a dynamic programming method from Polzin and Vahdati to solve Steiner tree problems to optimality.
 * See also "Practical Partitioning-Based Methods for the Steiner Problem", WEA 2006
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_DPBORDER_H_
#define APPLICATIONS_STP_SRC_DPBORDER_H_

#include "scip/scip.h"
#include "graph.h"



#ifdef __cplusplus
extern "C" {
#endif

typedef struct dynamic_programming_border DPBORDER;


extern SCIP_RETCODE     dpborder_init(SCIP*, const GRAPH*, DPBORDER**);
extern SCIP_RETCODE     dpborder_probePotential(SCIP*, GRAPH*, DPBORDER*, SCIP_Bool*);
extern SCIP_RETCODE     dpborder_solve(SCIP*, GRAPH*, DPBORDER*, int*, SCIP_Bool*);
extern void             dpborder_free(SCIP*, DPBORDER**);



#ifdef __cplusplus
}
#endif

#endif /* APPLICATIONS_STP_SRC_DPBORDER_H_ */
