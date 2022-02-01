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

/**@file   solhistory.h
 * @brief  includes methods working on the (reduction) history of solutions to Steiner tree problems
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef APPLICATIONS_STP_SRC_SOLHISTORY_H_
#define APPLICATIONS_STP_SRC_SOLHISTORY_H_

#include "scip/scip.h"
#include "graph.h"


#ifdef __cplusplus
extern "C" {
#endif


/** bi-decomposition reduction parameters */
typedef struct solution_history
{
   STP_Bool*             orgnodes_isInSol;   /**< node contained? */
   STP_Bool*             orgedges_isInSol;   /**< edge contained? */
   int                   nsolnodes;          /**< number */
   int                   nsoledges;          /**< number */
   int                   norgnodes;          /**< number */
   int                   norgedges;          /**< number */
} SOLHISTORY;


extern SCIP_RETCODE solhistory_init(SCIP*, const GRAPH*, SOLHISTORY**);
extern void solhistory_free(SCIP*, SOLHISTORY**);
extern SCIP_RETCODE solhistory_computeHistory(SCIP*, SCIP_SOL*, const GRAPH*, SOLHISTORY*);



#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_SOLHISTORY_H_ */
