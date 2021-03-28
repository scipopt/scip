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

/**@file   dpborder_base.c
 * @brief  Dynamic programming solver for Steiner tree (sub-) problems with small border
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dpborder.h"
#include "dpborderinterns.h"
#include "stpvector.h"


/*
 * Local methods
 */

/*
 * Interface methods
 */


/** initializes */
SCIP_RETCODE dpborder_init(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< original graph */
   DPBORDER**            dpborder            /**< to initialize */
)
{
   DPBORDER* dpb;
   SCIP_CALL( SCIPallocMemory(scip, dpborder) );
   dpb = *dpborder;

   assert(graph);

   dpb->nnodes = graph->knots;

   return SCIP_OKAY;
}


/** frees */
void dpborder_free(
   SCIP*                 scip,               /**< SCIP data structure */
   DPBORDER**            dpborder            /**< to be freed */
)
{
   SCIPfreeMemory(scip, dpborder);
}
