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

/**@file   dpborder_util.c
 * @brief  Utility methods for dynamic programming solver for Steiner tree (sub-) problems with small border
 * @author Daniel Rehfeldt
 *
 * Implements utility methods.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include "dpborder.h"
#include "dpborderinterns.h"

/*
 * Local methods
 */


/*
 * Interface methods
 */


/** builds map between old and new border char representation */
void dpborder_buildBorderMap(
   DPBORDER*             dpborder            /**< border */
)
{
   const int extnode = dpborder_getTopLevel(dpborder)->extnode;
   const SCIP_Bool* const nodes_isBorder = dpborder->nodes_isBorder;
   int* RESTRICT bordercharmap = dpborder->bordercharmap;
   int nbordernew = 0;

   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      const int bordernode = dpborder->bordernodes[i];

      if( nodes_isBorder[bordernode] )
         bordercharmap[i] = nbordernew++;
      else
         bordercharmap[i] = -1;
   }

   if( dpborder->nodes_outdeg[extnode] != 0 )
   {
      nbordernew++;
   }

   /* now we set the delimiter */
   bordercharmap[StpVecGetSize(dpborder->bordernodes)] = nbordernew;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("char border map, old to new: \n");
   for( int i = 0; i < StpVecGetSize(dpborder->bordernodes); i++ )
   {
      SCIPdebugMessage("%d->%d \n", i, bordercharmap[i]);
   }
   SCIPdebugMessage("delimiter: %d->%d \n", StpVecGetSize(dpborder->bordernodes),
         bordercharmap[StpVecGetSize(dpborder->bordernodes)]);

#endif
}

/** builds distances to extension node */
void dpborder_buildBorderDists(
   const GRAPH*          graph,              /**< graph */
   DPBORDER*             dpborder            /**< border */
)
{
   SCIP_Real* RESTRICT borderchardists = dpborder->borderchardists;
   const SCIP_Bool* const nodes_isBorder = dpborder->nodes_isBorder;
   const int* const bordernodes = dpborder->bordernodes;
   const CSR* const csr = graph->csr_storage;
   const int* const start_csr = csr->start;
   const int extnode = dpborder_getTopLevel(dpborder)->extnode;
   const int nbordernodes = StpVecGetSize(bordernodes);

   SCIPdebugMessage("setting up border distances for extnode=%d \n", extnode);

   for( int i = 0; i < nbordernodes; i++ )
      borderchardists[i] = FARAWAY;

   for( int e = start_csr[extnode]; e != start_csr[extnode + 1]; e++ )
   {
      const int head = csr->head[e];
      if( nodes_isBorder[head] )
      {
         int i;
         for( i = 0; i < BPBORDER_MAXBORDERSIZE; i++ )
         {
            const int bordernode = bordernodes[i];
            assert(bordernode != extnode);

            if( head == bordernode )
            {
               SCIPdebugMessage("setting edge distance for border node %d to %f \n", head, csr->cost[e]);
               assert(EQ(borderchardists[i], FARAWAY));
               borderchardists[i] = csr->cost[e];
               break;
            }
         }
         assert(i != BPBORDER_MAXBORDERSIZE);
      }
   }
}
