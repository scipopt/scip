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


/** is partition valid? */
SCIP_Bool dpborder_partIsValid(
   const DPBPART*        borderpartition     /**< partition */
)
{
   const DPB_Ptype* const partitionchars = borderpartition->partchars;
   const DPB_Ptype delimiter = borderpartition->delimiter;
   const int partsize = borderpartition->partsize;

   for( int i = 0; i < partsize; i++ )
   {
      const DPB_Ptype borderchar = partitionchars[i];
      if( borderchar > delimiter )
      {
         SCIPdebugMessage("char %d to large \n", i);
         return FALSE;
      }

      if( borderchar == delimiter )
         continue;

      for( int j = 0; j < partsize; j++ )
      {
         const DPB_Ptype borderchar2 = partitionchars[j];

         if( i != j && borderchar == borderchar2 )
         {
            SCIPdebugMessage("duplicate char, positions %d %d \n", i, j);
            return FALSE;
         }
      }
   }

   return TRUE;
}


/** prints partition */
void dpborder_partPrint(
   const DPBPART*        borderpartition     /**< partition */
)
{
   const DPB_Ptype* const partitionchars = borderpartition->partchars;
   const DPB_Ptype delimiter = borderpartition->delimiter;
   const int partsize = borderpartition->partsize;

   for( int i = 0; i < partsize; i++ )
   {
      const DPB_Ptype borderchar = partitionchars[i];

      if( borderchar == delimiter )
      {
         printf("X ");
         continue;
      }

      printf("%d ", borderchar);
   }

   printf(" \n");

   assert(dpborder_partIsValid(borderpartition));
}


/** gets candidates start for given partition */
STP_Vectype(int) dpborder_partGetCandstarts(
   SCIP*                 scip,               /**< SCIP data structure */
   const DPBPART*        borderpartition,    /**< partition */
   const DPBORDER*       dpborder            /**< border */
)
{
   STP_Vectype(int) candstarts = NULL;
   const DPB_Ptype* const partitionchars = borderpartition->partchars;
   const DPB_Ptype delimiter = borderpartition->delimiter;
   const int partsize = borderpartition->partsize;
   const SCIP_Real* const borderchardists = dpborder->borderchardists;

   assert(dpborder_partIsValid(borderpartition));

   for( int i = 0; i < partsize; i++ )
   {
      const DPB_Ptype borderchar = partitionchars[i];
      assert(0 <= borderchar && borderchar <= delimiter);

      if( borderchar == delimiter )
         continue;

      if( LT(borderchardists[borderchar], FARAWAY) )
      {
         int startpos;
         for( startpos = i; startpos > 0; startpos-- )
         {
            if( partitionchars[startpos] == delimiter )
               break;
         }

         if( partitionchars[startpos] == delimiter )
            startpos++;

         StpVecPushBack(scip, candstarts, startpos);

         /* move to next set of the partition */
         for( ; i < partsize && partitionchars[i] != delimiter; i++ )
         {
            assert(partitionchars[i] < delimiter);
         }
      }
   }

   return candstarts;
}


/** gets cardinality from global index of new global partition. */
int dpborder_partglobalGetCard(
   int                   globalindex,        /**< global index */
   int                   delimiter,          /**< delimiter */
   const DPBORDER*       dpborder            /**< border */
)
{
   int card = 1;
   const int globalstart = dpborder->global_partstarts[globalindex];
   const int globalend = dpborder->global_partstarts[globalindex + 1];
   const DPB_Ptype* const global_partitions = dpborder->global_partitions;

   assert(0 <= globalindex && globalindex < dpborder->global_npartitions);
   assert(0 <= delimiter);
   assert(globalstart < globalend);

   for( int i = globalstart; i != globalend; i++ )
   {
      if( global_partitions[i] == delimiter )
         card++;
   }

   return card;
}


/** gets minimum connection cost of connection selected sets of partition to extension vertex */
SCIP_Real dpborder_partGetConnectionCost(
   const DPBORDER*       dpborder,           /**< border */
   const DPBPART*        borderpartition,    /**< base partition */
   const int*            candstarts_sub,     /**< candidate starts from which to construct new partition */
   int                   ncandstarts_sub     /**< number of candidate starts */
)
{
   SCIP_Real costsum = 0.0;
   const SCIP_Real* const borderchardists = dpborder->borderchardists;
   const DPB_Ptype* const partitionchars = borderpartition->partchars;
   const DPB_Ptype delimiter_prev = borderpartition->delimiter;
   const int partsize = borderpartition->partsize;

   assert(dpborder_partIsValid(borderpartition));

   for( int i = 0; i < ncandstarts_sub; i++ )
   {
      SCIP_Real minedgecost = FARAWAY;
      const int candstart = candstarts_sub[i];
      assert(0 <= candstart && candstart < partsize);

      for( int j = candstart; j < partsize; j++ )
      {
         const DPB_Ptype partchar = partitionchars[j];
         assert(0 <= partchar && partchar <= delimiter_prev);

         if( partchar == delimiter_prev )
            break;

         if( LT(borderchardists[partchar], minedgecost) )
            minedgecost = borderchardists[partchar];
      }

      costsum += minedgecost;

      if( GE(costsum, FARAWAY) )
         break;
   }

   return costsum;
}


/** Gets global index of new global partition.
 *  Returns -1 if no valid partition could be built. */
int dpborder_partGetIdxNew(
   SCIP*                 scip,               /**< SCIP data structure */
   const DPBPART*        borderpartition,    /**< base partition */
   const int*            candstarts_sub,     /**< candidate starts from which to construct new partition */
   int                   ncandstarts_sub,    /**< number of candidate starts */
   DPBORDER*             dpborder            /**< border */
)
{
   int i;
   int globalstart = dpborder->global_partstarts[dpborder->global_npartitions];
   int globalend = globalstart;
   DPB_Ptype* RESTRICT global_partitions = dpborder->global_partitions;
   const int* const bordercharmap = dpborder->bordercharmap;
   DPB_Ptype* RESTRICT partitionchars = borderpartition->partchars;
   const DPB_Ptype delimiter_prev = borderpartition->delimiter;
   const DPB_Ptype delimiter_new = dpborder_getTopDelimiter(dpborder);
   const int partsize = borderpartition->partsize;
   SCIP_Bool doCopy;

   assert(dpborder_partIsValid(borderpartition));

   for( i = 0; i < ncandstarts_sub; i++ )
   {
      const int candstart = candstarts_sub[i];
      assert(0 <= candstart && candstart < partsize);

      for( int j = candstart; j < partsize; j++ )
      {
         const DPB_Ptype partchar = partitionchars[j];
         assert(0 <= partchar && partchar <= delimiter_prev);

         if( partchar == delimiter_prev )
            break;

         if( bordercharmap[partchar] != -1 )
            global_partitions[globalend++] = bordercharmap[partchar];
      }

      assert(partitionchars[candstart] < delimiter_prev);

      /* we mark the starts to skip them later on */
      partitionchars[candstart] = -partitionchars[candstart] - 1;
      assert(partitionchars[candstart] < 0);
   }

   /* adds char for extension node... */
   global_partitions[globalend++] = dpborder_getTopLevel(dpborder)->nbordernodes - 1;
   assert(dpborder_getTopLevel(dpborder)->extnode
      == dpborder->bordernodes[dpborder_getTopLevel(dpborder)->nbordernodes - 1]);

   doCopy = TRUE;

   if( partitionchars[0] >= 0 )
   {
      assert(delimiter_prev != partitionchars[0]);
      global_partitions[globalend++] = delimiter_new;
   }

   /* now we add the remaining sets of the partition */
   for( i = 0; i < partsize; i++ )
   {
      const DPB_Ptype partchar = partitionchars[i];
      if( partchar < 0 )
      {
         partitionchars[i] = -(partitionchars[i] + 1);
         doCopy = FALSE;
         continue;
      }

      if( partchar == delimiter_prev )
      {
         assert(i < partsize);

         if( partitionchars[i + 1] >= 0 )
         {
            /* empty subset? */
            if( delimiter_new == global_partitions[globalend - 1] )
            {
               globalstart = -1;
               break;
            }

            global_partitions[globalend++] = delimiter_new;
            doCopy = TRUE;
         }
         continue;
      }

      if( !doCopy )
         continue;

      if( bordercharmap[partchar] != -1 )
         global_partitions[globalend++] = bordercharmap[partchar];
   }

   if( delimiter_new == global_partitions[globalend - 1] )
      globalstart = -1;

   if( globalstart == -1 )
   {
      for( int j = i + 1; j < partsize; j++ )
      {
         if( partitionchars[j] < 0 )
            partitionchars[j] = -(partitionchars[j] + 1);
      }
   }
   else
   {
      StpVecPushBack(scip, dpborder->global_partstarts, globalend);
      StpVecPushBack(scip, dpborder->global_partcosts, FARAWAY);
      dpborder->global_npartitions++;
   }


#ifndef NDEBUG
   for( int j = 0; j < partsize; j++ )
      assert(0 <= partitionchars[j] && partitionchars[j] <= delimiter_prev);
#endif

#ifdef SCIP_DEBUG
   if( globalstart != -1 )
   {
      DPBPART partition;
      partition.partchars = &(global_partitions[globalstart]);
      partition.partsize = (globalend - globalstart);
      partition.delimiter = delimiter_new;
      printf("new (sub) partition: \n");
      dpborder_partPrint(&partition);
   }
#endif

   if( globalstart != -1 )
   {
      return dpborder->global_npartitions - 1;
   }

   return -1;
}


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
