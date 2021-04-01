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

   if( partitionchars[0] == delimiter )
   {
      SCIPdebugMessage("partition starts with delimiter\n");
      return FALSE;
   }

   for( int i = 1; i < partsize; i++ )
   {
      if( partitionchars[i] == delimiter && partitionchars[i - 1] == delimiter )
      {
         SCIPdebugMessage("empty subset at %d %d \n", i - 1, i);
         return FALSE;
      }
   }

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
   assert(globalstart + partsize + 2 < dpborder->global_partcap);

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

   if( dpborder->extborderchar >= 0 )
   {
      assert(dpborder_getTopLevel(dpborder)->extnode == dpborder->bordernodes[dpborder->extborderchar]);
      global_partitions[globalend++] = dpborder->extborderchar;
   }

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
      StpVecPushBack(scip, dpborder->global_predparts, -1);
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
      printf("new (sub) partition (range %d-%d): \n", globalstart, globalend);
      dpborder_partPrint(&partition);
   }
#endif

   if( globalstart != -1 )
   {
      return dpborder->global_npartitions - 1;
   }

   return -1;
}



/** Gets global index of new global partition, similar to above, but merely removes prev. border
 *  nodes.
 *  Returns -1 if no valid partition could be built. */
int dpborder_partGetIdxNewExclusive(
   SCIP*                 scip,               /**< SCIP data structure */
   const DPBPART*        borderpartition,    /**< base partition */
   DPBORDER*             dpborder            /**< border */
)
{
   int globalstart = dpborder->global_partstarts[dpborder->global_npartitions];
   int globalend = globalstart;
   DPB_Ptype* RESTRICT global_partitions = dpborder->global_partitions;
   const int* const bordercharmap = dpborder->bordercharmap;
   const DPB_Ptype* const partitionchars = borderpartition->partchars;
   const DPB_Ptype delimiter_new = dpborder_getTopDelimiter(dpborder);
   const int partsize = borderpartition->partsize;

   assert(dpborder_partIsValid(borderpartition));
   assert(globalstart + partsize < dpborder->global_partcap);

   for( int i = 0; i < partsize; i++ )
   {
      const DPB_Ptype partchar = partitionchars[i];
      assert(0 <= partchar && partchar <= borderpartition->delimiter);
      assert(partchar != borderpartition->delimiter || bordercharmap[partchar] == delimiter_new);

      if( bordercharmap[partchar] != -1 )
         global_partitions[globalend++] = bordercharmap[partchar];
   }

   if( globalstart == globalend || global_partitions[globalstart] == delimiter_new )
   {
      assert(globalend - globalstart <= 1 || global_partitions[globalstart + 1] == delimiter_new);
      SCIPdebugMessage("exlusive sub-partition is invalid (empty)... \n");
      return -1;
   }

   for( int i = globalstart + 1; i != globalend; i++ )
   {
      if( global_partitions[i] == delimiter_new && global_partitions[i - 1] == delimiter_new )
      {
         SCIPdebugMessage("exlusive sub-partition is invalid (empty subset)... \n");
         return -1;
      }
   }

#ifdef SCIP_DEBUG
   {
      DPBPART partition;
      partition.partchars = &(global_partitions[globalstart]);
      partition.partsize = (globalend - globalstart);
      partition.delimiter = delimiter_new;
      printf("new (exclusive sub) partition (range %d-%d): \n", globalstart, globalend);
      dpborder_partPrint(&partition);
   }
#endif

   StpVecPushBack(scip, dpborder->global_partstarts, globalend);
   StpVecPushBack(scip, dpborder->global_partcosts, FARAWAY);
   StpVecPushBack(scip, dpborder->global_predparts, -1);
   dpborder->global_npartitions++;

   return dpborder->global_npartitions - 1;
}


/** marks optimal solution nodes */
void dpborder_markSolNodes(
   const DPBORDER*       dpborder,           /**< border */
   STP_Bool* RESTRICT    nodes_isSol         /**< solution nodes */
)
{

   const DPB_Ptype* const global_partitions = dpborder->global_partitions;
   const STP_Vectype(int) global_partstarts = dpborder->global_partstarts;
   const int nnodes = dpborder->nnodes;
   const int optposition = dpborder->global_optposition;
   int nlevels = 0;

   assert(dpborder && nodes_isSol);

   BMSclearMemoryArray(nodes_isSol, nnodes);
   SCIPdebugMessage("marking solution nodes: \n");

   for( int pos = optposition; pos != 0; pos = dpborder->global_predparts[pos] )
      nlevels++;

   for( int pos = optposition, level = nlevels; pos != 0; pos = dpborder->global_predparts[pos], level-- )
   {
      const int globalend = global_partstarts[pos + 1];
      const STP_Vectype(int) nodemap = dpborder->borderlevels[level]->bordernodesMapToOrg;
      const DPB_Ptype delimiter = dpborder_getDelimiter(dpborder, level);

      assert(level > 0);
      assert(0 <= pos && pos < dpborder->global_npartitions);
      assert(global_partstarts[pos + 1] > global_partstarts[pos]);
      assert(nodemap);
      assert(delimiter == StpVecGetSize(nodemap));

      SCIPdebugMessage("pos=%d, size %d range: %d-%d \n", pos,
            global_partstarts[pos + 1] - global_partstarts[pos], global_partstarts[pos], globalend);

      for( int i = global_partstarts[pos]; i != globalend; i++ )
      {
         const DPB_Ptype borderchar = global_partitions[i];
         assert(0 <= borderchar && borderchar <= delimiter);

         if( borderchar != delimiter )
         {
            const int node = nodemap[borderchar];
            SCIPdebugMessage("solnode=%d \n", node);
            assert(0 <= node && node < nnodes);

            nodes_isSol[node] = TRUE;
         }
      }
   }

   SCIPdebugMessage("final solnode=%d \n", dpborder->borderlevels[0]->extnode);
   nodes_isSol[dpborder->borderlevels[0]->extnode] = TRUE;
}
