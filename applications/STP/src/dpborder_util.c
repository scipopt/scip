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

/**@file   dpborder_util.c
 * @brief  Utility methods for dynamic programming solver for Steiner tree (sub-) problems with small border
 * @author Daniel Rehfeldt
 *
 * Implements utility methods.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include "dpborder.h"
#include "dpborderinterns.h"

/*
 * Local methods
 */

#ifndef NDEBUG
static
SCIP_Bool partitionIsIncluded(
   const DPBORDER*       dpborder,           /**< border */
   const DPB_Ptype*      partition,          /**< array of size 'size' */
   int                   size                /**< size */
)
{
   const int levelstart = dpborder_getTopLevel(dpborder)->globalstartidx;
   const int levelend = dpborder->global_npartitions - 1;

   assert(levelstart <= levelend);

   for( int i = levelstart; i != levelend; i++ )
   {
      const int start = dpborder->global_partstarts[i];
      const int end = dpborder->global_partstarts[i + 1];
      assert(start < end);

      if( size != end - start )
         continue;

      if( memcmp(partition, &(dpborder->global_partitions[start]), size) == 0 )
      {
         SCIPdebugMessage("included at pos=%d \n", i);
         return TRUE;
      }
   }

   return FALSE;
}



/** fully sorted? */
static
SCIP_Bool partitionIsSorted(
   const DPB_Ptype*      partition,          /**< array of size 'size' with delimiter at [size] */
   DPB_Ptype             delimiter,          /**< delimiter */
   int                   size                /**< size */
)
{
   int nsubsets = 0;
   int substarts[BPBORDER_MAXBORDERSIZE + 1];

   assert(size > 0);
   assert(partition[size] == delimiter);

   substarts[0] = 0;

   for( int iter = 1; iter < size; iter++ )
   {
      if( partition[iter - 1] == delimiter )
      {
         substarts[++nsubsets] = iter;
         continue;
      }

      if( partition[iter - 1] >= partition[iter] )
      {
         SCIPdebugMessage("unsorted subset \n");
         return FALSE;
      }
   }
   substarts[++nsubsets] = size + 1;

   for( int iter = 1; iter < nsubsets; iter++ )
   {
      const int start_curr = substarts[iter];
      const int start_prev = substarts[iter - 1];
      const int start_next = substarts[iter + 1];

      assert(start_prev < start_curr && start_curr < start_next);

      if( (start_curr - start_prev) < (start_next - start_curr) )
         continue;

      if( (start_curr - start_prev) > (start_next - start_curr) )
      {
         SCIPdebugMessage("wrongly ordered sizes \n");
         return FALSE;
      }

      if( memcmp(&partition[start_prev], &partition[start_curr], start_curr - start_prev) >= 0 )
      {
         SCIPdebugMessage("wrongly ordered substrings \n");
         return FALSE;
      }
   }

   return TRUE;
}
#endif


/** sorts all subsets
 * NOTE: partition needs to have allocated entry [-1] and [size] */
static inline
void partitionSortSubsets(
   DPB_Ptype* RESTRICT   partition,          /**< array of size 'size' */
   DPB_Ptype             delimiter,          /**< delimiter */
   int                   size                /**< size */
)
{
   int nsubsets = 0;
   int dummy[BPBORDER_MAXBORDERSIZE + 1];
   DPB_Ptype subbuffer[BPBORDER_MAXBORDERSIZE + 1];
   int* subsizes = &dummy[1];
   assert(size >= 1);

   /* sentinel */
   dummy[0] = 0;

   for( int iter = 0; iter < size; iter++ )
   {
      int iter2;
      for( iter2 = iter + 1; iter2 < size; iter2++ )
      {
         if( partition[iter2] == delimiter )
            break;
      }

      subsizes[nsubsets++] = (iter2 - iter) + 1;

      /* sentinel */
      for( int i = iter + 1; i < iter2; i++ )
      {
         int j;
         const DPB_Ptype curr = partition[i];

         for( j = i - 1; curr < partition[j] && j >= 0; j-- )
         {
            assert(j >= 0);
            partition[j + 1] = partition[j];
         }
         partition[j + 1] = curr;
      }

      iter = iter2;
   }

#ifndef NDEBUG
   for( int iter = 1; iter < size; iter++ )
   {
      if( partition[iter - 1] == delimiter )
         continue;

      assert(partition[iter - 1] < partition[iter]);
   }
#endif

   assert(nsubsets >= 1);

   // todo extra method

   /* add delimiter so that each subset ends with a delimiter */
   partition[size] = delimiter;

   for( int i = 1, curr_pos = subsizes[0]; i < nsubsets; i++ )
   {
      const int curr_size = subsizes[i];
      int j;
      int j_pos;
      assert(curr_size > 0);

      memcpy(subbuffer, &partition[curr_pos], curr_size * sizeof(partition[0]));

      for( j = i - 1, j_pos = curr_pos - subsizes[i - 1];
         curr_size < subsizes[j] ||
        (curr_size == subsizes[j] && memcmp(&partition[curr_pos], &partition[j_pos], curr_size) < 0 );
         j-- )
      {
         assert(j >= 0 && j_pos >= 0);
         j_pos -= subsizes[j - 1];
         subsizes[j + 1] = subsizes[j];
      }
      j_pos += subsizes[j];
      memmove(&partition[j_pos + curr_size], &partition[j_pos], (curr_pos - j_pos) * sizeof(partition[0]));
      memcpy(&partition[j_pos], subbuffer, curr_size * sizeof(partition[0]));
      subsizes[j + 1] = curr_size;
      curr_pos += curr_size;
   }

   assert(partition[size] == delimiter);
   assert(partitionIsSorted(partition, delimiter, size));
}


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

   assert(partsize > 0);

   if( partitionchars[0] == delimiter )
   {
      SCIPdebugMessage("partition starts with delimiter\n");
      return FALSE;
   }

   if( partitionchars[partsize - 1] == delimiter )
   {
      SCIPdebugMessage("partition ends with delimiter\n");
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

      if( LT(borderchardists[(unsigned char)borderchar], FARAWAY) )
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

         if( LT(borderchardists[(unsigned char)partchar], minedgecost) )
            minedgecost = borderchardists[(unsigned char)partchar];
      }

      costsum += minedgecost;

      if( GE(costsum, FARAWAY) )
         break;
   }
   assert(GE(costsum, 0.0));

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
#ifdef SCIP_DEBUG
   DPBPART partition;
#endif

   assert(dpborder_partIsValid(borderpartition));
   assert(globalstart + partsize + 2 < dpborder->global_partcap);
   assert(globalstart >= 1);

   /* form the union of marked subsets, as well as of extension node (if in border) */
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

         if( bordercharmap[(unsigned char)partchar] != -1 )
            global_partitions[globalend++] = bordercharmap[(unsigned char)partchar];
      }

      assert(partitionchars[candstart] < delimiter_prev);

      /* we mark the starts to skip them later on */
      partitionchars[candstart] = -partitionchars[candstart] - 1;
      assert(partitionchars[candstart] < 0);
   }

   if( dpborder->extborderchar >= 0 )
   {
      assert(dpborder_getTopLevel(dpborder)->extnode == dpborder->bordernodes[(unsigned char)dpborder->extborderchar]);
      global_partitions[globalend++] = dpborder->extborderchar;
   }

   if( globalend == globalstart )
   {
      SCIPdebugMessage("...empty first subset... \n");
      for( int j = 0; j < partsize; j++ )
      {
         if( partitionchars[j] < 0 )
            partitionchars[j] = -(partitionchars[j] + 1);
      }

#ifndef NDEBUG
      for( int j = 0; j < partsize; j++ )
         assert(0 <= partitionchars[j] && partitionchars[j] <= delimiter_prev);
#endif

      return -1;
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

      if( bordercharmap[(unsigned char)partchar] != -1 )
         global_partitions[globalend++] = bordercharmap[(unsigned char)partchar];
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

#ifndef NDEBUG
   for( int j = 0; j < partsize; j++ )
      assert(0 <= partitionchars[j] && partitionchars[j] <= delimiter_prev);
#endif

   if( globalstart != -1 )
   {
      int position;
      const int partition_length = globalend - globalstart;

#ifdef SCIP_DEBUG
      partition.partchars = &(global_partitions[globalstart]);
      partition.partsize = (globalend - globalstart);
      partition.delimiter = delimiter_new;
      printf("new (sub) partition \n");
      dpborder_partPrint(&partition);
#endif

      partitionSortSubsets(&global_partitions[globalstart], delimiter_new, partition_length);

#ifdef SCIP_DEBUG
      printf("sorted: \n");
      dpborder_partPrint(&partition);
#endif

      position = hashmap_get(&dpborder->hashmap, globalstart, partition_length);

      /* not found? */
      if( -1 == position )
      {
         SCIPdebugMessage("...partition is new \n");
         StpVecPushBack(scip, dpborder->global_partstarts, globalend);
         StpVecPushBack(scip, dpborder->global_partcosts, FARAWAY);
         StpVecPushBack(scip, dpborder->global_partsUseExt, TRUE);
         StpVecPushBack(scip, dpborder->global_predparts, -1);
         position = dpborder->global_npartitions;
         dpborder->global_npartitions++;

         assert(!partitionIsIncluded(dpborder, &global_partitions[globalstart], partition_length));
         hashmap_put(&dpborder->hashmap, globalstart, partition_length, position);
      }
      else
      {
         assert(LT(dpborder->global_partcosts[position], FARAWAY));
         assert(0 == memcmp(&global_partitions[globalstart], &global_partitions[dpborder->global_partstarts[position]], partition_length));
      }

#ifdef SCIP_DEBUG
      printf("final new (sub) partition glbpos=%d \n", position);
#endif

      assert(1 <= position && position < dpborder->global_npartitions);
      assert(dpborder_getTopLevel(dpborder)->globalstartidx <= position);

      return position;
   }

   SCIPdebugMessage("invalid partition... \n");

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
   int position;
   int partition_length;
   int globalstart = dpborder->global_partstarts[dpborder->global_npartitions];
   int globalend = globalstart;
   DPB_Ptype* RESTRICT global_partitions = dpborder->global_partitions;
   const int* const bordercharmap = dpborder->bordercharmap;
   const DPB_Ptype* const partitionchars = borderpartition->partchars;
   const DPB_Ptype delimiter_new = dpborder_getTopDelimiter(dpborder);
   const int partsize = borderpartition->partsize;
#ifdef SCIP_DEBUG
   DPBPART partition;
#endif

   assert(dpborder_partIsValid(borderpartition));
   assert(globalstart + partsize < dpborder->global_partcap);
   assert(globalstart >= 1);

   for( int i = 0; i < partsize; i++ )
   {
      const DPB_Ptype partchar = partitionchars[i];
      assert(0 <= partchar && partchar <= borderpartition->delimiter);
      assert(partchar != borderpartition->delimiter || bordercharmap[(unsigned char)partchar] == delimiter_new);

      if( bordercharmap[(unsigned char)partchar] != -1 )
         global_partitions[globalend++] = bordercharmap[(unsigned char)partchar];
   }

   if( globalstart == globalend
      || global_partitions[globalstart] == delimiter_new
      || global_partitions[globalend - 1] == delimiter_new )
   {
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

   partition_length = (globalend - globalstart);

#ifdef SCIP_DEBUG
   partition.partchars = &(global_partitions[globalstart]);
   partition.partsize = partition_length;
   partition.delimiter = delimiter_new;
   printf("new (exclusive sub) partition \n");
   dpborder_partPrint(&partition);
#endif

   partitionSortSubsets(&global_partitions[globalstart], delimiter_new, partition_length);
#ifdef SCIP_DEBUG
   printf("sorted: \n");
   dpborder_partPrint(&partition);
#endif

   position = hashmap_get(&dpborder->hashmap, globalstart, partition_length);

   /* not found? */
   if( -1 == position )
   {
      SCIPdebugMessage("...partition is new \n");
      StpVecPushBack(scip, dpborder->global_partstarts, globalend);
      StpVecPushBack(scip, dpborder->global_partcosts, FARAWAY);
      StpVecPushBack(scip, dpborder->global_predparts, -1);
      StpVecPushBack(scip, dpborder->global_partsUseExt, FALSE);
      position = dpborder->global_npartitions;
      dpborder->global_npartitions++;

      assert(!partitionIsIncluded(dpborder, &global_partitions[globalstart], partition_length));
      hashmap_put(&dpborder->hashmap, globalstart, partition_length, position);
   }
   else
   {
      assert(LT(dpborder->global_partcosts[position], FARAWAY));
      assert(0 == memcmp(&global_partitions[globalstart], &global_partitions[dpborder->global_partstarts[position]], partition_length));
   }

#ifdef SCIP_DEBUG
      printf("final new (sub) partition glbpos=%d \n", position);
#endif
   assert(1 <= position && position < dpborder->global_npartitions);
   assert(dpborder_getTopLevel(dpborder)->globalstartidx <= position);
   assert(dpborder->global_npartitions == StpVecGetSize(dpborder->global_predparts));
   assert(dpborder->global_npartitions == StpVecGetSize(dpborder->global_partcosts));
   assert(dpborder->global_npartitions == StpVecGetSize(dpborder->global_partsUseExt));
   assert(dpborder->global_npartitions + 1 == StpVecGetSize(dpborder->global_partstarts));

   return position;
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

      if( dpborder->global_partsUseExt[pos] )
      {
         const int extnode = dpborder->borderlevels[level]->extnode;
         assert(0 <= extnode && extnode < nnodes);
         SCIPdebugMessage("solnode=%d (ext) \n", extnode);

         nodes_isSol[extnode] = TRUE;
      }

      for( int i = global_partstarts[pos]; i != globalend; i++ )
      {
         const DPB_Ptype borderchar = global_partitions[i];
         assert(0 <= borderchar && borderchar <= delimiter);

         if( borderchar != delimiter )
         {
            const int node = nodemap[(unsigned char)borderchar];
            SCIPdebugMessage("solnode=%d \n", node);
            assert(0 <= node && node < nnodes);

            nodes_isSol[node] = TRUE;
         }
      }
   }

   SCIPdebugMessage("final solnode=%d \n", dpborder->borderlevels[0]->extnode);
   nodes_isSol[dpborder->borderlevels[0]->extnode] = TRUE;
}
