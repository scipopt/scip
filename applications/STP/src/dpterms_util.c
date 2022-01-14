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

/**@file   dpterms_util.c
 * @brief  Utility methods for dynamic programming solver for Steiner tree (sub-) problems
 * @author Daniel Rehfeldt
 *
 * Implements two implementations for finding valid intersections of sub-trees during DP.
 * One naive one, and one based on a search tree (DPS tree). Performance is slightly better for the latter one.
 * NOTE: DPS tree design is mostly taken from "Separator-Based Pruned Dynamic Programming for Steiner Tree"
 * by Iwata and Shigemura.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include "scip/scipdefplugins.h"
#include "dpterms.h"
#include "dptermsinterns.h"

#define CHILD_NONE -1
#define CHILD_LEFT  0
#define CHILD_RIGHT 1

/*
 * Data structures
 */

/** sub tree */
typedef struct dynamic_programming_search_tree_node
{
   int                   children_id[2];     /**< indices of children */
   STP_Bitset            terms_intersection; /**< intersection of all terminals of leaves of subtree */
   STP_Bitset            roots_union;        /**< union of all roots of leaves of subtree */
   int64_t               nsubsets;           /**< number of subsets */
   int                   split_pos;          /**< split position */
} DPSNODE;


/** search tree */
struct dynamic_programming_search_tree
{
   STP_Vectype(DPSNODE)  treenodes;          /**< nodes */
   STP_Vectype(int)      intersects_tmp;     /**< helper */
   int                   treeroot;           /**< tree root */
   int                   nterms;             /**< number of terminals of underlying graph */
   int                   nnodes;             /**< number of nodes of underlying graph */
};


/*
 * Local methods
 */

#ifndef NDEBUG

/** valid sizes? for debug checks */
static
SCIP_Bool bitsetsizesAreValid(
   STP_Bitset            termsmark,          /**< terminal mark */
   STP_Bitset            rootsmark,          /**< marks roots of extension trees */
   const DPSTREE*        dpstree             /**< to check for */
)
{
   assert(dpstree && termsmark && rootsmark);

   if( stpbitset_getCapacity(termsmark) != (((dpstree->nterms + 63) / 64) * 64) )
   {
      printf("termsmark size is wrong %d!=%d \n", stpbitset_getCapacity(termsmark), (((dpstree->nterms + 63) / 64) * 64));

      return FALSE;
   }

   if( stpbitset_getCapacity(rootsmark) != (((dpstree->nnodes + 63) / 64) * 64) )
   {
      printf("rootsmark size is wrong %d!=%d \n", stpbitset_getCapacity(rootsmark), (((dpstree->nnodes + 63) / 64) * 64));

      return FALSE;
   }

   return TRUE;
}


/** valid node? for debug checks */
static
SCIP_Bool treenodeIsInRange(
   int                   node_pos,           /**< position of node */
   const DPSTREE*        dpstree             /**< to check for */
)
{
   assert(dpstree);

   if( node_pos < 0 )
   {
      return FALSE;
   }

   if( node_pos >= StpVecGetSize(dpstree->treenodes) )
   {
      return FALSE;
   }


   return TRUE;
}
#endif


/** NOTE: only debug prints */
static
void printNodeDebugInfo(
   int                   node_pos,           /**< position of node */
   const DPSTREE*        dpstree             /**< to check for */
)
{
   assert(treenodeIsInRange(node_pos, dpstree));
#ifdef SCIP_DEBUG
   SCIPdebugMessage("node_pos=%d; split_pos=%d child1=%d, child2=%d \n",
         node_pos,
         dpstree->treenodes[node_pos].split_pos,
         dpstree->treenodes[node_pos].children_id[CHILD_LEFT],
         dpstree->treenodes[node_pos].children_id[CHILD_RIGHT]);

   SCIPdebugMessage("intersection/union: \n");
   stpbitset_print(dpstree->treenodes[node_pos].terms_intersection);
   stpbitset_print(dpstree->treenodes[node_pos].roots_union);
#endif
}



/** adds leaf
 * NOTE: takes ownership of termsmark and rootsmark*/
static
void addLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Bitset            termsmark,          /**< terminal mark */
   STP_Bitset            rootsmark,          /**< marks roots of extension trees */
   int64_t               nsubsets,           /**< number of subsets */
   int*                  leafpos,            /**< position of new leaf */
   DPSTREE*              dpstree             /**< to insert to */
)
{
   DPSNODE leaf = { .children_id = {CHILD_NONE, CHILD_NONE},
                    .terms_intersection = termsmark,
                    .roots_union = rootsmark,
                    .nsubsets = nsubsets,
                    .split_pos = dpstree->nterms /* mark as leaf */
                  };

   assert(bitsetsizesAreValid(termsmark, rootsmark, dpstree));

   *leafpos = StpVecGetSize(dpstree->treenodes);
   StpVecPushBack(scip, dpstree->treenodes, leaf);
}


/** inserts (recursively) */
static
void insertData(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Bitset            termsmark,          /**< terminal mark */
   STP_Bitset            rootsmark,          /**< marks roots of extension trees */
   int64_t               nsubsets,           /**< number of subsets */
   int                   node_pos,           /**< node position */
   int                   split_pos,          /**< current split position */
   int*                  subrootpos,         /**< position of new sub-root */
   DPSTREE*              dpstree             /**< to insert to */
)
{
   STP_Vectype(DPSNODE) tnodes = dpstree->treenodes;
   STP_Bitset terms_intersection;

   assert(0 <= split_pos && split_pos <= dpstree->nterms);
   assert(treenodeIsInRange(node_pos, dpstree));
   assert(bitsetsizesAreValid(termsmark, rootsmark, dpstree));

   terms_intersection = tnodes[node_pos].terms_intersection;
   assert(terms_intersection);

   /* find the correct position */
   while( tnodes[node_pos].split_pos != split_pos &&
      stpbitset_bitIsTrue(terms_intersection, split_pos) == stpbitset_bitIsTrue(termsmark, split_pos) )
   {
      split_pos++;
      assert(split_pos <= dpstree->nterms);
   }

   if( tnodes[node_pos].split_pos == split_pos )
   {
      /* split position is greater equal than split position of the current node */
      int root;
      const int down_direction =
         stpbitset_bitIsTrue(termsmark, split_pos) ? CHILD_RIGHT : CHILD_LEFT;
      const int child_position = tnodes[node_pos].children_id[down_direction];
      STP_Bitset roots_union = tnodes[node_pos].roots_union;

      assert(child_position >= 0 && child_position != node_pos);

      /* update node data */
      stpbitset_and(scip, termsmark, terms_intersection, terms_intersection);
      stpbitset_or(scip, rootsmark, roots_union, roots_union);

      SCIPdebugMessage("update internal node:\n ");
      printNodeDebugInfo(node_pos, dpstree);
      SCIPdebugMessage("continue with child %d \n", down_direction);

      /* continue from child */
      insertData(scip, termsmark, rootsmark, nsubsets, child_position, split_pos + 1,
            &(root), dpstree);
      dpstree->treenodes[node_pos].children_id[down_direction] = root;

      *subrootpos = node_pos;
   }
   else
   {
      /* split position is smaller than split position of the current node */
      int leaf_pos;
      STP_Bitset roots_union = tnodes[node_pos].roots_union;
      const SCIP_Bool leafIsRight = stpbitset_bitIsTrue(termsmark, split_pos);

      assert(stpbitset_bitIsTrue(terms_intersection, split_pos) != stpbitset_bitIsTrue(termsmark, split_pos));

      /* we add a new internal node with current node and new leaf as children */

      addLeaf(scip,
            stpbitset_newCopy(scip, termsmark),
            stpbitset_newCopy(scip, rootsmark),
            nsubsets, &leaf_pos, dpstree);

      stpbitset_and(scip, terms_intersection, termsmark, termsmark);
      stpbitset_or(scip, roots_union, rootsmark, rootsmark);

      assert(stpbitset_getPopcount(rootsmark) > 0);

      {
         DPSNODE parent = {
            .terms_intersection = termsmark,
            .roots_union = rootsmark,
            .nsubsets = -1, /* mark as internal node */
            .split_pos = split_pos };

         if( leafIsRight )
         {
            parent.children_id[CHILD_LEFT] = node_pos;
            parent.children_id[CHILD_RIGHT] = leaf_pos;
         }
         else
         {
            parent.children_id[CHILD_LEFT] = leaf_pos;
            parent.children_id[CHILD_RIGHT] = node_pos;
         }

         *subrootpos = StpVecGetSize(dpstree->treenodes);
         StpVecPushBack(scip, dpstree->treenodes, parent);

         SCIPdebugMessage("created new internal node: \n");
         printNodeDebugInfo(*subrootpos, dpstree);
      }
   }
}


/** gets indices of intersections */
static
void streeCollectIntersects(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Bitset            termsmark,          /**< terminal mark */
   STP_Bitset            rootsmark,          /**< marks roots of extension trees */
   int                   node_pos,           /**< node position */
   DPSTREE*              dpstree             /**< to insert to */
)
{
   STP_Vectype(DPSNODE) tnodes = dpstree->treenodes;
   STP_Bitset terms_intersection = tnodes[node_pos].terms_intersection;
   STP_Bitset roots_union = tnodes[node_pos].roots_union;

   assert(scip && dpstree && termsmark && rootsmark && tnodes);
   assert(bitsetsizesAreValid(termsmark, rootsmark, dpstree));
   assert(treenodeIsInRange(node_pos, dpstree));

   SCIPdebugMessage("at node: \n");
   printNodeDebugInfo(node_pos, dpstree);

   if( stpbitset_haveIntersection(terms_intersection, termsmark) )
   {
      SCIPdebugMessage("common terminal with subtree, returning \n");
      return;
   }
   else if( !stpbitset_haveIntersection(roots_union, rootsmark) )
   {
      SCIPdebugMessage("no common extension-root with subtree, returning \n");
      return;
   }
   /* at leaf? */
   else if( tnodes[node_pos].split_pos == dpstree->nterms )
   {
      SCIPdebugMessage("...is leaf \n");
      assert(tnodes[node_pos].nsubsets >= 0); /* double check that really leaf */

      StpVecPushBack(scip, dpstree->intersects_tmp, tnodes[node_pos].nsubsets);
   }
   else
   {
      SCIPdebugMessage("check left child \n");
      streeCollectIntersects(scip, termsmark, rootsmark,
            tnodes[node_pos].children_id[CHILD_LEFT],
            dpstree);

      /* any chance to find something in the right child? */
      if( !stpbitset_bitIsTrue(termsmark, tnodes[node_pos].split_pos) )
      {
         SCIPdebugMessage("check right child \n");
         streeCollectIntersects(scip, termsmark, rootsmark,
               tnodes[node_pos].children_id[CHILD_RIGHT],
               dpstree);
      }
   }
}


/*
 * Interface methods
 */


/** for debugging: checks whether given intersection is equal to naively computed one */
SCIP_Bool dpterms_intersectsEqualNaive(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Bitset            termsmark,          /**< terminal mark */
   STP_Bitset            rootsmark,          /**< marks roots of extension trees */
   STP_Vectype(int)      intersects,         /**< intersection indices */
   DPMISC*               dpmisc              /**< MISC DP data */
)
{
   STP_Vectype(int) intersects_naive = dpterms_collectIntersectsNaive(scip, termsmark, rootsmark, dpmisc);
   SCIP_Bool isEqual = TRUE;

   if( StpVecGetSize(intersects_naive) != StpVecGetSize(intersects) )
   {
#ifdef SCIP_DEBUG
      SCIPdebugMessage("wrong sizes: %d %d\n", StpVecGetSize(intersects_naive), StpVecGetSize(intersects));
      printf("naive: \n");

      for( int i = 0; i < StpVecGetSize(intersects_naive); i++ )
         printf("%d \n", intersects_naive[i]);

      printf("org: \n");
      for( int i = 0; i < StpVecGetSize(intersects); i++ )
         printf("%d \n", intersects[i]);
#endif
      isEqual = FALSE;
   }
   else
   {
      const int ninds = StpVecGetSize(intersects);
      int* intersects_org;

      SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &intersects_org, ninds) );
      BMScopyMemoryArray(intersects_org, intersects, ninds);

      SCIPsortDownInt(intersects_org, ninds);
      SCIPsortDownInt(intersects_naive, ninds);

      for( int i = 0; i < ninds; i++ )
      {
         if( intersects_org[i] != intersects_naive[i] )
         {
            SCIPdebugMessage("wrong index %d: %d!=%d \n", i, intersects_org[i], intersects_naive[i]);

            isEqual = FALSE;
            break;
         }
      }

      SCIPfreeMemoryArray(scip, &intersects_org);
   }


   StpVecFree(scip, intersects_naive);

   return isEqual;
}


/** gets indices of intersections by using naive computation */
STP_Vectype(int) dpterms_collectIntersectsNaive(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Bitset            termsmark,          /**< terminal mark */
   STP_Bitset            rootsmark,          /**< marks roots of extension trees */
   DPMISC*               dpmisc              /**< MISC DP data */
)
{
   STP_Vectype(int) intersects = NULL;
   STP_Vectype(SOLTRACE) global_traces = dpmisc->global_traces;
   STP_Vectype(STP_Bitset) global_termbits = dpmisc->global_termbits;
   STP_Vectype(int) global_starts = dpmisc->global_starts;
   const int nsets = StpVecGetSize(global_termbits);

   for( int i = 0; i < nsets; i++ )
   {
      STP_Bitset termsmark_i = global_termbits[i];
      assert(stpbitset_setsAreCompatible(termsmark_i, termsmark));

      if( !stpbitset_haveIntersection(termsmark_i, termsmark) )
      {
         SOLTRACE* soltraces_i = &(global_traces[global_starts[i]]);
         const int nsoltraces = global_starts[i + 1] - global_starts[i];
         SCIP_Bool hasIntersection = FALSE;

         assert(nsoltraces >= 0);

         for( int j = 0; j < nsoltraces; j++ )
         {
            const int root = soltraces_i[j].root;

            if( stpbitset_bitIsTrue(rootsmark, root) )
            {
               hasIntersection = TRUE;
               break;
            }
         }

         if( hasIntersection )
         {
            StpVecPushBack(scip, intersects, i);
         }
      }
   }

   return intersects;
}


/** inserts
 *  NOTE: dps-tree takes ownership of bitsets! */
SCIP_RETCODE dpterms_streeInsert(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Bitset            termsmark,          /**< terminal mark; will become OWNED */
   STP_Bitset            rootsmark,          /**< marks roots of extension trees; will become OWNED */
   int64_t               nsubsets,           /**< number of subsets */
   DPSTREE*              dpstree             /**< to insert to */
)
{
   assert(scip && dpstree && termsmark && rootsmark);
   assert(nsubsets >= 0);
   assert(bitsetsizesAreValid(termsmark, rootsmark, dpstree));
   assert(stpbitset_getPopcount(rootsmark) > 0);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("inserting: \n");
   SCIPdebugMessage("termsmark: ");
   stpbitset_print(termsmark);
   SCIPdebugMessage("rootsmark: ");
   stpbitset_print(rootsmark);
#endif

   // todo if ever too big, need to change some data to int64_t
   if( nsubsets >= (INT_MAX - 1) )
   {
      SCIPerrorMessage("too many subsets in terminal-based DP \n");
      return SCIP_ERROR;
   }

   if( dpstree->treeroot == -1 )
   {
      int leafpos;
      SCIPdebugMessage("...is first element \n");

      addLeaf(scip, termsmark, rootsmark, nsubsets, &leafpos, dpstree);
      assert(leafpos == 0);
      dpstree->treeroot = 0;
   }
   else
   {
      int subroot;
      const int node_pos = dpstree->treeroot;
      const int split_pos = 0;
      SCIPdebugMessage("insert from root %d \n", dpstree->treeroot);

      insertData(scip, termsmark, rootsmark, nsubsets, node_pos, split_pos, &subroot, dpstree);
      assert(subroot >= 0);
      SCIPdebugMessage("new root=%d \n", subroot);
      dpstree->treeroot = subroot;
   }

   return SCIP_OKAY;
}


/** gets indices of intersections */
STP_Vectype(int) dpterms_streeCollectIntersects(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Bitset            termsmark,          /**< terminal mark */
   STP_Bitset            rootsmark,          /**< marks roots of extension trees */
   DPSTREE*              dpstree             /**< to insert to */
)
{
   STP_Vectype(int) intersects;

   assert(scip && dpstree && termsmark && rootsmark);
   assert(bitsetsizesAreValid(termsmark, rootsmark, dpstree));
   assert(!dpstree->intersects_tmp);

   SCIPdebugMessage("collect intersections \n");

   if( dpstree->treeroot != -1 )
   {
      assert(dpstree->treeroot >= 0);
      streeCollectIntersects(scip, termsmark, rootsmark, dpstree->treeroot, dpstree);
   }

   intersects = dpstree->intersects_tmp;
   dpstree->intersects_tmp = NULL;

   return intersects;
}


/** initializes */
SCIP_RETCODE dpterms_streeInit(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nterms,             /**< number of terminals */
   int                   nnodes,             /**< number of nodes */
   DPSTREE**             dpstree             /**< initialize */
)
{
   DPSTREE* tree;

   assert(scip);
   assert(nterms > 0);
   assert(nnodes >= nterms);

   SCIP_CALL( SCIPallocMemory(scip, dpstree) );
   tree = *dpstree;

   tree->treeroot = -1;
   tree->nterms = nterms;
   tree->nnodes = nnodes;
   tree->treenodes = NULL;
   tree->intersects_tmp = NULL;

   return SCIP_OKAY;
}


/** frees */
void dpterms_streeFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSTREE**             dpstree             /**< initialize */
)
{
   DPSTREE* tree;
   int ntreenodes;

   assert(scip);
   assert(!(*dpstree)->intersects_tmp);

   tree = *dpstree;
   ntreenodes = StpVecGetSize(tree->treenodes);

   for( int i = ntreenodes - 1; i >= 0; i-- )
   {
      stpbitset_free(scip, &(tree->treenodes[i].terms_intersection));
      stpbitset_free(scip, &(tree->treenodes[i].roots_union));
   }

   StpVecFree(scip, tree->treenodes);

   SCIPfreeMemory(scip, dpstree);
}
