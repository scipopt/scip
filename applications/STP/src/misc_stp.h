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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc_stp.h
 * @brief  miscellaneous methods used for solving Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file includes miscellaneous methods used for solving Steiner problems. For more details see \ref STP_MISC page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MISC_STP_H__
#define __SCIP_MISC_STP_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SWAP_INTS(int1, int2)  \
   do                                  \
   {                                   \
      const int _tmp_ = int1;           \
      int1 = int2;                     \
      int2 = _tmp_;                     \
   } while(0)


#define SWAP_REALS(real1, real2)  \
   do                                  \
   {                                   \
      const SCIP_Real _tmp_ = real1;     \
      real1 = real2;                     \
      real2 = _tmp_;                     \
   } while(0)


/** graph node structure storing number and distance */
typedef struct Graph_Node
{
   int                   number;             /**< node number   */
   SCIP_Real             dist;               /**< node distance */
} GNODE;


/** Voronoi list node structure storing distance, incoming edge,base and pointer to next list node */
typedef struct Vnoi_List_Node
{
   double                dist;               /**< Distance to the end of the path             */
   signed int            edge;               /**< First edge to go                            */
   signed int            base;               /**< Voronoi base                                */
   struct Vnoi_List_Node *next;
} VLIST;

/** link-cut_node */
typedef struct link_cut_node
{
   int                   edge;               /**< edge to the node       */
   int                   parent;             /**< index of parent node */
} LCNODE;


/** a  weighted-quick-union-path-compression union find structure */
typedef struct UnionFind_Structure
{
   int*                  parent;             /**< parent[i] stores the parent of i                       */
   int*                  size;               /**< size[i] stores number of nodes in the tree rooted at i */
   int                   nComponents;        /**< number of components                                   */
   int                   nElements;          /**< number of elements                                     */
} UF;


/** integer list node */
typedef struct Int_List_Node
{
   int                   index;              /**< int value to store     */
   struct Int_List_Node  *parent;            /**< pointer to parent node */
} IDX;


/** Pairing heap node */
typedef struct PHeap_Node
{
   struct PHeap_Node*    child;              /**< pointer to child node */
   struct PHeap_Node*    sibling;            /**< pointer to right sibling */
   struct PHeap_Node*    prev;               /**< pointer to to previous node */
   SCIP_Real             key;                /**< key value        */
   int                   element;            /**< integer data value        */
} PHNODE;


/** returns maximum of given SCIP_Real values */
SCIP_EXPORT
SCIP_Real miscstp_maxReal(
   const SCIP_Real*      realarr,            /**< array of reals */
   unsigned              nreals              /**< size of array of reals */
  );

/**
 * Int List
 */

/** append copy of list pertaining to node2 to node1 */
SCIP_EXPORT
SCIP_RETCODE SCIPintListNodeAppendCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX**                 node1,              /**< pointer to the last node of list to be enlarged */
   IDX*                  node2,              /**< pointer to the last node of source list */
   SCIP_Bool*            conflict            /**< pointer to store whether a conflict has been detected by the method */
   );

/** append list pertaining to node2 to (non-empty!) node1 */
SCIP_EXPORT
void SCIPintListNodeAppend(
   IDX*                  node1,              /**< pointer to the last node of non-empty list to be enlarged */
   IDX*                  node2               /**< pointer to the last node of source list */
   );

/** insert a new node */
SCIP_EXPORT
SCIP_RETCODE SCIPintListNodeInsert(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX**                 node,               /**< pointer to the last list node */
   int                   nodeval             /**< data of the new node */
   );

/** free list */
SCIP_EXPORT
void SCIPintListNodeFree(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX**                 node                /**< pointer to the last list node */
   );

/** compares distances of two GNODE structures */
SCIP_EXPORT
int GNODECmpByDist(
   void                  *first_arg,         /**< first argument */
   void                  *second_arg         /**< second argument */
   );

/**
 * Linear Link Cut Tree
 */

/** inits a node, setting 'parent' and 'edge' to its default values */
SCIP_EXPORT
void SCIPlinkcuttreeInitNode(
   LCNODE*               v                   /**< pointer to node representing the tree */
   );

/** renders w a child of v; v has to be the root of its tree */
SCIP_EXPORT
void SCIPlinkcuttreeLink(
   LCNODE*               tree,               /**< the tree */
   int                   v,                  /**< pointer to node representing the tree */
   int                   w,                  /**< pointer to node of another tree */
   int                   edge                /**< link edge */
   );

/** cut tree at given node */
SCIP_EXPORT
void SCIPlinkcuttreeCutNode(
   LCNODE*               v                   /**< node to cut at */
   );

/** finds minimum weight chain between node 'start' and distinct root node (for maximum-weight connected subgraph) **/
SCIP_EXPORT
SCIP_Real SCIPlinkcuttreeFindMinChainMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const LCNODE*         tree,               /**< tree */
   const SCIP_Real*      nodeweight,         /**< node weight array */
   const int*            heads,              /**< head of an arc */
   const int*            stdeg,              /**< degree in Steiner tree */
   const SCIP_Bool*      nodeIsBlocked,      /**< has node been blocked? */
   int                   start,              /**< the node to start at */
   int*                  first,              /**< first node of chain */
   int*                  last                /**< last node of chain */
   );


/** finds maximum cost chain between node 'start' and distinct root node **/
SCIP_EXPORT
SCIP_Real SCIPlinkcuttreeFindMaxChain(
   SCIP*                 scip,               /**< SCIP data structure */
   const LCNODE*         tree,               /**< tree */
   const SCIP_Real*      edgecosts,          /**< edge cost array */
   const SCIP_Real*      prizes,             /**< node weight array for PC/RPC */
   const int*            heads,              /**< head of an arc */
   const int*            nonTermDeg,         /**< degree in Steiner tree, or UNKNOWN if vertex is terminal */
   const SCIP_Bool*      nodeIsBlocked,      /**< has node been blocked? */
   int                   start,              /**< the node to start at (NOT the root!) */
   int*                  first,              /**< first node of chain */
   int*                  last                /**< last node of chain */
   );


/** finds the max value between node 'v' and the root of the tree **/
SCIP_EXPORT
int SCIPlinkcuttreeFindMax(
   const LCNODE*         tree,                /**< tree */
   const SCIP_Real*      cost,                /**< edge cost array */
   int                   v                    /**< the node */
   );


/** makes vertex v the root of the link cut tree */
SCIP_EXPORT
void SCIPlinkcuttreeEvert(
   LCNODE* RESTRICT      tree,                /**< tree */
   int                   root_new             /**< the vertex to become the root  */
   );


/*
 * Pairing Heap
 */

/** links nodes 'root1' and 'root2' together */
SCIP_EXPORT
PHNODE* SCIPpairheapMergeheaps(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE                *root1,             /**< pointer to root of first heap */
   PHNODE                *root2              /**< pointer to root of second heap */
   );

/** inserts a new node into the pairing heap */
SCIP_EXPORT
SCIP_RETCODE SCIPpairheapInsert(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              root,
   PHNODE*               pheapelems,         /**< data array */
   int                   element,
   SCIP_Real             key,
   int*                  size
   );

/** deletes the root of the paring heap, concomitantly storing its data and key in '*element' and '*key' respectively */
SCIP_EXPORT
SCIP_RETCODE SCIPpairheapDeletemin(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  element,            /**< data of the root */
   SCIP_Real*            key,                /**< key of the root */
   PHNODE**              root,               /**< pointer to root of the heap */
   int*                  size                /**< pointer to size of the heap */
   );

/** links nodes 'root1' and 'root2' together, roots the resulting tree at root1 and sets root2 to NULL */
SCIP_EXPORT
void SCIPpairheapMeldheaps(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              root1,              /**< pointer to root of first heap */
   PHNODE**              root2,              /**< pointer to root of second heap */
   int*                  sizeroot1,          /**< pointer to size of first heap */
   int*                  sizeroot2           /**< pointer to size of second heap */
   );

/** frees the paring heap with root 'p' */
SCIP_EXPORT
void SCIPpairheapFree(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              root                /**< root of heap to be freed */
   );

/** stores all elements of the pairing heap in an array */
SCIP_EXPORT
SCIP_RETCODE SCIPpairheapBuffarr(
   SCIP*                 scip,               /**< SCIP data structure */
   const PHNODE*         root,               /**< root of the heap */
   int                   size,               /**< size of the array */
   int**                 elements            /**< pointer to array (will be allocated) */
   );



/*
 * Union-Find data structure
 */

/** initializes the union-find structure 'uf' with 'length' many components (of size one) */
SCIP_EXPORT
SCIP_RETCODE SCIPStpunionfindInit(
   SCIP*                 scip,               /**< SCIP data structure */
   UF*                   uf,                 /**< union find data structure */
   int                   length              /**< number of components */
   );

/** clears the union-find structure 'uf'*/
SCIP_EXPORT
void SCIPStpunionfindClear(
   SCIP*                 scip,               /**< SCIP data structure */
   UF*                   uf                  /**< union find data structure */
   );


/** is the union-find structure 'uf' clear? */
SCIP_EXPORT
SCIP_Bool SCIPStpunionfindIsClear(
   SCIP*                 scip,               /**< SCIP data structure */
   const UF*             uf                  /**< union find data structure */
   );


/** finds and returns the component identifier */
SCIP_EXPORT
int SCIPStpunionfindFind(
   UF*                   uf,                 /**< union find data structure */
   int                   element             /**< element to be found */
   );

/** Merges the components containing p and q respectively.
 *  Identifier of 'p' will always be used if 'compress' is FALSE. */
SCIP_EXPORT
void SCIPStpunionfindUnion(
   UF*                   uf,                 /**< union find data structure */
   int                   p,                  /**< first component */
   int                   q,                  /**< second component*/
   SCIP_Bool             compress            /**< compress union find structure? */
   );

/** frees the data fields of the union-find structure */
SCIP_EXPORT
void SCIPStpunionfindFreeMembers(
   SCIP*                 scip,               /**< SCIP data structure */
   UF*                   uf                  /**< union find data structure */
   );


#ifdef __cplusplus
}
#endif

#endif
