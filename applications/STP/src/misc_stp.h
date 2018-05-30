/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc_stp.h
 * @brief  miscellaneous methods used for solving Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file includes miscellaneous methods used for solving Steiner problems. For more details see \ref MISCSTP page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MISC_STP_H__
#define __SCIP_MISC_STP_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns maximum of given SCIP_Real values */
extern
SCIP_Real misc_stp_maxReal(
   SCIP_Real*            realarr,            /**< array of reals */
   unsigned              nreals              /**< size of array of reals */
  );

/** graph node structure storing number and distance */
typedef struct Graph_Node
{
   int                   number;             /**< node number   */
   SCIP_Real             dist;               /**< node distance */
}GNODE;

/** voronoi list node structure storing distance, incoming edge,base and pointer to next list node */
typedef struct Vnoi_List_Node
{
   double                dist;               /**< Distance to the end of the path             */
   signed int            edge;               /**< First edge to go                            */
   signed int            base;               /**< Voronoi base                                */
   struct Vnoi_List_Node *next;
}VLIST;

/** node */
typedef struct ST_Node
{
   int                   edge;               /**< edge to the node       */
   struct ST_Node        *parent;            /**< pointer to parent node */
}NODE;


/** a  weighted-quick-union-path-compression union find structure */
typedef struct UnionFind_Structure
{
   int*                  parent;             /**< parent[i] stores the parent of i                       */
   int*                  size;               /**< size[i] stores number of nodes in the tree rooted at i */
   int                   count;              /**< number of components                                   */
}UF;

/** integer list node */
typedef struct Int_List_Node
{
   int                   index;              /**< int value to store     */
   struct Int_List_Node  *parent;            /**< pointer to parent node */
}IDX;

/** Pairing heap node */
typedef struct PHeap_Node
{
   int                   element;            /**< int data value        */
   SCIP_Real             key;                /**< key value        */
   struct PHeap_Node*    child;              /**< pointer to child node */
   struct PHeap_Node*    sibling;            /**< pointer to right sibling */
   struct PHeap_Node*    prev;               /**< pointer to to previous node */
}PHNODE;


/**
 * Int List
 */

/** append copy of list pertaining to node2 to node1 */
extern
SCIP_RETCODE SCIPintListNodeAppendCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX**                 node1,              /**< pointer to the last node of list to be enlarged */
   IDX*                  node2,              /**< pointer to the last node of source list */
   SCIP_Bool*            conflict            /**< pointer to store whether a conflict has been detected by the method */
   );

/** insert a new node */
extern
SCIP_RETCODE SCIPintListNodeInsert(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX**                 node,               /**< pointer to the last list node */
   int                   nodeval             /**< data of the new node */
   );

/** free list */
extern
void SCIPintListNodeFree(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX**                 node                /**< pointer to the last list node */
   );

/** compares distances of two GNODE structures */
extern
int GNODECmpByDist(
   void                  *first_arg,         /**< first argument */
   void                  *second_arg         /**< second argument */
   );

/**
 * Linear Link Cut Tree
 */

/** inits a node, setting 'parent' and 'edge' to its default values */
extern
void SCIPlinkcuttreeInit(
   NODE*                 v                   /**< pointer to node representing the tree */
   );

/** renders w a child of v; v has to be the root of its tree */
extern
void SCIPlinkcuttreeLink(
   NODE*                 v,                  /**< pointer to node representing the tree */
   NODE*                 w,                  /**< pointer to the child */
   int                   edge                /**< link edge */
   );

/** cut tree at given node */
extern
void SCIPlinkcuttreeCut(
   NODE*                 v                   /**< node to cut at */
   );

/** finds minimum weight chain between node 'start' and distinct root node **/
SCIP_Real SCIPlinkcuttreeFindMinChain(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      nodeweight,         /**< node weight array */
   const int*            head,               /**< head of an arc */
   const int*            stdeg,              /**< degree in Steiner tree */
   const NODE*           start,              /**< the node to start at */
   NODE**                first,              /**< first node of chain */
   NODE**                last                /**< last node of chain */
   );

/** finds the max value between node 'v' and the root of the tree **/
extern
NODE* SCIPlinkcuttreeFindMax(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      cost,               /**< edge cost array */
   NODE*                 v                   /**< the node */
   );

/** makes vertex v the root of the link cut tree */
extern
void SCIPlinkcuttreeEvert(
   NODE*                 v                   /**< the vertex to become the root */
   );


/*
 * Pairing Heap
 */

/** links nodes 'root1' and 'root2' together */
extern
PHNODE* SCIPpairheapMergeheaps(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE                *root1,             /**< pointer to root of first heap */
   PHNODE                *root2              /**< pointer to root of second heap */
   );

extern
PHNODE* SCIPpairheapAddtoheap(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE*               root1,              /**< pointer to root of first heap */
   PHNODE*               root2               /**< pointer to root of second heap */
   );

/** inserts a new node into the pairing heap */
extern
SCIP_RETCODE SCIPpairheapInsert(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              root,
   int                   element,
   SCIP_Real             key,
   int*                  size
   );

/** deletes the root of the paring heap, concomitantly storing its data and key in '*element' and '*key' respectively */
extern
SCIP_RETCODE SCIPpairheapDeletemin(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  element,            /**< data of the root */
   SCIP_Real*            key,                /**< key of the root */
   PHNODE**              root,               /**< pointer to root of the heap */
   int*                  size                /**< pointer to size of the heap */
   );

/** links nodes 'root1' and 'root2' together, roots the resulting tree at root1 and sets root2 to NULL */
extern
void SCIPpairheapMeldheaps(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              root1,              /**< pointer to root of first heap */
   PHNODE**              root2,              /**< pointer to root of second heap */
   int*                  sizeroot1,          /**< pointer to size of first heap */
   int*                  sizeroot2           /**< pointer to size of second heap */
   );

/** frees the paring heap with root 'p' */
extern
void SCIPpairheapFree(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              root                /**< root of heap to be freed */
   );

/** stores all elements of the pairing heap in an array */
extern
SCIP_RETCODE SCIPpairheapBuffarr(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE*               root,               /**< root of the heap */
   int                   size,               /**< size of the array */
   int**                 elements            /**< pointer to array */
   );

/*
 * Union-Find data structure
 */

/** initializes the union-find structure 'uf' with 'length' many components (of size one) */
extern
SCIP_RETCODE SCIPStpunionfindInit(
   SCIP*                 scip,               /**< SCIP data structure */
   UF*                   uf,                 /**< union find data structure */
   int                   length              /**< number of components */
   );

/** clears the union-find structure 'uf'*/
extern
void SCIPStpunionfindClear(
   SCIP*                 scip,               /**< SCIP data structure */
   UF*                   uf,                 /**< union find data structure */
   int                   length              /**< number of components */
   );

/** finds and returns the component identifier */
extern
int SCIPStpunionfindFind(
   UF*                   uf,                 /**< union find data structure */
   int                   element             /**< element to be found */
   );

/** merges the components containing p and q respectively */
extern
void SCIPStpunionfindUnion(
   UF*                   uf,                 /**< union find data structure */
   int                   p,                  /**< first component */
   int                   q,                  /**< second component*/
   SCIP_Bool             compress            /**< compress union find structure? */
   );

/** frees the data fields of the union-find structure */
extern
void SCIPStpunionfindFree(
   SCIP*                 scip,               /**< SCIP data structure */
   UF*                   uf                  /**< union find data structure */
   );


#ifdef __cplusplus
}
#endif

#endif
