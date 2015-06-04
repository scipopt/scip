/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc_stp.h
 * @brief  miscellaneous data structures and methods
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MISC_STP_H__
#define __SCIP_MISC_STP_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Graph_Node
{
   int number;
   SCIP_Real dist;
}GNODE;

/* a  weighted-quick-union-path-compression union find structure */
typedef struct UnionFind_Structure
{
   int* parent;    /* parent[i] stores the parent of i */
   int* size;      /* size[i] stores number of nodes in the tree rooted at i */
   int count;      /* number of components */
}UF;


typedef struct Vnoi_List_Node
{
   double       dist;         /* Distance to the end of the path             */
   signed int   edge;         /* First edge to go                            */
   signed int   base;         /* Voronoi base                            */
   struct Vnoi_List_Node *next;
} VLIST;

typedef struct ST_Node
{
   int edge;
   struct ST_Node *parent;
} NODE;

typedef struct Index_List_Node
{
   int index;
   struct Index_List_Node *parent;
} IDX;


typedef struct PairingHeap_Node
{
   int element;
   SCIP_Real key;
   struct PairingHeap_Node* child;
   struct PairingHeap_Node* sibling;
   struct PairingHeap_Node* prev;

}PHNODE;

extern
SCIP_RETCODE SCIPindexListNodeAppendCopy(
   SCIP* scip,
   IDX** node1,
   IDX* node2
   );

extern
SCIP_RETCODE SCIPindexListNodeInsert(
   SCIP* scip,
   IDX** node,
   int   index
   );

extern
void SCIPindexListNodeFree(
   SCIP* scip,
   IDX** node
   );

/***  ***/
extern
int GNODECmpByDist(
   void *first_arg,
   void *second_arg
   );

/**
 * Linear Link Cut Tree
 */

/** inits a node, setting 'parent' and 'edge' to its default values */
extern
void SCIPlinkcuttreeInit(
   NODE* v
   );

/** renders w a child of v; v has to be the root of its tree */
extern
void SCIPlinkcuttreeLink(
   NODE* v,
   NODE* w,
   int edge
   );

extern
void SCIPlinkcuttreeCut(
   NODE* v
   );

/** finds the max value between node 'v' and the root of the tree **/
extern
NODE* SCIPlinkcuttreeFindMax(
   SCIP* scip,
   const SCIP_Real* cost,
   NODE* v
   );

/** makes vertex v the root of the link cut tree */
extern
void SCIPlinkcuttreeEvert(
   NODE* v
   );


/*
 * Pairing Heap
 */

/** links nodes 'root1' and 'root2' together */
extern
PHNODE* SCIPpairheapMergeheaps(
   SCIP* scip,
   PHNODE *root1,
   PHNODE *root2
   );

extern
PHNODE* SCIPpairheapAddtoheap(
   SCIP* scip,
   PHNODE *root1,
   PHNODE *root2
   );

/** inserts a new node into the pairing heap */
extern
SCIP_RETCODE SCIPpairheapInsert(
   SCIP* scip,
   PHNODE** root,
   int element,
   SCIP_Real key,
   int* size
   );

/** deletes the root of the paring heap, concomitantly storing its data and key in '*element' and '*key' respectively */
extern
void SCIPpairheapDeletemin(
   SCIP* scip,
   int* element,
   SCIP_Real *key,
   PHNODE** root,
   int* size
   );

/** links nodes 'root1' and 'root2' together, roots the resulting tree at root1 and sets root2 to NULL */
extern
void SCIPpairheapMeldheaps(
   SCIP* scip,
   PHNODE** root1,
   PHNODE** root2,
   int* sizeroot1,
   int* sizeroot2
   );

/** frees the paring heap with root 'p' */
extern
void SCIPpairheapFree(
   SCIP* scip,
   PHNODE** root
   );

/** stores all elements of the pairing heap in an array */
extern
SCIP_RETCODE SCIPpairheapBuffarr(
   SCIP* scip,
   PHNODE* root,
   int size,
   int** elements
   );

/*
 * Union-Find data structure
 */

/** initializes the union-find structure 'uf' with 'length' many components (of size one) */
extern
SCIP_RETCODE SCIPunionfindInit(
   SCIP* scip,
   UF* uf,
   int length
   );

/** finds and returns the component identifier */
extern
int SCIPunionfindFind(
   UF* uf,
   int element
   );

/** merges the components containing p and q respectively */
extern
void SCIPunionfindUnion(
   UF* uf,
   int p,
   int q,
   SCIP_Bool compress
   );

/** frees the data fields of the union-find structure */
extern
void SCIPunionfindFree(
   SCIP* scip,
   UF* uf
   );


#ifdef __cplusplus
}
#endif

#endif
