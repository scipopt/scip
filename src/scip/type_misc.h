/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_misc.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for miscellaneous datastructures
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_MISC_H__
#define __SCIP_TYPE_MISC_H__

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Sparse solution data structure
 *
 *  - \ref SparseSol "List of all available methods"
 */
typedef struct SCIP_SparseSol SCIP_SPARSESOL;

/** Priority queue data structure
 *
 *  - \ref PriorityQueue "List of all available methods"
 */
typedef struct SCIP_PQueue SCIP_PQUEUE;

/** Hash table data structure
 *
 *  - \ref HashTable "List of all available methods"
 */
typedef struct SCIP_HashTable SCIP_HASHTABLE;

/** Hash table element list to store single elements of a hash table */
typedef struct SCIP_HashTableList SCIP_HASHTABLELIST;

/** Hash map data structure
 *
 *  - \ref HashMap "List of all available methods"
 */
typedef struct SCIP_HashMap SCIP_HASHMAP;

/** Hash map  element list to store single mappings of a hash map */
typedef struct SCIP_HashMapList SCIP_HASHMAPLIST;

/** dynamic array for storing SCIP_Real values */
typedef struct SCIP_RealArray SCIP_REALARRAY;

/** dynamic array for storing int values */
typedef struct SCIP_IntArray SCIP_INTARRAY;

/** dynamic array for storing SCIP_Bool values */
typedef struct SCIP_BoolArray SCIP_BOOLARRAY;

/** dynamic array for storing pointers */
typedef struct SCIP_PtrArray SCIP_PTRARRAY;

/** Resource profile data structure
 *
 *  - \ref ResourceProfile "List of all available methods"
 */
typedef struct SCIP_Profile SCIP_PROFILE;

/** Directed graph data structure (stored as adjacency list)
 *
 *  - \ref DirectedGraph "List of all available methods"
 */
typedef struct SCIP_Digraph SCIP_DIGRAPH;

/** Binary tree data structure
 *
 *  - \ref BinaryTree "List of all available methods"
 */
typedef struct SCIP_Bt SCIP_BT;

/** search node of \ref SCIP_BT "binary tree" */
typedef struct SCIP_BtNode SCIP_BTNODE;

/** compares two element indices
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
#define SCIP_DECL_SORTINDCOMP(x) int x (void* dataptr, int ind1, int ind2)

/** compares two data element pointers
 *  result:
 *    < 0: elem1 comes before (is better than) elem2
 *    = 0: both elements have the same value
 *    > 0: elem2 comes after (is worse than) elem2
 */
#define SCIP_DECL_SORTPTRCOMP(x) int x (void* elem1, void* elem2)

/** gets the key of the given element */
#define SCIP_DECL_HASHGETKEY(x) void* x (void* userptr, void* elem)

/** returns TRUE iff both keys are equal */
#define SCIP_DECL_HASHKEYEQ(x) SCIP_Bool x (void* userptr, void* key1, void* key2)

/** returns the hash value of the key */
#define SCIP_DECL_HASHKEYVAL(x) unsigned int x (void* userptr, void* key)

#ifdef __cplusplus
}
#endif

#endif
