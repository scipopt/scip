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

/**@file   struct_misc.h
 * @brief  miscellaneous datastructures
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_MISC_H__
#define __SCIP_STRUCT_MISC_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data structure for sparse solutions */
struct SCIP_SparseSol
{
   SCIP_VAR**            vars;               /**< variables */
   SCIP_Longint*         lbvalues;           /**< array of lower bounds */
   SCIP_Longint*         ubvalues;           /**< array of upper bounds */
   int                   nvars;              /**< number of variables */
};

/** priority queue data structure
 *  Elements are stored in an array, which grows dynamically in size as new elements are added to the queue.
 *  The ordering is done through a pointer comparison function.
 *  The array is organized as follows. The root element (that is the "best" element $r$ with $r <= x$ for all $x$)
 *  is stored in position 0. The children of an element at position $p$ are stored at positions $q_1 = 2*p+1$ and
 *  $q_2 = 2*p+2$. That means, the parent of the element at position $q$ is at position $p = (q-1)/2$.
 *  At any time, the condition holds that $p <= q$ for each parent $p$ and its children $q$.
 *  Insertion and removal of single elements needs time $O(log n)$.
 */
struct SCIP_PQueue
{
   SCIP_Real             sizefac;            /**< memory growing factor */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp));        /**< compares two data elements */
   void**                slots;              /**< array of element slots */
   int                   len;                /**< number of used element slots */
   int                   size;               /**< total number of available element slots */
};

/** element list to store single elements of a hash table */
struct SCIP_HashTableList
{
   void*                 element;            /**< this element */
   SCIP_HASHTABLELIST*   next;               /**< rest of the hash table list */
};

/** hash table data structure */
struct SCIP_HashTable
{
   SCIP_DECL_HASHGETKEY((*hashgetkey));      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq));       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval));      /**< returns the hash value of the key */
   BMS_BLKMEM*           blkmem;             /**< block memory used to store hash map entries */
   SCIP_HASHTABLELIST**  lists;              /**< hash table lists of the hash table */
   int                   nlists;             /**< number of lists stored in the hash table */
   void*                 userptr;            /**< user pointer */
};

/** element list to store single mappings of a hash map */
struct SCIP_HashMapList
{
   void*                 origin;             /**< origin of the mapping origin -> image */
   void*                 image;              /**< image of the mapping origin -> image */
   SCIP_HASHMAPLIST*     next;               /**< rest of the hash map list */
};

/** hash map data structure to map pointers on pointers */
struct SCIP_HashMap
{
   BMS_BLKMEM*           blkmem;             /**< block memory used to store hash map entries */
   SCIP_HASHMAPLIST**    lists;              /**< hash map lists of the hash map */
   int                   nlists;             /**< number of lists stored in the hash map */
};

/** dynamic array for storing real values */
struct SCIP_RealArray
{
   BMS_BLKMEM*           blkmem;             /**< block memory that stores the vals array */
   SCIP_Real*            vals;               /**< array values */
   int                   valssize;           /**< size of vals array */
   int                   firstidx;           /**< index of first element in vals array */
   int                   minusedidx;         /**< index of first non zero element in vals array */
   int                   maxusedidx;         /**< index of last non zero element in vals array */
};

/** dynamic array for storing int values */
struct SCIP_IntArray
{
   BMS_BLKMEM*           blkmem;             /**< block memory that stores the vals array */
   int*                  vals;               /**< array values */
   int                   valssize;           /**< size of vals array */
   int                   firstidx;           /**< index of first element in vals array */
   int                   minusedidx;         /**< index of first non zero element in vals array */
   int                   maxusedidx;         /**< index of last non zero element in vals array */
};

/** dynamic array for storing bool values */
struct SCIP_BoolArray
{
   BMS_BLKMEM*           blkmem;             /**< block memory that stores the vals array */
   SCIP_Bool*            vals;               /**< array values */
   int                   valssize;           /**< size of vals array */
   int                   firstidx;           /**< index of first element in vals array */
   int                   minusedidx;         /**< index of first non zero element in vals array */
   int                   maxusedidx;         /**< index of last non zero element in vals array */
};

/** dynamic array for storing pointers */
struct SCIP_PtrArray
{
   BMS_BLKMEM*           blkmem;             /**< block memory that stores the vals array */
   void**                vals;               /**< array values */
   int                   valssize;           /**< size of vals array */
   int                   firstidx;           /**< index of first element in vals array */
   int                   minusedidx;         /**< index of first non zero element in vals array */
   int                   maxusedidx;         /**< index of last non zero element in vals array */
};

/** resource profile */
struct SCIP_Profile
{
   int*                  timepoints;         /**< time point array */
   int*                  loads;              /**< array holding the load for each time point */
   int                   capacity;           /**< capacity of the resource profile */
   int                   ntimepoints;        /**< current number of entries */
   int                   arraysize;          /**< current array size */
};

/** digraph structure to store and handle graphs */
struct SCIP_Digraph
{
   int**                 successors;         /**< adjacency list: for each node (first dimension) list of all successors */
   void***               arcdatas;           /**< arc datas corresponding to the arcs to successors given by the successors array  */
   int*                  successorssize;     /**< sizes of the successor lists for the nodes */
   int*                  nsuccessors;        /**< number of successors stored in the adjacency lists of the nodes */
   int*                  components;         /**< array to store the node indices of the components, one component after the other */
   int*                  componentstarts;    /**< array to store the start indices of the components in the components array */
   int                   ncomponents;        /**< number of undirected components stored */
   int                   componentstartsize; /**< size of array componentstarts */
   int                   nnodes;             /**< number of nodes, nodes should be numbered from 0 to nnodes-1 */
};

/** binary node data structure for binary tree */
struct SCIP_BtNode
{
   SCIP_BTNODE*          parent;             /**< pointer to the parent node */
   SCIP_BTNODE*          left;               /**< pointer to the left child node */
   SCIP_BTNODE*          right;              /**< pointer to the right child node */
   void*                 dataptr;            /**< user pointer */
};

/** binary search tree data structure */
struct SCIP_Bt
{
   SCIP_BTNODE*          root;               /**< pointer to the dummy root node; root is left child */
   BMS_BLKMEM*           blkmem;             /**< block memory used to store tree nodes */
};

#ifdef __cplusplus
}
#endif

#endif
