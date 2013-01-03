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

/**@file   pub_misc.h
 * @brief  public data structures and miscellaneous methods
 * @author Tobias Achterberg
 *
 * This file contains a bunch of data structures and miscellaneous methods:
 *
 * - \ref DataStructures "Data structures"
 * - \ref MiscellaneousMethods "Miscellaneous Methods"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MISC_H__
#define __SCIP_PUB_MISC_H__



#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_misc.h"
#include "scip/type_message.h"
#include "scip/type_var.h"

/* in optimized mode some of the function are handled via defines, for that the structs are needed */
#ifdef NDEBUG
#include "scip/struct_misc.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * GML graphical printing methods
 * For a detailed format decription see http://docs.yworks.com/yfiles/doc/developers-guide/gml.html
 */

/**@defgroup GMLgraph GML graphical printing
 *
 * @{
 */

/** writes a node section to the given graph file */
EXTERN
void SCIPgmlWriteNode(
   FILE*                 file,               /**< file to write to */
   unsigned int          id,                 /**< id of the node */
   const char*           label,              /**< label of the node */
   const char*           nodetype,           /**< type of the node, or NULL */
   const char*           fillcolor,          /**< color of the node's interior, or NULL */
   const char*           bordercolor         /**< color of the node's border, or NULL */
   );

/** writes an edge section to the given graph file */
EXTERN
void SCIPgmlWriteEdge(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   );

/** writes an arc section to the given graph file */
EXTERN
void SCIPgmlWriteArc(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   );

/** writes the starting line to a GML graph file, does not open a file */
EXTERN
void SCIPgmlWriteOpening(
   FILE*                 file,               /**< file to write to */
   SCIP_Bool             directed            /**< is the graph directed */
   );

/** writes the ending lines to a GML graph file, does not close a file */
EXTERN
void SCIPgmlWriteClosing(
   FILE*                 file                /**< file to close */
   );

/**@} */


/** @defgroup DataStructures Data Structures
 *
 *  Below you find a list of available data structures
 *
 * @{
 */

/*
 * Sparse solution
 */

/**@defgroup SparseSol Sparse solution
 *
 * @{
 */

/** creates a sparse solution */
EXTERN
SCIP_RETCODE SCIPsparseSolCreate(
   SCIP_SPARSESOL**      sparsesol,          /**< pointer to store the created sparse solution */
   SCIP_VAR**            vars,               /**< variables in the sparse solution, must not contain continuous
					      *   variables
					      */
   int                   nvars,              /**< number of variables to store, size of the lower and upper bound
					      *   arrays
					      */
   SCIP_Bool             cleared             /**< should the lower and upper bound arrays be cleared (entries set to
					      *	  0)
					      */
   );

/** frees priority queue, but not the data elements themselves */
EXTERN
void SCIPsparseSolFree(
   SCIP_SPARSESOL**      sparsesol           /**< pointer to a sparse solution */
   );

/** returns the variables in the given sparse solution */
EXTERN
SCIP_VAR** SCIPsparseSolGetVars(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   );

/** returns the number of variables in the given sparse solution */
EXTERN
int SCIPsparseSolGetNVars(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   );

/** returns the the lower bound array for all variables for a given sparse solution */
EXTERN
SCIP_Longint* SCIPsparseSolGetLbs(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   );

/** returns the the upper bound array for all variables for a given sparse solution */
EXTERN
SCIP_Longint* SCIPsparseSolGetUbs(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   );


/**@} */

/*
 * Priority Queue
 */

/**@defgroup PriorityQueue Priority Queue
 *
 * @{
 */

/** creates priority queue */
EXTERN
SCIP_RETCODE SCIPpqueueCreate(
   SCIP_PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** frees priority queue, but not the data elements themselves */
EXTERN
void SCIPpqueueFree(
   SCIP_PQUEUE**         pqueue              /**< pointer to a priority queue */
   );

/** clears the priority queue, but doesn't free the data elements themselves */
EXTERN
void SCIPpqueueClear(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** inserts element into priority queue */
EXTERN
SCIP_RETCODE SCIPpqueueInsert(
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   void*                 elem                /**< element to be inserted */
   );

/** removes and returns best element from the priority queue */
EXTERN
void* SCIPpqueueRemove(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the best element of the queue without removing it */
EXTERN
void* SCIPpqueueFirst(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the number of elements in the queue */
EXTERN
int SCIPpqueueNElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
EXTERN
void** SCIPpqueueElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/**@} */


/*
 * Hash Table
 */

/**@defgroup HashTable Hash Table
 *
 *@{
 */

/** returns a reasonable hash table size (a prime number) that is at least as large as the specified value */
EXTERN
int SCIPcalcHashtableSize(
   int                   minsize             /**< minimal size of the hash table */
   );

/** creates a hash table */
EXTERN
SCIP_RETCODE SCIPhashtableCreate(
   SCIP_HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash table entries */
   int                   tablesize,          /**< size of the hash table */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr             /**< user pointer */
   );

/** frees the hash table */
EXTERN
void SCIPhashtableFree(
   SCIP_HASHTABLE**      hashtable           /**< pointer to the hash table */
   );

/** removes all elements of the hash table
 *
 *  @note From a performance point of view you should not fill and clear a hash table too often since the clearing can
 *        be expensive. Clearing is done by looping over all buckets and removing the hash table lists one-by-one.
 */
EXTERN
void SCIPhashtableClear(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   );

/** inserts element in hash table (multiple inserts of same element possible) */
EXTERN
SCIP_RETCODE SCIPhashtableInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   );

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
EXTERN
SCIP_RETCODE SCIPhashtableSafeInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   );

/** retrieve element with key from hash table, returns NULL if not existing */
EXTERN
void* SCIPhashtableRetrieve(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 key                 /**< key to retrieve */
   );

/** retrieve element with key from hash table, returns NULL if not existing
 * can be used to retrieve all entries with the same key (one-by-one) */
EXTERN
void* SCIPhashtableRetrieveNext(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_HASHTABLELIST**  hashtablelist,      /**< input: entry in hash table list from which to start searching, or NULL; output: entry in hash table list corresponding to element after retrieved one, or NULL */
   void*                 key                 /**< key to retrieve */
   );

/** returns whether the given element exists in the table */
EXTERN
SCIP_Bool SCIPhashtableExists(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to search in the table */
   );

/** removes element from the hash table, if it exists */
EXTERN
SCIP_RETCODE SCIPhashtableRemove(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to remove from the table */
   );

/** prints statistics about hash table usage */
EXTERN
void SCIPhashtablePrintStatistics(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** standard hash key comparator for string keys */
EXTERN
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqString);

/** standard hashing function for string keys */
EXTERN
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValString);

/** gets the element as the key */
EXTERN
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyStandard);

/** returns TRUE iff both keys(pointer) are equal */
EXTERN
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqPtr);

/** returns the hash value of the key */
EXTERN
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValPtr);

/**@} */


/*
 * Hash Map
 */

/**@defgroup HashMap Hash Map
 *
 *@{
 */

/** creates a hash map mapping pointers to pointers */
EXTERN
SCIP_RETCODE SCIPhashmapCreate(
   SCIP_HASHMAP**        hashmap,            /**< pointer to store the created hash map */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash map entries */
   int                   mapsize             /**< size of the hash map */
   );

/** frees the hash map */
EXTERN
void SCIPhashmapFree(
   SCIP_HASHMAP**        hashmap             /**< pointer to the hash map */
   );

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
EXTERN
SCIP_RETCODE SCIPhashmapInsert(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   );

/** retrieves image of given origin from the hash map, or NULL if no image exists */
EXTERN
void* SCIPhashmapGetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   );

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
EXTERN
SCIP_RETCODE SCIPhashmapSetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   );

/** checks whether an image to the given origin exists in the hash map */
EXTERN
SCIP_Bool SCIPhashmapExists(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to search for */
   );

/** removes origin->image pair from the hash map, if it exists */
EXTERN
SCIP_RETCODE SCIPhashmapRemove(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to remove from the list */
   );

/** prints statistics about hash map usage */
EXTERN
void SCIPhashmapPrintStatistics(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** indicates whether a hash map has no entries */
EXTERN
SCIP_Bool SCIPhashmapIsEmpty(
   SCIP_HASHMAP*         hashmap             /**< hash map */
);

/** gives the number of entries in a hash map */ 
EXTERN
int SCIPhashmapGetNEntries(
   SCIP_HASHMAP*         hashmap             /**< hash map */
);

/** gives the number of lists (buckets) in a hash map */ 
EXTERN
int SCIPhashmapGetNLists(
   SCIP_HASHMAP*         hashmap             /**< hash map */
);

/** gives a specific list (bucket) in a hash map */
EXTERN
SCIP_HASHMAPLIST* SCIPhashmapGetList(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   int                   listindex           /**< index of hash map list */
);

/** gives the number of entries in a list of a hash map */ 
EXTERN
int SCIPhashmapListGetNEntries(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list, can be NULL */
);

/** retrieves origin of given entry in a hash map */ 
EXTERN
void* SCIPhashmapListGetOrigin(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list */
);

/** retrieves image of given entry in a hash map */ 
EXTERN
void* SCIPhashmapListGetImage(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list */
);

/** retrieves next entry from given entry in a hash map list, or NULL if at end of list. */ 
EXTERN
SCIP_HASHMAPLIST* SCIPhashmapListGetNext(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list */
);

/** removes all entries in a hash map. */ 
EXTERN
SCIP_RETCODE SCIPhashmapRemoveAll(
   SCIP_HASHMAP*         hashmap             /**< hash map */
);

/**@} */

/*
 * Resource Profile
 */

/**@defgroup ResourceProfile Resource Profile
 *
 *@{
 */

/** creates resource profile */
EXTERN
SCIP_RETCODE SCIPprofileCreate(
   SCIP_PROFILE**        profile,            /**< pointer to store the resource profile */
   int                   capacity            /**< resource capacity */
   );

/** frees given resource profile */
EXTERN
void SCIPprofileFree(
   SCIP_PROFILE**        profile             /**< pointer to the resource profile */
   );

/** output of the given resource profile */
EXTERN
void SCIPprofilePrint(
   SCIP_PROFILE*         profile,            /**< resource profile to output */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** returns the capacity of the resource profile */
EXTERN
int SCIPprofileGetCapacity(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the number time points of the resource profile */
EXTERN
int SCIPprofileGetNTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the time points of the resource profile */
EXTERN
int* SCIPprofileGetTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the loads of the resource profile */
EXTERN
int* SCIPprofileGetLoads(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the time point for given position of the resource profile */
EXTERN
int SCIPprofileGetTime(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   pos                 /**< position */
   );

/** returns the loads of the resource profile at the given position */
EXTERN
int SCIPprofileGetLoad(
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   pos                 /**< position */
   );

/** returns if the given time point exists in the resource profile and stores the position of the given time point if it
 *  exists; otherwise the position of the next smaller existing time point is stored
 */
EXTERN
SCIP_Bool SCIPprofileFindLeft(
   SCIP_PROFILE*         profile,            /**< resource profile to search */
   int                   timepoint,          /**< time point to search for */
   int*                  pos                 /**< pointer to store the position */
   );

/** insert a core into resource profile; if the core is non-empty the resource profile will be updated otherwise nothing
 *  happens
 */
EXTERN
SCIP_RETCODE SCIPprofileInsertCore(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height,             /**< height of the core */
   int*                  pos,                /**< pointer to store the first position were it gets infeasible */
   SCIP_Bool*            infeasible          /**< pointer to store if the core does not fit due to capacity */
   );

/** subtracts the height from the resource profile during core time */
EXTERN
SCIP_RETCODE SCIPprofileDeleteCore(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height              /**< height of the core */
   );

/** return the earliest possible starting point within the time interval [lb,ub] for a given core (given by its height
 *  and duration)
 */
EXTERN
int SCIPprofileGetEarliestFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   est,                /**< earliest starting time of the given core */
   int                   lst,                /**< latest starting time of the given core */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the corer cannot be inserted */
   );

/** return the latest possible starting point within the time interval [lb,ub] for a given core (given by its height and
 *  duration)
 */
EXTERN
int SCIPprofileGetLatestFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   lb,                 /**< earliest possible start point */
   int                   ub,                 /**< latest possible start point */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the core cannot be inserted */
   );

/**@} */

/*
 * Directed graph
 */

/**@defgroup DirectedGraph Directed Graph
 *
 *@{
 */

/** creates directed graph structure */
EXTERN
SCIP_RETCODE SCIPdigraphCreate(
   SCIP_DIGRAPH**        digraph,            /**< pointer to store the created directed graph */
   int                   nnodes              /**< number of nodes */
   );

/** copies directed graph structure */
EXTERN
SCIP_RETCODE SCIPdigraphCopy(
   SCIP_DIGRAPH**        targetdigraph,      /**< pointer to store the copied directed graph */
   SCIP_DIGRAPH*         sourcedigraph       /**< source directed graph */
   );

/** sets the sizes of the successor lists for the nodes in a directed graph and allocates memory for the lists */
EXTERN
SCIP_RETCODE SCIPdigraphSetSizes(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int*                  sizes               /**< sizes of the successor lists */
   );

/** frees given directed graph structure */
EXTERN
void SCIPdigraphFree(
   SCIP_DIGRAPH**        digraph             /**< pointer to the directed graph */
   );

/** add (directed) arc and a related data to the directed graph structure
 *
 *  @note if the arc is already contained, it is added a second time
 */
EXTERN
SCIP_RETCODE SCIPdigraphAddArc(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the arc */
   int                   endnode,            /**< start node of the arc */
   void*                 data                /**< data that should be stored for the arc; or NULL */
   );

/** add (directed) arc to the directed graph structure, if it is not contained, yet
 *
 * @note if there already exists an arc from startnode to endnode, the new arc is not added,
 *       even if its data is different
 */
EXTERN
SCIP_RETCODE SCIPdigraphAddArcSafe(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the arc */
   int                   endnode,            /**< start node of the arc */
   void*                 data                /**< data that should be stored for the arc; or NULL */
   );

/** returns the number of nodes of the given digraph */
EXTERN
int SCIPdigraphGetNNodes(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the total number of arcs in the given digraph */
EXTERN
int SCIPdigraphGetNArcs(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the number of successor nodes of the given node */
EXTERN
int SCIPdigraphGetNSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the number of outgoing arcs is returned */
   );

/** returns the array of indices of the successor nodes; this array must not be changed from outside */
EXTERN
int* SCIPdigraphGetSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the array of outgoing arcs is returned */
   );

/** returns the array of datas corresponding to the arcs originating at the given node, or NULL if no data exist; this
 *  array must not be changed from outside
 */
EXTERN
void** SCIPdigraphGetSuccessorsDatas(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the data corresponding to the outgoing arcs is returned */
   );

/** Compute undirected connected components on the given graph.
 *
 *  @note For each arc, its reverse is added, so the graph does not need to be the directed representation of an
 *        undirected graph.
 */
EXTERN
SCIP_RETCODE SCIPdigraphComputeUndirectedComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   minsize,            /**< all components with less nodes are ignored */
   int*                  components,         /**< array with as many slots as there are nodes in the directed graph
                                              *   to store for each node the component to which it belongs
                                              *   (components are numbered 0 to ncomponents - 1); or NULL, if components
                                              *   are accessed one-by-one using SCIPdigraphGetComponent() */
   int*                  ncomponents         /**< pointer to store the number of components; or NULL, if the
                                              *   number of components is accessed by SCIPdigraphGetNComponents() */
   );

/** Performes an (almost) topological sort on the undirected components of the given directed graph. The undirected
 *  components should be computed before using SCIPdigraphComputeUndirectedComponents().
 *
 *  @note In general a topological sort is not unique.  Note, that there might be directed cycles, that are randomly
 *        broken, which is the reason for having only almost topologically sorted arrays.
 */
EXTERN
SCIP_RETCODE SCIPdigraphTopoSortComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the number of previously computed undirected components for the given directed graph */
EXTERN
int SCIPdigraphGetNComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** Returns the previously computed undirected component of the given number for the given directed graph.
 *  If the components were sorted using SCIPdigraphTopoSortComponents(), the component is (almost) topologically sorted.
 */
EXTERN
void SCIPdigraphGetComponent(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   compidx,            /**< number of the component to return */
   int**                 nodes,              /**< pointer to store the nodes in the component; or NULL, if not needed */
   int*                  nnodes              /**< pointer to store the number of nodes in the component;
                                              *   or NULL, if not needed */
   );

/** frees the component information for the given directed graph */
EXTERN
void SCIPdigraphFreeComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** output of the given directed graph via the given message handler */
EXTERN
void SCIPdigraphPrint(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** prints the given directed graph structure in GML format into the given file */
EXTERN
void SCIPdigraphPrintGml(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   FILE*                 file                /**< file to write to */
   );


/** output of the given directed graph via the given message handler */
EXTERN
void SCIPdigraphPrintComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/**@} */

/*
 * Binary search tree
 */

/**@defgroup BinaryTree Binary Tree
 *
 *@{
 */

/** creates a binary tree node with sorting value and user data */
EXTERN
SCIP_RETCODE SCIPbtnodeCreate(
   SCIP_BT*              tree,               /**< binary search tree */
   SCIP_BTNODE**         node,               /**< pointer to store the created search node */
   void*                 dataptr             /**< user node data pointer, or NULL */
   );

/** frees the binary node including the rooted subtree
 *
 *  @note The user pointer (object) is not freed. If needed, it has to be done by the user.
 */
EXTERN
void SCIPbtnodeFree(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node                /**< node to be freed */
   );

/** returns the user data pointer stored in that node */
EXTERN
void* SCIPbtnodeGetData(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns the parent which can be NULL if the given node is the root */
EXTERN
SCIP_BTNODE* SCIPbtnodeGetParent(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns left child which can be NULL if the given node is a leaf */
EXTERN
SCIP_BTNODE* SCIPbtnodeGetLeftchild(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns right child which can be NULL if the given node is a leaf */
EXTERN
SCIP_BTNODE* SCIPbtnodeGetRightchild(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns the sibling of the node or NULL if does not exist */
EXTERN
SCIP_BTNODE* SCIPbtnodeGetSibling(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns whether the node is a root node */
EXTERN
SCIP_Bool SCIPbtnodeIsRoot(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns whether the node is a leaf */
EXTERN
SCIP_Bool SCIPbtnodeIsLeaf(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns TRUE if the given node is left child */
EXTERN
SCIP_Bool SCIPbtnodeIsLeftchild(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns TRUE if the given node is right child */
EXTERN
SCIP_Bool SCIPbtnodeIsRightchild(
   SCIP_BTNODE*          node                /**< node */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPbtnodeGetData(node)               ((node)->dataptr)
#define SCIPbtnodeGetParent(node)             ((node)->parent)
#define SCIPbtnodeGetLeftchild(node)          ((node)->left)
#define SCIPbtnodeGetRightchild(node)         ((node)->right)
#define SCIPbtnodeGetSibling(node)            ((node)->parent == NULL ? NULL : \
                                               (node)->parent->left == (node) ? (node)->parent->right : (node)->parent->left)
#define SCIPbtnodeIsRoot(node)                ((node)->parent == NULL)
#define SCIPbtnodeIsLeaf(node)                ((node)->left == NULL && (node)->right == NULL)
#define SCIPbtnodeIsLeftchild(node)           ((node)->parent == NULL ? FALSE : (node)->parent->left == (node) ? TRUE : FALSE)
#define SCIPbtnodeIsRightchild(node)          ((node)->parent == NULL ? FALSE : (node)->parent->right == (node) ? TRUE : FALSE)

#endif

/** sets the give node data
 *
 *  @note The old user pointer is not freed.
 */
EXTERN
void SCIPbtnodeSetData(
   SCIP_BTNODE*          node,               /**< node */
   void*                 dataptr             /**< node user data pointer */
   );

/** sets parent node
 *
 *  @note The old parent including the rooted subtree is not delete.
 */
EXTERN
void SCIPbtnodeSetParent(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          parent              /**< new parent node, or NULL */
   );

/** sets left child
 *
 *  @note The old left child including the rooted subtree is not delete.
 */
EXTERN
void SCIPbtnodeSetLeftchild(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          left                /**< new left child, or NULL */
   );

/** sets right child
 *
 *  @note The old right child including the rooted subtree is not delete.
 */
EXTERN
void SCIPbtnodeSetRightchild(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          right               /**< new right child, or NULL */
   );

/** creates an binary tree */
EXTERN
SCIP_RETCODE SCIPbtCreate(
   SCIP_BT**             tree,               /**< pointer to store the created binary tree */
   BMS_BLKMEM*           blkmem              /**< block memory used to create nodes */
   );

/** frees binary tree
 *
 *  @note The user pointers (object) of the search nodes are not freed. If needed, it has to be done by the user.
 */
EXTERN
void SCIPbtFree(
   SCIP_BT**             tree                /**< pointer to binary tree */
   );

/** prints the binary tree in GML format into the given file */
EXTERN
void SCIPbtPrintGml(
   SCIP_BT*              tree,               /**< binary tree */
   FILE*                 file                /**< file to write to */
   );

/** returns whether the binary tree is empty (has no nodes) */
EXTERN
SCIP_Bool SCIPbtIsEmpty(
   SCIP_BT *             tree                /**< binary tree */
   );

/** returns the root node of the binary tree or NULL if the binary tree is empty */
EXTERN
SCIP_BTNODE* SCIPbtGetRoot(
   SCIP_BT*              tree                /**< tree to be evaluated */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPbtIsEmpty(tree) (tree->root == NULL)
#define SCIPbtGetRoot(tree) (tree->root)

#endif

/** sets root node
 *
 *  @note The old root including the rooted subtree is not delete.
 */
EXTERN
void SCIPbtSetRoot(
   SCIP_BT*              tree,               /**< tree to be evaluated */
   SCIP_BTNODE*          root                /**< new root, or NULL */
   );

/**@} */

/**@} */

/*
 * Sorting algorithms
 */

/**@defgroup SortingAlgorithms Sorting Algorithms
 *
 * @{
 */

/** default comparer for integers */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPsortCompInt);

/* first all upwards-sorting methods */

/** sort an indexed element set in non-decreasing order, resulting in a permutation index array */
EXTERN
void SCIPsort(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-decreasing order */
EXTERN
void SCIPsortInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-decreasing order */
EXTERN
void SCIPsortPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );


/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Reals in non-decreasing order */
EXTERN
void SCIPsortReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort array of ints in non-decreasing order */
EXTERN
void SCIPsortInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Longints in non-decreasing order */
EXTERN
void SCIPsortLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/* now all downwards-sorting methods */

/** sort an indexed element set in non-increasing order, resulting in a permutation index array */
EXTERN
void SCIPsortDown(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-increasing order */
EXTERN
void SCIPsortDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-increasing order */
EXTERN
void SCIPsortDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Reals in non-increasing order */
EXTERN
void SCIPsortDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );


/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort array of ints in non-increasing order */
EXTERN
void SCIPsortDownInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Longints in non-increasing order */
EXTERN
void SCIPsortDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two three arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/*
 * Sorted vectors
 */

/* upwards insertion */

/** insert a new element into an index array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtr(
   void**                ptrarray,           /**< pointer to the pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Bools/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,          /**< pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   intval,             /**<  additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an arrays of Reals, sorted in non-decreasing order */
EXTERN
void SCIPsortedvecInsertReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted  */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Longint          field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealRealIntInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   void*                 field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of ints in non-decreasing order */
EXTERN
void SCIPsortedvecInsertInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int*                  intarray3,          /**< third int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntIntPtr(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntIntReal(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntPtrReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntIntIntPtr(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int*                  intarray3,          /**< second int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntPtrIntReal(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of Longints, sorted in non-decreasing order */
EXTERN
void SCIPsortedvecInsertLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   SCIP_Bool             field5val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );


/* downwards insertion */

/** insert a new element into an index array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of Reals, sorted in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array  where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Longint          field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   void*                 field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of ints in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int*                  intarray3,          /**< third int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   int*                  intarray3,          /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );


/** insert a new element into four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of Longints, sorted in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increased order */
EXTERN
void SCIPsortedvecInsertDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   SCIP_Bool             field5val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/* upwards position deletion */

/** delete the element at the given position from an index array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Bools/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an arrays of Reals, sorted in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealPtrPtrInt(
   SCIP_Real*            realarray,          /**<  first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array  where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of ints in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int*                  intarray3,          /**< third int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntPtrReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted by in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/* downwards position deletion */

/** delete the element at the given position from an index array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Reals, sorted in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of ints in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int*                  intarray3,          /**< third int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosDownIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three two arrays of Long/pointer, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/* upwards binary search */

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindInd(
   int*                  indarray,           /**< index array to be searched */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindPtr(
   void**                ptrarray,           /**< pointer array to be searched */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be searched */
   SCIP_Real             val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindInt(
   int*                  intarray,           /**< int array to be searched */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be searched */
   SCIP_Longint          val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );


/* downwards binary search */

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindDownInd(
   int*                  indarray,           /**< index array to be searched */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindDownPtr(
   void**                ptrarray,           /**< pointer array to be searched */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be searched */
   SCIP_Real             val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindDownInt(
   int*                  intarray,           /**< int array to be searched */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
EXTERN
SCIP_Bool SCIPsortedvecFindDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be searched */
   SCIP_Longint          val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/**@} */

/**@defgroup MiscellaneousMethods Miscellaneous Methods
 *
 * Below you find a list of miscellaneous methods grouped by different categories
 *@{
 */

/*
 * Numerical methods
 */

/**@defgroup NumericalMethods Numerical Methods
 *
 *@{
 */

/** returns the machine epsilon: the smallest number eps > 0, for which 1.0 + eps > 1.0 */
EXTERN
SCIP_Real SCIPcalcMachineEpsilon(
   void
   );

/** calculates the greatest common divisor of the two given values */
EXTERN
SCIP_Longint SCIPcalcGreComDiv(
   SCIP_Longint          val1,               /**< first value of greatest common devisor calculation */
   SCIP_Longint          val2                /**< second value of greatest common devisor calculation */
   );

/** calculates the smallest common multiple of the two given values */
EXTERN
SCIP_Longint SCIPcalcSmaComMul(
   SCIP_Longint          val1,               /**< first value of smallest common multiple calculation */
   SCIP_Longint          val2                /**< second value of smallest common multiple calculation */
   );

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
EXTERN
SCIP_Bool SCIPrealToRational(
   SCIP_Real             val,                /**< real value r to convert into rational number */
   SCIP_Real             mindelta,           /**< minimal allowed difference r - q of real r and rational q = n/d */
   SCIP_Real             maxdelta,           /**< maximal allowed difference r - q of real r and rational q = n/d */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   );

/** tries to find a value, such that all given values, if scaled with this value become integral */
EXTERN
SCIP_RETCODE SCIPcalcIntegralScalar(
   SCIP_Real*            vals,               /**< values to scale */
   int                   nvals,              /**< number of values to scale */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   );

/** given a (usually very small) interval, tries to find a rational number with simple denominator (i.e. a small
 *  number, probably multiplied with powers of 10) out of this interval; returns TRUE iff a valid rational
 *  number inside the interval was found
 */
EXTERN
SCIP_Bool SCIPfindSimpleRational(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed for resulting rational number */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   );

/** given a (usually very small) interval, selects a value inside this interval; it is tried to select a rational number
 *  with simple denominator (i.e. a small number, probably multiplied with powers of 10);
 *  if no valid rational number inside the interval was found, selects the central value of the interval
 */
EXTERN
SCIP_Real SCIPselectSimpleValue(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom             /**< maximal denominator allowed for resulting rational number */
   );

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
EXTERN
SCIP_Real SCIPrelDiff(
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPrelDiff(val1, val2)         ( ((val1)-(val2))/(MAX3(1.0,REALABS(val1),REALABS(val2))) )

#endif

/**@} */


/*
 * Random Numbers
 */

/**@defgroup RandomNumbers Random Numbers
 *
 *@{
 */

/** returns a random integer between minrandval and maxrandval */
EXTERN
int SCIPgetRandomInt(
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   );

/** returns a random real between minrandval and maxrandval */
EXTERN
SCIP_Real SCIPgetRandomReal(
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   );

/**@} */


/*
 * Additional math functions
 */

/**@defgroup AdditionalMathFunctions Additional math functions
 *
 *@{
 */

/** calculates a binomial coefficient n over m, choose m elements out of n, maximal value will be 33 over 16 (because
 *  the n=33 is the last line in the Pascal's triangle where each entry fits in a 4 byte value), an error occurs due to
 *  big numbers or an negative value m (and m < n) and -1 will be returned
 */
EXTERN
SCIP_Longint SCIPcalcBinomCoef(
   int                   n,                  /**< number of different elements */
   int                   m                   /**< number to choose out of the above */
   );

/**@} */

/*
 * Permutations / Shuffling
 */

/**@defgroup PermutationsShuffling Permutations Shuffling
 *
 *@{
 */

/** swaps two ints */
EXTERN
void SCIPswapInts(
   int*                  value1,             /**< pointer to first integer */
   int*                  value2              /**< pointer ti second integer */
   );

/** swaps the addresses of two pointers */
EXTERN
void SCIPswapPointers(
   void**                pointer1,           /**< first pointer */
   void**                pointer2            /**< second pointer */
   );

/** randomly shuffles parts of an integer array using the Fisher-Yates algorithm */
EXTERN
void SCIPpermuteIntArray(
   int*                  array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end,                /**< last index that should be subject to shuffling (array size for whole
					      *   array)
					      */
   unsigned int*         randseed            /**< seed value for the random generator */
   );

/** randomly shuffles parts of an array using the Fisher-Yates algorithm */
EXTERN
void SCIPpermuteArray(
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end,                /**< last index that should be subject to shuffling (array size for whole
					      *   array)
					      */
   unsigned int*         randseed            /**< pointer to seed value for the random generator */
   );

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 */
EXTERN
SCIP_RETCODE SCIPgetRandomSubset(
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems,          /**< number of elements that should be drawn and stored */
   unsigned int          randseed            /**< seed value for random generator */
   );

/**@} */

/*
 * Strings
 */

/**@defgroup StringMethods String Methods
 *
 *@{
 */

/** copies characters from 'src' to 'dest', copying is stopped when either the 'stop' character is reached or after
 *  'cnt' characters have been copied, whichever comes first.
 *
 *  @note undefined behaviuor on overlapping arrays
 */
EXTERN
int SCIPmemccpy(
   char*                 dest,               /**< destination pointer to copy to */
   const char*           src,                /**< source pointer to copy to */
   char                  stop,               /**< character when found stop copying */
   unsigned int          cnt                 /**< maximal number of characters to copy too */
   );

/** prints an error message containing of the given string followed by a string describing the current system error;
 *  prefers to use the strerror_r method, which is threadsafe; on systems where this method does not exist,
 *  NO_STRERROR_R should be defined (see INSTALL), in this case, srerror is used which is not guaranteed to be
 *  threadsafe (on SUN-systems, it actually is) 
 */
EXTERN
void SCIPprintSysError(
   const char*           message             /**< first part of the error message, e.g. the filename */
   );

/** extracts tokens from strings - wrapper method for strtok_r() */
EXTERN
char* SCIPstrtok(
   char*                 s,                  /**< string to parse */
   const char*           delim,              /**< delimiters for parsing */
   char**                ptrptr              /**< pointer to working char pointer - must stay the same while parsing */
   );

/** translates the given string into a string where symbols ", ', and spaces are escaped with a \ prefix */
EXTERN
void SCIPescapeString(
   char*                 t,                  /**< target buffer to store escaped string */
   int                   bufsize,            /**< size of buffer t */
   const char*           s                   /**< string to transform into escaped string */
   );

/** safe version of snprintf */
EXTERN
int SCIPsnprintf(
   char*                 t,                  /**< target string */
   int                   len,                /**< length of the string to copy */
   const char*           s,                  /**< source string */
   ...                                       /**< further parameters */
   );

/** extract the next token as a integer value if it is one; in case no value is parsed the endptr is set to str */
EXTERN
SCIP_Bool SCIPstrToIntValue(
   const char*           str,                /**< string to search */
   int*                  value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed */
   );

/** extract the next token as a double value if it is one; in case a value is parsed the endptr is set to NULL */
EXTERN
SCIP_Bool SCIPstrToRealValue(
   const char*           str,                /**< string to search */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed
                                              *   otherwise NULL */
   );

/** copies the string between a start and end character */
EXTERN
void SCIPstrCopySection(
   const char*           str,                /**< string to search */
   char                  startchar,          /**< character which defines the beginning */
   char                  endchar,            /**< character which defines the ending */
   char*                 token,              /**< string to store the copy */
   int                   size,               /**< size of the token char array */
   char**                endptr              /**< pointer to store the final string position if successfully parsed
                                              *   otherwise NULL */
   );

/**@} */

/*
 * File methods
 */

/**@defgroup FileMethods File Methods
 *
 *@{
 */

/** returns, whether the given file exists */
EXTERN
SCIP_Bool SCIPfileExists(
   const char*           filename            /**< file name */
   );

/** splits filename into path, name, and extension */
EXTERN
void SCIPsplitFilename(
   char*                 filename,           /**< filename to split; is destroyed (but not freed) during process */
   char**                path,               /**< pointer to store path, or NULL if not needed */
   char**                name,               /**< pointer to store name, or NULL if not needed */
   char**                extension,          /**< pointer to store extension, or NULL if not needed */
   char**                compression         /**< pointer to store compression extension, or NULL if not needed */
   );

/**@} */

/**@} */

#ifdef __cplusplus
}
#endif

#endif
