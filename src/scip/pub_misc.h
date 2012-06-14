/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
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
extern
void SCIPgmlWriteNode(
   FILE*                 file,               /**< file to write to */
   unsigned int          id,                 /**< id of the node */
   const char*           label,              /**< label of the node */
   const char*           nodetype,           /**< type of the node, or NULL */
   const char*           fillcolor,          /**< color of the node's interior, or NULL */
   const char*           bordercolor         /**< color of the node's border, or NULL */
   );

/** writes an edge section to the given graph file */
extern
void SCIPgmlWriteEdge(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   );

/** writes an arc section to the given graph file */
extern
void SCIPgmlWriteArc(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   );

/** writes the starting line to a GML graph file, does not open a file */
extern
void SCIPgmlWriteOpening(
   FILE*                 file,               /**< file to write to */
   SCIP_Bool             directed            /**< is the graph directed */
   );

/** writes the ending lines to a GML graph file, does not close a file */
extern
void SCIPgmlWriteCosing(
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
 * Priority Queue
 */

/**@defgroup PriorityQueue Priority Queue
 *
 * @{
 */

/** creates priority queue */
extern
SCIP_RETCODE SCIPpqueueCreate(
   SCIP_PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** frees priority queue, but not the data elements themselves */
extern
void SCIPpqueueFree(
   SCIP_PQUEUE**         pqueue              /**< pointer to a priority queue */
   );

/** clears the priority queue, but doesn't free the data elements themselves */
extern
void SCIPpqueueClear(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** inserts element into priority queue */
extern
SCIP_RETCODE SCIPpqueueInsert(
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   void*                 elem                /**< element to be inserted */
   );

/** removes and returns best element from the priority queue */
extern
void* SCIPpqueueRemove(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the best element of the queue without removing it */
extern
void* SCIPpqueueFirst(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the number of elements in the queue */
extern
int SCIPpqueueNElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
extern
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
extern
int SCIPcalcHashtableSize(
   int                   minsize             /**< minimal size of the hash table */
   );

/** creates a hash table */
extern
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
extern
void SCIPhashtableFree(
   SCIP_HASHTABLE**      hashtable           /**< pointer to the hash table */
   );

/** removes all elements of the hash table
 *
 *  @note From a performance point of view you should not fill and clear a hash table too often since the clearing can
 *        be expensive. Clearing is done by looping over all buckets and removing the hash table lists one-by-one.
 */
extern
void SCIPhashtableClear(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   );

/** inserts element in hash table (multiple inserts of same element possible) */
extern
SCIP_RETCODE SCIPhashtableInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   );

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
extern
SCIP_RETCODE SCIPhashtableSafeInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   );

/** retrieve element with key from hash table, returns NULL if not existing */
extern
void* SCIPhashtableRetrieve(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 key                 /**< key to retrieve */
   );

/** retrieve element with key from hash table, returns NULL if not existing
 * can be used to retrieve all entries with the same key (one-by-one) */
extern
void* SCIPhashtableRetrieveNext(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_HASHTABLELIST**  hashtablelist,      /**< input: entry in hash table list from which to start searching, or NULL; output: entry in hash table list corresponding to element after retrieved one, or NULL */
   void*                 key                 /**< key to retrieve */
   );

/** returns whether the given element exists in the table */
extern
SCIP_Bool SCIPhashtableExists(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to search in the table */
   );

/** removes element from the hash table, if it exists */
extern
SCIP_RETCODE SCIPhashtableRemove(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to remove from the table */
   );

/** prints statistics about hash table usage */
extern
void SCIPhashtablePrintStatistics(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** standard hash key comparator for string keys */
extern
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqString);

/** standard hashing function for string keys */
extern
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValString);

/** gets the element as the key */
extern
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyStandard);

/** returns TRUE iff both keys(pointer) are equal */
extern
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqPtr);

/** returns the hash value of the key */
extern
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
extern
SCIP_RETCODE SCIPhashmapCreate(
   SCIP_HASHMAP**        hashmap,            /**< pointer to store the created hash map */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash map entries */
   int                   mapsize             /**< size of the hash map */
   );

/** frees the hash map */
extern
void SCIPhashmapFree(
   SCIP_HASHMAP**        hashmap             /**< pointer to the hash map */
   );

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
extern
SCIP_RETCODE SCIPhashmapInsert(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   );

/** retrieves image of given origin from the hash map, or NULL if no image exists */
extern
void* SCIPhashmapGetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   );

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
extern
SCIP_RETCODE SCIPhashmapSetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   );

/** checks whether an image to the given origin exists in the hash map */
extern
SCIP_Bool SCIPhashmapExists(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to search for */
   );

/** removes origin->image pair from the hash map, if it exists */
extern
SCIP_RETCODE SCIPhashmapRemove(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to remove from the list */
   );

/** prints statistics about hash map usage */
extern
void SCIPhashmapPrintStatistics(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** indicates whether a hash map has no entries */
extern
SCIP_Bool SCIPhashmapIsEmpty(
   SCIP_HASHMAP*         hashmap             /**< hash map */
);

/** gives the number of entries in a hash map */ 
extern
int SCIPhashmapGetNEntries(
   SCIP_HASHMAP*         hashmap             /**< hash map */
);

/** gives the number of lists (buckets) in a hash map */ 
extern
int SCIPhashmapGetNLists(
   SCIP_HASHMAP*         hashmap             /**< hash map */
);

/** gives a specific list (bucket) in a hash map */
extern
SCIP_HASHMAPLIST* SCIPhashmapGetList(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   int                   listindex           /**< index of hash map list */
);

/** gives the number of entries in a list of a hash map */ 
extern
int SCIPhashmapListGetNEntries(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list, can be NULL */
);

/** retrieves origin of given entry in a hash map */ 
extern
void* SCIPhashmapListGetOrigin(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list */
);

/** retrieves image of given entry in a hash map */ 
extern
void* SCIPhashmapListGetImage(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list */
);

/** retrieves next entry from given entry in a hash map list, or NULL if at end of list. */ 
extern
SCIP_HASHMAPLIST* SCIPhashmapListGetNext(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list */
);

/** removes all entries in a hash map. */ 
extern
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
extern
SCIP_RETCODE SCIPprofileCreate(
   SCIP_PROFILE**        profile,            /**< pointer to store the resource profile */
   int                   capacity            /**< resource capacity */
   );

/** frees given resource profile */
extern
void SCIPprofileFree(
   SCIP_PROFILE**        profile             /**< pointer to the resource profile */
   );

/** output of the given resource profile */
extern
void SCIPprofilePrint(
   SCIP_PROFILE*         profile,            /**< resource profile to output */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** returns the capacity of the resource profile */
extern
int SCIPprofileGetCapacity(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the number time points of the resource profile */
extern
int SCIPprofileGetNTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the time points of the resource profile */
extern
int* SCIPprofileGetTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the loads of the resource profile */
extern
int* SCIPprofileGetLoads(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the time point for given position of the resource profile */
extern
int SCIPprofileGetTime(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   pos                 /**< position */
   );

/** returns the loads of the resource profile at the given position */
extern
int SCIPprofileGetLoad(
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   pos                 /**< position */
   );

/** returns if the given time point exists in the resource profile and stores the position of the given time point if it
 *  exists; otherwise the position of the next smaller existing time point is stored
 */
extern
SCIP_Bool SCIPprofileFindLeft(
   SCIP_PROFILE*         profile,            /**< resource profile to search */
   int                   timepoint,          /**< time point to search for */
   int*                  pos                 /**< pointer to store the position */
   );

/** insert a core into resource profile; if the core is non-empty the resource profile will be updated otherwise nothing
 *  happens
 */
extern
SCIP_RETCODE SCIPprofileInsertCore(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height,             /**< height of the core */
   int*                  pos,                /**< pointer to store the first position were it gets infeasible */
   SCIP_Bool*            infeasible          /**< pointer to store if the core does not fit due to capacity */
   );

/** subtracts the height from the resource profile during core time */
extern
SCIP_RETCODE SCIPprofileDeleteCore(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height              /**< height of the core */
   );

/** return the earliest possible starting point within the time interval [lb,ub] for a given core (given by its height
 *  and duration)
 */
extern
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
extern
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
extern
SCIP_RETCODE SCIPdigraphCreate(
   SCIP_DIGRAPH**        digraph,            /**< pointer to store the created directed graph */
   int                   nnodes              /**< number of nodes */
   );

/** copies directed graph structure */
extern
SCIP_RETCODE SCIPdigraphCopy(
   SCIP_DIGRAPH**        targetdigraph,      /**< pointer to store the copied directed graph */
   SCIP_DIGRAPH*         sourcedigraph       /**< source directed graph */
   );

/** sets the sizes of the successor lists for the nodes in a directed graph and allocates memory for the lists */
extern
SCIP_RETCODE SCIPdigraphSetSizes(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int*                  sizes               /**< sizes of the successor lists */
   );

/** frees given directed graph structure */
extern
void SCIPdigraphFree(
   SCIP_DIGRAPH**        digraph             /**< pointer to the directed graph */
   );

/** add (directed) arc and a related data to the directed graph structure
 *
 *  @note if the arc is already contained, it is added a second time
 */
extern
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
extern
SCIP_RETCODE SCIPdigraphAddArcSafe(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the arc */
   int                   endnode,            /**< start node of the arc */
   void*                 data                /**< data that should be stored for the arc; or NULL */
   );

/** returns the number of nodes of the given digraph */
extern
int SCIPdigraphGetNNodes(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the total number of arcs in the given digraph */
extern
int SCIPdigraphGetNArcs(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the number of successor nodes of the given node */
extern
int SCIPdigraphGetNSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the number of outgoing arcs is returned */
   );

/** returns the array of indices of the successor nodes; this array must not be changed from outside */
extern
int* SCIPdigraphGetSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the array of outgoing arcs is returned */
   );

/** returns the array of datas corresponding to the arcs originating at the given node, or NULL if no data exist; this
 *  array must not be changed from outside
 */
extern
void** SCIPdigraphGetSuccessorsDatas(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the data corresponding to the outgoing arcs is returned */
   );

/** Compute undirected connected components on the given graph.
 *
 *  @note For each arc, its reverse is added, so the graph does not need to be the directed representation of an
 *        undirected graph.
 */
extern
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
extern
SCIP_RETCODE SCIPdigraphTopoSortComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the number of previously computed undirected components for the given directed graph */
extern
int SCIPdigraphGetNComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** Returns the previously computed undirected component of the given number for the given directed graph.
 *  If the components were sorted using SCIPdigraphTopoSortComponents(), the component is (almost) topologically sorted.
 */
extern
void SCIPdigraphGetComponent(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   compidx,            /**< number of the component to return */
   int**                 nodes,              /**< pointer to store the nodes in the component; or NULL, if not needed */
   int*                  nnodes              /**< pointer to store the number of nodes in the component;
                                              *   or NULL, if not needed */
   );

/** frees the component information for the given directed graph */
extern
void SCIPdigraphFreeComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** output of the given directed graph via the given message handler */
extern
void SCIPdigraphPrint(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** prints the given directed graph structure in GML format into the given file */
extern
void SCIPdigraphPrintGml(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   FILE*                 file                /**< file to write to */
   );


/** output of the given directed graph via the given message handler */
extern
void SCIPdigraphPrintComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/**@} */

/*
 * Binary search tree
 */

/**@defgroup BinarySearchTree Binary Search Tree
 *
 *@{
 */

/** creates a search tree node with sorting value and user data */
extern
SCIP_RETCODE SCIPbstnodeCreate(
   SCIP_BST*             tree,               /**< binary search tree */
   SCIP_BSTNODE**        node,               /**< pointer to store the created search node */
   void*                 key,                /**< sorting key, or NULL */
   void*                 dataptr             /**< user node data pointer, or NULL */
   );

/** frees the search node including the rooted subtree
 *
 *  @note The user pointer (object) is not freed. If needed, it has to be done by the user.
 */
extern
void SCIPbstnodeFree(
   SCIP_BST*             tree,               /**< binary search tree */
   SCIP_BSTNODE**        node                /**< search node to be freed */
   );

#ifndef NDEBUG

/** returns whether the search node is a leaf */
extern
SCIP_Bool SCIPbstnodeIsLeaf(
   SCIP_BSTNODE*         node                /**< search node */
   );

/** returns the user data pointer stored in that search node */
extern
void* SCIPbstnodeGetData(
   SCIP_BSTNODE*         node                /**< search node */
   );

/** returns the key of the search node */
extern
void* SCIPbstnodeGetKey(
   SCIP_BSTNODE*         node                /**< search node */
   );

/** returns the parent which can be NULL if the given node is the root */
extern
SCIP_BSTNODE* SCIPbstnodeGetParent(
   SCIP_BSTNODE*         node                /**< search node */
   );

/** returns left child which can be NULL if the given node is a leaf */
extern
SCIP_BSTNODE* SCIPbstnodeGetLeftchild(
   SCIP_BSTNODE*         node                /**< search node */
   );

/** returns right child which can be NULL if the given node is a leaf */
extern
SCIP_BSTNODE* SCIPbstnodeGetRightchild(
   SCIP_BSTNODE*         node                /**< search node */
   );

#else

#define SCIPbstnodeIsLeaf(node)               (node->left == NULL && node->right == NULL)
#define SCIPbstnodeGetData(node)              (node->dataptr)
#define SCIPbstnodeGetKey(node)               (node->key)
#define SCIPbstnodeGetParent(node)            (node->parent)
#define SCIPbstnodeGetLeftchild(node)         (node->left)
#define SCIPbstnodeGetRightchild(node)        (node->right)

#endif

/** sets the give node data
 *
 *  @note The old user pointer is not freed.
 */
extern
void SCIPbstnodeSetData(
   SCIP_BSTNODE*         node,               /**< search node */
   void*                 dataptr             /**< node user data pointer */
   );

/** sets the key to the search node
 *
 *  @note The old key pointer is not freed.
 */
extern
void SCIPbstnodeSetKey(
   SCIP_BSTNODE*         node,               /**< search node */
   void*                 key                 /**< key value */
   );

/** sets parent node
 *
 *  @note The old parent including the rooted subtree is not delete.
 */
extern
void SCIPbstnodeSetParent(
   SCIP_BSTNODE*         node,               /**< search node */
   SCIP_BSTNODE*         parent              /**< new parent node, or NULL */
   );

/** sets left child
 *
 *  @note The old left child including the rooted subtree is not delete.
 */
extern
void SCIPbstnodeSetLeftchild(
   SCIP_BSTNODE*         node,               /**< search node */
   SCIP_BSTNODE*         left                /**< new left child, or NULL */
   );

/** sets right child
 *
 *  @note The old right child including the rooted subtree is not delete.
 */
extern
void SCIPbstnodeSetRightchild(
   SCIP_BSTNODE*         node,               /**< search node */
   SCIP_BSTNODE*         right               /**< new right child, or NULL */
   );

/** creates an binary search tree */
extern
SCIP_RETCODE SCIPbstCreate(
   SCIP_BST**            tree,               /**< pointer to store the created binary search tree */
   BMS_BLKMEM*           blkmem,             /**< block memory used to create search node */
   SCIP_DECL_BSTINSERT   ((*inserter)),      /**< inserter used to insert a new search node */
   SCIP_DECL_BSTDELETE   ((*deleter)),       /**< deleter used to delete new search node */
   SCIP_DECL_SORTPTRCOMP ((*comparer))       /**< comparer used to compares two search keys */
   );

/** frees binary search tree
 *
 *  @note The user pointers (object) of the search nodes are not freed. If needed, it has to be done by the user.
 */
extern
void SCIPbstFree(
   SCIP_BST**            tree                /**< pointer to binary search tree */
   );

/** prints the binary search tree in GML format into the given file */
extern
void SCIPbstPrintGml(
   SCIP_BST*             tree,               /**< binary search tree */
   FILE*                 file                /**< file to write to */
   );

#ifndef NDEBUG

/** returns whether the binary search tree is empty (has no nodes) */
extern
SCIP_Bool SCIPbstIsEmpty(
   SCIP_BST*             tree                /**< binary search tree */
   );

/** returns the the root node of the binary search or NULL if the binary search tree is empty */
extern
SCIP_BSTNODE* SCIPbstGetRoot(
   SCIP_BST*             tree                /**< tree to be evaluated */
   );

#else


#define SCIPbstIsEmpty(tree) (tree->root == NULL)
#define SCIPbstGetRoot(tree) (tree->root)

#endif

/** sets root node
 *
 *  @note The old root including the rooted subtree is not delete.
 */
extern
void SCIPbstSetRoot(
   SCIP_BST*             tree,               /**< tree to be evaluated */
   SCIP_BSTNODE*         root                /**< new root, or NULL */
   );

/** inserts the given node into the binary search tree; uses the given SCIP_DECL_BSTINSERT() method */
extern
SCIP_RETCODE SCIPbstInsert(
   SCIP_BST*             tree,               /**< binary search tree */
   SCIP_BSTNODE*         node,               /**< search node */
   SCIP_Bool*            inserted            /**< pointer to store whether the node was inserted */
   );

/** deletes the given node from the binary search tree; uses the given SCIP_DECL_BSTDELETE() method */
extern
SCIP_RETCODE SCIPbstDelete(
   SCIP_BST*             tree,               /**< binary search tree */
   SCIP_BSTNODE*         node,               /**< search node */
   SCIP_Bool*            deleted             /**< pointer to store whether the node was deleted */
   );

/** compares to search nodes using the search tree comparer */
extern
int SCIPbstComp(
   SCIP_BST*             tree,               /**< binary search tree */
   SCIP_BSTNODE*         node1,              /**< search node 1 */
   SCIP_BSTNODE*         node2               /**< search node 2 */
   );

/** Finds the position at which the given node is located in the search tree or has to be inserted. If the search tree
 *  is empty NULL is return. If a search node with the same node key exists, the method returns the last search node and
 *  sets the found pointer to TRUE. If the element does not exist, the method returns the search node with the last
 *  highest node key value which is smaller than the given one and sets the found pointer to FALSE.
 */
extern
SCIP_BSTNODE* SCIPbstFindInsertNode(
   SCIP_BST*             tree,               /**< binary search tree */
   SCIP_BSTNODE*         node,               /**< search node to find */
   SCIP_Bool*            found               /**< pointer to store if a search node with the given key was found */
   );

/** Finds the position at which the given key is located in the search tree. If the search tree is empty NULL is
 *  return. If a search node with the same key exists, the method returns the last search node and sets the found
 *  pointer to TRUE. If the element does not exist, the method returns the search node with the last highest key value
 *  which is smaller than the given one and sets the found pointer to FALSE.
 */
extern
SCIP_BSTNODE* SCIPbstFindKey(
   SCIP_BST*             tree,               /**< binary search tree */
   void*                 key,                /**< key value */
   SCIP_Bool*            found               /**< pointer to store if a search node with the given key was found */
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
extern
SCIP_DECL_SORTPTRCOMP(SCIPsortCompInt);

/* first all upwards-sorting methods */

/** sort an indexed element set in non-decreasing order, resulting in a permutation index array */
extern
void SCIPsort(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-decreasing order */
extern
void SCIPsortInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-decreasing order */
extern
void SCIPsortPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );


/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
extern
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
extern
void SCIPsortReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
extern
void SCIPsortRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort array of ints in non-decreasing order */
extern
void SCIPsortInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order */
extern
void SCIPsortIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Longints in non-decreasing order */
extern
void SCIPsortLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
extern
void SCIPsortLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
extern
void SCIPsortLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
extern
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
extern
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
extern
void SCIPsortDown(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-increasing order */
extern
void SCIPsortDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-increasing order */
extern
void SCIPsortDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
extern
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
extern
void SCIPsortDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );


/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort array of ints in non-increasing order */
extern
void SCIPsortDownInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Longints in non-increasing order */
extern
void SCIPsortDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
extern
void SCIPsortDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two three arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
extern
void SCIPsortDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
extern
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
extern
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
extern
void SCIPsortedvecInsertInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-decreasing order */
extern
void SCIPsortedvecInsertPtr(
   void**                ptrarray,           /**< pointer to the pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers sorted by first array in non-decreasing order */
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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

/** insert a new element into four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
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
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecInsertRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
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
extern
void SCIPsortedvecInsertReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecInsertInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecInsertLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
extern
void SCIPsortedvecInsertLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecInsertDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-increasing order */
extern
void SCIPsortedvecInsertDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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

/** insert a new element into four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecInsertDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
extern
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
extern
void SCIPsortedvecInsertDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/pointers, sorted by first array in non-increasing order */
extern
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
extern
void SCIPsortedvecInsertDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecInsertDownInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecInsertDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
extern
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
extern
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
extern
void SCIPsortedvecInsertDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
extern
void SCIPsortedvecInsertDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecDelPosInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-decreasing order */
extern
void SCIPsortedvecDelPosPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecDelPosRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an arrays of Reals, sorted in non-decreasing order */
extern
void SCIPsortedvecDelPosReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealPtrPtrInt(
   SCIP_Real*            realarray,          /**<  first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array  where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
extern
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
extern
void SCIPsortedvecDelPosInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int*                  intarray3,          /**< third int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntPtrReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted by in non-decreasing order */
extern
void SCIPsortedvecDelPosLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
extern
void SCIPsortedvecDelPosLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-decreasing order */
extern
void SCIPsortedvecDelPosLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecDelPosDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
extern
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
extern
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
extern
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
extern
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
extern
void SCIPsortedvecDelPosDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
extern
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
extern
void SCIPsortedvecDelPosDownInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int*                  intarray3,          /**< third int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosDownIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosDownIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted in non-increasing order */
extern
void SCIPsortedvecDelPosDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three two arrays of Long/pointer, sorted by the first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
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
extern
SCIP_Real SCIPcalcMachineEpsilon(
   void
   );

/** calculates the greatest common divisor of the two given values */
extern
SCIP_Longint SCIPcalcGreComDiv(
   SCIP_Longint          val1,               /**< first value of greatest common devisor calculation */
   SCIP_Longint          val2                /**< second value of greatest common devisor calculation */
   );

/** calculates the smallest common multiple of the two given values */
extern
SCIP_Longint SCIPcalcSmaComMul(
   SCIP_Longint          val1,               /**< first value of smallest common multiple calculation */
   SCIP_Longint          val2                /**< second value of smallest common multiple calculation */
   );

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
extern
SCIP_Bool SCIPrealToRational(
   SCIP_Real             val,                /**< real value r to convert into rational number */
   SCIP_Real             mindelta,           /**< minimal allowed difference r - q of real r and rational q = n/d */
   SCIP_Real             maxdelta,           /**< maximal allowed difference r - q of real r and rational q = n/d */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   );

/** tries to find a value, such that all given values, if scaled with this value become integral */
extern
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
extern
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
extern
SCIP_Real SCIPselectSimpleValue(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom             /**< maximal denominator allowed for resulting rational number */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
extern
SCIP_Real SCIPrelDiff(
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
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
extern
int SCIPgetRandomInt(
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   );

/** returns a random real between minrandval and maxrandval */
extern
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
extern
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

/** swaps the addresses of two pointers */
extern
void SCIPswapPointers(
   void**                pointer1,           /**< first pointer */
   void**                pointer2            /**< second pointer */
   );

/** randomly shuffles parts of an array using the Fisher-Yates algorithm */
extern
void SCIPpermuteArray(
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end,                /**< last index that should be subject to shuffling (array size for whole array) */
   unsigned int*         randseed            /**< pointer to seed value for the random generator */
   );

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 */
extern
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
extern
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
extern
void SCIPprintSysError(
   const char*           message             /**< first part of the error message, e.g. the filename */
   );

/** extracts tokens from strings - wrapper method for strtok_r() */
extern
char* SCIPstrtok(
   char*                 s,                  /**< string to parse */
   const char*           delim,              /**< delimiters for parsing */
   char**                ptrptr              /**< pointer to working char pointer - must stay the same while parsing */
   );

/** translates the given string into a string where symbols ", ', and spaces are escaped with a \ prefix */
extern
void SCIPescapeString(
   char*                 t,                  /**< target buffer to store escaped string */
   int                   bufsize,            /**< size of buffer t */
   const char*           s                   /**< string to transform into escaped string */
   );

/** safe version of snprintf */
extern
int SCIPsnprintf(
   char*                 t,                  /**< target string */
   int                   len,                /**< length of the string to copy */
   const char*           s,                  /**< source string */
   ...                                       /**< further parameters */
   );

/** extract the next token as a double value if it is one; in case a value is parsed the endptr is set to NULL */
extern 
SCIP_Bool SCIPstrToRealValue(
   const char*           str,                /**< string to search */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed
                                              *   otherwise NULL */
   );

/** copies the string between a start and end character */
extern
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
extern
SCIP_Bool SCIPfileExists(
   const char*           filename            /**< file name */
   );

/** splits filename into path, name, and extension */
extern
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
