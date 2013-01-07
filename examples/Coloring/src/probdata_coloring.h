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

/**@file   probdata_coloring.h
 * @brief  problem data for vertex coloring algorithm
 * @author Gerald Gamrath
 *
 * This file implements the problem data for the coloring algorithm.
 *
 * The problem data contains the original graph, preprocessing information, the preprocessed graph,
 * the array with the node-constraints, and an array with all stable sets and corresponding
 * variables.
 *
 * The preprocessing deletes nodes that have a lower degree than the size of a maximum clique.
 * Additionally, it also deletes nodes that have a dominated neighborhood. For further information,
 * look at the documentation for the method preprocessGraph().
 *
 * The deleted nodes and the relation between the nodes of the original graph and the nodes of the
 * preprocessed graph are stored in order to convert a solution of the preprocessed problem to a
 * solution for the original graph and vice versa.
 *
 * Each variable has a pointer of type SCIP_VARDATA* that is used in this case to store an integer
 * representing the number of the stable set. With the aid of this, the corresponding stable set can
 * be found in the array returned by COLORprobGetStableSets().  This array contains all stable sets
 * and is also used to check whether a stable set found by the pricer is really new. This can be
 * done by calling COLORprobStableSetIsNew(). All sets are sorted decreasingly with respect to the
 * indices of the nodes. New candidates should also be sorted that way.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_COLORING__
#define __SCIP_PROBDATA_COLORING__

#include <assert.h>

#include "scip/scip.h"
#include "tclique/tclique.h"   /* def. of clique data structures */
#include "scip/cons_setppc.h"
#include "reader_col.h"

#ifdef __cplusplus
extern "C" {
#endif

/* stable set / variable functions */

/** returns the number of stable sets / variables */
extern
int COLORprobGetNStableSets(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the stable set with the given index */
extern
void COLORprobGetStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,           /**< index of the stable set */
   int**                 stableset,          /**< return value: pointer to the stable set */
   int*                  nelements           /**< return value: number of elements in the stable set */
   );

/** returns all stable sets */
extern
void COLORprobGetStableSets(
   SCIP*                 scip,               /**< SCIP data structure */
   int***                stablesets,         /**< return value: pointer to the stable sets */
   int**                 nelements,          /**< return value: number of elements in the stable sets */
   int*                  nstablesets         /**< return value: number of sets */
   );

/** adds a new stable set, the set must be sorted descendingly, 
    attention: you need to check whether it is new before adding it*/
extern
SCIP_RETCODE COLORprobAddNewStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  cliquenodes,        /**< array of nodes in the stable set */
   int                   ncliquenodes,       /**< number of nodes in the stable set */
   int*                  setindex            /**< return value: index of the stable set, -i-1 if set was not new 
                                              *   and is already saved as set i */
   );

/** adds a variable that belongs to a given stable set */
extern
SCIP_RETCODE COLORprobAddVarForStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,              /**< index of the stable set */
   SCIP_VAR*             var                 /**< pointer to the variable */
   );

/** gets the variable belonging to a given stable set */
extern
SCIP_VAR* COLORprobGetVarForStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex            /**< index of the stable set */
   );

/** checks whether the given stable set is new, returns TRUE if the stable is new and 
 *  FALSE if it is part of or equal to an already existing stable set 
 */
extern
SCIP_Bool COLORprobStableSetIsNew(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  stablesetnodes,     /**< array of nodes in the stable set */
   int                   nstablesetnodes     /**< number of nodes in the stable set */
   );

/** checks whether the first set is equal to the second set, both sets have to be sorted in a decreasing way */
extern
SCIP_Bool COLORprobStableSetsAreEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  set1,               /**< array of nodes in the first set */ 
   int                   nset1nodes,         /**< number of nodes in the first set */
   int*                  set2,               /**< array of nodes in the second set */ 
   int                   nset2nodes          /**< number of nodes in the second set */
   );

/** prints all stable sets to standart output */
extern
void COLORprobPrintStableSets(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** prints the requested stable set to standart output */
extern
void COLORprobPrintStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setnumber           /**< the number of the requested set */
   );

/** checks whether a node is in a given stable set, returns true iff it is */
extern
SCIP_Bool COLORprobIsNodeInStableSet( 
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,           /**< index of the stable set */
   int                   node                /**< number of the node */
   );

/* constraint functions */

/** creates all node-constraints and saves them in an array */
extern
SCIP_RETCODE COLORprobSetUpArrayOfCons(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns all node-constraints */
extern
SCIP_CONS** COLORprobGetConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the node-constraint belonging to a given node */
extern
SCIP_CONS* COLORprobGetConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node                /**< number of the node, for which this constraint assures coloring */
   );



/* graph and preprocessing functions */

/** returns the (preprocessed) graph */
extern
TCLIQUE_GRAPH* COLORprobGetGraph(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the original graph */
extern
TCLIQUE_GRAPH* COLORprobGetOriginalGraph(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** computes the complementary graph for a given graph and stores it in the given pointer */
extern
void COLORprobGetComplementaryGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        graph,              /**< the given graph */
   TCLIQUE_GRAPH*        cgraph              /**< the complementary graph is saved in here */
   );

/** returns the number of nodes in the (preprocessed) graph */
extern
int COLORprobGetNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of nodes in the original graph */
extern
int COLORprobGetOriginalNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array of nodes deleted during preprocessing, length = COLORprobGetOriginalNNodes(), 
 *  filled with -1 at the end 
 */
extern
int* COLORprobGetDeletedNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the array in which for every node in the preprocessed graph, the related node in the original graph is saved */
extern
int* COLORprobGetOriginalNodesForNewNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the node in the preprocessed graph, that belongs to the given node, returns -1 if node was deleted */
extern
int COLORprobGetNewNodeForOriginalNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node                /**< a node in the original graph */
   );



/* miscellaneous functions */

/** checks whether the given node is in the given array */ 
extern
SCIP_Bool COLORprobIsNodeInArray(
   int                   node,               /**< the number of the node */
   int*                  arrayNodes,         /**< the nodes of the maximum clique */
   int                   narraynodes         /**< number of nodes in the maximum clique */
   );

/** checks whether the given two given arrays are equal */
extern
SCIP_Bool COLORprobEqualSortedArrays(
   int*                  array1nodes,         /**< the nodes of the first set */
   int                   narray1nodes,        /**< number of nodes in the first set */
   int*                  array2nodes,         /**< the nodes of the second set */
   int                   narray2nodes         /**< number of nodes in the second set */
   );


/* create probdate */

/** sets up the problem data */
extern
SCIP_RETCODE SCIPcreateProbColoring(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */           
   int                   nnodes,             /**< number of nodes */
   int                   nedges,             /**< number of edges */
   int**                 edges               /**< array with start- and endpoints of the edges */
   );

#ifdef __cplusplus
}
#endif

#endif
