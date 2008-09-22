/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: probdata_coloring.h,v 1.2 2008/09/22 16:21:32 bzfgamra Exp $"

/**@file   probdata_coloring.h
 * @brief  problem data for coloring algorithm
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_COLORING__
#define __SCIP_PROBDATA_COLORING__

#include "scip/scip.h"
#include "tclique/tclique.h"   /* def. of clique data structures */


/* stable set / variable functions */

extern
/** returns the number of stable sets / variables */
int COLORprobGetNStableSets(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** returns the stable set with the given index */
void COLORprobGetStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   index,              /**< index of the stable set */
   int**                 stableset,          /**< return value: pointer to the stable set */
   int*                  nelements           /**< return value: number of elements in the stable set */
   );

extern
/** returns all stable sets */
void COLORprobGetStableSets(
   SCIP*                 scip,               /**< SCIP data structure */
   int***                stablesets,         /**< return value: pointer to the stable sets */
   int**                 nelements,          /**< return value: number of elements in the stable sets */
   int*                  nstablesets               /**< return value: number of sets */
   );

extern
/** adds a new stable set, the set must be sorted descendingly, 
    attention: you need to check whether it is new before adding it*/
SCIP_RETCODE COLORprobAddNewStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  cliquenodes,        /**< array of nodes in the stable set */
   int                   ncliquenodes,       /**< number of nodes in the stable set */
   int*                  index               /**< return value: index of the stable set, -i-1 if set was not new 
                                                and is already saved as set i */
   );

extern
/** adds a variable that belongs to a given stable set */
void COLORprobAddVarForStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   index,              /**< index of the stable set */
   SCIP_VAR*             var                 /**< pointer to the variable */
   );

extern
/** gets the variable belonging to a given stable set */
SCIP_VAR* COLORprobGetVarForStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   index               /**< index of the stable set */
   );

extern
/** checks whether the given stable set is new
    returns TRUE if the stable is new, 
            FALSE if it is part of or equal to an already existing stable set */
SCIP_Bool COLORprobStableSetIsNew(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  stablesetnodes,     /**< array of nodes in the stable set */
   int                   nstablesetnodes     /**< number of nodes in the stable set */
   );

extern
/** checks whether the first set is equal to the second set, both sets have to be sorted in a decreasing way */
SCIP_Bool COLORprobStableSetsAreEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  set1,               /**< array of nodes in the first set */ 
   int                   nset1nodes,         /**< number of nodes in the first set */
   int*                  set2,               /**< array of nodes in the second set */ 
   int                   nset2nodes          /**< number of nodes in the second set */
   );

extern
/** prints all stable sets to standart output */
void COLORprobPrintStableSets(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** prints the requested stable set to standart output */
void COLORprobPrintStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setnumber           /**< the number of the requested set */
   );

extern
/** checks whether a node is in a given stable set, returns true iff it is */
SCIP_Bool COLORprobIsNodeInStableSet( 
   SCIP*                 scip,               /**< SCIP data structure */
   int                   index,              /**< index of the stable set */
   int                   node                /**< number of the node */
   );

/* constraint functions */

extern
/** creates all node-constraints and saves them in an array */
SCIP_RETCODE COLORprobSetUpArrayOfCons(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** returns all node-constraints */
SCIP_CONS** COLORprobGetConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** returns the node-constraint belonging to a given node */
SCIP_CONS* COLORprobGetConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node                /**< number of the node, for which this constraint assures coloring */
   );



/* graph and preprocessing functions */

extern
/** returns the (preprocessed) graph */
TCLIQUE_GRAPH* COLORprobGetGraph(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** returns the original graph */
TCLIQUE_GRAPH* COLORprobGetOriginalGraph(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** computes the complementary graph for a given graph and stores it in the given pointer */
void COLORprobGetComplementaryGraph(
   SCIP*                 scip,                /**< SCIP data structure */
   TCLIQUE_GRAPH*        graph,               /**< the given graph */
   TCLIQUE_GRAPH*        cgraph               /**< the complementary graph is saved in here */
   );

extern
/** returns the number of nodes in the (preprocessed) graph */
int COLORprobGetNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** returns the number of nodes in the original graph */
int COLORprobGetOriginalNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** returns the array of nodes deleted during preprocessing, length = COLORprobGetOriginalNNodes(), filled with -1 at the end */
int* COLORprobGetDeletedNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** returns the array in which for every node in the preprocessed graph, the related node in the original graph is saved */
int* COLORprobGetOriginalNodesForNewNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** returns the node in the preprocessed graph, that belongs to the given node, returns -1 if node was deleted */
int COLORprobGetNewNodeForOriginalNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node                /**< a node in the original graph */
   );



/* miscellaneous functions */

extern
/** checks whether the given node is in the given array */ 
SCIP_Bool COLORprobIsNodeInArray(
   int                   node,               /**< the number of the node */
   int*                  arrayNodes,         /**< the nodes of the maximum clique */
   int                   narraynodes         /**< number of nodes in the maximum clique */
   );

extern
/** checks whether the given two given arrays are equal */
SCIP_Bool COLORprobEqualSortedArrays(
   int*                  array1nodes,         /**< the nodes of the first set */
   int                   narray1nodes,        /**< number of nodes in the first set */
   int*                  array2nodes,         /**< the nodes of the second set */
   int                   narray2nodes         /**< number of nodes in the second set */
   );


/* create probdate */

extern
/** sets up the problem data */
SCIP_RETCODE SCIPcreateProbColoring(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */           
   int                   nnodes,             /**< number of nodes */
   int                   nedges,             /**< number of edges */
   int**                 edges               /**< array with start- and endpoints of the edges */
   );

#endif
