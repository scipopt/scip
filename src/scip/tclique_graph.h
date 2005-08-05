/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                        This file is part of the program                   */
/*                    TCLIQUE --- Algorithm for Maximum Cliques              */
/*                                                                           */
/*    Copyright (C) 1996-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: tclique_graph.h,v 1.8 2005/08/05 16:04:32 bzfpfend Exp $"

/**@file   tclique_graph.h
 * @brief  graph interface used in the tclique methods
 * @author Tobias Achterberg
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TCLIQUE_GRAPH_H__
#define __SCIP_TCLIQUE_GRAPH_H__

#include "scip/def.h"



/*
 * Data Types and Structures
 */

typedef int  WEIGHT;

typedef struct _TCLIQUEDATA TCLIQUEDATA;




/*
 * Interface Methods used by the TClique algorithm
 */

/** gets number of nodes in the graph */
extern
int tcliqueGetNNodes(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets weight of nodes in the graph */
extern
WEIGHT* tcliqueGetWeights(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** returns, whether the edge (node1, node2) is in the graph */
extern
int tcliqueIsEdge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node1,              /**< start node of edge */
   int              node2               /**< end node of edge */
   );

/* selects all nodes from a given set of nodes which are adjacent to a given node
 * and returns the number of selected nodes
 */
extern
int tcliqueSelectAdjnodes( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node,               /**< given node */
   int*             nodes,              /**< given set of nodes (ordered by node index) */ 
   int              nnodes,             /**< number of nodes in given set nodes */ 
   int*             adjnodes            /**< pointer to store adjacent nodes of node, contained in nodes */ 
   );




/*
 * External Interface Methods to access the graph (this can be changed without affecting the TClique algorithm)
 */

/** creates tclique data structure */
extern
Bool tcliqueCreate(
   TCLIQUEDATA**    tcliquedata         /**< pointer to store tclique data structure */
   );

/** frees tclique data structure */
extern
void tcliqueFree(
   TCLIQUEDATA**    tcliquedata         /**< pointer to tclique data structure */
   );

/** adds nodes up to the given node number to tclique data structure (intermediate nodes have weight 0) */
extern
Bool tcliqueAddNode(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              node,               /**< node number to add */
   WEIGHT           weight              /**< weight of node to add */
   );

/** changes weight of node in tclique data structure */
extern
void tcliqueChangeWeight(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              node,               /**< node to set new weight */
   WEIGHT           weight              /**< new weight of node (allready scaled) */
   );

/** adds edge (node1, node2) to tclique data structure (node1 and node2 have to be contained in 
 *  tclique data structure);
 *  new edges are cached, s.t. the graph data structures are not correct until a call to tcliqueFlush();
 *  you have to make sure, that no double edges are inserted
 */
extern
Bool tcliqueAddEdge(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              node1,              /**< start node of edge to add */
   int              node2               /**< end node of edge to add */
   );

/** inserts all cached edges into the data structures */
extern
Bool tcliqueFlush(
   TCLIQUEDATA*     tcliquedata         /**< tclique data structure */
   );

/** loads tclique data structure from file */
extern
Bool tcliqueLoadFile(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to store tclique data structure */
   const char*      filename,           /**< name of file with graph data */
   double           scaleval,           /**< value to scale weights (only integral part of scaled weights is considered) */
   char*            probname            /**< buffer to store the name of the problem */
   );

/** saves tclique data structure to file */
extern
Bool tcliqueSaveFile(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   const char*      filename,           /**< name of file to create */
   double           scaleval,           /**< value to unscale weights with */
   const char*      probname            /**< name of the problem */
   );

/** gets number of edges in the graph */
extern
int tcliqueGetNEdges(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets degree of nodes in graph */
extern
int* tcliqueGetDegrees(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets adjacent nodes of edges in graph */
extern
int* tcliqueGetAdjnodes(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets pointer to first adjacent edge of given node in graph */
extern
int* tcliqueGetFirstAdjedge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node                /**< given node */
   );

/** gets pointer to last adjacent edge of given node in graph */
extern
int* tcliqueGetLastAdjedge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node                /**< given node */
   );

/* prints tclique data structure */
extern
void tcliquePrintData(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );


#endif
