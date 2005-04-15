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
#pragma ident "@(#) $Id: tclique_graph.h,v 1.1 2005/04/15 11:46:55 bzfpfend Exp $"

/**@file   tclique_graph.h
 * @brief  tclique data part of algorithm for maximum cliques
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TCLIQUE_GRAPH_H__
#define __TCLIQUE_GRAPH_H__



typedef enum {
   NO, 
   YES
} BOOL;

typedef int  WEIGHT;

typedef struct _HEAD_ADJ{
	int first;
	int last;
} HEAD_ADJ;

typedef struct _TCLIQUEDATA{
   int              nnodes;		/**< number of nodes in graph */
   int              nedges;		/**< number of edges in graph */
   WEIGHT*          weights;	        /**< weight of nodes */
   int*             degrees;	        /**< degree of nodes */
   int*             adjnodes;	        /**< adjacent nodes of edges */
   HEAD_ADJ*        adjedges;           /**< pointer to first and one after last adjacent edge of nodes */
   int              sizenodes;		/**< size of arrays concerning nodes (weights, degrees and adjedges) */
   int              sizeedges;		/**< size of arrays concerning edges (adjnodes) */
   int*             cacheddegrees;      /**< number of adjacent cached edges for each node */
   int*             cachedorigs;        /**< origin nodes of cached edges */
   int*             cacheddests;        /**< destination nodes of cached edges */
   int              ncachededges;       /**< number of cached edges (not yet inserted in all data structures) */
   int              sizecachededges;    /**< size of arrays concerning cached edges */
} TCLIQUEDATA; 




/** creates tclique data structure */
extern
BOOL tcliqueCreate(
   TCLIQUEDATA**    tcliquedata         /**< pointer to store tclique data structure */
   );

/** frees tclique data structure */
extern
void tcliqueFree(
   TCLIQUEDATA**    tcliquedata         /**< pointer to tclique data structure */
   );

/** adds nodes up to the given node number to tclique data structure (intermediate nodes have weight 0) */
extern
BOOL tcliqueAddNode(
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
BOOL tcliqueAddEdge(
   TCLIQUEDATA*     tcliquedata,        /**< tclique data structure */
   int              node1,              /**< start node of edge to add */
   int              node2               /**< end node of edge to add */
   );

/** inserts all cached edges into the data structures */
extern
BOOL tcliqueFlush(
   TCLIQUEDATA*     tcliquedata         /**< tclique data structure */
   );

/** loads tclique data structure from file */
extern
BOOL tcliqueLoadFile(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to store tclique data structure */
   const char*      filename,           /**< name of file with graph data */
   char*            probname            /**< buffer to store the name of the problem */
   );

/** gets number of nodes in the graph */
extern
int tcliqueGetNNodes(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets number of edges in the graph */
extern
int tcliqueGetNEdges(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets weight of nodes in the graph */
extern
WEIGHT* tcliqueGetWeights(
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

/** gets pointer to first and one after last adjacent edge of nodes in graph */
extern
HEAD_ADJ* tcliqueGetAdjedges(
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

/* prints tclique data structure */
extern
void tcliquePrintData(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );


#endif
