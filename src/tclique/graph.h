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
#pragma ident "@(#) $Id: graph.h,v 1.1 2005/03/10 17:11:18 bzfpfend Exp $"

/**@file   graph.h
 * @brief  tclique data part of algorithm for maximum cliques
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GRAPH_H__
#define __GRAPH_H__



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
} TCLIQUEDATA; 




/** creates tclique data structure */
void tcliqueCreate(
   TCLIQUEDATA**    tcliquedata         /**< pointer to store tclique data structure */
   );

/** frees tclique data structure */
void tcliqueFree(
   TCLIQUEDATA**    tcliquedata         /**< pointer to tclique data structure */
   );

/** adds node to tclique data structure */
void tcliqueAddNode(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to tclique data structure */
   int              node,               /**< node to add */
   WEIGHT           weight              /**< weight of node to add (allready scaled) */
   );

/** adds edge (node1, node2) to tclique data structure (node1 and node2 have to be contained in 
 *  tclique data structure) */
void tcliqueAddEdge(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to tclique data structure */
   int              node1,              /**< start node of edge to add */
   int              node2               /**< end node of edge to add */
   );

/** changes weight of node in tclique data structure */
void tcliqueNodeChangeWeight(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to tclique data structure */
   int              node,               /**< node to set new weight */
   WEIGHT           weight              /**< new weight of node (allready scaled) */
   );

/** loads tclique data structure */
BOOL tcliqueLoadComplete(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to store tclique data structure */
   char*            filename,           /**< name of file with graph data */
   char*            probname            /**< name of problem */
   );

/** loads tclique data structure */
BOOL tcliqueLoadStepwise(
   TCLIQUEDATA**    tcliquedata,        /**< pointer to store tclique data structure */
   char*            filename,           /**< name of file with graph data */
   char*            probname            /**< name of problem */
   );

/** gets number of nodes in the graph */
int getNnodes(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets number of edges in the graph */
int getNedges(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets weight of nodes in the graph */
WEIGHT* getWeights(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets degree of nodes in graph */
int* getDegrees(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets adjacent nodes of edges in graph */
int* getAdjnodes(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets pointer to first and one after last adjacent edge of nodes in graph */
HEAD_ADJ* getAdjedges(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );

/** gets pointer to first adjacent edge of given node in graph */
int* getFirstAdjedge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node                /**< given node */
   );

/** gets pointer to last adjacent edge of given node in graph */
int* getLastAdjedge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node                /**< given node */
   );

/** returns, whether the edge (node1, node2) is in the graph */
int isEdge(
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node1,              /**< start node of edge */
   int              node2               /**< end node of edge */
   );

/* selects all nodes from a given set of nodes which are adjacent to a given node
 * and returns the number of selected nodes */
int selectAdjnodes( 
   TCLIQUEDATA*     tcliquedata,        /**< pointer to tclique data structure */
   int              node,               /**< given node */
   int*             nodes,              /**< given set of nodes (ordered by node index) */ 
   int              nnodes,             /**< number of nodes in given set nodes */ 
   int*             adjnodes            /**< pointer to store adjacent nodes of node, contained in nodes */ 
   );

/* prints tclique data structure */
void printTcliquedata(
   TCLIQUEDATA*     tcliquedata         /**< pointer to tclique data structure */
   );


#endif
