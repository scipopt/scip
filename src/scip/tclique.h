/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                        This file is part of the program                   */
/*                    TCLIQUE --- Algorithm for Maximum Cliques              */
/*                                                                           */
/*    Copyright (C) 1996-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  TCLIQUE is distributed under the terms of the ZIB Academic License.      */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with TCLIQUE; see the file COPYING.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: tclique.h,v 1.1 2005/08/08 13:20:36 bzfpfend Exp $"

/**@file   tclique.h
 * @brief  tclique user interface
 * @author Tobias Achterberg
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TCLIQUE_H__
#define __TCLIQUE_H__


/*
 * Data Types and Structures
 */

typedef int  WEIGHT;                    /**< type used for node weights in the graph */
typedef struct TcliqueGraph TCLIQUEGRAPH; /**< user defined structure for storing the graph, passed to graph callbacks */
typedef struct TcliqueData TCLIQUEDATA; /**< user defined data to pass to new solution callback method */

#ifndef Bool
#define Bool unsigned int               /**< type used for boolean values */
#endif
#ifndef TRUE
#define TRUE  1                         /**< boolean value TRUE */
#define FALSE 0                         /**< boolean value FALSE */
#endif



/*
 * User Callback Methods
 */

/** user callback method which is called whenever a feasible clique was found
 *  input:
 *   - tcliquedata  : user data given to tcliqueMaxClique()
 *   - cliquenodes  : array with nodes of the clique
 *   - ncliquenodes : number of nodes in the clique
 *   - cliqueweight : weight of the clique
 *  output:
 *   - minweight    : new minimal weight for feasible cliques
 *   - acceptsol    : setting TRUE makes clique the new best clique, and updates minweight
 *   - stopsolving  : setting TRUE aborts the search for cliques
 */
#define TCLIQUE_NEWSOL(x) void x (TCLIQUEDATA* tcliquedata, int* cliquenodes, int ncliquenodes, \
      WEIGHT cliqueweight, WEIGHT* minweight, Bool* acceptsol, Bool* stopsolving)

/** user callback method to get number of nodes in the graph
 *  input:
 *   - tcliquegraph : user defined graph data structure given to tcliqueMaxClique()
 *  returns:
 *   number of nodes in the graph
 */
#define TCLIQUE_GETNNODES(x) int x (TCLIQUEGRAPH* tcliquegraph)

/** user callback method to get weights of nodes in the graph
 *  input:
 *   - tcliquegraph : user defined graph data structure given to tcliqueMaxClique()
 *  returns:
 *   array of node weights (of length at least equal to the number of nodes in the graph)
 */
#define TCLIQUE_GETWEIGHTS(x) const WEIGHT* x (TCLIQUEGRAPH* tcliquegraph)

/** user callback method to return whether the edge (node1, node2) is in the graph
 *  input:
 *   - tcliquegraph : user defined graph data structure given to tcliqueMaxClique()
 *   - node1        : start node of edge (tail)
 *   - node2        : end node of edge (head)
 *  returns:
 *   TRUE if edge is in the graph, FALSE otherwise
 */
#define TCLIQUE_ISEDGE(x) Bool x (TCLIQUEGRAPH* tcliquegraph, int node1, int node2)

/** user callback method to select all nodes from a given set of nodes which are adjacent to a given node
 *  input:
 *   - tcliquegraph : user defined graph data structure given to tcliqueMaxClique()
 *   - node         : node to select adjacent nodes from
 *   - nodes        : array of nodes to select nodes from
 *   - nnodes       : number of nodes in 'nodes' array
 *   - adjnodes     : pointer to memory to store the resulting nodes
 *                    'adjnodes' and 'nodes' may point to the same memory location
 *  output:
 *   - adjnodes     : array of nodes that are contained in 'nodes' and that are adjacent to 'node'
 *  returns:
 *   number of nodes in 'adjnodes'
 */
#define TCLIQUE_SELECTADJNODES(x) int x (TCLIQUEGRAPH* tcliquegraph, int node, int* nodes, int nnodes, int* adjnodes)




/*
 * Default Graph Implementation: Interface Methods used by the TClique algorithm
 */

/** gets number of nodes in the graph */
extern
TCLIQUE_GETNNODES(tcliqueGetNNodes);

/** gets weight of nodes in the graph */
extern
TCLIQUE_GETWEIGHTS(tcliqueGetWeights);

/** returns, whether the edge (node1, node2) is in the graph */
extern
TCLIQUE_ISEDGE(tcliqueIsEdge);

/* selects all nodes from a given set of nodes which are adjacent to a given node
 * and returns the number of selected nodes
 */
extern
TCLIQUE_SELECTADJNODES(tcliqueSelectAdjnodes);




/*
 * Default Graph Implementation: External Interface Methods to access the graph
 */

/** creates graph data structure */
extern
Bool tcliqueCreate(
   TCLIQUEGRAPH**   tcliquegraph        /**< pointer to store graph data structure */
   );

/** frees graph data structure */
extern
void tcliqueFree(
   TCLIQUEGRAPH**   tcliquegraph        /**< pointer to graph data structure */
   );

/** adds nodes up to the given node number to graph data structure (intermediate nodes have weight 0) */
extern
Bool tcliqueAddNode(
   TCLIQUEGRAPH*    tcliquegraph,       /**< graph data structure */
   int              node,               /**< node number to add */
   WEIGHT           weight              /**< weight of node to add */
   );

/** changes weight of node in graph data structure */
extern
void tcliqueChangeWeight(
   TCLIQUEGRAPH*    tcliquegraph,       /**< graph data structure */
   int              node,               /**< node to set new weight */
   WEIGHT           weight              /**< new weight of node (allready scaled) */
   );

/** adds edge (node1, node2) to graph data structure (node1 and node2 have to be contained in 
 *  graph data structure);
 *  new edges are cached, s.t. the graph data structures are not correct until a call to tcliqueFlush();
 *  you have to make sure, that no double edges are inserted
 */
extern
Bool tcliqueAddEdge(
   TCLIQUEGRAPH*    tcliquegraph,       /**< graph data structure */
   int              node1,              /**< start node of edge to add */
   int              node2               /**< end node of edge to add */
   );

/** inserts all cached edges into the data structures */
extern
Bool tcliqueFlush(
   TCLIQUEGRAPH*    tcliquegraph        /**< graph data structure */
   );

/** loads graph data structure from file */
extern
Bool tcliqueLoadFile(
   TCLIQUEGRAPH**   tcliquegraph,       /**< pointer to store graph data structure */
   const char*      filename,           /**< name of file with graph data */
   double           scaleval,           /**< value to scale weights (only integral part of scaled weights is considered) */
   char*            probname            /**< buffer to store the name of the problem */
   );

/** saves graph data structure to file */
extern
Bool tcliqueSaveFile(
   TCLIQUEGRAPH*    tcliquegraph,       /**< graph data structure */
   const char*      filename,           /**< name of file to create */
   double           scaleval,           /**< value to unscale weights with */
   const char*      probname            /**< name of the problem */
   );

/** gets number of edges in the graph */
extern
int tcliqueGetNEdges(
   TCLIQUEGRAPH*    tcliquegraph        /**< pointer to graph data structure */
   );

/** gets degree of nodes in graph */
extern
int* tcliqueGetDegrees(
   TCLIQUEGRAPH*    tcliquegraph        /**< pointer to graph data structure */
   );

/** gets adjacent nodes of edges in graph */
extern
int* tcliqueGetAdjnodes(
   TCLIQUEGRAPH*    tcliquegraph        /**< pointer to graph data structure */
   );

/** gets pointer to first adjacent edge of given node in graph */
extern
int* tcliqueGetFirstAdjedge(
   TCLIQUEGRAPH*    tcliquegraph,       /**< pointer to graph data structure */
   int              node                /**< given node */
   );

/** gets pointer to last adjacent edge of given node in graph */
extern
int* tcliqueGetLastAdjedge(
   TCLIQUEGRAPH*    tcliquegraph,       /**< pointer to graph data structure */
   int              node                /**< given node */
   );

/* prints graph data structure */
extern
void tcliquePrintGraph(
   TCLIQUEGRAPH*    tcliquegraph        /**< pointer to graph data structure */
   );




/*
 * Interface Methods
 */

/** finds maximum weight clique */
extern
void tcliqueMaxClique(
   TCLIQUE_GETNNODES((*getnnodes)),     /**< user function to get the number of nodes */
   TCLIQUE_GETWEIGHTS((*getweights)),   /**< user function to get the node weights */
   TCLIQUE_ISEDGE   ((*isedge)),        /**< user function to check for existence of an edge */
   TCLIQUE_SELECTADJNODES((*selectadjnodes)), /**< user function to select adjacent edges */
   TCLIQUEGRAPH*    tcliquegraph,       /**< pointer to graph data structure that is passed to graph callbacks */
   TCLIQUE_NEWSOL   ((*newsol)),        /**< user function to call on every new solution */
   TCLIQUEDATA*     tcliquedata,        /**< user data to pass to new solution callback function */
   int*             maxcliquenodes,     /**< pointer to store nodes of the maximum weight clique */
   int*             nmaxcliquenodes,    /**< pointer to store number of nodes in the maximum weight clique */
   WEIGHT*          maxcliqueweight,    /**< pointer to store weight of the maximum weight clique */
   WEIGHT           maxfirstnodeweight, /**< maximum weight of branching nodes in level 0; 0 if not used 
                                         *   for cliques with at least one fractional node) */
   WEIGHT           minweight,          /**< lower bound for weight of generated cliques */
   int              maxtreenodes 	/**< maximum number of nodes of b&b tree */
   );

#endif
