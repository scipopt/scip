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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   gminucut.h
 * @brief  generator for global minimum cuts in undirected graphs
 * @author Georg Skorobohatyj
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GMINUCUT_H__
#define __GMINUCUT_H__

#include "objscip/objscip.h"


typedef struct GraphNode
{
   int         id;
   double      x;
   double      y;
   struct GraphEdge  *first_edge;   /**< in cyclic list of incident edges */
   struct GraphEdge  *scan_ptr;     /**< next edge to be scanned when node will be visited again */
   
   /* subsequent entries for use by gmincut */
   int         dist;
   double       excess;
   struct GraphNode  *bfs_link;     /**< for one way BFS working queue */
   struct GraphNode  *stack_link;   /**< for stack of active node */
   struct GraphNode  *left_link;    /**< for doubly linked lists */
   struct GraphNode  *right_link;   /**< of dormant and W nodes  */
   Bool         unmarked;      /**< while BFS in progress */
   Bool         in_S;          /**< in set of source nodes  */
   Bool         in_W;          /**< in set of W-valid nodes */
   Bool         partition;        /**< final mark for placement of node in one of the cut shores */
} GRAPHNODE;

typedef struct GraphEdge
{ 
   GRAPHNODE         *adjac; /**< pointer to adjacent node */
   struct GraphEdge  *next;  /**< in incidence list of node from which edge is emanating */
   struct GraphEdge  *back;  /**< pointer to reverse edge */
   double       cap;
   double       rcap;        /**< residual capacity used in gmincut */
   double       length;      /**< length of the edge measured by some fixed metric */
   VAR*          var;
} GRAPHEDGE;

typedef struct Graph
{ 
   int nuses;
   int    nnodes;
   GRAPHNODE    *nodes;
   int    nedges;
   GRAPHEDGE    *edges;
} GRAPH;


extern
Bool create_graph(int n, int m, GRAPH** gr);

extern
void capture_graph(GRAPH* gr);
   
extern
void release_graph(GRAPH** gr);
   
extern
Bool gmincut(GRAPH *gr, double *mincap, long *n_shore);


#endif
