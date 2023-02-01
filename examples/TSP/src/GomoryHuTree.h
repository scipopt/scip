/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   GomoryHuTree.h
 * @brief  generator for global cuts in undirected graphs
 * @author Georg Skorobohatyj
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GOMORYHUTREE_H__
#define __GOMORYHUTREE_H__

#include "objscip/objscip.h"


/** a node in the graph */
typedef struct GraphNode
{
   int                   id;                 /**< number of the node*/
   int                   dist;               /**< distances used in push-relabel */

   double                x;                  /**< 2D-coordinate in some metric */
   double                y;                  /**< second coordinate */
   double                excess;             /**< excess of node used in push-relabel */
   double                mincap;             /**< capacity of minimum cut between node and parent in GH cut tree */

   SCIP_Bool             unmarked;           /**< while BFS in progress */
   SCIP_Bool             alive;              /**< marks alive (active) nodes in push-relabel */

   struct GraphEdge*     first_edge;         /**< in list of incident edges */
   struct GraphEdge*     scan_ptr;           /**< next edge to be scanned when node will be visited again */

   struct GraphNode*     bfs_link;           /**< for one way BFS working queue */
   struct GraphNode*     stack_link;         /**< for stack of active node */
   struct GraphNode*     parent;             /**< pointer of Gomory-Hu cut tree */
} GRAPHNODE;


/** an edge in the graph */
typedef struct GraphEdge
{
   double                cap;                /**< capacity used in maxflow */
   double                rcap;               /**< residual capacity used in maxflow */
   double                length;             /**< length of the edge measured by some fixed metric */

   struct GraphEdge*     next;               /**< in incidence list of node from which edge is emanating */
   struct GraphEdge*     back;               /**< pointer to reverse edge */

   GRAPHNODE*            adjac;              /**< pointer to adjacent node */

   SCIP_VAR*             var;                /**< variable associated to edge */
} GRAPHEDGE;


/** undirected graph */
typedef struct Graph
{
   int                   nuses;              /**< usage counter */
   int                   nnodes;             /**< number of nodes of the graph */
   int                   nedges;             /**< number of edges */
   int                   nedgesnonzero;      /**< nonzero edges (not currently used) */

   GRAPHNODE*            nodes;              /**< array containing the nodes of the graph */
   GRAPHEDGE*            edges;              /**< array containing all halfedges (thus, it's size is two times nedges) */
} GRAPH;


/** create a graph */
SCIP_Bool create_graph(
   int                   n,                  /**< number of nodes */
   int                   m,                  /**< number of edges */
   GRAPH**               gr                  /**< pointer to store graph */
   );

/** capture graph */
void capture_graph(
   GRAPH*                gr                  /**< graph */
   );

/** release graph */
void release_graph(
   GRAPH**               gr                  /**< graph */
   );

/** Determines Gomory/Hu cut tree for input graph with capacitated edges */
SCIP_Bool ghc_tree(
   GRAPH*                gr,                 /**< graph */
   SCIP_Bool**           cuts,               /**< array of arrays to store cuts */
   int*                  ncuts,              /**< pointer to store number of cuts */
   double                minviol             /**< minimal violation of a cut to be returned */
   );

#endif
