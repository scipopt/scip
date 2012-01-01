/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dijkstra.h
 * @brief  Definitions for Disjkstra's shortest path algorithm using binary heaps
 * @author Thorsten Koch
 * @author Marc Pfetsch
 */

#ifndef DIJSKSTRA_H
#define DIJSKSTRA_H

#ifdef __cplusplus
extern "C" {
#endif

/* declare own bools, if necessary */
#ifndef DIJKSTRA_Bool
#define DIJKSTRA_Bool unsigned int           /**< type used for boolean values */
#endif
#ifndef TRUE
#define TRUE  1                              /**< boolean value TRUE */
#define FALSE 0                              /**< boolean value FALSE */
#endif


#define DIJKSTRA_FARAWAY 0xffffffffffffffffLL   /**< node has distance 'infinity' */
#define DIJKSTRA_UNUSED  0xffffffff             /**< node is unused */

/** graph structure - use consecutive storage for arcs */
struct Dijkstra_graph
{
   unsigned int          nodes;              /**< number of nodes */
   unsigned int*         outbeg;             /**< indices of out-arcs for each node in arcs array */
   unsigned int*         outcnt;             /**< number of out-arcs for each node */
   unsigned int          arcs;               /**< consecutive storage for all arcs */
   unsigned int*         weight;             /**< corresponding weights for all arcs */
   unsigned int*         head;               /**< target nodes for all arcs */
   unsigned long long    min_weight;         /**< total minimal weight */
   unsigned long long    max_weight;         /**< total maximal weight */
};

/** graph structure - use consecutive storage for arcs */
typedef struct Dijkstra_graph Dijkstra_Graph;


/** Check whether the data structures of the graph are valid */
extern
DIJKSTRA_Bool Dijsktra_graphIsValid(
  const Dijkstra_Graph*  G                   /**< directed graph */
  );

/** Dijkstra's algorithm using binary heaps */
extern
unsigned int graph_dijkstra_bh(
  const Dijkstra_Graph*  G,                  /**< directed graph */
  unsigned int           start,              /**< start node */
  unsigned long long*    dist,               /**< node distances */
  unsigned int*          pred,               /**< node predecessors in final shortest path tree */
  unsigned int*          entry,              /**< temporary storage (for each node - must be allocated by user) */
  unsigned int*          order               /**< temporary storage (for each node - must be allocated by user) */
  );

#ifdef __cplusplus
}
#endif

#endif
