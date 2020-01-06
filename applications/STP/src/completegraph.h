/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   completegraph.c
 * @brief  includes complete graph methods, in particular for MST calculation
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_COMPLETEGRAPH_H_
#define APPLICATIONS_STP_SRC_COMPLETEGRAPH_H_

#define CGRAPH_EDGECOST_UNDEFINED_VALUE -1.0

#include "scip/scip.h"
#include "graph.h"

/** complete, undirected graph storage */
typedef struct complete_graph
{
   SCIP_Real*            edgecosts;          /**< edge cost array; of size nnodes_max * (nnodes_max - 1) */
   SCIP_Real*            adjedgecosts;       /**< adjacency cost array of size nnodes_max + 1, to be filled in by user */
   int*                  nodeids;            /**< node ids; of size nnodes_max */
   SCIP_Bool*            node_has_adjcosts;  /**< does node have adjacency costs? */
   int                   nnodes_max;         /**< maximum number of nodes */
   int                   nnodes_curr;        /**< current number of nodes */
   int                   nnodes_active;      /**< current number of nodes for which adjacency costs have been added */
} CGRAPH;


/** buffer that allows to (partially) restore a complete graph node
 * todo need something more special here that takes dfs depth and a node identifier (which) as input applies corresponding
 * part to graph...
 * maybe whenever you store something you get an identifier back??? */
typedef struct complete_graph_node_buffer
{
   SCIP_Real*            adjedgecosts;       /**< adjacency cost array for nodes */
   int*                  nodeids;            /**< corresponding node ids (only needed for debugging) */
   int                   nnodes;             /**< number of nodes that are currently stored */
   int                   dfsdepth;           /**< current number of nodes */
} CNBUFF;


/** complete, undirected graph MST storage; for computing an MST on a CGRAPH */
typedef struct complete_mst
{
   DHEAP*                heap;               /**< heap needed for MST computation */
   SCIP_Real*            dist;               /**< distance array of size nnodes_max */
   int*                  predecessors;       /**< predecessor array of size nnodes_max */
   SCIP_Real             mstobj;             /**< objective of current MST */
   int                   nnodes_max;         /**< maximum number of nodes */
} CMST;


/* methods for the complete graph */
SCIP_Bool cgraph_valid(const CGRAPH*);
SCIP_Bool cgraph_isEmpty(const CGRAPH*);
SCIP_Bool cgraph_idsInSync(const CGRAPH*, const int*, int);
SCIP_Bool cgraph_idIsContained(const CGRAPH*, int);
SCIP_RETCODE cgraph_init(SCIP*, CGRAPH**, int);
void cgraph_free(SCIP*, CGRAPH**);
void cgraph_clean(CGRAPH*);
SCIP_Bool cgraph_node_hasAdjCosts(const CGRAPH*, int);
void cgraph_node_append(CGRAPH*, int);
void cgraph_node_applyMinAdjCosts(CGRAPH*, int, int);
void cgraph_node_repositionTop(CGRAPH*, int);
void cgraph_node_deleteTop(CGRAPH*);
void cgraph_node_delete(CGRAPH*, int);
int cgraph_node_getTopId(const CGRAPH*);
SCIP_Real cgraph_edge_getCost(const CGRAPH*, int, int);

/*  methods for the corresponding MST structure */
SCIP_Bool cmst_isSync(const CGRAPH*, const CMST*);
SCIP_RETCODE cmst_init(SCIP*, CMST**, int);
void cmst_free(SCIP*, CMST**);
void cmst_computeMst(const CGRAPH*, int, CMST*);


#endif /* APPLICATIONS_STP_SRC_COMPLETEGRAPH_H_ */
