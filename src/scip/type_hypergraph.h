/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   type_hypergraph.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for hypergraphs
 * @author Matthias Walter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_HYPERGRAPH_H__
#define __SCIP_TYPE_HYPERGRAPH_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "scip/def.h"

/** a hypergraph with vertices, edges and overlaps of edge pairs */
typedef struct SCIP_Hypergraph SCIP_HYPERGRAPH;

/** data for iterating over adjacent edges. */
typedef struct SCIP_Hypergraph_Iter SCIP_HYPERGRAPH_ITER;

/** locally defined data for each vertex in a hypergraph */
typedef struct SCIP_Hypergraph_NodeData SCIP_HYPERGRAPH_VERTEXDATA;

/** locally defined data for each edge in a hypergraph */
typedef struct SCIP_Hypergraph_EdgeData SCIP_HYPERGRAPH_EDGEDATA;

/** locally defined data for each overlap set in a hypergraph */
typedef struct SCIP_Hypergraph_OverlapData SCIP_HYPERGRAPH_OVERLAPDATA;

/** vertex in a hypergraph */
typedef int SCIP_HYPERGRAPH_VERTEX;

/** edge in a hypergraph */
typedef int SCIP_HYPERGRAPH_EDGE;

/** overlap set in a hypergraph */
typedef int SCIP_HYPERGRAPH_OVERLAP;


/** Called by \ref SCIPhypergraphOverlapFind, \ref SCIPhypergraphIntersectEdges and \ref SCIPhypergraphComputeOverlaps
 * whenever a new overlap set is created or an existing overlap is found.
 */
#define SCIP_DECL_HYPERGRAPH_OVERLAP(x) SCIP_RETCODE x (SCIP_HYPERGRAPH* hypergraph, SCIP_HYPERGRAPH_OVERLAP overlap, \
  SCIP_HYPERGRAPH_OVERLAPDATA* data, SCIP_HYPERGRAPH_EDGE first, SCIP_HYPERGRAPH_EDGE second, SCIP_Bool created, \
  void* userdata)


/** masks to control the iteration over adjacent edges. */
enum SCIP_Hypergraph_IterCtrl
{
   SCIP_HYPERGRAPH_ITERCTRL_MINOVERLAP = 255,    /**< Mask for minimum required size of edge intersections.  */
   SCIP_HYPERGRAPH_ITERCTRL_ONLYLATER = 256,     /**< Whether to only consider edges with larger index than the base. */
   SCIP_HYPERGRAPH_ITERCTRL_FINDOVERLAPS = 512   /**< Whether to compute the corresponding overlaps. */
};
typedef enum SCIP_Hypergraph_IterCtrl SCIP_HYPERGRAPH_ITERCTRL; /**< controls the iteration over adjacent edges. */

/** data for iterating over adjacent edges. */
struct SCIP_Hypergraph_Iter
{
   SCIP_HYPERGRAPH_EDGE  base;               /**< Base edge for iteration. */
   int                   vertexidx;          /**< Index of incident vertex w.r.t. to base; initially -1; invalid if -2 */
   SCIP_HYPERGRAPH_VERTEX minvertex;         /**< Incident vertex. */
   int                   edgeidx;            /**< Index of the adjacent edge w.r.t. to the vertex. */
   SCIP_HYPERGRAPH_EDGE  adjacent;           /**< Adjacent edge. */
   int                   ncommonvertices;    /**< Number of common vertices. */
   int                   sizecommonvertices; /**< Memory allocated for common vertices. */
   SCIP_HYPERGRAPH_VERTEX* commonvertices;   /**< Array of common vertices. */
   SCIP_HYPERGRAPH_OVERLAP overlap;          /**< Overlap set of the intersection, if available; -2 if disabled. */
   unsigned int          minoverlapsize : 6; /**< Minimum size of the overlap. */
   unsigned int          onlylater : 1;      /**< Whether to only consider edges greater than the base edge. */
   unsigned int          findoverlaps : 1;   /**< Whether to compute the overlap sets. */
};

#ifdef __cplusplus
}
#endif

#endif

