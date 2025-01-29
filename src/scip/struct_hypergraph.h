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

/**@file   struct_hypergraph.h
 * @ingroup INTERNALAPI
 * @brief  datastructures hypergraphs
 * @author Matthias Walter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_HYPERGRAPH_H__
#define __SCIP_STRUCT_HYPERGRAPH_H__

#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_hypergraph.h"

#ifdef __cplusplus
extern "C" {
#endif

/** a hypergraph with vertices, edges and overlaps of edge pairs */
struct SCIP_Hypergraph
{
   BMS_BLKMEM*           blkmem;             /**< Block memory for storage. */

   size_t                sizevertexdata;     /**< Size (in bytes) of additional vertex data. */
   size_t                sizeedgedata;       /**< Size (in bytes) of additional edge data. */
   size_t                sizeoverlapdata;    /**< Size (in bytes) of additional overlap data. */

   int                   nvertices;          /**< Number of vertices. */
   int                   nedges;             /**< Number of edges. */
   int                   noverlaps;          /**< Number of overlaps. */

   int                   memvertices;        /**< Number of vertices for which memory is allocated. */
   size_t*               verticesdata;       /**< Array with vertex data. */
   int                   memedges;           /**< Number of edges for which memory is allocated. */

   size_t*               edgesdata;          /**< Array with vertex data. */
   int*                  edgesverticesbeg;   /**< Array with indices of edges' incident vertices. */
   int                   memedgesvertices;   /**< Number of edges' vertices for which memory is allocated. */
   SCIP_HYPERGRAPH_EDGE* edgesvertices;      /**< Array with all edges' vertices. */

   SCIP_Bool             hasvertexedges;     /**< Whether there is a mapping from vertices to incident edges. */
   int                   memverticesedgesbeg;/**< Number of vertices for which memory is allocated for vertices' edges. */
   int*                  verticesedgesbeg;   /**< Array with indices of vertices' incident edges. */
   int                   memverticesedges;   /**< Number of incidences for which memory is allocated. */
   SCIP_HYPERGRAPH_VERTEX* verticesedges;    /**< Array with all vertices' incident edges. */

   SCIP_Bool             hasoverlaps;        /**< Whether overlap sets are known. */
   SCIP_HASHTABLE*       overlaphashtable;   /**< Hashtable for overlap sets. */
   int                   memoverlaps;        /**< Number of overlaps for which memory is allocated. */
   int*                  overlapsverticesbeg;/**< Array with indices of overlaps' vertices. */
   int                   memoverlapsvertices;/**< Number of overlaps' vertices. */
   SCIP_HYPERGRAPH_VERTEX* overlapsvertices; /**< Array with all overlaps' vertices. */
   size_t*               overlapsdata;       /**< Array with overlaps' data. */
   int                   memedgesoverlapsbeg;/**< Memory allocated for \p edgesoverlapsbeg minus 1. */
   int*                  edgesoverlapsbeg;   /**< Array with indices of edges' incident overlaps. */
   int                   memedgesoverlaps;   /**< Number of edges' overlaps for which memory is allocated. */
   SCIP_HYPERGRAPH_OVERLAP* edgesoverlaps;   /**< Array with edges' incident overlaps. */

   SCIP_Bool             hasoverlapsedges;   /**< Whether overlaps' edges are known. */
   int                   memoverlapsedgesbeg;/**< Memory allocated for \p edgesOverlapsSlice minus 1. */
   int*                  overlapsedgesbeg;   /**< Array with indices of overlaps' incident edges. */
   int                   memoverlapsedges;   /**< Memory allocated for \p overlapsEdges. */
   SCIP_HYPERGRAPH_EDGE* overlapsedges;      /**< Array with overlaps' incident edges. */

   SCIP_Bool             hasverticesoverlaps;/**< Whether vertices' overlaps are known. */
   int                   memverticesoverlapsbeg;/**< Memory allocated for \p verticesOverlapsSlice minus 1. */
   int*                  verticesoverlapsbeg;/**< Array with indices of vertices' incident overlaps. */
   int                   memverticesoverlaps;/**< Memory allocated for \p verticesOverlaps. */
   SCIP_HYPERGRAPH_OVERLAP* verticesoverlaps;/**< Array with vertices' incident overlaps. */
};

#ifdef __cplusplus
}
#endif

#endif
