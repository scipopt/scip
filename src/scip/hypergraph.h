/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   hypergraph.h
 * @ingroup INTERNALAPI
 * @brief  Internal methods for dealing with hypergraphs
 * @author Matthias Walter
 *
 * A <em>hypergraph</em> \f$ G \f$ is a pair \f$ (V,E) \f$ of <em>vertices</em> \f$ v \in V \f$ and
 * <em>(hyper-)edges</em> \f$ e \in E \f$ such that \f$ e \subseteq V \f$ holds. We say that a vertex \f$ v \f$ and an
 * edge \f$ e \f$ are <em>incident</em> if \f$ v \in e \f$ holds.
 * A subset \f$ U \subseteq V \f$ is called an <em>overlap set</em> if it is a nontrivial intersection of two distinct
 * edges, i.e., if \f$ U = e \cap f \f$ for \f$ e,f \in E \f$ and \f$ |U| \geq 2 \f$. In this case, we say that
 * \f$ U \f$ and \f$ e \f$ are <em>incident</em> (same for \f$ f \f$. We denote the set of overlap sets by
 * \f$ \mathcal{U} \f$.
 *
 * A \ref SCIP_HYPERGRAPH struct defines such a hypergraph \f$ G = (V,E) \f$, potentially with the set
 * \f$ \mathcal{U} \f$ of overlaps. Efficient access for retrieving elements is provided, e.g., the vertices incident
 * to an edge, the edges incident to a vertex, or the edges incident to an overlap.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HYPERGRAPH_H__
#define __SCIP_HYPERGRAPH_H__

#include "scip/type_retcode.h"
#include "blockmemshell/memory.h"
#include "scip/type_hypergraph.h"

#ifdef NDEBUG
#include "scip/struct_hypergraph.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** @brief creates a hypergraph */
SCIP_RETCODE SCIPhypergraphCreate(
   SCIP_HYPERGRAPH**     phypergraph,        /**< Pointer for storing the hypergraph. */
   BMS_BLKMEM*           blkmem,             /**< Block memory for storage. */
   int                   memvertices,        /**< Upper bound on expected number of vertices. */
   int                   memedges,           /**< Upper bound on expected number of edges. */
   int                   memoverlaps,        /**< Upper bound on expected number of overlaps. */
   int                   memedgesvertices,   /**< Upper bound on expected average size of edges. */
   size_t                sizevertexdata,     /**< Size (in bytes) of additional vertex data. */
   size_t                sizeedgedata,       /**< Size (in bytes) of additional edge data. */
   size_t                sizeoverlapdata     /**< Size (in bytes) of additional overlap data. */
   );

/** @brief frees a hypergraph */
SCIP_RETCODE SCIPhypergraphFree(
   SCIP_HYPERGRAPH**     phypergraph         /**< Pointer to the hypergraph. */
   );

/** @brief clears a hypergraph, deleting all vertices and edges */
SCIP_RETCODE SCIPhypergraphClear(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/** @brief adds a new vertex to the hypergraph */
SCIP_RETCODE SCIPhypergraphAddVertex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX* pvertex,          /**< Pointer for storing the new vertex. */
   SCIP_HYPERGRAPH_VERTEXDATA** pvertexdata  /**< Pointer for returning the vertex's data (may be \c NULL). */
   );

/**
 * @brief adds a new edge to the hypergraph
 *
 * Does not update the vertices' or overlaps' lists of incident edges.
 */
SCIP_RETCODE SCIPhypergraphAddEdge(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   nvertices,          /**< Number of vertices of this edge. */
   SCIP_HYPERGRAPH_VERTEX* vertices,         /**< Array with vertices. */
   SCIP_HYPERGRAPH_EDGE* pedge,              /**< Pointer for storing the new edge. */
   SCIP_HYPERGRAPH_EDGEDATA** pedgedata      /**< Pointer for returning the edge's data (may be \c NULL). */
   );

/** @brief computes each vertex' list of incident edges */
SCIP_RETCODE SCIPhypergraphComputeVerticesEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/**
 * @brief computes all overlaps and stores overlaps' vertices and all edges' overlaps
 *
 * Requires \ref SCIPhypergraphHasVertexEdges to be \c TRUE which results from
 * \ref SCIPhypergraphComputeVerticesEdges.
 */
SCIP_RETCODE SCIPhypergraphComputeOverlaps(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_DECL_HYPERGRAPH_OVERLAP((*handler)), /**< Function to be called once the overlap is found. */
   void*                 userdata            /**< Pointer passed to \p handler. */
   );

/**
 * @brief finds the overlap corresponding to vertex set \p vertices; returns -1 if it is not found
 *
 * Requires \ref SCIPhypergraphHasOverlaps to be \c TRUE which results from \ref SCIPhypergraphComputeOverlaps.
 */
SCIP_RETCODE SCIPhypergraphOverlapFind(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   nvertices,          /**< Number of vertices. */
   SCIP_HYPERGRAPH_VERTEX* vertices,         /**< Sorted array of vertices. */
   SCIP_HYPERGRAPH_OVERLAP* poverlap         /**< Pointer for storing the overlap. */
   );

/**
 * @brief finds the overlap or singleton vertex corresponding to the intersection of edges \p first and \p second
 *
 * Requires \ref SCIPhypergraphHasOverlaps to be \c TRUE which results from \ref SCIPhypergraphComputeOverlaps.
 */
SCIP_RETCODE SCIPhypergraphIntersectEdges(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  first,              /**< First edge to intersect. */
   SCIP_HYPERGRAPH_EDGE  second,             /**< Second edge to intersect. */
   SCIP_HYPERGRAPH_OVERLAP* poverlap,        /**< Pointer for storing the overlap. */
   SCIP_HYPERGRAPH_VERTEX* pvertex           /**< Pointer for storing the vertex. */
   );

/**
 * @brief computes all overlaps' lists of incident edges
 *
 * Requires \ref SCIPhypergraphHasOverlaps to be \c TRUE which results from \ref SCIPhypergraphComputeOverlaps.
 */
SCIP_RETCODE SCIPhypergraphComputeOverlapsEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/**
 * @brief computes all vertices' lists of incident overlaps
 *
 * Requires \ref SCIPhypergraphHasOverlaps.
 */
SCIP_RETCODE SCIPhypergraphComputeVerticesOverlaps(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/** @brief returns whether overlaps \p overlap1 and \p overlap2 are disjoint */
SCIP_Bool SCIPhypergraphOverlapsDisjoint(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap1,         /**< First overlap. */
   SCIP_HYPERGRAPH_OVERLAP overlap2          /**< Second overlap. */
   );

/** @brief asserts that the hypergraph data structures are valid */
SCIP_Bool SCIPhypergraphIsValid(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   FILE*                 file                /**< A file stream to print problems to (can be \c NULL). */
   );

/** @brief initializes a hypergraph iterator's internal memory */
SCIP_RETCODE SCIPhypergraphIterInit(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph of this iterator. */
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator to be initialized. */
   );

/** @brief frees a hypergraph iterator's internal memory */
void SCIPhypergraphIterClear(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph of this iterator. */
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator to be cleared. */
   );

/**
 * @brief initializes the iterator to the first adjacent edge of \p base
 *
 * \p findoverlaps indicates whether the overlap corresponding to each edge intersection shall be determined. This
 * requires \ref SCIPhypergraphHasOverlaps.
 */
void SCIPhypergraphIterStart(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph of this iterator. */
   SCIP_HYPERGRAPH_ITER* iterator,           /**< Iterator to be initialized. */
   SCIP_HYPERGRAPH_EDGE  base,               /**< Base edge. */
   unsigned int          minoverlapsize,     /**< Minimum size of the overlap. */
   SCIP_Bool             onlylater,          /**< Whether to only consider edges greater than \p base. */
   SCIP_Bool             findoverlaps        /**< Whether to compute the overlap sets. */
   );

/** @brief returns whether the iterator is valid */
SCIP_Bool SCIPhypergraphIterValid(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   );

/** @brief initializes the iterator to the first adjacent edge of \p base */
void SCIPhypergraphIterNext(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph of this iterator. */
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   );

/** @brief returns the base edge */
SCIP_HYPERGRAPH_EDGE SCIPhypergraphIterBase(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   );

/** @brief returns the current adjacent edge */
SCIP_HYPERGRAPH_EDGE SCIPhypergraphIterAdjacent(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   );

/** @brief returns the minimum vertex in the intersection of the base and the current adjacent edge */
SCIP_HYPERGRAPH_VERTEX SCIPhypergraphIterMinVertex(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   );

/**
 * @brief returns the overlap for the intersection of the base and the current adjacent edge
 *
 * If the intersection of the two edges has only one element, then -1 is returned, and if overlap information was not
 * requested in \ref SCIPhypergraphIterStart, then -2 is returned.
 */
SCIP_HYPERGRAPH_OVERLAP SCIPhypergraphIterOverlap(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   );

/**
 * @brief returns whether vertices' incident edges are known.
 *
 * Use \ref SCIPhypergraphComputeVerticesEdges to compute them.
 */

SCIP_Bool SCIPhypergraphHasVertexEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/**
 * @brief returns whether edges' overlaps and overlaps' vertices are known.
 *
 * Use \ref SCIPhypergraphComputeOverlaps to compute them.
 */
SCIP_Bool SCIPhypergraphHasOverlaps(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/**
 * @brief returns whether overlaps' incident edges are known.
 *
 * Use \ref SCIPhypergraphComputeOverlapsEdges to compute them.
 */
SCIP_Bool SCIPhypergraphHasOverlapsEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/**
 * @brief returns whether vertices' incident overlaps are known
 *
 * Use \ref SCIPhypergraphComputeOverlaps to compute them.
 */
SCIP_Bool SCIPhypergraphHasVertexOverlaps(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/** @brief returns the number of vertices */
int SCIPhypergraphGetNVertices(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/** @brief returns the number of edges */
int SCIPhypergraphGetNEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/** @brief returns the hypergraph's block memory structure */
BMS_BLKMEM* SCIPhypergraphBlkmem(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/** @brief returns the number of overlaps */
int SCIPhypergraphGetNOverlaps(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/** @brief returns the number of vertex-edge incidences. */
int SCIPhypergraphGetNIncidences(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   );

/** @brief returns additional data of \p vertex */
SCIP_HYPERGRAPH_VERTEXDATA* SCIPhypergraphVertexData (
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex             /**< A vertex. */
   );

/** @brief returns additional data of \p edge */
SCIP_HYPERGRAPH_EDGEDATA* SCIPhypergraphEdgeData(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< Edge. */
   );

/** @brief returns the number of vertices of \p edge */
int SCIPhypergraphEdgeSize(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< Edge. */
   );

/**
 * @brief returns the array of vertices of \p edge
 *
 * The length of the array is \ref SCIPhypergraphEdgeSize.
 */

SCIP_HYPERGRAPH_VERTEX* SCIPhypergraphEdgeVertices(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< Edge. */
   );

/** @brief returns an index for the first edge incident to \p vertex */
int SCIPhypergraphVertexEdgesFirst(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex             /**< A vertex. */
   );

/** @brief returns an index beyond the last edge incident to \p vertex */
int SCIPhypergraphVertexEdgesBeyond(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex             /**< A vertex. */
   );

/**
 * @brief returns the edge corresponding to \p index that is incident to a vertex
 *
 * See \ref SCIPhypergraphVertexEdgesFirst and \ref SCIPhypergraphVertexEdgesBeyond to obtain such indices for a vertex.
 */
SCIP_HYPERGRAPH_EDGE SCIPhypergraphVertexEdgesGetAtIndex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   index               /**< Index. */
   );

/** @brief returns additional data of \p overlap */
SCIP_HYPERGRAPH_OVERLAPDATA* SCIPhypergraphOverlapData(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< Overlap. */
   );

/**
 * @brief returns the number of vertices of \p overlap
 *
 * Requires \ref SCIPhypergraphHasOverlaps to be \c TRUE which results from \ref SCIPhypergraphComputeOverlaps.
 */
int SCIPhypergraphOverlapSize(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< Overlap. */
   );

/**
 * @brief returns the array of sorted vertices of \p overlap
 *
 * The length of the array is \ref SCIPhypergraphOverlapSize.
 * Requires \ref SCIPhypergraphHasOverlaps to be \c TRUE which results from \ref SCIPhypergraphComputeOverlaps.
 */
SCIP_HYPERGRAPH_VERTEX* SCIPhypergraphOverlapVertices(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< Overlap. */
   );

/** @brief returns an index for the first overlap incident to \p edge */
int SCIPhypergraphEdgesOverlapsFirst(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< Edge. */
   );

/** @brief returns an index beyond the last overlap incident to \p edge */
int SCIPhypergraphEdgesOverlapsBeyond(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< Edge. */
   );

/**
 * @brief returns the overlap corresponding to \p index that is incident to an edge
 *
 * See \ref SCIPhypergraphEdgesOverlapsFirst and \ref SCIPhypergraphEdgesOverlapsBeyond to obtain such indices for an
 * edge.
 */
SCIP_HYPERGRAPH_OVERLAP SCIPhypergraphEdgesOverlapsGetAtIndex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   index               /**< Index. */
   );

/** @brief returns an index for the first edge incident to \p overlap */
int SCIPhypergraphOverlapsEdgesFirst(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< Overlap. */
   );

/** @brief returns an index beyond the last edge incident to \p overlap */
int SCIPhypergraphOverlapsEdgesBeyond(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< Overlap. */
   );

/**
 * @brief returns the edge corresponding to \p index that is incident to an overlap
 *
 * See \ref SCIPhypergraphOverlapsEdgesFirst and \ref SCIPhypergraphOverlapsEdgesBeyond to obtain such indices for an
 * overlap.
 */
SCIP_HYPERGRAPH_OVERLAP SCIPhypergraphOverlapsEdgesGetAtIndex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   index               /**< Index. */
   );

/** @brief returns an index for the first overlap containing \p vertex */
int SCIPhypergraphVertexOverlapsFirst(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex            /**< A vertex. */
   );

/** @brief returns an index beyond the last overlap incident to \p vertex */
int SCIPhypergraphVertexOverlapsBeyond(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex             /**< A vertex. */
   );

/**
 * @brief returns the overlap corresponding to \p index that is incident to a vertex
 *
 * See \ref SCIPhypergraphVertexOverlapsFirst and \ref SCIPhypergraphVertexOverlapsBeyond to obtain such indices for a
 * vertex.
 */
SCIP_HYPERGRAPH_OVERLAP SCIPhypergraphVertexOverlapsGetAtIndex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   index               /**< Index. */
   );


#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and speed up
 * the algorithms.
 */

#define SCIPhypergraphGetNVertices(hypergraph) ((hypergraph)->nvertices)
#define SCIPhypergraphGetNEdges(hypergraph) ((hypergraph)->nedges)
#define SCIPhypergraphBlkmem(hypergraph) ((hypergraph)->blkmem)
#define SCIPhypergraphGetNOverlaps(hypergraph) ((hypergraph)->noverlaps)
#define SCIPhypergraphGetNIncidences(hypergraph) ((hypergraph)->edgesverticesbeg[(hypergraph)->nedges])
#define SCIPhypergraphVertexData(hypergraph,vertex) \
   ((SCIP_HYPERGRAPH_VERTEXDATA*)((hypergraph)->verticesdata + ((vertex) * ((hypergraph)->sizevertexdata) / sizeof(*((hypergraph)->verticesdata))  )))
#define SCIPhypergraphEdgeData(hypergraph,edge) \
   ((SCIP_HYPERGRAPH_EDGEDATA*)((hypergraph)->edgesdata + ((edge) * (hypergraph)->sizeedgedata / sizeof(*((hypergraph)->edgesdata)) )))
#define SCIPhypergraphOverlapData(hypergraph,overlap) ((SCIP_HYPERGRAPH_OVERLAPDATA*)((hypergraph)->overlapsdata \
   + ((overlap) * (hypergraph)->sizeoverlapdata / sizeof(*((hypergraph)->overlapsdata)) )))
#define SCIPhypergraphEdgeSize(hypergraph,edge) ((hypergraph)->edgesverticesbeg[(edge) + 1] \
   - (hypergraph)->edgesverticesbeg[edge])
#define SCIPhypergraphEdgeVertices(hypergraph, edge) (&(hypergraph)->edgesvertices[(hypergraph)->edgesverticesbeg[edge]])
#define SCIPhypergraphHasVertexEdges(hypergraph) ((hypergraph)->hasvertexedges)
#define SCIPhypergraphVertexEdgesFirst(hypergraph,vertex) ((hypergraph)->verticesedgesbeg[vertex])
#define SCIPhypergraphVertexEdgesBeyond(hypergraph,vertex) ((hypergraph)->verticesedgesbeg[(vertex) + 1])
#define SCIPhypergraphVertexEdgesGetAtIndex(hypergraph,index) ((hypergraph)->verticesedges[index])
#define SCIPhypergraphHasOverlaps(hypergraph) ((hypergraph)->hasoverlaps)
#define SCIPhypergraphOverlapSize(hypergraph,overlap) ((hypergraph)->overlapsverticesbeg[(overlap) + 1] \
   - (hypergraph)->overlapsverticesbeg[overlap])
#define SCIPhypergraphOverlapVertices(hypergraph,overlap) \
   (&(hypergraph)->overlapsvertices[(hypergraph)->overlapsverticesbeg[overlap]])
#define SCIPhypergraphEdgesOverlapsFirst(hypergraph,edge) ((hypergraph)->edgesoverlapsbeg[edge])
#define SCIPhypergraphEdgesOverlapsBeyond(hypergraph,edge) ((hypergraph)->edgesoverlapsbeg[(edge) + 1])
#define SCIPhypergraphEdgesOverlapsGetAtIndex(hypergraph,index) ((hypergraph)->edgesoverlaps[index])
#define SCIPhypergraphHasOverlapsEdges(hypergraph) ((hypergraph)->hasoverlapsedges)
#define SCIPhypergraphOverlapsEdgesFirst(hypergraph,overlap) ((hypergraph)->overlapsedgesbeg[overlap])
#define SCIPhypergraphOverlapsEdgesBeyond(hypergraph,overlap) ((hypergraph)->overlapsedgesbeg[(overlap) + 1])
#define SCIPhypergraphOverlapsEdgesGetAtIndex(hypergraph,index) ((hypergraph)->overlapsedges[index])
#define SCIPhypergraphHasVertexOverlaps(hypergraph) ((hypergraph)->hasverticesoverlaps)
#define SCIPhypergraphVertexOverlapsFirst(hypergraph,vertex) ((hypergraph)->verticesoverlapsbeg[vertex])
#define SCIPhypergraphVertexOverlapsBeyond(hypergraph,vertex) ((hypergraph)->verticesoverlapsbeg[(vertex) + 1])
#define SCIPhypergraphVertexOverlapsGet(hypergraph,index) ((hypergraph)->verticesoverlaps[index])

#endif /* NDEBUG */

#ifdef __cplusplus
}
#endif

#endif
