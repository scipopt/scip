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

/**@file   hypergraph.c
 * @ingroup OTHER_CFILES
 * @brief  internal methods for dealing with hypergraphs
 * @author Matthias Walter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/hypergraph.h"
#include "scip/misc.h"
#include "scip/pub_message.h" /* For SCIPdebugMessage */

#ifndef NDEBUG
#include "scip/struct_hypergraph.h"
#endif

/** @brief enlarges vertex arrays to have capacity for at least \p required vertices */
static
SCIP_RETCODE ensureNumVertices(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   required            /**< Required vertex capacity. */
   )
{
   assert(hypergraph);
   assert(required >= 0);

   if( hypergraph->memvertices < required )
   {
      int newcapacity;
      newcapacity = MAX(256, hypergraph->memvertices);
      while( newcapacity < required )
         newcapacity *= 2;

      assert( hypergraph->memvertices >= 0 );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesdata,
            hypergraph->memvertices * hypergraph->sizevertexdata / sizeof(*(hypergraph->verticesdata)),
            newcapacity * hypergraph->sizevertexdata / sizeof(*(hypergraph->verticesdata))) ); /*lint --e{737}*/
      hypergraph->memvertices = newcapacity;
      SCIPdebugMessage("ensuring memory for %d vertices.\n", newcapacity);
   }

   return SCIP_OKAY;
}

/** @brief enlarges edge arrays to have capacity for at least \p required edges */
static
SCIP_RETCODE ensureNumEdges(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   required            /**< Required edge capacity. */
   )
{
   assert(hypergraph);
   assert(required >= 0);

   if( hypergraph->memedges < required )
   {
      int newcapacity;
      newcapacity = MAX(256, hypergraph->memedges);
      while( newcapacity < required )
         newcapacity *= 2;

      assert( hypergraph->memedges >= 0 );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesdata,
            hypergraph->memedges * hypergraph->sizeedgedata / sizeof(*(hypergraph->edgesdata)),
            newcapacity * hypergraph->sizeedgedata / sizeof(*(hypergraph->edgesdata))) ); /*lint --e{737}*/

      if( hypergraph->edgesverticesbeg )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesverticesbeg,
               hypergraph->memedges + 1, newcapacity + 1) );
      }
      else
      {
         assert(hypergraph->memedges == 0);
         SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesverticesbeg, newcapacity + 1) );
         hypergraph->edgesverticesbeg[0] = 0;
      }
      hypergraph->memedges = newcapacity;
      SCIPdebugMessage("ensuring memory for %d edges.\n", newcapacity);
   }

   return SCIP_OKAY;
}

/** @brief enlarges arrays for incident vertex/edge pairs to have capacity for at least \p required */
static
SCIP_RETCODE ensureNumEdgesVertices(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   required            /**< Required incident vertex/edge capacity. */
   )
{
   assert(hypergraph);
   assert(required >= 0);

   if( hypergraph->memedgesvertices < required )
   {
      int newcapacity;
      newcapacity = MAX(256, hypergraph->memedgesvertices);
      while( newcapacity < required )
         newcapacity *= 2;

      SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesvertices,
            hypergraph->memedgesvertices, newcapacity) );
      hypergraph->memedgesvertices = newcapacity;
      SCIPdebugMessage("ensuring memory for %d vertex/edge pairs.\n", newcapacity);
   }

   return SCIP_OKAY;
}

/** @brief enlarges overlap arrays to have capacity for at least \p required overlaps */
static
SCIP_RETCODE ensureNumOverlaps(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   required            /**< Required overlap capacity. */
   )
{
   assert(hypergraph);
   assert(required >= 0);

   if( hypergraph->memoverlaps < required )
   {
      int newcapacity;
      newcapacity = MAX(256, hypergraph->memoverlaps);
      while( newcapacity < required )
         newcapacity *= 2;

      if( hypergraph->overlapsverticesbeg )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsverticesbeg,
               hypergraph->memoverlaps + 1, newcapacity + 1) );
      }
      else
      {
         assert(hypergraph->memoverlaps == 0);
         SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsverticesbeg, newcapacity + 1) );
      }
      if( hypergraph->overlapsdata )
      {
         assert( hypergraph->memoverlaps >= 0 );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsdata,
               hypergraph->memoverlaps * hypergraph->sizeoverlapdata / sizeof(*(hypergraph->overlapsdata)),
               newcapacity * hypergraph->sizeoverlapdata / sizeof(*(hypergraph->overlapsdata))) ); /*lint --e{737}*/
      }
      else
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsdata,
               newcapacity * hypergraph->sizeoverlapdata / sizeof(*(hypergraph->overlapsdata))) ); /*lint --e{737}*/
      }
      hypergraph->memoverlaps = newcapacity;
      SCIPdebugMessage("ensuring memory for %d overlaps.\n", newcapacity);
   }

   return SCIP_OKAY;
}

/** @brief enlarges overlapVertices array to have capacity for at least \p required overlaps' vertices */
static
SCIP_RETCODE ensureNumOverlapsVertices(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   required            /**< Required overlapVertices capacity. */
   )
{
   assert(hypergraph);
   assert(required >= 0);

   if( hypergraph->memoverlapsvertices < required )
   {
      int newcapacity;
      newcapacity = MAX(256, hypergraph->memoverlapsvertices);
      while( newcapacity < required )
         newcapacity *= 2;

      SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsvertices,
            hypergraph->memoverlapsvertices, newcapacity) );
      hypergraph->memoverlapsvertices = newcapacity;
      SCIPdebugMessage("ensuring memory for %d overlaps' vertices.\n", newcapacity);
   }

   return SCIP_OKAY;
}

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
   )
{
   SCIP_HYPERGRAPH* hypergraph;

   assert(phypergraph);
   assert(blkmem);

   assert(memvertices >= 0);
   assert(memedges >= 0);
   assert(memoverlaps >= 0);
   assert(memedgesvertices > 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, phypergraph) );
   hypergraph = *phypergraph;
   hypergraph->blkmem = blkmem;
   hypergraph->memvertices = 0;
   hypergraph->memedges = 0;
   hypergraph->memoverlaps = 0;
   hypergraph->memedgesvertices = 0;
   hypergraph->memverticesedges = 0;
   hypergraph->memverticesedgesbeg = 0;
   hypergraph->memoverlapsvertices = 0;
   hypergraph->memedgesoverlaps = 0;
   hypergraph->memedgesoverlapsbeg = 0;
   hypergraph->sizevertexdata = sizevertexdata;
   hypergraph->sizeedgedata = sizeedgedata;
   hypergraph->sizeoverlapdata = sizeoverlapdata;

   hypergraph->verticesdata = NULL;
   hypergraph->verticesedgesbeg = NULL;
   SCIP_CALL( ensureNumVertices(hypergraph, memvertices) );
   hypergraph->edgesdata = NULL;
   hypergraph->edgesverticesbeg = NULL;
   SCIP_CALL( ensureNumEdges(hypergraph, memedges) );
   hypergraph->overlapsdata = NULL;
   hypergraph->overlapsverticesbeg = NULL;
   hypergraph->overlapsvertices = NULL;
   SCIP_CALL( ensureNumOverlaps(hypergraph, memoverlaps) );
   hypergraph->edgesvertices = NULL;
   hypergraph->verticesedges = NULL;
   hypergraph->edgesoverlaps = NULL;
   hypergraph->edgesoverlapsbeg = NULL;
   SCIP_CALL( ensureNumEdgesVertices(hypergraph, memedgesvertices * memedges) );
   hypergraph->overlaphashtable = NULL;

   hypergraph->memoverlapsedgesbeg = 0;
   hypergraph->overlapsedgesbeg = NULL;
   hypergraph->memoverlapsedges = 0;
   hypergraph->overlapsedges = NULL;

   hypergraph->memverticesoverlapsbeg = 0;
   hypergraph->verticesoverlapsbeg = NULL;
   hypergraph->memverticesoverlaps = 0;
   hypergraph->verticesoverlaps = NULL;

   SCIP_CALL( SCIPhypergraphClear(hypergraph) );

   return SCIP_OKAY;
}

/** @brief frees a hypergraph */
SCIP_RETCODE SCIPhypergraphFree(
   SCIP_HYPERGRAPH**     phypergraph         /**< Pointer to the hypergraph. */
   )
{
   SCIP_HYPERGRAPH* hypergraph;

   assert(phypergraph);

   hypergraph = *phypergraph;
   if( hypergraph == NULL )
      return SCIP_OKAY;

   BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesdata,
      hypergraph->memvertices * hypergraph->sizevertexdata / sizeof(*hypergraph->verticesdata) );
   BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesdata,
      hypergraph->memedges * hypergraph->sizeedgedata / sizeof(*hypergraph->edgesdata));
   if( hypergraph->edgesverticesbeg )
      BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesverticesbeg, hypergraph->memedges + 1);
   if( hypergraph->edgesvertices )
      BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesvertices, hypergraph->memedgesvertices);

   if( hypergraph->hasvertexedges )
   {
      if( hypergraph->verticesedges )
         BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesedges, hypergraph->memverticesedges);
      if( hypergraph->verticesedgesbeg )
         BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesedgesbeg, hypergraph->memvertices + 1);
   }

   if( hypergraph->overlapsdata )
   {
      BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsdata,
         hypergraph->memoverlaps * hypergraph->sizeoverlapdata / sizeof(*hypergraph->overlapsdata));
   }

   if( hypergraph->overlapsverticesbeg )
      BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsverticesbeg, hypergraph->memoverlaps + 1);

   if( hypergraph->hasoverlaps )
   {
      if( hypergraph->overlapsvertices )
         BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsvertices, hypergraph->memoverlapsvertices);
      if( hypergraph->overlaphashtable )
         SCIPhashtableFree(&hypergraph->overlaphashtable);
      if( hypergraph->edgesoverlapsbeg )
      {
         BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesoverlapsbeg,
            hypergraph->memedgesoverlapsbeg + 1);
      }
      if( hypergraph->edgesoverlaps )
         BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesoverlaps, hypergraph->memedgesoverlaps);
   }

   if( hypergraph->hasoverlapsedges )
   {
      if( hypergraph->overlapsedges )
         BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsedges, hypergraph->memoverlapsedges);
      if( hypergraph->overlapsedgesbeg )
      {
         BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsedgesbeg,
            hypergraph->memoverlapsedgesbeg + 1);
      }
   }

   if( hypergraph->hasverticesoverlaps )
   {
      BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesoverlaps, hypergraph->memverticesoverlaps);
      if( hypergraph->verticesoverlapsbeg )
      {
         BMSfreeBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesoverlapsbeg,
            hypergraph->memverticesoverlapsbeg + 1);
      }
   }

   BMSfreeBlockMemory(hypergraph->blkmem, phypergraph);

   return SCIP_OKAY;
}

/** @brief clears a hypergraph, deleting all vertices and edges */
SCIP_RETCODE SCIPhypergraphClear(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph);

   hypergraph->nvertices = 0;
   hypergraph->nedges = 0;
   hypergraph->noverlaps = 0;
   if( hypergraph->edgesverticesbeg )
      hypergraph->edgesverticesbeg[0] = 0;
   hypergraph->hasvertexedges = FALSE;
   hypergraph->hasoverlaps = FALSE;

   if( hypergraph->overlaphashtable )
      SCIPhashtableFree(&hypergraph->overlaphashtable);

   hypergraph->hasoverlapsedges = FALSE;
   hypergraph->hasverticesoverlaps = FALSE;

   return SCIP_OKAY;
}


/** @brief adds a new vertex to the hypergraph */
SCIP_RETCODE SCIPhypergraphAddVertex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX* pvertex,          /**< Pointer for storing the new vertex. */
   SCIP_HYPERGRAPH_VERTEXDATA** pvertexdata  /**< Pointer for returning the vertex's data (may be \c NULL). */
   )
{
   assert(hypergraph);
   assert(pvertex);

   SCIP_CALL( ensureNumVertices(hypergraph, hypergraph->nvertices + 1) );

   *pvertex = hypergraph->nvertices;
   if( pvertexdata != NULL )
   {
      assert( hypergraph->nvertices >= 0 );
      *pvertexdata = (SCIP_HYPERGRAPH_VERTEXDATA*)(hypergraph->verticesdata
         + hypergraph->nvertices * hypergraph->sizevertexdata / sizeof(*(hypergraph->verticesdata)) ); /*lint --e{737}*/
   }

   hypergraph->nvertices++;
   hypergraph->hasvertexedges = FALSE;

   return SCIP_OKAY;
}

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
   )
{
   int first;
   int beyond;

   assert(hypergraph);
   assert(nvertices > 0);
   assert(vertices);
   assert(pedge);

   SCIP_CALL( ensureNumEdges(hypergraph, hypergraph->nedges + 1) );
   first = hypergraph->edgesverticesbeg[hypergraph->nedges];
   beyond = first + nvertices;
   SCIP_CALL( ensureNumEdgesVertices(hypergraph, beyond) );
   hypergraph->edgesverticesbeg[hypergraph->nedges + 1] = beyond;
   for( int i = 0; i < nvertices; ++i )
      hypergraph->edgesvertices[first + i] = vertices[i];
   SCIPsortInt(&hypergraph->edgesvertices[first], nvertices);

   *pedge = hypergraph->nedges;
   if( pedgedata != NULL )
   {
      assert( hypergraph->nedges >= 0 );
      *pedgedata = (SCIP_HYPERGRAPH_EDGEDATA*)(hypergraph->edgesdata +
         hypergraph->nedges * hypergraph->sizeedgedata / sizeof(*(hypergraph->edgesdata))); /*lint --e{737}*/
   }
   hypergraph->nedges++;
   hypergraph->hasvertexedges = FALSE;

   return SCIP_OKAY;
}

/** @brief get-key function for overlap hashtable */
static
SCIP_DECL_HASHGETKEY(overlapHashGetKey)
{  /*lint --e{715}*/
   return elem;
}

/** @brief equality function for overlap hashtable */
static
SCIP_DECL_HASHKEYEQ(overlapHashKeyEqPtr)
{
   size_t o1;
   size_t o2;
   SCIP_HYPERGRAPH* hypergraph;
   int begin1;
   int beyond1;
   int begin2;
   int beyond2;
   int i2;

   o1 = (size_t)key1 - 1;
   o2 = (size_t)key2 - 1;
   hypergraph = (SCIP_HYPERGRAPH*) userptr;
   begin1 = hypergraph->overlapsverticesbeg[o1];
   beyond1 = hypergraph->overlapsverticesbeg[o1 + 1];
   begin2 = hypergraph->overlapsverticesbeg[o2];
   beyond2 = hypergraph->overlapsverticesbeg[o2 + 1];

   if( beyond1 - begin1 != beyond2 - begin2 )
      return FALSE;

   i2 = begin2;
   for( int i1 = begin1; i1 < beyond1; ++i1 )
   {
      if( hypergraph->overlapsvertices[i1] != hypergraph->overlapsvertices[i2] )
         return FALSE;
      ++i2;
   }

   return TRUE;
}

/** @brief hash function for overlap hashtable */
static
SCIP_DECL_HASHKEYVAL(overlapHashKeyValPtr)
{
   size_t o;
   SCIP_HYPERGRAPH* hypergraph;
   int begin;
   int beyond;
   int i;
   uint32_t hash;

   o = (size_t)key - 1;
   hypergraph = (SCIP_HYPERGRAPH*) userptr;
   begin = hypergraph->overlapsverticesbeg[o];
   beyond = hypergraph->overlapsverticesbeg[o + 1];

   assert(beyond > begin);

   hash = (uint32_t) (beyond - begin);
   for( i = begin + 2; i < beyond; i += 3 )
   {
      hash = SCIPhashFour(hash, hypergraph->overlapsvertices[i - 2], hypergraph->overlapsvertices[i - 1],
         hypergraph->overlapsvertices[i]);
   }
   if( i - 1 < beyond )
      hash = SCIPhashThree(hash, hypergraph->overlapsvertices[i - 2], hypergraph->overlapsvertices[i - 1]);
   else if( i - 2 < beyond )
      hash = SCIPhashTwo(hash, hypergraph->overlapsvertices[i - 2]);

   return hash;
}

/** @brief finds an overlap specified by a sorted vertex set, potentially adding it */
static
SCIP_RETCODE findOverlap(
   SCIP_HYPERGRAPH*      hypergraph,         /**< Hypergraph. */
   int                   nvertices,          /**< Number of vertices. */
   int*                  vertices,           /**< Sorted array of vertices. */
   int*                  poverlap,           /**< Pointer for storing the overlap. */
   SCIP_Bool             add,                /**< Whether a missing overlap set should be added. */
   SCIP_Bool*            padded              /**< Whether a missing overlap set was added. */
   )
{
   int first;
   int beyond;
   void* element;
   size_t nextoverlap;

   assert(hypergraph);
   assert(nvertices > 0);
   assert(vertices);
   assert(poverlap);

   /* In order to find a set of vertices we need to add it to the overlaps' vertex storage. */

   SCIP_CALL( ensureNumOverlaps(hypergraph, hypergraph->noverlaps + 1) );
   first = hypergraph->overlapsverticesbeg[hypergraph->noverlaps];
   beyond = first + nvertices;
   SCIP_CALL( ensureNumOverlapsVertices(hypergraph, beyond) );
   for( int i = 0; i < nvertices; ++i )
      hypergraph->overlapsvertices[first + i] = vertices[i];
   hypergraph->overlapsverticesbeg[hypergraph->noverlaps + 1] = beyond;

   nextoverlap = (size_t) hypergraph->noverlaps + 1UL; /*lint !e571 */
   element = SCIPhashtableRetrieve(hypergraph->overlaphashtable, (void*) nextoverlap);
   if( element != NULL )
   {
      *poverlap = (size_t)element - 1; /*lint !e712*/
      if( padded )
         *padded = FALSE;
   }
   else if( add )
   {
      SCIP_CALL_ABORT( SCIPhashtableInsert(hypergraph->overlaphashtable, (void*) nextoverlap) );

      *poverlap = hypergraph->noverlaps;
      hypergraph->noverlaps++;
      if( padded )
         *padded = TRUE;
   }
   else
      *poverlap = -1;

   return SCIP_OKAY;
}

/**
 * @brief computes all overlaps and stores overlaps' vertices and all edges' overlaps
 *
 * Requires \ref SCIPhypergraphHasVertexEdges to be \c TRUE which results from \ref SCIPhypergraphComputeVerticesEdges.
 */
SCIP_RETCODE SCIPhypergraphComputeOverlaps(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_DECL_HYPERGRAPH_OVERLAP((*handler)), /**< Function to be called once the overlap is found. */
   void*                 userdata            /**< Pointer passed to \p handler. */
   )
{
   int memCommonVertices = 32;
   int* commonVertices = NULL;
   int numCommonVertices;
   int memEdgeOverlapPairs;
   int numEdgeOverlapPairs = 0;
   SCIP_Longint* edgeOverlapPairs = NULL;
   int o;
   int lastEdge;

   assert(hypergraph);

   if( !hypergraph->hasvertexedges )
      return SCIP_ERROR;

   if( hypergraph->hasoverlaps )
      return SCIP_OKAY;

   hypergraph->noverlaps = 0;
   hypergraph->overlapsverticesbeg[0] = 0;

   hypergraph->hasoverlaps = TRUE;
   assert(hypergraph->overlaphashtable == NULL);
   SCIP_CALL( SCIPhashtableCreate(&hypergraph->overlaphashtable, hypergraph->blkmem, 4 * hypergraph->nedges + 4,
      overlapHashGetKey, overlapHashKeyEqPtr, overlapHashKeyValPtr, hypergraph) );

   /* Allocate memory for the vertices in common to two edges. */
   SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &commonVertices, memCommonVertices) );

   /* Allocate memory for an array consisting of pairs (e,o) where o is an overlap incident to edge e. */
   memEdgeOverlapPairs = 2 * (hypergraph->nedges + hypergraph->noverlaps) + 2;
   SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &edgeOverlapPairs, memEdgeOverlapPairs) );

   /* We iterate over all edges e. */
   for( SCIP_HYPERGRAPH_EDGE e = 0; e < hypergraph->nedges; ++e )
   {
      int eBeyond;
      int eSize;

      eBeyond = hypergraph->edgesverticesbeg[e + 1];
      eSize = eBeyond - hypergraph->edgesverticesbeg[e];
      if( eSize > memCommonVertices )
      {
         int newcapacity;
         newcapacity = memCommonVertices;
         while( memCommonVertices < eSize )
            memCommonVertices *= 2;
         SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &commonVertices, newcapacity, memCommonVertices) );
      }

      /* We iterate over all vertices v of e and then over all edges f incident to v to avoid scanning all pairs
       * of edges because most of them will be disjoint. */
      for( int i = 0; i < eSize; ++i )
      {
         SCIP_HYPERGRAPH_VERTEX v;
         int vFirst, vBeyond;

         v = hypergraph->edgesvertices[hypergraph->edgesverticesbeg[e] + i];
         vFirst = hypergraph->verticesedgesbeg[v];
         vBeyond = hypergraph->verticesedgesbeg[v + 1];
         for( int j = vFirst; j < vBeyond; ++j )
         {
            SCIP_HYPERGRAPH_EDGE f;
            int fBeyond;
            int ie;
            int je;
            int u;
            int w;

            f = hypergraph->verticesedges[j];

            /* Avoid considering (e,f) and (f,e). */
            if( f <= e )
               continue;

            /* We now traverse the sorted lists of vertices of e and f simultaneously to find common vertices. */
            fBeyond = hypergraph->edgesverticesbeg[f + 1];
            numCommonVertices = 0;
            ie = hypergraph->edgesverticesbeg[e];
            je = hypergraph->edgesverticesbeg[f];
            while( ie < eBeyond && je < fBeyond )
            {
               u = hypergraph->edgesvertices[ie];
               w = hypergraph->edgesvertices[je];
               if( u < w )
                  ++ie;
               else if( u > w )
                  ++je;
               else
               {
                  /* Vertex u = w is part of e and of f.  */
                  if( numCommonVertices == 0 && u != v )
                  {
                     /* v is not the minimum vertex in e \cap f, so we skip f since we have considered it already
                      * for a smaller v. */
                     break;
                  }
                  commonVertices[numCommonVertices++] = u;
                  ++ie;
                  ++je;
               }
            }

            /* We only consider overlap sets that are not just vertices. */
            if( numCommonVertices >= 2 )
            {
               int overlap;
               SCIP_Bool added;
               uint64_t ebits;
               uint64_t fbits;
               uint64_t obits;
               uint64_t ocardinalitybits;

               /* Find out if that overlap (as a set) already exists. */
               SCIP_CALL( findOverlap(hypergraph, numCommonVertices, commonVertices, &overlap, TRUE, &added) );
               if( handler != NULL )
               {
                  SCIP_CALL( handler(hypergraph, overlap, SCIPhypergraphOverlapData(hypergraph, overlap), e, f, added,
                     userdata) );
               }

               if( numEdgeOverlapPairs >= memEdgeOverlapPairs - 1 )
               {
                  SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &edgeOverlapPairs, memEdgeOverlapPairs,
                     2 * memEdgeOverlapPairs) ); /*lint !e647 */
                  memEdgeOverlapPairs *= 2;
               }

               /* For the pair (overlap,e) we store (e << 40) | (|overlap| << 24) | overlap and sort by it to find
                * duplicates and the get overlaps incident to e in order of increasing cardinality. */
               assert(overlap >= 0);
               assert(numCommonVertices >= 0);
               assert(e >= 0);
               assert(f >= 0);
               obits = (uint64_t) overlap; /*lint !e571*/
               ocardinalitybits = ((uint64_t) numCommonVertices) << 24; /*lint !e571*/
               ebits = ((uint64_t) e) << 40; /*lint !e571*/
               ebits |= ocardinalitybits | obits;
               fbits = ((uint64_t) f) << 40; /*lint !e571*/
               fbits |= ocardinalitybits | obits;
               edgeOverlapPairs[numEdgeOverlapPairs] = (long long) ebits;
               edgeOverlapPairs[numEdgeOverlapPairs + 1] = (long long) fbits;
               numEdgeOverlapPairs += 2;
            }
         }
      }
   }

   SCIPdebugMessage("found %d edge-overlap incidences, including possible duplicates.\n", numEdgeOverlapPairs);
   SCIPsortLong(edgeOverlapPairs, numEdgeOverlapPairs);

   o = -1; /* Last pair written. */
   for( int i = 0; i < numEdgeOverlapPairs; ++i)
   {
      if( o < 0 || edgeOverlapPairs[i] != edgeOverlapPairs[o] )
         edgeOverlapPairs[++o] = edgeOverlapPairs[i];
   }

   numEdgeOverlapPairs = o + 1;
   SCIPdebugMessage("found %d unique edge-overlap incidences.\n", numEdgeOverlapPairs);

   /* Compute edges' incident overlaps. */
   if( hypergraph->memedgesoverlapsbeg == 0 || hypergraph->memedgesoverlapsbeg < hypergraph->nedges )
   {
      int newcapacity;

      newcapacity = MAX(256, hypergraph->memedgesoverlapsbeg);
      while( newcapacity < hypergraph->nedges )
         newcapacity *= 2;
      if( hypergraph->edgesoverlapsbeg )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesoverlapsbeg,
               hypergraph->memedgesoverlapsbeg + 1, newcapacity + 1) );
      }
      else
      {
         assert(hypergraph->memedgesoverlapsbeg == 0);
         SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesoverlapsbeg, newcapacity + 1) );
      }
      hypergraph->memedgesoverlapsbeg = newcapacity;
   }

   if( hypergraph->memedgesoverlaps < numEdgeOverlapPairs )
   {
      int newcapacity;

      newcapacity = MAX(256, hypergraph->memedgesoverlaps);
      while( newcapacity < numEdgeOverlapPairs )
         newcapacity *= 2;
      SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->edgesoverlaps,
            hypergraph->memedgesoverlaps, newcapacity) );
      hypergraph->memedgesoverlaps = newcapacity;
   }

   lastEdge = -1;
   for( int i = 0; i < numEdgeOverlapPairs; ++i )
   {
      uint64_t bits;
      int edge;
      int overlap;

      bits = (uint64_t) edgeOverlapPairs[i];
      edge = (int) (bits >> 40); /*lint !e704 */
      overlap = bits & ((1L << 24) - 1L); /*lint !e712 */

      while( lastEdge < edge )
         hypergraph->edgesoverlapsbeg[++lastEdge] = i;
      hypergraph->edgesoverlaps[i] = overlap;
   }
   while( lastEdge < hypergraph->nedges )
      hypergraph->edgesoverlapsbeg[++lastEdge] = numEdgeOverlapPairs;

   BMSfreeBlockMemoryArray(hypergraph->blkmem, &edgeOverlapPairs, memEdgeOverlapPairs);
   BMSfreeBlockMemoryArray(hypergraph->blkmem, &commonVertices, memCommonVertices);

   return SCIP_OKAY;
}

/** @brief computes each vertex' list of incident edges */
SCIP_RETCODE SCIPhypergraphComputeVerticesEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   int nincidences;

   assert(hypergraph);

   if( hypergraph->hasvertexedges )
      return SCIP_OKAY;

   if( hypergraph->memverticesedgesbeg == 0 || hypergraph->memverticesedgesbeg < hypergraph->nvertices )
   {
      int newcapacity = MAX(256, hypergraph->memvertices);
      while( newcapacity < hypergraph->nvertices )
         newcapacity *= 2;

      if( hypergraph->verticesedgesbeg )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesedgesbeg,
            hypergraph->memverticesedgesbeg + 1, newcapacity + 1) );
      }
      else
      {
         assert(hypergraph->memverticesedgesbeg == 0);
         SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesedgesbeg, newcapacity + 1) );
      }
      hypergraph->memverticesedgesbeg = newcapacity;
      SCIPdebugMessage("ensuring memory for %d vertices' edges slice.\n", newcapacity);
   }

   /* We know the total number of vertex-edge incidences from the edges. */
   nincidences = hypergraph->edgesverticesbeg[hypergraph->nedges];
   if( hypergraph->memverticesedges < nincidences )
   {
      int newcapacity;

      newcapacity = MAX(256, hypergraph->memverticesedges);
      while( newcapacity < nincidences )
         newcapacity *= 2;

      SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesedges,
            hypergraph->memverticesedges, newcapacity) );
      hypergraph->memverticesedges = newcapacity;
      SCIPdebugMessage("ensuring memory for %d vertices' edges.\n", newcapacity);
   }

   /* We initialize beg[v] = 0. */
   for( SCIP_HYPERGRAPH_VERTEX v = 0; v < hypergraph->nvertices; ++v)
      hypergraph->verticesedgesbeg[v] = 0;

   /* We traverse through all incidences (using the edges' vertex arrays) and count the degree deg(v) in beg[v]. */
   for( int i = 0; i < nincidences; ++i)
      hypergraph->verticesedgesbeg[hypergraph->edgesvertices[i]]++;

   /* This loop sets beg[0] = 0, beg[1] = deg(0), beg[2] = deg(0) + deg(1), beg[v] = deg(0) + ... + deg(v-1). */
   int current = 0;
   for( SCIP_HYPERGRAPH_VERTEX v = 0; v < hypergraph->nvertices; ++v )
   {
      int temp;

      temp = current;
      current += hypergraph->verticesedgesbeg[v];
      hypergraph->verticesedgesbeg[v] = temp;
   }

   /* We now traverse through all edges e and their incident vertices v, storing e at beg[v], which is incremented. */
   for( SCIP_HYPERGRAPH_EDGE e = 0; e < hypergraph->nedges; ++e )
   {
      int first;
      int beyond;

      first = hypergraph->edgesverticesbeg[e];
      beyond = hypergraph->edgesverticesbeg[e + 1];
      for( int i = first; i < beyond; ++i )
      {
         SCIP_HYPERGRAPH_VERTEX v = hypergraph->edgesvertices[i];
         hypergraph->verticesedges[hypergraph->verticesedgesbeg[v]++] = e;
      }
   }

   /* We now revert the incrementation of the beg[v]-values again. */
   for( SCIP_HYPERGRAPH_VERTEX v = hypergraph->nvertices; v > 0; --v )
      hypergraph->verticesedgesbeg[v] = hypergraph->verticesedgesbeg[v-1];
   hypergraph->verticesedgesbeg[0] = 0;
   hypergraph->hasvertexedges = TRUE;

   return SCIP_OKAY;
}

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
   )
{
   assert(hypergraph);
   assert(nvertices > 0);
   assert(vertices);
   assert(poverlap);
   assert(hypergraph->hasoverlaps);

   SCIP_CALL( findOverlap(hypergraph, nvertices, vertices, poverlap, FALSE, NULL) );

   return SCIP_OKAY;
}

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
   )
{
   int* commonVertices = NULL;
   int memCommonVertices;
   int firstBeyond;
   int secondBeyond;
   int numCommonVertices = 0;
   int i;
   int j;

   assert(hypergraph);
   assert(first >= 0 && first < hypergraph->nedges);
   assert(second >= 0 && second < hypergraph->nedges);
   assert(poverlap);
   assert(hypergraph->hasoverlaps);

   memCommonVertices = SCIPhypergraphEdgeSize(hypergraph, first);
   SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &commonVertices, memCommonVertices) );
   firstBeyond = hypergraph->edgesverticesbeg[first + 1];
   secondBeyond = hypergraph->edgesverticesbeg[second + 1];
   i = hypergraph->edgesverticesbeg[first];
   j = hypergraph->edgesverticesbeg[second];
   while( i < firstBeyond && j < secondBeyond )
   {
      int u;
      int w;

      u = hypergraph->edgesvertices[i];
      w = hypergraph->edgesvertices[j];
      if( u < w )
         ++i;
      else if( u > w )
         ++j;
      else
      {
         commonVertices[numCommonVertices++] = u;
         ++i;
         ++j;
      }
   }

   SCIP_CALL( findOverlap(hypergraph, numCommonVertices, commonVertices, poverlap, FALSE, NULL) );
   if( (*poverlap) < 0 && pvertex != NULL )
      *pvertex = numCommonVertices > 0 ? commonVertices[0] : -1;

   BMSfreeBlockMemoryArray(hypergraph->blkmem, &commonVertices, memCommonVertices);

   return SCIP_OKAY;
}

/**
 * @brief computes all overlaps' lists of incident edges
 *
 * Requires \ref SCIPhypergraphHasOverlaps to be \c TRUE which results from \ref SCIPhypergraphComputeOverlaps.
 */
SCIP_RETCODE SCIPhypergraphComputeOverlapsEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   int nincidences;
   int current;

   assert(hypergraph);
   assert(hypergraph->hasoverlaps);

   if( hypergraph->hasoverlapsedges )
      return SCIP_OKAY;

   if( hypergraph->memoverlapsedgesbeg < hypergraph->noverlaps + 1 )
   {
      int newcapacity;

      newcapacity = MAX(256, hypergraph->memoverlapsedgesbeg);
      while( newcapacity < hypergraph->noverlaps )
         newcapacity *= 2;

      if( hypergraph->overlapsedgesbeg )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsedgesbeg,
            hypergraph->memoverlapsedgesbeg + 1, newcapacity + 1) );
      }
      else
      {
         assert(hypergraph->memoverlapsedgesbeg == 0);
         SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsedgesbeg, newcapacity + 1) );
      }
      hypergraph->memoverlapsedgesbeg = newcapacity;
      SCIPdebugMessage("ensuring memory for %d overlaps' edges slice.\n", newcapacity);
   }

   /* We know the total number of overlap-edge incidences from the edges. */
   nincidences = hypergraph->edgesoverlapsbeg[hypergraph->nedges];
   if( hypergraph->memoverlapsedges < nincidences )
   {
      int newcapacity;

      newcapacity = MAX(256, hypergraph->memoverlaps);
      while( newcapacity < nincidences )
         newcapacity *= 2;

      SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->overlapsedges,
         hypergraph->memoverlapsedges, newcapacity) );
      hypergraph->memoverlapsedges = newcapacity;
      SCIPdebugMessage("ensuring memory for %d overlaps' edges.\n", newcapacity);
   }

   /* We initialize beg[v] = 0. */
   for( SCIP_HYPERGRAPH_OVERLAP o = 0; o < hypergraph->noverlaps; ++o )
      hypergraph->overlapsedgesbeg[o] = 0;

   /* We traverse through all incidences (using the edges' overlap arrays) and count the degree deg(o) in beg[o],
    * where deg(o) shall denote the number of edges incident to overlap o. */
   for( int i = 0; i < nincidences; ++i )
      hypergraph->overlapsedgesbeg[hypergraph->edgesoverlaps[i]]++;

   /* This loop sets beg[0] = 0, beg[1] = deg(0), beg[2] = deg(0) + deg(1), beg[o] = deg(0) + ... + deg(o-1) */
   current = 0;
   for( SCIP_HYPERGRAPH_OVERLAP o = 0; o < hypergraph->noverlaps; ++o )
   {
      int temp;

      temp = current;
      current += hypergraph->overlapsedgesbeg[o];
      hypergraph->overlapsedgesbeg[o] = temp;
   }

   /* We now traverse through all edges e and their incident overlaps o, storing e at beg[o], which is incremented. */
   for( SCIP_HYPERGRAPH_EDGE e = 0; e < hypergraph->nedges; ++e )
   {
      int first;
      int beyond;

      first = hypergraph->edgesoverlapsbeg[e];
      beyond = hypergraph->edgesoverlapsbeg[e + 1];
      for( int i = first; i < beyond; ++i )
      {
         SCIP_HYPERGRAPH_OVERLAP o = hypergraph->edgesoverlaps[i];
         hypergraph->overlapsedges[hypergraph->overlapsedgesbeg[o]++] = e;
      }
   }

   /* We now revert the incrementation of the beg[o]-values again. */
   for( SCIP_HYPERGRAPH_OVERLAP o = hypergraph->noverlaps; o > 0; --o )
      hypergraph->overlapsedgesbeg[o] = hypergraph->overlapsedgesbeg[o-1];
   hypergraph->overlapsedgesbeg[0] = 0;
   hypergraph->hasoverlapsedges = TRUE;

   return SCIP_OKAY;
}


/**
 * @brief computes all vertices' lists of incident overlaps
 *
 * Requires \ref SCIPhypergraphHasOverlaps.
 */
SCIP_RETCODE SCIPhypergraphComputeVerticesOverlaps(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   int nincidences;
   int current;

   assert(hypergraph);
   assert(hypergraph->hasoverlaps);

   if( hypergraph->hasverticesoverlaps )
      return SCIP_OKAY;

   if( hypergraph->memverticesoverlapsbeg < hypergraph->nvertices )
   {
      int newcapacity;

      newcapacity = MAX(256, hypergraph->memverticesoverlapsbeg);
      while( newcapacity < hypergraph->nvertices )
         newcapacity *= 2;

      if( hypergraph->verticesoverlapsbeg )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesoverlapsbeg,
               hypergraph->memverticesoverlapsbeg, newcapacity) );
      }
      else
      {
         assert(hypergraph->memverticesoverlapsbeg == 0);
         SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesoverlapsbeg, newcapacity) );
      }
      hypergraph->memverticesoverlapsbeg = newcapacity;
      SCIPdebugMessage("ensuring memory for %d vertices' overlaps slice.\n", newcapacity);
   }

   /* We know the total number of vertex-overlap incidences from the over. */
   nincidences = hypergraph->overlapsverticesbeg[hypergraph->noverlaps];
   if( hypergraph->memverticesoverlaps < nincidences )
   {
      int newcapacity;

      newcapacity = MAX(256, hypergraph->memverticesoverlaps);
      while( newcapacity < nincidences )
         newcapacity *= 2;

      SCIP_ALLOC( BMSreallocBlockMemoryArray(hypergraph->blkmem, &hypergraph->verticesoverlaps,
            hypergraph->memverticesoverlaps, newcapacity) );
      hypergraph->memverticesoverlaps = newcapacity;
      SCIPdebugMessage("ensuring memory for %d vertices' overlaps.\n", newcapacity);
   }

   /* We initialize beg[v] = 0. */
   for( SCIP_HYPERGRAPH_VERTEX v = 0; v < hypergraph->nvertices; ++v )
      hypergraph->verticesoverlapsbeg[v] = 0;

   /* We traverse through all incidences (using the overlaps' vertex arrays) and count the degree deg(v) in beg[v],
    * where deg(v) shall denote the number of overlaps incident to vertex v. */
   for( int i = 0; i < nincidences; ++i )
      hypergraph->verticesoverlapsbeg[hypergraph->overlapsvertices[i]]++;

   /* This loop sets beg[0] = 0, beg[1] = deg(0), beg[2] = deg(0) + deg(1), beg[v] = deg(0) + ... + deg(v-1) */
   current = 0;
   for( SCIP_HYPERGRAPH_VERTEX v = 0; v < hypergraph->nvertices; ++v )
   {
      int temp;

      temp = current;
      current += hypergraph->verticesoverlapsbeg[v];
      hypergraph->verticesoverlapsbeg[v] = temp;
   }

   /* We traverse through all overlaps o and their incident vertices v, storing o at beg[v], which is incremented. */
   for( SCIP_HYPERGRAPH_OVERLAP o = 0; o < hypergraph->noverlaps; ++o )
   {
      int first;
      int beyond;

      first = hypergraph->overlapsverticesbeg[o];
      beyond = hypergraph->overlapsverticesbeg[o + 1];
      for( int i = first; i < beyond; ++i )
      {
         SCIP_HYPERGRAPH_VERTEX v = hypergraph->overlapsvertices[i];
         hypergraph->verticesoverlaps[hypergraph->verticesoverlapsbeg[v]++] = o;
      }
   }

   /* We now revert the incrementation of the beg[o]-values again. */
   for( SCIP_HYPERGRAPH_VERTEX v = hypergraph->nvertices; v > 0; --v )
      hypergraph->verticesoverlapsbeg[v] = hypergraph->verticesoverlapsbeg[v-1];
   hypergraph->verticesoverlapsbeg[0] = 0;
   hypergraph->hasverticesoverlaps = TRUE;

   return SCIP_OKAY;
}

/** @brief returns whether overlaps \p overlap1 and \p overlap2 are disjoint */
SCIP_Bool SCIPhypergraphOverlapsDisjoint(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap1,         /**< First overlap. */
   SCIP_HYPERGRAPH_OVERLAP overlap2          /**< Second overlap. */
   )
{
   int size1;
   SCIP_HYPERGRAPH_VERTEX* vertices1;
   int size2;
   SCIP_HYPERGRAPH_VERTEX* vertices2;
   int i = 0;
   int j = 0;

   assert(hypergraph);

   /* We loop over the vertices in parallel because they are sorted. */
   size1 = SCIPhypergraphOverlapSize(hypergraph, overlap1);
   vertices1 = SCIPhypergraphOverlapVertices(hypergraph, overlap1);
   size2 = SCIPhypergraphOverlapSize(hypergraph, overlap2);
   vertices2 = SCIPhypergraphOverlapVertices(hypergraph, overlap2);
   while( i < size1 && j < size2 )
   {
      if( vertices1[i] < vertices2[j] )
         ++i;
      else if( vertices1[i] > vertices2[j] )
         ++j;
      else
         return FALSE;
   }

   return TRUE;
}

/** @brief asserts that the hypergraph data structures are valid */
SCIP_Bool SCIPhypergraphIsValid(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   FILE*                 file                /**< A file stream to print problems to (can be \c NULL). */
   )
{
   SCIP_Bool valid = TRUE;

   assert(hypergraph);

   if( SCIPhypergraphHasVertexEdges(hypergraph) )
   {
      int nincidences;
      long long* incidences1 = NULL;
      long long* incidences2 = NULL;
      int i1 = 0;
      int i2 = 0;

      nincidences = SCIPhypergraphGetNIncidences(hypergraph);
      SCIP_ALLOC_ABORT( BMSallocBlockMemoryArray(hypergraph->blkmem, &incidences1, nincidences) );
      SCIP_ALLOC_ABORT( BMSallocBlockMemoryArray(hypergraph->blkmem, &incidences2, nincidences) );
      for( SCIP_HYPERGRAPH_EDGE e = 0; e < SCIPhypergraphGetNEdges(hypergraph); ++e )
      {
         int esize = SCIPhypergraphEdgeSize(hypergraph, e);
         for( int i = 0; i < esize; ++i )
         {
            SCIP_HYPERGRAPH_VERTEX v = SCIPhypergraphEdgeVertices(hypergraph, e)[i];
            long long pair = v + (long long) SCIPhypergraphGetNVertices(hypergraph) * (long long) e;
            incidences1[i1++] = pair;
         }
      }

      if ( i1 != nincidences )
      {
         valid = FALSE;
         if( file )
         {
            fprintf(file, "SCIPhypergraphIsValid detected inconsistency: "
               "number of incidences is claimed to be %d but counting via edges yields %d!\n", nincidences, i1);
            fflush(file);
         }
      }

      for( SCIP_HYPERGRAPH_VERTEX v = 0; v < SCIPhypergraphGetNVertices(hypergraph); ++v )
      {
         int first = SCIPhypergraphVertexEdgesFirst(hypergraph, v);
         int beyond = SCIPhypergraphVertexEdgesBeyond(hypergraph, v);
         for( int j = first; j < beyond; ++j )
         {
            SCIP_HYPERGRAPH_EDGE e = SCIPhypergraphVertexEdgesGetAtIndex(hypergraph, j);
            long long pair = v + (long long) SCIPhypergraphGetNVertices(hypergraph) * (long long) e;
            incidences2[i2++] = pair;
         }
      }

      if ( i2 != nincidences )
      {
         valid = FALSE;
         if( file )
         {
            fprintf(file, "SCIPhypergraphIsValid detected inconsistency: "
               "number of incidences is claimed to be %d but counting via vertices yields %d!\n", nincidences, i2);
            fflush(file);
         }
      }

      if( valid )
      {
         SCIPsortLong(incidences1, nincidences);
         SCIPsortLong(incidences2, nincidences);

         for( int i = 0; i < nincidences; ++i )
         {
            if( incidences1[i] != incidences2[i] )
            {
               valid = FALSE;
               if( file )
               {
                  fprintf(file, "SCIPhypergraphIsValid detected inconsistency: "
                     "incidence #%d is %lld (via edges) and %lld (via vertices)!\n", i, incidences1[i], incidences2[i]);
                  fflush(file);
               }
            }
         }
      }

      BMSfreeBlockMemoryArray(hypergraph->blkmem, &incidences2, nincidences);
      BMSfreeBlockMemoryArray(hypergraph->blkmem, &incidences1, nincidences);
   }

   if( SCIPhypergraphHasVertexEdges(hypergraph) )
   {
      int nincidences;
      long long* incidences1 = NULL;
      long long* incidences2 = NULL;
      int i1 = 0;
      int i2 = 0;

      nincidences = SCIPhypergraphGetNIncidences(hypergraph);
      SCIP_ALLOC_ABORT( BMSallocBlockMemoryArray(hypergraph->blkmem, &incidences1, nincidences) );
      SCIP_ALLOC_ABORT( BMSallocBlockMemoryArray(hypergraph->blkmem, &incidences2, nincidences) );
      for( SCIP_HYPERGRAPH_EDGE e = 0; e < SCIPhypergraphGetNEdges(hypergraph); ++e )
      {
         int esize = SCIPhypergraphEdgeSize(hypergraph, e);
         for( int i = 0; i < esize; ++i )
         {
            SCIP_HYPERGRAPH_VERTEX v = SCIPhypergraphEdgeVertices(hypergraph, e)[i];
            long long pair = v + (long long) SCIPhypergraphGetNVertices(hypergraph) * (long long) e;
            incidences1[i1++] = pair;
         }
      }

      if ( i1 != nincidences )
      {
         valid = FALSE;
         if( file )
         {
            fprintf(file, "SCIPhypergraphIsValid detected inconsistency: "
               "number of incidences is claimed to be %d but counting via edges yields %d!\n", nincidences, i1);
            fflush(file);
         }
      }

      for( SCIP_HYPERGRAPH_VERTEX v = 0; v < SCIPhypergraphGetNVertices(hypergraph); ++v )
      {
         int first = SCIPhypergraphVertexEdgesFirst(hypergraph, v);
         int beyond = SCIPhypergraphVertexEdgesBeyond(hypergraph, v);
         for( int j = first; j < beyond; ++j )
         {
            SCIP_HYPERGRAPH_EDGE e = SCIPhypergraphVertexEdgesGetAtIndex(hypergraph, j);
            long long pair = v + (long long) SCIPhypergraphGetNVertices(hypergraph) * (long long) e;
            incidences2[i2++] = pair;
         }
      }

      if ( i2 != nincidences )
      {
         valid = FALSE;
         if( file )
         {
            fprintf(file, "SCIPhypergraphIsValid detected inconsistency: "
               "number of incidences is claimed to be %d but counting via vertices yields %d!\n", nincidences, i2);
            fflush(file);
         }
      }

      if( valid )
      {
         SCIPsortLong(incidences1, nincidences);
         SCIPsortLong(incidences2, nincidences);

         for( int i = 0; i < nincidences; ++i )
         {
            if( incidences1[i] != incidences2[i] )
            {
               valid = FALSE;
               if( file )
               {
                  fprintf(file, "SCIPhypergraphIsValid detected inconsistency: "
                     "incidence #%d is %lld (via edges) and %lld (via vertices)!\n", i, incidences1[i], incidences2[i]);
                  fflush(file);
               }
            }
         }
      }

      BMSfreeBlockMemoryArray(hypergraph->blkmem, &incidences2, nincidences);
      BMSfreeBlockMemoryArray(hypergraph->blkmem, &incidences1, nincidences);
   }

   if( SCIPhypergraphHasOverlaps(hypergraph) )
   {
      int n;
      SCIP_Bool* marked = NULL;

      /* Verify edge-overlap incidences. */
      n = SCIPhypergraphGetNVertices(hypergraph);
      SCIP_ALLOC_ABORT( BMSallocBlockMemoryArray(hypergraph->blkmem, &marked, n) );
      for( SCIP_HYPERGRAPH_VERTEX v = 0; v < n; ++v )
         marked[v] = FALSE;
      for( SCIP_HYPERGRAPH_EDGE e = 0; e < hypergraph->nedges && valid; ++e )
      {
         int first;
         int beyond;
         int esize;
         int lastcardinality = 0;

         esize = SCIPhypergraphEdgeSize(hypergraph, e);
         SCIP_HYPERGRAPH_VERTEX* evertices = SCIPhypergraphEdgeVertices(hypergraph, e);
         for( int i = 0; i < esize; ++i )
            marked[evertices[i]] = TRUE;
         first = SCIPhypergraphEdgesOverlapsFirst(hypergraph, e);
         beyond = SCIPhypergraphEdgesOverlapsBeyond(hypergraph, e);
         for( int i = first; i < beyond && valid; ++i )
         {
            SCIP_HYPERGRAPH_OVERLAP o;
            int cardinality;

            o = SCIPhypergraphEdgesOverlapsGetAtIndex(hypergraph, i);
            cardinality = SCIPhypergraphOverlapSize(hypergraph, o);

            /* Check if the cardinalities are nontrivial and form a non-decreasing sequences for each edge. */
            if (cardinality < 1 || cardinality < lastcardinality)
            {
               valid = FALSE;
               if( file )
               {
                  fprintf(file, "SCIPhypergraphIsValid detected inconsistency: edge e%d = {", e);
                  for( int k = 0; k < esize; ++k )
                     fprintf(file, "%sv%d", k ? "," : "", SCIPhypergraphEdgeVertices(hypergraph, e)[k]);
                  fprintf(file, "} has incidence #%d of [%d,%d) to overlap o%d = {", i, first, beyond, o);
                  for( int k = 0; k < SCIPhypergraphOverlapSize(hypergraph, o); ++k )
                     fprintf(file, "%sv%d", k ? "," : "", SCIPhypergraphOverlapVertices(hypergraph, o)[k]);
                  fprintf(file, "} of cardinality |o%d| = %d, but previous cardinality is %d.\n", o, cardinality,
                     lastcardinality);
                  fflush(file);
               }
            }
            lastcardinality = cardinality;

            for( int j = 0; j < SCIPhypergraphOverlapSize(hypergraph, o) && valid; ++j )
            {
               if( !marked[SCIPhypergraphOverlapVertices(hypergraph, o)[j]])
               {
                  valid = FALSE;
                  if( file )
                  {
                     fprintf(file, "SCIPhypergraphIsValid detected inconsistency: edge e%d = {", e);
                     for( int k = 0; k < esize; ++k )
                        fprintf(file, "%sv%d", k ? "," : "", SCIPhypergraphEdgeVertices(hypergraph, e)[k]);
                     fprintf(file, "} has incidence #%d of [%d,%d) to overlap o%d = {", i, first, beyond, o);
                     for( int k = 0; k < SCIPhypergraphOverlapSize(hypergraph, o); ++k )
                        fprintf(file, "%sv%d", k ? "," : "", SCIPhypergraphOverlapVertices(hypergraph, o)[k]);
                     fprintf(file, "} which is not a subset!\n");
                     fflush(file);
                  }
               }
            }
         }
         for( int i = 0; i < esize; ++i )
            marked[evertices[i]] = FALSE;
      }
      BMSfreeBlockMemoryArray(hypergraph->blkmem, &marked, n);
   }

   return valid;
}

/** @brief initializes a hypergraph iterator's internal memory */
SCIP_RETCODE SCIPhypergraphIterInit(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph of this iterator. */
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator to be initialized. */
   )
{
   assert(hypergraph);
   assert(iterator);

   iterator->base = -1;
   iterator->vertexidx = -1;
   iterator->minvertex = -1;
   iterator->edgeidx = -1;
   iterator->adjacent = -1;
   iterator->sizecommonvertices = 32;
   iterator->commonvertices = NULL;
   iterator->minoverlapsize = 2;
   iterator->onlylater = FALSE;
   iterator->findoverlaps = FALSE;

   SCIP_ALLOC( BMSallocBlockMemoryArray(hypergraph->blkmem, &iterator->commonvertices, 32) );

   return SCIP_OKAY;
}

/** @brief frees a hypergraph iterator's internal memory */
void SCIPhypergraphIterClear(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph of this iterator. */
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator to be cleared. */
   )
{
   assert(hypergraph);
   assert(iterator);

   BMSfreeBlockMemoryArray(hypergraph->blkmem, &iterator->commonvertices, iterator->sizecommonvertices);
}

/** @brief initializes the iterator to the first adjacent edge of \p base */
void SCIPhypergraphIterStart(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph of this iterator. */
   SCIP_HYPERGRAPH_ITER* iterator,           /**< Iterator to be initialized. */
   SCIP_HYPERGRAPH_EDGE  base,               /**< Base edge. */
   unsigned int          minoverlapsize,     /**< Minimum size of the overlap. */
   SCIP_Bool             onlylater,          /**< Whether to only consider edges greater than \p base. */
   SCIP_Bool             findoverlaps        /**< Whether to compute the overlap sets. */
   )
{
   assert(hypergraph);
   assert(iterator);
   assert(base >= 0);
   assert(base < hypergraph->nedges);
   assert(hypergraph->hasvertexedges);
   assert(hypergraph->hasoverlaps || !findoverlaps);

   iterator->base = base;
   iterator->vertexidx = -1;
   iterator->edgeidx = -1;
   iterator->minvertex = -1;
   iterator->adjacent = -1;
   iterator->ncommonvertices = 0;
   iterator->overlap = -1;
   iterator->minoverlapsize = minoverlapsize;
   iterator->onlylater = onlylater;
   iterator->findoverlaps = findoverlaps;

   SCIPhypergraphIterNext(hypergraph, iterator);
}

/** @brief returns whether the iterator is valid */
SCIP_Bool SCIPhypergraphIterValid(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   )
{
   assert(iterator);

   return iterator->vertexidx != -2;
}

/** @brief initializes the iterator to the first adjacent edge of \p base */
void SCIPhypergraphIterNext(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph of this iterator. */
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   )
{
   assert(hypergraph);
   assert(iterator);
   assert(iterator->base >= 0);
   assert(iterator->base < hypergraph->nedges);
   assert(iterator->vertexidx == -1 || iterator->vertexidx >= hypergraph->edgesverticesbeg[iterator->base]);
   assert(iterator->vertexidx < hypergraph->edgesverticesbeg[iterator->base + 1]);

   if( iterator->vertexidx < 0 )
   {
      /* Freshly initialized iterator. */
      iterator->vertexidx = hypergraph->edgesverticesbeg[iterator->base];
      if( iterator->vertexidx == hypergraph->edgesverticesbeg[iterator->base + 1] )
      {
         /* No vertices in this edge. */
         iterator->vertexidx = -2;
         return;
      }

      /* Edge index is set to one before the first edge and incremented in next loop. */
      iterator->minvertex = hypergraph->edgesvertices[iterator->vertexidx];
      iterator->edgeidx = hypergraph->verticesedgesbeg[iterator->minvertex] - 1;
   }

   while( iterator->vertexidx >= 0 )
   {
      int i, iend;
      int j, jend;

      assert(iterator->minvertex >= 0);
      assert(iterator->minvertex < hypergraph->nvertices);
      assert(iterator->edgeidx >= hypergraph->verticesedgesbeg[iterator->minvertex] - 1);
      assert(iterator->edgeidx < hypergraph->verticesedgesbeg[iterator->minvertex + 1]);

      /* Go to next edge of vertex. */

      iterator->edgeidx++;
      if( iterator->edgeidx < hypergraph->verticesedgesbeg[iterator->minvertex + 1] )
      {
         /* There is a next edge. */
         iterator->adjacent = hypergraph->verticesedges[iterator->edgeidx];
      }
      else
      {
         /* There is no next edge, so we consider the next vertex. */
         iterator->vertexidx++;
         if( iterator->vertexidx == hypergraph->edgesverticesbeg[iterator->base + 1] )
         {
            /* There is no next vertex either. */
            iterator->vertexidx = -2;
            return;
         }

         /* Edge index is set to one before the first edge of the next vertex and incremented in the next iteration. */
         iterator->minvertex = hypergraph->edgesvertices[iterator->vertexidx];
         iterator->edgeidx = hypergraph->verticesedgesbeg[iterator->minvertex] - 1;
         continue;
      }

      /* Skip if we consider the base edge. */
      if( iterator->adjacent == iterator->base )
         continue;

      /* Check for adjacent < base and discard if not desired. */
      if( iterator->onlylater && (iterator->adjacent < iterator->base) )
         continue;

      /* We traverse the common vertices. */
      iterator->ncommonvertices = 0;
      i = hypergraph->edgesverticesbeg[iterator->base];
      iend = hypergraph->edgesverticesbeg[iterator->base + 1];
      j = hypergraph->edgesverticesbeg[iterator->adjacent];
      jend = hypergraph->edgesverticesbeg[iterator->adjacent + 1];
      while( i < iend && j < jend )
      {
         SCIP_HYPERGRAPH_VERTEX ivertex, jvertex;

         ivertex = hypergraph->edgesvertices[i];
         jvertex = hypergraph->edgesvertices[j];

         if( ivertex < jvertex )
            ++i;
         else if( ivertex > jvertex )
            ++j;
         else
         {
            if( iterator->ncommonvertices == 0 )
            {
               /* First common vertex must be the iterator vertex. Otherwise we abort. */
               if( ivertex != iterator->minvertex )
               {
                  iterator->ncommonvertices = -1;
                  break;
               }
            }
            if( iterator->ncommonvertices == iterator->sizecommonvertices )
            {
               int newsize = 2 * iterator->sizecommonvertices;
               SCIP_ALLOC_ABORT( BMSreallocBlockMemoryArray(hypergraph->blkmem, &iterator->commonvertices,
                  iterator->sizecommonvertices, newsize) );
               iterator->sizecommonvertices = newsize;
            }
            iterator->commonvertices[iterator->ncommonvertices] = ivertex;
            iterator->ncommonvertices++;
            ++i;
            ++j;
         }
      }
      assert(iterator->ncommonvertices != 0);

      if( iterator->ncommonvertices < iterator->minoverlapsize )
      {
         /* Abort since either ncommonvertices = -1 (minvertex not minimum vertex in intersection)
          * or if it was smaller than requested. */
         continue;
      }

      if( iterator->findoverlaps )
      {
         /* If requested, find corresponding overlap set. */
         if( iterator->ncommonvertices >= 2 )
         {
            SCIP_CALL_ABORT( findOverlap(hypergraph, iterator->ncommonvertices, iterator->commonvertices,
               &iterator->overlap, FALSE, NULL) );
         }
         else
            iterator->overlap = -1;
      }

      break;
   }
}

/** @brief returns the base edge */
SCIP_HYPERGRAPH_EDGE SCIPhypergraphIterBase(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   )
{
   return iterator->base;
}

/** @brief returns the current adjacent edge */
SCIP_HYPERGRAPH_EDGE SCIPhypergraphIterAdjacent(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   )
{
   return iterator->adjacent;
}

/** @brief returns the minimum vertex in the intersection of the base and the current adjacent edge */
SCIP_HYPERGRAPH_VERTEX SCIPhypergraphIterMinVertex(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   )
{
   return iterator->minvertex;
}

/**
 * @brief returns the overlap for the intersection of the base and the current adjacent edge
 *
 * If the intersection of the two edges has only one element, then -1 is returned, and if overlap information was not
 * requested in \ref SCIPhypergraphIterStart, then -2 is returned.
 */
SCIP_HYPERGRAPH_OVERLAP SCIPhypergraphIterOverlap(
   SCIP_HYPERGRAPH_ITER* iterator            /**< Iterator. */
   )
{
   return iterator->overlap;
}


/*
 * simple functions implemented as defines
 */

#ifndef NDEBUG

/*
 * In debug mode, the following methods are implemented as function calls to ensure type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPhypergraphNumVertices
#undef SCIPhypergraphNumEdges
#undef SCIPhypergraphBlkmem
#undef SCIPhypergraphNumOverlaps
#undef SCIPhypergraphVertexData
#undef SCIPhypergraphEdgeData
#undef SCIPhypergraphOverlapData
#undef SCIPhypergraphEdgeSize
#undef SCIPhypergraphEdgeVertices
#undef SCIPhypergraphHasVertexEdges
#undef SCIPhypergraphVerticesEdgesFirst
#undef SCIPhypergraphVerticesEdgesBeyond
#undef SCIPhypergraphVerticesEdgesGet
#undef SCIPhypergraphHasOverlaps
#undef SCIPhypergraphOverlapSize
#undef SCIPhypergraphOverlapVertices
#undef SCIPhypergraphEdgesOverlapsFirst
#undef SCIPhypergraphEdgesOverlapsBeyond
#undef SCIPhypergraphEdgesOverlapsGet
#undef SCIPhypergraphHasOverlapsEdges
#undef SCIPhypergraphOverlapsEdgesFirst
#undef SCIPhypergraphOverlapsEdgesBeyond
#undef SCIPhypergraphOverlapsEdgesGet
#undef SCIPhypergraphHasVerticesOverlaps
#undef SCIPhypergraphVerticesOverlapsFirst
#undef SCIPhypergraphVerticesOverlapsBeyond
#undef SCIPhypergraphVerticesOverlapsGet


/** @brief returns the number of vertices */
int SCIPhypergraphGetNVertices(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph);

   return hypergraph->nvertices;
}

/** @brief returns the number of edges */
int SCIPhypergraphGetNEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph);

   return hypergraph->nedges;
}

/** @brief returns the hypergraph's block memory structure */
BMS_BLKMEM* SCIPhypergraphBlkmem(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph != NULL);

   return hypergraph->blkmem;
}

/** @brief returns the number of overlaps */
int SCIPhypergraphGetNOverlaps(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph);

   return hypergraph->noverlaps;
}

/** @brief returns the number of edge-vertex incidences */
int SCIPhypergraphGetNIncidences(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph);

   return hypergraph->edgesverticesbeg[hypergraph->nedges];
}

/** @brief returns additional data of \p vertex */
SCIP_HYPERGRAPH_VERTEXDATA* SCIPhypergraphVertexData(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex             /**< A vertex. */
   )
{
   assert(hypergraph != NULL);

   return (SCIP_HYPERGRAPH_VERTEXDATA*)(hypergraph->verticesdata
      + vertex * hypergraph->sizevertexdata / sizeof(*(hypergraph->verticesdata))); /*lint --e{737}*/
}

/** @brief returns additional data of \p edge */
SCIP_HYPERGRAPH_EDGEDATA* SCIPhypergraphEdgeData(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< An edge. */
   )
{
   assert(hypergraph != NULL);

   return (SCIP_HYPERGRAPH_EDGEDATA*)(hypergraph->edgesdata
      + edge * hypergraph->sizeedgedata / sizeof(*(hypergraph->edgesdata))); /*lint --e{737}*/
}

/** @brief returns additional data of \p overlap */
SCIP_HYPERGRAPH_OVERLAPDATA* SCIPhypergraphOverlapData(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< An overlap. */
   )
{
   assert(hypergraph != NULL);

   return (SCIP_HYPERGRAPH_OVERLAPDATA*)(hypergraph->overlapsdata
      + overlap * hypergraph->sizeoverlapdata / sizeof(*(hypergraph->overlapsdata))); /*lint --e{737}*/
}

/** @brief returns the number of vertices of \p edge */
int SCIPhypergraphEdgeSize(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< An edge. */
   )
{
   assert(hypergraph);
   assert(edge >= 0 && edge < hypergraph->nedges);

   return hypergraph->edgesverticesbeg[edge + 1] - hypergraph->edgesverticesbeg[edge];
}

/**
 * @brief returns the array of vertices of \p edge
 *
 * The length of the array is \ref SCIPhypergraphEdgeSize.
 */

SCIP_HYPERGRAPH_VERTEX* SCIPhypergraphEdgeVertices(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< An edge. */
   )
{
   assert(hypergraph);
   assert(edge >= 0 && edge < hypergraph->nedges);

   return &hypergraph->edgesvertices[hypergraph->edgesverticesbeg[edge]];
}

/**
 * @brief returns whether vertices' incident edges are known.
 *
 * Use \ref SCIPhypergraphComputeVerticesEdges to compute them.
 */
SCIP_Bool SCIPhypergraphHasVertexEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph);

   return hypergraph->hasvertexedges;
}

/** @brief returns an index for the first edge incident to \p vertex */
int SCIPhypergraphVertexEdgesFirst(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex             /**< A vertex. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasvertexedges);

   return hypergraph->verticesedgesbeg[vertex];
}

/** @brief returns an index beyond the last edge incident to \p vertex */
int SCIPhypergraphVertexEdgesBeyond(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex             /**< A vertex. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasvertexedges);

   return hypergraph->verticesedgesbeg[vertex + 1];
}

/**
 * @brief returns the edge corresponding to \p index that is incident to a vertex
 *
 * See \ref SCIPhypergraphVertexEdgesFirst and \ref SCIPhypergraphVertexEdgesBeyond to obtain such indices for a vertex.
 */
SCIP_HYPERGRAPH_EDGE SCIPhypergraphVertexEdgesGetAtIndex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   idx                 /**< Index. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasvertexedges);

   return hypergraph->verticesedges[idx];
}

/**
 * @brief returns whether edges' overlaps and overlaps' vertices are known.
 *
 * Use \ref SCIPhypergraphComputeOverlaps to compute them.
 */
SCIP_Bool SCIPhypergraphHasOverlaps(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph);

   return hypergraph->hasoverlaps;
}

/**
 * @brief returns the number of vertices of \p overlap
 *
 * Requires \ref SCIPhypergraphHasOverlaps to be \c TRUE which results from \ref SCIPhypergraphComputeOverlaps.
 */
int SCIPhypergraphOverlapSize(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< Overlap. */
   )
{
   assert(hypergraph);
   assert(overlap >= 0 && overlap < hypergraph->noverlaps);
   assert(hypergraph->hasoverlaps);

   return hypergraph->overlapsverticesbeg[overlap + 1] - hypergraph->overlapsverticesbeg[overlap];
}

/**
 * @brief returns the array of sorted vertices of \p overlap
 *
 * The length of the array is \ref SCIPhypergraphOverlapSize.
 * Requires \ref SCIPhypergraphHasOverlaps to be \c TRUE which results from \ref SCIPhypergraphComputeOverlaps.
 */
SCIP_HYPERGRAPH_VERTEX* SCIPhypergraphOverlapVertices(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< Overlap. */
   )
{
   assert(hypergraph);
   assert(overlap >= 0 && overlap < hypergraph->noverlaps);
   assert(hypergraph->hasoverlaps);

   return &hypergraph->overlapsvertices[hypergraph->overlapsverticesbeg[overlap]];
}

/** @brief returns an index for the first overlap incident to \p edge */
int SCIPhypergraphEdgesOverlapsFirst(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< Edge. */
   )
{
   assert(hypergraph);
   assert(edge >= 0 && edge < hypergraph->nedges);
   assert(hypergraph->hasoverlaps);

   return hypergraph->edgesoverlapsbeg[edge];
}

/** @brief returns an index beyond the last overlap incident to \p edge */
int SCIPhypergraphEdgesOverlapsBeyond(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_EDGE  edge                /**< Edge. */
   )
{
   assert(hypergraph);
   assert(edge >= 0 && edge < hypergraph->nedges);
   assert(hypergraph->hasoverlaps);

   return hypergraph->edgesoverlapsbeg[edge + 1];
}

/**
 * @brief returns the overlap corresponding to \p idx that is incident to an edge
 *
 * See \ref SCIPhypergraphEdgesOverlapsFirst and \ref SCIPhypergraphEdgesOverlapsBeyond to obtain such indices for an
 * edge.
 */
SCIP_HYPERGRAPH_OVERLAP SCIPhypergraphEdgesOverlapsGetAtIndex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   idx                 /**< Index. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasoverlaps);

   return hypergraph->edgesoverlaps[idx];
}

/**
 * @brief returns whether overlaps' incident edges are known.
 *
 * Use \ref SCIPhypergraphComputeOverlapsEdges to compute them.
 */
SCIP_Bool SCIPhypergraphHasOverlapsEdges(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph);

   return hypergraph->hasoverlapsedges;
}

/** @brief returns an index for the first edge incident to \p overlap */
int SCIPhypergraphOverlapsEdgesFirst(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< Overlap. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasoverlapsedges);

   return hypergraph->overlapsedgesbeg[overlap];
}

/** @brief returns an index beyond the last edge incident to \p overlap */
int SCIPhypergraphOverlapsEdgesBeyond(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_OVERLAP overlap           /**< Overlap. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasoverlapsedges);

   return hypergraph->overlapsedgesbeg[overlap + 1];
}

/**
 * @brief returns the edge corresponding to \p idx that is incident to an overlap
 *
 * See \ref SCIPhypergraphOverlapsEdgesFirst and \ref SCIPhypergraphOverlapsEdgesBeyond to obtain such indices for an
 * overlap.
 */
SCIP_HYPERGRAPH_OVERLAP SCIPhypergraphOverlapsEdgesGetAtIndex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   idx                 /**< Index. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasoverlapsedges);

   return hypergraph->overlapsedges[idx];
}


/**
 * @brief returns whether vertices' incident overlaps are known
 *
 * Use \ref SCIPhypergraphComputeOverlaps to compute them.
 */
SCIP_Bool SCIPhypergraphHasVertexOverlaps(
   SCIP_HYPERGRAPH*      hypergraph          /**< The hypergraph. */
   )
{
   assert(hypergraph);

   return hypergraph->hasverticesoverlaps;
}

/** @brief returns an index for the first overlap containing \p vertex */
int SCIPhypergraphVertexOverlapsFirst(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex             /**< A vertex. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasverticesoverlaps);

   return hypergraph->verticesoverlapsbeg[vertex];
}

/** @brief returns an index beyond the last overlap incident to \p vertex */
int SCIPhypergraphVertexOverlapsBeyond(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   SCIP_HYPERGRAPH_VERTEX vertex             /**< A vertex. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasverticesoverlaps);

   return hypergraph->verticesoverlapsbeg[vertex + 1];
}

/**
 * @brief returns the overlap corresponding to \p idx that is incident to a vertex
 *
 * See \ref SCIPhypergraphVertexOverlapsFirst and \ref SCIPhypergraphVertexOverlapsBeyond to obtain such indices for a
 * vertex.
 */
SCIP_HYPERGRAPH_OVERLAP SCIPhypergraphVertexOverlapsGetAtIndex(
   SCIP_HYPERGRAPH*      hypergraph,         /**< The hypergraph. */
   int                   idx                 /**< Index. */
   )
{
   assert(hypergraph);
   assert(hypergraph->hasverticesoverlaps);

   return hypergraph->verticesoverlaps[idx];
}

#endif /* NDEBUG */
