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

/**@file   graph_path.c
 * @brief  Shortest path based graph algorithms for Steiner problems
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file encompasses various (heap-based) shortest path based algorithms including
 * Dijkstra's algorithm and Voronoi diagram algorithms
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include "portab.h"
#include "graph.h"



inline static SCIP_Real sdwalk_getdistnewEdge(
   const int*            prevedges,          /**< previous edges per node */
   const int*            nprevedges,         /**< number of previous edges per node */
   const SCIP_Real*      cost,               /**< cost */
   const SCIP_Real*      dist,               /**< distance */
   int                   k,                  /**< previous node  */
   int                   e,                  /**< current outgoing edge  */
   int                   maxnprevs           /**< maximum number of previous  */
)
{
   const int nprevs = nprevedges[k];
   SCIP_Real dist_e;

   /* ancestor list not full? */
   if( nprevs != maxnprevs + 1 )
   {
      int i;
      const int e2 = e / 2;
      assert(nprevs <= maxnprevs);

      /* check whether m is contained in ancestor list */
      for( i = 0; i < nprevs; i++ )
      {
         const int prevedge = prevedges[maxnprevs * k + i];

         if( e2 == prevedge )
            break;
      }

      /* e2 in list? */
      if( i != nprevs )
      {
         assert(e2 == prevedges[maxnprevs * k + i]);
         dist_e = dist[k];
      }
      else
         dist_e = dist[k] + cost[e];
   }
   else
   {
      dist_e = dist[k] + cost[e];
   }

   return dist_e;
}


inline static SCIP_Real sdwalk_getdistnewPrize(
   const int*            prevNPterms,        /**< previous np terminals per node */
   const int*            nprevNPterms,       /**< number of previous np terminals per node */
   const int*            termmark,           /**< terminal mark */
   const STP_Bool*       visited,            /**< visited */
   const SCIP_Real*      prize,              /**< prize */
   int                   k,                  /**< current node  */
   int                   m,                  /**< next node  */
   SCIP_Real             distnew,            /**< distance of m */
   int                   maxnprevs           /**< maximum number of previous  */
)
{
   SCIP_Real distnewP = distnew;

   assert(termmark[m] == 1 || termmark[m] == 2 );

   if( termmark[m] == 2 || !visited[m] )
   {
      distnewP = MAX(0.0, distnewP - prize[m]);
   }
   else
   {
      const int nprevs = nprevNPterms[k];

      /* ancestor list not full? */
      if( nprevs != maxnprevs + 1 )
      {
         int i;
         assert(nprevs <= maxnprevs);

         /* check whether m is contained in ancestor list */
         for( i = 0; i < nprevs; i++ )
         {
            const int prevterm = prevNPterms[maxnprevs * k + i];

            if( m == prevterm )
               break;
         }

         /* m not in list? */
         if( i == nprevs )
            distnewP = MAX(0.0, distnewP - prize[m]);
      }
   }

   return distnewP;
}



inline static SCIP_Bool sdwalk_conflict(
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   const int*            prevterms,          /**< previous terminals */
   const int*            nprevterms,         /**< number of previous terminals */
   const STP_Bool        nodevisited         /**< visited node already? */
   )
{
   const int nprevs = nprevterms[prednode];
   SCIP_Bool conflict = FALSE;

   assert(Is_term(g->term[node]));

   if( !nodevisited )
      return FALSE;

   if( nprevs > maxnprevs )
   {
      assert(nprevs == maxnprevs + 1);
      return TRUE;
   }

   for( int i = 0; i < nprevs; i++ )
   {
      const int prevterm = prevterms[maxnprevs * prednode + i];
      assert(prevterm >= 0);

      if( prevterm == node )
      {
         conflict = TRUE;
         break;
      }
   }

   return conflict;
}

inline static void sdwalk_update(
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms          /**< number of previous terminals */
   )
{
   const int predsize = nprevterms[prednode];
   const SCIP_Bool isterm = Is_term(g->term[node]);

   assert(predsize <= maxnprevs + 1);

   if( predsize == maxnprevs + 1 || (isterm && predsize == maxnprevs) )
   {
      nprevterms[node] = maxnprevs + 1;
   }
   else
   {
#ifndef NDEBUG
      for( int j = 0; j < predsize; j++ )
         assert(prevterms[maxnprevs * prednode + j] != node);
#endif

      for( int i = 0; i < predsize; i++ )
         prevterms[maxnprevs * node + i] = prevterms[maxnprevs * prednode + i];

      nprevterms[node] = predsize;

      if( isterm )
      {
         assert(predsize < maxnprevs);
         prevterms[maxnprevs * node + predsize] = node;
         nprevterms[node]++;
      }

      assert(nprevterms[node] <= maxnprevs);
   }
}

inline
static void sdwalk_updateCopy(
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   int*                  prev,               /**< previous data elements */
   int*                  nprev               /**< number of previous data elements */
   )
{
   const int predsize = nprev[prednode];

   assert(predsize <= maxnprevs);

   /* copy data from predecesseor */
   for( int i = 0; i < predsize; i++ )
      prev[maxnprevs * node + i] = prev[maxnprevs * prednode + i];

   nprev[node] = predsize;
}

static void sdwalk_update2(
   const int*            termmark,           /**< terminal mark */
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   edge,               /**< the edge to be updated */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   SCIP_Bool             clear,              /**< clear arrays */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms,         /**< number of previous terminals */
   int*                  prevNPterms,        /**< previous non-proper terminals */
   int*                  nprevNPterms,       /**< number of previous non-proper terminals */
   int*                  prevedges,          /**< previous edges */
   int*                  nprevedges          /**< number of previous edges */
   )
{
   int predsize = nprevterms[prednode];

   /*** 1. proper terminals ***/

   /* not enough space? */
   if( predsize == maxnprevs + 1 || (termmark[node] == 2 && predsize == maxnprevs) )
   {
      nprevterms[node] = maxnprevs + 1;
   }
   else
   {
#ifndef NDEBUG
      for( int j = 0; j < predsize; j++ )
         assert(prevterms[maxnprevs * prednode + j] != node);
#endif

      sdwalk_updateCopy(node, prednode, maxnprevs, prevterms, nprevterms);

      if( termmark[node] == 2 )
      {
         assert(predsize < maxnprevs);
         prevterms[maxnprevs * node + predsize] = node;
         nprevterms[node]++;
      }

      assert(nprevterms[node] <= maxnprevs);
   }


   /*** 2. edges ***/

   if( clear )
   {
      nprevNPterms[node] = 0;
      nprevedges[node] = 0;
      return;
   }

   predsize = nprevedges[prednode];

   if( predsize >= maxnprevs )
   {
      assert(predsize == maxnprevs || predsize == maxnprevs + 1);

      nprevedges[node] = maxnprevs + 1;
      nprevNPterms[node] = maxnprevs + 1;
      return;
   }
   assert(predsize < maxnprevs);

   sdwalk_updateCopy(node, prednode, maxnprevs, prevedges, nprevedges);

   prevedges[maxnprevs * node + predsize] = edge / 2;
   nprevedges[node]++;

   assert(nprevedges[node] <= maxnprevs);


   /*** 3. non-proper terminals ***/

   predsize = nprevNPterms[prednode];

   if( predsize == maxnprevs + 1 || (termmark[node] == 1 && predsize == maxnprevs) )
   {
      nprevNPterms[node] = maxnprevs + 1;
   }
   else
   {
      sdwalk_updateCopy(node, prednode, maxnprevs, prevNPterms, nprevNPterms);

      if( termmark[node] == 1 )
      {
         assert(predsize < maxnprevs);
         prevNPterms[maxnprevs * node + predsize] = node;
         nprevNPterms[node]++;
      }

      assert(nprevNPterms[node] <= maxnprevs);
   }
}

inline static void sdwalk_reset(
   int                   nvisits,            /**< number of visited nodes */
   const int*            visitlist,          /**< stores all visited nodes */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   STP_Bool*             visited             /**< stores whether a node has been visited */
)
{
   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      assert(node >= 0);

      visited[node] = FALSE;
      dist[node] = FARAWAY;
      state[node] = UNKNOWN;
   }
}


/*---------------------------------------------------------------------------*/
/*--- Name     : get NEAREST knot                                         ---*/
/*--- Function : Holt das oberste Element vom Heap (den Knoten mit der    ---*/
/*---            geringsten Entfernung zu den bereits Verbundenen)        ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzte Kanten              ---*/
/*--- Returns  : Nummer des bewussten Knotens                             ---*/
/*---------------------------------------------------------------------------*/
inline static int nearest(
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   const PATH* path)
{
   int   k;
   int   t;
   int   c;
   int   j;

   /* Heap shift down
    * (Oberstes Element runter und korrigieren)
    */
   k              = heap[1];
   j              = 1;
   c              = 2;
   heap[1]        = heap[(*count)--];
   state[heap[1]] = 1;

   if ((*count) > 2)
      if (LT(path[heap[3]].dist, path[heap[2]].dist))
         c++;

   while((c <= (*count)) && GT(path[heap[j]].dist, path[heap[c]].dist))
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c             += c;

      if ((c + 1) <= (*count))
         if (LT(path[heap[c + 1]].dist, path[heap[c]].dist))
            c++;
   }
   return(k);
}


/*---------------------------------------------------------------------------*/
/*--- Name     : get NEAREST knot                                         ---*/
/*--- Function : Holt das oberste Element vom Heap (den Knoten mit der    ---*/
/*---            geringsten Entfernung zu den bereits Verbundenen)        ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzte Kanten              ---*/
/*--- Returns  : Nummer des bewussten Knotens                             ---*/
/*---------------------------------------------------------------------------*/
inline static int nearestX(
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   const SCIP_Real* pathdist)
{
   int   k;
   int   t;
   int   c;
   int   j;
   int   dcount;

   /* Heap shift down
    * (Oberstes Element runter und korrigieren)
    */
   k              = heap[1];
   j              = 1;
   c              = 2;
   heap[1]        = heap[(*count)--];
   state[heap[1]] = 1;

   dcount = *count;

   if (dcount > 2)
      if (LT(pathdist[heap[3]], pathdist[heap[2]]))
         c++;

   while((c <= dcount) && GT(pathdist[heap[j]], pathdist[heap[c]]))
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c             += c;

      if ((c + 1) <= dcount)
         if (LT(pathdist[heap[c + 1]], pathdist[heap[c]]))
            c++;
   }
   return(k);
}


/*---------------------------------------------------------------------------*/
/*--- Name     : CORRECT heap                                             ---*/
/*--- Function : Setzt ein neues Element auf den Heap, bzw. korrigiert    ---*/
/*---            die Position eines vorhandenen Elementes                 ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzten Kanten,            ---*/
/*---            Neuer Knoten, Vorgaengerknoten, Kante von der man aus    ---*/
/*---            den neuen Knoten erreicht, Kosten der Kante,             ---*/
/*---            sowie Betriebsmodus                                      ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
inline static void correct(
   SCIP* scip,
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   PATH*  path,
   int    l,
   int    k,
   int    e,
   SCIP_Real cost,
   int    mode)
{
   int    t;
   int    c;
   int    j;

   path[l].dist = (mode == MST_MODE) ? cost : (path[k].dist + cost);
   path[l].edge = e;

   /* new node? */
   if( state[l] == UNKNOWN )
   {
      heap[++(*count)] = l;
      state[l]      = (*count);
   }

   /* Heap shift up */
   j = state[l];
   c = j / 2;
   while( (j > 1) && path[heap[c]].dist > path[heap[j]].dist )
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c              = j / 2;
   }
}


/*---------------------------------------------------------------------------*/
/*--- Name     : CORRECT heap                                             ---*/
/*--- Function : Setzt ein neues Element auf den Heap, bzw. korrigiert    ---*/
/*---            die Position eines vorhandenen Elementes                 ---*/
/*--- Parameter: Derzeitige Entfernungen und benutzten Kanten,            ---*/
/*---            Neuer Knoten, Vorgaengerknoten, Kante von der man aus    ---*/
/*---            den neuen Knoten erreicht, Kosten der Kante,             ---*/
/*---            sowie Betriebsmodus                                      ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
inline static void correctX(
   SCIP* scip,
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   SCIP_Real*  pathdist,
   int*   pathedge,
   int    l,
   int    k,
   int    e,
   SCIP_Real cost
   )
{
   int    t;
   int    c;
   int    j;

   pathdist[l] = (pathdist[k] + cost);

   if( pathedge != NULL )
      pathedge[l] = e;

   if (state[l] == UNKNOWN)
   {
      heap[++(*count)] = l;
      state[l]      = (*count);
   }

   /* Heap shift up
    */
   j = state[l];
   c = j / 2;

   while( (j > 1) && pathdist[heap[c]] > pathdist[heap[j]] )
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c              = j / 2;
   }
}


inline static void correctXwalk(
   SCIP* scip,
   int* heap,
   int* state,
   int* count,    /* pointer to store the number of elements on the heap */
   SCIP_Real*  pathdist,
   int    l,
   SCIP_Real newcost
   )
{
   int    t;
   int    c;
   int    j;

   pathdist[l] = newcost;

   if (state[l] == UNKNOWN)
   {
      heap[++(*count)] = l;
      state[l]      = (*count);
   }

   /* Heap shift up
    */
   j = state[l];
   c = j / 2;

   while( (j > 1) && pathdist[heap[c]] > pathdist[heap[j]] )
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c              = j / 2;
   }
}

static
void updatmaxprize(
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      orderedprizes,      /**< ordered prizes for (pseudo) terminals */
   const int*            orderedprizes_id,   /**< ordered prizes ids */
   const STP_Bool*       connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int                   node,               /**< current node */
   int*                  maxprizeidx_p,      /**< pointer to */
   SCIP_Real*            maxprizeval_p       /**< pointer to */
)
{
   int maxprizeidx = *maxprizeidx_p;
   assert(maxprizeidx >= 0 && maxprizeidx <= g->terms);

   /* sentinel? */
   if( orderedprizes_id[maxprizeidx] < 0 )
   {
      assert(orderedprizes_id[maxprizeidx] == -1);
      assert(*maxprizeval_p == 0.0);
      return;
   }

   assert(maxprizeidx < g->terms);

   /* is current node at the maximum? */
   if( node == orderedprizes_id[maxprizeidx] )
   {
      while( orderedprizes_id[maxprizeidx] >= 0 && connected[orderedprizes_id[maxprizeidx]] )
      {
         maxprizeidx++;
         assert(maxprizeidx <= g->terms);
      }

      *maxprizeidx_p = maxprizeidx;

      if( orderedprizes_id[maxprizeidx] < 0 )
         *maxprizeval_p = 0.0;
      else
         *maxprizeval_p = orderedprizes[maxprizeidx];
   }
}

/** adds element 'node' to heap */
void heap_add(
   const PATH* path,
   int node,               /* the node to be added */
   int* RESTRICT heap,     /* heap array */
   int* RESTRICT state,
   int* count             /* pointer to store the number of elements on the heap */
   )
{
   int c;
   int j;

   heap[++(*count)] = node;
   state[node]      = (*count);

   /* Heap shift up */
   j = state[node];
   c = j / 2;

   while((j > 1) && GT(path[heap[c]].dist, path[heap[j]].dist))
   {
      const int t = heap[c];
      heap[c] = heap[j];
      heap[j] = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j = c;
      c = j / 2;
   }
}


inline static void resetX(
   SCIP* scip,
   SCIP_Real* pathdist,
   int* heap,
   int* state,
   int* count,
   int node,
   SCIP_Real distnew
   )
{
   int c;
   int j;

   pathdist[node] = distnew;

   heap[++(*count)] = node;
   state[node] = (*count);

   /* heap shift up */
   j = state[node];
   c = j / 2;

   while( (j > 1) && pathdist[heap[c]] > pathdist[heap[j]] )
   {
      const int t = heap[c];
      heap[c] = heap[j];
      heap[j] = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j = c;
      c = j / 2;
   }
}

inline static void reset(
   SCIP* scip,
   PATH* path,
   int* heap,
   int* state,
   int* count,
   int    node
   )
{
   int    t;
   int    c;
   int    j;

   path[node].dist = 0.0;

   heap[++(*count)] = node;
   state[node]      = (*count);

   /* heap shift up */
   j = state[node];
   c = j / 2;

   while( (j > 1) && SCIPisGT(scip, path[heap[c]].dist, path[heap[j]].dist) )
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c              = j / 2;
   }
}

inline static void utdist(
   SCIP*            scip,
   const GRAPH*  g,
   PATH* path,
   SCIP_Real ecost,
   int* vbase,
   int k,
   int l,
   int k2,
   int shift,
   int nnodes
   )
{
   SCIP_Real dist;
   int vbk;
   int vbk2;

   if( Is_term(g->term[k]) )
      vbk = k;
   else
      vbk = vbase[k];

   if( l == 0 )
   {
      assert(shift == 0);

      dist = ecost;
      if( !Is_term(g->term[k]) )
         dist += path[k].dist;

      if( !Is_term(g->term[k2]) )
      {
         dist += path[k2].dist;
         vbk2 = vbase[k2];
      }
      else
      {
         vbk2 = k2;
      }

      if( SCIPisLT(scip, dist, path[vbk].dist) )
      {
         path[vbk].dist = dist;
         vbase[vbk] = vbk2;
         return;
      }
   }
   else
   {
      int max;
      int pos;
      int r;
      int s;
      int t;

      pos = vbk + shift;
      max = MIN((l + 1), 3);

      for( r = 0; r <= max; r++ )
      {
         if( Is_term(g->term[k2]) )
         {
            if( r == 0 )
               t = k2;
            else
               break;
         }
         else
         {
            t = vbase[k2 + (r * nnodes)];
         }
         for( s = 0; s < l; s++ )
            if( vbase[vbk + s * nnodes] == t )
               break;
         if( s < l || vbk == t )
            continue;

         dist = ecost;
         if( !Is_term(g->term[k]) )
            dist += path[k].dist;
         if( !Is_term(g->term[k2]) )
            dist += path[k2 + (r * nnodes)].dist;

         if( SCIPisLT(scip, dist, path[pos].dist) )
         {
            path[pos].dist = dist;
            vbase[pos] = t;
            return;
         }
      }
   }
   return;
}


/** connect given node to tree */
static inline
void stPcmwConnectNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   k,                  /**< the vertex */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      orderedprizes,      /**< ordered prizes for (pseudo) terminals */
   const int*            orderedprizes_id,   /**< ordered prizes ids */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int*                  maxprizeidx,        /**< index */
   SCIP_Real*            maxprizeval,        /**< value */
   int*                  count,              /**< for the heap */
   int*                  nterms              /**< terminal count */
)
{
   assert(k >= 0);

   connected[k] = TRUE;
   pathdist[k] = 0.0;
   updatmaxprize(g, orderedprizes, orderedprizes_id, connected, k, maxprizeidx, maxprizeval);
   (*nterms)++;

   assert(pathedge[k] != -1);

   /* connect k to current subtree */
   for( int node = g->tail[pathedge[k]]; !connected[node]; node = g->tail[pathedge[node]] )
   {
      connected[node] = TRUE;
      resetX(scip, pathdist, g->path_heap, g->path_state, count, node, 0.0);

      if( Is_pseudoTerm(g->term[node]) )
      {
         updatmaxprize(g, orderedprizes, orderedprizes_id, connected, node, maxprizeidx, maxprizeval);
         (*nterms)++;
      }

      assert(pathedge[node] != -1);
   }
}


/** initialize */
static inline
void stPcmwInit(
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int*                  npseudoterms        /**< number of pseudo terminals */
)
{
   const int nnodes = graph_get_nNodes(g);
   int ntermspos = 0;
   int* const RESTRICT state = g->path_state;

   for( int k = 0; k < nnodes; k++ )
   {
      g->mark[k] = ((g->grad[k] > 0) && !Is_term(g->term[k]));
      state[k] = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;

      if( Is_pseudoTerm(g->term[k]) )
         ntermspos++;
   }

   *npseudoterms = ntermspos;
}


/** initialize */
static inline
void stRpcmwInit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   int*                  nrealterms          /**< number of real terminals */
)
{
   const int nnodes = graph_get_nNodes(g);
   int nrterms = 0;
   int* const RESTRICT state = g->path_state;

   /* unmark dummy terminals */
   graph_pc_markOrgGraph(scip, g);
   assert(graph_pc_knotIsFixedTerm(g, g->source));

   for( int k = 0; k < nnodes; k++ )
   {
      state[k] = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;

      if( graph_pc_knotIsFixedTerm(g, k) )
      {
         assert(g->mark[k]);
         nrterms++;
         assert(!Is_pseudoTerm(g->term[k]));
      }
   }

   *nrealterms = nrterms;
}


/*---------------------------------------------------------------------------*/
/*--- Name     : INIT shortest PATH algorithm                             ---*/
/*--- Function : Initialisiert den benoetigten Speicher fuer die          ---*/
/*---            kuerzeste Wege Berechnung                                ---*/
/*--- Parameter: Graph in dem berechnet werden soll.                      ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
SCIP_RETCODE graph_path_init(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
   )
{
   assert(g != NULL);
   assert(g->path_heap  == NULL);
   assert(g->path_state == NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->path_heap), g->knots + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->path_state), g->knots) );

   return SCIP_OKAY;
}

/*---------------------------------------------------------------------------*/
/*--- Name     : EXIT shortest PATH algorithm                             ---*/
/*--- Function : Gibt den bei der initialisierung angeforderten Speicher  ---*/
/*---            wieder frei                                              ---*/
/*--- Parameter: Keine                                                    ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
void graph_path_exit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
   )
{
   assert(g != NULL);
   assert(g->path_heap  != NULL);
   assert(g->path_state != NULL);

   SCIPfreeMemoryArray(scip, &(g->path_state));
   SCIPfreeMemoryArray(scip, &(g->path_heap));
}

/*---------------------------------------------------------------------------*/
/*--- Name     : FIND shortest PATH / minimum spanning tree               ---*/
/*--- Function : Dijkstras Algorithmus zur bestimmung eines kuerzesten    ---*/
/*---            Weges in einem gerichteten Graphen.                      ---*/
/*--- Parameter: Graph, Nummer des Startknotens und Betriebsmodus         ---*/
/*---            (Kuerzeste Wege oder minimal spannender Baum)            ---*/
/*--- Returns  : Liefert einen Vektor mit einem PATH Element je Knotem.   ---*/
/*----           Setzt das .dist Feld zu jedem Knoten auf die Entfernung  ---*/
/*---            zum Startknoten, sowie das .prev Feld auf den Knoten von ---*/
/*---            dem man zum Knoten gekommen ist. Im MST Modus steht im   ---*/
/*---            .dist Feld die Entfernung zum naechsten Knoten.          ---*/
/*---------------------------------------------------------------------------*/
void graph_path_exec(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          p,                  /**< graph data structure */
   const int             mode,               /**< shortest path (FSP_MODE) or minimum spanning tree (MST_MODE)? */
   int                   start,              /**< start vertex */
   const SCIP_Real*      cost,               /**< edge costs */
   PATH*                 path                /**< shortest paths data structure */
   )
{
   int* heap;
   int* state;
   int k;
   const int nnodes = p->knots;
   int count;

   assert(scip      != NULL);
   assert(p      != NULL);
   assert(start  >= 0);
   assert(start  <  p->knots);
   assert((mode  == FSP_MODE) || (mode == MST_MODE));
   assert(p->path_heap   != NULL);
   assert(p->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);

   /* no nodes?, return*/
   if( nnodes == 0 )
      return;

   heap = p->path_heap;
   state = p->path_state;

   /* initialize */
   for( int i = 0; i < nnodes; i++ )
   {
      state[i]     = UNKNOWN;
      path[i].dist = FARAWAY + 1;
      path[i].edge = UNKNOWN;
   }
   /* add first node to heap */
   k            = start;
   path[k].dist = 0.0;

   if( nnodes > 1 )
   {
      count       = 1;
      heap[count] = k;
      state[k]    = count;

      while( count > 0 )
      {
         /* get nearest labeled node */
         k = nearest(heap, state, &count, path);

         /* mark as scanned */
         state[k] = CONNECT;

         for( int i = p->outbeg[k]; i >= 0; i = p->oeat[i] )
         {
            const int m = p->head[i];

            assert(i != EAT_LAST);

            /* node not scanned and valid? */
            if( state[m] )
            {
               /* closer than previously and valid? */
               if( path[m].dist > ((mode == MST_MODE) ? cost[i] : (path[k].dist + cost[i])) && p->mark[m] )
                  correct(scip, heap, state, &count, path, m, k, i, cost[i], mode);
            }
         }
      }
   }
}

/** limited Dijkstra, stopping at terminals */
void graph_sdPaths(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PATH*                 path,               /**< shortest paths data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned during SP calculation */
   int*                  memlbl,             /**< array to save labelled nodes */
   int*                  nlbl,               /**< number of labelled nodes */
   int                   tail,               /**< tail of the edge */
   int                   head,               /**< head of the edge */
   int                   limit               /**< maximum number of edges to consider during execution */
   )
{
   int count;
   int nchecks;
   const int limit1 = limit / 2;

   assert(g      != NULL);
   assert(heap   != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(nlbl   != NULL);
   assert(memlbl != NULL);
   assert(limit1 >= 0);

   *nlbl = 0;

   if( g->grad[tail] == 0 || g->grad[head] == 0 )
      return;

   assert(g->mark[head] && g->mark[tail]);

   count = 0;
   nchecks = 0;
   path[tail].dist = 0.0;
   state[tail] = CONNECT;
   memlbl[(*nlbl)++] = tail;

   if( g->stp_type != STP_MWCSP )
      g->mark[head] = FALSE;

   for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && (distlimit >= cost[e]) )
      {
         assert(SCIPisGT(scip, path[m].dist, path[tail].dist + cost[e]));

         /* m labelled the first time */
         memlbl[(*nlbl)++] = m;
         correct(scip, heap, state, &count, path, m, tail, e, cost[e], FSP_MODE);

         if( nchecks++ > limit1 )
            break;
      }
   }
   g->mark[head] = TRUE;

   while( count > 0 && nchecks <= limit )
   {
      /* get nearest labelled node */
      const int k = nearest(heap, state, &count, path);

      /* scanned */
      state[k] = CONNECT;

      /* distance limit reached? */
      if( SCIPisGT(scip, path[k].dist, distlimit) )
         break;

      /* stop at terminals */
      if( Is_term(g->term[k]) || k == head )
         continue;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];
         if( state[m] && g->mark[m] && (distlimit >= cost[e]) && (path[m].dist > path[k].dist + cost[e]) )
         {
            /* m labelled for the first time? */
            if( state[m] == UNKNOWN )
               memlbl[(*nlbl)++] = m;
            correct(scip, heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
         }
         if( nchecks++ > limit )
            break;
      }
   }
}


/** limited Dijkstra, stopping at terminals */
void graph_sdStar(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             with_zero_edges,    /**< telling name */
   int                   star_root,          /**< root of the start */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   int*                  star_base,          /**< star base node, must be initially set to SDSTAR_BASE_UNSET */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool*             visited,            /**< stores whether a node has been visited */
   SCIP_Bool*            success             /**< will be set to TRUE iff at least one edge can be deleted */
   )
{
   int nchecks;
   int nstarhits;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const RESTRICT range_csr = dcsr->range;
   const int* const RESTRICT head_csr = dcsr->head;
   const SCIP_Real* const RESTRICT cost_csr = dcsr->cost;
   const int star_degree = range_csr[star_root].end - range_csr[star_root].start;
   SCIP_Real distlimit;

   assert(dcsr && g && dist && visitlist && nvisits && visited && dheap && success);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(g->mark[star_root] && star_degree >= 1);
   assert(dheap->size == 0);
   assert(edgelimit >= 1);

   *nvisits = 0;
   *success = FALSE;

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      assert(dist[k] == FARAWAY);
      assert(star_base[k] == SDSTAR_BASE_UNSET);
      assert(state[k] == UNKNOWN);
   }
#endif

   distlimit = 0.0;
   dist[star_root] = 0.0;
   state[star_root] = CONNECT;
   visitlist[(*nvisits)++] = star_root;

   for( int e = range_csr[star_root].start, end = range_csr[star_root].end; e < end; e++ )
   {
      const int m = head_csr[e];

      assert(g->mark[m]);
      assert(!visited[m]);

      visitlist[(*nvisits)++] = m;
      visited[m] = TRUE;
      dist[m] = cost_csr[e];
      star_base[m] = m;
      graph_heap_correct(m, cost_csr[e], dheap);

      if( cost_csr[e] > distlimit )
         distlimit = cost_csr[e];
   }


   nchecks = 0;
   nstarhits = 0;

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labeled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      assert(k != star_root);
      assert(state[k] == CONNECT);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( with_zero_edges && star_base[k] == k )
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];
         assert(g->mark[m] && star_base[k] >= 0);

         if( state[m] != CONNECT )
         {
            const SCIP_Real distnew = dist[k] + cost_csr[e];

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               if( star_base[m] == m )
                  nstarhits++;

               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m);
            }
            else if( SCIPisEQ(scip, distnew, dist[m]) && star_base[m] == m )
            {
               if( with_zero_edges && star_base[k] == star_base[m] )
                  continue;

               assert(visited[m]);
               nstarhits++;

               assert(star_base[m] != star_base[k]);

               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m);
            }

            /* all star nodes hit already? */
            if( nstarhits == star_degree )
            {
               nchecks = edgelimit + 1;
               break;
            }
         }
         nchecks++;
      }
   }

  *success = (nstarhits > 0);
}


/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalks(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise)*/
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int count;
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;

   assert(g != NULL);
   assert(heap != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);

   *nvisits = 0;

   if( g->grad[start] == 0 || g->grad[end] == 0 )
      return FALSE;

   assert(g->mark[start] && g->mark[end]);

   count = 0;
   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   g->mark[start] = FALSE;
   g->mark[end] = FALSE;

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
         {
            const SCIP_Real newcost = MAX(cost[e] - g->prize[m], 0.0);
            correctXwalk(scip, heap, state, &count, dist, m, newcost);
         }
         else
         {
            correctXwalk(scip, heap, state, &count, dist, m, cost[e]);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }
   g->mark[end] = TRUE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( (state[m] != CONNECT) && g->mark[m] )
         {
            SCIP_Real distnew = dist[k] + cost[e];

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( termmark[m] != 0 )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               correctXwalk(scip, heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;
   return success;
}



/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalks_csr(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise)*/
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   RANGE* const RESTRICT range_csr = dcsr->range;
   int* const RESTRICT head_csr = dcsr->head;
   SCIP_Real* const RESTRICT cost_csr = dcsr->cost;

   assert(dcsr && g && dist && visitlist && nvisits && visited && dheap);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(g->grad[start] != 0 && g->grad[end] != 0);
   assert(g->mark[start] && g->mark[end]);
   assert(dheap->size == 0);

   *nvisits = 0;

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
      assert(state[k] == UNKNOWN);
#endif

   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   for( int e = range_csr[start].start; e < range_csr[start].end; e++ )
   {
      const int m = head_csr[e];
      assert(g->mark[m]);

      if( SCIPisLE(scip, cost_csr[e], distlimit) && m != end )
      {
         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
         {
            const SCIP_Real newcost = MAX(cost_csr[e] - g->prize[m], 0.0);

            dist[m] = newcost;
            graph_heap_correct(m, newcost, dheap);
         }
         else
         {
            dist[m] = cost_csr[e];
            graph_heap_correct(m, cost_csr[e], dheap);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];

         if( state[m] != CONNECT && m != start )
         {
            SCIP_Real distnew = dist[k] + cost_csr[e];

            assert(g->mark[m]);

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( termmark[m] != 0 )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
         nchecks++;
      }
   }

   return success;
}


/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalksTriangle(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   const int*            stateprev,          /**< state of previous run or NULL */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            prizeoffset,        /**< array for storing prize offset or NULL */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   RANGE* const RESTRICT range_csr = dcsr->range;
   int* const RESTRICT head_csr = dcsr->head;
   SCIP_Real* const RESTRICT cost_csr = dcsr->cost;

   assert(dcsr && g && dist && visitlist && visited && dheap);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(g->grad[start] != 0 && g->grad[end] != 0);
   assert(g->mark[start] && g->mark[end]);
   assert(dheap->size == 0);

   *nvisits = 0;

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
      assert(state[k] == UNKNOWN);
#endif

   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   for( int e = range_csr[start].start; e < range_csr[start].end; e++ )
   {
      const int m = head_csr[e];
      assert(g->mark[m]);

      if( SCIPisLE(scip, cost_csr[e], distlimit) && m != end )
      {
         assert(!visited[m]);

         if( stateprev && stateprev[m] == CONNECT )
            continue;

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
         {
            const SCIP_Real newcost = MAX(cost_csr[e] - g->prize[m], 0.0);

            dist[m] = newcost;
            graph_heap_correct(m, newcost, dheap);

            if( prizeoffset )
            {
               if( g->prize[m] > cost_csr[e] )
               {
                  prizeoffset[m] = cost_csr[e];
                  assert(SCIPisZero(scip, newcost));
               }
               else
                  prizeoffset[m] = g->prize[m];
            }
         }
         else
         {
            dist[m] = cost_csr[e];
            graph_heap_correct(m, cost_csr[e], dheap);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];

         if( state[m] != CONNECT )
         {
            SCIP_Real distnew;

            assert(m != start);

            if( stateprev && stateprev[m] == CONNECT )
                continue;

            distnew = dist[k] + cost_csr[e];

            assert(g->mark[m]);

            if( distnew > distlimit )
               continue;

            if( termmark[m] != 0 )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               if( prizeoffset && termmark[m] != 0 )
               {
                  const SCIP_Real distnew0 = dist[k] + cost_csr[e];

                  if( g->prize[m] > distnew0 )
                  {
                     prizeoffset[m] = distnew0;
                     assert(SCIPisZero(scip, distnew));
                  }
                  else
                     prizeoffset[m] = g->prize[m];
               }

               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
         nchecks++;
      }
   }

   return success;
}

/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalksExt(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms,         /**< number of previous terminals */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int count;
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;

   assert(g != NULL);
   assert(heap != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);

   *nvisits = 0;

   if( g->grad[start] == 0 || g->grad[end] == 0 )
      return FALSE;

   assert(g->mark[start] && g->mark[end]);

   count = 0;
   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   g->mark[start] = FALSE;
   g->mark[end] = FALSE;

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;
         sdwalk_update(g, m, start, maxnprevs, prevterms, nprevterms);

         if( Is_term(g->term[m]) )
         {
            const SCIP_Real newcost = MAX(cost[e] - g->prize[m], 0.0);
            correctXwalk(scip, heap, state, &count, dist, m, newcost);
         }
         else
         {
            correctXwalk(scip, heap, state, &count, dist, m, cost[e]);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }
   assert(nprevterms[start] == 0);

   g->mark[end] = TRUE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( g->mark[m] )
         {
            SCIP_Real distnew = dist[k] + cost[e];

            assert(state[m] != CONNECT);

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( Is_term(g->term[m]) )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               const SCIP_Bool mvisited = visited[m];
               if( !mvisited )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               if( Is_term(g->term[m]) && sdwalk_conflict(g, m, k, maxnprevs, prevterms, nprevterms, mvisited) )
                  continue;

               sdwalk_update(g, m, k, maxnprevs, prevterms, nprevterms);
               correctXwalk(scip, heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;
   return success;
}



/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalksExt2(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms,         /**< number of previous terminals */
   int*                  prevNPterms,        /**< previous non-proper terminals */
   int*                  nprevNPterms,       /**< number of previous non-proper terminals */
   int*                  prevedges,          /**< previous edges */
   int*                  nprevedges,         /**< number of previous edges */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int count;
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;

   assert(g != NULL);
   assert(heap != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);

   *nvisits = 0;

   if( g->grad[start] == 0 || g->grad[end] == 0 )
      return FALSE;

   assert(g->mark[start] && g->mark[end]);

   count = 0;
   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   g->mark[start] = FALSE;
   g->mark[end] = FALSE;

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         SCIP_Real distnew = cost[e];

         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
            distnew = MAX(distnew - g->prize[m], 0.0);

         sdwalk_update2(termmark, m, start, e, maxnprevs, SCIPisZero(scip, distnew),
               prevterms, nprevterms, prevNPterms, nprevNPterms, prevedges, nprevedges);
         correctXwalk(scip, heap, state, &count, dist, m, distnew);

         if( ++nchecks > edgelimit1 )
            break;
      }
   }
   assert(nprevterms[start] == 0);

   g->mark[end] = TRUE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( g->mark[m] )
         {
            SCIP_Real distnew = sdwalk_getdistnewEdge(prevedges, nprevedges, cost, dist, k, e, maxnprevs);

            assert(state[m] != CONNECT);

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( termmark[m] != 0 )
               distnew = sdwalk_getdistnewPrize(prevNPterms, nprevNPterms, termmark, visited, g->prize, k, m, distnew, maxnprevs);

            if( SCIPisLT(scip, distnew, dist[m]) )
            {
               const SCIP_Bool mvisited = visited[m];
               if( !mvisited )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               /* continue if m is proper terminals and is on the walk to k */
               if( termmark[m] == 2 && sdwalk_conflict(g, m, k, maxnprevs, prevterms, nprevterms, mvisited) )
                  continue;

               sdwalk_update2(termmark, m, k, e, maxnprevs, SCIPisZero(scip, distnew),
                     prevterms, nprevterms, prevNPterms, nprevNPterms, prevedges, nprevedges);
               correctXwalk(scip, heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;
   return success;
}


/** modified Dijkstra along walks for PcMw */
SCIP_Bool graph_sdWalksConnected(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   const SCIP_Real*      cost,               /**< edge costs */
   const STP_Bool*       endpoint,           /**< stores whether search should be ended at vertex */
   int                   start,              /**< start vertex */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited,            /**< stores whether a node has been visited */
   SCIP_Bool             resetarrays         /**< should arrays be reset? */
   )
{
   int* const heap = g->path_heap;
   int* const state = g->path_state;
   int count;
   int nchecks;
   const SCIP_Real prize = g->prize[start];

   assert(heap != NULL);
   assert(state != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(Is_term(g->term[start]));
   assert(g->grad[start] > 0);
   assert(g->mark[start]);

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      assert(state[k] == UNKNOWN);
      assert(visited[k] == FALSE);
      assert(dist[k] == FARAWAY);
   }
#endif

   *nvisits = 0;
   nchecks = 0;
   count = 1;
   heap[count] = start;
   state[start] = count;
   dist[start] = 0.0;
   visitlist[(*nvisits)++] = start;
   g->mark[start] = FALSE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(SCIPisLE(scip, dist[k], prize));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( (state[m] != CONNECT) && g->mark[m] )
         {
            SCIP_Real distnew = dist[k] + cost[e];

            if( SCIPisGT(scip, distnew, prize) )
               continue;

            if( termmark[m] != 0 )
               distnew -= g->prize[m];

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( endpoint != NULL && endpoint[m] )
               {
                  g->mark[start] = TRUE;
                  if( resetarrays )
                     sdwalk_reset(*nvisits, visitlist, dist, state, visited);

                  return TRUE;
               }

               correctXwalk(scip, heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;

   if( resetarrays )
      sdwalk_reset(*nvisits, visitlist, dist, state, visited);

   return FALSE;
}



/** limited Dijkstra for PCSTP, taking terminals into account */
void graph_path_PcMwSd(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PATH*                 path,               /**< shortest paths data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int*                  pathmaxnode,        /**< maximum weight node on path */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned during SP calculation */
   int*                  stateblock,         /**< array to indicate whether a node has been scanned during previous SP calculation */
   int*                  memlbl,             /**< array to save labelled nodes */
   int*                  nlbl,               /**< number of labelled nodes */
   int                   tail,               /**< tail of the edge */
   int                   head,               /**< head of the edge */
   int                   limit               /**< maximum number of edges to consider during execution */
   )
{
   const int limit1 = limit / 2;
   int count;
   int nchecks;
   const SCIP_Bool block = (stateblock != NULL);

   assert(limit > 0);
   assert(g      != NULL);
   assert(heap   != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(nlbl   != NULL);
   assert(memlbl != NULL);
   assert(g->prize != NULL);
   assert(pathmaxnode != NULL);

   *nlbl = 0;

   if( g->grad[tail] == 0 || g->grad[head] == 0 )
      return;

   assert(g->mark[head] && g->mark[tail]);

   nchecks = 0;
   count = 0;
   path[tail].dist = 0.0;
   state[tail] = CONNECT;
   memlbl[(*nlbl)++] = tail;

   if( g->stp_type != STP_MWCSP )
      g->mark[head] = FALSE;

   for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         assert(SCIPisGT(scip, path[m].dist, path[tail].dist + cost[e]));

         /* m labelled the first time */
         memlbl[(*nlbl)++] = m;
         correct(scip, heap, state, &count, path, m, tail, e, cost[e], FSP_MODE);

         if( nchecks++ > limit1 )
            break;
      }
   }

   g->mark[head] = TRUE;

   /* main loop */
   while( count > 0 )
   {
      const int k = nearest(heap, state, &count, path);
      SCIP_Real maxweight = pathmaxnode[k] >= 0 ? g->prize[pathmaxnode[k]] : 0.0;

      assert(k != tail);
      assert(maxweight >= 0);
      assert(SCIPisLE(scip, path[k].dist - maxweight, distlimit));

      /* scanned */
      state[k] = CONNECT;

      /* stop at other end */
      if( k == head )
         continue;

      if( Is_term(g->term[k]) && g->prize[k] > maxweight && distlimit >= path[k].dist )
      {
         pathmaxnode[k] = k;
         maxweight = g->prize[k];
      }

      /* stop at node scanned in first run */
      if( block && stateblock[k] == CONNECT )
         continue;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( state[m] && g->mark[m] && path[m].dist > (path[k].dist + cost[e])
               && distlimit >= (path[k].dist + cost[e] - maxweight) )
         {
            if( state[m] == UNKNOWN ) /* m labeled for the first time? */
               memlbl[(*nlbl)++] = m;

            pathmaxnode[m] = pathmaxnode[k];
            correct(scip, heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
         }
         if( nchecks++ > limit )
            break;
      }
   }
}



/** Dijkstra's algorithm starting from node 'start' */
void graph_path_execX(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   start,              /**< start vertex */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge            /**< predecessor edge array (on vertices) */
   )
{
   int   k;
   int   m;
   int   i;
   int   count;
   int   nnodes;
   int* heap;
   int* state;

   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(g->path_heap   != NULL);
   assert(g->path_state  != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(cost   != NULL);

   nnodes = g->knots;

   if( nnodes == 0 )
      return;

   heap = g->path_heap;
   state = g->path_state;

   for( i = nnodes - 1; i >= 0; --i )
   {
      state[i]     = UNKNOWN;
      pathdist[i] = FARAWAY;
      pathedge[i] = -1;
   }

   k = start;
   pathdist[k] = 0.0;

   if( nnodes > 1 )
   {
      count       = 1;
      heap[count] = k;
      state[k]    = count;

      while( count > 0 )
      {
         k = nearestX(heap, state, &count, pathdist);

         state[k] = CONNECT;

         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            if( state[m] && g->mark[m] && SCIPisGT(scip, pathdist[m], (pathdist[k] + cost[i])) )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, i, cost[i]);
         }
      }
   }
}

/** Dijkstra on incoming edges until root is reached */
void graph_path_invroot(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   start,              /**< start vertex */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge            /**< predecessor edge array (on vertices) */
   )
{
   SCIP_Real rootdist;
   int   k;
   int   m;
   int   i;
   int   count;
   int   nnodes;
   int* heap;
   int* state;

   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(g->path_heap   != NULL);
   assert(g->path_state  != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(cost   != NULL);

   nnodes = g->knots;

   if( nnodes == 0 )
      return;

   heap = g->path_heap;
   state = g->path_state;
   rootdist = FARAWAY;

   for( i = nnodes - 1; i >= 0; --i )
   {
      state[i]     = UNKNOWN;
      pathdist[i] = FARAWAY;
      pathedge[i] = -1;
   }

   k = start;
   pathdist[k] = 0.0;

   if( nnodes > 1 )
   {
      int root = g->source;

      count       = 1;
      heap[count] = k;
      state[k]    = count;

      while( count > 0 )
      {
         k = nearestX(heap, state, &count, pathdist);

         state[k] = CONNECT;

         if( k == root )
            rootdist = pathdist[k];
         else if( SCIPisGT(scip, pathdist[k], rootdist) )
            break;

         for( i = g->inpbeg[k]; i != EAT_LAST; i = g->ieat[i] )
         {
            m = g->tail[i];

            if( state[m] && g->mark[m] && SCIPisGT(scip, pathdist[m], (pathdist[k] + cost[i])) )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, i, cost[i]);
         }
      }
   }
}


/** extension heuristic */
void graph_path_st_pcmw_extendOut(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   start,              /**< start */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   SCIP_Real*            dist,               /**< distances array */
   int*                  pred,               /**< predecessor node */
   STP_Bool*             connected_out,      /**< array for internal stuff */
   DHEAP*                dheap,              /**< Dijkstra heap */
   SCIP_Bool*            success             /**< extension successful? */
   )
{
   int* const state = dheap->position;
   const CSR* const csr = g->csr_storage;
   const int* const start_csr = csr->start;
   const int* const head_csr = csr->head;
   const SCIP_Real* const cost_csr = csr->cost;
   SCIP_Real outerprofit;
   SCIP_Real outermaxprize;
   const int nnodes = g->knots;

   assert(csr && g && dist && dheap && pred && connected && connected_out);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(!connected[start]);

   *success = FALSE;
   outermaxprize = 0.0;

   for( int k = 0; k < nnodes; k++ )
   {
      state[k] = UNKNOWN;
      dist[k] = FARAWAY;
      connected_out[k] = FALSE;
#ifndef NDEBUG
      pred[k] = -1;
#endif

      if( !connected[k] && Is_term(g->term[k]) && g->prize[k] > outermaxprize && k != start )
         outermaxprize = g->prize[k];
   }

   graph_heap_clean(FALSE, dheap);

   dist[start] = 0.0;
   graph_heap_correct(start, 0.0, dheap);
   connected_out[start] = TRUE;
   outerprofit = g->prize[start];

   for( int rounds = 0; rounds < 2 && *success == FALSE; rounds++  )
   {
      if( rounds == 1 )
      {
         /* no improvement in last round? */
         if( !SCIPisGT(scip, outerprofit, g->prize[start]) )
            break;

         if( dheap->size > 0 )
            graph_heap_clean(TRUE, dheap);

         assert(dheap->size == 0);

         /* insert outer tree vertices into heap */
         for( int k = 0; k < nnodes; k++ )
         {
            if( connected_out[k] )
            {
               dist[k] = 0.0;
               graph_heap_correct(k, 0.0, dheap);
            }
            else
               dist[k] = FARAWAY;
         }
      }

      while( dheap->size > 0 )
      {
         /* get nearest labelled node */
         const int k = graph_heap_deleteMinReturnNode(dheap);
         const int k_start = start_csr[k];
         const int k_end = start_csr[k + 1];
         state[k] = UNKNOWN;

         /* if k is positive vertex and close enough, connect k to current subtree */
         if( (connected[k] && SCIPisGT(scip, outerprofit, dist[k]))
               || (!connected[k] && !connected_out[k] && Is_term(g->term[k]) && SCIPisGE(scip, g->prize[k], dist[k])) )
         {
            assert(k != start);
            assert(pred[k] != -1);
            assert(!connected_out[k] || !connected[k]);

            outerprofit += g->prize[k] - dist[k];
            connected_out[k] = TRUE;
            dist[k] = 0.0;

            assert(SCIPisGE(scip, outerprofit, g->prize[start]) || connected[k]);

            /* connect k to current subtree */
            for( int node = pred[k]; !connected_out[node]; node = pred[node] )
            {
               connected_out[node] = TRUE;
               dist[node] = 0.0;
               graph_heap_correct(node, 0.0, dheap);

               if( Is_term(g->term[node]) )
                  outerprofit += g->prize[node];

               assert(state[node]);
               assert(pred[node] >= 0);
            }

            if( connected[k] )
            {
               *success = TRUE;
               break;
            }
         }
         else if( outerprofit + outermaxprize < dist[k] )
         {
            assert(*success == FALSE);
            break;
         }

         /* correct incident nodes */
         for( int e = k_start; e < k_end; e++ )
         {
            const int m = head_csr[e];
            const SCIP_Real distnew = dist[k] + cost_csr[e];

            if( distnew < dist[m] )
            {
               dist[m] = distnew;
               pred[m] = k;
               graph_heap_correct(m, distnew, dheap);
            }
         }
      }
   }

   if( *success )
   {
      for( int k = 0; k < nnodes; k++ )
         if( connected_out[k] )
            connected[k] = TRUE;
   }
}

/** Find a directed tree rooted in node 'start' and spanning all terminals */
void graph_path_st(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edgecosts */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   int   count;
   int   nnodes;
   int*  heap;
   int*  state;

   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(cost   != NULL);
   assert(connected != NULL);

   count = 0;
   nnodes = g->knots;
   heap = g->path_heap;
   state = g->path_state;

   /* initialize */
   for( int k = 0; k < nnodes; k++ )
   {
      state[k]     = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;
   }

   pathdist[start] = 0.0;
   connected[start] = TRUE;

   if( nnodes > 1 )
   {
      int nterms = 0;

      if( Is_term(g->term[start]) )
         nterms++;

      /* add start vertex to heap */
      count = 1;
      heap[count] = start;
      state[start] = count;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         const int k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* k is terminal an not connected yet? */
         if( Is_term(g->term[k]) && k != start )
         {
            assert(pathedge[k] >= 0 && !connected[k]);

            connected[k] = TRUE;
            pathdist[k] = 0.0;

            /* connect k to current solution */
            for( int node = g->tail[pathedge[k]]; !connected[node]; node = g->tail[pathedge[node]] )
            {
               assert(pathedge[node] != -1);
               assert(!Is_term(g->term[node]));

               connected[node] = TRUE;
               resetX(scip, pathdist, heap, state, &count, node, 0.0);
            }

            /* have all terminals been reached? */
            if( ++nterms == g->terms )
               break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int m = g->head[e];

            assert(state[m]);

            /* is m not connected, allowed and closer (as close)? */
            if( !connected[m] && pathdist[m] > (pathdist[k] + cost[e]) && g->mark[m] )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}



/** Find a tree rooted in node 'start' and connecting
 *  positive vertices as long as this is profitable.
 *  Note that this function overwrites g->mark.
 *  */
void graph_path_st_pcmw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      orderedprizes,      /**< ordered prizes for (pseudo) terminals */
   const int*            orderedprizes_id,   /**< ordered prizes ids */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   SCIP_Bool             costIsBiased,       /**< is cost biased? */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   const int nnodes = g->knots;
   int* const RESTRICT heap = g->path_heap;
   int* const RESTRICT state = g->path_state;
   int ntermspos = -1;

   assert(scip && g && orderedprizes && orderedprizes_id && cost && prize && pathdist && pathedge && connected);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(g->extended);
   assert(graph_pc_isPcMw(g) && !graph_pc_isRootedPcMw(g));

   /* initialize */
   stPcmwInit(g, pathdist, pathedge, connected, &ntermspos);

   assert(g->mark[start]);
   assert(ntermspos >= 0);

   pathdist[start] = 0.0;
   connected[start] = TRUE;

   if( nnodes > 1 )
   {
      int count = 1;
      int maxprizeidx = 0;
      SCIP_Real maxprizeval = orderedprizes[0];
      int nterms = 0;
      const SCIP_Bool isPc = graph_pc_isPc(g);

      if( Is_pseudoTerm(g->term[start]) )
      {
         nterms++;
         updatmaxprize(g, orderedprizes, orderedprizes_id, connected, start, &maxprizeidx, &maxprizeval);
      }

      /* add start vertex to heap */
      heap[count] = start;
      state[start] = count;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         SCIP_Bool connectK = FALSE;
         const int k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* if k is positive vertex and close enough, connect k to current subtree */
         if( !connected[k] && Is_pseudoTerm(g->term[k]) )
         {
            connectK = (prize[k] >= pathdist[k]);

            assert(k != start);

            if( !connectK )
            {
               SCIP_Real prizesum = 0.0;

               for( int node = g->tail[pathedge[k]]; !connected[node]; node = g->tail[pathedge[node]] )
               {
                  if( Is_pseudoTerm(g->term[node]) )
                  {
                     prizesum += prize[node];
                  }
                  else if( isPc && !costIsBiased && graph_pc_knotIsNonLeafTerm(g, node) )
                  {
                     prizesum += prize[node];
                  }
               }

               assert(prizesum >= 0.0 && SCIPisLT(scip, prizesum, FARAWAY));

               connectK = (prize[k] + prizesum >= pathdist[k]);
            }

            if( connectK )
            {
               stPcmwConnectNode(scip, k, g, orderedprizes, orderedprizes_id, pathdist, pathedge, connected,
                     &maxprizeidx, &maxprizeval, &count, &nterms);

               assert(nterms <= ntermspos);

               /* have all biased terminals been connected? */
               if( nterms == ntermspos )
               {
                  SCIPdebugMessage("all terms reached \n");
                  break;
               }
            }
         }

         if( !connectK && pathdist[k] > maxprizeval )
         {
            break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int m = g->head[e];

            assert(state[m] && e != EAT_LAST);

            /* is m not connected, allowed and closer? */
            if( g->mark[m] && !connected[m] && pathdist[m] > (pathdist[k] + cost[e]) )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      } /* while( count > 0 ) */
   } /* nnodes > 1*/
}

/** Reduce given solution
 *  Note that this function overwrites g->mark.
 *  */
void graph_path_st_pcmw_reduce(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            tmpnodeweight,      /**< node weight array */
   int*                  result,             /**< incoming arc array */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   assert(tmpnodeweight != NULL);
   assert(result   != NULL);
   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(cost   != NULL);
   assert(connected != NULL);

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      if( result[e] == CONNECT )
      {
         const int head = g->head[e];

         if( Is_term(g->term[head]) )
            continue;

         graph_path_st_pcmw_reduce(scip, g, cost, tmpnodeweight, result, head, connected);

         assert(connected[head]);

         if( SCIPisGE(scip, cost[e], tmpnodeweight[head]) )
         {
            connected[head] = FALSE;
            result[e] = UNKNOWN;
         }
         else
         {
            tmpnodeweight[start] += tmpnodeweight[head] - cost[e];
         }
      }
   }

   SCIPfreeBufferArray(scip, &tmpnodeweight);
}



/** Find a tree rooted in node 'start' and connecting
 *  all positive vertices.
 *  Note that this function overwrites g->mark.
 *  */
void graph_path_st_pcmw_full(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   int* const heap = g->path_heap;
   int* const state = g->path_state;
   const int nnodes = g->knots;

   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   assert(g      != NULL);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(cost   != NULL);
   assert(connected != NULL);
   assert(graph_pc_isPcMw(g) && !graph_pc_isRootedPcMw(g));

   /* initialize */
   for( int k = nnodes - 1; k >= 0; --k )
   {
      g->mark[k] = ((g->grad[k] > 0) && !Is_term(g->term[k]));
      state[k]     = UNKNOWN;
      pathdist[k] = FARAWAY;
      pathedge[k] = -1;
      connected[k] = FALSE;
   }

   pathdist[start] = 0.0;
   connected[start] = TRUE;

   if( nnodes > 1 )
   {
      int count = 1;
      int nterms = 0;

      /* add start vertex to heap */
      heap[count] = start;
      state[start] = count;

      if( Is_pseudoTerm(g->term[start]) )
         nterms++;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         const int k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* if k is unconnected positive vertex, connect its path to current subtree */
         if( !connected[k] && Is_pseudoTerm(g->term[k]) )
         {
            connected[k] = TRUE;
            pathdist[k] = 0.0;

            assert(k != start);
            assert(pathedge[k] != -1);

            for( int node = g->tail[pathedge[k]]; !connected[node]; node = g->tail[pathedge[node]] )
            {
               assert(!Is_anyTerm(g->term[node]));

               connected[node] = TRUE;
               resetX(scip, pathdist, heap, state, &count, node, 0.0);

               assert(pathedge[node] != -1);
            }

            /* have all terminals been reached? */
            if( ++nterms == g->terms - 1 )
               break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            int m = g->head[e];

            assert(state[m]);

            /* is m not connected, allowed and closer (as close)? */
            if( !connected[m] && g->mark[m] && SCIPisGT(scip, pathdist[m], (pathdist[k] + cost[e])) )
               correctX(scip, heap, state, &count, pathdist, pathedge, m, k, e, cost[e]);
         }
      }
   }
}


/** greedy extension of a given tree for PC or MW problems */
void graph_path_st_pcmw_extend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Bool             breakearly,         /**< finish computation early if no profitable extension possible? */
   PATH*                 path,               /**< shortest paths data structure */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   SCIP_Bool*            extensions          /**< extensions performed? */
   )
{
   SCIP_Real maxprize;
   int   count;
   int   nnodes;
   int   nstnodes;
   int   outerterms;
   int*  heap;
   int*  state;

   assert(path   != NULL);
   assert(g      != NULL);
   assert(cost   != NULL);
   assert(connected != NULL);
   assert(g->extended);

   maxprize = 0.0;
   count = 0;
   nstnodes = 0;
   nnodes = g->knots;
   heap = g->path_heap;
   state = g->path_state;
   *extensions = FALSE;
   outerterms = 0;

   /* initialize */
   for( int k = 0; k < nnodes; k++ )
   {
      g->mark[k] = ((g->grad[k] > 0) && !Is_term(g->term[k]));
      if( connected[k] && g->mark[k] )
      {
         /* add node to heap */
         nstnodes++;
         if( nnodes > 1 )
            heap[++count] = k;

         state[k]     = count;
         path[k].dist = 0.0;
         assert(path[k].edge != UNKNOWN || k == g->source);
      }
      else
      {
         state[k]     = UNKNOWN;
         path[k].dist = FARAWAY;

         if( Is_pseudoTerm(g->term[k]) && g->mark[k] )
         {
            outerterms++;
            if( g->prize[k] > maxprize )
            {
               maxprize = g->prize[k];
            }
         }
      }

      if( !connected[k] )
         path[k].edge = UNKNOWN;
   }

   /* nothing to extend? */
   if( nstnodes == 0 )
      return;

   if( nnodes > 1 )
   {
      int nterms = 0;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         const int k = nearest(heap, state, &count, path);
         state[k] = UNKNOWN;

         /* if k is positive vertex and close enough (or fixnode), connect its path to current subtree */
         if( !connected[k] && Is_pseudoTerm(g->term[k]) && SCIPisGE(scip, g->prize[k], path[k].dist) )
         {
            int node;

            nterms++;
            *extensions = TRUE;
            connected[k] = TRUE;
            path[k].dist = 0.0;

            assert(path[k].edge >= 0);
            node = g->tail[path[k].edge];

            while( !connected[node] )
            {
               assert(path[node].edge != UNKNOWN);
               connected[node] = TRUE;
               reset(scip, path, heap, state, &count, node);
               assert(state[node]);

               if( Is_pseudoTerm(g->term[node]) )
                  nterms++;

               node = g->tail[path[node].edge];
            }

            assert(path[node].dist == 0.0);
            assert(nterms <= outerterms);

            /* have all terminals been reached? */
            if( nterms == outerterms )
               break;
         }
         else if( breakearly && SCIPisGT(scip, path[k].dist, maxprize) )
         {
            break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int m = g->head[e];

            assert(state[m]);
            assert(e != EAT_LAST);

            /* is m not connected, allowed and closer (as close)? */

            if( !connected[m] && path[m].dist > (path[k].dist + cost[e]) && g->mark[m] )
               correct(scip, heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
         }
      }
   }
}


/** greedy extension of a given tree for PC or MW problems; path[i].edge needs to be initialized */
void graph_path_st_pcmw_extendBiased(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   PATH*                 path,               /**< shortest paths data structure with .edge initialized */
   STP_Bool*             connected,          /**< array to mark whether a vertex is part of computed Steiner tree */
   SCIP_Bool*            extensions          /**< extensions performed? */
   )
{
   SCIP_Real maxprize;
   int count;
   int nstnodes;
   int outermscount;
   const int nnodes = g->knots;
   int* const heap = g->path_heap;
   int* const state = g->path_state;

   assert(path   != NULL);
   assert(g      != NULL);
   assert(cost   != NULL);
   assert(connected != NULL);
   assert(g->extended);

   maxprize = 0.0;
   count = 0;
   nstnodes = 0;
   outermscount = 0;

   *extensions = FALSE;

   /* unmark dummy terminals */
   graph_pc_markOrgGraph(scip, g);
   assert(graph_pc_knotIsFixedTerm(g, g->source));

   /* initialize */
   for( int k = 0; k < nnodes; k++ )
   {
      state[k]     = UNKNOWN;
      path[k].dist = FARAWAY;

      if( !g->mark[k] )
         continue;

      if( connected[k] )
      {
         /* add node to heap */
         nstnodes++;
         if( nnodes > 1 )
            heap[++count] = k;

         state[k]     = count;
         path[k].dist = 0.0;
         assert(path[k].edge != UNKNOWN || k == g->source);
      }
      else if( Is_pseudoTerm(g->term[k]) )
      {
         assert(g->mark[k]);
         outermscount++;

         if( prize[k] > maxprize )
            maxprize = prize[k];

      }
   }

   /* with at least two nodes and at least one in the solution? */
   if( nnodes > 1 && nstnodes > 0 )
   {
      int nterms = 0;

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         const int k = nearest(heap, state, &count, path);
         state[k] = UNKNOWN;

         assert(g->mark[k]);

         /* if k is positive vertex and close enough (or fixnode), connect its path to current subtree */
         if( !connected[k] && Is_pseudoTerm(g->term[k]) && SCIPisGE(scip, prize[k], path[k].dist) )
         {
            int node;

            nterms++;
            *extensions = TRUE;
            connected[k] = TRUE;
            path[k].dist = 0.0;

            assert(path[k].edge >= 0);
            node = g->tail[path[k].edge];

            while( !connected[node] )
            {
               assert(g->mark[node]);
               assert(path[node].edge >= 0);
               connected[node] = TRUE;
               reset(scip, path, heap, state, &count, node);
               assert(state[node]);

               if( Is_pseudoTerm(g->term[node]) )
                  nterms++;

               node = g->tail[path[node].edge];
            }

            assert(path[k].dist == 0.0);

            assert(nterms <= outermscount);

            /* have all terminals been reached? */
            if( nterms == outermscount )
               break;
         }
         else if( path[k].dist > maxprize )
         {
            break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int m = g->head[e];
            assert(state[m]);

            if( connected[m] )
               continue;

            /* is m allowed and closer? */
            if( path[m].dist > (path[k].dist + cost[e]) && g->mark[m] )
               correct(scip, heap, state, &count, path, m, k, e, cost[e], FSP_MODE);
         }
      }
   }
}


/** Shortest path heuristic for the RMWCSP and RPCSPG
 * Find a directed tree rooted in node 'start' and connecting all terminals as well as all
 *  positive vertices (as long as this is profitable).
 *  */
void graph_path_st_rpcmw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   const SCIP_Real*      orderedprizes,      /**< ordered prizes for (pseudo) terminals */
   const int*            orderedprizes_id,   /**< ordered prizes ids */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      prize,              /**< (possibly biased) prize */
   SCIP_Real*            pathdist,           /**< distance array (on vertices) */
   int*                  pathedge,           /**< predecessor edge array (on vertices) */
   int                   start,              /**< start vertex */
   STP_Bool*             connected           /**< array to mark whether a vertex is part of computed Steiner tree */
   )
{
   const int nnodes = g->knots;
   int nrterms = -1;
   int* const heap = g->path_heap;
   int* const state = g->path_state;

   assert(scip && g && orderedprizes && orderedprizes_id && cost && prize && pathdist && pathedge && connected);
   assert(start  >= 0);
   assert(start  <  g->knots);
   assert(g->extended);
   assert(graph_pc_isRootedPcMw(g));

   stRpcmwInit(scip, g, pathdist, pathedge, connected, &nrterms);

   assert(nrterms >= 1);
   pathdist[start] = 0.0;
   connected[start] = TRUE;

   if( nnodes > 1 )
   {
      int count;
      const int nterms = g->terms;
      int termscount = 0;
      int rtermscount = 0;
      int maxprizeidx = 0;
      SCIP_Real maxprizeval = orderedprizes[0];

      /* add start vertex to heap */
      count = 1;
      heap[count] = start;
      state[start] = count;

      if( Is_anyTerm(g->term[start]) )
      {
         updatmaxprize(g, orderedprizes, orderedprizes_id, connected, start, &maxprizeidx, &maxprizeval);

         termscount++;
      }

      if( Is_term(g->term[start]) )
      {
         assert(graph_pc_knotIsFixedTerm(g, start));
         rtermscount++;
      }

      /* repeat until heap is empty */
      while( count > 0 )
      {
         /* get closest node */
         const int k = nearestX(heap, state, &count, pathdist);
         state[k] = UNKNOWN;

         /* if k is fixed terminal positive vertex and close enough, connect its path to current subtree */
         if( Is_anyTerm(g->term[k]) && (Is_term(g->term[k]) || SCIPisGE(scip, prize[k], pathdist[k]))
            && !connected[k] )
         {
            int node;

            assert(k != start);
            assert(pathedge[k] != -1);
            assert(!graph_pc_knotIsDummyTerm(g, k));
            assert(graph_pc_knotIsFixedTerm(g, k) || SCIPisGE(scip, prize[k], pathdist[k]));

            if( !graph_pc_knotIsNonLeafTerm(g, k) )
               termscount++;

            if( Is_term(g->term[k]) )
            {
               assert(graph_pc_knotIsFixedTerm(g, k));
               rtermscount++;
            }
            else if( Is_pseudoTerm(g->term[k]) )
            {
               updatmaxprize(g, orderedprizes, orderedprizes_id, connected, k, &maxprizeidx, &maxprizeval);
            }

            connected[k] = TRUE;
            pathdist[k] = 0.0;

            node = k;

            while( !connected[node = g->tail[pathedge[node]]] )
            {
               assert(pathedge[node] != -1);
               assert(!Is_term(g->term[node]));
               assert(!graph_pc_knotIsFixedTerm(g, node));
               assert(g->mark[node]);

               connected[node] = TRUE;
               resetX(scip, pathdist, heap, state, &count, node, 0.0);

               if( Is_pseudoTerm(g->term[node]) )
               {
                  termscount++;
                  updatmaxprize(g, orderedprizes, orderedprizes_id, connected, node, &maxprizeidx, &maxprizeval);
               }
            }

            assert(termscount <= nterms);

            /* have all terminals been reached? */
            if( termscount == nterms )
            {
               SCIPdebugMessage("all terminals reached \n");
               break;
            }
         }
         else if( rtermscount >= nrterms && pathdist[k] > maxprizeval )
         {
            SCIPdebugMessage("all fixed terminals reached \n");

            assert(rtermscount == nrterms);
            break;
         }

         /* update adjacent vertices */
         for( int e = g->outbeg[k]; e >= 0; e = g->oeat[e] )
         {
            const int head = g->head[e];

            assert(state[head]);

            /* is m not connected, allowed and closer (as close)? */
            if( !connected[head] && g->mark[head] && pathdist[head] > (pathdist[k] + cost[e]) )
               correctX(scip, heap, state, &count, pathdist, pathedge, head, k, e, cost[e]);
         }
      }
   }

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
   {
      if( graph_pc_knotIsFixedTerm(g, k) )
      {
         assert(connected[k]);
      }
   }
#endif
}


/** extend a Voronoi region until all neighbouring terminals are spanned */
SCIP_RETCODE graph_voronoiExtend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edgecosts */
   PATH*                 path,               /**< shortest paths data structure */
   SCIP_Real**           distarr,            /**< array to store distance from each node to its base */
   int**                 basearr,            /**< array to store the bases */
   int**                 edgearr,            /**< array to store the ancestor edge */
   STP_Bool*             termsmark,          /**< array to mark terminal */
   int*                  reachednodes,       /**< array to mark reached nodes */
   int*                  nreachednodes,      /**< pointer to number of reached nodes */
   int*                  nodenterms,         /**< array to store number of terimals to each node */
   int                   nneighbterms,       /**< number of neighbouring terminals */
   int                   base,               /**< voronoi base */
   int                   countex             /**< count of heap elements */
   )
{
   int   k;
   int   m;
   int   i;
   int* heap;
   int* state;
   int count = countex;

   assert(g      != NULL);
   assert(g->path_heap   != NULL);
   assert(g->path_state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);

   if( g->knots == 0 || nneighbterms <= 0 )
      return SCIP_OKAY;

   heap = g->path_heap;
   state = g->path_state;

   if( g->knots > 1 )
   {
      while( count > 0 && nneighbterms > 0 )
      {
         k = nearest(heap, state, &count, path);
         state[k] = CONNECT;
	 reachednodes[(*nreachednodes)++] = k;
         distarr[k][nodenterms[k]] = path[k].dist;
         edgearr[k][nodenterms[k]] = path[k].edge;
         basearr[k][nodenterms[k]] = base;

         nodenterms[k]++;

         if( termsmark[k] == TRUE )
         {
            termsmark[k] = FALSE;
            if( --nneighbterms == 0 )
            {
               while( count > 0 )
                  reachednodes[(*nreachednodes)++] = heap[count--];

               return SCIP_OKAY;
            }
         }
         else if( !Is_term(g->term[k]) )
         {
            /* iterate over all outgoing edges of vertex k */
            for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
            {
               m = g->head[i];

               /* check whether the path (to m) including (k, m) is shorter than the so far best known */
               if( (state[m]) && (g->mark[m]) && (GT(path[m].dist, (path[k].dist + cost[i]))) )
                  correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
            }
         }
      }
      assert(nneighbterms == 0);
   }
   return SCIP_OKAY;
}


/** build a Voronoi region, w.r.t. shortest paths, for a given set of bases */
void graph_voronoi(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   const STP_Bool*       basemark,           /**< array to indicate whether a vertex is a Voronoi base */
   int*                  vbase,              /**< voronoi base to each vertex */
   PATH*                 path                /**< path data struture (leading to respective Voronoi base) */
   )
{
   int* heap;
   int* state;
   int root;
   int count = 0;
   int nbases = 0;

   assert(g && scip && cost && costrev && basemark && vbase && path);
   assert(g->path_heap != NULL);
   assert(g->path_state != NULL);

   if( g->knots == 0 )
      return;

   root = g->source;
   heap = g->path_heap;
   state = g->path_state;

   /* initialize */
   for( int i = 0; i < g->knots; i++ )
   {
      /* set the base of vertex i */
      if( basemark[i] )
      {
         nbases++;
         if( g->knots > 1 )
            heap[++count] = i;
         vbase[i] = i;
         path[i].dist = 0.0;
         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i] = UNKNOWN;
      }
   }
   assert(nbases > 0);

   if( g->knots > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         const int k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all outgoing edges of vertex k */
         for( int i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            const int m = g->head[i];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + ((vbase[k] == root)? cost[i] : costrev[i])) )
            {
               assert(!basemark[m]);
               correct(scip, heap, state, &count, path, m, k, i, ((vbase[k] == root)? cost[i] : costrev[i]), FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
}

/** 2th next terminal to all non terminal nodes */
void graph_get2next(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data structure (leading to first and second nearest terminal) */
   int*                  vbase,              /**< first and second nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int j;
   int i;
   int e;
   int root;
   int count;
   int nnodes;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);
   assert(costrev   != NULL);

   root = g->source;
   count = 0;
   nnodes = g->knots;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      state[i] = CONNECT;

      /* copy of node i */
      k = i + nnodes;
      vbase[k] = UNKNOWN;
      state[k] = UNKNOWN;
      path[k].edge = UNKNOWN;
      path[k].dist = FARAWAY;
   }

   /* scan original nodes */
   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
         continue;

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         j = g->head[e];
         k = j + nnodes;
         if( !Is_term(g->term[j]) && SCIPisGT(scip, path[k].dist, path[i].dist + ((root == vbase[i])? cost[e] : costrev[e])) &&
            vbase[i] != vbase[j] && g->mark[j] )
         {
            correct(scip, heap, state, &count, path, k, i, e, ((root == vbase[i])? cost[e] : costrev[e]), FSP_MODE);
            vbase[k] = vbase[i];
         }
      }
   }

   if( nnodes > 1 )
   {
      int jc;
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as removed from heap */
         state[k] = UNKNOWN;

         assert(k - nnodes >= 0);
         /* iterate over all outgoing edges of vertex k */
         for( e = g->outbeg[k - nnodes]; e != EAT_LAST; e = g->oeat[e] )
         {
            j = g->head[e];
            if( Is_term(g->term[j]) || !g->mark[j] )
               continue;
            jc = j + nnodes;

            /* check whether the path (to j) including k is shorter than the so far best known */
            if( vbase[j] != vbase[k] && SCIPisGT(scip, path[jc].dist, path[k].dist + ((root == vbase[k])? cost[e] : costrev[e])) )
            {
               correct(scip, heap, state, &count, path, jc, k, e, (root == vbase[k])? cost[e] : costrev[e], FSP_MODE);
               vbase[jc] = vbase[k];
            }
         }
      }
   }
   return;
}

/* 3th next terminal to all non terminal nodes */
void graph_get3next(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*       costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data structure (leading to first, second and third nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int j;
   int i;
   int l;
   int v;
   int e;
   int root;
   int count = 0;
   int nnodes;
   int dnnodes;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);
   assert(costrev   != NULL);

   root = g->source;
   nnodes = g->knots;
   dnnodes = 2 * nnodes;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      /* copy of node i */
      k = i + dnnodes;
      vbase[k] = UNKNOWN;
      state[k] = UNKNOWN;
      path[k].edge = UNKNOWN;
      path[k].dist = FARAWAY;
   }

   /* scan original nodes */
   for( i = 0; i < nnodes; i++ )
   {
      state[i] = CONNECT;
      state[i + nnodes] = CONNECT;
      if( !g->mark[i] )
         continue;

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         j = g->head[e];
         k = j + dnnodes;
         if( !Is_term(g->term[j]) && g->mark[j] )
         {
            v = i;
            for( l = 0; l < 2; l++ )
            {
               if( SCIPisGT(scip, path[k].dist, path[v].dist + ((root == vbase[v])? cost[e] : costrev[e])) &&
                  vbase[v] != vbase[j] && vbase[v] != vbase[j + nnodes] )
               {
                  correct(scip, heap, state, &count, path, k, v, e, ((root == vbase[v])? cost[e] : costrev[e]), FSP_MODE);
                  vbase[k] = vbase[v];
               }
               v += nnodes;
            }
         }
      }
   }
   if( nnodes > 1 )
   {
      int jc;
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as removed from heap */
         state[k] = UNKNOWN;

         assert(k - dnnodes >= 0);
         /* iterate over all outgoing edges of vertex k */
         for( e = g->outbeg[k - dnnodes]; e != EAT_LAST; e = g->oeat[e] )
         {
            j = g->head[e];
            if( Is_term(g->term[j]) || !g->mark[j] )
               continue;
            jc = j + dnnodes;
            /* check whether the path (to j) including k is shorter than the so far best known */
            if( vbase[j] != vbase[k] && vbase[j + nnodes] != vbase[k]
               && SCIPisGT(scip, path[jc].dist, path[k].dist + ((root == vbase[k])? cost[e] : costrev[e])) ) /*TODO(state[jc])??*/
            {
               correct(scip, heap, state, &count, path, jc, k, e, (root == vbase[k])? cost[e] : costrev[e], FSP_MODE);
               vbase[jc] = vbase[k];
            }
         }
      }
   }
   return;
}


/* 4th next terminal to all non terminal nodes */
void graph_get4next(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data struture (leading to first, second and third nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int j;
   int i;
   int l;
   int v;
   int e;
   int root;
   int count = 0;
   int nnodes;
   int dnnodes;
   int tnnodes;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);
   assert(costrev   != NULL);

   root = g->source;
   nnodes = g->knots;
   dnnodes = 2 * nnodes;
   tnnodes = 3 * nnodes;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      /* copy of node i */
      k = i + tnnodes;
      vbase[k] = UNKNOWN;
      state[k] = UNKNOWN;
      path[k].edge = UNKNOWN;
      path[k].dist = FARAWAY;
   }

   /* scan original nodes */
   for( i = 0; i < nnodes; i++ )
   {
      state[i] = CONNECT;
      state[i + nnodes] = CONNECT;
      state[i + dnnodes] = CONNECT;
      if( !g->mark[i] )
         continue;

      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         j = g->head[e];
         k = j + tnnodes;
         if( !Is_term(g->term[j]) && g->mark[j] )
         {
            v = i;
            for( l = 0; l < 3; l++ )
            {
               if( SCIPisGT(scip, path[k].dist, path[v].dist + ((root == vbase[v])? cost[e] : costrev[e])) &&
                  vbase[v] != vbase[j] && vbase[v] != vbase[j + nnodes] && vbase[v] != vbase[j + dnnodes] )
               {
                  correct(scip, heap, state, &count, path, k, v, e, ((root == vbase[v])? cost[e] : costrev[e]), FSP_MODE);
                  vbase[k] = vbase[v];
               }
               v += nnodes;
            }
         }
      }
   }
   if( nnodes > 1 )
   {
      int jc;
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as removed from heap */
         state[k] = UNKNOWN;

         assert(k - tnnodes >= 0);
         /* iterate over all outgoing edges of vertex k */
         for( e = g->outbeg[k - tnnodes]; e != EAT_LAST; e = g->oeat[e] )
         {
            j = g->head[e];
            if( Is_term(g->term[j]) || !g->mark[j] )
               continue;
            jc = j + tnnodes;
            /* check whether the path (to j) including k is shorter than the so far best known */
            if( vbase[j] != vbase[k] && vbase[j + nnodes] != vbase[k] && vbase[j + dnnodes] != vbase[k]
               && SCIPisGT(scip, path[jc].dist, path[k].dist + ((root == vbase[k])? cost[e] : costrev[e])) ) /*TODO(state[jc])??*/
            {
               correct(scip, heap, state, &count, path, jc, k, e, (root == vbase[k])? cost[e] : costrev[e], FSP_MODE);
               vbase[jc] = vbase[k];
            }
         }
      }
   }
   return;
}

/** build a voronoi region in presolving, w.r.t. shortest paths, for all terminals*/
void graph_get3nextTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   PATH*                 path3,              /**< path data struture (leading to first, second and third nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   assert(g      != NULL);
   assert(path3   != NULL);
   assert(cost   != NULL);
   assert(costrev   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   if( !graph_pc_isPcMw(g) )
      for( int k = 0; k < g->knots; k++ )
         g->mark[k] = (g->grad[k] > 0);

   /* build voronoi diagram */
   graph_voronoiTerms(scip, g, cost, path3, vbase, heap, state);

   /* get 2nd nearest terms */
   graph_get2next(scip, g, cost, costrev, path3, vbase, heap, state);

   /* get 3th nearest terms */
   graph_get3next(scip, g, cost, costrev, path3, vbase, heap, state);

   return;
}

/** build a voronoi region in presolving, w.r.t. shortest paths, for all terminals*/
void graph_get4nextTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data struture (leading to first, second, third and fouth nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   assert(g         != NULL);
   assert(path      != NULL);
   assert(cost      != NULL);
   assert(heap      != NULL);
   assert(state     != NULL);
   assert(costrev   != NULL);

   if( !graph_pc_isPcMw(g) )
      for( int k = 0; k < g->knots; k++ )
         g->mark[k] = (g->grad[k] > 0);

   /* build voronoi diagram */
   graph_voronoiTerms(scip, g, cost, path, vbase, heap, state);

   /* get 2nd nearest terms */
   graph_get2next(scip, g, cost, costrev, path, vbase, heap, state);

   /* get 3th nearest terms */
   graph_get3next(scip, g, cost, costrev, path, vbase, heap, state);

   /* get 4th nearest terms */
   graph_get4next(scip, g, cost, costrev, path, vbase, heap, state);

   return;
}


/** get 4 close terminals to each terminal */
SCIP_RETCODE graph_get4nextTTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   PATH*                 path,               /**< path data structure (leading to first, second, third and fourth nearest terminal) */
   int*                  vbase,              /**< first, second and third nearest terminal to each non terminal */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int* boundedges;
   int k;
   int e;
   int l;
   int k2;
   int bedge;
   int shift;
   int nnodes;
   int nboundedges;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &boundedges, g->edges) );

   shift = 0;
   nnodes = g->knots;

   if( g->stp_type != STP_PCSPG && g->stp_type != STP_RPCSPG )
      for( k = 0; k < g->knots; k++ )
         g->mark[k] = (g->grad[k] > 0);

   nboundedges = 0;
   for( k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] )
         continue;

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         k2 = g->head[e];
         if( !g->mark[k2] || k2 < k )
            continue;

         /* is e a boundary edge? */
         if( vbase[k] != vbase[k2] )
         {
            boundedges[nboundedges++] = e;
         }
      }
      if( Is_term(g->term[k]) )
      {
         path[k].dist = FARAWAY;
         vbase[k] = UNKNOWN;
      }

   }

   for( l = 0; l < 4; l++ )
   {
      for( e = 0; e < nboundedges; e++ )
      {
         bedge = boundedges[e];
         k = g->tail[bedge];
         k2 = g->head[bedge];
         utdist(scip, g, path, cost[bedge], vbase, k, l, k2, shift, nnodes);
         utdist(scip, g, path, cost[bedge], vbase, k2, l, k, shift, nnodes);
      }
      shift += nnodes;
   }

   SCIPfreeBufferArray(scip, &boundedges);

   return SCIP_OKAY;
}

/** build a Voronoi region in presolving, w.r.t. shortest paths, for all terminals */
void graph_voronoiTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   PATH*                 path,               /**< path data structure (leading to respective Voronoi base) */
   int*                  vbase,              /**< Voronoi base to each vertex */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int m;
   int count = 0;
   int nbases = 0;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   /* initialize */
   for( int i = 0; i < g->knots; i++ )
   {
      /* set the base of vertex i */
      if( Is_term(g->term[i]) && g->mark[i] )
      {
         nbases++;
         if( g->knots > 1 )
            heap[++count] = i;
         vbase[i] = i;
         path[i].dist = 0.0;
         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
      }
   }
   if( nbases == 0 )
      return;
   if( g->knots > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all outgoing edges of vertex k */
         for( int i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) && g->mark[m] )
            {
               correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
   return;
}


/** build a Voronoi region, w.r.t. shortest paths, for all positive vertices */
void graph_voronoiMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   PATH*                 path,               /**< path data structure (leading to respective Voronoi base) */
   int*                  vbase,              /**< Voronoi base to each vertex */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark the state of each node during calculation */
   )
{
   int k;
   int m;
   int i;
   int count = 0;
   int nbases = 0;
   int nnodes;

   assert(g      != NULL);
   assert(path   != NULL);
   assert(costrev   != NULL);
   assert(heap   != NULL);
   assert(state   != NULL);

   nnodes = g->knots;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      /* set the base of vertex i */
      if( Is_term(g->term[i]) && g->mark[i] )
      {
         nbases++;
         if( nnodes > 1 )
            heap[++count] = i;
         vbase[i] = i;
         path[i].dist = -g->prize[i];
         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
      }
   }
   if( nbases == 0 )
      return;
   if( nnodes > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + costrev[i]) && g->mark[m] && !Is_term(g->term[m]) )
            {
               assert(vbase[m] != m);
               correct(scip, heap, state, &count, path, m, k, i, costrev[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }
   return;
}


/** build a voronoi region, w.r.t. shortest paths, for all terminal and the distance */
SCIP_RETCODE graph_voronoiWithDist(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            distance,           /**< array storing path from a terminal over shortest
                                                incident edge to nearest terminal */
   int*                  minedgepred,        /**< shortest edge predecessor array */
   int*                  vbase,              /**< array containing Voronoi base to each node */
   int*                  minarc,             /**< array to mark whether an edge is one a path corresponding to 'distance' */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array indicating state of each vertex during calculation of Voronoi regions */
   int*                  distnode,           /**< array to store terminal corresponding to distance stored in distance array */
   PATH*                 path                /**< array containing Voronoi paths data */
   )
{
   SCIP_Real new;
   int e;
   int k;
   int m;
   int i;
   int pred;
   int count;
   int nbases;
   int nnodes;
   int prededge;

   assert(g        != NULL);
   assert(path     != NULL);
   assert(cost     != NULL);
   assert(heap     != NULL);
   assert(state    != NULL);
   assert(distance != NULL);

   count = 0;
   nbases = 0;
   nnodes = g->knots;

   /* initialize */
   for( i = 0; i < nnodes; i++ )
   {
      distance[i] = FARAWAY;
      if( distnode != NULL )
         distnode[i] = UNKNOWN;

      /* set the base of vertex i */
      if( Is_term(g->term[i]) && g->mark[i] )
      {
         nbases++;
         if( nnodes > 1 )
            heap[++count] = i;
         vbase[i] = i;
         state[i] = count;
         path[i].dist = 0.0;
      }
      else
      {
         vbase[i] = UNKNOWN;
         state[i] = UNKNOWN;
         path[i].dist = FARAWAY;
      }
      path[i].edge = UNKNOWN;
   }

   /* initialize predecessor array */

   for( e = 0; e < g->edges; e++ )
      minedgepred[e] = FALSE;

   for( k = 0; k < nbases; k++ )
      if( minarc[k] != UNKNOWN )
         minedgepred[minarc[k]] = TRUE;

   /* if graph contains more than one vertices, start main loop */
   if( nnodes > 1 )
   {
      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         prededge = path[k].edge;

         if( prededge != UNKNOWN )
         {
            pred = g->tail[prededge];

            assert(vbase[k] != UNKNOWN);
            assert(vbase[pred] != UNKNOWN);
            assert(vbase[pred] == vbase[k]);
            assert(g->head[prededge] == k);

            if( !Is_term(g->term[pred]) && g->mark[pred] )
            {
               assert(path[pred].edge != UNKNOWN);
               minedgepred[prededge] = minedgepred[path[pred].edge];
            }

         }

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            if( state[m] == CONNECT && vbase[m] != vbase[k] && g->mark[m] )
            {
               if( minedgepred[i] || (prededge != UNKNOWN && minedgepred[prededge] ) )
               {
                  new = path[k].dist + g->cost[i] + path[m].dist;
                  if( SCIPisLT(scip, new, distance[vbase[k]]) )
                  {
                     if( distnode != NULL )
                        distnode[vbase[k]] = vbase[m];

                     distance[vbase[k]] = new;
                  }
               }
               if( minedgepred[Edge_anti(i)] || (path[m].edge != UNKNOWN && minedgepred[path[m].edge] ) )
               {
                  new = path[m].dist + g->cost[i] + path[k].dist;
                  if( SCIPisLT(scip, new, distance[vbase[m]]) )
                  {
                     if( distnode != NULL )
                        distnode[vbase[m]] = vbase[k];

                     distance[vbase[m]] = new;
                  }
               }
            }

            /* Check whether the path (to m) including k is shorter than the so far best known.
               In case of equality, also update if k is sucessor of minedge. */
            if( state[m] && g->mark[m] &&
               (SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) ||
                  (prededge != UNKNOWN && minedgepred[prededge] && SCIPisEQ(scip, path[m].dist, path[k].dist + cost[i]))) )
            {
               correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** build voronoi regions, w.r.t. shortest paths, for all terminals and compute the radii */
SCIP_RETCODE graph_voronoiWithRadius(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   GRAPH*                adjgraph,           /**< graph data structure */
   PATH*                 path,               /**< array containing Voronoi paths data */
   SCIP_Real*            rad,                /**< array storing shortest way from a terminal out of its Voronoi region */
   SCIP_Real*            cost,               /**< edge costs */
   SCIP_Real*            costrev,            /**< reversed edge costs */
   int*                  vbase,              /**< array containing Voronoi base of each node */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark state of each node during calculation */
   )
{
   int* nodesid;
   int count = 0;
   int nterms = 0;
   const int root = graph->source;
   const int nnodes = graph->knots;
   const STP_Bool mw = (graph->stp_type == STP_MWCSP);
   const STP_Bool pc = ((graph->stp_type == STP_PCSPG) || (graph->stp_type == STP_RPCSPG));

   assert(graph != NULL);
   assert(heap   != NULL);
   assert(state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(costrev   != NULL);
   assert(rad != NULL);
   assert(vbase != NULL);
   assert(!graph->extended);

   if( graph->terms == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodesid, nnodes) );

   /* initialize */
   for( int i = 0; i < nnodes; i++ )
   {
      rad[i] = FARAWAY;

      /* set the base of vertex i */
      if( Is_term(graph->term[i]) && graph->mark[i] )
      {
         if( nnodes > 1 )
            heap[++count] = i;

         if( !mw && adjgraph != NULL )
         {
            adjgraph->mark[adjgraph->knots] = TRUE;
            graph_knot_add(adjgraph, -1);
         }

         nodesid[i] = nterms++;
         vbase[i] = i;

         if( mw )
            assert(SCIPisGE(scip, graph->prize[i], 0.0));
         if( mw )
            path[i].dist = -graph->prize[i];
         else
            path[i].dist = 0.0;

         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         nodesid[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
      }
   }

   if( nnodes > 1 )
   {
      SCIP_Real ecost;

      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         const int k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate over all outgoing edges of vertex k */
         for( int i = graph->outbeg[k]; i != EAT_LAST; i = graph->oeat[i] )
         {
            const int m = graph->head[i];
            const int vbm = vbase[m];
            const int vbk = vbase[k];

            if( state[m] == CONNECT && vbm != vbk && graph->mark[m] )
            {
               assert(graph->mark[vbm]);
               assert(graph->mark[vbk]);
               if( mw )
               {
                  if( SCIPisGT(scip, rad[vbk], path[k].dist) )
                     rad[vbk] = path[k].dist;
                  if( SCIPisGT(scip, rad[vbm], path[m].dist) )
                     rad[vbm] = path[m].dist;
               }
               else
               {

                  if( SCIPisGT(scip, rad[vbk], path[k].dist + ((vbk == root) ? cost[i] : costrev[i])) )
                     rad[vbk] = path[k].dist + ((vbk == root) ? cost[i] : costrev[i]);

                  if( SCIPisGT(scip, rad[vbm], path[m].dist + ((vbm == root) ? costrev[i] : cost[i])) )
                     rad[vbm] = path[m].dist + ((vbm == root) ? costrev[i] : cost[i]);

                  if( adjgraph != NULL )
                  {
                     int ne;

                     if( graph->stp_type == STP_DHCSTP )
                     {
                        SCIP_Real c1;
                        SCIP_Real c2;

                        if( m == root )
                           c1 = path[m].dist + costrev[i];
                        else
                           c1 = path[m].dist + cost[i];
                        if( k == root )
                           c2 = path[k].dist + cost[i];
                        else
                           c2 = path[k].dist + costrev[i];

                        if( SCIPisGT(scip, c1, c2) )
                           ecost = c2;
                        else
                           ecost = c1;
                     }
                     else
                     {
                        if( SCIPisGT(scip, path[m].dist, path[k].dist) )
                           ecost = path[k].dist + cost[i];
                        else
                           ecost = path[m].dist + cost[i];
                     }

                     if( pc )
                     {
                        if( SCIPisGT(scip, ecost, graph->prize[vbm]) && root != vbm && !graph_pc_knotIsFixedTerm(graph, vbm) )
                           ecost = graph->prize[vbm];

                        if( SCIPisGT(scip, ecost, graph->prize[vbk]) && root != vbk && !graph_pc_knotIsFixedTerm(graph, vbk) )
                           ecost = graph->prize[vbk];
                     }

                     /* find edge in adjgraph */
                     for( ne = adjgraph->outbeg[nodesid[vbk]]; ne != EAT_LAST; ne = adjgraph->oeat[ne] )
                        if( adjgraph->head[ne] == nodesid[vbm] )
                           break;

                     /* edge exists? */
                     if( ne != EAT_LAST )
                     {
                        assert(ne >= 0);
                        assert(adjgraph->head[ne] == nodesid[vbm]);
                        assert(adjgraph->tail[ne] == nodesid[vbk]);
                        if( SCIPisGT(scip, adjgraph->cost[ne], ecost) )
                        {
                           adjgraph->cost[ne] = ecost;
                           adjgraph->cost[Edge_anti(ne)] = ecost;
                        }
                     }
                     else
                     {
                        graph_edge_add(scip, adjgraph, nodesid[vbm], nodesid[vbk], ecost, ecost);
                     }
                  }
               }
            }

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( state[m] && graph->mark[m] && !Is_term(graph->term[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + ((vbk == root)? cost[i] : costrev[i])) )
            {
               correct(scip, heap, state, &count, path, m, k, i, ((vbk == root)? cost[i] : costrev[i]), FSP_MODE);
               vbase[m] = vbk;
            }
         }
      }
   }
   SCIPfreeBufferArray(scip, &nodesid);

   return SCIP_OKAY;
}

/** Build partition of graph for MWCSP into regions around the positive vertices.
 * Store the weight of a minimum weight center-boundary path for each region
 * in the radius array (has to be reverted to obtain the final r() value).
 */
void graph_voronoiWithRadiusMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   PATH*                 path,               /**< array containing graph decomposition data */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            radius,             /**< array storing shortest way from a positive vertex out of its region */
   int*                  vbase,              /**< array containing base of each node */
   int*                  heap,               /**< array representing a heap */
   int*                  state               /**< array to mark state of each node during calculation */
   )
{
   SCIP_Real* nodeweight;
   int i;
   int k;
   int m;
   int count;
   int nbases;
   int nnodes;
   int nterms;

   assert(g != NULL);
   assert(heap   != NULL);
   assert(state  != NULL);
   assert(path   != NULL);
   assert(cost   != NULL);
   assert(radius != NULL);
   assert(vbase != NULL);

   count = 0;
   nbases = 0;
   nnodes = g->knots;
   nterms = g->terms - 1;
   nodeweight = g->prize;

   assert(nodeweight != NULL);

   if( nnodes == 0 || nterms <= 0 )
      return;

   assert(!g->mark[g->source]);

   /* initialize data */
   for( i = 0; i < nnodes; i++ )
   {
      radius[i] = FARAWAY;

      /* set the base of vertex i */
      if( Is_term(g->term[i]) && g->mark[i] )
      {
         nbases++;
         if( nnodes > 1 )
            heap[++count] = i;
         vbase[i] = i;
         path[i].dist = 0.0;
         path[i].edge = UNKNOWN;
         state[i] = count;
      }
      else
      {
         vbase[i] = UNKNOWN;
         path[i].dist = FARAWAY;
         path[i].edge = UNKNOWN;
         state[i]     = UNKNOWN;
      }
   }

   if( nbases == 0 )
      return;

   if( nnodes > 1 )
   {
      int basem;
      int basek;

      /* until the heap is empty */
      while( count > 0 )
      {
         /* get the next (i.e. a nearest) vertex of the heap */
         k = nearest(heap, state, &count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;
         basek = vbase[k];

         assert(g->mark[basek]);

         /* iterate over all outgoing edges of vertex k */
         for( i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            m = g->head[i];

            if( !(g->mark[m]) )
               continue;

            basem = vbase[m];

            /* are both m and k finally in different regions? */
            if( basem != basek && (state[m] == CONNECT || vbase[m] == m ||
               SCIPisGE(scip, path[k].dist, nodeweight[basek])) )
            {
               if( SCIPisGT(scip, radius[basek], path[k].dist) )
                  radius[basek] = path[k].dist;
               if( (state[m] == CONNECT || vbase[m] == m) && SCIPisGT(scip, radius[basem], path[m].dist) )
               {
                  assert(g->mark[basem]);
                  radius[basem] = path[m].dist;
               }
            }

            /* is the distance of vertex k smaller than the weight of its base node and smaller than the temp. radius? */
            if( SCIPisLT(scip, path[k].dist, nodeweight[basek]) && SCIPisLT(scip, path[k].dist, radius[basek]) )
            {
               /* check whether the path (to m) including k is shorter than the so far best known */
               if( (state[m]) && SCIPisGT(scip, path[m].dist, path[k].dist + cost[i]) )
               {
                  assert(vbase[m] != m);
                  correct(scip, heap, state, &count, path, m, k, i, cost[i], FSP_MODE);
                  vbase[m] = basek;
               }
            }
         }
      }
   }
   return;
}

/** repair a Voronoi diagram for a given set of base nodes */
void graph_voronoiRepair(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs (possibly biased) */
   const SCIP_Real*      cost_org,           /**< original edge costs (only needed for PC) */
   int*                  nheapelems,         /**< pointer to number of heap elements */
   int*                  vbase,              /**< array containing Voronoi base of each node */
   PATH*                 path,               /**< Voronoi paths data struture */
   int*                  newedge,            /**< the new edge */
   int                   crucnode,           /**< the current crucial node */
   UF*                   uf                  /**< union find data structure */
   )
{
   const int nnodes = graph_get_nNodes(g);
   int boundaryedge = UNKNOWN;
   SCIP_Real boundarydist = FARAWAY;

   assert(cost && scip && nheapelems && newedge && vbase && uf);

   *newedge = UNKNOWN;

   if( nnodes > 1 )
   {
      int* const heap = g->path_heap;
      int* const state = g->path_state;
      const SCIP_Bool isPcMw = graph_pc_isPcMw(g);
      const SCIP_Bool isPc = graph_pc_isPc(g);

      assert(heap && state);

      /* until the heap is empty */
      while( *nheapelems > 0 )
      {
         /* get the next (i.e. a nearest) vertex from the heap */
         const int k = nearest(heap, state, nheapelems, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         assert(vbase[k] >= 0);
         assert(g->mark[k] && g->mark[vbase[k]]);

         /* iterate over all outgoing edges of vertex k */
         for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int m = g->head[e];

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( (state[m]) && (SCIPisGT(scip, path[m].dist, path[k].dist + cost[e])) )
            {
               assert(g->mark[m]);
               correct(scip, heap, state, nheapelems, path, m, k, e, cost[e], FSP_MODE);
               vbase[m] = vbase[k];
            }
            /* check whether there is a better new boundary edge adjacent to vertex k */
            else if( state[m] == CONNECT && g->mark[m] )
            {
               const int node1 = SCIPStpunionfindFind(uf, vbase[m]);
               const int node2 = SCIPStpunionfindFind(uf, vbase[k]);

               if( (node1 == crucnode) == (node2 == crucnode) )
                  continue;

               if( !g->mark[node1] || !g->mark[node2] || !g->mark[vbase[m]] )
                  continue;

               /* now to the actual checks */

               if( boundaryedge == UNKNOWN )
               {
                  boundaryedge = e;
               }
               else
               {
                  SCIP_Real dist_new = path[k].dist + path[m].dist;

                  if( isPcMw )
                  {
                     if( isPc )
                     {
                        assert(cost_org);
                        dist_new += cost_org[e];
                     }
                  }
                  else
                  {
                     assert(!cost_org);
                     dist_new += cost[e];
                  }

                  if( dist_new < boundarydist )
                  {
                     boundarydist = dist_new;
                     boundaryedge = e;
                  }
               }
            } /* check for boundary edge */
         }
      }
   }

   *newedge = boundaryedge;
}


/** repair the Voronoi diagram for a given set nodes */
void graph_voronoiRepairMult(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const STP_Bool*       nodesmark,          /**< array to mark temporarily discarded nodes */
   int* RESTRICT         count,              /**< pointer to number of heap elements */
   int* RESTRICT         vbase,              /**< array containing Voronoi base of each node */
   int* RESTRICT         boundedges,         /**< boundary edges */
   int* RESTRICT         nboundedges,        /**< number of boundary edges */
   UF* RESTRICT          uf,                 /**< union find data structure */
   PATH* RESTRICT        path                /**< Voronoi paths data structure */
   )
{
   assert(scip && g && cost && nodesmark && count && vbase && boundedges && nboundedges && uf && path);
   assert(g->path_heap && g->path_state);

   if( g->knots > 1 )
   {
      int node1;
      int node2;
      int* const heap = g->path_heap;
      int* const state = g->path_state;

      /* until the heap is empty */
      while( (*count) > 0 )
      {
         /* get the next (i.e. a nearest) vertex from the heap */
         const int k = nearest(heap, state, count, path);

         /* mark vertex k as scanned */
         state[k] = CONNECT;

         /* iterate all outgoing edges of vertex k */
         for( int i = g->outbeg[k]; i != EAT_LAST; i = g->oeat[i] )
         {
            const int m = g->head[i];

            if( !(g->mark[m]) )
               continue;

            /* check whether the path (to m) including k is shorter than the so far best known */
            if( state[m] && (SCIPisGT(scip, path[m].dist, path[k].dist + cost[i])) )
            {
               correct(scip, heap, state, count, path, m, k, i, cost[i], FSP_MODE);
               vbase[m] = vbase[k];
            }
            /* check whether there is a new boundary edge adjacent to vertex k */
            else if( (state[m] == CONNECT) && ((node1 = SCIPStpunionfindFind(uf, vbase[m])) != (node2 = SCIPStpunionfindFind(uf, vbase[k])))
               && g->mark[node1] && g->mark[node2] && (nodesmark[node1] || nodesmark[node2]) )
            {
               boundedges[(*nboundedges)++] = i;
            }
         }
      }
   }
}
