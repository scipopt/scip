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

/**@file   dijkstra_bh.c
 * @brief  C implementation of Dijkstra's algorithm using binary heaps
 * @author Thorsten Koch
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "dijkstra_bh.h"


/** Check whether the data structures of the graph are valid */
DIJKSTRA_Bool Dijsktra_graphIsValid(
  const Dijkstra_Graph*  G                   /**< directed graph */
  )
{
   unsigned int i;
   unsigned int k;
   unsigned int count = 0;

   if ( G == NULL || G->outbeg == NULL || G->outcnt == NULL || G->weight == NULL || G->head == NULL )
      abort();

   for (i = 0; i < G->nodes; ++i)
   {
      for (k = G->outbeg[i]; k < G->outbeg[i] + G->outcnt[i]; ++k)
      {
         if (G->head[k] >= G->nodes)
            abort();

         if (G->weight[k] > G->max_weight || G->weight[k] < G->min_weight)
            abort();

         ++count;
      }
      if (G->head[k] != DIJKSTRA_UNUSED)
         abort();

      ++count;
   }
   if (count > G->arcs)
      abort();

   return TRUE;
}



#ifndef NDEBUG
/** Check whether heap is valid
 *
 *  @note Sift up/down does not use the swap, only for the last the changed one is entered.
 */
static
DIJKSTRA_Bool heap_is_valid(
   const unsigned int*        entry,
   const unsigned long long*  value,
   const unsigned int*        order,
   const unsigned int         used,
   const unsigned int         size
   )
{
   unsigned int i;

   if ( entry == NULL || value == NULL || order == NULL || used  >  size )
      return FALSE;

   /* check heap property */
   for (i = 0; i < used / 2; ++i)
   {
      if ( value[entry[i]] > value[entry[i + i]] )
         return FALSE;
      if ( i + i + 1 < used && value[entry[i]] > value[entry[i + i + 1]] )
         return FALSE;
   }

   return TRUE;
}
#endif



/** Moves an entry down in the vector until the sorting is valid again. */
static
void sift_down(
   unsigned int*              entry,
   const unsigned long long*  value,
   unsigned int*              order,
   unsigned int               used,
   unsigned int               current
   )
{
   unsigned int        child = current + current;
   unsigned int        ent   = entry[current];
   unsigned long long  val   = value[ent];
   unsigned int        e;

   while (child < used)
   {
      e = entry[child];

      if (child + 1 < used)
      {
         if (value[entry[child + 1]] < value[e])
         {
            ++child;
            e = entry[child];
         }
      }
      if (value[e] >= val)
         break;

      entry[current] = e;
      order[e] = current;

      current = child;
      child += child;
   }
   entry[current] = ent;
   order[ent] = current;
}


/** Moves an entry up in the vector until the sorting is valid again. */
static
void sift_up(
   unsigned int*              entry,
   const unsigned long long*  value,
   unsigned int*              order,
   unsigned int               current
   )
{
   unsigned int       parent;
   unsigned int       ent = entry[current];
   unsigned long long val = value[ent];
   unsigned int       e;

   while (current > 0)
   {
      parent = current / 2;
      e = entry[parent];

      if (value[e] <= val)
         break;

      entry[current] = e;
      order[e] = current;
      current = parent;
   }
   entry[current] = ent;
   order[ent] = current;
}



/** Dijkstra's algorithm using binary heaps */
unsigned int
graph_dijkstra_bh(
  const Dijkstra_Graph*  G,                  /**< directed graph */
  unsigned int           start,              /**< start node */
  unsigned long long*    dist,               /**< node distances */
  unsigned int*          pred,               /**< node predecessors in final shortest path tree */
  unsigned int*          entry,              /**< temporary storage (for each node - must be allocated by user) */
  unsigned int*          order               /**< temporary storage (for each node - must be allocated by user) */
  )
{
   unsigned int       used;
   unsigned int       head;
   unsigned int       tail;
   unsigned long long weight;
   unsigned int       i;
   unsigned int       e;
   unsigned int       iters = 0;

   assert( Dijsktra_graphIsValid(G) );
   assert( start < G->nodes );
   assert( dist != NULL );
   assert( pred != NULL );
   used = 0;

   assert( heap_is_valid(entry, dist, order, used, G->nodes) );

   /* initialize nodes */
   for (i = 0; i < G->nodes; ++i)
   {
      dist[i] = DIJKSTRA_FARAWAY;
      order[i] = DIJKSTRA_UNUSED;
      pred[i] = DIJKSTRA_UNUSED;
   }

   /* enter start node into heap */
   entry[0]     = start;
   order[start] = 0;
   pred[start]  = DIJKSTRA_UNUSED;
   dist[start]  = 0;

   ++used;

   /* loop while heap is not empty */
   while (used > 0)
   {
      /* get next node */
      --used;

      tail = entry[0];

      entry[0]        = entry[used];
      order[entry[0]] = 0;
      order[tail]     = DIJKSTRA_UNUSED;

      sift_down(entry, dist, order, used, 0);

      assert( heap_is_valid(entry, dist, order, used, G->nodes) );
      assert( entry[used] < G->nodes );

      /* check adjacent nodes */
      for (e = G->outbeg[tail]; G->head[e] != DIJKSTRA_UNUSED; ++e)
      {
         head   = G->head[e];
         weight = G->weight[e] + dist[tail];

	 /* Can we improve the current shortest path? */
         if (dist[head] > weight)
         {
            assert( heap_is_valid(entry, dist, order, used, G->nodes) );
            assert( used < G->nodes );
            assert( head <= G->nodes );

            pred[head] = tail;
            dist[head] = weight;

            if (order[head] == DIJKSTRA_UNUSED)
            {
               assert( head < G->nodes );

               entry[used] = head;
               order[head] = used;

               ++used;

               sift_up(entry, dist, order, used - 1);
            }
            else
            {
               sift_up(entry, dist, order, order[head]);
            }
            assert( heap_is_valid(entry, dist, order, used, G->nodes) );

            ++iters;
         }
      }
   }
   assert( heap_is_valid(entry, dist, order, used, G->nodes) );

   return iters;
}
