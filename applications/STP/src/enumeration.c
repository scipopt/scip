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

/**@file   enumeration.c
 * @brief  includes enumeration algorithms for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * Methods for enumerating solutions (i.e. trees) to Steiner tree problems.
 *
 * A list of all interface methods can be found in enumeration.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */
//#define SCIP_DEBUG

#include <assert.h>
#include "portab.h"
#include "solstp.h"
#include "enumeration.h"




/** Finds maximum two terminals. Terminal for rooted problems is always the maximum. */
static
void pcmwFindMax2Terms(
   const GRAPH*          g,                  /**< graph data structure */
   int*                  term_max,
   int*                  term_max2
   )
{
   const int pseudoroot = graph_pc_isRootedPcMw(g) ? -1 : g->source;
   const int nnodes = graph_get_nNodes(g);

   SCIP_Real max = -1.0;
   SCIP_Real max2 = -1.0;
   *term_max = UNKNOWN;
   *term_max2 = UNKNOWN;

   for( int i = 0; i < nnodes; ++i )
   {
      if( !Is_term(g->term[i]) || i == pseudoroot )
         continue;

      if( g->prize[i] > max || i == g->source )
      {
         max2 = max;
          *term_max2 = *term_max;

         max = g->prize[i];
         *term_max = i;
      }
      else if( g->prize[i] > max2 )
      {
         max2 = g->prize[i];
         *term_max2 = i;
      }
   }
}

/** builds path and connect if possible */
static
SCIP_RETCODE tryPathPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   term_start,
   int                   term_end,
   STP_Bool* RESTRICT    nodes_inPath
)
{
   PATH* path;
   const int nnodes = graph_get_nNodes(g);
   assert(term_start != term_end);
   assert(term_end != g->source);

   SCIP_CALL(SCIPallocBufferArray(scip, &path, nnodes));
   BMSclearMemoryArray(nodes_inPath, nnodes);
   nodes_inPath[term_start] = TRUE;

   graph_path_exec(scip, g, FSP_MODE, term_start, g->cost, path);

   /* profitable to connect 2nd terminal? */
   if( LT(path[term_end].dist, g->prize[term_end]) )
   {
      for( int node = term_end; node != term_start; node = g->tail[path[node].edge] )
      {
         nodes_inPath[node] = TRUE;
         assert(graph_edge_isInRange(g, path[node].edge));
      }
   }

   SCIPfreeBufferArray(scip, &path);

   return SCIP_OKAY;
}


/** enumeration of rooted PC/MW */
static
SCIP_RETCODE findSolRPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int* RESTRICT         result              /**< solution array, indicating whether an edge is in the solution */
)
{
   int term_max;
   int term_min;
   const int nnodes = graph_get_nNodes(g);
   STP_Bool* RESTRICT connected;

   assert(graph_pc_isRootedPcMw(g));
   assert(graph_isMarked(g));

   assert(g->terms <= 2);

   if( g->terms == 1 )
      return SCIP_OKAY;

   pcmwFindMax2Terms(g, &term_max, &term_min);
   assert(term_max != UNKNOWN && term_min != UNKNOWN);
   assert(term_max == g->source);
   assert(GE(g->prize[term_max], g->prize[term_min]));

   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );

   SCIP_CALL( tryPathPcMw(scip, g, term_max, term_min, connected) );

   assert(connected[term_max]);
   assert(connected[term_min] || !graph_pc_knotIsFixedTerm(g, term_min));

   graph_pc_2trans(scip, g);
   solstp_pruneFromNodes(scip, g, result, connected);
   graph_pc_2org(scip, g);

   SCIPfreeBufferArray(scip, &connected);

   return SCIP_OKAY;
}


/** enumeration of unrooted PC/MW with one terminal other than the root  */
static
void findSolPcMw1Term(
   const GRAPH*          g,                  /**< graph data structure */
   int* RESTRICT         result              /**< solution array, indicating whether an edge is in the solution */
)
{
   int a;
   const int root = g->source;

   assert(g->terms == 2);

   for( a = g->outbeg[root]; a != EAT_LAST; a = g->oeat[a] )
   {
      if( !graph_pc_knotIsDummyTerm(g, g->head[a]) )
      {
         result[a] = CONNECT;
         break;
      }
   }

   assert(a != EAT_LAST);
}


/** enumeration of unrooted PC/MW with one terminal other than the root  */
static
SCIP_RETCODE findSolPcMw2Term(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int* RESTRICT         result              /**< solution array, indicating whether an edge is in the solution */
)
{
   int term_max;
   int term_min;
   const int nnodes = graph_get_nNodes(g);
   STP_Bool* RESTRICT connected;

   assert(g->terms == 3);
   assert(!g->extended);

   pcmwFindMax2Terms(g, &term_max, &term_min);
   assert(term_max != UNKNOWN && term_min != UNKNOWN);
   assert(term_max != g->source && term_min != g->source);
   assert(GE(g->prize[term_max], g->prize[term_min]));

   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );

   SCIP_CALL( tryPathPcMw(scip, g, term_max, term_min, connected) );

   graph_pc_2trans(scip, g);
   solstp_pruneFromNodes(scip, g, result, connected);
   graph_pc_2org(scip, g);

   SCIPfreeBufferArray(scip, &connected);

   return SCIP_OKAY;
}


/** enumeration of unrooted PC/MW */
static
SCIP_RETCODE findSolPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int* RESTRICT         result              /**< solution array, indicating whether an edge is in the solution */
)
{
   assert(!graph_pc_isRootedPcMw(g));
   assert(g->terms <= 3);
   assert(g->terms > 1);
   assert(graph_isMarked(g));

   if( g->terms == 2 )
   {
      SCIPdebugMessage("finding solution for 1 terminal \n");

      findSolPcMw1Term(g, result);
   }
   else
   {
      SCIPdebugMessage("finding solution for 2 terminals \n");

      SCIP_CALL( findSolPcMw2Term(scip, g, result) );
   }

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


/** enumeration of PC/MW */
SCIP_RETCODE enumeration_findSolPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int* RESTRICT         result              /**< solution array, indicating whether an edge is in the solution */
)
{
   const int nedges = graph_get_nEdges(g);

   assert(scip && g && result);
   assert(graph_pc_isPcMw(g));
   assert(enumeration_isPossible(g));

   for( int e = 0; e < nedges; e++ )
      result[e] = UNKNOWN;

   graph_mark(g);

   if( graph_pc_isRootedPcMw(g) )
   {
      SCIP_CALL( findSolRPcMw(scip, g, result) );
   }
   else
   {
      SCIP_CALL( findSolPcMw(scip, g, result) );
   }

   return TRUE;
}


/** enumeration possible? */
SCIP_Bool enumeration_isPossible(
   const GRAPH*          g                  /**< graph data structure */
)
{
   assert(g);

   if( graph_pc_isRootedPcMw(g) )
   {
      return (g->terms <= 2);
   }
   else if( graph_pc_isPcMw(g) )
   {
      return (g->terms <= 3);
   }

   return FALSE;

}
