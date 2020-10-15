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

/**@file   graph_stats.c
 * @brief  includes several graph statistic methods for Steiner problem graphs
 * @author Daniel Rehfeldt
 *
 * This file contains methods to obtain statistics on the given Steiner problem instance.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h)                             */
//#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "portab.h"
#include "graph.h"



/** is the edge blocked? */
SCIP_Bool graph_edge_isBlocked(
   const GRAPH*          g,                  /**< the graph */
   int                   e                   /**< the edge */
   )
{
   assert(g);
   assert(e >= 0 && e < g->edges);

   if( EQ(g->cost[e], BLOCKED_MINOR) || EQ(g->cost[e], BLOCKED) )
      return TRUE;

   return FALSE;
}


/** is edge in range? */
SCIP_Bool graph_edge_isInRange(
   const GRAPH*          g,                  /**< the graph */
   int                   e                   /**< the edge */
   )
{
   assert(g);

   return (0 <= e && e < g->edges);
}

/** has edge been deleted? */
SCIP_Bool graph_edge_isDeleted(
   const GRAPH*          g,                  /**< the graph */
   int                   e                   /**< the edge */
   )
{
   assert(g);

   return (g->oeat[e] == EAT_FREE);
}


/** print information on edge */
void graph_edge_printInfo(
   const GRAPH*          g,                  /**< the graph */
   int                   e                   /**< the edge */
   )
{
   const int t = g->tail[e];
   const int h = g->head[e];

   if( graph_pc_isPcMw(g) )
      printf("e: %d   %d->%d (%d->%d) cost:=%f costrev=%f \n", e, t, h, g->term[t], g->term[h], g->cost[e], g->cost[flipedge(e)]);
   else
      printf("e: %d   %d->%d (%d->%d) cost:=%f \n", e, t, h, g->term[t], g->term[h], g->cost[e]);
}


/** print information on node */
void graph_knot_printInfo(
   const GRAPH*          g,                  /**< the graph */
   int                   k                   /**< the node */
   )
{
   assert(!graph_pc_isPcMw(g) || g->term2edge != NULL);

   if( graph_pc_isPcMw(g) )
   {
      assert(g->prize);

      if( graph_pc_knotIsNonLeafTerm(g, k) )
         printf("node %d: term=%d grad=%d prize=%f (non-leaf terminal) \n", k, g->term[k], g->grad[k], g->prize[k]);
      else if( graph_pc_knotIsFixedTerm(g, k) )
         printf("node %d: term=%d grad=%d prize=%f (fixed terminal) \n", k, g->term[k], g->grad[k], g->prize[k]);
      else if( Is_term(g->term[k]) )
         printf("node %d: term=%d grad=%d prize=%f (standard terminal) \n", k, g->term[k], g->grad[k], g->prize[k]);
      else
         printf("node %d: term=%d grad=%d prize=%f \n", k, g->term[k], g->grad[k], g->prize[k]);
   }
   else
   {
      printf("node %d: term=%d grad=%d  \n", k, g->term[k], g->grad[k]);
   }

   if( k == g->source )
      printf("...%d is the root! \n", k);

}


/** print information on graph  */
void graph_printInfo(
   const GRAPH*          g                   /**< the graph */
   )
{
   char type[64];

   switch( g->stp_type )
      {
   case 0:
      strcpy(type, "STP_SPG");
      break;
   case 1:
      strcpy(type, "STP_SAP");
      break;
   case 2:
      strcpy(type, "STP_PCSPG");
      break;
   case 3:
      strcpy(type, "STP_RPCSPG");
      break;
   case 4:
      strcpy(type, "STP_NWSPG");
      break;
   case 5:
      strcpy(type, "STP_DCSTP");
      break;
   case 6:
      strcpy(type, "STP_NWPTSPG");
      break;
   case 7:
      strcpy(type, "STP_RSMT");
      break;
   case 8:
      strcpy(type, "STP_OARSMT");
      break;
   case 9:
      strcpy(type, "STP_MWCSP");
      break;
   case 10:
      strcpy(type, "STP_DHCSTP");
      break;
   case 11:
      strcpy(type, "STP_GSTP");
      break;
   case 12:
      strcpy(type, "STP_RMWCSP");
      break;
   case 13:
      strcpy(type, "STP_BRMWCSP");
      break;
   default:
      strcpy(type, "UNKNOWN");
   }

   assert(g != NULL);
   if( graph_pc_isPcMw(g) )
   {
      printf("nodes=%d, edges=%d, terminals=%d, root=%d, type=%s, isExtended=%d\n", g->knots, g->edges, g->terms, g->source, type, g->extended);
      if( graph_pc_isPc(g) )
      {
         printf("non-leaf terminals=%d, ", graph_pc_nNonLeafTerms(g));
         printf("fixed terminals=%d, ", graph_pc_nFixedTerms(g));
         printf("proper terminals=%d \n", graph_pc_nProperPotentialTerms(g));
      }
   }
   else
   {
      printf("nodes=%d, edges=%d, terminals=%d, root=%d, type=%s \n", g->knots, g->edges, g->terms, g->source, type);
   }

   if( g->stp_type == STP_BRMWCSP )
      printf("budget=%f \n", g->budget);
}


/** print information on graph that has been subject to reductions  */
void graph_printInfoReduced(
   const GRAPH*          g                   /**< the graph */
)
{
   int nnodes;
   int nedges;
   int nterms;

   graph_get_nVET(g, &nnodes, &nedges, &nterms);

   if( graph_pc_isPcMw(g) )
   {
      printf("nodes=%d, edges=%d, terminals=%d, root=%d, type=%d, isExtended=%d \n",
            nnodes, nedges, nterms, g->source, g->stp_type, g->extended);

      if( graph_pc_isPc(g) )
      {
         printf("non-leaf terminals=%d, ", graph_pc_nNonLeafTerms(g));
         printf("fixed terminals=%d, ", graph_pc_nFixedTerms(g));
         printf("proper terminals=%d \n", graph_pc_nProperPotentialTerms(g));
      }
   }
   else
   {
      printf("nodes=%d, edges=%d, terminals=%d, root=%d, type=%d \n", nnodes,
            nedges, nterms, g->source, g->stp_type);
   }
}


/** get number of nodes */
int graph_get_nNodes(
   const GRAPH*    graph               /**< the graph */
)
{
   assert(graph);
   assert(graph->knots >= 1);

   return graph->knots;
}


/** get number of edges */
int graph_get_nEdges(
   const GRAPH*    graph               /**< the graph */
)
{
   assert(graph);
   assert(graph->edges >= 0);

   return graph->edges;
}


/** get number of terminals */
int graph_get_nTerms(
   const GRAPH*    graph               /**< the graph */
)
{
   assert(graph);

   return graph->terms;
}


/** get (real) number of nodes, edges, terminals */
void graph_get_nVET(
   const GRAPH*    graph,              /**< the graph */
   int*            nnodes,             /**< number of nodes */
   int*            nedges,             /**< number of edges */
   int*            nterms              /**< number of terminals */
   )
{
   int v = 0;
   int e = 0;
   int t = 0;
   int vorg;

   assert(graph != NULL);

   vorg = graph->knots;

   for( int k = 0; k < vorg; k++ )
   {
      if( graph->grad[k] > 0 )
      {
         v++;
         e += graph->grad[k];
         if( Is_term(graph->term[k]) )
            t++;
      }
   }

   if( nnodes )
      *nnodes = v;

   if( nedges )
      *nedges = e;

   if( nterms )
      *nterms = t;

   return;
}


/** get (real) average degree of graph */
SCIP_Real graph_get_avgDeg(
   const GRAPH*    graph               /**< the graph */
   )
{
   int v = 0;
   int e = 0;
   const int vorg = graph->knots;

   assert(graph);

   for( int k = 0; k < vorg; k++ )
   {
      if( graph->grad[k] > 0 )
      {
         v++;
         e += graph->grad[k];
      }
   }

   if( v == 0 )
      return 0.0;

   return ((double) e / v );
}
