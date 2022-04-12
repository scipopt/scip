/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
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


#define STP_UNIFORM_MINRATIO 0.9
#define STP_UNIFORM_RANGEMIN 0.9
#define STP_UNIFORM_RANGEMAX 1.1


/**@name Local methods
 *
 * @{
 */


/**@} */

/**@name Interface methods
 *
 * @{
 */


/** is the given graph a variant that is effectively an STP?? */
SCIP_Bool graph_typeIsSpgLike(
   const GRAPH*          g                   /**< the graph */
   )
{
   const int type = g->stp_type;
   assert(g);

   return (type == STP_SPG || type == STP_RSMT || type == STP_GSTP || type == STP_OARSMT);
}


/** is the given graph (originally) undirected? */
SCIP_Bool graph_typeIsUndirected(
   const GRAPH*          g                   /**< the graph */
   )
{
   assert(g);

   if( g->stp_type == STP_SAP || g->stp_type == STP_DHCSTP || g->stp_type == STP_NWPTSPG )
      return FALSE;

   return TRUE;
}


/** is the given graph (originally) directed? */
SCIP_Bool graph_typeIsDirected(
   const GRAPH*          g                   /**< the graph */
   )
{
   assert(g);

   return !(graph_typeIsUndirected(g));
}


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

   if( !graph_typeIsUndirected(g) )
      printf("e: %d   %d->%d (%d->%d) cost:=%f costrev=%f \n", e, t, h, g->term[t], g->term[h], g->cost[e], g->cost[flipedge(e)]);
   else if( graph_pc_isPcMw(g) )
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

   if( g->stp_type == STP_NWPTSPG && graph_knotIsNWLeaf(g, k) )
      printf("...%d is a leaf-terminal \n", k);
}


/** graph with multi-edges? */
SCIP_Bool graph_hasMultiEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             verbose             /**< be verbose? */
)
{
   const int nnodes = graph_get_nNodes(g);
   int* count;
   SCIP_Bool hasMultiEdges = FALSE;

   assert(scip);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &count, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      count[k] = 0;

   for( int k = 0; k < nnodes && !hasMultiEdges; k++ )
   {
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( count[head] > 0 )
         {
            SCIPdebugMessage("problem for edge %d->%d \n", k, head);

            if( verbose )
            {
               SCIPerrorMessage("parallel edge %d->%d found \n", k + 1, head + 1);
            }

            hasMultiEdges = TRUE;
         }

         count[head]++;
      }

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         count[head]--;
      }
   }

   SCIPfreeBufferArray(scip, &count);

   return hasMultiEdges;
}


/** has the graph almost uniform edge weights? */
SCIP_Bool graph_isAlmostUniform(
   const GRAPH*          g                   /**< the graph */
   )
{
   const int nedges = graph_get_nEdges(g);
   SCIP_Real avg = 0.0;
   SCIP_Real avg_lower;
   SCIP_Real avg_upper;
   int edgecount = 0;
   int nrangeedges = 0;

   for( int e = 0; e < nedges; e+= 2 )
   {
      SCIP_Real edgecost;
      if( g->oeat[e] == EAT_FREE )
         continue;

      edgecost = g->cost[e];

      if( EQ(edgecost, 0.0) || EQ(edgecost, FARAWAY) )
         continue;

      edgecount++;
      avg += edgecost;
   }

   if( edgecount == 0 )
      return TRUE;

   avg /= (SCIP_Real) edgecount;

   assert(STP_UNIFORM_RANGEMIN < 1.0);
   assert(STP_UNIFORM_RANGEMAX > 1.0);

   avg_lower = avg * STP_UNIFORM_RANGEMIN;
   avg_upper = avg * STP_UNIFORM_RANGEMAX;

   for( int e = 0; e < nedges; e += 2 )
   {
      SCIP_Real edgecost;
      if( g->oeat[e] == EAT_FREE )
         continue;

      edgecost = g->cost[e];

      if( EQ(edgecost, 0.0) || EQ(edgecost, FARAWAY) )
         continue;

      if( edgecost < avg_lower || edgecost > avg_upper )
         continue;

      nrangeedges++;
   }

   SCIPdebugMessage("number of unreduced edges=%d \n", edgecount);
   SCIPdebugMessage("number of in-range edges=%d \n", nrangeedges);

   return (((SCIP_Real) nrangeedges / (SCIP_Real) edgecount) > STP_UNIFORM_MINRATIO);
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

   printf("reduced graph stats: ");

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


/** prints edge conflicts */
SCIP_RETCODE graph_printEdgeConflicts(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< the graph */
)
{
   int* childcount;
   int* pseudonodecount;
   int nconflicts;
   int npseudoconflicts;
   int npseudofixed;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const int nedgesorg = MAX(g->orgedges, g->edges);

   assert(scip != NULL && g != NULL);
   assert(g->ancestors != NULL);
   assert(nedgesorg % 2 == 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &childcount, nedgesorg / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pseudonodecount, nnodes) );

   for( int e = 0; e < nedgesorg / 2; e++ )
      childcount[e] = 0;

   for( int e = 0; e < nnodes; e++ )
      pseudonodecount[e] = 0;

   npseudofixed = graph_getNfixpseudonodes(g);
   nconflicts = 0;
   npseudoconflicts = 0;

   for( int e = 0; e < nedges; e += 2 )
   {
      SCIP_Bool hasConflict = FALSE;
      SCIP_Bool hasPseudoConflict = FALSE;

      const int nPseudoAncestors = graph_edge_nPseudoAncestors(g, e);
      const int* pseudoAncestors = graph_edge_getPseudoAncestors(g, e);
      for( IDX* curr = g->ancestors[e]; curr != NULL; curr = curr->parent )
      {
         assert(curr->index >= 0 && curr->index / 2 < nedgesorg / 2);
         if( childcount[curr->index / 2] > 0 )
            hasConflict = TRUE;

         childcount[curr->index / 2]++;
      }

      for( int k = 0; k < nPseudoAncestors; k++ )
      {
         const int a = pseudoAncestors[k];
         assert(a >= 0 && a < nnodes);

         if( pseudonodecount[a] > 0 )
            hasPseudoConflict = TRUE;

         pseudonodecount[a]++;
      }

      if( hasConflict )
      {
         nconflicts++;
         assert(hasPseudoConflict);
      }

      if( hasPseudoConflict )
         npseudoconflicts++;

   }

   if( graph_pc_isPcMw(g) )
   {
      assert(g->extended);

      for( int e = 0; e < nnodes; e++ )
      {
         if( !Is_term(g->term[e]) )
         {
            SCIP_Bool hasConflict = FALSE;
            SCIP_Bool hasPseudoConflict = FALSE;

            const int nPseudoAncestors = graph_knot_nPseudoAncestors(g, e);
            const int* pseudoAncestors = graph_knot_getPseudoAncestors(g, e);
            for( IDX* curr = g->pcancestors[e]; curr != NULL; curr = curr->parent )
            {
               assert(curr->index >= 0 && curr->index / 2 < nedgesorg / 2);
               if( childcount[curr->index / 2] > 0 )
                  hasConflict = TRUE;

               childcount[curr->index / 2]++;
            }

            for( int k = 0; k < nPseudoAncestors; k++ )
            {
               const int a = pseudoAncestors[k];
               assert(a >= 0 && a < nnodes);

               if( pseudonodecount[a] > 0 )
                  hasPseudoConflict = TRUE;

               pseudonodecount[a]++;
            }

            if( hasConflict )
            {
               nconflicts++;
               assert(hasPseudoConflict);
            }

            if( hasPseudoConflict )
               npseudoconflicts++;
         }
      }
   }

   printf("number of edge conflicts %d \n", nconflicts);
   printf("number of pseudo-ancestor conflicts %d \n", npseudoconflicts);
   printf("number of fixed pseudo-ancestors %d \n", npseudofixed);

   SCIPfreeBufferArray(scip, &pseudonodecount);
   SCIPfreeBufferArray(scip, &childcount);

   return SCIP_OKAY;
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
      }

      if( nterms && Is_term(graph->term[k]) )
         t++;
   }

   if( nnodes )
      *nnodes = MAX(v, 1);

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

/**@} */
