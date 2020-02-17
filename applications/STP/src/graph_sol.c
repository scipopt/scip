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

/**@file   graph_sol.c
 * @brief  includes methods working on solutions (i.e. trees) to Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * Methods for manipulating solutions (i.e. trees) to Steiner tree problems, such as pruning.
 * Also includes methods for obtaining information about solutions.
 *
 * A list of all interface methods can be found in graph.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */


#include "graph.h"
#include "probdata_stp.h"
#include "portab.h"


/** prune a Steiner tree in such a way that all leaves are terminals */
static
SCIP_RETCODE pruneSteinerTreeStp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int*                  result,             /**< ST edges, which need to be set to UNKNOWN */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   PATH* mst;
   int count;
   const int nnodes = graph_get_nNodes(g);
#ifndef NEDBUG
   int nconnected = 0;
#endif

   assert(scip != NULL);
   assert(cost != NULL);
   assert(result != NULL);
   assert(connected != NULL);

#ifndef NEDBUG
   for( int i = 0; i < g->edges; i++ )
      assert(UNKNOWN == result[i]);

   for( int i = nnodes - 1; i >= 0; --i )
      if( connected[i] )
         nconnected++;

   assert(nconnected >= g->terms);
   assert(g->source >= 0);
   assert(g->source < nnodes);
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );

   /* compute the MST */
   for( int i = nnodes - 1; i >= 0; --i )
      g->mark[i] = connected[i];

   graph_path_exec(scip, g, MST_MODE, g->source, cost, mst);

   for( int i = nnodes - 1; i >= 0; --i )
   {
      if( connected[i] && (mst[i].edge != -1) )
      {
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == UNKNOWN);

         result[mst[i].edge] = 0;
      }
   }

   /* prune */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for( int i = nnodes - 1; i >= 0; --i )
      {
         int j;

         if( !g->mark[i] )
            continue;

         if( g->term[i] == 0 )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == 0 )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == 0 )
               {
                  result[j]    = -1;
                  g->mark[i]   = FALSE;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
         }
      }
   }
   while( count > 0 );

   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}


/* prune the (rooted) prize collecting Steiner tree in such a way that all leaves are terminals */
static
SCIP_RETCODE pruneSteinerTreePc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int*                  result,             /**< ST edges (need to be set to UNKNOWN) */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   PATH* mst;
   int count;
   int root = g->source;
   const int nnodes = g->knots;
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(g);

   assert(g != NULL && cost != NULL && result != NULL && connected != NULL);
   assert(g->extended);

#ifndef NEDBUG
   for( int i = 0; i < g->edges; i++ )
      assert(UNKNOWN == result[i]);
#endif

   if( rpcmw )
   {
      for( int i = 0; i < nnodes; i++ )
      {
         if( connected[i] && !graph_pc_knotIsDummyTerm(g, i) )
            g->mark[i] = TRUE;
         else
            g->mark[i] = FALSE;

         assert(g->mark[i] || !graph_pc_knotIsFixedTerm(g, i));
      }

      if( !g->mark[root] )
      {
         printf("FAIL in SCIPStpHeurTMPrunePc, root not connected \n");
         return SCIP_ERROR;
      }
   }
   else
   {
      int proot;
      for( int i = 0; i < nnodes; i++ )
      {
         if( connected[i] && !Is_term(g->term[i]) )
            g->mark[i] = TRUE;
         else
            g->mark[i] = FALSE;
      }

      proot = -1;
      if( SCIPprobdataGetNTerms(scip) == g->terms && SCIPprobdataGetNNodes(scip) == nnodes )
      {
         int min = nnodes;
         const int* termsorder = SCIPprobdataGetPctermsorder(scip);

         for( int k = 0; k < nnodes; k++ )
         {
            if( termsorder[k] < min && connected[k] )
            {
               assert(Is_pseudoTerm(g->term[k]));

               min = termsorder[k];
               proot = k;
            }
         }

         assert(min >= 0);
         assert(proot == -1 || min < nnodes);
      }
      else
      {
         for( int a = g->outbeg[root]; a != EAT_LAST; a = g->oeat[a] )
         {
            const int head = g->head[a];
            if( !Is_term(g->term[head]) && connected[head] )
            {
               proot = head;
               break;
            }
         }
      }

      /* trivial solution? */
      if( proot == -1 )
      {
         printf("trivial solution in pruning \n");
         for( int a = g->outbeg[g->source]; a != EAT_LAST; a = g->oeat[a] )
         {
            const int head = g->head[a];
            if( Is_term(g->term[head]) )
            {
               assert(connected[head]);
               result[a] = CONNECT;
            }
         }
         return SCIP_OKAY;
      }

      assert(g->mark[proot]);
      root = proot;
   }
   assert(root >= 0);
   assert(root < nnodes);

   SCIPdebugMessage("(non-artificial) root=%d \n", root);

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );
   graph_path_exec(scip, g, MST_MODE, root, cost, mst);

   for( int i = 0; i < nnodes; i++ )
   {
      if( g->mark[i] && (mst[i].edge != UNKNOWN) )
      {
         assert(g->path_state[i] == CONNECT);  assert(g->head[mst[i].edge] == i);  assert(result[mst[i].edge] == -1);
         result[mst[i].edge] = CONNECT;
      }
   }

   /* connect all terminals */
   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && i != g->source )
      {
         int e1;
         int e2;

         if( rpcmw && g->mark[i] )
         {
            assert(g->prize[i] == FARAWAY && connected[i]);
            continue;
         }

         assert(!graph_pc_knotIsFixedTerm(g, i));
         connected[i] = TRUE;

         e1 = g->inpbeg[i];
         assert(e1 >= 0);
         e2 = g->ieat[e1];

         if( e2 == EAT_LAST )
         {
            result[e1] = CONNECT;
         }
         else
         {
            const int k1 = g->tail[e1];
            const int k2 = g->tail[e2];

            assert(e2 >= 0);
            assert(g->ieat[e2] == EAT_LAST);
            assert(k1 == g->source || k2 == g->source);

            if( k1 != g->source && g->path_state[k1] == CONNECT )
               result[e1] = CONNECT;
            else if( k2 != g->source && g->path_state[k2] == CONNECT )
               result[e2] = CONNECT;
            else if( k1 == g->source )
               result[e1] = CONNECT;
            else if( k2 == g->source )
               result[e2] = CONNECT;
         }
      }
      else if( i == root && !rpcmw )
      {
         int e1;
         for( e1 = g->inpbeg[i]; e1 != EAT_LAST; e1 = g->ieat[e1] )
            if( g->tail[e1] == g->source )
               break;
         assert(e1 != EAT_LAST);
         result[e1] = CONNECT;
      }
   }

   /* prune */
   do
   {
      count = 0;

      for( int i = nnodes - 1; i >= 0; --i )
      {
         int j;
         if( !g->mark[i] || g->path_state[i] != CONNECT || Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == CONNECT )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            assert(!Is_term(g->term[i]) && !Is_pseudoTerm(g->term[i]));

            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == CONNECT )
               {
                  result[j]    = -1;
                  g->mark[i]   = FALSE;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
            assert(j != EAT_LAST);
         }
      }
   }
   while( count > 0 );

#ifndef NDEBUG
   /* make sure there is no unconnected vertex */
   for( int i = 0; i < nnodes; i++ )
   {
      if( connected[i] && i != g->source )
      {
         int j;
         for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            if( result[j] == CONNECT )
               break;

         assert(j != EAT_LAST);
      }
   }
#endif

   assert(graph_solIsValid(scip, g, result));
   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}

/*
 * Interface methods
 */

/** Prune solution given by included nodes.
 *  NOTE: For PC/RPC this method will get the original edge costs before pruning! */
SCIP_RETCODE graph_solPrune(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   const int nedges = graph_get_nEdges(g);

   assert(scip && result && connected);
   assert(g->stp_type != STP_DHCSTP);

   for( int e = 0; e < nedges; e++ )
      result[e] = UNKNOWN;

   if( graph_pc_isPcMw(g) )
   {
      SCIP_Real* edgecosts = NULL;
      assert(g->extended);

      /* do we have biased edge costs? */
      if( graph_pc_isPc(g) )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &edgecosts, nedges) );

         graph_pc_getOrgCosts(scip, g, edgecosts);
      }
      else
      {
         edgecosts = g->cost;
      }

      SCIP_CALL( pruneSteinerTreePc(scip, g, edgecosts, result, connected) );

      if( graph_pc_isPc(g) )
         SCIPfreeBufferArray(scip, &edgecosts);
   }
   else
   {
      SCIP_CALL( pruneSteinerTreeStp(scip, g, g->cost, result, connected) );
   }

   return SCIP_OKAY;
}


/** prune solution given by included nodes */
SCIP_RETCODE graph_solPruneFromNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   assert(scip && g && result && connected);
   assert(g->stp_type != STP_DHCSTP);

   SCIP_CALL( graph_solPrune(scip, g, result, connected) );

   return SCIP_OKAY;
}


/** prune solution given by included edges */
SCIP_RETCODE graph_solPruneFromEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int*                  result              /**< ST edges */
   )
{
   STP_Bool* connected;
   const int nedges = g->edges;
   const int nnodes = g->knots;

   assert(scip && g && result);

   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      connected[k] = FALSE;

   for( int e = 0; e < nedges; e++ )
     if( result[e] == CONNECT )
     {
        connected[g->tail[e]] = TRUE;
        connected[g->head[e]] = TRUE;
     }

   SCIP_CALL( graph_solPruneFromNodes(scip, g, result, connected) );

   SCIPfreeBufferArray(scip, &connected);

   return SCIP_OKAY;
}


/** Prunes solution with respect to the provided edges costs. */
SCIP_RETCODE graph_solPruneOnGivenCosts(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs for DHCSTP and PC */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   const int nedges = graph_get_nEdges(g);

   assert(scip && cost && result && connected);

   if( g->stp_type != STP_DHCSTP )
   {
      for( int e = 0; e < nedges; e++ )
         result[e] = UNKNOWN;
   }

   if( graph_pc_isPcMw(g) )
   {
      if( graph_pc_isPc(g) )
      {
         assert(cost);
         SCIP_CALL( pruneSteinerTreePc(scip, g, cost, result, connected) );
      }
      else
      {
         assert(!cost);
         SCIP_CALL( pruneSteinerTreePc(scip, g, g->cost, result, connected) );
      }
   }
   else
      SCIP_CALL( pruneSteinerTreeStp(scip, g, (g->stp_type != STP_DHCSTP) ? g->cost : cost, result, connected) );

   return SCIP_OKAY;
}


/** changes solution according to given root */
SCIP_RETCODE graph_solReroot(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int*                  result,             /**< solution array (CONNECT/UNKNOWN) */
   int                   newroot             /**< the new root */
   )
{
   int* queue;
   int* gmark;
   int size;
   const int nnodes = graph_get_nNodes(g);

   assert(scip != NULL);
   assert(g != NULL);
   assert(result != NULL);
   assert(newroot >= 0 && newroot < nnodes);
   assert(Is_term(g->term[newroot]));

   if( g->grad[newroot] == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &gmark, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &queue, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      gmark[k] = FALSE;

   gmark[newroot] = TRUE;
   size = 0;
   queue[size++] = newroot;

   /* BFS loop */
   while( size )
   {
      const int node = queue[--size];

      /* traverse outgoing arcs */
      for( int a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
      {
         const int head = g->head[a];

         if( !gmark[head] && (result[a] == CONNECT || result[flipedge(a)] == CONNECT ) )
         {
            if( result[flipedge(a)] == CONNECT  )
            {
               result[a] = CONNECT;
               result[flipedge(a)] = UNKNOWN;
            }
            gmark[head] = TRUE;
            queue[size++] = head;
         }
      }
   }

   /* adjust solution if infeasible */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !gmark[k] )
      {
         for( int a = g->outbeg[k]; a != EAT_LAST; a = g->oeat[a] )
         {
            result[a] = UNKNOWN;
            result[flipedge(a)] = UNKNOWN;
         }

         /* not yet connected terminal? */
         if( Is_term(g->term[k]) )
         {
            int a;
            assert(g->stp_type != STP_SPG);

            for( a = g->inpbeg[k]; a != EAT_LAST; a = g->ieat[a] )
            {
               const int node = g->tail[a];
               if( gmark[node] && node != newroot )
               {
                  result[a] = CONNECT;
                  break;
               }
            }
            if( a == EAT_LAST )
            {
               for( a = g->inpbeg[k]; a != EAT_LAST; a = g->ieat[a] )
               {
                  const int node = g->tail[a];
                  if( node == newroot )
                  {
                     result[a] = CONNECT;
                     break;
                  }
               }
            }
            else
               gmark[k] = TRUE;
         }
      }
   }

   SCIPfreeBufferArray(scip, &queue);
   SCIPfreeBufferArray(scip, &gmark);

#ifndef NDEBUG
   {
      const int realroot = g->source;
      g->source = newroot;
      assert(graph_solIsValid(scip, g, result));
      g->source = realroot;
   }
#endif

   return SCIP_OKAY;
}


/** checks whether edge(s) of given primal solution have been deleted */
SCIP_Bool graph_solIsUnreduced(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const int*            result              /**< solution array, indicating whether an edge is in the solution */
   )
{
   const int nedges = graph_get_nEdges(graph);

   assert(scip != NULL);
   assert(result != NULL);

   for( int i = 0; i < nedges; i++ )
   {
      if( result[i] == CONNECT && graph->oeat[i] == EAT_FREE )
         return FALSE;
   }

   return TRUE;
}

/** verifies whether a given primal solution is feasible */
SCIP_Bool graph_solIsValid(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const int*            result              /**< solution array, indicating whether an edge is in the solution */
   )
{
   int* queue = NULL;
   STP_Bool* reached = NULL;
   int size;
   int nterms;
   int termcount;
   const int nnodes = graph_get_nNodes(graph);
   const int root = graph->source;
   SCIP_Bool countpseudo;

   assert(scip && result);
   assert(root >= 0);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &reached, nnodes) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &queue, nnodes) );

   if( graph_pc_isPcMw(graph) && !graph->extended )
   {
      countpseudo = TRUE;
      nterms = graph_pc_nProperPotentialTerms(graph);

      if( !graph_pc_isRootedPcMw(graph) )
         nterms++;
   }
   else
   {
      countpseudo = FALSE;
      nterms = graph->terms;
   }

   for( int i = 0; i < nnodes; i++ )
      reached[i] = FALSE;

   /* BFS until all terminals are reached */

   termcount = 1;
   size = 0;
   reached[root] = TRUE;
   queue[size++] = root;

   while( size )
   {
      const int node = queue[--size];

      for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( result[e] == CONNECT )
         {
            const int i = graph->head[e];

            /* cycle? */
            if( reached[i] )
            {
               SCIPfreeBufferArray(scip, &queue);
               SCIPfreeBufferArray(scip, &reached);

               SCIPdebugMessage("solution contains a cycle ... \n");
               return FALSE;
            }

            if( countpseudo )
            {
               if( Is_pseudoTerm(graph->term[i]) || graph_pc_knotIsFixedTerm(graph, i) )
                  termcount++;
            }
            else
            {
               if( Is_term(graph->term[i]) )
                  termcount++;
            }

            reached[i] = TRUE;
            queue[size++] = i;
         }
      }
   }

#ifdef SCIP_DEBUG
   if( termcount != nterms )
   {
      printf("termcount %d graph->terms %d \n", termcount, nterms);
      printf("root %d \n", root);

      for( int i = 0; i < nnodes; i++ )
      {
         const int isMandatoryTerm = countpseudo?
               (Is_pseudoTerm(graph->term[i]) || graph_pc_knotIsFixedTerm(graph, i)) : Is_term(graph->term[i]);

         if( !reached[i] && isMandatoryTerm )
         {
            if( graph_pc_isPc(graph) && graph_pc_termIsNonLeafTerm(graph, i) )
               continue;

            printf("fail: ");
            graph_knot_printInfo(graph, i);

            for( int e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
            {
               printf("...neighbor: ");
               graph_knot_printInfo(graph, graph->tail[e]);
            }
         }
      }
   }
#endif

   SCIPfreeBufferArray(scip, &queue);
   SCIPfreeBufferArray(scip, &reached);

   return (termcount == nterms);
}

/** mark endpoints of edges in given list */
void graph_solSetNodeList(
   const GRAPH*          g,              /**< graph data structure */
   STP_Bool*             solnode,        /**< solution nodes array (TRUE/FALSE) */
   IDX*                  listnode        /**< edge list */
   )
{
   int i;
   IDX* curr;

   assert(g != NULL);
   assert(solnode != NULL);

   curr = listnode;

   while( curr != NULL )
   {
      i = curr->index;

      solnode[g->head[i]] = TRUE;
      solnode[g->tail[i]] = TRUE;

      curr = curr->parent;
   }
}

/** compute solution value for given edge-solution array (CONNECT/UNKNOWN) and offset */
SCIP_Real graph_solGetObj(
   const GRAPH*          g,                  /**< the graph */
   const int*            soledge,            /**< solution */
   SCIP_Real             offset,             /**< offset */
   int                   nedges              /**< number of edges todo delete */
   )
{
   SCIP_Real obj = offset;
   const SCIP_Real* const edgecost = g->cost;

   assert(nedges == g->edges);
   assert(!graph_pc_isPcMw(g) || g->extended);

   for( int e = 0; e < nedges; e++ )
      if( soledge[e] == CONNECT )
         obj += edgecost[e];

   return obj;
}


/** computes number of edges in solution value */
int graph_solGetNedges(
   const GRAPH*          g,                  /**< the graph */
   const int*            soledge             /**< solution */
   )
{
   const int nedges = graph_get_nEdges(g);
   int edgecount = 0;

   assert(soledge);

   for( int e = 0; e < nedges; e++ )
      if( soledge[e] == CONNECT )
         edgecount++;

   return edgecount;
}


/** marks vertices for given edge-solution array (CONNECT/UNKNOWN) */
void graph_solSetVertexFromEdge(
   const GRAPH*          g,                  /**< the graph */
   const int*            result,             /**< solution array (CONNECT/UNKNOWN) */
   STP_Bool*             solnode             /**< marks whether node is in solution */
)
{
   const int nedges = g->edges;
   const int nnodes = g->knots;

   assert(g && result && solnode);

   for( int i = 0; i < nnodes; i++ )
      solnode[i] = FALSE;

   solnode[g->source] = TRUE;

   for( int e = 0; e < nedges; e++ )
   {
      if( result[e] == CONNECT )
      {
         assert(g->oeat[e] != EAT_FREE);

         solnode[g->head[e]] = TRUE;
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < nedges; e++ )
      if( result[e] == CONNECT )
         assert(solnode[g->head[e]] && solnode[g->tail[e]]);
#endif
}

/** get original solution */
SCIP_RETCODE graph_solGetOrg(
   SCIP*           scip,               /**< SCIP data structure */
   const GRAPH*    transgraph,         /**< the transformed graph */
   const GRAPH*    orggraph,           /**< the original graph */
   const int*      transsoledge,       /**< solution for transformed problem */
   int*            orgsoledge          /**< new retransformed solution */
)
{
   STP_Bool* orgnodearr;
   IDX** const ancestors = transgraph->ancestors;

   const int transnedges = transgraph->edges;
   const int orgnnodes = orggraph->knots;
   const SCIP_Bool pcmw = graph_pc_isPcMw(transgraph);

   assert(transgraph != NULL && orggraph != NULL && transsoledge != NULL && orgsoledge != NULL);
   assert(transgraph->ancestors != NULL);
   assert(transgraph->stp_type == orggraph->stp_type);

   SCIP_CALL( SCIPallocBufferArray(scip, &orgnodearr, orgnnodes) );

   for( int k = 0; k < orgnnodes; k++ )
      orgnodearr[k] = FALSE;

   for( int e = 0; e < transnedges; e++ )
      if( transsoledge[e] == CONNECT )
         graph_solSetNodeList(orggraph, orgnodearr, ancestors[e]);

   /* retransform edges fixed during graph reduction */
   graph_solSetNodeList(orggraph, orgnodearr, graph_get_fixedges(transgraph));

   if( pcmw )
   {
      // potentially single-vertex solution?
      if( graph_pc_isRootedPcMw(transgraph) && transgraph->terms == 1 && graph_pc_nFixedTerms(orggraph) == 1 )
         orgnodearr[orggraph->source] = TRUE;

      SCIP_CALL( graph_solMarkPcancestors(scip, transgraph->pcancestors, orggraph->tail, orggraph->head, orgnnodes,
            orgnodearr, NULL, NULL, NULL, NULL ) );
   }

   /* prune solution (in original graph) */
   SCIP_CALL( graph_solPrune(scip, orggraph, orgsoledge, orgnodearr) );

   SCIPfreeBufferArray(scip, &orgnodearr);

   assert(graph_solIsValid(scip, orggraph, orgsoledge));

   return SCIP_OKAY;
}



/** mark original solution */
SCIP_RETCODE graph_solMarkPcancestors(
   SCIP*           scip,               /**< SCIP data structure */
   IDX**           pcancestors,        /**< the ancestors */
   const int*      tails,              /**< tails array */
   const int*      heads,              /**< heads array */
   int             orgnnodes,          /**< original number of nodes */
   STP_Bool*       solnodemark,        /**< solution nodes mark array */
   STP_Bool*       soledgemark,        /**< solution edges mark array or NULL */
   int*            solnodequeue,       /**< solution nodes queue or NULL  */
   int*            nsolnodes,          /**< number of solution nodes or NULL */
   int*            nsoledges           /**< number of solution edges or NULL */
)
{
   int* queue;
   int nnodes;
   int nedges = (nsoledges != NULL)? *nsoledges : 0;
   int qstart;
   int qend;

   assert(scip != NULL && tails != NULL && heads != NULL && pcancestors != NULL && solnodemark != NULL);

   if( solnodequeue != NULL )
      queue = solnodequeue;
   else
      SCIP_CALL( SCIPallocBufferArray(scip, &queue, orgnnodes) );

   if( nsolnodes == NULL )
   {
      assert(solnodequeue == NULL);
      nnodes = 0;
      for( int k = 0; k < orgnnodes; k++ )
         if( solnodemark[k] )
            queue[nnodes++] = k;
   }
   else
   {
      nnodes = *nsolnodes;
      assert(solnodequeue != NULL);
   }

   qstart = 0;
   qend = nnodes;

   while( qend != qstart )
   {
      int k = qstart;

      assert(qstart < qend);
      qstart = qend;

      for( ; k < qend; k++ )
      {
         const int ancestornode = queue[k];

         assert(solnodemark[ancestornode]);

         for( IDX* curr = pcancestors[ancestornode]; curr != NULL; curr = curr->parent )
         {
            const int ancestoredge = curr->index;
            assert(tails[ancestoredge] < orgnnodes && heads[ancestoredge] < orgnnodes);

            if( soledgemark != NULL && !soledgemark[ancestoredge] )
            {
               soledgemark[ancestoredge] = TRUE;
               nedges++;
            }
            if( !solnodemark[tails[ancestoredge]] )
            {
               solnodemark[tails[ancestoredge]] = TRUE;
               queue[nnodes++] = tails[ancestoredge];
            }
            if( !solnodemark[heads[ancestoredge]] )
            {
               solnodemark[heads[ancestoredge]] = TRUE;
               queue[nnodes++] = heads[ancestoredge];
            }
         }
      }
      qend = nnodes;
   }

   if( nsolnodes != NULL )
      *nsolnodes = nnodes;

   if( nsoledges != NULL )
      *nsoledges = nedges;

   if( solnodequeue == NULL )
      SCIPfreeBufferArray(scip, &queue);

   return SCIP_OKAY;
}
