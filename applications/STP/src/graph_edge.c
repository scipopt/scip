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

/**@file   graph_edge.c
 * @brief  includes graph edge methods for Steiner problems
 * @author Daniel Rehfeldt
 *
 * Graph edge methods for Steiner problems
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */

#include "graph.h"
#include "portab.h"

/*
 * Local methods
 */

inline static
void removeEdge(
   GRAPH*                g,                  /**< the graph */
   int                   e                   /**< the edge to be removed */
   )
{
   int i;
   const int head = g->head[e];
   const int tail = g->tail[e];

   assert(graph_edge_isInRange(g, e));

   if( g->inpbeg[head] == e )
   {
      g->inpbeg[head] = g->ieat[e];
   }
   else
   {
      if( g->rootedgeprevs != NULL && head == g->source )
      {
         i = g->rootedgeprevs[e];
         assert(g->ieat[i] == e);
         if( g->ieat[e] >= 0 )
            g->rootedgeprevs[g->ieat[e]] = i;
      }
      else
      {
         for( i = g->inpbeg[head]; g->ieat[i] != e; i = g->ieat[i] )
         {
            assert(i >= 0);
         }
      }

      g->ieat[i] = g->ieat[e];
   }

   if( g->outbeg[tail] == e )
   {
      g->outbeg[tail] = g->oeat[e];
   }
   else
   {
      if( g->rootedgeprevs != NULL && tail == g->source )
      {
         i = g->rootedgeprevs[e];
         assert(g->oeat[i] == e);
         if( g->oeat[e] >= 0 )
            g->rootedgeprevs[g->oeat[e]] = i;
      }
      else
      {
         for( i = g->outbeg[tail]; g->oeat[i] != e; i = g->oeat[i] )
         {
            assert(i >= 0);
         }
      }

      g->oeat[i] = g->oeat[e];
   }
}

/*
 * Interface methods
 */


/** redirects given edge eki */
int graph_edge_redirect(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   eki,                /**< the edge */
   int                   k,                  /**< new tail */
   int                   j,                  /**< new head */
   SCIP_Real             cost,               /**< new cost */
   SCIP_Bool             forcedelete,        /**< delete edge eki if it is not used? */
   SCIP_Bool             checkexist          /**< check if there is already an edge kj */
   )
{
   int e;

   assert(scip && g);
   assert(SCIPisGE(scip, cost, 0.0));

   if( forcedelete )
      graph_edge_del(NULL, g, eki, FALSE);

   if( checkexist )
   {
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(g->tail[e] == k);

         if( g->head[e] == j )
            break;
      }
   }
   else
   {
#ifndef NDEBUG
      /* make sure that the edge does not exist! */
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(g->tail[e] == k);

         if( g->head[e] == j )
            break;
      }
      assert(e == EAT_LAST);
#endif

      e = EAT_LAST;
   }

   /* does edge already exist? */
   if( e != EAT_LAST )
   {
      /* correct cost */
      if( SCIPisGT(scip, g->cost[e], cost) )
      {
         g->cost[e]            = cost;
         g->cost[Edge_anti(e)] = cost;
      }
      else
      {
         e = -1;
      }
   }
   else
   {
      assert(graph_edge_isInRange(g, eki));

      if( !forcedelete && g->ieat[eki] != EAT_FREE )
         graph_edge_del(NULL, g, eki, FALSE);

      assert(g->oeat[eki] == EAT_FREE);

      e = eki;

      g->grad[k]++;
      g->grad[j]++;

      g->cost[e]   = cost;
      g->head[e]   = j;
      g->tail[e]   = k;
      g->ieat[e]   = g->inpbeg[j];
      g->oeat[e]   = g->outbeg[k];
      g->inpbeg[j] = e;
      g->outbeg[k] = e;

      e = Edge_anti(eki);

      g->cost[e]   = cost;
      g->head[e]   = k;
      g->tail[e]   = j;
      g->ieat[e]   = g->inpbeg[k];
      g->oeat[e]   = g->outbeg[j];
      g->inpbeg[k] = e;
      g->outbeg[j] = e;
      return eki;
   }
   return e;
}

/** reinsert an edge to replace two other edges */
SCIP_RETCODE graph_edge_reinsert(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   edge,               /**< edge to reinsert */
   int                   tail,               /**< new tail */
   int                   head,               /**< new head */
   SCIP_Real             cost,               /**< new edge cost */
   int                   ancestornode,       /**< node to copy ancestors from or -1 */
   SINGLETONANS*         ancestorsBackward,  /**< backwards (edge) ancestors */
   SINGLETONANS*         ancestorsForward,   /**< forward (edge) ancestors */
   int*                  insertedge,         /**< pointer to inserted edge or -1 */
   SCIP_Bool*            conflict            /**< does the newly edge contain conflicts? (i.e. is  redundant) */
   )
{
   const int newedge = graph_edge_redirect(scip, g, edge, tail, head, cost, FALSE, TRUE);

   assert(ancestorsBackward && ancestorsForward && insertedge && conflict);
   assert(ancestornode >= 0 || ancestornode == -1);

   *conflict = FALSE;

   /* is there a new edge? */
   if( newedge >= 0 )
   {
      IDX* const ancestorsBack = ancestorsBackward->ancestors;
      IDX* const revancestorsBack = ancestorsBackward->revancestors;
      IDX* const ancestorsFor = ancestorsForward->ancestors;
      IDX* const revancestorsFor = ancestorsForward->revancestors;

      assert(ancestorsBack && ancestorsFor);
      assert((revancestorsBack == NULL) == graph_typeIsUndirected(g));
      assert((revancestorsFor == NULL) == graph_typeIsUndirected(g));

      graph_edge_delHistory(scip, g, newedge);

      if( !revancestorsBack )
      {
         const int edge_even = Edge_even(newedge);
         assert(graph_typeIsUndirected(g));
         assert(edge_even == flipedge(newedge) || edge_even == newedge);

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[edge_even]), ancestorsBack, NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[edge_even]), ancestorsFor, NULL) );
      }
      else
      {
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[newedge]), revancestorsBack, NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[newedge]), ancestorsFor, NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(newedge)]), ancestorsBack, NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(newedge)]), revancestorsFor, NULL) );
      }

      SCIP_CALL( graph_pseudoAncestors_appendCopySingToEdge(scip, newedge, ancestorsBackward, TRUE, g, conflict) );
      assert(!(*conflict));

      SCIP_CALL( graph_pseudoAncestors_appendCopySingToEdge(scip, newedge, ancestorsForward, TRUE, g, conflict) );

      /* ancestor node given?*/
      if( ancestornode >= 0 )
      {
         const int edge_even = Edge_even(newedge);

         assert(graph_pc_isPcMw(g));

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[edge_even]), g->pcancestors[ancestornode], NULL) );

         if( !(*conflict) )
            SCIP_CALL( graph_pseudoAncestors_appendCopyNodeToEdge(scip, newedge, ancestornode, TRUE, g, conflict) );
      }
   }

   *insertedge = newedge;

   return SCIP_OKAY;
}


/** add a new edge to the graph */
void graph_edge_add(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   tail,               /**< tail of the new edge */
   int                   head,               /**< head of the new edge*/
   SCIP_Real             cost1,              /**< tail to head cost */
   SCIP_Real             cost2               /**< head to tail cost */
   )
{
   int    e;

   assert(g      != NULL);
   assert(SCIPisGE(scip, cost1, 0.0) || SCIPisEQ(scip, cost1, (double) UNKNOWN));
   assert(SCIPisGE(scip, cost2, 0.0) || SCIPisEQ(scip, cost2, (double) UNKNOWN));
   assert(tail   >= 0);
   assert(tail   <  g->knots);
   assert(head   >= 0);
   assert(head   <  g->knots);

   assert(g->esize >= g->edges + 2);

   e = g->edges;

   g->grad[head]++;
   g->grad[tail]++;

   if( cost1 != UNKNOWN )
      g->cost[e]           = cost1;
   g->tail[e]           = tail;
   g->head[e]           = head;
   g->ieat[e]           = g->inpbeg[head];
   g->oeat[e]           = g->outbeg[tail];
   g->inpbeg[head]      = e;
   g->outbeg[tail]      = e;

   e++;

   if( cost2 != UNKNOWN )
      g->cost[e]           = cost2;
   g->tail[e]           = head;
   g->head[e]           = tail;
   g->ieat[e]           = g->inpbeg[tail];
   g->oeat[e]           = g->outbeg[head];
   g->inpbeg[tail]      = e;
   g->outbeg[head]      = e;

   g->edges += 2;
}

/** add a bi-edge to the graph */
void graph_edge_addBi(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   tail,               /**< tail of the new edge */
   int                   head,               /**< head of the new edge */
   SCIP_Real             cost                /**< head to tail cost */
   )
{
   graph_edge_add(scip, g, tail, head, cost, cost);
}


/** add a new edge to a subgraph */
void graph_edge_addSubgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< original graph */
   const int*            nodemapOrg2sub,     /**< node mapping from original to subgraph */
   int                   orgedge,            /**< original edge */
   GRAPH*                subgraph            /**< the subgraph */
   )
{
   const int e = orgedge;
   const int orgtail = orggraph->tail[e];
   const int orghead = orggraph->head[e];
   const int newtail = nodemapOrg2sub[orgtail];
   const int newhead = nodemapOrg2sub[orghead];

   assert(e >= 0 && e < orggraph->edges);

   if( graph_pc_isPcMw(orggraph) )
   {
      assert(orggraph->extended);
      graph_pc_updateSubgraphEdge(orggraph, nodemapOrg2sub, e, subgraph);
   }

   graph_edge_add(scip, subgraph, newtail, newhead, orggraph->cost[e], orggraph->cost[flipedge(e)]);
}


/** delete an edge from standard data structure */
void graph_edge_del(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   e,                  /**< the edge */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   assert(g);
   assert(e >= 0 && e < g->edges);

   if( freeancestors && g->ancestors )
   {
      graph_edge_delHistory(scip, g, e);
   }

   /* delete first arc */
   e -= e % 2;
   assert(g->head[e] == g->tail[e + 1]);
   assert(g->tail[e] == g->head[e + 1]);

   g->grad[g->head[e]]--;
   g->grad[g->tail[e]]--;

   removeEdge(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_FREE;
   g->oeat[e] = EAT_FREE;

   /* delete second arc */
   e++;
   removeEdge(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_FREE;
   g->oeat[e] = EAT_FREE;

   assert(g->tail[e] >= 0 && g->tail[e] < g->knots);
   assert(g->head[e] >= 0 && g->head[e] < g->knots);
}

/** delete an edge from standard, DCSR (if existent) and CSR (if existent) data structures */
void graph_edge_delFull(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   e,                  /**< the edge */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   assert(scip && g);
   assert(e >= 0 && e < g->edges);

   if( g->dcsr_storage )
   {
      int csredge;
      assert(graph_valid_dcsr(g, FALSE));

      csredge = g->dcsr_storage->id2csredge[e];
      graph_dcsr_deleteEdgeBi(scip, g->dcsr_storage, csredge);

      assert(graph_valid_dcsr(g, FALSE));
   }

   if( g->csr_storage )
   {
      assert(0 && "not yet supported");
   }

   graph_edge_del(scip, g, e, freeancestors);
}

/** deletes edges marked by given array */
void graph_edge_delBlocked(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   const SCIP_Bool*      edge_deletable,     /**< marks edges to delete (of size nedges / 2) */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   const int nedges = g->edges;

   assert(g && edge_deletable);

   for( int e = 0; e < nedges / 2; e++ )
      if( edge_deletable[e] )
         graph_edge_del(scip, g, e * 2, freeancestors);
}


/** free the history of an edge */
void graph_edge_delHistory(
   SCIP*                 scip,               /**< SCIP data */
   GRAPH*                g,                  /**< graph data */
   int                   edge                /**< edge for which to free the history */
   )
{
   assert(scip && g);
   assert(edge >= 0 && edge < g->edges);

   SCIPintListNodeFree(scip, &(g->ancestors[edge]));
   SCIPintListNodeFree(scip, &(g->ancestors[flipedge_Uint(edge)]));

   graph_edge_delPseudoAncestors(scip, edge, g);
}


/** hide edge */
void graph_edge_hide(
   GRAPH*                g,                  /**< the graph */
   int                   e                   /**< the edge */
   )
{
   assert(g          != NULL);
   assert(e          >= 0);
   assert(e          <  g->edges);

   /* Immer mit der ersten von beiden Anfangen
    */
   e -= e % 2;

   assert(g->head[e] == g->tail[e + 1]);
   assert(g->tail[e] == g->head[e + 1]);

   g->grad[g->head[e]]--;
   g->grad[g->tail[e]]--;

   removeEdge(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_HIDE;
   g->oeat[e] = EAT_HIDE;

   e++;

   removeEdge(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_HIDE;
   g->oeat[e] = EAT_HIDE;
}
