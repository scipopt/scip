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

/**@file   graph_node.c
 * @brief  includes graph node methods for Steiner problems
 * @author Daniel Rehfeldt
 *
 * Graph node methods for Steiner problems
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */

#include "graph.h"
#include "portab.h"



/* is the vertex a leaf (for NWPTSPG) */
SCIP_Bool graph_knotIsNWLeaf(
   const GRAPH*          g,                  /**< the graph */
   int                   vertex              /**< the vertex  */
)
{
   int e;
   assert(g != NULL && g->stp_type == STP_NWPTSPG);

   for( e = g->outbeg[vertex]; e != EAT_LAST; e = g->oeat[e] )
      if( LT(g->cost[e], FARAWAY) )
         break;

   return (e == EAT_LAST );
}


/** is node in range? */
SCIP_Bool graph_knot_isInRange(
   const GRAPH*          g,                  /**< the graph */
   int                   k                   /**< the node */
   )
{
   assert(g);

   return (0 <= k && k < g->knots);
}


/** add a vertex */
void graph_knot_add(
   GRAPH*                p,                  /**< the graph */
   int                   term                /**< terminal property */
   )
{
   assert(p        != NULL);
   assert(p->ksize >  p->knots);
   assert(term     <  p->layers);

   p->term  [p->knots] = term;
   p->mark  [p->knots] = TRUE;
   p->grad  [p->knots] = 0;
   p->inpbeg[p->knots] = EAT_LAST;
   p->outbeg[p->knots] = EAT_LAST;

   if( Is_term(term) )
      p->terms++;

   p->knots++;
}

/** change terminal property of a vertex */
void graph_knot_chg(
   GRAPH*                p,                  /**< the graph */
   int                   node,               /**< node to be changed */
   int                   term                /**< terminal property */
   )
{
   assert(p      != NULL);
   assert(node   >= 0);
   assert(node   < p->knots);
   assert(term == STP_TERM || term == STP_TERM_NONE || term == STP_TERM_NONLEAF || term == STP_TERM_PSEUDO);

   if( term != p->term[node] )
   {
      if( Is_term(p->term[node]) )
         p->terms--;

      p->term[node] = term;

      if( Is_term(p->term[node]) )
         p->terms++;
   }
}


/** deletes node */
void graph_knot_del(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   k,                  /**< the node */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   assert(scip && g);
   assert(graph_knot_isInRange(g, k));

   while( g->outbeg[k] != EAT_LAST )
      graph_edge_del(scip, g, g->outbeg[k], freeancestors);

   assert(g->grad[k] == 0);
   assert(g->outbeg[k] == EAT_LAST && g->inpbeg[k] == EAT_LAST);
}



/** deletes node, and also adapts DCSR storage */
void graph_knot_delFull(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   k,                  /**< the node */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   assert(scip && g);
   assert(graph_knot_isInRange(g, k));

   if( g->dcsr_storage )
   {
      while( g->outbeg[k] != EAT_LAST )
         graph_edge_delFull(scip, g, g->outbeg[k], freeancestors);
   }
   else
   {
      while( g->outbeg[k] != EAT_LAST )
         graph_edge_del(scip, g, g->outbeg[k], freeancestors);
   }

   assert(g->grad[k] == 0);
   assert(g->outbeg[k] == EAT_LAST && g->inpbeg[k] == EAT_LAST);
}


/** pseudo deletes non-terminal of degree 2 */
SCIP_RETCODE graph_knot_replaceDeg2(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   vertex,             /**< the vertex to replace */
   SCIP_Real             edgecost,           /**< new edge cost */
   int                   ancestornode,       /**< node to copy ancestors from or -1 */
   GRAPH*                g,                  /**< the graph */
   SCIP_Bool*            edgeEliminated      /**< edge eliminated? (due to conflict) */
   )
{
   IDX* ancestors1;
   IDX* ancestors2;
   int* pseudoancestors1 = NULL;
   int* pseudoancestors2 = NULL;
   const int e1 = g->outbeg[vertex];
   const int e2 = g->oeat[e1];
   const int i1 = g->head[e1];
   const int i2 = g->head[e2];
   const int npseudoancestors1 = graph_edge_nPseudoAncestors(g, e1);
   const int npseudoancestors2 = graph_edge_nPseudoAncestors(g, e2);
   int newedge;
   SCIP_Bool conflict = FALSE;

   assert(scip && g);
   assert(vertex >= 0 && vertex < g->knots);
   assert(!Is_term(g->term[vertex]));
   assert(g->grad[vertex] == 2);
   assert(e1 >= 0 && e2 >= 0);
   assert(SCIPisEQ(scip, g->cost[e1], g->cost[flipedge(e1)]));
   assert(SCIPisEQ(scip, g->cost[e2], g->cost[flipedge(e2)]));
   assert(graph_valid_pseudoAncestors(scip, g));
   assert(graph_typeIsUndirected(g));

   *edgeEliminated = FALSE;

   if( npseudoancestors1 > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(pseudoancestors1), npseudoancestors1) );
      BMScopyMemoryArray(pseudoancestors1, graph_edge_getPseudoAncestors(g, e1), npseudoancestors1);
   }

   if( npseudoancestors2 > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(pseudoancestors2), npseudoancestors2) );
      BMScopyMemoryArray(pseudoancestors2, graph_edge_getPseudoAncestors(g, e2), npseudoancestors2);
   }

   ancestors1 = graph_edge_getAncestors(g, e1);
   ancestors2 = graph_edge_getAncestors(g, e2);

   g->ancestors[e1] = g->ancestors[flipedge(e1)] = NULL;
   g->ancestors[e2] = g->ancestors[flipedge(e2)] = NULL;

   newedge = graph_edge_redirect(scip, g, e1, i2, i1, edgecost, FALSE, TRUE);

   assert(ancestornode >= 0 || ancestornode == -1);

   /* is there a new edge? */
   if( newedge >= 0 )
   {
      const int edge_even = Edge_even(newedge);

      graph_edge_delHistory(scip, g, newedge);

      g->ancestors[edge_even] = ancestors1;
      SCIPintListNodeAppend(g->ancestors[edge_even], ancestors2);

      SCIP_CALL( graph_pseudoAncestors_appendCopyArrayToEdge(scip, newedge, pseudoancestors1, npseudoancestors1, g, &conflict) );
      assert(!conflict);
      SCIP_CALL( graph_pseudoAncestors_appendCopyArrayToEdge(scip, newedge, pseudoancestors2, npseudoancestors2, g, &conflict) );

      /* ancestor node given?*/
      if( ancestornode >= 0 )
      {
         assert(graph_pc_isPcMw(g));

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[edge_even]), g->pcancestors[ancestornode], NULL) );

         if( !conflict )
            SCIP_CALL( graph_pseudoAncestors_appendCopyNodeToEdge(scip, newedge, ancestornode, TRUE, g, &conflict) );
      }
   }
   else
   {
      SCIPintListNodeFree(scip, &ancestors1);
      SCIPintListNodeFree(scip, &ancestors2);
   }

   graph_knot_del(scip, g, vertex, TRUE);

   SCIPfreeBlockMemoryArrayNull(scip, &pseudoancestors1, npseudoancestors1);
   SCIPfreeBlockMemoryArrayNull(scip, &pseudoancestors2, npseudoancestors2);

   if( conflict )
   {
      SCIPdebugMessage("conflict in graph_knot_replaceDeg2 \n");
      assert(newedge >= 0);
      graph_edge_del(scip, g, newedge, TRUE);
      *edgeEliminated = TRUE;
   }

   return SCIP_OKAY;
}



/** contracts node s into node t */
SCIP_RETCODE graph_knot_contract(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                or NULL */
   int                   t,                  /**< tail node to be contracted (surviving node) */
   int                   s                   /**< head node to be contracted */
   )
{
   SCIP_Real* incost = NULL;
   SCIP_Real* outcost = NULL;
   SINGLETONANS* ancestors = NULL;
   int* mark = NULL;
   int* edge = NULL;
   int* knot = NULL;
   int slc = 0;
   int sgrad;
   const SCIP_Bool isUndirected = graph_typeIsUndirected(g);

   assert(g && scip);
   assert(graph_knot_isInRange(g, t));
   assert(graph_knot_isInRange(g, s));
   assert(t != s);

   /* save solution nodes? */
   if( solnode )
   {
      if( solnode[s] == CONNECT )
         solnode[t] = CONNECT;
   }

   /* trace contractions? */
   if( g->contracttrace )
   {
      assert(g->contracttrace[s] == -1);
      g->contracttrace[s] = t;
      SCIPdebugMessage("contract trace %d->%d \n", s, t);
   }

   /* change terminal property */
   if( Is_term(g->term[s]) )
   {
      graph_knot_chg(g, t, g->term[s]);
      graph_knot_chg(g, s, -1);
   }

   /* retain root */
   if( g->source == s )
      g->source = t;

   sgrad = g->grad[s];
   if( sgrad >= 1 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &incost, sgrad) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &outcost, sgrad) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &mark, sgrad) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &edge, sgrad) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &knot, sgrad) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ancestors, sgrad) );
   }

   /* store edges to be moved/removed */
   for( int es = g->outbeg[s]; es != EAT_LAST; es = g->oeat[es] )
   {
      assert(g->tail[es] == s);

      if( g->head[es] != t )
      {
         assert(ancestors && mark && incost && outcost && edge && knot);
         assert(slc < sgrad);

         SCIP_CALL( graph_singletonAncestors_init(scip, g, es, &(ancestors[slc])) );
         mark[slc] = FALSE;
         edge[slc] = es;
         knot[slc] = g->head[es];
         outcost[slc] = g->cost[es];
         incost[slc] = g->cost[Edge_anti(es)];
         slc++;
      }
   }

   assert(slc == sgrad - 1 || slc == sgrad);

   /* traverse edges */
   for( int i = 0; i < slc; i++ )
   {
      int et;
      const int ihead = knot[i];

      assert(knot != NULL && outcost != NULL && incost != NULL && mark != NULL);

      /* search for an edge out of t with same head as current edge */

      if( g->grad[ihead] >= g->grad[t] )
      {
         for( et = g->outbeg[t]; et >= 0; et = g->oeat[et] )
            if( g->head[et] == ihead )
               break;
      }
      else
      {
         for( et = g->inpbeg[ihead]; et >= 0; et = g->ieat[et] )
            if( g->tail[et] == t )
               break;
      }

      /* does such an edge not exist? */
      if( et == EAT_LAST )
      {
         mark[i] = TRUE;
      }
      else
      {
         const int anti = flipedge_Uint(et);
         SCIP_Bool copyPseudoancestors = FALSE;
         assert(et != EAT_LAST);

         /* This is for nodes with edges to s and t.
          * Need to adjust the out and in costs of the edge
          */

         if( isUndirected && SCIPisGT(scip, g->cost[et], outcost[i]) && SCIPisGT(scip, g->cost[anti], incost[i]) )
            copyPseudoancestors = TRUE;

         if( copyPseudoancestors )
            graph_edge_delPseudoAncestors(scip, et, g);

         /* NOTE: even in the undirected case there might be different anti-parallel weights for PC/MW */
         if( isUndirected )
         {
            assert((SCIPisGT(scip, g->cost[et], outcost[i]) && SCIPisGT(scip, g->cost[anti], incost[i]))
                || (SCIPisLE(scip, g->cost[et], outcost[i]) && SCIPisLE(scip, g->cost[anti], incost[i])));

            if( SCIPisGT(scip, g->cost[et], outcost[i]) )
            {
               const int even = Edge_even(et);

               assert(!g->ancestors[et] || !g->ancestors[anti]);
               assert(g->ancestors[even]);
               assert(SCIPisGT(scip, g->cost[anti], incost[i]));

               SCIPintListNodeFree(scip, &(g->ancestors[even]));
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[even]), ancestors[i].ancestors, NULL) );

               g->cost[et] = outcost[i];
               g->cost[anti] = incost[i];
            }
         }
         else
         {
            if( SCIPisGT(scip, g->cost[et], outcost[i]) )
            {
               SCIPintListNodeFree(scip, &((g->ancestors)[et]));
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((g->ancestors)[et]), ancestors[i].ancestors, NULL) );

               assert(graph_edge_nPseudoAncestors(g, et) == 0);
               g->cost[et] = outcost[i];
            }

            if( SCIPisGT(scip, g->cost[anti], incost[i]) )
            {
               SCIPintListNodeFree(scip, &(g->ancestors[anti]));
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((g->ancestors)[anti]), ancestors[i].revancestors, NULL) );

               assert(graph_edge_nPseudoAncestors(g, anti) == 0);
               g->cost[anti] = incost[i];
            }
         }

         if( copyPseudoancestors )
         {
            SCIP_Bool conflict;
            SCIP_CALL( graph_pseudoAncestors_appendCopySingToEdge(scip, et, &(ancestors[i]), FALSE, g, &conflict) );
            assert(!conflict);
         }
      }
   }

   /* insert edges */
   for( int i = 0; i < slc; i++ )
   {
      assert(mark && ancestors && knot && outcost && incost);

      if( mark[i] )
      {
         int es = g->outbeg[s];
         const int head = knot[i];
         const int tail = t;
         SCIP_Bool conflict;

         assert(es != EAT_LAST);

         graph_edge_del(scip, g, es, TRUE);
         assert(!g->ancestors[es] && !g->ancestors[flipedge(es)]);

         if( isUndirected )
         {
            const int even = Edge_even(es);
            assert(ancestors[i].ancestors && !ancestors[i].revancestors);
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[even]), ancestors[i].ancestors, NULL) );
         }

         if( !isUndirected )
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[es]), ancestors[i].ancestors, NULL) );

         SCIP_CALL( graph_pseudoAncestors_appendCopySingToEdge(scip, es, &(ancestors[i]), FALSE, g, &conflict) );
         assert(!conflict);

         g->grad[head]++;
         g->grad[tail]++;

         g->cost[es]     = outcost[i];
         g->tail[es]     = tail;
         g->head[es]     = head;
         g->ieat[es]     = g->inpbeg[head];
         g->oeat[es]     = g->outbeg[tail];
         g->inpbeg[head] = es;
         g->outbeg[tail] = es;

         es = Edge_anti(es);

         if( !isUndirected )
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[es]), ancestors[i].revancestors, NULL) );

         g->cost[es]     = incost[i];
         g->tail[es]     = head;
         g->head[es]     = tail;
         g->ieat[es]     = g->inpbeg[tail];
         g->oeat[es]     = g->outbeg[head];
         g->inpbeg[tail] = es;
         g->outbeg[head] = es;
      }
   }

   /* delete remaining edges */
   graph_knot_del(scip, g, s, TRUE);

   if( sgrad >= 1 )
   {
      assert(ancestors && knot && outcost && edge && mark && incost);

      for( int i = 0; i < slc; i++ )
         graph_singletonAncestors_freeMembers(scip, &(ancestors[i]));

      SCIPfreeBlockMemoryArray(scip, &ancestors, sgrad);
      SCIPfreeBlockMemoryArray(scip, &knot, sgrad);
      SCIPfreeBlockMemoryArray(scip, &edge, sgrad);
      SCIPfreeBlockMemoryArray(scip, &mark, sgrad);
      SCIPfreeBlockMemoryArray(scip, &outcost, sgrad);
      SCIPfreeBlockMemoryArray(scip, &incost, sgrad);
   }

   return SCIP_OKAY;
}


/** contract an edge, given index and by its endpoints, which is to be fixed */
SCIP_RETCODE graph_knot_contractFixed(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                or NULL */
   int                   edge,               /**< the edge */
   int                   t,                  /**< tail node to be contracted (surviving node) */
   int                   s                   /**< head node to be contracted */
   )
{
   SCIP_CALL( graph_fixed_addEdge(scip, edge, g) );
   SCIP_CALL( graph_knot_contract(scip, g, solnode, t, s) );

   return SCIP_OKAY;
}


/** contract endpoint of lower degree into endpoint of higher degree */
SCIP_RETCODE graph_knot_contractLowdeg2High(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                or NULL */
   int                   t,                  /**< tail node to be contracted */
   int                   s                   /**< head node to be contracted */
   )
{
   assert(g != NULL);

   if( g->grad[t] >= g->grad[s] )
      SCIP_CALL( graph_knot_contract(scip, g, solnode, t, s) );
   else
      SCIP_CALL( graph_knot_contract(scip, g, solnode, s, t) );

   return SCIP_OKAY;
}
