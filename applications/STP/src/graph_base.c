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

/**@file   graph_base.c
 * @brief  includes several methods for Steiner problem graphs
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file contains several basic methods to process Steiner problem graphs.
 * A graph can not be reduced in terms of edge or node size, but edges can be marked as
 * EAT_FREE (to not be used anymore) and nodes may have degree one.
 * The method 'graph_pack()' can be used to build a new graph that discards these nodes and edges.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h)                             */
//#define SCIP_DEBUG

#include "scip/misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "portab.h"
#include "misc_stp.h"
#include "graph.h"
#include "reduce.h"
#include "stpvector.h"
#include "heur_tm.h"


/** is the PC/MW part of the given graph valid? */
static
SCIP_Bool graphisValidPcMw(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the graph */
   SCIP_Bool*            nodevisited         /**< array */
   )
{
   int nFixedTerms = 0;
   int nProperTerms = 0;
   int nPseudoTerms = 0;
   int nNonLeafTerms = 0;
   const int nnodes = g->knots;
   const int root = g->source;
   const SCIP_Bool isExtended = g->extended;
   const SCIP_Bool isRooted = graph_pc_isRootedPcMw(g);

   assert(graph_pc_isPcMw(g));
   assert(g->prize != NULL);
   assert(g->term2edge != NULL);

   assert(graph_pc_term2edgeIsConsistent(scip, g));

   for( int k = 0; k < nnodes; k++ )
   {
      if( graph_pc_knotIsFixedTerm(g, k) )
      {
         assert(isRooted || k == root);

         nFixedTerms++;

         continue;
      }

      assert(k != root);

      if( graph_pc_knotIsNonLeafTerm(g, k) )
      {
         nNonLeafTerms++;
      }
      /* is k a terminal in the original graph? */
      else if( (isExtended ? Is_pseudoTerm(g->term[k]) : Is_term(g->term[k])) )
      {
         nPseudoTerms++;
      }
      /* is k an artificial (newly added) terminal? */
      else if( (isExtended ? Is_term(g->term[k]) : Is_pseudoTerm(g->term[k])) )
      {
         int e;
         int e2;
         int pterm = -1;
         const int term = k;

         nProperTerms++;

         if( g->grad[k] != 2 )
         {
            SCIPdebugMessage("terminal degree != 2 for %d \n", k);
            return FALSE;
         }

         for( e = g->inpbeg[term]; e != EAT_LAST; e = g->ieat[e] )
            if( g->tail[e] == root )
               break;

         if( e == EAT_LAST )
         {
            SCIPdebugMessage("no edge to root for term %d \n", term);
            return FALSE;
         }

         for( e2 = g->outbeg[term]; e2 != EAT_LAST; e2 = g->oeat[e2] )
         {
            pterm = g->head[e2];
            if( (isExtended ? Is_pseudoTerm(g->term[pterm]) : Is_term(g->term[pterm])) && pterm != root  )
               break;
         }

         if( e2 == EAT_LAST)
         {
            SCIPdebugMessage("no terminal for dummy %d \n", g->head[e2]);
            return FALSE;
         }

         assert(pterm >= 0);
         assert(pterm != root);

         if( e2 != g->term2edge[term] )
         {
            SCIPdebugMessage("term2edge for node %d faulty \n", term);
            return FALSE;
         }

         if( !EQ(g->cost[e], g->prize[pterm]) )
         {
            SCIPdebugMessage("prize mismatch for node %d: %f!=%f \n", pterm, g->cost[e], g->prize[pterm]);
            return FALSE;
         }
      }
   } /* for k = 0; k < nnodes; k++ */

   if( nProperTerms != nPseudoTerms )
   {
      SCIPdebugMessage("wrong nPseudoTerms terminal count: %d != %d \n", nProperTerms, nPseudoTerms);

      return FALSE;
   }

   if( isExtended )
   {
      /* non-leaf terms are not part of g->terms */
      if( (nProperTerms + nFixedTerms) != g->terms )
      {
         SCIPdebugMessage("wrong overall terminal count(1): %d != %d \n", nProperTerms + nFixedTerms, g->terms);

         return FALSE;
      }
   }
   else
   {
      if( (nProperTerms + nNonLeafTerms + nFixedTerms) != g->terms )
      {
         SCIPdebugMessage("wrong overall terminal count(2): %d != %d \n", nProperTerms + nNonLeafTerms + nFixedTerms - 1, g->terms);

         return FALSE;
      }
   }

   if( isRooted )
   {
      SCIP_CALL_ABORT( graph_trail_costAware(scip, g, g->source, nodevisited) );

      for( int k = 0; k < nnodes; k++ )
      {
         if( graph_pc_knotIsFixedTerm(g, k) && !nodevisited[k] )
         {
            SCIPdebugMessage("disconnected fixed terminal %d \n", k);
            return FALSE;
         }
      }
   }

   if( isRooted && graph_pc_isMw(g) )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( !graph_pc_knotIsFixedTerm(g, k))
            continue;

         for( int e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         {
            if( !EQ(g->cost[e], 0.0) )
            {
               if( k == g->source && graph_pc_knotIsDummyTerm(g, g->tail[e]) )
                  continue;

               SCIPdebugMessage("non-zero incoming arc for fixed MW terminal %d \n", k);
               return FALSE;
            }
         }
      }
   }

   return TRUE;
}


/** do changes for Pc/Mw variants for vanished graph */
static
SCIP_RETCODE packPcMwVanished(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g_old               /**< the old graph */
   )
{
   int term = -1;
   const int nnodes_old = g_old->knots;

   if( graph_pc_isRootedPcMw(g_old) )
      return SCIP_OKAY;

   for( int i = 0; i < nnodes_old; ++i )
   {
      if( graph_pc_knotIsNonLeafTerm(g_old, i) )
      {
         assert(-1 == term);
         assert(g_old->source != i);

         term = i;
         break;
      }
   }

   assert(term >= 0);
   assert(graph_pc_termIsNonLeafTerm(g_old, term));

   SCIP_CALL( graph_fixed_addNodePc(scip, term, g_old) );

   return SCIP_OKAY;
}


/** do changes for Pc/Mw variants */
static
SCIP_RETCODE packPcMwInit(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nnodes_new,         /**< number of nodes */
   GRAPH*                g_old,              /**< the old graph */
   GRAPH*                g_new               /**< the new graph */
   )
{
   assert(nnodes_new > 0);
   assert(graph_pc_isPcMw(g_old));

   assert(g_new->ksize == nnodes_new);

   SCIP_CALL( graph_pc_initSubgraph(scip, g_new) );

   if( g_old->stp_type == STP_BRMWCSP )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g_new->costbudget), nnodes_new) );

   return SCIP_OKAY;
}

/** add nodes to new graph during graph packing */
static
void packNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g_old,              /**< the old graph */
   GRAPH*                g_new               /**< the new graph */
   )
{
   const int oldnnodes = graph_get_nNodes(g_old);
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(g_old);
   const SCIP_Bool pcmw = graph_pc_isPcMw(g_old);

   if( rpcmw )
   {
      for( int i = 0; i < oldnnodes; i++ )
         g_old->mark[i] = (g_old->grad[i] > 0);

      for( int e = g_old->outbeg[g_old->source]; e != EAT_LAST; e =
            g_old->oeat[e] )
      {
         const int i = g_old->head[e];

         if( SCIPisGT(scip, g_old->cost[e], 0.0) && Is_term(g_old->term[i])
               && g_old->term2edge[i] >= 0 )
         {
            g_old->mark[i] = FALSE;
            assert(g_old->grad[i] == 2);
         }
      }
   }

   for( int i = 0; i < oldnnodes; i++ )
   {
      if( g_old->grad[i] > 0 )
      {
         if( pcmw )
         {
            if( !Is_term(g_old->term[i]) || (rpcmw && g_old->mark[i]) )
               g_new->prize[g_new->knots] = g_old->prize[i];
            else
               g_new->prize[g_new->knots] = 0.0;
         }

         if( g_old->stp_type == STP_BRMWCSP )
         {
            if( !Is_term(g_old->term[i]) || (g_old->mark[i]) )
               g_new->costbudget[g_new->knots] = g_old->costbudget[i];
            else
               g_new->costbudget[g_new->knots] = 0.0;
         }
         graph_knot_add(g_new, g_old->term[i]);
      }
   }
}


/** adds pseudo-ancestor to new graph during graph packing */
static
SCIP_RETCODE packPseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g_old,              /**< the old graph */
   GRAPH*                g_new               /**< the new graph */
   )
{
   const int oldnacestors = graph_getNpseudoAncestors(g_old);
   const int oldnedges = graph_get_nEdges(g_old);
   int e_new = 0;

   SCIP_CALL( graph_initPseudoAncestors(scip, g_new) );

   if( oldnacestors == 0 )
      return SCIP_OKAY;

   graph_addPseudoAncestors(oldnacestors, g_new);

   for( int e_old = 0; e_old < oldnedges; e_old += 2 )
   {
      if( g_old->ieat[e_old] != EAT_FREE )
      {
         const int* ancestors;
         const int nancestors = graph_edge_nPseudoAncestors(g_old, e_old);
         SCIP_Bool conflict = FALSE;

         if( nancestors == 0 )
         {
            e_new += 2;
            continue;
         }

         ancestors = graph_edge_getPseudoAncestors(g_old, e_old);
         assert(ancestors);

         SCIP_CALL( graph_pseudoAncestors_appendCopyArrayToEdge(scip, e_new, ancestors, nancestors, g_new, &conflict) );
         assert(!conflict);
         assert(graph_edge_nPseudoAncestors(g_new, e_new) == nancestors);
         assert(graph_edge_nPseudoAncestors(g_new, e_new + 1) == nancestors);

         e_new += 2;
      }
   }

   assert(e_new == g_new->edges);

 //  printf("added %d pseudo ancestors \n", oldnacestors);

   return SCIP_OKAY;
}


/** adds edges to new graph during graph packing */
static
SCIP_RETCODE packEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            old2newNode,        /**< node mapping */
   GRAPH*                g_old,              /**< the old graph */
   int                   nnodes,             /**< number of nodes for new graph */
   int                   nedges,             /**< number of edges for new graph */
   GRAPH*                g_new               /**< the new graph */
   )
{
   const int oldnedges = graph_get_nEdges(g_old);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_new->ancestors), nedges) );

   /* add edges */
   for( int e_old = 0; e_old < oldnedges; e_old += 2 )
   {
      int e_new;
      if( g_old->ieat[e_old] == EAT_FREE )
      {
         assert(g_old->oeat[e_old]     == EAT_FREE);
         assert(g_old->ieat[e_old + 1] == EAT_FREE);
         assert(g_old->oeat[e_old + 1] == EAT_FREE);

         graph_edge_delHistory(scip, g_old, e_old);
         continue;
      }

      assert(g_old->oeat[e_old]      != EAT_FREE);
      assert(g_old->ieat[e_old + 1]  != EAT_FREE);
      assert(g_old->oeat[e_old + 1]  != EAT_FREE);
      assert(old2newNode[g_old->tail[e_old]] >= 0);
      assert(old2newNode[g_old->head[e_old]] >= 0);

      e_new = g_new->edges;
      assert(e_new % 2 == 0);

      g_new->ancestors[e_new] = g_old->ancestors[e_old];
      g_new->ancestors[e_new + 1] = g_old->ancestors[e_old + 1];

      assert(old2newNode[g_old->tail[e_old]] < nnodes && old2newNode[g_old->head[e_old]] < nnodes);

      graph_edge_addSubgraph(scip, g_old, old2newNode, e_old, g_new);

      g_old->ancestors[e_old] = g_old->ancestors[e_old + 1] = NULL;
   }

   assert(nedges == g_new->edges);

   SCIP_CALL( packPseudoAncestors(scip, g_old, g_new) );

   return SCIP_OKAY;
}

/*
 * global functions
 */


/** builds complete graph (Kn) with unit edge weights and no terminals */
SCIP_RETCODE graph_buildCompleteGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               g,                  /**< new graph */
   int                   nnodes              /**< number of nodes */
   )
{
   const int nedges = nnodes * (nnodes - 1);
   GRAPH* graph;

   assert(scip && g);
   assert(nnodes >= 1);

   SCIP_CALL( graph_init(scip, g, nnodes, nedges, 1) );

   graph = *g;

   for( int k = 0; k < nnodes; k++ )
      graph_knot_add(graph, -1);

   for( int k = 0; k < nnodes - 1; k++ )
      for( int k2 = nnodes - 1; k2 >= k + 1; k2-- )
         graph_edge_add(scip, graph, k, k2, 1.0, 1.0);

   return SCIP_OKAY;
}


/* get compressed sparse row arrays representing current graph */
void graph_getCsr(
   const GRAPH*          g,                  /**< the graph */
   int* RESTRICT         edgearr,            /**< original edge array [0,...,nedges - 1] */
   int* RESTRICT         tailarr,            /**< tail of csr edge [0,...,nedges - 1]  */
   int* RESTRICT         start,              /**< start array [0,...,nnodes] */
   int*                  nnewedges           /**< pointer to store number of new edges */
      )
{
   int i = 0;
   const int nnodes = g->knots;

   assert(g != NULL);
   assert(tailarr != NULL);
   assert(edgearr != NULL);
   assert(start != NULL);

   for( int k = 0; k < nnodes; k++ )
   {
      start[k] = i;
      for( int e = g->inpbeg[k]; e >= 0; e = g->ieat[e] )
      {
         edgearr[i] = e;
         tailarr[i++] = g->tail[e] + 1;
      }
   }

   *nnewedges = i;
   start[nnodes] = i;
}


/** get edge costs */
void graph_getEdgeCosts(
   const GRAPH*          graph,              /**< the graph */
   SCIP_Real* RESTRICT   cost,               /**< reduced edge costs */
   SCIP_Real* RESTRICT   costrev             /**< reduced reverse edge costs */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real* const gcost = graph->cost;

   assert(cost && costrev);

   for( int e = 0; e < nedges; e++ )
   {
      cost[e] = gcost[e];
      costrev[e] = gcost[flipedge(e)];
      assert(GE(cost[e], 0.0));
   }
}


/** gets reversed edge costs */
void graph_getEdgeRevCosts(
   const GRAPH*          graph,              /**< the graph */
   SCIP_Real* RESTRICT   costrev             /**< reduced reverse edge costs */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real* const gcost = graph->cost;

   assert(costrev);

   for( int e = 0; e < nedges; e++ )
   {
      costrev[e] = gcost[flipedge(e)];
      assert(GE(costrev[e], 0.0));
   }
}


/* modifies 'isterm' to mark whether node is a terminal (or proper terminal for PC) */
void graph_getIsTermArray(
   const GRAPH*          g,                  /**< the graph */
   SCIP_Bool*            isterm              /**< marks whether node is a terminal (or proper terminal for PC) */
)
{
   const int nnodes = g->knots;
   const SCIP_Bool pcmw = graph_pc_isPcMw(g);

   assert(g && isterm);
   assert(!pcmw || !g->extended);

   for( int i = 0; i < nnodes; i++ )
   {
      isterm[i] = FALSE;

      if( Is_term(g->term[i]) )
      {
         if( pcmw && graph_pc_termIsNonLeafTerm(g, i) )
           continue;

         isterm[i] = TRUE;
      }
   }
}


/** gets terminals */
void graph_getTerms(
   const GRAPH*          g,                  /**< the graph */
   int*                  terms               /**< array of size g->terms */
)
{
   int nterms = 0;
   const int nnodes = graph_get_nNodes(g);

   assert(terms);

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) )
         terms[nterms++] = i;
   }

   assert(g->terms == nterms);
}


/** gets randomly permuted terminals */
SCIP_RETCODE graph_getTermsRandom(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   int*                  terms               /**< array of size g->terms */
)
{
   SCIP_RANDNUMGEN* rand;
   SCIP_CALL( SCIPcreateRandom(scip, &rand, g->terms, TRUE) );

   graph_getTerms(g, terms);
   SCIPrandomPermuteIntArray(rand, terms, 0, g->terms);

   SCIPfreeRandom(scip, &rand);

   return SCIP_OKAY;
}


/** initialize graph */
SCIP_RETCODE graph_init(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               g,                  /**< new graph */
   int                   ksize,              /**< slots for nodes */
   int                   esize,              /**< slots for edges */
   int                   layers              /**< number of layers (only needed for packing, otherwise 1) */
   )
{
   GRAPH* p;

   assert(ksize > 0);
   assert(ksize < INT_MAX);
   assert(esize >= 0);
   assert(esize < INT_MAX);
   assert(layers > 0);
   assert(layers < SHRT_MAX);

   SCIP_CALL( SCIPallocMemory(scip, g) );
   p = *g;
   assert(p != NULL);

   p->ancestors = NULL;
   p->pcancestors = NULL;
   p->fixedcomponents = NULL;
   p->orgtail = NULL;
   p->orghead = NULL;
   p->rootedgeprevs = NULL;
   p->norgmodelterms = 0;
   p->norgmodelknots = 0;
   p->norgmodeledges = 0;
   p->ksize  = ksize;
   p->orgknots = 0;
   p->orgedges = 0;
   p->knots  = 0;
   p->terms  = 0;
   p->grid_dim = -1;
   p->orgsource = UNKNOWN;
   p->stp_type = UNKNOWN;
   p->layers = layers;
   p->hoplimit = UNKNOWN;
   p->extended = FALSE;
   p->source = -1;
   p->is_packed = FALSE;
   p->withInexactReductions = FALSE;
   p->cost_org_pc = NULL;
   p->contracttrace = NULL;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->term), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mark), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->grad), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->inpbeg), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->outbeg), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->cost), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->tail), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->head), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->ieat), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->oeat), esize) );

   p->esize = esize;
   p->edges = 0;
   p->mincut_nnodes = 0;
   p->mincut_nedges = 0;
   p->prize  = NULL;
   p->maxdeg = NULL;
   p->grid_coordinates = NULL;
   p->grid_ncoords = NULL;
   p->mincut_dist = NULL;
   p->mincut_head = NULL;
   p->mincut_numb = NULL;
   p->mincut_prev = NULL;
   p->mincut_next = NULL;
   p->mincut_temp = NULL;
   p->mincut_e = NULL;
   p->mincut_x = NULL;
   p->mincut_r = NULL;
   p->path_heap = NULL;
   p->path_state = NULL;
   p->term2edge = NULL;
   p->budget = -1.0;
   p->costbudget = NULL;
   p->csr_storage = NULL;
   p->dcsr_storage = NULL;
   p->pseudoancestors = NULL;

   return SCIP_OKAY;
}

/** initialize data structures required to keep track of reductions */
SCIP_RETCODE graph_initHistory(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< graph */
   )
{
   IDX** pcancestors;        /* ancestor lists array (over all nodes) */
   const int nedges = graph_get_nEdges(graph);
   const int* tail = graph->tail;
   const int* head = graph->head;
   int* orgtail;
   int* orghead;
   const SCIP_Bool isPcMw = graph_pc_isPcMw(graph);

   assert(scip);
   assert(!graph->orgtail && !graph->orghead);
   assert(!graph->ancestors);
   assert(nedges > 0);

   SCIP_CALL( graph_initPseudoAncestors(scip, graph) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->orgtail), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->orghead), nedges) );

   orgtail = graph->orgtail;
   orghead = graph->orghead;

   BMScopyMemoryArray(orgtail, tail, nedges);
   BMScopyMemoryArray(orghead, head, nedges);

   if( isPcMw )
   {
      const int nnodes = graph->knots;

      SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->pcancestors), nnodes) );

      pcancestors = graph->pcancestors;

      for( int k = 0; k < nnodes; k++ )
         pcancestors[k] = NULL;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->ancestors), nedges) );
   SCIP_CALL( graph_initAncestors(scip, graph) );

   return SCIP_OKAY;
}


/** enlarge given graph */
SCIP_RETCODE graph_resize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph to be resized */
   int                   ksize,              /**< new node slots */
   int                   esize,              /**< new edge slots */
   int                   layers              /**< layers (set to -1 by default) */
   )
{
   assert(scip      != NULL);
   assert(g      != NULL);
   assert((ksize  < 0) || (ksize  >= g->knots));
   assert((esize  < 0) || (esize  >= g->edges));
   assert((layers < 0) || (layers >= g->layers));

   if( (layers > 0) && (layers != g->layers) )
      g->layers = layers;

   if( (ksize > 0) && (ksize != g->ksize) )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->term), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->mark), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->grad), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->inpbeg), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->outbeg), ksize) );

      if( g->prize != NULL)
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->prize), ksize) );

      if( g->costbudget != NULL)
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->costbudget), ksize) );

      if( g->term2edge != NULL)
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->term2edge), ksize) );

      g->ksize  = ksize;
   }
   if( (esize > 0) && (esize != g->esize) )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->cost), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->tail), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->head), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->ieat), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->oeat), esize) );

      g->esize = esize;
   }

   return SCIP_OKAY;
}


/** free the graph */
void graph_free(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph to be freed */
   SCIP_Bool             final               /**< delete ancestor data structures? */
   )
{
   GRAPH* p;

   assert(scip != NULL);
   assert(graph != NULL);

   p = *graph;
   assert(p != NULL);

   if( graph_mincut_isInitialized(p) )
      graph_mincut_exit(scip, p);

   if( graph_path_exists(p) )
      graph_path_exit(scip, p);

   graph_freeHistory(scip, p);

   if( final )
      graph_freeHistoryDeep(scip, p);

   SCIPfreeMemoryArrayNull(scip, &(p->prize));
   SCIPfreeMemoryArrayNull(scip, &(p->costbudget));
   SCIPfreeMemoryArrayNull(scip, &(p->term2edge));

   if( p->csr_storage != NULL )
      graph_free_csr(scip, p);

   if( p->dcsr_storage != NULL )
      graph_free_dcsr(scip, p);

   if( p->stp_type == STP_DCSTP )
   {
      SCIPfreeMemoryArray(scip, &(p->maxdeg));
   }
   else if( p->stp_type == STP_RSMT )
   {
      if( p->grid_coordinates != NULL )
      {
         assert(p->grid_coordinates != NULL);
         for( int i = p->grid_dim - 1; i >= 0;  i-- )
            SCIPfreeMemoryArray(scip, &(p->grid_coordinates[i]));

         SCIPfreeMemoryArray(scip, &(p->grid_coordinates));
      }

      if( p->grid_ncoords != NULL )
         SCIPfreeMemoryArray(scip, &(p->grid_ncoords));
   }

   SCIPfreeMemoryArray(scip, &(p->oeat));
   SCIPfreeMemoryArray(scip, &(p->ieat));
   SCIPfreeMemoryArray(scip, &(p->head));
   SCIPfreeMemoryArray(scip, &(p->tail));
   SCIPfreeMemoryArray(scip, &(p->cost));
   SCIPfreeMemoryArray(scip, &(p->outbeg));
   SCIPfreeMemoryArray(scip, &(p->inpbeg));
   SCIPfreeMemoryArray(scip, &(p->grad));
   SCIPfreeMemoryArray(scip, &(p->mark));
   SCIPfreeMemoryArray(scip, &(p->term));
   SCIPfreeMemoryArrayNull(scip, &(p->rootedgeprevs));
   SCIPfreeMemoryArrayNull(scip, &(p->cost_org_pc));

   SCIPfreeMemory(scip, graph);
}


/** frees the history */
void graph_freeHistory(
   SCIP*                 scip,               /**< SCIP data */
   GRAPH*                p                   /**< graph data */
   )
{
   SCIPfreeMemoryArrayNull(scip, &(p->contracttrace));

   if( p->pseudoancestors != NULL )
      graph_freePseudoAncestors(scip, p);

   if( p->ancestors != NULL )
   {
      const int nedges = p->edges;

      for( int e = nedges - 1; e >= 0; e-- )
      {
         IDX* curr = p->ancestors[e];
         while( curr != NULL )
         {
            p->ancestors[e] = curr->parent;
            SCIPfreeBlockMemory(scip, &(curr));
            curr = p->ancestors[e];
         }
      }
      SCIPfreeMemoryArray(scip, &(p->ancestors));
   }

   assert(p->pseudoancestors == NULL && p->ancestors == NULL);
}


/** frees the deep history */
void graph_freeHistoryDeep(
   SCIP*                 scip,               /**< SCIP data */
   GRAPH*                p                   /**< graph data */
   )
{
   assert(scip != NULL);
   assert(p != NULL);
   assert(p->path_heap == NULL);
   assert(p->path_state == NULL);

   if( p->pcancestors != NULL )
   {
      for( int e = p->norgmodelknots - 1; e >= 0; e-- )
      {
         IDX* curr = p->pcancestors[e];
         while( curr != NULL )
         {
            p->pcancestors[e] = curr->parent;
            SCIPfreeBlockMemory(scip, &(curr));
            curr = p->pcancestors[e];
         }
      }
      SCIPfreeMemoryArray(scip, &(p->pcancestors));
   }

   if( p->orgtail != NULL )
   {
      assert(p->orghead != NULL);

      SCIPfreeMemoryArray(scip, &(p->orghead));
      SCIPfreeMemoryArray(scip, &(p->orgtail));
   }

   if( p->fixedcomponents )
      graph_free_fixed(scip, p);
}


/** copy the graph */
SCIP_RETCODE graph_copy(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< original graph */
   GRAPH**               copygraph           /**< graph to be created */
   )
{
   const GRAPH* p = orggraph;
   assert(p != NULL);

   SCIP_CALL( graph_init(scip, copygraph, p->ksize, p->esize, p->layers) );

   SCIP_CALL( graph_copyData(scip, orggraph, *copygraph) );

   return SCIP_OKAY;
}


/** copy the data of the graph */
SCIP_RETCODE graph_copyData(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orgraph,            /**< original graph */
   GRAPH*                copygraph           /**< graph to be copied to */
   )
{
   GRAPH* g_copy = copygraph;
   const GRAPH* g_org = orgraph;
   const int ksize = g_org->ksize;
   const int esize = g_org->esize;

   assert(scip != NULL);
   assert(orgraph != NULL);
   assert(copygraph != NULL);
   assert(ksize == g_copy->ksize && ksize > 0);
   assert(esize == g_copy->esize && esize >= 0);
   assert(g_org->dcsr_storage == NULL && g_org->csr_storage == NULL);

   g_copy->norgmodeledges = g_org->norgmodeledges;
   g_copy->norgmodelknots = g_org->norgmodelknots;
   g_copy->norgmodelterms = g_org->norgmodelterms;
   g_copy->knots = g_org->knots;
   g_copy->terms = g_org->terms;
   g_copy->edges = g_org->edges;
   g_copy->source = g_org->source;
   g_copy->orgsource = g_org->orgsource;
   g_copy->orgedges = g_org->orgedges;
   g_copy->orgknots = g_org->orgknots;
   g_copy->grid_dim = g_org->grid_dim;
   g_copy->stp_type = g_org->stp_type;
   g_copy->hoplimit = g_org->hoplimit;
   g_copy->extended = g_org->extended;
   g_copy->budget = g_org->budget;
   g_copy->is_packed = g_org->is_packed;
   g_copy->withInexactReductions = g_org->withInexactReductions;

   BMScopyMemoryArray(g_copy->term, g_org->term, ksize);
   BMScopyMemoryArray(g_copy->mark, g_org->mark, ksize);
   BMScopyMemoryArray(g_copy->grad, g_org->grad, ksize);
   BMScopyMemoryArray(g_copy->inpbeg, g_org->inpbeg, ksize);
   BMScopyMemoryArray(g_copy->outbeg, g_org->outbeg, ksize);
   BMScopyMemoryArray(g_copy->cost, g_org->cost, esize);
   BMScopyMemoryArray(g_copy->tail, g_org->tail, esize);
   BMScopyMemoryArray(g_copy->head, g_org->head, esize);
   BMScopyMemoryArray(g_copy->ieat, g_org->ieat, esize);
   BMScopyMemoryArray(g_copy->oeat, g_org->oeat, esize);

   if( graph_pc_isPcMw(g_copy) )
   {
      const SCIP_Bool brpcmw = (g_copy->stp_type == STP_BRMWCSP);

      assert(g_copy->extended && g_org->extended);

      if( g_copy->costbudget != NULL )
      {
         assert(brpcmw);
         SCIPfreeMemoryArray(scip, &(g_copy->costbudget));
      }

      if( brpcmw )
      {
         SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->costbudget), g_copy->knots));
         BMScopyMemoryArray(g_copy->costbudget, g_org->costbudget, g_copy->knots);
      }

      SCIPfreeMemoryArrayNull(scip, &(g_copy->prize));
      SCIPfreeMemoryArrayNull(scip, &(g_copy->term2edge));
      SCIPfreeMemoryArrayNull(scip, &(g_copy->cost_org_pc));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->prize), g_copy->knots));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->term2edge), g_copy->knots));

      if( graph_pc_isPc(g_org) )
      {
         assert(g_org->cost_org_pc != NULL);
         SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->cost_org_pc), g_copy->edges));
         BMScopyMemoryArray(g_copy->cost_org_pc, g_org->cost_org_pc, g_copy->edges);
      }
      else
      {
         g_copy->cost_org_pc = NULL;
      }

      assert(g_org->prize != NULL);
      BMScopyMemoryArray(g_copy->prize, g_org->prize, g_copy->knots);

#ifndef NDEBUG
      for( int k = 0; k < g_copy->knots; k++ )
      {
         if( Is_term(g_org->term[k]) && (!graph_pc_isRootedPcMw(g_copy) || !graph_pc_knotIsFixedTerm(g_org, k)) )
            assert(EQ(g_copy->prize[k], 0.0));
      }
#endif

      assert(g_org->term2edge != NULL);
      BMScopyMemoryArray(g_copy->term2edge, g_org->term2edge, g_copy->knots);
   }
   else if( g_copy->stp_type == STP_DCSTP )
   {
      assert(g_org->maxdeg != NULL);

      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->maxdeg), g_copy->knots));

      for( int k = 0; k < g_copy->knots; k++ )
         g_copy->maxdeg[k] = g_org->maxdeg[k];
   }
   else if( g_org->stp_type == STP_RSMT )
   {
      assert(g_org->grid_ncoords != NULL);
      assert(g_org->grid_coordinates != NULL);
      assert(g_org->grid_dim > 0);

      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->grid_coordinates), g_org->grid_dim));

      BMScopyMemoryArray(g_copy->grid_coordinates, g_org->grid_coordinates, g_org->grid_dim);
      for( int k = 0; k < g_org->grid_dim; k++ )
      {
         SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->grid_coordinates[k]), g_org->terms)); /*lint !e866*/
         BMScopyMemoryArray(g_copy->grid_coordinates[k], g_org->grid_coordinates[k], g_org->terms); /*lint !e866*/
      }
      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->grid_ncoords), g_org->grid_dim));

      BMScopyMemoryArray(g_copy->grid_ncoords, g_org->grid_ncoords, g_org->grid_dim);
   }

   assert(graph_valid(scip, g_copy));

   return SCIP_OKAY;
}


/** copies the pseudo-ancestors */
SCIP_RETCODE graph_copyPseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< original graph */
   GRAPH*                copygraph           /**< graph to be created */
   )
{
   const int nedges = graph_get_nEdges(orggraph);
   const int nancestors_all = graph_getNpseudoAncestors(orggraph);
   assert(orggraph && copygraph);
   assert(orggraph->edges == copygraph->edges);
   assert(graph_getNpseudoAncestors(copygraph) == 0);

   if( nancestors_all == 0 )
      return SCIP_OKAY;

   graph_addPseudoAncestors(nancestors_all, copygraph);

   for( int e = 0; e < nedges; e += 2 )
   {
      if( orggraph->ieat[e] != EAT_FREE )
      {
         const int* ancestors;
         const int nancestors = graph_edge_nPseudoAncestors(orggraph, e);
         SCIP_Bool conflict = FALSE;

         if( nancestors == 0 )
            continue;

         ancestors = graph_edge_getPseudoAncestors(orggraph, e);
         assert(ancestors);

         SCIP_CALL( graph_pseudoAncestors_appendCopyArrayToEdge(scip, e, ancestors, nancestors, copygraph, &conflict) );
         assert(!conflict);
      }
   }

   return SCIP_OKAY;
}



/** marks the current graph */
void graph_mark(
   GRAPH*                g                   /**< the graph */
   )
{
   const int nnodes = graph_get_nNodes(g);
   const int root = g->source;
   const int* const grad = g->grad;
   int* RESTRICT isMarked = g->mark;

   for( int k = 0; k < nnodes; k++ )
      isMarked[k] = (grad[k] > 0);

   isMarked[root] = TRUE;

   if( graph_pc_isPcMw(g) && !g->extended )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( Is_pseudoTerm(g->term[k]) )
            isMarked[k] = FALSE;
         else if( Is_term(g->term[k]) )
            isMarked[k] = TRUE;
      }

      if( graph_pc_isRootedPcMw(g) )
         isMarked[root] = TRUE;
      else
         isMarked[root] = FALSE;
   }

   assert(graph_isMarked(g));
}


/** is the current graph properly marked? */
SCIP_Bool graph_isMarked(
   const GRAPH*          g                   /**< the graph */
   )
{
   assert(g);

   if( graph_pc_isPcMw(g) && !g->extended )
   {
      const int nnodes = g->knots;
      const int root = g->source;
      const int rooted = graph_pc_isRootedPcMw(g);

      for( int k = 0; k < nnodes; k++ )
      {
         if( Is_pseudoTerm(g->term[k]) || (!rooted && k == root) )
         {
            assert(g->grad[k] == 2 || k == root);

            if( g->mark[k] )
            {
               graph_knot_printInfo(g, k);
               printf("pseudo-terminal %d is marked \n", k);
               return FALSE;
            }
         }
         else
         {
            if( g->mark[k] != (g->grad[k] > 0) )
            {
               if( k == root && g->mark[k] )
                  continue;

               if( Is_term(g->term[k]) && g->mark[k] )
                  continue;

               graph_knot_printInfo(g, k);
               printf("node %d: mark=%d grad=%d \n", k, g->mark[k], g->grad[k]);
               return FALSE;
            }
         }
      }

      if( rooted )
      {
         if( !g->mark[root] )
         {
            printf("root not marked \n");
            return FALSE;
         }
      }
   }
   else
   {
      const int nnodes = g->knots;

      for( int k = 0; k < nnodes; k++ )
      {
         if( g->mark[k] != (g->grad[k] > 0) )
         {
            if( k == g->source && g->mark[k] )
               continue;

            graph_knot_printInfo(g, k);
            printf("node %d: mark=%d grad=%d \n", k, g->mark[k], g->grad[k]);
            return FALSE;
         }
      }
   }

   return TRUE;
}



/** is the current graph already set up? (with history and path) */
SCIP_Bool graph_isSetUp(
   const GRAPH*          g                   /**< the graph */
   )
{
   assert(g);
   assert((g->orgtail == NULL) == (g->ancestors == NULL));

   return (g->ancestors != NULL);
}


void graph_show(
   const GRAPH*          p                   /**< the graph */
   )
{
   int i;

   assert(p != NULL);

   for(i = 0; i < p->knots; i++)
      if (p->grad[i] > 0)
         (void)printf("Knot %d, term=%d, grad=%d, inpbeg=%d, outbeg=%d\n",
            i, p->term[i], p->grad[i], p->inpbeg[i], p->outbeg[i]);

   (void)fputc('\n', stdout);

   for(i = 0; i < p->edges; i++)
      if (p->ieat[i] != EAT_FREE)
         (void)printf("Edge %d, cost=%g, tail=%d, head=%d, ieat=%d, oeat=%d\n",
            i, p->cost[i], p->tail[i], p->head[i], p->ieat[i], p->oeat[i]);

   (void)fputc('\n', stdout);
}


/** reinsert all hidden edges */
void graph_uncover(
   GRAPH*                g                   /**< the graph */
   )
{/*lint --e{850}*/
   int head;
   int tail;
   int e;

   assert(g      != NULL);

   for( e = 0; e < g->edges; e++ )
   {
      if( g->ieat[e] == EAT_HIDE )
      {
         assert(e % 2 == 0);
         assert(g->oeat[e] == EAT_HIDE);

         head            = g->head[e];
         tail            = g->tail[e];

         g->grad[head]++;
         g->grad[tail]++;

         g->ieat[e]      = g->inpbeg[head];
         g->oeat[e]      = g->outbeg[tail];
         g->inpbeg[head] = e;
         g->outbeg[tail] = e;

         e++;

         assert(g->ieat[e] == EAT_HIDE);
         assert(g->oeat[e] == EAT_HIDE);
         assert(g->head[e] == tail);
         assert(g->tail[e] == head);

         head            = g->head[e];
         tail            = g->tail[e];
         g->ieat[e]      = g->inpbeg[head];
         g->oeat[e]      = g->outbeg[tail];
         g->inpbeg[head] = e;
         g->outbeg[tail] = e;
      }
   }
}


/** Packs the graph, i.e. builds a new graph that discards deleted edges and nodes.
 *  The original graph is deleted. */
SCIP_RETCODE graph_pack(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   REDSOL*               redsol,             /**< reduce solution or NULL */
   SCIP_Bool             verbose             /**< verbose? */
   )
{
   GRAPH *g_old;
   GRAPH *g_new;
   int* old2newNode;
   const int oldnnodes = graph_get_nNodes(graph);
   const int oldnedges = graph_get_nEdges(graph);
   int nnodes = 0;
   int nedges = 0;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
   SCIP_Bool graphHasVanished = FALSE;

   assert(scip && newgraph);
   assert(graph_valid(scip, graph));
   assert(!graph->is_packed);

   g_old = graph;

   SCIP_CALL( SCIPallocBufferArray(scip, &old2newNode, oldnnodes) );

   if( verbose )
      printf("Reduced graph: ");

   /* build node mapping */
   for( int i = 0; i < oldnnodes; i++ )
   {
      if( g_old->grad[i] > 0 )
         old2newNode[i] = nnodes++;
      else
         old2newNode[i] = -1;
   }

   /* graph vanished? */
   if( nnodes == 0 )
   {
      SCIPfreeBufferArray(scip, &old2newNode);
      old2newNode = NULL;
      if( verbose )
         printf(" graph vanished!\n");

      nnodes = 1;
      graphHasVanished = TRUE;

      if( pcmw )
      {
         SCIP_CALL( packPcMwVanished(scip, g_old) );
      }
   }

   /* count surviving edges */
   for( int i = 0; i < oldnedges; i++ )
   {
      if( g_old->oeat[i] != EAT_FREE )
      {
         assert(g_old->ieat[i] != EAT_FREE);
         nedges++;
      }
   }

   assert(nnodes > 1 || nedges == 0);
   SCIP_CALL( graph_init(scip, newgraph, nnodes, nedges, g_old->layers) );
   g_new = *newgraph;
   g_new->norgmodelknots = g_old->norgmodelknots;
   g_new->norgmodelterms = g_old->norgmodelterms;
   g_new->norgmodeledges = g_old->norgmodeledges;
   g_new->orgsource = g_old->orgsource;
   g_new->orgtail = g_old->orgtail;
   g_new->orghead = g_old->orghead;
   g_new->orgknots = g_old->knots;
   g_new->orgedges = g_old->edges;
   g_new->stp_type = g_old->stp_type;
   g_new->maxdeg = g_old->maxdeg;
   g_new->grid_dim = g_old->grid_dim;
   g_new->grid_ncoords = g_old->grid_ncoords;
   g_new->grid_coordinates = g_old->grid_coordinates;
   g_new->pcancestors = g_old->pcancestors;
   g_new->fixedcomponents = g_old->fixedcomponents;
   g_new->hoplimit = g_old->hoplimit;
   g_new->extended = g_old->extended;
   g_new->budget = g_old->budget;
   g_new->is_packed = TRUE;

   if( graphHasVanished )
   {
      assert(!old2newNode);

      g_new->ancestors = NULL;
      graph_free(scip, &g_old, FALSE);

      if( g_new->stp_type == STP_RSMT )
      {
         g_new->grid_ncoords = NULL;
         g_new->grid_coordinates = NULL;
      }

      graph_knot_add(g_new, 0);
      g_new->source = 0;

      return SCIP_OKAY;
   }

   assert(nnodes >= 2 && nedges >= 1);

   if( pcmw )
   {
      SCIP_CALL( packPcMwInit(scip, nnodes, g_old, g_new) );
   }

   if( redsol  )
   {
      reduce_solPack(g_old, old2newNode, nnodes, redsol);
   }

   /* add nodes (of positive degree) to new graph */
   packNodes(scip, g_old, g_new);

   /* add root */
   assert(Is_term(g_new->term[old2newNode[g_old->source]]));
   g_new->source = old2newNode[g_old->source];
   assert(!graph_pc_isRootedPcMw(graph) || FARAWAY == g_new->prize[g_new->source]);

   /* NOTE: also handles ancestors */
   SCIP_CALL( packEdges(scip, old2newNode, g_old, nnodes, nedges, g_new) );

   SCIPfreeBufferArray(scip, &old2newNode);

   SCIP_CALL( graph_pc_finalizeSubgraph(scip, g_new) );

   if( graph_pc_isPc(g_new) )
   {
      SCIP_Real* offset = reduce_solGetOffsetPointer(redsol);
      *offset += graph_pc_getNonLeafTermOffset(scip, g_new);
   }

   if( g_old->path_heap != NULL )
      graph_path_exit(scip, g_old);

   g_old->stp_type = UNKNOWN;
   graph_free(scip, &g_old, FALSE);

   assert(graph_valid(scip, g_new));

   if( verbose )
      printf("Nodes: %d  Edges: %d  Terminals: %d\n", g_new->knots, g_new->edges, g_new->terms);

   return SCIP_OKAY;
}


/** is the given graph valid? */
SCIP_Bool graph_valid(
   SCIP*                 scip,               /**< scip struct */
   const GRAPH*          g                   /**< the graph */
   )
{
   int nterms;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   SCIP_Bool isValid = TRUE;
   SCIP_Bool* nodevisited = NULL;

   assert(scip && g);
   nterms = g->terms;

   for( int k = 0; k < nnodes; k++ )
   {
      int e;

      if( Is_term(g->term[k]) )
         nterms--;

      for( e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         if( g->head[e] != k )
            break;

      if( e != EAT_LAST )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: Head invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n",
               k, e, g->tail[e], g->head[e]);

         goto EXIT;
      }

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         if( g->tail[e] != k )
            break;

      if( e != EAT_LAST )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: Tail invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n",
               k, e, g->tail[e], g->head[e]);

         goto EXIT;
      }
   }

   if( nterms != 0 )
   {
      isValid = FALSE;
      SCIPdebugMessage("*** Graph invalid: Wrong Terminal count, count is %d, should be %d\n",
            g->terms, g->terms - nterms);

      goto EXIT;
   }

   if( (g->source < 0 ) || (g->source >= g->knots) || (g->term[g->source] != 0) )
   {
      isValid = FALSE;
      SCIPdebugMessage("*** Graph invalid: Source invalid, Layer %d, Source %d, Terminal %d\n",
            0, g->source, g->term[g->source]);

      goto EXIT;
   }

   for( int e = 0; e < nedges; e += 2 )
   {
      if( (g->ieat[e] == EAT_FREE) && (g->oeat[e] == EAT_FREE)
         && (g->ieat[e + 1] == EAT_FREE) && (g->oeat[e + 1] == EAT_FREE) )
      {
         continue;
      }

      if( (g->ieat[e] == EAT_FREE) || (g->oeat[e] == EAT_FREE)
         || (g->ieat[e + 1] == EAT_FREE) || (g->oeat[e + 1] == EAT_FREE) )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: FREE invalid, Edge %d/%d\n",
               e, e + 1);

         goto EXIT;
      }

      if( (g->head[e] != g->tail[e + 1]) || (g->tail[e] != g->head[e + 1]) )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: Anti invalid, Edge %d/%d, Tail=%d/%d, Head=%d/%d\n",
               e, e + 1, g->head[e], g->tail[e + 1],  g->tail[e], g->head[e + 1]);

         goto EXIT;
      }
   }

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nodevisited, nnodes) );
   SCIP_CALL_ABORT( graph_trail_arr(scip, g, g->source, nodevisited) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( (g->grad[k] == 0) && ((g->inpbeg[k] != EAT_LAST) || (g->outbeg[k] != EAT_LAST)) )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: Knot %d with Grad 0 has Edges\n", k);

         goto EXIT;
      }

      if( !nodevisited[k] && ((g->grad[k] > 0) || (Is_term(g->term[k]))) )
      {
         if( graph_pc_isPcMw(g) && !graph_pc_knotIsFixedTerm(g, k) )
         {
            continue;
         }

         isValid = FALSE;
#ifdef SCIP_DEBUG
         graph_knot_printInfo(g, k);
         SCIPdebugMessage("*** Graph invalid: Knot %d not connected\n", k);
#endif

         goto EXIT;
      }
   }

   if( isValid && graph_pc_isPcMw(g) )
   {
      isValid = graphisValidPcMw(scip, g, nodevisited);
   }


   EXIT:

   SCIPfreeBufferArrayNull(scip, &nodevisited);
   return isValid;
}



/** is the given input graph valid? */
SCIP_Bool graph_validInput(
   SCIP*                 scip,               /**< scip struct */
   const GRAPH*          g                   /**< the graph */
   )
{
   const int nnodes = graph_get_nNodes(g);
   int nterms = g->terms;
   SCIP_Bool isValid = TRUE;
   SCIP_Bool* nodevisited = NULL;

   assert(scip && g);

   for( int k = 0; k < nnodes; k++ )
   {
      int e;

      if( Is_term(g->term[k]) )
         nterms--;

      for( e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         if( g->head[e] != k )
            break;

      if( e != EAT_LAST )
      {
         isValid = FALSE;
         SCIPerrorMessage("Input graph invalid: Head invalid, node %d, edge %d->%d \n",
               k, g->tail[e] + 1, g->head[e] + 1);

         goto EXIT;
      }

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         if( g->tail[e] != k )
            break;

      if( e != EAT_LAST )
      {
         isValid = FALSE;
         SCIPerrorMessage("Input graph invalid: Tail invalid, node %d, edge %d->%d \n",
               k, g->tail[e] + 1, g->head[e] + 1);

         goto EXIT;
      }
   }

   if( nterms != 0 )
   {
      isValid = FALSE;
      SCIPerrorMessage("Input graph invalid: Wrong number of terminals, count is %d, should be %d\n",
            g->terms, g->terms - nterms);

      goto EXIT;
   }

   if( (g->source < 0 ) || (g->source >= g->knots) || (g->term[g->source] != 0) )
   {
      isValid = FALSE;
      SCIPerrorMessage("Input graph invalid: Root invalid, root %d, terminal state %d\n",
            g->source + 1, g->term[g->source]);

      goto EXIT;
   }

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nodevisited, nnodes) );
   SCIP_CALL_ABORT( graph_trail_arr(scip, g, g->source, nodevisited) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( (g->grad[k] == 0) && ((g->inpbeg[k] != EAT_LAST) || (g->outbeg[k] != EAT_LAST)) )
      {
         isValid = FALSE;
         SCIPerrorMessage("Input graph invalid: Node %d of degree 0 has edges \n", k + 1);

         goto EXIT;
      }

      if( !nodevisited[k] && ((g->grad[k] > 0) || (Is_term(g->term[k]))) )
      {
         if( graph_pc_isPcMw(g) && !graph_pc_knotIsFixedTerm(g, k) )
         {
            continue;
         }

         isValid = FALSE;
         SCIPerrorMessage("Input graph invalid: Node %d not connected to terminal node %d \n", k + 1, g->source + 1);

         goto EXIT;
      }
   }

   if( isValid && graph_pc_isPcMw(g) )
   {
      isValid = graphisValidPcMw(scip, g, nodevisited);

      if( !isValid )
      {
         SCIPerrorMessage("Input graph invalid (non-standard Steiner tree specific error) \n");
      }
   }


   EXIT:

   SCIPfreeBufferArrayNull(scip, &nodevisited);
   return isValid;
}
