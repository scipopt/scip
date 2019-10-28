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

/**@file   graph_pcbase.c
 * @brief  includes several methods for prize-collecting Steiner problem graphs
 * @author Daniel Rehfeldt
 *
 * This file contains several basic methods to process prize-collecting Steiner problem graphs and kinsmen
 * such as the maximum-weight connected subgraph problem.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scip/misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "portab.h"
#include "misc_stp.h"
#include "graph.h"



/*
 * local functions
 */


/** is terminal a non-leaf? Checked by evaluation of the current graph */
static
SCIP_Bool termIsNonLeafTerm_evaluate(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the graph */
   int                   term                /**< terminal to be checked */
   )
{
   SCIP_Bool isNonLeaf = TRUE;

   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(Is_term(g->term[term]));

   for( int e = g->inpbeg[term]; e != EAT_LAST; e = g->ieat[e] )
   {
      const int neighbor = g->tail[e];

      if( !graph_pc_knotIsDummyTerm(g, neighbor) && SCIPisLT(scip, g->cost[e], g->prize[term]) )
      {
         assert(neighbor != g->source || graph_pc_isRootedPcMw(g));
         isNonLeaf = FALSE;
         break;
      }
   }

   return isNonLeaf;
}


/** is vertex a non-leaf (call before graph transformation was performed)  */
static inline
SCIP_Bool isNonLeaf_pretrans(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the graph */
   int                   vertex              /**< node check */
)
{
   const SCIP_Real prize = g->prize[vertex];
   const SCIP_Real* const cost = g->cost;

   for( int e = g->inpbeg[vertex]; e != EAT_LAST; e = g->ieat[e] )
   {
      if( SCIPisGT(scip, prize, cost[e]) )
         return FALSE;
   }

   return TRUE;
}


/** remove non-leaf terminals by edge weight shifting (call before graph transformation was performed,
 *  call only from graph transformation method!) */
static
void markNonLeafTerms_pretrans(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g                   /**< the graph */
)
{
   const int nnodes = g->knots;

   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(g->term[k]) )
         continue;

      if( isNonLeaf_pretrans(scip, g, k) )
      {
         graph_knot_chg(g, k, STP_TERM_NONLEAF);
      }
   }
}


/** remove non-leaf terminals by edge weight shifting (call before graph transformation was performed)  */
static
void markNonLeafTerms_2trans(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g                   /**< the graph */
)
{
   const int nnodes = g->knots;

   assert(!g->extended);

   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(g->term[k]) )
         continue;

      if( graph_pc_termIsNonLeafTerm(g, k) )
      {
         graph_knot_chg(g, k, STP_TERM_NONLEAF);
      }
   }
}


/** shift costs of non-leaf terminals (call right after transformation to extended has been performed)  */
static
void shiftNonLeafCosts_2trans(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g                   /**< the graph */
)
{
   const int nnodes = g->knots;
   SCIP_Real* const cost = g->cost;

#ifndef NDEBUG
   assert(g->cost_org_pc);
   assert(graph_pc_isPc(g));
   assert(g->extended);

   for( int e = 0; e < g->edges; e++ )
   {
      if( g->oeat[e] == EAT_LAST )
         continue;

      assert(SCIPisEQ(scip, cost[e], cost[flipedge(e)]) || SCIPisGE(scip, cost[e], FARAWAY)
               || SCIPisGE(scip, cost[flipedge(e)], FARAWAY));
   }
#endif

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_nonleafTerm(g->term[k]) )
      {
         const SCIP_Real prize = g->prize[k];

         assert(SCIPisPositive(scip, prize));

         for( int e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         {
            assert(SCIPisLT(scip, cost[e], FARAWAY));
            assert(SCIPisLT(scip, cost[flipedge(e)], FARAWAY));

            cost[e] -= prize;
            assert(SCIPisGE(scip, cost[e], 0.0));

            if( cost[e] < 0.0 )
               cost[e] = 0.0;
         }
      }
   }
}


/** initializes cost_org_pc array (call right after transformation to extended has been performed)  */
static
SCIP_RETCODE initCostOrgPc(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g                   /**< the graph */
)
{
   const int nedges = g->edges;
   SCIP_Real* const cost = g->cost;

   assert(!g->cost_org_pc);
   assert(graph_pc_isPc(g));
   assert(g->extended);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->cost_org_pc), nedges) );
   BMScopyMemoryArray(g->cost_org_pc, cost, nedges);

   return SCIP_OKAY;
}


/** contract an edge of rooted prize-collecting Steiner tree problem or maximum-weight connected subgraph problem
 *  such that this edge is incident to least one fixed terminal */
static
SCIP_RETCODE contractEdgeWithFixedEnd(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int*                  solnode,            /**< solution nodes or NULL */
   int                   t,                  /**< tail node to be contracted (surviving node) */
   int                   s,                  /**< head node to be contracted */
   int                   ets                 /**< edge from t to s */
)
{
   assert(ets >= 0);
   assert(g->tail[ets] == t && g->head[ets] == s);
   assert(graph_pc_isRootedPcMw(g));
   assert(graph_pc_knotIsFixedTerm(g, s) || graph_pc_knotIsFixedTerm(g, t));

   SCIP_CALL(graph_fixed_moveNodePc(scip, s, g));
   SCIP_CALL(graph_fixed_addEdge(scip, ets, g));

   if( !graph_pc_knotIsFixedTerm(g, s) )
   {
      assert(graph_pc_knotIsFixedTerm(g, t));
      if( Is_term(g->term[s]) )
         graph_pc_deleteTermExtension(scip, g, s);
   }
   else if( !graph_pc_knotIsFixedTerm(g, t) )
   {
      /* need to make t a fixed term */
      assert(g->source != t);
      graph_pc_deleteTermExtension(scip, g, t);
      assert(g->prize[s] == FARAWAY);
      g->prize[t] = FARAWAY;
   }

   /* contract s into t */
   SCIP_CALL(graph_knot_contract(scip, g, solnode, t, s));
   g->term2edge[t] = TERM2EDGE_FIXEDTERM;

   return SCIP_OKAY;
}



/** contract an edge of (rooted) prize-collecting Steiner tree problem or maximum-weight connected subgraph problem
 *  such that this edge is NOT incident to least one fixed terminal */
static
SCIP_RETCODE contractEdgeNoFixedEnd(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int*                  solnode,            /**< solution nodes or NULL */
   int                   t,                  /**< tail node to be contracted (surviving node) */
   int                   s,                  /**< head node to be contracted */
   int                   ets,                /**< edge from t to s */
   int                   term4offset         /**< terminal to add offset to */
)
{
   assert(ets >= 0);
   assert(g->tail[ets] == t && g->head[ets] == s);
   assert(Is_term(g->term[term4offset]));
   assert(!graph_pc_knotIsFixedTerm(g, s) && !graph_pc_knotIsFixedTerm(g, t));
   assert(!(g->term[s] && !g->term[t])); // todo might lead to issues with term2edge

   SCIP_CALL( graph_pc_contractNodeAncestors(scip, g, t, s, ets) );

   /* are both end-points of the edge to be contracted terminals? */
   if( Is_term(g->term[t]) && Is_term(g->term[s]) )
   {
      if( !graph_pc_termIsNonLeafTerm(g, s) )
         graph_pc_deleteTermExtension(scip, g, s);

      if( !graph_pc_knotIsFixedTerm(g, term4offset) )
         graph_pc_subtractPrize(scip, g, g->cost[ets] - g->prize[s], term4offset);
   }
   else
   {
      if( !graph_pc_knotIsFixedTerm(g, term4offset) )
      {
         if( g->stp_type != STP_MWCSP && g->stp_type != STP_RMWCSP )
            graph_pc_subtractPrize(scip, g, g->cost[ets], term4offset);
         else
            graph_pc_subtractPrize(scip, g, -(g->prize[s]), term4offset);
      }
      g->term2edge[s] = TERM2EDGE_NOTERM;
   }

   /* contract s into t */
   SCIP_CALL( graph_knot_contract(scip, g, solnode, t, s) );

   if( Is_term(g->term[t]) )
   {
      if( termIsNonLeafTerm_evaluate(scip, g, t) )
      {
         if( g->term2edge[t] != TERM2EDGE_NONLEAFTERM )
         {
            assert(g->term2edge[t] >= 0);

            graph_pc_deleteTermExtension(scip, g, t);
            g->term2edge[t] = TERM2EDGE_NONLEAFTERM;
         }
      }
      else
      {
         assert(TERM2EDGE_NONLEAFTERM != g->term2edge[t] && "currently not supported");
         /* NOTE: this case could happen, but not with current contraction calls; would need to keep edges for extension in this case */
      }
   }
   else
   {
      assert(g->term2edge[t] == TERM2EDGE_NOTERM);
   }

   return SCIP_OKAY;
}

/*
 * global functions
 */


#if 0
/** transforms an MWCSP to an SAP */
SCIP_RETCODE graph_MwcsToSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real*            maxweights          /**< array containing the weight of each node */
   )
{
   int e;
   int i;
   int nnodes;
   int nterms = 0;

   assert(maxweights != NULL);
   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->cost != NULL);
   assert(graph->terms == 0);

   nnodes = graph->knots;

   /* count number of terminals, modify incoming edges for non-terminals */
   for( i = 0; i < nnodes; i++ )
   {
      if( SCIPisLT(scip, maxweights[i], 0.0) )
      {
         for( e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
         {
            graph->cost[e] -= maxweights[i];
         }
      }
      else
      {
         graph_knot_chg(graph, i, 0);
         nterms++;
      }
   }
   nterms = 0;
   for( i = 0; i < nnodes; i++ )
   {
      graph->prize[i] = maxweights[i];
      if( Is_term(graph->term[i]) )
      {
         assert(SCIPisGE(scip, maxweights[i], 0.0));
         nterms++;
      }
      else
      {
         assert(SCIPisLT(scip, maxweights[i], 0.0));
      }
   }
   assert(nterms == graph->terms);
   graph->stp_type = STP_MWCSP;

   SCIP_CALL( graph_PcToSap(scip, graph) );
   assert(graph->stp_type == STP_MWCSP);
   return SCIP_OKAY;
}


/** alters the graph for prize collecting problems */
SCIP_RETCODE graph_PcToSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   SCIP_Real* prize;
   int k;
   int root;
   int node;
   int nnodes;
   int nterms;
   int pseudoroot;

   assert(graph != NULL);
   assert(graph->prize != NULL);
   assert(graph->knots == graph->ksize);
   assert(graph->edges == graph->esize);

   prize = graph->prize;
   nnodes = graph->knots;
   nterms = graph->terms;
   graph->norgmodelknots = nnodes;
   graph->norgmodeledges = graph->edges;

   /* for each terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + graph->terms + 2), (graph->esize + graph->terms * 8) , -1) );

   /* create a new nodes */
   for( k = 0; k < nterms; ++k )
      graph_knot_add(graph, -1);

   /* new pseudo-root */
   pseudoroot = graph->knots;
   graph_knot_add(graph, -1);

   /* new root */
   root = graph->knots;
   graph_knot_add(graph, 0);

   nterms = 0;
   for( k = 0; k < nnodes; ++k )
   {
      /* is the kth node a terminal other than the root? */
      if( Is_term(graph->term[k]) )
      {
         /* the copied node */
         node = nnodes + nterms;
         nterms++;

         /* switch the terminal property, mark k */
         graph_knot_chg(graph, k, -2);
         graph_knot_chg(graph, node, 0);
         assert(SCIPisGE(scip, prize[k], 0.0));

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, k, BLOCKED, FARAWAY);
         graph_edge_add(scip, graph, pseudoroot, node, prize[k], FARAWAY);
         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);
         graph_edge_add(scip, graph, k, pseudoroot, 0.0, FARAWAY);
      }
      else if( graph->stp_type != STP_MWCSP )
      {
         prize[k] = 0;
      }
   }
   graph->source = root;
   graph->extended = TRUE;
   assert((nterms + 1) == graph->terms);
   if( graph->stp_type != STP_MWCSP )
      graph->stp_type = STP_PCSPG;

   return SCIP_OKAY;
}
#endif


/** allocates (first and second) and initializes (only second) arrays for PC and MW problems */
SCIP_RETCODE graph_pc_init(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   sizeprize,          /**< size of prize array to allocate (or -1) */
   int                   sizeterm2edge       /**< size of term2edge array to allocate and initialize to -1 (or -1) */
   )
{
   assert(scip != NULL);
   assert(g != NULL);

   if( sizeprize > 0 )
   {
      assert(NULL == g->prize);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->prize), sizeprize) );
   }

   if( sizeterm2edge > 0 )
   {
      assert(NULL == g->term2edge);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->term2edge), sizeterm2edge) );
      for( int i = 0; i < sizeterm2edge; i++ )
         g->term2edge[i] = TERM2EDGE_NOTERM;
   }

   return SCIP_OKAY;
}

/** changes graph of PC and MW problems needed for presolving routines */
SCIP_RETCODE graph_pc_presolInit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   int prev;
   const int root = g->source;
   const int nedges = g->edges;

   if( g->stp_type == STP_RPCSPG )
      return SCIP_OKAY;

   assert(scip != NULL && g != NULL);
   assert(g->rootedgeprevs == NULL);
   assert(nedges > 0 && g->grad[root] > 0);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->rootedgeprevs), nedges) );

   for( int e = 0; e < nedges; e++ )
      g->rootedgeprevs[e] = -1;

   prev = g->outbeg[root];
   assert(prev != EAT_LAST);

   for( int e = g->oeat[prev]; e != EAT_LAST; e = g->oeat[e] )
   {
      g->rootedgeprevs[e] = prev;
      prev = e;
   }

   prev = g->inpbeg[root];
   assert(prev != EAT_LAST);

   for( int e = g->ieat[prev]; e != EAT_LAST; e = g->ieat[e] )
   {
      g->rootedgeprevs[e] = prev;
      prev = e;
   }

   return SCIP_OKAY;
}

/** changes graph of PC and MW problems needed after exiting presolving routines */
void graph_pc_presolExit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   assert(scip != NULL && g != NULL);

   if( g->stp_type == STP_RPCSPG )
      return;

   assert(g->rootedgeprevs != NULL);

   SCIPfreeMemoryArray(scip, &(g->rootedgeprevs));
}

/** checks consistency of term2edge array  */
SCIP_Bool graph_pc_term2edgeIsConsistent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< the graph */
)
{
   const int root = g->source;
   const SCIP_Bool isExtended = g->extended;

   assert(scip && g && g->term2edge);
   assert(graph_pc_isPcMw(g));

   if( g->term2edge[root] != TERM2EDGE_FIXEDTERM )
   {
      SCIPdebugMessage("term2edge root consistency \n");
      return FALSE;
   }

   for( int i = 0; i < g->knots; i++ )
   {
      if( i == root )
         continue;

      if( !isExtended && Is_term(g->term[i]) )
      {
          const SCIP_Bool isNonLeaf = (g->term2edge[i] == TERM2EDGE_NONLEAFTERM);
          const SCIP_Bool isNonLeaf_evaluate = termIsNonLeafTerm_evaluate(scip, g, i);

          if( isNonLeaf != isNonLeaf_evaluate )
          {
             SCIPdebugMessage("erm2edge consistency fail0 %d \n", i);
             return FALSE;
          }
      }

      if( Is_anyTerm(g->term[i])
            && !graph_pc_knotIsFixedTerm(g, i)
            && !graph_pc_knotIsNonLeafTerm(g, i)
            && g->term2edge[i] < 0 )
      {
         SCIPdebugMessage("term2edge consistency fail1 %d \n", i);
         return FALSE;
      }

      if( !Is_anyTerm(g->term[i]) && g->term2edge[i] != TERM2EDGE_NOTERM )
      {
         SCIPdebugMessage("term2edge consistency fail2 %d \n", i);
         return FALSE;
      }

      if( isExtended && Is_nonleafTerm(g->term[i]) && g->term2edge[i] != TERM2EDGE_NONLEAFTERM )
      {
         SCIPdebugMessage("term2edge consistency fail2b %d \n", i);
         return FALSE;
      }

      if( !isExtended && g->term2edge[i] == TERM2EDGE_NONLEAFTERM && !Is_term(g->term[i]) )
      {
         SCIPdebugMessage("term2edge consistency fail2c %d \n", i);
         return FALSE;
      }

      if( Is_pseudoTerm(g->term[i]) )
      {
         int k = -1;
         int e;

         for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         {
            k = g->head[e];
            if( Is_term(g->term[k]) && k != root )
               break;
         }
         assert(e != EAT_LAST);
         assert(k >= 0);

         if( g->term2edge[i] != e )
         {
            SCIPdebugMessage("term2edge consistency fail3 %d \n", i);
            return FALSE;
         }

         if( g->term2edge[k] != flipedge(e) )
         {
            SCIPdebugMessage("term2edge consistency fail4 %d \n", i);
            return FALSE;
         }
      }
   }
   return TRUE;
}


/** transformed problem consistent to original one? Call only for extended graph */
SCIP_Bool graph_pc_transOrgAreConistent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< the graph */
   SCIP_Bool             verbose             /**< be verbose? */
   )
{
   const int nedges = graph->edges;
   const SCIP_Real* const cost = graph->cost;
   const SCIP_Real* const cost_org = graph->cost_org_pc;

   assert(graph->cost_org_pc && graph->prize);
   assert(graph->extended);
   assert(graph_pc_isPcMw(graph));

   for( int e = 0; e < nedges; e++ )
   {
      if( graph->oeat[e] != EAT_FREE )
      {
         const int head = graph->head[e];

         if( Is_nonleafTerm(graph->term[head]) )
         {
            const SCIP_Real prize = graph->prize[head];
            assert(prize > 0.0);

            if( !SCIPisEQ(scip, cost_org[e], cost[e] + prize) )
            {
               if( verbose )
               {
                  graph_edge_printInfo(graph, e);
                  printf("cost_org=%f cost=%f prize=%f \n", cost_org[e], cost[e], prize);
               }

               return FALSE;
            }
         }
         else
         {
            if( !SCIPisEQ(scip, cost_org[e], cost[e]) )
            {
               if( verbose )
               {
                  graph_edge_printInfo(graph, e);
                  printf("cost_org=%f cost=%f \n", cost_org[e], cost[e]);
               }

               return FALSE;
            }
         }
      }
   }

   return TRUE;
}


/** change property of node to be a non-terminal; prize is not changed! */
void graph_pc_knotToNonTerm(
   GRAPH*                g,                  /**< the graph */
   int                   node                /**< node to be changed */
   )
{
   assert(g);
   assert(node >= 0 && node < g->knots);

   assert(graph_pc_isPcMw(g));
   assert(g->term2edge);

   g->term2edge[node] = TERM2EDGE_NOTERM;

   graph_knot_chg(g, node, STP_TERM_NONE);
}


/** change property of node to be a terminal; prize is changed! */
void graph_pc_knotToFixedTerm(
   GRAPH*                g,                  /**< the graph */
   int                   node                /**< node to be changed */
   )
{
   assert(g);
   assert(node >= 0 && node < g->knots);
   assert(g->term2edge && g->prize);

   g->term2edge[node] = TERM2EDGE_FIXEDTERM;
   g->prize[node] = FARAWAY;

   graph_knot_chg(g, node, STP_TERM);
}


/** check whether node is fixed terminal */
SCIP_Bool graph_pc_knotIsFixedTerm(
   const GRAPH*          g,                  /**< the graph */
   int                   node                /**< node to be checked */
   )
{
   assert(g);
   assert(node >= 0 && node < g->knots);
   assert(graph_pc_isPcMw(g));
   assert(g->term2edge);

#ifndef NDEBUG
   if( TERM2EDGE_FIXEDTERM == g->term2edge[node] )
   {
      assert(Is_term(g->term[node]));
      assert(node == g->source || g->prize[node] == FARAWAY);
   }
#endif

   return (TERM2EDGE_FIXEDTERM == g->term2edge[node]);
}


/** check whether node is a dummy (pseudo) terminal */
SCIP_Bool graph_pc_knotIsDummyTerm(
   const GRAPH*          g,                  /**< the graph */
   int                   node                /**< node to be checked */
   )
{
   assert(g      != NULL);
   assert(node   >= 0);
   assert(node   < g->knots);
   assert(graph_pc_isPcMw(g));
   assert(graph_pc_knotIsFixedTerm(g, g->source));

   if( node == g->source && !graph_pc_isRootedPcMw(g) )
      return TRUE;

   if( g->extended )
   {
      if( Is_term(g->term[node]) && !graph_pc_knotIsFixedTerm(g, node) )
      {
         assert(g->grad[node] == 2 );

         return TRUE;
      }
   }
   else
   {
      if( Is_pseudoTerm(g->term[node]) )
      {
         assert(g->grad[node] == 2 );
         assert(!graph_pc_knotIsFixedTerm(g, node));

         return TRUE;
      }
   }

   return FALSE;
}

/** check whether terminal is not a leaf in at least one optimal tree */
void graph_pc_termMarkProper(
   const GRAPH*          g,                  /**< the graph */
   int*                  termmark            /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
)
{
   const int nnodes = g->knots;

   assert(!g->extended);

   assert(g && termmark);

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) )
      {
         if( graph_pc_termIsNonLeafTerm(g, i) )
            termmark[i] = 1;
         else
            termmark[i] = 2;
      }
      else
      {
         termmark[i] = 0;
      }
   }
}


/** check whether node is a terminal AND is not a leaf (or not contained) in at least one optimal tree */
SCIP_Bool graph_pc_knotIsNonLeafTerm(
   const GRAPH*          g,                  /**< the graph */
   int                   node                /**< node to be checked */
   )
{
   assert(g);
   assert(node >= 0 && node < g->knots);

   if( !Is_anyTerm(g->term[node]) )
      return FALSE;

   return graph_pc_termIsNonLeafTerm(g, node);
}


/** check whether terminal is not a leaf (or not contained) in at least one optimal tree */
SCIP_Bool graph_pc_termIsNonLeafTerm(
   const GRAPH*          g,                  /**< the graph */
   int                   term                /**< terminal to be checked */
   )
{
   SCIP_Bool isNonLeafTerm = FALSE;

   assert(g && g->term2edge);
   assert(term >= 0 && term < g->knots);
   assert(Is_anyTerm(g->term[term]));
   assert(graph_pc_isPc(g));

   if( graph_pc_knotIsFixedTerm(g, term) )
   {
      isNonLeafTerm = FALSE;
   }
   else if( g->extended )
   {
      isNonLeafTerm = Is_nonleafTerm(g->term[term]);
   }
   else
   {
      /* original graph: */
      assert(Is_term(g->term[term]) || Is_pseudoTerm(g->term[term]));

      isNonLeafTerm = (g->term2edge[term] == TERM2EDGE_NONLEAFTERM);

      assert(!(Is_pseudoTerm(g->term[term]) && isNonLeafTerm));
   }

   assert(!isNonLeafTerm || g->term2edge[term] == TERM2EDGE_NONLEAFTERM);

   return isNonLeafTerm;
}


/** set high costs for not including given pseudo-terminal */
void graph_pc_enforcePseudoTerm(
   GRAPH*          graph,              /**< graph */
   int             pterm               /**< the pseudo-terminal */
)
{
   const int term = graph_pc_getTwinTerm(graph, pterm);
   const int root2term = graph_pc_getRoot2PtermEdge(graph, term);

   assert(graph != NULL && Is_pseudoTerm(graph->term[pterm]));
   assert(Is_term(graph->term[term]));
   assert(graph->cost[root2term] == graph->prize[pterm]);

   /* don't change because of weird prize sum in reduce.c */
   if( graph->prize[pterm] < BLOCKED )
   {
      graph->prize[pterm] = BLOCKED - 1.0;
      graph->cost[root2term] = BLOCKED - 1.0;
   }
}


/** enforces non-leaf terminal */
void graph_pc_enforceNonLeafTerm(
   GRAPH*          graph,              /**< graph */
   int             nonleafterm         /**< the terminal */
)
{
   const SCIP_Bool prize = graph->prize[nonleafterm];

   assert(graph_pc_isPcMw(graph) && graph->extended);
   assert(graph_pc_knotIsNonLeafTerm(graph, nonleafterm));
   assert(prize > 0.0 && prize < FARAWAY);

   for( int e = graph->inpbeg[nonleafterm]; e != EAT_LAST; e = graph->ieat[e] )
   {
      graph->cost[e] += prize;
   }

   if( graph_pc_isRootedPcMw(graph) )
   {
      /* make it a proper fixed terminal */
      graph_pc_knotToFixedTerm(graph, nonleafterm);
   }
   else if( graph->prize[nonleafterm] < BLOCKED )
   {
      /* don't change because of weird prize sum in reduce.c */
      graph->prize[nonleafterm] = BLOCKED - 1.0; // todo quite hacky, because it destroys the invariant of non-leaf terms!
   }
}


/** updates term2edge array for new graph that is a subgraph of an old graph
 *  needs to be called right before corresponding edge is added */
void graph_pc_updateTerm2edge(
   GRAPH*                newgraph,           /**< the new graph */
   const GRAPH*          oldgraph,           /**< the old graph */
   int                   newtail,            /**< tail in new graph */
   int                   newhead,            /**< head in new graph */
   int                   oldtail,            /**< tail in old graph */
   int                   oldhead             /**< head in old graph */
)
{
   assert(newgraph != NULL);
   assert(oldgraph != NULL);
   assert(newgraph->term2edge != NULL);
   assert(oldgraph->term2edge != NULL);
   assert(newtail >= 0);
   assert(newhead >= 0);
   assert(oldtail >= 0);
   assert(oldhead >= 0);
   assert(oldgraph->extended);
   assert(newgraph->extended);

   if( oldgraph->term2edge[oldtail] >= 0 && oldgraph->term2edge[oldhead] >= 0 && oldgraph->term[oldtail] != oldgraph->term[oldhead] )
   {
      assert(Is_anyTerm(newgraph->term[newtail]) && Is_anyTerm(newgraph->term[newhead]));
      assert(Is_anyTerm(oldgraph->term[oldtail]) && Is_anyTerm(oldgraph->term[oldhead]));
      assert(oldgraph->source != oldtail && oldgraph->source != oldhead);
      assert(flipedge(newgraph->edges) == newgraph->edges + 1);

      newgraph->term2edge[newtail] = newgraph->edges;
      newgraph->term2edge[newhead] = newgraph->edges + 1;
   }

#ifndef NDEBUG
   if( TERM2EDGE_NOTERM == oldgraph->term2edge[oldtail] )
      assert(TERM2EDGE_NOTERM == newgraph->term2edge[newtail]);

   if( TERM2EDGE_NOTERM == oldgraph->term2edge[oldhead] )
      assert(TERM2EDGE_NOTERM == newgraph->term2edge[newhead]);
#endif

   /* now save the terminal states if there are any */

   if( oldgraph->term2edge[oldtail] < 0 )
      newgraph->term2edge[newtail] = oldgraph->term2edge[oldtail];

   if( oldgraph->term2edge[oldhead] < 0 )
      newgraph->term2edge[newhead] = oldgraph->term2edge[oldhead];
}


/** mark original graph (without dummy terminals) */
void graph_pc_markOrgGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
)
{
   const int root = g->source;
   const int nnodes = g->knots;

   assert(g != NULL);
   assert(graph_pc_isPcMw(g));
   assert(g->extended);

   for( int k = 0; k < nnodes; k++ )
      g->mark[k] = (g->grad[k] > 0);

   for( int e = g->outbeg[root]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int head = g->head[e];

      if( Is_term(g->term[head]) && !graph_pc_knotIsFixedTerm(g, head) )
      {
         g->mark[head] = FALSE;
         assert(g->grad[head] == 2);
         assert(SCIPisGT(scip, g->cost[e], 0.0));
      }
   }

   if( !graph_pc_isRootedPcMw(g) )
      g->mark[root] = FALSE;
}


/** mark terminals and switch terminal property to original terminals */
void graph_pc_2org(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   const int root = graph->source;
   const int nnodes = graph->knots;

   assert(scip && graph);
   assert(graph->extended && graph_pc_isPcMw(graph));

   /* restore original edge weights */
   if( graph_pc_isPc(graph) )
   {
      assert(graph_pc_transOrgAreConistent(scip, graph, FALSE));

      BMScopyMemoryArray(graph->cost, graph->cost_org_pc, graph->edges);
   }

   /* swap terminal properties and mark original graph */
   for( int k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);

      if( Is_pseudoTerm(graph->term[k]) || Is_nonleafTerm(graph->term[k]) )
      {
         graph_knot_chg(graph, k, 0);
      }
      else if( Is_term(graph->term[k]) && !graph_pc_knotIsFixedTerm(graph, k) )
      {
         assert(k != root);

         graph->mark[k] = FALSE;
         graph_knot_chg(graph, k, -2);
      }
   }

   if( graph->stp_type == STP_RPCSPG || graph->stp_type == STP_RMWCSP )
      graph->mark[root] = TRUE;
   else
      graph->mark[root] = FALSE;

   graph->extended = FALSE;
}

/** mark transformed graph and adapt terminal properties to transformed graph */
void graph_pc_2trans(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   const int nnodes = graph->knots;;

   assert(scip);
   assert(!graph->extended);
   assert(graph_pc_isPcMw(graph));

   /* adapt terminal properties for non-leaf terminals (so far not explicitly marked) */
   if( graph_pc_isPc(graph) )
      markNonLeafTerms_2trans(scip, graph);

   /* adapt terminal properties and mark transformed graph */
   for( int k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);

      if( Is_pseudoTerm(graph->term[k]) )
      {
         graph_knot_chg(graph, k, 0);
      }
      else if( Is_term(graph->term[k]) && !graph_pc_knotIsFixedTerm(graph, k) )
      {
         assert(k != graph->source);
         graph_knot_chg(graph, k, -2);
      }
   }

   graph->extended = TRUE;

   /* restore transformed edge weights (shift) and store original ones */
   if( graph_pc_isPc(graph) )
   {
      assert(graph->cost_org_pc);
      BMScopyMemoryArray(graph->cost_org_pc, graph->cost, graph->edges);

      shiftNonLeafCosts_2trans(scip, graph);
   }
}

/** graph_pc_2org if extended */
void graph_pc_2orgcheck(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   assert(graph && scip);

   if( !graph->extended )
      return;

   graph_pc_2org(scip, graph);
}

/** graph_pc_2trans if not extended */
void graph_pc_2transcheck(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   assert(graph && scip);

   if( graph->extended )
      return;

   graph_pc_2trans(scip, graph);
}

/* returns sum of positive vertex weights */
SCIP_Real graph_pc_getPosPrizeSum(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph               /**< the graph */
   )
{
   SCIP_Real prizesum = 0.0;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->prize != NULL);
   assert(!graph->extended);

   for( int i = 0; i < graph->knots; i++ )
      if( Is_term(graph->term[i]) && i != graph->source && graph->prize[i] < BLOCKED )
         prizesum += graph->prize[i];

   return prizesum;
}


/* returns degree of non-root vertex in non-extended graph */
int graph_pc_realDegree(
   const GRAPH*          g,                  /**< graph data structure */
   int                   i,                  /**< the vertex to be checked */
   SCIP_Bool             fixedterm           /**< fixed terminal? */
)
{
   const int ggrad = g->grad[i];
   const SCIP_Bool rpc = (g->stp_type == STP_RPCSPG);
   int newgrad;

   assert(g != NULL);
   assert(!g->extended);
   assert(!Is_pseudoTerm(g->term[i]));
   assert(i != g->source);
   assert(rpc || g->stp_type == STP_PCSPG);

   if( !Is_term(g->term[i]) || fixedterm )
      newgrad = ggrad;
   else if( rpc )
      newgrad = ggrad - 1;
   else
      newgrad = ggrad - 2;

   return newgrad;
}


/** obtains an SAP from prize collecting problems */
SCIP_RETCODE graph_pc_getSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Real*            offset              /**< offset */
   )
{
   SCIP_Real prizesum = 0.0;
   int e;
   int maxpvert;
   const int root = graph->source;
   const int nnodes = graph->knots;
   const int nterms = graph->terms;
   const int stp_type = graph->stp_type;
   int pseudoroot;

   assert(scip && graph && graph->prize);
   assert(graph->knots == graph->ksize);
   assert(graph->edges == graph->esize);
   assert(graph->extended);

   graph->stp_type = STP_SAP;
   SCIP_CALL( graph_copy(scip, graph, newgraph) );
   graph->stp_type = stp_type;

   /* for each terminal, except for the root, three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, (*newgraph), ((*newgraph)->ksize + 1), ((*newgraph)->esize + 2 * (nterms - 1)) , -1) );

   assert((*newgraph)->source == root);

   /* new pseudo-root */
   pseudoroot = (*newgraph)->knots;
   graph_knot_add((*newgraph), -1);

   maxpvert = -1;

   for( int k = 0; k < nnodes; k++ )
      if( Is_pseudoTerm(graph->term[k]) && (maxpvert == -1 || graph->prize[k] > graph->prize[maxpvert]) )
         maxpvert = k;

   /* compute upper bound on best prize sum */
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_pseudoTerm(graph->term[k]) )
      {
         prizesum += graph->prize[k];

         if( stp_type == STP_PCSPG && k != maxpvert )
         {
            SCIP_Real minin = FARAWAY;
            for( e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
               if( !graph_pc_knotIsDummyTerm(graph, graph->tail[e]) && graph->cost[e] < minin )
                  minin = graph->cost[e];

            assert(!SCIPisZero(scip, minin));

            prizesum -= MIN(minin, graph->prize[k]);
         }
      }
   }

   if( stp_type != STP_PCSPG && maxpvert >= 0 )
      prizesum -= graph->prize[maxpvert];

   assert(SCIPisLT(scip, prizesum, FARAWAY));

   *offset -= prizesum;

   SCIP_CALL( graph_pc_presolInit(scip, *newgraph) );

   e = (*newgraph)->outbeg[root];

   while( e != EAT_LAST )
   {
      const int enext = (*newgraph)->oeat[e];
      const int head = (*newgraph)->head[e];

      if( Is_term((*newgraph)->term[head]) )
      {
         (void) graph_edge_redirect(scip, (*newgraph), e, pseudoroot, head, graph->cost[e], TRUE, FALSE);
         (*newgraph)->cost[flipedge(e)] = FARAWAY;
         assert((*newgraph)->head[e] == head);
         assert((*newgraph)->tail[e] == pseudoroot);
      }
      else
      {
         (*newgraph)->cost[e] = prizesum;
      }

      e = enext;
   }

   graph_pc_presolExit(scip, *newgraph);

   for( int k = 0; k < nnodes; k++ )
   {
      /* is the kth node a terminal other than the root? */
      if( Is_pseudoTerm((*newgraph)->term[k]) )
         graph_edge_add(scip, (*newgraph), k, pseudoroot, 0.0, FARAWAY);
   }

   return SCIP_OKAY;
}


/** builds new rooted SAP graph for prize-collecting problems (with given root for SAP) */
SCIP_RETCODE graph_pc_getRsap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   const int*            rootcands,          /**< array containing all vertices that could be used as root */
   int                   nrootcands,         /**< number of all vertices that could be used as root */
   int                   saproot             /**< the root of the new SAP */
   )
{
   GRAPH* p;
   int twinterm;
   int nnodes;
   const int oldroot = graph->source;
   const int stp_type = graph->stp_type;

   assert(scip && graph && graph->prize && rootcands);
   assert(graph->knots == graph->ksize);
   assert(graph->edges == graph->esize);
   assert(graph_pc_isPcMw(graph) && !graph_pc_isRootedPcMw(graph));

   graph_pc_2transcheck(scip, graph);

   assert(Is_pseudoTerm(graph->term[saproot]));

   /* copy graph to obtain an SAP */
   graph->stp_type = STP_SAP;
   SCIP_CALL( graph_copy(scip, graph, newgraph) );
   graph->stp_type = stp_type;

   p = *newgraph;
   twinterm = -1;

   for( int e = p->outbeg[saproot]; e != EAT_LAST; e = p->oeat[e] )
   {
      const int head = p->head[e];
      if( Is_term(p->term[head]) && head != oldroot )
      {
         graph_knot_chg(p, head, STP_TERM_NONE);
         twinterm = head;
         graph_edge_del(scip, p, e, FALSE);
         break;
      }
   }

   assert(twinterm >= 0);

   SCIP_CALL( graph_pc_presolInit(scip, p) );

   for( int e = graph->outbeg[oldroot]; e != EAT_LAST; e = graph->oeat[e] )
   {
      const int head = graph->head[e];

      assert(graph->head[e] == p->head[e]);
      assert(graph->tail[e] == p->tail[e]);

      if( Is_term(graph->term[head]) && head != twinterm )
      {
         assert(Is_term(p->term[head]));

         (void) graph_edge_redirect(scip, p, e, saproot, head, graph->cost[e], TRUE, FALSE);
         p->cost[flipedge(e)] = FARAWAY;

#ifndef NDEBUG
         assert(p->grad[head] == 2);
         for( int e2 = p->outbeg[head]; e2 != EAT_LAST; e2 = p->oeat[e2] )
            assert(p->head[e2] == saproot || Is_pseudoTerm(p->term[p->head[e2]]));
#endif
      }
      else
      {
         graph_edge_del(scip, p, e, FALSE);
      }
   }

   assert(p->grad[twinterm] == 0 && p->grad[oldroot] == 0);

   graph_pc_presolExit(scip, p);

   nnodes = p->knots;
   p->source = saproot;
   graph_knot_chg(p, saproot, STP_TERM);

   for( int k = 0; k < nnodes; k++ )
      p->mark[k] = (p->grad[k] > 0);

   SCIP_CALL( graph_pc_init(scip, p, nnodes, nnodes) );

   for( int k = 0; k < nnodes; k++)
   {
      p->term2edge[k] = graph->term2edge[k];
      if( k < graph->norgmodelknots )
         p->prize[k] = graph->prize[k];
      else
         p->prize[k] = 0.0;
   }

   p->term2edge[saproot] = TERM2EDGE_FIXEDTERM;
   p->term2edge[twinterm] = TERM2EDGE_NOTERM;

   if( nrootcands > 0 )
   {
      SCIP_CALL( graph_pc_presolInit(scip, p) );
      for( int k = 0; k < nrootcands; k++ )
      {
         int e;
         int head = -1;
         const int rootcand = rootcands[k];

         if( rootcand == saproot )
            continue;

         assert(Is_pseudoTerm(p->term[rootcand]));

         for( e = p->outbeg[rootcand]; e != EAT_LAST; e = p->oeat[e] )
         {
            head = p->head[e];

            if( Is_term(p->term[head]) && p->term2edge[head] >= 0 )
            {
               assert(p->grad[head] == 2 && head != saproot);

               graph_knot_chg(p, head, STP_TERM_NONE);
               p->term2edge[head] = TERM2EDGE_NOTERM;
               graph_knot_del(scip, p, head, FALSE);
               break;
            }
         }
         assert(e != EAT_LAST && head >= 0);

         graph_pc_knotToFixedTerm(p, rootcand);
      }
      graph_pc_presolExit(scip, p);
   }

   graph_knot_chg(p, oldroot, STP_TERM_NONE);
   p->prize[saproot] = 0.0;

   return SCIP_OKAY;
}


/** adapts SAP deriving from PCST or MWCS problem with new big M */
void graph_pc_adaptSap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bigM,               /**< new big M value */
   GRAPH*                graph,              /**< the SAP graph */
   SCIP_Real*            offset              /**< the offset */
   )
{
   SCIP_Real oldbigM;
   const int root = graph->source;

   assert(bigM > 0.0);
   assert(scip != NULL && graph != NULL && offset != NULL);
   assert(graph->outbeg[root] >= 0);

   oldbigM = graph->cost[graph->outbeg[root]];
   assert(oldbigM > 0.0);

   *offset += (oldbigM - bigM);

   SCIPdebugMessage("new vs old %f, %f \n", bigM, oldbigM);

   for( int e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
   {
      assert(graph->cost[e] == oldbigM);
      graph->cost[e] = bigM;
   }
}


/** initializes biased data structure */
void graph_pc_getBiased(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Bool             fullbias,           /**< fully bias or only dominated terminals? */
   SCIP_Real*            costbiased,         /**< biased costs */
   SCIP_Real*            prizebiased         /**< biased prizes */
)
{
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(graph);
   const int root = graph->source;

   assert(scip && graph && costbiased && prizebiased);
   assert(graph->extended);
   assert(!rpcmw || graph_pc_knotIsFixedTerm(graph, root));

   BMScopyMemoryArray(costbiased, graph->cost, nedges);
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_pseudoTerm(graph->term[k]) && graph->grad[k] != 0 )
      {
         SCIP_Real mincost = FARAWAY;

         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            const int head = graph->head[e];
            if( Is_term(graph->term[head]) && !graph_pc_knotIsFixedTerm(graph, head) )
               continue;

            if( graph->cost[e] < mincost )
               mincost = graph->cost[e];
         }

#ifndef NDEBUG
         {
            SCIP_Real mincostt = FARAWAY;

            for( int e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
            {
               if( !rpcmw && graph->tail[e] == root )
                  continue;

               if( graph->cost[e] < mincostt )
                  mincostt = graph->cost[e];
            }
            assert((graph->stp_type != STP_PCSPG && graph->stp_type != STP_RPCSPG )
                  || SCIPisEQ(scip, mincostt, mincost) || graph->grad[k] <= 2 );
         }
#endif

         if( fullbias )
         {
            mincost = MIN(mincost, graph->prize[k]);
         }
         else
         {
            if( mincost < graph->prize[k] )
            {
               prizebiased[k] = graph->prize[k];
               continue;
            }

            mincost = graph->prize[k];
         }

         for( int e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
         {
            if( !rpcmw && graph->tail[e] == root )
               continue;

            if( SCIPisGE(scip, graph->cost[e], FARAWAY) )
               continue;

            costbiased[e] -= mincost;
            assert(!SCIPisNegative(scip, costbiased[e]) || (graph->stp_type != STP_PCSPG && graph->stp_type != STP_RPCSPG));
         }

         prizebiased[k] = graph->prize[k] - mincost;
         assert(!SCIPisNegative(scip, prizebiased[k]));
      }
      else
      {
         if( rpcmw && graph_pc_knotIsFixedTerm(graph, k) )
            prizebiased[k] = graph->prize[k];
         else
            prizebiased[k] = 0.0;
      }
   }
}


/** alters the graph for prize collecting problems */
SCIP_RETCODE graph_pc_2pc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   int root;
   const int nnodes = graph->knots;
   int termscount;

   assert(scip && graph->prize);
   assert(!graph->extended);
   assert(graph->edges == graph->esize && nnodes == graph->ksize);

   graph->norgmodeledges = graph->edges;
   graph->norgmodelknots = nnodes;

   /* for PC remove terminal property from non-leaves and reduce terminal count */
   if( graph->stp_type != STP_MWCSP )
      markNonLeafTerms_pretrans(scip, graph);

   /* for each proper terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + graph->terms + 1), (graph->esize + graph->terms * 6) , -1) );

   /* add future terminals */
   for( int k = 0; k < graph->terms; ++k )
      graph_knot_add(graph, -1);

   /* add new root */
   root = graph->knots;
   graph_knot_add(graph, 0);
   graph->prize[root] = 0.0;

   graph_pc_init(scip, graph, -1, graph->knots);
   assert(graph->term2edge);

   graph->term2edge[root] = TERM2EDGE_FIXEDTERM;

   termscount = 0;
   for( int k = 0; k < nnodes; ++k )
   {
      if( Is_nonleafTerm(graph->term[k]) )
      {
         assert(graph->stp_type != STP_MWCSP);
         graph->term2edge[k] = TERM2EDGE_NONLEAFTERM;

         continue;
      }

      if( Is_term(graph->term[k]) )
      {
         /* get the new terminal */
         const int node = nnodes + termscount;
         termscount++;

         assert(node != root && k != root);

         /* switch the terminal property, mark k */
         graph_knot_chg(graph, k, STP_TERM_PSEUDO);
         graph_knot_chg(graph, node, STP_TERM);
         graph->prize[node] = 0.0;
         assert(SCIPisGE(scip, graph->prize[k], 0.0));

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, k, 0.0, FARAWAY);
         graph_edge_add(scip, graph, root, node, graph->prize[k], FARAWAY);

         graph->term2edge[k] = graph->edges;
         graph->term2edge[node] = graph->edges + 1;
         assert(graph->edges + 1 == flipedge(graph->edges));

         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);

         assert(graph->head[graph->term2edge[k]] == node);
         assert(graph->head[graph->term2edge[node]] == k);
      }
      else if( graph->stp_type != STP_MWCSP )
      {
         graph->prize[k] = 0.0;
      }
   }

   assert((termscount + 1) == graph->terms);

   graph->source = root;
   graph->extended = TRUE;

   if( graph->stp_type != STP_MWCSP )
   {
      graph->stp_type = STP_PCSPG;

      SCIP_CALL( initCostOrgPc(scip, graph) );
      shiftNonLeafCosts_2trans(scip, graph);
      SCIPdebugMessage("Transformed to PC \n");
   }

   assert(graph_pc_term2edgeIsConsistent(scip, graph));
   assert(graph_valid(scip, graph));
   assert(graph->orgsource == -1);

   return SCIP_OKAY;
}


/** alters the graph for rooted prize collecting problems */
SCIP_RETCODE graph_pc_2rpc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   const int root = graph->source;
   const int nnodes = graph->knots;
   int nterms;
   int nfixterms;
   int npotterms;

   assert(graph && graph->prize);
   assert(graph->edges == graph->esize);
   assert(graph->knots == graph->ksize);
   assert(!graph->extended);
   assert(root >= 0 && root < graph->knots) ;

   graph->norgmodeledges = graph->edges;
   graph->norgmodelknots = nnodes;
   nfixterms = 0;
   npotterms = 0;

   /* remove terminal property from non-leaves and reduce terminal count */
   markNonLeafTerms_pretrans(scip, graph);

   /* count number of fixed and potential terminals and adapt prizes */
   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_nonleafTerm(graph->term[i]) )
         continue;

      if( !Is_term(graph->term[i]) )
      {
         graph->prize[i] = 0.0;
         assert(graph->term[i] == STP_TERM_NONE);
      }
      else if( SCIPisGE(scip, graph->prize[i], FARAWAY) )
      {
         assert(SCIPisEQ(scip, graph->prize[i], FARAWAY));
         nfixterms++;
      }
      else if( SCIPisGT(scip, graph->prize[i], 0.0) )
      {
         assert(i != root);
         assert(Is_term(graph->term[i]));
         npotterms++;
      }
      else
      {
         assert(i != root);
         assert(SCIPisEQ(scip, graph->prize[i], 0.0));
         graph->prize[i] = 0.0;
         graph_knot_chg(graph, i, STP_TERM_NONE);
      }
   }

   assert(npotterms + nfixterms == graph->terms);

   /* for each terminal, except for the root, one node and two edges (i.e. four arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + npotterms), (graph->esize + npotterms * 4) , -1) );

   /* create new nodes corresponding to potential terminals */
   for( int k = 0; k < npotterms; ++k )
      graph_knot_add(graph, STP_TERM_NONE);

   /* allocate and default initialize term2edge array */
   graph_pc_init(scip, graph, -1, graph->knots);
   assert(graph->term2edge);

   nterms = 0;

   for( int k = 0; k < nnodes; ++k )
   {
      if( Is_nonleafTerm(graph->term[k]) )
      {
         graph->term2edge[k] = TERM2EDGE_NONLEAFTERM;
         continue;
      }

      /* is the kth node a potential terminal? */
      if( Is_term(graph->term[k]) && SCIPisLT(scip, graph->prize[k], FARAWAY) )
      {
         /* the future terminal */
         const int node = nnodes + nterms;
         nterms++;

         assert(k != root && node != root);

         /* switch the terminal property, mark k as former terminal */
         graph_knot_chg(graph, k, STP_TERM_PSEUDO);
         graph_knot_chg(graph, node, STP_TERM);
         assert(SCIPisGE(scip, graph->prize[k], 0.0));
         graph->prize[node] = 0.0;

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, node, graph->prize[k], FARAWAY);

         graph->term2edge[k] = graph->edges;
         graph->term2edge[node] = graph->edges + 1;
         assert(graph->edges + 1 == flipedge(graph->edges));

         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);

         assert(graph->head[graph->term2edge[k]] == node);
         assert(graph->head[graph->term2edge[node]] == k);
      }
      else
      {
         assert(graph->prize[k] == FARAWAY || graph->prize[k] == 0.0);
      }
   }

   graph->stp_type = STP_RPCSPG;
   graph->orgsource = graph->source;
   graph->extended = TRUE;

   SCIP_CALL( initCostOrgPc(scip, graph) );
   shiftNonLeafCosts_2trans(scip, graph);

   assert(nterms == npotterms);
   assert(graph->prize[graph->source] == FARAWAY);
   SCIPdebugMessage("Transformed problem to (RPC) SAP \n");

   return SCIP_OKAY;
}

/** alters the graph for MWCS problems */
SCIP_RETCODE graph_pc_2mw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real*            maxweights          /**< array containing the weight of each node */
   )
{
   int nnodes;
   int nterms = 0;

   assert(maxweights != NULL);
   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->cost != NULL);
   assert(graph->terms == 0);

   nnodes = graph->knots;

   /* count number of terminals, modify incoming edges for non-terminals */
   for( int i = 0; i < nnodes; i++ )
   {
      if( SCIPisLT(scip, maxweights[i], 0.0) )
      {
         for( int e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
            graph->cost[e] -= maxweights[i];
      }
      else if( SCIPisGT(scip, maxweights[i], 0.0) )
      {
         graph_knot_chg(graph, i, 0);
         nterms++;
      }
   }
   nterms = 0;
   for( int i = 0; i < nnodes; i++ )
   {
      graph->prize[i] = maxweights[i];
      if( Is_term(graph->term[i]) )
      {
         assert(SCIPisGE(scip, maxweights[i], 0.0));
         nterms++;
      }
      else
      {
         assert(SCIPisLE(scip, maxweights[i], 0.0));
      }
   }
   assert(nterms == graph->terms);
   graph->stp_type = STP_MWCSP;

   SCIP_CALL( graph_pc_2pc(scip, graph) );
   assert(graph->stp_type == STP_MWCSP);

   SCIPdebugMessage("Transformed to MW \n");

   return SCIP_OKAY;
}



/** alters the graph for RMWCS problems */
SCIP_RETCODE graph_pc_2rmw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   SCIP_Real* maxweights;
   int i;
   int root;
   int nnodes;
   int npterms;
   int nrterms;
   int maxgrad;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->cost != NULL);

   root = -1;
   maxgrad = -1;
   npterms = 0;
   nrterms = 0;
   nnodes = graph->knots;
   maxweights = graph->prize;

   assert(maxweights != NULL);

   /* count number of terminals, modify incoming edges for non-terminals */
   for( i = 0; i < nnodes; i++ )
   {
      if( SCIPisLT(scip, maxweights[i], 0.0) )
      {
         for( int e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
            graph->cost[e] -= maxweights[i];
      }
      else if( SCIPisGE(scip, maxweights[i], FARAWAY) )
      {
         assert(Is_term(graph->term[i]));
         if( graph->grad[i] > maxgrad )
         {
            root = i;
            maxgrad = graph->grad[i];
         }

         nrterms++;
      }
      else if( SCIPisGT(scip, maxweights[i], 0.0) )
      {
         graph_knot_chg(graph, i, 0);
         npterms++;
      }
   }

   assert(root >= 0);
   assert(graph->terms == (npterms + nrterms));

   graph->norgmodeledges = graph->edges;
   graph->norgmodelknots = nnodes;
   graph->source = root;

   /* for each terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + npterms), (graph->esize + npterms * 4) , -1) );
   maxweights = graph->prize;

   /* create a new nodes */
   for( int k = 0; k < npterms; k++ )
      graph_knot_add(graph, -1);

   /* allocate and initialize term2edge array */
   graph_pc_init(scip, graph, -1, graph->knots);
   assert(graph->term2edge);

   i = 0;
   for( int k = 0; k < nnodes; ++k )
   {
      /* is the kth node a non-fixed terminal */
      if( Is_term(graph->term[k]) && SCIPisLT(scip, maxweights[k], FARAWAY) )
      {
         /* the copied node */
         const int node = nnodes + i;
         i++;

         /* switch the terminal property, mark k */
         graph_knot_chg(graph, k, -2);
         graph_knot_chg(graph, node, 0);
         graph->prize[node] = 0.0;
         assert(SCIPisGE(scip, maxweights[k], 0.0));

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, node, maxweights[k], FARAWAY);

         graph->term2edge[k] = graph->edges;
         graph->term2edge[node] = graph->edges + 1;
         assert(graph->edges + 1 == flipedge(graph->edges));

         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);

         assert(graph->head[graph->term2edge[k]] == node);
         assert(graph->head[graph->term2edge[node]] == k);
      }
   }

   assert(i == npterms);
   graph->extended = TRUE;
   graph->stp_type = STP_RMWCSP;
   graph->orgsource = graph->source;

   SCIPdebugMessage("Transformed to RMW \n");

   return SCIP_OKAY;
}


/** transforms PCSPG or MWCSP to RPCSPG or RMWCSP if possible */
SCIP_RETCODE graph_pc_pcmw2rooted(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real             prizesum            /**< sum of positive prizes */
   )
{
   int e;
   int p;
   int newroot;
   int maxgrad;
   int nfixedterms;
   const int orgnterms = graph->terms;
   const int root = graph->source;
   const int pc = (graph->stp_type == STP_PCSPG);

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->term2edge != NULL);
   assert(graph->extended);
   assert(pc || graph->stp_type == STP_MWCSP);

   newroot = -1;
   maxgrad = -1;

#ifndef WITH_UG
   printf("attempt transformation to rooted problem \n");
#endif

   nfixedterms = 0;
   e = graph->outbeg[root];
   while( e != EAT_LAST )
   {
      const int enext = graph->oeat[e];
      if( SCIPisGE(scip, graph->cost[e], prizesum) )
      {
         int e2;
         const int k = graph->head[e];

         assert(Is_term(graph->term[k]));
         assert(graph->grad[k] == 2);

         graph->cost[e] = FARAWAY;

         for( e2 = graph->outbeg[k]; e2 != EAT_LAST; e2 = graph->oeat[e2] )
            if( graph->head[e2] != root )
               break;

         p = graph->head[e2];
         assert(e2 == graph->term2edge[k]);
         assert(Is_pseudoTerm(graph->term[p]));
         assert(SCIPisGE(scip, graph->prize[p], prizesum));

         graph->prize[p] = FARAWAY;

         /* delete terminal */
         graph_knot_chg(graph, k, -1);
         while( graph->outbeg[k] != EAT_LAST )
            graph_edge_del(scip, graph, graph->outbeg[k], TRUE);

         graph->term2edge[k] = TERM2EDGE_NOTERM;
         graph->term2edge[p] = TERM2EDGE_NOTERM;

         graph_knot_chg(graph, p, 0);

         if( graph->grad[p] > maxgrad )
         {
            newroot = p;
            maxgrad = graph->grad[p];
         }
         nfixedterms++;
      }
      e = enext;
   }

   /* is there a new root? */
   if( newroot >= 0 )
   {
      graph->source = newroot;

      if( graph->rootedgeprevs != NULL )
         graph_pc_presolExit(scip, graph);

      e = graph->outbeg[root];
      while( e != EAT_LAST )
      {
         const int enext = graph->oeat[e];
         const int k = graph->head[e];
         if( Is_term(graph->term[k]) && !SCIPisZero(scip, graph->cost[e]) )
         {
            (void) graph_edge_redirect(scip, graph, e, newroot, k, graph->cost[e], TRUE, TRUE);
            graph->cost[flipedge(e)] = FARAWAY;
         }
         e = enext;
      }

      /* delete old root */
      graph_knot_chg(graph, root, -1);
      while( graph->outbeg[root] != EAT_LAST )
         graph_edge_del(scip, graph, graph->outbeg[root], TRUE);

      if( pc )
      {
         graph->stp_type = STP_RPCSPG;
#ifndef WITH_UG
         printf("...transformed PC to RPC; fixed %d out of %d terminals \n", nfixedterms, orgnterms - 1);
#endif

      }
      else
      {
         graph->stp_type = STP_RMWCSP;
#ifndef WITH_UG
         printf("...transformed MW to RMW; fixed %d out of %d terminals \n", nfixedterms, orgnterms - 1);
#endif
      }

      assert(orgnterms - 1 == graph->terms);
   }

#ifndef WITH_UG
   if( !graph_pc_isRootedPcMw(graph) )
      printf("...failed \n");
#endif

   return SCIP_OKAY;
}

/* deletes dummy terminal to given node and edge from pseudo-root (if existent) */
void graph_pc_deleteTermExtension(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< index of the terminal */
   )
{
   int e;
   int dummyterm;
   const SCIP_Bool has_pseudoroot = !graph_pc_isRootedPcMw(g);

   assert(g != NULL && g->pcancestors != NULL && g->term2edge != NULL);
   assert(graph_pc_isPcMw(g) && !g->extended);
   assert(Is_term(g->term[i]));
   assert(!graph_pc_knotIsFixedTerm(g, i) && !graph_pc_knotIsNonLeafTerm(g, i));
   assert(i != g->source);

   /* get edge from i to its artificial terminal */
   e = g->term2edge[i];
   assert(e >= 0);

   dummyterm = g->head[e];
   assert(dummyterm != g->source && !g->mark[dummyterm]);

#ifndef NDEBUG
   {
      int e2;
      for( e2 = g->outbeg[i]; e2 != EAT_LAST; e2 = g->oeat[e2] )
         if( Is_pseudoTerm(g->term[g->head[e2]]) )
            break;
      assert(e2 == e);
   }
#endif

   /* delete edge and unmark artificial terminal */
   graph_knot_chg(g, dummyterm, STP_TERM_NONE);
   graph_edge_del(scip, g, e, TRUE);
   g->term2edge[dummyterm] = TERM2EDGE_NOTERM;

   /* delete remaining incident edge of artificial terminal */
   e = g->inpbeg[dummyterm];

   assert(e != EAT_LAST);
   assert(g->source == g->tail[e]);
   assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));

   graph_edge_del(scip, g, e, TRUE);

   g->term2edge[i] = TERM2EDGE_NOTERM;
   assert(g->inpbeg[dummyterm] == EAT_LAST && g->grad[dummyterm] == 0);

   if( has_pseudoroot )
   {
      const int edgeRoot2i = graph_pc_getRoot2PtermEdge(g, i);

      assert(g->cost[edgeRoot2i] == 0.0);
      graph_edge_del(scip, g, edgeRoot2i, TRUE);
   }
}


/** delete a terminal for a (rooted) prize-collecting problem */
int graph_pc_deleteTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< index of the terminal */
   )
{
   int e;
   int t;
   int grad = g->grad[i];

   assert(g != NULL && scip != NULL && g->term2edge != NULL);
   assert(graph_pc_isPcMw(g));
   assert(Is_term(g->term[i]));
   assert(g->term2edge[i] >= 0);
   assert(i != g->source);

   t = UNKNOWN;

   /* delete terminal */

   graph_pc_knotToNonTerm(g, i);
   g->mark[i] = FALSE;

   while( (e = g->outbeg[i]) != EAT_LAST )
   {
      const int i1 = g->head[e];

      if( Is_pseudoTerm(g->term[i1]) && g->source != i1 )
         t = g->head[e];
      graph_edge_del(scip, g, e, TRUE);
   }

   assert(g->grad[i] == 0);
   assert(t != UNKNOWN);

   /* delete artificial terminal */

   graph_pc_knotToNonTerm(g, t);
   g->mark[t] = FALSE;
   grad += g->grad[t] - 1;

   graph_knot_del(scip, g, t, TRUE);

   return grad;
}


/** subtract a given sum from the prize of a terminal */
void graph_pc_subtractPrize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   SCIP_Real             cost,               /**< cost to be subtracted */
   int                   i                   /**< the terminal */
   )
{
   int e;
   int j;

   assert(scip != NULL);
   assert(g != NULL);
   assert(!g->extended);
   assert(!graph_pc_knotIsFixedTerm(g, i));
   assert(i != g->source);

   g->prize[i] -= cost;
   for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      if( Is_pseudoTerm(g->term[g->head[e]]) )
         break;

   assert(e != EAT_LAST);

   j = g->head[e];

   assert(j != g->source);
   assert(!g->mark[j]);

   for( e = g->inpbeg[j]; e != EAT_LAST; e = g->ieat[e] )
      if( g->source == g->tail[e] )
         break;

   assert(e != EAT_LAST);
   assert(!g->mark[g->tail[e]] || g->stp_type == STP_RPCSPG);

   g->cost[e] -= cost;

   assert(g->stp_type == STP_MWCSP  || g->stp_type == STP_RMWCSP || SCIPisGE(scip, g->prize[i], 0.0));
   assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));
   assert(SCIPisGE(scip, g->prize[i], 0.0) || g->stp_type == STP_MWCSP);
}

/** change prize of a terminal */
void graph_pc_chgPrize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   SCIP_Real             newprize,           /**< prize to be subtracted */
   int                   i                   /**< the terminal */
   )
{
   int e;
   int j;

   assert(scip != NULL);
   assert(g != NULL);
   assert(newprize > 0.0);

   if( g->stp_type == STP_RPCSPG && i == g->source )
      return;

   g->prize[i] = newprize;
   for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      if( Is_pseudoTerm(g->term[g->head[e]]) )
         break;

   assert(e != EAT_LAST);

   j = g->head[e];

   assert(j != g->source);
   assert(!g->mark[j]);

   for( e = g->inpbeg[j]; e != EAT_LAST; e = g->ieat[e] )
      if( g->source == g->tail[e] )
         break;

   assert(e != EAT_LAST);
   assert(!g->mark[g->tail[e]] || g->stp_type == STP_RPCSPG);

   g->cost[e] = newprize;

   assert(g->stp_type == STP_MWCSP  || g->stp_type == STP_RMWCSP || SCIPisGE(scip, g->prize[i], 0.0));
   assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));
   assert(SCIPisGE(scip, g->prize[i], 0.0) || g->stp_type == STP_MWCSP);
}


/** contract ancestors of an edge of (rooted) prize-collecting Steiner tree problem or maximum-weight connected subgraph problem */
SCIP_RETCODE graph_pc_contractNodeAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   t,                  /**< tail node to be contracted (surviving node) */
   int                   s,                  /**< head node to be contracted */
   int                   ets                 /**< edge from t to s or -1 */
   )
{
   assert(g != NULL);
   assert(scip != NULL);

   if( ets < 0 )
   {
      assert(ets == -1);

      for( ets = g->outbeg[t]; ets != EAT_LAST; ets = g->oeat[ets] )
         if( g->head[ets] == s )
            break;
   }

   assert(ets >= 0 && ets < g->edges);

   SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[t]), g->pcancestors[s], NULL) );
   SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[s]), g->pcancestors[t], NULL) );

   SCIP_CALL(SCIPintListNodeAppendCopy(scip, &(g->pcancestors[s]), g->ancestors[ets], NULL));
   SCIP_CALL(SCIPintListNodeAppendCopy(scip, &(g->pcancestors[t]), g->ancestors[ets], NULL));

#if 0
   SCIP_Bool conflict;

   SCIP_CALL( graph_pseudoAncestors_appendCopyEdgeToNode(scip, t, ets, FALSE, g, &conflict) );
   assert(!conflict);

   SCIP_CALL( graph_pseudoAncestors_appendCopyNode(scip, t, s, FALSE, g, &conflict) );
   assert(!conflict);
#endif

   return SCIP_OKAY;
}


/** contract an edge of (rooted) prize-collecting Steiner tree problem or maximum-weight connected subgraph problem */
SCIP_RETCODE graph_pc_contractEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int*                  solnode,            /**< solution nodes or NULL */
   int                   t,                  /**< tail node to be contracted (surviving node) */
   int                   s,                  /**< head node to be contracted */
   int                   term4offset         /**< terminal to add offset to */
   )
{
   int ets;

   assert(scip && g);
   assert(Is_term(g->term[term4offset]));
   assert(!g->extended);
   assert(g->source != s);

   /* get edge from t to s */
   for( ets = g->outbeg[t]; ets != EAT_LAST; ets = g->oeat[ets] )
      if( g->head[ets] == s )
         break;

   assert(ets != EAT_LAST);

   if( graph_pc_knotIsFixedTerm(g, s) || graph_pc_knotIsFixedTerm(g, t) )
   {
      SCIP_CALL( contractEdgeWithFixedEnd(scip, g, solnode, t, s, ets) );
   }
   else
   {
      assert(!graph_pc_isRootedPcMw(g) || (s != g->source && t != g->source));
      SCIP_CALL( contractEdgeNoFixedEnd(scip, g, solnode, t, s, ets, term4offset) );
   }

   assert(g->grad[s] == 0);
   assert(TERM2EDGE_NOTERM == g->term2edge[s]);

   SCIPdebugMessage("PcMw contraction: %d into %d \n", s, t);

   return SCIP_OKAY;
}


/** mark original solution */
SCIP_RETCODE graph_sol_markPcancestors(
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


/** is this graph a prize-collecting variant? */
SCIP_Bool graph_pc_isPc(
   const GRAPH*          g                   /**< the graph */
)
{
   const int type = g->stp_type;
   assert(g != NULL);

   return (type == STP_PCSPG || type == STP_RPCSPG);
}


/** is this graph a maximum-weight variant? */
SCIP_Bool graph_pc_isMw(
   const GRAPH*          g                   /**< the graph */
)
{
   const int type = g->stp_type;
   assert(g != NULL);

   return (type == STP_MWCSP || type == STP_RMWCSP ||type == STP_BRMWCSP);
}


/** is this graph a prize-collecting or maximum-weight variant? */
SCIP_Bool graph_pc_isPcMw(
   const GRAPH*          g                   /**< the graph */
)
{
   const int type = g->stp_type;
   assert(g != NULL);

   return (type == STP_PCSPG || type == STP_RPCSPG || type == STP_MWCSP || type == STP_RMWCSP || type == STP_BRMWCSP);
}


/** get edge from root to (pseudo) terminal */
int graph_pc_getRoot2PtermEdge(
   const GRAPH*          graph,               /**< the graph */
   int                   pseudoterm           /**< the pseudo terminal  */
)
{
   int rootedge = -1;
   const int root = graph->source;

   assert(graph != NULL);
   assert(pseudoterm >= 0 && pseudoterm < graph->knots);

   for( int e = graph->inpbeg[pseudoterm]; e != EAT_LAST; e = graph->ieat[e] )
   {
      if( graph->tail[e] == root )
      {
         rootedge = e;
         break;
      }
   }

   assert(rootedge >= 0);

   return rootedge;
}


/** get number of fixed terminals */
int graph_pc_nFixedTerms(
   const GRAPH*          graph                /**< the graph */
)
{
   int nfixterms = 0;
   const int nnodes = graph->knots;
   assert(graph != NULL);
   assert(graph_pc_isPcMw(graph));

   if( !graph_pc_isRootedPcMw(graph) )
      return 0;

   for( int k = 0; k < nnodes; k++ )
      if( graph_pc_knotIsFixedTerm(graph, k) )
         nfixterms++;

   return nfixterms;
}

/** get number of potential terminals */
int graph_pc_nPotentialTerms(
   const GRAPH*          graph                /**< the graph */
)
{
   assert(graph != NULL);
   assert(graph_pc_isPcMw(graph));

   if( !graph_pc_isRootedPcMw(graph) )
      return (graph->terms - 1);

   return (graph->terms - graph_pc_nFixedTerms(graph));
}

/** get twin-terminal */
int graph_pc_getTwinTerm(
   const GRAPH*          g,                  /**< the graph */
   int                   vertex              /**< the vertex  */
)
{
   assert(g != NULL);
   assert(graph_pc_isPcMw(g));
   assert(g->term2edge != NULL && g->term2edge[vertex] >= 0);
   assert(Is_anyTerm(g->term[vertex]));

   return g->head[g->term2edge[vertex]];
}


/** is this graph a rooted prize-collecting or rooted maximum-weight variant? */
SCIP_Bool graph_pc_isRootedPcMw(
   const GRAPH*          g                   /**< the graph */
)
{
   const int type = g->stp_type;
   assert(g != NULL);

   return (type == STP_RPCSPG || type == STP_RMWCSP || type == STP_BRMWCSP);
}
