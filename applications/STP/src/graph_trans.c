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

/**@file   graph_trans.c
 * @brief  Transformation algorithms for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This includes method for transformations (or reductions) between the different Steiner problem classes.
 *
 * A list of all interface methods can be found in graph.h.
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "portab.h"
#include "graph.h"


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


/** is vertex a non-leaf (call before graph transformation was performed)  */
static inline
SCIP_Bool isNonLeaf_pretransPc(
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
void markNonLeafTerms_pretransPc(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g                   /**< the graph */
)
{
   const int nnodes = g->knots;

   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(g->term[k]) )
         continue;

      if( isNonLeaf_pretransPc(scip, g, k) )
      {
         graph_knot_chg(g, k, STP_TERM_NONLEAF);
      }
   }
}


/*
 * Interface methods
 */


/** alters the graph for prize collecting problems */
SCIP_RETCODE graph_transPc(
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
      markNonLeafTerms_pretransPc(scip, graph);

   /* for each proper terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + graph->terms + 1), (graph->esize + graph->terms * 6) , -1) );

   /* add future terminals */
   for( int k = 0; k < graph->terms; ++k )
      graph_knot_add(graph, -1);

   /* add new root */
   root = graph->knots;
   graph_knot_add(graph, 0);
   graph->prize[root] = 0.0;

   SCIP_CALL( graph_pc_initTerm2Edge(scip, graph, graph->knots) );

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
      graph_pc_shiftNonLeafCosts(scip, graph);
      SCIPdebugMessage("Transformed to PC \n");
   }

   assert(graph_pc_term2edgeIsConsistent(scip, graph));
   assert(graph_valid(scip, graph));
   assert(graph->orgsource == -1);

   return SCIP_OKAY;
}


/** alters the graph for prize collecting problems and transforms it to an SPG */
SCIP_RETCODE graph_transPc2Spg(
   SCIP*                 scip,               /**< SCIP data structure */
   PRESOL*               presol,             /**< presolving struct */
   GRAPH*                graph               /**< the graph */
   )
{
   SCIP_Real baseprize = FARAWAY;
   int root;
   int baseterm = -1;
   int termscount;
   const int nnodes_org = graph->knots;
   SCIP_Real prizesum;
   int naddedarcs;
   int naddednodes;
   const int nterms_org = graph->terms;

   assert(graph && graph->prize);
   assert(graph->edges == graph->esize);
   assert(graph->knots == graph->ksize);
   assert(!graph->extended);
   assert(nterms_org > 0);

   prizesum = 1.0;
   for( int i = 0; i < nnodes_org; i++ )
   {
      if( graph->prize[i] > 0.0 &&  SCIPisLT(scip, graph->prize[i], FARAWAY) )
      {
         prizesum += graph->prize[i];
      }
   }

   printf("prizesum=%f \n", prizesum);
   assert(prizesum < BLOCKED);

   presol->fixed -= prizesum * (nterms_org + 1);

   printf("fixed=%f \n", prizesum * (nterms_org + 1));

   naddednodes = nterms_org + 1;
   naddedarcs = 4 * nterms_org + 4 * (nterms_org - 1);

   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + naddednodes), (graph->esize + naddedarcs) , -1) );

   /* create new nodes corresponding to potential terminals */
   for( int k = 0; k < nterms_org; ++k )
      graph_knot_add(graph, STP_TERM_NONE);

   /* add new root */
   root = graph->knots;
   graph_knot_add(graph, STP_TERM);

   graph->source = root;

   /* select the base term */
   termscount = 0;
   for( int k = 0; k < nnodes_org; ++k )
   {
      if( Is_term(graph->term[k]) )
      {
         if( graph->prize[k] < baseprize )
         {
            baseterm = nnodes_org + termscount;
            baseprize = graph->prize[k];
         }

         termscount++;
      }
   }
   assert(termscount == nterms_org);
   assert(baseterm >= 0);

   termscount = 0;
   for( int k = 0; k < nnodes_org; ++k )
   {
      /* is the kth node a potential terminal? */
      if( Is_term(graph->term[k])  )
      {
         /* the future terminal */
         const int newterm = nnodes_org + termscount;
         termscount++;

         assert(k != root && newterm != root);
         assert(!Is_term(graph->term[newterm]));

         if( baseterm == newterm )
         {
            assert(EQ(graph->prize[k], baseprize));
         }
         else
         {
            graph_edge_addBi(scip, graph, k, baseterm, prizesum + baseprize);
            graph_edge_addBi(scip, graph, newterm, baseterm, prizesum + graph->prize[k]);
         }

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_addBi(scip, graph, root, k, prizesum);
         graph_edge_addBi(scip, graph, k, newterm, prizesum);

         /* switch the terminal property, mark k as former terminal */
         graph_knot_chg(graph, k, STP_TERM_NONE);
         graph_knot_chg(graph, newterm, STP_TERM);
         assert(SCIPisGE(scip, graph->prize[k], 0.0));
      }
   }

   SCIPfreeMemoryArrayNull(scip, &(graph->prize));
   SCIPfreeMemoryArrayNull(scip, &(graph->term2edge));

   assert(termscount == nterms_org);
   graph->stp_type = STP_SPG;
   graph->extended = FALSE;

   return SCIP_OKAY;
}

/** alters the graph for rooted prize collecting problems */
SCIP_RETCODE graph_transRpc(
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
   markNonLeafTerms_pretransPc(scip, graph);

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

   SCIP_CALL( graph_pc_initTerm2Edge(scip, graph, graph->knots) );

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
      else if( Is_term(graph->term[k]) )
      {
         assert(EQ(graph->prize[k], FARAWAY));
         graph_pc_knotToFixedTermProperty(graph, k);
      }
      else
      {
         assert(EQ(graph->prize[k], 0.0));
      }
   }

   graph->stp_type = STP_RPCSPG;
   graph->orgsource = graph->source;
   graph->extended = TRUE;

   SCIP_CALL( initCostOrgPc(scip, graph) );
   graph_pc_shiftNonLeafCosts(scip, graph);

   assert(nterms == npotterms);
   assert(graph->prize[graph->source] == FARAWAY);
   SCIPdebugMessage("Transformed problem to (RPC) SAP \n");

   return SCIP_OKAY;
}


/** changes the graph for rooted prize collecting problems such that no proper potential terminal are fixed */
SCIP_RETCODE graph_transRpc2FixedProper(
   SCIP*                 scip,               /**< SCIP data structure */
   PRESOL*               presol,             /**< presolving struct */
   GRAPH*                graph               /**< the graph */
   )
{
   const int root = graph->source;
   const int nnodes = graph->knots;
   SCIP_Real prizesum;
   int nterms;
   int nfixterms = 0;
   int npropterms = 0;

   assert(graph && graph->prize);
   assert(graph->edges == graph->esize);
   assert(graph->knots == graph->ksize);
   assert(!graph->extended);

   /* count number of fixed and potential terminals */
   for( int i = 0; i < nnodes; i++ )
   {
      if( isNonLeaf_pretransPc(scip, graph, i) )
         continue;

      if( SCIPisGE(scip, graph->prize[i], FARAWAY) )
      {
         assert(SCIPisEQ(scip, graph->prize[i], FARAWAY));
         nfixterms++;
      }
      else if( SCIPisGT(scip, graph->prize[i], 0.0) )
      {
         assert(i != root);
         assert(Is_term(graph->term[i]));
         npropterms++;
      }
   }

   prizesum = 0.0;
   for( int i = 0; i < nnodes; i++ )
   {
      if( graph->prize[i] > 0.0 &&  SCIPisLT(scip, graph->prize[i], FARAWAY) )
      {
         prizesum += graph->prize[i];
      }
   }

   printf("prizesum=%f \n", prizesum);
   assert(prizesum < BLOCKED);

   presol->fixed -= prizesum * npropterms;

   printf("number of proper potential terminals %d \n", npropterms);
   printf("fixed=%f \n", prizesum * npropterms);

   /* for each terminal, except for the root, one node and two edges (i.e. four arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + npropterms), (graph->esize + npropterms * 4) , -1) );

   /* create new nodes corresponding to potential terminals */
   for( int k = 0; k < npropterms; ++k )
      graph_knot_add(graph, STP_TERM_NONE);

   nterms = 0;

   for( int k = 0; k < nnodes; ++k )
   {
      if( isNonLeaf_pretransPc(scip, graph, k) )
         continue;

      /* is the kth node a potential terminal? */
      if( Is_term(graph->term[k]) && SCIPisLT(scip, graph->prize[k], FARAWAY) )
      {
         /* the future terminal */
         const int newterm = nnodes + nterms;
         nterms++;

         assert(k != root && newterm != root);
         assert(!Is_term(graph->term[newterm]));

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_addBi(scip, graph, root, newterm, prizesum + graph->prize[k]);
         graph_edge_addBi(scip, graph, k, newterm, prizesum);

         /* switch the terminal property, mark k as former terminal */
         graph_knot_chg(graph, k, STP_TERM_NONE);
         graph_knot_chg(graph, newterm, STP_TERM);
         assert(SCIPisGE(scip, graph->prize[k], 0.0));
         graph->prize[k] = 0.0;
         graph->prize[newterm] = FARAWAY;
      }
   }

   assert(nterms == npropterms);
   assert(graph->prize[root] == FARAWAY);

   SCIP_CALL( graph_transRpc(scip, graph) );

   return SCIP_OKAY;
}


/** alters the graph for MWCS problems */
SCIP_RETCODE graph_transMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   int nterms = 0;
   const int nnodes = graph_get_nNodes(graph);
   SCIP_Real* const maxweights = graph->prize;

   assert(maxweights != NULL);
   assert(scip != NULL);
   assert(graph->cost != NULL);
   assert(graph->terms == 0);

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

   SCIP_CALL( graph_transPc(scip, graph) );
   assert(graph->stp_type == STP_MWCSP);

   SCIPdebugMessage("Transformed to MW \n");

   return SCIP_OKAY;
}



/** alters the graph for RMWCS problems */
SCIP_RETCODE graph_transRmw(
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
            graph->cost[e] = -maxweights[i];
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
         graph_knot_chg(graph, i, STP_TERM);
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

   /* create new nodes */
   for( int k = 0; k < npterms; k++ )
      graph_knot_add(graph, -1);

   SCIP_CALL( graph_pc_initTerm2Edge(scip, graph, graph->knots) );

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
      else if( Is_term(graph->term[k]) )
      {
         assert(EQ(graph->prize[k], FARAWAY));
         graph_pc_knotToFixedTermProperty(graph, k);
      }
      else
      {
         assert(LE(graph->prize[k], 0.0));
      }
   }

   graph->extended = TRUE;
   graph->stp_type = STP_RMWCSP;
   graph->orgsource = graph->source;

   assert(i == npterms);
   assert(graph_valid(scip, graph));

   SCIPdebugMessage("Transformed to RMW \n");

   return SCIP_OKAY;
}


/** transforms PCSPG or MWCSP to RPCSPG or RMWCSP if possible */
SCIP_RETCODE graph_transPcmw2rooted(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real             fixprize,           /**< prize at which to make terminals */
   SCIP_Bool             verbose             /**< be verbose? */
   )
{
   int e;
   int newroot;
   int maxgrad;
   int nfixedterms;
   const int orgnterms = graph->terms;
   const int root = graph->source;
   const int pc = (graph->stp_type == STP_PCSPG);

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->term2edge != NULL);
   assert(!graph->extended);
   assert(pc || graph->stp_type == STP_MWCSP);

   newroot = -1;
   maxgrad = -1;

   if( verbose )
      printf("attempt transformation to rooted problem \n");

   nfixedterms = 0;
   e = graph->outbeg[root];
   while( e != EAT_LAST )
   {
      const int enext = graph->oeat[e];
      if( SCIPisGE(scip, graph->cost[e], fixprize) )
      {
         const int dummyterm = graph->head[e];
         const int pseudoterm = graph_pc_getTwinTerm(graph, dummyterm);

         assert(Is_pseudoTerm(graph->term[dummyterm]));
         assert(graph->grad[dummyterm] == 2);
         assert(Is_term(graph->term[pseudoterm]));
         assert(SCIPisGE(scip, graph->prize[pseudoterm], fixprize));

         graph_pc_knotToNonTermProperty(graph, dummyterm);

         graph_knot_del(scip, graph, dummyterm, TRUE);

         graph_pc_knotToFixedTermProperty(graph, pseudoterm);

         SCIPdebugMessage("fix terminal %d (delete %d)\n", pseudoterm, dummyterm);

         if( graph->grad[pseudoterm] > maxgrad )
         {
            newroot = pseudoterm;
            maxgrad = graph->grad[pseudoterm];
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

         if( Is_pseudoTerm(graph->term[k]) )
         {
            (void) graph_edge_redirect(scip, graph, e, newroot, k, graph->cost[e], TRUE, TRUE);
            graph->cost[flipedge(e)] = FARAWAY;
         }
         e = enext;
      }

      /* delete old root (cannot use default function) */
      graph_knot_chg(graph, root, STP_TERM_NONE);
      graph->term2edge[root] = TERM2EDGE_NOTERM;
      graph->prize[root] = 0.0;
      graph_knot_del(scip, graph, root, TRUE);

      if( pc )
         graph->stp_type = STP_RPCSPG;
      else
         graph->stp_type = STP_RMWCSP;

      assert(graph_valid(scip, graph));

      if( verbose )
      {
         if( pc )
            printf("...transformed PC to RPC; fixed %d out of %d terminals \n", nfixedterms, orgnterms - 1);
         else
            printf("...transformed MW to RMW; fixed %d out of %d terminals \n", nfixedterms, orgnterms - 1);
      }

      assert(orgnterms - 1 == graph->terms);
   }

   if( verbose )
   {
      if( !graph_pc_isRootedPcMw(graph) )
         printf("...failed \n");
   }

   assert(graph_valid(scip, graph));

   return SCIP_OKAY;
}




/** obtains an SAP from prize collecting problems */
SCIP_RETCODE graph_transPcGetSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Real*            offset              /**< offset (in/out) */
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

   assert(scip && graph && graph->prize && offset);
   assert(graph->knots == graph->ksize);
   assert(graph->edges == graph->esize);
   assert(graph->extended);
   assert(*offset >= 0.0);

   if( graph_pc_isPc(graph) )
   {
      *offset += graph_pc_getNonLeafTermOffset(scip, graph);
   }

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
   {
      if( Is_pseudoTerm(graph->term[k]) && (maxpvert == -1 || graph->prize[k] > graph->prize[maxpvert]) )
         maxpvert = k;
   }

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
SCIP_RETCODE graph_transPcGetRsap(
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

   SCIP_CALL( graph_pc_initPrizes(scip, p, nnodes) );
   SCIP_CALL( graph_pc_initTerm2Edge(scip, p, nnodes) );

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

   assert(p->stp_type == STP_SAP);

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

         graph_pc_knotToFixedTermProperty(p, rootcand);
      }
      graph_pc_presolExit(scip, p);
   }

   graph_knot_chg(p, oldroot, STP_TERM_NONE);
   p->prize[saproot] = 0.0;

   return SCIP_OKAY;
}


/** alters the graph for node-weighted Steiner tree problems */
SCIP_RETCODE graph_transNw2pc(
   SCIP*                 scip,               /**< SCIP data structure */
   PRESOL*               presol,             /**< presolving struct */
   GRAPH*                g                   /**< the graph */
   )
{
   SCIP_Real maxweight = 0.0;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const int nterms_org = g->terms;

   assert(scip && presol);
   assert(g->prize);

   for( int k = 0; k < nnodes; k++ )
   {
      if( maxweight < g->prize[k] )
         maxweight = g->prize[k];
   }

   assert(LT(maxweight, BLOCKED));

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) )
      {
         presol->fixed += g->prize[k];
         g->prize[k] = FARAWAY;
         g->source = k;
      }
      else
      {
         g->prize[k] = maxweight - g->prize[k];

         assert(GE(g->prize[k], 0.0));

         if( GT(g->prize[k], 0.0) )
            graph_knot_chg(g, k, STP_TERM);
      }
   }

   SCIPdebugMessage("maxweight=%f \n", maxweight);

   for( int e = 0; e < nedges; e++ )
   {
      g->cost[e] += maxweight;
   }

   assert(nterms_org >= 1);

   presol->fixed -= (nterms_org - 1) * maxweight;

   for( int k = 0; k < nnodes; k++ )
   {
      if( !EQ(g->prize[k], FARAWAY) )
         presol->fixed -= g->prize[k];
   }

   SCIPdebugMessage("presol->fixed=%f \n", presol->fixed);

   SCIP_CALL( graph_transRpc(scip, g) );

   return SCIP_OKAY;
}


/** alters the graph for node-weighted Steiner tree problems */
SCIP_RETCODE graph_transNw2sap(
   SCIP*                 scip,               /**< SCIP data structure */
   PRESOL*               presol,             /**< presolving struct */
   GRAPH*                g                   /**< the graph */
   )
{
   const int nnodes = graph_get_nNodes(g);

   assert(scip && presol);
   assert(g->prize);

   for( int k = 0; k < nnodes; k++ )
   {
      const SCIP_Real nodeweight = g->prize[k];

      if( Is_term(g->term[k]) )
      {
         presol->fixed += nodeweight;
      }
      else
      {
         /* add node-weight to edge-weights of all incoming edges */
         for( int i = g->inpbeg[i]; i != EAT_LAST; i = g->ieat[i] )
            g->cost[i] += nodeweight;
      }
   }

   return SCIP_OKAY;
}


/** alters the graph for node-weighted Steiner tree problems */
SCIP_RETCODE graph_transNw(
   SCIP*                 scip,               /**< SCIP data structure */
   PRESOL*               presol,             /**< presolving struct */
   GRAPH*                g                   /**< the graph */
   )
{
#ifdef NW_TRANS_SIMPLE
   SCIP_CALL( graph_transNw2sap(scip, presol, g) );
#else
   SCIP_CALL( graph_transNw2pc(scip, presol, g) );
#endif

   return SCIP_OKAY;
}
