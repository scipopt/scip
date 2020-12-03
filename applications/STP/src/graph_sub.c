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

/**@file   graph_sub.c
 * @brief  includes several methods for Steiner problem sub-graphs
 * @author Daniel Rehfeldt
 *
 * This file contains several basic methods to process subgraph of Steiner problem graphs.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h)                             */

//#define SCIP_DEBUG

#include "scip/misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "portab.h"
#include "graph.h"
#include "stpvector.h"



/** edge transfer */
typedef struct edge_transfer
{
   int                   source_edge;
   int                   target_tail;
   int                   target_head;
   int                   target_edge;
   SCIP_Bool             redirectEdge;
} EDGETRANS;



/** in/out */
struct subgraph_extraction_insertion
{
   STP_Vectype(int)      org_bordernodes;
   STP_Vectype(int)      org_spareedges;
   int*                  org_contractRecord;
   int*                  nodemap_subToOrg;    /**< node map */
   int*                  nodemap_orgToSub;    /**< node map */
   int                   org_nnodes;
   int                   sub_nnodes;
   SCIP_Bool             rootIsTransfered;
};


/*
 * local methods
 */



/** initializes subgraph */
static
void extractSubgraphGetSizeAndMap(
   const GRAPH*          orggraph,           /**< original graph */
   int*                  orgToSub,           /**< node map to be filled */
   int*                  subToOrg,           /**< node map to be filled */
   int*                  sub_nnodes,         /**< number of nodes */
   int*                  sub_nedges         /**< number of edges */
   )
{
   const int nnodes = graph_get_nNodes(orggraph);
   const int* const isMarked = orggraph->mark;
   int sub_n = 0;
   int sub_m = 0;
#ifdef SCIP_DEBUG
   int sub_t = 0;
#endif

   for( int i = 0; i < nnodes; i++ )
   {
      if( !isMarked[i] )
      {
         orgToSub[i] = -1;
         continue;
      }

      subToOrg[sub_n] = i;
      orgToSub[i] = sub_n;

      sub_n++;

#ifdef SCIP_DEBUG
      if( Is_term(orggraph->term[i]) )
         sub_t++;
#endif

      for( int e = orggraph->outbeg[i]; e != EAT_LAST; e = orggraph->oeat[e] )
      {
         const int head = orggraph->head[e];
         if( isMarked[head] )
            sub_m++;
      }
   }

   assert(sub_m % 2 == 0);

   *sub_nnodes = sub_n;
   *sub_nedges = sub_m;

#ifdef SCIP_DEBUG
   printf("subgraph: nodes=%d, edges=%d, terms=%d \n", sub_n, sub_m, sub_t);
#endif
}


/** helper */
static
void extractSubgraphAddNodes(
   const GRAPH*          orggraph,           /**< original graph */
   SUBINOUT*             subinout,
   GRAPH*                subgraph            /**< graph to fill */
   )
{
   const int* const nodemap_subToOrg = subinout->nodemap_subToOrg;
   const int ksize = subgraph->ksize;
   assert(subgraph->source == -1);

   for( int i = 0; i < ksize; i++ )
   {
      const int orgnode = nodemap_subToOrg[i];
      assert(graph_knot_isInRange(orggraph, orgnode));

      if( Is_term(orggraph->term[orgnode]) )
      {
         graph_knot_add(subgraph, STP_TERM);
         if( subgraph->source == -1 )
            subgraph->source = i;

         if( orggraph->source == orgnode )
         {
            subgraph->source = i;
            subinout->rootIsTransfered = TRUE;
         }
      }
      else
      {
         graph_knot_add(subgraph, STP_TERM_NONE);
      }
   }

   assert(subgraph->knots == ksize);
   assert(subgraph->source != -1);
}


/** helper */
static
SCIP_RETCODE extractSubgraphInitHistory(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< original graph */
   GRAPH*                subgraph            /**< graph to fill */
   )
{
   assert(subgraph->esize >= 2);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(subgraph->ancestors), subgraph->esize) );

   SCIP_CALL( graph_initPseudoAncestorsSized(scip, subgraph->esize, subgraph) );

   if( graph_getNpseudoAncestors(orggraph) > 0 )
      graph_addPseudoAncestors(graph_getNpseudoAncestors(orggraph), subgraph);

   assert(graph_getNpseudoAncestors(orggraph) == graph_getNpseudoAncestors(subgraph));

   SCIP_CALL( graph_copyFixed(scip, orggraph, subgraph) );

   return SCIP_OKAY;
}


/** helper */
static
SCIP_RETCODE extractSubgraphAddEdge(
   SCIP*                 scip,                /**< SCIP data structure */
   const GRAPH*          source_graph,        /**< source graph */
   const EDGETRANS*      edgetrans,
   GRAPH*                target_graph         /**< graph to fill */
   )
{
   IDX** organcestors = source_graph->ancestors;
   IDX** ancestors = target_graph->ancestors;
   const int source_edge = edgetrans->source_edge;
   const int target_edge = edgetrans->target_edge;
   const int target_tail = edgetrans->target_tail;
   const int target_head = edgetrans->target_head;
   const int target_edgeRev = flipedge(target_edge);
   const int npseudoancestors = graph_edge_nPseudoAncestors(source_graph, source_edge);
   SCIP_Bool conflict = FALSE;

   assert(graph_edge_isInRange(source_graph, source_edge));
   assert(graph_knot_isInRange(target_graph, target_tail));
   assert(graph_knot_isInRange(target_graph, target_head));

   if( npseudoancestors != 0 )
   {
      const int* const pseudoancestors = graph_edge_getPseudoAncestors(source_graph, source_edge);
      assert(pseudoancestors);

      SCIP_CALL( graph_pseudoAncestors_appendCopyArrayToEdge(scip, target_edge, pseudoancestors, npseudoancestors, target_graph, &conflict));
      assert(!conflict);
   }
   assert(graph_edge_nPseudoAncestors(target_graph, target_edge) == npseudoancestors);
   assert(graph_edge_nPseudoAncestors(target_graph, target_edgeRev) == npseudoancestors);

   ancestors[target_edge] = NULL;
   ancestors[target_edgeRev] = NULL;

   SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[target_edge]), organcestors[source_edge], &conflict) );
   assert(!conflict);
   SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[target_edgeRev]), organcestors[flipedge(source_edge)], &conflict) );
   assert(!conflict);

   assert(EQ(source_graph->cost[source_edge], source_graph->cost[flipedge(source_edge)]));

   if( edgetrans->redirectEdge )
   {
      (void) graph_edge_redirect(scip, target_graph, target_edge, target_tail, target_head,
            source_graph->cost[source_edge], FALSE, FALSE);
   }
   else
   {
      graph_edge_add(scip, target_graph, target_tail, target_head, source_graph->cost[source_edge], source_graph->cost[source_edge]);
   }

   return SCIP_OKAY;
}


/** helper */
static
SCIP_RETCODE extractSubgraphAddEdgesWithHistory(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< original graph */
   SUBINOUT*             subinout,
   GRAPH*                subgraph            /**< graph to fill */
   )
{
   const int* const isMarked = orggraph->mark;
   const int subnnodes = graph_get_nNodes(subgraph);
   const int* const nodemap_subToOrg = subinout->nodemap_subToOrg;
   const int* const nodemap_orgToSub = subinout->nodemap_orgToSub;

   assert(subgraph->ancestors);

   for( int i = 0; i < subnnodes; i++ )
   {
      const int orgtail = nodemap_subToOrg[i];
      assert(graph_knot_isInRange(orggraph, orgtail));

      for( int e = orggraph->outbeg[orgtail]; e != EAT_LAST; e = orggraph->oeat[e] )
      {
         const int orghead = orggraph->head[e];

         assert(isMarked[orgtail]);
         assert(orgtail == orggraph->tail[e]);

         /* NOTE: avoid double inclusion of edge */
         if( isMarked[orghead] && orghead > orgtail )
         {
            const int subhead = nodemap_orgToSub[orghead];
            const EDGETRANS edgetransfer = { .source_edge = e, .target_edge = subgraph->edges,
                                             .target_tail = i, .target_head = subhead, .redirectEdge = FALSE };
            assert(graph_knot_isInRange(subgraph, subhead));
            assert(nodemap_orgToSub[orgtail] == i);

            SCIP_CALL( extractSubgraphAddEdge(scip, orggraph, &edgetransfer, subgraph) );
         }
      }
   }

   subgraph->orgedges = MAX(orggraph->edges, orggraph->orgedges);
   assert(subgraph->edges == subgraph->esize);

   return SCIP_OKAY;
}


/** collects border nodes (cut vertices) of subgraph */
static
void borderNodesCollect(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< original graph */
   const GRAPH*          subgraph,           /**< sub graph */
   SUBINOUT*             subinout
   )
{
   const int* const isMarked = orggraph->mark;
   const int subnnodes = graph_get_nNodes(subgraph);
   const int* const nodemap_subToOrg = subinout->nodemap_subToOrg;

   assert(StpVecGetSize(subinout->org_bordernodes) == 0);

   for( int i = 0; i < subnnodes; i++ )
   {
      const int orgnode = nodemap_subToOrg[i];
      assert(graph_knot_isInRange(orggraph, orgnode));
      assert(isMarked[orgnode]);

      for( int e = orggraph->outbeg[orgnode]; e != EAT_LAST; e = orggraph->oeat[e] )
      {
         const int orghead = orggraph->head[e];

         if( !isMarked[orghead] )
         {
            StpVecPushBack(scip, subinout->org_bordernodes, orgnode);
            SCIPdebugMessage("added border node %d \n", orgnode);

            assert(Is_term(subgraph->term[i]));
            break;
         }
      }
   }
}


/** contracts border nodes of subgraph with remaining */
static
SCIP_RETCODE borderNodesContract(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          subgraph,           /**< sub graph */
   SUBINOUT*             subinout,           /**< helper */
   GRAPH*                orggraph            /**< original graph */
   )
{
   const int* const nodemap_subToOrg = subinout->nodemap_subToOrg;
   const int* const nodemap_orgToSub = subinout->nodemap_orgToSub;
   int* const contractRecord = subinout->org_contractRecord;
   STP_Vectype(int) bordernodes = subinout->org_bordernodes;

   for( int i = 0; i < StpVecGetSize(bordernodes); i++ )
   {
      const int orgnode = bordernodes[i];
      const int subnode = nodemap_orgToSub[orgnode];

      assert(graph_knot_isInRange(orggraph, orgnode));
      assert(graph_knot_isInRange(subgraph, subnode));
      assert(orggraph->grad[orgnode] > 0);

      if( subgraph->grad[subnode] == 0 )
      {
         int subnode_traced;
         int orgnode_traced;

         if( Is_term(subgraph->term[subnode]) )
         {
            assert(!graph_knot_hasContractTrace(subnode, subgraph));
            assert(subgraph->terms == 1);
            subnode_traced = subnode;
            continue;
         }

         subnode_traced = graph_contractTrace(subnode, subgraph);
         orgnode_traced = nodemap_subToOrg[subnode_traced];

         SCIPdebugMessage("traced %d to %d \n", orgnode, orgnode_traced);
         SCIPdebugMessage("contracting %d->%d \n", orgnode, orgnode_traced);

         /* NOTE: orgnode_traced survives */
         SCIP_CALL( graph_knot_contract(scip, orggraph, NULL, orgnode_traced, orgnode) );
         assert(contractRecord[orgnode] == -1 && contractRecord[orgnode_traced] == -1);
         contractRecord[orgnode] = orgnode_traced;
      }
   }

   if( bordernodes )
      StpVecClear(bordernodes);

   return SCIP_OKAY;
}


/** builds subgraph */
static
SCIP_RETCODE extractSubgraphBuild(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< original graph */
   SUBINOUT*             subinout,
   GRAPH**               subgraph            /**< graph to be created */
   )
{
   GRAPH* subg;
   int sub_nnodes;
   int sub_nedges;
   int* nodemap_subToOrg = subinout->nodemap_subToOrg;
   int* nodemap_orgToSub = subinout->nodemap_orgToSub;

   extractSubgraphGetSizeAndMap(orggraph, nodemap_orgToSub, nodemap_subToOrg, &sub_nnodes, &sub_nedges);

   SCIP_CALL( graph_init(scip, subgraph, sub_nnodes, sub_nedges, 1) );
   subg = *subgraph;
   assert(subg->ksize == sub_nnodes);
   subinout->sub_nnodes = sub_nnodes;

   if( graph_typeIsSpgLike(orggraph) )
      subg->stp_type = STP_SPG;

   extractSubgraphAddNodes(orggraph, subinout, subg);
   SCIP_CALL( extractSubgraphInitHistory(scip, orggraph, subg) );
   SCIP_CALL( extractSubgraphAddEdgesWithHistory(scip, orggraph, subinout, subg) );

   SCIP_CALL( graph_path_init(scip, subg) );
   SCIP_CALL( graph_initContractTracing(scip, subg) );

   return SCIP_OKAY;
}


/** helper */
static
void reinsertSubgraphTransferFixedHistory(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          subgraph,           /**< graph to be created */
   GRAPH*                orggraph            /**< original graph */
   )
{
   const int npseudoans_sub = graph_getNpseudoAncestors(subgraph);
   const int npseudoans_org = graph_getNpseudoAncestors(orggraph);

   assert(npseudoans_org <= npseudoans_sub);

   if( npseudoans_org < npseudoans_sub )
      graph_addPseudoAncestors(npseudoans_sub - npseudoans_org, orggraph);

   assert(graph_getNpseudoAncestors(subgraph) == graph_getNpseudoAncestors(orggraph));

   if( orggraph->fixedcomponents )
      graph_free_fixed(scip, orggraph);
   assert(orggraph->fixedcomponents == NULL);
   orggraph->fixedcomponents = subgraph->fixedcomponents;

   assert(graph_get_nFixpseudonodes(subgraph) == graph_get_nFixpseudonodes(orggraph));
}


/** helper */
static
SCIP_RETCODE reinsertSubgraphTransferEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          subgraph,           /**< graph to be inserted */
   SUBINOUT*             subinsertion,
   GRAPH*                orggraph            /**< original graph */
   )
{
   STP_Vectype(int) spareedges = subinsertion->org_spareedges;
   int sparecount = StpVecGetSize(spareedges) - 1;
   const int* nodemap_subToOrg = subinsertion->nodemap_subToOrg;
   const int sub_nedges = graph_get_nEdges(subgraph);

   assert(sparecount >= 0);

   for( int e = 0; e < sub_nedges; e += 2 )
   {
      assert(sparecount >= 0 || subgraph->oeat[e] == EAT_FREE);
      assert(sparecount == StpVecGetSize(spareedges) - 1);

      if( subgraph->oeat[e] != EAT_FREE )
      {
         const int orgtail = nodemap_subToOrg[subgraph->tail[e]];
         const int orghead = nodemap_subToOrg[subgraph->head[e]];
         const EDGETRANS edgetransfer = { .source_edge = e,
                                          .target_edge = spareedges[sparecount],
                                          .target_tail = orgtail, .target_head = orghead,
                                          .redirectEdge = TRUE };
         assert(graph_knot_isInRange(orggraph, orgtail));
         assert(graph_knot_isInRange(orggraph, orghead));
         assert(graph_edge_isDeleted(orggraph, spareedges[sparecount]));

         SCIP_CALL( extractSubgraphAddEdge(scip, subgraph, &edgetransfer, orggraph) );
         StpVecPopBack(spareedges);
         sparecount--;
      }
   }

   return SCIP_OKAY;
}


/** helper */
static
SCIP_RETCODE reinsertSubgraphTransferTerminals(
   const GRAPH*          subgraph,           /**< graph to be inserted */
   SUBINOUT*             subinout,
   GRAPH*                orggraph            /**< original graph */
   )
{
   const int* nodemap_subToOrg = subinout->nodemap_subToOrg;
   const int sub_nnodes = graph_get_nNodes(subgraph);

   assert(nodemap_subToOrg);

   for( int i = 0; i < sub_nnodes; i++ )
   {
      const int orgnode = nodemap_subToOrg[i];

      graph_knot_chg(orggraph, orgnode, subgraph->term[i]);
   }

   if( subinout->rootIsTransfered )
   {
      assert(graph_knot_isInRange(subgraph, subgraph->source));

      orggraph->source = nodemap_subToOrg[subgraph->source];
      assert(Is_term(orggraph->term[orggraph->source]));
   }

   return SCIP_OKAY;
}


/** helper */
static
SCIP_RETCODE reinsertSubgraphDeleteOldEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          subgraph,           /**< graph to be inserted */
   SUBINOUT*             subinout,
   GRAPH*                orggraph            /**< original graph */
   )
{
   const int* nodemap_orgToSub = subinout->nodemap_orgToSub;
   const int* nodemap_subToOrg = subinout->nodemap_subToOrg;
   const int sub_nnodes = graph_get_nNodes(subgraph);

   assert(nodemap_orgToSub && nodemap_subToOrg);

   for( int i = 0; i < sub_nnodes; i++ )
   {
      const int orgnode = nodemap_subToOrg[i];
      int e = orggraph->outbeg[orgnode];

      while( e != EAT_LAST )
      {
         const int enext = orggraph->oeat[e];
         const int orghead = orggraph->head[e];

         /* does the head correspond to subgraph? */
         if( nodemap_orgToSub[orghead] >= 0 )
         {
            assert(graph_knot_isInRange(subgraph, nodemap_orgToSub[orghead]));
            StpVecPushBack(scip, subinout->org_spareedges, e);
            graph_edge_del(scip, orggraph, e, TRUE);
         }

         e = enext;
      }
   }

   return SCIP_OKAY;
}


/** builds subgraph */
static
SCIP_RETCODE reinsertSubgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          subgraph,           /**< graph to be inserted */
   SUBINOUT*             subinout,
   GRAPH*                orggraph            /**< original graph */
   )
{
   assert(StpVecGetSize(subinout->org_spareedges) == 0);
   assert(subgraph->edges <= orggraph->edges);
   assert(subgraph->knots <= orggraph->knots);
   assert(subinout->sub_nnodes == subgraph->knots);
   assert(subinout->org_nnodes == orggraph->knots);

   reinsertSubgraphTransferFixedHistory(scip, subgraph, orggraph);
   SCIP_CALL( reinsertSubgraphDeleteOldEdges(scip, subgraph, subinout, orggraph) );
   SCIP_CALL( reinsertSubgraphTransferEdges(scip, subgraph, subinout, orggraph) );
   reinsertSubgraphTransferTerminals(subgraph, subinout, orggraph);

   SCIP_CALL( borderNodesContract(scip, subgraph, subinout, orggraph) );

   StpVecClear(subinout->org_spareedges);

   return SCIP_OKAY;
}


/*
 * Interface methods
 */



/** Obtains a new subgraph corresponding to marked nodes.
 *  Also fills a node map from the new to the original nodes.
 *  Creates a shallow copy wrt ancestor information:
 *  ancestors and pseudo-ancestors are kept and will be modified also for original graph! */
SCIP_RETCODE graph_subgraphExtract(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< original graph */
   SUBINOUT*             subinout,
   GRAPH**               subgraph            /**< graph to be created */
   )
{
   assert(scip && orggraph && subinout);
   assert(subinout->org_nnodes == orggraph->knots);
   assert(!graph_pc_isPcMw(orggraph) && "not yet supported");

   SCIP_CALL( extractSubgraphBuild(scip, orggraph, subinout, subgraph) );

   borderNodesCollect(scip, orggraph, *subgraph, subinout);

   return SCIP_OKAY;
}


/** initializes */
SCIP_RETCODE graph_subinoutInit(
  SCIP*                 scip,               /**< SCIP data structure */
  const GRAPH*          orggraph,           /**< original graph */
  SUBINOUT**            subinout
  )
{
   SUBINOUT* sub;
   const int nnodes = graph_get_nNodes(orggraph);

   assert(scip && orggraph);

   SCIP_CALL( SCIPallocMemory(scip, subinout) );
   sub = *subinout;

   sub->org_nnodes = nnodes;
   sub->sub_nnodes = -1;
   sub->rootIsTransfered = FALSE;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(sub->org_contractRecord), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sub->nodemap_subToOrg), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sub->nodemap_orgToSub), nnodes) );
   sub->org_spareedges = NULL;
   sub->org_bordernodes = NULL;

   for( int i = 0; i < nnodes; i++ )
      sub->org_contractRecord[i] = -1;

   return SCIP_OKAY;
}


/** get nodes map */
const int* graph_subinoutGetSubToOrgNodeMap(
  const SUBINOUT*       subinout
  )
{
   assert(subinout);

   return subinout->nodemap_subToOrg;
}



/** get contraction record for cut nodes (-1 if no contraction) */
const int* graph_subinoutGetContractionRecord(
  const SUBINOUT*       subinout
  )
{
   assert(subinout);

   return subinout->org_contractRecord;
}



/** gets ancesotr */
int graph_knot_getContractionRecordAncestor(
  int                   node,
  const SUBINOUT*       subinout
  )
{
   int ancestor;
   const int* record;

   assert(subinout);
   assert(subinout->org_contractRecord);
   assert(0 <= node && node < subinout->org_nnodes);

   record = subinout->org_contractRecord;


   for( ancestor = node; record[ancestor] != -1; ancestor = record[ancestor]  )
   {
      assert(0 <= ancestor && ancestor < subinout->org_nnodes);
   }

   return ancestor;
}



/** frees */
void graph_subinoutFree(
  SCIP*                 scip,               /**< SCIP data structure */
  SUBINOUT**            subinout
  )
{
   SUBINOUT* sub;

   assert(scip && subinout);

   sub = *subinout;

   StpVecFree(scip, sub->org_bordernodes);
   StpVecFree(scip, sub->org_spareedges);
   SCIPfreeMemoryArray(scip, &(sub->nodemap_orgToSub));
   SCIPfreeMemoryArray(scip, &(sub->nodemap_subToOrg));
   SCIPfreeMemoryArray(scip, &(sub->org_contractRecord));

   SCIPfreeMemory(scip, subinout);
}


/** reinserts extracted (and modified) subgraph; also deletes subgraph.
 *  Dual to "graph_subgraphExtract"  */
SCIP_RETCODE graph_subgraphReinsert(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBINOUT*             subinout,
   GRAPH*                orggraph,           /**< original graph */
   GRAPH**               subgraph            /**< graph to be reinserted (and freed) */
   )
{
   assert(scip && subinout && orggraph && subgraph);
   assert(*subgraph);
   assert(!graph_pc_isPcMw(orggraph) && "not yet supported");
   assert(graph_valid(scip, orggraph));

   SCIP_CALL( reinsertSubgraph(scip, *subgraph, subinout, orggraph) );

   graph_path_exit(scip, *subgraph);

   /* NOTE: fixed components are not deleted */
   graph_free(scip, subgraph, FALSE);

   assert(graph_valid(scip, orggraph));

   return SCIP_OKAY;
}


/** frees extracted subgraph */
SCIP_RETCODE graph_subgraphFree(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               subgraph            /**< graph to be freed */
   )
{
   assert(scip && subgraph);

   graph_path_exit(scip, *subgraph);
   graph_free(scip, subgraph, TRUE);

   return SCIP_OKAY;
}
