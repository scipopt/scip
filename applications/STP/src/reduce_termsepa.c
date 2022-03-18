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

/**@file   reduce_termsepa.c
 * @brief  base class for terminal-separator reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements base classes for terminal-separator reduction techniques Steiner tree problems.
 * The actual algorithms can be found in reduce_termsepada.c and reduce_termsepafull.c
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "extreduce.h"
#include "solstp.h"
#include "mincut.h"
#include "portab.h"
#include "stpvector.h"
#include "scip/scip.h"



#define MARK_NONACTIVE 0
#define MARK_SUBNODE   1
#define MARK_SEPARATOR 2

/*
 * Data structures
 */


/** tree bottleneck node */
typedef struct tree_bottleneck_node
{
   SCIP_Real             edgecost;           /**< cost of outgoing edge */
   SCIP_Real             bottleneck;         /**< (temporary) bottleneck up to this node; -1.0 for UNSET */
   int                   parent;             /**< parent node index */
   int                   degree;             /**< degree of node */
} TBNODE;



/** terminal separator tree bottleneck structure */
struct terminal_separator_tree_bottleneck
{
   TBNODE*               tree;               /**< tree for representation; NOTE: oriented towards root */
   int                   nnodes;             /**< number of nodes of underlying graph */
};



/*
 * Local methods
 */



/** initializes */
static
SCIP_RETCODE tbottleneckInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            soledges,           /**< solution tree to be represented */
   TBOTTLENECK**         tbottleneck         /**< to initialize */
   )
{
   TBOTTLENECK* tbneck;
   TBNODE* tbtree;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);

   assert(scip && soledges);
   assert(solstp_isValid(scip, g, soledges));

   SCIP_CALL( SCIPallocMemory(scip, tbottleneck) );
   tbneck = *tbottleneck;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(tbneck->tree), nnodes) );
   tbtree = tbneck->tree;
   tbneck->nnodes = nnodes;

   for( int i = 0; i < nnodes; i++ )
   {
      tbtree[i].edgecost = -FARAWAY;
      tbtree[i].parent = UNKNOWN;
      tbtree[i].bottleneck = -1.0;

      /* NOTE: slight hack to make sure that terminals are not taken as part of key path */
      if( Is_term(g->term[i]) )
      {
         tbtree[i].degree = nnodes;
      }
      else
      {
         tbtree[i].degree = 0.0;
      }
   }

#ifdef SCIP_DEBUG
   graph_printInfo(g);
#endif

   for( int e = 0; e < nedges; e++ )
   {
      if( soledges[e] == CONNECT )
      {
         const int tail = g->tail[e];
         const int head = g->head[e];

#ifdef SCIP_DEBUG
         SCIPdebugMessage("soledge: ");
         graph_edge_printInfo(g, e);
#endif

         assert(tbtree[head].parent == UNKNOWN);

         tbtree[head].parent = tail;
         tbtree[head].edgecost = g->cost[e];

         tbtree[tail].degree++;
         tbtree[head].degree++;
      }
   }

   assert(tbtree[g->source].parent == UNKNOWN);

   return SCIP_OKAY;
}

// todo really needed? There should be not need for the key-paths to be disjoint!
#ifdef SCIP_DISABLED
/** helper */
static
void tbottleneckCut(
   int                   startnode,          /**< start node */
   int                   endnode,            /**< end node */
   SCIP_Real             maxlength,          /**< length of bottleneck */
   TBOTTLENECK*          tbottleneck         /**< tree bottleneck structure */
   )
{
   TBNODE* tbtree = tbottleneck->tree;
   int k;
   int cutstart = startnode;
   int cutend;
   SCIP_Real max_local = 0.0;
#ifndef NDEBUG
   SCIP_Real max_debug = 0.0;
#endif

   assert(startnode != endnode);
   assert(tbottleneck->nnodes >= 2);

   /* go up to the root, stop at end node  */
   for( k = startnode; k != endnode; k = tbtree[k].parent )
   {
      assert(0 <= k && k < tbottleneck->nnodes);

      if( GE(max_local, maxlength) )
      {
         break;
      }

      assert(tbtree[k].parent != UNKNOWN);

      if( tbtree[k].degree > 2 )
      {
         cutstart = k;
         max_local = 0.0;
      }

      max_local += tbtree[k].edgecost;
   }

   assert(EQ(max_local, maxlength));
   cutend = k;

   for( k = cutstart; k != cutend; k = tbtree[k].parent )
   {
      assert(k == cutstart || tbtree[k].degree <= 2);
#ifndef NDEBUG
      max_debug += tbtree[k].edgecost;
#endif

      SCIPdebugMessage("removing tree bottleneck node=%d, degree=%d edgecost=%f \n", k, tbtree[k].degree, tbtree[k].edgecost);
      tbtree[k].edgecost = -FARAWAY;
   }

#ifndef NDEBUG
   assert(EQ(max_debug, maxlength));
#endif

}
#endif

/** gets tree bottleneck cost  */
static
void tbottleneckGetMax(
   int                   node1,              /**< first node */
   int                   node2,              /**< second node */
   TBOTTLENECK*          tbottleneck,        /**< tree bottleneck structure */
   SCIP_Real*            maxlength           /**< IN/OUT! user needs to provide a maximum that is to be surpassed */
   )
{
   int k;
   TBNODE* tbtree = tbottleneck->tree;
   SCIP_Real max_local = 0.0;
   SCIP_Real max1 = 0.0;
   SCIP_Real max2 = 0.0;
   SCIP_Real maxtotal;

   assert(node1 >= 0 && node2 >= 0);
   assert(node1 != node2);
   assert(GT(*maxlength, 0.0));

   /* go up to the root and mark the tree on the way */
   for( k = node1; k != UNKNOWN; k = tbtree[k].parent )
   {
      assert(0 <= k && k < tbottleneck->nnodes);
      assert(EQ(tbtree[k].bottleneck, -1.0));

      tbtree[k].bottleneck = max1;

      if( tbtree[k].degree > 2 )
         max_local = 0.0;

      if( tbtree[k].parent != UNKNOWN )
      {
         max_local += tbtree[k].edgecost;
      }

      if( max_local > max1 )
         max1 = max_local;
   }

   max_local = 0.0;

   /* go up from second node */
   for( k = node2; tbtree[k].bottleneck < -0.5; k = tbtree[k].parent )
   {
      assert(0 <= k && k < tbottleneck->nnodes);
      /* NOTE: we should not reach the root */
      assert(tbtree[k].parent != UNKNOWN);

      if( tbtree[k].degree > 2 )
         max_local = 0.0;

      max_local += tbtree[k].edgecost;

      if( max_local > max2 )
         max2 = max_local;
   }

   assert(GE(tbtree[k].bottleneck, 0.0));
   max1 = tbtree[k].bottleneck;

   SCIPdebugMessage("joint node: %d \n", k);
   SCIPdebugMessage("tree bdist from node %d: %f \n", node1, max1);
   SCIPdebugMessage("tree bdist from node %d: %f \n", node2, max2);

   maxtotal = MAX(max1, max2);

   /* not yet removed sub-path and improvement over given distance? */
   if( GT(maxtotal, *maxlength)  )
   {
      SCIPdebugMessage("improved tree bottleneck  %f -> %f \n", *maxlength, maxtotal);

      *maxlength = maxtotal;

#ifdef SCIP_DISABELD
      if( GT(max1, max2) )
         tbottleneckCut(node1, k, max1, tbottleneck);
      else
         tbottleneckCut(node2, k, max2, tbottleneck);
#endif
   }

   /* reset */
   for( k = node1; k != UNKNOWN; k = tbtree[k].parent )
   {
      assert(0 <= k && k < tbottleneck->nnodes);
      tbtree[k].bottleneck = -1.0;
   }
}


/** frees */
static
void tbottleneckFree(
   SCIP*                 scip,               /**< SCIP data structure */
   TBOTTLENECK**         tbottleneck         /**< to initialize */
   )
{
   TBOTTLENECK* tbneck;
   tbneck = *tbottleneck;

   SCIPfreeMemoryArray(scip, &(tbneck->tree));
   SCIPfreeMemory(scip, tbottleneck);
}


/** identifies subgraph */
static
void subgraphIdentify(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   TERMCOMP*             termcomp            /**< component */
   )
{
   const COMPBUILDER* const builder = termcomp->builder;
   const int* const sepaterms = builder->sepaterms;
   int* const nodemap_orgToSub = termcomp->nodemap_orgToSub;
   STP_Vectype(int) bfsqueue = NULL;
   int* RESTRICT gmark = termcomp->nodes_mark;
   int sub_e = 0;
   int sub_n = 0;
   const int nsepaterms = builder->nsepatterms;
   const int sourceterm = builder->sourceterm;

   assert(graph_knot_isInRange(g, sourceterm) && Is_term(g->term[sourceterm]));
   assert(sepaterms && nodemap_orgToSub);
   assert(nsepaterms >= 2 && nsepaterms < g->terms);

   BMSclearMemoryArray(gmark, g->knots);

   /* mark separator */
   for( int i = 0; i < nsepaterms; i++ )
   {
      const int term = sepaterms[i];
      assert(graph_knot_isInRange(g, term) && Is_term(g->term[term]));
      assert(!gmark[term]);
      assert(nodemap_orgToSub[term] == UNKNOWN);

      SCIPdebugMessage("mapping separator %d to %d \n", term, sub_n);

      nodemap_orgToSub[term] = sub_n;
      sub_n++;
      sub_e += g->grad[term];
      gmark[term] = MARK_SEPARATOR;
   }

   assert(!gmark[sourceterm]);
   StpVecReserve(scip, bfsqueue, 16);
   StpVecPushBack(scip, bfsqueue, sourceterm);
   gmark[sourceterm] = MARK_SUBNODE;

   /* BFS loop */
   for( int i = 0; i < StpVecGetSize(bfsqueue); i++ )
   {
      const int k = bfsqueue[i];

      assert(gmark[k]);
      assert(nodemap_orgToSub[k] == UNKNOWN);

      nodemap_orgToSub[k] = sub_n;

      SCIPdebugMessage("mapping component node %d to %d \n", k, sub_n);

      sub_n++;
      sub_e += g->grad[k];

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         if( !gmark[head] )
         {
            StpVecPushBack(scip, bfsqueue, head);
            gmark[head] = MARK_SUBNODE;
         }
      }
   }

   assert(termcomp->bfsqueue == NULL);
   assert(termcomp->subnnodes == -1);
   assert(termcomp->subnedges == -1);

   /* reserve space for the separator terminal clique */
   sub_e += (nsepaterms) * (nsepaterms - 1);

   termcomp->bfsqueue = bfsqueue;
   termcomp->subnnodes = sub_n;
   termcomp->subnedges = sub_e;
}


/**  creates subgraph */
static
SCIP_RETCODE subgraphBuild(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< graph data structure */
   SCIP_Bool             useSd,              /**< use special distances */
   EXTPERMA*             extperma,           /**< extension data or NULL */
   TERMCOMP*             termcomp            /**< component */
   )
{
   GRAPH* subgraph;
   DISTDATA* const distdata = useSd ? extperma->distdata_default : NULL;
   const COMPBUILDER* const builder = termcomp->builder;
   const int* const nodemap_orgToSub = termcomp->nodemap_orgToSub;
   int* nodemap_subToOrg;
   int* edgemap_subToOrg;
   const int* const orgmark = termcomp->nodes_mark;
   const int* const sepaterms = builder->sepaterms;
   const int nsepaterms = builder->nsepatterms;
   const int nnodes_sub = termcomp->subnnodes;
   const int nedges_sub = termcomp->subnedges;
   STP_Vectype(int) bfsqueue = termcomp->bfsqueue;

   assert(nsepaterms >= 2 && nsepaterms < orggraph->terms);
   assert(nnodes_sub > 0 && nedges_sub > 0);
   assert(nodemap_orgToSub && sepaterms);
   assert(!termcomp->subgraph);
   assert((extperma != NULL) == useSd);

   SCIP_CALL( graph_init(scip, &(termcomp->subgraph), nnodes_sub, nedges_sub, 1) );
   subgraph = termcomp->subgraph;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(nodemap_subToOrg), nnodes_sub) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(edgemap_subToOrg), nedges_sub) );

   for( int i = 0; i < nsepaterms; i++ )
   {
      assert(graph_knot_isInRange(orggraph, sepaterms[i]) && Is_term(orggraph->term[sepaterms[i]]));
      assert(nodemap_orgToSub[sepaterms[i]] == subgraph->knots);

      graph_knot_add(subgraph, STP_TERM);
      nodemap_subToOrg[i] = sepaterms[i];
   }

   for( int i = 0; i < StpVecGetSize(bfsqueue); i++ )
   {
      const int orgnode = bfsqueue[i];
      assert(nodemap_orgToSub[orgnode] == subgraph->knots);

      nodemap_subToOrg[subgraph->knots] = orgnode;
      graph_knot_add(subgraph, orggraph->term[orgnode]);
   }

#ifndef NDEBUG
   assert(nnodes_sub == subgraph->knots);
   for( int i = 0; i < nnodes_sub; i++ )
   {
      assert(i == nodemap_orgToSub[nodemap_subToOrg[i]]);
   }
#endif

   for( int subsepaterm = 0; subsepaterm < nsepaterms; subsepaterm++ )
   {
      const int orgsepaterm = sepaterms[subsepaterm];

      assert(graph_knot_isInRange(orggraph, orgsepaterm) && Is_term(orggraph->term[orgsepaterm]));
      assert(orgmark[orgsepaterm]);
      assert(nodemap_subToOrg[subsepaterm] == orgsepaterm);

      for( int e = orggraph->outbeg[orgsepaterm]; e != EAT_LAST; e = orggraph->oeat[e] )
      {
         const int orghead = orggraph->head[e];

         /* NOTE: add only small to big edges, to avoid adding an edge twice */
         if( orghead < orgsepaterm )
            continue;

         /* NOTE: ignore separator nodes, to not add edge twice */
         if( orgmark[orghead] == MARK_SUBNODE )
         {
            const int subnode = nodemap_orgToSub[orghead];
            assert(subnode >= 0);
            assert(EQ(orggraph->cost[e], orggraph->cost[flipedge(e)]));

            edgemap_subToOrg[subgraph->edges] = e;
            edgemap_subToOrg[subgraph->edges + 1] = flipedge(e);
            graph_edge_addBi(scip, subgraph, subsepaterm, subnode, orggraph->cost[e]);

#ifdef SCIP_DEBUG
            SCIPdebugMessage("adding separator-to-default edge: ");
            graph_edge_printInfo(subgraph, subgraph->edges - 2);
#endif
         }
      }

      for( int subsepaterm2 = subsepaterm + 1; subsepaterm2 < nsepaterms; subsepaterm2++ )
      {
         const int orgsepaterm2 = sepaterms[subsepaterm2];
         const SCIP_Real sd = useSd ?
               extreduce_distDataGetSdDouble(scip, orggraph, orgsepaterm, orgsepaterm2, distdata)
             : BLOCKED;

         assert(GE(sd, 0.0));

         edgemap_subToOrg[subgraph->edges] = -1;
         edgemap_subToOrg[subgraph->edges + 1] = -1;
         graph_edge_addBi(scip, subgraph, subsepaterm, subsepaterm2, sd);

#ifdef SCIP_DEBUG
         SCIPdebugMessage("adding separator-to-separator edge: ");
         graph_edge_printInfo(subgraph, subgraph->edges - 2);
#endif
      }
   }

   for( int i = 0; i < StpVecGetSize(bfsqueue); i++ )
   {
      const int orgnode = bfsqueue[i];
      const int subnode = nodemap_orgToSub[orgnode];

      for( int e = orggraph->outbeg[orgnode]; e != EAT_LAST; e = orggraph->oeat[e] )
      {
         const int orghead = orggraph->head[e];

         if( orghead < orgnode )
            continue;

         if( orgmark[orghead] )
         {
            const int subhead = nodemap_orgToSub[orghead];
            assert(subhead >= 0);
            assert(EQ(orggraph->cost[e], orggraph->cost[flipedge(e)]));

            edgemap_subToOrg[subgraph->edges] = e;
            edgemap_subToOrg[subgraph->edges + 1] = flipedge(e);
            graph_edge_addBi(scip, subgraph, subnode, subhead, orggraph->cost[e]);

#ifdef SCIP_DEBUG
            SCIPdebugMessage("adding default-to-? edge: ");
            graph_edge_printInfo(subgraph, subgraph->edges - 2);
#endif
         }
      }
   }

   /* finally, we need to update the edge map on sepa-sepa edges */
   for( int subsepaterm = 0; subsepaterm < nsepaterms; subsepaterm++ )
   {
      const int orgsepaterm = sepaterms[subsepaterm];

      assert(orgmark[orgsepaterm]);
      assert(nodemap_subToOrg[subsepaterm] == orgsepaterm);

      for( int orge = orggraph->outbeg[orgsepaterm]; orge != EAT_LAST; orge = orggraph->oeat[orge] )
      {
         const int orghead = orggraph->head[orge];
         if( orghead < orgsepaterm )
            continue;

         if( orgmark[orghead] == MARK_SEPARATOR )
         {
            const int subhead = nodemap_orgToSub[orghead];
            assert(subhead >= 0);

            /* find and update the corresponding edge in subgraph */
            for( int e = subgraph->outbeg[subsepaterm]; e != EAT_LAST; e = subgraph->oeat[e] )
            {
               if( subgraph->head[e] == subhead )
               {
                  const SCIP_Real newcost = MIN(orggraph->cost[orge], subgraph->cost[e]);
                  assert(edgemap_subToOrg[e] == -1);

                  edgemap_subToOrg[e] = orge;
                  edgemap_subToOrg[flipedge(e)] = flipedge(orge);

                  subgraph->cost[e] = newcost;
                  subgraph->cost[flipedge(e)] = newcost;

                  break;
               }
            }
         }
      }
   }

   assert(graph_knot_isInRange(orggraph, builder->sourceterm));
   assert(graph_knot_isInRange(orggraph, nodemap_orgToSub[builder->sourceterm]));
   assert(orggraph->stp_type != UNKNOWN);

   subgraph->source = nodemap_orgToSub[builder->sourceterm];
   subgraph->stp_type = orggraph->stp_type;

   SCIP_CALL( graph_path_init(scip, subgraph) );

   assert(!termcomp->nodemap_subToOrg && !termcomp->edgemap_subToOrg);

   termcomp->nodemap_subToOrg = nodemap_subToOrg;
   termcomp->edgemap_subToOrg = edgemap_subToOrg;

   return SCIP_OKAY;
}



/** gets extended bottleneck distances */
static
SCIP_Real getExtBottleneckDist(
   const GRAPH*          g,                  /**< graph data structure */
   int                   term1,              /**< first terminal */
   int                   term2,              /**< second terminal */
   int                   term1_sub,          /**< corresponding terminal in subgraph */
   int                   term2_sub,          /**< corresponding terminal in subgraph */
   int                   mstsource,          /**< source node of minimum spanning tree */
   PATH*                 mst,                /**< minimum spanning tree */
   TBOTTLENECK*          subsolbottleneck,   /**< tree bottleneck */
   SCIP_Real* RESTRICT   mstsdist            /**< node distance helper */
   )
{
   SCIP_Real bdist = 0.0;
   int tempnode = term1;

   assert(mst && mstsdist);
   assert(Is_term(g->term[term1]));
   assert(Is_term(g->term[term2]));
   assert(term1 != term2);

   mstsdist[tempnode] = 0.0;

   while( tempnode != mstsource )
   {
      const int ne = mst[tempnode].edge;

      assert(ne >= 0);
      assert(g->head[ne] == tempnode);
      tempnode = g->tail[ne];

      if( g->cost[ne] > bdist )
         bdist = g->cost[ne];

      mstsdist[tempnode] = bdist;
      if( tempnode == term2 )
         break;
   }

   /* already finished? */
   if( tempnode == term2 )
   {
      tempnode = mstsource;
   }
   else
   {
      tempnode = term2;
      bdist = 0.0;
   }

   while( tempnode != mstsource )
   {
      const int ne = mst[tempnode].edge;
      assert(ne >= 0);

      tempnode = g->tail[ne];

      if( g->cost[ne] > bdist )
         bdist = g->cost[ne];

      /* already visited? */
      if( mstsdist[tempnode] > -0.5 )
      {
         if( mstsdist[tempnode] > bdist )
            bdist = mstsdist[tempnode];
         break;
      }

      assert(EQ(mstsdist[tempnode], -1.0));
   }

   /* restore */
   tempnode = term1;
   mstsdist[tempnode] = -1.0;
   while( tempnode != mstsource )
   {
      const int ne = mst[tempnode].edge;
      tempnode = g->tail[ne];
      mstsdist[tempnode] = -1.0;
      if( tempnode == term2 )
         break;
   }

#ifndef NDEBUG
   for( int i = 0; i < g->knots; i++ )
      assert(EQ(mstsdist[i], -1.0));
#endif

   /* updates with tree bottleneck of sub-solution if available */
   if( subsolbottleneck )
      tbottleneckGetMax(term1_sub, term2_sub, subsolbottleneck, &bdist);

   return bdist;
}



/*
 * Interface methods
 */



/** initializes */
SCIP_RETCODE reduce_compbuilderInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   COMPBUILDER**         compbuilder         /**< to initialize */
   )
{
   COMPBUILDER* builder;
   const int nnodes = graph_get_nNodes(g);

   SCIP_CALL( SCIPallocMemory(scip, compbuilder) );
   builder = *compbuilder;

   builder->sepaterms = NULL;
   builder->sourceterm = -1;
   builder->nsepatterms = 2;
   builder->componentnumber = 0;
   builder->ncomponentnodes = -1;
   builder->ngraphnodes = -1;
   builder->rootcompIsProcessed = FALSE;
   builder->maxncompchecks = -1;
   builder->maxsepasize = -1;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(builder->nodes_bdist), nnodes) );
   graph_get_nVET(g, &(builder->ngraphnodes), NULL, NULL);

   for( int i = 0; i < nnodes; i++ )
   {
      builder->nodes_bdist[i] = -1.0;
   }

   return SCIP_OKAY;
}


/** frees */
void reduce_compbuilderFree(
   SCIP*                 scip,               /**< SCIP data structure */
   COMPBUILDER**         compbuilder         /**< to free */
   )
{
   COMPBUILDER* builder;
   builder = *compbuilder;

   SCIPfreeMemoryArray(scip, &(builder->nodes_bdist));
   SCIPfreeMemory(scip, compbuilder);
}


/** gets nodes ratio of subgraph  */
SCIP_Real reduce_compbuilderGetSubNodesRatio(
   const COMPBUILDER*    compbuilder         /**< builder */
   )
{
   const SCIP_Real ratio = (SCIP_Real) compbuilder->ncomponentnodes / (SCIP_Real) compbuilder->ngraphnodes;

   assert(GT(ratio, 0.0));
   assert(LE(ratio, 1.0));

   return ratio;
}


/** prints */
void reduce_compbuilderPrintSeparators(
   const GRAPH*          orggraph,           /**< graph data structure */
   const COMPBUILDER*    compbuilder         /**< builder */
   )
{
   assert(orggraph && compbuilder);

   printf("number of separator nodes=%d, nodes: \n", compbuilder->nsepatterms);

   for( int i = 0; i < compbuilder->nsepatterms; i++ )
   {
      const int sepaterm = compbuilder->sepaterms[i];
      assert(graph_knot_isInRange(orggraph, sepaterm));

      graph_knot_printInfo(orggraph, sepaterm);
   }
}


/** builds extended subgraph with SD weighted edges between terminal-separator nodes */
SCIP_RETCODE reduce_termcompBuildSubgraphWithSds(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                orggraph,           /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   TERMCOMP*             termcomp            /**< component */
   )
{
   const SCIP_Bool useSd = TRUE;

   assert(scip && orggraph && extperma && termcomp);

   subgraphIdentify(scip, orggraph, termcomp);
   SCIP_CALL( subgraphBuild(scip, orggraph, useSd, extperma, termcomp) );

   return SCIP_OKAY;
}


/** builds extended subgraph with BLOCKED weighted edges between terminal-separator nodes */
SCIP_RETCODE reduce_termcompBuildSubgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                orggraph,           /**< graph data structure */
   TERMCOMP*             termcomp            /**< component */
   )
{
   const SCIP_Bool useSd = FALSE;

   assert(scip && orggraph && termcomp);

   subgraphIdentify(scip, orggraph, termcomp);
   SCIP_CALL( subgraphBuild(scip, orggraph, useSd, NULL, termcomp) );

   return SCIP_OKAY;
}


/** Changes weights between terminal separator nodes to bottlenecks.
 *  NOTE: also computes the bottleneck distances, so potentially expensive! */
SCIP_RETCODE reduce_termcompChangeSubgraphToBottleneck(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   TERMCOMP*             termcomp,           /**< component */
   SCIP_Bool*            success             /**< pointer to store whether computation was successful (OUT)  */
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   PATH* mst;
   COMPBUILDER* const builder = termcomp->builder;
   const int* const sepaterms = builder->sepaterms;
   const int* const nodemap_orgToSub = termcomp->nodemap_orgToSub;
   const int* const nodemap_subToOrg = termcomp->nodemap_subToOrg;
   const int* const nodes_mark = termcomp->nodes_mark;
   const int nsepaterms = builder->nsepatterms;
   const int nnodes = graph_get_nNodes(g);
   const int mstroot = sepaterms[0];

   *success = TRUE;

   /* mark the anti-component */
   for( int i = 0; i < nnodes; i++ )
      g->mark[i] = (nodes_mark[i] != MARK_SUBNODE);

   SCIPallocBufferArray(scip, &mst, nnodes);
   graph_path_exec(scip, g, MST_MODE, mstroot, g->cost, mst);

   for( int i = 0; i < nsepaterms; i++ )
   {
      const int term = sepaterms[i];
      if( term != mstroot && mst[term].edge == -1 )
      {
         SCIPfreeBufferArray(scip, &mst);
         *success = FALSE;
         graph_mark(g);
         return SCIP_OKAY;
      }
   }

#ifndef NDEBUG
   for( int i = 0; i < g->knots; i++ )
   {
      if( Is_term(g->term[i]) && g->mark[i] && i != mstroot )
      {
         assert(mst[i].edge != -1);
      }
   }
#endif

   for( int i = 0; i < subgraph->knots; i++ )
      subgraph->mark[i] = MARK_NONACTIVE;

   for( int i = 0; i < nsepaterms; i++ )
   {
      const int term = sepaterms[i];
      const int term_sub = nodemap_orgToSub[term];
      assert(graph_knot_isInRange(subgraph, term_sub));
      assert(subgraph->mark[term_sub] == MARK_NONACTIVE);

      subgraph->mark[term_sub] = MARK_SEPARATOR;
   }

   for( int i = 0; i < nsepaterms; i++ )
   {
      const int term = sepaterms[i];
      const int term_sub = nodemap_orgToSub[term];

      for( int e = subgraph->outbeg[term_sub]; e != EAT_LAST; e = subgraph->oeat[e] )
      {
         const int head_sub = subgraph->head[e];

         if( head_sub > term_sub && subgraph->mark[head_sub] == MARK_SEPARATOR )
         {
            const int head = nodemap_subToOrg[head_sub];
            const SCIP_Real b = getExtBottleneckDist(g, term, head, term_sub, head_sub, mstroot, mst,
                  termcomp->subsolbottleneck, builder->nodes_bdist);

            subgraph->cost[e] = b;
            subgraph->cost[flipedge(e)] = b;

            SCIPdebugMessage("%d->%d setting b=%f \n", term_sub, head_sub, b);
         }
      }
   }

   graph_mark(g);
   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}


/** changes weights between terminal separator nodes to original costs (or BLOCKED) */
void reduce_termcompChangeSubgraphToOrgCosts(
   const GRAPH*          orggraph,           /**< graph data structure */
   TERMCOMP*             termcomp            /**< component */
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   const COMPBUILDER* const builder = termcomp->builder;
   const int* const edgemap_subToOrg = termcomp->edgemap_subToOrg;
   const int nsepaterms = builder->nsepatterms;

   assert(subgraph && orggraph);
   assert(nsepaterms > 0);

   /*  NOTE: we exploit the fact that the first nodes of subgraph are the separator terminals */

   for( int subsepaterm = 0; subsepaterm < nsepaterms; subsepaterm++ )
   {
      for( int sube = subgraph->outbeg[subsepaterm]; sube != EAT_LAST; sube = subgraph->oeat[sube] )
      {
         const int subhead = subgraph->head[sube];
         if( subhead < nsepaterms )
         {
            const int orge = edgemap_subToOrg[sube];

            assert(Is_term(subgraph->term[subsepaterm]));
            assert(Is_term(subgraph->term[subhead]));

            if( orge >= 0 )
            {
               assert(Is_term(orggraph->term[orggraph->head[orge]]));
               assert(Is_term(orggraph->term[orggraph->tail[orge]]));

               subgraph->cost[sube] = orggraph->cost[orge];
            }
            else
            {
               subgraph->cost[sube] = BLOCKED;
            }
         }
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < subgraph->edges; e += 2 )
   {
      assert(EQ(subgraph->cost[e], subgraph->cost[e + 1]));
   }
#endif
}


/** initializes */
SCIP_RETCODE reduce_termcompInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   COMPBUILDER*          builder,            /**< initializer */
   TERMCOMP**            termcomp            /**< to initialize */
   )
{
   TERMCOMP* comp;

   SCIP_CALL( SCIPallocMemory(scip, termcomp) );
   comp = *termcomp;

   comp->builder = builder;
   comp->subgraph = NULL;
   comp->subsolution = NULL;
   comp->edgemap_subToOrg = NULL;
   comp->nodemap_subToOrg = NULL;
   comp->bfsqueue = NULL;
   comp->subnedges = -1;
   comp->subnnodes = -1;
   comp->subprimalobj = -FARAWAY;
   comp->subsolbottleneck = NULL;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(comp->nodes_mark), g->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(comp->nodemap_orgToSub), g->knots) );

#ifndef NDEBUG
   for( int i = 0; i < g->knots; i++ )
      comp->nodemap_orgToSub[i] = UNKNOWN;
#endif

   return SCIP_OKAY;
}



/** initializes tree bottleneck for terminal component */
SCIP_RETCODE reduce_termcompInitTbottleneck(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            subsoledges,        /**< solution tree to be represented */
   TERMCOMP*             termcomp            /**< to initialize for */
   )
{
   const GRAPH* subgraph;

   assert(scip && subsoledges && termcomp);
   assert(!termcomp->subsolbottleneck);

   subgraph = termcomp->subgraph;

   SCIP_CALL( tbottleneckInit(scip, subgraph, subsoledges, &(termcomp->subsolbottleneck)) );

   return SCIP_OKAY;
}


/** frees */
void reduce_termcompFree(
   SCIP*                 scip,               /**< SCIP data structure */
   TERMCOMP**            termcomp            /**< to initialize */
   )
{
   TERMCOMP* comp;
   comp = *termcomp;

   assert(comp->nodemap_subToOrg && comp->edgemap_subToOrg);


   if( comp->subsolbottleneck )
      tbottleneckFree(scip, &(comp->subsolbottleneck));
   StpVecFree(scip, comp->bfsqueue);

   SCIPfreeMemoryArrayNull(scip, &(comp->subsolution));
   SCIPfreeMemoryArray(scip, &(comp->nodemap_subToOrg));
   SCIPfreeMemoryArray(scip, &(comp->edgemap_subToOrg));

   if( comp->subgraph )
   {
      graph_free(scip, &(comp->subgraph), TRUE);
   }

   SCIPfreeMemoryArray(scip, &(comp->nodemap_orgToSub));
   SCIPfreeMemoryArray(scip, &(comp->nodes_mark));

   SCIPfreeMemory(scip, termcomp);
}



/** gets next separator component */
void reduce_termsepaGetNextComp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   TERMSEPAS*            termsepas,          /**< terminal separator store */
   COMPBUILDER*          builder,            /**< builder */
   SCIP_Bool*            compWasFound        /**< was new component found? */
   )
{
   const int* sepaterms = NULL;
   int sinkterm;
   int nsinknodes;
   int nsepaterms = builder->nsepatterms;

   *compWasFound = FALSE;

   assert(nsepaterms >= 2);
   assert(builder->componentnumber >= 0);

   if( builder->componentnumber > builder->maxncompchecks )
   {
      SCIPdebugMessage("maximum number of components reached, stopping sepaDualAscent reductions \n");
      return;
   }

   for( ; nsepaterms <= builder->maxsepasize; nsepaterms++  )
   {
      sepaterms = mincut_termsepasGetNext(nsepaterms, termsepas, &sinkterm, &nsinknodes);

      if( sepaterms != NULL )
      {
         *compWasFound = TRUE;
         break;
      }
   }

   if( *compWasFound )
   {
      const int nsourcenodes = builder->ngraphnodes - nsinknodes;

      assert(nsinknodes > 0);
      assert(nsourcenodes > 0);

      builder->nsepatterms = nsepaterms;
      builder->sepaterms = sepaterms;

      /* NOTE: we want to take the smaller component if possible */
      if( nsinknodes > nsourcenodes && !builder->rootcompIsProcessed )
      {
         const int sourceterm = mincut_termsepasGetSource(termsepas);
         assert(graph_knot_isInRange(g, sourceterm));

         builder->sourceterm = sourceterm;
         builder->ncomponentnodes = nsourcenodes;
         /* NOTE: even if the component turns out to not be promising, we only want to check it once */
         builder->rootcompIsProcessed = TRUE;
      }
      else
      {
         builder->sourceterm = sinkterm;
         builder->ncomponentnodes = nsinknodes;
      }

      SCIPdebugMessage("selecting component: source=%d, |S|=%d, |V'|=%d \n", builder->sourceterm, nsepaterms, builder->ncomponentnodes);
   }
   else
   {
      SCIPdebugMessage("no further components available, stopping sepaDualAscent reductions \n");
   }
}
