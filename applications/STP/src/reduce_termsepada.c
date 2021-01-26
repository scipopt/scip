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

/**@file   reduce_sepada.c
 * @brief  several terminal-separator/dual-ascent based reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements terminal-separator based reduction techniques for several Steiner problems that
 * use dual-ascent for reductions
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "dualascent.h"
#include "reduce.h"
#include "extreduce.h"
#include "solstp.h"
#include "mincut.h"
#include "heur_ascendprune.h"
#include "portab.h"
#include "stpvector.h"
#include "scip/scip.h"

#define MARK_NONACTIVE 0
#define MARK_SUBNODE   1
#define MARK_SEPARATOR 2
#define COMPONENT_MINNODESRATIO 0.01
#define COMPONENT_MAXNODESRATIO 0.5
#define SEPARATOR_MAXSIZE 5
#define SEPARATOR_MAXNCHECKS 50


/** separator data needed to build component */
typedef struct terminial_component_initializes
{
   SCIP_Real*            nodes_bdist;        /**< bottleneck computation distance for each node, always reset to -1.0 */
   const int*            sepaterms;          /**< separator terminals NON OWNED */
   int                   sourceterm;         /**< source terminal NOTE: we eliminate the associated sub-graph! */
   int                   nsepatterms;        /**< size of separator */
   int                   ncomponentnodes;    /**< NOTE: possibly overestimate */
   int                   componentnumber;    /**< number of component (0,1,...)*/
   int                   ngraphnodes;        /**< number of nodes of underlying graph, not counting degree 0 nodes */
   SCIP_Bool             rootcompIsProcessed;/**< already processed root component? */
} COMPBUILDER;



/** (extended) terminal component */
typedef struct terminal_separator_component
{
   COMPBUILDER*          builder;            /**< initializer; NON-OWNED */
   GRAPH*                subgraph;           /**< graph for (extended) component */
   int*                  subsolution;        /**< primal solution for (extended) component (CONNECTED/UNKNOWN) */
   int*                  nodemap_orgToSub;   /**< map */
   int*                  nodemap_subToOrg;   /**< map */
   int*                  edgemap_subToOrg;   /**< map */
   STP_Vectype(int)      bfsqueue;           /**< queue for BFS */
   int                   subnnodes;
   int                   subnedges;
} TERMCOMP;


/*
 * Local methods
 */


/** initializes */
static
SCIP_RETCODE compbuilderInit(
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

   SCIP_CALL( SCIPallocMemoryArray(scip, &(builder->nodes_bdist), nnodes) );
   graph_get_nVET(g, &(builder->ngraphnodes), NULL, NULL);

   for( int i = 0; i < nnodes; i++ )
   {
      builder->nodes_bdist[i] = -1.0;
   }

   return SCIP_OKAY;
}


/** frees */
static
void compbuilderFree(
   SCIP*                 scip,               /**< SCIP data structure */
   COMPBUILDER**         compbuilder         /**< to free */
   )
{
   COMPBUILDER* builder;
   builder = *compbuilder;

   SCIPfreeMemoryArray(scip, &(builder->nodes_bdist));
   SCIPfreeMemory(scip, compbuilder);
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
   int* RESTRICT gmark = g->mark;
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
   EXTPERMA*             extperma,           /**< extension data */
   TERMCOMP*             termcomp            /**< component */
   )
{
   GRAPH* subgraph;
   DISTDATA* const distdata = extperma->distdata_default;
   const COMPBUILDER* const builder = termcomp->builder;
   const int* const nodemap_orgToSub = termcomp->nodemap_orgToSub;
   int* nodemap_subToOrg;
   int* edgemap_subToOrg;
   const int* const orgmark = orggraph->mark;
   const int* const sepaterms = builder->sepaterms;
   const int nsepaterms = builder->nsepatterms;
   const int nnodes_sub = termcomp->subnnodes;
   const int nedges_sub = termcomp->subnedges;
   STP_Vectype(int) bfsqueue = termcomp->bfsqueue;

   assert(nsepaterms >= 2 && nsepaterms < orggraph->terms);
   assert(nnodes_sub > 0 && nedges_sub > 0);
   assert(nodemap_orgToSub && sepaterms && distdata);
   assert(!termcomp->subgraph);

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
         const SCIP_Real sd = extreduce_distDataGetSdDouble(scip, orggraph, orgsepaterm, orgsepaterm2, distdata);
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

   assert(graph_knot_isInRange(orggraph, builder->sourceterm));
   assert(graph_knot_isInRange(orggraph, nodemap_orgToSub[builder->sourceterm]));

   subgraph->source = nodemap_orgToSub[builder->sourceterm];
   subgraph->stp_type = orggraph->stp_type;

   SCIP_CALL( graph_path_init(scip, subgraph) );

   assert(!termcomp->nodemap_subToOrg && !termcomp->edgemap_subToOrg);

   termcomp->nodemap_subToOrg = nodemap_subToOrg;
   termcomp->edgemap_subToOrg = edgemap_subToOrg;

   return SCIP_OKAY;
}




/** initializes */
static
SCIP_RETCODE termcompInit(
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

   SCIP_CALL( SCIPallocMemoryArray(scip, &(comp->nodemap_orgToSub), g->knots) );

#ifndef NDEBUG
   for( int i = 0; i < g->knots; i++ )
      comp->nodemap_orgToSub[i] = UNKNOWN;
#endif

   return SCIP_OKAY;
}


/** frees */
static
void termcompFree(
   SCIP*                 scip,               /**< SCIP data structure */
   TERMCOMP**            termcomp            /**< to initialize */
   )
{
   TERMCOMP* comp;
   comp = *termcomp;

   StpVecFree(scip, comp->bfsqueue);

   SCIPfreeMemoryArrayNull(scip, &(comp->subsolution));
   SCIPfreeMemoryArray(scip, &(comp->nodemap_subToOrg));
   SCIPfreeMemoryArray(scip, &(comp->edgemap_subToOrg));

   if( comp->subgraph )
   {
      graph_free(scip, &(comp->subgraph), TRUE);
   }

   SCIPfreeMemoryArray(scip, &(comp->nodemap_orgToSub));

   SCIPfreeMemory(scip, termcomp);
}



/** builds extended subgraph with SD weighted edges between terminal-separator nodes */
static
SCIP_RETCODE termcompBuildSubgraphWithSds(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   TERMCOMP*             termcomp            /**< component */
   )
{
   subgraphIdentify(scip, g, termcomp);
   SCIP_CALL( subgraphBuild(scip, g, extperma, termcomp) );

   return SCIP_OKAY;
}


/** gets extended bottleneck distances */
static
SCIP_Real termcompGetExtBottleneckDist(
   const GRAPH*          g,                  /**< graph data structure */
   int                   term1,              /**< terminal */
   int                   term2,              /**< terminal */
   int                   mstsource,          /**< source node of minimum spanning tree */
   PATH*                 mst,                /**< minimum spanning tree */
   SCIP_Real* RESTRICT   mstsdist            /**< node distance helper */
   )
{
   SCIP_Real sdist = 0.0;
   int tempnode = term1;

   assert(Is_term(g->term[term1]));
   assert(Is_term(g->term[term2]));
   assert(term1 != term2);

   mstsdist[tempnode] = 0.0;

   while( tempnode != mstsource )
   {
      const int ne = mst[tempnode].edge;

      assert(g->head[ne] == tempnode);
      tempnode = g->tail[ne];

      if( g->cost[ne] > sdist )
         sdist = g->cost[ne];

      mstsdist[tempnode] = sdist;
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
      sdist = 0.0;
   }

   while( tempnode != mstsource )
   {
      const int ne = mst[tempnode].edge;
      tempnode = g->tail[ne];

      if( g->cost[ne] > sdist )
         sdist = g->cost[ne];

      /* already visited? */
      if( mstsdist[tempnode] > -0.5 )
      {
         if( mstsdist[tempnode] > sdist )
            sdist = mstsdist[tempnode];
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

   return sdist;
}


/** changes weights between terminal separator nodes to bottlenecks */
static
SCIP_RETCODE termcompChangeSubgraphToBottleneck(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   TERMCOMP*             termcomp            /**< component */
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   PATH* mst;
   COMPBUILDER* const builder = termcomp->builder;
   const int* const subsol = termcomp->subsolution; // maybe also have a path for the subsolution? and modify it?
   const int* const sepaterms = builder->sepaterms;
   const int* const nodemap_orgToSub = termcomp->nodemap_orgToSub;
   const int* const nodemap_subToOrg = termcomp->nodemap_subToOrg;
   const int nsepaterms = builder->nsepatterms;
   const int nnodes = graph_get_nNodes(g);
   const int mstroot = sepaterms[0];

   /* mark the anti-component */
   for( int i = 0; i < nnodes; i++ )
      g->mark[i] = (g->mark[i] != MARK_SUBNODE);

   SCIPallocBufferArray(scip, &mst, nnodes);
   graph_path_exec(scip, g, MST_MODE, mstroot, g->cost, mst);

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
            const SCIP_Real b = termcompGetExtBottleneckDist(g, term, head, mstroot, mst, builder->nodes_bdist);

            subgraph->cost[e] = b;
            subgraph->cost[flipedge(e)] = b;

            SCIPdebugMessage("%d->%d setting b=%f \n", term_sub, head_sub, b);

            // using also primal solution! and unmark primal solution
         }
      }
   }

   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}


/** deletes */
static
void termcompDeleteEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   REDCOST*              redcostdata,        /**< reduced costs */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   TERMCOMP*             termcomp            /**< component */
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   const int subnnodes = graph_get_nNodes(subgraph);
   const PATH* vnoi = redcosts_getNodeToTermsPathsTop(redcostdata);
   const SCIP_Real* cost = redcosts_getEdgeCostsTop(redcostdata);
   const SCIP_Real* pathdist = redcosts_getRootToNodeDistTop(redcostdata);
   const SCIP_Real cutoffbound = redcosts_getCutoffTop(redcostdata);
   const int* const edgemap_subToOrg = termcomp->edgemap_subToOrg;

   graph_mark(g);

   for( int k = 0; k < subnnodes; k++ )
   {
      for( int e = subgraph->outbeg[k]; e != EAT_LAST; e = subgraph->oeat[e] )
      {
         const int subhead = subgraph->head[e];

         /* NOTE: we avoid double checking and deletion of artificial edges */
         if( subhead > k && edgemap_subToOrg[e] >= 0)
         {
            SCIP_Real redcost = pathdist[k] + cost[e] + vnoi[subhead].dist;

            if( !SCIPisGT(scip, redcost, cutoffbound) )
               continue;

           // graph_edge_printInfo(g, edgemap_subToOrg[e]);
          //  printf("%f %f %f \n", pathdist[k], cost[e], vnoi[subhead].dist);
            redcost = pathdist[subhead] + cost[flipedge(e)] + vnoi[k].dist;

            if( SCIPisGT(scip, redcost, cutoffbound) )
            {
#ifdef SCIP_DEBUG
               SCIPdebugMessage("deleting original edge (%f>%f): ", redcost, cutoffbound);
               graph_edge_printInfo(g, edgemap_subToOrg[e]);
#endif
               extreduce_edgeRemove(scip, edgemap_subToOrg[e], g, extperma->distdata_default, extperma);
            }
         }
      }
   }
}


/** perform reductions */
static
SCIP_RETCODE termcompReduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   TERMCOMP*             termcomp            /**< component */
   )
{
   DAPARAMS daparams = { .addcuts = FALSE, .ascendandprune = FALSE, .root = -1,
           .is_pseudoroot = FALSE, .damaxdeviation = -1.0 };
   GRAPH* subgraph = termcomp->subgraph;
   RCPARAMS rcparams = { .cutoff = -1.0, .nLevels = 1, .nCloseTerms = 2, .nnodes = subgraph->knots,
                       .nedges = subgraph->edges, .redCostRoot = subgraph->source };
   REDCOST* redcostdata;
   const int* const subsol = termcomp->subsolution; // maybe also have a path for the subsolution? and modify it?
   const SCIP_Real subprimal = solstp_getObj(subgraph, subsol, 0.0);
   SCIP_Real subdual;

   SCIP_CALL( redcosts_initFromParams(scip, &rcparams, &redcostdata) );

   daparams.root = subgraph->source;
   SCIP_CALL( dualascent_exec(scip, subgraph, NULL, &daparams, redcosts_getEdgeCostsTop(redcostdata), &subdual) );

   SCIPdebugMessage("subdual=%f, subprimal=%f \n", subdual, subprimal);

   redcosts_setDualBoundTop(subdual, redcostdata);
   graph_mark(subgraph);
   SCIP_CALL( redcosts_initializeDistancesTop(scip, subgraph, redcostdata) );
   redcosts_setCutoffFromBoundTop(subprimal, redcostdata);

   termcompDeleteEdges(scip, redcostdata, g, extperma, termcomp);

   // todo check for nodes...don't delete
   {
      int todo; // check nodes and also run dual ascent twice, once without guiding solution, once with
                //...or maybe run with different root if component is not big!
   }

   redcosts_free(scip, &redcostdata);

   return SCIP_OKAY;
}


/** computes primal solution on subgraph */
static
SCIP_RETCODE termcompComputeSubgraphSol(
   SCIP*                 scip,               /**< SCIP data structure */
   TERMCOMP*             termcomp            /**< component */
   )
{
   DAPARAMS daparams = { .addcuts = FALSE, .ascendandprune = FALSE, .root = -1,
           .is_pseudoroot = FALSE, .damaxdeviation = -1.0 };

   GRAPH* subgraph = termcomp->subgraph;
   SCIP_Real* redcosts;
   int* subsol;
   const int* const nodemap_orgToSub = termcomp->nodemap_orgToSub;
   SCIP_Real dualobjval;
   SCIP_Bool success;
   const int nedges = graph_get_nEdges(subgraph);
   const int sourceterm_org = termcomp->builder->sourceterm;

   assert(nedges >= 2);
   assert(!termcomp->subsolution);
   assert(nodemap_orgToSub);
   assert(sourceterm_org >= 0);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(subsol), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(redcosts), nedges) );

   daparams.root = nodemap_orgToSub[sourceterm_org];
   SCIP_CALL( dualascent_exec(scip, subgraph, NULL, &daparams, redcosts, &dualobjval) );

   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, NULL, subgraph, redcosts, subsol,
         daparams.root, &success, FALSE));
   assert(success);

   termcomp->subsolution = subsol;
   SCIPfreeBufferArray(scip, &redcosts);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("primal sub-sol value: %f \n", solstp_getObj(subgraph, subsol, 0.0));
#endif

   return SCIP_OKAY;
}


/** promising to perform reductions on given component? */
static
SCIP_Bool termcompIsPromising(
   const GRAPH*          g,                  /**< graph data structure */
   const COMPBUILDER*    builder             /**< terminal separator component initializer */
   )
{
   /* NOTE: we allow a few components regardless of their size */
   if( builder->componentnumber < 10 )
   {
      SCIPdebugMessage("...component is promising \n");
      return TRUE;
   }
   else
   {
      const SCIP_Real noderatio = builder->ncomponentnodes / builder->ngraphnodes;
      assert(GT(noderatio, 0.0));

      SCIPdebugMessage(" noderatio=%f \n", noderatio);

      if( COMPONENT_MINNODESRATIO < noderatio && noderatio < COMPONENT_MAXNODESRATIO )
      {
         SCIPdebugMessage("...component is promising \n");
         return TRUE;
      }
   }

   SCIPdebugMessage("...component is NOT promising! \n");

   return TRUE;
}


/** processes subgraph associated with SEPADA */
static
SCIP_RETCODE processComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   COMPBUILDER*          builder,            /**< terminal separator component initializer */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   int*                  nelims              /**< number of eliminations*/
   )
{
   TERMCOMP* termcomp;

   SCIP_CALL( termcompInit(scip, g, builder, &termcomp) );
   SCIP_CALL( termcompBuildSubgraphWithSds(scip, g, extperma, termcomp) );
   SCIP_CALL( termcompComputeSubgraphSol(scip, termcomp) );
   SCIP_CALL( termcompChangeSubgraphToBottleneck(scip, g, termcomp) );

   // todo, also mark pseudo-eliminaitons
   SCIP_CALL( termcompReduce(scip, g, extperma, termcomp) );

   termcompFree(scip, &termcomp);

   return SCIP_OKAY;
}


/** gets next separator component */
static
void getNextComponent(
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

   if( builder->componentnumber > SEPARATOR_MAXNCHECKS )
   {
      SCIPdebugMessage("maximum number of components reached, stopping sepaDualAscent reductions \n");
      return;
   }

   for( ; nsepaterms <= SEPARATOR_MAXSIZE; nsepaterms++  )
   {
      // todo sinknodes should be computed correctly by the routine!
      // todo!!!
      // todo!!!
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


/*
 * Interface methods
 */


/** terminal-separator/dual-ascent reduction method given extperma
 * todo maybe we want to also give terminal separators to be able to reuse them? */
SCIP_RETCODE reduce_sepaDualAscentWithExperma(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   int*                  nelims              /**< number of eliminations*/
   )
{
   COMPBUILDER* builder;
   TERMSEPAS* termsepas;

   assert(scip && g && nelims && extperma);
   *nelims = 0;

   // todo probably we want to have an array (parameter) to keep the nodes for pseudo-elimination

   SCIP_CALL( mincut_termsepasInit(scip, g, &termsepas) );
   SCIP_CALL( compbuilderInit(scip, g, &builder) );

   // todo different random seed! g->terms maybe
   SCIP_CALL( mincut_findTerminalSeparators(scip, 1, g, termsepas) );

   for( ;; )
   {
      SCIP_Bool compWasFound;
      getNextComponent(scip, g, termsepas, builder, &compWasFound);

      if( !compWasFound )
         break;

      if( !termcompIsPromising(g, builder) )
         continue;

      SCIP_CALL( processComponent(scip, builder, g, extperma, nelims) );

      builder->componentnumber++;
   }

   compbuilderFree(scip, &builder);
   mincut_termsepasFree(scip, &termsepas);


   return SCIP_OKAY;
}


/** terminal-separator/dual-ascent reduction method
 *  NOTE: intended to be used for debugging and unit tests, since creation of SD data is expensive
 *  and is best combined with extended reduction techniques */
SCIP_RETCODE reduce_sepaDualAscent(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  nelims              /**< number of eliminations*/
   )
{
   EXTPERMA* extperma;
   const SCIP_Bool useSd = TRUE;

   assert(scip && g && nelims);

   SCIP_CALL( graph_init_dcsr(scip, g) );
   SCIP_CALL( extreduce_extPermaInit(scip, extred_fast, g, NULL, &extperma) );
   SCIP_CALL( extreduce_distDataInit(scip, g, 50, useSd, FALSE, &(extperma->distdata_default)) );
   SCIP_CALL( reduce_sdRepairSetUp(scip, g, extperma->distdata_default->sdistdata) );

   SCIP_CALL( reduce_sepaDualAscentWithExperma(scip, g, extperma, nelims) );

   /* NOTE: also frees DCSR */
   extreduce_exit(scip, g, &extperma);

   return SCIP_OKAY;
}
