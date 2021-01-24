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
#include "reduce.h"
#include "extreduce.h"
#include "mincut.h"
#include "portab.h"
#include "stpvector.h"
#include "scip/scip.h"


#define MARK_SUBNODE 1
#define MARK_SEPARATOR 2
#define SEPARATOR_MAXSIZE 5
#define SEPARATOR_MAXNCHECKS 50


/** separator data needed to build component */
typedef struct terminial_component_initializes
{
   const int*            sepaterms;          /**< separator terminals */
   int                   sourceterm;         /**< source terminal NOTE: we eliminate the associated sub-graph! */
   int                   nsepatterms;        /**< size of separator */
   int                   componentnumber;    /**< */
   int                   ngraphnodes;
} COMPINIT;



/** (extended) terminal component */
typedef struct terminal_separator_component
{
   const COMPINIT*       sepainitializer;    /**< initializer; NON-OWNED */
   GRAPH*                subgraph;           /**< graph for (extended) component */
   int*                  subsolution;        /**< primal solution for (extended) component (CONNECTED/UNKNOWN) */
   int*                  nodemap_orgToSub;   /**< map */
   int*                  nodemap_subToOrg;   /**< map */
   int*                  edgemap_subToOrg;   /**< map */
   STP_Vectype(int)      dfsstack;           /**< DFS */
   STP_Vectype(int)      elim_edges;         /**< edges to eliminate */
   STP_Vectype(int)      elim_nodes;         /**< nodes to pseudo-eliminate */
   int                   subnnodes;
   int                   subnedges;
} TERMCOMP;



/*
 * Local methods
 */



/** initializes */
static
SCIP_RETCODE termcompInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const COMPINIT*       sepainitializer,    /**< initializer */
   TERMCOMP**            termcomp            /**< to initialize */
   )
{
   TERMCOMP* comp;

   SCIP_CALL( SCIPallocMemory(scip, termcomp) );
   comp = *termcomp;

   comp->sepainitializer = sepainitializer;
   comp->subgraph = NULL;
   comp->subsolution = NULL;
   comp->edgemap_subToOrg = NULL;
   comp->nodemap_subToOrg = NULL;
   comp->elim_edges = NULL;
   comp->elim_nodes = NULL;
   comp->dfsstack = NULL;
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

   StpVecFree(scip, comp->elim_nodes);
   StpVecFree(scip, comp->elim_edges);
   StpVecFree(scip, comp->dfsstack);

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


/**  identifies subgraph */
static
void subgraphIdentify(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   TERMCOMP*             termcomp            /**< component */
   )
{
   const COMPINIT* const sepainitializer = termcomp->sepainitializer;
   const int* const sepaterms = sepainitializer->sepaterms;
   int* const nodemap_orgToSub = termcomp->nodemap_orgToSub;
   STP_Vectype(int) dfsstack = NULL;
   int* RESTRICT gmark = g->mark;
   int sub_e = 0;
   int sub_n = 0;
   const int nsepaterms = sepainitializer->nsepatterms;
   const int sourceterm = sepainitializer->sourceterm;

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

      nodemap_orgToSub[term] = sub_n;
      sub_n++;
      sub_e += g->grad[term];
      gmark[term] = MARK_SEPARATOR;
   }

   assert(!gmark[sourceterm]);
   StpVecReserve(scip, dfsstack, 16);
   StpVecPushBack(scip, dfsstack, sourceterm);
   gmark[sourceterm] = MARK_SUBNODE;

   while( StpVecGetSize(dfsstack)  )
   {
      const int k = dfsstack[StpVecGetSize(dfsstack) - 1];
      StpVecPopBack(dfsstack);

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
            StpVecPushBack(scip, dfsstack, head);
            gmark[head] = MARK_SUBNODE;
         }
      }
   }

   assert(termcomp->dfsstack == NULL);
   assert(termcomp->subnnodes == -1);
   assert(termcomp->subnedges == -1);

   /* reserve space for the separator terminal clique */
   sub_e += (nsepaterms) * (nsepaterms - 1);

   termcomp->dfsstack = dfsstack;
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
   const COMPINIT* const sepainitializer = termcomp->sepainitializer;
   const int* const nodemap_orgToSub = termcomp->nodemap_orgToSub;
   int* nodemap_subToOrg;
   int* edgemap_subToOrg;
   const int* const orgmark = orggraph->mark;
   const int* const sepaterms = sepainitializer->sepaterms;
   const int nsepaterms = sepainitializer->nsepatterms;
   const int nnodes_sub = termcomp->subnnodes;
   const int nedges_sub = termcomp->subnedges;
   STP_Vectype(int) dfsstack = termcomp->dfsstack;

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
      assert(nodemap_orgToSub[sepaterms[i]] == orggraph->knots);

      graph_knot_add(subgraph, STP_TERM);
      nodemap_subToOrg[i] = sepaterms[i];
   }

   for( int i = 0; i < StpVecGetSize(dfsstack); i++ )
   {
      const int orgnode = dfsstack[i];
      assert(nodemap_orgToSub[orgnode] == orggraph->knots);

      nodemap_subToOrg[orggraph->knots] = orgnode;
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

   for( int i = 0; i < StpVecGetSize(dfsstack); i++ )
   {
      const int orgnode = dfsstack[i];
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

   assert(!termcomp->nodemap_subToOrg && !termcomp->edgemap_subToOrg);

   termcomp->nodemap_subToOrg = nodemap_subToOrg;
   termcomp->edgemap_subToOrg = edgemap_subToOrg;

   return SCIP_OKAY;
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

/** processes subgraph associated with SEPADA */
static
SCIP_RETCODE processComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   const COMPINIT*       sepainitializer,    /**< terminal separator component initializer */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   int*                  nelims              /**< number of eliminations*/
   )
{
   TERMCOMP* component;

   SCIP_CALL( termcompInit(scip, g, sepainitializer, &component) );
   SCIP_CALL( termcompBuildSubgraphWithSds(scip, g, extperma, component) );

   // compute solution on graph (A & P)

   // compute b on anti-component and add to sub-graph

   // dual-scent, and mark edges, nodes (need nodes, edges)

   // delete edges, save nodes

   termcompFree(scip, &component);

   return SCIP_OKAY;
}


/** gets next separator component */
static
void getNextComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   TERMSEPAS*            termsepas,
   COMPINIT*             sepacomp,           /**< component identifier */
   SCIP_Bool*            compWasFound
   )
{
   const int* sepaterms = NULL;
   int sinkterm;
   int nsinknodes;
   int nsepaterms = sepacomp->nsepatterms;

   *compWasFound = FALSE;

   assert(nsepaterms >= 2);
   assert(sepacomp->componentnumber >= 0);

   if( sepacomp->componentnumber++ > SEPARATOR_MAXNCHECKS )
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
     // int todo; // ADD extra methd to get source from temrinal sepa
      const int nsourcenodes = sepacomp->ngraphnodes - nsinknodes;

      assert(nsinknodes > 0);
      assert(nsourcenodes > 0);

      sepacomp->nsepatterms = nsepaterms;
      sepacomp->sepaterms = sepaterms;

      /* NOTE: we want to take the smaller component */
      if( nsinknodes > nsourcenodes )
      {
         const int sourceterm = g->source;
         assert(0);

         sepacomp->sourceterm = sourceterm;
      }
      else
      {
         sepacomp->sourceterm = sinkterm;
      }

      SCIPdebugMessage("selecting component with source %d and %d separator terminals \n", sepacomp->sourceterm, nsepaterms);
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
   COMPINIT sepacomp = { .sepaterms = NULL, .sourceterm = -1,
                         .nsepatterms = 2, .componentnumber = 0, .ngraphnodes = -1 };
   TERMSEPAS* termsepas;

   assert(scip && g && nelims && extperma);
   *nelims = 0;

   graph_get_nVET(g, &(sepacomp.ngraphnodes), NULL, NULL);

   // todo probably we want to have an array (parameter) to keep the nodes for pseudo-elimination

   SCIP_CALL( mincut_termsepasInit(scip, g, &termsepas) );
   // todo different random seed! g->terms maybe
   SCIP_CALL( mincut_findTerminalSeparators(scip, 1, g, termsepas) );

   for( ;; )
   {
      SCIP_Bool compWasFound;
      getNextComponent(scip, g, termsepas, &sepacomp, &compWasFound);

      if( !compWasFound )
         break;

      SCIP_CALL( processComponent(scip, &sepacomp, g, extperma, nelims) );
   }

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

   SCIP_CALL( extreduce_extPermaInit(scip, extred_fast, g, NULL, &extperma) );
   SCIP_CALL( extreduce_distDataInit(scip, g, 50, useSd, FALSE, &(extperma->distdata_default)) );

   SCIP_CALL( reduce_sepaDualAscentWithExperma(scip, g, extperma, nelims) );

   extreduce_exit(scip, g, &extperma);

   return SCIP_OKAY;
}
