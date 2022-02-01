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

/**@file   reduce_termsepada.c
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
//#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "dualascent.h"
#include "reduce.h"
#include "extreduce.h"
#include "solstp.h"
#include "heur_local.h"
#include "mincut.h"
#include "heur_ascendprune.h"
#include "portab.h"
#include "stpvector.h"
#include "scip/scip.h"

// todo tune values, more or less random

#define COMPONENT_NODESRATIO_MIN 0.01
#define COMPONENT_NODESRATIO_SMALL 0.1
#define COMPONENT_NODESRATIO_MAX 0.5

#define SEPARATOR_MAXSIZE 5
#define SEPARATOR_MAXNCHECKS 75
#define SEPARATOR_MAXSIZE_FAST 4
#define SEPARATOR_MAXNCHECKS_FAST 40
#define SEPARATOR_MINTERMRATIO 0.12
#define SEPARATOR_MINTERMRATIO_FAST 0.005


/*
 * Local methods
 */


/** deletes */
static
void termcompDeleteEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced costs */
   const TERMCOMP*       termcomp,           /**< component */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   int*                  nelims              /**< number of eliminations*/
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   const int subnnodes = graph_get_nNodes(subgraph);
   const PATH* vnoi = redcosts_getNodeToTermsPathsTop(redcostdata);
   const SCIP_Real* cost = redcosts_getEdgeCostsTop(redcostdata);
   const SCIP_Real* pathdist = redcosts_getRootToNodeDistTop(redcostdata);
   const SCIP_Real cutoffbound = redcosts_getCutoffTop(redcostdata);
   const int* const edgemap_subToOrg = termcomp->edgemap_subToOrg;
   const int* result = extperma->result;

   assert(graph_isMarked(g));
   assert(!extperma->solIsValid || result);

   for( int k = 0; k < subnnodes; k++ )
   {
      int e = subgraph->outbeg[k];
      while( e != EAT_LAST )
      {
         const int enext = subgraph->oeat[e];
         const int subhead = subgraph->head[e];
         const int eorg = edgemap_subToOrg[e];

         /* NOTE: we avoid double checking and deletion of artificial edges */
         if( subhead > k && eorg >= 0)
         {
            SCIP_Real redcost = pathdist[k] + cost[e] + vnoi[subhead].dist;

            if( !SCIPisGT(scip, redcost, cutoffbound) )
            {
               e = enext;
               continue;
            }

            redcost = pathdist[subhead] + cost[flipedge(e)] + vnoi[k].dist;

            if( SCIPisGT(scip, redcost, cutoffbound) && !graph_edge_isDeleted(g, eorg) )
            {
#ifdef SCIP_DEBUG
               SCIPdebugMessage("deleting original edge (%f>%f): ", redcost, cutoffbound);
               graph_edge_printInfo(g, eorg);
#endif
               extreduce_edgeRemove(scip, eorg, g, extperma->distdata_default, extperma);
               // todo  todo is that ok??? might also change the primal bound
               graph_edge_del(scip, subgraph, e, FALSE);
               if( extperma->solIsValid && (result[e] == CONNECT || result[flipedge(e)] == CONNECT) )
                  extperma->solIsValid = FALSE;

               (*nelims)++;
            }
         }

         e = enext;
      }
   }



   assert(graph_valid(scip, g));
}



/** marks nodes for pseudo-deletion */
static
void termcompMarkPseudoDelNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced costs */
   const TERMCOMP*       termcomp,           /**< component */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   SCIP_Bool*            pseudoDelNodes      /**< array to mark pseudo deletable nodes  */
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   const int subnnodes = graph_get_nNodes(subgraph);
   const PATH* const vnoi = redcosts_getNodeToTermsPathsTop(redcostdata);
   const SCIP_Real* const pathdist = redcosts_getRootToNodeDistTop(redcostdata);
   const int* const nodemap_subToOrg = termcomp->nodemap_subToOrg;
   const SCIP_Real cutoffbound = redcosts_getCutoffTop(redcostdata);

   for( int k = 0; k < subnnodes; k++ )
   {
      if( !Is_term(subgraph->term[k]) )
      {
         const SCIP_Real redcost = pathdist[k] + vnoi[k].dist + vnoi[k + subnnodes].dist;

         if( SCIPisGT(scip, redcost, cutoffbound) )
         {
            const int orgnode = nodemap_subToOrg[k];
            assert(graph_knot_isInRange(g, orgnode));
            pseudoDelNodes[orgnode] = TRUE;

#ifdef SCIP_DEBUG
            SCIPdebugMessage("marking original node (%f>%f): ", redcost, cutoffbound);
            graph_knot_printInfo(g, orgnode);
#endif
         }
      }
   }
}


/** perform reductions */
static
SCIP_RETCODE termcompReduceWithParams(
   SCIP*                 scip,               /**< SCIP data structure */
   DAPARAMS*             daparams,           /**< dual-ascent parameters */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   TERMCOMP*             termcomp,           /**< component */
   SCIP_Bool*            pseudoDelNodes,     /**< array to mark pseudo deletable nodes or NULL (OUT) */
   int*                  nelims              /**< number of eliminations*/
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   RCPARAMS rcparams = { .cutoff = -1.0, .nLevels = 1, .nCloseTerms = 2, .nnodes = subgraph->knots,
                       .nedges = subgraph->edges, .redCostRoot = daparams->root };
   REDCOST* redcostdata;
   const SCIP_Real subprimal = termcomp->subprimalobj;
   SCIP_Real subdual;

   assert(GE(subprimal, 0.0));

   SCIP_CALL( redcosts_initFromParams(scip, &rcparams, &redcostdata) );

   SCIP_CALL( dualascent_exec(scip, subgraph, NULL, daparams, redcosts_getEdgeCostsTop(redcostdata), &subdual) );

   SCIPdebugMessage("subdual=%f, subprimal=%f \n", subdual, subprimal);

   redcosts_setDualBoundTop(subdual, redcostdata);
   graph_mark(subgraph);
   SCIP_CALL( redcosts_initializeDistancesTop(scip, subgraph, redcostdata) );
   redcosts_setCutoffFromBoundTop(subprimal, redcostdata);

   termcompDeleteEdges(scip, redcostdata, termcomp, g, extperma, nelims);

   if( pseudoDelNodes )
   {
      termcompMarkPseudoDelNodes(scip, redcostdata, termcomp, g, extperma, pseudoDelNodes);
   }

   redcosts_free(scip, &redcostdata);

   return SCIP_OKAY;
}


/** perform reductions */
static
SCIP_RETCODE termcompReduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   TERMCOMP*             termcomp,           /**< component */
   SCIP_Bool*            pseudoDelNodes,     /**< array to mark pseudo deletable nodes or NULL (OUT) */
   int*                  nelims              /**< number of eliminations*/
   )
{
   DAPARAMS daparams = { .addcuts = FALSE, .ascendandprune = FALSE, .root = -1,
           .is_pseudoroot = FALSE, .damaxdeviation = -1.0 };
   const COMPBUILDER* builder = termcomp->builder;

   daparams.root = termcomp->subgraph->source;
   SCIP_CALL( termcompReduceWithParams(scip, &daparams, g, extperma, termcomp, pseudoDelNodes, nelims) );

   if( reduce_compbuilderGetSubNodesRatio(builder) <= COMPONENT_NODESRATIO_SMALL && *nelims > 0 )
   {
      const int sepaterm = builder->sepaterms[0];
      assert(daparams.root != termcomp->nodemap_orgToSub[sepaterm]);
      assert(graph_knot_isInRange(termcomp->subgraph, termcomp->nodemap_orgToSub[sepaterm]));

      daparams.root = termcomp->nodemap_orgToSub[sepaterm];

      SCIP_CALL( termcompReduceWithParams(scip, &daparams, g, extperma, termcomp, pseudoDelNodes, nelims) );
   }

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

   // todo maybe deactivate if still too expensive...
   if( reduce_compbuilderGetSubNodesRatio(termcomp->builder) <= COMPONENT_NODESRATIO_SMALL )
      SCIP_CALL( SCIPStpHeurLocalRun(scip, subgraph, subsol) );

   termcomp->subsolution = subsol;
   termcomp->subprimalobj = solstp_getObj(subgraph, subsol, 0.0);
   SCIPfreeBufferArray(scip, &redcosts);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("primal sub-sol value: %f \n", solstp_getObj(subgraph, subsol, 0.0));
#endif

   SCIP_CALL( reduce_termcompInitTbottleneck(scip, subsol, termcomp) );

   return SCIP_OKAY;
}


/** promising to perform reductions on given component? */
static
SCIP_Bool termcompIsPromising(
   const GRAPH*          g,                  /**< graph data structure */
   const COMPBUILDER*    builder             /**< terminal separator component initializer */
   )
{
   const SCIP_Real noderatio = reduce_compbuilderGetSubNodesRatio(builder);

   assert(GT(noderatio, 0.0));
   SCIPdebugMessage(" noderatio=%f \n", noderatio);

   /* NOTE: we allow a few components regardless of their size */
   if( builder->componentnumber < 10 && noderatio < COMPONENT_NODESRATIO_MAX )
   {

      SCIPdebugMessage("...component is promising \n");

      return TRUE;
   }
   else
   {
      if( noderatio < COMPONENT_NODESRATIO_MAX )
      {
         SCIPdebugMessage("...component is promising \n");
         return TRUE;
      }
   }

   // todo correct currently, we always give true, but seems to be no problem....

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
   SCIP_Bool*            pseudoDelNodes,     /**< array to mark pseudo deletable nodes or NULL (OUT) */
   int*                  nelims              /**< number of eliminations*/
   )
{
   TERMCOMP* termcomp;
   SCIP_Bool success;

   SCIPdebugMessage("component nodes=%d \n", builder->ncomponentnodes);

   SCIP_CALL( reduce_termcompInit(scip, g, builder, &termcomp) );
   SCIP_CALL( reduce_termcompBuildSubgraphWithSds(scip, g, extperma, termcomp) );
   SCIP_CALL( termcompComputeSubgraphSol(scip, termcomp) );
   SCIP_CALL( reduce_termcompChangeSubgraphToBottleneck(scip, g, termcomp, &success) );

   assert(success || *nelims > 0);

   /* NOTE: we might fail because the separator is not connected anymore one one side */
   if( success )
   {
      SCIP_CALL( termcompReduce(scip, g, extperma, termcomp, pseudoDelNodes, nelims) );
   }

   reduce_termcompFree(scip, &termcomp);

   return SCIP_OKAY;
}


/** initializes helpers */
static
SCIP_RETCODE initHelpers(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const EXTPERMA*       extperma,           /**< extension data */
   COMPBUILDER**         builder,            /**< to initialize */
   TERMSEPAS**           termsepas           /**< to initialize */
   )
{
   int mincheckbound;
   int maxsepasize;
   int maxncompchecks;

   if( extperma->mode == extred_fast )
   {
      mincheckbound = (int) (SEPARATOR_MINTERMRATIO_FAST * g->terms);
      maxsepasize = SEPARATOR_MAXSIZE_FAST;
      maxncompchecks = SEPARATOR_MAXNCHECKS_FAST;
   }
   else
   {
      mincheckbound = (int) (SEPARATOR_MINTERMRATIO * g->terms);
      maxsepasize = SEPARATOR_MAXSIZE;
      maxncompchecks = SEPARATOR_MAXNCHECKS;
   }

   if( maxncompchecks < mincheckbound )
   {
      SCIPdebugMessage("update nChecks %d->%d \n", maxncompchecks, mincheckbound);
      maxncompchecks = mincheckbound;
   }

   /* NOTE: we want to allow a few more terminal separators to be able to choose small ones */
   SCIP_CALL( mincut_termsepasInit(scip, g, (int) (1.5 * maxncompchecks), maxsepasize, termsepas) );
   SCIP_CALL( reduce_compbuilderInit(scip, g, builder) );

   (*builder)->maxncompchecks = maxncompchecks;
   (*builder)->maxsepasize = maxsepasize;

   return SCIP_OKAY;
}

/** frees helper */
static
void freeHelpers(
   SCIP*                 scip,               /**< SCIP data structure */
   COMPBUILDER**         builder,            /**< to initialize */
   TERMSEPAS**           termsepas           /**< to initialize */
   )
{
   reduce_compbuilderFree(scip, builder);
   mincut_termsepasFree(scip, termsepas);
}


/*
 * Interface methods
 */


/** terminal-separator/dual-ascent reduction method given extperma
 * todo maybe we want to also give terminal separators to be able to reuse them? */
SCIP_RETCODE reduce_termsepaDaWithExperma(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   EXTPERMA*             extperma,           /**< extension data */
   SCIP_Bool*            pseudoDelNodes,     /**< array to mark pseudo deletable nodes or NULL (OUT) */
   int*                  nelims              /**< number of eliminations (OUT)*/
   )
{
   COMPBUILDER* builder;
   TERMSEPAS* termsepas;

   assert(scip && g && nelims && extperma);
   *nelims = 0;

   if( pseudoDelNodes )
      BMSclearMemoryArray(pseudoDelNodes, g->knots);

   if( g->terms == 1 )
   {
      return SCIP_OKAY;
   }

   SCIP_CALL( initHelpers(scip, g, extperma, &builder, &termsepas) );
   SCIP_CALL( mincut_findTerminalSeparators(scip, extperma->randnumgen, g, termsepas) );

   for( ;; )
   {
      SCIP_Bool compWasFound;

      reduce_termsepaGetNextComp(scip, g, termsepas, builder, &compWasFound);

      if( !compWasFound )
         break;

      if( !termcompIsPromising(g, builder) )
         continue;

      SCIP_CALL( processComponent(scip, builder, g, extperma, pseudoDelNodes, nelims) );

      builder->componentnumber++;
   }

   freeHelpers(scip, &builder, &termsepas);

   return SCIP_OKAY;
}


/** terminal-separator/dual-ascent reduction method
 *  NOTE: intended to be used for debugging and unit tests, since creation of SD data is expensive
 *  and is best combined with extended reduction techniques */
SCIP_RETCODE reduce_termsepaDa(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  nelims              /**< number of eliminations*/
   )
{
   EXTPERMA* extperma;
   const SCIP_Bool useSd = TRUE;

   assert(scip && g && nelims);

   if( g->terms == 1 )
   {
      *nelims = 0;
      return SCIP_OKAY;
   }

   SCIP_CALL( graph_init_dcsr(scip, g) );
   SCIP_CALL( extreduce_extPermaInit(scip, extred_fast, g, NULL, &extperma) );
   SCIP_CALL( extreduce_distDataInit(scip, g, 50, useSd, FALSE, &(extperma->distdata_default)) );
   SCIP_CALL( reduce_sdRepairSetUp(scip, g, extperma->distdata_default->sdistdata) );

   SCIP_CALL( reduce_termsepaDaWithExperma(scip, g, extperma, NULL, nelims) );

   /* NOTE: also frees DCSR */
   extreduce_exit(scip, g, &extperma);

   return SCIP_OKAY;
}
