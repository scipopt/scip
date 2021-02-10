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

/**@file   reduce_termsepafull.c
 * @brief  several terminal-separator/exact solution reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements terminal-separator based reduction techniques for Steiner tree problems.
 * These techniques require to solve subproblems to optimality. See reduce_termsepada.c for related,
 * but heuristic methods (based on dual-ascent).
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
#include "solstp.h"
#include "reduce.h"
#include "mincut.h"
#include "portab.h"
#include "substpsolver.h"
#include "extreduce.h"
#include "stpvector.h"
#include "scip/scip.h"

#define SEPARATOR_MAXSIZE 4
#define SEPARATOR_MAXNCHECKS 75
#define SEPARATOR_MINTERMRATIO 0.1
#define COMPONENT_MAXNODESRATIO_2SEPA 0.5
#define COMPONENT_MAXNODESRATIO_3SEPA 0.05
#define COMPONENT_MAXNODESRATIO_4SEPA 0.025

#define COMPONENT_MAXNODESRATIO_1CANDS 0.8
#define COMPONENT_MAXNODESRATIO_2CANDS 0.4
#define COMPONENT_MAXNODESRATIO_3CANDS 0.2
#define COMPONENT_MAXNODESRATIO_4CANDS 0.1
#define COMPONENT_MAXNODESRATIO_5PLUSCANDS 0.02


/*
 * Data structures
 */

/** full terminal separator reduction data structure */
typedef struct terminal_separator_full
{
   STP_Vectype(int*)     subsols;            /**< stores all required solutions to the current component */
   DISTDATA*             subdistdata;        /**< distance data on sub-problem */
   GRAPH*                subgraph;           /**< component graph NON-OWNED */
   SCIP_Real*            subgraph_bcosts;    /**< edge costs with bottleneck distances */
   SCIP_Real*            subgraph_orgcosts;  /**< original edge costs, with artificial edges having cost FARAWAY */
   STP_Vectype(int*)     solcands_sepaedges; /**< bottleneck edges (w.r.t. subgraph) for each solution candidate */
   int                   nsolcands;          /**< number of solution candidates */
} TSEPAFULL;



/*
 * Local methods
 */


/** sets up distance data for subgraph */
static
SCIP_RETCODE sepafullInitDistdata(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                subgraph,           /**< graph data structure */
   TSEPAFULL*            tsepafull           /**< to initialize for */
   )
{
   const SCIP_Bool useSd = TRUE;
   const SCIP_Bool useBias = TRUE;

   assert(tsepafull);
   assert(!tsepafull->subdistdata);

   SCIP_CALL( graph_init_dcsr(scip, subgraph) );
   SCIP_CALL( extreduce_distDataInit(scip, subgraph, 50, useSd, useBias, &(tsepafull->subdistdata)) );
   graph_free_dcsr(scip, subgraph);

   return SCIP_OKAY;
}


/** initializes */
static
SCIP_RETCODE sepafullInit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                subgraph,           /**< component graph */
   TSEPAFULL**           tsepafull           /**< to initialize */
   )
{
   TSEPAFULL* tfull;

   assert(scip && subgraph);

   SCIP_CALL( SCIPallocMemory(scip, tsepafull) );
   tfull = *tsepafull;

   tfull->subsols = NULL;
   tfull->subgraph = subgraph;
   tfull->solcands_sepaedges = NULL;
   tfull->nsolcands = 0;
   tfull->subgraph_bcosts = NULL;
   tfull->subgraph_orgcosts = NULL;
   tfull->subdistdata = NULL;

   return SCIP_OKAY;
}


/** frees */
static
void sepafullFree(
   SCIP*                 scip,               /**< SCIP data structure */
   TSEPAFULL**           tsepafull           /**< to initialize */
   )
{
   TSEPAFULL* tfull;
   tfull = *tsepafull;

   for( int i = StpVecGetSize(tfull->subsols) - 1; i >= 0; i-- )
   {
      assert(tfull->subsols[i]);
      SCIPfreeMemoryArray(scip, &(tfull->subsols[i]));
   }
   StpVecFree(scip, tfull->subsols);

   for( int i = StpVecGetSize(tfull->solcands_sepaedges) - 1; i >= 0; i-- )
   {
      StpVecFree(scip, tfull->solcands_sepaedges[i]);
   }
   StpVecFree(scip, tfull->solcands_sepaedges);

   SCIPfreeMemoryArrayNull(scip, &(tfull->subgraph_bcosts));
   SCIPfreeMemoryArrayNull(scip, &(tfull->subgraph_orgcosts));

   if( tfull->subdistdata )
   {
      extreduce_distDataFree(scip, tfull->subgraph, &(tfull->subdistdata));
   }

   SCIPfreeMemory(scip, tsepafull);
}


/** builds solution candidates (and tries to rule-out some already) */
static
SCIP_RETCODE sepafullBuildSolcands(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                orggraph,           /**< graph data structure */
   TERMCOMP*             termcomp,           /**< component */
   TSEPAFULL*            tsepafull,          /**< sepa */
   SCIP_Bool*            success             /**< built? */
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   const COMPBUILDER* const builder = termcomp->builder;
  // const int* const edgemap_subToOrg = termcomp->edgemap_subToOrg;
   const int nsepaterms = builder->nsepatterms;
   const int nsubedges = graph_get_nEdges(subgraph);

   assert(!tsepafull->solcands_sepaedges);
   assert(tsepafull->nsolcands == 0);
   assert(nsubedges >= 2);

   SCIP_CALL( reduce_termcompChangeSubgraphToBottleneck(scip, orggraph, termcomp, success) );

   if( !(*success) )
   {
      printf("problem while building bottleneck distances, aborting \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(tsepafull->subgraph_bcosts), nsubedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tsepafull->subgraph_orgcosts), nsubedges) );

   BMScopyMemoryArray(tsepafull->subgraph_bcosts, subgraph->cost, nsubedges);

   /* NOTE: we need to reset the costs for SDs to be valid */
   reduce_termcompChangeSubgraphToOrgCosts(orggraph, termcomp);
   BMScopyMemoryArray(tsepafull->subgraph_orgcosts, subgraph->cost, nsubedges);

   SCIP_CALL( sepafullInitDistdata(scip, subgraph, tsepafull) );

   // todo here we want to look at each possible combination and rule out
   // for testing we just add one and remove the other
   // todo: also test for ruleout with SDs!

   assert(nsepaterms == 2);

   StpVecPushBack(scip, tsepafull->solcands_sepaedges, NULL);
   tsepafull->nsolcands++;

   {
      SCIP_Bool isRuledOut = FALSE;

      StpVecPushBack(scip, tsepafull->solcands_sepaedges, NULL);

      for( int subsepaterm = 0; subsepaterm < nsepaterms && !isRuledOut; subsepaterm++ )
      {
         for( int e = subgraph->outbeg[subsepaterm]; e != EAT_LAST; e = subgraph->oeat[e] )
         {
            const int head = subgraph->head[e];
            if( head < subsepaterm )
            {
               const SCIP_Real bdist = tsepafull->subgraph_bcosts[e];
               const SCIP_Real sd =
                  extreduce_distDataGetSdDouble(scip, subgraph, subsepaterm, head, tsepafull->subdistdata);
               StpVecPushBack(scip, tsepafull->solcands_sepaedges[tsepafull->nsolcands], e);

               printf("%d->%d bdist=%f, sd=%f \n", subsepaterm, head, bdist, sd);

               if( SCIPisGE(scip, bdist, sd) )
               {
                  isRuledOut = TRUE;
                  break;
               }

               assert(LE(sd, BLOCKED));
            }
         }
      }

      if( isRuledOut )
      {
         StpVecFree(scip, tsepafull->solcands_sepaedges[tsepafull->nsolcands]);
         StpVecPopBack(tsepafull->solcands_sepaedges);
      }
      else
      {
         tsepafull->nsolcands++;
      }

   }


   return SCIP_OKAY;
}


/** promising to do reductions? */
static
SCIP_Bool sepafullSolcandsArePromising(
   const TERMCOMP*       termcomp,           /**< component */
   const TSEPAFULL*      tsepafull           /**< sepa */
   )
{
   const COMPBUILDER* const builder = termcomp->builder;
   SCIP_Real maxratio;
   const SCIP_Real noderatio = reduce_compbuilderGetSubNodesRatio(builder);
   const int nsolcands = tsepafull->nsolcands;

   assert(nsolcands > 0);
   assert(GT(noderatio, 0.0));

   if( nsolcands == 1 )
   {
      maxratio = COMPONENT_MAXNODESRATIO_1CANDS;
   }
   else if( nsolcands == 2 )
   {
      maxratio = COMPONENT_MAXNODESRATIO_2CANDS;
   }
   else if( nsolcands == 3 )
   {
      maxratio = COMPONENT_MAXNODESRATIO_3CANDS;
   }
   else if( nsolcands == 4 )
   {
      maxratio = COMPONENT_MAXNODESRATIO_4SEPA;
   }
   else
   {
      assert(builder->nsepatterms >= 5);
      maxratio = COMPONENT_MAXNODESRATIO_5PLUSCANDS;
   }

   printf("FINAL CHECK: noderatio=%f, maxratio=%f\n", noderatio, maxratio);

   if( noderatio > maxratio )
   {
      printf("...component is finally NOT promising! \n");
      return FALSE;
   }

   printf("...component is finally promising! \n");
   return TRUE;
}


/** solves subproblem */
static
SCIP_RETCODE solveSub(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          subgraph,           /**< graph */
   int*                  subedges_sol        /**< to store solution */
   )
{
   GRAPH* subgraph_copy;
   SUBSTP* substp;
   SCIP_Bool success;

   assert(subgraph);

   SCIP_CALL( graph_copy(scip, subgraph, &subgraph_copy) );

   /* NOTE: subgraph will be moved into substp! */
   SCIP_CALL( substpsolver_init(scip, subgraph_copy, &substp) );
   SCIP_CALL( substpsolver_initHistory(substp) );

   SCIP_CALL( substpsolver_setMute(substp) );
   substpsolver_setProbIsIndependent(substp);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("starting to solve subgraph: ");
   graph_printInfo(subgraph);
#endif

   SCIP_CALL( substpsolver_solve(scip, substp, &success) );
   assert(success);

   /* fix solution in original graph */
   SCIP_CALL( substpsolver_getSolution(substp, subedges_sol) );

   substpsolver_free(scip, &substp);

   return SCIP_OKAY;
}


/** adds solution for candidate, or rules the solution out */
static
SCIP_RETCODE sepafullAddSolForCand(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   cand,               /**< candidate number */
   TERMCOMP*             termcomp,           /**< component */
   TSEPAFULL*            tsepafull           /**< full separation data  */
   )
{
   int* subsol;
   GRAPH* subgraph = termcomp->subgraph;
   STP_Vectype(int) solcands_sepaedges = tsepafull->solcands_sepaedges[cand];
   const SCIP_Real* const subgraph_bcosts = tsepafull->subgraph_bcosts;
   const SCIP_Real* const subgraph_orgcosts = tsepafull->subgraph_orgcosts;

   assert(0 <= cand && cand < tsepafull->nsolcands);
   assert(subgraph_bcosts && subgraph_orgcosts);

   BMScopyMemoryArray(subgraph->cost, subgraph_orgcosts, subgraph->edges);

   SCIP_CALL( SCIPallocMemoryArray(scip, &subsol, subgraph->edges) );

   for( int i = 0; i < StpVecGetSize(solcands_sepaedges); i++ )
   {
      const int subedge = solcands_sepaedges[i];
      assert(graph_edge_isInRange(subgraph, subedge));

      subgraph->cost[subedge] = subgraph_bcosts[subedge];
      subgraph->cost[flipedge(subedge)] = subgraph_bcosts[flipedge(subedge)];

#ifdef SCIP_DEBUG
      SCIPdebugMessage("adding artifical subedge (for solcand %d): ", cand);
      graph_edge_printInfo(subgraph, subedge);
#endif
   }

   SCIP_CALL( solveSub(scip, subgraph, subsol) );
   assert(solstp_isValid(scip, subgraph, subsol));

   for( int i = 0; i < StpVecGetSize(solcands_sepaedges); i++ )
   {
      const int subedge = solcands_sepaedges[i];

      if( subsol[subedge] == UNKNOWN && subsol[flipedge(subedge)] == UNKNOWN )
      {
         SCIPdebugMessage("bottleneck edge %d not used, discarding solution! \n", subedge);
         SCIPfreeMemoryArray(scip, &subsol);
         assert(subsol == NULL);
         break;
      }
   }

   if( subsol )
   {
      StpVecPushBack(scip, tsepafull->subsols, subsol);
   }

   return SCIP_OKAY;
}


/** performs reductions based on solutions to feasible subproblem */
static
SCIP_RETCODE sepafullReduceFromSols(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                orggraph,           /**< graph data structure */
   REDBASE*              redbase,            /**< reduction base */
   TSEPAFULL*            tsepafull,          /**< full separation data  */
   TERMCOMP*             termcomp,           /**< component */
   int*                  nelims              /**< number of eliminations */
   )
{
   const GRAPH* const subgraph = termcomp->subgraph;
   STP_Vectype(int*) subsols = tsepafull->subsols;
   const int* const edgemap_subToOrg = termcomp->edgemap_subToOrg;
   int8_t* edgesolcount;
   SCIP_Real* offset = reduce_solGetOffsetPointer(redbase->redsol);
   const int nsubedges = graph_get_nEdges(subgraph);
   const int nvalidsols = StpVecGetSize(subsols);

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &edgesolcount, nsubedges / 2) );

   assert(nvalidsols > 0);
   assert(edgemap_subToOrg && subsols);

   SCIPdebugMessage("number of valid solution: %d \n", nvalidsols);

   for( int i = 0; i < StpVecGetSize(subsols); i++ )
   {
      const int* const subsol = subsols[i];
      assert(solstp_isValid(scip, subgraph, subsol));

      for( int e = 0; e < nsubedges; e += 2 )
      {
         if( subsol[e] == CONNECT || subsol[e + 1] == CONNECT )
         {
#ifdef SCIP_DEBUG
            SCIPdebugMessage("solution %d (original) edge: ", i);
            if( edgemap_subToOrg[e] != -1  )
               graph_edge_printInfo(orggraph, edgemap_subToOrg[e]);
            else
               printf("artifical \n");
#endif
            edgesolcount[e / 2]++;
         }
         else
         {
            assert(subsol[e] == UNKNOWN && subsol[e + 1] == UNKNOWN);
            edgesolcount[e / 2]--;
         }
      }
   }

   for( int e = 0; e < nsubedges / 2; e++ )
   {
      const int orgedge = edgemap_subToOrg[e * 2];

      if( orgedge == -1 )
         continue;

      assert(graph_edge_isInRange(orggraph, orgedge));

      if( edgesolcount[e] == nvalidsols )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMessage("fix edge ");
         graph_edge_printInfo(orggraph, orgedge);
#endif
         /* NOTE: edge will be automatically contracted later; we avoid trouble other terminal separators */
         *offset += orggraph->cost[orgedge];
         orggraph->cost[orgedge] = 0.0;
         orggraph->cost[flipedge(orgedge)] = 0.0;
         graph_knot_chg(orggraph, orggraph->tail[orgedge], STP_TERM);
         graph_knot_chg(orggraph, orggraph->head[orgedge], STP_TERM);
         (*nelims)++;
      }
      else if( edgesolcount[e] == -nvalidsols )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMessage("delete edge ");
         graph_edge_printInfo(orggraph, orgedge);
#endif
         graph_edge_del(scip, orggraph, orgedge, TRUE);
         (*nelims)++;
      }
   }

   SCIPfreeMemoryArray(scip, &edgesolcount);

   return SCIP_OKAY;
}

/** performs reductions */
static
SCIP_RETCODE sepafullReduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                orggraph,           /**< graph data structure */
   REDBASE*              redbase,            /**< reduction base */
   TSEPAFULL*            tsepafull,          /**< full separation data  */
   TERMCOMP*             termcomp,           /**< component */
   int*                  nelims              /**< number of eliminations*/
   )
{
   const int nsolcands = tsepafull->nsolcands;

   assert(nsolcands > 0);

   for( int i = 0; i < nsolcands; i++ )
   {
      SCIP_CALL( sepafullAddSolForCand(scip, i, termcomp, tsepafull) );
   }

   SCIP_CALL( sepafullReduceFromSols(scip, orggraph, redbase, tsepafull, termcomp, nelims) );

   return SCIP_OKAY;
}


/** promising to perform reductions on given component? */
static
SCIP_Bool termcompIsPromising(
   const GRAPH*          g,                  /**< graph data structure */
   const COMPBUILDER*    builder             /**< terminal separator component initializer */
   )
{
   SCIP_Real maxratio;
   const SCIP_Real noderatio = reduce_compbuilderGetSubNodesRatio(builder);

   assert(GT(noderatio, 0.0));

   if( builder->nsepatterms == 2 )
   {
      maxratio = COMPONENT_MAXNODESRATIO_2SEPA;
   }
   else if( builder->nsepatterms == 3 )
   {
      maxratio = COMPONENT_MAXNODESRATIO_3SEPA;
   }
   else
   {
      assert(builder->nsepatterms == 4);
      maxratio = COMPONENT_MAXNODESRATIO_4SEPA;
   }

   if( builder->nsepatterms != 2 )
   {
      int todo; // deletme
      printf("not promising, sepa != 2 \n");

      return FALSE;
   }

   printf("noderatio=%f, maxratio=%f\n", noderatio, maxratio);

   if( noderatio > maxratio )
   {
      printf("...component is NOT promising! \n");
      return FALSE;
   }

   printf("...component is promising! \n");

   return TRUE;
}


/** processes subgraph associated with builder */
static
SCIP_RETCODE processComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   COMPBUILDER*          builder,            /**< terminal separator component initializer */
   GRAPH*                g,                  /**< graph data structure */
   REDBASE*              redbase,            /**< reduction base */
   int*                  nelims              /**< number of eliminations*/
   )
{
   TSEPAFULL* tsepafull;
   TERMCOMP* termcomp;
   SCIP_Bool success;

#ifdef SCIP_DEBUG
   reduce_compbuilderPrintSeparators(g, builder);
   SCIPdebugMessage("number of component nodes: %d \n", builder->ncomponentnodes);
#endif

   SCIP_CALL( reduce_termcompInit(scip, g, builder, &termcomp) );
   SCIP_CALL( reduce_termcompBuildSubgraph(scip, g, termcomp) );

   /* NOTE: the subgraph does not yet have the bottleneck short-cuts */
   SCIP_CALL( sepafullInit(scip, termcomp->subgraph, &tsepafull) );
   SCIP_CALL( sepafullBuildSolcands(scip, g, termcomp, tsepafull, &success) );

   /* NOTE: we might fail because the separator is not connected anymore on one side */

   if( success && sepafullSolcandsArePromising(termcomp, tsepafull) )
   {
      SCIP_CALL( sepafullReduce(scip, g, redbase, tsepafull, termcomp, nelims) );
   }

   sepafullFree(scip, &tsepafull);
   reduce_termcompFree(scip, &termcomp);

   return SCIP_OKAY;
}


/** initializes helpers */
static
SCIP_RETCODE initHelpers(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   COMPBUILDER**         builder,            /**< to initialize */
   TERMSEPAS**           termsepas           /**< to initialize */
   )
{
   const int mincheckbound = (int) (SEPARATOR_MINTERMRATIO * g->terms);
   const int maxsepasize = SEPARATOR_MAXSIZE;
   int maxncompchecks = SEPARATOR_MAXNCHECKS;

   if( maxncompchecks < mincheckbound )
   {
      SCIPdebugMessage("update nChecks %d->%d \n", maxncompchecks, mincheckbound);
      maxncompchecks = mincheckbound;
   }

   SCIPdebugMessage("max. number of components to be checked: %d \n", maxncompchecks);

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


/** terminal-separator reduction method using optimal (full) solution of sub-problems */
SCIP_RETCODE reduce_termsepaFull(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< solution nodes mark or NULL */
   REDBASE*              redbase,            /**< reduction base */
   int*                  nelims              /**< number of eliminations*/
   )
{
   SCIP_RANDNUMGEN* randnumgen;
   COMPBUILDER* builder;
   TERMSEPAS* termsepas;

   assert(scip && g && nelims && redbase);
   *nelims = 0;

   if( g->terms == 1 )
   {
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, (unsigned) g->terms, TRUE) );
   SCIP_CALL( initHelpers(scip, g, &builder, &termsepas) );
   SCIP_CALL( mincut_findTerminalSeparators(scip, randnumgen, g, termsepas) );

   SCIPdebugMessage("starting termsepFull, with %d separators \n", mincut_termsepasGetNall(termsepas));

   for( ;; )
   {
      SCIP_Bool compWasFound;

      reduce_termsepaGetNextComp(scip, g, termsepas, builder, &compWasFound);

      if( !compWasFound )
         break;

      if( !termcompIsPromising(g, builder) )
         continue;

      SCIP_CALL( processComponent(scip, builder, g, redbase, nelims) );

      builder->componentnumber++;
   }

   /* NOTE: solution edges have been fixed to 0 before */
   if( *nelims > 0 )
      SCIP_CALL( reduce_contract0Edges(scip, g, solnode, TRUE) );

   SCIPfreeRandom(scip, &randnumgen);
   freeHelpers(scip, &builder, &termsepas);

   return SCIP_OKAY;
}
