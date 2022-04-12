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
//#define SCIP_DEBUG

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
#define SEPARATOR_MINEFFTERMRATIO 0.3
#define COMPONENT_MAXNODESRATIO_2SEPA 0.5
#define COMPONENT_MAXNODESRATIO_3SEPA 0.1
#define COMPONENT_MAXNODESRATIO_4SEPA 0.025

#define COMPONENT_MAXNODESRATIO_1CANDS 0.8
#define COMPONENT_MAXNODESRATIO_2CANDS 0.33
#define COMPONENT_MAXNODESRATIO_3CANDS 0.15
#define COMPONENT_MAXNODESRATIO_4CANDS 0.08
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
   STP_Vectype(int*)     solcands_sepapart;  /**< partition of the separator terminals for each solution candidate */
   int                   nsolcands;          /**< number of solution candidates */
} TSEPAFULL;


/** for computing all partitions of a given ground set */
typedef struct bell_parititioner
{
   STP_Vectype(int*)     partitions;         /**< all partitions; each stores as array with mark 0,1,..., according to subset membership
                                                  of each element 0,...,nelems - 1*/
   int                   nelems;             /**< number of elements of ground set */
   int                   npartitions;        /**< number of partitions */
} BPARTITIONS;


/*
 * Local methods
 */


/* gets number of partition groups */
static
int getPartitionNgroups(
   int                   nelems,             /**< number of elements */
   const int*            partition           /**< group id for each element; 0,1,... */
   )
{
   int ngroups = 0;

   assert(nelems > 0);
   assert(partition);

   for( int i = 0; i < nelems; i++ )
   {
      const int group = partition[i];
      assert(0 <= group && group < nelems);

      if( group > ngroups )
         ngroups = group;
   }

   ngroups++;

   SCIPdebugMessage("number of partition groups: %d \n", ngroups);

   return ngroups;
}


/* gets groups sizes */
static
void getPartitionGroupsSizes(
   int                   nelems,             /**< number of elements */
   int                   ngroups,            /**< number of groups */
   const int*            partition,          /**< group id for each element; 0,1,... */
   int*                  groups_sizes        /**< to be filled: size of each group */
   )
{
   assert(nelems > 0);
   assert(partition && groups_sizes);
   assert(ngroups == getPartitionNgroups(nelems, partition));

   BMSclearMemoryArray(groups_sizes, ngroups);

   for( int i = 0; i < nelems; i++ )
   {
      const int group = partition[i];
      assert(0 <= group && group < ngroups);

      groups_sizes[group]++;
   }
}


/** prints top partition if SCIP_DEBUG is defined */
static
void bpartitionsDebugPrintTop(
   const BPARTITIONS*    bparitions          /**< to print for */
   )
{
   assert(bparitions);
   assert(StpVecGetSize(bparitions->partitions) > 0);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("partition: \n");
   for( int i = 0; i < bparitions->nelems; i++ )
   {
      const int* const partition = bparitions->partitions[StpVecGetSize(bparitions->partitions) - 1];
      printf("%d  ", partition[i]);
   }
   printf("\n");
#endif
}


/** recursive add
 *  NOTE: not extremely efficient, only meant for small sizes */
static
SCIP_RETCODE bsubpartAdd(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   subsize,            /**< size of subgroup */
   int                   maxgroupid,         /**< maximum id for partition subset, also marks subset
                                                  which will be further parititioned */
   BPARTITIONS*          bparts              /**< bell partitions */
   )
{

   assert(0 < subsize && subsize <= bparts->nelems);
   assert(maxgroupid >= 0);
   assert(StpVecGetSize(bparts->partitions) > 0);

   if( subsize >= 2 )
   {
      const int* const basepartition = bparts->partitions[StpVecGetSize(bparts->partitions) - 1];
      const int length = subsize - 1;
      const uint32_t powsize = (uint32_t) pow(2.0, length);
      const int nelems = bparts->nelems;
      const int id0 = maxgroupid;
      const int id1 = maxgroupid + 1;

      assert(basepartition);

      for( uint32_t counter = powsize - 1; counter >= 1; counter-- )
      {
         int* partition;
         int pos;
         int subsubsize = 0;

         /* we get the previous partition, and split the group with the highest ID further */
         SCIP_CALL( SCIPallocMemoryArray(scip, &partition, nelems) );
         BMScopyMemoryArray(partition, basepartition, nelems);

         for( pos = 0; pos < nelems; pos++ )
         {
            if( partition[pos] == maxgroupid )
            {
               partition[pos] = id0;
               break;
            }
         }
         assert(pos != nelems);

         for( uint32_t j = 0; j < (uint32_t) length; j++ )
         {
            pos++;
            for( ; pos < nelems; pos++ )
               if( partition[pos] == maxgroupid )
                  break;

            assert(pos != nelems);

            /* Check if jth bit in counter is set */
            if( counter & ((uint32_t) 1 << j) )
            {
               partition[pos] = id1;
               subsubsize++;
            }
            else
            {
               partition[pos] = id0;
            }
         }

         StpVecPushBack(scip, bparts->partitions, partition);
         bpartitionsDebugPrintTop(bparts);

         SCIP_CALL( bsubpartAdd(scip, subsubsize, id1, bparts) );
      }
   }

   return SCIP_OKAY;
}


/** computes partitions */
static
SCIP_RETCODE bpartitionsCompute(
   SCIP*                 scip,               /**< SCIP data structure */
   BPARTITIONS*          bparts              /**< bell partitions */
   )
{
   int* part0;
   const int maxgroupid = 0;
   const int nelems = bparts->nelems;

   assert(nelems >= 2);
   assert(StpVecGetSize(bparts->partitions) == 0);

   SCIP_CALL( SCIPallocMemoryArray(scip, &part0, nelems) );
   for( int i = 0; i < bparts->nelems; i++ )
      part0[i] = maxgroupid;

   StpVecPushBack(scip, bparts->partitions, part0);
   bpartitionsDebugPrintTop(bparts);

   /* start the recursion */
   SCIP_CALL( bsubpartAdd(scip, nelems, maxgroupid, bparts) );

   bparts->npartitions = StpVecGetSize(bparts->partitions);
   assert(bparts->npartitions >= 2);

#ifndef NDEBUG
   if( bparts->nelems == 2 )
   {
      assert(bparts->npartitions == 2);
   }
   else if( bparts->nelems == 3 )
   {
      assert(bparts->npartitions == 5);
   }
   else if( bparts->nelems == 4 )
   {
      assert(bparts->npartitions == 15);
   }
#endif


   return SCIP_OKAY;
}


/** initializes */
static
SCIP_RETCODE bpartitionsInit(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nelems,             /**< number of elements of ground set */
   BPARTITIONS**         bparitions          /**< to initialize */
   )
{
   BPARTITIONS* bparts;

   assert(scip);
   assert(nelems >= 2);

   SCIP_CALL( SCIPallocMemory(scip, bparitions) );
   bparts = *bparitions;

   bparts->partitions = NULL;
   bparts->nelems = nelems;
   bparts->npartitions = 0;

   return SCIP_OKAY;
}


/** frees */
static
void bpartitionsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BPARTITIONS**         bparitions          /**< to initialize */
   )
{
   BPARTITIONS* bparts;
   bparts = *bparitions;

   assert(scip);

   for( int i = StpVecGetSize(bparts->partitions) - 1; i >= 0; i-- )
   {
      SCIPfreeMemoryArrayNull(scip, &(bparts->partitions[i]));
   }
   StpVecFree(scip, bparts->partitions);

   SCIPfreeMemory(scip, bparitions);
}



/** sets up distance data for subgraph */
static
SCIP_RETCODE sepafullInitDistdata(
   SCIP*                 scip,               /**< SCIP data structure */
   TERMCOMP*             termcomp,           /**< component */
   TSEPAFULL*            tsepafull           /**< to initialize for */
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   const SCIP_Bool useSd = TRUE;
   // todo test the impact of bias. Could be activated again
   // if the separator terminals are not used for computing profits
   const SCIP_Bool useBias = FALSE;

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
   tfull->solcands_sepapart = NULL;
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

   for( int i = StpVecGetSize(tfull->solcands_sepapart) - 1; i >= 0; i-- )
   {
      SCIPfreeMemoryArray(scip, &(tfull->solcands_sepapart[i]));
   }
   StpVecFree(scip, tfull->solcands_sepapart);

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

/** add solution candidate edges from given partition */
static
void sepafullAddSingleSolcandEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   partitionidx,       /**< partition index for bpartitions */
   const TERMCOMP*       termcomp,           /**< component */
   BPARTITIONS*          bpartitions,        /**< partitions of separation terminals */
   TSEPAFULL*            tsepafull           /**< full separator */
   )
{
   const GRAPH* subgraph = termcomp->subgraph;
   const COMPBUILDER* const builder = termcomp->builder;
   SCIP_Bool isRuledOut = FALSE;
   const int nsepaterms = builder->nsepatterms;
   const int nsolcands = tsepafull->nsolcands;
   const int* const partition = bpartitions->partitions[partitionidx];

   SCIPdebugMessage("checking next solution candidate (number %d) \n", nsolcands);

   assert(partition);
   assert(StpVecGetSize(tsepafull->solcands_sepaedges) == StpVecGetSize(tsepafull->solcands_sepapart));

   StpVecPushBack(scip, tsepafull->solcands_sepaedges, NULL);
   assert(nsolcands == StpVecGetSize(tsepafull->solcands_sepaedges) - 1);

   for( int subsepaterm = 0; subsepaterm < nsepaterms && !isRuledOut; subsepaterm++ )
   {
      const int groupid = partition[subsepaterm];

      for( int e = subgraph->outbeg[subsepaterm]; e != EAT_LAST; e = subgraph->oeat[e] )
      {
         const int head = subgraph->head[e];

         if( head > subsepaterm )
            continue;

         if( groupid == partition[head] )
         {
            const SCIP_Real bdist = tsepafull->subgraph_bcosts[e];
            const SCIP_Real sd =
               extreduce_distDataGetSdDouble(scip, subgraph, subsepaterm, head, tsepafull->subdistdata);

            StpVecPushBack(scip, tsepafull->solcands_sepaedges[nsolcands], e);

#ifdef SCIP_DEBUG
            SCIPdebugMessage("%d->%d bdist=%f, sd=%f \n",
                  termcomp->nodemap_subToOrg[subsepaterm],
                  termcomp->nodemap_subToOrg[head],
                  bdist, sd);
#endif

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

      SCIPdebugMessage("...ruled-out! \n");
   }
   else
   {
      tsepafull->nsolcands++;

      StpVecPushBack(scip, tsepafull->solcands_sepapart, bpartitions->partitions[partitionidx]);
      bpartitions->partitions[partitionidx] = NULL;

      SCIPdebugMessage("...added! \n");
   }
}

/** computes artificial edges for each solution candidate (and tries to rule-out some candidates already) */
static
SCIP_RETCODE sepafullBuildSolcandsEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   TERMCOMP*             termcomp,           /**< component */
   TSEPAFULL*            tsepafull           /**< full separator */
   )
{
   BPARTITIONS* bpartitions;
   const COMPBUILDER* const builder = termcomp->builder;
   const int nsepaterms = builder->nsepatterms;

   assert(nsepaterms <= 4);
   assert(tsepafull->nsolcands == 0);
   assert(!tsepafull->solcands_sepaedges);

#ifdef SCIP_DEBUG
   for( int i = 0; i < nsepaterms; i++ )
      SCIPdebugMessage("sepa sub to org: %d->%d \n", i, termcomp->nodemap_subToOrg[i]);
#endif

   SCIP_CALL( bpartitionsInit(scip, nsepaterms, &bpartitions) );
   SCIP_CALL( bpartitionsCompute(scip, bpartitions) );

   for( int i = 0; i < bpartitions->npartitions; i++ )
   {
      sepafullAddSingleSolcandEdges(scip, i, termcomp, bpartitions, tsepafull);
   }

   bpartitionsFree(scip, &bpartitions);

   return SCIP_OKAY;
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
   const int nsubedges = graph_get_nEdges(subgraph);

   assert(!tsepafull->solcands_sepaedges);
   assert(tsepafull->nsolcands == 0);
   assert(nsubedges >= 2);

   SCIP_CALL( reduce_termcompChangeSubgraphToBottleneck(scip, orggraph, termcomp, success) );

   if( !(*success) )
   {
      SCIPdebugMessage("problem while building bottleneck distances, aborting \n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(tsepafull->subgraph_bcosts), nsubedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(tsepafull->subgraph_orgcosts), nsubedges) );

   BMScopyMemoryArray(tsepafull->subgraph_bcosts, subgraph->cost, nsubedges);

   /* NOTE: we need to reset the costs for SDs to be valid */
   reduce_termcompChangeSubgraphToOrgCosts(orggraph, termcomp);
   BMScopyMemoryArray(tsepafull->subgraph_orgcosts, subgraph->cost, nsubedges);

   SCIP_CALL( sepafullInitDistdata(scip, termcomp, tsepafull) );

   /* now do the actual work */
   SCIP_CALL( sepafullBuildSolcandsEdges(scip, termcomp, tsepafull) );

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
      assert(builder->nsepatterms >= 3);
      maxratio = COMPONENT_MAXNODESRATIO_5PLUSCANDS;
   }

   SCIPdebugMessage("FINAL CHECK: noderatio=%f, maxratio=%f\n", noderatio, maxratio);

   if( noderatio > maxratio )
   {
      SCIPdebugMessage("...component is finally NOT promising! \n");
      return FALSE;
   }

   SCIPdebugMessage("...component is finally promising! \n");

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
   SCIP_CALL( substpsolver_initBC(scip, subgraph_copy, &substp) );
   SCIP_CALL( substpsolver_initHistory(substp) );

   SCIP_CALL( substpsolver_setMute(substp) );
   //SCIP_CALL( substpsolver_setProbIsIndependent(substp) );
   SCIP_CALL( substpsolver_setProbFullPresolve(substp) );

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


/** returns TRUE if sub-solution can be ruled out */
static
SCIP_Bool subSolIsRedundant(
   SCIP*                 scip,               /**< SCIP data structure */
   const TERMCOMP*       termcomp,           /**< component */
   const int*            subsol,             /**< solution */
   const int*            sepapartition,      /**< partition */
   STP_Vectype(int)      solcand_sepaedges   /**< bottleneck edges */
   )
{
   const GRAPH* subgraph = termcomp->subgraph;
   int* groups_size;
   int* groups_nsolbedges;
   const int nsepaterms = termcomp->builder->nsepatterms;
   const int ngroups = getPartitionNgroups(nsepaterms, sepapartition);
   SCIP_Bool isRedundant = FALSE;

   assert(nsepaterms >= 2);
   assert(subgraph);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &groups_size, ngroups) );
   SCIP_CALL_ABORT( SCIPallocClearBufferArray(scip, &groups_nsolbedges, ngroups) );

   getPartitionGroupsSizes(nsepaterms, ngroups, sepapartition, groups_size);

   for( int i = 0; i < StpVecGetSize(solcand_sepaedges); i++ )
   {
      const int subedge = solcand_sepaedges[i];
      const int tail = subgraph->tail[subedge];
      const int group = sepapartition[tail];

      assert(graph_edge_isInRange(termcomp->subgraph, subedge));
      assert(tail < nsepaterms);
      assert(subgraph->head[subedge] < nsepaterms);
      assert(sepapartition[tail] == sepapartition[subgraph->head[subedge]]);
      assert(group < ngroups);
      assert(groups_size[group] >= 2);

      if( subsol[subedge] == CONNECT || subsol[flipedge(subedge)] == CONNECT )
      {
         groups_nsolbedges[group]++;
      }
   }

   for( int i = 0; i < ngroups; i++ )
   {
      const int ngroupnodes = groups_size[i];
      const int ngroupsoledges = groups_nsolbedges[i];

      assert(ngroupsoledges < ngroupnodes && ngroupnodes <= nsepaterms);
      assert(0 <= ngroupsoledges);

      /* is there at least one unconnected group node? */
      if( ngroupsoledges < ngroupnodes - 1 )
      {
         SCIPdebugMessage("ruled out: ngroupsoledges=%d ngroupnodes=%d  \n", ngroupsoledges, ngroupnodes);
         isRedundant = TRUE;
         break;
      }
   }

   SCIPfreeBufferArray(scip, &groups_nsolbedges);
   SCIPfreeBufferArray(scip, &groups_size);

   return isRedundant;
}



/** sets weight for artificial bottleneck edges between separator terminals */
static
void setSubBottleneckEdges(
   int                   cand,               /**< candidate number */
   TERMCOMP*             termcomp,           /**< component */
   TSEPAFULL*            tsepafull           /**< full separation data  */
   )
{
   GRAPH* subgraph = termcomp->subgraph;
   const SCIP_Real* const subgraph_bcosts = tsepafull->subgraph_bcosts;
   STP_Vectype(int) solcand_sepaedges = tsepafull->solcands_sepaedges[cand];

   assert(subgraph_bcosts);

   for( int i = 0; i < StpVecGetSize(solcand_sepaedges); i++ )
   {
      const int subedge = solcand_sepaedges[i];
      assert(graph_edge_isInRange(subgraph, subedge));

      subgraph->cost[subedge] = subgraph_bcosts[subedge];
      subgraph->cost[flipedge(subedge)] = subgraph_bcosts[flipedge(subedge)];

#ifdef SCIP_DEBUG
      SCIPdebugMessage("adding artifical subedge (for solcand %d): ", cand);
      graph_edge_printInfo(subgraph, subedge);
#endif
   }
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
   STP_Vectype(int) solcand_sepaedges = tsepafull->solcands_sepaedges[cand];
   const int* const solcand_sepapart = tsepafull->solcands_sepapart[cand];
   const SCIP_Real* const subgraph_orgcosts = tsepafull->subgraph_orgcosts;

   assert(0 <= cand && cand < tsepafull->nsolcands);
   assert(StpVecGetSize(tsepafull->solcands_sepaedges) == StpVecGetSize(tsepafull->solcands_sepapart));
   assert(subgraph_orgcosts && solcand_sepapart);

   BMScopyMemoryArray(subgraph->cost, subgraph_orgcosts, subgraph->edges);
   setSubBottleneckEdges(cand, termcomp, tsepafull);

   SCIP_CALL( SCIPallocMemoryArray(scip, &subsol, subgraph->edges) );
   SCIP_CALL( solveSub(scip, subgraph, subsol) );
   assert(solstp_isValid(scip, subgraph, subsol));

   if( subSolIsRedundant(scip, termcomp, subsol, solcand_sepapart, solcand_sepaedges) )
   {
      SCIPdebugMessage("DISCARDING sub-solution \n");
      SCIPfreeMemoryArray(scip, &subsol);
   }
   else
   {
      SCIPdebugMessage("ADDING sub-solution \n");
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
   const SCIP_Real* const subgraph_bcosts = tsepafull->subgraph_bcosts;
   const SCIP_Real* const subgraph_orgcosts = tsepafull->subgraph_orgcosts;
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
            const int* const subToOrg = termcomp->nodemap_subToOrg;
            SCIPdebugMessage("solution %d (original) edge: ", i);

            if( edgemap_subToOrg[e] != -1  )
               graph_edge_printInfo(orggraph, edgemap_subToOrg[e]);
            else
               printf("artificial %d->%d \n", subToOrg[subgraph->tail[e]], subToOrg[subgraph->head[e]]);
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
      const int subedge = e * 2;
      const int orgedge = edgemap_subToOrg[subedge];

      if( orgedge == -1 )
         continue;

      assert(graph_edge_isInRange(orggraph, orgedge));

      /* we want to make sure that we don't modify sepa->sepa edges that got a smaller
       * bottleneck distance */
      if( !EQ(subgraph_bcosts[subedge], subgraph_orgcosts[subedge]) )
      {
         assert(LT(subgraph_bcosts[subedge], subgraph_orgcosts[subedge]));
         assert(subgraph->tail[subedge] < termcomp->builder->nsepatterms);
         assert(subgraph->head[subedge] < termcomp->builder->nsepatterms);

         continue;
      }

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
      SCIPdebugMessage("---Checking solution candidate %d: \n", i);

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
   assert(builder->nsepatterms <= SEPARATOR_MAXSIZE);

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

   SCIPdebugMessage("noderatio=%f, maxratio=%f\n", noderatio, maxratio);

   if( noderatio > maxratio )
   {
      SCIPdebugMessage("...component is NOT promising! \n");
      return FALSE;
   }

   SCIPdebugMessage("...component is promising! \n");

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
   int nnodes_real;
   SCIP_Real termratio;

   graph_get_nVET(g, &nnodes_real, NULL, NULL);
   termratio = (SCIP_Real) g->terms / (SCIP_Real) nnodes_real;



   if( maxncompchecks < mincheckbound )
   {
      SCIPdebugMessage("update nChecks %d->%d \n", maxncompchecks, mincheckbound);
      maxncompchecks = mincheckbound;
   }

   if( termratio >= SEPARATOR_MINEFFTERMRATIO )
   {
      maxncompchecks *= 1.5;
   }

#ifdef SCIP_DEBUG
   graph_printInfoReduced(g);
   printf("max. number of components to be checked: %d \n", maxncompchecks);
#endif

   /* NOTE: we want to allow a few more terminal separators to be able to choose small ones */
   SCIP_CALL( mincut_termsepasInit(scip, g, (int) (2.0 * maxncompchecks), maxsepasize, termsepas) );
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

  // int todo;
   /* first, we want to get rid of medium sized connected components */
   SCIP_CALL( reduce_bidecompositionExact(scip, g, redbase, solnode, nelims) );

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
