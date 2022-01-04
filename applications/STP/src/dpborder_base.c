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

/**@file   dpborder_base.c
 * @brief  Dynamic programming solver for Steiner tree (sub-) problems with small border
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dpborder.h"
#include "dpborderinterns.h"
#include "stpvector.h"
#include "solstp.h"


#define DPB_MINTERMRATIO 0.4
#define DPB_MAXNNODES    1000
#define DPB_MINNNODES    100
#define DPB_MAXAVGDEG    4.5



/*
 * Local methods
 */


/** non-promising? quick check */
static
SCIP_Bool dpborderIsNonPromising(
   const GRAPH*          graph               /**< original graph */
)
{
   if( graph->knots > DPB_MAXNNODES )
      return TRUE;

   if( graph->knots < DPB_MINNNODES )
      return TRUE;

   if( ((SCIP_Real) graph->terms / (SCIP_Real) graph->knots) < DPB_MINTERMRATIO )
      return TRUE;

   if( ((SCIP_Real) graph->edges / (SCIP_Real) graph->knots) > DPB_MAXAVGDEG )
      return TRUE;

   return FALSE;
}


/** initializes helper */
static
SCIP_RETCODE dpborderInitHelper(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph of sub-problem */
   DPBORDER*             dpborder           /**< border */
   )
{
   const int nnodes = graph_get_nNodes(graph);

   assert(nnodes == dpborder->nnodes);
   assert(!dpborder->nodes_isBorder);
   assert(!dpborder->nodes_outdeg);

   SCIP_CALL( dpborder_dpbsequenceInit(scip, graph, &(dpborder->dpbsequence)) );

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &(dpborder->nodes_isBorder), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(dpborder->nodes_outdeg), nnodes) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(dpborder->bordercharmap), BPBORDER_MAXBORDERSIZE + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(dpborder->borderchardists), BPBORDER_MAXBORDERSIZE) );

   BMScopyMemoryArray(dpborder->nodes_outdeg, graph->grad, nnodes);

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


/** initializes */
SCIP_RETCODE dpborder_dpbsequenceInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< original graph */
   DPBSEQUENCE**         dpbsequence         /**< to initialize */
)
{
   DPBSEQUENCE* seq;
   const int nnodes = graph_get_nNodes(graph);

   SCIP_CALL( SCIPallocMemory(scip, dpbsequence) );
   seq = *dpbsequence;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(seq->nodessquence), nnodes) );
   seq->maxnpartitions = 0;
   seq->maxbordersize = nnodes + 1;
   seq->nnodes = nnodes;

   return SCIP_OKAY;
}


/** copies */
void dpborder_dpbsequenceCopy(
   const DPBSEQUENCE*    dpbsequence_source, /**< to copy */
   DPBSEQUENCE*          dpbsequence_target  /**< to copy to */
)
{
   const int nnodes = dpbsequence_source->nnodes;
   assert(dpbsequence_source->nnodes == dpbsequence_target->nnodes);

   BMScopyMemoryArray(dpbsequence_target->nodessquence, dpbsequence_source->nodessquence, nnodes);
   dpbsequence_target->maxbordersize = dpbsequence_source->maxbordersize;
   dpbsequence_target->maxnpartitions = dpbsequence_source->maxnpartitions;
}


/** frees */
void dpborder_dpbsequenceFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPBSEQUENCE**         dpbsequence         /**< to free */
)
{
   DPBSEQUENCE* seq =  *dpbsequence;

   assert(seq);

   SCIPfreeMemoryArray(scip, &(seq->nodessquence));
   SCIPfreeMemory(scip, dpbsequence);
}


/** initializes */
SCIP_RETCODE dpborder_dpblevelInit(
   SCIP*                 scip,               /**< SCIP data structure */
   DPBLEVEL**            dpblevel            /**< to initialize */
)
{
   DPBLEVEL* level;

   SCIP_CALL( SCIPallocMemory(scip, dpblevel) );
   level = *dpblevel;

   level->bordernodesMapToOrg = NULL;
   level->nbordernodes = 0;
   level->extnode = -1;
   level->globalstartidx = -1;
   level->exnodeIsTerm = FALSE;

   return SCIP_OKAY;
}


/** frees */
void dpborder_dpblevelFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPBLEVEL**            dpblevel            /**< to be freed */
)
{
   DPBLEVEL* level = *dpblevel;

   assert(level);

   StpVecFree(scip, level->bordernodesMapToOrg);
   SCIPfreeMemory(scip, dpblevel);
}


/** checks whether DP border has potential
 * NOTE: needs to be called before dpborder_solve! */
SCIP_RETCODE dpborder_probePotential(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph of sub-problem */
   DPBORDER*             dpborder,           /**< border */
   SCIP_Bool*            hasPotential        /**< was problem solved to optimality? */
)
{
   DPBSEQUENCE* dpbsequence;
   assert(scip && graph && dpborder && hasPotential);
   assert(!dpborder->dpbsequence);

   *hasPotential = FALSE;

   /* make rough check */
   if( dpborderIsNonPromising(graph) )
   {
      return SCIP_OKAY;
   }

   SCIP_CALL( graph_init_csrWithEdgeId(scip, graph) );

   SCIP_CALL( dpborderInitHelper(scip, graph, dpborder) );
   SCIP_CALL( dpborder_coreComputeOrderingSimple(scip, graph, dpborder) );

   dpbsequence = dpborder->dpbsequence;
   *hasPotential = TRUE;

   if( dpbsequence->maxnpartitions >= BPBORDER_MAXNPARTITIONS )
   {
      *hasPotential = FALSE;
   }

   if( dpbsequence->maxbordersize >= BPBORDER_MAXBORDERSIZE )
   {
      *hasPotential = FALSE;
   }

   graph_free_csr(scip, graph);

   return SCIP_OKAY;
}


/** solves problem given by graph */
SCIP_RETCODE dpborder_solve(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph of sub-problem */
   DPBORDER*             dpborder,           /**< border */
   int*                  solution,           /**< optimal solution (out) */
   SCIP_Bool*            wasSolved           /**< was problem solved to optimality? */
)
{
   uint64_t maxnpartitions;
   assert(scip && graph && solution);
   assert(dpborder->dpbsequence);

   SCIP_CALL( graph_init_csr(scip, graph) );

   SCIP_CALL( dpborder_coreUpdateOrdering(scip, graph, dpborder) );
   maxnpartitions = pow(2, ceil(log(dpborder->dpbsequence->maxnpartitions)/log(2)));
   assert(maxnpartitions >= dpborder->dpbsequence->maxnpartitions);

   SCIP_CALL( hashmap_create(DPBORDER_GROWTH_FACTOR * maxnpartitions, NULL, &dpborder->hashmap) );
   SCIP_CALL( dpborder_coreSolve(scip, graph, dpborder, wasSolved) );

   if( *wasSolved )
   {
      STP_Bool* connected;
      SCIP_CALL( SCIPallocBufferArray(scip, &connected, graph->knots) );

      dpborder_markSolNodes(dpborder, connected);
      SCIP_CALL( solstp_pruneFromNodes(scip, graph, solution, connected) );
      assert(solstp_isValid(scip, graph, solution));

      SCIPfreeBufferArray(scip, &connected);
   }

#ifndef NDEBUG
   if( *wasSolved )
   {
      const SCIP_Real solval =  solstp_getObj(graph, solution, 0.0);
      assert(EQ(solval, dpborder->global_obj));
   }
#endif

   graph_free_csr(scip, graph);

   return SCIP_OKAY;
}


/** initializes */
SCIP_RETCODE dpborder_init(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< original graph */
   DPBORDER**            dpborder            /**< to initialize */
)
{
   DPBORDER* dpb;
   SCIP_CALL( SCIPallocMemory(scip, dpborder) );
   dpb = *dpborder;

   assert(graph);

   dpb->global_partsUseExt = NULL;
   dpb->bordercharmap = NULL;
   dpb->borderchardists = NULL;
   dpb->dpbsequence = NULL;
   dpb->borderlevels = NULL;
   dpb->bordernodes = NULL;
   dpb->prevbordernodes = NULL;
   dpb->global_partitions = NULL;
   dpb->global_partstarts = NULL;
   dpb->global_predparts = NULL;
   dpb->global_partcosts = NULL;
   dpb->nodes_isBorder = NULL;
   dpb->nodes_outdeg = NULL;
   dpb->global_obj = FARAWAY;
   dpb->global_npartitions = 0;
   dpb->global_partcap = 0;
   dpb->nnodes = graph->knots;
   dpb->nterms = graph->terms;
   dpb->ntermsvisited = 0;
   dpb->global_optposition = -1;
   dpb->hashmap.elements = NULL;

   return SCIP_OKAY;
}


/** frees */
void dpborder_free(
   SCIP*                 scip,               /**< SCIP data structure */
   DPBORDER**            dpborder            /**< to be freed */
)
{
   DPBORDER* dpb = *dpborder;

   StpVecFree(scip, dpb->global_partcosts);
   StpVecFree(scip, dpb->global_predparts);
   StpVecFree(scip, dpb->global_partstarts);
   StpVecFree(scip, dpb->bordernodes);
   StpVecFree(scip, dpb->prevbordernodes);
   StpVecFree(scip, dpb->global_partsUseExt);

   SCIPfreeMemoryArrayNull(scip, &(dpb->bordercharmap));
   SCIPfreeMemoryArrayNull(scip, &(dpb->borderchardists));
   SCIPfreeMemoryArrayNull(scip, &(dpb->global_partitions));
   SCIPfreeMemoryArrayNull(scip, &(dpb->nodes_isBorder));
   SCIPfreeMemoryArrayNull(scip, &(dpb->nodes_outdeg));

   // todo make clean
   if( dpb->hashmap.elements )
   {
      hashmap_destroy(&dpb->hashmap);
   }

   if( dpb->dpbsequence )
   {
      dpborder_dpbsequenceFree(scip, &(dpb->dpbsequence));
   }

   if( dpb->borderlevels )
   {
      for( int i = 0; i < StpVecGetSize(dpb->borderlevels); i++ )
      {
         dpborder_dpblevelFree(scip, &(dpb->borderlevels[i]));
         assert(!dpb->borderlevels[i]);
      }

      StpVecFree(scip, dpb->borderlevels);
   }

   SCIPfreeMemory(scip, dpborder);
}
