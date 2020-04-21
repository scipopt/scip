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

/**@file   heur_tm.c
 * @brief  Shortest paths based primal heuristics for Steiner problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 *
 * This file implements several shortest paths based primal heuristics for Steiner problems, see
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 *
 * A list of all interface methods can be found in heur_tm.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>
#include "heur_tm.h"
#include "probdata_stp.h"
#include "portab.h"
#include "solstp.h"
#include "scip/misc.h"
#include "shortestpath.h"
#include <math.h>

#define HEUR_NAME             "TM"
#define HEUR_DESC             "shortest path based primal heuristics for Steiner trees"
#define HEUR_DISPCHAR         '+'
#define HEUR_PRIORITY         10000000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_EVALRUNS 20                  /**< number of runs */
#define DEFAULT_INITRUNS 100                 /**< number of initial runs */
#define DEFAULT_LEAFRUNS 25                  /**< number of runs at leafs */
#define DEFAULT_ROOTRUNS 50                  /**< number of runs at the root */
#define DEFAULT_DURINGLPFREQ 5               /**< frequency during LP solving */
#define DEFAULT_TYPE  0                      /**< heuristic to execute */
#define DEFAULT_PCMODE  4                    /**< solving mode for PC/MW */

#define DEFAULT_RANDSEED 5                   /**< seed for pseudo-random functions */

#define AUTO        0
#define TM_SP       1
#define TM_VORONOI  2
#define TM_DIJKSTRA 3

#define TM_USE_CSR
#define TM_USE_CSR_PCMW

#ifdef WITH_UG
int getUgRank(void);
#endif

/*
 * Data structures
 */


typedef
struct TM_base_data
{
   int*                  heap_position;      /**< heap position array or  */
   DENTRY*               heap_entries;       /**< entries array or NULL */
   CSR*                  csr;                /**< CSR with possible biased costs.
                                                  NOTE: shares memory with csr_orgcosts! */
   CSR*                  csr_orgcosts;       /**< CSR with original costs.
                                                  NOTE: shares memory with csr! */
   DHEAP*                dheap;              /**< Dijkstra heap */
   const SCIP_Real*      cost;               /**< arc costs */
   const SCIP_Real*      costrev;            /**< reversed arc costs */
   SCIP_Real*            nodes_dist;         /**< distance values for each node */
   int*                  nodes_pred;         /**< predecessor for each node */
   int*                  startnodes;         /**< array containing start vertices (NULL to not provide any) */
   int*                  result;             /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   int*                  best_result;        /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   STP_Bool*             connected;          /**< array indicating whether a node is part of solution (TRUE/FALSE) */
   SCIP_Real             best_obj;           /**< objective */
   int                   nruns;              /**< number of runs */
} TMBASE;


/** All shortest paths data
 *  NOTE: DEPRECATED */
typedef
struct TM_all_shortestpath
{
   int**                 pathedge;           /**< node predecessors for each terminal */
   SCIP_Real**           pathdist;           /**< node distances for each terminal */
} TMALLSP;


/** Voronoi TM heuristic data
 *  NOTE: DEPRECATED */
typedef
struct TM_voronoi
{
   SCIP_PQUEUE* pqueue;
   SCIP_Real** node_dist;
   GNODE** gnodearr;
   int* nodenterms;
   int** node_base;
   int** node_edge;
} TMVNOI;


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          nlpiterations;      /**< number of total LP iterations*/
   SCIP_Longint          ncalls;             /**< number of total calls (of TM) */
   SCIP_Longint          nexecs;             /**< number of total executions (of TM) */
   SCIP_Real             hopfactor;          /**< edge multiplication factor for hop constrained problems */
   int                   stp_type;           /**< problem type */
   int                   evalruns;           /**< number of runs */
   int                   initruns;           /**< number of initial runs */
   int                   leafruns;           /**< number of runs at leafs */
   int                   rootruns;           /**< number of runs at the root */
   int                   duringlpfreq;       /**< frequency during LP solving */
   int                   type;               /**< Heuristic type: 0 automatic, 1 TM_SP, 2 TM_VORONOI, 3 TM_DIJKSTRA */
   int                   beststartnode;      /**< start node of the so far best found solution */
   unsigned int          randseed;           /**< seed value for random number generator */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   unsigned int          timing;             /**< timing for timing mask */
   enum PCMW_TmMode      pcmw_mode;          /**< solving mode for PCMW */
};

/*
 * Static methods
 */

#if defined(TM_USE_CSR_PCMW) && !defined(NDEBUG)
/** debug method; checks whether CSR solution objective is actually computed correctly */
static
SCIP_Bool pcmwUpdateBestSol_csrInSync(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const TMBASE*         tmbase              /**< data */
)
{
   const int *const result_csr = tmbase->result;
   const STP_Bool *const connected = tmbase->connected;
   const SCIP_Real obj = solstp_pcGetObjCsr(graph, tmbase->csr_orgcosts, result_csr, connected);
   const int nedges = graph_get_nEdges(graph);
   const int nnodes = graph_get_nNodes(graph);
   int *result_dummy;
   STP_Bool *connected_dummy;
   SCIP_Real obj_dummy;
   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip, &result_dummy, nedges));
   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip, &connected_dummy, nnodes));
   BMScopyMemoryArray(connected_dummy, connected, nnodes);

   solstp_convertCsrToGraph(scip, graph, tmbase->csr_orgcosts, result_csr, connected_dummy, result_dummy);

   assert(!stpsol_pruningIsPossible(graph, result_dummy, connected_dummy));

   obj_dummy = solstp_getObj(graph, result_dummy, 0.0, nedges);
   obj_dummy += graph_pc_getNonLeafTermOffset(scip, graph);

   SCIPfreeMemoryArray(scip, &connected_dummy);
   SCIPfreeMemoryArray(scip, &result_dummy);

   return SCIPisEQ(scip, obj_dummy, obj);
}
#endif

/** information method for a parameter change of random seed */
static
SCIP_DECL_PARAMCHGD(paramChgdRandomseed)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int newrandseed;

   newrandseed = SCIPparamGetInt(param);

   heurdata = (SCIP_HEURDATA*)SCIPparamGetData(param);
   assert(heurdata != NULL);

   heurdata->randseed = (unsigned int)newrandseed;

   return SCIP_OKAY;
}


/*
 *  local functions
 */

/** initializes */
static
SCIP_RETCODE tmBaseInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   TMBASE*               tmbase              /**< data */
)
{
   SCIP_Real* costorg_csr;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);

   assert(tmbase && scip);
   assert(tmbase->nodes_dist == NULL);
   assert(tmbase->nodes_pred == NULL);
   assert(tmbase->startnodes == NULL);
   assert(tmbase->connected == NULL);
   assert(tmbase->result == NULL);
   assert(tmbase->cost);

   if( tmbase->nruns < 1 )
      tmbase->nruns = 1;

   SCIP_CALL( SCIPallocBufferArray(scip, &(tmbase->startnodes), MIN(tmbase->nruns, nnodes)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmbase->result), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmbase->nodes_dist), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmbase->nodes_pred), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmbase->connected), nnodes) );

#ifdef TM_USE_CSR
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmbase->heap_position), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmbase->heap_entries), nnodes + 2) );
   SCIP_CALL( graph_heap_create(scip, nnodes, tmbase->heap_position, tmbase->heap_entries, &(tmbase->dheap)) );
   SCIP_CALL( graph_csr_allocWithEdgeId(scip, nnodes, nedges, &(tmbase->csr)) );
   graph_csr_build(graph, tmbase->cost, tmbase->csr);

   /* build the original cost CSR.
    * NOTE: it consists mostly of alias, except for the edge costs */
   SCIP_CALL( SCIPallocBufferArray(scip, &(costorg_csr), nedges) );
   SCIP_CALL( SCIPallocMemory(scip, &(tmbase->csr_orgcosts)) );
   *(tmbase->csr_orgcosts) = *(tmbase->csr);
   tmbase->csr_orgcosts->cost = costorg_csr;

   if( graph_pc_isPc(graph) )
   {
      graph_pc_getOrgCostsCsr(scip, graph, tmbase->csr_orgcosts);
   }
   else
   {
      graph_csr_buildCosts(graph, tmbase->csr_orgcosts, graph->cost, costorg_csr);
   }
#else
   tmbase->heap_position = NULL;
   tmbase->heap_entries = NULL;
   tmbase->dheap = NULL;
   tmbase->csr = NULL;
   tmbase->csr_orgcosts = NULL;
#endif

   return SCIP_OKAY;
}


/** frees */
static
void tmBaseFree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   TMBASE*               tmbase              /**< data */
)
{
#ifdef TM_USE_CSR
   SCIPfreeBufferArray(scip, &(tmbase->csr_orgcosts->cost));
   SCIPfreeMemory(scip, &(tmbase->csr_orgcosts));
   graph_csr_free(scip, &(tmbase->csr));
   graph_heap_free(scip, FALSE, FALSE, &(tmbase->dheap));
   SCIPfreeBufferArray(scip, &(tmbase->heap_entries));
   SCIPfreeBufferArray(scip, &(tmbase->heap_position));
#endif

   assert(tmbase->dheap == NULL);
   assert(tmbase->csr == NULL);
   assert(tmbase->heap_position == NULL);
   assert(tmbase->heap_entries == NULL);
   SCIPfreeBufferArray(scip, &(tmbase->connected));
   SCIPfreeBufferArray(scip, &(tmbase->nodes_pred));
   SCIPfreeBufferArray(scip, &(tmbase->nodes_dist));
   SCIPfreeBufferArray(scip, &(tmbase->result));
   SCIPfreeBufferArray(scip, &(tmbase->startnodes));
}


/** initializes */
static
SCIP_RETCODE tmAllspInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   TMALLSP*              tmallsp             /**< data */
)
{
   const int nnodes = graph_get_nNodes(graph);

   SCIP_CALL(SCIPallocBufferArray(scip, &(tmallsp->pathdist), nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &(tmallsp->pathedge), nnodes));
   BMSclearMemoryArray((tmallsp->pathdist), nnodes);
   BMSclearMemoryArray((tmallsp->pathedge), nnodes);

   for( int k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);

      if( Is_term(graph->term[k]) )
      {
         SCIP_CALL(SCIPallocBufferArray(scip, &((tmallsp->pathdist[k])), nnodes)); /*lint !e866*/
         SCIP_CALL(SCIPallocBufferArray(scip, &((tmallsp->pathedge[k])), nnodes)); /*lint !e866*/
      }
      else
      {
         tmallsp->pathdist[k] = NULL;
         tmallsp->pathedge[k] = NULL;
      }
   }

   return SCIP_OKAY;
}


/** frees */
static
void tmAllspFree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   TMALLSP*              tmallsp             /**< data */
)
{
   const int nnodes = graph_get_nNodes(graph);

   assert(tmallsp->pathedge);
   assert(tmallsp->pathdist);

   for( int k = nnodes - 1; k >= 0; k-- )
   {
      SCIPfreeBufferArrayNull(scip, &(tmallsp->pathedge[k]));
      SCIPfreeBufferArrayNull(scip, &(tmallsp->pathdist[k]));
   }

   SCIPfreeBufferArray(scip, &(tmallsp->pathedge));
   SCIPfreeBufferArray(scip, &(tmallsp->pathdist));
}


/** initializes */
static
SCIP_RETCODE tmVnoiInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   TMVNOI*               tmvnoi              /**< data */
)
{
   const int nnodes = graph_get_nNodes(graph);
   const int nterms = graph->terms;

   assert(nnodes > 0 && nterms > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &(tmvnoi->nodenterms), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmvnoi->gnodearr), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmvnoi->node_base), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmvnoi->node_dist), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(tmvnoi->node_edge), nnodes) );

   for( int k = 0; k < nnodes; k++ )
   {
      SCIP_CALL( SCIPallocBuffer(scip, &(tmvnoi->gnodearr[k])) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &(tmvnoi->node_base[k]), nterms) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &(tmvnoi->node_dist[k]), nterms) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &(tmvnoi->node_edge[k]), nterms) ); /*lint !e866*/
   }

   SCIP_CALL( SCIPpqueueCreate( &(tmvnoi->pqueue), nnodes, 2.0, GNODECmpByDist) );

   return SCIP_OKAY;
}


/** frees */
static
void tmVnoiFree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   TMVNOI*               tmvnoi              /**< data */
)
{
   const int nnodes = graph_get_nNodes(graph);

   SCIPpqueueFree(&(tmvnoi->pqueue));

   assert(tmvnoi->node_edge != NULL);
   assert(tmvnoi->node_dist != NULL);
   assert(tmvnoi->node_base != NULL);
   assert(tmvnoi->gnodearr != NULL);

   for( int k = nnodes - 1; k >= 0; k-- )
   {
      SCIPfreeBufferArray(scip, &(tmvnoi->node_edge[k]));
      SCIPfreeBufferArray(scip, &(tmvnoi->node_dist[k]));
      SCIPfreeBufferArray(scip, &(tmvnoi->node_base[k]));
      SCIPfreeBuffer(scip, &(tmvnoi->gnodearr[k]));
   }

   SCIPfreeBufferArray(scip, &(tmvnoi->node_edge));
   SCIPfreeBufferArray(scip, &(tmvnoi->node_dist));
   SCIPfreeBufferArray(scip, &(tmvnoi->node_base));
   SCIPfreeBufferArray(scip, &(tmvnoi->gnodearr));
   SCIPfreeBufferArray(scip, &(tmvnoi->nodenterms));
}


/** returns TM heur data */
static inline
SCIP_HEURDATA* getTMheurData(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_HEURDATA* heurdata;

   assert(scip);

   heurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));
   assert(heurdata);

   return heurdata;
}


/** returns mode */
static
int getTmMode(
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   const GRAPH*          graph               /**< graph data structure */
)
{
   /* get user parameter */
   int mode = heurdata->type;
   assert(mode == AUTO || mode == TM_SP || mode == TM_VORONOI || mode == TM_DIJKSTRA);

   if( graph_pc_isPcMw(graph) )
   {
      mode = TM_DIJKSTRA;
   }
   else if( graph->stp_type == STP_DHCSTP )
   {
      mode = TM_SP;
   }
   else if( graph->stp_type == STP_DCSTP )
   {
      mode = TM_SP;
   }
   else if( graph->stp_type == STP_SAP )
   {
      mode = TM_SP;
   }
   else
   {
      if( mode == AUTO )
         mode = TM_DIJKSTRA;
   }

   return mode;
}


/** updates best PC/MW solution if possible */
static inline
void pcmwUpdateBestSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   TMBASE*               tmbase,             /**< data */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
)
{
#ifdef TM_USE_CSR_PCMW
   const int* const result_csr = tmbase->result;
   STP_Bool* const connected = tmbase->connected;
   /* compute objective value w.r.t. original costs! */
   const SCIP_Real obj = solstp_pcGetObjCsr(graph, tmbase->csr_orgcosts, result_csr, connected);

   assert(pcmwUpdateBestSol_csrInSync(scip, graph, tmbase));

   if( LT(obj, tmbase->best_obj) )
   {
      SCIPdebugMessage("\n improved obj=%f ", obj);

      tmbase->best_obj = obj;
      solstp_convertCsrToGraph(scip, graph, tmbase->csr_orgcosts, result_csr, connected, tmbase->best_result);

      (*success) = TRUE;
   }

#else
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real obj = solstp_getObj(graph, tmbase->result, 0.0, nedges);

   if( LT(obj, tmbase->best_obj) )
   {
      SCIPdebugMessage("\n improved obj=%f ", obj);

      tmbase->best_obj = obj;
      BMScopyMemoryArray(tmbase->best_result, tmbase->result, nedges);
      (*success) = TRUE;
   }
#endif
}


/** submethod for SCIPStpHeurTMRun in PC or MW mode */
static
void pcmwAdaptStarts(
   SCIP_HEURDATA*        heurdata,           /**< heurdata */
   const GRAPH*          graph,              /**< graph data structure */
   int                   maxtmruns,          /**> number of TM runs */
   int                   bestincstart,       /**< best incumbent start vertex */
   int*                  terminalperm        /**< terminal permutation */
   )
{
   const int nnodes = graph_get_nNodes(graph);

   assert(maxtmruns <= graph->terms);

   if( maxtmruns > 0 && bestincstart >= 0 && bestincstart < nnodes && Is_pseudoTerm(graph->term[bestincstart])
        && SCIPrandomGetInt(heurdata->randnumgen, 0, 2) == 1 )
   {
      int r;
      for( r = 0; r < maxtmruns; r++ )
         if( terminalperm[r] == bestincstart )
            break;

      if( r == maxtmruns )
      {
         assert(maxtmruns > 1);
         terminalperm[1] = bestincstart;
      }
   }
}


/** submethod for runPCMW */
static
SCIP_RETCODE pcmwGetStartNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      nodepriority,       /**< vertex priorities for vertices to be starting points (NULL for no priorities) */
   int                   maxtmruns,          /**> number of TM runs */
   int                   bestincstart,       /**< best incumbent start vertex */
   GRAPH*                graph,              /**< graph data structure */
   int*                  terminalperm        /**< terminal permutation */
   )
{
   SCIP_HEURDATA* heurdata = getTMheurData(scip);
   SCIP_Real *terminalprio; /**< terminal priority */
   const int root = graph->source;
   const int nnodes = graph->knots;
   const int nterms = graph->terms;
   int termcount = 0;

   assert(scip && heurdata && terminalperm);
   assert(graph->extended);

   SCIP_CALL( SCIPallocBufferArray(scip, &terminalprio, nterms) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_pseudoTerm(graph->term[k]) && graph->grad[k] != 0 )
      {
         assert(!nodepriority || nodepriority[k] >= 0.0);
         assert(graph->term2edge[k] < 0 || SCIPisGT(scip, graph->prize[k], 0.0));
         terminalperm[termcount] = k;

         if( nodepriority == NULL )
            terminalprio[termcount++] = -SCIPrandomGetReal(heurdata->randnumgen, 0.0, graph->prize[k]);
         else
            terminalprio[termcount++] = -SCIPrandomGetReal(heurdata->randnumgen, nodepriority[k] / 2.0, nodepriority[k]);
      }
   }

   if( graph_pc_isRootedPcMw(graph) )
   {
      SCIP_Real minprio = 0.0;

      for( int k = 0; k < termcount; k++ )
         minprio = MIN(terminalprio[k], minprio);

      // todo why only marking for this case???
      graph_pc_markOrgGraph(graph);

      graph->mark[graph->source] = TRUE;

      for( int k = 0; k < nnodes; k++ )
      {
         if( graph_pc_knotIsFixedTerm(graph, k) )
         {
            assert(graph->mark[k]);

            if( k == root )
               terminalprio[termcount] = -FARAWAY;

            terminalperm[termcount] = k;
            terminalprio[termcount++] = -SCIPrandomGetReal(heurdata->randnumgen, 0.0,
                  fabs(minprio) * (1.0 + (SCIP_Real) graph->grad[k] / (SCIP_Real) nnodes));
         }
         /* isolated vertex? */
         else if( Is_term(graph->term[k]) && graph->grad[k] == 1 )
            terminalprio[termcount] = FARAWAY;
      }

      assert(nterms == termcount);
      SCIPsortRealInt(terminalprio, terminalperm, nterms);
   }
   else
   {
      SCIPsortRealInt(terminalprio, terminalperm, nterms - 1);
   }

   SCIPfreeBufferArray(scip, &terminalprio);

   pcmwAdaptStarts(heurdata, graph, maxtmruns, bestincstart, terminalperm);

   return SCIP_OKAY;
}


/** submethod for SCIPStpHeurTMRun in PC or MW mode */
static
void pcmwSetEdgeCosts(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< arc costs */
   const SCIP_Real*      costorg,            /**< arc costs */
   const SCIP_Real*      costfullbiased,     /**< arc costs */
   enum PCMW_TmMode      pcmwmode,           /**< mode */
   TMBASE*               tmbase              /**< data */
)
{
   if( pcmwmode == pcmode_biasfull )
   {
      tmbase->cost = costfullbiased;
   }
   else if( pcmwmode == pcmode_simple )
   {
      assert(graph_pc_isPc(graph));
      tmbase->cost = costorg;
   }
   else
   {
      assert(pcmwmode == pcmode_bias || pcmwmode == pcmode_fulltree);
      tmbase->cost = cost;
   }

   assert(tmbase->cost);

   graph_csr_chgCosts(graph, tmbase->cost, tmbase->csr);
}


/* compute starting vertices and number of runs */
static
SCIP_RETCODE computeStarts(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   const int*            orgstarts,          /**< original start vertices */
   SCIP_Bool             startsgiven,        /**< start vertices given? */
   SCIP_Real*            nodepriority,       /**< nodepriority */
   TMBASE*               tmbase,             /**< data */
   int*                  bestp               /**< pointer to best start vertex */
   )
{
   SCIP_HEURDATA* heurdata = getTMheurData(scip);
   int runs = tmbase->nruns;
   const int root = graph->source;
   const int nnodes = graph->knots;
   const int nterms = graph->terms;
   const SCIP_Real* cost = tmbase->cost;
   int* const start = tmbase->startnodes;
   int* const dijkedge = tmbase->nodes_pred;
   SCIP_Real* const dijkdist = tmbase->nodes_dist;

   /* compute number of iterations and starting points for SPH */
   if( graph_pc_isPcMw(graph) )
   {
      if( runs > (nterms - 1) && !graph_pc_isRootedPcMw(graph) )
         runs = nterms - 1;
      else if( runs > (nterms) && graph_pc_isRootedPcMw(graph) )
         runs = nterms;
   }
   else if( graph->stp_type == STP_NWPTSPG )
   {
      int count = 0;
      const int shift = SCIPrandomGetInt(heurdata->randnumgen, 0, nnodes);

      assert(!startsgiven);

      /* add terminals */
      for( int k = 0; k < nnodes && count < runs; k++ )
      {
         const int v = (k + shift) % nnodes;
         if( !Is_term(graph->term[v]) || graph_nw_knotIsLeaf(graph, v) )
            continue;
         start[count++] = v;
      }

      /* add non-terminals */
      for( int k = 0; k < nnodes && count < runs; k++ )
      {
         const int v = (k + shift) % nnodes;
         if( Is_term(graph->term[v]) )
            continue;
         start[count++] = v;
      }
      runs = count;
   }
   else if( startsgiven )
   {
      for( int k = 0; k < MIN(runs, nnodes); k++ )
         start[k] = orgstarts[k];
   }
   else if( runs < nnodes || graph->stp_type == STP_DHCSTP )
   {
      int r = 0;
      int* perm;

      if( *bestp < 0 )
      {
         *bestp = root;
      }
      else
      {
         const int randint = SCIPrandomGetInt(heurdata->randnumgen, 0, 2);
         if( randint == 0 )
            *bestp = -1;
         else if( randint == 1 )
            *bestp = root;
      }

      if( graph->stp_type == STP_DHCSTP )
         graph_path_execX(scip, graph, root, cost, dijkdist, dijkedge);

      SCIP_CALL( SCIPallocBufferArray(scip, &perm, nnodes) );

      /* are there no nodes are to be prioritized a priori? */
      if( nodepriority == NULL )
      {
         for( int k = 0; k < nnodes; k++ )
            perm[k] = k;
         SCIPrandomPermuteIntArray(heurdata->randnumgen, perm, 0, nnodes);

         /* use terminals (randomly permutated) as starting points for TM heuristic */
         for( int k = 0; k < nnodes; k++ )
         {
            if( r >= runs || r >= nterms )
               break;

            if( Is_term(graph->term[perm[k]]) )
               start[r++] = perm[k];
         }

         /* fill empty slots */
         for( int k = 0; k < nnodes && r < runs; k++ )
         {
            if( graph->stp_type == STP_DHCSTP )
            {
               assert(dijkdist != NULL);
               if( SCIPisGE(scip, dijkdist[perm[k]], BLOCKED) )
                  continue;
            }
            if( !Is_term(graph->term[perm[k]]) && graph->mark[k] )
               start[r++] = perm[k];
         }
      }
      else
      {
         SCIP_Real max = 0.0;
         const int bbound = runs - runs / 3;

         for( int k = 0; k < nnodes; k++ )
         {
            perm[k] = k;
            if( SCIPisLT(scip, max, nodepriority[k]) && Is_term(graph->term[k]) )
               max = nodepriority[k];
         }
         for( int k = 0; k < nnodes; k++ )
         {
            if( Is_term(graph->term[k]) )
               nodepriority[k] += SCIPrandomGetReal(heurdata->randnumgen, 0.0, max);
            else if( SCIPisLE(scip, 1.0, nodepriority[k]) )
               nodepriority[k] = nodepriority[k] * SCIPrandomGetReal(heurdata->randnumgen, 1.0, 2.0);
         }

         SCIPsortRealInt(nodepriority, perm, nnodes);

         for( int k = nnodes - 1; k >= 0; k-- )
         {
            if( r >= nterms || r >= bbound )
               break;

            if( Is_term(graph->term[perm[k]]) )
            {
               start[r++] = perm[k];
               perm[k] = -1;
            }
         }

         /* fill empty slots */
         for( int k = nnodes - 1; k >= 0 && r < runs; k-- )
         {
            if( perm[k] == -1 )
               continue;
            if( graph->stp_type == STP_DHCSTP )
            {
               assert(dijkdist != NULL);
               if( SCIPisGE(scip, dijkdist[perm[k]], BLOCKED) )
                  continue;
            }
            if( graph->mark[k] )
               start[r++] = perm[k];
         }
      }
      /* not all slots filled? */
      if( r < runs )
         runs = r;

      if( *bestp >= 0 )
      {
         /* check whether we have a already selected the best starting node */
         for( r = 0; r < runs; r++ )
            if( start[r] == *bestp )
               break;

         /* no, we still have to */
         if( r == runs && runs > 0 )
            start[(int) heurdata->nexecs % runs] = *bestp;
      }

      SCIPfreeBufferArray(scip, &perm);
   } /* runs < nnodes || STP_DHCSTP */
   else
   {
      runs = nnodes;
      for( int k = 0; k < nnodes; k++ )
         start[k] = k;
   } /* runs > nnodes */

   tmbase->nruns = runs;

   return SCIP_OKAY;
}

#ifdef TM_USE_CSR
/** CSR based shortest paths heuristic */
static inline
SCIP_RETCODE computeSteinerTreeCsr(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int                   startnode,          /**< start vertex*/
   TMBASE*               tmbase              /**< (in/out) */
)
{
   int* RESTRICT result = tmbase->result;
   SPATHS spaths = { .csr = tmbase->csr, .csr_orgcosts = tmbase->csr_orgcosts, .nodes_dist = tmbase->nodes_dist,
                     .nodes_pred = tmbase->nodes_pred, .dheap = tmbase->dheap, .nodes_isConnected = tmbase->connected };

   assert(g->stp_type != STP_DHCSTP);

   shortestpath_computeSteinerTree(g, startnode, &spaths);

   // SCIP_CALL( solstp_pruneFromTmHeur(scip, g, NULL, result, connected) );
   SCIP_CALL( solstp_pruneFromTmHeur_csr(scip, g, &spaths, result));

   return SCIP_OKAY;
}
#endif


/** Dijkstra based shortest paths heuristic */
static inline
SCIP_RETCODE computeSteinerTreeDijk(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   int                   start,              /**< start vertex */
   TMBASE*               tmbase              /**< (in/out) */
   )
{
   const SCIP_Real* cost = tmbase->cost;
   int* dijkedge = tmbase->nodes_pred;
   SCIP_Real* dijkdist = tmbase->nodes_dist;
   int* RESTRICT result = tmbase->result;

   assert(!graph_pc_isPcMw(g));

   graph_mark(g);

   graph_path_st(g, cost, dijkdist, dijkedge, start, tmbase->connected);

   /* cost will actually only be used for hop-constrained problem */
   SCIP_CALL( solstp_pruneFromTmHeur(scip, g, cost, result, tmbase->connected) );

   return SCIP_OKAY;
}

/** Dijkstra based shortest paths heuristic for PCSTP and MWCSP */
static
SCIP_RETCODE computeSteinerTreeDijkPcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   const SCIP_Real*      cost_org,           /**< (un-biased) edge costs, only needed for PC/RPC */
   const SCIP_Real*      prize,              /**< (possibly biased) vertex prizes */
   SCIP_Bool             costsAreBiased,     /**< biased? */
   int                   start,              /**< start vertex */
   SPATHSPC*             sppc,               /**< PC/MW shortest path data */
   TMBASE*               tmbase              /**< data */
   )
{
#ifdef TM_USE_CSR_PCMW
   SPATHS spaths = { .csr = tmbase->csr, .csr_orgcosts = tmbase->csr_orgcosts, .nodes_dist = tmbase->nodes_dist,
                     .nodes_pred = tmbase->nodes_pred, .dheap = tmbase->dheap, .nodes_isConnected = tmbase->connected };

   assert(graph_pc_isPcMw(g));

   if( g->stp_type == STP_RPCSPG || g->stp_type == STP_RMWCSP )
   {
      shortestpath_computeSteinerTreeRpcMw(g, start, prize, sppc, &spaths);
   }
   else
   {
      shortestpath_computeSteinerTreePcMw(g, start, prize, costsAreBiased, sppc, &spaths);
   }

   SCIP_CALL( solstp_pruneFromTmHeur_csr(scip, g, &spaths, tmbase->result));

#else
   assert(graph_pc_isPcMw(g));

   if( g->stp_type == STP_RPCSPG || g->stp_type == STP_RMWCSP )
   {
      graph_path_st_rpcmw(g, sppc->orderedprizes, sppc->orderedprizes_id,
            tmbase->cost, prize, tmbase->nodes_dist, tmbase->nodes_pred, start, connected);
   }
   else
   {
      graph_path_st_pcmw(g, sppc->orderedprizes, sppc->orderedprizes_id,
            tmbase->cost, prize, costsAreBiased, tmbase->nodes_dist, tmbase->nodes_pred, start, connected);
   }

   SCIP_CALL( solstp_pruneFromTmHeur(scip, g, cost_org, tmbase->result, connected));

#endif


   return SCIP_OKAY;
}



/** Dijkstra based shortest paths heuristic for PCSTP and MWCSP that computes tree spanning all positive
 * vertex weights and subsequently prunes this tree */
static
SCIP_RETCODE computeSteinerTreeDijkPcMwFull(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   const SCIP_Real*      cost_org,           /**< (un-biased) edge costs, only needed for PC/RPC */
   int                   start,              /**< start vertex */
   TMBASE*               tmbase              /**< data */
   )
{
#ifdef TM_USE_CSR_PCMW
   SPATHS spaths = { .csr = tmbase->csr, .csr_orgcosts = tmbase->csr_orgcosts, .nodes_dist = tmbase->nodes_dist,
                     .nodes_pred = tmbase->nodes_pred,
                     .dheap = tmbase->dheap, .nodes_isConnected = tmbase->connected };

   shortestpath_computeSteinerTreePcMwFull(g, start, &spaths);

   SCIP_CALL( solstp_pruneFromTmHeur_csr(scip, g, &spaths, tmbase->result));
#else
   graph_path_st_pcmw_full(g, tmbase->cost, tmbase->nodes_dist, tmbase->nodes_pred, start, connected);
   SCIP_CALL( solstp_pruneFromTmHeur(scip, g, cost_org, tmbase->result, tmbase->connected) );
#endif

   return SCIP_OKAY;
}


/** shortest paths based heuristic
 *  NOTE: DEPRECATED */
static
SCIP_RETCODE computeSteinerTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   int                   start,              /**< start vertex */
   int*                  result,             /**< solution array (on edges) */
   TMALLSP*              tmallsp,            /**< all SP */
   STP_Bool*             connected,          /**< array marking all solution vertices */
   SCIP_RANDNUMGEN*      randnumgen          /**< random number generator */
   )
{
   SCIP_Real min;
   SCIP_Bool directed = (g->stp_type == STP_SAP || g->stp_type == STP_DHCSTP);
   int* perm = NULL;
   int* cluster = NULL;
   int    k;
   int    e;
   int    i;
   int    j;
   int    l;
   int    z;
   int    old;
   int    root;
   int    csize;
   int    newval;
   int    nnodes;
   SCIP_Real** pathdist = tmallsp->pathdist;
   int** pathedge = tmallsp->pathedge;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);

   root = g->source;
   csize = 0;
   nnodes = g->knots;

   assert(0 <= start && start < nnodes);

   SCIPdebugMessage("Heuristic: Start=%5d \n", start);

   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nnodes) );

   cluster[csize++] = start;

   for( e = 0; e < g->edges; e++ )
      result[e] = UNKNOWN;

   for( i = 0; i < nnodes; i++ )
   {
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
      perm[i] = i;
   }

   connected[start] = TRUE;

   SCIPrandomPermuteIntArray(randnumgen, perm, 0, nnodes - 1);

   assert(graph_valid(scip, g));

   /* CONSTCOND */
   for( ;; )
   {
      /* Find a terminal with minimal distance to current ST */
      min = FARAWAY;
      old = -1;
      newval = -1;

      /* time limit exceeded? */
      if( SCIPisStopped(scip) )
         break;

      /* find shortest path from current tree to unconnected terminal */
      for( l = nnodes - 1; l >= 0; --l )
      {
         i = perm[l];
         if( !Is_term(g->term[i]) || connected[i] || !g->mark[i] || (directed && !connected[root] && i != root) )
            continue;

         z = SCIPrandomGetInt(randnumgen, 0, nnodes - 1);

         for( k = csize - 1; k >= 0; k-- )
         {
            j = cluster[(k + z) % csize];
            assert(i != j);
            assert(connected[j]);

            if( SCIPisLT(scip, pathdist[i][j], min) )
            {
               min = pathdist[i][j];
               newval = i;
               old = j;
            }
         }
      }

      /* no new path found? */
      if( newval == -1 )
         break;

      /* mark the new path */
      assert(old > -1);
      assert(newval > -1);
      assert(pathdist[newval] != NULL);
      assert(pathdist[newval][old] < FARAWAY);
      assert(g->term[newval] == 0);
      assert(!connected[newval]);
      assert(connected[old]);

      SCIPdebug(fputc('R', stdout));
      SCIPdebug(fflush(stdout));

      /* start from current tree */
      k = old;

      while( k != newval )
      {
         e = pathedge[newval][k];
         k = g->tail[e];
         if (!connected[k])
         {
            connected[k] = TRUE;
            cluster[csize++] = k;
         }
      }
   }

   /* prune the tree */
   SCIP_CALL( solstp_pruneFromTmHeur(scip, g, cost, result, connected) );

   SCIPfreeBufferArrayNull(scip, &perm);
   SCIPfreeBufferArray(scip, &cluster);

   return SCIP_OKAY;
}

/** cumputes single node Steiner tree; only for rooted PC/MW */
static inline
void computeSteinerTreeSingleNode(
   const GRAPH*          graph,              /**< graph data structure */
   int                   node,               /**< the nodes */
   TMBASE*               tmbase,             /**< data */
   SCIP_Bool*            success             /**< telling name */
)
{
   assert(graph_pc_isPcMw(graph));

   /* only root remaining? */
   if( node == graph->source )
   {
      const int nedges = graph_get_nEdges(graph);
      int *RESTRICT best_result = tmbase->best_result;
      assert(graph_pc_isRootedPcMw(graph));

      for( int e = 0; e < nedges; e++ )
         best_result[e] = UNKNOWN;

      (*success) = TRUE;
   }
}


/** heuristic for degree constrained STPs
 *  NOTE: DEPRECATED */
static
SCIP_RETCODE computeDegConsTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   int                   start,              /**< start vertex */
   int*                  result,             /**< array to indicate whether an edge is in the solution */
   TMALLSP*              tmallsp,            /**< all SP */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   STP_Bool*             connected,          /**< array to indicate whether a vertex is in the solution */
   STP_Bool*             solfound            /**< pointer to store whether a solution has been found */
   )
{
   SCIP_Real min;
   int    csize = 0;
   int    k;
   int    e;
   int    l;
   int    i;
   int    j;
   int    t;
   int    u;
   int    z;
   int    n;
   int    old;
   int    tldegcount;
   int    degcount;
   int    mindegsum;
   int    degmax;
   int    newval;
   int    nnodes;
   int    nterms;
   int    termcount;
   int*   degs;
   int*   maxdegs;
   int* cluster;
   int* perm;
   SCIP_Real** pathdist = tmallsp->pathdist;
   int** pathedge = tmallsp->pathedge;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(g->maxdeg != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);

   z = 0;
   nterms = g->terms;
   nnodes = g->knots;
   mindegsum = 2;
   maxdegs = g->maxdeg;

   SCIPdebugMessage("Heuristic: Start=%5d ", start);

   SCIP_CALL( SCIPallocBufferArray(scip, &degs, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nnodes) );

   cluster[csize++] = start;

   for( i = 0; i < nnodes; i++ )
   {
      degs[i] = 0;
      g->mark[i] = (g->grad[i] > 0);
      connected[i] = FALSE;
      perm[i] = i;
   }

   for( e = 0; e < g->edges; e++ )
   {
      assert(SCIPisGT(scip, cost[e], 0.0));
      result[e] = UNKNOWN;
   }
   connected[start] = TRUE;
   tldegcount = MIN(g->grad[start], maxdegs[start]);
   SCIPrandomPermuteIntArray(randnumgen, perm, 0, nnodes - 1);

   if( Is_term(g->term[start]) )
      termcount = 1;
   else
      termcount = 0;

   for( n = 0; n < nnodes; n++ )
   {
      /* Find a terminal with minimal distance to the current ST */
      min = FARAWAY - 1;
      /* is the free degree sum at most one and are we not connecting the last terminal? */
      if( tldegcount <= 1 && termcount < nterms - 1)
         degmax = 1;
      else
         degmax = 0;
      old = -1;
      newval = -1;
      for( t = 0; t < nnodes; t++ )
      {
         i = perm[t];
         if( !Is_term(g->term[i]) || connected[i] || g->grad[i] == 0 )
            continue;

         z = SCIPrandomGetInt(randnumgen, 0, nnodes - 1);

         for( k = 0; k < csize; k++ )
         {
            j = cluster[(k + z) % csize];
            assert(i != j);
            assert(connected[j]);

            if( SCIPisLE(scip, pathdist[i][j], min) && degs[j] < maxdegs[j])
            {
               u = j;
               degcount = 1;
               while( u != i )
               {
                  u = g->tail[pathedge[i][u]];
                  if( !connected[u] )
                  {
                     if( (MIN(g->grad[u], maxdegs[u]) < 2 || Is_term(g->term[u])) && u != i )
                     {
                        degcount = -2;
                        break;
                     }
                     degcount += MIN(g->grad[u] - 2, maxdegs[u] - 2);
                  }
                  else
                  {
                     assert(u != i);
                     l = g->tail[pathedge[i][u]];
                     if( !connected[l] && degs[u] >= maxdegs[u] )
                     {
                        degcount = -2;
                        break;
                     }
                  }
               }
               if( degcount >= degmax || (degcount >= mindegsum && SCIPisLT(scip, pathdist[i][j], min)) )
               {
                  degmax = degcount;
                  min = pathdist[i][j];
                  newval = i;
                  old = j;
               }
            }
         }
      }

      if( newval == -1 )
      {
         j = UNKNOWN;
         for( k = 0; k < csize; k++ )
         {
            j = cluster[(k + z) % csize];
            if( degs[j] < maxdegs[j] )
               break;
         }
         if( j != UNKNOWN )
         {
            assert(k != csize);

            min = FARAWAY + 1;
            newval = UNKNOWN;
            for( e = g->outbeg[j]; e != EAT_LAST; e = g->oeat[e] )
            {
               u = g->head[e];
               if( !Is_term(g->term[u]) && !connected[u] && SCIPisGE(scip, min, cost[e]) )
               {
                  min = cost[e];
                  k = e;
                  newval = u;
               }
            }
            if( newval != UNKNOWN )
            {
               result[flipedge(k)] = CONNECT;
               degs[newval]++;
               degs[j]++;
               connected[newval] = TRUE;
               cluster[csize++] = newval;
               tldegcount += MIN(maxdegs[newval], g->grad[newval]) - 2;
               continue;
            }

         }
      }
      tldegcount += degmax - 1;
      /* break? */
      if( newval == -1 )
         break;
      /* Weg setzten
       */
      assert(old > -1);
      assert(newval > -1);
      assert(pathdist[newval] != NULL);
      assert(g->term[newval] == 0);
      assert(!connected[newval]);
      assert(connected[old]);

      SCIPdebug(fputc('R', stdout));
      SCIPdebug(fflush(stdout));

      /* mark new tree nodes/edges */
      k = old;
      if( Is_term(g->term[newval]) )
         termcount++;

      while(k != newval)
      {
         e = pathedge[newval][k];
         u = k;
         k = g->tail[e];

         if( !connected[k])
         {
            result[flipedge(e)] = CONNECT;
            degs[u]++;
            degs[k]++;
            connected[k] = TRUE;
            cluster[csize++] = k;
            if( k != newval )
               assert(!Is_term(g->term[k]));
         }
      }
      if( termcount == nterms )
         break;
      assert(degs[newval] == 1);
   }

   *solfound = TRUE;

   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && !connected[i] )
      {
         *solfound = FALSE;
         break;
      }
   }

   if( *solfound )
   {
      /* prune the solution */
      SCIP_CALL( SCIPStpHeurTMBuildTreeDc(scip, g, result, connected) );

      for( t = 0; t < nnodes; t++ )
         if( degs[t] > maxdegs[t] )
            *solfound = FALSE;
   }
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &cluster);
   SCIPfreeBufferArray(scip, &degs);
   return SCIP_OKAY;
}

/** Voronoi based shortest path heuristic.
 *  NOTE: DEPRECATED */
static
SCIP_RETCODE computeSteinerTreeVnoi(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   int                   start,              /**< start vertex */
   STP_Bool              firstrun,           /**< method called for the first time? (during one heuristic round) */
   TMVNOI*               tmvnoi,             /**< TM data */
   int*                  result,             /**< array to indicate whether an edge is in the solution */
   STP_Bool*             connected           /**< array to indicate whether a vertex is in the solution */
   )
{
   int *vcount;
   int k;
   int i;
   int j;
   int best;
   int term;
   int count;
   int nnodes;
   int nterms;
   SCIP_Real **node_dist = tmvnoi->node_dist;
   SCIP_PQUEUE *pqueue = tmvnoi->pqueue;
   GNODE **gnodearr = tmvnoi->gnodearr;
   int *nodenterms = tmvnoi->nodenterms;
   int **node_base = tmvnoi->node_base;
   int **node_edge = tmvnoi->node_edge;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   nnodes = g->knots;
   nterms = g->terms;

   SCIP_CALL( SCIPallocBufferArray(scip, &vcount, nnodes) );

   SCIPdebugMessage("TM_Polzin Heuristic: Start=%5d ", start);

   /* if the heuristic is called for the first time several data structures have to be set up */
   if( firstrun )
   {
      PATH* vnoi;
      SCIP_Real* vcost;
      int old;
      int oedge;
      int root = g->source;
      int   ntovisit;
      int   nneighbnodes;
      int   nneighbterms;
      int   nreachednodes;
      int*  state;
      int*  vbase;
      int*  terms;
      int*  tovisit;
      int*  reachednodes;
      STP_Bool* termsmark;
      STP_Bool* visited;
      int e;
      /* PHASE I: */
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

      /* allocate memory needed in PHASE I */
      SCIP_CALL( SCIPallocBufferArray(scip, &terms, nterms) );
      SCIP_CALL( SCIPallocBufferArray(scip, &termsmark, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &visited, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &reachednodes, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tovisit, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vcost, nnodes) );

      j = 0;
      for( i = 0; i < nnodes; i++ )
      {
         visited[i] = FALSE;
         if( Is_term(g->term[i]) )
         {
            termsmark[i] = TRUE;
            terms[j++] = i;
         }
         else
         {
            termsmark[i] = FALSE;
         }
      }

      for( e = 0; e < g->edges; e++)
      {
         assert(SCIPisGE(scip, cost[e], 0.0));
         assert(SCIPisGE(scip, costrev[e], 0.0));
      }

      assert(j == nterms);
      graph_voronoi(scip, g, cost, costrev, termsmark, vbase, vnoi);
      state = g->path_state;

      for( i = 0; i < nnodes; i++ )
         if( Is_term(g->term[i]) )
            assert(vbase[i] == i);

      for( k = 0; k < nnodes; k++ )
      {
         connected[k] = FALSE;
         vcount[k] = 0;
         gnodearr[k]->number = k;
         if( !Is_term(g->term[k]) )
         {
            node_dist[k][0] = vnoi[k].dist;
            node_edge[k][0] = vnoi[k].edge;
            node_base[k][0] = vbase[k];
            nodenterms[k] = 1;
         }
         else
         {
            nodenterms[k] = 0;
            node_edge[k][0] = UNKNOWN;
            termsmark[k] = FALSE;
         }
         state[k] = UNKNOWN;
         vcost[k] = vnoi[k].dist;
         vnoi[k].dist = FARAWAY;
      }

      /* for each terminal: extend the voronoi regions until all neighbouring terminals have been visited */
      for( i = 0; i < nterms; i++ )
      {
         term = terms[i];
         nneighbterms = 0;
         nneighbnodes = 0;
         nreachednodes = 0;
         for( k = 0; k < nnodes; k++ )
            assert(termsmark[k] == FALSE);
         /* DFS (starting from terminal i) until the entire voronoi region has been visited */
         tovisit[0] = term;
         ntovisit = 1;
         visited[term] = TRUE;
         state[term] = CONNECT;
         while( ntovisit > 0 )
         {
            /* iterate all incident edges */
            old = tovisit[--ntovisit];

            for( oedge = g->outbeg[old]; oedge != EAT_LAST; oedge = g->oeat[oedge] )
            {
               k = g->head[oedge];

               /* is node k in the voronoi region of the i-th terminal ? */
               if( vbase[k] == term )
               {
                  if( !visited[k] )
                  {
                     state[k] = CONNECT;
                     assert(nnodes - (nneighbnodes + 1) > ntovisit);
                     tovisit[ntovisit++] = k;
                     visited[k] = TRUE;
                     reachednodes[nreachednodes++] = k;
                  }
               }
               else
               {
                  if( !visited[k] )
                  {
                     visited[k] = TRUE;
                     vnoi[k].dist = vcost[old] + ((vbase[k] == root)? cost[oedge] : costrev[oedge]);
                     vnoi[k].edge = oedge;

                     if( termsmark[vbase[k]] == FALSE )
                     {
                        termsmark[vbase[k]] = TRUE;
                        nneighbterms++;
                     }
                     assert(nnodes - (nneighbnodes + 1) > ntovisit - 1);
                     tovisit[nnodes - (++nneighbnodes)] = k;
                  }
                  else
                  {
                     /* if edge 'oedge' allows a shorter connection of node k, update */
                     if( SCIPisGT(scip, vnoi[k].dist, vcost[old] + ((vbase[k] == root)? cost[oedge] : costrev[oedge])) )
                     {
                        vnoi[k].dist = vcost[old] + ((vbase[k] == root)? cost[oedge] : costrev[oedge]);
                        vnoi[k].edge = oedge;
                     }
                  }
               }
            }
         }

         count = 0;
         for( j = 0; j < nneighbnodes; j++ )
         {
            assert(termsmark[vbase[tovisit[nnodes - j - 1]]]);
            graph_pathHeapAdd(vnoi, tovisit[nnodes - j - 1], g->path_heap, state, &count);
         }
         SCIP_CALL( graph_voronoiExtend(scip, g, ((term == root)? cost : costrev), vnoi, node_dist, node_base, node_edge, termsmark, reachednodes, &nreachednodes, nodenterms,
               nneighbterms, term, nneighbnodes) );

         reachednodes[nreachednodes++] = term;

         for( j = 0; j < nreachednodes; j++ )
         {
            vnoi[reachednodes[j]].dist = FARAWAY;
            state[reachednodes[j]] = UNKNOWN;
            visited[reachednodes[j]] = FALSE;
         }

         for( j = 0; j < nneighbnodes; j++ )
         {
            vnoi[tovisit[nnodes - j - 1]].dist = FARAWAY;
            state[tovisit[nnodes - j - 1]] = UNKNOWN;
            visited[tovisit[nnodes - j - 1]] = FALSE;
         }
      }

      /* for each node v: sort the terminal arrays according to their distance to v */
      for( i = 0; i < nnodes && !SCIPisStopped(scip); i++ )
         SCIPsortRealIntInt(node_dist[i], node_base[i], node_edge[i], nodenterms[i]);

      /* free memory */
      SCIPfreeBufferArray(scip, &vcost);
      SCIPfreeBufferArray(scip, &tovisit);
      SCIPfreeBufferArray(scip, &vbase);
      SCIPfreeBufferArray(scip, &reachednodes);
      SCIPfreeBufferArray(scip, &visited);
      SCIPfreeBufferArray(scip, &vnoi);
      SCIPfreeBufferArray(scip, &termsmark);
      SCIPfreeBufferArray(scip, &terms);
   }

   /* PHASE II */
   else
   {
      for( k = 0; k < nnodes; k++ )
      {
         connected[k] = FALSE;
         vcount[k] = 0;
      }
   }

   connected[start] = TRUE;
   gnodearr[start]->dist = node_dist[start][0];
   SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[start]) );

   while( SCIPpqueueNElems(pqueue) > 0 && !SCIPisStopped(scip) )
   {
      best = ((GNODE*) SCIPpqueueRemove(pqueue))->number;

      term = node_base[best][vcount[best]];
      assert( Is_term(g->term[term]) );
      /* has the terminal already been connected? */
      if( !connected[term] )
      {
         /* connect the terminal */
         k = g->tail[node_edge[best][vcount[best]]];
         while( k != term )
         {
            j = 0;

            while( node_base[k][vcount[k] + j] != term )
               j++;

            assert(vcount[k] + j < nodenterms[k]);

            if( !connected[k] )
            {
               assert(vcount[k] == 0);

               connected[k] = TRUE;
               while( vcount[k] < nodenterms[k] && connected[node_base[k][vcount[k]]] )
               {
                  vcount[k]++;
                  j--;
               }

               if( vcount[k] < nodenterms[k] )
               {
                  gnodearr[k]->dist = node_dist[k][vcount[k]];
                  SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
               }
            }

            assert( vcount[k] + j < nodenterms[k] );
            assert(node_base[k][vcount[k] + j] == term);
            k = g->tail[node_edge[k][vcount[k] + j]];
         }

         /* finally, connected the terminal */
         assert( k == term );
         assert( !connected[k] );
         connected[k] = TRUE;

         assert( vcount[k] == 0 );
         while( vcount[k] < nodenterms[k] && connected[node_base[k][vcount[k]]] )
         {
            vcount[k]++;
         }
         if( vcount[k] < nodenterms[k] )
         {
            gnodearr[k]->dist = node_dist[k][vcount[k]];
            SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
         }
      }

      while( vcount[best] + 1 < nodenterms[best] )
      {
         if( !connected[node_base[best][++vcount[best]]] )
         {
            gnodearr[best]->dist = node_dist[best][vcount[best]];
            SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[best]) );
            break;
         }
      }
   }

   /* prune the ST, so that all leaves are terminals */
   SCIP_CALL( solstp_pruneFromTmHeur(scip, g, cost, result, connected) );

   SCIPfreeBufferArray(scip, &vcount);

   return SCIP_OKAY;
}


static
void initCostsAndPrioLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heurdata */
   SCIP_VAR**            vars,               /**< variables */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Real             randupper,          /**< random value bound */
   SCIP_Real             randlower,          /**< random value bound */
   const SCIP_Real*      xval,               /**< xval */
   SCIP_Real*            nodepriority,       /**< node priority (uninitialized) */
   SCIP_Real*            prize,              /**< prize (uninitialized) or NULL */
   SCIP_Real*            cost                /**< arc costs (uninitialized) */
)
{
   SCIP_Bool partrand = FALSE;
   SCIP_Bool totalrand = FALSE;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;

   if( (heurdata->nlpiterations == SCIPgetNLPIterations(scip) && SCIPrandomGetInt(heurdata->randnumgen, 0, 5) != 1)
         || SCIPrandomGetInt(heurdata->randnumgen, 0, 15) == 5 )
      partrand = TRUE;

   if( !partrand && (heurdata->nlpiterations == SCIPgetNLPIterations(scip)) )
      totalrand = TRUE;
   else if( graph->stp_type == STP_DCSTP && heurdata->ncalls != 1 && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 1
         && (graph->maxdeg[graph->source] == 1 || SCIPrandomGetInt(heurdata->randnumgen, 0, 5) == 5) )
   {
      totalrand = TRUE;
      partrand = FALSE;
   }

   for( int k = 0; k < nnodes; k++ )
      nodepriority[k] = 0.0;

   if( graph->stp_type != STP_MWCSP && graph->stp_type != STP_RMWCSP )
   {
      for( int e = 0; e < nedges; e++ )
      {
         nodepriority[graph->head[e]] += xval[e];
         nodepriority[graph->tail[e]] += xval[e];
      }
   }

   if( graph->stp_type == STP_DHCSTP )
   {
      for( int e = 0; e < nedges; e++ )
      {
         if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[flipedge_Uint(e)]) < 0.5 )
         {
            cost[e] = BLOCKED;
         }
         else
         {
            if( totalrand )
            {
               const SCIP_Real randval = SCIPrandomGetReal(heurdata->randnumgen, randlower, randupper);
               cost[e] = graph->cost[e] * randval;
            }
            else
            {
               cost[e] = ((1.0 - xval[e]) * graph->cost[e]);
            }
         }
         if( partrand )
         {
            const SCIP_Real randval = SCIPrandomGetReal(heurdata->randnumgen, randlower, randupper);
            cost[e] = cost[e] * randval;
         }

         assert(SCIPisGE(scip, cost[e], 0.0));
      }
   }
   else
   {
      /* swap costs; set a high cost if the variable is fixed to 0 */
      if( graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP )
      {
         for( int e = 0; e < nedges; e++ )
            nodepriority[graph->head[e]] += xval[e];

         for( int e = 0; e < nedges; e++ )
         {
            if( graph->cost[e] >= FARAWAY )
               cost[e] = graph->cost[e];
            else if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[flipedge_Uint(e)]) < 0.5 )
               cost[e] = MAX(graph->cost[e], BLOCKED);
            else
               cost[e] = graph->cost[e] * (1.0 - MIN(1.0, nodepriority[graph->head[e]]));
         }

         for( int e = 0; e < nedges; e++ )
            nodepriority[graph->tail[e]] += xval[e];
      }
      else
      {
         for( int e = 0; e < nedges; e++ )
         {
            SCIP_Real randval = 1.0;
            if( totalrand || partrand )
               randval = SCIPrandomGetReal(heurdata->randnumgen, randlower, randupper);

            if( SCIPvarGetUbLocal(vars[e]) < 0.5 && SCIPvarGetUbLocal(vars[flipedge_Uint(e)]) < 0.5 )
            {
               cost[e] = MAX(graph->cost[e], BLOCKED);
               continue;
            }

            if( totalrand )
               cost[e] = graph->cost[e] * randval;
            else
               cost[e] = ((1.0 - xval[e]) * graph->cost[e]);

            if( partrand )
               cost[e] = cost[e] * randval;
         }
#if 0
         for( int e = 0; e < nedges; e++ )
         if( graph->cost[e] == graph->cost[flipedge_Uint(e)] )
         cost[e] = MIN(cost[e], cost[flipedge_Uint(e)]);
#endif
      } /* graph->stp_type != STP_MWCSP && graph->stp_type != STP_RMWCSP */
   } /* graph->stp_type != STP_DHCSTP */


   if( graph_pc_isPcMw(graph) && prize )
   {
      assert(graph->extended);

      for( int k = 0; k < nnodes; k++ )
         prize[k] = graph->prize[k];

      for( int k = 0; k < nnodes; k++ )
      {
         if( Is_pseudoTerm(graph->term[k]) )
         {
            const int term = graph_pc_getTwinTerm(graph, k);
            const int rootedge = graph_pc_getRoot2PtermEdge(graph, term);

            prize[k] = cost[rootedge];

            assert(prize[k] >= 0.0);
         }
      }
   }
}


/** initializes for TM SP runs (computes shortest paths from all terminals) */
static
void buildTmAllSp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< arc costs */
   const SCIP_Real*      costrev,            /**< reversed arc costs */
   TMALLSP*              tmallsp             /**< TM data */
)
{
   SCIP_Real** pathdist = tmallsp->pathdist;
   int** pathedge = tmallsp->pathedge;
   const int nnodes = graph_get_nNodes(graph);
   const int root = graph->source;

   assert(pathdist != NULL);
   assert(pathedge != NULL);

   for( int k = 0; k < nnodes; k++ )
      graph->mark[k] = (graph->grad[k] > 0);

   /* initialize shortest paths from all terminals */
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(graph->term[k]) )
      {
         if( root == k )
            graph_path_execX(scip, graph, k, cost,  pathdist[k], pathedge[k]);
         else
            graph_path_execX(scip, graph, k, costrev, pathdist[k], pathedge[k]);
      }
   }
}

/** initial runs for DHCSTP to figure out a good hop factor */
static
SCIP_RETCODE dhcstpWarmUp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_Real*            cost,               /**< hacky */
   SCIP_Real*            costrev,            /**< hacky */
   SCIP_Real*            hopfactor,          /**< (in/out) */
   TMBASE*               tmbase,             /**< (in/out) */
   SCIP_Bool*            success             /**< success? */
)
{
   STP_Bool* connected;
   SCIP_Real* RESTRICT orgcost;
   int* RESTRICT result = tmbase->result;
   int* RESTRICT best_result = tmbase->best_result;
   SCIP_Real hopfactor_local;
   SCIP_Real hopfactor_best = -1.0;
   SCIP_Real maxcost = 0.0;
   const int nedges = graph_get_nEdges(graph);
   const int root = graph->source;
   const int nnodes = graph_get_nNodes(graph);

   assert(hopfactor && success);
   assert(*success == FALSE);
   assert(graph->stp_type == STP_DHCSTP);

   if( LE(*hopfactor, 0.0 ) )
      *hopfactor = DEFAULT_HOPFACTOR;

   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orgcost, nedges) );

   for( int e = 0; e < nedges; e++ )
   {
      if( GT(cost[e], maxcost) && LT(cost[e], BLOCKED) )
         maxcost = cost[e];
   }

   assert(GT(maxcost, 0.0));

   hopfactor_local = *hopfactor;

   BMScopyMemoryArray(orgcost, cost, nedges);

   /* do a warm-up run */
   for( int r = 0; r < 10; r++ )
   {
      const SCIP_Real* const gcost = graph->cost;
      SCIP_Real obj = 0.0;
      int edgecount = 0;
      SCIP_Bool lsuccess = FALSE;

      assert(GT(hopfactor_local, 0.0));

      for( int e = 0; e < nedges; e++ )
      {
         if( LT(cost[e], BLOCKED) )
            cost[e] = 1.0 + orgcost[e] / (hopfactor_local * maxcost);

         result[e] = UNKNOWN;
      }

      SCIP_CALL( computeSteinerTreeDijk(scip, graph, root, tmbase) );

      for( int e = 0; e < nedges; e++)
      {
         if( result[e] == CONNECT )
         {
            obj += gcost[e];
            edgecount++;
         }
      }

      if( SCIPisLT(scip, obj, tmbase->best_obj) && edgecount <= graph->hoplimit )
      {
         tmbase->best_obj = obj;

         for( int e = 0; e < nedges; e++ )
            best_result[e] = result[e];

         (*success) = TRUE;
         lsuccess = TRUE;
         hopfactor_best = hopfactor_local;
      }

      if( !lsuccess || SCIPisGT(scip, fabs((double) edgecount - graph->hoplimit) / (double) graph->hoplimit, 0.05) )
      {
         if( !lsuccess )
         {
            if( (*success) )
            {
               hopfactor_local = hopfactor_local * (1.0 + fabs((double) edgecount - graph->hoplimit) / (double) graph->hoplimit);
            }
            else
            {
               hopfactor_local = hopfactor_local * (1.0 + 3 * fabs((double) edgecount - graph->hoplimit) / (double) graph->hoplimit);
               hopfactor_best = hopfactor_local;
            }
         }
         else
         {
            hopfactor_local = hopfactor_local / (1.0 + fabs((double) edgecount - graph->hoplimit) / (double) graph->hoplimit);
         }

         assert(SCIPisGT(scip, hopfactor_local, 0.0));
      }
      else
      {
         break;
      }
   }

   assert(SCIPisGT(scip, hopfactor_best, 0.0));

   (*hopfactor) = hopfactor_best;

   for( int e = 0; e < nedges; e++ )
      if( (LT(cost[e], BLOCKED) ) )
         cost[e] = 1.0 + orgcost[e] / (hopfactor_best * maxcost);

   for( int e = 0; e < nedges; e++)
      costrev[e] = cost[flipedge(e)];

   SCIPfreeBufferArray(scip, &orgcost);
   SCIPfreeBufferArray(scip, &connected);

   return SCIP_OKAY;
}


/** submethod for SCIPStpHeurTMRun */
static
SCIP_RETCODE runTm(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   TMBASE*               tmbase,             /**< TM base data */
   TMALLSP*              tmallsp,            /**< TM data */
   TMVNOI*               tmvnoi,             /**< TM data */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
)
{
   SCIP_HEURDATA* heurdata = getTMheurData(scip);
   const int* const start = tmbase->startnodes;
   int* const result = tmbase->result;
   const SCIP_Real* cost = tmbase->cost;
   const SCIP_Real* costrev = tmbase->costrev;
   const int mode = getTmMode(heurdata, graph);
   const int nedges = graph_get_nEdges(graph);
   const int runs = tmbase->nruns;
   STP_Bool solfound = FALSE;

   assert(!graph_pc_isPcMw(graph));
   assert(runs >= 1);

   if( graph->stp_type == STP_DCSTP || mode == TM_SP  )
   {
      buildTmAllSp(scip, graph, cost, costrev, tmallsp);
   }

   /* main loop */
   for( int r = 0; r < runs; r++ )
   {
      SCIP_Real obj;

      assert(start[r] >= 0 && start[r] < graph->knots);
      assert(graph->stp_type != STP_NWPTSPG || !graph_nw_knotIsLeaf(graph, start[r]));

      if( graph->stp_type == STP_DCSTP )
      {
         SCIP_CALL( computeDegConsTree(scip, graph, cost, costrev, start[r], result, tmallsp,  heurdata->randnumgen, tmbase->connected, &solfound) );
      }
      else if( mode == TM_DIJKSTRA )
      {
#if 0
         SCIP_CALL( computeSteinerTreeDijk(scip, graph, start[r], tmbase) );
#else
         SCIP_CALL( computeSteinerTreeCsr(scip, graph, start[r], tmbase) );
#endif
      }
      else if( mode == TM_SP )
      {
         SCIP_CALL( computeSteinerTree(scip, graph, cost, costrev, start[r], result, tmallsp, tmbase->connected, heurdata->randnumgen) );
      }
      else
      {
         SCIP_CALL( computeSteinerTreeVnoi(scip, graph, cost, costrev, start[r], (r == 0), tmvnoi, result, tmbase->connected) );
      }

      /* here another measure than in the do_(...) heuristics is being used */
      obj = solstp_getObj(graph, result, 0.0, nedges);

      SCIPdebugMessage("run=%d, obj=%.12e\n", r, obj);

      if( SCIPisLT(scip, obj, tmbase->best_obj) && (graph->stp_type != STP_DCSTP || solfound) )
      {
         if( graph->stp_type != STP_DHCSTP || solstp_getNedges(graph, result) <= graph->hoplimit )
         {
            tmbase->best_obj = obj;
            BMScopyMemoryArray(tmbase->best_result, result, nedges);
            (*success) = TRUE;
         }
      }

      /* stop early? */
      if( SCIPisStopped(scip) )
         break;
   }

   return SCIP_OKAY;
}


/** submethod for SCIPStpHeurTMRun in PC or MW mode */
static
SCIP_RETCODE runTmPcMW_mode(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< arc costs */
   const SCIP_Real*      prize,              /**< prizes (for PCMW) or NULL */
   enum PCMW_TmMode      pcmwmode,           /**< mode */
   int                   bestincstart,       /**< best incumbent start vertex */
   SCIP_Real*            nodepriority,       /**< vertex priorities for vertices to be starting points (NULL for no priorities) */
   TMBASE*               tmbase,             /**< data */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
)
{
   SPATHSPC* sppc;
   SCIP_Real* costorg = NULL;
   SCIP_Real* costfullbiased = NULL;
   SCIP_Real* prizefullbiased = NULL;
   const SCIP_Real* const prize_in = (prize) ? prize : graph->prize;
   int* startnodes = NULL;
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   int maxruns = tmbase->nruns;

   assert(graph->extended);
   assert(pcmwmode != pcmode_all);
   assert(graph_valid(scip, graph));

#ifndef TM_USE_CSR_PCMW
   if( graph_pc_isPc(graph) )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &costorg, nedges) );
      graph_pc_getOrgCosts(scip, graph, costorg);
   }
#endif

   if( pcmwmode == pcmode_biasfull )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &costfullbiased, nedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &prizefullbiased, nnodes) );
      graph_pc_getBiased(scip, graph, costfullbiased, prizefullbiased);
      SCIP_CALL( shortestpath_pcInit(scip, graph, costfullbiased, prizefullbiased, &sppc) );
   }
   else
   {
      SCIP_CALL( shortestpath_pcInit(scip, graph, NULL, prize_in, &sppc) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &startnodes, graph->terms) );
   pcmwGetStartNodes(scip, nodepriority, maxruns, bestincstart, graph, startnodes);

   pcmwSetEdgeCosts(graph, cost, costorg, costfullbiased, pcmwmode, tmbase);

   /* main loop */
   for( int r = 0; r < maxruns; r++ )
   {
      const int start = startnodes[r];
      SCIPdebugMessage("TM run=%d start=%d\n", r, start);

      if( graph->grad[start] == 0 )
      {
         computeSteinerTreeSingleNode(graph, start, tmbase, success);
         continue;
      }
      else if( graph->grad[start] == 1 )
      {
         if( !graph_pc_isRootedPcMw(graph) || !graph_pc_knotIsFixedTerm(graph, start) )
            continue;
      }

      if( pcmwmode == pcmode_fulltree )
      {
         SCIP_CALL( computeSteinerTreeDijkPcMwFull(scip, graph, costorg, start, tmbase) );
      }
      else if( pcmwmode == pcmode_biasfull )
      {
         SCIP_CALL( computeSteinerTreeDijkPcMw(scip, graph,
               costorg, prizefullbiased, TRUE, start, sppc, tmbase) );
      }
      else if( pcmwmode == pcmode_simple )
      {
         assert(graph_pc_isPc(graph));
         SCIP_CALL( computeSteinerTreeDijkPcMw(scip, graph,
               costorg, prize_in, FALSE, start, sppc, tmbase) );
      }
      else
      {
         assert(pcmwmode == pcmode_bias);
         SCIP_CALL( computeSteinerTreeDijkPcMw(scip, graph,
               costorg, prize_in, TRUE, start, sppc, tmbase) );
      }

      /* update the best known solution if possible */
      pcmwUpdateBestSol(scip, graph, tmbase, success);

      if( SCIPisStopped(scip) )
         break;
   }

   shortestpath_pcFree(scip, &sppc);

   SCIPfreeBufferArray(scip, &startnodes);
   SCIPfreeBufferArrayNull(scip, &prizefullbiased);
   SCIPfreeBufferArrayNull(scip, &costfullbiased);
   SCIPfreeBufferArrayNull(scip, &costorg);

   return SCIP_OKAY;
}


/** SCIPStpHeurTMRun in PC or MW mode */
static inline
SCIP_RETCODE runTmPcMW(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< arc costs */
   const SCIP_Real*      prize,              /**< prizes (for PCMW) or NULL */
   enum PCMW_TmMode      pcmw_tmmode,        /**< mode */
   int                   beststart,          /**< best incumbent start vertex */
   SCIP_Real*            nodepriority,       /**< vertex priorities for vertices to be starting points (NULL for no priorities) */
   TMBASE*               tmbase,             /**< data */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
)
{
   SCIP_HEURDATA* heurdata = getTMheurData(scip);
   const enum PCMW_TmMode pcmw_mode = (pcmw_tmmode == pcmode_fromheurdata) ? heurdata->pcmw_mode : pcmw_tmmode;
   const SCIP_Bool run_all = (pcmw_mode == pcmode_all);

   assert(graph_pc_isPcMw(graph));
   assert(graph_pc_isPc(graph) || pcmw_mode != pcmode_simple);
   assert(pcmw_mode != pcmode_fromheurdata);

   if( graph_pc_isPc(graph) && (pcmw_mode == pcmode_simple || run_all) )
   {
      SCIP_CALL( runTmPcMW_mode(scip, graph, cost, prize, pcmode_simple, beststart, nodepriority, tmbase, success));
   }
   if( pcmw_mode == pcmode_bias || pcmw_mode == pcmode_biasAndFulltree || run_all )
   {
      SCIP_CALL( runTmPcMW_mode(scip, graph, cost, prize, pcmode_bias, beststart, nodepriority, tmbase, success));
   }
   if( pcmw_mode == pcmode_biasfull || run_all )
   {
      SCIP_CALL( runTmPcMW_mode(scip, graph, cost, prize, pcmode_biasfull, beststart, nodepriority, tmbase, success));
   }
   if( pcmw_mode == pcmode_fulltree || pcmw_mode == pcmode_biasAndFulltree || run_all )
   {
      SCIP_CALL( runTmPcMW_mode(scip, graph, cost, prize, pcmode_fulltree, beststart, nodepriority, tmbase, success));
   }

   return SCIP_OKAY;
}


/** submethod for SCIPStpHeurTMRun */
static
SCIP_RETCODE runTmDhcstp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   TMBASE*               tmbase,             /**< TM base data */
   TMALLSP*              tmallsp,            /**< TM data */
   TMVNOI*               tmvnoi,             /**< TM data */
   SCIP_Real*            hopfactor,          /**< hopfactor */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
)
{
   SCIP_Real* cost;
   SCIP_Real* costrev;
   const SCIP_Real* cost_org;
   const SCIP_Real* costrev_org;
   const int nedges = graph_get_nEdges(graph);

   assert(graph->stp_type == STP_DHCSTP);

   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   BMScopyMemoryArray(cost, tmbase->cost, nedges);
   BMScopyMemoryArray(costrev, tmbase->costrev, nedges);

   cost_org = tmbase->cost;
   costrev_org = tmbase->costrev;
   tmbase->cost = cost;
   tmbase->costrev = costrev;

   SCIP_CALL( dhcstpWarmUp(scip, graph, cost, costrev, hopfactor, tmbase, success) );
   SCIP_CALL( runTm(scip, graph, tmbase, tmallsp, tmvnoi, success) );

   tmbase->cost = cost_org;
   tmbase->costrev = costrev_org;
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTM)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPStpIncludeHeurTM(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeRandom(scip, &heurdata->randnumgen);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_PROBDATA* probdata;
   GRAPH* graph;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   if( graph == NULL )
   {
      heurdata->stp_type = STP_SPG;
      return SCIP_OKAY;
   }
   heurdata->stp_type = graph->stp_type;
   heurdata->beststartnode = -1;
   heurdata->ncalls = 0;
   heurdata->nlpiterations = -1;
   heurdata->nexecs = 0;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTM)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_PROBDATA* probdata;
   SCIP_HEURDATA* heurdata;
   GRAPH* graph;
   int* soledges;
   int runs;
   int nedges;
   SCIP_Bool success = FALSE;

   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   if( graph->stp_type == STP_BRMWCSP )
      return SCIP_OKAY;

   runs = 0;

   /* set the runs, i.e. number of different starting points for the heuristic */
   if( heurtiming & SCIP_HEURTIMING_BEFORENODE )
   {
      if( SCIPgetDepth(scip) > 0 )
         return SCIP_OKAY;

      runs = heurdata->initruns;
   }
   else if( ((heurtiming & SCIP_HEURTIMING_DURINGLPLOOP) && (heurdata->ncalls % heurdata->duringlpfreq == 0)) || (heurtiming & SCIP_HEURTIMING_AFTERLPLOOP) )
   {
      runs = heurdata->evalruns;

      if( graph_pc_isPcMw(graph) )
         runs *= 2;
   }
   else if( heurtiming & SCIP_HEURTIMING_AFTERNODE )
   {
      if( SCIPgetDepth(scip) == 0 )
         runs = heurdata->rootruns;
      else
         runs = heurdata->leafruns;
   }

   /* increase counter for number of (TM) calls */
   heurdata->ncalls++;

   if( runs == 0 )
      return SCIP_OKAY;

   heurdata->nexecs++;

   SCIPdebugMessage("Heuristic Start\n");

   /* get all variables (corresponding to the edges) */
   vars = SCIPprobdataGetVars(scip);
   if( vars == NULL )
      return SCIP_OKAY;

   assert(vars[0] != NULL);

   nedges = graph->edges;

   /* allocate memory */
   SCIP_CALL(SCIPallocBufferArray(scip, &soledges, nedges));

   *result = SCIP_DIDNOTFIND;

   /* call the actual heuristic */
   SCIP_CALL( SCIPStpHeurTMRunLP(scip, graph, heur, soledges, runs, &success) );

   if( success )
   {
      SCIP_Real* nval;
      const int nvars = SCIPprobdataGetNVars(scip);

      assert(nvars == nedges);

      SCIP_CALL(SCIPallocBufferArray(scip, &nval, nvars));

      for( int v = 0; v < nvars; v++ )
      {
         nval[v] = (soledges[v] == CONNECT) ? 1.0 : 0.0;
      }

      SCIP_CALL( SCIPStpValidateSol(scip, graph, nval, FALSE, &success) );

      assert(success);

      if( success )
      {
         SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, heur, &success) );

         if( success )
         {
            SCIPdebugMessage("TM solution added, value %f \n",
                  solstp_getObj(graph, soledges, SCIPprobdataGetOffset(scip), nedges));

            *result = SCIP_FOUNDSOL;
         }
      }
      SCIPfreeBufferArray(scip, &nval);
   }

   heurdata->nlpiterations = SCIPgetNLPIterations(scip);
   SCIPfreeBufferArray(scip, &soledges);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */



/** compute starting points among marked (w.r.t. g->mark) vertices for constructive heuristics */
void SCIPStpHeurTMCompStarts(
   GRAPH*                graph,              /**< graph data structure */
   int*                  starts,             /**< starting points array */
   int*                  runs                /**< pointer to number of runs */
   )
{
   int r;
   int l;
   int root;
   int nruns;
   int nnodes;
   int nterms;
   int randval;

   assert(runs != NULL);
   assert(graph != NULL);
   assert(starts != NULL);

   nruns = *runs;
   root = graph->source;
   nnodes = graph->knots;
   nterms = graph->terms;

   r = 0;
   if( graph->mark[root] && nruns > 0 )
      starts[r++] = root;

   randval = nnodes - nterms;

   /* use non-isolated terminals as starting points for TM heuristic */
   for( int k = 0; k < nnodes; k++ )
   {
      if( r >= nruns || r >= nterms )
         break;

      l = (k + randval) % nnodes;
      if( Is_term(graph->term[l]) && graph->mark[l] && l != root )
         starts[r++] = l;
   }

   /* fill empty slots randomly */
   for( int k = 0; k < nnodes && r < nruns; k++ )
   {
      l = (k + randval) % nnodes;
      if( !Is_term(graph->term[l]) && graph->mark[l] )
         starts[r++] = l;
   }

   *runs = r;
}


/** build Steiner tree in such a way that all leaves are terminals */
SCIP_RETCODE SCIPStpHeurTMBuildTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   PATH*                 mst,                /**< path data structure array */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            objresult,          /**< pointer to store objective value of result */
   int*                  connected           /**< CONNECT/UNKNOWN */
   )
{
   SCIP_Real obj;
   int i;
   int j;
   int e1;
   int root;
   int count;
   int nnodes;

   assert(g != NULL);
   assert(mst != NULL);
   assert(scip != NULL);
   assert(cost != NULL);
   assert(connected != NULL);
   assert(!graph_pc_isPcMw(g));

   obj = 0.0;
   root = g->source;
   nnodes = g->knots;

   /* compute the MST */
   for( i = nnodes - 1; i >= 0; --i )
      g->mark[i] = (connected[i] == CONNECT);

   graph_path_exec(scip, g, MST_MODE, root, cost, mst);

   /* prune */
   do
   {
      count = 0;

      for( i = nnodes - 1; i >= 0; --i )
      {
         if( !g->mark[i] || Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
         {
            e1 = mst[g->head[j]].edge;
            if( e1 == j )
               break;
         }

         if( j == EAT_LAST )
         {
            mst[i].edge = UNKNOWN;
            g->mark[i] = FALSE;
            connected[i] = UNKNOWN;
            count++;
            break;
         }
      }
   }
   while( count > 0 );

   for( i = nnodes - 1; i >= 0; --i )
   {
      if( mst[i].edge >= 0 )
         obj += cost[mst[i].edge];
      else if( Is_term(g->term[i]) && i != root )
      {
         obj = FARAWAY;
         break;
      }
   }

   *objresult = obj;

   return SCIP_OKAY;
}



/** build (rooted) prize collecting Steiner tree in such a way that all leaves are terminals; objresult is set FARAWAY if infeasible */
SCIP_RETCODE SCIPStpHeurTMBuildTreePcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   PATH*                 mst,                /**< path data structure array */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            objresult,          /**< pointer to store objective value of result */
   int*                  connected           /**< CONNECT/UNKNOWN */
   )
{
   SCIP_Real obj;
   int mstroot;
   int count;
   const int nnodes = g->knots;
   const int orgroot = g->source;
   const SCIP_Bool isrooted = graph_pc_isRootedPcMw(g);

   assert(g != NULL);
   assert(mst != NULL);
   assert(scip != NULL);
   assert(cost != NULL);
   assert(connected != NULL);
   assert(g->extended);

   obj = 0.0;

   /* unmark all dummy terminals and unconnected nodes */
   for( int i = nnodes - 1; i >= 0; --i )
   {
      if( connected[i] == CONNECT && !Is_term(g->term[i]) )
         g->mark[i] = TRUE;
      else
         g->mark[i] = FALSE;

      if( isrooted && graph_pc_knotIsFixedTerm(g, i) )
         g->mark[i] = TRUE;
   }

   if( isrooted )
   {
      mstroot = orgroot;
      assert(g->mark[mstroot]);
   }
   else
   {
      int a;
      int i;
      for( a = g->outbeg[orgroot]; a != EAT_LAST; a = g->oeat[a] )
      {
         i = g->head[a];
         if( g->mark[i] )
         {
            assert(Is_pseudoTerm(g->term[i]) && connected[i] == CONNECT);
            break;
         }
      }

      /* trivial solution? */
      if( a == EAT_LAST )
      {
         for( int k = 0; k < nnodes; k++ )
            mst[k].edge = UNKNOWN;

         printf("trivial solution in buildPcMwTree \n");
         for( a = g->outbeg[orgroot]; a != EAT_LAST; a = g->oeat[a] )
         {
            const int head = g->head[a];
            if( Is_term(g->term[head]) )
            {
               obj += cost[a];
               mst[head].edge = a;
            }
         }
         (*objresult) = obj;
         return SCIP_OKAY;
      }

      assert(g->mark[i]);
      mstroot = i;
   }
   assert(mstroot >= 0);
   assert(mstroot < nnodes);

   graph_path_exec(scip, g, MST_MODE, mstroot, cost, mst);

   assert(g->extended);

   /* connect all potential terminals */
   for( int i = nnodes - 1; i >= 0; --i )
   {
      if( graph_pc_knotIsDummyTerm(g, i) && i != orgroot )
      {
         int k1;
         int k2;
         const int e1 = g->inpbeg[i];
         const int e2 = g->ieat[e1];

         assert(e1 >= 0);
         assert(e2 != EAT_LAST);

#if 0
         if( e2 == EAT_LAST )
         {
            mst[i].edge = e1;
         }
         else
         {
#endif

         k1 = g->tail[e1];
         k2 = g->tail[e2];

         assert(e2 >= 0);
         assert(g->ieat[e2] == EAT_LAST);
         assert(k1 == orgroot || k2 == orgroot);

         if( k1 != orgroot && g->path_state[k1] == CONNECT )
         {
            mst[i].edge = e1;
         }
         else if( k2 != orgroot && g->path_state[k2] == CONNECT )
         {
            mst[i].edge = e2;
         }
         else if( k1 == orgroot )
         {
            mst[i].edge = e1;
         }
         else if( k2 == orgroot )
         {
            mst[i].edge = e2;
         }
         else
         {
            assert(0 && "should not happen");
         }
      }
   }    /* connect all potential terminals */

   if( !isrooted )
   {
      int e1;
      for( e1 = g->inpbeg[mstroot]; e1 != EAT_LAST; e1 = g->ieat[e1] )
         if( g->tail[e1] == orgroot )
            break;
      assert(e1 != EAT_LAST);
      mst[mstroot].edge = e1;
   }

   /* prune */
   do
   {
      int j;
      count = 0;

      for( int i = nnodes - 1; i >= 0; --i )
      {
         if( !g->mark[i] || g->path_state[i] != CONNECT || Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
         {
            const int e1 = mst[g->head[j]].edge;
            if( e1 == j )
               break;
         }

         if( j == EAT_LAST )
         {
            assert(!Is_pseudoTerm(g->term[i]));
            mst[i].edge = UNKNOWN;
            g->mark[i] = FALSE;
            connected[i] = UNKNOWN;
            count++;
            break;
         }
      }
   }
   while( count > 0 );

   if( isrooted )
   {
      /* check for feasibility */
      for( int i = 0; i < nnodes; i++ )
         if( graph_pc_knotIsFixedTerm(g, i) && i != orgroot && mst[i].edge == UNKNOWN )
         {
            *objresult = FARAWAY;
            return SCIP_OKAY;
         }
   }

   for( int i = nnodes - 1; i >= 0; --i )
      if( mst[i].edge >= 0 )
         obj += cost[mst[i].edge];

   *objresult = obj;
   return SCIP_OKAY;
}


/** prune a degree constrained Steiner tree in such a way that all leaves are terminals */
SCIP_RETCODE SCIPStpHeurTMBuildTreeDc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes (to be set) */
   )
{
   int* queue;
   int i;
   int j;
   int qsize;
   int count;
   int nnodes;
   assert(scip != NULL);
   assert(g != NULL);
   assert(result != NULL);
   assert(connected != NULL);

   nnodes = g->knots;

   /*
    * DFS until all terminals are reached
    */

   for( i = 0; i < nnodes; i++ )
      connected[i] = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &queue, nnodes) );

   qsize = 0;
   queue[qsize++] = g->source;
   connected[g->source] = TRUE;

   while( qsize )
   {
      const int node = queue[--qsize];
      for( int e = g->outbeg[node]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( (result[e] == CONNECT || result[flipedge(e)] == CONNECT) && !(connected[g->head[e]]) )
         {
            i = g->head[e];
            result[e] = CONNECT;
            result[flipedge(e)] = UNKNOWN;
            connected[i] = TRUE;
            queue[qsize++] = i;
         }
      }
   }

   SCIPfreeBufferArray(scip, &queue);

   /* prune */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for( i = 0; i < nnodes; i++ )
      {
         if( !connected[i] )
            continue;

         if( Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == CONNECT )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == CONNECT )
               {
                  result[j]    = UNKNOWN;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
            assert(j != EAT_LAST);
         }
      }
   }
   while( count > 0 );

   return SCIP_OKAY;
}


/** execute shortest paths heuristic to obtain a Steiner tree */
SCIP_RETCODE SCIPStpHeurTMRun(
   SCIP*                 scip,               /**< SCIP data structure */
   enum PCMW_TmMode      pcmw_tmmode,        /**< mode for PC/MW */
   GRAPH*                graph,              /**< graph data structure */
   int*                  starts,             /**< array containing start vertices (NULL to not provide any) */
   const SCIP_Real*      prize,              /**< prizes (for PCMW) or NULL */
   int*                  best_result,        /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   int                   runs,               /**< number of runs */
   int                   bestincstart,       /**< best incumbent start vertex */
   SCIP_Real*            cost,               /**< arc costs */
   SCIP_Real*            costrev,            /**< reversed arc costs */
   SCIP_Real*            hopfactor,          /**< edge cost multiplicator for HC problems */
   SCIP_Real*            nodepriority,       /**< vertex priorities for vertices to be starting points (NULL for no priorities) */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
   )
{
   int beststart;
   const int mode = getTmMode(getTMheurData(scip), graph);
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Bool startsgiven = (runs >= 1 && (starts != NULL));
   TMBASE tmbase = { .dheap = NULL, .cost = cost, .costrev = costrev, .nodes_dist = NULL,
                     .startnodes = NULL, .result = NULL, .best_result = best_result,
                     .nodes_pred = NULL, .connected = NULL,
                     .best_obj = FARAWAY, .nruns = runs };
   TMVNOI tmvnoi = {NULL, NULL, NULL, NULL, NULL, NULL};
   TMALLSP tmallsp = {NULL, NULL};

#ifndef NDEBUG
   assert(scip && cost && costrev && best_result);
   assert(graph->source >= 0 && nedges > 0);
   for( int e = 0; e < nedges; e++) assert(SCIPisGE(scip, cost[e], 0.0) && SCIPisGE(scip, costrev[e], 0.0));
#endif

   beststart = bestincstart;
   (*success) = FALSE;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   if( graph_pc_isPcMw(graph) )
      graph_pc_2transcheck(scip, graph);

   SCIP_CALL( tmBaseInit(scip, graph, &tmbase) );

   if( mode == TM_VORONOI )
   {
      SCIP_CALL( tmVnoiInit(scip, graph, &tmvnoi) );
   }
   else if( mode == TM_SP )
   {
      SCIP_CALL( tmAllspInit(scip, graph, &tmallsp) );
   }

   SCIP_CALL( computeStarts(scip, graph, starts, startsgiven, nodepriority, &tmbase, &beststart) );

   /* call the main routines for SPH computations, differentiate between STP variants */
   if( graph_pc_isPcMw(graph) )
   {
      SCIP_CALL( runTmPcMW(scip, graph, cost, prize, pcmw_tmmode, beststart, nodepriority, &tmbase, success));
   }
   else if( graph->stp_type == STP_DHCSTP )
   {
      SCIP_CALL( runTmDhcstp(scip, graph, &tmbase, &tmallsp, &tmvnoi, hopfactor, success) );
   }
   else
   {
      SCIP_CALL( runTm(scip, graph, &tmbase, &tmallsp, &tmvnoi, success) );
   }

   if( mode == TM_SP )
   {
      tmAllspFree(scip, graph, &tmallsp);
   }
   else if( mode == TM_VORONOI )
   {
      tmVnoiFree(scip, graph, &tmvnoi);
   }

   tmBaseFree(scip, graph, &tmbase);

   SCIPdebugMessage("final objective: %f \n", tmbase.best_obj);

  // printf("final objective: %f \n", tmbase.best_obj);

  // assert(0);


#if 0
   printf("final objective: %f \n", tmbase.best_obj);

   FILE *fp;
              fp = fopen("/nfs/optimi/kombadon/bzfrehfe/projects/scip/applications/STP/tm_csr.txt", "a+");
              fprintf(fp, "%s %f \n", SCIPgetProbName(scip), tmbase.best_obj);
              fclose(fp);
              exit(1);
#endif

   return SCIP_OKAY;
}

/** run shortest path heuristic, but bias edge costs towards best current LP solution */
SCIP_RETCODE SCIPStpHeurTMRunLP(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_HEUR*            heur,               /**< heuristic or NULL */
   int*                  result,             /**< (uninitialized) array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   int                   runs,               /**< number of runs */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
   )
{
   SCIP_VAR** vars = NULL;
   SCIP_HEURDATA* heurdata = NULL;
   SCIP_Real* xval = NULL;
   SCIP_Real* nodepriority = NULL;
   SCIP_Real* prize = NULL;
   SCIP_Real* cost = NULL;
   SCIP_Real* costrev = NULL;
   SCIP_Real randupper;
   SCIP_Real randlower;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);

   assert(scip != NULL);
   assert(result != NULL);
   assert(success != NULL);

   assert(SCIPfindHeur(scip, "TM") != NULL);
   heurdata = SCIPheurGetData(SCIPfindHeur(scip, "TM"));

   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);

   randupper = SCIPrandomGetReal(heurdata->randnumgen, 1.1, 2.5);
   randlower = SCIPrandomGetReal(heurdata->randnumgen, 1.1, randupper);

   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );

   /* LP was not solved */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      xval = NULL;
   }
   else
   {
      SCIP_SOL* sol = NULL;
      SCIP_CALL(SCIPcreateSol(scip, &sol, heur));

      /* copy the current LP solution to the working solution */
      SCIP_CALL(SCIPlinkLPSol(scip, sol));

      xval = SCIPprobdataGetXval(scip, sol);

      SCIP_CALL(SCIPfreeSol(scip, &sol));
   }

#if 0
   if( graph_pc_isPcMw(graph) )
      SCIP_CALL(SCIPallocBufferArray(scip, &prize, nnodes));
#endif

   if( xval == NULL )
   {
      BMScopyMemoryArray(cost, graph->cost, nedges);

      /* hop constraint problem? */
      if( graph->stp_type == STP_DHCSTP )
      {
         for( int e = 0; e < nedges; e++ )
         {
            if( SCIPvarGetUbGlobal(vars[e]) < 0.5 && SCIPvarGetUbGlobal(vars[flipedge_Uint(e)]) < 0.5 )
               cost[e] = BLOCKED;
         }
      }
      else
      {
         if( !graph_pc_isPcMw(graph) )
         {
            SCIP_CALL(SCIPallocBufferArray(scip, &nodepriority, nnodes));

            for( int k = 0; k < nnodes; k++ )
            {
               if( Is_term(graph->term[k]) )
                  nodepriority[k] = (SCIP_Real) graph->grad[k];
               else
                  nodepriority[k] = SCIPrandomGetReal(heurdata->randnumgen, 0.0, 1.0);
            }
         }

         for( int e = 0; e < nedges; e++ )
            if( SCIPvarGetUbGlobal(vars[e]) < 0.5 && SCIPvarGetUbGlobal(vars[flipedge_Uint(e)]) < 0.5 )
               cost[e] = MAX(cost[e], BLOCKED);
      }
   }
   else
   {
      SCIP_CALL(SCIPallocBufferArray(scip, &nodepriority, nnodes));
      initCostsAndPrioLP(scip, heurdata, vars, graph, randupper, randlower, xval, nodepriority, prize, cost);
   } /* xval != NULL */


   for( int e = 0; e < nedges; e++ )
   {
      const SCIP_Real eps = SCIPepsilon(scip);
      const SCIP_Real double_eps = 2.0 * eps;

      assert(cost[e] >= -eps);

      if( cost[e] <= eps )
      {
         assert(SCIPisZero(scip, cost[e]));

         cost[e] = double_eps;
      }

      assert(!SCIPisZero(scip, cost[e]));
   }

   for( int e = 0; e < nedges; e++ )
      costrev[e] = cost[flipedge_Uint(e)];

   /* set (edge) result array to default */
   for( int e = 0; e < nedges; e++ )
      result[e] = UNKNOWN;

   /* build a Steiner tree */
   SCIP_CALL( SCIPStpHeurTMRun(scip, pcmode_bias,
      graph, NULL, prize, result, runs, heurdata->beststartnode, cost, costrev, &(heurdata->hopfactor), nodepriority, success));

   SCIPfreeBufferArrayNull(scip, &prize);
   SCIPfreeBufferArrayNull(scip, &nodepriority);

   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &cost);

   return SCIP_OKAY;
}


/** creates the TM primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurTM(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;
   char paramdesc[SCIP_MAXSTRLEN];
   int pcmwmode = -1;

   /* create TM primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heur = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTM, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTM) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTM) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitTM) );

   heurdata->ncalls = 0;
   heurdata->nlpiterations = -1;
   heurdata->nexecs = 0;
   heurdata->randseed = DEFAULT_RANDSEED;

   /* add TM primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/evalruns",
         "number of runs for eval",
         &heurdata->evalruns, FALSE, DEFAULT_EVALRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/randseed",
         "random seed for heuristic",
         NULL, FALSE, DEFAULT_RANDSEED, 1, INT_MAX, paramChgdRandomseed, (SCIP_PARAMDATA*)heurdata) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/initruns",
         "number of runs for init",
         &heurdata->initruns, FALSE, DEFAULT_INITRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/leafruns",
         "number of runs for leaf",
         &heurdata->leafruns, FALSE, DEFAULT_LEAFRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/rootruns",
         "number of runs for root",
         &heurdata->rootruns, FALSE, DEFAULT_ROOTRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/duringlpfreq",
         "frequency for calling heuristic during LP loop",
         &heurdata->duringlpfreq, FALSE, DEFAULT_DURINGLPFREQ, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/type",
         "Heuristic: 0 automatic, 1 TM_SP, 2 TM_VORONOI, 3 TM_DIJKSTRA",
         &heurdata->type, FALSE, DEFAULT_TYPE, 0, 3, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/pcmwmode",
         "PC/MW solving mode: 0 simple, 1 bias, 2 full bias, 3 full tree, 4 all",
         &pcmwmode, FALSE, DEFAULT_PCMODE, 0, 4, NULL, NULL) );

   assert(pcmwmode >= 0 && pcmwmode <= 4);

   heurdata->pcmw_mode = pcmwmode;
   heurdata->hopfactor = DEFAULT_HOPFACTOR;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, heurdata->randseed, TRUE) );

   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing when heuristic should be called (%u:BEFORENODE, %u:DURINGLPLOOP, %u:AFTERLPLOOP, %u:AFTERNODE)", SCIP_HEURTIMING_BEFORENODE, SCIP_HEURTIMING_DURINGLPLOOP, SCIP_HEURTIMING_AFTERLPLOOP, SCIP_HEURTIMING_AFTERNODE);
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/timing", paramdesc,
         (int*) &heurdata->timing, TRUE, (int) HEUR_TIMING, (int) SCIP_HEURTIMING_BEFORENODE, 2 * (int) SCIP_HEURTIMING_AFTERPSEUDONODE - 1, NULL, NULL) ); /*lint !e713*/

   return SCIP_OKAY;
}
