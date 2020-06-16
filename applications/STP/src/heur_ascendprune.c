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

/**@file   heur_ascendprune.c
 * @brief  reduction-based primal heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a reduction and dual-cost based heuristic for Steiner problems. See
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" (2016) by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 *
 * A list of all interface methods can be found in heur_ascendprune.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define DEBUG_ASCENDPRUNE
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "heur_ascendprune.h"
#include "heur_local.h"
#include "heur_prune.h"
#include "graph.h"
#include "reduce.h"
#include "heur_tm.h"
#include "solstp.h"
#include "prop_stp.h"
#include "cons_stp.h"
#include "probdata_stp.h"

#define HEUR_NAME             "ascendprune"
#define HEUR_DESC             "Dual-cost reduction heuristic for Steiner problems"
#define HEUR_DISPCHAR         'A'
#define HEUR_PRIORITY         2
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE           /**< does the heuristic use a secondary SCIP instance?                                 */

#define DEFAULT_MAXFREQPRUNE     FALSE         /**< executions of the heuristic at maximum frequency?                               */
#define ASCENPRUNE_MINLPIMPROVE     0.2          /**< minimum percentual improvement of dual bound (wrt to gap) mandatory to execute heuristic */

#ifdef WITH_UG
int getUgRank(void);
#endif

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Real             lastdualbound;      /**< dual bound after the previous run                                 */
   int                   bestsolindex;       /**< best solution during the previous run                             */
   int                   nfailures;          /**< number of failures since last successful call                     */
   SCIP_Bool             maxfreq;            /**< should the heuristic be called at maximum frequency?              */
};


/** subgraph data */
typedef struct redcost0_graph
{
   const SCIP_Real*      redcosts;           /**< the reduced costs */
   GRAPH*                newgraph;           /**< graph                                 */
   int*                  edgelist;           /**< edge list for new graph                            */
   int*                  nodeOrg2NewMap;
   int*                  edgeNew2OrgMap;
   int                   nnodes;             /**< nodes                    */
   int                   nedges_half;        /**< edges                    */
   int                   root;               /**< red cost root             */
} RCGRAPH;


/*
 * Local methods
 */


#ifdef DEBUG_ASCENDPRUNE
/** debug method */
static
SCIP_RETCODE checkRedCostGraph(
   const GRAPH*          g,                  /**< the graph */
   const RCGRAPH*        redcostgraph        /**< the redcost graph */
)
{
   GRAPH* redgraph = redcostgraph->newgraph;
   const SCIP_Bool pcmw = graph_pc_isPcMw(g);
   const int nnodes = graph_get_nNodes(g);
   const int *nodechild = redcostgraph->nodeOrg2NewMap;

   if( graph_pc_isRootedPcMw(g) )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( graph_pc_knotIsFixedTerm(g, k) )
         {
            if( nodechild[k] < 0 || !graph_pc_knotIsFixedTerm(redgraph, nodechild[k]) )
            {
               printf("RPCMW child FAIL in AP \n\n\n");
               return SCIP_ERROR;
            }
         }
      }
   }

   for( int k = 0; k < nnodes && !pcmw; k++ )
   {
      if( Is_term(g->term[k]) )
      {
         const int i = nodechild[k];
         if( i < 0 )
         {
            printf("FAIL in AP \n\n\n");
            return SCIP_ERROR;
         }

         if( redgraph->grad[i] == 0 && redgraph->knots > 1 )
         {
            printf("FAIL GRAD \n\n\n");
            return SCIP_ERROR;
         }
      }
   }

   return SCIP_OKAY;
}


static
SCIP_RETCODE checkRedCostEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   const RCGRAPH*        redcostgraph        /**< the redcost graph */
)
{
   const int* const subgraphedges = redcostgraph->edgelist;
   SCIP_Bool* halfedge_mark;
   const int nedges = graph_get_nEdges(g);
   const int nedges_redocst = redcostgraph->nedges_half;

   assert(subgraphedges);

   SCIP_CALL( SCIPallocMemoryArray(scip, &halfedge_mark, nedges / 2) );
   BMSclearMemoryArray(halfedge_mark, nedges / 2);

   for( int e = 0; e < nedges_redocst; e++ )
   {
      const int halfedge = subgraphedges[e] / 2;

      if( halfedge_mark[halfedge] )
      {
         printf("checkRedCostEdges FAIL at %d \n", e);
         return SCIP_ERROR;
      }

      halfedge_mark[halfedge] = TRUE;
   }

   SCIPfreeMemoryArray(scip, &halfedge_mark);

   return SCIP_OKAY;
}
#endif



/** marks part of graph corresponding to zero cost paths from the root to all terminals */
static
SCIP_RETCODE redcostGraphMark(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   RCGRAPH*              redcostgraph
   )
{
   int* RESTRICT mark = g->mark;
   int* RESTRICT queue;
   STP_Bool* RESTRICT scanned;
   int* RESTRICT subgraphedges;
   const SCIP_Real* redcosts = redcostgraph->redcosts;
   int qsize;
   int nnewnodes = 0;
   int nnewedges = 0;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const SCIP_Bool pcmw = graph_pc_isPcMw(g);
   const int root_redcost = redcostgraph->root;

   assert(graph_knot_isInRange(g, root_redcost));

   /* perform DFS to identify 0-redcost subgraph.
    */

   SCIP_CALL( SCIPallocBufferArray(scip, &queue, nnodes ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scanned, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subgraphedges, nedges) );
   redcostgraph->edgelist = subgraphedges;

   BMSclearMemoryArray(mark, nnodes);
   BMSclearMemoryArray(scanned, nnodes);

   qsize = 0;
   mark[root_redcost] = TRUE;
   queue[qsize++] = root_redcost;
   nnewnodes++;

   /* DFS */
   while( qsize > 0 )
   {
      const int k = queue[--qsize];
      scanned[k] = TRUE;

      for( int a = g->outbeg[k]; a != EAT_LAST; a = g->oeat[a] )
      {
         const int head = g->head[a];

         if( SCIPisZero(scip, redcosts[a]) )
         {
            SCIP_Bool isRootEdge;
            if( pcmw && k == root_redcost && graph_pc_knotIsDummyTerm(g, head) )
               continue;

            /* vertex not visited yet? */
            if( !mark[head] )
            {
               mark[head] = TRUE;
               nnewnodes++;
               queue[qsize++] = head;
            }

            isRootEdge = (k == root_redcost || head == root_redcost);

            /* NOTE: we need to be careful to not mark anti-parallel edges. We also do not mark out-going edges from the root */
            if( !isRootEdge && (!scanned[head] || !SCIPisZero(scip, redcosts[flipedge(a)])) )
            {
               assert(g->tail[a] != root_redcost && g->head[a] != root_redcost);

               subgraphedges[nnewedges++] = a;
            }
         }
      }
   }

#ifndef NDEBUG
   for( int k = 0; k < nnewedges && pcmw; k++ )
   {
      const int e = subgraphedges[k];
      assert(g->tail[e] != root_redcost && g->head[e] != root_redcost);
   }
#endif

   for( int a = g->outbeg[root_redcost]; a != EAT_LAST; a = g->oeat[a] )
   {
      const int head = g->head[a];
      if( mark[head] )
      {
         subgraphedges[nnewedges++] = a;
      }
   }

   /* we need to make sure that the dummy terminals are connected from the root */
   if( pcmw && root_redcost != g->source )
   {
      assert(graph_pc_isRootedPcMw(g));
      assert(mark[g->source]);

      for( int a = g->outbeg[g->source]; a != EAT_LAST; a = g->oeat[a] )
      {
         const int head = g->head[a];
         if( mark[head] && graph_pc_knotIsDummyTerm(g, head) && !SCIPisZero(scip, redcosts[a]) )
         {
            subgraphedges[nnewedges++] = a;
         }
      }
   }

   SCIPfreeBufferArray(scip, &scanned);
   SCIPfreeBufferArray(scip, &queue);

   redcostgraph->nnodes = nnewnodes;
   redcostgraph->nedges_half = nnewedges;

#ifdef DEBUG_ASCENDPRUNE
   SCIP_CALL( checkRedCostEdges(scip, g, redcostgraph) );
#endif

   return SCIP_OKAY;
}


/** gets number of terms that are marked */
static
int redcostGetNTermsMarked(
   const GRAPH*          g,                  /**< the graph */
   const RCGRAPH*        redcostgraph        /**< reduced cost graph data */
   )
{
   int nterms = 0;
   const int nnodes = graph_get_nNodes(g);

   for( int i = 0; i < nnodes; i++ )
   {
      if( g->mark[i] && Is_term(g->term[i]) )
         nterms++;
   }

   return nterms;
}


/** builds graphs */
static
SCIP_RETCODE redcostGraphBuild(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   RCGRAPH*              redcostgraph
   )
{
   GRAPH* newgraph;
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool pcmw = graph_pc_isPcMw(g);
   const int probtype = g->stp_type;
   const int* const mark = g->mark;
   const int* const newedges = redcostgraph->edgelist;
   const int nnewnodes = redcostgraph->nnodes;
   const int nnewedges = redcostgraph->nedges_half;
   int* edgeNew2OrgMap;
   int* nodeOrg2NewMap;

   assert(newedges);
   assert(nnewnodes > 0 && nnewedges > 0);
   assert(graph_knot_isInRange(g, redcostgraph->root));

   SCIP_CALL( SCIPallocBufferArray(scip, &nodeOrg2NewMap, nnodes ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgeNew2OrgMap, 2 * nnewedges) );
   redcostgraph->edgeNew2OrgMap = edgeNew2OrgMap;
   redcostgraph->nodeOrg2NewMap = nodeOrg2NewMap;

   /* initialize new graph */
   SCIP_CALL( graph_init(scip, &(redcostgraph->newgraph), nnewnodes, 2 * nnewedges, 1) );
   newgraph = redcostgraph->newgraph;

   if( probtype == STP_RSMT || probtype == STP_OARSMT || probtype == STP_GSTP )
      newgraph->stp_type = STP_SPG;
   else
      newgraph->stp_type = probtype;

   if( pcmw )
      SCIP_CALL( graph_pc_initSubgraph(scip, newgraph) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( mark[k] )
      {
         if( pcmw )
         {
            assert(g->extended);
            newgraph->prize[newgraph->knots] = g->prize[k];

            assert(!graph_pc_knotIsDummyTerm(g, k) || 0.0 == newgraph->prize[newgraph->knots]);
         }

         nodeOrg2NewMap[k] = newgraph->knots;
         graph_knot_add(newgraph, g->term[k]);
      }
      else
      {
         nodeOrg2NewMap[k] = UNKNOWN;
      }
   }

   if( pcmw )
   {
      newgraph->norgmodelknots = nnewnodes;
      newgraph->extended = TRUE;
   }

   assert(nnewnodes == newgraph->knots);

   /* set root of new graph */
   if( pcmw )
      newgraph->source = nodeOrg2NewMap[g->source];
   else
      newgraph->source = nodeOrg2NewMap[redcostgraph->root];

   assert(newgraph->source >= 0);
   assert(!graph_pc_isRootedPcMw(g) || newgraph->prize[newgraph->source] == FARAWAY);

   /* add edges to new graph */
   for( int a = 0; a < nnewedges; a++ )
   {
      int i;
      const int e = newedges[a];
      const int tail = nodeOrg2NewMap[g->tail[e]];
      const int head = nodeOrg2NewMap[g->head[e]];

      assert(tail >= 0);
      assert(head >= 0);

      for( i = newgraph->outbeg[tail]; i != EAT_LAST; i = newgraph->oeat[i] )
      {
         if( newgraph->head[i] == head )
            break;
      }

      /* edge not added yet? */
      if( i == EAT_LAST )
      {
         edgeNew2OrgMap[newgraph->edges] = e;
         edgeNew2OrgMap[newgraph->edges + 1] = flipedge(e);

         graph_edge_addSubgraph(scip, g, nodeOrg2NewMap, e, newgraph);
      }
   }

   newgraph->norgmodeledges = newgraph->edges;

   assert(!pcmw || TERM2EDGE_FIXEDTERM == newgraph->term2edge[newgraph->source]);

   SCIP_CALL( graph_pc_finalizeSubgraph(scip, newgraph) );

   /* initialize ancestors of new graph edges */
   SCIP_CALL( graph_init_history(scip, newgraph) );

   return SCIP_OKAY;
}


/** frees members */
static
void redcostGraphFree(
   SCIP*                 scip,               /**< SCIP data structure */
   RCGRAPH*              redcostgraph
   )
{
   if( redcostgraph->newgraph )
   {
      graph_free(scip, &(redcostgraph->newgraph), TRUE);
   }
   SCIPfreeBufferArrayNull(scip, &(redcostgraph->edgeNew2OrgMap));
   SCIPfreeBufferArrayNull(scip, &(redcostgraph->nodeOrg2NewMap));
   SCIPfreeBufferArrayNull(scip, &(redcostgraph->edgelist));
}



/** builds Steiner tree on subgraph */
static
SCIP_RETCODE redcostGraphComputeSteinerTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   RCGRAPH*              redcostgraph,
   int*                  result
   )
{
   GRAPH* newgraph = redcostgraph->newgraph;
   STP_Bool* RESTRICT nodearrchar;
   const int* const edgeNew2OrgMap = redcostgraph->edgeNew2OrgMap;
   const int nnodes = graph_get_nNodes(g);
   const int nnewedges = graph_get_nEdges(newgraph);
   SCIP_Bool success;

   SCIP_CALL( graph_path_init(scip, newgraph) );
   SCIP_CALL( reduceLevel0(scip, newgraph) );
   assert(graph_valid(scip, newgraph));

#ifdef DEBUG_ASCENDPRUNE
   SCIP_CALL( checkRedCostGraph(g, redcostgraph) );
#endif

   /* get solution on new graph by PRUNE heuristic */
   SCIP_CALL( SCIPStpHeurPruneRun(scip, NULL, newgraph, result, &success, FALSE, TRUE) );

#ifndef NDEBUG
   for( int k = 0; k < newgraph->knots; ++k )
   {
      assert( !(Is_term(newgraph->term[k]) && newgraph->grad[k] == 0 && k != newgraph->source) );
   }
   assert(solstp_isValid(scip, newgraph, result));
#endif

   assert(success && solstp_isValid(scip, newgraph, result));

   SCIPdebugMessage("obj after prune %f \n", solstp_getObjBounded(newgraph, result, 0.0, newgraph->edges));

   SCIP_CALL( SCIPStpHeurLocalRun(scip, newgraph, result) );

   SCIPdebugMessage("obj after local %f \n", solstp_getObjBounded(newgraph, result, 0.0, newgraph->edges));

   assert(solstp_isValid(scip, newgraph, result));
   graph_path_exit(scip, newgraph);

    /*
    * prune solution (in the original graph)
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes ) );
   BMSclearMemoryArray(nodearrchar, nnodes);

   for( int e = 0; e < nnewedges; e++ )
   {
      if( result[e] == CONNECT )
      {
         const int eorg = edgeNew2OrgMap[e];
         assert(graph_edge_isInRange(g, eorg));

         nodearrchar[g->tail[eorg]] = TRUE;
         nodearrchar[g->head[eorg]] = TRUE;
      }
   }

   if( newgraph->knots == 1 )
      nodearrchar[g->source] = TRUE;

   SCIP_CALL( solstp_prune(scip, g, result, nodearrchar) );

   assert(solstp_isValid(scip, g, result));
   SCIPfreeBufferArray(scip, &nodearrchar);

   return SCIP_OKAY;
}



/*
 * Callback methods of primal heuristic
 */


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyAscendPrune)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPStpIncludeHeurAscendPrune(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeAscendPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */

static
SCIP_DECL_HEURINIT(heurInitAscendPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   /* initialize data */
   heurdata->nfailures = 0;
   heurdata->bestsolindex = -1;
   heurdata->lastdualbound = 0.0;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecAscendPrune)
{  /*lint --e{715}*/
   SCIP_HEURDATA*    heurdata;
   SCIP_PROBDATA*    probdata;
   SCIP_VAR**        vars;
   SCIP_SOL*         bestsol;                        /* best known solution */
   GRAPH*            graph;
   SCIP_Real         dualbound;
   SCIP_Real         gap;
   SCIP_Real*        redcosts;
   SCIP_Bool         solAdded;
   int       nedges;
   int*      edgearrint;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);

   *result = SCIP_DIDNOTRUN;

   if( !SCIPStpRedcostAvailable(scip) )
      return SCIP_OKAY;

   nedges = graph->edges;
   solAdded = FALSE;

   bestsol = SCIPgetBestSol(scip);

   /* no solution available? */
   if( bestsol == NULL )
      return SCIP_OKAY;

   /* get dual bound */
   dualbound = SCIPgetDualbound(scip);

   /* no new best solution available? */
   if( heurdata->bestsolindex == SCIPsolGetIndex(SCIPgetBestSol(scip)) && !(heurdata->maxfreq) )
   {
      /* current optimality gap */
      gap = SCIPgetSolOrigObj(scip, bestsol) - dualbound;

      if( SCIPisLT(scip, dualbound - heurdata->lastdualbound, gap * ASCENPRUNE_MINLPIMPROVE ) )
         return SCIP_OKAY;
   }

   printf("CALL ascend-prune \n");

   heurdata->lastdualbound = dualbound;

   /* allocate memory for ascend-and-prune */
   SCIP_CALL( SCIPallocBufferArray(scip, &redcosts, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );

   SCIPStpGetRedcosts(scip, vars, graph, redcosts);

   /* perform ascent and prune */
   SCIP_CALL( SCIPStpHeurAscendPruneRun(scip, heur, graph, redcosts, edgearrint, graph->source, &solAdded, TRUE) );

   if( solAdded )
   {
      heurdata->nfailures = 0;
      *result = SCIP_FOUNDSOL;
   }
   else
   {
      heurdata->nfailures++;
   }

   heurdata->bestsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   SCIPfreeBufferArray(scip, &edgearrint);
   SCIPfreeBufferArray(scip, &redcosts);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */


/** ascent and prune */
SCIP_RETCODE SCIPStpHeurAscendPruneRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure or NULL */
   const GRAPH*          g,                  /**< the graph */
   const SCIP_Real*      redcosts,           /**< the reduced costs */
   int*                  result,             /**< int edges array to store solution */
   int                   root,               /**< the root (used for dual ascent) */
   SCIP_Bool*            solfound,           /**< has a solution been found? And added, if requested? */
   SCIP_Bool             addsol              /**< should the solution be added to SCIP by this method? */
   )
{
   RCGRAPH redcostgraph = { .newgraph = NULL, .redcosts = redcosts, .edgelist = NULL, .nodeOrg2NewMap = NULL, .edgeNew2OrgMap = NULL,
         .nnodes = -1, . nedges_half = -1, .root = root  };

   assert(scip && redcosts && result && solfound);
   assert(!graph_pc_isPcMw(g) || g->extended);

   *solfound = TRUE;

   if( root < 0 )
   {
      redcostgraph.root = g->source;
   }

   assert(Is_term(g->term[root]));
   assert(graph_valid(scip, g));
   assert(!graph_pc_isPcMw(g) || graph_pc_knotIsFixedTerm(g, root));

   SCIP_CALL( redcostGraphMark(scip, g, &redcostgraph) );

   if( redcostgraph.nedges_half == 0 )
   {
      assert(g->stp_type == STP_RMWCSP);
      solstp_getTrivialSol(g, result);
      redcostGraphFree(scip, &redcostgraph);

      return SCIP_OKAY;
   }

   if( graph_typeIsSpgLike(g) && g->terms != redcostGetNTermsMarked(g, &redcostgraph) )
   {
      printf("not all terminals could be reached (%d < %d) abort ascend-and-prune \n",  g->terms,  redcostGetNTermsMarked(g, &redcostgraph) );
      redcostGraphFree(scip, &redcostgraph);
      *solfound = FALSE;

      return SCIP_OKAY;
   }

   assert(graph_pc_isPcMw(g) || redcostGetNTermsMarked(g, &redcostgraph) == g->terms);

   SCIP_CALL( redcostGraphBuild(scip, g, &redcostgraph) );
   SCIP_CALL( redcostGraphComputeSteinerTree(scip, g, &redcostgraph, result) );

   assert(solstp_isValid(scip, g, result));

#ifdef SCIP_DEBUG
   {
      const SCIP_Real solval = solstp_getObj(g, result, 0.0);
      graph_printInfo(g);
      printf("ascend-prune obj=%f \n", solval);
   }
#endif

   if( addsol )
   {
      const int nedges = graph_get_nEdges(g);
      const int nvars = SCIPprobdataGetNVars(scip);
      SCIP_Real* nval;
      SCIP_CALL( SCIPallocBufferArray(scip, &nval, nvars) );

      for( int e = 0; e < nedges; e++ )
      {
         if( result[e] == CONNECT )
            nval[e] = 1.0;
         else
            nval[e] = 0.0;
      }

      SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, heur, solfound) );
      SCIPdebugMessage("Ascend-and-prune added solution \n");

      SCIPfreeBufferArray(scip, &nval);
   }

   for( int k = 0; k < g->knots; k++ )
      g->mark[k] = (g->grad[k] > 0);

   redcostGraphFree(scip, &redcostgraph);

   return SCIP_OKAY;
}


/** creates the prune primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurAscendPrune(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create prune primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecAscendPrune, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyAscendPrune) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeAscendPrune) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitAscendPrune) );

   /* add ascend and prune primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/maxfreq",
         "should the heuristic be executed at maximum frequeny?",
         &heurdata->maxfreq, FALSE, DEFAULT_MAXFREQPRUNE, NULL, NULL) );

   return SCIP_OKAY;
}
