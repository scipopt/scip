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

/**@file   solstp.c
 * @brief  includes methods working on solutions (i.e. trees) to Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * Methods for manipulating solutions (i.e. trees) to Steiner tree problems, such as pruning.
 * Also includes methods for obtaining information about solutions.
 *
 * A list of all interface methods can be found in solstp.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */
//#define SCIP_DEBUG

#include "probdata_stp.h"
#include "portab.h"
#include "solstp.h"
#include "mst.h"
#include "shortestpath.h"

/** Deletes subtree from given node, marked by dfspos.
 *  NOTE: recursive method. */
static
void pcsubtreeDelete(
   const GRAPH*          g,                  /**< graph structure */
   int                   subtree_root,       /**< root of the subtree */
   int                   dfspos[],           /**< array to mark DFS positions of nodes */
   int                   result[],           /**< ST edges */
   STP_Bool              connected[]         /**< ST nodes */
)
{
   const int dfspos_root = dfspos[subtree_root];

   assert(dfspos_root > 0);
   assert(connected[subtree_root]);
   assert(g->mark[subtree_root]);

   connected[subtree_root] = FALSE;
   g->mark[subtree_root] = FALSE;

   SCIPdebugMessage("strong prune deletes tree vertex %d \n", subtree_root);

   for( int e = g->outbeg[subtree_root]; e != EAT_LAST; e = g->oeat[e] )
   {
      if( result[e] == CONNECT )
      {
         const int neighbor = g->head[e];

         assert(dfspos[neighbor] >= 0);
         assert(!graph_pc_knotIsDummyTerm(g, neighbor));
         assert(dfspos[neighbor] != dfspos_root);

         /* is neighbor a DFS child of the root?  */
         if( dfspos[neighbor] > dfspos_root)
         {
            result[e] = UNKNOWN;
#ifdef SCIP_DEBUG
            SCIPdebugMessage("strong prune deletes tree edge ");
            graph_edge_printInfo(g, e);
#endif
            pcsubtreeDelete(g, neighbor, dfspos, result, connected);
         }
      }
   }
}


/** Deletes subtree from given node, marked by dfspos.
 *  NOTE: recursive method. */
static
void pcsubtreeDelete_csr(
   const CSR*            csr_orgcosts,       /**< CSR */
   int                   subtree_root,       /**< root of the subtree */
   int* RESTRICT         dfspos,             /**< array to mark DFS positions of nodes */
   int* RESTRICT         result,             /**< ST edges */
   STP_Bool* RESTRICT    connected           /**< ST nodes */
)
{
   const int dfspos_root = dfspos[subtree_root];
   const int rootedges_start = csr_orgcosts->start[subtree_root];
   const int rootedges_end = csr_orgcosts->start[subtree_root + 1];

   assert(dfspos_root > 0);
   assert(connected[subtree_root]);

   connected[subtree_root] = FALSE;

   SCIPdebugMessage("strong prune deletes tree vertex %d \n", subtree_root);

   for( int e = rootedges_start; e != rootedges_end; e++ )
   {
      if( result[e] == CONNECT )
      {
         const int neighbor = csr_orgcosts->head[e];

         assert(dfspos[neighbor] >= 0);
         assert(dfspos[neighbor] != dfspos_root);

         /* is neighbor a DFS child of the root?  */
         if( dfspos[neighbor] > dfspos_root)
         {
            result[e] = UNKNOWN;
            pcsubtreeDelete_csr(csr_orgcosts, neighbor, dfspos, result, connected);
         }
      }
   }
}


/** Prunes subtree from given node such that it becomes most profitable and returns the profit.
 *  NOTE: recursive method. */
static
SCIP_Real pcsubtreePruneForProfit(
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int                   subtree_root,       /**< root of the subtree */
   int* RESTRICT         dfspos,             /**< array to mark DFS positions of nodes */
   int* RESTRICT         result,             /**< ST edges */
   STP_Bool* RESTRICT    connected,        /**< ST nodes */
   int* RESTRICT         dfscount            /**< counter */
)
{
   SCIP_Real profit = g->prize[subtree_root];

   if( !graph_pc_isPc(g) )
   {
      assert(graph_pc_isMw(g));

      if( LT(profit, 0.0) )
         profit = 0.0;
   }

   assert(0 <= *dfscount && *dfscount < g->knots);

   dfspos[subtree_root] = ++(*dfscount);

   SCIPdebugMessage("strong-prune from root %d \n", subtree_root);

   for( int e = g->outbeg[subtree_root]; e >= 0; e = g->oeat[e] )
   {
      if( result[e] == CONNECT )
      {
         const int neighbor = g->head[e];

         assert(dfspos[neighbor] >= 0);
         assert(!graph_pc_knotIsDummyTerm(g, neighbor));

         /* not visited yet? */
         if( dfspos[neighbor] == 0 )
         {
            const SCIP_Real neighbor_profit = pcsubtreePruneForProfit(g, cost, neighbor, dfspos, result, connected, dfscount);
            const SCIP_Real extension_profit = neighbor_profit - cost[e];

            if( LT(extension_profit, 0.0) )
            {
               result[e] = UNKNOWN;
#ifdef SCIP_DEBUG
               SCIPdebugMessage("strong prune deletes tree edge ");
               graph_edge_printInfo(g, e);
#endif
               pcsubtreeDelete(g, neighbor, dfspos, result, connected);
            }
            else
            {
               profit += extension_profit;
            }
         }
      }
   }

   return profit;
}



/** Prunes subtree from given node such that it becomes most profitable and returns the profit.
 *  NOTE: recursive method. */
static
SCIP_Real pcsubtreePruneForProfit_csr(
   const CSR*            csr_orgcosts,       /**< CSR */
   const SCIP_Real*      prize,              /**< the prize */
   SCIP_Bool             isPc,               /**< is PC? */
   int                   subtree_root,       /**< root of the subtree */
   int* RESTRICT         dfspos,             /**< array to mark DFS positions of nodes */
   int* RESTRICT         result,             /**< ST edges */
   STP_Bool* RESTRICT    connected,          /**< ST nodes */
   int* RESTRICT         dfscount            /**< counter */
)
{
   const SCIP_Real* const orgcosts_csr = csr_orgcosts->cost;
   const int rootedges_start = csr_orgcosts->start[subtree_root];
   const int rootedges_end = csr_orgcosts->start[subtree_root + 1];
   const int* const heads_csr = csr_orgcosts->head;
   SCIP_Real profit = prize[subtree_root];

   if( !isPc )
   {
      /* NOTE: for MW any negative prize is already counted with the edge costs! */
      if( LT(profit, 0.0) )
         profit = 0.0;
   }

   assert(0 <= *dfscount && *dfscount < csr_orgcosts->nnodes);
   assert(rootedges_start <= rootedges_end);

   dfspos[subtree_root] = ++(*dfscount);

   SCIPdebugMessage("strong-prune from root %d \n", subtree_root);

   for( int e = rootedges_start; e != rootedges_end; e++ )
   {
      if( result[e] == CONNECT )
      {
         const int neighbor = heads_csr[e];

         assert(dfspos[neighbor] >= 0);

         /* not visited yet? */
         if( dfspos[neighbor] == 0 )
         {
            const SCIP_Real neighbor_profit = pcsubtreePruneForProfit_csr(csr_orgcosts, prize, isPc, neighbor, dfspos,
                  result, connected, dfscount);
            const SCIP_Real extension_profit = neighbor_profit - orgcosts_csr[e];

            if( LE(extension_profit, 0.0) )
            {
               result[e] = UNKNOWN;
               pcsubtreeDelete_csr(csr_orgcosts, neighbor, dfspos, result, connected);
            }
            else
            {
               profit += extension_profit;
            }
         }
      }
   }

   return profit;
}


/** computes trivial solution and sets result edges */
static inline
void pcsolGetTrivialEdges(
   const GRAPH*          g,                  /**< graph structure */
   const STP_Bool*       connected,          /**< ST nodes */
   int* RESTRICT         result              /**< MST solution, which does not include artificial terminals */
)
{
   const int root = g->source;
   const int* const gOeat = g->oeat;
   const int* const gHead = g->head;

#ifndef NEDBUG
   for( int i = 0; i < g->edges; i++ )
      assert(UNKNOWN == result[i]);
#endif

   for( int a = g->outbeg[root]; a >= 0; a = gOeat[a] )
   {
      const int head = gHead[a];
      if( graph_pc_knotIsDummyTerm(g, head) )
      {
         assert(!connected || connected[head]);
         result[a] = CONNECT;
      }
   }
}


/** computes MST on marked graph and sets result edges */
static inline
SCIP_RETCODE pcsolGetMstEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int                   root,               /**< root of solution */
   int*                  result              /**< MST solution, which does not include artificial terminals */
)
{
   PATH* mst;
   const int nnodes = graph_get_nNodes(g);
   const int* const gmark = g->mark;

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );
   graph_path_exec(scip, g, MST_MODE, root, cost, mst);

   for( int i = 0; i < nnodes; i++ )
   {
      if( gmark[i] && (mst[i].edge != UNKNOWN) )
      {
         assert(g->path_state[i] == CONNECT);
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == -1);

         result[mst[i].edge] = CONNECT;
      }
   }

   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}


/** mark nodes of the solution in the graph */
static inline
void pcsolMarkGraphNodes(
   const STP_Bool*       connected,          /**< ST nodes */
   const GRAPH*          g                   /**< graph structure */
   )
{
   int* RESTRICT gmark = g->mark;
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(g);

   assert(g->extended);

   if( rpcmw )
   {
      for( int i = 0; i < nnodes; i++ )
      {
         if( connected[i] && !graph_pc_knotIsDummyTerm(g, i) )
            gmark[i] = TRUE;
         else
            gmark[i] = FALSE;

         assert(gmark[i] || !graph_pc_knotIsFixedTerm(g, i));
      }

      assert(gmark[g->source]);
   }
   else
   {
      const int* const gterm = g->term;

      for( int i = 0; i < nnodes; i++ )
      {
         if( connected[i] && !Is_term(gterm[i]) )
            gmark[i] = TRUE;
         else
            gmark[i] = FALSE;
      }
   }
}


/** mark nodes of the solution in the graph */
static inline
void pcsolMarkGraphNodes_csr(
   const GRAPH*          g,                  /**< graph structure */
   STP_Bool* RESTRICT    connected           /**< ST nodes */
   )
{
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(g);

   assert(g->extended);

   if( rpcmw )
   {
      const int* const gGrad = g->grad;

      for( int i = 0; i < nnodes; i++ )
      {
         if( gGrad[i] == 2 )
         {
            if( graph_pc_knotIsDummyTerm(g, i) )
               connected[i] = FALSE;
         }
         else
         {
            assert(!graph_pc_knotIsDummyTerm(g, i));
         }

         assert(connected[i] || !graph_pc_knotIsFixedTerm(g, i));
      }

      assert(connected[g->source]);
   }
   else
   {
      const int* const gterm = g->term;

      for( int i = 0; i < nnodes; i++ )
      {
         if( Is_term(gterm[i]) )
            connected[i] = FALSE;
      }
   }
}


/** prune a Steiner tree in such a way that all leaves are terminals */
static inline
void pcsolPrune(
   const GRAPH*          g,                  /**< graph structure */
   int*                  result,             /**< MST solution, which does not include artificial terminals */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   const int nnodes = graph_get_nNodes(g);
   int count;

   SCIPdebugMessage("starting (simple) pruning \n");

   do
   {
      count = 0;

      for( int i = nnodes - 1; i >= 0; --i )
      {
         int j;
         if( !g->mark[i] || g->path_state[i] != CONNECT || Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == CONNECT )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            assert(!Is_term(g->term[i]) && !Is_pseudoTerm(g->term[i]));

            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == CONNECT )
               {
                  SCIPdebugMessage("prune delete vertex %d \n", i);
#ifdef SCIP_DEBUG
                  SCIPdebugMessage("prune delete edge ");
                  graph_edge_printInfo(g, j);
#endif

                  result[j] = UNKNOWN;
                  g->mark[i] = FALSE;
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

#ifndef NDEBUG
   /* make sure there is no unconnected vertex */
   for( int i = 0; i < nnodes; i++ )
   {
      if( connected[i] && i != g->source )
      {
         int j;
         for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
         {
            if( result[j] == CONNECT )
               break;
         }

         assert(j != EAT_LAST);
      }
   }
#endif
}


/** prunes the Steiner tree in such a way that all leaves are terminals */
static
void stpsolPrune_csr(
   const GRAPH*          g,                  /**< graph structure */
   const MST*            mst,                /**< the MST */
   int* RESTRICT         result,             /**< ST edges (need to be set to UNKNOWN) */
   STP_Bool* RESTRICT    connected           /**< ST nodes */
   )
{
   const CSR* const csr = mst->csr;
   const int* const csr_start = csr->start;
   const int nnodes = graph_get_nNodes(g);
   int count;

   SCIPdebugMessage("starting (simple, CSR) pruning \n");

   do
   {
      count = 0;

      for( int i = 0; i < nnodes; i++ )
      {
         int outedge;

         if( !connected[i] || Is_term(g->term[i]) )
            continue;

         for( outedge = csr_start[i]; outedge != csr_start[i + 1]; outedge++ )
         {
            if( result[outedge] == CONNECT )
               break;
         }

         /* no outgoing edge? */
         if( outedge == csr_start[i + 1] )
         {
            /* there has to be exactly one incoming edge -> remove it */

            const int inedge = mst->nodes_predEdge[i];
            assert(inedge != UNKNOWN);
            assert(result[inedge] == CONNECT);

            SCIPdebugMessage("prune delete vertex %d \n", i);

            result[inedge] = UNKNOWN;
            connected[i] = FALSE;
            count++;
         }
      }
   }
   while( count > 0 );
}


/** computes MST on marked graph and sets result edges */
static inline
void solGetMstEdges_csr(
   const GRAPH*          g,                  /**< graph structure */
   const STP_Bool*       connected,          /**< ST nodes */
   int                   root,               /**< root of solution */
   MST*                  mst,                /**< the MST */
   int* RESTRICT         result              /**< MST solution, which does not include artificial terminals */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int* predEdge;

   mst_computeOnMarked(g, connected, root, mst);
   predEdge = mst->nodes_predEdge;

   for( int i = 0; i < nnodes; i++ )
   {
      if( connected[i] && (predEdge[i] != UNKNOWN) )
      {
         assert(mst->csr->head[predEdge[i]] == i);
         assert(result[predEdge[i]] == UNKNOWN);

         result[predEdge[i]] = CONNECT;
      }
   }
}

/** Finds optimal prize-collecting Steiner tree on given tree. */
static
SCIP_RETCODE strongPruneSteinerTreePc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int                   solroot,            /**< root of the solution */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
)
{
   int* dfspos;
   const int nnodes = graph_get_nNodes(g);
   int dfscount = 0;
   SCIP_Real profit;
#ifndef NDEBUG
   const int nsoledges = solstp_getNedges(g, result);
#endif

   assert(solroot >= 0);
   assert(connected[solroot]);
   assert(graph_pc_isPcMw(g));
   assert(!graph_pc_knotIsDummyTerm(g, solroot));
   assert(graph_pc_isMw(g) || graph_pc_costsEqualOrgCosts(scip, g, cost));
   assert(g->extended);

   /* todo find best root? */

   SCIP_CALL( SCIPallocBufferArray(scip, &dfspos, nnodes) );

   BMSclearMemoryArray(dfspos, nnodes);

   /* compute the subtree */
   profit = pcsubtreePruneForProfit(g, cost, solroot, dfspos, result, connected, &dfscount);

   assert(nsoledges + 1 == dfscount);

   if( LT(profit, 0.0) )
   {
      assert(!graph_pc_isRootedPcMw(g));
      assert(!Is_anyTerm(g->term[solroot]));
      assert(EQ(g->prize[solroot], 0.0));

      // todo can this ever happen?
      // if so, better have a flag here, because we dont wannt set edges to dummies here...
      return SCIP_ERROR;
//      SCIPdebugMessage("Best subtree is negative! Take empty solution \n");
//      pcsolGetTrivialEdges(g, connected, result);
   }

   SCIPfreeBufferArray(scip, &dfspos);

   return SCIP_OKAY;
}



/** Finds optimal prize-collecting Steiner tree on given tree. */
static
SCIP_RETCODE strongPruneSteinerTreePc_csr(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const CSR*            csr_orgcosts,       /**< CSR */
   int                   solroot,            /**< root of the solution */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
)
{
   int* dfspos;
   SCIP_Real profit;
   const int nnodes = graph_get_nNodes(g);
   int dfscount = 0;
#ifndef NDEBUG
   const int nsoledges = solstp_getNedgesBounded(g, result, graph_csr_getNedges(csr_orgcosts));
#endif
   const SCIP_Bool isPc = graph_pc_isPc(g);

   assert(solroot >= 0);
   assert(connected[solroot]);
   assert(graph_pc_isPcMw(g));
   assert(!graph_pc_knotIsDummyTerm(g, solroot));
   assert(g->extended);
   assert(nnodes == csr_orgcosts->nnodes);


   SCIP_CALL( SCIPallocBufferArray(scip, &dfspos, nnodes) );

   BMSclearMemoryArray(dfspos, nnodes);

   /* compute the subtree */

   profit = pcsubtreePruneForProfit_csr(csr_orgcosts, g->prize, isPc, solroot, dfspos, result, connected, &dfscount);

   assert(nsoledges + 1 == dfscount);

   SCIPfreeBufferArray(scip, &dfspos);

   if( LT(profit, 0.0) )
   {
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}



/** prune a Steiner tree in such a way that all leaves are terminals */
static
SCIP_RETCODE pruneSteinerTreeStp(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int*                  result,             /**< ST edges, which need to be set to UNKNOWN */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   PATH* mst;
   int count;
   const int nnodes = graph_get_nNodes(g);
#ifndef NEDBUG
   int nconnected = 0;
#endif

   assert(scip != NULL);
   assert(cost != NULL);
   assert(result != NULL);
   assert(connected != NULL);

#ifndef NEDBUG
   for( int i = 0; i < g->edges; i++ )
      assert(UNKNOWN == result[i]);

   for( int i = nnodes - 1; i >= 0; --i )
      if( connected[i] )
         nconnected++;

   assert(nconnected >= g->terms);
   assert(g->source >= 0 && g->source < nnodes);
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );

   /* compute the MST */
   for( int i = nnodes - 1; i >= 0; --i )
      g->mark[i] = connected[i];

   graph_path_exec(scip, g, MST_MODE, g->source, cost, mst);

   for( int i = nnodes - 1; i >= 0; --i )
   {
      if( connected[i] && (mst[i].edge != -1) )
      {
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == UNKNOWN);

         result[mst[i].edge] = 0;
      }
   }

   /* prune */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for( int i = nnodes - 1; i >= 0; --i )
      {
         int j;

         if( !g->mark[i] )
            continue;

         if( g->term[i] == 0 )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == 0 )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == 0 )
               {
                  result[j]    = -1;
                  g->mark[i]   = FALSE;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
         }
      }
   }
   while( count > 0 );

   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}


/** Prunes the Steiner tree in such a way that all leaves are terminals:
 *  1. Builds MST
 *  2. Removes non-terminal leaves repeatedly */
static
SCIP_RETCODE pruneSteinerTreeStp_csr(
   const GRAPH*          g,                  /**< graph structure */
   MST*                  mst,                /**< the MST */
   int* RESTRICT         result,             /**< ST edges (need to be set to UNKNOWN) */
   STP_Bool* RESTRICT    connected           /**< ST nodes */
   )
{
   assert(g && mst && result && connected);

#ifndef NDEBUG
   {
      const int nedges_csr = graph_csr_getNedges(mst->csr);
      for( int i = 0; i < nedges_csr; i++ )
         assert(UNKNOWN == result[i]);
   }
#endif

   /* 1. build MST on solution nodes */
   solGetMstEdges_csr(g, connected, g->source, mst, result);

   /* 2. prune MST */
   stpsolPrune_csr(g, mst, result, connected);

   return SCIP_OKAY;
}


/* prune the (rooted) prize collecting Steiner tree in such a way that all leaves are terminals
 * NOTE: graph is not really const, mark is changed! todo */
static
SCIP_RETCODE pruneSteinerTreePc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int*                  result,             /**< ST edges (need to be set to UNKNOWN) */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(g);
   int solroot = g->source;

#ifndef NDEBUG
   int* result_dbg;
   STP_Bool* connected_dbg;
   const int nedges = graph_get_nEdges(g);
   for( int i = 0; i < nedges; i++ )
      assert(UNKNOWN == result[i]);
#endif

   assert(scip && cost && result && connected);
   assert(g->extended);
   assert(graph_pc_isPcMw(g));

   pcsolMarkGraphNodes(connected, g);

   if( !rpcmw )
   {
      solroot = solstp_pcGetSolRoot(scip, g, connected);

      /* trivial solution? */
      if( solroot == -1 )
      {
         printf("trivial solution in pruning \n");

         pcsolGetTrivialEdges(g, connected, result);

         return SCIP_OKAY;
      }
   }

   assert(0 <= solroot && solroot < g->knots);
   assert(g->mark[solroot]);
   SCIPdebugMessage("(non-artificial) solution root=%d \n", solroot);

   SCIP_CALL( pcsolGetMstEdges(scip, g, cost, solroot, result) );

#ifndef NDEBUG
   for( int i = 0; i < g->knots; ++i )
      assert((g->path_state[i] == CONNECT) == g->mark[i]);

   SCIP_CALL( SCIPallocBufferArray(scip, &result_dbg, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &connected_dbg, g->knots) );

   BMScopyMemoryArray(result_dbg, result, nedges);
   BMScopyMemoryArray(connected_dbg, connected, g->knots);
#endif

   SCIP_CALL( strongPruneSteinerTreePc(scip, g, cost, solroot, result, connected) );

   solstp_pcConnectDummies(g, solroot, result, connected);

   /* simple pruning */
   pcsolPrune(g, result, connected);
   assert(!stpsol_pruningIsPossible(g, result, connected));

#ifndef NDEBUG
   solstp_pcConnectDummies(g, solroot, result_dbg, connected_dbg);
   pcsolPrune(g, result_dbg, connected_dbg);

   assert(LE(solstp_getObjBounded(g, result, 0.0, nedges), solstp_getObjBounded(g, result_dbg, 0.0, nedges)));

   SCIPfreeBufferArray(scip, &connected_dbg);
   SCIPfreeBufferArray(scip, &result_dbg);
#endif

   assert(solstp_isValid(scip, g, result));

   return SCIP_OKAY;
}



/* prune the (rooted) prize collecting Steiner tree in such a way that all leaves are terminals
 * NOTE: graph is not really const, mark is changed! todo */
static
SCIP_RETCODE pruneSteinerTreePc_csr(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   MST*                  mst,                /**< the MST */
   int* RESTRICT         result,             /**< ST edges (need to be set to UNKNOWN) */
   STP_Bool* RESTRICT    connected           /**< ST nodes */
   )
{
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(g);
   int solroot = g->source;

#ifndef NDEBUG
   {
      const int nedges_csr = graph_csr_getNedges(mst->csr);
      for( int i = 0; i < nedges_csr; i++ )
         assert(UNKNOWN == result[i]);

      assert(scip && result && connected);
      assert(g->extended);
      assert(graph_pc_isPcMw(g));
   }
#endif

   pcsolMarkGraphNodes_csr(g, connected);

   if( !rpcmw )
   {
      solroot = solstp_pcGetSolRoot(scip, g, connected);

      /* trivial solution? */
      if( solroot == -1 )
      {
         printf("trivial solution in pruning \n");
        // pcsolGetTrivialEdges(g, connected, result); // does not work for CSR!
         return SCIP_ERROR;
      }
   }

   assert(0 <= solroot && solroot < g->knots);
   assert(connected[solroot]);
   SCIPdebugMessage("(non-artificial) solution root=%d \n", solroot);

   solGetMstEdges_csr(g, connected, solroot, mst, result);

   SCIP_CALL( strongPruneSteinerTreePc_csr(scip, g, mst->csr, solroot, result, connected) );

   return SCIP_OKAY;
}

/*
 * Interface methods
 */



/** Gets root of solution for unrooted PC/MW.
 *  Returns -1 if the solution is empty. */
int solstp_pcGetSolRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const STP_Bool*       connected           /**< ST nodes */
   )
{
   int proot = -1;
   const int nnodes = graph_get_nNodes(g);

   assert(g && connected);
   assert(graph_pc_isPcMw(g));

   if( graph_pc_isRootedPcMw(g) )
   {
      return g->source;
   }

   /* todo remove this hack, better ask for the SCIP stage */
   if( SCIPprobdataGetNTerms(scip) == g->terms && SCIPprobdataGetNNodes(scip) == nnodes )
   {
      int min = nnodes;
      const int* termsorder = SCIPprobdataGetPctermsorder(scip);

      for( int k = 0; k < nnodes; k++ )
      {
         if( termsorder[k] < min && connected[k] )
         {
            assert(Is_pseudoTerm(g->term[k]));

            min = termsorder[k];
            proot = k;
         }
      }

      assert(min >= 0);
      assert(proot == -1 || min < nnodes);
   }
   else
   {
      const int* const gOeat = g->oeat;
      const int* const gHead = g->head;
      const int* const gTerm = g->term;
      const int groot = g->source;

      for( int a = g->outbeg[groot]; a >= 0; a = gOeat[a] )
      {
         const int head = gHead[a];
         if( !Is_term(gTerm[head]) && connected[head] )
         {
            proot = head;
            break;
         }
      }
   }

   return proot;
}


/** connects dummy terminals to given (pre-) PC solution */
void solstp_pcConnectDummies(
   const GRAPH*          g,                  /**< graph structure */
   int                   solroot,            /**< root of solution */
   int* RESTRICT         result,             /**< MST solution, which does not include artificial terminals */
   STP_Bool* RESTRICT    connected           /**< ST nodes */
   )
{
   const int* const gTerm = g->term;
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(g);
   const int gRoot = g->source;

   assert(graph_pc_isPcMw(g));

   /* connect all terminals */
   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(gTerm[i]) && i != gRoot )
      {
         const SCIP_Bool isFixedTerm = graph_pc_knotIsFixedTerm(g, i);

         assert(isFixedTerm || g->grad[i] == 2);
         assert(isFixedTerm || g->inpbeg[i] >= 0);

         if( isFixedTerm )
         {
            assert(rpcmw);
            assert(connected[i]);
            continue;
         }
         else
         {
            const int e1 = g->inpbeg[i];
            const int e2 = g->ieat[e1];
            const int k1 = g->tail[e1];
            const int k2 = g->tail[e2];

            connected[i] = TRUE;

            assert(graph_pc_knotIsDummyTerm(g, i));
            assert(g->ieat[e2] == EAT_LAST);
            assert(g->grad[i] == 2);
            assert(k1 == gRoot || k2 == gRoot);

            if( k1 != gRoot && connected[k1] )
               result[e1] = CONNECT;
            else if( k2 != gRoot && connected[k2] )
               result[e2] = CONNECT;
            else if( k1 == gRoot )
               result[e1] = CONNECT;
            else if( k2 == gRoot )
               result[e2] = CONNECT;

            /* xor: exactly one of e1 and e2 is used */
            assert((result[e1] != CONNECT) != (result[e2] != CONNECT));
         }
      }
      else if( i == solroot && !rpcmw )
      {
         int e1;
         for( e1 = g->inpbeg[i]; e1 != EAT_LAST; e1 = g->ieat[e1] )
         {
            if( g->tail[e1] == gRoot )
               break;
         }

         assert(e1 != EAT_LAST);
         result[e1] = CONNECT;
      }
   }

   if( !rpcmw )
      connected[gRoot] = TRUE;

   assert(connected[gRoot]);
}


/** (simple) pruning of given solution possible? */
SCIP_Bool stpsol_pruningIsPossible(
   const GRAPH*          g,                  /**< graph structure */
   const int*            result,             /**< ST edges */
   const STP_Bool*       connected           /**< ST nodes */
   )
{
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool isPc = graph_pc_isPc(g);

   for( int i = 0; i < nnodes; i++ )
   {
      int outedge;

      if( !connected[i] )
         continue;

      if( Is_term(g->term[i]) || Is_pseudoTerm(g->term[i]) )
         continue;

      for( outedge = g->outbeg[i]; outedge != EAT_LAST; outedge = g->oeat[outedge] )
         if( result[outedge] == CONNECT )
            break;

      if( outedge == EAT_LAST )
      {
         int e;

         if( isPc && graph_pc_knotIsNonLeafTerm(g, i) )
         {
            assert(g->prize && g->cost_org_pc);

            for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
               if( result[e] == CONNECT )
                  break;

            assert(e != EAT_LAST);

            if( GE(g->prize[i], g->cost_org_pc[e]) )
               continue;
         }

         for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
            if( result[e] == CONNECT )
               break;

         if( e != EAT_LAST && EQ(g->cost[e], 0.0) )
         {
            continue;
         }

         return TRUE;
      }
   }

   return FALSE;
}

/** Prune solution given by included nodes.
 *  NOTE: For PC/RPC this method will get the original edge costs before pruning! */
SCIP_RETCODE solstp_prune(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int*                  result,             /**< ST edges (out) */
   STP_Bool*             connected           /**< ST nodes (in/out) */
   )
{
   const int nedges = graph_get_nEdges(g);

   assert(scip && result && connected);
   assert(g->stp_type != STP_DHCSTP);

   for( int e = 0; e < nedges; e++ )
      result[e] = UNKNOWN;

   if( graph_pc_isPcMw(g) )
   {
      SCIP_Real* edgecosts = NULL;
      assert(g->extended);

      /* do we have biased edge costs? */
      if( graph_pc_isPc(g) )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &edgecosts, nedges) );

         graph_pc_getOrgCosts(scip, g, edgecosts);
      }
      else
      {
         edgecosts = g->cost;
      }

      SCIP_CALL( pruneSteinerTreePc(scip, g, edgecosts, result, connected) );

      if( graph_pc_isPc(g) )
         SCIPfreeBufferArray(scip, &edgecosts);
   }
   else
   {
      SCIP_CALL( pruneSteinerTreeStp(scip, g, g->cost, result, connected) );
   }

   assert(solstp_isValid(scip, g, result));

   return SCIP_OKAY;
}


/** prune solution given by included nodes */
SCIP_RETCODE solstp_pruneFromNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
   )
{
   assert(scip && g && result && connected);
   assert(g->stp_type != STP_DHCSTP);

   SCIP_CALL( solstp_prune(scip, g, result, connected) );

   return SCIP_OKAY;
}


/** prune solution given by included edges */
SCIP_RETCODE solstp_pruneFromEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int*                  result              /**< ST edges */
   )
{
   STP_Bool* connected;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);

   assert(scip && result);
   assert(solstp_isValid(scip, g, result));

   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      connected[k] = FALSE;

   for( int e = 0; e < nedges; e++ )
   {
      if( CONNECT == result[e] )
      {
         connected[g->head[e]] = TRUE;
         connected[g->tail[e]] = TRUE;
      }
   }

#ifdef SCIP_DEBUG
   SCIPdebugMessage("prune from edges: \n");
   solstp_print(g, result);
#endif

   SCIP_CALL( solstp_pruneFromNodes(scip, g, result, connected) );

   SCIPfreeBufferArray(scip, &connected);

   return SCIP_OKAY;
}


/** Prunes solution with respect to the provided edges costs.
 *  NOTE: method exists purely for optimization, so that unbiased costs for PC do not have to computed again! */
SCIP_RETCODE solstp_pruneFromTmHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs (original edge costs for PC!) */
   int* RESTRICT         result,             /**< ST edges */
   STP_Bool* RESTRICT    connected           /**< ST nodes */
   )
{
   const int nedges = graph_get_nEdges(g);

   assert(scip && result && connected);

   if( g->stp_type != STP_DHCSTP )
   {
      for( int e = 0; e < nedges; e++ )
         result[e] = UNKNOWN;
   }

   if( graph_pc_isPcMw(g) )
   {
      if( graph_pc_isPc(g) )
      {
         assert(cost);
         assert(graph_pc_costsEqualOrgCosts(scip, g, cost));
         SCIP_CALL( pruneSteinerTreePc(scip, g, cost, result, connected) );
      }
      else
      {
         assert(!cost);
         SCIP_CALL( pruneSteinerTreePc(scip, g, g->cost, result, connected) );
      }
   }
   else
      SCIP_CALL( pruneSteinerTreeStp(scip, g, (g->stp_type != STP_DHCSTP) ? g->cost : cost, result, connected) );

   return SCIP_OKAY;
}


/** Prunes solution with respect to the provided edges costs.
 *  CSR version! */
SCIP_RETCODE solstp_pruneFromTmHeur_csr(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SPATHS*               spaths,             /**< shortest paths;
                                                  NOTE: distances and preds not valid afterwards!
                                                  hacky, but improves cache-efficiency */
   int* RESTRICT         result              /**< ST edges */
   )
{
   const int nnodes = spaths->csr->nnodes;
   const int nedges = spaths->csr->start[nnodes];

   assert(scip && result);
   assert(nedges <= g->edges);

   assert(g->stp_type != STP_DHCSTP && "does not work because DHCSTP uses different edge costs");

   for( int e = 0; e < nedges; e++ )
      result[e] = UNKNOWN;

   if( graph_pc_isPcMw(g) )
   {
      MST mst = { .csr = spaths->csr_orgcosts, .dheap = spaths->dheap, .nodes_dist = spaths->nodes_dist,
                  .nodes_predEdge = spaths->nodes_pred};

      assert(graph_pc_isPc(g) || graph_csr_costsAreInSync(g, mst.csr, g->cost));

      SCIP_CALL( pruneSteinerTreePc_csr(scip, g, &mst, result, spaths->nodes_isConnected) );
   }
   else
   {
      MST mst = { .csr = spaths->csr_orgcosts, .dheap = spaths->dheap,
                  .nodes_dist = spaths->nodes_dist, .nodes_predEdge = spaths->nodes_pred };

      assert(graph_csr_costsAreInSync(g, mst.csr, g->cost));

      pruneSteinerTreeStp_csr(g, &mst, result, spaths->nodes_isConnected);
   }

   return SCIP_OKAY;
}


/** changes solution according to given root */
SCIP_RETCODE solstp_reroot(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int*                  result,             /**< solution array (CONNECT/UNKNOWN) */
   int                   newroot             /**< the new root */
   )
{
   int* queue;
   int* gmark;
   int size;
   const int nnodes = graph_get_nNodes(g);

   assert(scip != NULL);
   assert(g != NULL);
   assert(result != NULL);
   assert(newroot >= 0 && newroot < nnodes);
   assert(Is_term(g->term[newroot]));

   if( g->grad[newroot] == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &gmark, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &queue, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      gmark[k] = FALSE;

   gmark[newroot] = TRUE;
   size = 0;
   queue[size++] = newroot;

   /* BFS loop */
   while( size )
   {
      const int node = queue[--size];

      /* traverse outgoing arcs */
      for( int a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
      {
         const int head = g->head[a];

         if( !gmark[head] && (result[a] == CONNECT || result[flipedge(a)] == CONNECT ) )
         {
            if( result[flipedge(a)] == CONNECT  )
            {
               result[a] = CONNECT;
               result[flipedge(a)] = UNKNOWN;
            }
            gmark[head] = TRUE;
            queue[size++] = head;
         }
      }
   }

   /* adjust solution if infeasible */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !gmark[k] )
      {
         for( int a = g->outbeg[k]; a != EAT_LAST; a = g->oeat[a] )
         {
            result[a] = UNKNOWN;
            result[flipedge(a)] = UNKNOWN;
         }

         /* not yet connected terminal? */
         if( Is_term(g->term[k]) )
         {
            int a;
            assert(g->stp_type != STP_SPG);

            for( a = g->inpbeg[k]; a != EAT_LAST; a = g->ieat[a] )
            {
               const int node = g->tail[a];
               if( gmark[node] && node != newroot )
               {
                  result[a] = CONNECT;
                  break;
               }
            }
            if( a == EAT_LAST )
            {
               for( a = g->inpbeg[k]; a != EAT_LAST; a = g->ieat[a] )
               {
                  const int node = g->tail[a];
                  if( node == newroot )
                  {
                     result[a] = CONNECT;
                     break;
                  }
               }
            }
            else
               gmark[k] = TRUE;
         }
      }
   }

   SCIPfreeBufferArray(scip, &queue);
   SCIPfreeBufferArray(scip, &gmark);

#ifndef NDEBUG
   {
      const int realroot = g->source;
      g->source = newroot;
      assert(solstp_isValid(scip, g, result));
      g->source = realroot;
   }
#endif

   return SCIP_OKAY;
}


/** checks whether edge(s) of given primal solution have been deleted */
SCIP_Bool solstp_isUnreduced(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const int*            result              /**< solution array, indicating whether an edge is in the solution */
   )
{
   const int nedges = graph_get_nEdges(graph);

   assert(scip != NULL);
   assert(result != NULL);

   for( int i = 0; i < nedges; i++ )
   {
      if( result[i] == CONNECT && graph->oeat[i] == EAT_FREE )
         return FALSE;
   }

   return TRUE;
}


/** is the node contained in the solution? */
SCIP_Bool solstp_containsNode(
   const GRAPH*          g,                  /**< graph data structure */
   const int*            result,             /**< solution array, indicating whether an edge is in the solution */
   int                   node                /**< node to check for */
   )
{
   assert(g && result);
   assert(node >= 0 && node < g->knots);

   for( int e = g->outbeg[node]; e != EAT_LAST; e = g->oeat[e] )
   {
      if( result[e] == CONNECT || result[flipedge(e)] == CONNECT )
      {
         return TRUE;
      }
   }

   return FALSE;
}


/** verifies whether a given primal solution is feasible */
SCIP_Bool solstp_isValid(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const int*            result              /**< solution array, indicating whether an edge is in the solution */
   )
{
   int* queue = NULL;
   STP_Bool* reached = NULL;
   int size;
   int nterms;
   int termcount;
   const int nnodes = graph_get_nNodes(graph);
   const int root = graph->source;
   SCIP_Bool countpseudo;

   assert(scip && result);
   assert(root >= 0);

#ifndef NDEBUG
   for( int e = 0; e < graph->edges; ++e )
      assert(result[e] == CONNECT || result[e] == UNKNOWN);
#endif

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &reached, nnodes) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &queue, nnodes) );

   if( graph_pc_isPcMw(graph) && !graph->extended )
   {
      countpseudo = TRUE;
      nterms = graph_pc_nProperPotentialTerms(graph);

      if( !graph_pc_isRootedPcMw(graph) )
         nterms++;
   }
   else
   {
      countpseudo = FALSE;
      nterms = graph->terms;
   }

   for( int i = 0; i < nnodes; i++ )
      reached[i] = FALSE;

   /* BFS until all terminals are reached */

   termcount = 1;
   size = 0;
   reached[root] = TRUE;
   queue[size++] = root;

   while( size )
   {
      const int node = queue[--size];

      for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( result[e] == CONNECT )
         {
            const int i = graph->head[e];

            /* cycle? */
            if( reached[i] )
            {
               SCIPfreeBufferArray(scip, &queue);
               SCIPfreeBufferArray(scip, &reached);

               SCIPdebugMessage("solution contains a cycle ... \n");
               return FALSE;
            }

            if( countpseudo )
            {
               if( Is_pseudoTerm(graph->term[i]) || graph_pc_knotIsFixedTerm(graph, i) )
                  termcount++;
            }
            else
            {
               if( Is_term(graph->term[i]) )
                  termcount++;
            }

            reached[i] = TRUE;
            queue[size++] = i;
         }
      }
   }

#ifdef SCIP_DEBUG
   if( termcount != nterms )
   {
      printf("termcount %d graph->terms %d \n", termcount, nterms);
      printf("root %d \n", root);

      for( int i = 0; i < nnodes; i++ )
      {
         const int isMandatoryTerm = countpseudo?
               (Is_pseudoTerm(graph->term[i]) || graph_pc_knotIsFixedTerm(graph, i)) : Is_term(graph->term[i]);

         if( !reached[i] && isMandatoryTerm )
         {
            if( graph_pc_isPc(graph) && graph_pc_termIsNonLeafTerm(graph, i) )
               continue;

            printf("fail: ");
            graph_knot_printInfo(graph, i);

            for( int e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
            {
               printf("...neighbor: ");
               graph_knot_printInfo(graph, graph->tail[e]);
            }
         }
      }

      solstp_print(graph, result);
   }
#endif

   SCIPfreeBufferArray(scip, &queue);
   SCIPfreeBufferArray(scip, &reached);

   return (termcount == nterms);
}


/** prints given solution */
void solstp_print(
   const GRAPH*          graph,              /**< graph data structure */
   const int*            result              /**< solution array, indicating whether an edge is in the solution */
   )
{
   const int nedges = graph_get_nEdges(graph);

   assert(result);

   printf("solution tree edges: \n");

   for( int e = 0; e < nedges; ++e )
   {
      assert(result[e] == CONNECT || result[e] == UNKNOWN);

      if( CONNECT == result[e] )
      {
         printf("   ");
         graph_edge_printInfo(graph, e);
      }
   }
}


/** mark endpoints of edges in given list */
void solstp_setNodeList(
   const GRAPH*          g,              /**< graph data structure */
   STP_Bool*             solnode,        /**< solution nodes array (TRUE/FALSE) */
   IDX*                  listnode        /**< edge list */
   )
{
   int i;
   IDX* curr;

   assert(g != NULL);
   assert(solnode != NULL);

   curr = listnode;

   while( curr != NULL )
   {
      i = curr->index;

      solnode[g->head[i]] = TRUE;
      solnode[g->tail[i]] = TRUE;

      curr = curr->parent;
   }
}

/** compute solution value for given edge-solution array (CONNECT/UNKNOWN) and offset */
SCIP_Real solstp_getObjBounded(
   const GRAPH*          g,                  /**< the graph */
   const int*            soledge,            /**< solution */
   SCIP_Real             offset,             /**< offset */
   int                   nedges              /**< number of edges */
   )
{
   register SCIP_Real obj = offset;
   const SCIP_Real* const edgecost = g->cost;

   assert(nedges == g->edges);
   assert(!graph_pc_isPcMw(g) || g->extended);

   for( int e = 0; e < nedges; e++ )
   {
      assert(soledge[e] == CONNECT || soledge[e] == UNKNOWN);

      if( soledge[e] == CONNECT )
         obj += edgecost[e];
   }

   return obj;
}


/** compute solution value for given edge-solution array (CONNECT/UNKNOWN) and offset */
SCIP_Real solstp_getObj(
   const GRAPH*          g,                  /**< the graph */
   const int*            soledge,            /**< solution */
   SCIP_Real             offset             /**< offset */
   )
{
   assert(g);

   return solstp_getObjBounded(g, soledge, offset, g->edges);
}


/** compute solution value for given edge-solution array */
SCIP_Real solstp_pcGetObjCsr(
   const GRAPH*          g,                  /**< the graph */
   const CSR*            csr,                /**< the csr */
   const int*            soledge_csr,        /**< solution (CONNECT/UNKNOWN)  */
   const STP_Bool*       solnode             /**< solution vertices (TRUE/FALSE) */
   )
{
   const int nnodes = graph_get_nNodes(g);
   const int nedges_csr = graph_csr_getNedges(csr);
   const SCIP_Real* const edgecost_csr = csr->cost;
   const SCIP_Real* const prize = g->prize;
   const int* const gTerm = g->term;
   register SCIP_Real obj = 0.0;
   const SCIP_Bool isPc = graph_pc_isPc(g);

   assert(graph_pc_isPcMw(g));
   assert(g->extended);
   assert(nnodes == csr->nnodes);
   assert(nedges_csr <= g->edges);

   for( int e = 0; e < nedges_csr; e++ )
   {
      assert(soledge_csr[e] == CONNECT || soledge_csr[e] == UNKNOWN);

      if( soledge_csr[e] == CONNECT )
         obj += edgecost_csr[e];
   }

   for( int k = 0; k < nnodes; k++ )
   {
      if( solnode[k] )
         continue;

      if( Is_pseudoTerm(gTerm[k]) )
      {
         assert(g->grad[k] >= 1);
         obj += prize[k];
      }
      else if( isPc && Is_nonleafTerm(gTerm[k]) )
      {
         obj += prize[k];
      }
   }

   return obj;
}


/** compute solution value for given edge-solution array */
SCIP_Real solstp_getObjCsr(
   const GRAPH*          g,                  /**< the graph */
   const CSR*            csr,                /**< the csr */
   const int*            soledge_csr,        /**< solution (CONNECT/UNKNOWN)  */
   const STP_Bool*       solnode             /**< solution vertices (TRUE/FALSE) */
   )
{
   const int nedges_csr = graph_csr_getNedges(csr);
   const SCIP_Real* const edgecost_csr = csr->cost;
   register SCIP_Real obj = 0.0;

   assert(!graph_pc_isPcMw(g));
   assert(graph_get_nNodes(g) == csr->nnodes);
   assert(nedges_csr <= g->edges);

   for( int e = 0; e < nedges_csr; e++ )
   {
      assert(soledge_csr[e] == CONNECT || soledge_csr[e] == UNKNOWN);

      if( soledge_csr[e] == CONNECT )
         obj += edgecost_csr[e];
   }

   return obj;
}



/** converts solution from CSR to graph based */
void solstp_convertCsrToGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   const CSR*            csr,                /**< the CSR */
   const int*            soledge_csr,        /**< CSR solution (CONNECT/UNKNOWN)  */
   STP_Bool* RESTRICT    solnode,            /**< solution vertices (TRUE/FALSE) in/out! */
   int* RESTRICT         soledge_g           /**< graph solution (CONNECT/UNKNOWN) out */
)
{
   const int nedges_g = graph_get_nEdges(g);
   const int nedges_csr = graph_csr_getNedges(csr);
   const int* const edgeid = csr->edge_id;

   assert(solnode && soledge_csr && soledge_g);
   assert(edgeid);
   assert(0 <= nedges_csr && nedges_csr <= nedges_g);

   for( int i = 0; i < nedges_g; i++ )
   {
      soledge_g[i] = UNKNOWN;
   }

   for( int i = 0; i < nedges_csr; i++ )
   {
      if( CONNECT == soledge_csr[i] )
      {
         const int edge_g = edgeid[i];

         assert(0 <= edge_g && edge_g < nedges_g);
         assert(UNKNOWN == soledge_g[edge_g]);

         soledge_g[edge_g] = CONNECT;
      }
   }

   if( graph_pc_isPcMw(g) )
   {
      const int solroot = solstp_pcGetSolRoot(scip, g, solnode);
      solstp_pcConnectDummies(g, solroot, soledge_g, solnode);
   }
}


/** sets trivial solution (all UNKNOWN) */
void solstp_getTrivialSol(
   const GRAPH*          g,                  /**< the graph */
   int*                  soledge             /**< solution */
   )
{
   const int nedges = graph_get_nEdges(g);

   assert(soledge);

   for( int e = 0; e < nedges; e++ )
   {
      soledge[e] = UNKNOWN;
   }

   if( graph_pc_isPcMw(g) )
      pcsolGetTrivialEdges(g, NULL, soledge);
}


/** computes number of edges in solution value */
int solstp_getNedges(
   const GRAPH*          g,                  /**< the graph */
   const int*            soledge             /**< solution */
   )
{
   const int nedges = graph_get_nEdges(g);
   int edgecount = 0;

   assert(soledge);

   for( int e = 0; e < nedges; e++ )
      if( soledge[e] == CONNECT )
         edgecount++;

   return edgecount;
}

/** computes number of edges in solution value */
int solstp_getNedgesBounded(
   const GRAPH*          g,                  /**< the graph */
   const int*            soledge,            /**< solution */
   int                   nedges              /**< the (first) number of edges to consider */
   )
{
   int edgecount = 0;

   assert(nedges <= graph_get_nEdges(g));
   assert(soledge);

   for( int e = 0; e < nedges; e++ )
      if( soledge[e] == CONNECT )
         edgecount++;

   return edgecount;
}


/** marks vertices for given edge-solution array (CONNECT/UNKNOWN) */
void solstp_setVertexFromEdge(
   const GRAPH*          g,                  /**< the graph */
   const int*            result,             /**< solution array (CONNECT/UNKNOWN) */
   STP_Bool*             solnode             /**< marks whether node is in solution */
)
{
   const int nedges = g->edges;
   const int nnodes = g->knots;

   assert(g && result && solnode);

   for( int i = 0; i < nnodes; i++ )
      solnode[i] = FALSE;

   solnode[g->source] = TRUE;

   for( int e = 0; e < nedges; e++ )
   {
      if( result[e] == CONNECT )
      {
         assert(g->oeat[e] != EAT_FREE);

         solnode[g->head[e]] = TRUE;
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < nedges; e++ )
      if( result[e] == CONNECT )
         assert(solnode[g->head[e]] && solnode[g->tail[e]]);
#endif
}

/** get original solution */
SCIP_RETCODE solstp_getOrg(
   SCIP*           scip,               /**< SCIP data structure */
   const GRAPH*    transgraph,         /**< the transformed graph */
   const GRAPH*    orggraph,           /**< the original graph */
   const int*      transsoledge,       /**< solution for transformed problem */
   int*            orgsoledge          /**< new retransformed solution */
)
{
   STP_Bool* orgnodearr;
   IDX** const ancestors = transgraph->ancestors;

   const int transnedges = transgraph->edges;
   const int orgnnodes = orggraph->knots;
   const SCIP_Bool pcmw = graph_pc_isPcMw(transgraph);

   assert(transgraph != NULL && orggraph != NULL && transsoledge != NULL && orgsoledge != NULL);
   assert(transgraph->ancestors != NULL);
   assert(transgraph->stp_type == orggraph->stp_type);

   SCIP_CALL( SCIPallocBufferArray(scip, &orgnodearr, orgnnodes) );

   for( int k = 0; k < orgnnodes; k++ )
      orgnodearr[k] = FALSE;

   for( int e = 0; e < transnedges; e++ )
      if( transsoledge[e] == CONNECT )
         solstp_setNodeList(orggraph, orgnodearr, ancestors[e]);

   /* retransform edges fixed during graph reduction */
   solstp_setNodeList(orggraph, orgnodearr, graph_get_fixedges(transgraph));

   if( pcmw )
   {
      // potentially single-vertex solution?
      if( graph_pc_isRootedPcMw(transgraph) && transgraph->terms == 1 && graph_pc_nFixedTerms(orggraph) == 1 )
         orgnodearr[orggraph->source] = TRUE;

      SCIP_CALL( solstp_markPcancestors(scip, transgraph->pcancestors, orggraph->tail, orggraph->head, orgnnodes,
            orgnodearr, NULL, NULL, NULL, NULL ) );
   }

   /* prune solution (in original graph) */
   SCIP_CALL( solstp_prune(scip, orggraph, orgsoledge, orgnodearr) );

   SCIPfreeBufferArray(scip, &orgnodearr);

   assert(solstp_isValid(scip, orggraph, orgsoledge));

   return SCIP_OKAY;
}



/** mark original solution */
SCIP_RETCODE solstp_markPcancestors(
   SCIP*           scip,               /**< SCIP data structure */
   IDX**           pcancestors,        /**< the ancestors */
   const int*      tails,              /**< tails array */
   const int*      heads,              /**< heads array */
   int             orgnnodes,          /**< original number of nodes */
   STP_Bool*       solnodemark,        /**< solution nodes mark array */
   STP_Bool*       soledgemark,        /**< solution edges mark array or NULL */
   int*            solnodequeue,       /**< solution nodes queue or NULL  */
   int*            nsolnodes,          /**< number of solution nodes or NULL */
   int*            nsoledges           /**< number of solution edges or NULL */
)
{
   int* queue;
   int nnodes;
   int nedges = (nsoledges != NULL)? *nsoledges : 0;
   int qstart;
   int qend;

   assert(scip != NULL && tails != NULL && heads != NULL && pcancestors != NULL && solnodemark != NULL);

   if( solnodequeue != NULL )
      queue = solnodequeue;
   else
      SCIP_CALL( SCIPallocBufferArray(scip, &queue, orgnnodes) );

   if( nsolnodes == NULL )
   {
      assert(solnodequeue == NULL);
      nnodes = 0;
      for( int k = 0; k < orgnnodes; k++ )
         if( solnodemark[k] )
            queue[nnodes++] = k;
   }
   else
   {
      nnodes = *nsolnodes;
      assert(solnodequeue != NULL);
   }

   qstart = 0;
   qend = nnodes;

   while( qend != qstart )
   {
      int k = qstart;

      assert(qstart < qend);
      qstart = qend;

      for( ; k < qend; k++ )
      {
         const int ancestornode = queue[k];

         assert(solnodemark[ancestornode]);

         for( IDX* curr = pcancestors[ancestornode]; curr != NULL; curr = curr->parent )
         {
            const int ancestoredge = curr->index;
            assert(tails[ancestoredge] < orgnnodes && heads[ancestoredge] < orgnnodes);

            if( soledgemark != NULL && !soledgemark[ancestoredge] )
            {
               soledgemark[ancestoredge] = TRUE;
               nedges++;
            }
            if( !solnodemark[tails[ancestoredge]] )
            {
               solnodemark[tails[ancestoredge]] = TRUE;
               queue[nnodes++] = tails[ancestoredge];
            }
            if( !solnodemark[heads[ancestoredge]] )
            {
               solnodemark[heads[ancestoredge]] = TRUE;
               queue[nnodes++] = heads[ancestoredge];
            }
         }
      }
      qend = nnodes;
   }

   if( nsolnodes != NULL )
      *nsolnodes = nnodes;

   if( nsoledges != NULL )
      *nsoledges = nedges;

   if( solnodequeue == NULL )
      SCIPfreeBufferArray(scip, &queue);

   return SCIP_OKAY;
}
