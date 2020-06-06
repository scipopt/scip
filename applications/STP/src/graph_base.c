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

/**@file   graph_base.c
 * @brief  includes several methods for Steiner problem graphs
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file contains several basic methods to process Steiner problem graphs.
 * A graph can not be reduced in terms of edge or node size, but edges can be marked as
 * EAT_FREE (to not be used anymore) and nodes may have degree one.
 * The method 'graph_pack()' can be used to build a new graph that discards these nodes and edges.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h)                             */
//#define SCIP_DEBUG
#include "scip/misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "portab.h"
#include "misc_stp.h"
#include "graph.h"
#include "heur_tm.h"

#define STP_DELPSEUDO_NOEDGE    -1
#define STP_DELPSEUDO_SKIPEDGE  -2



/** internal data for pseudo-deletion */
typedef struct pseudo_deletion
{
   SCIP_Real*            ecost;              /**< edge cost */
   SCIP_Real*            ecostrev;           /**< reverse edge cost */
   SCIP_Real*            ecostreal;          /**< reverse edge cost */
   int*                  incedge;            /**< incident edges */
   int*                  adjvert;            /**< adjacent vertices */
   int*                  neigbedge;          /**< neighboring edges array */
   SCIP_Real             vertexprize;        /**< prize for PC */
   int                   degree;             /**< degree of vertex to be deleted */
   int                   ancestorsnode;      /**< ancestor node for PC */
} DELPSEUDO;


/*
 * local functions
 */


/** traverses the graph from vertex 'start' and marks all reached nodes */
static
SCIP_RETCODE trail(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the new graph */
   int                   start,              /**< node to start from */
   SCIP_Bool             costAware,          /**< be cost aware? */
   SCIP_Bool*            nodevisited         /**< marks which node has been visited */
   )
{
   int* stackarr;
   int stacksize;
   const int nnodes = g->knots;

   for( int i = 0; i < nnodes; i++ )
      nodevisited[i] = FALSE;

   nodevisited[start] = TRUE;

   if( g->grad[start] == 0 )
      return SCIP_OKAY;

   stacksize = 0;

   SCIP_CALL(SCIPallocBufferArray(scip, &stackarr, nnodes));

   stackarr[stacksize++] = start;

   /* DFS loop */
   while( stacksize != 0 )
   {
      const int node = stackarr[--stacksize];

      /* traverse outgoing arcs */
      for( int a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
      {
         const int head = g->head[a];

         if( !nodevisited[head] )
         {
            if( costAware && SCIPisGE(scip, g->cost[a], FARAWAY) )
               continue;

            nodevisited[head] = TRUE;
            stackarr[stacksize++] = head;
         }
      }
   }
   SCIPfreeBufferArray(scip, &stackarr);

   return SCIP_OKAY;
}


/** can edge in pseudo-elimination method be cut off? */
inline static
SCIP_Bool isCutoffEdge(
   SCIP*                 scip,               /**< SCIP data */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge */
   const SCIP_Real*      cutoffsrev,         /**< revere cutoff values (or NULL if undirected) */
   const SCIP_Real*      ecost,              /**< edge cost*/
   const SCIP_Real*      ecostrev,           /**< reverse edge cost */
   SCIP_Real             prize,              /**< prize if PcMw */
   int                   edgeidx1,           /**< index of first edge to be checked (wrt provided arrays) */
   int                   edgeidx2,           /**< index of second edge to be checked (wrt provided arrays) */
   int                   cutoffidx           /**< index for cutoff array */
   )
{
   SCIP_Real newcost;

   assert(edgeidx1 != edgeidx2);

   if( cutoffs == NULL )
      return FALSE;

   newcost = ecostrev[edgeidx1] + ecost[edgeidx2] - prize;

   if( !SCIPisGT(scip, newcost, cutoffs[cutoffidx]) )
      return FALSE;

   if( cutoffsrev != NULL )
   {
      const SCIP_Real newcostrev = ecost[edgeidx1] + ecostrev[edgeidx2];

      if( !SCIPisGT(scip, newcostrev, cutoffsrev[cutoffidx]) )
         return FALSE;
   }

   return TRUE;
}


inline static
void removeEdge(
   GRAPH*                g,                  /**< the graph */
   int                   e                   /**< the edge to be removed */
   )
{
   int    i;
   int    head;
   int    tail;

   assert(g          != NULL);
   assert(e          >= 0);
   assert(e          <  g->edges);

   head = g->head[e];
   tail = g->tail[e];

   if( g->inpbeg[head] == e )
      g->inpbeg[head] = g->ieat[e];
   else
   {
      if( g->rootedgeprevs != NULL && head == g->source )
      {
         i = g->rootedgeprevs[e];
         assert(g->ieat[i] == e);
         if( g->ieat[e] >= 0 )
            g->rootedgeprevs[g->ieat[e]] = i;
      }
      else
         for( i = g->inpbeg[head]; g->ieat[i] != e; i = g->ieat[i] )
            assert(i >= 0);

      g->ieat[i] = g->ieat[e];
   }
   if( g->outbeg[tail] == e )
      g->outbeg[tail] = g->oeat[e];
   else
   {
      if( g->rootedgeprevs != NULL && tail == g->source )
      {
         i = g->rootedgeprevs[e];
         assert(g->oeat[i] == e);
         if( g->oeat[e] >= 0 )
            g->rootedgeprevs[g->oeat[e]] = i;
      }
      else
         for( i = g->outbeg[tail]; g->oeat[i] != e; i = g->oeat[i] )
            assert(i >= 0);

      g->oeat[i] = g->oeat[e];
   }
}

/** used by graph_grid_create */
static
int getNodeNumber(
   int  grid_dim,
   int  shiftcoord,
   int* ncoords,
   int* currcoord
   )
{
   int number = 0;
   int tmp;
   int i;
   int j;
   for( i = 0; i < grid_dim; i++ )
   {
      tmp = 1;
      for( j = i + 1; j < grid_dim; j++ )
      {
         tmp = tmp * ncoords[j];
      }
      if( shiftcoord == i )
         number += (currcoord[i] + 1) * tmp;
      else
         number += currcoord[i] * tmp;
   }
   number++;
   return number;
}

/** used by graph_obstgrid_create */
static
void compEdgesObst(
   int   coord,
   int   grid_dim,
   int   nobstacles,
   int*  ncoords,
   int*  currcoord,
   int*  edgecosts,
   int*  gridedgecount,
   int** coords,
   int** gridedges,
   int** obst_coords,
   char* inobstacle
   )
{
   char inobst;
   int i;
   int j;
   int z;
   int x;
   int y;
   int node;
   i = 0;
   while( i < ncoords[coord] )
   {
      currcoord[coord] = i;
      if( coord < grid_dim - 1 )
         compEdgesObst(coord + 1, grid_dim, nobstacles, ncoords, currcoord, edgecosts, gridedgecount, coords, gridedges, obst_coords, inobstacle);
      else
      {
         x = coords[0][currcoord[0]];
         y = coords[1][currcoord[1]];
         inobst = FALSE;
         node = getNodeNumber(grid_dim, -1, ncoords, currcoord);
         for( z = 0; z < nobstacles; z++ )
         {
            assert(obst_coords[0][z] < obst_coords[2][z]);
            assert(obst_coords[1][z] < obst_coords[3][z]);
            if( x > obst_coords[0][z] && x < obst_coords[2][z] &&
               y > obst_coords[1][z] && y < obst_coords[3][z] )
            {
               inobst = TRUE;
               inobstacle[node-1] = TRUE;
               break;
            }
         }
         for( j = 0; j < grid_dim; j++ )
         {
            if( currcoord[j] + 1 < ncoords[j] )
            {
               if( inobst == FALSE )
               {
                  gridedges[0][*gridedgecount] = node;
                  gridedges[1][*gridedgecount] = getNodeNumber(grid_dim, j, ncoords, currcoord);
                  edgecosts[*gridedgecount] = coords[j][currcoord[j] + 1] - coords[j][currcoord[j]];
                  (*gridedgecount)++;
               }
            }
         }
      }
      i++;
   }
}

/** used by graph_grid_create */
static
void compEdges(
   int   coord,
   int   grid_dim,
   int*  ncoords,
   int*  currcoord,
   int*  edgecosts,
   int*  gridedgecount,
   int** coords,
   int** gridedges
   )
{
   int j;
   int i = 0;
   while( i < ncoords[coord] )
   {
      currcoord[coord] = i;
      if( coord < grid_dim - 1 )
         compEdges(coord + 1, grid_dim, ncoords, currcoord, edgecosts, gridedgecount, coords, gridedges);
      else
      {
         for( j = 0; j < grid_dim; j++ )
         {
            if( currcoord[j] + 1 < ncoords[j] )
            {
               gridedges[0][*gridedgecount] = getNodeNumber(grid_dim, -1, ncoords, currcoord);
               gridedges[1][*gridedgecount] = getNodeNumber(grid_dim, j, ncoords, currcoord);
               edgecosts[*gridedgecount] = coords[j][currcoord[j] + 1] - coords[j][currcoord[j]];
               (*gridedgecount)++;
            }
         }
      }
      i++;
   }
}

/** is the PC/MW part of the given graph valid? */
static
SCIP_Bool graphisValidPcMw(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the graph */
   SCIP_Bool*            nodevisited         /**< array */
   )
{
   int nFixedTerms = 0;
   int nProperTerms = 0;
   int nPseudoTerms = 0;
   int nNonLeafTerms = 0;
   const int nnodes = g->knots;
   const int root = g->source;
   const SCIP_Bool isExtended = g->extended;
   const SCIP_Bool isRooted = graph_pc_isRootedPcMw(g);

   assert(graph_pc_isPcMw(g));
   assert(g->prize != NULL);
   assert(g->term2edge != NULL);

   assert(graph_pc_term2edgeIsConsistent(scip, g));

#if 0
   if( !isRooted && g->grad[g->source] < 2 )
   {
      SCIPdebugMessage("only artificial root left \n");
      return FALSE;
   }
#endif

   for( int k = 0; k < nnodes; k++ )
   {
      if( graph_pc_knotIsFixedTerm(g, k) )
      {
         assert(isRooted || k == root);

         nFixedTerms++;

         continue;
      }

      assert(k != root);

      if( graph_pc_knotIsNonLeafTerm(g, k) )
      {
         nNonLeafTerms++;
      }
      /* is k a terminal in the original graph? */
      else if( (isExtended ? Is_pseudoTerm(g->term[k]) : Is_term(g->term[k])) )
      {
         nPseudoTerms++;
      }
      /* is k an artificial (newly added) terminal? */
      else if( (isExtended ? Is_term(g->term[k]) : Is_pseudoTerm(g->term[k])) )
      {
         int e;
         int e2;
         int pterm = -1;
         const int term = k;

         nProperTerms++;

         if( g->grad[k] != 2 )
         {
            SCIPdebugMessage("terminal degree != 2 for %d \n", k);
            return FALSE;
         }

         for( e = g->inpbeg[term]; e != EAT_LAST; e = g->ieat[e] )
            if( g->tail[e] == root )
               break;

         if( e == EAT_LAST )
         {
            SCIPdebugMessage("no edge to root for term %d \n", term);
            return FALSE;
         }

         for( e2 = g->outbeg[term]; e2 != EAT_LAST; e2 = g->oeat[e2] )
         {
            pterm = g->head[e2];
            if( (isExtended ? Is_pseudoTerm(g->term[pterm]) : Is_term(g->term[pterm])) && pterm != root  )
               break;
         }

         if( e2 == EAT_LAST)
         {
            SCIPdebugMessage("no terminal for dummy %d \n", g->head[e2]);
            return FALSE;
         }

         assert(pterm >= 0);
         assert(pterm != root);

         if( e2 != g->term2edge[term] )
         {
            SCIPdebugMessage("term2edge for node %d faulty \n", term);
            return FALSE;
         }

         if( !EQ(g->cost[e], g->prize[pterm]) )
         {
            SCIPdebugMessage("prize mismatch for node %d: %f!=%f \n", pterm, g->cost[e], g->prize[pterm]);
            return FALSE;
         }
      }
   } /* for k = 0; k < nnodes; k++ */

   if( nProperTerms != nPseudoTerms )
   {
      SCIPdebugMessage("wrong nPseudoTerms terminal count: %d != %d \n", nProperTerms, nPseudoTerms);

      return FALSE;
   }

   if( isExtended )
   {
      /* non-leaf terms are not part of g->terms */
      if( (nProperTerms + nFixedTerms) != g->terms )
      {
         SCIPdebugMessage("wrong overall terminal count(1): %d != %d \n", nProperTerms + nFixedTerms, g->terms);

         return FALSE;
      }
   }
   else
   {
      if( (nProperTerms + nNonLeafTerms + nFixedTerms) != g->terms )
      {
         SCIPdebugMessage("wrong overall terminal count(2): %d != %d \n", nProperTerms + nNonLeafTerms + nFixedTerms - 1, g->terms);

         return FALSE;
      }
   }

   if( isRooted )
   {
      SCIP_CALL_ABORT( graph_trail_costAware(scip, g, g->source, nodevisited) );

      for( int k = 0; k < nnodes; k++ )
      {
         if( graph_pc_knotIsFixedTerm(g, k) && !nodevisited[k] )
         {
            SCIPdebugMessage("disconnected fixed terminal %d \n", k);
            return FALSE;
         }
      }
   }

   if( isRooted && graph_pc_isMw(g) )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( !graph_pc_knotIsFixedTerm(g, k))
            continue;

         for( int e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         {
            if( !EQ(g->cost[e], 0.0) )
            {
               if( k == g->source && graph_pc_knotIsDummyTerm(g, g->tail[e]) )
                  continue;

               SCIPdebugMessage("non-zero incoming arc for fixed MW terminal %d \n", k);
               return FALSE;
            }
         }
      }
   }

   return TRUE;
}


/** do changes for Pc/Mw variants for vanished graph */
static
SCIP_RETCODE packPcMwVanished(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g_old,              /**< the old graph */
   GRAPH*                g_new               /**< the new graph */
   )
{
   int term = -1;
   const int nnodes_old = g_old->knots;

   if( graph_pc_isRootedPcMw(g_old) )
      return SCIP_OKAY;

   for( int i = 0; i < nnodes_old; ++i )
   {
      if( graph_pc_knotIsNonLeafTerm(g_old, i) )
      {
         assert(-1 == term);
         assert(g_old->source != i);

         term = i;
         break;
      }
   }

   assert(term >= 0);
   assert(graph_pc_termIsNonLeafTerm(g_old, term));

   SCIP_CALL( graph_fixed_addNodePc(scip, term, g_old) );

   return SCIP_OKAY;
}


/** do changes for Pc/Mw variants */
static
SCIP_RETCODE packPcMwInit(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nnodes_new,         /**< number of nodes */
   GRAPH*                g_old,              /**< the old graph */
   GRAPH*                g_new               /**< the new graph */
   )
{
   assert(nnodes_new > 0);
   assert(graph_pc_isPcMw(g_old));

   assert(g_new->ksize == nnodes_new);

   SCIP_CALL( graph_pc_initSubgraph(scip, g_new) );

   if( g_old->stp_type == STP_BRMWCSP )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g_new->costbudget), nnodes_new) );

   return SCIP_OKAY;
}

/** add nodes to new graph during graph packing */
static
void packNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g_old,              /**< the old graph */
   GRAPH*                g_new               /**< the new graph */
   )
{
   const int oldnnodes = graph_get_nNodes(g_old);
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(g_old);
   const SCIP_Bool pcmw = graph_pc_isPcMw(g_old);

   if( rpcmw )
   {
      for( int i = 0; i < oldnnodes; i++ )
         g_old->mark[i] = (g_old->grad[i] > 0);

      for( int e = g_old->outbeg[g_old->source]; e != EAT_LAST; e =
            g_old->oeat[e] )
      {
         const int i = g_old->head[e];

         if( SCIPisGT(scip, g_old->cost[e], 0.0) && Is_term(g_old->term[i])
               && g_old->term2edge[i] >= 0 )
         {
            g_old->mark[i] = FALSE;
            assert(g_old->grad[i] == 2);
         }
      }
   }

   for( int i = 0; i < oldnnodes; i++ )
   {
      if( g_old->grad[i] > 0 )
      {
         if( pcmw )
         {
            if( !Is_term(g_old->term[i]) || (rpcmw && g_old->mark[i]) )
               g_new->prize[g_new->knots] = g_old->prize[i];
            else
               g_new->prize[g_new->knots] = 0.0;
         }

         if( g_old->stp_type == STP_BRMWCSP )
         {
            if( !Is_term(g_old->term[i]) || (g_old->mark[i]) )
               g_new->costbudget[g_new->knots] = g_old->costbudget[i];
            else
               g_new->costbudget[g_new->knots] = 0.0;
         }
         graph_knot_add(g_new, g_old->term[i]);
      }
   }
}


/** add edges to new graph during graph packing */
static
SCIP_RETCODE packEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            old2newNode,        /**< node mapping */
   GRAPH*                g_old,              /**< the old graph */
   int                   nnodes,             /**< number of nodes for new graph */
   GRAPH*                g_new               /**< the new graph */
   )
{
   const int oldnedges = graph_get_nEdges(g_old);

   /* add edges */
   for( int i = 0; i < oldnedges; i += 2 )
   {
      int e;

      if( g_old->ieat[i] == EAT_FREE )
      {
         assert(g_old->oeat[i]     == EAT_FREE);
         assert(g_old->ieat[i + 1] == EAT_FREE);
         assert(g_old->oeat[i + 1] == EAT_FREE);

         graph_edge_delHistory(scip, g_old, i);
         continue;
      }

      assert(g_old->oeat[i]      != EAT_FREE);
      assert(g_old->ieat[i + 1]  != EAT_FREE);
      assert(g_old->oeat[i + 1]  != EAT_FREE);
      assert(old2newNode[g_old->tail[i]] >= 0);
      assert(old2newNode[g_old->head[i]] >= 0);

      e = g_new->edges;

      g_new->ancestors[e] = NULL;
      g_new->ancestors[e + 1] = NULL;
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g_new->ancestors[e]), g_old->ancestors[i], NULL) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g_new->ancestors[e + 1]), g_old->ancestors[i + 1], NULL) );

      assert(old2newNode[g_old->tail[i]] < nnodes && old2newNode[g_old->head[i]] < nnodes);

      graph_edge_addSubgraph(scip, g_old, old2newNode, i, g_new);
   }

   return SCIP_OKAY;
}



/** initializes */
static
SCIP_RETCODE delPseudoInit(
   SCIP*                 scip,               /**< SCIP data */
   const SCIP_Real*      edgecosts,          /**< edge costs for cutoff */
   int                   vertex,             /**< the vertex */
   GRAPH*                g,                  /**< graph */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   SCIP_Real* RESTRICT ecost;
   SCIP_Real* RESTRICT ecostrev;
   SCIP_Real* RESTRICT ecostreal;
   int* RESTRICT incedge;
   int* RESTRICT adjvert;
   int* RESTRICT neigbedge;
   int edgecount = 0;

#ifndef NDEBUG
   {
      int sum = 0;
      for( int i = 1; i < STP_DELPSEUDO_MAXGRAD; i++ )
         sum += i;
      assert(sum == STP_DELPSEUDO_MAXNEDGES);
   }
#endif

   if( Is_term(g->term[vertex]) )
   {
      assert(graph_pc_isPcMw(g));
      assert(!graph_pc_knotHasMaxPrize(g, vertex) );

      delpseudo->vertexprize = g->prize[vertex];
      delpseudo->ancestorsnode = vertex;

      /* NOTE: for degree 3 the deletion is always possible.
       * Thus, we can already manipulate the graph here */
      graph_pc_termToNonTerm(scip, g, vertex);

      delpseudo->degree = g->grad[vertex];

      assert(delpseudo->vertexprize > 0.0);
      assert(delpseudo->degree == 3);
   }
   else
   {
      delpseudo->vertexprize = 0.0;
      delpseudo->ancestorsnode = -1;
      delpseudo->degree = g->grad[vertex];

      if( g->pcancestors && g->pcancestors[vertex] )
      {
         assert(graph_pc_isPcMw(g));
         assert(SCIPisZero(scip, g->prize[vertex]));

         delpseudo->ancestorsnode = vertex;
      }
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ecost), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ecostrev), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(ecostreal), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(incedge), STP_DELPSEUDO_MAXGRAD) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(adjvert), STP_DELPSEUDO_MAXGRAD) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(neigbedge), STP_DELPSEUDO_MAXNEDGES) );

   /* save the state of all incident edges */
   for( int e = g->outbeg[vertex]; e != EAT_LAST; e = g->oeat[e] )
   {
      assert(e >= 0);

      incedge[edgecount] = e;
      ecostreal[edgecount] = g->cost[e];
      ecost[edgecount] = edgecosts[e];
      ecostrev[edgecount] = edgecosts[flipedge(e)];
      adjvert[edgecount++] = g->head[e];

      assert(edgecount <= STP_DELPSEUDO_MAXGRAD);
   }
   assert(edgecount == delpseudo->degree);

   delpseudo->ecost = ecost;
   delpseudo->ecostrev = ecostrev;
   delpseudo->ecostreal = ecostreal;
   delpseudo->incedge = incedge;
   delpseudo->adjvert = adjvert;
   delpseudo->neigbedge = neigbedge;

   return SCIP_OKAY;
}


/** tries to pseudo eliminate */
static
SCIP_RETCODE delPseudoDeleteVertex(
   SCIP*                 scip,               /**< SCIP data */
   int                   vertex,             /**< the vertex */
   GRAPH*                g,                  /**< graph */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   SINGLETONANS ancestors[STP_DELPSEUDO_MAXGRAD];
   const SCIP_Real* ecostreal = delpseudo->ecostreal;
   const int* incedge = delpseudo->incedge;
   const int* neigbedge = delpseudo->neigbedge;
   const int* adjvert = delpseudo->adjvert;
   const int degree = delpseudo->degree;
   const int nspareedges = degree; /* todo we might want to allow additional edges to be inserted */
   int edgecount = 0;
   int replacecount = 0;

   for( int i = 0; i < degree; i++ )
   {
      const int e = incedge[i];
      SCIP_CALL( graph_singletonAncestors_init(scip, g, e, &(ancestors[i])) );
   }

   assert(EQ(delpseudo->vertexprize, 0.0) || graph_pc_isPc(g));

   for( int i = 0; i < degree - 1; i++ )
   {
      for( int j = i + 1; j < degree; j++ )
      {
         const SCIP_Bool skipedge = (neigbedge[edgecount] == STP_DELPSEUDO_SKIPEDGE);

         /* do we need to insert edge at all? */
         if( !skipedge )
         {
            SCIP_Bool conflict;
            int newijedge;
            const SCIP_Real newijcost = ecostreal[i] + ecostreal[j] - delpseudo->vertexprize;
            const int oldincedge = incedge[(replacecount == nspareedges) ? replacecount - 1 : replacecount];
            const int oldijedge = neigbedge[edgecount];
#ifndef NDEBUG
            const int oldinctail = g->tail[oldincedge];
            const int oldinchead = g->head[oldincedge];
#endif
            assert(replacecount <= nspareedges);
            assert(replacecount < nspareedges || neigbedge[edgecount] != STP_DELPSEUDO_NOEDGE);

            SCIP_CALL( graph_edge_reinsert(scip, g, oldincedge, adjvert[i], adjvert[j], newijcost,
                  delpseudo->ancestorsnode, &(ancestors[i]), &(ancestors[j]), &newijedge, &conflict) );

            assert(!conflict);

            /* has a new edge been inserted (or the existing one been updated)? */
            if( newijedge >= 0 )
               SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, newijedge, vertex, g) );

            /* does no original edge exist? */
            if( oldijedge == STP_DELPSEUDO_NOEDGE )
            {
               replacecount++;
               assert(newijedge >= 0);
               assert(oldinctail != g->tail[oldincedge] || oldinchead != g->head[oldincedge]);
            }
            else
            {
               assert(newijedge == oldijedge || newijedge == - 1);
               assert(newijedge != oldincedge);
               assert(oldinctail == g->tail[oldincedge] && oldinchead == g->head[oldincedge]);
            }
         }

         edgecount++;
         assert(edgecount <= STP_DELPSEUDO_MAXNEDGES);
      }
   }

   /* delete remaining edges */
   graph_knot_del(scip, g, vertex, TRUE);

   for( int i = 0; i < degree; i++ )
      graph_singletonAncestors_freeMembers(scip, &(ancestors[i]));

   return SCIP_OKAY;
}


/** frees data */
static
void delPseudoFreeData(
   SCIP*                 scip,               /**< SCIP data */
   DELPSEUDO*            delpseudo           /**< data */
)
{
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->neigbedge), STP_DELPSEUDO_MAXNEDGES);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->adjvert), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->incedge), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->ecostreal), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->ecostrev), STP_DELPSEUDO_MAXGRAD);
   SCIPfreeBlockMemoryArray(scip, &(delpseudo->ecost), STP_DELPSEUDO_MAXGRAD);
}


/** gets replacement edge; helper function for pseudo-elimination */
static
SCIP_RETCODE delPseudoGetReplaceEdges(
   SCIP*                 scip,               /**< SCIP data */
   const GRAPH*          g,                  /**< graph */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge (or NULL) */
   const SCIP_Real*      cutoffsrev,         /**< reverse edge cutoff values (or NULL if undirected or non-existent) */
   DELPSEUDO*            delpseudo,          /**< data */
   SCIP_Bool*            success             /**< enough replace edges available?  */
)
{
   int* hasharr;
   int edgecount = 0;
   int replacecount = 0;
   const int degree = delpseudo->degree;
   const int nspareedges = degree;
   const int *incedge = delpseudo->incedge;
   const SCIP_Real *ecost = delpseudo->ecost;
   const SCIP_Real *ecostrev = delpseudo->ecostrev;
   const int *adjvert = delpseudo->adjvert;
   int* RESTRICT neigbedge = delpseudo->neigbedge;
   const SCIP_Real vertexprize = delpseudo->vertexprize;

   assert(scip && g && ecost && neigbedge);
   assert(degree >= 0 && degree <= STP_DELPSEUDO_MAXGRAD);

   *success = TRUE;

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, g->knots) );

   for( int i = 0; i < STP_DELPSEUDO_MAXNEDGES; i++ )
      neigbedge[i] = STP_DELPSEUDO_NOEDGE;

   for( int i = 0; i < degree - 1 && *success; i++ )
   {
      const int adjvertex = adjvert[i];
      const int iedge = incedge[i];

      graph_pseudoAncestors_hashEdge(g->pseudoancestors, iedge, hasharr);

      for( int j = i + 1; j < degree; j++ )
      {
         SCIP_Bool skipedge = isCutoffEdge(scip, cutoffs, cutoffsrev, ecost, ecostrev, vertexprize, i, j, edgecount);

         if( !skipedge )
         {
            const int jedge = incedge[j];
            skipedge = graph_pseudoAncestors_edgeIsHashed(g->pseudoancestors, jedge, hasharr);
         }

         edgecount++;
         assert(edgecount <= STP_DELPSEUDO_MAXNEDGES);

         /* can edge be discarded? */
         if( skipedge )
         {
            neigbedge[edgecount - 1] = STP_DELPSEUDO_SKIPEDGE;
         }
         else
         {
            int e;

            /* check whether edge already exists */
            for( e = g->outbeg[adjvertex]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( g->head[e] == adjvert[j] )
               {
                  assert(e >= 0);
                  neigbedge[edgecount - 1] = e;
                  break;
               }
            }

            if( e != EAT_LAST )
               continue;

            /* not enough spare edges? */
            if( ++replacecount > nspareedges )
            {
               *success = FALSE;
               break;
            }
         }
      } /* inner neighbor loop */

      graph_pseudoAncestors_unhashEdge(g->pseudoancestors, iedge, hasharr);
   }  /* outer neighbor loop */

   SCIPfreeCleanBufferArray(scip, &hasharr);

   return SCIP_OKAY;
}


/*
 * global functions
 */


/** creates a graph out of a given grid */
SCIP_RETCODE graph_obstgrid_create(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               gridgraph,          /**< the (obstacle) grid graph to be constructed */
   int**                 coords,             /**< coordinates of all points  */
   int**                 obst_coords,        /**< coordinates of obstacles */
   int                   nterms,             /**< number of terminals */
   int                   grid_dim,           /**< dimension of the problem */
   int                   nobstacles,         /**< number of obstacles*/
   int                   scale_order         /**< scale factor */
   )
{
   GRAPH* g;
   GRAPH* gp;
   double cost;
   int    i;
   int    j;
   int    k;
   int    tmp;
   int    shift;
   int    nnodes;
   int    nedges;
   double  scale_factor;
   int    gridedgecount;
   int*   ncoords;
   int*   currcoord;
   int*   edgecosts;
   int**  termcoords;
   int**  gridedges;
   char*  inobstacle;
   assert(coords != NULL);
   assert(grid_dim > 1);
   assert(nterms > 0);
   assert(grid_dim == 2);
   scale_factor = pow(10.0, (double) scale_order);

   /* initialize the terminal-coordinates array */
   SCIP_CALL( SCIPallocBufferArray(scip, &termcoords, grid_dim) );

   for( i = 0; i < grid_dim; i++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(termcoords[i]), nterms) ); /*lint !e866*/
      for( j = 0; j < nterms; j++ )
         termcoords[i][j] = coords[i][j];
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &ncoords, grid_dim) );
   SCIP_CALL( SCIPallocBufferArray(scip, &currcoord, grid_dim) );

   /* sort the coordinates and delete multiples */
   for( i = 0; i < grid_dim; i++ )
   {
      ncoords[i] = 1;
      SCIPsortInt(coords[i], nterms);
      shift = 0;
      for( j = 0; j < nterms - 1; j++ )
      {
         if( coords[i][j] == coords[i][j + 1] )
         {
            shift++;
         }
         else
         {
            coords[i][j + 1 - shift] = coords[i][j + 1];
            ncoords[i]++;
         }
      }
   }

   nnodes = 1;

   for( i = 0; i < grid_dim; i++ )
      nnodes = nnodes * ncoords[i];

   tmp = 0;

   for( i = 0; i < grid_dim; i++ )
      tmp = tmp + nnodes / ncoords[i];

   nedges = grid_dim * nnodes - tmp;
   SCIP_CALL( SCIPallocBufferArray(scip, &gridedges, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgecosts, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(gridedges[0]), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(gridedges[1]), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(inobstacle), nnodes) );
   gridedgecount = 0;
   for( i = 0; i < nnodes; i++ )
      inobstacle[i] = FALSE;
   compEdgesObst(0, grid_dim, nobstacles, ncoords, currcoord, edgecosts, &gridedgecount, coords, gridedges, obst_coords, inobstacle);
   nedges = gridedgecount;
   /* initialize empty g with allocated slots for nodes and edges */
   SCIP_CALL( graph_init(scip, gridgraph, nnodes, 2 * nedges, 1) );

   g = *gridgraph;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->grid_ncoords), grid_dim) );
   for( i = 0; i < grid_dim; i++ )
      g->grid_ncoords[i] = ncoords[i];

   g->grid_dim = grid_dim;
   g->grid_coordinates = coords;

   /* add nodes */
   for( i = 0; i < nnodes; i++ )
      graph_knot_add(g, -1);

   /* add edges */
   for( i = 0; i < nedges; i++ )
   {
      /* (re) scale edge costs */
      if( inobstacle[gridedges[1][i] - 1] == FALSE )
      {
         cost = ((double) edgecosts[i]) / scale_factor;
         graph_edge_add(scip, g, gridedges[0][i] - 1, gridedges[1][i] - 1, cost, cost);
      }
   }

   /* add terminals */
   for( i = 0; i < nterms; i++ )
   {
      for( j = 0; j < grid_dim; j++ )
      {
         for( k = 0; k <= ncoords[j]; k++ )
         {
            assert(k != ncoords[j]);
            if( coords[j][k] == termcoords[j][i] )
            {
               currcoord[j] = k;
               break;
            }
         }
      }
      /* the position of the (future) terminal */
      k = getNodeNumber(grid_dim, -1, ncoords, currcoord) - 1;

      if( i == 0 )
         g->source = k;

      /* make a terminal out of the node */
      graph_knot_chg(g, k, 0);
   }

   SCIP_CALL( graph_pack(scip, g, &gp, NULL, TRUE) );
   g = gp;
   g->stp_type = STP_OARSMT;

   SCIPfreeBufferArray(scip, &inobstacle);
   SCIPfreeBufferArray(scip, &(gridedges[1]));
   SCIPfreeBufferArray(scip, &(gridedges[0]));
   SCIPfreeBufferArray(scip, &edgecosts);
   SCIPfreeBufferArray(scip, &gridedges);
   SCIPfreeBufferArray(scip, &currcoord);
   SCIPfreeBufferArray(scip, &ncoords);

   for( i = grid_dim - 1; i >= 0 ; --i )
      SCIPfreeBufferArray(scip, &(termcoords[i]));

   SCIPfreeBufferArray(scip, &termcoords);

   return SCIP_OKAY;
}



/** creates a graph out of a given grid */
SCIP_RETCODE graph_grid_create(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               gridgraph,          /**< the grid graph to be constructed */
   int**                 coords,             /**< coordinates */
   int                   nterms,             /**< number of terminals*/
   int                   grid_dim,           /**< problem dimension */
   int                   scale_order         /**< scale order */
   )
{
   GRAPH* g;
   double cost;
   int    i;
   int    j;
   int    k;
   int    tmp;
   int    shift;
   int    nnodes;
   int    nedges;
   double  scale_factor;
   int    gridedgecount;
   int*   ncoords;
   int*   currcoord;
   int*   edgecosts;
   int**  termcoords;
   int**  gridedges;
   assert(coords != NULL);
   assert(grid_dim > 1);
   assert(nterms > 0);

   scale_factor = pow(10.0, (double) scale_order);

   /* initialize the terminal-coordinates array */
   SCIP_CALL( SCIPallocBufferArray(scip, &termcoords, grid_dim) );
   for( i = 0; i < grid_dim; i++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(termcoords[i]), nterms) ); /*lint !e866*/
      for( j = 0; j < nterms; j++ )
         termcoords[i][j] = coords[i][j];
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &ncoords, grid_dim) );
   SCIP_CALL( SCIPallocBufferArray(scip, &currcoord, grid_dim) );

   /* sort the coordinates and delete multiples */
   for( i = 0; i < grid_dim; i++ )
   {
      ncoords[i] = 1;
      SCIPsortInt(coords[i], nterms);
      shift = 0;
      for( j = 0; j < nterms - 1; j++ )
      {
         if( coords[i][j] == coords[i][j + 1] )
         {
            shift++;
         }
         else
         {
            coords[i][j + 1 - shift] = coords[i][j + 1];
            ncoords[i]++;
         }
      }
   }

   nnodes = 1;
   for( i = 0; i < grid_dim; i++ )
      nnodes = nnodes * ncoords[i];

   tmp = 0;
   for( i = 0; i < grid_dim; i++ )
   {
      tmp = tmp + nnodes / ncoords[i];
   }

   nedges = grid_dim * nnodes - tmp;

   SCIP_CALL( SCIPallocBufferArray(scip, &gridedges, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgecosts, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(gridedges[0]), nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(gridedges[1]), nedges) );

   gridedgecount = 0;

   compEdges(0, grid_dim, ncoords, currcoord, edgecosts, &gridedgecount, coords, gridedges);

   /* initialize empty graph with allocated slots for nodes and edges */
   SCIP_CALL( graph_init(scip, gridgraph, nnodes, 2 * nedges, 1) );

   g = *gridgraph;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->grid_ncoords), grid_dim) );
   for( i = 0; i < grid_dim; i++ )
      g->grid_ncoords[i] = ncoords[i];

   g->grid_dim = grid_dim;
   g->grid_coordinates = coords;

   /* add nodes */
   for( i = 0; i < nnodes; i++ )
      graph_knot_add(g, -1);

   /* add edges */
   for( i = 0; i < nedges; i++ )
   {
      /* (re) scale edge costs */
      cost = (double) edgecosts[i] / scale_factor;
      graph_edge_add(scip, g, gridedges[0][i] - 1, gridedges[1][i] - 1, cost, cost);
   }

   /* add terminals */
   for( i = 0; i < nterms; i++ )
   {
      for( j = 0; j < grid_dim; j++ )
      {
         for( k = 0; k <= ncoords[j]; k++ )
         {
            assert(k != ncoords[j]);
            if( coords[j][k] == termcoords[j][i] )
            {
               currcoord[j] = k;
               break;
            }
         }
      }
      /* the position of the (future) terminal */
      k = getNodeNumber(grid_dim, -1, ncoords, currcoord) - 1;

      /* make a terminal out of the node */
      graph_knot_chg(g, k, 0);
   }

   g->stp_type = STP_RSMT;

   SCIPfreeBufferArray(scip, &(gridedges[1]));
   SCIPfreeBufferArray(scip, &(gridedges[0]));
   SCIPfreeBufferArray(scip, &edgecosts);
   SCIPfreeBufferArray(scip, &gridedges);
   SCIPfreeBufferArray(scip, &currcoord);
   SCIPfreeBufferArray(scip, &ncoords);

   for( i = grid_dim - 1; i >= 0 ; i-- )
      SCIPfreeBufferArray(scip, &(termcoords[i]));

   SCIPfreeBufferArray(scip, &termcoords);

   return SCIP_OKAY;
}


/** computes coordinates of node 'node' */
SCIP_RETCODE graph_grid_coordinates(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 coords,             /**< coordinates */
   int**                 nodecoords,         /**< coordinates of the node (to be computed) */
   int*                  ncoords,            /**< array with number of coordinate */
   int                   node,               /**< the node */
   int                   grid_dim            /**< problem dimension */
   )
{
   int i;
   int j;
   int tmp;
   int coord;
   assert(grid_dim > 1);
   assert(node >= 0);
   assert(coords != NULL);
   assert(ncoords != NULL);
   if( *nodecoords == NULL )
      SCIP_CALL( SCIPallocMemoryArray(scip, nodecoords, grid_dim) );

   for( i = 0; i < grid_dim; i++ )
   {
      tmp = 1;
      for( j = i; j < grid_dim; j++ )
         tmp = tmp * ncoords[j];

      coord = node % tmp;
      tmp = tmp / ncoords[i];
      coord = coord / tmp;
      (*nodecoords)[i] = coords[i][coord];
   }
   return SCIP_OKAY;
}


/* is the vertex a leaf (for NWPTSPG) */
SCIP_Bool graph_nw_knotIsLeaf(
   const GRAPH*          g,                  /**< the graph */
   int                   vertex              /**< the vertex  */
)
{
   int e;
   assert(g != NULL && g->stp_type == STP_NWPTSPG);

   for( e = g->outbeg[vertex]; e != EAT_LAST; e = g->oeat[e] )
      if( g->cost[e] < FARAWAY )
         break;

   return (e == EAT_LAST );
}


/** is node in range? */
SCIP_Bool graph_knot_isInRange(
   const GRAPH*          g,                  /**< the graph */
   int                   k                   /**< the node */
   )
{
   assert(g);

   return (0 <= k && k < g->knots);
}


/** add a vertex */
void graph_knot_add(
   GRAPH*                p,                  /**< the graph */
   int                   term                /**< terminal property */
   )
{
   assert(p        != NULL);
   assert(p->ksize >  p->knots);
   assert(term     <  p->layers);

   p->term  [p->knots] = term;
   p->mark  [p->knots] = TRUE;
   p->grad  [p->knots] = 0;
   p->inpbeg[p->knots] = EAT_LAST;
   p->outbeg[p->knots] = EAT_LAST;

   if( Is_term(term) )
      p->terms++;

   p->knots++;
}

/** change terminal property of a vertex */
void graph_knot_chg(
   GRAPH*                p,                  /**< the graph */
   int                   node,               /**< node to be changed */
   int                   term                /**< terminal property */
   )
{
   assert(p      != NULL);
   assert(node   >= 0);
   assert(node   < p->knots);
   assert(term == STP_TERM || term == STP_TERM_NONE || term == STP_TERM_NONLEAF || term == STP_TERM_PSEUDO);

   if( term != p->term[node] )
   {
      if( Is_term(p->term[node]) )
         p->terms--;

      p->term[node] = term;

      if( Is_term(p->term[node]) )
         p->terms++;
   }
}


/** delete node */
void graph_knot_del(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   k,                  /**< the node */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   assert(scip && g);
   assert(k >= 0 && k < g->knots);

   while( g->outbeg[k] != EAT_LAST )
      graph_edge_del(scip, g, g->outbeg[k], freeancestors);

   assert(g->grad[k] == 0);
   assert(g->outbeg[k] == EAT_LAST);
   assert(g->inpbeg[k] == EAT_LAST);
}


/** pseudo deletes non-terminal of degree 2 */
SCIP_RETCODE graph_knot_replaceDeg2(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   vertex,             /**< the vertex */
   GRAPH*                g,                  /**< the graph */
   int*                  solnode,            /**< marks whether an node is part of a given solution (CONNECT), or is NULL */
   SCIP_Bool*            edgeEliminated      /**< edge eliminated? (due to conflict) */
   )
{
   SINGLETONANS ancestors1;
   SINGLETONANS ancestors2;
   const int e1 = g->outbeg[vertex];
   const int e2 = g->oeat[e1];
   const int i1 = g->head[e1];
   const int i2 = g->head[e2];
   int newedge;
   SCIP_Bool conflict;

   assert(scip && g);
   assert(vertex >= 0 && vertex < g->knots);
   assert(!Is_term(g->term[vertex]));
   assert(g->grad[vertex] == 2);
   assert(e1 >= 0 && e2 >= 0);
   assert(SCIPisEQ(scip, g->cost[e1], g->cost[flipedge(e1)]));
   assert(SCIPisEQ(scip, g->cost[e2], g->cost[flipedge(e2)]));
   assert(graph_valid_pseudoAncestors(scip, g));

   *edgeEliminated = FALSE;

   SCIP_CALL( graph_singletonAncestors_init(scip, g, e1, &(ancestors1)) );
   SCIP_CALL( graph_singletonAncestors_init(scip, g, e2, &(ancestors2)) );

   SCIP_CALL( graph_edge_reinsert(scip, g, e1, i2, i1, g->cost[e1] + g->cost[e2],
         -1, &ancestors2, &ancestors1, &newedge, &conflict) );

   graph_singletonAncestors_freeMembers(scip, &(ancestors1));
   graph_singletonAncestors_freeMembers(scip, &(ancestors2));

   graph_knot_del(scip, g, vertex, TRUE);

   if( conflict )
   {
      SCIPdebugMessage("conflict in graph_knot_replaceDeg2 \n");
      assert(newedge >= 0);
      graph_edge_del(scip, g, newedge, TRUE);
      *edgeEliminated = TRUE;
   }

   return SCIP_OKAY;
}


/** pseudo delete node, i.e. reconnect neighbors; maximum degree of STP_DELPSEUDO_MAXGRAD! */
SCIP_RETCODE graph_knot_delPseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   const SCIP_Real*      edgecosts,          /**< edge costs for cutoff */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge (or NULL) */
   const SCIP_Real*      cutoffsrev,         /**< reverse edge cutoff values (or NULL if undirected or non-existent) */
   int                   vertex,             /**< the vertex */
   SCIP_Bool*            success             /**< has node been pseudo-eliminated? */
   )
{
   DELPSEUDO delpseudo = { NULL, NULL, NULL, NULL, NULL, NULL, -1.0, -1, -1 };

   assert(scip && success && g);
   assert(vertex >= 0 && vertex < g->knots);
   assert(g->grad[vertex] <= STP_DELPSEUDO_MAXGRAD);

   *success = TRUE;

   if( g->grad[vertex] <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( delPseudoInit(scip, edgecosts, vertex, g, &delpseudo) );
   SCIP_CALL( delPseudoGetReplaceEdges(scip, g, cutoffs, cutoffsrev, &delpseudo, success) );

   /* enough spare edges? */
   if( (*success) )
   {
      SCIP_CALL( delPseudoDeleteVertex(scip, vertex, g, &delpseudo) );
   }

   delPseudoFreeData(scip, &delpseudo);

   return SCIP_OKAY;
}


/** contracts node s into node t */
SCIP_RETCODE graph_knot_contract(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                p,                  /**< graph data structure */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                or NULL */
   int                   t,                  /**< tail node to be contracted (surviving node) */
   int                   s                   /**< head node to be contracted */
   )
{
   SCIP_Real* incost = NULL;
   SCIP_Real* outcost = NULL;
   SINGLETONANS* ancestors = NULL;
   int* mark = NULL;
   int* edge = NULL;
   int* knot = NULL;
   int slc = 0;
   int sgrad;

   assert(p && scip);
   assert(t >= 0 && t < p->knots);
   assert(s >= 0 && s < p->knots);
   assert(p->grad[s] > 0 && p->grad[t] > 0);
   assert(p->layers == 1);

   /* save solution */
   if( solnode != NULL )
      if( solnode[s] == CONNECT )
         solnode[t] = CONNECT;

   /* change terminal property */
   if( Is_term(p->term[s]) )
   {
      graph_knot_chg(p, t, p->term[s]);
      graph_knot_chg(p, s, -1);
   }

   /* retain root */
   if( p->source == s )
      p->source = t;

   sgrad = p->grad[s];
   if( sgrad >= 2 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &incost, sgrad - 1) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &outcost, sgrad - 1) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &mark, sgrad - 1) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &edge, sgrad - 1) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &knot, sgrad - 1) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ancestors, sgrad - 1) );
   }

   /* store edges to be moved/removed */
   for( int es = p->outbeg[s]; es != EAT_LAST; es = p->oeat[es] )
   {
      assert(p->tail[es] == s);

      if( p->head[es] != t )
      {
         assert(ancestors && mark && incost && outcost && edge && knot);

         SCIP_CALL( graph_singletonAncestors_init(scip, p, es, &(ancestors[slc])) );
         mark[slc] = FALSE;
         edge[slc] = es;
         knot[slc] = p->head[es];
         outcost[slc] = p->cost[es];
         incost[slc] = p->cost[Edge_anti(es)];
         slc++;

         assert(slc < sgrad);
      }
   }

   assert(slc == sgrad - 1);

   /* traverse edges */
   for( int i = 0; i < slc; i++ )
   {
      int et;
      const int ihead = knot[i];

      assert(knot != NULL && outcost != NULL && incost != NULL && mark != NULL);

      /* search for an edge out of t with same head as current edge */

      if( p->grad[ihead] >= p->grad[t] )
      {
         for( et = p->outbeg[t]; et >= 0; et = p->oeat[et] )
            if( p->head[et] == ihead )
               break;
      }
      else
      {
         for( et = p->inpbeg[ihead]; et >= 0; et = p->ieat[et] )
            if( p->tail[et] == t )
               break;
      }

      /* does such an edge not exist? */
      if( et == EAT_LAST )
      {
         mark[i] = TRUE;
      }
      else
      {
         SCIP_Bool copyPseudoancestors = FALSE;
         assert(et != EAT_LAST);

         /* This is for nodes with edges to s and t.
          * Need to adjust the out and in costs of the edge
          */

         if( graph_typeIsUndirected(p) && SCIPisGT(scip, p->cost[et], outcost[i]) && SCIPisGT(scip, p->cost[Edge_anti(et)], incost[i]) )
            copyPseudoancestors = TRUE;

         if( copyPseudoancestors )
            graph_edge_delPseudoAncestors(scip, et, p);

         if( SCIPisGT(scip, p->cost[et], outcost[i]) )
         {
            SCIPintListNodeFree(scip, &((p->ancestors)[et]));
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[et]), ancestors[i].ancestors, NULL) );

            assert(graph_edge_nPseudoAncestors(p, et) == 0);

            p->cost[et] = outcost[i];
         }
         if( SCIPisGT(scip, p->cost[Edge_anti(et)], incost[i]) )
         {
            const int anti = Edge_anti(et);

            SCIPintListNodeFree(scip, &(p->ancestors[anti]));
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[anti]), ancestors[i].revancestors, NULL) );

            assert(graph_edge_nPseudoAncestors(p, anti) == 0);

            p->cost[anti] = incost[i];
         }

         if( copyPseudoancestors )
         {
            SCIP_Bool conflict;
            SCIP_CALL( graph_pseudoAncestors_appendCopySingToEdge(scip, et, &(ancestors[i]), FALSE, p, &conflict) );
            assert(!conflict);
         }
      }
   }

   /* insert edges */
   for( int i = 0; i < slc; i++ )
   {
      assert(mark && ancestors && knot && outcost && incost);

      if( mark[i] )
      {
         int es = p->outbeg[s];
         const int head = knot[i];
         const int tail = t;
         SCIP_Bool conflict;

         assert(es != EAT_LAST);

         graph_edge_del(scip, p, es, TRUE);

         assert(p->ancestors[es] == NULL);
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), ancestors[i].ancestors, NULL) );

         SCIP_CALL( graph_pseudoAncestors_appendCopySingToEdge(scip, es, &(ancestors[i]), FALSE, p, &conflict) );
         assert(!conflict);

         p->grad[head]++;
         p->grad[tail]++;

         p->cost[es]     = outcost[i];
         p->tail[es]     = tail;
         p->head[es]     = head;
         p->ieat[es]     = p->inpbeg[head];
         p->oeat[es]     = p->outbeg[tail];
         p->inpbeg[head] = es;
         p->outbeg[tail] = es;

         es = Edge_anti(es);

         assert(p->ancestors[es] == NULL);
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), ancestors[i].revancestors, NULL) );

         p->cost[es]     = incost[i];
         p->tail[es]     = head;
         p->head[es]     = tail;
         p->ieat[es]     = p->inpbeg[tail];
         p->oeat[es]     = p->outbeg[head];
         p->inpbeg[tail] = es;
         p->outbeg[head] = es;
      }
   }

   /* delete remaining edges */
   graph_knot_del(scip, p, s, TRUE);

   if( sgrad >= 2 )
   {
      assert(ancestors && knot && outcost && edge && mark && incost);

      for( int i = 0; i < slc; i++ )
         graph_singletonAncestors_freeMembers(scip, &(ancestors[i]));

      SCIPfreeBlockMemoryArray(scip, &ancestors, sgrad - 1);
      SCIPfreeBlockMemoryArray(scip, &knot, sgrad - 1);
      SCIPfreeBlockMemoryArray(scip, &edge, sgrad - 1);
      SCIPfreeBlockMemoryArray(scip, &mark, sgrad - 1);
      SCIPfreeBlockMemoryArray(scip, &outcost, sgrad - 1);
      SCIPfreeBlockMemoryArray(scip, &incost, sgrad - 1);
   }

   return SCIP_OKAY;
}


/** contract an edge, given index and by its endpoints, which is to be fixed */
SCIP_RETCODE graph_knot_contractFixed(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                or NULL */
   int                   edge,               /**< the edge */
   int                   t,                  /**< tail node to be contracted */
   int                   s                   /**< head node to be contracted */
   )
{
   SCIP_CALL( graph_fixed_addEdge(scip, edge, g) );
   SCIP_CALL( graph_knot_contract(scip, g, solnode, t, s) );

   return SCIP_OKAY;
}


/** contract endpoint of lower degree into endpoint of higher degree */
SCIP_RETCODE graph_knot_contractLowdeg2High(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                or NULL */
   int                   t,                  /**< tail node to be contracted */
   int                   s                   /**< head node to be contracted */
   )
{
   assert(g != NULL);

   if( g->grad[t] >= g->grad[s] )
      SCIP_CALL( graph_knot_contract(scip, g, solnode, t, s) );
   else
      SCIP_CALL( graph_knot_contract(scip, g, solnode, s, t) );

   return SCIP_OKAY;
}


/** redirects given edge eki */
int graph_edge_redirect(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   eki,                /**< the edge */
   int                   k,                  /**< new tail */
   int                   j,                  /**< new head */
   SCIP_Real             cost,               /**< new cost */
   SCIP_Bool             forcedelete,        /**< delete edge eki if it is not used? */
   SCIP_Bool             checkexist          /**< check if there is already an edge kj */
   )
{
   int e;

   assert(scip && g);
   assert(SCIPisGE(scip, cost, 0.0));

   if( forcedelete )
      graph_edge_del(NULL, g, eki, FALSE);

   if( checkexist )
   {
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(g->tail[e] == k);

         if( g->head[e] == j )
            break;
      }
   }
   else
   {
#ifndef NDEBUG
      /* make sure that the edge does not exist! */
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(g->tail[e] == k);

         if( g->head[e] == j )
            break;
      }
      assert(e == EAT_LAST);
#endif

      e = EAT_LAST;
   }

   /* does edge already exist? */
   if( e != EAT_LAST )
   {
      /* correct cost */
      if( SCIPisGT(scip, g->cost[e], cost) )
      {
         g->cost[e]            = cost;
         g->cost[Edge_anti(e)] = cost;
      }
      else
      {
         e = -1;
      }
   }
   else
   {
      if( !forcedelete )
         graph_edge_del(NULL, g, eki, FALSE);

      assert(g->oeat[eki] == EAT_FREE);

      e = eki;

      g->grad[k]++;
      g->grad[j]++;

      g->cost[e]   = cost;
      g->head[e]   = j;
      g->tail[e]   = k;
      g->ieat[e]   = g->inpbeg[j];
      g->oeat[e]   = g->outbeg[k];
      g->inpbeg[j] = e;
      g->outbeg[k] = e;

      e = Edge_anti(eki);

      g->cost[e]   = cost;
      g->head[e]   = k;
      g->tail[e]   = j;
      g->ieat[e]   = g->inpbeg[k];
      g->oeat[e]   = g->outbeg[j];
      g->inpbeg[k] = e;
      g->outbeg[j] = e;
      return eki;
   }
   return e;
}

/** reinsert an edge to replace two other edges */
SCIP_RETCODE graph_edge_reinsert(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   edge,               /**< edge to reinsert */
   int                   tail,               /**< new tail */
   int                   head,               /**< new head */
   SCIP_Real             cost,               /**< new edge cost */
   int                   ancestornode,       /**< node to copy ancestors from or -1 */
   SINGLETONANS*         ancestorsBackward,  /**< backwards (edge) ancestors */
   SINGLETONANS*         ancestorsForward,   /**< forward (edge) ancestors */
   int*                  insertedge,         /**< pointer to inserted edge or -1 */
   SCIP_Bool*            conflict            /**< does the newly edge contain conflicts? (i.e. is  redundant) */
   )
{
   const int newedge = graph_edge_redirect(scip, g, edge, tail, head, cost, FALSE, TRUE);

   assert(ancestorsBackward && ancestorsForward && insertedge && conflict);
   assert(ancestornode >= 0 || ancestornode == -1);

   *conflict = FALSE;

   /* is there a new edge? */
   if( newedge >= 0 )
   {
      IDX* const ancestorsBack = ancestorsBackward->ancestors;
      IDX* const revancestorsBack = ancestorsBackward->revancestors;
      IDX* const ancestorsFor = ancestorsForward->ancestors;
      IDX* const revancestorsFor = ancestorsForward->revancestors;

      graph_edge_delHistory(scip, g, newedge);

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[newedge]), revancestorsBack, NULL) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[newedge]), ancestorsFor, NULL) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(newedge)]), ancestorsBack, NULL) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(newedge)]), revancestorsFor, NULL) );

      SCIP_CALL( graph_pseudoAncestors_appendCopySingToEdge(scip, newedge, ancestorsBackward, TRUE, g, conflict) );
      assert(!(*conflict));

      SCIP_CALL( graph_pseudoAncestors_appendCopySingToEdge(scip, newedge, ancestorsForward, TRUE, g, conflict) );

      /* ancestor node given?*/
      if( ancestornode >= 0 )
      {
         IDX* const nodeans = g->pcancestors[ancestornode];

         assert(graph_pc_isPcMw(g));

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[newedge]), nodeans, NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(newedge)]), nodeans, NULL) );

         if( !(*conflict) )
            SCIP_CALL( graph_pseudoAncestors_appendCopyNodeToEdge(scip, newedge, ancestornode, TRUE, g, conflict) );
      }
   }

   *insertedge = newedge;

   return SCIP_OKAY;
}


/** add a new edge to the graph */
void graph_edge_add(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   tail,               /**< tail of the new edge */
   int                   head,               /**< head of the new edge*/
   SCIP_Real             cost1,              /**< tail to head cost */
   SCIP_Real             cost2               /**< head to tail cost */
   )
{
   int    e;

   assert(g      != NULL);
   assert(SCIPisGE(scip, cost1, 0.0) || SCIPisEQ(scip, cost1, (double) UNKNOWN));
   assert(SCIPisGE(scip, cost2, 0.0) || SCIPisEQ(scip, cost2, (double) UNKNOWN));
   assert(tail   >= 0);
   assert(tail   <  g->knots);
   assert(head   >= 0);
   assert(head   <  g->knots);

   assert(g->esize >= g->edges + 2);

   e = g->edges;

   g->grad[head]++;
   g->grad[tail]++;

   if( cost1 != UNKNOWN )
      g->cost[e]           = cost1;
   g->tail[e]           = tail;
   g->head[e]           = head;
   g->ieat[e]           = g->inpbeg[head];
   g->oeat[e]           = g->outbeg[tail];
   g->inpbeg[head]      = e;
   g->outbeg[tail]      = e;

   e++;

   if( cost2 != UNKNOWN )
      g->cost[e]           = cost2;
   g->tail[e]           = head;
   g->head[e]           = tail;
   g->ieat[e]           = g->inpbeg[tail];
   g->oeat[e]           = g->outbeg[head];
   g->inpbeg[tail]      = e;
   g->outbeg[head]      = e;

   g->edges += 2;
}

/** add a bi-edge to the graph */
void graph_edge_addBi(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   tail,               /**< tail of the new edge */
   int                   head,               /**< head of the new edge */
   SCIP_Real             cost                /**< head to tail cost */
   )
{
   graph_edge_add(scip, g, tail, head, cost, cost);
}


/** add a new edge to a subgraph */
void graph_edge_addSubgraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orggraph,           /**< original graph */
   const int*            nodemapOrg2sub,     /**< node mapping from original to subgraph */
   int                   orgedge,            /**< original edge */
   GRAPH*                subgraph            /**< the subgraph */
   )
{
   const int e = orgedge;
   const int orgtail = orggraph->tail[e];
   const int orghead = orggraph->head[e];
   const int newtail = nodemapOrg2sub[orgtail];
   const int newhead = nodemapOrg2sub[orghead];

   assert(e >= 0 && e < orggraph->edges);

   if( graph_pc_isPcMw(orggraph) )
   {
      assert(orggraph->extended);
      graph_pc_updateSubgraphEdge(orggraph, nodemapOrg2sub, e, subgraph);
   }

   graph_edge_add(scip, subgraph, newtail, newhead, orggraph->cost[e], orggraph->cost[flipedge(e)]);
}


/** delete an edge from standard data structure */
void graph_edge_del(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   e,                  /**< the edge */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   assert(g);
   assert(e >= 0 && e < g->edges);

   if( freeancestors )
   {
      graph_edge_delHistory(scip, g, e);
   }

   /* delete first arc */
   e -= e % 2;
   assert(g->head[e] == g->tail[e + 1]);
   assert(g->tail[e] == g->head[e + 1]);

   g->grad[g->head[e]]--;
   g->grad[g->tail[e]]--;

   removeEdge(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_FREE;
   g->oeat[e] = EAT_FREE;

   /* delete second arc */
   e++;
   removeEdge(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_FREE;
   g->oeat[e] = EAT_FREE;

   assert(g->tail[e] >= 0 && g->tail[e] < g->knots);
   assert(g->head[e] >= 0 && g->head[e] < g->knots);
}

/** delete an edge from standard, DCSR (if existent) and CSR (if existent) data structures */
void graph_edge_delFull(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   e,                  /**< the edge */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   assert(scip && g);
   assert(e >= 0 && e < g->edges);

   if( g->dcsr_storage )
   {
      int csredge;
      assert(graph_valid_dcsr(g, FALSE));

      csredge = g->dcsr_storage->id2csredge[e];
      graph_dcsr_deleteEdgeBi(scip, g->dcsr_storage, csredge);

      assert(graph_valid_dcsr(g, FALSE));
   }

   if( g->csr_storage )
   {
      assert(0 && "not yet supported");
   }

   graph_edge_del(scip, g, e, freeancestors);
}

/** deletes edges marked by given array */
void graph_edge_delBlocked(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   const SCIP_Bool*      edge_deletable,     /**< marks edges to delete (of size nedges / 2) */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   const int nedges = g->edges;

   assert(g && edge_deletable);

   for( int e = 0; e < nedges / 2; e++ )
      if( edge_deletable[e] )
         graph_edge_del(scip, g, e * 2, freeancestors);
}


/** hide edge */
void graph_edge_hide(
   GRAPH*                g,                  /**< the graph */
   int                   e                   /**< the edge */
   )
{
   assert(g          != NULL);
   assert(e          >= 0);
   assert(e          <  g->edges);

   /* Immer mit der ersten von beiden Anfangen
    */
   e -= e % 2;

   assert(g->head[e] == g->tail[e + 1]);
   assert(g->tail[e] == g->head[e + 1]);

   g->grad[g->head[e]]--;
   g->grad[g->tail[e]]--;

   removeEdge(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_HIDE;
   g->oeat[e] = EAT_HIDE;

   e++;

   removeEdge(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_HIDE;
   g->oeat[e] = EAT_HIDE;
}


/** is the edge blocked? */
SCIP_Bool graph_edge_isBlocked(
   const GRAPH*          g,                  /**< the graph */
   int                   e                   /**< the edge */
   )
{
   assert(g);
   assert(e >= 0 && e < g->edges);

   if( EQ(g->cost[e], BLOCKED_MINOR) || EQ(g->cost[e], BLOCKED) )
      return TRUE;

   return FALSE;
}


/** is edge in range? */
SCIP_Bool graph_edge_isInRange(
   const GRAPH*          g,                  /**< the graph */
   int                   e                   /**< the edge */
   )
{
   assert(g);

   return (0 <= e && e < g->edges);
}


/** print information on edge */
void graph_edge_printInfo(
   const GRAPH*          g,                  /**< the graph */
   int                   e                   /**< the edge */
   )
{
   const int t = g->tail[e];
   const int h = g->head[e];

   if( graph_pc_isPcMw(g) )
      printf("e: %d   %d->%d (%d->%d) cost:=%f costrev=%f \n", e, t, h, g->term[t], g->term[h], g->cost[e], g->cost[flipedge(e)]);
   else
      printf("e: %d   %d->%d (%d->%d) cost:=%f \n", e, t, h, g->term[t], g->term[h], g->cost[e]);
}


/** print information on node */
void graph_knot_printInfo(
   const GRAPH*          g,                  /**< the graph */
   int                   k                   /**< the node */
   )
{
   assert(!graph_pc_isPcMw(g) || g->term2edge != NULL);

   if( graph_pc_isPcMw(g) )
   {
      assert(g->prize);

      if( graph_pc_knotIsNonLeafTerm(g, k) )
         printf("node %d: term=%d grad=%d prize=%f (non-leaf terminal) \n", k, g->term[k], g->grad[k], g->prize[k]);
      else if( graph_pc_knotIsFixedTerm(g, k) )
         printf("node %d: term=%d grad=%d prize=%f (fixed terminal) \n", k, g->term[k], g->grad[k], g->prize[k]);
      else if( Is_term(g->term[k]) )
         printf("node %d: term=%d grad=%d prize=%f (standard terminal) \n", k, g->term[k], g->grad[k], g->prize[k]);
      else
         printf("node %d: term=%d grad=%d prize=%f \n", k, g->term[k], g->grad[k], g->prize[k]);
   }
   else
   {
      printf("node %d: term=%d grad=%d  \n", k, g->term[k], g->grad[k]);
   }

   if( k == g->source )
      printf("...%d is the root! \n", k);

}


/** print information on graph  */
void graph_printInfo(
   const GRAPH*          g                   /**< the graph */
   )
{
   char type[64];

   switch( g->stp_type )
      {
   case 0:
      strcpy(type, "STP_SPG");
      break;
   case 1:
      strcpy(type, "STP_SAP");
      break;
   case 2:
      strcpy(type, "STP_PCSPG");
      break;
   case 3:
      strcpy(type, "STP_RPCSPG");
      break;
   case 4:
      strcpy(type, "STP_NWSPG");
      break;
   case 5:
      strcpy(type, "STP_DCSTP");
      break;
   case 6:
      strcpy(type, "STP_NWPTSPG");
      break;
   case 7:
      strcpy(type, "STP_RSMT");
      break;
   case 8:
      strcpy(type, "STP_OARSMT");
      break;
   case 9:
      strcpy(type, "STP_MWCSP");
      break;
   case 10:
      strcpy(type, "STP_DHCSTP");
      break;
   case 11:
      strcpy(type, "STP_GSTP");
      break;
   case 12:
      strcpy(type, "STP_RMWCSP");
      break;
   case 13:
      strcpy(type, "STP_BRMWCSP");
      break;
   default:
      strcpy(type, "UNKNOWN");
   }

   assert(g != NULL);
   if( graph_pc_isPcMw(g) )
   {
      printf("nodes=%d, edges=%d, terminals=%d, root=%d, type=%s, isExtended=%d\n", g->knots, g->edges, g->terms, g->source, type, g->extended);
      if( graph_pc_isPc(g) )
      {
         printf("non-leaf terminals=%d, ", graph_pc_nNonLeafTerms(g));
         printf("fixed terminals=%d, ", graph_pc_nFixedTerms(g));
         printf("proper terminals=%d \n", graph_pc_nProperPotentialTerms(g));
      }
   }
   else
   {
      printf("nodes=%d, edges=%d, terminals=%d, root=%d, type=%s \n", g->knots, g->edges, g->terms, g->source, type);
   }

   if( g->stp_type == STP_BRMWCSP )
      printf("budget=%f \n", g->budget);
}


/** get number of nodes */
int graph_get_nNodes(
   const GRAPH*    graph               /**< the graph */
)
{
   assert(graph);
   assert(graph->knots >= 1);

   return graph->knots;
}


/** get number of edges */
int graph_get_nEdges(
   const GRAPH*    graph               /**< the graph */
)
{
   assert(graph);
   assert(graph->edges >= 0);

   return graph->edges;
}


/** get number of terminals */
int graph_get_nTerms(
   const GRAPH*    graph               /**< the graph */
)
{
   assert(graph);

   return graph->terms;
}


/** get (real) number of nodes, edges, terminals */
void graph_get_nVET(
   const GRAPH*    graph,              /**< the graph */
   int*            nnodes,             /**< number of nodes */
   int*            nedges,             /**< number of edges */
   int*            nterms              /**< number of terminals */
   )
{
   int v = 0;
   int e = 0;
   int t = 0;
   int vorg;

   assert(graph != NULL);

   vorg = graph->knots;

   for( int k = 0; k < vorg; k++ )
   {
      if( graph->grad[k] > 0 )
      {
         v++;
         e += graph->grad[k];
         if( Is_term(graph->term[k]) )
            t++;
      }
   }

   if( nnodes )
      *nnodes = v;

   if( nedges )
      *nedges = e;

   if( nterms )
      *nterms = t;

   return;
}


/** get (real) average degree of graph */
SCIP_Real graph_get_avgDeg(
   const GRAPH*    graph               /**< the graph */
   )
{
   int v = 0;
   int e = 0;
   const int vorg = graph->knots;

   assert(graph);

   for( int k = 0; k < vorg; k++ )
   {
      if( graph->grad[k] > 0 )
      {
         v++;
         e += graph->grad[k];
      }
   }

   if( v == 0 )
      return 0.0;

   return ((double) e / v );
}


/* get compressed sparse row arrays representing current graph */
void graph_get_csr(
   const GRAPH*          g,                  /**< the graph */
   int* RESTRICT         edgearr,            /**< original edge array [0,...,nedges - 1] */
   int* RESTRICT         tailarr,            /**< tail of csr edge [0,...,nedges - 1]  */
   int* RESTRICT         start,              /**< start array [0,...,nnodes] */
   int*                  nnewedges           /**< pointer to store number of new edges */
      )
{
   int i = 0;
   const int nnodes = g->knots;

   assert(g != NULL);
   assert(tailarr != NULL);
   assert(edgearr != NULL);
   assert(start != NULL);

   for( int k = 0; k < nnodes; k++ )
   {
      start[k] = i;
      for( int e = g->inpbeg[k]; e >= 0; e = g->ieat[e] )
      {
         edgearr[i] = e;
         tailarr[i++] = g->tail[e] + 1;
      }
   }

   *nnewedges = i;
   start[nnodes] = i;
}


/** get edge costs */
void graph_get_edgeCosts(
   const GRAPH*          graph,              /**< the graph */
   SCIP_Real* RESTRICT   cost,               /**< reduced edge costs */
   SCIP_Real* RESTRICT   costrev             /**< reduced reverse edge costs */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real* const gcost = graph->cost;

   assert(cost && costrev);

   for( int e = 0; e < nedges; e++ )
   {
      cost[e] = gcost[e];
      costrev[e] = gcost[flipedge(e)];
      assert(GE(cost[e], 0.0));
   }
}


/* modifies 'isterm' to mark whether node is a terminal (or proper terminal for PC) */
void graph_get_isTerm(
   const GRAPH*          g,                  /**< the graph */
   SCIP_Bool*            isterm              /**< marks whether node is a terminal (or proper terminal for PC) */
)
{
   const int nnodes = g->knots;
   const SCIP_Bool pcmw = graph_pc_isPcMw(g);

   assert(g && isterm);
   assert(!pcmw || !g->extended);

   for( int i = 0; i < nnodes; i++ )
   {
      isterm[i] = FALSE;

      if( Is_term(g->term[i]) )
      {
         if( pcmw && graph_pc_termIsNonLeafTerm(g, i) )
           continue;

         isterm[i] = TRUE;
      }
   }
}

/** gets edge conflicts */
SCIP_RETCODE graph_get_edgeConflicts(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< the graph */
)
{
   int* childcount;
   int* pseudonodecount;
   int nconflicts;
   int npseudoconflicts;
   int npseudofixed;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const int nedgesorg = MAX(g->orgedges, g->edges);

   assert(scip != NULL && g != NULL);
   assert(g->ancestors != NULL);
   assert(nedgesorg % 2 == 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &childcount, nedgesorg / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pseudonodecount, nnodes) );

   for( int e = 0; e < nedgesorg / 2; e++ )
      childcount[e] = 0;

   for( int e = 0; e < nnodes; e++ )
      pseudonodecount[e] = 0;

   npseudofixed = graph_get_nFixpseudonodes(scip, g);
   nconflicts = 0;
   npseudoconflicts = 0;

   for( int e = 0; e < nedges; e += 2 )
   {
      SCIP_Bool hasConflict = FALSE;
      SCIP_Bool hasPseudoConflict = FALSE;

      const int nPseudoAncestors = graph_edge_nPseudoAncestors(g, e);
      const int* pseudoAncestors = graph_edge_getPseudoAncestors(g, e);
      for( IDX* curr = g->ancestors[e]; curr != NULL; curr = curr->parent )
      {
         assert(curr->index >= 0 && curr->index / 2 < nedgesorg / 2);
         if( childcount[curr->index / 2] > 0 )
            hasConflict = TRUE;

         childcount[curr->index / 2]++;
      }

      for( int k = 0; k < nPseudoAncestors; k++ )
      {
         const int a = pseudoAncestors[k];
         assert(a >= 0 && a < nnodes);

         if( pseudonodecount[a] > 0 )
            hasPseudoConflict = TRUE;

         pseudonodecount[a]++;
      }

      if( hasConflict )
      {
         nconflicts++;
         assert(hasPseudoConflict);
      }

      if( hasPseudoConflict )
         npseudoconflicts++;

   }

   if( graph_pc_isPcMw(g) )
   {
      assert(g->extended);

      for( int e = 0; e < nnodes; e++ )
      {
         if( !Is_term(g->term[e]) )
         {
            SCIP_Bool hasConflict = FALSE;
            SCIP_Bool hasPseudoConflict = FALSE;

            const int nPseudoAncestors = graph_knot_nPseudoAncestors(g, e);
            const int* pseudoAncestors = graph_knot_getPseudoAncestors(g, e);
            for( IDX* curr = g->pcancestors[e]; curr != NULL; curr = curr->parent )
            {
               assert(curr->index >= 0 && curr->index / 2 < nedgesorg / 2);
               if( childcount[curr->index / 2] > 0 )
                  hasConflict = TRUE;

               childcount[curr->index / 2]++;
            }

            for( int k = 0; k < nPseudoAncestors; k++ )
            {
               const int a = pseudoAncestors[k];
               assert(a >= 0 && a < nnodes);

               if( pseudonodecount[a] > 0 )
                  hasPseudoConflict = TRUE;

               pseudonodecount[a]++;
            }

            if( hasConflict )
            {
               nconflicts++;
               assert(hasPseudoConflict);
            }

            if( hasPseudoConflict )
               npseudoconflicts++;
         }
      }
   }

   printf("number of edge conflicts %d \n", nconflicts);
   printf("number of pseudo-ancestor conflicts %d \n", npseudoconflicts);
   printf("number of fixed pseudo-ancestors %d \n", npseudofixed);

   SCIPfreeBufferArray(scip, &pseudonodecount);
   SCIPfreeBufferArray(scip, &childcount);

   return SCIP_OKAY;
}


/** initialize graph */
SCIP_RETCODE graph_init(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               g,                  /**< new graph */
   int                   ksize,              /**< slots for nodes */
   int                   esize,              /**< slots for edges */
   int                   layers              /**< number of layers (only needed for packing, otherwise 1) */
   )
{
   GRAPH* p;

   assert(ksize > 0);
   assert(ksize < INT_MAX);
   assert(esize >= 0);
   assert(esize < INT_MAX);
   assert(layers > 0);
   assert(layers < SHRT_MAX);

   SCIP_CALL( SCIPallocMemory(scip, g) );
   p = *g;
   assert(p != NULL);

   p->ancestors = NULL;
   p->pcancestors = NULL;
   p->fixedcomponents = NULL;
   p->orgtail = NULL;
   p->orghead = NULL;
   p->rootedgeprevs = NULL;
   p->norgmodelknots = 0;
   p->norgmodeledges = 0;
   p->ksize  = ksize;
   p->orgknots = 0;
   p->orgedges = 0;
   p->knots  = 0;
   p->terms  = 0;
   p->orgsource = UNKNOWN;
   p->stp_type = UNKNOWN;
   p->layers = layers;
   p->hoplimit = UNKNOWN;
   p->extended = FALSE;
   p->source = -1;
   p->is_packed = FALSE;
   p->cost_org_pc = NULL;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->term), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mark), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->grad), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->inpbeg), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->outbeg), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->cost), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->tail), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->head), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->ieat), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->oeat), esize) );

   p->esize = esize;
   p->edges = 0;
   p->prize  = NULL;
   p->maxdeg = NULL;
   p->grid_coordinates = NULL;
   p->grid_ncoords = NULL;
   p->mincut_dist = NULL;
   p->mincut_head = NULL;
   p->mincut_numb = NULL;
   p->mincut_prev = NULL;
   p->mincut_next = NULL;
   p->mincut_temp = NULL;
   p->mincut_e = NULL;
   p->mincut_x = NULL;
   p->mincut_r = NULL;
   p->path_heap = NULL;
   p->path_state = NULL;
   p->term2edge = NULL;
   p->budget = -1.0;
   p->costbudget = NULL;
   p->csr_storage = NULL;
   p->dcsr_storage = NULL;
   p->pseudoancestors = NULL;

   return SCIP_OKAY;
}

/** initialize data structures required to keep track of reductions */
SCIP_RETCODE graph_init_history(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< graph */
   )
{
   IDX** ancestors;          /* ancestor lists array (over all edges) */
   IDX** pcancestors;        /* ancestor lists array (over all nodes) */
   int* tail;                /* tail of all edges  */
   int* head;                /* head of all edges  */
   int* orgtail;             /* (original) tail of all original edges  */
   int* orghead;             /* (original) head of all original edges  */
   const int nedges = graph->edges;
   SCIP_Bool pcmw;

   assert(scip != NULL);
   assert(graph != NULL);

   pcmw = graph_pc_isPcMw(graph);

   SCIP_CALL( graph_init_pseudoAncestors(scip, graph) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->orgtail), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->orghead), nedges) );

   tail = graph->tail;
   head = graph->head;
   orgtail = graph->orgtail;
   orghead = graph->orghead;

   for( int e = 0; e < nedges; e++ )
   {
      orgtail[e] = tail[e];
      orghead[e] = head[e];
   }

   if( pcmw )
   {
      const int nnodes = graph->knots;

      SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->pcancestors), nnodes) );

      pcancestors = graph->pcancestors;

      for( int k = 0; k < nnodes; k++ )
         pcancestors[k] = NULL;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->ancestors), nedges) );

   ancestors = graph->ancestors;

   for( int e = 0; e < nedges; e++ )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &(ancestors[e])) ); /*lint !e866*/
      (ancestors)[e]->index = e;
      (ancestors)[e]->parent = NULL;
   }

   return SCIP_OKAY;
}


/** builds complete graph (Kn) with unit edge weights and no terminals */
SCIP_RETCODE graph_buildCompleteGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               g,                  /**< new graph */
   int                   nnodes              /**< number of nodes */
   )
{
   const int nedges = nnodes * (nnodes - 1);
   GRAPH* graph;

   assert(scip && g);
   assert(nnodes >= 1);

   SCIP_CALL( graph_init(scip, g, nnodes, nedges, 1) );

   graph = *g;

   for( int k = 0; k < nnodes; k++ )
      graph_knot_add(graph, -1);

   for( int k = 0; k < nnodes - 1; k++ )
      for( int k2 = nnodes - 1; k2 >= k + 1; k2-- )
         graph_edge_add(scip, graph, k, k2, 1.0, 1.0);

   return SCIP_OKAY;
}


/** enlarge given graph */
SCIP_RETCODE graph_resize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph to be resized */
   int                   ksize,              /**< new node slots */
   int                   esize,              /**< new edge slots */
   int                   layers              /**< layers (set to -1 by default) */
   )
{
   assert(scip      != NULL);
   assert(g      != NULL);
   assert((ksize  < 0) || (ksize  >= g->knots));
   assert((esize  < 0) || (esize  >= g->edges));
   assert((layers < 0) || (layers >= g->layers));

   if( (layers > 0) && (layers != g->layers) )
      g->layers = layers;

   if( (ksize > 0) && (ksize != g->ksize) )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->term), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->mark), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->grad), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->inpbeg), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->outbeg), ksize) );

      if( g->prize != NULL)
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->prize), ksize) );

      if( g->costbudget != NULL)
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->costbudget), ksize) );

      if( g->term2edge != NULL)
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->term2edge), ksize) );

      g->ksize  = ksize;
   }
   if( (esize > 0) && (esize != g->esize) )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->cost), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->tail), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->head), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->ieat), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->oeat), esize) );

      g->esize = esize;
   }

   return SCIP_OKAY;
}


/** free the graph */
void graph_free(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph to be freed */
   SCIP_Bool             final               /**< delete ancestor data structures? */
   )
{
   GRAPH* p;

   assert(scip != NULL);
   assert(graph != NULL);

   p = *graph;
   assert(p != NULL);

   graph_free_history(scip, p);

   if( final )
      graph_free_historyDeep(scip, p);

   SCIPfreeMemoryArrayNull(scip, &(p->prize));
   SCIPfreeMemoryArrayNull(scip, &(p->costbudget));
   SCIPfreeMemoryArrayNull(scip, &(p->term2edge));

   if( p->csr_storage != NULL )
      graph_free_csr(scip, p);

   if( p->dcsr_storage != NULL )
      graph_free_dcsr(scip, p);

   if( p->stp_type == STP_DCSTP )
   {
      SCIPfreeMemoryArray(scip, &(p->maxdeg));
   }
   else if( p->stp_type == STP_RSMT )
   {
      if( p->grid_coordinates != NULL )
      {
         assert(p->grid_coordinates != NULL);
         for( int i = p->grid_dim - 1; i >= 0;  i-- )
            SCIPfreeMemoryArray(scip, &(p->grid_coordinates[i]));

         SCIPfreeMemoryArray(scip, &(p->grid_coordinates));
      }

      if( p->grid_ncoords != NULL )
         SCIPfreeMemoryArray(scip, &(p->grid_ncoords));
   }

   SCIPfreeMemoryArray(scip, &(p->oeat));
   SCIPfreeMemoryArray(scip, &(p->ieat));
   SCIPfreeMemoryArray(scip, &(p->head));
   SCIPfreeMemoryArray(scip, &(p->tail));
   SCIPfreeMemoryArray(scip, &(p->cost));
   SCIPfreeMemoryArray(scip, &(p->outbeg));
   SCIPfreeMemoryArray(scip, &(p->inpbeg));
   SCIPfreeMemoryArray(scip, &(p->grad));
   SCIPfreeMemoryArray(scip, &(p->mark));
   SCIPfreeMemoryArray(scip, &(p->term));
   SCIPfreeMemoryArrayNull(scip, &(p->rootedgeprevs));
   SCIPfreeMemoryArrayNull(scip, &(p->cost_org_pc));

   SCIPfreeMemory(scip, graph);
}


/** frees the history */
void graph_free_history(
   SCIP*                 scip,               /**< SCIP data */
   GRAPH*                p                   /**< graph data */
   )
{
   if( p->pseudoancestors != NULL )
      graph_free_pseudoAncestors(scip, p);

   if( p->ancestors != NULL )
   {
      const int nedges = p->edges;

      for( int e = nedges - 1; e >= 0; e-- )
      {
         IDX* curr = p->ancestors[e];
         while( curr != NULL )
         {
            p->ancestors[e] = curr->parent;
            SCIPfreeBlockMemory(scip, &(curr));
            curr = p->ancestors[e];
         }
      }
      SCIPfreeMemoryArray(scip, &(p->ancestors));
   }

   assert(p->pseudoancestors == NULL && p->ancestors == NULL);
}


/** free the history of an edge */
void graph_edge_delHistory(
   SCIP*                 scip,               /**< SCIP data */
   GRAPH*                g,                  /**< graph data */
   int                   edge                /**< edge for which to free the history */
   )
{
   assert(scip && g);
   assert(edge >= 0 && edge < g->edges);

   SCIPintListNodeFree(scip, &(g->ancestors[edge]));
   SCIPintListNodeFree(scip, &(g->ancestors[flipedge_Uint(edge)]));

   graph_edge_delPseudoAncestors(scip, edge, g);
}

/** frees the deep history */
void graph_free_historyDeep(
   SCIP*                 scip,               /**< SCIP data */
   GRAPH*                p                   /**< graph data */
   )
{
   IDX* curr;

   assert(scip != NULL);
   assert(p != NULL);
   assert(p->path_heap == NULL);
   assert(p->path_state == NULL);

   if( p->pcancestors != NULL )
   {
      for( int e = p->norgmodelknots - 1; e >= 0; e-- )
      {
         curr = p->pcancestors[e];
         while( curr != NULL )
         {
            p->pcancestors[e] = curr->parent;
            SCIPfreeBlockMemory(scip, &(curr));
            curr = p->pcancestors[e];
         }
      }
      SCIPfreeMemoryArray(scip, &(p->pcancestors));
   }

   if( p->orgtail != NULL )
   {
      assert(p->orghead != NULL);

      SCIPfreeMemoryArray(scip, &(p->orghead));
      SCIPfreeMemoryArray(scip, &(p->orgtail));
   }

   if( p->fixedcomponents )
      graph_free_fixed(scip, p);
}

/** copy the data of the graph */
SCIP_RETCODE graph_copy_data(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orgraph,            /**< original graph */
   GRAPH*                copygraph           /**< graph to be copied to */
   )
{
   GRAPH* g_copy = copygraph;
   const GRAPH* g_org = orgraph;
   const int ksize = g_org->ksize;
   const int esize = g_org->esize;

   assert(scip != NULL);
   assert(orgraph != NULL);
   assert(copygraph != NULL);
   assert(ksize == g_copy->ksize && ksize > 0);
   assert(esize == g_copy->esize && esize >= 0);
   assert(g_org->dcsr_storage == NULL && g_org->csr_storage == NULL);

   g_copy->norgmodeledges = g_org->norgmodeledges;
   g_copy->norgmodelknots = g_org->norgmodelknots;
   g_copy->knots = g_org->knots;
   g_copy->terms = g_org->terms;
   g_copy->edges = g_org->edges;
   g_copy->source = g_org->source;
   g_copy->orgsource = g_org->orgsource;
   g_copy->orgedges = g_org->orgedges;
   g_copy->orgknots = g_org->orgknots;
   g_copy->grid_dim = g_org->grid_dim;
   g_copy->stp_type = g_org->stp_type;
   g_copy->hoplimit = g_org->hoplimit;
   g_copy->extended = g_org->extended;
   g_copy->budget = g_org->budget;
   g_copy->is_packed = g_org->is_packed;

   BMScopyMemoryArray(g_copy->term, g_org->term, ksize);
   BMScopyMemoryArray(g_copy->mark, g_org->mark, ksize);
   BMScopyMemoryArray(g_copy->grad, g_org->grad, ksize);
   BMScopyMemoryArray(g_copy->inpbeg, g_org->inpbeg, ksize);
   BMScopyMemoryArray(g_copy->outbeg, g_org->outbeg, ksize);
   BMScopyMemoryArray(g_copy->cost, g_org->cost, esize);
   BMScopyMemoryArray(g_copy->tail, g_org->tail, esize);
   BMScopyMemoryArray(g_copy->head, g_org->head, esize);
   BMScopyMemoryArray(g_copy->ieat, g_org->ieat, esize);
   BMScopyMemoryArray(g_copy->oeat, g_org->oeat, esize);

   if( graph_pc_isPcMw(g_copy) )
   {
      const SCIP_Bool brpcmw = (g_copy->stp_type == STP_BRMWCSP);

      assert(g_copy->extended && g_org->extended);

      if( g_copy->costbudget != NULL )
      {
         assert(brpcmw);
         SCIPfreeMemoryArray(scip, &(g_copy->costbudget));
      }

      if( brpcmw )
      {
         SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->costbudget), g_copy->knots));
         BMScopyMemoryArray(g_copy->costbudget, g_org->costbudget, g_copy->knots);
      }

      SCIPfreeMemoryArrayNull(scip, &(g_copy->prize));
      SCIPfreeMemoryArrayNull(scip, &(g_copy->term2edge));
      SCIPfreeMemoryArrayNull(scip, &(g_copy->cost_org_pc));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->prize), g_copy->knots));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->term2edge), g_copy->knots));

      if( graph_pc_isPc(g_org) )
      {
         assert(g_org->cost_org_pc != NULL);
         SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->cost_org_pc), g_copy->edges));
         BMScopyMemoryArray(g_copy->cost_org_pc, g_org->cost_org_pc, g_copy->edges);
      }
      else
      {
         g_copy->cost_org_pc = NULL;
      }

      assert(g_org->prize != NULL);
      BMScopyMemoryArray(g_copy->prize, g_org->prize, g_copy->knots);

#ifndef NDEBUG
      for( int k = 0; k < g_copy->knots; k++ )
      {
         if( Is_term(g_org->term[k]) && (!graph_pc_isRootedPcMw(g_copy) || !graph_pc_knotIsFixedTerm(g_org, k)) )
            assert(EQ(g_copy->prize[k], 0.0));
      }
#endif

      assert(g_org->term2edge != NULL);
      BMScopyMemoryArray(g_copy->term2edge, g_org->term2edge, g_copy->knots);
   }
   else if( g_copy->stp_type == STP_DCSTP )
   {
      assert(g_org->maxdeg != NULL);

      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->maxdeg), g_copy->knots));

      for( int k = 0; k < g_copy->knots; k++ )
         g_copy->maxdeg[k] = g_org->maxdeg[k];
   }
   else if( g_org->stp_type == STP_RSMT )
   {
      assert(g_org->grid_ncoords != NULL);
      assert(g_org->grid_coordinates != NULL);

      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->grid_coordinates), g_org->grid_dim));

      BMScopyMemoryArray(g_copy->grid_coordinates, g_org->grid_coordinates, g_org->grid_dim);
      for( int k = 0; k < g_org->grid_dim; k++ )
      {
         SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->grid_coordinates[k]), g_org->terms)); /*lint !e866*/
         BMScopyMemoryArray(g_copy->grid_coordinates[k], g_org->grid_coordinates[k], g_org->terms); /*lint !e866*/
      }
      SCIP_CALL(SCIPallocMemoryArray(scip, &(g_copy->grid_ncoords), g_org->grid_dim));

      BMScopyMemoryArray(g_copy->grid_ncoords, g_org->grid_ncoords, g_org->grid_dim);
   }

   assert(graph_valid(scip, g_copy));

   return SCIP_OKAY;
}

/** copy the graph */
SCIP_RETCODE graph_copy(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orgraph,            /**< original graph */
   GRAPH**               copygraph           /**< graph to be created */
   )
{
   const GRAPH* p = orgraph;
   assert(p != NULL);

   SCIP_CALL( graph_init(scip, copygraph, p->ksize, p->esize, p->layers) );

   SCIP_CALL( graph_copy_data(scip, orgraph, *copygraph) );

   return SCIP_OKAY;
}


/** marks the current graph */
void graph_mark(
   GRAPH*                g                   /**< the graph */
   )
{
   const int nnodes = graph_get_nNodes(g);
   const int root = g->source;

   for( int k = 0; k < nnodes; k++ )
      g->mark[k] = (g->grad[k] > 0);

   g->mark[root] = TRUE;

   if( graph_pc_isPcMw(g) && !g->extended )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( Is_pseudoTerm(g->term[k]) )
            g->mark[k] = FALSE;
         else if( Is_term(g->term[k]) )
            g->mark[k] = TRUE;
      }

      if( graph_pc_isRootedPcMw(g) )
         g->mark[root] = TRUE;
      else
         g->mark[root] = FALSE;
   }

   assert(graph_isMarked(g));
}


/** is the current graph properly marked? */
SCIP_Bool graph_isMarked(
   const GRAPH*          g                   /**< the graph */
   )
{
   assert(g);

   if( graph_pc_isPcMw(g) && !g->extended )
   {
      const int nnodes = g->knots;
      const int root = g->source;
      const int rooted = graph_pc_isRootedPcMw(g);

      for( int k = 0; k < nnodes; k++ )
      {
         if( Is_pseudoTerm(g->term[k]) || (!rooted && k == root) )
         {
            assert(g->grad[k] == 2 || k == root);

            if( g->mark[k] )
            {
               graph_knot_printInfo(g, k);
               printf("pseudo-terminal %d is marked \n", k);
               return FALSE;
            }
         }
         else
         {
            if( g->mark[k] != (g->grad[k] > 0) )
            {
               if( k == root && g->mark[k] )
                  continue;

               if( Is_term(g->term[k]) && g->mark[k] )
                  continue;

               graph_knot_printInfo(g, k);
               printf("node %d: mark=%d grad=%d \n", k, g->mark[k], g->grad[k]);
               return FALSE;
            }
         }
      }

      if( rooted )
      {
         if( !g->mark[root] )
         {
            printf("root not marked \n");
            return FALSE;
         }
      }
   }
   else
   {
      const int nnodes = g->knots;

      for( int k = 0; k < nnodes; k++ )
      {
         if( g->mark[k] != (g->grad[k] > 0) )
         {
            if( k == g->source && g->mark[k] )
               continue;

            graph_knot_printInfo(g, k);
            printf("node %d: mark=%d grad=%d \n", k, g->mark[k], g->grad[k]);
            return FALSE;
         }
      }
   }

   return TRUE;
}



void graph_show(
   const GRAPH*          p                   /**< the graph */
   )
{
   int i;

   assert(p != NULL);

   for(i = 0; i < p->knots; i++)
      if (p->grad[i] > 0)
         (void)printf("Knot %d, term=%d, grad=%d, inpbeg=%d, outbeg=%d\n",
            i, p->term[i], p->grad[i], p->inpbeg[i], p->outbeg[i]);

   (void)fputc('\n', stdout);

   for(i = 0; i < p->edges; i++)
      if (p->ieat[i] != EAT_FREE)
         (void)printf("Edge %d, cost=%g, tail=%d, head=%d, ieat=%d, oeat=%d\n",
            i, p->cost[i], p->tail[i], p->head[i], p->ieat[i], p->oeat[i]);

   (void)fputc('\n', stdout);
}


/** reinsert all hidden edges */
void graph_uncover(
   GRAPH*                g                   /**< the graph */
   )
{/*lint --e{850}*/
   int head;
   int tail;
   int e;

   assert(g      != NULL);

   for( e = 0; e < g->edges; e++ )
   {
      if( g->ieat[e] == EAT_HIDE )
      {
         assert(e % 2 == 0);
         assert(g->oeat[e] == EAT_HIDE);

         head            = g->head[e];
         tail            = g->tail[e];

         g->grad[head]++;
         g->grad[tail]++;

         g->ieat[e]      = g->inpbeg[head];
         g->oeat[e]      = g->outbeg[tail];
         g->inpbeg[head] = e;
         g->outbeg[tail] = e;

         e++;

         assert(g->ieat[e] == EAT_HIDE);
         assert(g->oeat[e] == EAT_HIDE);
         assert(g->head[e] == tail);
         assert(g->tail[e] == head);

         head            = g->head[e];
         tail            = g->tail[e];
         g->ieat[e]      = g->inpbeg[head];
         g->oeat[e]      = g->outbeg[tail];
         g->inpbeg[head] = e;
         g->outbeg[tail] = e;
      }
   }
}


/** Packs the graph, i.e. builds a new graph that discards deleted edges and nodes.
 *  The original graph is deleted. */
SCIP_RETCODE graph_pack(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Real*            offset,             /**< pointer to add offset from non-leaf terminals to (only PC) */
   SCIP_Bool             verbose             /**< verbose? */
   )
{
   GRAPH *g_old;
   GRAPH *g_new;
   int* old2newNode;
   const int oldnnodes = graph_get_nNodes(graph);
   const int oldnedges = graph_get_nEdges(graph);
   int nnodes = 0;
   int nedges = 0;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
   SCIP_Bool graphHasVanished = FALSE;

   assert(scip && newgraph);
   assert(graph_valid(scip, graph));
   assert(!graph->is_packed);

   g_old = graph;

   SCIP_CALL( SCIPallocBufferArray(scip, &old2newNode, oldnnodes) );

   if( verbose )
      printf("Reduced graph: ");

   /* build node mapping */
   for( int i = 0; i < oldnnodes; i++ )
   {
      if( g_old->grad[i] > 0 )
         old2newNode[i] = nnodes++;
      else
         old2newNode[i] = -1;
   }

   /* graph vanished? */
   if( nnodes == 0 )
   {
      SCIPfreeBufferArray(scip, &old2newNode);
      old2newNode = NULL;
      if( verbose )
         printf(" graph vanished!\n");

      nnodes = 1;
      graphHasVanished = TRUE;

      if( pcmw )
      {
         SCIP_CALL( packPcMwVanished(scip, g_old, g_new) );
      }
   }

   /* count surviving edges */
   for( int i = 0; i < oldnedges; i++ )
   {
      if( g_old->oeat[i] != EAT_FREE )
      {
         assert(g_old->ieat[i] != EAT_FREE);
         nedges++;
      }
   }

   assert(nnodes > 1 || nedges == 0);
   SCIP_CALL( graph_init(scip, newgraph, nnodes, nedges, g_old->layers) );
   g_new = *newgraph;
   g_new->norgmodelknots = g_old->norgmodelknots;
   g_new->norgmodeledges = g_old->norgmodeledges;
   g_new->orgsource = g_old->orgsource;
   g_new->orgtail = g_old->orgtail;
   g_new->orghead = g_old->orghead;
   g_new->orgknots = g_old->knots;
   g_new->orgedges = g_old->edges;
   g_new->stp_type = g_old->stp_type;
   g_new->maxdeg = g_old->maxdeg;
   g_new->grid_dim = g_old->grid_dim;
   g_new->grid_ncoords = g_old->grid_ncoords;
   g_new->grid_coordinates = g_old->grid_coordinates;
   g_new->pcancestors = g_old->pcancestors;
   g_new->fixedcomponents = g_old->fixedcomponents;
   g_new->hoplimit = g_old->hoplimit;
   g_new->extended = g_old->extended;
   g_new->budget = g_old->budget;
   g_new->is_packed = TRUE;

   if( graphHasVanished )
   {
      assert(!old2newNode);

      g_new->ancestors = NULL;
      graph_free(scip, &g_old, FALSE);

      if( g_new->stp_type == STP_RSMT )
      {
         g_new->grid_ncoords = NULL;
         g_new->grid_coordinates = NULL;
      }

      graph_knot_add(g_new, 0);
      g_new->source = 0;

      return SCIP_OKAY;
   }

   assert(nnodes >= 2 && nedges >= 1);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_new->ancestors), nedges) );

   if( pcmw )
   {
      SCIP_CALL( packPcMwInit(scip, nnodes, g_old, g_new) );
   }

   /* add nodes (of positive degree) to new graph */
   packNodes(scip, g_old, g_new);

   /* add root */
   assert(Is_term(g_new->term[old2newNode[g_old->source]]));

   g_new->source = old2newNode[g_old->source];

   assert(!graph_pc_isRootedPcMw(graph) || FARAWAY == g_new->prize[g_new->source]);

   SCIP_CALL( packEdges(scip, old2newNode, g_old, nnodes, g_new) );

   SCIPfreeBufferArray(scip, &old2newNode);

   SCIP_CALL( graph_pc_finalizeSubgraph(scip, g_new) );

   if( graph_pc_isPc(g_new) )
   {
      assert(offset);
      *offset += graph_pc_getNonLeafTermOffset(scip, g_new);
   }

   if( g_old->path_heap != NULL )
      graph_path_exit(scip, g_old);

   g_old->stp_type = UNKNOWN;
   graph_free(scip, &g_old, FALSE);

   assert(graph_valid(scip, g_new));

   if( verbose )
      printf("Nodes: %d  Edges: %d  Terminals: %d\n", g_new->knots, g_new->edges, g_new->terms);

   return SCIP_OKAY;
}


/** traverses the graph from vertex 'start' and marks all reached nodes */
SCIP_RETCODE graph_trail_arr(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the new graph */
   int                   start,              /**< node to start from */
   SCIP_Bool*            nodevisited         /**< marks which node has been visited */
   )
{
   assert(g && scip && nodevisited);
   assert(start >= 0 && start < g->knots);

   trail(scip, g, start, FALSE, nodevisited);

   return SCIP_OKAY;
}


/** traverses the graph and marks all reached nodes, does not take edge of cost >= FARAWAY */
SCIP_RETCODE graph_trail_costAware(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the new graph */
   int                   start,              /**< node to start from */
   SCIP_Bool*            nodevisited         /**< marks which node has been visited */
   )
{
   assert(g && scip && nodevisited);
   assert(start >= 0 && start < g->knots);

   trail(scip, g, start, TRUE, nodevisited);

   return SCIP_OKAY;
}


/** checks whether all terminals are reachable from root */
SCIP_RETCODE graph_termsReachable(
   SCIP*                 scip,               /**< scip struct */
   const GRAPH*          g,                  /**< the new graph */
   SCIP_Bool*            reachable           /**< are they reachable? */
   )
{
   const int nnodes = g->knots;
   SCIP_Bool* nodevisited;

   assert(g != NULL);
   assert(reachable != NULL);

   *reachable = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodevisited, nnodes) );
   SCIP_CALL( graph_trail_arr(scip, g, g->source, nodevisited) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && !nodevisited[k] )
      {
         *reachable = FALSE;
         break;
      }
   }

   SCIPfreeBufferArray(scip, &nodevisited);

   return SCIP_OKAY;
}


/*
 * distinguishes a terminal as the root; with centertype
 *      = CENTER_OK  : Do nothing
 *      = CENTER_DEG : find maximum degree
 *      = CENTER_SUM : find the minimum distance sum
 *      = CENTER_MIN : find the minimum largest distance
 *      = CENTER_ALL : find the minimum distance sum to all knots
 */
SCIP_RETCODE graph_findCentralTerminal(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   centertype,         /**< type of root selection */
   int*                  central_term        /**< pointer to store the selected (terminal) vertex */
   )
{
   PATH*   path;
   double* cost;
   int     i;
   int     k;
   int     center  = -1;
   int     degree  = 0;
   double  sum;
   double  max;
   double  minimum = FARAWAY;
   double  maximum = 0.0;
   double  oldval  = 0.0;

   assert(g         != NULL);
   assert(g->layers == 1);
   assert(centertype == STP_CENTER_OK || centertype == STP_CENTER_DEG ||
          centertype == STP_CENTER_SUM  || centertype == STP_CENTER_MIN || centertype == STP_CENTER_ALL );

   *central_term = g->source;

   if( centertype == STP_CENTER_OK || g->grad[g->source] == 0)
   {
      assert(Is_term(g->term[*central_term]));

      return SCIP_OKAY;
   }

   /* Find knot of maximum degree.
    */
   if( centertype == STP_CENTER_DEG )
   {
      degree = 0;

      for( i = 0; i < g->knots; i++ )
      {
         if( g->stp_type == STP_NWPTSPG && graph_nw_knotIsLeaf(g, i) )
            continue;

         if( Is_term(g->term[i]) && (g->grad[i] > degree) )
         {
            degree = g->grad[i];
            center = i;
         }
      }

      assert(degree > 0);
      assert(Is_term(g->term[center]));

      *central_term = center;

      return SCIP_OKAY;
   }

   /* For the other methods we need the shortest paths */
   SCIP_CALL( SCIPallocBufferArray(scip, &path, g->knots) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, g->edges) );

   assert(path != NULL);
   assert(cost != NULL);

   for( i = 0; i < g->knots; i++ )
      g->mark[i] = TRUE;

   for( i = 0; i < g->edges; i++ )
      cost[i] = 1.0;

   for( i = 0; i < g->knots; i++ )
   {
      if (!Is_term(g->term[i]))
         continue;

      if( g->stp_type == STP_NWPTSPG && graph_nw_knotIsLeaf(g, i) )
         continue;

      graph_path_exec(scip, g, FSP_MODE, i, cost, path);

      sum = 0.0;
      max = 0.0;

      for( k = 0; k < g->knots; k++ )
      {
         assert((path[k].edge >= 0) || (k == i));
         assert((path[k].edge >= 0) || (path[k].dist == 0));

         if( Is_term(g->term[k]) || (centertype == STP_CENTER_ALL) )
         {
            sum += path[k].dist;

            if( path[k].dist > max )
               max = path[k].dist;
         }
      }

      if( (centertype == STP_CENTER_SUM) || (centertype == STP_CENTER_ALL) )
      {
         if( sum < minimum )
         {
            minimum = sum;
            center  = i;
         }
         if( sum > maximum )
            maximum = sum;

         if( i == g->source )
            oldval = sum;
      }
      else
      {
         assert(centertype == STP_CENTER_MIN);

         /* If the maximum distance to terminal ist shorter or if
          * it is of the same length but the degree of the knot is
          * higher, we change the center.
          */
         if( SCIPisLT(scip, max, minimum) || (SCIPisEQ(scip, max, minimum) && (g->grad[i] > degree)) )
         {
            minimum = max;
            center  = i;
            degree  = g->grad[i];
         }
         if( max > maximum )
            maximum = max;

         if( i == g->source )
            oldval = max;
      }
   }
   assert(center >= 0);
   assert(Is_term(g->term[center]));

   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &path);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Central Terminal is %d (min=%g, max=%g, old=%g)\n",
      center, minimum, maximum, oldval);

   assert(Is_term(g->term[center]));
   *central_term = center;

   return SCIP_OKAY;
}


/** is the given graph valid? */
SCIP_Bool graph_valid(
   SCIP*                 scip,               /**< scip struct */
   const GRAPH*          g                   /**< the graph */
   )
{
   int nterms;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   SCIP_Bool isValid = TRUE;
   SCIP_Bool* nodevisited = NULL;

   assert(scip && g);
   nterms = g->terms;

   for( int k = 0; k < nnodes; k++ )
   {
      int e;

      if( Is_term(g->term[k]) )
         nterms--;

      for( e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         if( g->head[e] != k )
            break;

      if( e != EAT_LAST )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: Head invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n",
               k, e, g->tail[e], g->head[e]);

         goto EXIT;
      }

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         if( g->tail[e] != k )
            break;

      if( e != EAT_LAST )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: Tail invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n",
               k, e, g->tail[e], g->head[e]);

         goto EXIT;
      }
   }

   if( nterms != 0 )
   {
      isValid = FALSE;
      SCIPdebugMessage("*** Graph invalid: Wrong Terminal count, count is %d, should be %d\n",
            g->terms, g->terms - nterms);

      goto EXIT;
   }

   if( (g->source < 0 ) || (g->source >= g->knots) || (g->term[g->source] != 0) )
   {
      isValid = FALSE;
      SCIPdebugMessage("*** Graph invalid: Source invalid, Layer %d, Source %d, Terminal %d\n",
            0, g->source, g->term[g->source]);

      goto EXIT;
   }

   for( int e = 0; e < nedges; e += 2 )
   {
      if( (g->ieat[e] == EAT_FREE) && (g->oeat[e] == EAT_FREE)
         && (g->ieat[e + 1] == EAT_FREE) && (g->oeat[e + 1] == EAT_FREE) )
      {
         continue;
      }

      if( (g->ieat[e] == EAT_FREE) || (g->oeat[e] == EAT_FREE)
         || (g->ieat[e + 1] == EAT_FREE) || (g->oeat[e + 1] == EAT_FREE) )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: FREE invalid, Edge %d/%d\n",
               e, e + 1);

         goto EXIT;
      }

      if( (g->head[e] != g->tail[e + 1]) || (g->tail[e] != g->head[e + 1]) )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: Anti invalid, Edge %d/%d, Tail=%d/%d, Head=%d/%d\n",
               e, e + 1, g->head[e], g->tail[e + 1],  g->tail[e], g->head[e + 1]);

         goto EXIT;
      }
   }

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nodevisited, nnodes) );
   SCIP_CALL_ABORT( graph_trail_arr(scip, g, g->source, nodevisited) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( (g->grad[k] == 0) && ((g->inpbeg[k] != EAT_LAST) || (g->outbeg[k] != EAT_LAST)) )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: Knot %d with Grad 0 has Edges\n", k);

         goto EXIT;
      }

      if( !nodevisited[k] && ((g->grad[k] > 0) || (Is_term(g->term[k]))) )
      {
         if( graph_pc_isPcMw(g) && !graph_pc_knotIsFixedTerm(g, k) )
         {
            continue;
         }

         isValid = FALSE;
#ifdef SCIP_DEBUG
         graph_knot_printInfo(g, k);
         SCIPdebugMessage("*** Graph invalid: Knot %d not connected\n", k);
#endif

         goto EXIT;
      }
   }

   if( isValid && graph_pc_isPcMw(g) )
   {
      isValid = graphisValidPcMw(scip, g, nodevisited);
   }


   EXIT:

   SCIPfreeBufferArrayNull(scip, &nodevisited);
   return isValid;
}

/** is the given graph a variant that is effectively an STP?? */
SCIP_Bool graph_typeIsSpgLike(
   const GRAPH*          g                   /**< the graph */
   )
{
   const int type = g->stp_type;
   assert(g);

   return (type == STP_SPG || type == STP_RSMT || type == STP_GSTP);
}

/** is the given graph undirected? */
SCIP_Bool graph_typeIsUndirected(
   const GRAPH*          g                   /**< the graph */
   )
{
   assert(g);

   if( g->stp_type == STP_SAP || g->stp_type == STP_DHCSTP )
      return FALSE;

   return TRUE;
}
