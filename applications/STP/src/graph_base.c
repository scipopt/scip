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
/*lint -esym(766,string.h)                                                   */
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

#define STP_DELPSEUDO_MAXGRAD   5
#define STP_DELPSEUDO_MAXNEDGES 10
#define STP_DELPSEUDO_NOEDGE    -1
#define STP_DELPSEUDO_SKIPEDGE  -2


/*
 * local functions
 */


/** is vertex a non-leaf (call before transformation was performed)  */
static inline
SCIP_Bool isNonLeaf_pretrans(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the graph */
   int                   vertex              /**< node check */
)
{
   const SCIP_Real prize = g->prize[vertex];
   const SCIP_Real* const cost = g->cost;

   for( int e = g->inpbeg[vertex]; e != EAT_LAST; e = g->ieat[e] )
   {
      if( SCIPisGT(scip, prize, cost[e]) )
         return FALSE;
   }

   return TRUE;
}


/** remove non-leaf terminals by edge weight shifting (call before transformation was performed)  */
static
void removeNonLeafTerms_pretrans(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g                   /**< the graph */
)
{
   const int nnodes = g->knots;

   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(g->term[k]) )
         continue;

      if( isNonLeaf_pretrans(scip, g, k) )
      {
         graph_knot_chg(g, k, STP_TERM_NONLEAF);
      }
   }
}


/** shift costs of non-leaf terminals (call before transformation was performed)  */
static
SCIP_RETCODE shiftNonLeafCosts(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g                   /**< the graph */
)
{
   const int nnodes = g->knots;
   const int nedges = g->edges;
   SCIP_Real* const cost = g->cost;

   assert(!g->cost_org_pc);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->cost_org_pc), nedges) );
   BMScopyMemoryArray(g->cost_org_pc, cost, nedges);

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_nonleafTerm(g->term[k]) )
      {
         const SCIP_Real prize = g->prize[k];

         assert(SCIPisPositive(scip, prize));

         for( int e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         {
            assert(SCIPisLT(scip, cost[e], FARAWAY));
            assert(SCIPisLT(scip, cost[flipedge(e)], FARAWAY));

            cost[e] -= prize;
            assert(SCIPisGE(scip, cost[e], 0.0));

            if( cost[e] < 0.0 )
               cost[e] = 0.0;
         }
      }
   }

   return SCIP_OKAY;
}


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

/** gets replacement edge; helper function for pseudo-elimination */
static
SCIP_RETCODE delPseudoGetReplaceEdges(
   SCIP*                 scip,               /**< SCIP data */
   const GRAPH*          g,                  /**< graph */
   const SCIP_Real*      cutoffs,            /**< cutoff values for each incident edge (or NULL) */
   const SCIP_Real*      cutoffsrev,         /**< reverse edge cutoff values (or NULL if undirected or non-existent) */
   const int*            incedge,            /**< incident edges */
   const SCIP_Real*      ecost,              /**< edge cost */
   const SCIP_Real*      ecostrev,           /**< reverse edge cost */
   const int*            adjvert,            /**< adjacent vertices */
   SCIP_Real             vertexprize,        /**< prize if Pc/Mw */
   int                   degree,             /**< degree of node to be replaced */
   int*                  neigbedge,          /**< neighboring edges array */
   SCIP_Bool*            success             /**< enough replace edges available?  */
)
{
   int* hasharr;
   int edgecount = 0;
   int replacecount = 0;
   const int nspareedges = degree;

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
               if( g->head[e] == adjvert[j] )
               {
                  assert(e >= 0);
                  neigbedge[edgecount - 1] = e;
                  break;
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


#if 0
/** transforms an MWCSP to an SAP */
SCIP_RETCODE graph_MwcsToSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real*            maxweights          /**< array containing the weight of each node */
   )
{
   int e;
   int i;
   int nnodes;
   int nterms = 0;

   assert(maxweights != NULL);
   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->cost != NULL);
   assert(graph->terms == 0);

   nnodes = graph->knots;

   /* count number of terminals, modify incoming edges for non-terminals */
   for( i = 0; i < nnodes; i++ )
   {
      if( SCIPisLT(scip, maxweights[i], 0.0) )
      {
         for( e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
         {
            graph->cost[e] -= maxweights[i];
         }
      }
      else
      {
         graph_knot_chg(graph, i, 0);
         nterms++;
      }
   }
   nterms = 0;
   for( i = 0; i < nnodes; i++ )
   {
      graph->prize[i] = maxweights[i];
      if( Is_term(graph->term[i]) )
      {
         assert(SCIPisGE(scip, maxweights[i], 0.0));
         nterms++;
      }
      else
      {
         assert(SCIPisLT(scip, maxweights[i], 0.0));
      }
   }
   assert(nterms == graph->terms);
   graph->stp_type = STP_MWCSP;

   SCIP_CALL( graph_PcToSap(scip, graph) );
   assert(graph->stp_type == STP_MWCSP);
   return SCIP_OKAY;
}


/** alters the graph for prize collecting problems */
SCIP_RETCODE graph_PcToSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   SCIP_Real* prize;
   int k;
   int root;
   int node;
   int nnodes;
   int nterms;
   int pseudoroot;

   assert(graph != NULL);
   assert(graph->prize != NULL);
   assert(graph->knots == graph->ksize);
   assert(graph->edges == graph->esize);

   prize = graph->prize;
   nnodes = graph->knots;
   nterms = graph->terms;
   graph->norgmodelknots = nnodes;
   graph->norgmodeledges = graph->edges;

   /* for each terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + graph->terms + 2), (graph->esize + graph->terms * 8) , -1) );

   /* create a new nodes */
   for( k = 0; k < nterms; ++k )
      graph_knot_add(graph, -1);

   /* new pseudo-root */
   pseudoroot = graph->knots;
   graph_knot_add(graph, -1);

   /* new root */
   root = graph->knots;
   graph_knot_add(graph, 0);

   nterms = 0;
   for( k = 0; k < nnodes; ++k )
   {
      /* is the kth node a terminal other than the root? */
      if( Is_term(graph->term[k]) )
      {
         /* the copied node */
         node = nnodes + nterms;
         nterms++;

         /* switch the terminal property, mark k */
         graph_knot_chg(graph, k, -2);
         graph_knot_chg(graph, node, 0);
         assert(SCIPisGE(scip, prize[k], 0.0));

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, k, BLOCKED, FARAWAY);
         graph_edge_add(scip, graph, pseudoroot, node, prize[k], FARAWAY);
         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);
         graph_edge_add(scip, graph, k, pseudoroot, 0.0, FARAWAY);
      }
      else if( graph->stp_type != STP_MWCSP )
      {
         prize[k] = 0;
      }
   }
   graph->source = root;
   graph->extended = TRUE;
   assert((nterms + 1) == graph->terms);
   if( graph->stp_type != STP_MWCSP )
      graph->stp_type = STP_PCSPG;

   return SCIP_OKAY;
}




#endif


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

   SCIP_CALL( graph_pack(scip, g, &gp, TRUE) );
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


/** allocates (first and second) and initializes (only second) arrays for PC and MW problems */
SCIP_RETCODE graph_pc_init(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   sizeprize,          /**< size of prize array to allocate (or -1) */
   int                   sizeterm2edge       /**< size of term2edge array to allocate and initialize to -1 (or -1) */
   )
{
   assert(scip != NULL);
   assert(g != NULL);

   if( sizeprize > 0 )
   {
      assert(NULL == g->prize);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->prize), sizeprize) );
   }

   if( sizeterm2edge > 0 )
   {
      assert(NULL == g->term2edge);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->term2edge), sizeterm2edge) );
      for( int i = 0; i < sizeterm2edge; i++ )
         g->term2edge[i] = -1;
   }

   return SCIP_OKAY;
}

/** changes graph of PC and MW problems needed for presolving routines */
SCIP_RETCODE graph_pc_presolInit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   int prev;
   const int root = g->source;
   const int nedges = g->edges;

   if( g->stp_type == STP_RPCSPG )
      return SCIP_OKAY;

   assert(scip != NULL && g != NULL);
   assert(g->rootedgeprevs == NULL);
   assert(nedges > 0 && g->grad[root] > 0);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->rootedgeprevs), nedges) );

   for( int e = 0; e < nedges; e++ )
      g->rootedgeprevs[e] = -1;

   prev = g->outbeg[root];
   assert(prev != EAT_LAST);

   for( int e = g->oeat[prev]; e != EAT_LAST; e = g->oeat[e] )
   {
      g->rootedgeprevs[e] = prev;
      prev = e;
   }

   prev = g->inpbeg[root];
   assert(prev != EAT_LAST);

   for( int e = g->ieat[prev]; e != EAT_LAST; e = g->ieat[e] )
   {
      g->rootedgeprevs[e] = prev;
      prev = e;
   }

   return SCIP_OKAY;
}

/** changes graph of PC and MW problems needed after exiting presolving routines */
void graph_pc_presolExit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   assert(scip != NULL && g != NULL);

   if( g->stp_type == STP_RPCSPG )
      return;

   assert(g->rootedgeprevs != NULL);

   SCIPfreeMemoryArray(scip, &(g->rootedgeprevs));
}

/** checks consistency of term2edge array ONLY for non-extended graphs! */
SCIP_Bool graph_pc_term2edgeConsistent(
   const GRAPH*          g                   /**< the graph */
)
{
   const int root = g->source;

   assert(g != NULL);
   assert(g->term2edge);
   assert(!g->extended);

   if( g->term2edge[g->source] != -1 )
      return FALSE;

   for( int i = 0; i < g->knots; i++ )
   {
      if( Is_anyTerm(g->term[i]) && !graph_pc_knotIsFixedTerm(g, i) && i != root && g->term2edge[i] < 0 )
      {
         SCIPdebugMessage("term2edge consistency fail1 %d \n", i);
         return FALSE;
      }

      if( !Is_anyTerm(g->term[i]) && g->term2edge[i] != -1 )
      {
         SCIPdebugMessage("term2edge consistency fail2 %d \n", i);
         return FALSE;
      }

      if( Is_pseudoTerm(g->term[i]) && i != root )
      {
         int k = -1;
         int e;

         for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         {
            k = g->head[e];
            if( Is_term(g->term[k]) && k != root )
               break;
         }
         assert(e != EAT_LAST);
         assert(k >= 0);

         if( g->term2edge[i] != e )
         {
            SCIPdebugMessage("term2edge consistency fail3 %d \n", i);
            return FALSE;
         }

         if( g->term2edge[k] != flipedge(e) )
         {
            SCIPdebugMessage("term2edge consistency fail4 %d \n", i);
            return FALSE;
         }
      }
   }
   return TRUE;
}

/** change property of node to non-terminal */
void graph_pc_knot2nonTerm(
   GRAPH*                g,                  /**< the graph */
   int                   node                /**< node to be changed */
   )
{
   assert(g      != NULL);
   assert(node   >= 0);
   assert(node   < g->knots);
   assert(g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG || g->stp_type == STP_MWCSP || g->stp_type == STP_RMWCSP);
   assert(g->term2edge);

   if( Is_term(g->term[node]) )
      g->terms--;

   g->term[node] = -1;
   g->term2edge[node] = -1;
}

/** check whether node is fixed terminal */
SCIP_Bool graph_pc_knotIsFixedTerm(
   const GRAPH*          g,                  /**< the graph */
   int                   node                /**< node to be checked */
   )
{
   assert(g      != NULL);
   assert(node   >= 0);
   assert(node   < g->knots);
   assert(graph_pc_isPcMw(g));
   assert(g->term2edge);

#ifndef NDEBUG
   if( Is_term(g->term[node]) && g->term2edge[node] < 0 )
      assert(node == g->source || g->prize[node] == FARAWAY);
#endif

   return (Is_term(g->term[node]) && g->term2edge[node] < 0);
}


/** check whether node is a dummy (pseudo) terminal */
SCIP_Bool graph_pc_knotIsDummyTerm(
   const GRAPH*          g,                  /**< the graph */
   int                   node                /**< node to be checked */
   )
{
   assert(g      != NULL);
   assert(node   >= 0);
   assert(node   < g->knots);
   assert(graph_pc_isPcMw(g));
   assert(g->term2edge);
   assert(graph_pc_knotIsFixedTerm(g, g->source));

   if( node == g->source && !graph_pc_isRootedPcMw(g) )
      return TRUE;

   if( g->extended )
   {
      if( Is_term(g->term[node]) && !graph_pc_knotIsFixedTerm(g, node) )
      {
         assert(g->grad[node] == 2 );

         return TRUE;
      }
   }
   else
   {
      if( Is_pseudoTerm(g->term[node]) )
      {
         assert(g->grad[node] == 2 );
         assert(!graph_pc_knotIsFixedTerm(g, node));

         return TRUE;
      }
   }

   return FALSE;
}

/** check whether terminal is not a leaf in at least one optimal tree */
void graph_pc_termMarkProper(
   const GRAPH*          g,                  /**< the graph */
   int*                  termmark            /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
)
{
   const int nnodes = g->knots;

   assert(!g->extended);

   assert(g && termmark);

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) )
      {
         if( graph_pc_termIsNonLeaf(g, i) )
            termmark[i] = 1;
         else
            termmark[i] = 2;
      }
      else
      {
         termmark[i] = 0;
      }
   }
}

/** check whether terminal is not a leaf (or not contained) in at least one optimal tree */
SCIP_Bool graph_pc_termIsNonLeaf(
   const GRAPH*          g,                  /**< the graph */
   int                   term                /**< terminal to be checked */
   )
{
   int e;

   assert(g      != NULL);
   assert(term   >= 0);
   assert(term   < g->knots);
   assert(graph_pc_isPcMw(g));
   assert(g->term2edge);

   if( graph_pc_knotIsFixedTerm(g, term) )
      return FALSE;

   if( g->extended )
      assert(Is_pseudoTerm(g->term[term]));
   else
      assert(Is_term(g->term[term]));

   for( e = g->inpbeg[term]; e != EAT_LAST; e = g->ieat[e] )
   {
      const int neighbor = g->tail[e];
      if( !graph_pc_knotIsDummyTerm(g, neighbor) && g->cost[e] < g->prize[term] )
      {
         assert(neighbor != g->source || graph_pc_isRootedPcMw(g));
         break;
      }
   }

   return (e == EAT_LAST);
}


/** set high costs for not including given pseudo-terminal */
void graph_pc_enforcePterm(
   GRAPH*          graph,              /**< graph */
   int             pterm               /**< the pseudo-terminal */
)
{
   const int term = graph->head[graph->term2edge[pterm]];
   const int root2term = graph_pc_getRoot2PtermEdge(graph, term);

   assert(graph != NULL && Is_pseudoTerm(graph->term[pterm]));
   assert(graph->term2edge[pterm] >= 0 && Is_term(graph->term[term]));
   assert(graph->cost[root2term] == graph->prize[pterm]);

   /* don't change because of weird prize sum in reduce.c */
   if( graph->prize[pterm] < BLOCKED )
   {
      graph->prize[pterm] = BLOCKED - 1.0;
      graph->cost[root2term] = BLOCKED - 1.0;
   }
}

/** updates term2edge array for new graph */
void graph_pc_updateTerm2edge(
   GRAPH*                newgraph,           /**< the new graph */
   const GRAPH*          oldgraph,           /**< the old graph */
   int                   newtail,            /**< tail in new graph */
   int                   newhead,            /**< head in new graph */
   int                   oldtail,            /**< tail in old graph */
   int                   oldhead             /**< head in old graph */
)
{
   assert(newgraph != NULL);
   assert(oldgraph != NULL);
   assert(newgraph->term2edge != NULL);
   assert(oldgraph->term2edge != NULL);
   assert(newtail >= 0);
   assert(newhead >= 0);
   assert(oldtail >= 0);
   assert(oldhead >= 0);
   assert(oldgraph->extended);
   assert(newgraph->extended);

   assert(newgraph->term2edge != NULL);
   if( oldgraph->term2edge[oldtail] >= 0 && oldgraph->term2edge[oldhead] >= 0 && oldgraph->term[oldtail] != oldgraph->term[oldhead] )
   {
      assert(Is_anyTerm(newgraph->term[newtail]) && Is_anyTerm(newgraph->term[newhead]));
      assert(Is_anyTerm(oldgraph->term[oldtail]) && Is_anyTerm(oldgraph->term[oldhead]));
      assert(oldgraph->source != oldtail && oldgraph->source != oldhead);
      assert(flipedge(newgraph->edges) == newgraph->edges + 1);

      newgraph->term2edge[newtail] = newgraph->edges;
      newgraph->term2edge[newhead] = newgraph->edges + 1;
   }

   assert(-1 == newgraph->term2edge[newgraph->source]);
}


/** mark original graph (without dummy terminals) */
void graph_pc_markOrgGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
)
{
   const int root = g->source;
   const int nnodes = g->knots;

   assert(g != NULL);
   assert(graph_pc_isPcMw(g));
   assert(g->term2edge != NULL);
   assert(g->extended);

   for( int k = 0; k < nnodes; k++ )
      g->mark[k] = (g->grad[k] > 0);

   for( int e = g->outbeg[root]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int head = g->head[e];

      if( Is_term(g->term[head]) && !graph_pc_knotIsFixedTerm(g, head) )
      {
         g->mark[head] = FALSE;
         assert(g->grad[head] == 2);
         assert(SCIPisGT(scip, g->cost[e], 0.0));
      }
   }

   if( !graph_pc_isRootedPcMw(g) )
      g->mark[root] = FALSE;
}

/** mark terminals and switch terminal property to original terminals */
SCIP_RETCODE graph_pc_2org(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   const int root = graph->source;
   const int nnodes = graph->knots;

   assert(graph != NULL);
   assert(graph->term2edge != NULL);
   assert(graph->extended);
   assert(graph_pc_isPcMw(graph));

   for( int k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);

      if( Is_pseudoTerm(graph->term[k]) )
      {
         graph_knot_chg(graph, k, 0);
      }
      else if( Is_term(graph->term[k]) && !graph_pc_knotIsFixedTerm(graph, k) )
      {
         assert(k != root);

         graph->mark[k] = FALSE;
         graph_knot_chg(graph, k, -2);
      }
   }

   // todo adapt edge costs

   // add prizes extra method? Assert that eq to anti


   if( graph_pc_isPc(graph) )
   {
      const int nedges = graph->edges;
      assert(graph->cost_org_pc);


      SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->cost_org_pc), nedges) );
      BMScopyMemoryArray(graph->cost_org_pc, graph->cost, nedges);
   }

   if( graph->stp_type == STP_RPCSPG || graph->stp_type == STP_RMWCSP )
      graph->mark[root] = TRUE;
   else
      graph->mark[root] = FALSE;

   graph->extended = FALSE;

   return SCIP_OKAY;
}

/** unmark terminals and switch terminal property to transformed terminals */
void graph_pc_2trans(
   GRAPH*                graph               /**< the graph */
   )
{
   const int nnodes = graph->knots;;

   assert(graph != NULL);
   assert(!(graph->extended));
   assert(graph_pc_isPcMw(graph));

   for( int k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);

      if( Is_pseudoTerm(graph->term[k]) )
         graph_knot_chg(graph, k, 0);
      else if( Is_term(graph->term[k]) && !graph_pc_knotIsFixedTerm(graph, k) )
      {
         assert(k != graph->source);
         graph_knot_chg(graph, k, -2);
      }
   }

   if( graph_pc_isPc(graph) )
   {
      assert(graph->cost_org_pc);
      BMScopyMemoryArray(graph->cost_org_pc, graph->cost, graph->edges);
      // subtract prizes extra method?
   }


   graph->extended = TRUE;

   return;
}

/** graph_pc_2org if extended */
SCIP_RETCODE graph_pc_2orgcheck(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   assert(graph != NULL);

   if( !graph->extended )
      return SCIP_OKAY;

   SCIP_CALL_ABORT( graph_pc_2org(scip, graph) );

   return SCIP_OKAY;
}

/** graph_pc_2trans if not extended */
void graph_pc_2transcheck(
   GRAPH*                graph               /**< the graph */
   )
{
   assert(graph != NULL);

   if( graph->extended )
      return;

   graph_pc_2trans(graph);
}

/* returns sum of positive vertex weights */
SCIP_Real graph_pc_getPosPrizeSum(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph               /**< the graph */
   )
{
   SCIP_Real prizesum = 0.0;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->prize != NULL);
   assert(!graph->extended);

   for( int i = 0; i < graph->knots; i++ )
      if( Is_term(graph->term[i]) && i != graph->source && graph->prize[i] < BLOCKED )
         prizesum += graph->prize[i];

   return prizesum;
}


/* returns degree of non-root vertex in non-extended graph */
int graph_pc_realDegree(
   const GRAPH*          g,                  /**< graph data structure */
   int                   i,                  /**< the vertex to be checked */
   SCIP_Bool             fixedterm           /**< fixed terminal? */
)
{
   const int ggrad = g->grad[i];
   const SCIP_Bool rpc = (g->stp_type == STP_RPCSPG);
   int newgrad;

   assert(g != NULL);
   assert(!g->extended);
   assert(!Is_pseudoTerm(g->term[i]));
   assert(i != g->source);
   assert(rpc || g->stp_type == STP_PCSPG);

   if( !Is_term(g->term[i]) || fixedterm )
      newgrad = ggrad;
   else if( rpc )
      newgrad = ggrad - 1;
   else
      newgrad = ggrad - 2;

   return newgrad;
}

/** alters the graph for prize collecting problems */
SCIP_RETCODE graph_pc_getSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Real*            offset              /**< offset */
   )
{
   SCIP_Real prizesum;
   int k;
   int e;
   int maxpvert;
   int enext;
   int root;
   int head;
   int nnodes;
   int nterms;
   int stp_type;
   int pseudoroot;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->prize != NULL);
   assert(graph->knots == graph->ksize);
   assert(graph->edges == graph->esize);
   assert(graph->extended);

   nnodes = graph->knots;
   nterms = graph->terms;
   prizesum = 0.0;

   stp_type = graph->stp_type;
   graph->stp_type = STP_SAP;

   /* for each terminal, except for the root, three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_copy(scip, graph, newgraph) );

   graph->stp_type = stp_type;

   SCIP_CALL( graph_resize(scip, (*newgraph), ((*newgraph)->ksize + 1), ((*newgraph)->esize + 2 * (nterms - 1)) , -1) );

   (*newgraph)->source = graph->source;
   root = (*newgraph)->source;

   /* new pseudo-root */
   pseudoroot = (*newgraph)->knots;
   graph_knot_add((*newgraph), -1);

   maxpvert = -1;

   for( k = 0; k < nnodes; k++ )
      if( Is_pseudoTerm(graph->term[k]) && (maxpvert == -1 || graph->prize[k] > graph->prize[maxpvert]) )
         maxpvert = k;

   for( k = 0; k < nnodes; k++ )
   {
      if( Is_pseudoTerm(graph->term[k]) )
      {
         prizesum += graph->prize[k];

         if( stp_type == STP_PCSPG && k != maxpvert )
         {
            SCIP_Real minin = FARAWAY;
            for( e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
            {
               if( !graph_pc_knotIsDummyTerm(graph, graph->tail[e]) && graph->cost[e] < minin )
                  minin = graph->cost[e];
            }
            assert(!SCIPisZero(scip, minin));

            prizesum -= MIN(minin, graph->prize[k]);
         }
      }
   }

   if( stp_type != STP_PCSPG && maxpvert >= 0 )
      prizesum -= graph->prize[maxpvert];

   *offset -= prizesum;

   SCIP_CALL( graph_pc_presolInit(scip, *newgraph) );

   e = (*newgraph)->outbeg[root];

   while( e != EAT_LAST )
   {
      enext = (*newgraph)->oeat[e];
      head = (*newgraph)->head[e];
      if( Is_term((*newgraph)->term[head]) )
      {
         (void) graph_edge_redirect(scip, (*newgraph), e, pseudoroot, head, graph->cost[e], TRUE, FALSE);
         (*newgraph)->cost[flipedge(e)] = FARAWAY;
         assert((*newgraph)->head[e] == head);
         assert((*newgraph)->tail[e] == pseudoroot);
      }
      else
      {
         (*newgraph)->cost[e] = prizesum;
      }

      e = enext;
   }

   graph_pc_presolExit(scip, *newgraph);

   for( k = 0; k < nnodes; k++ )
   {
      /* is the kth node a terminal other than the root? */
      if( Is_pseudoTerm((*newgraph)->term[k]) )
      {
         graph_edge_add(scip, (*newgraph), k, pseudoroot, 0.0, FARAWAY);
      }
   }

   return SCIP_OKAY;
}

/** adapts SAP deriving from PCST or MWCS problem with new big M */
void graph_pc_adaptSap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bigM,               /**< new big M value */
   GRAPH*                graph,              /**< the SAP graph */
   SCIP_Real*            offset              /**< the offset */
   )
{
   SCIP_Real oldbigM;
   const int root = graph->source;

   assert(bigM > 0.0);
   assert(scip != NULL && graph != NULL && offset != NULL);
   assert(graph->outbeg[root] >= 0);

   oldbigM = graph->cost[graph->outbeg[root]];
   assert(oldbigM > 0.0);

   *offset += (oldbigM - bigM);

   printf("new vs old %f, %f \n", bigM, oldbigM);

   for( int e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
   {
      assert(graph->cost[e] == oldbigM);
      graph->cost[e] = bigM;
   }
}


/** alters the graph for prize collecting problems and shifts weights to reduce number of terminal */
SCIP_RETCODE graph_pc_getSapShift(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Real*            offset              /**< offset */
   )
{
   GRAPH* newg;
   SCIP_Real maxp;
   SCIP_Real prizesum;
   int e;
   int root;
   const int nnodes = graph->knots;
   const int stp_type = graph->stp_type;
   int maxvert;
   int pseudoroot;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->prize != NULL);
   assert(graph->knots == graph->ksize);
   assert(graph->edges == graph->esize);
   assert(stp_type == STP_MWCSP || stp_type == STP_PCSPG);
   assert(graph->extended);
   graph->stp_type = STP_SAP;

   /* for each terminal, except for the root, three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_copy(scip, graph, newgraph) );

   graph->stp_type = stp_type;
   newg = *newgraph;

   /* get max prize and max vertex */
   maxvert = -1;
   maxp = -1.0;
   for( int k = 0; k < nnodes; k++ )
      if( Is_pseudoTerm(graph->term[k]) && SCIPisGT(scip, graph->prize[k], maxp) )
      {
         assert(graph->grad[k] > 0);
         maxp = graph->prize[k];
         maxvert = k;
      }

   assert(maxvert >= 0);

   /* shift the costs */
   for( int k = 0; k < nnodes; k++ )
   {
      newg->mark[k] = (newg->grad[k] > 0);
      if( Is_pseudoTerm(graph->term[k]) && k != maxvert )
      {
         SCIP_Real p;

         assert(newg->mark[k]);

         p = graph->prize[k];
         for( e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
            if( SCIPisLT(scip, graph->cost[e], p) && !Is_term(graph->term[graph->tail[e]]) )
               break;

         /* if there is no incoming arc of lower cost than prize[k], make k a common node */
         if( e == EAT_LAST )
         {
            int e2 = -1;
            int term = -1;

            newg->term[k] = -1;

            for( e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
            {
               const int tail = graph->tail[e];

               if( Is_term(graph->term[tail]) )
               {
                  if( tail == graph->source )
                     e2 = e;
                  else
                  {
                     assert(term == -1);
                     term = tail;
                  }
               }
               else
               {
                  newg->cost[e] -= p;
                  assert(SCIPisGE(scip, newg->cost[e], 0.0));
               }
            }
            graph->prize[k] = 0.0;

            (*offset) += p;
            assert(e2 != -1);
            assert(term != -1);

            while( newg->inpbeg[term] != EAT_LAST )
               graph_edge_del(scip, newg, newg->inpbeg[term], FALSE);

            newg->mark[term] = FALSE;
            graph_knot_chg(newg, k, -1);
            graph_knot_chg(newg, term, -1);
            graph_edge_del(scip, newg, e2, FALSE);
         }
      }
   }

   SCIP_CALL( graph_resize(scip, newg, (newg->ksize + 1), (newg->esize + 2 * (newg->terms - 1)) , -1) );

   assert(newg->source == graph->source);
   root = newg->source;

   /* new pseudo-root */
   pseudoroot = newg->knots;
   graph_knot_add(newg, -1);

   prizesum = 0.0;
   for( int k = 0; k < nnodes; k++ )
      if( Is_pseudoTerm(graph->term[k]) )
         prizesum += graph->prize[k];

   prizesum += 1;

   *offset -= prizesum;

   /* move edges to terminal from root to pseudo-root */
   e = newg->outbeg[root];
   while( e != EAT_LAST )
   {
      const int head = newg->head[e];
      const int enext = newg->oeat[e];

      if( Is_term(newg->term[head]) )
      {
         (void) graph_edge_redirect(scip, newg, e, pseudoroot, head, graph->cost[e], TRUE, TRUE);
         newg->cost[flipedge(e)] = FARAWAY;
         assert(newg->head[e] == head);
         assert(newg->tail[e] == pseudoroot);
      }
      else
      {
         newg->cost[e] = prizesum;
      }

      e = enext;
   }

   /* add edges from pterminals to pseudo-root */
   for( int k = 0; k < nnodes; k++ )
   {
      /* is the kth node a terminal other than the root? */
      if( Is_pseudoTerm(newg->term[k]) )
      {
         assert(newg->mark[k]);
         graph_edge_add(scip, newg, k, pseudoroot, 0.0, FARAWAY);
      }
   }

   return SCIP_OKAY;
}



/** initializes biased data structure */
void graph_pc_getBiased(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Bool             fullbias,           /**< fully bias or only dominated terminals? */
   SCIP_Real*            costbiased,         /**< biased costs */
   SCIP_Real*            prizebiased         /**< biased prizes */
)
{
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(graph);
   const int root = graph->source;

   assert(scip && graph && costbiased && prizebiased);
   assert(graph->extended);
   assert(!rpcmw || graph_pc_knotIsFixedTerm(graph, root));

   BMScopyMemoryArray(costbiased, graph->cost, nedges);
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_pseudoTerm(graph->term[k]) && graph->grad[k] != 0 )
      {
         SCIP_Real mincost = FARAWAY;

         for( int e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
         {
            const int head = graph->head[e];
            if( Is_term(graph->term[head]) && !graph_pc_knotIsFixedTerm(graph, head) )
               continue;

            if( graph->cost[e] < mincost )
               mincost = graph->cost[e];
         }

#ifndef NDEBUG
         {
            SCIP_Real mincostt = FARAWAY;

            for( int e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
            {
               if( !rpcmw && graph->tail[e] == root )
                  continue;

               if( graph->cost[e] < mincostt )
                  mincostt = graph->cost[e];
            }
            assert((graph->stp_type != STP_PCSPG && graph->stp_type != STP_RPCSPG )
                  || SCIPisEQ(scip, mincostt, mincost) || graph->grad[k] <= 2 );
         }
#endif

         if( fullbias )
         {
            mincost = MIN(mincost, graph->prize[k]);
         }
         else
         {
            if( mincost < graph->prize[k] )
            {
               prizebiased[k] = graph->prize[k];
               continue;
            }

            mincost = graph->prize[k];
         }

         for( int e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
         {
            if( !rpcmw && graph->tail[e] == root )
               continue;

            if( SCIPisGE(scip, graph->cost[e], FARAWAY) )
               continue;

            costbiased[e] -= mincost;
            assert(!SCIPisNegative(scip, costbiased[e]) || (graph->stp_type != STP_PCSPG && graph->stp_type != STP_RPCSPG));
         }

         prizebiased[k] = graph->prize[k] - mincost;
         assert(!SCIPisNegative(scip, prizebiased[k]));
      }
      else
      {
         if( rpcmw && graph_pc_knotIsFixedTerm(graph, k) )
            prizebiased[k] = graph->prize[k];
         else
            prizebiased[k] = 0.0;
      }
   }
}

/** alters the graph for prize-collecting problems with given root */
SCIP_RETCODE graph_pc_getRsap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   int*                  rootcands,          /**< array containing all vertices that could be used as root */
   int                   nrootcands,         /**< number of all vertices could be used as root */
   int                   root                /**< the root of the new SAP */
   )
{
   GRAPH* p;
   int k;
   int e;
   int e2;
   int head;
   int aterm;
   int proot;
   int nnodes;
   int stp_type;

   assert(graph != NULL);
   assert(graph->prize != NULL);
   assert(graph->knots == graph->ksize);
   assert(graph->edges == graph->esize);

   graph_pc_2transcheck(graph);

   aterm = -1;
   proot = graph->source;
   stp_type = graph->stp_type;
   graph->stp_type = STP_SAP;

   /* copy graph */
   SCIP_CALL( graph_copy(scip, graph, newgraph) );

   p = *newgraph;

   graph->stp_type = stp_type;

   assert(Is_pseudoTerm(graph->term[root]));

   for( e = p->outbeg[root]; e != EAT_LAST; e = p->oeat[e] )
   {
      head = p->head[e];
      if( Is_term(p->term[head]) && head != proot )
      {
         graph_knot_chg(p, head, -1);
         aterm = head;
         graph_edge_del(scip, p, e, FALSE);
         break;
      }
   }

   assert(aterm != -1);
   SCIP_CALL( graph_pc_presolInit(scip, p) );

   for( e = graph->outbeg[proot]; e != EAT_LAST; e = graph->oeat[e] )
   {
      head = graph->head[e];

      assert(graph->head[e] == p->head[e]);
      assert(graph->tail[e] == p->tail[e]);

      if( Is_term(graph->term[head]) && head != aterm )
      {
         assert(Is_term(p->term[head]));

         (void) graph_edge_redirect(scip, p, e, root, head, graph->cost[e], TRUE, FALSE);
         p->cost[flipedge(e)] = FARAWAY;

         for( e2 = p->outbeg[head]; e2 != EAT_LAST; e2 = p->oeat[e2] )
            if( p->head[e2] != root )
               assert(p->term[p->head[e2]] == -2);
      }
      else
      {
         graph_edge_del(scip, p, e, FALSE);
      }
   }

   graph_pc_presolExit(scip, p);

   assert(p->grad[aterm] == 0);

   nnodes = p->knots;
   p->source = root;
   graph_knot_chg(p, root, 0);

   for( k = 0; k < nnodes; k++ )
      p->mark[k] = (p->grad[k] > 0);

   assert(p->grad[graph->source] == 0);

   SCIP_CALL( graph_pc_init(scip, p, nnodes, nnodes) );

   assert(graph->term2edge != NULL);

   for( k = 0; k < nnodes; k++)
   {
      p->term2edge[k] = graph->term2edge[k];
      if( k < graph->norgmodelknots )
         p->prize[k] = graph->prize[k];
      else
         p->prize[k] = 0.0;
   }
   p->term2edge[root] = -1;
   p->term2edge[aterm] = -1;

   if( nrootcands > 0 )
   {
      SCIP_CALL( graph_pc_presolInit(scip, p) );
      for( k = 0; k < nrootcands; k++ )
      {
         aterm = rootcands[k];
         if( aterm == root )
            continue;

         for( e = p->outbeg[aterm]; e != EAT_LAST; e = p->oeat[e] )
         {
            head = p->head[e];

            if( Is_term(p->term[head]) && p->term2edge[head] >= 0 )
            {
               assert(p->grad[head] == 2);
               assert(head != root);

               while( p->outbeg[head] != EAT_LAST )
                  graph_edge_del(scip, p, p->outbeg[head], FALSE);

               graph_knot_chg(p, head, -1);
               break;
            }
         }

         p->term2edge[head] = -1;
         p->term2edge[aterm] = -1;

         assert(e != EAT_LAST);
         graph_knot_chg(p, aterm, 0);
      }
      graph_pc_presolExit(scip, p);
   }

   graph_knot_chg(p, proot, -1);
   p->prize[root] = 0.0;

   return SCIP_OKAY;
}


/** alters the graph for prize collecting problems */
SCIP_RETCODE graph_pc_2pc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   int root;
   const int nnodes = graph->knots;
   int termscount;

   assert(scip && graph->prize);
   assert(graph->edges == graph->esize && nnodes == graph->ksize);

   graph->norgmodeledges = graph->edges;
   graph->norgmodelknots = nnodes;

   /* for PC remove terminal property from non-leaves and reduce terminal count */
   if( graph->stp_type != STP_MWCSP )
      removeNonLeafTerms_pretrans(scip, graph);

   /* for each proper terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + graph->terms + 1), (graph->esize + graph->terms * 6) , -1) );

   /* add future terminals */
   for( int k = 0; k < graph->terms; ++k )
      graph_knot_add(graph, -1);

   /* add new root */
   root = graph->knots;
   graph_knot_add(graph, 0);
   graph->prize[root] = 0.0;

   /* allocate and initialize term2edge array */
   graph_pc_init(scip, graph, -1, graph->knots);
   assert(NULL != graph->term2edge);

   termscount = 0;
   for( int k = 0; k < nnodes; ++k )
   {
      assert(k != graph->source);

      /* is the kth node a proper terminal? */
      if( Is_term(graph->term[k]) )
      {
         /* get the new terminal */
         const int node = nnodes + termscount;
         termscount++;

         /* switch the terminal property, mark k */
         graph_knot_chg(graph, k, -2);
         graph_knot_chg(graph, node, 0);
         graph->prize[node] = 0.0;
         assert(SCIPisGE(scip, graph->prize[k], 0.0));

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, k, 0.0, FARAWAY);
         graph_edge_add(scip, graph, root, node, graph->prize[k], FARAWAY);

         graph->term2edge[k] = graph->edges;
         graph->term2edge[node] = graph->edges + 1;
         assert(graph->edges + 1 == flipedge(graph->edges));

         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);

         assert(graph->head[graph->term2edge[k]] == node);
         assert(graph->head[graph->term2edge[node]] == k);
      }
      else if( graph->stp_type != STP_MWCSP && !Is_nonleafTerm(graph->term[k]) )
      {
         graph->prize[k] = 0.0;
      }
   }

   graph->source = root;
   graph->extended = TRUE;
   assert((termscount + 1) == graph->terms);

   if( graph->stp_type != STP_MWCSP )
   {
      SCIP_CALL( shiftNonLeafCosts(scip, graph) );
      graph->stp_type = STP_PCSPG;
   }

   assert(graph->orgsource == -1);
   SCIPdebugMessage("Transformed to PC \n");

   return SCIP_OKAY;
}


/** alters the graph for rooted prize collecting problems */
SCIP_RETCODE graph_pc_2rpc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   const int root = graph->source;
   const int nnodes = graph->knots;
   int nterms;
   int nfixterms;
   int npotterms;

   assert(graph && graph->prize);
   assert(graph->edges == graph->esize);
   assert(nnodes == graph->ksize);
   assert(root >= 0);

   graph->norgmodeledges = graph->edges;
   graph->norgmodelknots = nnodes;

   nfixterms = 0;
   npotterms = 0;

   removeNonLeafTerms_pretrans(scip, graph);

   /* count number of fixed and potential terminals and adapt prizes */
   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_nonleafTerm(graph->term[i]) )
         continue;

      if( !Is_term(graph->term[i]) )
      {
         graph->prize[i] = 0.0;
      }
      else if( SCIPisGE(scip, graph->prize[i], FARAWAY) )
      {
         assert(graph->prize[i] == FARAWAY);
         nfixterms++;
      }
      else if( SCIPisGT(scip, graph->prize[i], 0.0) )
      {
         assert(i != root);
         assert(Is_term(graph->term[i]));
         npotterms++;
      }
      else
      {
         assert(i != root);
         assert(SCIPisEQ(scip, graph->prize[i], 0.0));
         graph->prize[i] = 0.0;
         graph_knot_chg(graph, i, -1);
      }
   }

   assert(npotterms + nfixterms == graph->terms);

   /* for each terminal, except for the root, one node and two edges (i.e. four arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + npotterms), (graph->esize + npotterms * 4) , -1) );

   /* create new nodes corresponding to potential terminals */
   for( int k = 0; k < npotterms; ++k )
      graph_knot_add(graph, -1);

   /* allocate and initialize term2edge array */
   graph_pc_init(scip, graph, -1, graph->knots);
   assert(graph->term2edge != NULL);

   nterms = 0;

   for( int k = 0; k < nnodes; ++k )
   {
      /* is the kth node a potential terminal? */
      if( Is_term(graph->term[k]) && SCIPisLT(scip, graph->prize[k], FARAWAY) )
      {
         /* the future terminal */
         const int node = nnodes + nterms;
         nterms++;

         assert(k != root);

         /* switch the terminal property, mark k as former terminal */
         graph_knot_chg(graph, k, -2);
         graph_knot_chg(graph, node, 0);
         assert(SCIPisGE(scip, graph->prize[k], 0.0));
         graph->prize[node] = 0.0;

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, node, graph->prize[k], FARAWAY);

         graph->term2edge[k] = graph->edges;
         graph->term2edge[node] = graph->edges + 1;
         assert(graph->edges + 1 == flipedge(graph->edges));

         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);

         assert(graph->head[graph->term2edge[k]] == node);
         assert(graph->head[graph->term2edge[node]] == k);
      }
      else
      {
         assert(graph->prize[k] == FARAWAY || graph->prize[k] == 0.0 || Is_nonleafTerm(graph->term[k]));
      }
   }

   SCIP_CALL( shiftNonLeafCosts(scip, graph) );

   graph->extended = TRUE;
   assert(nterms == npotterms);
   assert(graph->prize[graph->source] == FARAWAY);
   graph->stp_type = STP_RPCSPG;
   graph->orgsource = graph->source;

   SCIPdebugMessage("Transformed problem to (RPC) SAP \n");

   return SCIP_OKAY;
}

/** alters the graph for MWCS problems */
SCIP_RETCODE graph_pc_2mw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real*            maxweights          /**< array containing the weight of each node */
   )
{
   int nnodes;
   int nterms = 0;

   assert(maxweights != NULL);
   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->cost != NULL);
   assert(graph->terms == 0);

   nnodes = graph->knots;

   /* count number of terminals, modify incoming edges for non-terminals */
   for( int i = 0; i < nnodes; i++ )
   {
      if( SCIPisLT(scip, maxweights[i], 0.0) )
      {
         for( int e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
            graph->cost[e] -= maxweights[i];
      }
      else if( SCIPisGT(scip, maxweights[i], 0.0) )
      {
         graph_knot_chg(graph, i, 0);
         nterms++;
      }
   }
   nterms = 0;
   for( int i = 0; i < nnodes; i++ )
   {
      graph->prize[i] = maxweights[i];
      if( Is_term(graph->term[i]) )
      {
         assert(SCIPisGE(scip, maxweights[i], 0.0));
         nterms++;
      }
      else
      {
         assert(SCIPisLE(scip, maxweights[i], 0.0));
      }
   }
   assert(nterms == graph->terms);
   graph->stp_type = STP_MWCSP;

   SCIP_CALL( graph_pc_2pc(scip, graph) );
   assert(graph->stp_type == STP_MWCSP);

   SCIPdebugMessage("Transformed to MW \n");

   return SCIP_OKAY;
}



/** alters the graph for RMWCS problems */
SCIP_RETCODE graph_pc_2rmw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   SCIP_Real* maxweights;
   int i;
   int root;
   int nnodes;
   int npterms;
   int nrterms;
   int maxgrad;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->cost != NULL);

   root = -1;
   maxgrad = -1;
   npterms = 0;
   nrterms = 0;
   nnodes = graph->knots;
   maxweights = graph->prize;

   assert(maxweights != NULL);

   /* count number of terminals, modify incoming edges for non-terminals */
   for( i = 0; i < nnodes; i++ )
   {
      if( SCIPisLT(scip, maxweights[i], 0.0) )
      {
         for( int e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
            graph->cost[e] -= maxweights[i];
      }
      else if( SCIPisGE(scip, maxweights[i], FARAWAY) )
      {
         assert(Is_term(graph->term[i]));
         if( graph->grad[i] > maxgrad )
         {
            root = i;
            maxgrad = graph->grad[i];
         }

         nrterms++;
      }
      else if( SCIPisGT(scip, maxweights[i], 0.0) )
      {
         graph_knot_chg(graph, i, 0);
         npterms++;
      }
   }

   assert(root >= 0);
   assert(graph->terms == (npterms + nrterms));

   graph->norgmodeledges = graph->edges;
   graph->norgmodelknots = nnodes;
   graph->source = root;

   /* for each terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + npterms), (graph->esize + npterms * 4) , -1) );
   maxweights = graph->prize;

   /* create a new nodes */
   for( int k = 0; k < npterms; k++ )
      graph_knot_add(graph, -1);

   /* allocate and initialize term2edge array */
   graph_pc_init(scip, graph, -1, graph->knots);
   assert(graph->term2edge != NULL);

   i = 0;
   for( int k = 0; k < nnodes; ++k )
   {
      /* is the kth node a non-fixed terminal */
      if( Is_term(graph->term[k]) && SCIPisLT(scip, maxweights[k], FARAWAY) )
      {
         /* the copied node */
         const int node = nnodes + i;
         i++;

         /* switch the terminal property, mark k */
         graph_knot_chg(graph, k, -2);
         graph_knot_chg(graph, node, 0);
         graph->prize[node] = 0.0;
         assert(SCIPisGE(scip, maxweights[k], 0.0));

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, node, maxweights[k], FARAWAY);

         graph->term2edge[k] = graph->edges;
         graph->term2edge[node] = graph->edges + 1;
         assert(graph->edges + 1 == flipedge(graph->edges));

         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);

         assert(graph->head[graph->term2edge[k]] == node);
         assert(graph->head[graph->term2edge[node]] == k);
      }
   }

   assert(i == npterms);
   graph->extended = TRUE;
   graph->stp_type = STP_RMWCSP;
   graph->orgsource = graph->source;

   SCIPdebugMessage("Transformed to RMW \n");

   return SCIP_OKAY;
}

/** transforms PCSPG or MWCSP to RPCSPG or RMWCSP if possible */
SCIP_RETCODE graph_pc_pcmw2rooted(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real             prizesum            /**< sum of positive prizes */
   )
{
   int e;
   int p;
   int newroot;
   int maxgrad;
   int nfixedterms;
   const int orgnterms = graph->terms;
   const int root = graph->source;
   const int pc = (graph->stp_type == STP_PCSPG);

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->term2edge != NULL);
   assert(graph->extended);
   assert(pc || graph->stp_type == STP_MWCSP);

   newroot = -1;
   maxgrad = -1;

#ifndef WITH_UG
   printf("attempt transformation to rooted problem \n");
#endif

   nfixedterms = 0;
   e = graph->outbeg[root];
   while( e != EAT_LAST )
   {
      const int enext = graph->oeat[e];
      if( SCIPisGE(scip, graph->cost[e], prizesum) )
      {
         int e2;
         const int k = graph->head[e];

         assert(Is_term(graph->term[k]));
         assert(graph->grad[k] == 2);

         graph->cost[e] = FARAWAY;

         for( e2 = graph->outbeg[k]; e2 != EAT_LAST; e2 = graph->oeat[e2] )
            if( graph->head[e2] != root )
               break;

         p = graph->head[e2];
         assert(e2 == graph->term2edge[k]);
         assert(Is_pseudoTerm(graph->term[p]));
         assert(SCIPisGE(scip, graph->prize[p], prizesum));

         graph->prize[p] = FARAWAY;

         /* delete terminal */
         graph_knot_chg(graph, k, -1);
         while( graph->outbeg[k] != EAT_LAST )
            graph_edge_del(scip, graph, graph->outbeg[k], TRUE);

         graph->term2edge[k] = -1;
         graph->term2edge[p] = -1;

         graph_knot_chg(graph, p, 0);

         if( graph->grad[p] > maxgrad )
         {
            newroot = p;
            maxgrad = graph->grad[p];
         }
         nfixedterms++;
      }
      e = enext;
   }

   /* is there a new root? */
   if( newroot >= 0 )
   {
      graph->source = newroot;

      if( graph->rootedgeprevs != NULL )
         graph_pc_presolExit(scip, graph);

      e = graph->outbeg[root];
      while( e != EAT_LAST )
      {
         const int enext = graph->oeat[e];
         const int k = graph->head[e];
         if( Is_term(graph->term[k]) && !SCIPisZero(scip, graph->cost[e]) )
         {
            (void) graph_edge_redirect(scip, graph, e, newroot, k, graph->cost[e], TRUE, TRUE);
            graph->cost[flipedge(e)] = FARAWAY;
         }
         e = enext;
      }

      /* delete old root */
      graph_knot_chg(graph, root, -1);
      while( graph->outbeg[root] != EAT_LAST )
         graph_edge_del(scip, graph, graph->outbeg[root], TRUE);

      if( pc )
      {
         graph->stp_type = STP_RPCSPG;
#ifndef WITH_UG
         printf("...transformed PC to RPC; fixed %d out of %d terminals \n", nfixedterms, orgnterms - 1);
#endif

      }
      else
      {
         graph->stp_type = STP_RMWCSP;
#ifndef WITH_UG
         printf("...transformed MW to RMW; fixed %d out of %d terminals \n", nfixedterms, orgnterms - 1);
#endif
      }

      assert(orgnterms - 1 == graph->terms);
   }

#ifndef WITH_UG
   if( !graph_pc_isRootedPcMw(graph) )
      printf("...failed \n");
#endif

   return SCIP_OKAY;
}

/* deletes dummy terminal to given node and edge from pseudo-root (if existent) */
void graph_pc_deleteTermExtension(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< index of the terminal */
   )
{
   int e;
   int dummyterm;
   const SCIP_Bool has_pseudoroot = !graph_pc_isRootedPcMw(g);

   assert(g != NULL && graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(Is_term(g->term[i]));
   assert(!graph_pc_knotIsFixedTerm(g, i));
   assert(i != g->source);
   assert(g->pcancestors != NULL && g->term2edge != NULL);

   /* get edge from s to its artificial terminal */
   e = g->term2edge[i];
   assert(e >= 0);

   dummyterm = g->head[e];
   assert(dummyterm != g->source && !g->mark[dummyterm]);

#ifndef NDEBUG
   {
      int e2;
      for( e2 = g->outbeg[i]; e2 != EAT_LAST; e2 = g->oeat[e2] )
         if( Is_pseudoTerm(g->term[g->head[e2]]) )
            break;
      assert(e2 == e);
   }
#endif

   /* delete edge and unmark artificial terminal */
   graph_knot_chg(g, dummyterm, -1);
   graph_edge_del(scip, g, e, TRUE);
   g->term2edge[dummyterm] = -1;

   /* delete remaining incident edge of artificial terminal */
   e = g->inpbeg[dummyterm];

   assert(e != EAT_LAST);
   assert(g->source == g->tail[e]); // || g->source == dummyterm ??? deleteme
   assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));

   graph_edge_del(scip, g, e, TRUE);

   g->term2edge[i] = -1;
   assert(g->inpbeg[dummyterm] == EAT_LAST);

   if( has_pseudoroot )
   {
      /* delete edge from i to pseudoroot */
      for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
         if( g->tail[e] == g->source )
            break;

      assert(e != EAT_LAST);
      assert(g->cost[e] == 0.0);
      graph_edge_del(scip, g, e, TRUE);
   }
}


/** delete a terminal for a (rooted) prize-collecting problem */
int graph_pc_deleteTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< index of the terminal */
   )
{
   int e;
   int t;
   int grad = g->grad[i];

   assert(g != NULL);
   assert(scip != NULL);
   assert(Is_term(g->term[i]));
   assert(i != g->source);
   assert(!graph_pc_knotIsFixedTerm(g, i));

   t = UNKNOWN;

   /* delete terminal */

   assert(g->term2edge[i] != -1);
   graph_pc_knot2nonTerm(g, i);
   g->mark[i] = FALSE;

   while( (e = g->outbeg[i]) != EAT_LAST )
   {
      const int i1 = g->head[e];

      if( Is_pseudoTerm(g->term[i1]) && g->source != i1 )
         t = g->head[e];
      graph_edge_del(scip, g, e, TRUE);
   }

   assert(g->grad[i] == 0);
   assert(t != UNKNOWN);
   assert(g->term2edge != NULL);

   /* delete artificial terminal */

   graph_pc_knot2nonTerm(g, t);
   g->mark[t] = FALSE;
   grad += g->grad[t] - 1;

   while( g->outbeg[t] != EAT_LAST )
      graph_edge_del(scip, g, g->outbeg[t], TRUE);

   return grad;
}


/** subtract a given sum from the prize of a terminal */
void graph_pc_subtractPrize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   SCIP_Real             cost,               /**< cost to be subtracted */
   int                   i                   /**< the terminal */
   )
{
   int e;
   int j;

   assert(scip != NULL);
   assert(g != NULL);
   assert(!g->extended);
   assert(!graph_pc_knotIsFixedTerm(g, i));
   assert(i != g->source);

   g->prize[i] -= cost;
   for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      if( Is_pseudoTerm(g->term[g->head[e]]) )
         break;

   assert(e != EAT_LAST);

   j = g->head[e];

   assert(j != g->source);
   assert(!g->mark[j]);

   for( e = g->inpbeg[j]; e != EAT_LAST; e = g->ieat[e] )
      if( g->source == g->tail[e] )
         break;

   assert(e != EAT_LAST);
   assert(!g->mark[g->tail[e]] || g->stp_type == STP_RPCSPG);

   g->cost[e] -= cost;

   assert(g->stp_type == STP_MWCSP  || g->stp_type == STP_RMWCSP || SCIPisGE(scip, g->prize[i], 0.0));
   assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));
   assert(SCIPisGE(scip, g->prize[i], 0.0) || g->stp_type == STP_MWCSP);
}

/** change prize of a terminal */
void graph_pc_chgPrize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   SCIP_Real             newprize,           /**< prize to be subtracted */
   int                   i                   /**< the terminal */
   )
{
   int e;
   int j;

   assert(scip != NULL);
   assert(g != NULL);
   assert(newprize > 0.0);

   if( g->stp_type == STP_RPCSPG && i == g->source )
      return;

   g->prize[i] = newprize;
   for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      if( Is_pseudoTerm(g->term[g->head[e]]) )
         break;

   assert(e != EAT_LAST);

   j = g->head[e];

   assert(j != g->source);
   assert(!g->mark[j]);

   for( e = g->inpbeg[j]; e != EAT_LAST; e = g->ieat[e] )
      if( g->source == g->tail[e] )
         break;

   assert(e != EAT_LAST);
   assert(!g->mark[g->tail[e]] || g->stp_type == STP_RPCSPG);

   g->cost[e] = newprize;

   assert(g->stp_type == STP_MWCSP  || g->stp_type == STP_RMWCSP || SCIPisGE(scip, g->prize[i], 0.0));
   assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));
   assert(SCIPisGE(scip, g->prize[i], 0.0) || g->stp_type == STP_MWCSP);
}


/** contract ancestors of an edge of (rooted) prize-collecting Steiner tree problem or maximum-weight connected subgraph problem */
SCIP_RETCODE graph_pc_contractNodeAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   t,                  /**< tail node to be contracted (surviving node) */
   int                   s,                  /**< head node to be contracted */
   int                   ets                 /**< edge from t to s or -1 */
   )
{
   assert(g != NULL);
   assert(scip != NULL);

   if( ets < 0 )
   {
      assert(ets == -1);

      for( ets = g->outbeg[t]; ets != EAT_LAST; ets = g->oeat[ets] )
         if( g->head[ets] == s )
            break;
   }

   assert(ets >= 0 && ets < g->edges);

   SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[t]), g->pcancestors[s], NULL) );
   SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[s]), g->pcancestors[t], NULL) );

   SCIP_CALL(SCIPintListNodeAppendCopy(scip, &(g->pcancestors[s]), g->ancestors[ets], NULL));
   SCIP_CALL(SCIPintListNodeAppendCopy(scip, &(g->pcancestors[t]), g->ancestors[ets], NULL));

#if 0
   SCIP_Bool conflict;

   SCIP_CALL( graph_pseudoAncestors_appendCopyEdgeToNode(scip, t, ets, FALSE, g, &conflict) );
   assert(!conflict);

   SCIP_CALL( graph_pseudoAncestors_appendCopyNode(scip, t, s, FALSE, g, &conflict) );
   assert(!conflict);
#endif

   return SCIP_OKAY;
}


/** contract an edge of (rooted) prize-collecting Steiner tree problem or maximum-weight connected subgraph problem */
SCIP_RETCODE graph_pc_contractEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int*                  solnode,            /**< solution nodes or NULL */
   int                   t,                  /**< tail node to be contracted (surviving node) */
   int                   s,                  /**< head node to be contracted */
   int                   i                   /**< terminal to add offset to */
   )
{
   int ets;

   assert(g != NULL);
   assert(scip != NULL);
   assert(Is_term(g->term[i]));
   assert(!g->extended);
   assert(g->source != s);

   /* get edge from t to s */
   for( ets = g->outbeg[t]; ets != EAT_LAST; ets = g->oeat[ets] )
      if( g->head[ets] == s )
         break;

   assert(ets != EAT_LAST);

   /* is either endpoint a fixed terminal? */
   if( graph_pc_knotIsFixedTerm(g, s) || graph_pc_knotIsFixedTerm(g, t) )
   {
      assert(graph_pc_isRootedPcMw(g));

      SCIP_CALL( graph_fixed_moveNodePc(scip, s, g) );
      SCIP_CALL( graph_fixed_addEdge(scip, ets, g) );

      if( !graph_pc_knotIsFixedTerm(g, s) )
      {
         assert(graph_pc_knotIsFixedTerm(g, t));
         if( Is_term(g->term[s]))
            graph_pc_deleteTermExtension(scip, g, s);
      }
      else if( !graph_pc_knotIsFixedTerm(g, t) )
      {
         /* need to make t a fixed term */
         assert(g->source != t);
         graph_pc_deleteTermExtension(scip, g, t);
         assert(g->prize[s] == FARAWAY);
         g->prize[t] = FARAWAY;
      }

      /* contract s into t */
      SCIP_CALL( graph_knot_contract(scip, g, solnode, t, s) );
      g->term2edge[t] = -1;
      assert(g->term2edge[s] == -1);

      return SCIP_OKAY;
   }

   assert(!graph_pc_isRootedPcMw(g) || (s != g->source && t != g->source));

   SCIP_CALL( graph_pc_contractNodeAncestors(scip, g, t, s, ets) );

   /* are both endpoints of the edge to be contracted terminals? */
   if( Is_term(g->term[t]) && Is_term(g->term[s]) )
   {
      graph_pc_deleteTermExtension(scip, g, s);

      if( !graph_pc_knotIsFixedTerm(g, i ) )
         graph_pc_subtractPrize(scip, g, g->cost[ets] - g->prize[s], i);
   }
   else
   {
      if( !graph_pc_knotIsFixedTerm(g, i ) )
      {
         if( g->stp_type != STP_MWCSP && g->stp_type != STP_RMWCSP )
            graph_pc_subtractPrize(scip, g, g->cost[ets], i);
         else
            graph_pc_subtractPrize(scip, g, -(g->prize[s]), i);
      }
   }

   /* contract s into t */
   SCIP_CALL( graph_knot_contract(scip, g, solnode, t, s) );
   assert(g->grad[s] == 0);
   SCIPdebugMessage("PcMw contract: %d into %d \n", s, t);

   return SCIP_OKAY;
}


/** is this graph a prize-collecting variant? */
SCIP_Bool graph_pc_isPc(
   const GRAPH*          g                   /**< the graph */
)
{
   const int type = g->stp_type;
   assert(g != NULL);

   return (type == STP_PCSPG || type == STP_RPCSPG);
}


/** is this graph a prize-collecting or maximum-weight variant? */
SCIP_Bool graph_pc_isPcMw(
   const GRAPH*          g                   /**< the graph */
)
{
   const int type = g->stp_type;
   assert(g != NULL);

   return (type == STP_PCSPG || type == STP_RPCSPG || type == STP_MWCSP || type == STP_RMWCSP || type == STP_BRMWCSP);
}

/** get edge from root to (pseudo) terminal */
int graph_pc_getRoot2PtermEdge(
   const GRAPH*          graph,               /**< the graph */
   int                   pseudoterm           /**< the pseudo terminal  */
)
{
   int rootedge = -1;
   assert(graph != NULL);
   assert(pseudoterm >= 0 && pseudoterm < graph->knots);

   for( int e = graph->inpbeg[pseudoterm]; e != EAT_LAST; e = graph->ieat[e] )
      if( graph->tail[e] == graph->source )
         rootedge = e;
   assert(rootedge >= 0);

   return rootedge;
}

/** get number of fixed terminals */
int graph_pc_nFixedTerms(
   const GRAPH*          graph                /**< the graph */
)
{
   int nfixterms = 0;
   const int nnodes = graph->knots;
   assert(graph != NULL);
   assert(graph_pc_isPcMw(graph));

   if( !graph_pc_isRootedPcMw(graph) )
      return 0;

   for( int k = 0; k < nnodes; k++ )
      if( graph_pc_knotIsFixedTerm(graph, k) )
         nfixterms++;

   return nfixterms;
}

/** get number of potential terminals */
int graph_pc_nPotentialTerms(
   const GRAPH*          graph                /**< the graph */
)
{
   assert(graph != NULL);
   assert(graph_pc_isPcMw(graph));

   if( !graph_pc_isRootedPcMw(graph) )
      return (graph->terms - 1);

   return (graph->terms - graph_pc_nFixedTerms(graph));
}

/** get twin-terminal */
int graph_pc_getTwinTerm(
   const GRAPH*          g,                  /**< the graph */
   int                   vertex              /**< the vertex  */
)
{
   assert(g != NULL);
   assert(graph_pc_isPcMw(g));
   assert(g->term2edge != NULL && g->term2edge[vertex] > 0);
   assert(Is_anyTerm(g->term[vertex]));

   return g->head[g->term2edge[vertex]];
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

/** is this graph a rooted prize-collecting or rooted maximum-weight variant? */
SCIP_Bool graph_pc_isRootedPcMw(
   const GRAPH*          g                   /**< the graph */
)
{
   const int type = g->stp_type;
   assert(g != NULL);

   return (type == STP_RPCSPG || type == STP_RMWCSP || type == STP_BRMWCSP);
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
   assert(term   < p->layers);

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
   int*                  solnode             /**< marks whether an node is part of a given solution (CONNECT), or is NULL */
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
   }

   return SCIP_OKAY;
}


/** pseudo delete node, i.e. reconnect neighbors; maximum degree of 4! */
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
   SINGLETONANS ancestors[STP_DELPSEUDO_MAXGRAD];
   SCIP_Real ecost[STP_DELPSEUDO_MAXGRAD];
   SCIP_Real ecostrev[STP_DELPSEUDO_MAXGRAD];
   SCIP_Real ecostreal[STP_DELPSEUDO_MAXGRAD];
   SCIP_Real vertexprize;
   int incedge[STP_DELPSEUDO_MAXGRAD];
   int adjvert[STP_DELPSEUDO_MAXGRAD];
   int neigbedge[STP_DELPSEUDO_MAXNEDGES];
   int edgecount;
   int nspareedges;
   int replacecount;
   int degree = g->grad[vertex];
   int ancestorsnode;

   assert(scip && success && g);
   assert(vertex >= 0 && vertex < g->knots);
   assert(degree <= STP_DELPSEUDO_MAXGRAD);

#ifndef NDEBUG
   {
      int sum = 0;
      for( int i = 1; i < STP_DELPSEUDO_MAXGRAD; i++ )
         sum += i;
      assert(sum == STP_DELPSEUDO_MAXNEDGES);
   }
#endif

   *success = TRUE;

   if( degree <= 1 )
      return SCIP_OKAY;

   if( Is_term(g->term[vertex]) )
   {
      assert(graph_pc_isPcMw(g));

      graph_pc_deleteTermExtension(scip, g, vertex);

      degree = g->grad[vertex];
      vertexprize = g->prize[vertex];
      ancestorsnode = vertex;
      g->prize[vertex] = 0.0;
      graph_knot_chg(g, vertex, -1);

      assert(vertexprize > 0.0);
      assert(degree == 3);
   }
   else
   {
      vertexprize = 0.0;
      ancestorsnode = -1;

      if( g->pcancestors && g->pcancestors[vertex] )
      {
         assert(graph_pc_isPcMw(g));
         assert(SCIPisZero(scip, g->prize[vertex]));

         ancestorsnode = vertex;
      }
   }

   edgecount = 0;

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

   assert(edgecount == degree);

   SCIP_CALL( delPseudoGetReplaceEdges(scip, g, cutoffs, cutoffsrev, incedge, ecost, ecostrev, adjvert,
         vertexprize, degree, neigbedge, success) );

   /* not enough spare edges? */
   if( !(*success) )
      return SCIP_OKAY;

   nspareedges = degree; /* todo */

   for( int i = 0; i < degree; i++ )
   {
      const int e = incedge[i];
      SCIP_CALL( graph_singletonAncestors_init(scip, g, e, &(ancestors[i])) );
   }

   /* replace edges */
   edgecount = 0;
   replacecount = 0;
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
            const SCIP_Real newijcost = ecostreal[i] + ecostreal[j] - vertexprize;
            const int oldincedge = incedge[(replacecount == nspareedges) ? replacecount - 1 : replacecount];
            const int oldijedge = neigbedge[edgecount];
#ifndef NDEBUG
            const int oldinctail = g->tail[oldincedge];
            const int oldinchead = g->head[oldincedge];
#endif
            assert(replacecount <= nspareedges);
            assert(replacecount < nspareedges || neigbedge[edgecount] != STP_DELPSEUDO_NOEDGE);

            SCIP_CALL( graph_edge_reinsert(scip, g, oldincedge, adjvert[i], adjvert[j], newijcost,
                  ancestorsnode, &(ancestors[i]), &(ancestors[j]), &newijedge, &conflict) );

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
      printf("nodes=%d, edges=%d, terminals=%d, root=%d, type=%s, isExtended=%d\n", g->knots, g->edges, g->terms, g->source, type, g->extended);
   else
      printf("nodes=%d, edges=%d, terminals=%d, root=%d, type=%s \n", g->knots, g->edges, g->terms, g->source, type);

   if( g->stp_type == STP_BRMWCSP )
      printf("budget=%f \n", g->budget);
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

   if( graph_pc_isPcMw(g) && (g->term2edge[k] < 0 || !Is_term(g->term[k])) )
   {
      assert(g->prize != NULL);
      printf("Xnode %d: term=%d grad=%d prize=%f \n", k, g->term[k], g->grad[k], g->prize[k]);
   }
   else
      printf("node %d: term=%d grad=%d  \n", k, g->term[k], g->grad[k]);
}

/** changes solution according to given root */
SCIP_RETCODE graph_sol_reroot(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int*                  result,             /**< solution array (CONNECT/UNKNOWN) */
   int                   newroot             /**< the new root */
   )
{
   int* queue;
   int* gmark;
   int size;
   const int nnodes = g->knots;

   assert(scip != NULL);
   assert(g != NULL);
   assert(result != NULL);
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

   return SCIP_OKAY;
}


/** checks whether edge(s) of given primal solution have been deleted */
SCIP_Bool graph_sol_unreduced(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const int*            result              /**< solution array, indicating whether an edge is in the solution */
   )
{
   const int nedges = graph->edges;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(result != NULL);

   for( int i = 0; i < nedges; i++ )
      if( result[i] == CONNECT && graph->oeat[i] == EAT_FREE )
         return FALSE;

   return TRUE;
}

/** verifies whether a given primal solution is feasible */
SCIP_Bool graph_sol_valid(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const int*            result              /**< solution array, indicating whether an edge is in the solution */
   )
{
   int* queue;
   STP_Bool* reached;
   int root;
   int size;
   int nnodes;
   int termcount;
   SCIP_Bool usepterms;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(result != NULL);

   reached = NULL;
   nnodes = graph->knots;
   root = graph->source;
   assert(root >= 0);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &reached, nnodes) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &queue, nnodes) );

   if( graph_pc_isPcMw(graph) && !graph->extended )
      usepterms = TRUE;
   else
      usepterms = FALSE;

   assert(reached != NULL);

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
               return FALSE;
            }

            if( usepterms)
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

#if 0
   if(termcount != graph->terms)
   {
      printf("termcount %d graph->terms %d \n", termcount, graph->terms);
      printf("root %d \n", root);

      for( int i = 0; i < nnodes && 0; i++ )
      {
         if( Is_term(graph->term[i]) && !reached[i] )
         {
            printf("fail %d grad %d\n", i, graph->grad[i]);
            for( int e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
            {
               printf("tail %d %d \n", graph->tail[e], graph->term[graph->tail[e]]);
            }
         }
      }
   }
#endif
   SCIPfreeBufferArray(scip, &queue);
   SCIPfreeBufferArray(scip, &reached);

   return (termcount == graph->terms);
}

/** mark endpoints of edges in given list */
void graph_sol_setNodeList(
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
SCIP_Real graph_sol_getObj(
   const SCIP_Real*      edgecost,
   const int*            soledge,
   SCIP_Real             offset,
   int                   nedges
   )
{
   SCIP_Real obj = offset;
   int e;

   for( e = 0; e < nedges; e++ )
      if( soledge[e] == CONNECT )
         obj += edgecost[e];

   return obj;
}

/** marks vertices for given edge-solution array (CONNECT/UNKNOWN) */
void graph_sol_setVertexFromEdge(
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
      if( result[e] == CONNECT )
         solnode[g->head[e]] = TRUE;

#ifndef NDEBUG
   for( int e = 0; e < nedges; e++ )
      if( result[e] == CONNECT )
         assert(solnode[g->head[e]] && solnode[g->tail[e]]);
#endif
}

/** get original solution */
SCIP_RETCODE graph_sol_getOrg(
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
         graph_sol_setNodeList(orggraph, orgnodearr, ancestors[e]);

   /* retransform edges fixed during graph reduction */
   graph_sol_setNodeList(orggraph, orgnodearr, graph_get_fixedges(scip, transgraph));

   if( pcmw )
   {
      // potentially single-vertex solution?
      if( graph_pc_isRootedPcMw(transgraph) && transgraph->terms == 1 && graph_pc_nFixedTerms(orggraph) == 1 )
         orgnodearr[orggraph->source] = TRUE;

      SCIP_CALL( graph_sol_markPcancestors(scip, transgraph->pcancestors, orggraph->tail, orggraph->head, orgnnodes,
            orgnodearr, NULL, NULL, NULL, NULL ) );
   }

   for( int e = 0; e < orggraph->edges; e++ )
      orgsoledge[e] = UNKNOWN;

   /* prune solution (in original graph) */
   if( pcmw )
      SCIP_CALL( SCIPStpHeurTMPrunePc(scip, orggraph, orggraph->cost, orgsoledge, orgnodearr) );
   else
      SCIP_CALL( SCIPStpHeurTMPrune(scip, orggraph, orggraph->cost, 0, orgsoledge, orgnodearr) );

   SCIPfreeBufferArray(scip, &orgnodearr);

   assert(graph_sol_valid(scip, orggraph, orgsoledge));

   return SCIP_OKAY;
}


/** mark original solution */
SCIP_RETCODE graph_sol_markPcancestors(
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

/** get (real) number of nodes , edges, terminals */
void graph_get_NVET(
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
      for( int e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
      {
         edgearr[i] = e;
         tailarr[i++] = g->tail[e] + 1;
      }
   }

   *nnewedges = i;
   start[nnodes] = i;
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
         if( pcmw && graph_pc_termIsNonLeaf(g, i) )
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

   SCIPdebugMessage("Initialized new graph \n");

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

   if( p->prize != NULL )
      SCIPfreeMemoryArray(scip, &(p->prize));

   if( p->costbudget != NULL )
      SCIPfreeMemoryArray(scip, &(p->costbudget));

   if( p->term2edge != NULL )
      SCIPfreeMemoryArray(scip, &(p->term2edge));

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
   GRAPH* g = copygraph;
   const GRAPH* p = orgraph;
   const int ksize = p->ksize;
   const int esize = p->esize;

   assert(scip != NULL);
   assert(orgraph != NULL);
   assert(copygraph != NULL);
   assert(ksize == g->ksize && ksize > 0);
   assert(esize == g->esize && esize >= 0);
   assert(p->dcsr_storage == NULL && p->csr_storage == NULL);

   g->norgmodeledges = p->norgmodeledges;
   g->norgmodelknots = p->norgmodelknots;
   g->knots = p->knots;
   g->terms = p->terms;
   g->edges = p->edges;
   g->source = p->source;
   g->orgsource = p->orgsource;
   g->orgedges = p->orgedges;
   g->orgknots = p->orgknots;
   g->grid_dim = p->grid_dim;
   g->stp_type = p->stp_type;
   g->hoplimit = p->hoplimit;
   g->extended = p->extended;
   g->budget = p->budget;
   g->is_packed = p->is_packed;

   BMScopyMemoryArray(g->term, p->term, ksize);
   BMScopyMemoryArray(g->mark, p->mark, ksize);
   BMScopyMemoryArray(g->grad, p->grad, ksize);
   BMScopyMemoryArray(g->inpbeg, p->inpbeg, ksize);
   BMScopyMemoryArray(g->outbeg, p->outbeg, ksize);
   BMScopyMemoryArray(g->cost, p->cost, esize);
   BMScopyMemoryArray(g->tail, p->tail, esize);
   BMScopyMemoryArray(g->head, p->head, esize);
   BMScopyMemoryArray(g->ieat, p->ieat, esize);
   BMScopyMemoryArray(g->oeat, p->oeat, esize);

   if( graph_pc_isPcMw(g) )
   {
      const SCIP_Bool rpcmw = graph_pc_isRootedPcMw(g);
      const SCIP_Bool brpcmw = (g->stp_type == STP_BRMWCSP);

      if( g->prize != NULL )
         SCIPfreeMemoryArray(scip, &(g->prize));

      if( g->costbudget != NULL )
      {
         assert(brpcmw);
         SCIPfreeMemoryArray(scip, &(g->costbudget));
      }

      if( g->term2edge != NULL )
         SCIPfreeMemoryArray(scip, &(g->term2edge));

      SCIP_CALL(SCIPallocMemoryArray(scip, &(g->prize), g->knots));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(g->term2edge), g->knots));

      if( brpcmw )
      {
         SCIP_CALL(SCIPallocMemoryArray(scip, &(g->costbudget), g->knots));
         BMScopyMemoryArray(g->costbudget, p->costbudget, g->knots);
      }

      for( int k = 0; k < g->knots; k++ )
         if( Is_term(p->term[k]) && (!rpcmw || !graph_pc_knotIsFixedTerm(p, k)) )
            g->prize[k] = 0.0;
         else
            g->prize[k] = p->prize[k];

      assert(p->term2edge != NULL);

      BMScopyMemoryArray(g->term2edge, p->term2edge, g->knots);
   }
   else if( g->stp_type == STP_DCSTP )
   {
      assert(p->maxdeg != NULL);

      SCIP_CALL(SCIPallocMemoryArray(scip, &(g->maxdeg), g->knots));

      for( int k = 0; k < g->knots; k++ )
         g->maxdeg[k] = p->maxdeg[k];
   }
   else if( p->stp_type == STP_RSMT )
   {
      assert(p->grid_ncoords != NULL);
      assert(p->grid_coordinates != NULL);

      SCIP_CALL(SCIPallocMemoryArray(scip, &(g->grid_coordinates), p->grid_dim));

      BMScopyMemoryArray(g->grid_coordinates, p->grid_coordinates, p->grid_dim);
      for( int k = 0; k < p->grid_dim; k++ )
      {
         SCIP_CALL(SCIPallocMemoryArray(scip, &(g->grid_coordinates[k]), p->terms)); /*lint !e866*/
         BMScopyMemoryArray(g->grid_coordinates[k], p->grid_coordinates[k], p->terms); /*lint !e866*/
      }
      SCIP_CALL(SCIPallocMemoryArray(scip, &(g->grid_ncoords), p->grid_dim));

      BMScopyMemoryArray(g->grid_ncoords, p->grid_ncoords, p->grid_dim);
   }
   assert(graph_valid(scip, g));

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
   const int nnodes = g->knots;

   assert(g);

   for( int k = 0; k < nnodes; k++ )
      g->mark[k] = (g->grad[k] > 0);

   if( graph_pc_isPcMw(g) && !g->extended )
   {
      const int root = g->source;

      for( int k = 0; k < nnodes; k++ )
      {
         if( (Is_pseudoTerm(g->term[k]) || k == root) )
            g->mark[k] = FALSE;
      }

      if( graph_pc_isRootedPcMw(g) )
         g->mark[root] = TRUE;
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


/** pack the graph, i.e. build a new graph that discards deleted edges and nodes */
SCIP_RETCODE graph_pack(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Bool             verbose             /**< verbose? */
   )
{
   GRAPH* g;
   GRAPH* q;
   int*   new;
   int    e;
   int    oldnnodes;
   int    oldnedges;
   int    nnodes;
   int    nedges;
   SCIP_Bool rpcmw;
   SCIP_Bool pcmw;

   assert(scip      != NULL);
   assert(graph     != NULL);
   assert(graph_valid(scip, graph));
   assert(!graph->is_packed);

   g = graph;
   nnodes = 0;
   nedges = 0;
   oldnnodes = g->knots;
   oldnedges = g->edges;
   SCIP_CALL( SCIPallocBufferArray(scip, &new, oldnnodes) );

   if( verbose )
      printf("Reduced graph: ");

   /* count nodes */
   for( int i = 0; i < oldnnodes; i++ )
   {
      /* are there incident edges to current node? */
      if( g->grad[i] > 0 )
         new[i] = nnodes++;
      else
         new[i] = -1;
   }

   /* graph vanished? */
   if( nnodes == 0 )
   {
      SCIPfreeBufferArray(scip, &new);
      new = NULL;
      if( verbose )
         printf(" graph vanished!\n");

      nnodes = 1;
   }

   /* count edges */
   for( int i = 0; i < oldnedges; i++ )
   {
      if( g->oeat[i] != EAT_FREE )
      {
         assert(g->ieat[i] != EAT_FREE);
         nedges++;
      }
   }

   assert(nnodes > 1 || nedges == 0);
   SCIP_CALL( graph_init(scip, newgraph, nnodes, nedges, g->layers) );
   q = *newgraph;
   q->norgmodelknots = g->norgmodelknots;
   q->norgmodeledges = g->norgmodeledges;
   q->orgsource = g->orgsource;
   q->orgtail = g->orgtail;
   q->orghead = g->orghead;
   q->orgknots = g->knots;
   q->orgedges = g->edges;
   q->stp_type = g->stp_type;
   q->maxdeg = g->maxdeg;
   q->grid_dim = g->grid_dim;
   q->grid_ncoords = g->grid_ncoords;
   q->grid_coordinates = g->grid_coordinates;
   q->pcancestors = g->pcancestors;
   q->fixedcomponents = g->fixedcomponents;
   q->hoplimit = g->hoplimit;
   q->extended = g->extended;
   q->budget = g->budget;
   q->is_packed = TRUE;

   if( new == NULL )
   {
      q->ancestors = NULL;
      graph_free(scip, &g, FALSE);

      if( q->stp_type == STP_RSMT )
      {
         q->grid_ncoords = NULL;
         q->grid_coordinates = NULL;
      }

      graph_knot_add(q, 0);
      q->source = 0;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(q->ancestors), nedges) );

   rpcmw = graph_pc_isRootedPcMw(g);
   pcmw = (rpcmw || g->stp_type == STP_MWCSP || g->stp_type == STP_PCSPG);

   if( pcmw )
   {
      SCIP_CALL( graph_pc_init(scip, q, nnodes, nnodes) );
      if( g->stp_type == STP_BRMWCSP )
         SCIP_CALL( SCIPallocMemoryArray(scip, &(q->costbudget), nnodes) );
   }

   /* add nodes (of positive degree) */

   if( rpcmw )
   {
      for( int i = 0; i < oldnnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

      for( e = g->outbeg[g->source]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int i = g->head[e];

         if( SCIPisGT(scip, g->cost[e], 0.0) && Is_term(g->term[i]) && g->term2edge[i] >= 0 )
         {
            g->mark[i] = FALSE;
            assert(g->grad[i] == 2);
         }
      }
   }

   for( int i = 0; i < oldnnodes; i++ )
   {
      assert(g->term[i] < g->layers);
      if( g->grad[i] > 0 )
      {
         if( pcmw )
         {
            if( !Is_term(g->term[i]) || (rpcmw && g->mark[i]) )
               q->prize[q->knots] = g->prize[i];
            else
               q->prize[q->knots] = 0.0;
         }
         if( g->stp_type == STP_BRMWCSP )
         {
            if( !Is_term(g->term[i]) || (g->mark[i]) )
               q->costbudget[q->knots] = g->costbudget[i];
            else
               q->costbudget[q->knots] = 0.0;
         }
         graph_knot_add(q, g->term[i]);
      }
   }

   /* add root */
   assert(q->term[new[g->source]] == 0);

   q->source = new[g->source];

   if( rpcmw )
      assert(FARAWAY == q->prize[q->source]);

   /* add edges */
   for( int i = 0; i < oldnedges; i += 2 )
   {
      if( g->ieat[i] == EAT_FREE )
      {
         assert(g->oeat[i]     == EAT_FREE);
         assert(g->ieat[i + 1] == EAT_FREE);
         assert(g->oeat[i + 1] == EAT_FREE);

         graph_edge_delHistory(scip, g, i);
         continue;
      }

      assert(g->oeat[i]      != EAT_FREE);
      assert(g->ieat[i + 1]  != EAT_FREE);
      assert(g->oeat[i + 1]  != EAT_FREE);
      assert(new[g->tail[i]] >= 0);
      assert(new[g->head[i]] >= 0);

      e = q->edges;

      q->ancestors[e] = NULL;
      q->ancestors[e + 1] = NULL;
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(q->ancestors[e]), g->ancestors[i], NULL) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(q->ancestors[e + 1]), g->ancestors[i + 1], NULL) );

      assert(new[g->tail[i]] < nnodes && new[g->head[i]] < nnodes);

      if( pcmw )
         graph_pc_updateTerm2edge(q, g, new[g->tail[i]], new[g->head[i]], g->tail[i], g->head[i]);

      graph_edge_add(scip, q, new[g->tail[i]], new[g->head[i]], g->cost[i], g->cost[Edge_anti(i)]);
   }

   SCIPfreeBufferArray(scip, &new);

   if( g->path_heap != NULL )
      graph_path_exit(scip, g);

   g->stp_type = UNKNOWN;
   graph_free(scip, &g, FALSE);

   assert(graph_valid(scip, q));

   if( verbose )
      printf("Nodes: %d  Edges: %d  Terminals: %d\n", q->knots, q->edges, q->terms);

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
         SCIPdebugMessage("*** Graph invalid: Knot %d with Grad 0 has Edges\n",
               k);

         goto EXIT;
      }

      if( !nodevisited[k] && ((g->grad[k] > 0) || (Is_term(g->term[k])))
         && g->stp_type != STP_PCSPG && g->stp_type != STP_MWCSP && g->stp_type != STP_RMWCSP )
      {
         isValid = FALSE;
         SCIPdebugMessage("*** Graph invalid: Knot %d not connected\n",
               k);

         goto EXIT;
      }
   }

   if( graph_pc_isPcMw(g) )
   {
      int npterms = 0;
      const int root = g->source;
      const SCIP_Bool extended = g->extended;
      const SCIP_Bool rooted = graph_pc_isRootedPcMw(g);
      nterms = 0;

      assert(g->prize != NULL);
      assert(g->term2edge != NULL);

      for( int k = 0; k < nnodes; k++ )
      {
         if( rooted && g->term2edge[k] < 0 && Is_pseudoTerm(g->term[k]) )
         {
            assert(k != root);

            SCIPdebugMessage("inconsistent term2edge for %d \n", k);
            isValid = FALSE;
            goto EXIT;
         }

         if( rooted && graph_pc_knotIsFixedTerm(g, k) && g->prize[k] != FARAWAY )
         {
            SCIPdebugMessage("inconsistent prize for %d \n", k);
            isValid = FALSE;
            goto EXIT;
         }

         if( k == root || (rooted && g->term2edge[k] < 0) )
            continue;

         if( (extended ? Is_term(g->term[k]) : Is_pseudoTerm(g->term[k])) )
         {
            int e;
            int e2;
            int pterm;
            const int term = k;
            nterms++;

            if( g->grad[k] != 2 )
            {
               SCIPdebugMessage("terminal degree != 2 for %d \n", k);
               isValid = FALSE;
               goto EXIT;
            }

            for( e = g->inpbeg[term]; e != EAT_LAST; e = g->ieat[e] )
               if( g->tail[e] == root )
                  break;

            if( e == EAT_LAST )
            {
               SCIPdebugMessage("no edge to root for term %d \n", term);
               isValid = FALSE;
               goto EXIT;
            }

            for( e2 = g->outbeg[term]; e2 != EAT_LAST; e2 = g->oeat[e2] )
            {
               pterm = g->head[e2];
               if( (extended ? Is_pseudoTerm(g->term[pterm]) : Is_term(g->term[pterm])) && pterm != root  )
                  break;
            }

            if( e2 == EAT_LAST)
            {
               SCIPdebugMessage("no terminal for dummy %d \n", g->head[e2]);
               isValid = FALSE;
               goto EXIT;
            }

            assert(pterm != root);

            if( e2 != g->term2edge[term] )
            {
               SCIPdebugMessage("term2edge for node %d faulty \n", term);
               isValid = FALSE;
               goto EXIT;
            }

            if( g->cost[e] != g->prize[pterm] )
            {
               SCIPdebugMessage("prize mismatch for node %d: \n", k);
               isValid = FALSE;
               goto EXIT;
            }
         }
         else if( (extended ? Is_pseudoTerm(g->term[k]) : Is_term(g->term[k])) )
         {
            npterms++;
         }
      }

      if( nterms != npterms || nterms != g->terms - 1 )
      {
         if( !rooted )
         {
            SCIPdebugMessage("wrong terminal count \n");
            isValid = FALSE;
            goto EXIT;
         }
      }

      if( rooted )
      {
         SCIP_CALL_ABORT( graph_trail_costAware(scip, g, g->source, nodevisited) );

         for( int k = 0; k < nnodes; k++ )
         {
            if( graph_pc_knotIsFixedTerm(g, k) && !nodevisited[k] )
            {
               SCIPdebugMessage("disconnected fixed terminal %d \n", k);
               isValid = FALSE;
               goto EXIT;
            }
         }
      }
   } /* is PC/MW */


   EXIT:

   SCIPfreeBufferArrayNull(scip, &nodevisited);
   return isValid;
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
