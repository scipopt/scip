/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   grphbase.c
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


#include "scip/misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "portab.h"
#include "misc_stp.h"
#include "scip/misc.h"
#include "grph.h"

#define BLK 0

/** initialize graph */
SCIP_RETCODE graph_init(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               g,                  /**< new graph */
   int                   ksize,              /**< slots for nodes */
   int                   esize,              /**< slots for edges */
   int                   layers,             /**< number of layers (only needed for packing, otherwise 1) */
   int                   flags               /**< flags */
   )
{
   GRAPH* p;
   int    i;

   assert(ksize > 0);
   assert(ksize < INT_MAX);
   assert(esize >= 0);
   assert(esize < INT_MAX);
   assert(layers > 0);
   assert(layers < SHRT_MAX);

#if BLK
   SCIP_CALL( SCIPallocBlockMemory(scip, g) );
#else
   SCIP_CALL( SCIPallocMemory(scip, g) );
#endif
   p = *g;
   assert(p != NULL);

   /* ancestor data for retransformation after reductions */
   p->fixedges = NULL;
   p->ancestors = NULL;
   p->pcancestors = NULL;
   p->orgtail = NULL;
   p->orghead = NULL;
   p->norgmodelknots = 0;
   p->norgmodeledges = 0;
   p->ksize  = ksize;
   p->orgknots = 0;
   p->orgedges = 0;
   p->knots  = 0;
   p->terms  = 0;
   p->orgsource = UNKNOWN;
   p->stp_type = UNKNOWN;
   p->flags  = flags;
   p->layers = layers;
   p->hoplimit = UNKNOWN;

#if BLK
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->locals), layers) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->source), layers) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->term), ksize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->mark), ksize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->grad), ksize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->inpbeg), ksize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->outbeg), ksize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->cost), esize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->tail), esize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->head), esize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->ieat), esize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(p->oeat), esize) );
#else
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->locals), layers) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->source), layers) );
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
#endif

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

   for( i = 0; i < p->layers; i++ )
   {
      p->source[i] = -1;
      p->locals[i] =  0;
   }
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
   int e;
   int nedges;
   SCIP_Bool pcmw;

   assert(scip != NULL);
   assert(graph != NULL);

   pcmw = (graph->stp_type == STP_PCSPG || graph->stp_type == STP_RPCSPG || graph->stp_type == STP_MWCSP || graph->stp_type == STP_RMWCSP);

   nedges = graph->edges;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->orgtail), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->orghead), nedges) );

   tail = graph->tail;
   head = graph->head;
   orgtail = graph->orgtail;
   orghead = graph->orghead;

   for( e = 0; e < nedges; e++ )
   {
      orgtail[e] = tail[e];
      orghead[e] = head[e];
   }

   if( pcmw )
   {
      int k;
      int nnodes = graph->knots;

      SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->pcancestors), nnodes) );

      pcancestors = graph->pcancestors;

      for( k = 0; k < nnodes; k++ )
         pcancestors[k] = NULL;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->ancestors), nedges) );

   ancestors = graph->ancestors;

   for( e = 0; e < nedges; e++ )
   {
      SCIP_CALL( SCIPallocMemory(scip, &(ancestors[e])) ); /*lint !e866*/
      (ancestors)[e]->index = e;
      (ancestors)[e]->parent = NULL;
   }

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
   int i;
   assert(scip      != NULL);
   assert(g      != NULL);
   assert((ksize  < 0) || (ksize  >= g->knots));
   assert((esize  < 0) || (esize  >= g->edges));
   assert((layers < 0) || (layers >= g->layers));

   if( (layers > 0) && (layers != g->layers) )
   {
#if BLK
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->locals), g->layers, layers) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->source), g->layers, layers) );
#else
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->locals), layers) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->source), layers) );
#endif
      for( i = g->layers; i < layers; i++ )
      {
         g->source[i] = -1;
         g->locals[i] =  0;
      }
      g->layers = layers;
   }
   if( (ksize > 0) && (ksize != g->ksize) )
   {
#if BLK
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->term), g->ksize, ksize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->mark), g->ksize, ksize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->grad), g->ksize, ksize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->inpbeg), g->ksize, ksize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->outbeg), g->ksize, ksize) );
#else
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->term), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->mark), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->grad), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->inpbeg), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->outbeg), ksize) );
#endif
      g->ksize  = ksize;
   }
   if( (esize > 0) && (esize != g->esize) )
   {

#if BLK
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->cost), g->esize, esize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->tail), g->esize, esize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->head), g->esize, esize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->ieat), g->esize, esize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(g->oeat), g->esize, esize) );
#else
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->cost), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->tail), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->head), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->ieat), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->oeat), esize) );
#endif
      g->esize = esize;
   }

   return SCIP_OKAY;
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
   SCIP_CALL( graph_init(scip, gridgraph, nnodes, 2 * nedges, 1, 0) );

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
         g->source[0] = k;

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
   SCIP_CALL( graph_init(scip, gridgraph, nnodes, 2 * nedges, 1, 0) );

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

/** alters the graph for prize collecting problems */
SCIP_RETCODE graph_prize_transform(
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
   assert(graph != NULL);
   assert(graph->edges == graph->esize);
   nnodes = graph->knots;
   nterms = graph->terms;
   prize = graph->prize;
   assert(prize != NULL);
   assert(nnodes == graph->ksize);
   graph->norgmodeledges = graph->edges;
   graph->norgmodelknots = nnodes;

   /* for each terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + graph->terms + 1), (graph->esize + graph->terms * 6) , -1) );

   /* create a new nodes */
   for( k = 0; k < nterms; ++k )
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
         graph_edge_add(scip, graph, root, k, 0.0, FARAWAY);
         graph_edge_add(scip, graph, root, node, prize[k], FARAWAY);
         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);
      }
      else if( graph->stp_type != STP_MWCSP )
      {
         prize[k] = 0;
      }
   }
   graph->source[0] = root;
   assert((nterms + 1) == graph->terms);
   if( graph->stp_type != STP_MWCSP )
      graph->stp_type = STP_PCSPG;

   return SCIP_OKAY;
}


/** changes solution according to given root */
SCIP_RETCODE graph_RerootSol(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int*                  result,             /**< solution array (CONNECT/UNKNOWN) */
   int                   newroot             /**< the new root */
   )
{
   SCIP_QUEUE* queue;
   int a;
   int k;
   int head;
   int node;
   int nnodes;
   int* pnode;
   int* gmark;

   assert(scip != NULL);
   assert(g != NULL);
   assert(result != NULL);
   assert(Is_term(g->term[newroot]));

   if( g->grad[newroot] == 0 )
      return SCIP_OKAY;

   gmark = g->mark;
   nnodes = g->knots;

   for( k = 0; k < nnodes; k++ )
      gmark[k] = FALSE;

   gmark[newroot] = TRUE;

   SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 1.1) );
   SCIP_CALL( SCIPqueueInsert(queue, &newroot));

   /* BFS loop */
   while( !SCIPqueueIsEmpty(queue) )
   {
      pnode = (SCIPqueueRemove(queue));
      node = *pnode;

      /* traverse outgoing arcs */
      for( a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
      {
         head = g->head[a];

         if( !gmark[head] && (result[a] == CONNECT || result[flipedge(a)] == CONNECT ) )
         {
            if( result[flipedge(a)] == CONNECT  )
            {
               result[a] = CONNECT;
               result[flipedge(a)] = UNKNOWN;
            }
            gmark[head] = TRUE;
            SCIP_CALL(SCIPqueueInsert(queue, &(g->head[a])));
         }
      }
   }
   SCIPqueueFree(&queue);

   /* adjust solution if infeasible */
   for( k = 0; k < nnodes; k++ )
   {
      if( !gmark[k] )
      {
         for( a = g->outbeg[k]; a != EAT_LAST; a = g->oeat[a] )
         {
            result[a] = UNKNOWN;
            result[flipedge(a)] = UNKNOWN;
         }

         /* not yet connected terminal? */
         if( Is_term(g->term[k]) )
         {
            for( a = g->inpbeg[k]; a != EAT_LAST; a = g->ieat[a] )
            {
               node = g->tail[a];
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
                  node = g->tail[a];
                  if( node == newroot )
                  {
                     result[a] = CONNECT;
                     break;
                  }
               }
            }
            if( a != EAT_LAST )
               gmark[k] = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}


/** alters the graph for prize collecting problems */
SCIP_RETCODE graph_PcSapCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Real*            offset              /**< offset */
   )
{
   SCIP_Real* prize;
   SCIP_Real prizesum;
   int k;
   int e;
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

   prize = graph->prize;
   nnodes = graph->knots;
   nterms = graph->terms;
   prizesum = 0.0;

   stp_type = graph->stp_type;
   graph->stp_type = STP_SAP;

   /* for each terminal, except for the root, three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_copy(scip, graph, newgraph) );

   graph->stp_type = stp_type;

   SCIP_CALL( graph_resize(scip, (*newgraph), ((*newgraph)->ksize + 1), ((*newgraph)->esize + 2 * (nterms - 1)) , -1) );

   (*newgraph)->source[0] = graph->source[0];
   root = (*newgraph)->source[0];

   /* new pseudo-root */
   pseudoroot = (*newgraph)->knots;
   graph_knot_add((*newgraph), -1);

   for( k = 0; k < nnodes; k++ )
      if( Is_pterm(graph->term[k]) )
         prizesum += prize[k];

   prizesum += 1;

   *offset -= prizesum;

   e = (*newgraph)->outbeg[root];

   while( e != EAT_LAST )
   {
      enext = (*newgraph)->oeat[e];
      head = (*newgraph)->head[e];
      if( Is_term((*newgraph)->term[head]) )
      {
         (void) graph_edge_redirect(scip, (*newgraph), e, pseudoroot, head, graph->cost[e]);
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

   for( k = 0; k < nnodes; k++ )
   {
      /* is the kth node a terminal other than the root? */
      if( Is_pterm((*newgraph)->term[k]) )
      {
         graph_edge_add(scip, (*newgraph), k, pseudoroot, 0.0, FARAWAY);
      }
   }

   return SCIP_OKAY;
}

#if 1
/** alters the graph for prize collecting problems and shifts weights to reduce number of terminal */
SCIP_RETCODE graph_PcSapCopyShift(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Real*            offset              /**< offset */
   )
{
   GRAPH* newg;
   SCIP_Real p;
   SCIP_Real maxp;
   SCIP_Real* prize;
   SCIP_Real prizesum;
   int k;
   int e;
   int e2;
   int i;
   int l;
   int enext;
   int root;
   int head;
   int nnodes;
   int nterms;
   int maxvert;
   int stp_type;
   int pseudoroot;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->prize != NULL);
   assert(graph->knots == graph->ksize);
   assert(graph->edges == graph->esize);

   prize = graph->prize;
   nnodes = graph->knots;
   nterms = graph->terms;
   prizesum = 0.0;

   stp_type = graph->stp_type;
   graph->stp_type = STP_SAP;

   /* for each terminal, except for the root, three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_copy(scip, graph, newgraph) );

   graph->stp_type = stp_type;
   newg = *newgraph;

   maxvert = -1;
   maxp = -1.0;
   for( k = 0; k < nnodes; k++ )
   {
      if( Is_pterm(graph->term[k]) && graph->grad[k] > 0 )
      {
         if( SCIPisGT(scip, graph->prize[k], maxp) )
         {
            maxp = graph->prize[k];
            maxvert = k;
         }
      }
   }

   i = 0;
   for( k = 0; k < nnodes; k++ )
   {
      newg->mark[k] = (newg->grad[k] > 0);
      if( Is_pterm(graph->term[k]) && newg->mark[k] && k != maxvert )
      {
         assert(newg->grad[k] > 0);

         p = prize[k];
         for( e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
         {
            if( SCIPisLE(scip, graph->cost[e], p) && !Is_term(graph->term[graph->tail[e]]) )
               break;
         }

         /* if there is no incoming arc of lower cost than prize[k], make k a common node */
         if( e == EAT_LAST && i < graph->terms - 1)
         {

            newg->mark[k] = FALSE;
            newg->term[k] = -1;
            i++;
            e2 = UNKNOWN;
            head = UNKNOWN;
            for( e = graph->inpbeg[k]; e != EAT_LAST; e = graph->ieat[e] )
            {
               l = graph->tail[e];

               if( Is_term(graph->term[l]) )
               {
                  if( l == graph->source[0] )
                     e2 = e;
                  else
                     head = l;
               }
               else
               {
                  newg->cost[e] -= p;
                  assert(SCIPisGE(scip, newg->cost[e], 0.0));
               }
            }
            (*offset) += p;
            assert(e2 != UNKNOWN);
            assert(head != UNKNOWN);

            while( newg->inpbeg[head] != EAT_LAST )
               graph_edge_del(scip, newg, newg->inpbeg[head], FALSE);

            newg->mark[head] = FALSE;
            graph_knot_chg(newg, head, -1);
            graph_edge_del(scip, newg, e2, FALSE);
         }
      }
   }

   SCIP_CALL( graph_resize(scip, newg, (newg->ksize + 1), (newg->esize + 2 * (nterms - 1)) , -1) );

   assert(newg->source[0] == graph->source[0]);
   root = newg->source[0];

   /* new pseudo-root */
   pseudoroot = newg->knots;
   graph_knot_add(newg, -1);

   for( k = 0; k < nnodes; k++ )
      if( Is_pterm(graph->term[k]) && newg->mark[k] )
         prizesum += prize[k];

   prizesum += 1;

   *offset -= prizesum;

   e = newg->outbeg[root];

   while( e != EAT_LAST )
   {
      enext = newg->oeat[e];
      head = newg->head[e];
      if( Is_term(newg->term[head]) )
      {
         (void) graph_edge_redirect(scip, newg, e, pseudoroot, head, graph->cost[e]);
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

   for( k = 0; k < nnodes; k++ )
   {
      /* is the kth node a terminal other than the root? */
      if( Is_pterm(newg->term[k]) )
      {
         assert(newg->mark[k]);
         graph_edge_add(scip, newg, k, pseudoroot, 0.0, FARAWAY);
      }
   }

   return SCIP_OKAY;
}
#endif

/** alters the graph for prize-collecting problems with given root */
SCIP_RETCODE graph_PcRSapCopy(
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

   aterm = -1;
   proot = graph->source[0];
   stp_type = graph->stp_type;
   graph->stp_type = STP_SAP;

   /* copy graph */
   SCIP_CALL( graph_copy(scip, graph, newgraph) );

   p = *newgraph;

   graph->stp_type = stp_type;

   assert(Is_pterm(graph->term[root]));

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

   for( e = graph->outbeg[proot]; e != EAT_LAST; e = graph->oeat[e] )
   {
      head = graph->head[e];

      assert(graph->head[e] == p->head[e]);
      assert(graph->tail[e] == p->tail[e]);

      if( Is_term(graph->term[head]) && head != aterm )
      {
         assert(Is_term(p->term[head]));

         (void) graph_edge_redirect(scip, p, e, root, head, graph->cost[e]);
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

   assert(p->grad[aterm] == 0);

   nnodes = p->knots;
   p->source[0] = root;
   graph_knot_chg(p, root, 0);

   for( k = 0; k < nnodes; k++ )
      p->mark[k] = (p->grad[k] > 0);

   assert(p->grad[graph->source[0]] == 0);

   SCIP_CALL( SCIPallocMemoryArray(scip, &((*newgraph)->prize), nnodes) );

   for( k = 0; k < nnodes; k++)
   {
      if( k < graph->norgmodelknots )
      {
         p->prize[k] = graph->prize[k];
      }
      else
         p->prize[k] = 0.0;
   }

   if( nrootcands > 0 )
   {
      for( k = 0; k < nrootcands; k++ )
      {
         aterm = rootcands[k];
         if( aterm == root )
            continue;

         for( e = p->outbeg[aterm]; e != EAT_LAST; e = p->oeat[e] )
         {
            if( SCIPisZero(scip, p->cost[e]) )
            {
               head = p->head[e];

               for( e2 = p->inpbeg[head]; e2 != EAT_LAST; e2 = p->ieat[e] )
               {
                  if( p->tail[e2] == root )
                  {
                     p->cost[e2] = FARAWAY;
                     break;
                  }
               }
               assert(e2 != EAT_LAST);

               break;
            }
         }
         assert(e != EAT_LAST);
         p->prize[aterm] = FARAWAY;
      }
   }

   graph_knot_chg(p, proot, -1);
   p->prize[root] = 0.0;

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
   graph->source[0] = root;
   assert((nterms + 1) == graph->terms);
   if( graph->stp_type != STP_MWCSP )
      graph->stp_type = STP_PCSPG;

   return SCIP_OKAY;
}


/** alters the graph for rooted prize collecting problems */
SCIP_RETCODE graph_rootprize_transform(
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

   assert(graph != NULL);
   assert(graph->edges == graph->esize);

   root = graph->source[0];
   nnodes = graph->knots;
   nterms = graph->terms;
   prize = graph->prize;
   graph->norgmodeledges = graph->edges;
   graph->norgmodelknots = nnodes;

   assert(prize != NULL);
   assert(nnodes == graph->ksize);
   assert(root >= 0);

   /* for each terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + graph->terms), (graph->esize + graph->terms * 4) , -1) );

   /* create a new nodes */
   for( k = 0; k < nterms - 1; ++k )
      graph_knot_add(graph, -1);

   nterms = 0;

   for( k = 0; k < nnodes; ++k )
   {
      /* is the kth node a terminal other than the root? */
      if( Is_term(graph->term[k]) && k != root )
      {
         /* the copied node */
         node = nnodes + nterms;
         nterms++;
         /* switch the terminal property, mark k as former terminal */
         graph_knot_chg(graph, k, -2);
         graph_knot_chg(graph, node, 0);
         assert(SCIPisGE(scip, prize[k], 0.0));

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, node, prize[k], FARAWAY);
         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);
      }
      else
      {
         prize[k] = 0.0;
      }
   }
   /* one for the root */
   nterms++;
   assert((nterms) == graph->terms);
   graph->stp_type = STP_RPCSPG;
   return SCIP_OKAY;
}

/** alters the graph for MWCS problems */
SCIP_RETCODE graph_maxweight_transform(
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
            graph->cost[e] -= maxweights[i];
      }
      else if( SCIPisGT(scip, maxweights[i], 0.0) )
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
         assert(SCIPisLE(scip, maxweights[i], 0.0));
      }
   }
   assert(nterms == graph->terms);
   graph->stp_type = STP_MWCSP;

   SCIP_CALL( graph_prize_transform(scip, graph) );
   assert(graph->stp_type == STP_MWCSP);
   return SCIP_OKAY;
}




/** transforms MWCSP to RMWCSP if possible */
SCIP_RETCODE graph_MwToRmw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{

   int e;
   int k;
   int p;
   int e2;
   int root;
   int enext;
   int newroot;
   int maxgrad;

   assert(scip != NULL);
   assert(graph != NULL);

   newroot = -1;
   maxgrad = -1;
   root = graph->source[0];

   e = graph->outbeg[root];
   while( e != EAT_LAST )
   {
      enext = graph->oeat[e];
      if( SCIPisGE(scip, graph->cost[e], FARAWAY) )
      {
         k = graph->head[e];

         assert(Is_term(graph->term[k]));
         assert(graph->grad[k] == 2);

         for( e2 = graph->outbeg[k]; e2 != EAT_LAST; e2 = graph->oeat[e2] )
            if( graph->head[e2] != root )
               break;

         p = graph->head[e2];

         assert(Is_pterm(graph->term[p]));
         assert(SCIPisGE(scip, graph->prize[p], FARAWAY));

         /* delete terminal */
         graph_knot_chg(graph, k, -1);
         while( graph->outbeg[k] != EAT_LAST )
            graph_edge_del(scip, graph, graph->outbeg[k], TRUE);

         graph_knot_chg(graph, p, 0);

         if( graph->grad[p] > maxgrad )
         {
            newroot = p;
            maxgrad = graph->grad[p];
         }
      }
      e = enext;
   }

   /* is there a new root? */
   if( newroot >= 0 )
   {
      graph->source[0] = newroot;

      e = graph->outbeg[root];
      while( e != EAT_LAST )
      {
         enext = graph->oeat[e];
         k = graph->head[e];
         if( Is_term(graph->term[k]) && !SCIPisZero(scip, graph->cost[e]) )
         {
            (void) graph_edge_redirect(scip, graph, e, newroot, k, graph->cost[e]);
            graph->cost[flipedge(e)] = FARAWAY;
         }
         e = enext;
      }

      /* delete old root */
      graph_knot_chg(graph, root, -1);
      while( graph->outbeg[root] != EAT_LAST )
         graph_edge_del(scip, graph, graph->outbeg[root], TRUE);

      graph->stp_type = STP_RMWCSP;

      // todo
      printf("new problem type: STP_RMWCSP \n \n \n");
   }

   return SCIP_OKAY;
}



/** alters the graph for RMWCS problems */
SCIP_RETCODE graph_rootmaxweight_transform(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   SCIP_Real* maxweights;
   int e;
   int i;
   int k;
   int node;
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
         for( e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
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
   graph->source[0] = root;

   /* for each terminal, except for the root, one node and three edges (i.e. six arcs) are to be added */
   SCIP_CALL( graph_resize(scip, graph, (graph->ksize + npterms), (graph->esize + npterms * 4) , -1) );

   /* create a new nodes */
   for( k = 0; k < npterms; k++ )
      graph_knot_add(graph, -1);

   i = 0;
   for( k = 0; k < nnodes; ++k )
   {
      /* is the kth node a terminal other than the root? */
      if( Is_term(graph->term[k]) && SCIPisLT(scip, maxweights[k], FARAWAY) )
      {
         /* the copied node */
         node = nnodes + i;
         i++;

         /* switch the terminal property, mark k */
         graph_knot_chg(graph, k, -2);
         graph_knot_chg(graph, node, 0);
         assert(SCIPisGE(scip, maxweights[k], 0.0));

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(scip, graph, root, node, maxweights[k], FARAWAY);
         graph_edge_add(scip, graph, k, node, 0.0, FARAWAY);
      }
   }

   assert(i == npterms);

   graph->stp_type = STP_RMWCSP;

   assert(graph->stp_type == STP_RMWCSP);
   return SCIP_OKAY;
}


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

/** free the graph */
void graph_free(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                p,                  /**< graph to be freed */
   SCIP_Bool             final               /**< delete ancestor data structures? */
   )
{
   IDX* curr;
   int e;
   int i;
   int nedges;

   assert(scip != NULL);
   assert(p != NULL);

   nedges = p->edges;

   if( p->ancestors != NULL )
   {
        for( e = nedges - 1; e >= 0; e-- )
        {
           curr = p->ancestors[e];
           while( curr != NULL )
           {
              p->ancestors[e] = curr->parent;
              SCIPfreeMemory(scip, &(curr));
              curr = p->ancestors[e];
           }
        }
        SCIPfreeMemoryArray(scip, &(p->ancestors));
   }

   if( final )
   {
      assert(p->path_heap == NULL);
      assert(p->path_state == NULL);

      if( p->pcancestors != NULL )
      {
         for( e = p->norgmodelknots - 1; e >= 0; e-- )
         {
            curr = p->pcancestors[e];
            while( curr != NULL )
            {
               p->pcancestors[e] = curr->parent;
               SCIPfreeMemory(scip, &(curr));
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
      curr = p->fixedges;
      while( curr != NULL )
      {
         p->fixedges = curr->parent;
         SCIPfreeMemory(scip, &(curr));

         curr = p->fixedges;
      }
   }

   if( p->prize != NULL )
      SCIPfreeMemoryArray(scip, &(p->prize));

   if( p->stp_type == STP_DCSTP )
   {
      SCIPfreeMemoryArray(scip, &(p->maxdeg));
   }
   else if( p->stp_type == STP_RSMT )
   {
      if( p->grid_coordinates != NULL )
      {
         assert(p->grid_coordinates != NULL);
         for( i = p->grid_dim - 1; i >= 0;  i-- )
            SCIPfreeMemoryArray(scip, &(p->grid_coordinates[i]));

         SCIPfreeMemoryArray(scip, &(p->grid_coordinates));
      }

      if( p->grid_ncoords != NULL )
         SCIPfreeMemoryArray(scip, &(p->grid_ncoords));
   }

#if BLK
   SCIPfreeBlockMemoryArray(scip, &(p->oeat), p->esize);
   SCIPfreeBlockMemoryArray(scip, &(p->ieat), p->esize);
   SCIPfreeBlockMemoryArray(scip, &(p->head), p->esize);
   SCIPfreeBlockMemoryArray(scip, &(p->tail), p->esize);
   SCIPfreeBlockMemoryArray(scip, &(p->cost), p->esize);
   SCIPfreeBlockMemoryArray(scip, &(p->outbeg), p->ksize);
   SCIPfreeBlockMemoryArray(scip, &(p->inpbeg), p->ksize);
   SCIPfreeBlockMemoryArray(scip, &(p->grad), p->ksize);
   SCIPfreeBlockMemoryArray(scip, &(p->mark), p->ksize);
   SCIPfreeBlockMemoryArray(scip, &(p->term), p->ksize);
   SCIPfreeBlockMemoryArray(scip, &(p->source), p->layers);
   SCIPfreeBlockMemoryArray(scip, &(p->locals), p->layers);

   SCIPfreeBlockMemory(scip, &(p));
#else
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
   SCIPfreeMemoryArray(scip, &(p->source));
   SCIPfreeMemoryArray(scip, &(p->locals));

   SCIPfreeMemory(scip, &(p));
#endif
}

/** copy the graph */
SCIP_RETCODE graph_copy(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          orgraph,            /**< original graph */
   GRAPH**               copygraph           /**< graph to be copied */
   )
{
   GRAPH* g;
   const GRAPH* p;
   int k;
   int ksize;

   p = orgraph;
   assert(p != NULL);

   ksize = p->ksize;

   SCIP_CALL( graph_init(scip, copygraph, p->ksize, p->esize, p->layers, p->flags) );
   g = *copygraph;

   g->norgmodeledges = p->norgmodeledges;
   g->norgmodelknots = p->norgmodelknots;
   g->knots = p->knots;
   g->terms = p->terms;
   g->edges = p->edges;
   g->orgsource = p->orgsource;
   g->orgedges = p->orgedges;
   g->orgknots = p->orgknots;
   g->grid_dim = p->grid_dim;
   g->stp_type = p->stp_type;
   g->hoplimit = p->hoplimit;

   BMScopyMemoryArray(g->locals, p->locals, p->layers);
   BMScopyMemoryArray(g->source, p->source, p->layers);
   BMScopyMemoryArray(g->term, p->term, ksize);
   BMScopyMemoryArray(g->mark, p->mark, ksize);
   BMScopyMemoryArray(g->grad, p->grad, ksize);
   BMScopyMemoryArray(g->inpbeg, p->inpbeg, ksize);
   BMScopyMemoryArray(g->outbeg, p->outbeg, ksize);
   BMScopyMemoryArray(g->cost, p->cost, p->esize);
   BMScopyMemoryArray(g->tail, p->tail, p->esize);
   BMScopyMemoryArray(g->head, p->head, p->esize);
   BMScopyMemoryArray(g->ieat, p->ieat, p->esize);
   BMScopyMemoryArray(g->oeat, p->oeat, p->esize);

   if( g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG || g->stp_type == STP_MWCSP || g->stp_type == STP_RMWCSP )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->prize), g->knots) );

      for( k = 0; k < g->knots; k++)
         if( Is_term(p->term[k]) )
            g->prize[k] = 0.0;
         else
            g->prize[k] = p->prize[k];
   }
   else if( g->stp_type == STP_DCSTP )
   {
      assert(p->maxdeg != NULL);

      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->maxdeg), g->knots) );

      for( k = 0; k < g->knots; k++)
         g->maxdeg[k] = p->maxdeg[k];
   }
   else if( p->stp_type == STP_RSMT )
   {
      assert(p->grid_ncoords != NULL);
      assert(p->grid_coordinates != NULL);

      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->grid_coordinates), p->grid_dim) );

      BMScopyMemoryArray(g->grid_coordinates, p->grid_coordinates, p->grid_dim);
      for( k = 0; k < p->grid_dim; k++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(g->grid_coordinates[k]), p->terms) ); /*lint !e866*/
         BMScopyMemoryArray(g->grid_coordinates[k], p->grid_coordinates[k], p->terms); /*lint !e866*/
      }
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->grid_ncoords), p->grid_dim) );

      BMScopyMemoryArray(g->grid_ncoords, p->grid_ncoords, p->grid_dim);
   }
   assert(graph_valid(g));

   return SCIP_OKAY;
}

/** set flags */
void graph_flags(
   GRAPH*                p,                  /**< the graph */
   int                   flags               /**< new flags */
   )
{
   assert(p     != NULL);
   assert(flags >= 0);

   p->flags |= flags;
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

void graph_ident(
   const GRAPH*          p                   /**< the graph */
   )
{
   int i;
   int ident = 0;

   assert(p != NULL);

   for(i = 0; i < p->knots; i++)
      ident += (i + 1) * (p->term[i] * 2 + p->grad[i] * 3
         + p->inpbeg[i] * 5 + p->outbeg[i] * 7);

   for(i = 0; i < p->edges; i++)
      ident += (i + 1) * ((int)p->cost[i] + p->tail[i]
         + p->head[i] + p->ieat[i] + p->oeat[i]);

   (void)printf("Graph Ident = %d\n", ident);
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

   if (Is_term(term))
   {
      p->terms++;
      p->locals[term]++;
   }
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
      {
         p->terms--;
         p->locals[p->term[node]]--;
      }
      p->term[node] = term;

      if( Is_term(p->term[node]) )
      {
         p->terms++;
         p->locals[p->term[node]]++;
      }
   }
}

/** contract an edge, given by its endpoints */
SCIP_RETCODE graph_knot_contract(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                p,                  /**< graph data structure */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                or NULL */
   int                   t,                  /**< tail node to be contracted */
   int                   s                   /**< head node to be contracted */
   )
{
   SCIP_Real* incost = NULL;
   SCIP_Real* outcost = NULL;
   IDX**   ancestors = NULL;
   IDX**   revancestors = NULL;
   IDX*   tsancestors = NULL;
   IDX*   stancestors = NULL;
   int* mark = NULL;
   int* edge = NULL;
   int* knot = NULL;
   int    slc = 0;
   int    i;
   int    et;
   int    anti;
   int    es;
   int    cedgeout = UNKNOWN;
   int    head;
   int    tail;
   int    sgrad;
   SCIP_Bool conflict;

   assert(p          != NULL);
   assert(t          >= 0);
   assert(t          <  p->knots);
   assert(s          >= 0);
   assert(s          <  p->knots);
   assert(s          != t);
   assert(scip          != NULL);
   assert(p->grad[s] >  0);
   assert(p->grad[t] >  0);
   assert(p->layers  == 1);

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
   if( p->source[0] == s )
      p->source[0] = t;

   sgrad = p->grad[s];
   if( sgrad >= 2 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &incost, sgrad - 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &outcost, sgrad - 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &mark, sgrad - 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &edge, sgrad - 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &knot, sgrad - 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ancestors, sgrad - 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &revancestors, sgrad - 1) );
   }

   /* store edges to be moved/removed */
   for( es = p->outbeg[s]; es != EAT_LAST; es = p->oeat[es] )
   {
      assert(p->tail[es] == s);

      if( p->head[es] != t )
      {
         assert(ancestors != NULL);
         assert(revancestors != NULL);
         assert(mark != NULL);
         assert(incost != NULL);
         assert(outcost != NULL);
         assert(edge != NULL);
         assert(knot != NULL);

         ancestors[slc] = NULL;
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[slc]), p->ancestors[es], NULL) );
         revancestors[slc] = NULL;
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[slc]), p->ancestors[Edge_anti(es)], NULL) );

         mark[slc] = FALSE;
         edge[slc] = es;
         knot[slc] = p->head[es];
         outcost[slc] = p->cost[es];
         incost[slc] = p->cost[Edge_anti(es)];
         slc++;

         assert(slc < sgrad);
      }
      else
      {
         cedgeout = Edge_anti(es); /* The edge out of t and into s. */
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &stancestors, p->ancestors[es], NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &tsancestors, p->ancestors[cedgeout], NULL) );
      }
   }

   assert(slc == sgrad - 1);
   assert(tsancestors != NULL);
   assert(stancestors != NULL);

   /* traverse edges */
   for( i = 0; i < slc; i++ )
   {
      assert(knot != NULL && outcost != NULL && incost != NULL && mark != NULL);

      /* search for an edge out of t with same head as current edge */
      for( et = p->outbeg[t]; et != EAT_LAST; et = p->oeat[et] )
         if( p->head[et] == knot[i] )
            break;

      /* does such an edge not exist? */
      if( et == EAT_LAST )
      {
         mark[i] = TRUE;
      }
      else
      {
         assert(et != EAT_LAST);

         /* This is for nodes with edges to s and t.
          * Need to adjust the out and in costs of the edge
          */
         if( SCIPisGT(scip, p->cost[et], outcost[i]) )
         {
            SCIPintListNodeFree(scip, &((p->ancestors)[et]));
            assert(ancestors != NULL);
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[et]), ancestors[i], NULL) );
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[et]), tsancestors, NULL) );
            p->cost[et] = outcost[i];
         }
         if( SCIPisGT(scip, p->cost[Edge_anti(et)], incost[i]) )
         {
            anti = Edge_anti(et);
            SCIPintListNodeFree(scip, &(p->ancestors[anti]));
            assert(revancestors != NULL);
            assert(incost != NULL);
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[anti]), revancestors[i], NULL) );
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[anti]), stancestors, NULL) );
            p->cost[anti] = incost[i];
         }
      }
   }

   /* insert edges */
   for( i = 0; i < slc; i++ )
   {
      assert(mark != NULL);
      if( mark[i] )
      {
         es = p->outbeg[s];

         assert(es != EAT_LAST);
         assert(ancestors != NULL);
         assert(revancestors != NULL);
         assert(ancestors[i] != NULL);
         assert(revancestors[i] != NULL);
         assert(knot != NULL);
         assert(outcost != NULL);
         assert(incost != NULL);
         SCIPintListNodeFree(scip, &(p->ancestors[es]));
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), ancestors[i], NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), tsancestors, &conflict) );

         if( conflict && 0 ) /* todo */
         {
            graph_edge_del(scip, p, es, TRUE);
            printf("conflict %d \n", es);
            continue;
         }
         else
         {
            graph_edge_del(scip, p, es, FALSE);
         }

         head = knot[i];
         tail = t;

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
         SCIPintListNodeFree(scip, &(p->ancestors[es]));

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), revancestors[i], NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), stancestors, NULL) );
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
   while( p->outbeg[s] != EAT_LAST )
   {
      es = p->outbeg[s];
      SCIPintListNodeFree(scip, &(p->ancestors[es]));
      SCIPintListNodeFree(scip, &(p->ancestors[Edge_anti(es)]));
      graph_edge_del(scip, p, es, FALSE);
   }

   SCIPintListNodeFree(scip, &stancestors);
   SCIPintListNodeFree(scip, &tsancestors);

   if( sgrad >= 2 )
   {
      assert(ancestors != NULL);
      assert(revancestors != NULL);
      for( i = 0; i < slc; i++ )
      {
         SCIPintListNodeFree(scip, &(ancestors[i]));
         SCIPintListNodeFree(scip, &(revancestors[i]));
      }
      SCIPfreeBufferArray(scip, &revancestors);
      SCIPfreeBufferArray(scip, &ancestors);
      SCIPfreeBufferArray(scip, &knot);
      SCIPfreeBufferArray(scip, &edge);
      SCIPfreeBufferArray(scip, &mark);
      SCIPfreeBufferArray(scip, &outcost);
      SCIPfreeBufferArray(scip, &incost);
   }
   assert(p->grad[s]   == 0);
   assert(p->outbeg[s] == EAT_LAST);
   assert(p->inpbeg[s] == EAT_LAST);
   return SCIP_OKAY;
}

/** subtract a given sum from the prize of a terminal */
void prize_subtract(
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

   if( g->stp_type == STP_RPCSPG && i == g->source[0] )
      return;

   g->prize[i] -= cost;
   for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      if( Is_pterm(g->term[g->head[e]]) )
         break;

   assert(e != EAT_LAST);

   j = g->head[e];

   assert(j != g->source[0]);
   assert(!g->mark[j]);

   for( e = g->inpbeg[j]; e != EAT_LAST; e = g->ieat[e] )
      if( g->source[0] == g->tail[e] )
         break;

   assert(e != EAT_LAST);
   assert(!g->mark[g->tail[e]] || g->stp_type == STP_RPCSPG);

   g->cost[e] -= cost;

   assert(g->stp_type == STP_MWCSP  || g->stp_type == STP_RMWCSP || SCIPisGE(scip, g->prize[i], 0.0));
   assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));
   assert(SCIPisGE(scip, g->prize[i], 0.0) || g->stp_type == STP_MWCSP);
}

/** contract an edge of (rooted) prize-collecting Steiner tree problem or maximum-weight connected subgraph problem */
SCIP_RETCODE graph_knot_contractpc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int*                  solnode,            /**< solution nodes or NULL */
   int                   t,                  /**< tail node to be contracted */
   int                   s,                  /**< head node to be contracted */
   int                   i                   /**< terminal to add offset to */
   )
{
   int ets;

   assert(g != NULL);
   assert(scip != NULL);
   assert(Is_term(g->term[i]));

   /* get edge from t to s */
   for( ets = g->outbeg[t]; ets != EAT_LAST; ets = g->oeat[ets] )
      if( g->head[ets] == s )
         break;

   assert(ets != EAT_LAST);

   if( g->pcancestors[s] != NULL )
   {
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[t]), g->pcancestors[s], NULL) );
      SCIPintListNodeFree(scip, &(g->pcancestors[s]));
   }

   SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[t]), g->ancestors[ets], NULL) );

   /* are both endpoints of the edge to be contracted terminals? */
   if( Is_term(g->term[t]) && Is_term(g->term[s]) )
   {
      int e;
      int j;

      /* get edge from s to its artificial terminal */
      for( e = g->outbeg[s]; e != EAT_LAST; e = g->oeat[e] )
         if( Is_pterm(g->term[g->head[e]]) )
            break;

      assert(e != EAT_LAST);
      assert(g->pcancestors != NULL);

      /* artificial terminal to s */
      j = g->head[e];

      assert(j != g->source[0]);
      assert(!g->mark[j]);

      /* delete edge and unmark artificial terminal */
      graph_knot_chg(g, j, -1);
      graph_edge_del(scip, g, e, TRUE);

      /* delete remaining incident edge of artificial terminal */
      e = g->inpbeg[j];

      assert(e != EAT_LAST);
      assert(g->source[0] == g->tail[e] || g->source[0] == j);
      assert(SCIPisEQ(scip, g->prize[s], g->cost[e]));

      prize_subtract(scip, g, g->cost[ets] - g->prize[s], i);
      graph_edge_del(scip, g, e, TRUE);

      assert(g->inpbeg[j] == EAT_LAST);

      /* contract s into t */
      SCIP_CALL( graph_knot_contract(scip, g, solnode, t, s) );
      graph_knot_chg(g, s, -1);

      assert(g->grad[s] == 0);

      SCIPdebugMessage("PC contract: %d, %d \n", t, s);
   }
   else
   {
      if( g->stp_type != STP_MWCSP && g->stp_type != STP_RMWCSP )
         prize_subtract(scip, g, g->cost[ets], i);
      else
         prize_subtract(scip, g, -(g->prize[s]), i);
      SCIP_CALL( graph_knot_contract(scip, g, solnode, t, s) );
   }
   return SCIP_OKAY;
}

int graph_edge_redirect(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   eki,                /**< the edge */
   int                   k,                  /**< new tail */
   int                   j,                  /**< new head */
   SCIP_Real             cost                /**< new cost */
   )
{
   int e;

   graph_edge_del(NULL, g, eki, FALSE);

   for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      if( (g->tail[e] == k) && (g->head[e] == j) )
         break;

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
   int                   e1,                 /**< edge to reinsert */
   int                   k1,                 /**< tail */
   int                   k2,                 /**< head */
   SCIP_Real             cost,               /**< edgecost */
   IDX*                  ancestors0,         /**< ancestors of first edge */
   IDX*                  ancestors1,         /**< ancestors of second edge */
   IDX*                  revancestors0,      /**< reverse ancestors of first edge */
   IDX*                  revancestors1       /**< reverse ancestors of first edge */
   )
{
   int n1;

   /* redirect; store new edge in n1 */
   n1 = graph_edge_redirect(scip, g, e1, k1, k2, cost);
   if( n1 >= 0 )
   {
      SCIPintListNodeFree(scip, &(g->ancestors[n1]));
      SCIPintListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));

      SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), revancestors0, NULL) );
      SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors1, NULL) );

      SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), ancestors0, NULL) );
      SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors1, NULL) );
   }
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

inline static void edge_remove(
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
      for( i = g->inpbeg[head]; g->ieat[i] != e; i = g->ieat[i] )
         assert(i >= 0);

      g->ieat[i] = g->ieat[e];
   }
   if( g->outbeg[tail] == e )
      g->outbeg[tail] = g->oeat[e];
   else
   {
      for( i = g->outbeg[tail]; g->oeat[i] != e; i = g->oeat[i] )
         assert(i >= 0);

      g->oeat[i] = g->oeat[e];
   }
}

/** delete an edge */
void graph_edge_del(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   e,                  /**< the edge */
   SCIP_Bool             freeancestors       /**< free edge ancestors? */
   )
{
   assert(g          != NULL);
   assert(e          >= 0);
   assert(e          <  g->edges);

   if( freeancestors )
   {
      assert(scip != NULL);
      SCIPintListNodeFree(scip, &((g->ancestors)[e]));
      SCIPintListNodeFree(scip, &((g->ancestors)[Edge_anti(e)]));
   }

   /* delete first arc */
   e -= e % 2;
   assert(g->head[e] == g->tail[e + 1]);
   assert(g->tail[e] == g->head[e + 1]);

   g->grad[g->head[e]]--;
   g->grad[g->tail[e]]--;

   edge_remove(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_FREE;
   g->oeat[e] = EAT_FREE;

   /* delete second arc */
   e++;
   edge_remove(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_FREE;
   g->oeat[e] = EAT_FREE;
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

   edge_remove(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_HIDE;
   g->oeat[e] = EAT_HIDE;

   e++;

   edge_remove(g, e);

   assert(g->ieat[e] != EAT_FREE);
   assert(g->ieat[e] != EAT_HIDE);
   assert(g->oeat[e] != EAT_FREE);
   assert(g->oeat[e] != EAT_HIDE);

   g->ieat[e] = EAT_HIDE;
   g->oeat[e] = EAT_HIDE;
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

/** mark terminals and switch terminal property to original terminals */
SCIP_RETCODE pcgraphorg(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   int k;
   int  root;
   int nnodes;

   assert(scip != NULL);
   assert(graph != NULL);

   root = graph->source[0];
   nnodes = graph->knots;

   for( k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);

      if( Is_pterm(graph->term[k]) )
      {
         graph_knot_chg(graph, k, 0);
      }
      else if( Is_term(graph->term[k]) )
      {
         graph->mark[k] = FALSE;
         if( k != root )
            graph_knot_chg(graph, k, -2);
      }
   }

   if( graph->stp_type == STP_RPCSPG )
      graph->mark[root] = TRUE;

   return SCIP_OKAY;
}

/** unmark terminals and switch terminal property to transformed terminals */
SCIP_RETCODE pcgraphtrans(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< the graph */
   )
{
   int k;
   int  root;
   int nnodes;

   assert(scip != NULL);
   assert(graph != NULL);

   root = graph->source[0];
   nnodes = graph->knots;

   for( k = 0; k < nnodes; k++ )
   {
      graph->mark[k] = (graph->grad[k] > 0);

      if( Is_pterm(graph->term[k]) )
         graph_knot_chg(graph, k, 0);
      else if( Is_term(graph->term[k]) && k != root )
         graph_knot_chg(graph, k, -2);
   }

   return SCIP_OKAY;
}

/** pack the graph, i.e. build a new graph that discards deleted edges and nodes */
SCIP_RETCODE graph_pack(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Bool             verbose             /**< verbose? */
   )
{
   const char* msg1 = "Nodes: %d  Edges: %d  Terminals: %d\n";
   GRAPH* g;
   GRAPH* q;
   int*   new;
   int    i;
   int    e;
   int    oldnnodes;
   int    oldnedges;
   int    nnodes;
   int    nedges;
   SCIP_Bool rmw;
   SCIP_Bool pcmw;

   assert(scip      != NULL);
   assert(graph     != NULL);
   assert(graph_valid(graph));

   g = graph;
   nnodes = 0;
   nedges = 0;
   oldnnodes = g->knots;
   oldnedges = g->edges;
   SCIP_CALL( SCIPallocBufferArray(scip, &new, oldnnodes) );

   if( verbose )
      printf("Reduced graph: ");

   /* count nodes */
   for( i = 0; i < oldnnodes; i++ )
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
   for( i = 0; i < oldnedges; i++ )
   {
      if( g->oeat[i] != EAT_FREE )
      {
         assert(g->ieat[i] != EAT_FREE);
         nedges++;
      }
   }

   assert(nnodes > 1 || nedges == 0);
   SCIP_CALL( graph_init(scip, newgraph, nnodes, nedges, g->layers, g->flags) );
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
   q->fixedges = g->fixedges;
   q->hoplimit = g->hoplimit;
   q->pcancestors = g->pcancestors;

   if( new == NULL )
   {
      q->ancestors = NULL;
      graph_free(scip, g, FALSE);

      if( q->stp_type == STP_RSMT )
      {
         q->grid_ncoords = NULL;
         q->grid_coordinates = NULL;
      }

      graph_knot_add(q, 0);
      q->source[0] = 0;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(q->ancestors), nedges) );

   rmw = g->stp_type == STP_RMWCSP;
   pcmw = (g->stp_type == STP_MWCSP || g->stp_type == STP_RPCSPG || g->stp_type == STP_PCSPG || g->stp_type == STP_RMWCSP);
   if( pcmw )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(q->prize), nnodes) );
   }

   /* add nodes (of positive degree) */

   if( rmw )
   {
      for( i = 0; i < oldnnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

      for( e = g->outbeg[g->source[0]]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( SCIPisGT(scip, g->cost[e], 0.0) && Is_term(g->term[g->head[e]]) )
         {
            i = g->head[e];
            g->mark[i] = FALSE;
            assert(g->grad[i] == 2);
         }
      }
   }

   for( i = 0; i < oldnnodes; i++ )
   {
      assert(g->term[i] < g->layers);
      if( g->grad[i] > 0 )
      {
         if( pcmw )
         {
            if( !Is_term(g->term[i]) || (rmw && g->mark[i]) )
               q->prize[q->knots] = g->prize[i];
            else
               q->prize[q->knots] = 0.0;
         }
         graph_knot_add(q, g->term[i]);
      }
   }

   /* add edges */
   for( i = 0; i < oldnedges; i += 2 )
   {
      if( g->ieat[i] == EAT_FREE )
      {
         assert(g->oeat[i]     == EAT_FREE);
         assert(g->ieat[i + 1] == EAT_FREE);
         assert(g->oeat[i + 1] == EAT_FREE);
         SCIPintListNodeFree(scip, &(g->ancestors[i]));
         SCIPintListNodeFree(scip, &(g->ancestors[i + 1]));
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

      graph_edge_add(scip, q, new[g->tail[i]], new[g->head[i]], g->cost[i], g->cost[Edge_anti(i)]);
   }

   /* add root */
   assert(q->term[new[g->source[0]]] == 0);

   q->source[0] = new[g->source[0]];

   if( g->stp_type == STP_RPCSPG )
      q->prize[q->source[0]] = FARAWAY;

   SCIPfreeBufferArray(scip, &new);

   g->stp_type = UNKNOWN;
   graph_free(scip, g, FALSE);

   assert(graph_valid(q));

   if( verbose )
      printf(msg1, q->knots, q->edges, q->terms);

   return SCIP_OKAY;
}


/** get (real) number of nodes , edges, terminals */
void graph_getNVET(
   const GRAPH*    graph,              /**< the graph */
   int*            nnodes,             /**< number of nodes */
   int*            nedges,             /**< number of edges */
   int*            nterms              /**< number of terminals */
   )
{
   int k;
   int v = 0;
   int e = 0;
   int t = 0;
   int vorg;

   assert(graph != NULL);

   vorg = graph->knots;

   for( k = 0; k < vorg; k++ )
   {
      if( graph->grad[k] > 0 )
      {
         v++;
         e += graph->grad[k];
         if( Is_term(graph->term[k]) )
            t++;
      }
   }

   *nnodes = v;
   *nedges = e;
   *nterms = t;

   return;
}

/** traverse the graph and mark all reached nodes (g->mark[i] has to be FALSE for all i) */
void graph_trail(
   const GRAPH*          g,                  /**< the new graph */
   int                   i                   /**< node to start from */
   )
{
   int* gmark;

   assert(g      != NULL);
   assert(i      >= 0);
   assert(i      <  g->knots);

   gmark = g->mark;

   if( !gmark[i] )
   {
      SCIP_QUEUE* queue;
      int a;
      int head;
      int node;
      int* pnode;

      gmark[i] = TRUE;

      if( g->grad[i] == 0 )
         return;

      SCIP_CALL_ABORT( SCIPqueueCreate(&queue, g->knots, 1.1) );
      SCIP_CALL_ABORT( SCIPqueueInsert(queue, &i));

      /* BFS loop */
      while( !SCIPqueueIsEmpty(queue) )
      {
         pnode = (SCIPqueueRemove(queue));
         node = *pnode;

         /* traverse outgoing arcs */
         for( a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
         {
            head = g->head[a];

            if( !gmark[head] )
            {
               gmark[head] = TRUE;
               SCIP_CALL_ABORT(SCIPqueueInsert(queue, &(g->head[a])));
            }
         }
      }
      SCIPqueueFree(&queue);
   }
}


/** traverse the graph and mark all reached nodes (g->mark[i] has to be FALSE for all i) .... uses an array and should be faster
 *  than graph_trail, but needs a scip */
SCIP_RETCODE graph_trail_arr(
   SCIP*                 scip,               /**< scip struct */
   const GRAPH*          g,                  /**< the new graph */
   int                   i                   /**< node to start from */
   )
{
   int* gmark;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(i      >= 0);
   assert(i      <  g->knots);

   gmark = g->mark;

   if( !gmark[i] )
   {
      int* stackarr;
      int a;
      int head;
      int node;
      int nnodes;
      int stacksize;

      gmark[i] = TRUE;

      if( g->grad[i] == 0 )
         return SCIP_OKAY;;

      nnodes = g->knots;
      stacksize = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &stackarr, nnodes) );

      stackarr[stacksize++] = i;

      /* BFS loop */
      while( stacksize != 0 )
      {
         node = stackarr[--stacksize];

         /* traverse outgoing arcs */
         for( a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
         {
            head = g->head[a];

            if( !gmark[head] )
            {
               gmark[head] = TRUE;
               stackarr[stacksize++] = head;
            }
         }
      }
      SCIPfreeBufferArray(scip, &stackarr);
   }
   return SCIP_OKAY;
}

/** is the given graph valid? */
int graph_valid(
   const GRAPH*          g                   /**< the new graph */
   )
{
   const char* fehler1  = "*** Graph Validation Error: Head invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n";
   const char* fehler2  = "*** Graph Validation Error: Tail invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n";
   const char* fehler3  = "*** Graph Validation Error: Source invalid, Layer %d, Source %d, Terminal %d\n";
   const char* fehler4  = "*** Graph Validation Error: FREE invalid, Edge %d/%d\n";
   const char* fehler5  = "*** Graph Validation Error: Anti invalid, Edge %d/%d, Tail=%d/%d, Head=%d/%d\n";
   const char* fehler6  = "*** Graph Validation Error: Knot %d with Grad 0 has Edges\n";
   const char* fehler7  = "*** Graph Validation Error: Knot %d not connected\n";
#if 0
   const char* fehler8  = "*** Graph Validation Error: Wrong locals count, Layer %d, count is %d, should be %d\n";
#endif
   const char* fehler9  = "*** Graph Validation Error: Wrong Terminal count, count is %d, should be %d\n";

   int    k;
   int    l;
   int    e;
   int    nterms;
   int    nnodes;
   int    nedges;

   assert(g != NULL);

   nterms = g->terms;
   nedges = g->edges;
   nnodes = g->knots;

   for( k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) )
      {
         nterms--;
      }
      for( e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         if( g->head[e] != k )
            break;

      if( e != EAT_LAST )
         return((void)fprintf(stderr, fehler1, k, e, g->tail[e], g->head[e]), FALSE);

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         if( g->tail[e] != k )
            break;

      if( e != EAT_LAST )
         return((void)fprintf(stderr, fehler2, k, e, g->tail[e], g->head[e]), FALSE);
   }
   if( nterms != 0 )
      return((void)fprintf(stderr, fehler9, g->terms, g->terms - nterms), FALSE);

   for( l = 0; l < g->layers; l++ )
   {
      if( (g->source[l] < 0 )
         || (g->source[l] >= g->knots)
         || (g->term[g->source[l]] != l))
         return((void)fprintf(stderr, fehler3,
               l, g->source[l], g->term[g->source[l]]), FALSE);
   }

   for( e = 0; e < nedges; e += 2 )
   {
      if( (g->ieat[e] == EAT_FREE) && (g->oeat[e] == EAT_FREE)
         && (g->ieat[e + 1] == EAT_FREE) && (g->oeat[e + 1] == EAT_FREE) )
         continue;

      if( (g->ieat[e] == EAT_FREE) || (g->oeat[e] == EAT_FREE)
         || (g->ieat[e + 1] == EAT_FREE) || (g->oeat[e + 1] == EAT_FREE) )
         return((void)fprintf(stderr, fehler4, e, e + 1), FALSE);

      if( (g->head[e] != g->tail[e + 1]) || (g->tail[e] != g->head[e + 1]) )
         return((void)fprintf(stderr, fehler5,
               e, e + 1, g->head[e], g->tail[e + 1],
               g->tail[e], g->head[e + 1]), FALSE);
   }

   for( k = 0; k < nnodes; k++ )
      g->mark[k] = FALSE;

   graph_trail(g, g->source[0]);

   for( k = 0; k < nnodes; k++ )
   {
      if( (g->grad[k] == 0)
         && ((g->inpbeg[k] != EAT_LAST) || (g->outbeg[k] != EAT_LAST)) )
         return((void)fprintf(stderr, fehler6, k), FALSE);

      if( !g->mark[k] && ((g->grad[k] > 0) || (Is_term(g->term[k])))
         && g->stp_type != STP_PCSPG && g->stp_type != STP_MWCSP && g->stp_type != STP_RMWCSP )
         return((void)fprintf(stderr, fehler7, k), FALSE);
   }
   return TRUE;
}

/** verifies whether a given primal solution is feasible */
SCIP_Bool graph_sol_valid(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int*                  result              /**< solution array, indicating whether an edge is in the solution */
   )
{
   SCIP_QUEUE* queue;

   STP_Bool* terminal;
   int* pnode;
   int e;
   int i;
   int root;
   int nnodes;
   int termcount;

   assert(graph != NULL);
   assert(result != NULL);

   terminal = NULL;
   nnodes = graph->knots;
   root = graph->source[0];
   assert(root >= 0);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &terminal, nnodes) );

   assert(terminal != NULL);

   for( i = 0; i < nnodes; i++ )
      terminal[i] = FALSE;

   /* BFS until all terminals are reached */
   SCIP_CALL_ABORT( SCIPqueueCreate(&queue, nnodes, 2.0) );

   SCIP_CALL_ABORT( SCIPqueueInsert(queue, &root) );
   termcount = 1;
   terminal[root] = TRUE;

   while( !SCIPqueueIsEmpty(queue) )
   {
      pnode = (SCIPqueueRemove(queue));
      for( e = graph->outbeg[*pnode]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( result[e] == CONNECT )
         {
            i = graph->head[e];
            if( Is_term(graph->term[i]) )
            {
               assert(!terminal[i]);
               terminal[i] = TRUE;
               termcount++;
            }
            SCIP_CALL_ABORT( SCIPqueueInsert(queue, &graph->head[e]) );
         }
      }
   }

#if 0
   if(termcount != graph->terms)
   {
      for( i = 0; i < nnodes; i++ )
      {
         if( Is_term(graph->term[i]) && !terminal[i] )
         {
            printf("root %d \n", root);
            printf("fail %d grad %d\n", i, graph->grad[i]);
            for( e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
            {
               printf("tail %d %d \n", graph->tail[e], graph->term[graph->tail[e]]);
            }
         }
      }
   }
#endif



   SCIPqueueFree(&queue);
   SCIPfreeBufferArray(scip, &terminal);

   return (termcount == graph->terms);
}

/** mark endpoints of edges in given list */
void graph_listSetSolNode(
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
SCIP_Real graph_computeSolVal(
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
