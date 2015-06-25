/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   grphbase.c
 * @brief  includes several methodes for Steiner problem graphs
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file contains several basic methods to process Steiner problem graphs.
 * A graph can not be reduced in terms of edge or node size, but edges can be marked as
 * EAT_FREE (to not be used anymore) and nodes may have degree one.
 * The method 'graph_pack()' can be used to build a new graph, discarding those nodes and edges
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
#include "grph.h"
#if 0

#define GRPHBASE_C -esym(750,GRPHBASE_C)
GRAPH* graph_init2(
   int ksize,
   int esize,
   int layers,
   int flags)
{
   GRAPH* p;
   int    i;

   assert(ksize > 0);
   assert(ksize < INT_MAX);
   assert(esize >= 0);
   assert(esize < INT_MAX);
   assert(layers > 0);
   assert(layers < SHRT_MAX);

   p = malloc(sizeof(*p));

   assert(p != NULL);

   p->fixedges = NULL;
   p->ancestors = NULL;

   p->norgmodelknots = 0;
   p->norgmodeledges = 0;
   p->ksize  = ksize;
   p->orgknots = 0;
   p->orgedges = 0;
   p->knots  = 0;
   p->terms  = 0;
   p->stp_type = UNKNOWN;
   p->flags  = flags;
   p->layers = layers;
   p->hoplimit = UNKNOWN;
   p->locals = malloc((size_t)layers * sizeof(int));
   p->source = malloc((size_t)layers * sizeof(int));

   p->term   = malloc((size_t)ksize * sizeof(int));
   p->mark   = malloc((size_t)ksize * sizeof(int));
   p->grad   = malloc((size_t)ksize * sizeof(int));
   p->inpbeg = malloc((size_t)ksize * sizeof(int));
   p->outbeg = malloc((size_t)ksize * sizeof(int));

   p->esize = esize;
   p->edges = 0;

   p->cost  = malloc((size_t)esize * sizeof(SCIP_Real));

   if( p->stp_type == STP_ROOTED_PRIZE_COLLECTING || p->stp_type == STP_PRIZE_COLLECTING || p->stp_type == STP_MAX_NODE_WEIGHT )
      p->prize = malloc((size_t)ksize * sizeof(SCIP_Real));
   else
      p->prize  = NULL;

   p->tail  = malloc((size_t)esize * sizeof(int));

   p->head  = malloc((size_t)esize * sizeof(int));

   p->orgtail  = NULL;

   p->orghead  = NULL;

   p->ieat  = malloc((size_t)esize * sizeof(int));
   p->oeat  = malloc((size_t)esize * sizeof(int));

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

   assert(p->locals != NULL);
   assert(p->source != NULL);
   assert(p->term   != NULL);
   assert(p->mark   != NULL);
   assert(p->grad   != NULL);
   assert(p->inpbeg != NULL);
   assert(p->outbeg != NULL);
   assert(p->cost   != NULL);
   assert(p->tail   != NULL);
   assert(p->head   != NULL);
   assert(p->ieat   != NULL);
   assert(p->oeat   != NULL);

   for(i = 0; i < p->layers; i++)
   {
      p->source[i] = -1;
      p->locals[i] =  0;
   }
   return(p);
}

#endif

/** initalize a graph */
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

   SCIP_CALL( SCIPallocMemory(scip, g) );
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
   p->stp_type = UNKNOWN;
   p->flags  = flags;
   p->layers = layers;
   p->hoplimit = UNKNOWN;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->locals), layers) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->source), layers) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->term), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mark), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->grad), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->inpbeg), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->outbeg), ksize) );

   p->esize = esize;
   p->edges = 0;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->cost), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->tail), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->head), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->ieat), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->oeat), esize) );

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
   int* orgtail;             /* (original) tail of all orginal edges  */
   int* orghead;             /* (original) head of all orginal edges  */
   int e;
   int nedges;
   SCIP_Bool pc;

   assert(scip != NULL);
   assert(graph != NULL);

   pc = graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_ROOTED_PRIZE_COLLECTING;

   nedges = graph->edges;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->orgtail), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->orghead), nedges) );

   orgtail = graph->orgtail;
   orghead = graph->orghead;

   for( e = 0; e < nedges; e++ )
   {
      orgtail[e] = graph->tail[e];
      orghead[e] = graph->head[e];
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->ancestors), nedges) );
   ancestors = graph->ancestors;

   for( e = 0; e < nedges; e++ )
   {
      SCIP_CALL( SCIPallocMemory(scip, &(ancestors[e])) ); /*lint !e866*/
      (ancestors)[e]->index = e;
      (ancestors)[e]->parent = NULL;

   }

   if( pc )
   {
      int k;
      int nnodes = graph->knots;
      SCIP_CALL( SCIPallocMemoryArray(scip, &(graph->pcancestors), nnodes) );
      pcancestors = graph->pcancestors;
      for( k = 0; k < nnodes; k++ )
         pcancestors[k] = NULL;
   }

   return SCIP_OKAY;
}

/** enlarge the graph */
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
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->locals), layers) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->source), layers) );
      for( i = g->layers; i < layers; i++ )
      {
         g->source[i] = -1;
         g->locals[i] =  0;
      }
      g->layers = layers;
   }
   if( (ksize > 0) && (ksize != g->ksize) )
   {
      g->ksize  = ksize;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->term), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->mark), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->grad), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->inpbeg), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->outbeg), ksize) );
   }
   if( (esize > 0) && (esize != g->esize) )
   {
      g->esize = esize;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->cost), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->tail), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->head), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->ieat), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->oeat), esize) );
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

            /*printf("curr x1, y1 (%d,%d)  \n ", obst_coords[0][z], obst_coords[1][z]);
              printf("curr x2, y2(%d,%d)  \n ", obst_coords[2][z], obst_coords[3][z]);
              printf("  \n ");*/
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
               /*   printf("edge %d_%d \n ", getNodeNumber(-1), getNodeNumber(j) );*/
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
               /*     printf("edgeXXX %d_%d %d \n ", coords[j][currcoord[j] + 1],  coords[j][currcoord[j]], gridedgecount );*/
               (*gridedgecount)++;
               /*   printf("edge %d_%d \n ", getNodeNumber(-1), getNodeNumber(j) );*/
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

   /* initalize the terminal-coordinates array */
   SCIP_CALL( SCIPallocMemoryArray(scip, &termcoords, grid_dim) );

   for( i = 0; i < grid_dim; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(termcoords[i]), nterms) ); /*lint !e866*/
      for( j = 0; j < nterms; j++ )
         termcoords[i][j] = coords[i][j];
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &ncoords, grid_dim) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &currcoord, grid_dim) );

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
   SCIP_CALL( SCIPallocMemoryArray(scip, &gridedges, 2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &edgecosts, nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(gridedges[0]), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(gridedges[1]), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(inobstacle), nnodes) );
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
   g->stp_type = STP_OBSTACLES_GRID;

   for( i = 0; i < grid_dim; i++ )
      SCIPfreeMemoryArray(scip, &(termcoords[i]));
   SCIPfreeMemoryArray(scip, &(gridedges[0]));
   SCIPfreeMemoryArray(scip, &(gridedges[1]));
   SCIPfreeMemoryArray(scip, &inobstacle);
   SCIPfreeMemoryArray(scip, &termcoords);
   SCIPfreeMemoryArray(scip, &edgecosts);
   SCIPfreeMemoryArray(scip, &gridedges);
   SCIPfreeMemoryArray(scip, &ncoords);
   SCIPfreeMemoryArray(scip, &currcoord);

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

   /* initalize the terminal-coordinates array */
   SCIP_CALL( SCIPallocMemoryArray(scip, &termcoords, grid_dim) );
   for( i = 0; i < grid_dim; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(termcoords[i]), nterms) ); /*lint !e866*/
      for( j = 0; j < nterms; j++ )
         termcoords[i][j] = coords[i][j];
   }
   SCIP_CALL( SCIPallocMemoryArray(scip, &ncoords, grid_dim) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &currcoord, grid_dim) );

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

   SCIP_CALL( SCIPallocMemoryArray(scip, &gridedges, 2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &edgecosts, nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(gridedges[0]), nedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(gridedges[1]), nedges) );

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

#if 0
      if( i == 0 )
      {
         g->source[0] = k;
         printf("root: (%d", termcoords[0][i]);
         for( j = 1; j < grid_dim; j++ )
            printf(", %d", termcoords[j][i]);
         printf(")\n");
      }
#endif
      /* make a terminal out of the node */
      graph_knot_chg(g, k, 0);
   }

   g->stp_type = STP_GRID;

   for( i = 0; i < grid_dim; i++ )
      SCIPfreeMemoryArray(scip, &(termcoords[i]));

   SCIPfreeMemoryArray(scip, &(gridedges[0]));
   SCIPfreeMemoryArray(scip, &(gridedges[1]));
   SCIPfreeMemoryArray(scip, &termcoords);
   SCIPfreeMemoryArray(scip, &edgecosts);
   SCIPfreeMemoryArray(scip, &gridedges);
   SCIPfreeMemoryArray(scip, &ncoords);
   SCIPfreeMemoryArray(scip, &currcoord);
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
   SCIP*                 scip ,              /**< SCIP data structure */
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
      else
      {
         prize[k] = 0;
      }
   }
   graph->source[0] = root;
   assert((nterms + 1) == graph->terms);
   graph->stp_type = STP_PRIZE_COLLECTING;

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
   graph->stp_type = STP_ROOTED_PRIZE_COLLECTING;
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
      if( Is_term(graph->term[i]) )
      {
         assert(SCIPisGE(scip, maxweights[i], 0.0));
         graph->prize[i] = maxweights[i];
         nterms++;
      }
      else
      {
         assert(SCIPisLT(scip, maxweights[i], 0.0));
         graph->prize[i] = 0.0;
      }
   }
   assert(nterms == graph->terms);

   SCIP_CALL( graph_prize_transform(scip, graph) );

   graph->stp_type = STP_MAX_NODE_WEIGHT;
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

   assert(scip != NULL);
   assert(p != NULL);

   SCIPfreeMemoryArray(scip, &(p->locals));
   SCIPfreeMemoryArray(scip, &(p->source));
   SCIPfreeMemoryArray(scip, &(p->term));
   SCIPfreeMemoryArray(scip, &(p->mark));
   SCIPfreeMemoryArray(scip, &(p->grad));
   SCIPfreeMemoryArray(scip, &(p->inpbeg));
   SCIPfreeMemoryArray(scip, &(p->outbeg));
   SCIPfreeMemoryArray(scip, &(p->cost));
   SCIPfreeMemoryArray(scip, &(p->tail));
   SCIPfreeMemoryArray(scip, &(p->head));
   SCIPfreeMemoryArray(scip, &(p->ieat));
   SCIPfreeMemoryArray(scip, &(p->oeat));

   if( p->prize != NULL )
      SCIPfreeMemoryArray(scip, &(p->prize));
   if( p->ancestors != NULL )
   {
      for( e = 0; e < p->edges; e++ )
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
      if( p->orgtail != NULL )
      {
         assert(p->orghead != NULL);
         SCIPfreeMemoryArray(scip, &(p->orgtail));
         SCIPfreeMemoryArray(scip, &(p->orghead));
      }
      curr = p->fixedges;
      while( curr != NULL )
      {
         p->fixedges = curr->parent;
         SCIPfreeMemory(scip, &(curr));
         curr = p->fixedges;
      }

      if( p->pcancestors != NULL )
      {
         for( e = 0; e < p->norgmodelknots; e++ )
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
   }

   if( p->stp_type == STP_DEG_CONS )
   {
      SCIPfreeMemoryArray(scip, &(p->maxdeg));
   }
   else if( p->stp_type == STP_GRID )
   {
      int i;
      for( i = 0; i < p->grid_dim; i++ )
         SCIPfreeMemoryArray(scip, &(p->grid_coordinates[i]));

      SCIPfreeMemoryArray(scip, &(p->grid_coordinates));
      SCIPfreeMemoryArray(scip, &(p->grid_ncoords));
   }
   SCIPfreeMemory(scip, &(p));
}

/** copy the graph */
SCIP_RETCODE graph_copy(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                orgraph,            /**< original graph */
   GRAPH**               copygraph           /**< graph to be copied */
   )
{
   GRAPH* g;
   GRAPH* p;

   p = orgraph;
   assert(p != NULL);

   SCIP_CALL( graph_init(scip, copygraph, p->ksize, p->esize, p->layers, p->flags) );
   g = *copygraph;

   g->norgmodeledges = p->norgmodeledges;
   g->norgmodelknots = p->norgmodelknots;
   g->knots = p->knots;
   g->terms = p->terms;
   g->edges = p->edges;
   g->orgedges = p->orgedges;
   g->orgknots = p->orgknots;
   g->grid_dim = p->grid_dim;
   g->stp_type = p->stp_type;
   g->hoplimit = p->hoplimit;

   BMScopyMemoryArray(g->locals, p->locals, p->layers);
   BMScopyMemoryArray(g->source, p->source, p->layers);
   BMScopyMemoryArray(g->term, p->term, p->ksize);
   BMScopyMemoryArray(g->mark, p->mark, p->ksize);
   BMScopyMemoryArray(g->grad, p->grad, p->ksize);
   BMScopyMemoryArray(g->inpbeg, p->inpbeg, p->ksize);
   BMScopyMemoryArray(g->outbeg, p->outbeg, p->ksize);
   BMScopyMemoryArray(g->term, p->term, p->ksize);
   BMScopyMemoryArray(g->cost, p->cost, p->esize);
   BMScopyMemoryArray(g->tail, p->tail, p->esize);
   BMScopyMemoryArray(g->head, p->head, p->esize);
   BMScopyMemoryArray(g->ieat, p->ieat, p->esize);
   BMScopyMemoryArray(g->oeat, p->oeat, p->esize);
   /*
     memcpy(g->locals, p->locals, p->layers * sizeof(*p->locals));
     memcpy(g->source, p->source, p->layers * sizeof(*p->source));
     memcpy(g->term,   p->term,   p->ksize  * sizeof(*p->term));
     memcpy(g->mark,   p->mark,   p->ksize  * sizeof(*p->mark));
     memcpy(g->grad,   p->grad,   p->ksize  * sizeof(*p->grad));
     memcpy(g->inpbeg, p->inpbeg, p->ksize  * sizeof(*p->inpbeg));
     memcpy(g->outbeg, p->outbeg, p->ksize  * sizeof(*p->outbeg));
     memcpy(g->cost,   p->cost,   p->esize  * sizeof(*p->cost));
     memcpy(g->tail,   p->tail,   p->esize  * sizeof(*p->tail));
     memcpy(g->head,   p->head,   p->esize  * sizeof(*p->head));
     memcpy(g->ieat,   p->ieat,   p->esize  * sizeof(*p->ieat));
     memcpy(g->oeat,   p->oeat,   p->esize  * sizeof(*p->oeat));
   */
   if( g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING
      || g->stp_type == STP_MAX_NODE_WEIGHT )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->prize), g->knots) );
      BMScopyMemoryArray(g->prize, p->prize, g->knots);
   }
   else if( g->stp_type == STP_DEG_CONS )
   {
      assert(p->maxdeg != NULL);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->maxdeg), g->knots) );
      BMScopyMemoryArray(&(g->maxdeg), p->maxdeg, g->knots);
   }
   else if( p->stp_type == STP_GRID )
   {
      int i;
      assert(p->grid_ncoords != NULL);
      assert(p->grid_coordinates != NULL);

      SCIP_CALL( SCIPallocMemoryArray(scip, &(g->grid_coordinates), p->grid_dim) );
      BMScopyMemoryArray(g->grid_coordinates, p->grid_coordinates, p->grid_dim);
      for( i = 0; i < p->grid_dim; i++ )
      {
	 SCIP_CALL( SCIPallocMemoryArray(scip, &(g->grid_coordinates[i]), p->terms) ); /*lint !e866*/
	 BMScopyMemoryArray(g->grid_coordinates[i], p->grid_coordinates[i], p->terms); /*lint !e866*/
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

   if (term != p->term[node])
   {
      if (Is_term(p->term[node]))
      {
         p->terms--;
         p->locals[p->term[node]]--;
      }
      p->term[node] = term;

      if (Is_term(p->term[node]))
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
   int                   t,                  /**< tail node to be contracted */
   int                   s                   /**< head node to be contracted */
   )
{
   typedef struct save_list
   {
      unsigned int mark;
      signed int   edge;
      signed int   knot;
      double       incost;
      double       outcost;
   } SLIST;

   SLIST* slp = NULL;
   IDX**   ancestors = NULL;
   IDX**   revancestors = NULL;
   IDX*   tsancestors = NULL;
   IDX*   stancestors = NULL;
   int    slc = 0;
   int    i;
   int    et;
   int    anti;
   int    es;
   int    cedgeout = UNKNOWN;
   int    head;
   int    tail;
   int    sgrad;

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

   /* change terminal property */
   if( Is_term(p->term[s]) )
   {
      graph_knot_chg(p, t, p->term[s]);
      graph_knot_chg(p, s, -1);
   }

   /* retain root */
   if( p->source[0] == s )
      p->source[0] = t;

   sgrad =  p->grad[s];
   if( sgrad >= 2 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &slp, sgrad - 1) );
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
         assert(slp != NULL);

         ancestors[slc] = NULL;
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors[slc]), p->ancestors[es]) );
         revancestors[slc] = NULL;
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors[slc]), p->ancestors[Edge_anti(es)]) );

         slp[slc].mark = FALSE;
         slp[slc].edge = es;
         slp[slc].knot = p->head[es];
         slp[slc].outcost = p->cost[es];
         slp[slc].incost = p->cost[Edge_anti(es)];
         slc++;

         assert(slc < sgrad);
      }
      else
      {
         cedgeout = Edge_anti(es); /* The edge out of t and into s. */
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &stancestors, p->ancestors[es]) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &tsancestors, p->ancestors[cedgeout]) );
      }
   }

   assert(slc == sgrad - 1);
   assert(tsancestors != NULL);
   assert(stancestors != NULL);
   /* Kantenliste durchgehen
    */
   for( i = 0; i < slc; i++ )
   {
      assert(slp != NULL);

      /* search for an edge out of t with same head as current edge */
      for(et = p->outbeg[t]; et != EAT_LAST; et = p->oeat[et])
         if( p->head[et] == slp[i].knot )
            break;

      /* does such an edge not exist? */
      if( et == EAT_LAST )
      {
         slp[i].mark = TRUE;
      }
      else
      {
         assert(et != EAT_LAST);

         /* This is for nodes with edges to s and t.
          * Need to adjust the out and in costs of the edge
          */
         if( SCIPisGT(scip, p->cost[et], slp[i].outcost) )
         {
            SCIPintListNodeFree(scip, &((p->ancestors)[et]));
	    assert(ancestors != NULL);
	    assert(slp != NULL);
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[et]), ancestors[i]) );
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[et]), tsancestors) );
            p->cost[et] = slp[i].outcost;
         }
         if( SCIPisGT(scip, p->cost[Edge_anti(et)], slp[i].incost) )
         {
            anti = Edge_anti(et);
            SCIPintListNodeFree(scip, &(p->ancestors[anti]));
	    assert(revancestors != NULL);
	    assert(slp != NULL);
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[anti]), revancestors[i]) );
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &((p->ancestors)[anti]), stancestors) );
            p->cost[anti] = slp[i].incost;
         }
      }
   }

   /* insert edges */
   for( i = 0; i < slc; i++ )
   {
      assert(slp != NULL);
      if( slp[i].mark )
      {
         es = p->outbeg[s];

         assert(es != EAT_LAST);
	 assert(ancestors != NULL);
         assert(revancestors != NULL);
         assert(ancestors[i] != NULL);
         assert(revancestors[i] != NULL);
	 assert(slp != NULL);
         SCIPintListNodeFree(scip, &(p->ancestors[es]));
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), ancestors[i]) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), tsancestors) );
         graph_edge_del(scip, p, es, FALSE);

         head = slp[i].knot;
         tail = t;

         p->grad[head]++;
         p->grad[tail]++;

         p->cost[es]     = slp[i].outcost;
         p->tail[es]     = tail;
         p->head[es]     = head;
         p->ieat[es]     = p->inpbeg[head];
         p->oeat[es]     = p->outbeg[tail];
         p->inpbeg[head] = es;
         p->outbeg[tail] = es;

         es = Edge_anti(es);
         SCIPintListNodeFree(scip, &(p->ancestors[es]));

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), revancestors[i]) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(p->ancestors[es]), stancestors) );
         p->cost[es]     = slp[i].incost;
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
      SCIPfreeBufferArray(scip, &slp);
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
   g->prize[i] -= cost;
   for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      if( Is_pterm(g->term[g->head[e]]) )
         break;

   assert(e != EAT_LAST);
   assert(!g->mark[g->head[e]]);

   j = g->head[e];

   assert(j != g->source[0]);

   for( e = g->inpbeg[j]; e != EAT_LAST; e = g->ieat[e] )
      if( g->source[0] == g->tail[e] )
         break;

   assert(e != EAT_LAST);
   assert(!g->mark[g->tail[e]] || g->stp_type == STP_ROOTED_PRIZE_COLLECTING);

   g->cost[e] -= cost;

   assert(SCIPisGE(scip, g->prize[i], 0.0));
   assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));
}

/** contract an edge in (rooted) price collecting */
SCIP_RETCODE graph_knot_contractpc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   t,                  /**< tail node to be contracted */
   int                   s,                  /**< head node to be contracted */
   int                   i                   /**< terminal to add offset to */
   )
{
   int ets;
   assert(g != NULL);
   assert(scip != NULL);
   assert(Is_term(g->term[i]));

   for( ets = g->outbeg[t]; ets != EAT_LAST; ets = g->oeat[ets] )
      if( g->head[ets] == s )
         break;
   assert(ets != EAT_LAST);

   if( Is_term(g->term[t]) && Is_term(g->term[s]) )
   {
      IDX*   etsancestors = NULL;
      int e;
      int j;

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &etsancestors, g->ancestors[ets]) );

      for( e = g->outbeg[s]; e != EAT_LAST; e = g->oeat[e] )
         if( Is_pterm(g->term[g->head[e]]) )
            break;
      assert(e != EAT_LAST);
      /*
        assert(g->pcancestors != NULL);
        if( g->pcancestors[s] != NULL )
        {
        SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[t]), g->pcancestors[s]) );
        SCIPintListNodeFree(scip, &(g->pcancestors[s]));
        }
        SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[t]), g->ancestors[e]) );
      */
      j = g->head[e];

      assert(j != g->source[0]);
      assert(!g->mark[j]);

      graph_knot_chg(g, j, -1);
      graph_edge_del(scip, g, e, TRUE);

      e = g->inpbeg[j];
      assert(e != EAT_LAST);
      assert(g->source[0] == g->tail[e]);

      assert(SCIPisEQ(scip, g->prize[s], g->cost[e]));

      prize_subtract(scip, g, g->cost[ets] - g->prize[s], i);
      graph_edge_del(scip, g, e, TRUE);

      SCIP_CALL( graph_knot_contract(scip, g, t, s) );

      SCIPdebugMessage("PC contract: %d, %d \n", t, s);

      for( e = g->inpbeg[t]; e != EAT_LAST; e = g->ieat[e] )
	 SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[e]), etsancestors) );
      for( e = g->outbeg[t]; e != EAT_LAST; e = g->oeat[e] )
	 SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[e]), etsancestors) );

      SCIPintListNodeFree(scip, &etsancestors);
   }
   else
   {
      prize_subtract(scip, g, g->cost[ets], i);
      SCIP_CALL( graph_knot_contract(scip, g, t, s) );
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

      SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), revancestors0) );
      SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors1) );

      SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), ancestors0) );
      SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors1) );
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
   GRAPH* p,
   int    e)
{
   int    i;
   int    head;
   int    tail;

   assert(p          != NULL);
   assert(e          >= 0);
   assert(e          <  p->edges);

   head = p->head[e];
   tail = p->tail[e];

   if (p->inpbeg[head] == e)
      p->inpbeg[head] = p->ieat[e];
   else
   {
      for(i = p->inpbeg[head]; p->ieat[i] != e; i = p->ieat[i])
         assert(i >= 0);

      p->ieat[i] = p->ieat[e];
   }
   if (p->outbeg[tail] == e)
      p->outbeg[tail] = p->oeat[e];
   else
   {
      for(i = p->outbeg[tail]; p->oeat[i] != e; i = p->oeat[i])
         assert(i >= 0);

      p->oeat[i] = p->oeat[e];
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

   if( graph->stp_type == STP_ROOTED_PRIZE_COLLECTING )
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
#if 0
GRAPH *graph_pack2(
   SCIP*  scip,
   GRAPH* p,
   SCIP_Bool verbose)
{
   const char* msg1   = "Knots: %d  Edges: %d  Terminals: %d\n";
   SCIP_RETCODE rcode;
   GRAPH* q;
   int*   new;
   int    knots = 0;
   int    edges = 0;
   int    i;
   int    l;

   assert(p      != NULL);
   assert(graph_valid(p));
   if( verbose )
      printf("Packing graph: ");
   printf("Packing graph: \n ");
   new = malloc((size_t)p->knots * sizeof(new[0]));

   assert(new != NULL);

   /* Knoten zaehlen
    */
   for(i = 0; i < p->knots; i++)
   {
      new[i] = knots;

      /* Hat der Knoten noch Kanten ?
       */
      if (p->grad[i] > 0)
         knots++;
      else
         new[i] = -1;
   }

   /* Ist ueberhaupt noch ein Graph vorhanden ?
    */
   if (knots == 0)
   {
      free(new);
      new = NULL;
      if( verbose )
         printf(" graph vanished!\n");

      knots = 1;
   }

   /* Kanten zaehlen
    */
   for(i = 0; i < p->edges; i++)
   {
      if (p->oeat[i] != EAT_FREE)
      {
         assert(p->ieat[i] != EAT_FREE);
         edges++;
      }
   }
   if( knots == 1 )
      assert(edges == 0);
   rcode = graph_init(scip, &q, knots, edges, p->layers, p->flags);
   q->norgmodelknots = p->norgmodelknots;
   q->norgmodeledges = p->norgmodeledges;
   q->orgtail = p->orgtail;
   q->orghead = p->orghead;
   q->orgknots = p->knots;
   q->orgedges = p->edges;
   q->stp_type = p->stp_type;
   q->maxdeg = p->maxdeg;
   q->grid_dim = p->grid_dim;
   q->grid_ncoords = p->grid_ncoords;
   q->grid_coordinates = p->grid_coordinates;
   q->fixedges = p->fixedges;
   /*SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(q->fixedges), p->fixedges);*/

   q->hoplimit = p->hoplimit;
   if( new == NULL )
   {
      q->ancestors = NULL;
      graph_free(scip, p, FALSE);
      graph_knot_add(q, 0);
      q->source[0] = 0;
      return q;
   }

   q->ancestors = malloc((size_t)(edges) * sizeof(IDX*));
   for( i = 0; i < edges; i++ )
      q->ancestors[i] = NULL;

   /* Knoten umladen
    */
   for(i = 0; i < p->knots; i++)
   {
      assert(p->term[i] < p->layers);
#if 0
      if ((i % 100) == 0)
      {
         (void)fputc('k', stdout);
         (void)fflush(stdout);
      }
#endif
      if( p->grad[i] > 0 )
         graph_knot_add(q, p->term[i]);
   }

   /* Kanten umladen
    */
   for( i = 0; i < p->edges; i += 2 )
   {
#if 0
      if ((i % 1000) == 0)
      {
         (void)fputc('e', stdout);
         (void)fflush(stdout);
      }
#endif
      if (p->ieat[i] == EAT_FREE)
      {
         assert(p->oeat[i]     == EAT_FREE);
         assert(p->ieat[i + 1] == EAT_FREE);
         assert(p->oeat[i + 1] == EAT_FREE);
         SCIPintListNodeFree(scip, &(p->ancestors[i]));
         SCIPintListNodeFree(scip, &(p->ancestors[i + 1]));
         continue;
      }
      assert(p->ieat[i]      != EAT_FREE);
      assert(p->oeat[i]      != EAT_FREE);
      assert(p->ieat[i + 1]  != EAT_FREE);
      assert(p->oeat[i + 1]  != EAT_FREE);
      assert(new[p->tail[i]] >= 0);
      assert(new[p->head[i]] >= 0);

      rcode = SCIPintListNodeAppendCopy(scip, &(q->ancestors[q->edges]), p->ancestors[i]);
      rcode = SCIPintListNodeAppendCopy(scip, &(q->ancestors[q->edges + 1]), p->ancestors[i + 1]);
      graph_edge_add(scip, q, new[p->tail[i]], new[p->head[i]],
         p->cost[i], p->cost[Edge_anti(i)]);


   }

   /* Wurzeln umladen
    */
   for(l = 0; l < q->layers; l++)
   {
      assert(q->term[new[p->source[l]]] == l);
      q->source[l] = new[p->source[l]];
   }

   free(new);

   p->stp_type = UNKNOWN;
   graph_free(scip, p, FALSE);

#if 0
   for(l = 0; l < q->layers; l++)
      q->source[l] = -1;

   for(i = 0; i < q->knots; i++)
      if ((q->term[i] >= 0) && ((q->source[q->term[i]] < 0)
            || (q->grad[i] > q->grad[q->source[q->term[i]]])))
         q->source[q->term[i]] = i;
#endif
   assert(q->source[0] >= 0);
   if( verbose )
      printf(msg1, q->knots, q->edges, q->terms);

   return(q);
}
#endif

/** pack the graph, i.e. build a new graph that discards deleted edges and nodes */
SCIP_RETCODE graph_pack(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   GRAPH**               newgraph,           /**< the new graph */
   SCIP_Bool             verbose             /**< verbose? */
   )
{
   const char* msg1   = "Nodes: %d  Edges: %d  Terminals: %d\n";
   GRAPH* g;
   GRAPH* q;
   int*   new;
   int    i;
   int    nnodes;
   int    nedges;

   assert(scip      != NULL);
   assert(graph      != NULL);
   assert(graph_valid(graph));

   g = graph;
   nnodes = 0;
   nedges = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &new, g->knots) );

   if( verbose )
      printf("Packing graph: ");

   /* count nodes */
   for( i = 0; i < g->knots; i++ )
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
   for( i = 0; i < g->edges; i++ )
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
      graph_knot_add(q, 0);
      q->source[0] = 0;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(q->ancestors), nedges) );
   for( i = 0; i < nedges; i++ )
      q->ancestors[i] = NULL;
#if 0
   if( pc )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(q->pcancestors), nedges) );
      for( i = 0; i < nnodes; i++ )
         q->pcancestors[i] = NULL;
   }
#endif

   /* add nodes (of positive degree) */
   for( i = 0; i < g->knots; i++ )
   {
      assert(g->term[i] < g->layers);
      if( g->grad[i] > 0 )
      {
#if 0
         if( pc )
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(q->pcancestors[q->knots]), g->pcancestors[i]) );
#endif
         graph_knot_add(q, g->term[i]);
      }
   }

   /* add edges */
   for( i = 0; i < g->edges; i += 2 )
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
      assert(g->ieat[i]      != EAT_FREE);
      assert(g->oeat[i]      != EAT_FREE);
      assert(g->ieat[i + 1]  != EAT_FREE);
      assert(g->oeat[i + 1]  != EAT_FREE);
      assert(new[g->tail[i]] >= 0);
      assert(new[g->head[i]] >= 0);

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(q->ancestors[q->edges]), g->ancestors[i]) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(q->ancestors[q->edges + 1]), g->ancestors[i + 1]) );

      graph_edge_add(scip, q, new[g->tail[i]], new[g->head[i]], g->cost[i], g->cost[Edge_anti(i)]);
   }

   /* add root */
   assert(q->term[new[g->source[0]]] == 0);
   q->source[0] = new[g->source[0]];

   SCIPfreeBufferArray(scip, &new);

   g->stp_type = UNKNOWN;
   graph_free(scip, g, FALSE);

   assert(q->source[0] >= 0);
   if( verbose )
      printf(msg1, q->knots, q->edges, q->terms);

   return SCIP_OKAY;
}

/** traverse the graph */
void graph_trail(
   const GRAPH*          p,                  /**< the new graph */
   int                   i                   /**< node to start from */
   )
{
   int   k;

   assert(p      != NULL);
   assert(i      >= 0);
   assert(i      <  p->knots);

   if( !p->mark[i] )
   {
      p->mark[i] = TRUE;

      for( k = p->outbeg[i]; k != EAT_LAST; k = p->oeat[k] )
         if (!p->mark[p->head[k]])
            graph_trail(p, p->head[k]);
   }
}

/** is the graph valid? */
int graph_valid(
   const GRAPH*          p                   /**< the new graph */
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
   int    terms;
#if 0
   int*   locals;
   locals = malloc((size_t)p->layers * sizeof(int));
   assert(locals != NULL);
   for( l = 0; l < p->layers; l++ )
      locals[l] = p->locals[l];
#endif
   assert(p      != NULL);

   terms  = p->terms;

   for( k = 0; k < p->knots; k++ )
   {
      if( Is_term(p->term[k]) )
      {
#if 0
         locals[p->term[k]]--;
#endif
         terms--;
      }
      for( e = p->inpbeg[k]; e != EAT_LAST; e = p->ieat[e] )
         if( p->head[e] != k )
            break;

      if( e != EAT_LAST )
         return((void)fprintf(stderr, fehler1, k, e, p->tail[e], p->head[e]), FALSE);

      for( e = p->outbeg[k]; e != EAT_LAST; e = p->oeat[e] )
         if( p->tail[e] != k )
            break;

      if( e != EAT_LAST )
         return((void)fprintf(stderr, fehler2, k, e, p->tail[e], p->head[e]), FALSE);
   }
   if( terms != 0 )
      return((void)fprintf(stderr, fehler9, p->terms, p->terms - terms), FALSE);

   for( l = 0; l < p->layers; l++ )
   {
#if 0
      if( locals[l] != 0 )
         return((void)fprintf(stderr, fehler8,
               l, p->locals[l], p->locals[l] - locals[l]), FALSE);
#endif
      if( (p->source[l] < 0 )
         || (p->source[l] >= p->knots)
         || (p->term[p->source[l]] != l))
         return((void)fprintf(stderr, fehler3,
               l, p->source[l], p->term[p->source[l]]), FALSE);
   }
#if 0
   free(locals);
#endif
   for( e = 0; e < p->edges; e += 2 )
   {
      if( (p->ieat[e    ] == EAT_FREE) && (p->oeat[e    ] == EAT_FREE)
         && (p->ieat[e + 1] == EAT_FREE) && (p->oeat[e + 1] == EAT_FREE) )
         continue;

      if( (p->ieat[e] == EAT_FREE) || (p->oeat[e] == EAT_FREE)
         || (p->ieat[e + 1] == EAT_FREE) || (p->oeat[e + 1] == EAT_FREE) )
         return((void)fprintf(stderr, fehler4, e, e + 1), FALSE);

      if( (p->head[e] != p->tail[e + 1]) || (p->tail[e] != p->head[e + 1]) )
         return((void)fprintf(stderr, fehler5,
               e, e + 1, p->head[e], p->tail[e + 1],
               p->tail[e], p->head[e + 1]), FALSE);

   }
   for( k = 0; k < p->knots; k++ )
      p->mark[k] = FALSE;

   graph_trail(p, p->source[0]);

   for( k = 0; k < p->knots; k++ )
   {
      if( (p->grad[k] == 0)
         && ((p->inpbeg[k] != EAT_LAST) || (p->outbeg[k] != EAT_LAST)) )
         return((void)fprintf(stderr, fehler6, k), FALSE);

      if( !p->mark[k] && (p->grad[k] > 0) && p->stp_type != STP_PRIZE_COLLECTING && p->stp_type != STP_MAX_NODE_WEIGHT )
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

   char* terminal;
   int* pnode;
   int e;
   int i;
   int root;
   int nnodes;
   int termcount;
   assert(graph != NULL);
   assert(result != NULL);
   nnodes = graph->knots;
   root = graph->source[0];
   assert(root >= 0);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &terminal, nnodes) );
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

   if( termcount != graph->terms )
   {
      for( i = 0; i < nnodes; i++ )
         if( Is_term(graph->term[i]) && !terminal[i] )
            printf("not reached, node: %d\n", i);
      printf("a: %d, b: %d: \n", termcount, graph->terms);
   }

   SCIPfreeBufferArray(scip, &terminal);
   SCIPqueueFree(&queue);

   return (termcount == graph->terms);
}


char graph_valid2(
   SCIP* scip,
   const GRAPH* graph,
   SCIP_Real* cost
   )
{

   SCIP_QUEUE* queue;

   char* terminal;
   char* reached;
   int* pnode;
   int e;
   int i;
   int root;
   int nnodes;
   int termcount;
   assert(graph != NULL);
   assert(cost != NULL);
   nnodes = graph->knots;
   root = graph->source[0];
   assert(root >= 0);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &terminal, nnodes) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &reached, nnodes) );
   for( i = 0; i < nnodes; i++ )
   {
      terminal[i] = FALSE;
      reached[i] = FALSE;
   }
   /* BFS until all terminals are reached */
   SCIP_CALL_ABORT( SCIPqueueCreate(&queue, nnodes, 2.0) );

   SCIP_CALL_ABORT( SCIPqueueInsert(queue, &root) );
   termcount = 1;
   terminal[root] = TRUE;
   reached[root] = TRUE;
   while( !SCIPqueueIsEmpty(queue) )
   {
      pnode = (SCIPqueueRemove(queue));
      for( e = graph->outbeg[*pnode]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( SCIPisLT(scip, cost[e], BLOCKED) && !reached[graph->head[e]] )
         {
            i = graph->head[e];
            reached[i] = TRUE;
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
   SCIPfreeBufferArray(scip, &reached);
   SCIPfreeBufferArray(scip, &terminal);
   SCIPqueueFree(&queue);
   if (termcount != graph->terms)
   {
      for( i = 0; i < nnodes; i++ )
         if( Is_term(graph->term[i]) && !terminal[i] )
            printf("not reached, node: %d\n", i);
      return FALSE;
   }

   return (termcount == graph->terms);
}
