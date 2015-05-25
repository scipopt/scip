/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Functions                                                     */
/*   File....: grphbase.c                                                    */
/*   Name....: Basic Graph Routines                                          */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*lint -esym(750,GRPHBASE_C) -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h)                                                   */

/* Ein Graph wird initialisiert, und dannach veraendert sich die Anzahl
 * seiner Knoten 'g->knots' und Kanten 'g->edges' nicht mehr kleiner.
 * Allerdings kann der Grad eines Knotens auf 0 zurueckgehen und eine
 * Kante als EAT_FREE gekennzeichnet werden. Wird dann 'graph_pack()'
 * aufgerufen, werden diese Knoten und Kanten nicht uebernommen.
 */
#define GRPHBASE_C
#include "scip/misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "portab.h"
#include "misc_stp.h"
#include "grph.h"
#if 0
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

   //(*g) = malloc(sizeof(GRAPH));
   SCIP_CALL( SCIPallocMemory(scip, g) );
   p = *g;
   assert(p != NULL);

   /* ancestor data for retransformation after reductions */
   p->fixedges = NULL;
   p->ancestors = NULL;
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
   /*
     p->locals = malloc((size_t)layers * sizeof(int));
     p->source = malloc((size_t)layers * sizeof(int));

     p->term   = malloc((size_t)ksize * sizeof(int));
     p->mark   = malloc((size_t)ksize * sizeof(int));
     p->grad   = malloc((size_t)ksize * sizeof(int));
     p->inpbeg = malloc((size_t)ksize * sizeof(int));
     p->outbeg = malloc((size_t)ksize * sizeof(int));
   */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->locals), layers) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->source), layers) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->term), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->mark), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->grad), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->inpbeg), ksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->outbeg), ksize) );

   if( p->stp_type == STP_ROOTED_PRIZE_COLLECTING || p->stp_type == STP_PRIZE_COLLECTING || p->stp_type == STP_MAX_NODE_WEIGHT )
      //p->prize = malloc((size_t)ksize * sizeof(SCIP_Real));
      SCIP_CALL( SCIPallocMemoryArray(scip, &(p->prize), ksize) );
   else
      p->prize  = NULL;

   p->esize = esize;
   p->edges = 0;

   //p->cost  = malloc((size_t)esize * sizeof(SCIP_Real));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->cost), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->tail), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->head), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->ieat), esize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(p->oeat), esize) );
   /*
     p->tail  = malloc((size_t)esize * sizeof(int));

     p->head  = malloc((size_t)esize * sizeof(int));
   */

   /*
     p->ieat  = malloc((size_t)esize * sizeof(int));
     p->oeat  = malloc((size_t)esize * sizeof(int));
   */
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

   for( i = 0; i < p->layers; i++ )
   {
      p->source[i] = -1;
      p->locals[i] =  0;
   }
   return SCIP_OKAY;
}

/* initialize data structures required to keep track of reductions */
SCIP_RETCODE graph_init_history(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph */
   int**                 orgtail,            /**< (original) tail of all orginal edges */
   int**                 orghead,            /**< (original) head of all orginal edges  */
   IDX***                ancestors           /**< pointer to ancestor list array (on all edges) */
   )
{
   int e;
   int nedges;
   assert(graph != NULL);

   nedges = graph->edges;

   (*orgtail) = malloc((size_t)nedges * sizeof(int));
   (*orghead) = malloc((size_t)nedges * sizeof(int));

   for( e = 0; e < nedges; e++ )
   {
      (*orgtail)[e] = graph->tail[e];
      (*orghead)[e] = graph->head[e];
   }

   *ancestors = malloc((size_t)(nedges) * sizeof(IDX*));
   for( e = 0; e < nedges; e++ )
   {
      //SCIP_CALL( SCIPallocMemory(scip, &(*ancestors)[e]) );
      (*ancestors)[e] = malloc((size_t)sizeof(IDX));
      (*ancestors)[e]->index = e;
      (*ancestors)[e]->parent = NULL;

   }
   return SCIP_OKAY;
}

/* enlarge graph */
SCIP_RETCODE graph_resize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph to be resized */
   int                   ksize,              /**< new node slots */
   int                   esize,              /**< new edge slots */
   int                   layers              /**< layers (set to 1 by default) */
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
      //g->locals = realloc(g->locals, (size_t)layers * sizeof(int));
      //g->source = realloc(g->source, (size_t)layers * sizeof(int));
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
      /*
        g->term   = realloc(g->term,   (size_t)ksize * sizeof(int));
        g->mark   = realloc(g->mark,   (size_t)ksize * sizeof(int));
        g->grad   = realloc(g->grad,   (size_t)ksize * sizeof(int));
        g->inpbeg = realloc(g->inpbeg, (size_t)ksize * sizeof(int));
        g->outbeg = realloc(g->outbeg, (size_t)ksize * sizeof(int));
      */
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->term), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->mark), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->grad), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->inpbeg), ksize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->outbeg), ksize) );
   }
   if( (esize > 0) && (esize != g->esize) )
   {
      g->esize = esize;
      /*
        g->cost  = realloc(g->cost, (size_t)esize * sizeof(SCIP_Real));
        g->tail  = realloc(g->tail, (size_t)esize * sizeof(int));
        g->head  = realloc(g->head, (size_t)esize * sizeof(int));
        g->ieat  = realloc(g->ieat, (size_t)esize * sizeof(int));
        g->oeat  = realloc(g->oeat, (size_t)esize * sizeof(int));
      */
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->cost), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->tail), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->head), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->ieat), esize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(g->oeat), esize) );
   }
   if( g->stp_type == STP_GRID )
   {
      g->grid_ncoords = realloc(g->grid_ncoords, (size_t)g->grid_dim * sizeof(int));
      assert(g->grid_ncoords != NULL);
   }
   /*
     assert(g->locals != NULL);
     assert(g->source != NULL);
     assert(g->term   != NULL);
     assert(g->mark   != NULL);
     assert(g->grad   != NULL);
     assert(g->inpbeg != NULL);
     assert(g->outbeg != NULL);
     assert(g->cost   != NULL);
     assert(g->tail   != NULL);
     assert(g->head   != NULL);
     assert(g->ieat   != NULL);
     assert(g->oeat   != NULL);
   */
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
   while( i < ncoords[coord] ) //for( i = 0; i < ncoords[coord]; i++ )
   {
      currcoord[coord] = i;
      if( coord < grid_dim - 1 )
         compEdgesObst(coord + 1, grid_dim, nobstacles, ncoords, currcoord, edgecosts, gridedgecount, coords, gridedges, obst_coords, inobstacle);
      else
      {
         x = coords[0][currcoord[0]];
         y = coords[1][currcoord[1]];
         // printf("curr cord (%d,%d)  \n ", x, y);
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
               //printf("obst node found %d, %d inobst:%d \n", x, y, node);
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
                  //  printf("edgeXXX %d_%d %d \n ", coords[j][currcoord[j] + 1],  coords[j][currcoord[j]], *gridedgecount );
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
   int    nnodes = 0;
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
      SCIP_CALL( SCIPallocMemoryArray(scip, &(termcoords[i]), nterms) );
      for( j = 0; j < nterms; j++ )
         termcoords[i][j] = coords[i][j];
   }

   SCIPallocMemoryArray(scip, &ncoords, grid_dim);
   SCIPallocMemoryArray(scip, &currcoord, grid_dim);

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
   g->grid_ncoords = (int*) malloc(grid_dim * sizeof(int));
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
         graph_edge_add(g, gridedges[0][i] - 1, gridedges[1][i] - 1, cost, cost);
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
      tmp = -1;

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
   int    nnodes = 0;
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
      SCIP_CALL( SCIPallocMemoryArray(scip, &(termcoords[i]), nterms) );
      for( j = 0; j < nterms; j++ )
         termcoords[i][j] = coords[i][j];
   }
   SCIPallocMemoryArray(scip, &ncoords, grid_dim);
   SCIPallocMemoryArray(scip, &currcoord, grid_dim);

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

   g->grid_ncoords = (int*) malloc(grid_dim * sizeof(int));
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
      graph_edge_add(g, gridedges[0][i] - 1, gridedges[1][i] - 1, cost, cost);
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
      tmp = -1;
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
void graph_grid_coordinates(
   int**                 coords,             /**< coordinates */
   int**                 nodecoords,         /**<  coordinates of the node (to be computed) */
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
      *nodecoords = (int*) malloc(grid_dim * sizeof(int));

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
   double tmpsum = 0.0;
   assert(graph != NULL);
   assert(graph->edges == graph->esize);
   root = graph->source[0];
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
         assert(GT(prize[k], 0));
         tmpsum += prize[k];

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(graph, root, k, 0, FARAWAY);
         graph_edge_add(graph, root, node, prize[k], FARAWAY);
         graph_edge_add(graph, k, node, 0, FARAWAY);
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
   SCIP_Real tmpsum = 0.0;

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
         assert(GT(prize[k], 0));
         tmpsum += prize[k];

         /* add one edge going from the root to the 'copied' terminal and one going from the former terminal to its copy */
         graph_edge_add(graph, root, node, prize[k], FARAWAY);
         graph_edge_add(graph, k, node, 0, FARAWAY);
      }
      else
      {
         prize[k] = 0;
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
   SCIP_Real*            maxweights
   )
{
   //double* prize;
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
            //printf("edgecost2:  %d-%d  %.12f  \n", graph->tail[e], graph->head[e], graph->cost[e]);
         }
      }
      else
      {
         //printf("term:  %d \n", i);
         graph_knot_chg(graph, i, 0);
         nterms++;
      }
   }
   nterms = 0;
   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(graph->term[i]) )
      {
         assert(!SCIPisLT(scip, maxweights[i], 0.0));
         graph->prize[i] = maxweights[i];
         //printf("termcost/prize2:  %.12f \n", maxweights[i]);
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


#if 0

/** transform the graph to enabls solving of group Steiner problems */
double graph_group_convert(
   GRAPH* g)
{
   double  max_cost = 0.0;
   double  fixed    = 0.0;
   int     real_knots;
   int     pseudo;
   int     i;

   assert(g != NULL);

   if (g->groups < 2)
      return 0.0;

   assert(graph_valid(g));

   graph_resize(g, g->knots + g->groups, g->edges + 2 * g->terms,
      NO_CHANGE, NO_CHANGE);

   assert(g != NULL);

   real_knots = g->knots;

   for(i = 0; i < g->edges; i++)
      if (g->cost[i] > max_cost)
         max_cost = g->cost[i];

   /* So hoch muss man wohl sein, um bei den Reduktionen keinen Probleme
    * zu bekommen.
    */
   max_cost += 1.0;

   for(i = 0; i < g->groups; i++)
      graph_knot_add(g, 0, 0, 0, 0);

   for(i = 0; i < real_knots; i++)
   {
      if (!Is_term(g->term[i]))
         continue;

      assert(g->group[i] >= 0);

      pseudo = real_knots + g->group[i];

      graph_edge_add(g, i, pseudo, max_cost, max_cost);
      graph_knot_chg(g, i, -1, -1, NO_CHANGE, NO_CHANGE);

      /* Das ist noch nicht richtig toll.
       */
      graph_knot_chg(g, pseudo, NO_CHANGE, NO_CHANGE,
         g->xpos[i] + 1, g->ypos[i] + 1);
   }
   g->source[0] = real_knots;

   fixed = max_cost * g->groups;

   printf("Group converter: Added %d knots, cost: %g\n",
      g->groups, fixed);

   return -fixed;
}



void graph_group_compsdcost(
   const GRAPH* g,
   double*      cost)
{
   PATH*  path;
   double max_dist;
   int    i;
   int    k;
   int    j;

   assert(g    != NULL);
   assert(cost != NULL);

   if (g->groups < 2)
      return;

   path = malloc(g->knots * sizeof(*path));

   assert(path != NULL);

   for(i = 0; i < g->knots; i++)
      g->mark[i] = TRUE;

   for(k = 0; k < g->edges; k++)
      cost[k] = g->cost[k];

   for(i = 0; i < g->knots; i++)
   {
      /* Wir gehen die Gruppen anhand der Pseudoterminal durch.
       */
      if (!Is_term(g->term[i]))
         continue;

      for(k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k])
      {
         graph_path_exec(g, FSP_MODE, g->head[k], g->cost, path);

         max_dist = 0.0;

         for(j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j])
            if (path[g->head[j]].dist > max_dist)
               max_dist = path[g->head[j]].dist;

         printf("Term %d max_dist: %g\n", g->head[k], max_dist);

         /* Entweder ist sonst niemand in der Gruppe, dann ist
          * die Distanz 0 oder es sind noch welche da, dann ist da
          * auch eine positive Distanz.
          */
         assert(((g->grad[i] == 1) && EQ(max_dist, 0.0)) || GT(max_dist, 0.0));
         assert(g->head[k] == g->tail[Edge_anti(k)]);

         cost[k]            = max_dist + 1.0;
         cost[Edge_anti(k)] = max_dist + 1.0;
      }
   }
   free(path);
}


#endif
void graph_free(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                p,                  /**< graph to be freed */
   SCIP_Bool             final               /**< delete ancestor data structures? */
   )
{
   IDX* curr;
   int e;
   assert(p != NULL);
   /*
     free(p->locals);
     free(p->source);
     free(p->term);
     free(p->mark);
     free(p->grad);
     free(p->inpbeg);
     free(p->outbeg);
     free(p->cost);
   */
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
      //free(p->prize);
      SCIPfreeMemoryArray(scip, &(p->prize));
   if( p->ancestors != NULL )
   {
      for( e = 0; e < p->edges; e++ )
      {
         curr = p->ancestors[e];
         while( curr != NULL )
         {
            p->ancestors[e] = curr->parent;
            free(curr);
            curr = p->ancestors[e];
         }
      }
      free(p->ancestors);
   }

   if( final )
   {
      if( p->orgtail != NULL )
      {
         assert(p->orghead != NULL);
         free(p->orgtail);
         free(p->orghead);
      }
      curr = p->fixedges;
      while( curr != NULL )
      {
         p->fixedges = curr->parent;
         free(curr);
         curr = p->fixedges;
      }
   }

   if( p->stp_type == STP_DEG_CONS )
      free(p->maxdeg);
   else if(p->stp_type == STP_GRID )
   {
      int i;
      for( i = 0; i < p->grid_dim; i++ )
         free(p->grid_coordinates[i]);
      free(p->grid_coordinates);
      free(p->grid_ncoords);
   }
   //free(p);
   SCIPfreeMemory(scip, &(p));
}

SCIP_RETCODE graph_copy(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                orgraph,            /**< original graph */
   GRAPH**               copygraph           /**< graph to be copied */
   )
{
   GRAPH* g;
   GRAPH* p;

   assert(p != NULL);

   p = orgraph;

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
      //memcpy(g->prize,   p->prize,   p->ksize  * sizeof(*p->prize));
      BMScopyMemoryArray(g->prize, p->prize, p->ksize);
   }
   else if( g->stp_type == STP_DEG_CONS )
   {
      assert(p->maxdeg != NULL);
      g->maxdeg = malloc((size_t)(g->knots) * sizeof(int));
      memcpy(g->maxdeg,   p->maxdeg,   p->knots  * sizeof(*p->maxdeg));
   }
   else if( p->stp_type == STP_GRID )
   {
      int i;
      assert(p->grid_ncoords != NULL);
      assert(p->grid_coordinates != NULL);

      g->grid_coordinates = malloc((size_t)(p->grid_dim) * sizeof(int*));
      memcpy(g->grid_coordinates,   p->grid_coordinates,   p->grid_dim * sizeof(*p->grid_coordinates));
      for( i = 0; i < p->grid_dim; i++)
      {
         g->grid_coordinates[i] = malloc((size_t)(p->terms) * sizeof(int));
         memcpy(g->grid_coordinates[i],   p->grid_coordinates[i],   p->terms * sizeof(*(p->grid_coordinates[i])));
      }
      g->grid_ncoords = malloc((size_t)(p->grid_dim) * sizeof(int));
      memcpy(g->grid_ncoords,   p->grid_ncoords,   p->grid_dim * sizeof(*p->grid_ncoords));
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

/* ARGSUSED */
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

/* ARGSUSED */
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


SCIP_RETCODE graph_knot_contract(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                p,                  /**< the graph */
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
   if (p->source[0] == s)
      p->source[0] = t;

   sgrad =  p->grad[s];
   if( sgrad >= 2 )
   {
      slp = malloc((size_t)(sgrad - 1) * sizeof(SLIST));
      ancestors = malloc((size_t)(sgrad - 1) * sizeof(IDX*));
      revancestors = malloc((size_t)(sgrad - 1) * sizeof(IDX*));
      assert(slp != NULL);
   }

   /* Liste mit Kanten des aufzuloesenden Knotens merken
    */
   for( es = p->outbeg[s]; es != EAT_LAST; es = p->oeat[es] )
   {
      assert(p->tail[es] == s);

      if( p->head[es] != t )
      {
         i = 0;
         ancestors[slc] = NULL;
         SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(ancestors[slc]), p->ancestors[es]) );
         revancestors[slc] = NULL;
         SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(revancestors[slc]), p->ancestors[Edge_anti(es)]) );

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
         SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &stancestors, p->ancestors[es]) );
         SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &tsancestors, p->ancestors[cedgeout]) );
      }
   }

   assert(slc == sgrad - 1);
   assert(tsancestors != NULL);
   assert(stancestors != NULL);
   /* Kantenliste durchgehen
    */
   for(i = 0; i < slc; i++)
   {
      /* Hat t schon eine Kante mit diesem Ziel ?
       */
      for(et = p->outbeg[t]; et != EAT_LAST; et = p->oeat[et])
         if( p->head[et] == slp[i].knot )
            break;

      /* Keine gefunden, Kante aus der Liste muss eingefuegt werden.
       */
      if (et == EAT_LAST)
      {
         slp[i].mark = TRUE;
      }
      else
      {
         /* Ja ist vorhanden !
          */
         assert(et != EAT_LAST);

         /* Muessen die Kosten korrigiert werden ?
          * This is for nodes with edges to s and t.
          * Need to adjust the out and in costs of the edge
          */
         if( SCIPisGT(scip, p->cost[et], slp[i].outcost) )
         {
            SCIPindexListNodeFree(scip, &((p->ancestors)[et]));
            assert(p->ancestors[et] == NULL);
            SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &((p->ancestors)[et]), ancestors[i]) );
            SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &((p->ancestors)[et]), tsancestors) );
            p->cost[et] = slp[i].outcost;
         }
         if( SCIPisGT(scip, p->cost[Edge_anti(et)], slp[i].incost) )
         {
            anti = Edge_anti(et);
            SCIPindexListNodeFree(scip, &(p->ancestors[anti]));
            assert(p->ancestors[anti] == NULL);
            SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &((p->ancestors)[anti]), revancestors[i]) );
            SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &((p->ancestors)[anti]), stancestors) );
            p->cost[anti] = slp[i].incost;
         }
      }
   }

   /* Einzufuegenden Kanten einfuegen
    */
   for( i = 0; i < slc; i++ )
   {
      if( slp[i].mark )
      {
         es = p->outbeg[s];

         assert(es != EAT_LAST);
         assert(ancestors[i] != NULL);
         assert(revancestors[i] != NULL);
         SCIPindexListNodeFree(scip, &(p->ancestors[es]));
         SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(p->ancestors[es]), ancestors[i]) );
         SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(p->ancestors[es]), tsancestors) );
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
         SCIPindexListNodeFree(scip, &(p->ancestors[es]));

         SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(p->ancestors[es]), revancestors[i]) );
         SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(p->ancestors[es]), stancestors) );
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
      SCIPindexListNodeFree(scip, &(p->ancestors[es]));
      SCIPindexListNodeFree(scip, &(p->ancestors[Edge_anti(es)]));
      assert(p->ancestors[es] == NULL);
      assert(p->ancestors[Edge_anti(es)] == NULL);
      graph_edge_del(scip, p, es, FALSE);
   }

   SCIPindexListNodeFree(scip, &stancestors);
   SCIPindexListNodeFree(scip, &tsancestors);

   if( sgrad >= 2 )
   {
      for( i = 0; i < slc; i++ )
      {
         SCIPindexListNodeFree(scip, &(ancestors[i]));
         SCIPindexListNodeFree(scip, &(revancestors[i]));
      }
      free(ancestors);
      free(revancestors);
      free(slp);
   }
   assert(p->grad[s]   == 0);
   assert(p->outbeg[s] == EAT_LAST);
   assert(p->inpbeg[s] == EAT_LAST);
   return SCIP_OKAY;
}


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
   assert(SCIPisGE(scip, g->prize[i], 0));
   assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));
}

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
      int e;
      int j;

      for( e = g->outbeg[s]; e != EAT_LAST; e = g->oeat[e] )
         if( Is_pterm(g->term[g->head[e]]) )
            break;
      assert(e != EAT_LAST);
      j = g->head[e];

      assert(j != g->source[0]);
      assert(!g->mark[j]);

      graph_knot_chg(g, j, -1);
      graph_edge_del(scip, g, e, TRUE);

      for( e = g->inpbeg[j]; e != EAT_LAST; e = g->ieat[e] )
         if( g->source[0] == g->tail[e] )
            break;
      assert(e != EAT_LAST);
      assert(!g->mark[g->tail[e]]);

      assert(SCIPisEQ(scip, g->prize[s], g->cost[e]));

      prize_subtract(scip, g, g->cost[ets] - g->prize[s], i);
      graph_edge_del(scip, g, e, TRUE);

      SCIP_CALL( graph_knot_contract(scip, g, t, s) );
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


SCIP_RETCODE graph_edge_reinsert(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< the graph */
   int                   e1,
   int                   k1,
   int                   k2,
   SCIP_Real             cost,
   IDX*                  ancestors0,
   IDX*                  ancestors1,
   IDX*                  revancestors0,
   IDX*                  revancestors1
   )
{
   int n1;

   /* redirect; store new edge in n1 */
   n1 = graph_edge_redirect(scip, g, e1, k1, k2, cost);
   if( n1 >= 0 )
   {
      SCIPindexListNodeFree(scip, &(g->ancestors[n1]));
      SCIPindexListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));

      SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[n1]), revancestors0) );
      SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors1) );

      SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), ancestors0) );
      SCIP_CALL(  SCIPindexListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors1) );
   }
   return SCIP_OKAY;
}


/* ARGSUSED */
void graph_edge_add(
   GRAPH*                g,                  /**< the graph */
   int                   tail,               /**< tail of new edge */
   int                   head,               /**< head of new edge*/
   SCIP_Real             cost1,              /**< tail->head cost */
   SCIP_Real             cost2               /**< head->tail cost */
   )
{
   int    e;

   assert(g      != NULL);
   assert(GE(cost1, 0) || cost1 == UNKNOWN );
   assert(GE(cost2, 0) || cost2 == UNKNOWN );
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
      SCIPindexListNodeFree(scip, &((g->ancestors)[e]));
      SCIPindexListNodeFree(scip, &((g->ancestors)[Edge_anti(e)]));
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
{
   int head;
   int tail;
   int e;

   assert(g      != NULL);

   for(e = 0; e < g->edges; e++)
   {
      if (g->ieat[e] == EAT_HIDE)
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

/* unmark terminals and switch terminal property to orgininal terminals */
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
   //q = graph_init(knots, edges, p->layers, p->flags);
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
   /*SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(q->fixedges), p->fixedges);*/

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
         SCIPindexListNodeFree(scip, &(p->ancestors[i]));
         SCIPindexListNodeFree(scip, &(p->ancestors[i + 1]));
         continue;
      }
      assert(p->ieat[i]      != EAT_FREE);
      assert(p->oeat[i]      != EAT_FREE);
      assert(p->ieat[i + 1]  != EAT_FREE);
      assert(p->oeat[i + 1]  != EAT_FREE);
      assert(new[p->tail[i]] >= 0);
      assert(new[p->head[i]] >= 0);

      rcode = SCIPindexListNodeAppendCopy(scip, &(q->ancestors[q->edges]), p->ancestors[i]);
      rcode = SCIPindexListNodeAppendCopy(scip, &(q->ancestors[q->edges + 1]), p->ancestors[i + 1]);
      graph_edge_add(q, new[p->tail[i]], new[p->head[i]],
         p->cost[i], p->cost[Edge_anti(i)]);

      // SCIPindexListNodeFree(scip, &(p->ancestors[i]));
      // SCIPindexListNodeFree(scip, &(p->ancestors[i + 1]));
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
   //q = graph_init(knots, edges, g->layers, g->flags);
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

   if( new == NULL )
   {
      q->ancestors = NULL;
      graph_free(scip, g, FALSE);
      graph_knot_add(q, 0);
      q->source[0] = 0;
      return SCIP_OKAY;
   }

   q->ancestors = malloc((size_t)(nedges) * sizeof(IDX*));
   for( i = 0; i < nedges; i++ )
      q->ancestors[i] = NULL;

   /* add nodes (of positive degree) */
   for( i = 0; i < g->knots; i++ )
   {
      assert(g->term[i] < g->layers);
      if( g->grad[i] > 0 )
         graph_knot_add(q, g->term[i]);
   }

   /* add edges */
   for( i = 0; i < g->edges; i += 2 )
   {
      if( g->ieat[i] == EAT_FREE )
      {
         assert(g->oeat[i]     == EAT_FREE);
         assert(g->ieat[i + 1] == EAT_FREE);
         assert(g->oeat[i + 1] == EAT_FREE);
         SCIPindexListNodeFree(scip, &(g->ancestors[i]));
         SCIPindexListNodeFree(scip, &(g->ancestors[i + 1]));
         continue;
      }
      assert(g->ieat[i]      != EAT_FREE);
      assert(g->oeat[i]      != EAT_FREE);
      assert(g->ieat[i + 1]  != EAT_FREE);
      assert(g->oeat[i + 1]  != EAT_FREE);
      assert(new[g->tail[i]] >= 0);
      assert(new[g->head[i]] >= 0);

      SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(q->ancestors[q->edges]), g->ancestors[i]) );
      SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(q->ancestors[q->edges + 1]), g->ancestors[i + 1]) );
      graph_edge_add(q, new[g->tail[i]], new[g->head[i]], g->cost[i], g->cost[Edge_anti(i)]);
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

int graph_valid(
   const GRAPH*          p                  /**< the new graph */
   )
{
   const char* fehler1  = "*** Graph Validation Error: Head invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n";
   const char* fehler2  = "*** Graph Validation Error: Tail invalid, Knot %d, Edge %d, Tail=%d, Head=%d\n";
   const char* fehler3  = "*** Graph Validation Error: Source invalid, Layer %d, Source %d, Terminal %d\n";
   const char* fehler4  = "*** Graph Validation Error: FREE invalid, Edge %d/%d\n";
   const char* fehler5  = "*** Graph Validation Error: Anti invalid, Edge %d/%d, Tail=%d/%d, Head=%d/%d\n";
   const char* fehler6  = "*** Graph Validation Error: Knot %d with Grad 0 has Edges\n";
   const char* fehler7  = "*** Graph Validation Error: Knot %d not connected\n";
   const char* fehler8  = "*** Graph Validation Error: Wrong locals count, Layer %d, count is %d, should be %d\n";
   const char* fehler9  = "*** Graph Validation Error: Wrong Terminal count, count is %d, should be %d\n";

   int    k;
   int    l;
   int    e;
   int    terms;
   int*   locals;

   assert(p      != NULL);

   terms  = p->terms;
   locals = malloc((size_t)p->layers * sizeof(int));

   assert(locals != NULL);

   for( l = 0; l < p->layers; l++ )
      locals[l] = p->locals[l];

   for( k = 0; k < p->knots; k++ )
   {
      if( Is_term(p->term[k]) )
      {
         locals[p->term[k]]--;
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
      if( locals[l] != 0 )
         return((void)fprintf(stderr, fehler8,
               l, p->locals[l], p->locals[l] - locals[l]), FALSE);

      if( (p->source[l] < 0 )
         || (p->source[l] >= p->knots)
         || (p->term[p->source[l]] != l))
         return((void)fprintf(stderr, fehler3,
               l, p->source[l], p->term[p->source[l]]), FALSE);
   }
   free(locals);

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

char graph_sol_valid(
   const GRAPH*          graph,                  /**< the new graph */
   int*                  result
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

   terminal = malloc((size_t)nnodes * sizeof(char));
   for( i = 0; i < nnodes; i++ )
      terminal[i] = FALSE;
   /* BFS until all terminals are reached */
   SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2) );

   SCIP_CALL( SCIPqueueInsert(queue, &root) );
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
            SCIP_CALL( SCIPqueueInsert(queue, &graph->head[e]) );
         }
      }

   }

   if (termcount != graph->terms)
   {
      for( i = 0; i < nnodes; i++ )
         if( Is_term(graph->term[i]) && !terminal[i] )
            printf("not reached, node: %d\n", i);
      printf("a: %d, b: %d: \n", termcount, graph->terms);
   }

   free(terminal);
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

   terminal = malloc((size_t)nnodes * sizeof(char));
   reached = malloc((size_t)nnodes * sizeof(char));
   for( i = 0; i < nnodes; i++ )
   {
      terminal[i] = FALSE;
      reached[i] = FALSE;
   }
   /* BFS until all terminals are reached */
   SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2) );

   SCIP_CALL( SCIPqueueInsert(queue, &root) );
   termcount = 1;
   terminal[root] = TRUE;
   reached[root] = TRUE;
   while( !SCIPqueueIsEmpty(queue) )
   {
      pnode = (SCIPqueueRemove(queue));
      for( e = graph->outbeg[*pnode]; e != EAT_LAST; e = graph->oeat[e] )
      {

         if( SCIPisLT(scip, cost[e], 1e+10 - 10 ) && !reached[graph->head[e]] )
         {
            i = graph->head[e];
            reached[i] = TRUE;
            if( Is_term(graph->term[i]) )
            {
               assert(  !terminal[i]   );
               terminal[i] = TRUE;
               termcount++;

            }
            SCIP_CALL( SCIPqueueInsert(queue, &graph->head[e]) );
         }
      }

   }
   free(terminal);
   free(reached);
   SCIPqueueFree(&queue);
   if (termcount != graph->terms)
   {
      for( i = 0; i < nnodes; i++ )
         if( Is_term(graph->term[i]) && !terminal[i] )
            printf("not reached, node: %d\n", i);
      printf("a: %d, b: %d: \n", termcount, graph->terms);
      assert(0);
   }

   return (termcount == graph->terms);
}
