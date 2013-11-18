/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_tm.c
 * @brief  TM primal heuristic
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_tm.h"
#include "probdata_stp.h"
#include "grph.h"
#include "portab.h"


#define HEUR_NAME             "TM"
#define HEUR_DESC             "takahashi matsuyama primal heuristic for steiner trees"
#define HEUR_DISPCHAR         '+'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_EVALRUNS 10
#define DEFAULT_INITRUNS 100 /* TODO CHG TO 100*/
#define DEFAULT_LEAFRUNS 10
#define DEFAULT_ROOTRUNS 50
#define DEFAULT_DURINGLPFREQ 10


/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint ncalls;
   int evalruns;
   int initruns;
   int leafruns;
   int rootruns;
   int duringlpfreq;
};

typedef struct ST_Node
{
   int edge;
   struct ST_Node *parent;
} NODE;

typedef struct Index_List_Node
{
   int index;
   struct Index_List_Node *parent;
} IDX;

#if 0
/** primal heuristic data */
typedef struct ST_Node
{
   int s;
   int ls;
   int id;
   SCIP_Bool flip;
   SCIP_Bool lflip;
   struct ST_Node *left;
   struct ST_Node *right;
   struct ST_Node *parent;
} NODE;


/*
 * Local methods
 */

/* Die Heuristic stoert sich nicht dran, wenn sie einzelne Wege nicht
 * routen kann, sondern erklaert das jeweilige Netz einfach fuer fertig.
 * Die Loesung muss also ueberprueft werden.
 */

static void createNode2(
   NODE**          k,
   int            status,
   int            id
   )
{
   (*k)->s = status;
   (*k)->ls = status;
   (*k)->id = id;
   (*k)->left = (*k)->right = (*k)->parent = NULL;
   (*k)->flip = (*k)->lflip = FALSE;
   //  return k;
}

static void createNode(
   NODE*          k,
   int            status,
   int            id
   )
{ /* TODO really create....*/
   (k)->s = status;
   (k)->ls = status;
   (k)->id = id;
   (k)->left = (k)->right = (k)->parent = NULL;
   (k)->flip = (k)->lflip = FALSE;
   //  return k;
}
static SCIP_Bool isRoot(
   NODE* v
   )
{
   NODE* p = k->parent;
   return (p == NULL || !(p->left == v || p->right == v));
}

static void adjustNode(
   NODE* k
   )
{
   if( k->flip )
   {
      k->flip = FALSE;
      k->lflip = !k->lflip;
      if( k->right != NULL )
      {
         k->right->flip = !k->right->flip;
      }
      if( k->left != NULL )
      {
         k->left->flip = !k->left->flip;
      }
   }
}

static void up(
   NODE* v
   )
{
   v->s = v->ls;
   if( v->left != NULL && v->left->s > v->s)
      v->s = v->left->s;
    if( v->right != NULL && v->right->s > v->s)
      v->s = v->right->s;
}

static void rotR(
   NODE* k
   )
{
   NODE* q = k->parent;
   NODE* r = q->parent;
   adjustNode(q);
   adjustNode(k);
   if( (q->left = k->right) != NULL )
      q->left->parent = q;
   k->right = q;
   q->parent = k;
   if( (k->parent = r) != NULL )
   {
      if( r->left == q )
         r->left =k;
      else if( r->right == q )
         r->right = k;
   }
   up(q);
}

static void rotL(
   NODE* k
   )
{
   NODE* q = k->parent;
   NODE* r = q->parent;
   adjustNode(q);
   adjustNode(k);
   if( (q->right = k->left) != NULL )
      q->right->parent = q;
   k->left = q;
   q->parent = k;
   if( (k->parent = r) != NULL )
   {
      if( r->left == q )
         r->left =k;
      else if( r->right == q )
         r->right = k;
   }
   up(q);
}
static void splay(
   NODE* k
   )
{
   NODE* q;
   NODE* r;
   while( !isRoot(k) )
   {
      q = k->parent;
      if( isRoot(q) )
      {
         if( q->left == k )
            rotR(k);
         else
            rotL(k);
      }
      else
      {
         r = q->parent;
         if( r->left == q)
         {
            if( q->left == k)
            {
               rotR(q);
               rotR(k);
            }
            else
            {
               rotL(k);
               rotR(k);
            }
         }
         else
         {
            if( q->right == k)
            {
               rotL(q);
               rotL(k);
            }
            else
            {
               rotR(k);
               rotL(k);
            }
         }
      }
   }
   up(k);
}

static void access(
   NODE* k
   )
{
   NODE* r = NULL;
   NODE* q;
   for( q = k; q != NULL; q = q->parent )
   {
      splay(q);
      q->left = r;
      up(q);
      r = q;
   }
   splay(k);
   assert( k->parent == NULL );
}

static void link(
   NODE* k,
   NODE* q
   )
{
   access(k);
   assert( k->right == NULL );
   k->parent = q;
}

static void cut(
   NODE* v
   )
{
   access(v);
   v->right = NULL;
}



static void update(
   NODE* k
   )
{
   k->s = k->ls;

   if( k->left != NULL )
   {
      k->s += k->left->s;
      if( k->left->flip )
      {

      }
      else
      {

      }
   }
}

static void rotRed(
         NODE* p
            )
         {
         TRUE;
            }

#endif

static int flipEdge(
    int edge
	   )
{
  assert(edge >= 0);
  return ((edge % 2) == 0) ? edge + 1 : edge - 1;
}
static void init(
   NODE* v
   )
{
   v->parent = NULL;
   v->edge = -1;
}

static void link(
   NODE* v,
   NODE* w,
   int edge
   )
{

   assert( v->parent == NULL );
   v->parent = w;
   v->edge = edge;
}

static void cut(
   NODE* v
   )
{
   v->edge = -1;
   v->parent = NULL;
}

static NODE* findMax(
   SCIP* scip,
   const GRAPH* graph,
   NODE* v
   )
{
   NODE* p = v;
   NODE* q = NULL;
   SCIP_Real max = -1;
   while( p != NULL )
   {
      if( SCIPisGE(scip, graph->cost[p->edge], max) )
      {
         max = graph->cost[p->edge];
         q = p;
      }
      p = p->parent;
   }
   return q;
}

static void evert(
   NODE* v
   )
{
   NODE* p = NULL;
   NODE* q = v;
   NODE* r;

   int val = -1;
   int tmpval;
   while( q != NULL )
   {
      r = q->parent;
      tmpval =  q->edge;
      q->edge = (val >= 0) ? flipEdge(val) : val;
      val = tmpval;
      q->parent = p;
      p = q;
      q = r;
   }
}


static
SCIP_RETCODE do_heuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*  g,
   int           layer,
   int*          result,
   int           start,
   char*         connected,
   const SCIP_Real* cost,
   PATH**        path
   )
{
   PATH*  path1;
   PATH*  mst;
   int*   cluster;
   int    csize = 0;
   int    k;
   int    e;
   int    count;
   SCIP_Real min;
   int    i;
   int    j;
   int    old;
   int    newval;
   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(path      != NULL);
   assert(layer >= 0 && layer < g->layers);


   SCIPdebugMessage("Heuristic: Start=%5d ", start);


   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, g->knots) );
   cluster[csize++] = start;
#if 0
   SCIP_CALL( SCIPallocBufferArray(scip, &path1, g->knots) );
    for( i = 0; i < g->knots; i++ )
   {
      g->mark[i]   = (g->grad[i] > 0);
   }
   graph_path_exec2(g, FSP_MODE, start, cost, path1, connected, cluster, &csize);
   SCIPfreeBufferArray(scip, &path1);

   printf("ff %d \n", csize);

#else


   for( i = 0; i < g->knots; i++ )
   {
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
   }
   connected[start] = TRUE;


   /* CONSTCOND */
   for(;;)
   {
      /* Suche das Terminal, das am dichtesten dran ist
       */
      min = FARAWAY;
      old = -1;
      newval = -1;

      for(i = 0; i < g->knots; i++)
      {
         if (g->grad[i] == 0)
            continue;

         if (g->term[i] != layer)
            continue;

         if (connected[i])
            continue;

         /* Jetzt brauchen wir die Entfernungen.
          */
         if (path[i] == NULL)
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(path[i]), g->knots) );

            assert(path[i] != NULL);

            /* ! Ob das die guenstiges Richtung ist die Wege zu berechnen, wenn
             * ! die Kosten fuer Hin und Rueckweg unterschiedlich sind ist doch
             * ! sehr fraglich.
             * ! Koennte aber sein, weil wir die Kosten unten umgedreht haben.
             */
            graph_path_exec(g, FSP_MODE, i, cost, path[i]);
         }
         for(k = 0; k < csize; k++)
         {
            j = cluster[k];

            assert(i != j);
            assert(connected[j]);

            if (LT(path[i][j].dist, min))
            {
               min = path[i][j].dist;
               newval = i;
               old = j;
            }
         }
      }
      /* Nichts mehr gefunden, also fertig
       */
      if (newval == -1)
         break;

      /* Weg setzten
       */
      assert((old > -1) && (newval > -1));
      assert(path[newval] != NULL);
      assert(path[newval][old].dist < FARAWAY);
      assert(g->term[newval] == layer);
      assert(!connected[newval]);
      assert(connected[old]);

      SCIPdebug(fputc('R', stdout));
      SCIPdebug(fflush(stdout));

      /*    printf("Connecting Knot %d-%d dist=%d\n", newval, old, path[newval][old].dist);
       */
      /* Gegen den Strom schwimmend alles markieren
       */
      k = old;

      while(k != newval)
      {
         e = path[newval][k].edge;
         k = g->tail[e];

         if (!connected[k])
         {
            connected[k] = TRUE;
            cluster[csize++] = k;
         }
      }
   }


   SCIPdebug(fputc('M', stdout));
   SCIPdebug(fflush(stdout));
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, g->knots) );

   /* MST berechnen
    */
   for(i = 0; i < g->knots; i++)
      g->mark[i] = connected[i];

   assert(g->source[layer] >= 0);
   assert(g->source[layer] <  g->knots);

   graph_path_exec(g, MST_MODE, g->source[layer], g->cost, mst);

   for(i = 0; i < g->knots; i++)
   {
      if (connected[i] && (mst[i].edge != -1))
      {
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == -1);

         result[mst[i].edge] = layer;
      }
   }

   /* Baum beschneiden
    */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for(i = 0; i < g->knots; i++)
      {
         if (!g->mark[i])
            continue;

         if (g->term[i] == layer)
            continue;

         for(j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j])
            if (result[j] == layer)
               break;

         if (j == EAT_LAST)
         {
            /* Es muss genau eine eingehende Kante geben
             */
            for(j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j])
            {
               if (result[j] == layer)
               {
                  result[j]    = -1;
                  g->mark[i]   = FALSE;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
            assert(j != EAT_LAST);
         }
      }
   }


   while(count > 0);


   /* TEST!!!!!!
      TRYING TO DO A VERTEX INSERTION
   */


   if( 0 )
   {
      NODE* nodes;
      IDX* cand_start;
      IDX* cand_curr;
      IDX* add_start;
      IDX* add_curr;
      IDX* cut_start;
      IDX* cut_curr;
      int oedge;
      NODE* v;
      NODE* w;
      NODE* max;
      SCIP_Real diff;
      SCIP_CALL( SCIPallocBufferArray(scip, &nodes, g->knots) );
      for( i = 0; i< g->knots; i++ )
         init(&nodes[i]);

      for( i = 0; i < g->edges; i++ )
         assert( g->head[i] == g->tail[flipEdge(i)] );
      /* create a link-cut tree representing the current Steiner tree */
      for( i = 0; i < g->knots; i++ )
      {
         if( connected[i] && mst[i].edge != -1  )
         {

            //printf( "IINNNN %d __ %d\n", i, g->head[mst[i].edge] );
            assert( connected[g->head[mst[i].edge]] );
            //nodes[i].edge = flipEdge(mst[i].edge);

            link(&nodes[i], &nodes[g->tail[mst[i].edge]], flipEdge(mst[i].edge));

            assert(g->tail[nodes[i].edge] == i );

            if ( !(result[mst[i].edge] != -1) )
               printf("damnff %d\n", mst[i].edge) ;
         }
      }

      for( i = 0; i < g->knots; i++ )
      {
         if( !connected[i] && g->grad[i] > 1 )
         {
            cand_start = NULL;
            /* if an outgoing edge of node i points to the current tree, link the edge to a list */
            oedge = g->outbeg[i];
            assert( g->tail[oedge] == i);
            while( oedge >= 0 )
            {
               assert( g->tail[oedge] == i);
               if( connected[g->head[oedge]] )
               {
                  cand_curr = (IDX *)malloc(sizeof(IDX)); //TODO extra method
                  cand_curr->index = oedge;
                  cand_curr->parent  = cand_start;
                  cand_start = cand_curr;
               }
               oedge = g->oeat[oedge];
            }


            /* if there are at least two edges between node i and the current tree, start the insertion process */
            if( cand_start != NULL && cand_start->parent != NULL )
            {
               diff = 0;
               //  printf("try to insert node %d \n", i);
               /* the node to insert */
               v = &nodes[i];

               link(v, &nodes[g->head[cand_start->index]], cand_start->index);
               // printf(" firstedge: %d->%d \n", g->tail[cand_start->index], g->head[cand_start->index]);

               diff +=  g->cost[v->edge];
               cand_curr = cand_start->parent;

               cut_start = NULL;
               add_start = NULL;

               while( cand_curr != NULL )
               {
                  evert(v);
                  w = &nodes[g->head[cand_curr->index]];


                  /* if there is an edge with cost greater than that of the current edge... */
                  max = findMax(scip, g, w);
                  if( SCIPisGT(scip, g->cost[max->edge], g->cost[cand_curr->index]) )
                  {
                     diff += g->cost[cand_curr->index];
                     diff -= g->cost[max->edge];
                     //	 printf(" addedge: %d->%d cutedge %d->%d \n", g->tail[cand_curr->index], g->head[cand_curr->index], g->tail[max->edge], g->head[max->edge]);
                     cut_curr = (IDX *)malloc(sizeof(IDX));
                     cut_curr->index = max->edge;
                     cut_curr->parent  = cut_start;
                     cut_start = cut_curr;

                     cut(max);
                     link(v, w, cand_curr->index);

                     add_curr = (IDX *)malloc(sizeof(IDX));
                     add_curr->index = v->edge;
                     add_curr->parent  = add_start;
                     add_start = add_curr;
                  }
                  cand_curr = cand_curr->parent;
               }

               /* if the new tree is more expensive than the old one, restore the latter */
               if( SCIPisPositive(scip, diff) )
               {
		  //   printf("positive, restore \n");
                  cut_curr = cut_start;
                  add_curr = add_start;

                  evert(v);
                  while( cut_curr != NULL )
                  {
                     cut(&nodes[g->head[add_curr->index]]);
                     evert(&nodes[g->tail[cut_curr->index]]);
                     link(&nodes[g->tail[cut_curr->index]], &nodes[g->head[cut_curr->index]], cut_curr->index);
                     cut_curr = cut_curr->parent;
                     add_curr = add_curr->parent;
                  }

                  /* finally, cut the edge added first */
                  evert(v);
                  //      printf(" firstedge: %d->%d \n", g->tail[cand_start->index], g->head[cand_start->index]);
		  //    printf(" 7 = ? %d \n",  g->head[cand_start->index] );
                  //   printf( "head: %d \n", g->head[nodes[g->head[cand_start->index]].edge] );
                  assert( g->head[nodes[g->head[cand_start->index]].edge] == i); // || (0 &&v->edge == cand_start->index));

                  cut(&nodes[g->head[cand_start->index]]);
               }
               else /* TODO adjust tree st we only have to adjust result for the new edges*/
               {
                  cluster[csize++] = i;
                  connected[i] = TRUE;
               }

               add_curr = add_start;
               while( add_curr != NULL ){
		  add_start = add_curr->parent;
		  free( add_curr );
                  add_curr = add_start;
               }
               cut_curr = cut_start;
               while( cut_curr != NULL ){
		  cut_start = cut_curr->parent;
		  free( cut_curr );
                  cut_curr = cut_start;
               }


            }
            cand_curr = cand_start;
            while( cand_curr != NULL ){
               cand_start = cand_curr->parent;
               free( cand_curr );
               cand_curr = cand_start;
            }



         }

      }

      evert(&nodes[ g->source[0]] );
      for(e = 0; e < g->edges; e++)
      {
         result[e] = -1;
      }
      for(e = 0; e < csize; e++)
      {

         if( nodes[cluster[e]].edge != -1){
            result[flipEdge(nodes[cluster[e]].edge)] = 0;
            //  printf("resulta: %d-->%d \n", g->tail[flipEdge(nodes[cluster[e]].edge)] , g->head[flipEdge(nodes[cluster[e]].edge)]);
         }

      }
      /*for(e = 0; e < g->edges; e++)
        printf("resultb: %d-->%d \n", g->tail[e] , g->head[e]); */
      SCIPfreeBufferArray(scip, &nodes);
   }

   //printf("durchgang \n");






   if( 0 )
   {

  }


# if 0
   for( i = 0; i < g->knots; i++ )
   {

      if( connected[i] == TRUE )
         printf(" %d \n", i);
   }
   printf("------------ \n");
#endif
   SCIPfreeBufferArray(scip, &mst);

   SCIPfreeBufferArray(scip, &cluster);
   return SCIP_OKAY;
}

static
SCIP_RETCODE do_layer(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*  graph,
   int           layer,
   int*          best_result,
   int           runs,
   const SCIP_Real* cost
   )
{
   PATH** path;
   char* connected;
   int* result;
   int* start;
   SCIP_Real obj;
   SCIP_Real min = FARAWAY;
   int best = -1;
   int k;
   int r;
   int e;
   int nnodes;
   assert(scip != NULL);
   assert(graph != NULL);
   assert(best_result != NULL);
   assert(cost != NULL);
   assert(layer >= 0 && layer < graph->layers);

   nnodes = graph->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &start, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &result, graph->edges) );
#if 1
   BMSclearMemoryArray(path, nnodes);
#else
   for( k = 0; k < nnodes; k++ )
      path[k] = NULL;
#endif

   /* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    * Patch um die heuristic nach einem restruct starten zu koennen
    * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    */
   if( best >= nnodes )
      best = -1;

   if( graph->layers > 1 )
   {
      SCIP_CALL( do_heuristic(scip, graph, layer, best_result, graph->source[layer], connected, cost, path) );
   }
   else
   {
      runs = (runs > nnodes) ? nnodes         : runs;
      best = (best < 0)            ? graph->source[layer] : best;

      for( k = 0; k < nnodes; k++ )
      {
         assert(graph->grad[k] > 0);

         start[k] = k;
      }

      /* if we run over all nodes, we do not need to do the following */
      if( runs < nnodes )
      {
         int random;
         int tmp;

#if 0
         /* swap the starting values randomly
         for( r = 0; r < runs; r++ )
         {
            random   = rand() % nnodes;
            tmp      = start[r];
            start[r] = start[random];
            start[random] = tmp;
         }*/
#else

	 int* realterms = SCIPprobdataGetRTerms(scip);


	 int nrealterms = SCIPprobdataGetRNTerms(scip);
	 start[0] = graph->source[0];
         for( r = 1; r < runs; r++ )
         {
	    if( r < nrealterms + 1 )
	    {
	       start[r] = realterms[r - 1];
               start[realterms[r - 1]] = r;
	    }
	    else
	    {
	       break;
	    }
	 }
 #endif
         /* check if we have a best starting value */
         for( r = 0; r < runs; r++ )
            if( start[r] == best )
               break;

         /* do we need to set the start by hand */
         if( r == runs )
            start[0] = best;
      }
      else
      {
	 runs = nnodes;
      }

      for( r = 0; r < runs; r++ )
      {
         /* incorrect if layers > 1 ! */
         assert(graph->layers == 1);

         for( e = 0; e < graph->edges; e++ )
            result[e] = -1;

         SCIP_CALL( do_heuristic(scip, graph, layer, result, start[r], connected, cost, path) );

         obj = 0.0;

         /* here we take another measure than in do_heuristic() */
         for( e = 0; e < graph->edges; e++)
            obj += (result[e] > -1) ? graph->cost[e] : 0.0;

         //SCIPdebugMessage(" Obj=%.12e\n", obj);

         if( LT(obj, min) )
         {
            min = obj;
            printf(" Obj=%.12e\n", obj);
            for( e = 0; e < graph->edges; e++ )
               best_result[e] = result[e];

            best = start[r];
         }
      }

   }




   /*** VERTEX  INSERTION *//////
   if( 1 )
   {
      int i;
      int newnode = 0;
      int counter;
      NODE* nodes;
      IDX* cand_start;
      IDX* cand_curr;
      IDX* add_start;
      IDX* add_curr;
      IDX* cut_start;
      IDX* cut_curr;
      int oedge;
      NODE* v;
      NODE* w;
      NODE* max;
      char* steinertree;
      SCIP_Real diff;
      SCIP_CALL( SCIPallocBufferArray(scip, &nodes, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &steinertree, nnodes) );
      for( i = 0; i< nnodes; i++ )
      {
         steinertree[i] = FALSE;
         init(&nodes[i]);
      }

      /* create a link-cut tree representing the current Steiner tree */
      for( e = 0; e < graph->edges; e++ )
      {
         assert(graph->head[e] == graph->tail[flipEdge(e)]);

         /* if edge e is in the tree, so are its incident vertices */
         if( best_result[e] != -1 )
         {
            steinertree[graph->tail[e]] = TRUE;
            steinertree[graph->head[e]] = TRUE;
            link(&nodes[graph->head[e]], &nodes[graph->tail[e]], flipEdge(e));
         }
      }


      /* vertex insertion */
    //  for( counter = 0; counter < 3; counter++ )
     // {
      i = 0;
      for( ;; )
      {
         if( !steinertree[i] && graph->grad[i] > 1 )
         {
            cand_start = NULL;

            /* if an outgoing edge of node i points to the current tree, link the edge to a list */
            oedge = graph->outbeg[i];
            assert( graph->tail[oedge] == i);
            while( oedge >= 0 )
            {
               assert( graph->tail[oedge] == i);
               if( steinertree[graph->head[oedge]] )
               {
                  cand_curr = (IDX *)malloc(sizeof(IDX)); //TODO extra method
                  cand_curr->index = oedge;
                  cand_curr->parent  = cand_start;
                  cand_start = cand_curr;
               }
               oedge = graph->oeat[oedge];
            }


            /* if there are at least two edges connecting node i and the current tree, start the insertion process */
            if( cand_start != NULL && cand_start->parent != NULL )
            {
               diff = 0;

               /* the node to insert */
               v = &nodes[i];

               link(v, &nodes[graph->head[cand_start->index]], cand_start->index);
               diff +=  graph->cost[v->edge];
               cand_curr = cand_start->parent;

               cut_start = NULL;
               add_start = NULL;

               while( cand_curr != NULL )
               {
                  evert(v);

                  /* next vertex in the current Steiner tree adjacent to vertex i (the one being scrutinized for possible insertion) */
                  w = &nodes[graph->head[cand_curr->index]];


                  /* if there is an edge with cost greater than that of the current edge... */
                  max = findMax(scip, graph, w);
                  if( SCIPisGT(scip, graph->cost[max->edge], graph->cost[cand_curr->index]) )
                  {
                     diff += graph->cost[cand_curr->index];
                     diff -= graph->cost[max->edge];

                     cut_curr = (IDX *)malloc(sizeof(IDX));
                     cut_curr->index = max->edge;
                     cut_curr->parent  = cut_start;
                     cut_start = cut_curr;

                     cut(max);
                     link(v, w, cand_curr->index);

                     add_curr = (IDX *)malloc(sizeof(IDX));
                     add_curr->index = v->edge;
                     add_curr->parent  = add_start;
                     add_start = add_curr;
                  }
                  cand_curr = cand_curr->parent;
               }

               add_curr = add_start;

               /* if the new tree is more expensive than the old one, restore the latter */
               if( SCIPisPositive(scip, diff) )
               {
                  cut_curr = cut_start;
                  evert(v);
                  while( cut_curr != NULL )
                  {
                     cut(&nodes[graph->head[add_curr->index]]);
                     evert(&nodes[graph->tail[cut_curr->index]]);
                     link(&nodes[graph->tail[cut_curr->index]], &nodes[graph->head[cut_curr->index]], cut_curr->index);
                     cut_curr = cut_curr->parent;
                     add_curr = add_curr->parent;
                  }

                  /* finally, cut the edge added first */
                  evert(v);
                  cut(&nodes[graph->head[cand_start->index]]);
               }
               else
               {
		  newnode = i;
                  steinertree[i] = TRUE;
		  /* TODO adjust tree st we only have to adjust best_result for the new edges*/
#if 0
                  assert(best_result[flipEdge(cand_start->index)] == - 1);
                  best_result[flipEdge(cand_start->index)] = 0;
                  while( add_curr != NULL ){
                     printf("inserted %d->%d \n", graph->tail[add_curr->index], graph->head[add_curr->index]);
                     assert(best_result[flipEdge(add_curr->index)] == - 1);
                     best_result[flipEdge(add_curr->index)] = 0;
                     add_curr = add_curr->parent;
                  }

                  cut_curr = cut_start;
                  while( cut_curr != NULL ){
                     printf("cut %d->%d \n", graph->tail[cut_curr->index], graph->head[cut_curr->index]);
                     assert(best_result[flipEdge(cut_curr->index)] != -1 || best_result[(cut_curr->index)] != -1);
                     if( best_result[flipEdge(cut_curr->index)] != -1)
                        best_result[flipEdge(cut_curr->index)] = -1;
                     else
                        best_result[cut_curr->index] = -1;
                     cut_curr = cut_curr->parent;
                  }
#endif
               }

               add_curr = add_start;
               while( add_curr != NULL ){
		  add_start = add_curr->parent;
		  free( add_curr );
                  add_curr = add_start;
               }
               cut_curr = cut_start;
               while( cut_curr != NULL ){
		  cut_start = cut_curr->parent;
		  free( cut_curr );
                  cut_curr = cut_start;
               }


            }
            cand_curr = cand_start;
            while( cand_curr != NULL ){
               cand_start = cand_curr->parent;
               free( cand_curr );
               cand_curr = cand_start;
            }


         }
         if( i < nnodes - 1 )
	 {
	    i++;
	 }
	 else
	 {
	    printf("VertInsert newrun \n");
	    i = 0;
	 }
	 if( newnode == i )
	 {
	    break;
	 }

      }


    //  } //counter
      evert(&nodes[ graph->source[0]] );

      for( e = 0; e < graph->edges; e++ )
      {
         best_result[e] = -1;
      }
      for( i = 0; i < nnodes; i++ )
      {

         if( steinertree[i] && nodes[i].edge != -1){
            best_result[flipEdge(nodes[i].edge)] = 0;

         }
      }

      SCIPfreeBufferArray(scip, &nodes);
      SCIPfreeBufferArray(scip, &steinertree);

      obj = 0.0;
         for( e = 0; e < graph->edges; e++)
            obj += (best_result[e] > -1) ? graph->cost[e] : 0.0;

      printf(" ObjAfterVertexInsertion=%.12e\n", obj);
   }



   for( k = 0; k < nnodes; k++ )
   {
      assert(path[k] == NULL || graph->term[k] == layer);
      SCIPfreeBufferArrayNull(scip, &(path[k]));
   }


   SCIPfreeBufferArray(scip, &result);
   SCIPfreeBufferArray(scip, &start);
   SCIPfreeBufferArray(scip, &connected);
   SCIPfreeBufferArray(scip, &path);

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

   /* @todo copy heuristic? (probdata needs to be copied as well) */
#if 0
   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTM(scip) );
#endif
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

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   heurdata->ncalls = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitTM)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of TM primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitTM NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolTM)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of TM primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolTM NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolTM)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of TM primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolTM NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTM)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_PROBDATA* probdata;
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol;
   GRAPH* graph;
   SCIP_Real* cost;
   SCIP_Real* nval;
   SCIP_Real* xval;
   int* results;
   SCIP_Real pobj;
   int nvars;
   int layer;
   int runs;
   int e;
   int v;

   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DELAYED;
   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);

   runs = 0;

   if( heurtiming & SCIP_HEURTIMING_BEFORENODE )
      runs = heurdata->initruns;
   else if( ((heurtiming & SCIP_HEURTIMING_DURINGLPLOOP) && (heurdata->ncalls % heurdata->duringlpfreq == 0)) || (heurtiming & SCIP_HEURTIMING_AFTERLPLOOP) )
      runs = heurdata->evalruns;
   else if( heurtiming & SCIP_HEURTIMING_AFTERNODE )
   {
      if( SCIPgetDepth(scip) == 0 )
         runs = heurdata->rootruns;
      else
         runs = heurdata->leafruns;
   }

   heurdata->ncalls++;

   if( runs == 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("Heuristic Start\n");

   nvars = SCIPprobdataGetNVars(scip);
   vars = SCIPprobdataGetVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &cost, graph->edges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &results, graph->edges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nval, nvars) );

   *result = SCIP_DIDNOTFIND;

   /* */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      sol = NULL;
      xval = NULL;
   }
   else
   {
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

      /* copy the current LP solution to the working solution */
      SCIP_CALL( SCIPlinkLPSol(scip, sol) );

      xval = SCIPprobdataGetXval(scip, sol);

      SCIPfreeSol(scip, &sol);
   }

   for( e = 0; e < graph->edges; e++ )
      results[e] = -1;

   for( layer = 0; layer < graph->layers; layer++ )
   {
      if( xval == NULL )
      {
         BMScopyMemoryArray(cost, graph->cost, graph->edges);
      }
      else
      {
         /* swap costs; set a high cost if the variable is fixed to 0 */
         for( e = 0; e < graph->edges; e += 2)
         {
            if( SCIPvarGetUbLocal(vars[layer * graph->edges + e + 1]) < 0.5 )
               cost[e] = 1e+10; /* ???? why does FARAWAY/2 not work? */
            else
               cost[e] = ((1.0 - xval[layer * graph->edges + e + 1]) * graph->cost[e + 1]);

            if( SCIPvarGetUbLocal(vars[layer * graph->edges + e]) < 0.5 )
               cost[e + 1] = 1e+10; /* ???? why does FARAWAY/2 not work? */
            else
               cost[e + 1] = ((1.0 - xval[layer * graph->edges + e]) * graph->cost[e]);
         }
      }
      /* can we connect the network */
      SCIP_CALL( do_layer(scip, graph, layer, results, runs, cost) );

      /* take the path */
      if( graph->layers > 1 )
      {
         for( e = 0; e < graph->edges; e += 2)
         {
            if( (results[e] == layer) || (results[e + 1] == layer) )
               graph_edge_hide(graph, e);
         }
      }
   }

   if( graph->layers > 1 )
      graph_uncover(graph);

   for( v = 0; v < nvars; v++ )
      nval[v] = (results[v % graph->edges] == (v / graph->edges)) ? 1.0 : 0.0;

   if( validate(graph, nval) )
   {
      pobj = 0.0;

      for( v = 0; v < nvars; v++ )
         pobj += graph->cost[v % graph->edges] * nval[v];

      if( SCIPisLT(scip, pobj, SCIPgetPrimalbound(scip)) )
      {
         SCIP_Bool success;

         SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

         if( success )
            *result = SCIP_FOUNDSOL;
      }
   }

   SCIPfreeBufferArray(scip, &nval);
   SCIPfreeBufferArray(scip, &results);
   SCIPfreeBufferArray(scip, &cost);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the TM primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTM(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

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
#if 0
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitTM) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolTM) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolTM) );
#endif

   /* add TM primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   /* add TM primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/evalruns",
         "number of runs for eval",
         &heurdata->evalruns, FALSE, DEFAULT_EVALRUNS, -1, INT_MAX, NULL, NULL) );
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

   return SCIP_OKAY;
}
