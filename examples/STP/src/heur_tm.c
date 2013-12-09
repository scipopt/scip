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

#define DEFAULT_EVALRUNS 10 /*10*/
#define DEFAULT_INITRUNS 100 /* TODO CHG TO 100*/
#define DEFAULT_LEAFRUNS 10 /*10*/
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

typedef struct Restore_List_Node
{
   int index;
   int pred;
   SCIP_Real dist;
   struct Index_List_Node *next;
} RES;

typedef struct PairingHeap_Node
{
   int size;
   int element;
   SCIP_Real key;
   struct PairingHeap_Node* child;
   struct PairingHeap_Node* sibling;
   struct PairingHeap_Node* prev;

}PHeap;


static int* dfstree;


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

static void IDX_freeList(
   SCIP* scip,
   IDX* start
   )
{
   IDX* curr = start;
   if( start == NULL )
      return;

   while( curr != NULL ){
      start = curr->parent;
      SCIPfreeMemory(scip, &curr);
      curr = start;
   }

}


static void IDX_add(
   SCIP* scip,
   IDX* start,
   IDX* curr,
   int index
   )
{

   /*IDX* curr = NULL;
     SCIP_CALL( SCIPallocMemory(scip, &curr) ); */
   curr->index = index;
   curr->parent  =  start;
   start = curr;

}
static void crucList(
   const GRAPH* g,
   int* edges,
   int* node,
   int* counter
   )
{
   int tmp;

   int oedge = g->outbeg[*node];

   while( oedge >= 0 )
   {
      if( edges[oedge] >= 0 )
      {

         crucList(g, edges, &(g->head[oedge]), counter);
      }
      oedge = g->oeat[oedge];
   }
   dfstree[*counter] = *node;
   tmp = *counter;
   tmp++;
   *counter = tmp;

}

static char nodeIsCrucial(
   const GRAPH* graph,
   char* steinertree,
   int node
   )
{
   assert( graph != NULL );
   assert( steinertree != NULL );

   if( graph->term[node] == -1 )
   {
      int counter = 0;
      int e = graph->outbeg[node];
      while( e >= 0 )
      {

         /* check if the adjacent node is in the ST */
         if( steinertree[graph->head[e]] == TRUE )
         {
            counter++;
         }
         e = graph->oeat[e];
      }

      if( counter < 3 )
      {
         return FALSE;
      }
   }
   return TRUE;

}
static int flipEdge(
   int edge
   )
{
   assert(edge >= 0);
   return ((edge % 2) == 0) ? edge + 1 : edge - 1;
}


/*** Linear Link Cut Tree ***/
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

   assert(v != NULL);

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



/*** Pairing Heap ***/
static PHeap** PHeap_doubleArr(
   SCIP* scip,
   PHeap** arr,
   int length
   )
{
   int i;
   // SCIP_RETCODE* ret;

   // PHeap** oldArr = arr;
   PHeap* newarr[2 * length];

   //SCIP_CALL( SCIPallocMemoryArray(scip, &arr, 2 * length) );
   for( i = 0; i < length; i++ )
   {
      newarr[i] = arr[i];
   }
   arr = newarr;
   return arr;
}



static PHeap* PHeap_mergeheaps(
   SCIP* scip,
   PHeap *p,
   PHeap *q
   )
{
   if( q == NULL )
      return p;


   if( p->key <= q->key )
   {
      q->prev = p;
      p->sibling = q->sibling;
      if(p->sibling != NULL)
         p->sibling->prev = p;

      q->sibling = p->child;
      if(q->sibling != NULL)
         q->sibling->prev = q;

      p->child = q;
      //q = NULL;
      return p;
   }
   else
   {
      q->prev = p->prev;
      p->prev = q;
      p->sibling = q->child;
      if( p->sibling != NULL )
         p->sibling->prev = p;

      q->child = p;
      //p = NULL;
      return q;
   }
}

static PHeap* PHeap_combine_siblings(
   SCIP* scip,
   PHeap *p
   )
{
   //SCIP_RETCODE retcode;
   //PHeap** treearray;
   PHeap* q;
   //retcode = SCIPallocBufferArray(scip, &treearray, p->size);
   PHeap* treearray[p->size];
   //assert(retcode == SCIP_OKAY);

   int i;
   int j;
   int nsiblings;

   if( p->sibling == NULL )
      return p;

   for( nsiblings = 0; p != NULL; nsiblings++ )
   {
      /* if the last entry is reached, double the size of the array */
      if( nsiblings == p->size - 1 )
      {
         //treearray = PHeap_doubleArr(scip, treearray, p->size);
         p->size = p->size * 2;
      }
      treearray[nsiblings] = p;
      p->prev->sibling = NULL;
      p = p->sibling;
   }
   treearray[nsiblings] = NULL;

#if 0
   // combine the subtrees (simple)
   for(i = 1; i < nsiblings; i++)
      treearray[i] = PHeap_mergeheaps(scip, treearray[i-1], treearray[i]);


   return treearray[nsiblings-1];
#endif

   // combine the subtrees (two at a time)
   for( i = 0; i + 1 < nsiblings; i += 2 )
   {
      treearray[i] = PHeap_mergeheaps(scip, treearray[i], treearray[i + 1]);
   }
   j = i - 2;

   /* if the number of trees is odd, get the last one */
   if( j == nsiblings - 3 )
   {
      treearray[j] = PHeap_mergeheaps(scip, treearray[j], treearray[j + 2]);
   }

   for( ; j >= 2; j -= 2 )
   {
      treearray[j - 2] = PHeap_mergeheaps(scip, treearray[j - 2], treearray[j]);
   }

   q = treearray[0];

   //SCIPfreeBufferArray(scip, &treearray);
   return q;
}

static PHeap* PHeap_insert(SCIP* scip,
   PHeap* p,
   int element,
   SCIP_Real key,
   int size
   )
{
   PHeap *node;
   //SCIP_RETCODE retcode;

   assert(size > 0);

   //retcode = SCIPallocMemory(scip, &node);
   //assert(retcode == SCIP_OKAY);

   node = (PHeap*)malloc(sizeof(PHeap));

   node->size = size;
   node->key = key;
   node->element = element;
   node->child = NULL;
   node->sibling = NULL;
   node->prev = NULL;

   if( p == NULL )
   {
      return node;
   }
   else
   {
      return PHeap_mergeheaps(scip, p, node);
   }

}


static PHeap* PHeap_deleteMin(SCIP* scip, int* element, SCIP_Real *key, PHeap *p)
{
   PHeap *new_root = NULL;

   if( p == NULL )
   {
      *element = -1;
      return NULL;
   }
   else
   {
      *element = p->element;
      *key = p->key;
      /* printf("keyD: %f \n", *key); */
      if( p->child != NULL )
         new_root = PHeap_combine_siblings(scip, p->child);

      //SCIPfreeMemory(scip, &p);
      //p = NULL;
      free(p);
   }
   return new_root;
}


/* initialize the union find structure uf with 'length' many components (of size one) */
static SCIP_RETCODE UF_init(
   SCIP* scip,
   UF* uf,
   int length
   )
{
   int i;
   uf->count = length;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(uf->parent), length) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(uf->size), length) );
   for( i = 0; i < length; i++ ) {
      uf->parent[i] = i;
      uf->size[i] = 1;
   }

   return SCIP_OKAY;
}


/* Find the component identifier */
int UF_find(
   UF* uf,
   int element
   )
{
   int newelement;
   int root = element;
   int* parent = uf->parent;

   while( root != parent[root] )
   {
      root = parent[root];
   }

   while( element != root )
   {
      newelement = parent[element];
      parent[element] = root;
      element = newelement;
   }
   return root;
}

/* merge the components containing p and q respectively */
static void UF_union(
   UF* uf,
   int p,
   int q
   )
{
   int idp;
   int idq;
   int* size = uf->size;
   int* parent = uf->parent;
   idp = UF_find(uf, p);
   idq = UF_find(uf, q);

   /* if p and q lie in the same component, there is nothing to be done */
   if( idp == idq )
      return;

   parent[idq] = idp;
   size[idp] += size[idq];
   /*  TODO FLAG: COMPRESS
       if( size[idp] < size[idq] )
       {
       parent[idp] = idq;
       size[idq] += size[idp];
       }
       else
       {
       parent[idq] = idp;
       size[idp] += size[idq];
       }
   */
   /* one less component */
   uf->count--;

}

static void UF_free(
   SCIP* scip,
   UF* uf
   )
{
   SCIPfreeMemoryArray(scip, &uf->parent);
   SCIPfreeMemoryArray(scip, &uf->size);
   uf = NULL;

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

#if 0
	 int random;
         int tmp;

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



   char* steinertree;
   int i;
   NODE* nodes;
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &steinertree, nnodes) );
   /*** VERTEX  INSERTION ***/
   if( !graph->rootisfixed )
   {

      int newnode = 0;

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
      nodes[graph->source[0]].edge = -1;

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

                  // SCIP_CALL( IDX_add(scip, cand_start, cand_curr, oedge) );

                  SCIP_CALL( SCIPallocMemory(scip, &cand_curr) );

                  // cand_curr = (IDX *)malloc(sizeof(IDX)); //TODO extra method and try to work without cand list
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

                  /* next vertex in the current Steiner tree adjacent to vertex i resp. v (the one being scrutinized for possible insertion) */
                  w = &nodes[graph->head[cand_curr->index]];


                  /* if there is an edge with cost greater than that of the current edge... */
                  max = findMax(scip, graph, w);
                  if( SCIPisGT(scip, graph->cost[max->edge], graph->cost[cand_curr->index]) )
                  {
                     diff += graph->cost[cand_curr->index];
                     diff -= graph->cost[max->edge];

                     //cut_curr = (IDX *)malloc(sizeof(IDX));
                     SCIP_CALL( SCIPallocMemory(scip, &cut_curr) );

                     cut_curr->index = max->edge;
                     cut_curr->parent  = cut_start;
                     cut_start = cut_curr;
                     //SCIP_CALL( IDX_add(scip, cut_start, max->edge) );


                     cut(max);
                     link(v, w, cand_curr->index);
                     SCIP_CALL( SCIPallocMemory(scip, &add_curr) );

                     //    add_curr = (IDX *)malloc(sizeof(IDX));
                     add_curr->index = v->edge;
                     add_curr->parent  = add_start;
                     add_start = add_curr;
                     // SCIP_CALL( IDX_add(scip, add_start, v->edge) );

                  }
                  cand_curr = cand_curr->parent;
               }

               add_curr = add_start;

               /* if the new tree is more expensive than the old one, restore the latter */
               if( !SCIPisNegative(scip, diff) )
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

                  /* finally, cut the edge added first (if it had been cut during the insertion process, it will have been restored above) */
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


               // IDX_freeList(scip, cut_start);
               add_curr = add_start;
               while( add_curr != NULL ){
		  add_start = add_curr->parent;
                  SCIPfreeMemory(scip, &add_curr);
                  //free(add_curr);
                  add_curr = add_start;
               }
               cut_curr = cut_start;
               while( cut_curr != NULL ){
		  cut_start = cut_curr->parent;
                  SCIPfreeMemory(scip, &cut_curr);
                  //free(cut_curr);
                  cut_curr = cut_start;
               }


            }
            cand_curr = cand_start;
            while( cand_curr != NULL ){
               cand_start = cand_curr->parent;
               SCIPfreeMemory(scip, &cand_curr);

               cand_curr = cand_start;
            }

            // IDX_freeList(scip, cand_start);

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

      obj = 0.0;
      for( e = 0; e < graph->edges; e++)
         obj += (best_result[e] > -1) ? graph->cost[e] : 0.0;

      printf(" ObjAfterVertexInsertion=%.12e\n", obj);
   }





   /*** Key Path Exchange***/

#if 0

   PHeap* x;
   PHeap *y;
   PHeap *z;
   int el;
   z = PHeap_insert(scip, NULL, 0, 0.1, 10);
   z = PHeap_insert(scip, z, 0, 6.5, 10) ;
   x = PHeap_insert(scip, NULL, 0, 3.1, 10) ;
   x = PHeap_insert(scip, x, 0, 1.5, 10) ;
   x = PHeap_insert(scip, x, 0, 3.5, 10) ;

   y = PHeap_insert(scip, NULL, 0, 1.6, 10)  ;
   y = PHeap_insert(scip, y, 0, 0.19, 10) ;


   SCIP_Real key = -1;
   x = PHeap_mergeheaps(scip, y,x);

   x = PHeap_deleteMin(scip, &el, &key, x);

   printf(   "key 0.19 ==%f \n", key);

   x = PHeap_deleteMin(scip, &el,&key, x);

   printf(   "key 1.5 ==%f \n", key);
   x = PHeap_mergeheaps(scip, x,z);
   x = PHeap_deleteMin(scip, &el,&key, x);

   printf(   "key==0.1 %f \n", key);


   x = PHeap_deleteMin(scip, &el, &key, x);

   printf(   "key ==1.6 %f \n", key);

#endif

   if( !graph->rootisfixed )
   {
      /* array [1,..,nnodes],
       * if node i is in the current ST, blists_start[i] points to a linked list of all nodes having i as their base */
      IDX* newpath_start;
      IDX* newpath_curr;
      IDX* blists_curr;
      IDX** blists_start;
      IDX* kpath_curr;
      IDX* kpath_start;
      PHeap** heaparr;
      UF uf;  /* union find*/
      int* memvbase;
      int* meminedges;
      SCIP_Real* memdist;
      int nresnodes;
      int* vbase;     /* array [1,...,nnodes] */
      int* state;
      int node;
      int kptailnode;  /* tail node of the current keypath*/
      int crucnode;   /* current crucial node*/
      int keynode = - 1; /* first node of the key-path to be replaced */
      int newedge;
      int oldedge;
      int adjnode;
      int edge;
      int count = 0;
      SCIP_Real bestdiff = 0;
      int newpathstart = -1;
      SCIP_Real kpathcost;
      SCIP_Real edgecost;
      PATH* vnoi;
      char dummy = FALSE;


      SCIP_CALL( SCIPallocBufferArray(scip, &blists_start, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &memvbase, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &memdist, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &meminedges, nnodes) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &heaparr, nnodes) );

      /* init data structures */
      UF_init(scip, &uf, nnodes);
      BMSclearMemoryArray(heaparr, nnodes);

      /* DFS tree TODO!! STACK based iteration! */
      SCIP_CALL( SCIPallocBufferArray(scip, &dfstree, nnodes) );
      crucList(graph, best_result, &(graph->source[0]), &count);

      /* compute a voronoi diagram with the current ST nodes as bases */
      voronoi(graph, graph->cost, steinertree, vbase, vnoi);
      state = graph->path_state;
      // for(  e = 0; e < nnodes; e++)
      //   printf("(completevoronoi)base[%d] = %d \n", e, vbase[e]);

      for(  e = 0; e < graph->edges; e += 2 )
      {
         /* check if e is a boundary edge */
         if( vbase[graph->tail[e]] != vbase[graph->head[e]] )
         {
            //printf("!= .. %d_%d \n ", graph->tail[e], graph->head[e] );
            edgecost = vnoi[graph->tail[e]].dist + graph->cost[e] + vnoi[graph->head[e]].dist;
            //printf("put in pairheap: %d_%d cost : %f bases: %d %d \n ", graph->tail[e], graph->head[e], edgecost, vbase[graph->tail[e]], vbase[graph->head[e]] );
            heaparr[vbase[graph->tail[e]]] = PHeap_insert( scip, heaparr[vbase[graph->tail[e]]], e, edgecost, graph->edges);
            heaparr[vbase[graph->head[e]]] = PHeap_insert( scip, heaparr[vbase[graph->head[e]]], e, edgecost, graph->edges);
         }
      }

      /* init the baselists */
      BMSclearMemoryArray(blists_start, nnodes);

      for( k = 0; k < nnodes; k++ )
      {
         assert( state[k] == CONNECT);
         if( steinertree[k] )
         {
            if( 0 || dummy )
               printf("node[%d].parent = %d \n", k, graph->head[nodes[k].edge] );

         }
         SCIP_CALL( SCIPallocMemory(scip, &blists_curr) ); /*TODO extra method  (inline) */
         blists_curr->index = k;
         blists_curr->parent = blists_start[vbase[k]];
         blists_start[vbase[k]] = blists_curr;

      }


      newpath_start = NULL;

      /* main loop visiting all crucial (i.e. terminal- and key-) vertices in post-order */
      for( i = 0; dfstree[i] != graph->source[0]; i++ )
      {
         dummy = FALSE;
         crucnode = dfstree[i];
         if( dummy )
            printf("INNXX %d\n", crucnode);
         /* is crucnode a crucial vertex? */
         if( !nodeIsCrucial(graph, steinertree, crucnode) )
         {
            continue;
         }


         for( k = 0; k < nnodes; k++ )
         {
            assert( state[k] == CONNECT);
         }
         /* find the (unique) key-path containing the parent of the current crucial node crucnode */
         kpath_start = NULL;
         kptailnode = graph->head[nodes[crucnode].edge];
         kpathcost = graph->cost[nodes[crucnode].edge];

         while( !nodeIsCrucial(graph, steinertree, kptailnode) )
         {

            /*dummy = TRUE; */
            kpathcost += graph->cost[nodes[kptailnode].edge];
            SCIP_CALL( SCIPallocMemory(scip, &kpath_curr) );
            kpath_curr->index = kptailnode;
            if( dummy )
            {
               printf("kpathintern: %d \n " , kptailnode);
            }
            kpath_curr->parent = kpath_start;
            kpath_start = kpath_curr;
            kptailnode = graph->head[nodes[kptailnode].edge];
         }
         if( dummy )
         {
            printf("kpathstart: %d \n " ,crucnode);
            printf("kpathend: %d \n " , kptailnode);
         }

         nresnodes = 0;


         /* reset all nodes (referred to as 'C' henceforth) whose bases are internal nodes of the current keypath */
         kpath_curr = kpath_start;
         while( kpath_curr != NULL )
         {
            /* reset all nodes having the current (internal) keypath node as their voronoi base */
            blists_curr = blists_start[kpath_curr->index];
            assert( blists_curr != NULL ); //TODO reset later!!
            while( blists_curr != NULL )
            {
               node = blists_curr->index;
               if( dummy )
                  printf("C-node %d \n", blists_curr->index);
               memvbase[nresnodes] = vbase[node];
               memdist[nresnodes] =  vnoi[node].dist;
               meminedges[nresnodes] = vnoi[node].edge;
               nresnodes++;
               vbase[node] = UNKNOWN;
               vnoi[node].dist = FARAWAY;
               vnoi[node].edge = UNKNOWN;
               state[node] = UNKNOWN;
               blists_curr = blists_curr->parent;
            }

            /* update the union find data structure */
            if( kpath_curr != kpath_start )
               UF_union(&uf, kpath_start->index, kpath_curr->index);

            kpath_curr = kpath_curr->parent;

         }

         e = UNKNOWN;
         while( heaparr[crucnode] != NULL )
         {
            heaparr[crucnode] = PHeap_deleteMin(scip, &e, &edgecost, heaparr[crucnode]);

            node = (vbase[graph->head[e]] == UNKNOWN)? UNKNOWN : UF_find(&uf, vbase[graph->head[e]]);
            adjnode = (vbase[graph->tail[e]] == UNKNOWN)? UNKNOWN : UF_find(&uf, vbase[graph->tail[e]]);

            if ( e != -1 && dummy )
            {
               printf("prenodes %d_%d  \n ", graph->head[e], graph->tail[e]);
               printf("basenodes:  %d_%d\n ", node, adjnode);
            }

            if( (( node != crucnode && node!= UNKNOWN ) || ( adjnode != crucnode && adjnode!= UNKNOWN ) ) ) /*&& (graph->head[e] != 9 && graph->tail[e] !=9) */
            {
               if( dummy )
                  printf("break %d %d\n", node, adjnode);
               heaparr[crucnode] = PHeap_insert( scip, heaparr[crucnode], e, edgecost, graph->edges);
               break;

            }

            if( heaparr[crucnode] == NULL )
            {
               e = UNKNOWN;
               break;
            }

         }
         oldedge = e;
         /*  for(  e = 0; e < nnodes; e++)
             printf("(befrepairvoronoi)base[%d] = %d \n", e, vbase[e]); */
         count = 0;
         /* try to connect the nodes of C (directly) to COMP(C), as a preprocessing for voronoi_repair */
         kpath_curr = kpath_start;
         while( kpath_curr != NULL )
         {
            blists_curr = blists_start[kpath_curr->index];
            assert( blists_curr != NULL );
            while( blists_curr != NULL ) /* TODO flip edge implicitly*/
            {
               node = blists_curr->index;
               edge = graph->outbeg[node];
               while( edge >= 0 )
               {
                  adjnode = graph->head[edge];

                  /* check if the adjacent node is not in C and allows a better voronoi assignment of the current node */
                  if( state[adjnode] == CONNECT && SCIPisGT(scip, vnoi[node].dist, vnoi[adjnode].dist + graph->cost[edge]) )
                  {
                     vnoi[node].dist = vnoi[graph->head[edge]].dist + graph->cost[edge];
                     vbase[node] = vbase[adjnode];
                     vnoi[node].edge = flipEdge(edge);

                  }
                  edge = graph->oeat[edge];
               }
               if( vbase[node] != UNKNOWN )
               {
                  if( 0 && SCIPisGT(scip, (newedge == -1) ? FARAWAY : graph->cost[newedge], graph->cost[vnoi[node].edge]) )
                  {
                     newedge = vnoi[node].edge;
                  }
                  heap_add(graph->path_heap, state, &count, node, vnoi);
               }
               blists_curr = blists_curr->parent;
            }
            kpath_curr = kpath_curr->parent;

         }

         newedge = UNKNOWN;

         /* if there is no key path, nothing has to be repaired */
         if( kpath_start != NULL )
         {
            voronoi_repair(scip, graph, graph->cost, &count, vbase, vnoi, &newedge, crucnode, &uf);
         }

         //for(  e = 0; e < nnodes; e++)
         //    printf("(completevoronoi)base[%d] = %d \n", e, vbase[e]);
         if( dummy ){
            printf("newedge  %d_%d\n ", graph->tail[newedge], graph->head[newedge]);

            printf("oldedge   %d_%d\n ", graph->tail[oldedge], graph->head[oldedge]);
            printf("newedge pred  %d_%d  \n ", graph->tail[vnoi[graph->tail[newedge]].edge], graph->tail[vnoi[graph->head[newedge]].edge]);
            printf("oldedge pred  %d_%d  \n ", graph->tail[vnoi[graph->tail[oldedge]].edge], graph->tail[vnoi[graph->head[oldedge]].edge]);

         }
         if( oldedge != UNKNOWN && newedge != UNKNOWN && SCIPisLT(scip, edgecost,
               vnoi[graph->tail[newedge]].dist + graph->cost[newedge] + vnoi[graph->head[newedge]].dist) )
            newedge = oldedge;
         if( oldedge != UNKNOWN && newedge == UNKNOWN )
            newedge = oldedge;
         if( dummy ){
            printf("KOSTENVERGLEICH old/new: %f_ %f \n ", kpathcost,  vnoi[graph->tail[newedge]].dist + graph->cost[newedge]
               + vnoi[graph->head[newedge]].dist );
            printf("final edge %d_%d \n ", graph->tail[newedge], graph->head[newedge]);
         }

         edgecost = vnoi[graph->tail[newedge]].dist + graph->cost[newedge] + vnoi[graph->head[newedge]].dist;
         if( SCIPisLT(scip, edgecost, kpathcost) && SCIPisLT(scip, edgecost - kpathcost, bestdiff) )
         {
            bestdiff = edgecost - kpathcost;
            newpath_curr = newpath_start;
            while( newpath_curr != NULL )
            {
               newpath_start = newpath_curr->parent;
               SCIPfreeMemory(scip, &newpath_curr);
               newpath_curr = newpath_start;
            }
            node = UF_find(&uf, vbase[graph->head[newedge]]);
            // node = UF_find(&uf, vbase[graph->tail[newedge]]);

            newpath_start = NULL;
            keynode = crucnode;

            k = (node == crucnode)? graph->head[newedge] : graph->tail[newedge];
            while( k != vbase[k] )
            {
               SCIP_CALL( SCIPallocMemory(scip, &newpath_curr) );
               newpath_curr->index = flipEdge(vnoi[k].edge);
               //printf("bestpath: %d->%d \n " , graph->tail[newpath_curr->index], graph->head[newpath_curr->index]);
               newpath_curr->parent = newpath_start;
               newpath_start = newpath_curr;
               k = graph->tail[vnoi[k].edge];
            }

            k = (node == crucnode)? graph->tail[newedge] : graph->head[newedge];
            while( k != vbase[k] )
            {
               SCIP_CALL( SCIPallocMemory(scip, &newpath_curr) );
               newpath_curr->index = vnoi[k].edge;
               //printf("bestpath: %d->%d \n " , graph->tail[newpath_curr->index], graph->head[newpath_curr->index]);
               newpath_curr->parent = newpath_start;
               newpath_start = newpath_curr;
               k = graph->tail[vnoi[k].edge];
            }

            //printf("newedge : %d \n" , graph->head[newedge]);
            SCIP_CALL( SCIPallocMemory(scip, &newpath_curr) );
            newpath_curr->index = (node == crucnode)? newedge : flipEdge(newedge);
            //printf("bestpathI: %d->%d \n " , graph->tail[newpath_curr->index], graph->head[newpath_curr->index]);
            newpath_curr->parent = newpath_start;
            newpath_start = newpath_curr;
            newpathstart = (node == crucnode)? vbase[graph->head[newedge]] : vbase[graph->tail[newedge]];


         }

         /* update union find data structure */
         if( kpath_start != NULL )
         {
            UF_union(&uf, crucnode, kpath_start->index);
         }
         UF_union(&uf, kptailnode, crucnode);


         /* restore the old voronoi digram TODO besser implementieren!!!!*/
         k = 0;
         kpath_curr = kpath_start;
         while( kpath_curr != NULL )
         {
            /* reset all nodes having the current (internal) keypath node as their voronoi base */
            blists_curr = blists_start[kpath_curr->index];
            while( blists_curr != NULL )
            {
               node = blists_curr->index;
               vbase[node] = memvbase[k];
               vnoi[node].dist = memdist[k];
               vnoi[node].edge = meminedges[k];
               k++;
               blists_curr = blists_curr->parent;
            }

            if( kpath_curr->parent != NULL )
               heaparr[kpath_curr->parent->index] = PHeap_mergeheaps(scip, heaparr[kpath_curr->parent->index], heaparr[kpath_curr->index]);
            else
               heaparr[kptailnode] = PHeap_mergeheaps(scip, heaparr[kptailnode], heaparr[kpath_curr->index]);



            /* free the current keypath list element */
            kpath_start = kpath_curr->parent;
            SCIPfreeMemory(scip, &kpath_curr);
            //free(kpath_curr);
            kpath_curr = kpath_start;

         }
         heaparr[kptailnode] = PHeap_mergeheaps(scip, heaparr[kptailnode], heaparr[crucnode]);

#if 0
         /* debug! */
         int* xvbase;
         PATH* xvnoi;
         SCIP_CALL( SCIPallocBufferArray(scip, &xvbase, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &xvnoi, nnodes) );

         voronoi(graph, graph->cost, steinertree, xvbase, xvnoi);

         for( e = 0; e < nnodes; e++ )
         {
            assert(vbase[e] == xvbase[e] && vnoi[e].dist == xvnoi[e].dist && vnoi[e].edge == xvnoi[e].edge);
         }
         //for(  e = 0; e < nnodes; e++)
         // printf("(repairedvoronoi)base[%d] = %d \n", e, vbase[e]);

         /* debug! */
         SCIPfreeBufferArray(scip, &xvbase);
         SCIPfreeBufferArray(scip, &xvnoi);
#endif


      }
#if 0
      SCIP_Bool* edgemark;
      SCIP_CALL( SCIPallocBufferArray(scip, &edgemark, graph->edges / 2) );
      for( e = 0; e < graph->edges / 2; e++ ){
         edgemark[e] = FALSE;
      }
#endif
      /* update the solution */
      if( bestdiff != 0 )
      {
         assert(bestdiff < 0 );
         printf( "ADDING NEW KEY PATH (%f )\n", bestdiff );
         /* remove old keypath */
         k = keynode;
         do
         {

            assert(  best_result[flipEdge(nodes[k].edge)] != -1 );
            best_result[flipEdge(nodes[k].edge)] = -1;
            /* steinertree[k] = FALSE;*/
            //printf("delete: : %d->%d \n " , graph->head[nodes[k].edge], graph->tail[nodes[k].edge]);


            /*edgemark[flipEdge(nodes[k].edge) / 2] = TRUE; */


            k = graph->head[nodes[k].edge];

         }
         while( !nodeIsCrucial(graph, steinertree, k) );

         /* add new keypath */
         newpath_curr = newpath_start;
         while( newpath_curr != NULL ){
            newpath_start = newpath_curr->parent;
            assert( best_result[newpath_curr->index] != 0 );
            /* edgemark[newpath_curr->index / 2] = TRUE; */
            best_result[newpath_curr->index] = 0;
            //printf("add: : %d->%d \n " , graph->tail[newpath_curr->index], graph->head[newpath_curr->index]);

            SCIPfreeMemory(scip, &newpath_curr);
            newpath_curr = newpath_start;
         }
         //printf("newpathstart %d \n", newpathstart);
         //printf("crucnode %d \n", crucnode);
         evert(&nodes[newpathstart]); /* TODO cut new path (wrt nodes[]) before for less comp. time */
	 k = keynode;
	 while( k != newpathstart )
	 {
            //printf("delete: %d ->%d\n ",graph->tail[nodes[k].edge], graph->head[nodes[k].edge]);
            //assert( best_result[nodes[k].edge] == 0);
	    best_result[nodes[k].edge] = -1;
	    //assert( best_result[nodes[k].edge] == -1);
	    best_result[flipEdge(nodes[k].edge)] = 0;
	    k = graph->head[nodes[k].edge];
	 }
      }




      /* SCIP_CALL( SCIPprobdataPrintGraph2(graph,"TESTG.gml", edgemark) );
         SCIPfreeBufferArray(scip, &edgemark);*/

      SCIPfreeBufferArray(scip, &dfstree);
      /**********************************************************/

      /* free data structures */
      UF_free(scip, &uf);


      SCIPfreeBufferArray(scip, &vbase);
      SCIPfreeBufferArray(scip, &vnoi);
      SCIPfreeBufferArray(scip, &memvbase);
      SCIPfreeBufferArray(scip, &memdist);
      SCIPfreeBufferArray(scip, &meminedges);


      for( k = 0; k < nnodes; k++ )
      {
         //SCIPfreeMemoryNull(scip, &heaparr[k]);
         blists_curr = blists_start[k];
         while( blists_curr != NULL ){
            blists_start[k] = blists_curr->parent;
            SCIPfreeMemory(scip, &blists_curr);
            blists_curr = blists_start[k];
         }


      }
      SCIPfreeMemoryArray(scip, &heaparr);
      SCIPfreeBufferArray(scip, &blists_start);

   }

   for( k = 0; k < nnodes; k++ )
   {
      assert(path[k] == NULL || graph->term[k] == layer);
      SCIPfreeBufferArrayNull(scip, &(path[k]));

   }

   SCIPfreeBufferArray(scip, &nodes);
   SCIPfreeBufferArray(scip, &steinertree);

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
