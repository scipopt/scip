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
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "heur_tm.h"
#include "probdata_stp.h"
#include "grph.h"
#include "portab.h"
#include "scip/misc.h"
#include <time.h>
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
#define VEL 1
#define AUTO 0
#define TM 1
#define TMPOLZIN 2

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


typedef struct PairingHeap_Node
{
   int element;
   SCIP_Real key;
   struct PairingHeap_Node* child;
   struct PairingHeap_Node* sibling;
   struct PairingHeap_Node* prev;

}PHNODE;

typedef struct PairingHeap
{
   int size;
   PHNODE* root;
}PAIRHEAP;

#if 0
static
int sorta(void *first_arg, void *second_arg )
{
   int first = *(int*)first_arg;
   int second = *(int*)second_arg;
   if ( first < second )
   {
      return -1;
   }
   else if ( first == second )
   {
      return 0;
   }
   else
   {
      return 1;
   }
}
#endif
int GNODECmpByDist(void *first_arg, void *second_arg )
{
   SCIP_Real first = ((GNODE*)first_arg)->dist;
   SCIP_Real second = ((GNODE*)second_arg)->dist;
   if( first < second )
   {
      return -1;
   }
   else if( first == second )
   {
      return 0;
   }
   else
   {
      return 1;
   }
}
#if 0

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



#if 0


static
void IDX_freeList(
   SCIP* scip,
   IDX** start
   )
{
   IDX* curr = *start;
   if( (*start) == NULL )
      return;

   while( (curr) != NULL ){
      *start = (curr)->parent;
      SCIPfreeMemory(scip, &curr);
      curr = *start;
   }

}


static
SCIP_RETCODE IDX_add(
   SCIP* scip,
   IDX** start,
   IDX** curr,
   int index
   )
{
   SCIP_CALL( SCIPallocMemory(scip, &curr) );
   (*curr)->index = index;
   (*curr)->parent  =  *start;
   (*start) = (*curr);
   return SCIP_OKAY;

}


static
void crucList(
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

static
void phnode_insert2(SCIP* scip,
   PHNODE* p,
   int element,
   SCIP_Real key,
   int size
   )
{
   PHNODE* node;

   assert(size > 0);

   node = (PHNODE*)malloc(sizeof(PHNODE));

   node->size = size;
   node->key = key;
   node->element = element;
   node->child = NULL;
   node->sibling = NULL;
   node->prev = NULL;

   if( p == NULL )
   {
      printf("return \n");



      p = (PHNODE*)malloc(sizeof(PHNODE));

      p->size = size;
      p->key = key;
      p->element = element;
      p->child = NULL;
      p->sibling = NULL;
      p->prev = NULL;

      printf("p->key %f return \n", p->key);
   }
   else
   {
      printf("merge \n");
      p = phnode_mergeheaps(scip, p, node);
   }

}



/** returns the opposite edge */
static
int flipedge(
   int edge
   )
{
   assert(edge >= 0);
   return ((edge % 2) == 0) ? edge + 1 : edge - 1;
}


/** reverse the edge orientation in the Steiner Tree between two nodes */
static
void revorient(
   const GRAPH* graph,
   NODE* stnodes,
   int* stedges,
   int start,
   int end
   )
{
   int e;
   int node = start;

   while( node != end && node != graph->source[0])
   {
      /* the ST edge pointing towards the root */
      e = stnodes[node].edge;
      assert(stedges[e] == -1 && stedges[flipedge(e)] != -1 );
      printf(" switch : %d->%d \n ", graph->tail[e], graph->head[e]);
      stedges[e] = 0;
      stedges[flipedge(e)] = -1;
      node = graph->head[e];
   }
}

#endif

/** recursive methode for a DFS ordering of graph 'g' */
static
void dfsorder(
   const GRAPH* graph,
   int* edges,
   int* node,
   int* counter,
   int* dfst
   )
{
   int oedge = graph->outbeg[*node];

   while( oedge >= 0 )
   {
      if( edges[oedge] >= 0 )
      {
         dfsorder(graph, edges, &(graph->head[oedge]), counter, dfst);
      }
      oedge = graph->oeat[oedge];
   }

   dfst[*counter] = *node;
   (*counter)++;
}


/** checks whether node is crucial, i.e. a terminal or a vertex with degree at least 3 (w.r.t. the steinertree) */
static
char nodeIsCrucial(
   const GRAPH* graph,
   int* steineredges,
   int node
   )
{
   assert(graph != NULL);
   assert(steineredges != NULL);

   if( graph->term[node] == -1 )
   {
      int counter = 0;
      int e = graph->outbeg[node];
      while( e >= 0 )
      {

         /* check if the adjacent node is in the ST */
         if( steineredges[e] > -1 || steineredges[flipedge(e)] > -1 )
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

/** node degree (w.r.t. the steinertree) */
static
int stdeg(
   const GRAPH* graph,
   int* steineredges,
   int node
   )
{

   int counter = 0;
   int e;

   assert(graph != NULL);
   assert(steineredges != NULL);

   e = graph->outbeg[node];

   while( e >= 0 )
   {
      /* check whether edge 'e' is in the ST */
      if( steineredges[e] > -1 || steineredges[flipedge(e)] > -1 )
      {
         counter++;
      }
      e = graph->oeat[e];
   }

   return counter;
}


/** Linear Link Cut Tree */

/** inits a node, setting 'parent' and 'edge' to its default values */
static
void init(
   NODE* v
   )
{
   v->parent = NULL;
   v->edge = -1;
}

/** renders w a child of v; v has to be the root of its tree */
static
void link(
   NODE* v,
   NODE* w,
   int edge
   )
{
   assert(v->parent == NULL);
   v->parent = w;
   v->edge = edge;
}

static
void cut(
   NODE* v
   )
{
   v->edge = -1;
   v->parent = NULL;
}

/** finds the max value between node 'v' and the root of the tree **/
static
NODE* findMax(
   SCIP* scip,
   const SCIP_Real* cost,
   NODE* v
   )
{
   NODE* p = v;
   NODE* q = NULL;
   SCIP_Real max = -1;
   while( p != NULL )
   {
      if( SCIPisGE(scip, cost[p->edge], max) )
      {
         max = cost[p->edge];
         q = p;
      }
      p = p->parent;
   }
   return q;
}

/** makes vertex v the root of the link cut tree */
static
void evert(
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
      q->edge = (val >= 0) ? flipedge(val) : val;
      val = tmpval;
      q->parent = p;
      p = q;
      q = r;
   }
}



/*** Pairing Heap ****/

/** links nodes 'root1' and 'root2' together */
static
PHNODE* phnode_mergeheaps(
   SCIP* scip,
   PHNODE *root1,
   PHNODE *root2
   )
{
   if( root2 == NULL )
      return root1;
   if( root1 == NULL )
      return root2;

   if( root1->key <= root2->key )
   {
      /* attach root2 as (the leftmost) child of root1 */
      root2->prev = root1;
      root1->sibling = root2->sibling;
      if( root1->sibling != NULL )
         root1->sibling->prev = root1;

      root2->sibling = root1->child;
      if( root2->sibling != NULL )
         root2->sibling->prev = root2;

      root1->child = root2;

      return root1;
   }
   else
   {
      /* attach root1 as (the leftmost) child of root2 */
      root2->prev = root1->prev;
      root1->prev = root2;
      root1->sibling = root2->child;
      if( root1->sibling != NULL )
         root1->sibling->prev = root1;

      root2->child = root1;

      return root2;
   }
}

static
PHNODE* phnode_addtoheap(
   SCIP* scip,
   PHNODE *root1,
   PHNODE *root2
   )
{
   assert(root2 != NULL);
   assert(root1 != NULL);

   if( root1->key <= root2->key )
   {
      /* attach root2 as (the leftmost) child of root1 */
      root2->prev = root1;
      root1->sibling = root2->sibling;
      if( root1->sibling != NULL )
      {
         root1->sibling->prev = root1;
	 assert(0);
      }

      root2->sibling = root1->child;
      if( root2->sibling != NULL )
         root2->sibling->prev = root2;

      root1->child = root2;

      return root1;
   }
   else
   {
      /* attach root1 as (the leftmost) child of root2 */
      root2->prev = root1->prev;
      root1->prev = root2;
      root1->sibling = root2->child;
      if( root1->sibling != NULL )
         root1->sibling->prev = root1;

      root2->child = root1;

      return root2;
   }
}

/** internal method for combining the siblings after the root has been deleted */
static
SCIP_RETCODE phnode_combine_siblings(
   SCIP*                 scip,              /**< SCIP data structure */
   PHNODE**              p,                  /**< the first sibling */
   int                   size                 /**< the size of the heap */
   )
{
   PHNODE** treearray;
   int i;
   int j;
   int nsiblings;
   if( (*p)->sibling == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &treearray, size) );

   /* store all siblings in an array */
   for( nsiblings = 0; (*p) != NULL; nsiblings++ )
   {
#if 0
      /* if the last entry is reached, double the size of the array */
      if( nsiblings == *size )
      {
         SCIP_CALL( phnode_double(scip, &treearray, *size) );
         *size = *size * 2;
      }
#endif
      treearray[nsiblings] = (*p);
      if( (*p)->prev != NULL )
         (*p)->prev->sibling = NULL;
      (*p) = (*p)->sibling;
   }
   treearray[nsiblings] = NULL;

#if 0
   /*combine the subtrees (simple) */
   for(i = 1; i < nsiblings; i++)
      treearray[i] = phnode_mergeheaps(scip, treearray[i-1], treearray[i]);


   return treearray[nsiblings-1];
#endif

   /* combine the subtrees (two at a time) */
   for( i = 0; i + 1 < nsiblings; i += 2 )
   {
      treearray[i] = phnode_mergeheaps(scip, treearray[i], treearray[i + 1]);
   }
   j = i - 2;

   /* if the number of trees is odd, get the last one */
   if( j == nsiblings - 3 )
   {
      treearray[j] = phnode_mergeheaps(scip, treearray[j], treearray[j + 2]);
   }

   for( ; j >= 2; j -= 2 )
   {
      treearray[j - 2] = phnode_mergeheaps(scip, treearray[j - 2], treearray[j]);
   }

   (*p) = treearray[0];

   SCIPfreeBufferArray(scip, &treearray);

   return SCIP_OKAY;
}


/** inserts a new node into the pairing heap */
static
SCIP_RETCODE phnode_insert(
   SCIP* scip,
   PHNODE** root,
   int element,
   SCIP_Real key,
   int* size
   )
{
   if( (*root) == NULL )
   {
      (*size) = 1;
      SCIP_CALL( SCIPallocBuffer(scip, root) );
      (*root)->key = key;
      (*root)->element = element;
      (*root)->child = NULL;
      (*root)->sibling = NULL;
      (*root)->prev = NULL;
   }
   else
   {
      PHNODE* node;
      (*size)++;
      SCIP_CALL( SCIPallocBuffer(scip, &node) );
      node->key = key;
      node->element = element;
      node->child = NULL;
      node->sibling = NULL;
      node->prev = NULL;
      (*root) = phnode_addtoheap(scip, (*root), node);
   }
   return SCIP_OKAY;
}

/** deletes the root of the paring heap, concomitantly storing its data and key in '*element' and '*key' respectively */
static
void phnode_deletemin(
   SCIP* scip,
   int* element,
   SCIP_Real *key,
   PHNODE** root,
   int* size
   )
{

   if( (*root) == NULL )
   {
      *element = -1;
      return;
   }
   else
   {
      PHNODE *newroot = NULL;
      *element = (*root)->element;
      *key = (*root)->key;
      if( (*root)->child != NULL )
      {
	 newroot = (*root)->child;
         phnode_combine_siblings(scip, &newroot, --(*size));
      }

      SCIPfreeBuffer(scip, root);
      (*root) = newroot;
   }
}


/** links nodes 'root1' and 'root2' together, roots the resulting tree at root1 and sets root2 to NULL */
static
void phnode_meldheaps(
   SCIP* scip,
   PHNODE** root1,
   PHNODE** root2,
   int* sizeroot1,
   int* sizeroot2
   )
{
   if( root1 == NULL && root2 == NULL )
   {
      assert(sizeroot1 == 0 && sizeroot2 == 0);
      return;
   }

   (*root1) = phnode_mergeheaps(scip, *root1, *root2);
   (*sizeroot1) += (*sizeroot2);
   (*root2) = NULL;
}

#if 0


/** SCIP_Reals the array 'arr', maintaining its original values */
static
SCIP_RETCODE phnode_double(
   SCIP* scip,
   PHNODE*** arr,
   int length
   )
{
   int i;
   PHNODE** newarr;
   SCIP_CALL( SCIPallocBufferArray(scip, &newarr, length * 2) );
   printf("doublingARR 2\n");

   for( i = 0; i < length; i++ )
   {
      newarr[i] = (*arr)[i];
   }
   printf("doublingARR3 \n");
   SCIPfreeBufferArray(scip, arr);
   arr = &newarr;
   return SCIP_OKAY;
}

/** links heaps 'p' and 'q' together */
static
PAIRHEAP* pairheap_mergeheaps(
   SCIP* scip,
   PAIRHEAP** p,
   PAIRHEAP** q
   )
{
   PHNODE* rootp;
   PHNODE* rootq;
   PAIRHEAP** newheap;
   PHNODE* newroot;

   assert( *p != NULL && *q != NULL);
   rootp = (*p)->root;
   rootq = (*q)->root;
   (*newheap)->size = (*p)->size + (*q)->size;
   if( rootp->key <= rootq->key )
   {
      newheap = p;
      SCIPfreeBuffer(scip, q);
      q = NULL;
   }
   else
   {
      newheap = q;
      SCIPfreeBuffer(scip, p);
      p = NULL;
   }
   newroot = phnode_mergeheaps(scip, rootp, rootq);
   (*newheap)->root = newroot;
   return *newheap;
}

static
SCIP_RETCODE pairheap_insert(
   SCIP* scip,
   PAIRHEAP** p,
   int element,
   SCIP_Real key
   )
{
   if( (*p) == NULL )
   {
      SCIP_CALL( SCIPallocBuffer(scip, p) );
      (*p)->root = NULL;
   }
   SCIP_CALL( phnode_insert(scip, &((*p)->root), element, key, &((*p)->size)) );
   return SCIP_OKAY;
}


/** deletes the root of the paring heap, concomitantly storing its data and key in '*element' and '*key' respectively */
static
void pairheap_deletemin(
   SCIP* scip,
   int* element,
   SCIP_Real *key,
   PAIRHEAP** p
   )
{
   if( (*p) == NULL )
   {
      *element = -1;
   }
   else if( (*p)->root == NULL )
   {
      (*p) = NULL;
      *element = -1;
   }
   else
   {
      PHNODE** root = &((*p)->root);
      PHNODE* newroot = NULL;
      *element = (*root)->element;
      *key = (*root)->key;
      if( (*root)->child != NULL )
      {
	 newroot = (*root)->child;
         phnode_combine_siblings(scip, &newroot, (*p)->size);
      }

      SCIPfreeBuffer(scip, root);

      if( newroot == NULL )
	 (*p) = NULL;
      else
         (*p)->root = newroot;
   }
}

/** stores all elements of the pairing heap in an array */
static
SCIP_RETCODE pairheap_buffarr(SCIP* scip, PAIRHEAP* p, int** elements)
{
   int n = 0;
   PHNODE* root = p->root;
   SCIP_CALL( SCIPallocBufferArray(scip, elements, p->size) );
   phnode_rec(root, elements, &n);
   return SCIP_OKAY;
}


/** frees the paring heap */
static
void pairheap_free(SCIP* scip, PAIRHEAP** p)
{
   PHNODE* root = (*p)->root;
   phnode_free(scip, &root);
   SCIPfreeBuffer(scip, p);
}

#endif

/** frees the paring heap with root 'p' */
static
void phnode_free(SCIP* scip, PHNODE** root)
{
   if( (*root) == NULL )
   {
      return;
   }
   if( (*root)->sibling != NULL )
   {
      phnode_free(scip, &((*root)->sibling));
   }
   if( (*root)->child != NULL )
   {
      phnode_free(scip, &((*root)->child));
   }

   SCIPfreeBuffer(scip, root);
   (*root) = NULL;

}


/** internal method used by 'pairheap_buffarr' */
static
void phnode_rec(PHNODE* p, int** arr, int* n)
{
   if( p == NULL )
   {
      return;
   }
   (*arr)[(*n)++] = p->element;
   phnode_rec(p->sibling, arr, n);
   phnode_rec(p->child, arr, n);
}


/** stores all elements of the pairing heap in an array */
static
SCIP_RETCODE phnode_buffarr(SCIP* scip, PHNODE* root, int size, int** elements)
{
   int n = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, elements, size) );
   phnode_rec(root, elements, &n);
   return SCIP_OKAY;
}



/*** Union-Find data structure ***/

/** initializes the union-find structure 'uf' with 'length' many components (of size one) */
static
SCIP_RETCODE uf_init(
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


/** finds and returns the component identifier */
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

/** merges the components containing p and q respectively */
static
void uf_union(
   UF* uf,
   int p,
   int q,
   SCIP_Bool compress
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

   if( !compress )
   {
      parent[idq] = idp;
      size[idp] += size[idq];
   }
   else
   {
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
   }

   /* one less component */
   uf->count--;

}

/** frees the data fields of the union-find structure */
static
void uf_free(
   SCIP* scip,
   UF* uf
   )
{
   SCIPfreeMemoryArray(scip, &uf->parent);
   SCIPfreeMemoryArray(scip, &uf->size);
   uf = NULL;
}

/** computes lowest common ancestors for all pairs {vbase(v), vbase(u)} such that {u,w} is a boundary edge,
 * first call should be with u := root TODO only steineredges needed? */
static
SCIP_RETCODE lca(
   SCIP* scip,
   const GRAPH* graph,
   int u,
   UF* uf, char* nodesmark,
   int* steineredges,
   IDX** lcalists,
   PHNODE** boundpaths,
   int* heapsize,
   int* vbase
   )
{
   int* uboundpaths; /* boundary-paths (each one represented by its boundary edge) having node 'u' as an endpoint */
   int ancestor;
   int v;
   int i;
   int oedge; /* outgoing edge */
   IDX* curr;
   uf->parent[u] = u;

   for( oedge = graph->outbeg[u]; oedge != EAT_LAST; oedge = graph->oeat[oedge] )
   {
      v = graph->head[oedge];
      if( steineredges[oedge] == 0 )
      {
         SCIP_CALL( lca(scip, graph, v, uf, nodesmark, steineredges, lcalists, boundpaths, heapsize, vbase) );
         uf_union(uf, u, v, FALSE);
         uf->parent[UF_find(uf, u)] = u;
      }
   }
   nodesmark[u] = TRUE;

   /* iterate through all boundary-paths having one endpoint in the voronoi region of node 'u' */
   phnode_buffarr(scip, boundpaths[u], heapsize[u], &uboundpaths);
   for( i = 0; i < heapsize[u]; i++ )
   {
      oedge = uboundpaths[i];
      v = vbase[graph->head[oedge]];
      if( nodesmark[v] )
      {
	 ancestor = uf->parent[UF_find(uf, v)];

         /* if the ancestor of 'u' and 'v' is one of the two, the boundary-edge is already in boundpaths[u] */
         if( ancestor != u && ancestor != v)
         {
            SCIP_CALL( SCIPallocMemory(scip, &curr) );
            curr->index = oedge;
            curr->parent = lcalists[ancestor];
            lcalists[ancestor] = curr;
         }
      }
   }

   /* free the boundary-paths array */
   SCIPfreeBufferArray(scip, &uboundpaths);

   return SCIP_OKAY;
}



/*
  int CmpByDist(void *arg1, void *arg2 )
  {
  int node1 = *((int*)arg1);
  int node2 = *((int*)arg2);
  SCIP_Real dist1 = distarr[node1][vcount[node1]];
  SCIP_Real dist2 = distarr[node2][vcount[node2]];
  if( dist1 < dist2 )
  {
  return -1;
  }
  else if( dist1 == dist2 )
  {
  return 0;
  }
  else
  {
  return 1;
  }
  }

*/
/** for debug purposes only */
static
SCIP_RETCODE printGraph(
   SCIP* scip,
   const GRAPH*          graph,              /**< Graph to be printed */
   const char*           filename,           /**< Name of the output file */
   int*                  result
   )
{
   char label[SCIP_MAXSTRLEN];
   FILE* file;
   int e;
   int n;
   int m;
   char* stnodes;
   SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, graph->knots ) );

   assert(graph != NULL);
   file = fopen((filename != NULL) ? filename : "graphX.gml", "w");

   for( e = 0; e < graph->knots; e++ )
   {
      stnodes[e] = FALSE;
   }
   for( e = 0; e < graph->edges; e++ )
   {
      if( result[e] == CONNECT )
      {
	 stnodes[graph->tail[e]] = TRUE;
	 stnodes[graph->head[e]] = TRUE;
      }
   }

   /* write GML format opening, undirected */
   SCIPgmlWriteOpening(file, FALSE);

   /* write all nodes, discriminate between root, terminals and the other nodes */
   e = 0;
   m = 0;
   for( n = 0; n < graph->knots; ++n )
   {
      if( stnodes[n] )
      {
         if( n == graph->source[0] )
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Root", n);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "rectangle", "#666666", NULL);
            m = 1;
         }
         else if( graph->term[n] == 0 )
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Terminal %d", n, e + 1);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#ff0000", NULL);
            e += 1;
         }
         else
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Node %d", n, n + 1 - e - m);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#336699", NULL);
         }

      }
   }

   /* write all edges (undirected) */
   for( e = 0; e < graph->edges; e ++ )
   {
      if( result[e] == CONNECT )
      {
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%8.2f", graph->cost[e]);

	 SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, "#ff0000");
      }
   }
   SCIPfreeBufferArray(scip, &stnodes);
   /* write GML format closing */
   SCIPgmlWriteClosing(file);

   return SCIP_OKAY;
}

#if 0
static
SCIP_RETCODE do_heuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_PQUEUE*          pqueue,
   PATH**                path,
   GNODE**               gnodearr,
   SCIP_Real**           distarr,
   const SCIP_Real*      cost,
   const SCIP_Real*      costrev,
   int                   layer,
   int                   start,
   int*                  result,
   int*                  vcount,
   int*                  nodesid,
   int*                  nodenterms,
   int**                 basearr,
   int**                 edgearr,
   char                  firstrun,
   char                  lastrun,
   char*                 connected
   )
{
   PATH*  mst;
   int*   cluster;
   int    csize = 0;
   int    k;
   int    e;
   int    count;
   int    i;
   int    j;
   int    old;
   int    newval;
   int    best;
   int    term;
   int nnodes;
   int nterms;
   int* state;
   int* vbase;
   char* termsmark;
   char* visited;
   SCIP_Real* vcost;
   SCIP_Real min;
   int* reachednodes;
   int nreachednodes;
   int* tovisit;
   int ntovisit;
   int oedge;
   int* terms;
   int nneighbterms;
   int nneighbnodes;

   PATH* vnoi;
   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   assert(path      != NULL);
   assert(layer >= 0 && layer < g->layers);
   nnodes = g->knots;

   SCIPdebugMessage("TM Heuristic: Start=%5d ", start);

   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, nnodes) );
   nterms = g->terms;

   cluster[csize++] = start;

#if TM
   /* PATH*  path1;
      SCIP_CALL( SCIPallocBufferArray(scip, &path1, nnodes) );
      graph_path_exec2(g, FSP_MODE, start, cost, path1, connected, cluster, &csize);
      SCIPfreeBufferArray(scip, &path1); */

   if( firstrun )
   {
      /* PHASE I: */
      for( i = 0; i < nnodes; i++ )
      {
         g->mark[i] = (g->grad[i] > 0);
      }

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
	 nodesid[i] = i;
         visited[i] = FALSE;
         if( g->term[i] != -1 )
         {
            termsmark[i] = TRUE;
            terms[j++] = i;
         }
         else
         {
            termsmark[i] = FALSE;
         }
      }
      assert(j == nterms);

      voronoi(g, cost, costrev, termsmark, vbase, vnoi); /* TODO: chg s.t. edgecosts for base[root] are reversed */
      state = g->path_state;
      for( k = 0; k < nnodes; k++ )
      {
	 assert(termsmark[vbase[k]]);
      }

      for( k = 0; k < nnodes; k++ )
      {

         connected[k] = FALSE;
	 vcount[k] = 0;
	 SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[k]) );
	 gnodearr[k]->number = k;
         if( g->term[k] == -1 )
         {
	    //SCIP_CALL( SCIPallocBuffer(scip, &nodeterms[k][0]) );
	    //nodeterms[k][0]->number = k;
	    //nodeterms[k][0]->dist = vnoi[k].dist;
	    distarr[k][0] = vnoi[k].dist;
            edgearr[k][0] = vnoi[k].edge;
            basearr[k][0] = vbase[k];
            nodenterms[k] = 1;
	    //printf(" add %d to list %d \n", vbase[k], k );
         }
         else
         {
            nodenterms[k] = 0;
	    edgearr[k][0] = UNKNOWN;
            termsmark[k] = FALSE;
         }
         state[k] = UNKNOWN;
         vcost[k] = vnoi[k].dist;
         vnoi[k].dist = FARAWAY;
         /* vnoi[k].edge = -1;*/
      }

      /* for each terminal: extend the voronoi regions until all neighbouring terminals have been visited */
      for( i = 0; i < nterms; i++ )
      {
         //printf("term: %d \n", terms[i]);
         /* for( k = 0; k < nnodes; k++ )
            {

            assert(!termsmark[k]);
            //termsmark[k] = FALSE;
            assert(!visited[k]);
            } */
         nneighbterms = 0;
         nneighbnodes = 0;
         nreachednodes = 0;

         /* DFS (starting from terminal i) until the entire voronoi region has been visited */
         tovisit[0] = terms[i];
         ntovisit = 1;
         visited[terms[i]] = TRUE;
         state[terms[i]] = CONNECT;
         while( ntovisit > 0 )
         {
            /* iterate all incident edges */
            old = tovisit[--ntovisit];
            oedge = g->outbeg[old];
            while( oedge >= 0 )
            {
               k = g->head[oedge];

               /* is node k in the voronoi region of the i-th terminal ? */
               if( vbase[k] == terms[i] )
               {
                  if( !visited[k] )
                  {
                     state[k] = CONNECT;
		     //printf("   vnoinode: %d ",k);
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
                     vnoi[k].dist = vcost[old] + ((vbase[k] == g->source[0])? cost[oedge] : costrev[oedge]);
                     vnoi[k].edge = oedge;

                     if( termsmark[vbase[k]] == FALSE )
                     {
                        termsmark[vbase[k]] = TRUE;
                        nneighbterms++;
                     }
                     tovisit[nnodes - (++nneighbnodes)] = k;
                  }
                  else
                  {
                     /* if edge 'oedge' allows a shorter connection of node k, update */
                     if( SCIPisGT(scip, vnoi[k].dist, vcost[old] + ((vbase[k] == g->source[0])? cost[oedge] : costrev[oedge])) )  /* TODO: chg s.t. edgecosts for base[root] are reversed */
                     {
                        vnoi[k].dist = vcost[old] + ((vbase[k] == g->source[0])? cost[oedge] : costrev[oedge]);
                        vnoi[k].edge = oedge;
                     }
                  }
               }
               oedge = g->oeat[oedge];
            }
         }

         count = 0;
         for( j = 0; j < nneighbnodes; j++ )
         {

            heap_add(g->path_heap, state, &count, tovisit[nnodes - j - 1], vnoi);
            //   printf( "heap add %d cost %e\n", tovisit[nnodes - j - 1], vnoi[tovisit[nnodes - j - 1]].dist);
         }
         //printf("   nneighbterms %d \n", nneighbterms);
         //SCIP_CALL( voronoi_extend3(scip, g, ((vbase[k] == g->source[0])? cost : costrev), vnoi, nodeterms, basearr, edgearr, termsmark, reachednodes, &nreachednodes, nodenterms,
         //    nneighbterms, terms[i], nneighbnodes) );

         SCIP_CALL( voronoi_extend2(scip, g, ((vbase[k] == g->source[0])? cost : costrev), vnoi, distarr, basearr, edgearr, termsmark, reachednodes, &nreachednodes, nodenterms,
               nneighbterms, terms[i], nneighbnodes) );

         // printf( "terminal: %d\n", terms[i]);
         reachednodes[nreachednodes++] = terms[i];

         for( j = 0; j < nreachednodes; j++ )
         {
            //reachednodes[j] = k; restore vnoi!! TODO
            //	 printf( "reachdnode: : %d, vnoibase: %d \n", reachednodes[j], vbase[reachednodes[j]]);

            vnoi[reachednodes[j]].dist = FARAWAY;
            state[reachednodes[j]] = UNKNOWN;
            visited[reachednodes[j]] = FALSE;
         }

         for( j = 0; j < nneighbnodes; j++ ) // TODO AVOID DOUBLE WORK
         {
            vnoi[tovisit[nnodes - j - 1]].dist = FARAWAY;
            state[tovisit[nnodes - j - 1]] = UNKNOWN;
            visited[tovisit[nnodes - j - 1]] = FALSE;
         }
      }


      /* for each node v: sort the terminal arrays according to their distance to v */
      for( i = 0; i < nnodes; i++ )
      {
         //SCIPsortPtrIntInt((void**)nodeterms[i], edgearr[i], basearr[i], GNODECmpByDist, nodenterms[i]);
	 SCIPsortRealIntInt(distarr[i], basearr[i], edgearr[i], nodenterms[i]);
	 //printf(" node %d terms: ", i);
         //for( j = 0; j < nodenterms[i]; j++ )
         //printf( " %d, ", basearr[i][j]);
         //printf(" \n ");
      }

   }

   /** PHASE II **/
   else
   {
      for( k = 0; k < nnodes; k++ )
      {
         connected[k] = FALSE;
	 vcount[k] = 0;
      }
   }

   connected[start] = TRUE;
   gnodearr[start]->dist = distarr[start][0];
   //SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[start][0]) );
   SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[start]) );

   //printf("inserted\n");
   while( SCIPpqueueNElems(pqueue) > 0 )
   {
      best = ((GNODE*) SCIPpqueueRemove(pqueue))->number;
      //printf("best: %d\n", best);
      term = basearr[best][vcount[best]];

      /* has the terminal already been connected? */
      if( !connected[term] )
      {
         /* connect the terminal */
         //k = g->tail[nodeterms[best][ncheckedterms[best]]->edge];
         k = g->tail[edgearr[best][vcount[best]]];
         while( k != term )
         {
            j = 0;
            // printf(" term  : %d != term %d \n", nodeterms[k][ncheckedterms[k] + j]->terminal, best->terminal);
	    while( basearr[k][vcount[k] + j] != term )
            {
               j++;
            }
            //printf("  move on node %d \n", k);

            if( !connected[k] )
            {
	       assert(vcount[k] == 0);

               connected[k] = TRUE;
	       while( vcount[k] < nodenterms[k] && connected[basearr[k][vcount[k]]] )
               {
	          vcount[k]++;
		  j--;
               }

               if( vcount[k] < nodenterms[k] )
	       {
		  gnodearr[k]->dist = distarr[k][vcount[k]];
		  SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
                  //SCIP_CALL( SCIPpqueueInsert(pqueue, &(nodesid[k])) );
                  //SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[k][vcount[k]]) );
	       }
	    }

            assert( vcount[k] + j < nodenterms[k] );
	    k = g->tail[edgearr[k][vcount[k] + j]];
         }
         /* finally, connected the terminal reached*/
         //assert( k == best->terminal );
         assert( k == term );
         if( !connected[k] )
         {
            connected[k] = TRUE;
	    //printf("  move on term: %d \n", k);
	    //assert( ncheckedterms[k] == 0 );
	    assert( vcount[k] == 0 );
            while( vcount[k] < nodenterms[k] && connected[basearr[k][vcount[k]]] )
            {
	       vcount[k]++;
            }
            if( vcount[k] < nodenterms[k] )
	    {
	       gnodearr[k]->dist = distarr[k][vcount[k]];
	       SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
               //SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[k][vcount[k]]) );
	    }
	 }

      }
      /*k = best;
        while( ncheckedterms[k] + 1 < nodenterms[k] )
        {
        if( !connected[nodeterms[k][++ncheckedterms[k]]->terminal] )
        {
        SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[k][ncheckedterms[k]]) );
        break;
        }
        } */
      while( vcount[best] + 1 < nodenterms[best] )
      {
         if( !connected[basearr[best][++vcount[best]]] )
         {
	    gnodearr[best]->dist = distarr[best][vcount[best]];
	    SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[best]) );
            //SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[best][vcount[best]]) );
            break;
         }
      }
      // printf("after reinsert, k : %d \n", k);
   }

   if( firstrun )
   {
      SCIPfreeBufferArray(scip, &vnoi);
      SCIPfreeBufferArray(scip, &vbase);
      SCIPfreeBufferArray(scip, &terms);
      SCIPfreeBufferArray(scip, &termsmark);
      SCIPfreeBufferArray(scip, &visited);
      SCIPfreeBufferArray(scip, &tovisit);
      SCIPfreeBufferArray(scip, &reachednodes);
      SCIPfreeBufferArray(scip, &vcost);

   }
   if( lastrun || DEFAULT_INITRUNS == 1 )
   {
      SCIPpqueueFree(&pqueue);
      for( i = nnodes - 1; i >= 0; i-- )
      {

         /* for( j = 0; j < nodenterms[i]; j++ )
            {
            SCIPfreeBuffer(scip, &nodeterms[i][j]);
            }*/

         SCIPfreeBuffer(scip, &gnodearr[i]);
	 //SCIPfreeBufferArray(scip, &nodeterms[i]);
	 SCIPfreeBufferArray(scip, &distarr[i]);
	 SCIPfreeBufferArray(scip, &edgearr[i]);
	 SCIPfreeBufferArray(scip, &basearr[i]);

      }
      // SCIPfreeBufferArray(scip, &nodeterms);
      SCIPfreeBufferArray(scip, &distarr);
      SCIPfreeBufferArray(scip, &edgearr);
      SCIPfreeBufferArray(scip, &basearr);
      SCIPfreeBufferArray(scip, &gnodearr);
      SCIPfreeBufferArray(scip, &vcount);
      SCIPfreeBufferArray(scip, &nodesid);
      //SCIPfreeBufferArray(scip, &nodeterms);
      SCIPfreeBufferArray(scip, &nodenterms);
      //SCIPfreeBufferArray(scip, &ncheckedterms);
   }

#else

   for( i = 0; i < nnodes; i++ )
   {
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
   }
   connected[start] = TRUE;


   /* CONSTCOND */
   for( ;; )
   {
      /* Find a terminal with minimal distance to the current ST
       */
      min = FARAWAY;
      old = -1;
      newval = -1;

      for( i = 0; i < nnodes; i++ )
      {
         if( g->grad[i] == 0 || g->term[i] != layer || connected[i] )
            continue;


         /*
          */
         if( path[i] == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(path[i]), nnodes) );

            assert(path[i] != NULL);
            if( g->source[0] == i )
               graph_path_exec(g, FSP_MODE, i, cost, path[i]);
            else
               graph_path_exec(g, FSP_MODE, i, costrev, path[i]);
         }
         for( k = 0; k < csize; k++ )
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

   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );

   /* compute MST */
   for(i = 0; i < nnodes; i++)
      g->mark[i] = connected[i];

   assert(g->source[layer] >= 0);
   assert(g->source[layer] <  nnodes);

   graph_path_exec(g, MST_MODE, g->source[layer], g->cost, mst);

   for(i = 0; i < nnodes; i++)
   {
      if (connected[i] && (mst[i].edge != -1))
      {
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == -1);

         result[mst[i].edge] = layer;
      }
   }

   /* prune */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for( i = 0; i < nnodes; i++ )
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

   SCIPfreeBufferArray(scip, &mst);

   SCIPfreeBufferArray(scip, &cluster);
   return SCIP_OKAY;
}
#endif

/* prune the Steiner Tree in such a way, that all leaves are terminals */
static
SCIP_RETCODE do_prune(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int                   layer,
   int*                  result,             /**< ST edges */
   char*                 connected           /**< ST nodes */
   )
{
   PATH*  mst;
   int i;
   int j;
   int count;
   int nnodes;
   nnodes = g->knots;
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );

   /* compute the MST */
   for( i = 0; i < nnodes; i++ )
      g->mark[i] = connected[i];

   assert(g->source[layer] >= 0);
   assert(g->source[layer] <  nnodes);

   graph_path_exec(g, MST_MODE, g->source[layer], g->cost, mst);

   for( i = 0; i < nnodes; i++ )
   {
      if( connected[i] && (mst[i].edge != -1) )
      {
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == -1);

         result[mst[i].edge] = layer;
      }
   }

   /* prune */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for( i = 0; i < nnodes; i++ )
      {
         if( !g->mark[i] )
            continue;

         if( g->term[i] == layer )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if (result[j] == layer)
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == layer )
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
   while( count > 0 );

   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;

}

/* pure TM heuristic */
static
SCIP_RETCODE do_tm(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   PATH**                path,
   const SCIP_Real*      cost,
   const SCIP_Real*      costrev,
   int                   layer,
   int                   start,
   int*                  result,
   char*                 connected
   )
{
   int*   cluster;
   int    csize = 0;
   int    k;
   int    e;
   int    i;
   int    j;
   int    old;
   int    newval;
   int nnodes;
   SCIP_Real min;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   assert(path      != NULL);
   assert(layer >= 0 && layer < g->layers);
   nnodes = g->knots;

   SCIPdebugMessage("Heuristic: Start=%5d ", start);

   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, nnodes) );

   cluster[csize++] = start;

   for( i = 0; i < nnodes; i++ )
   {
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
   }
   connected[start] = TRUE;


   /* CONSTCOND */
   for( ;; )
   {
      /* Find a terminal with minimal distance to the current ST
       */
      min = FARAWAY;
      old = -1;
      newval = -1;

      for( i = 0; i < nnodes; i++ )
      {
         if( g->grad[i] == 0 || g->term[i] != layer || connected[i] )
            continue;

         /*
          */
         if( path[i] == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(path[i]), nnodes) );

            assert(path[i] != NULL);
            if( g->source[0] == i )
               graph_path_exec(g, FSP_MODE, i, cost, path[i]);
            else
               graph_path_exec(g, FSP_MODE, i, costrev, path[i]);
         }
         for( k = 0; k < csize; k++ )
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
   SCIPfreeBufferArray(scip, &cluster);

   SCIP_CALL( do_prune(scip, g, cost, layer, result, connected) );

   return SCIP_OKAY;
}


static
SCIP_RETCODE do_tm_polzin(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_PQUEUE*          pqueue,
   GNODE**               gnodearr,
   const SCIP_Real*      cost,
   const SCIP_Real*      costrev,
   int                   layer,
   SCIP_Real**           distarr,
   int                   start,
   int*                  result,
   int*                  vcount,
   int*                  nodenterms,
   int**                 basearr,
   int**                 edgearr,
   char                  firstrun,
   char*                 connected
   )
{
   PATH* vnoi;
   SCIP_Real* vcost;
   int    k;
   int    i;
   int    j;
   int    old;
   int    best;
   int    term;
   int    count;
   int   oedge;
   int   nnodes;
   int   nterms;
   int   ntovisit;
   int   nneighbnodes;
   int   nneighbterms;
   int   nreachednodes;
   int*  state;
   int*  vbase;
   int*  terms;
   int*  tovisit;
   int*  reachednodes;
   char* termsmark;
   char* visited;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   nnodes = g->knots;
   nterms = g->terms;

   SCIPdebugMessage("TM_Polzin Heuristic: Start=%5d ", start);

   /* if the heuristic is called for the first time several data structures have to be set up */
   if( firstrun )
   {
      /* PHASE I: */
      for( i = 0; i < nnodes; i++ )
      {
         g->mark[i] = (g->grad[i] > 0);
      }

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
         if( g->term[i] != -1 )
         {
            termsmark[i] = TRUE;
            terms[j++] = i;
         }
         else
         {
            termsmark[i] = FALSE;
         }
      }
      assert(j == nterms);

      voronoi(g, cost, costrev, termsmark, vbase, vnoi); /* TODO: chg s.t. edgecosts for base[root] are reversed */
      state = g->path_state;
      for( k = 0; k < nnodes; k++ )
      {
	 assert(termsmark[vbase[k]]);
      }

      for( k = 0; k < nnodes; k++ )
      {
         connected[k] = FALSE;
	 vcount[k] = 0;
	 SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[k]) );
	 gnodearr[k]->number = k;
         if( g->term[k] == -1 )
         {
	    distarr[k][0] = vnoi[k].dist;
            edgearr[k][0] = vnoi[k].edge;
            basearr[k][0] = vbase[k];
            nodenterms[k] = 1;
	    //printf(" add %d to list %d \n", vbase[k], k );
         }
         else
         {
            nodenterms[k] = 0;
	    edgearr[k][0] = UNKNOWN;
            termsmark[k] = FALSE;
         }
         state[k] = UNKNOWN;
         vcost[k] = vnoi[k].dist;
         vnoi[k].dist = FARAWAY;
         /* vnoi[k].edge = -1;*/
      }

      /* for each terminal: extend the voronoi regions until all neighbouring terminals have been visited */
      for( i = 0; i < nterms; i++ )
      {
         //printf("term: %d \n", terms[i]);
         nneighbterms = 0;
         nneighbnodes = 0;
         nreachednodes = 0;

         /* DFS (starting from terminal i) until the entire voronoi region has been visited */
         tovisit[0] = terms[i];
         ntovisit = 1;
         visited[terms[i]] = TRUE;
         state[terms[i]] = CONNECT;
         while( ntovisit > 0 )
         {
            /* iterate all incident edges */
            old = tovisit[--ntovisit];
            oedge = g->outbeg[old];
            while( oedge >= 0 )
            {
               k = g->head[oedge];

               /* is node k in the voronoi region of the i-th terminal ? */
               if( vbase[k] == terms[i] )
               {
                  if( !visited[k] )
                  {
                     state[k] = CONNECT;
		     //printf("   vnoinode: %d ",k);
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
                     vnoi[k].dist = vcost[old] + ((vbase[k] == g->source[0])? cost[oedge] : costrev[oedge]);
                     vnoi[k].edge = oedge;

                     if( termsmark[vbase[k]] == FALSE )
                     {
                        termsmark[vbase[k]] = TRUE;
                        nneighbterms++;
                     }
                     tovisit[nnodes - (++nneighbnodes)] = k;
                  }
                  else
                  {
                     /* if edge 'oedge' allows a shorter connection of node k, update */
                     if( SCIPisGT(scip, vnoi[k].dist, vcost[old] + ((vbase[k] == g->source[0])? cost[oedge] : costrev[oedge])) )  /* TODO: chg s.t. edgecosts for base[root] are reversed */
                     {
                        vnoi[k].dist = vcost[old] + ((vbase[k] == g->source[0])? cost[oedge] : costrev[oedge]);
                        vnoi[k].edge = oedge;
                     }
                  }
               }
               oedge = g->oeat[oedge];
            }
         }

         count = 0;
         for( j = 0; j < nneighbnodes; j++ )
         {
            heap_add(g->path_heap, state, &count, tovisit[nnodes - j - 1], vnoi);
            //   printf( "heap add %d cost %e\n", tovisit[nnodes - j - 1], vnoi[tovisit[nnodes - j - 1]].dist);
         }
         SCIP_CALL( voronoi_extend2(scip, g, ((vbase[k] == g->source[0])? cost : costrev), vnoi, distarr, basearr, edgearr, termsmark, reachednodes, &nreachednodes, nodenterms,
               nneighbterms, terms[i], nneighbnodes) );
         // printf( "terminal: %d\n", terms[i]);
         reachednodes[nreachednodes++] = terms[i];

         for( j = 0; j < nreachednodes; j++ )
         {
            //	 printf( "reachdnode: : %d, vnoibase: %d \n", reachednodes[j], vbase[reachednodes[j]]);
            vnoi[reachednodes[j]].dist = FARAWAY;
            state[reachednodes[j]] = UNKNOWN;
            visited[reachednodes[j]] = FALSE;
         }

         for( j = 0; j < nneighbnodes; j++ ) // TODO AVOID DOUBLE WORK
         {
            vnoi[tovisit[nnodes - j - 1]].dist = FARAWAY;
            state[tovisit[nnodes - j - 1]] = UNKNOWN;
            visited[tovisit[nnodes - j - 1]] = FALSE;
         }
      }

      /* for each node v: sort the terminal arrays according to their distance to v */
      for( i = 0; i < nnodes; i++ )
      {
	 SCIPsortRealIntInt(distarr[i], basearr[i], edgearr[i], nodenterms[i]);
	 //printf(" node %d terms: ", i);
         //for( j = 0; j < nodenterms[i]; j++ )
         //printf( " %d, ", basearr[i][j]);
         //printf(" \n ");
      }

   }

   /** PHASE II **/
   else
   {
      for( k = 0; k < nnodes; k++ )
      {
         connected[k] = FALSE;
	 vcount[k] = 0;
      }
   }

   connected[start] = TRUE;
   gnodearr[start]->dist = distarr[start][0];
   SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[start]) );

   //printf("inserted\n");
   while( SCIPpqueueNElems(pqueue) > 0 )
   {
      best = ((GNODE*) SCIPpqueueRemove(pqueue))->number;
      //printf("best: %d\n", best);
      term = basearr[best][vcount[best]];

      /* has the terminal already been connected? */
      if( !connected[term] )
      {
         /* connect the terminal */
         k = g->tail[edgearr[best][vcount[best]]];
         while( k != term )
         {
            j = 0;
            // printf(" term  : %d != term %d \n", nodeterms[k][ncheckedterms[k] + j]->terminal, best->terminal);
	    while( basearr[k][vcount[k] + j] != term )
            {
               j++;
            }
            //printf("  move on node %d \n", k);

            if( !connected[k] )
            {
	       assert(vcount[k] == 0);

               connected[k] = TRUE;
	       while( vcount[k] < nodenterms[k] && connected[basearr[k][vcount[k]]] )
               {
	          vcount[k]++;
		  j--;
               }

               if( vcount[k] < nodenterms[k] )
	       {
		  gnodearr[k]->dist = distarr[k][vcount[k]];
		  SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
                  //SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[k][vcount[k]]) );
	       }
	    }

            assert( vcount[k] + j < nodenterms[k] );
	    k = g->tail[edgearr[k][vcount[k] + j]];
         }
         /* finally, connected the terminal reached*/
         //assert( k == best->terminal );
         assert( k == term );
         if( !connected[k] )
         {
            connected[k] = TRUE;
	    //printf("  move on term: %d \n", k);
	    //assert( ncheckedterms[k] == 0 );
	    assert( vcount[k] == 0 );
            while( vcount[k] < nodenterms[k] && connected[basearr[k][vcount[k]]] )
            {
	       vcount[k]++;
            }
            if( vcount[k] < nodenterms[k] )
	    {
	       gnodearr[k]->dist = distarr[k][vcount[k]];
	       SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
               //SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[k][vcount[k]]) );
	    }
	 }
      }

      while( vcount[best] + 1 < nodenterms[best] )
      {
         if( !connected[basearr[best][++vcount[best]]] )
         {
	    gnodearr[best]->dist = distarr[best][vcount[best]];
	    SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[best]) );
            //SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[best][vcount[best]]) );
            break;
         }
      }
      // printf("after reinsert, k : %d \n", k);
   }

   if( firstrun )
   {
      SCIPfreeBufferArray(scip, &vnoi);
      SCIPfreeBufferArray(scip, &vbase);
      SCIPfreeBufferArray(scip, &terms);
      SCIPfreeBufferArray(scip, &termsmark);
      SCIPfreeBufferArray(scip, &visited);
      SCIPfreeBufferArray(scip, &tovisit);
      SCIPfreeBufferArray(scip, &reachednodes);
      SCIPfreeBufferArray(scip, &vcost);
   }

   /* prune the ST, so that all leaves are terminals */
   SCIP_CALL( do_prune(scip, g, cost, layer, result, connected) );

   return SCIP_OKAY;
}


/* local heuristics */
static
SCIP_RETCODE do_local(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*  graph,
   const SCIP_Real* cost,
   const SCIP_Real* costrev,
   int*          best_result
   )
{
   NODE* nodes;
   SCIP_Real obj;
   int e;
   int i;
   int k;
   int root;
   int nnodes;
   int nedges;
   char* steinertree;

   root = graph->source[0];
   nnodes = graph->knots;
   nedges = graph->edges;
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &steinertree, nnodes) );
   for( i = 0; i < nnodes; i++ )
   {
      steinertree[i] = FALSE;
      init(&nodes[i]);
   }

   /* create a link-cut tree representing the current Steiner tree */
   for( e = 0; e < nedges; e++ )
   {
      assert(graph->head[e] == graph->tail[flipedge(e)]);

      /* if edge e is in the tree, so are its incident vertices */
      if( best_result[e] != -1 )
      {
         steinertree[graph->tail[e]] = TRUE;
         steinertree[graph->head[e]] = TRUE;
         link(&nodes[graph->head[e]], &nodes[graph->tail[e]], flipedge(e));
      }
   }
   assert( nodes[root].edge == -1 );
   nodes[root].edge = -1;

   /** VERTEX  INSERTION */
   if( 0 && !graph->rootisfixed )  /* TODO adapt function to directed graph */
   {
      int newnode = 0;
      int oedge;
      int* insert;
      int* adds;
      int* cuts;
      int counter;
      int insertcount;
      NODE* v;
      NODE* w;
      NODE* max;
      SCIP_Real diff;

      SCIP_CALL( SCIPallocBufferArray(scip, &insert, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &adds, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cuts, nnodes) );

      i = 0;
      for( ;; )
      {
         /* if vertex i is not in the current ST and has at least two adjacent nodes, it might be added */
         if( !steinertree[i] && graph->grad[i] > 1 )
         {
            insertcount = 0;

            /* if an outgoing edge of vertex i points to the current ST, link the edge to a list */
            oedge = graph->outbeg[i];
            while( oedge >= 0 )
            {
               if( steinertree[graph->head[oedge]] )
               {
                  insert[insertcount++] = oedge;
               }
               oedge = graph->oeat[oedge];
            }

            /* if there are at least two edges connecting node i and the current tree, start the insertion process */
            if( insertcount > 1 )
            {
#if 0
	       int* tmpcosts;

	       SCIP_CALL( SCIPallocBufferArray(scip, &tmpcosts, insertcount) );
               for( k = 0; k < insertcount; k++ )
		  tmpcosts[k] = graph->cost[insert[k]];


               SCIPsortIntInt( tmpcosts, insert, insertcount );
               for( k = 0; k < insertcount -1 ; k++ )
                  assert( graph->cost[insert[k]] <= graph->cost[insert[k+1] ]);
               SCIPfreeBufferArray(scip, &tmpcosts);
#endif
               diff = 0.0;

               /* the node to insert */
               v = &nodes[i];

               link(v, &nodes[graph->head[insert[0]]], insert[0]);
               diff +=  graph->cost[v->edge];
	       /*
                 if( SCIPisGE(scip, cost[v->edge], 1e+10 - 1) || SCIPisGE(scip, costrev[v->edge], 1e+10 - 1) )
                 {
                 printf(" XRX ");
                 }*/

               counter = 0;
               for( k = 1; k < insertcount; k++ )
               {
                  evert(v);

                  /* next vertex in the current Steiner tree adjacent to vertex i resp. v (the one being scrutinized for possible insertion) */
                  w = &nodes[graph->head[insert[k]]];

                  /* if there is an edge with cost greater than that of the current edge... */
                  max = findMax(scip, graph->cost, w);
                  if( SCIPisGT(scip, graph->cost[max->edge], graph->cost[insert[k]]) )
                  {
                     diff += graph->cost[insert[k]];
                     diff -= graph->cost[max->edge];
                     /*    if( SCIPisGE(scip, cost[insert[k]], 1e+10 - 1) || SCIPisGE(scip, costrev[insert[k]], 1e+10 - 1) )
                           {
                           printf(" OFB  \n");
                           } */

                     cuts[counter] = max->edge;
                     cut(max);
                     link(v, w, insert[k]);
                     adds[counter++] = v->edge;
                  }
               }

               /* if the new tree is more expensive than the old one, restore the latter */
               if( !SCIPisNegative(scip, diff) )
               {

                  evert(v);
                  for( k = counter - 1; k >= 0; k-- )
                  {
                     cut(&nodes[graph->head[adds[k]]]);
                     evert(&nodes[graph->tail[cuts[k]]]);
                     link(&nodes[graph->tail[cuts[k]]], &nodes[graph->head[cuts[k]]], cuts[k]);
                  }

                  /* finally, cut the edge added first (if it had been cut during the insertion process, it will have been restored above) */
                  evert(v);
                  cut(&nodes[graph->head[insert[0]]]);
               }
               else
               {
                  /* check if a fixed edge has been added */
                  evert(&nodes[root]);
                  adds[counter] = insert[0];
                  for( k = 0; k <= counter; k++ )
                  {
                     if( nodes[i].edge == adds[k] ){
                        if( SCIPisGE(scip, costrev[adds[k]], 1e10 -1) )
                           break;
                     }
                     if( nodes[i].edge == flipedge(adds[k]) ){
                        if( SCIPisGE(scip, cost[adds[k]], 1e10 -1) )
                           break;
                     }
                     if( i == graph->head[adds[k]] ){
                        if( SCIPisGE(scip, costrev[adds[k]], 1e10 -1) )
                           break;
                     }
                     if( i == graph->head[flipedge(adds[k])] ){
                        if( SCIPisGE(scip, cost[adds[k]], 1e10 -1) )
                           break;
                     }
                  }

                  if( k != counter + 1 )
                  {
                     printf("RSTORING OFB \n\n");
                     evert(v);
                     for( k = counter - 1; k >= 0; k-- )
                     {
                        cut(&nodes[graph->head[adds[k]]]);
                        evert(&nodes[graph->tail[cuts[k]]]);
                        link(&nodes[graph->tail[cuts[k]]], &nodes[graph->head[cuts[k]]], cuts[k]);
                     }

                     /* finally, cut the edge added first (if it had been cut during the insertion process, it will have been restored above) */
                     evert(v);
                     cut(&nodes[graph->head[insert[0]]]);
                  }
                  else
                  {
                     newnode = i;
                     steinertree[i] = TRUE;
                     printf("ADDED VERTEX \n");
                  }
                  /* TODO adjust tree st we only have to adjust best_result for the new edges*/
               }
            }
         }

         if( i < nnodes - 1 )
         {
            i++;
         }
         else
         {
            i = 0;
         }
         if( newnode == i )
         {
            break;
         }
         if( i == 0 )
            printf("VertInsert newrun \n");
      }

      evert(&nodes[root]);

      for( e = 0; e < nedges; e++ )
      {
         best_result[e] = -1;
      }
      for( i = 0; i < nnodes; i++ )
      {
         if( steinertree[i] && nodes[i].edge != -1)
            best_result[flipedge(nodes[i].edge)] = 0;
      }
      SCIPfreeBufferArray(scip, &insert);
      SCIPfreeBufferArray(scip, &cuts);
      SCIPfreeBufferArray(scip, &adds);

      obj = 0.0;
      for( e = 0; e < nedges; e++)
         obj += (best_result[e] > -1) ? graph->cost[e] : 0.0;

      printf(" ObjAfterVertexInsertion=%.12e\n", obj);

   }


   /** Key-Vertex Elimination & Key-Path Exchange*/
   if( 1 && !graph->rootisfixed )
   {
      /* TODO declare variable needed only once, seperately*/
      IDX* blists_curr;
      IDX** blists_start;  /* array [1,..,nnodes],
                            * if node i is in the current ST, blists_start[i] points to a linked list of all nodes having i as their base */
      PATH* mst;           /* minimal spanning tree structure */
      GRAPH* supergraph;
      IDX** lvledges_start;  /* horizontal edges */
      IDX* lvledges_curr;
      PHNODE** boundpaths;
      UF uf;  /* union-find*/

      SCIP_Real bestdiff = 0;
      SCIP_Real* memdist;
      int* supernodes;
      int* supernodesid;
      int* heapsize;
      int* boundedges;
      int* memvbase;
      int* meminedges;
      int* kpedges;
      int* kpnodes;
      int* newedges;
      int* vbase;     /* array [1,...,nnodes] */
      int* state;
      int* graphmark;

      int node;
      int nresnodes;
      int kptailnode;  /* tail node of the current keypath*/
      int crucnode;   /* current crucial node*/
      int adjnode;
      int edge;
      int count;
      int newedge;
      int oldedge;
      int nsupernodes;
      int nkpedges;
      int nstnodes;
      int nkpnodes;
      int nboundedges;
      int rootpathstart;
      int l;
      SCIP_Real kpcost;
      SCIP_Real mstcost;
      SCIP_Real edgecost;
      PATH* vnoi;
      int* dfstree;
      char* nodesmark;
      char* pinned;
      char* scanned;
      int localmoves = 2;
      char debg = FALSE; //TRUE;
      int nruns;
      int newpathend = -1;
      SCIP_Real kpathcost;

      obj = 0.0;

      for( e = 0; e < nedges; e++)
      {
         // if(best_result[e] > -1)
	 //  printf("st edge: %d->%d \n", graph->tail[e], graph->head[e]);
         obj += (best_result[e] > -1) ? graph->cost[e] : 0.0;
      }
      if( debg )
         printf(" ObjBEFKEYVertexELimination=%.12e\n", obj);
      graphmark = graph->mark;
      //  SCIP_CALL( printGraph(scip, graph, "abef.gml", best_result) );

      /* allocate memory */

      /* only needed for Key-Path Elimination */
      SCIP_CALL( SCIPallocBufferArray(scip, &newedges, nedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lvledges_start, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &boundedges, nedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &supernodesid, nnodes) );

      /* only needed for Key-Path Exchange */

      /* memory needed for both Key-Path Elimination and Exchange */
      SCIP_CALL( SCIPallocBufferArray(scip, &scanned, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &heapsize, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &blists_start, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );

      SCIP_CALL( SCIPallocBufferArray(scip, &memvbase, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &memdist, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &meminedges, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &boundpaths, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pinned, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &dfstree, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nodesmark, nnodes) );

      /*  const char baseo[] = "X/org";
          char filenameo [ FILENAME_MAX ];
          sprintf(filenameo, "%s%d.gml", baseo, i);
          SCIP_CALL( printGraph(scip, graph, filenameo, best_result) );*/

      for( nruns = 0; nruns < 3 && localmoves > 0; nruns++ )
      {
         localmoves = 0;

         /* initialize data structures */
         uf_init(scip, &uf, nnodes);

         //BMSclearMemoryArray(lvledges_start, nnodes); FASTER?? TODO
         BMSclearMemoryArray(blists_start, nnodes);

         /* find a DFS order of the ST nodes */
         nstnodes = 0;
         dfsorder(graph, best_result, &(root), &nstnodes, dfstree);

         SCIP_CALL( SCIPallocBufferArray(scip, &supernodes, nstnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &kpnodes, nstnodes) );
         /* TODO warum nich nstnodes? e07 fehler:::::delete kpedges */
         SCIP_CALL( SCIPallocBufferArray(scip, &kpedges, nstnodes) );

         /* compute a voronoi diagram with the ST nodes as bases */
         voronoi(graph, graph->cost, graph->cost, steinertree, vbase, vnoi);
         state = graph->path_state;

         //printf( "root %d \n \n ", root);

         /* initialize data structures  */
         for( k = 0; k < nnodes; k++ )
         {
            assert(state[k] == CONNECT);
	    graphmark[k] = TRUE;
            pinned[k] = FALSE;
            nodesmark[k] = FALSE;
            scanned[k] = FALSE;

            /* initialize pairing heaps */
            heapsize[k] = 0;
            boundpaths[k] = NULL;

	    lvledges_start[k] = NULL;

            /* link all nodes to their (respective) voronoi base */
            SCIP_CALL( SCIPallocMemory(scip, &blists_curr) ); /* TODO extra method (inline) block memory? */
            blists_curr->index = k;
            blists_curr->parent = blists_start[vbase[k]];
            blists_start[vbase[k]] = blists_curr;
         }

         /* for each node, store all of its outgoing boundary-edges in a (respective) heap*/
         for( e = 0; e < nedges; e += 2 )
         {
            node = graph->tail[e];
            adjnode = graph->head[e];
            newedges[e] = UNKNOWN;
            newedges[e + 1] = UNKNOWN;

            /* is edge 'e' a boundary-edge? */
            if( vbase[node] != vbase[adjnode] )
            {
               edgecost = vnoi[node].dist + graph->cost[e] + vnoi[adjnode].dist;
               //printf("put in pairheap[%d]: %d_%d cost : %f bases: %d %d \n ", vbase[node], node, adjnode, edgecost, vbase[graph->tail[e]], vbase[graph->head[e]] );
               //printf("put in pairheap[%d]: %d_%d cost : %f bases: %d %d \n ", vbase[adjnode], graph->tail[flipedge(e)], graph->head[flipedge(e)], edgecost, vbase[graph->tail[flipedge(e)]], vbase[graph->head[flipedge(e)]] );

               /* add the boundary-edge 'e' and its reversed to the corresponding heaps */
               phnode_insert(scip, &boundpaths[vbase[node]], e, edgecost, &(heapsize[vbase[node]]));
               phnode_insert(scip, &boundpaths[vbase[adjnode]], flipedge(e), edgecost, &(heapsize[vbase[adjnode]]));
            }
         }

         /* find LCAs for all edges */
         SCIP_CALL( lca(scip, graph, root, &uf, nodesmark, best_result, lvledges_start, boundpaths, heapsize, vbase) );

         /* henceforth, the union-find structure will be used on the ST */
         uf_free(scip, &uf);
         uf_init(scip, &uf, nnodes);

         /* henceforth, nodesmark will be used to mark the current supervertices (except for the one representing the root-component) */
         for( i = 0; dfstree[i] != root; i++ )
         {
            nodesmark[dfstree[i]] = FALSE;
         }
         nodesmark[dfstree[i]] = FALSE;

         /* debug test; to be deleted later on TODO */
         assert(dfstree[i] == root);
         for( k = 0; k < nnodes; k++ )
         {
            assert( !nodesmark[k] );
         }

         /* main loop visiting all nodes of the current ST in post-order */
         for( i = 0; dfstree[i] != root; i++ )
         {
            crucnode = dfstree[i];
	    scanned[crucnode] = TRUE;
	    if( debg )
               printf("iteration %d (%d) \n", i, crucnode);

	    /*  has the node been temporarily removed from the ST? */
	    if( !graphmark[crucnode] )
	    {
	       //printf("^ is removed \n");
	       continue;
	    }

	    /* is node 'crucnode' a removable crucial node? (i.e. not pinned or a terminal) */
            if( !pinned[crucnode] && !Is_term(graph->term[crucnode]) && nodeIsCrucial(graph, best_result, crucnode) )
	    {
               /* if current ST node is a terminal or pinned, update union-find structure and heaps before continuing */
               if( !(1) && Is_term(graph->term[crucnode]) )//|| pinned[crucnode] )
               {
                  /* update union-find and pairing heaps: unite terminal 'crucnode' with all of its ancestor key-paths */
                  //printf("^ is pinned or term \n");
                  for( edge = graph->outbeg[crucnode]; edge != EAT_LAST; edge = graph->oeat[edge] )
                  {
                     /* check whether edge 'edge' leads to an ancestor of terminal 'crucnode' */
                     if( best_result[edge] != -1 && steinertree[graph->head[edge]] )
                     {
                        adjnode = graph->head[edge];

                        /* meld the heaps */
                        phnode_meldheaps(scip, &boundpaths[crucnode], &boundpaths[adjnode], &heapsize[crucnode], &heapsize[adjnode]);

                        if( debg )
                           printf( "unite 0 (%d) (%d) \n ",  crucnode, adjnode);
                        /* update the union-find data structure */
                        uf_union(&uf, crucnode, adjnode, FALSE);

                        /* move along the key-path until its end (i.e. until a crucial node is reached) */
                        while( !nodeIsCrucial(graph, best_result, adjnode) && !pinned[adjnode] )
                        {
                           for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                           {
                              if( best_result[e] != -1 )
                                 break;
                           }

                           /* assert that each leaf of the ST is a terminal */
                           /* TODO mustn be true after vertex insertion!!) */
                           assert( e != EAT_LAST );
                           adjnode = graph->head[e];
                           if( !steinertree[adjnode] )
                              break;

                           if( debg )
                              printf( "unite 1 (%d) (%d) \n ",  crucnode, adjnode);
                           /* update the union-find data structure */
                           uf_union(&uf, crucnode, adjnode, FALSE);

                           /* meld the heaps */
                           phnode_meldheaps(scip, &boundpaths[crucnode], &boundpaths[adjnode], &heapsize[crucnode], &heapsize[adjnode]);
                        }
                     }
                  }
                  continue;
               }
               if( !(1) && !nodeIsCrucial(graph, best_result, crucnode) )
               {
                  /* TODO replace*/
                  //printf("^ is not crucial \n");
                  assert( stdeg(graph, best_result, crucnode) <= 2 );
                  continue;
               }
               if (debg)
                  printf("Elimination: %d \n", crucnode);

               /* debug, TODO delete*/
               for( k = 0; k < nnodes; k++ )
               {
                  assert( state[k] == CONNECT );
               }

               /* find all (unique) key-paths starting in node 'crucnode' */
               k = UNKNOWN;
               kpcost = 0.0;
               nkpnodes = 0;
               nkpedges = 0;
               nsupernodes = 0;
               for( edge = graph->outbeg[crucnode]; edge != EAT_LAST; edge = graph->oeat[edge] )
               {
                  /* check whether the outgoing edge is in the ST */
                  if( (best_result[edge] > -1 && steinertree[graph->head[edge]]) || (best_result[flipedge(edge)] > -1 && steinertree[graph->tail[edge]]) )
                  {
                     kpcost += graph->cost[edge];
                     //printf(" kpcost1 +  %f \n", graph->cost[edge]);
                     /* check whether the current edge leads to the ST root*/
                     if( best_result[flipedge(edge)] > -1 )
                     {
                        k = flipedge(edge);
                        kpedges[nkpedges++] = k;
                        assert( edge == nodes[crucnode].edge );
                     }
                     else
                     {
                        kpedges[nkpedges++] = edge;
                        adjnode = graph->head[edge];
                        e = edge;

                        /* move along the key-path until its end (i.e. a crucial or pinned node) is reached */
                        while( !pinned[adjnode] && !nodeIsCrucial(graph, best_result, adjnode) && steinertree[adjnode] )
                        {
                           if( debg )
                              printf( "unite 2 (%d) (%d) \n ",  crucnode, adjnode);
                           /* update the union-find data structure */
                           uf_union(&uf, crucnode, adjnode, FALSE);

                           kpnodes[nkpnodes++] = adjnode;
                           //printf(" kp node: %d \n", adjnode);

                           for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                           {
                              if( best_result[e] > -1 )
                              {
                                 kpcost += graph->cost[e];
                                 //printf(" kpcost %f \n", graph->cost[e]);
                                 kpedges[nkpedges++] = e;
                                 break;
                              }
                           }
                           /* assert that each leaf of the ST is a terminal */
                           assert( e != EAT_LAST );

                           adjnode = graph->head[e];
                        }
                        /* does the last node on the path belong to a removed component? */
                        if( !steinertree[adjnode] )
                        {
                           kpcost -= graph->cost[e];
                           nkpedges--;
                           adjnode = graph->tail[e];
                           if( adjnode != crucnode )
                           {
                              supernodes[nsupernodes++] = adjnode;
                              if( debg )
                                 printf(" (art) supernode: %d \n", adjnode);
                              nodesmark[adjnode] = TRUE;
                           }
                        }
                        else
                        {
                           supernodes[nsupernodes++] = adjnode;
                           if( debg )
                              printf(" supernode: %d \n", adjnode);
                           nodesmark[adjnode] = TRUE;
                        }
                     }
                  }
               }

               /* traverse the key-path leading to the root-component */
               rootpathstart = nkpnodes;
               if( k != -1 )
               {
                  //printf(" rootedge: %d_%d \n", graph->tail[nodes[crucnode].edge], graph->head[nodes[crucnode].edge]);
                  /* begin with the edge starting in the root-component of node 'crucnode' */
                  e = k;
                  adjnode = graph->tail[e];
                  while( !pinned[adjnode] && !nodeIsCrucial(graph, best_result, adjnode) && steinertree[adjnode] )
                  {
                     /* update the union-find data structure */
                     //printf(" kp node root: %d \n", adjnode);
                     kpnodes[nkpnodes++] = adjnode;

                     for( e = graph->inpbeg[adjnode]; e != EAT_LAST; e = graph->ieat[e] )
                     {
                        if( best_result[e] > -1 )//&& steinertree[graph->tail[e]] )
                        {
                           assert(steinertree[graph->tail[e]]);
                           kpcost += graph->cost[e];
                           //printf(" kpcost %f \n", graph->cost[e]);
                           kpedges[nkpedges++] = e;
                           break;
                        }
                     }

                     assert( e != EAT_LAST );
                     adjnode = graph->tail[e];
                  }
                  supernodes[nsupernodes++] = adjnode;
                  if( debg )
                     printf("root supernode: %d \n", graph->tail[e]);
               }

               /* the last of the key-path nodes to be stored is the current key-node */
               kpnodes[nkpnodes++] = crucnode;

               /* number of reset nodes */
               nresnodes = 0;

               /* reset all nodes (referred to as 'C' henceforth) whose bases are internal nodes of the current key-paths */
               for( k = 0; k < nkpnodes; k++ )
               {
                  /* reset all nodes having the current (internal) keypath node as their voronoi base */
                  blists_curr = blists_start[kpnodes[k]];
                  while( blists_curr != NULL )
                  {
                     node = blists_curr->index;
		     assert(graphmark[node]);
                     //printf("C-node %d \n", blists_curr->index);
                     /* store all relevant data */
                     memvbase[nresnodes] = vbase[node];
                     memdist[nresnodes] =  vnoi[node].dist;
                     meminedges[nresnodes] = vnoi[node].edge;
                     nresnodes++;

                     /* reset data */
                     vbase[node] = UNKNOWN;
                     vnoi[node].dist = FARAWAY;
                     vnoi[node].edge = UNKNOWN;
                     state[node] = UNKNOWN;
                     blists_curr = blists_curr->parent;
                  }
               }

               /* add vertical boundary-paths between the child components and the root-component (wrt node 'crucnode') */
               nboundedges = 0;
               for( k = 0; k < nsupernodes - 1; k++ )
               {
                  l = supernodes[k];
                  edge = UNKNOWN;
                  while( boundpaths[l] != NULL )
                  {
                     phnode_deletemin(scip, &edge, &edgecost, &boundpaths[l], &heapsize[l]);

                     node = (vbase[graph->head[edge]] == UNKNOWN)? UNKNOWN : UF_find(&uf, vbase[graph->head[edge]]);
                     adjnode = (vbase[graph->tail[edge]] == UNKNOWN)? UNKNOWN : UF_find(&uf, vbase[graph->tail[edge]]);
                     assert(adjnode == l);
                     /*                    if ( edge != UNKNOWN )
                                           {
                                           printf("min edge from heap[%d]: %d_%d  |  ", l, graph->head[edge], graph->tail[edge]);
                                           printf("vorbases %d_%d  |  ", vbase[graph->head[edge]], vbase[graph->tail[edge]]);
                                           printf("basenodes:  %d_%d\n ", node, adjnode);
                                           }        */
                     /* check whether edge 'edge' represents a boundary-path having an endpoint in the kth-component and in the root-component respectively */
                     if( node != UNKNOWN && !nodesmark[node] && graphmark[node] )//&& !pinned[vbase[graph->tail[edge]]] && !pinned[vbase[graph->tail[edge]]] )
                     {
                        boundedges[nboundedges++] = edge;
                        if( debg )
                           printf("ADD vertical edge: %d_%d  \n", graph->tail[edge], graph->head[edge]);
                        phnode_insert(scip, &boundpaths[l], edge, edgecost, &heapsize[l]);
                        break;
                     }
                  }
               }

               /* add horizontal boundary-paths (between the  child-components) */
               lvledges_curr = lvledges_start[crucnode];
               while( lvledges_curr != NULL )
               {
		  edge = lvledges_curr->index;
		  k = vbase[graph->tail[edge]];
		  l = vbase[graph->head[edge]];
                  node = (l == UNKNOWN)? UNKNOWN : UF_find(&uf, l);
                  adjnode = (k == UNKNOWN)? UNKNOWN : UF_find(&uf, k);

                  /* check whether the current boundary-path connects two child components */
                  if( node != UNKNOWN && nodesmark[node] && adjnode != UNKNOWN && nodesmark[adjnode] )
                     //  && !pinned[k] && !pinned[l] )
                  {
		     assert(graphmark[node]);
		     assert(graphmark[adjnode]);
                     if( debg )
		     {
                        printf("ADD horizontal edge: %d_%d  \n ", graph->tail[edge], graph->head[edge]);
			printf("ADD horizontal edge vbases: %d_%d  \n ", vbase[graph->tail[edge]], vbase[graph->head[edge]]);
			printf("ADD horizontal edge ident: %d_%d  \n ", node, adjnode);

		     }
                     boundedges[nboundedges++] = edge;
                  }
                  lvledges_curr = lvledges_curr->parent;
               }

               /* try to connect the nodes of C (directly) to COMP(C), as a preprocessing for voronoi_repair */
               count = 0;
               for( k = 0; k < nkpnodes; k++ )
               {
                  blists_curr = blists_start[kpnodes[k]];
                  assert( blists_curr != NULL );
                  while( blists_curr != NULL )
                  {
                     node = blists_curr->index;

                     /* iterate through all outgoing edges of 'node' */
                     for( edge = graph->inpbeg[node]; edge != EAT_LAST; edge = graph->ieat[edge] )
                     {
                        adjnode = graph->tail[edge];

                        /* check whether the adjacent node is not in C and allows a better voronoi assignment of the current node */
                        if( state[adjnode] == CONNECT && SCIPisGT(scip, vnoi[node].dist, vnoi[adjnode].dist + graph->cost[edge])
			   && graphmark[vbase[adjnode]] && graphmark[adjnode] )
                        {
                           vnoi[node].dist = vnoi[adjnode].dist + graph->cost[edge];
                           vbase[node] = vbase[adjnode];
                           vnoi[node].edge = edge;
                        }
                     }
                     if( vbase[node] != UNKNOWN )
                     {
                        //printf("add to heap %d \n", node );
                        heap_add(graph->path_heap, state, &count, node, vnoi);
                     }
                     blists_curr = blists_curr->parent;
                  }
               }

               /* if there are no key-path nodes, something has gone wrong */
               assert( nkpnodes != 0 );

               voronoi_repair_mult(scip, graph, graph->cost, &count, vbase, boundedges, &nboundedges, nodesmark, &uf, vnoi);

               /* create a supergraph, having the endpoints of the key-paths incident to the current crucial node as (super-) vertices */
               supergraph = graph_init(nsupernodes, nboundedges * 2, 1, 0);

               /* add vertices to the supergraph */
               for( k = 0; k < nsupernodes; k++ )
               {
                  supernodesid[supernodes[k]] = k;
                  //printf("adding node %d (org: %d) \n ", k , supernodes[k]);
                  graph_knot_add(supergraph, graph->term[supernodes[k]], 0, 0);
               }

               /* the (super-) vertex representing the current root-component of the ST */
               k = supernodes[nsupernodes - 1];

               /* add edges to the supergraph */
               for( l = 0; l < nboundedges; l++ )
               {
                  edge = boundedges[l];
                  if( debg )
                     printf("boundedgeALL: %d_%d  vbases: %d_%d \n ", graph->tail[edge], graph->head[edge],  vbase[graph->tail[edge]], vbase[graph->head[edge]]);
                  node = UF_find(&uf, vbase[graph->tail[edge]]);
                  adjnode = UF_find(&uf, vbase[graph->head[edge]]);

                  /* if node 'node' or 'adjnode' belongs to the root-component, take the (temporary) root-component identifier instead */
                  node = ((nodesmark[node])? node : k);
                  adjnode = ((nodesmark[adjnode])? adjnode : k);

                  /* compute the cost of the boundary-path pertaining to the boundary-edge 'edge' */
                  edgecost = vnoi[graph->tail[edge]].dist + graph->cost[edge] + vnoi[graph->head[edge]].dist;
                  graph_edge_add(supergraph, supernodesid[node], supernodesid[adjnode], edgecost, edgecost);
               }

               /* compute a MST on the supergraph */
               SCIP_CALL( SCIPallocBufferArray(scip, &mst, nsupernodes) );
               graph_path_init(supergraph);
               graph_path_exec(supergraph, MST_MODE, nsupernodes - 1, supergraph->cost, mst);

               /* compute the cost of the MST */
               mstcost = 0.0;
               /*
                 for( l = 0; l < nsupernodes - 1; l++ )
                 printf(" SUPERGRAPH edge: : %d -> %d \n", supergraph->tail[mst[l].edge], supergraph->head[mst[l].edge] );
               */
               /* compute the cost of the MST */
               for( l = 0; l < nsupernodes - 1; l++ )
               {
                  /* compute the edge in the original graph corresponding to the current MST edge */
                  if( mst[l].edge % 2  == 0 )
                     edge = boundedges[mst[l].edge / 2 ];
                  else
                     edge = flipedge(boundedges[mst[l].edge / 2 ]);

                  //printf(" MST egde: : %d -> %d \n", graph->tail[edge], graph->head[edge]);
                  mstcost += graph->cost[edge];
                  assert( newedges[edge] != crucnode && newedges[flipedge(edge)] != crucnode );

                  /* mark the edge (in the original graph) as visited */
                  newedges[edge] = crucnode;
                  //printf(" ADD edge: : %d -> %d \n", graph->tail[(edge)], graph->head[(edge)]);

                  /* traverse along the boundary-path belonging to the boundary-edge 'edge' */
                  for( node = graph->tail[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                  {
                     e = vnoi[node].edge;
                     //printf(" edge: : %d -> %d \n", graph->tail[(e)], graph->head[(e)]);

                     /* if edge 'e' and its reversed have not been visited yet */
                     if( newedges[e] != crucnode && newedges[flipedge(e)] != crucnode )
                     {
                        //printf(" ADD edge: : %d -> %d \n", graph->tail[(e)], graph->head[(e)]);
                        newedges[e] = crucnode;
                        mstcost += graph->cost[e];
                     }
                  }
                  for( node = graph->head[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                  {
                     e = flipedge(vnoi[node].edge);
                     //printf(" edge: : %d -> %d \n", graph->tail[e], graph->head[e]);

                     /* if edge 'e' and its reversed have not been visited yet */
                     if( newedges[vnoi[node].edge] != crucnode && newedges[e] != crucnode )
                     {
                        //printf(" ADD edge: : %d -> %d \n", graph->tail[e], graph->head[e]);
                        newedges[e] = crucnode;
                        mstcost += graph->cost[e];
                     }
                  }
               }
               //printf(" mstcost: %f \n", mstcost );
               //printf(" kpcost: %f \n", kpcost );

               if( SCIPisLT(scip, mstcost, kpcost) )
               {
                  int added = 0;
                  int removed = 0;
		  /*const char base2[] = "X/impro";
                    char filename2 [ FILENAME_MAX ];
                    sprintf(filename2, "%s%d.gml", base2, i);
                    SCIP_CALL( printGraph(scip, graph, filename2, best_result) );
                    printf("(run %d) \n ", i); */

                  localmoves++;
                  if( debg )
                     printf("found improving solution in KEY VERTEX ELIMINATION (round: %d) \n ", nruns);

                  /* unmark the original edges spanning the supergraph */
                  for( e = 0; e < nkpedges; e++ )
                  {
                     assert(best_result[kpedges[e]] != -1);
                     best_result[kpedges[e]] = -1;
                     removed += graph->cost[kpedges[e]];
                     if( debg )
		     {
                        printf(" unmark: : %d -> %d \n", graph->tail[kpedges[e]], graph->head[kpedges[e]]);
			printf(" unmarkidentif: :  %d -> %d \n", UF_find(&uf, graph->tail[kpedges[e]]), UF_find(&uf, graph->head[kpedges[e]]) );
		     }
                  }

                  /* mark all ST nodes except for those belonging to the root-component as forbidden */
                  for( k = rootpathstart; k < nkpnodes; k++ )
                  {
                     graphmark[kpnodes[k]] = FALSE;
                     steinertree[kpnodes[k]] = FALSE;
                     if( debg )
                        printf("ungraphmark(rootcomp) %d \n", kpnodes[k]);
                  }

                  for( k = 0; k < i; k++ )
                  {
		     node = UF_find(&uf, dfstree[k]);
                     if( nodesmark[node] || node == crucnode )
                     {
                        graphmark[dfstree[k]] = FALSE;
                        steinertree[dfstree[k]] = FALSE;
                        if( debg )
                           printf("ungraphmark %d \n", dfstree[k]);
                     }
                  }

                  /* add the new edges reconnecting the (super-) components */
                  for( l = 0; l < nsupernodes - 1; l++ )
                  {
                     if( mst[l].edge % 2  == 0 )
                        edge = boundedges[mst[l].edge / 2 ];
                     else
                        edge = flipedge(boundedges[mst[l].edge / 2 ]);
                     if( debg )
                        printf("MST edge vbase tail %d vbase head: %d \n",vbase[graph->tail[edge]],  vbase[graph->head[edge]] );

                     /* change the orientation within the target-component if necessary */
                     if( !nodesmark[vbase[graph->head[edge]]] )
                     {
                        node = vbase[graph->head[edge]];
                        k = UF_find(&uf, node);
                        assert(nodesmark[k]);
                        while( node != k )
                        {
                           /* the ST edge pointing towards the root */
                           e = nodes[node].edge;

                           assert(best_result[e] == -1 && best_result[flipedge(e)] != -1 );
                           if( debg )
                              printf(" switch : %d->%d \n ", graph->tail[e], graph->head[e]);
                           best_result[e] = CONNECT;
                           best_result[flipedge(e)] = UNKNOWN;
                           node = graph->head[e];
                        }
                     }

                     /* is the vbase of the current boundary-edge tail in the root-component? */
                     if( !nodesmark[UF_find(&uf, vbase[graph->tail[edge]])] )
                     {
                        if( debg )
                           printf(" FINAL ADD root edgee: : %d -> %d \n", graph->tail[edge], graph->head[edge]);
                        //                     assert( best_result[edge] != 0 && best_result[flipedge(edge)] != 0 );
                        best_result[edge] = CONNECT;
                        added += graph->cost[edge];

                        if ( !graphmark[vbase[graph->tail[edge]]])
                        {
                           const char base[] = "X/debug";
                           char filename [ FILENAME_MAX ];
                           sprintf(filename, "%s%d.gml", base, i);
                           SCIP_CALL( printGraph(scip, graph, filename, best_result) );
                           printf("nodenumber: %d \n", vbase[graph->tail[edge]] );
                           printf("nodenumberidentifier: %d \n", UF_find(&uf,  vbase[graph->tail[edge]] ) );
			   if( pinned[vbase[graph->tail[edge]]] )
                              printf("vbase pinned \n");
                           else
                              printf("vbase not pinned \n");
                           assert(0);
                        }

                        for( node = graph->tail[edge], adjnode = graph->head[edge]; node != vbase[node]; adjnode = node, node = graph->tail[vnoi[node].edge] )
                        {
                           graphmark[node] = FALSE;
			   if( debg )
                              printf("ungraphmark %d \n", node);
                           if( best_result[flipedge(vnoi[node].edge)] == CONNECT )
                           {
			      assert(0); /* should never happen (?) */

                              best_result[flipedge(vnoi[node].edge)] = UNKNOWN;
                              removed += graph->cost[flipedge(vnoi[node].edge)];
                              if( debg )
                                 printf(" FINAL delete reverse1 of : %d -> %d \n", graph->tail[(vnoi[node].edge)], graph->head[(vnoi[node].edge)]);
                           }
                           if( debg )
                              printf("FINAL ADD rootedge: : %d -> %d \n", graph->tail[(vnoi[node].edge)], graph->head[(vnoi[node].edge)]);
                           best_result[vnoi[node].edge] = CONNECT;
                           added += graph->cost[vnoi[node].edge];
                        }

                        assert(!nodesmark[node] && vbase[node] == node);
                        //assert( graph->tail[(vnoi[adjnode].edge)] == node );
                        assert( graphmark[node] == TRUE );

                        /* is the pinned node its own component identifier? */
                        if( !Is_term(graph->term[node]) && !pinned[node] && !scanned[node] )
                        {
                           oldedge = edge;
			   if( debg )
                              printf("A DAY IN THE LIFE \n \n");
                           /* move down the ST until a crucial or pinned node has been reached
                            * and check whether it has already been scanned */
                           for( edge = graph->outbeg[node]; edge != EAT_LAST; edge = graph->oeat[edge] )
                           {

                              /* check whether edge 'edge' leads to an ancestor of node 'node' */
                              if( best_result[edge] == CONNECT && steinertree[adjnode = graph->head[edge]] )
                              {
                                 /* move along the key-path until it ends (i.e. until a crucial or pinned node has been reached) */
                                 while( !nodeIsCrucial(graph, best_result, adjnode) && !pinned[adjnode] )
                                 {
                                    for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                                    {
                                       if( best_result[e] == CONNECT )
                                          break;
                                    }
                                    adjnode = graph->head[e];
                                    if( !steinertree[adjnode] )
                                    {
                                       adjnode = graph->tail[e];
                                       break;
                                    }
                                 }
                                 assert( adjnode != node );

                                 /*  has node 'adjnode' already been scanned? */
                                 if( scanned[adjnode] )
                                 {
                                    /* move up again, comcomitantly updating the data structures */
                                    while( adjnode != node )
                                    {
                                       phnode_meldheaps(scip, &boundpaths[node], &boundpaths[adjnode], &heapsize[node], &heapsize[adjnode]);
                                       if( debg )
                                          printf( "unite 3 (%d) (%d) \n ",  node, adjnode);
                                       uf_union(&uf, node, adjnode, FALSE);
                                       adjnode = graph->head[nodes[adjnode].edge];
                                    }
                                 }
                              }
                           }
#if 0
                           /* update union-find and pairing heaps: unite terminal 'node' with all of its ancestor key-paths */
                           for( edge = graph->outbeg[node]; edge != EAT_LAST; edge = graph->oeat[edge] )
                           {

                              /* check whether edge 'edge' leads to an ancestor of terminal 'node' */
                              if( best_result[edge] != -1 && steinertree[graph->head[edge]] )
                              {
                                 adjnode = graph->head[edge];

                                 /* meld the heaps */
                                 phnode_meldheaps(scip, &boundpaths[node], &boundpaths[adjnode], &heapsize[node], &heapsize[adjnode]);
                                 if( debg )
                                    printf( "unite 3 (%d) (%d) \n ",  node, adjnode);
                                 /* update the union-find data structure */
                                 uf_union(&uf, node, adjnode, FALSE);

                                 //printf("unite1 %d + %d \n", node, adjnode);


                                 /* move along the key-path until its end (i.e. until a crucial node is reached) */
                                 while( !nodeIsCrucial(graph, best_result, adjnode) && !pinned[adjnode] )
                                 {
                                    for( e = graph->outbeg[adjnode]; e != EAT_LAST; e = graph->oeat[e] )
                                    {
                                       if( best_result[e] != -1 )
                                          break;
                                    }

                                    /* assert that each leaf of the ST is a terminal */
                                    assert( e != EAT_LAST );
                                    adjnode = graph->head[e];
                                    if( !steinertree[adjnode] )
                                       break;
                                    if( debg )
                                       printf( "unite 4 (%d) (%d) \n ",  node, adjnode);
                                    /* update the union-find data structure */
                                    uf_union(&uf, node, adjnode, FALSE);

                                    //printf("unite1x %d + %d \n", node, adjnode);

                                    /* meld the heaps */
                                    phnode_meldheaps(scip, &boundpaths[node], &boundpaths[adjnode], &heapsize[node], &heapsize[adjnode]);
                                 }
                              }
                           }
#endif
                           edge = oldedge;
                        }


                        /* mark the start node (lying in the root-component of the ST) of the current boundary-path as pinned,
                         * so that it may not be removed later on */
                        pinned[node] = TRUE;
                        if( debg )
                           printf("pinned node: %d \n", node);

                        for( node = graph->head[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                        {
                           //printf("ungraphmark %d \n", node);
                           graphmark[node] = FALSE;
                           if( best_result[vnoi[node].edge] == CONNECT )
                           {
                              best_result[vnoi[node].edge] = -1;
                              removed += graph->cost[vnoi[node].edge];
                              if( debg )
                                 printf(" FINAL delete reverse2 of : %d -> %d \n", graph->head[(vnoi[node].edge)], graph->tail[(vnoi[node].edge)]);
                           }
                           if( debg )
                              printf("FINAL ADD rootedge: : %d -> %d \n", graph->tail[flipedge(vnoi[node].edge)], graph->head[flipedge(vnoi[node].edge)]);
                           best_result[flipedge(vnoi[node].edge)] = CONNECT;
                           added += graph->cost[flipedge(vnoi[node].edge)];

                        }
                     }
                     else
                     {
                        if( debg )
                           printf(" FINAL ADD egde: : %d -> %d \n", graph->tail[edge], graph->head[edge]);

                        //                     assert( best_result[edge] != 0 && best_result[flipedge(edge)] != 0 );
                        best_result[edge] = CONNECT;
                        added += graph->cost[edge];
                        for( node = graph->tail[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                        {
                           graphmark[node] = FALSE;
                           //printf("ungraphmark %d \n", node);
                           //printf(" FINAL edge: : %d -> %d \n", graph->tail[(vnoi[node].edge)], graph->head[(vnoi[node].edge)]);
                           if( best_result[vnoi[node].edge] != CONNECT && best_result[flipedge(vnoi[node].edge)] != CONNECT )
                           {
                              if( debg )
                                 printf("FINAL ADD edge: : %d -> %d \n", graph->tail[(vnoi[node].edge)], graph->head[(vnoi[node].edge)]);
                              best_result[vnoi[node].edge] = CONNECT;
                              added+= graph->cost[(vnoi[node].edge)];
                           }
                        }

                        for( node = graph->head[edge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                        {
                           graphmark[node] = FALSE;
                           //printf("FINAL edge: : %d -> %d \n", graph->tail[flipedge(vnoi[node].edge)], graph->head[flipedge(vnoi[node].edge)]);
                           if( 1 )//if( best_result[vnoi[node].edge] != CONNECT && best_result[flipedge(vnoi[node].edge)] != CONNECT )
                           {
                              if( debg )
                                 printf("FINAL ADD edge: : %d -> %d \n", graph->tail[flipedge(vnoi[node].edge)], graph->head[flipedge(vnoi[node].edge)]);
                              best_result[flipedge(vnoi[node].edge)] = CONNECT;
			      best_result[vnoi[node].edge] = UNKNOWN;
                              added += graph->cost[flipedge(vnoi[node].edge)];
                           }
                        }
                     }
                  }

                  for( k = 0; k < nkpnodes; k++ )
                  {
                     if (graphmark[kpnodes[k]] != FALSE)
		     {

		        const char base[] = "X/debugMark";
                        char filename [ FILENAME_MAX ];
                        sprintf(filename, "%s%d.gml", base, i);
                        SCIP_CALL( printGraph(scip, graph, filename, best_result) );
                        printf(" node: %d \n", kpnodes[k]);
                        assert(0);
		     }
                     assert(steinertree[kpnodes[k]] == FALSE);
                     //printf("ungraphmark %d \n", kpnodes[k]);
                  }
                  assert(!graphmark[crucnode]);
                  //printf(" added : %d \n", added);
                  //printf("deleted %d \n", removed);
               }
               else
               {
                  /* no improving solution has been found during the move */

                  /* meld the heap pertaining to 'crucnode' and all heaps pertaining to descendant key-paths of node 'crucnode' */
                  for( k = 0; k < rootpathstart; k++ )
                  {
                     phnode_meldheaps(scip, &boundpaths[crucnode], &boundpaths[kpnodes[k]], &heapsize[crucnode], &heapsize[kpnodes[k]]);
                  }
                  for( k = 0; k < nsupernodes - 1; k++ )
                  {
                     phnode_meldheaps(scip, &boundpaths[crucnode], &boundpaths[supernodes[k]], &heapsize[crucnode], &heapsize[supernodes[k]]);
                     if( debg )
                        printf( "unite 5 (%d) (%d) \n ",  crucnode, supernodes[k]);
                     /* update the union-find data structure */
                     uf_union(&uf, crucnode, supernodes[k], FALSE);
                  }
               }

               /* free the supergraph and the MST data structure */
               graph_path_exit(supergraph);
               graph_free(supergraph);
               SCIPfreeBufferArray(scip, &mst);

               /* unmark the descendant supervertices */
               for( k = 0; k < nsupernodes - 1; k++ )
               {
                  nodesmark[supernodes[k]] = FALSE;
               }

               /* debug test; to be deleted later on */
               for( k = 0; k < nnodes; k++ )
               {
                  assert( !nodesmark[k] );
               }

               /* restore the original voronoi diagram */
               l = 0;
               for( k = 0; k < nkpnodes; k++ )
               {
                  /* restore data of all nodes having the current (internal) key-path node as their voronoi base */
                  blists_curr = blists_start[kpnodes[k]];
                  while( blists_curr != NULL )
                  {
                     node = blists_curr->index;
		     /*if( !graphmark[node] )
                       {
                       TODO? dont reset then?
                       }*/
                     vbase[node] = memvbase[l];
                     vnoi[node].dist = memdist[l];
                     vnoi[node].edge = meminedges[l];
                     l++;
                     blists_curr = blists_curr->parent;
                  }
               }

               /* debug TODO delete*/
               assert(l == nresnodes);
	    }

	    /** Key-Path Exchange */
	    if( 1 )
	    {
               //printf("ST NODE: %d\n", crucnode);
               /* if the node has just been eliminated, skip Key-Path Exchange */
               if( !graphmark[crucnode] )
                  continue;

               /* is crucnode a crucial or pinned vertex? */
               if( (!nodeIsCrucial(graph, best_result, crucnode) && !pinned[crucnode]) )
                  continue;

               /* counts the internal nodes of the keypath */
               nkpnodes = 0;

               for( k = 0; k < nnodes; k++ )
               {
                  assert( state[k] == CONNECT);
               }

               /* find the (unique) key-path containing the parent of the current crucial node 'crucnode' */
               kptailnode = graph->head[nodes[crucnode].edge];
               kpathcost = graph->cost[nodes[crucnode].edge];
               if( debg )
                  printf("kpathhead: %d \n " ,crucnode);

               while( !nodeIsCrucial(graph, best_result, kptailnode) && !pinned[kptailnode] )
               {
                  kpathcost += graph->cost[nodes[kptailnode].edge];
		  if( debg )
                     printf("kpathinternal: %d \n " , kptailnode);
                  kpnodes[nkpnodes++] = kptailnode;
                  kptailnode = graph->head[nodes[kptailnode].edge];
               }
               if( debg )
                  printf("kpathtail: %d \n " , kptailnode);

               /* counts the reset nodes during voronoi repair */
               nresnodes = 0;

               /* reset all nodes (henceforth referred to as 'C') whose bases are internal nodes of the current keypath */
               for( k = 0; k < nkpnodes; k++ )
               {
                  /* reset all nodes having the current (internal) keypath node as their voronoi base */
                  blists_curr = blists_start[kpnodes[k]];
                  while( blists_curr != NULL )
                  {
                     node = blists_curr->index;
                     //printf("C-node %d \n", blists_curr->index);
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

                  /* update the union-find data structure
                     if( k != 0 )
                     {
                     uf_union(&uf, kpnodes[0], kpnodes[k], FALSE);
                     } */
               }

               edgecost = UNKNOWN;
               e = UNKNOWN;
               while( boundpaths[crucnode] != NULL )
               {
                  phnode_deletemin(scip, &e, &edgecost, &boundpaths[crucnode], &(heapsize[crucnode]));
                  assert( e != UNKNOWN );
		  k = vbase[graph->tail[e]];
                  l = vbase[graph->head[e]];
                  if( !graphmark[k] )
                  {
                     assert(graphmark[graph->tail[e]]);
                     //printf(" unmarked: %d \n", graph->tail[e]);
                     //printf(" unmarkedhead: %d \n", graph->head[e]);
                     //printf(" vbase unmarked: %d \n", k);
                  }
                  assert(graphmark[k]);
                  node = (l == UNKNOWN || !graphmark[l] )? UNKNOWN : UF_find(&uf, l);
                  adjnode = (k == UNKNOWN)? UNKNOWN : UF_find(&uf, k);

                  if ( 0 && e != -1 &&  debg )
                  {
                     printf("prenodes %d_%d  \n ", graph->head[e], graph->tail[e]);
                     printf("basenodes:  %d_%d\n ", node, adjnode);
                  }
                  if(  adjnode != crucnode && graphmark[adjnode] )
                  {
                     const char base[] = "X/debugX";
                     char filename [ FILENAME_MAX ];
                     sprintf(filename, "%s%d.gml", base, i);
                     printf( "adjnode: %d \n ", adjnode);
                     printf("vnoi %d_%d  \n ", k, l);
                     printf("vnoi %d_%d  \n ", UF_find(&uf,  k), UF_find( &uf, l) );

                     SCIP_CALL( printGraph(scip, graph, filename, best_result) );
                     printf("nodenumber: %d \n", vbase[graph->tail[edge]] );
                     printf("nodenumberidentifier: %d \n", UF_find(&uf,  vbase[graph->tail[edge]] ) );

                     assert(0);
                  }
                  assert( graphmark[adjnode] );

                  /* does the boundary-path end in the root component? */
                  if( node != UNKNOWN && node != crucnode && graphmark[l] ) //&& !pinned[k] && !pinned[l] )
                  {
                     if( debg )
                        printf("add exg vbase : %d %d \n", k,  l);
                     phnode_insert(scip, &boundpaths[crucnode], e, edgecost, &(heapsize[crucnode]));
                     break;
                  }
               }

               if( boundpaths[crucnode] == NULL )
               {
                  oldedge = UNKNOWN;
               }
               else
               {
                  oldedge = e;
               }

               /* counts the nodes connected during the following 'preprocessing' */
               count = 0;
               /* try to connect the nodes of C (directly) to COMP(C), as a preprocessing for voronoi-repair */
               for( k = 0; k < nkpnodes; k++ )
               {
                  blists_curr = blists_start[kpnodes[k]];
                  assert( blists_curr != NULL );
                  while( blists_curr != NULL )
                  {
                     node = blists_curr->index;

                     /* iterate through all outgoing edges of 'node' */
                     for( edge = graph->inpbeg[node]; edge != EAT_LAST; edge = graph->ieat[edge] )
                     {
                        adjnode = graph->tail[edge];

                        /* check whether the adjacent node is not in C and allows a better voronoi assignment of the current node */
                        if( state[adjnode] == CONNECT && SCIPisGT(scip, vnoi[node].dist, vnoi[adjnode].dist + graph->cost[edge])
			   && graphmark[vbase[adjnode]] && graphmark[adjnode] )
                        {
                           vnoi[node].dist = vnoi[adjnode].dist + graph->cost[edge];
                           vbase[node] = vbase[adjnode];
                           vnoi[node].edge = edge;
                        }
                     }
                     if( vbase[node] != UNKNOWN )
                     {
                        //printf("add to heap %d \n", node );
                        heap_add(graph->path_heap, state, &count, node, vnoi);
                     }
                     blists_curr = blists_curr->parent;
                  }
               }
               if( nkpnodes > 0 )
                  assert(count > 0);
               newedge = UNKNOWN;

               /* if there is no key path, nothing has to be repaired */
               if( nkpnodes > 0 )
                  voronoi_repair(scip, graph, graph->cost, &count, vbase, vnoi, &newedge, crucnode, &uf);
               else
                  newedge = nodes[crucnode].edge;
               if( 0 && debg ){
                  for(  e = 0; e < nnodes; e++)
                     printf("(completevoronoi)base[%d] = %d \n", e, vbase[e]);
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

               //   printf("KOSTENVERGLEICH old/new: %f_ %f \n ", kpathcost,  vnoi[graph->tail[newedge]].dist + graph->cost[newedge]
               //    + vnoi[graph->head[newedge]].dist );
               //printf("final edge %d_%d \n ", graph->tail[newedge], graph->head[newedge]);
	       if( debg )
                  printf("final edge vronoi  %d_%d \n ", vbase[graph->tail[newedge]], vbase[graph->head[newedge]]);
               assert( newedge != UNKNOWN );
               edgecost = vnoi[graph->tail[newedge]].dist + graph->cost[newedge] + vnoi[graph->head[newedge]].dist;
               if( SCIPisLT(scip, edgecost, kpathcost) )
               {

                  bestdiff = edgecost - kpathcost;
                  node = UF_find(&uf, vbase[graph->head[newedge]]);
                  obj += bestdiff;
                  /*   const char bas[] = "X/exchange";
                       char filenam [ FILENAME_MAX ];
                       sprintf(filenam, "%s%d.gml", bas, i);

                       SCIP_CALL( printGraph(scip, graph, filenam, best_result) ); */
		  if( debg )
                     printf( "ADDING NEW KEY PATH (%f )\n", bestdiff );
                  localmoves++;

                  /* remove old keypath */
                  assert(  best_result[flipedge(nodes[crucnode].edge)] != UNKNOWN );
                  best_result[flipedge(nodes[crucnode].edge)] = UNKNOWN;
                  steinertree[crucnode] = FALSE;
                  graphmark[crucnode] = FALSE;
		  if( debg )
                     printf("unmarkcruc %d \n", crucnode);

		  if( debg )
                     printf("delete: %d->%d \n", graph->tail[ flipedge(nodes[crucnode].edge) ], graph->head[ flipedge(nodes[crucnode].edge) ]);
                  for( k = 0; k < nkpnodes; k++ )
                  {
                     assert(  best_result[flipedge(nodes[kpnodes[k]].edge)] != UNKNOWN );
                     best_result[flipedge(nodes[kpnodes[k]].edge)] = UNKNOWN;
                     steinertree[kpnodes[k]] = FALSE;
                     graphmark[kpnodes[k]] = FALSE;
                     if( debg )
                        printf("unmarkkp %d \n", kpnodes[k]);
		     if( debg )
                        printf("delete: %d->%d \n", graph->tail[ flipedge(nodes[kpnodes[k]].edge) ], graph->head[ flipedge(nodes[kpnodes[k]].edge)]);
                  }
                  assert(graphmark[kptailnode]);

                  /* add new key-path */
		  /*
                    k = (node == crucnode)? graph->head[newedge] : graph->tail[newedge];
                    while( k != vbase[k] )
                    {
                    best_result[flipedge(vnoi[k].edge)] = CONNECT;
                    if( debg )
                    printf("add %d->%d \n", graph->tail[ flipedge(vnoi[k].edge) ], graph->head[ flipedge(vnoi[k].edge) ]);
                    k = graph->tail[vnoi[k].edge];
                    }

                    k = (node == crucnode)? graph->tail[newedge] : graph->head[newedge];

                    while( k != vbase[k] )
                    {
                    best_result[vnoi[k].edge] = CONNECT;
                    if( debg )
                    printf("add %d->%d \n", graph->tail[ (vnoi[k].edge) ], graph->head[ (vnoi[k].edge) ]);
                    k = graph->tail[vnoi[k].edge];

                    }
                    best_result[(node == crucnode)? newedge : flipedge(newedge)] = CONNECT;
		  */
		  if( node == crucnode )
		  {
		     if( debg )
                        printf("whoaa \n \n");
		     newedge = flipedge(newedge);
		  }
		  if( debg )
                     printf("vbases newedge %d %d \n", vbase[graph->tail[newedge]], vbase[graph->head[newedge]] );
                  for( node = graph->tail[newedge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                  {
                     //printf("ungraphmark %d \n", node);
                     if( debg )
                        printf("unmarknew %d \n", node);
                     graphmark[node] = FALSE;

                     best_result[flipedge(vnoi[node].edge)] = CONNECT;
                     best_result[vnoi[node].edge] = UNKNOWN;
                     if( debg ){
                        printf("add(Tail) %d->%d \n", graph->tail[ flipedge(vnoi[node].edge) ], graph->head[ flipedge(vnoi[node].edge) ]);
                        printf("(->X)vbase %d  \n", vbase[graph->head[ flipedge(vnoi[node].edge)] ]);
                     }
                  }

		  for( node = graph->head[newedge]; node != vbase[node]; node = graph->tail[vnoi[node].edge] )
                  {
                     //printf("ungraphmark %d \n", node);
                     if( debg )
                        printf("unmarknew %d \n", node);
                     graphmark[node] = FALSE;

                     best_result[vnoi[node].edge] = CONNECT;
                     if( debg )
                        printf("add(head) %d->%d \n", graph->tail[ (vnoi[node].edge) ], graph->head[ (vnoi[node].edge) ]);
                  }


		  if( debg )
                     printf("add %d->%d \n", graph->tail[ (node == crucnode)? newedge : flipedge(newedge) ], graph->head[ (node == crucnode)? newedge : flipedge(newedge) ]);
		  best_result[flipedge(newedge)] = CONNECT;
                  //printf("bestpathI: %d->%d \n " , graph->tail[newpath_curr->index], graph->head[newpath_curr->index]);
                  /* TODO could be solved in a better way */
                  // pinned[kptailnode] = TRUE;

		  newpathend = vbase[graph->tail[newedge]];
		  assert(node == vbase[graph->head[newedge]] );
		  pinned[node] = TRUE;
                  /* if( node == crucnode )
                     {
                     newpathend = vbase[graph->head[newedge]];
                     pinned[vbase[graph->tail[newedge]]] = TRUE;
                     }
                     else
                     {
                     newpathend = vbase[graph->tail[newedge]];
                     pinned[vbase[graph->head[newedge]]] = TRUE;
                     }
                  */
                  /* flip all edges on the ST path between the endnode of the new key-path and the current crucial node */
                  k = newpathend;
                  //printf(" root: %d \n ", graph->source[0]);
                  if( UF_find(&uf, newpathend) != crucnode )
                  {
                     printf(" newpath: %d crucnode: %d \n ", newpathend, crucnode);
                     assert(0);
                  }
                  while( k != crucnode )
                  {
                     //printf("k %d, \n", k);
                     assert(graphmark[k]);
                     assert( best_result[flipedge(nodes[k].edge)] != -1);
                     best_result[flipedge(nodes[k].edge)] = UNKNOWN;

                     best_result[nodes[k].edge] = CONNECT;
		     if( debg )
                        printf("flipedge:  %d->%d \n", graph->tail[nodes[k].edge ], graph->head[nodes[k].edge ]);
                     k = graph->head[nodes[k].edge];
                  }


                  for( k = 0; k < i; k++ )
                  {
                     if( crucnode == UF_find(&uf, dfstree[k]) )
                     {
                        graphmark[dfstree[k]] = FALSE;
			steinertree[dfstree[k]] = FALSE;
			if( debg )
                           printf("unmarkEx %d \n", dfstree[k]);
                        //printf("graphmark node CHILDREDN: %d = FALSE; \n", dfstree[k]);
                     }
                  }

                  /*   const char base[] = "X/KEYPATHEXCHG";
                       char filename [ FILENAME_MAX ];

                       sprintf(filename, "%s%d%s", base, i, ".gml");

                       SCIP_CALL( printGraph(scip, graph, filename, best_result) );*/
               }
               else
               {
                  if( Is_term(graph->term[kptailnode]) || pinned[kptailnode] )
                  {
                     /* update union-find data structure
                        if( nkpnodes > 0 )
                        {
                        uf_union(&uf, crucnode, kpnodes[0], FALSE);
                        } */

                     /* merge the heaps pertaining to the current key-path */
                     for( k = 0; k < nkpnodes - 1; k++ )
                     {
                        phnode_meldheaps(scip, &boundpaths[kpnodes[k + 1]], &boundpaths[kpnodes[k]], &heapsize[kpnodes[k + 1]], &heapsize[kpnodes[k]]);
                        //printf("1merge %d, %d \n", kpnodes[k + 1], kpnodes[k ]);
                        /*boundpaths[kpnodes[k + 1]] = phnode_mergeheaps(scip, boundpaths[kpnodes[k + 1]], boundpaths[kpnodes[k]]);
                          heapsize[kpnodes[k + 1]] += heapsize[kpnodes[k]];
                          boundpaths[kpnodes[k]] = NULL;*/

                        /* update union-find data structure */
                        uf_union(&uf, crucnode, kpnodes[k], FALSE);
			if( debg )
			   printf("uniteA, %d, %d \n", crucnode, kpnodes[k]);
                     }

                     if( nkpnodes > 0 )
                     {
                        phnode_meldheaps(scip, &boundpaths[kptailnode], &boundpaths[kpnodes[nkpnodes - 1]], &heapsize[kptailnode], &heapsize[kpnodes[nkpnodes - 1]]);
                        /*boundpaths[kptailnode] = phnode_mergeheaps(scip, boundpaths[kptailnode], boundpaths[kpnodes[nkpnodes - 1]]);
                          heapsize[kptailnode] += heapsize[kpnodes[nkpnodes - 1]];
                          boundpaths[kpnodes[nkpnodes - 1]] = NULL;*/

                        uf_union(&uf, crucnode, kpnodes[nkpnodes - 1], FALSE);
                        //printf("2merge %d, %d \n", kptailnode, kpnodes[nkpnodes - 1]);
                     }

                     phnode_meldheaps(scip, &boundpaths[kptailnode], &boundpaths[crucnode], &heapsize[kptailnode], &heapsize[crucnode]);
                     /*boundpaths[kptailnode] = phnode_mergeheaps(scip, boundpaths[kptailnode], boundpaths[crucnode]);
                       heapsize[kptailnode] += heapsize[crucnode];
                       boundpaths[crucnode] = NULL;*/

                     uf_union(&uf, kptailnode, crucnode, FALSE);
		     if( debg )
                        printf("uniteB, %d, %d \n", kptailnode, crucnode);
                  }
                  //	    printf("3merge %d, %d \n", kptailnode, crucnode);
               }

               /* restore the original voronoi digram */
               l = 0;
               for( k = 0; k < nkpnodes; k++ )
               {
                  /* reset all nodes having the current (internal) keypath node as their voronoi base */
                  blists_curr = blists_start[kpnodes[k]];
                  while( blists_curr != NULL )
                  {
                     node = blists_curr->index;
                     vbase[node] = memvbase[l];
                     vnoi[node].dist = memdist[l];
                     vnoi[node].edge = meminedges[l];
                     l++;
                     blists_curr = blists_curr->parent;
                  }
               }
               assert(l == nresnodes);
            }
         }
#if 0
         /* debug! */
         int* xvbase;
         PATH* xvnoi;
         SCIP_CALL( SCIPallocBufferArray(scip, &xvbase, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &xvnoi, nnodes) );

         voronoi(graph, cost, steinertree, xvbase, xvnoi);

         for( e = 0; e < nnodes; e++ )
         {
            assert(vbase[e] == xvbase[e] && vnoi[e].dist == xvnoi[e].dist && vnoi[e].edge == xvnoi[e].edge);
         }

         /* debug! */
         SCIPfreeBufferArray(scip, &xvbase);
         SCIPfreeBufferArray(scip, &xvnoi);
         SCIP_Bool* edgemark;
         SCIP_CALL( SCIPallocBufferArray(scip, &edgemark, nedges / 2) );
         for( e = 0; e < nedges / 2; e++ ){
            if( best_result[2*e] == 0 || best_result[flipedge(2*e)] == 0)
               edgemark[e] = TRUE;
            else
               edgemark[e] = FALSE;
         }

         SCIP_CALL( SCIPprobdataPrintGraph2(graph,"TESTXX.gml", edgemark) );
         SCIPfreeBufferArray(scip, &edgemark);
         assert(0);
#endif

#if 0
         SCIP_Bool* edgemark;
         SCIP_CALL( SCIPallocBufferArray(scip, &edgemark, nedges / 2) );
         for( e = 0; e < nedges / 2; e++ ){
            edgemark[e] = FALSE;
         }
#endif

         /*
           const char base[] = "X/GRAPHLL";
           char filename [ FILENAME_MAX ];

           sprintf(filename, "%s%d.gml", base, i);


           SCIP_CALL( printGraph(scip, graph, filename, best_result) );


           for( e = 0; e < nedges; e++ )
           {
           if( best_result[e] == CONNECT )
           printf(" EDGE: %d->%d \n", graph->tail[e], graph->head[e]);
           }*/

         /* SCIP_CALL( SCIPprobdataPrintGraph2(graph,"TESTG.gml", edgemark) );
            SCIPfreeBufferArray(scip, &edgemark);*/


         /**********************************************************/

         /* free data structures */
         uf_free(scip, &uf);
         SCIPfreeBufferArray(scip, &supernodes);
         SCIPfreeBufferArray(scip, &kpedges);
         SCIPfreeBufferArray(scip, &kpnodes);

         for( k = 0; k < nnodes; k++ )
         {
            if( boundpaths[k] != NULL )
            {
               phnode_free(scip, &boundpaths[k]);
            }

            blists_curr = blists_start[k];
            lvledges_curr = lvledges_start[k];
            while( lvledges_curr != NULL )
            {
               lvledges_start[k] = lvledges_curr->parent;
               SCIPfreeMemory(scip, &lvledges_curr);
               lvledges_curr = lvledges_start[k];
            }

            while( blists_curr != NULL )
            {
               blists_start[k] = blists_curr->parent;
               SCIPfreeMemory(scip, &blists_curr);
               blists_curr = blists_start[k];
            }
         }

	 if( localmoves > 0 )
	 {

            /*    const char base[] = "X/chged";
                  char filename [ FILENAME_MAX ];

                  sprintf(filename, "%s%d.gml", base, i);


                  SCIP_CALL( printGraph(scip, graph, filename, best_result) ); */

            for( i = 0; i < nnodes; i++ )
            {
               steinertree[i] = FALSE;
               init(&nodes[i]);
            }

            /* create a link-cut tree representing the current Steiner tree */
            for( e = 0; e < nedges; e++ )
            {
               assert(graph->head[e] == graph->tail[flipedge(e)]);

               /* if edge e is in the tree, so are its incident vertices */
               if( best_result[e] != -1 )
               {
                  steinertree[graph->tail[e]] = TRUE;
                  steinertree[graph->head[e]] = TRUE;
                  link(&nodes[graph->head[e]], &nodes[graph->tail[e]], flipedge(e));
               }
            }
            assert( nodes[root].edge == -1 );
            nodes[root].edge = -1;


	 }

      }

      /* free data structures */
      SCIPfreeBufferArray(scip, &vnoi);
      SCIPfreeBufferArray(scip, &dfstree);
      SCIPfreeBufferArray(scip, &supernodesid);
      SCIPfreeBufferArray(scip, &scanned);
      SCIPfreeBufferArray(scip, &heapsize);
      SCIPfreeBufferArray(scip, &boundedges);
      SCIPfreeBufferArray(scip, &newedges);

      SCIPfreeBufferArray(scip, &vbase);
      SCIPfreeBufferArray(scip, &memvbase);
      SCIPfreeBufferArray(scip, &memdist);
      SCIPfreeBufferArray(scip, &meminedges);
      SCIPfreeBufferArray(scip, &nodesmark);
      SCIPfreeBufferArray(scip, &pinned);


      SCIPfreeBufferArray(scip, &lvledges_start);
      SCIPfreeBufferArray(scip, &boundpaths);
      SCIPfreeBufferArray(scip, &blists_start);

      /******/
   }

   obj = 0.0;
   for( e = 0; e < nedges; e++)
      obj += (best_result[e] > -1) ? graph->cost[e] : 0.0;
   if( 0 )
      printf(" ObjAfterLocal=%.12e\n", obj);

   SCIPfreeBufferArray(scip, &nodes);
   SCIPfreeBufferArray(scip, &steinertree);
   return SCIP_OKAY;
}

static
SCIP_RETCODE do_layer(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*  graph,
   int           layer,
   int*          best_result,
   int           runs,
   const SCIP_Real* cost,
   const SCIP_Real* costrev
   )
{
   PATH** path;
   char* connected;
   int* result;
   int* start;
   SCIP_Real obj;
   SCIP_Real objt;
   SCIP_Real min = FARAWAY;
   int best = -1;
   int k;
   int r;
   int e;
   int nnodes;
   int nedges;
   int mode;

   GNODE** gnodearr;
   int* nodenterms;
   int nterms;
   int** basearr;
   SCIP_Real** distarr;
   int** edgearr;
   int* vcount;
   SCIP_PQUEUE* pqueue;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(best_result != NULL);
   assert(cost != NULL);
   assert(costrev != NULL);
   assert(layer >= 0 && layer < graph->layers);

   nnodes = graph->knots;
   nedges = graph->edges;
   nterms = graph->terms;

   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &start, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );

   /* get user parameter */
   SCIP_CALL( SCIPgetIntParam(scip, "stp/tmheuristic", &mode) );
   if( 1 )
      printf(" tmmode: %d ->", mode);
   assert(mode == AUTO || mode == TM || mode == TMPOLZIN);

   if( mode == AUTO )
   {
      /* are there enough terminals for the TM Polzin variant to (expectably) be advantageous? */
      if( SCIPisGE(scip, ((double) nterms) / ((double) nnodes ), 0.1) )
         mode = TMPOLZIN;
      else
         mode = TM;
   }
   if( 1 )
      printf(" %d \n", mode);
   /* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    * Patch um die heuristic nach einem restruct starten zu koennen
    * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    */
   if( best >= nnodes )
      best = -1;

   if( graph->layers > 1 )
   {
      /*  SCIP_CALL( do_heuristic(scip, graph, layer, best_result, graph->source[layer], connected, cost, costrev, path) );*/
      assert(0);
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

         /* do we need to set the start by hand? */
         if( r == runs )
            start[0] = best;
      }
      else
      {
         runs = nnodes;
      }

      if( mode == TM )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );
         BMSclearMemoryArray(path, nnodes);

         /* for( k = 0; k < nnodes; k++ )
            path[k] = NULL;  TODO why???*/
      }
      else
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &nodenterms, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &basearr, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &distarr, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &edgearr, nnodes) );

         for( k = 0; k < nnodes; k++ )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &basearr[k], nterms) );
            SCIP_CALL( SCIPallocBufferArray(scip, &distarr[k], nterms) );
            SCIP_CALL( SCIPallocBufferArray(scip, &edgearr[k], nterms) );
         }

         SCIP_CALL( SCIPallocBufferArray(scip, &vcount, nnodes) );
         SCIP_CALL( SCIPpqueueCreate( &pqueue, nnodes, 2, GNODECmpByDist) );
      }
      /*graph_path_exit(graph);
        graph_path_init(graph); */

      for( r = 0; r < runs; r++ )
      {
         /* incorrect if layers > 1 ! */
         assert(graph->layers == 1);

         for( e = 0; e < nedges; e++ )
            result[e] = -1;

         /* SCIP_CALL( do_heuristic(scip, graph, layer, result, start[r], connected, cost, costrev, path, nodeterms, nodenterms, ncheckedterms, pqueue, (r == 0) ,
            (r == runs - 1)) );
            SCIP_CALL( do_heuristic(scip, graph, pqueue, path, gnodearr, distarr, cost, costrev, layer, start[r], result,
            vcount, nodesid, nodenterms, basearr, edgearr, (r == 0), (r == runs - 1), connected) );*/


	 if( mode == TM )
            SCIP_CALL( do_tm(scip, graph, path, cost, costrev, layer, start[r], result, connected) );
	 else
            SCIP_CALL( do_tm_polzin(scip, graph, pqueue, gnodearr, cost, costrev, layer, distarr, start[r], result, vcount,
                  nodenterms, basearr, edgearr, (r == 0), connected) );


         obj = 0.0;
         objt = 0.0;
         /* here another measure than in the do_(...) heuristics is being used*/
         for( e = 0; e < nedges; e++)
            obj += (result[e] > -1) ? graph->cost[e] : 0.0;
         for( e = 0; e < nedges; e++)
            objt += (result[e] > -1) ? cost[e] : 0.0;
         //SCIPdebugMessage(" Obj=%.12e\n", obj);

         if( SCIPisLT(scip, obj, min) )
         {
            min = obj;

            SCIPdebugMessage(" Objt=%.12e    ", objt);
            printf(" Obj(run: %d)=%.12e\n", r, obj);

            for( e = 0; e < nedges; e++ )
               best_result[e] = result[e];

            best = start[r];
         }
      }
   }

   /* free allocated memory */
   if( mode == TM )
   {
      for( k = 0; k < nnodes; k++ )
      {
         assert(path[k] == NULL || graph->term[k] == layer);
         SCIPfreeBufferArrayNull(scip, &(path[k]));
      }
      SCIPfreeBufferArray(scip, &path);
   }
   else if( mode == TMPOLZIN )
   {
      SCIPpqueueFree(&pqueue);
      for( k = nnodes - 1; k >= 0; k-- )
      {
         SCIPfreeBuffer(scip, &gnodearr[k]);
	 SCIPfreeBufferArray(scip, &distarr[k]);
	 SCIPfreeBufferArray(scip, &edgearr[k]);
	 SCIPfreeBufferArray(scip, &basearr[k]);
      }
      SCIPfreeBufferArray(scip, &distarr);
      SCIPfreeBufferArray(scip, &edgearr);
      SCIPfreeBufferArray(scip, &basearr);
      SCIPfreeBufferArray(scip, &gnodearr);
      SCIPfreeBufferArray(scip, &vcount);
      SCIPfreeBufferArray(scip, &nodenterms);
   }

   /* NOW IN EXTRA FILE
    *SCIP_CALL( do_local(scip, graph, cost, costrev, best_result) ); */
   SCIPfreeBufferArray(scip, &result);
   SCIPfreeBufferArray(scip, &start);
   SCIPfreeBufferArray(scip, &connected);

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
   SCIP_Real* costrev;
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
   {
      if( SCIPgetDepth(scip) > 0 )
         return SCIP_OKAY;

      runs = heurdata->initruns;
   }
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
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, graph->edges) );
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
         /* TODO chg. for asymmetric graphs */
         for( e = 0; e < graph->edges; e += 2 )
         {
            costrev[e] = cost[e + 1];
            costrev[e + 1] = cost[e];
         }
      }
      else
      {
         /* swap costs; set a high cost if the variable is fixed to 0 */
         for( e = 0; e < graph->edges; e += 2)
         {
            if( SCIPvarGetUbLocal(vars[layer * graph->edges + e + 1]) < 0.5 )
            {
               costrev[e] = 1e+10; /* ???? why does FARAWAY/2 not work? */
               cost[e + 1] = 1e+10;
            }
            else
            {
               costrev[e] = ((1.0 - xval[layer * graph->edges + e + 1]) * graph->cost[e + 1]);
               cost[e + 1] = costrev[e];
            }

            if( SCIPvarGetUbLocal(vars[layer * graph->edges + e]) < 0.5 )
            {
               costrev[e + 1] = 1e+10; /* ???? why does FARAWAY/2 not work? */
               cost[e] = 1e+10;
            }
            else
            {
               costrev[e + 1] = ((1.0 - xval[layer * graph->edges + e]) * graph->cost[e]);
               cost[e] = costrev[e + 1];
            }
         }
      }
      /* can we connect the network */
      SCIP_CALL( do_layer(scip, graph, layer, results, runs, cost, costrev) );

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
   SCIPfreeBufferArray(scip, &costrev);
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
