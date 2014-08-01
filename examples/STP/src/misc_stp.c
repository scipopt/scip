/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc_stp.c
 * @brief  miscellaneous methods
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "heur_tm.h"
#include "probdata_stp.h"
#include "portab.h"
#include "scip/misc.h"


/* compares distances of two GNODE structures */
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


/*
 * Linear Link Cut Tree
 */

/** inits a node, setting 'parent' and 'edge' to its default values */
void SCIPlinkcuttreeInit(
   NODE* v
   )
{
   v->parent = NULL;
   v->edge = -1;
}

/** renders w a child of v; v has to be the root of its tree */
void SCIPlinkcuttreeLink(
   NODE* v,
   NODE* w,
   int edge
   )
{
   assert(v->parent == NULL);
   v->parent = w;
   v->edge = edge;
}

void SCIPlinkcuttreeCut(
   NODE* v
   )
{
   v->edge = -1;
   v->parent = NULL;
}

/** finds the max value between node 'v' and the root of the tree **/
NODE* SCIPlinkcuttreeFindMax(
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
void SCIPlinkcuttreeEvert(
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



/*
 * Pairing Heap
 */

/** links nodes 'root1' and 'root2' together */
PHNODE* SCIPpairheapMergeheaps(
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

PHNODE* SCIPpairheapAddtoheap(
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
SCIP_RETCODE pairheapCombineSiblings(
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
      treearray[i] = SCIPpairheapMergeheaps(scip, treearray[i-1], treearray[i]);


   return treearray[nsiblings-1];
#endif

   /* combine the subtrees (two at a time) */
   for( i = 0; i + 1 < nsiblings; i += 2 )
   {
      treearray[i] = SCIPpairheapMergeheaps(scip, treearray[i], treearray[i + 1]);
   }
   j = i - 2;

   /* if the number of trees is odd, get the last one */
   if( j == nsiblings - 3 )
   {
      treearray[j] = SCIPpairheapMergeheaps(scip, treearray[j], treearray[j + 2]);
   }

   for( ; j >= 2; j -= 2 )
   {
      treearray[j - 2] = SCIPpairheapMergeheaps(scip, treearray[j - 2], treearray[j]);
   }

   (*p) = treearray[0];

   SCIPfreeBufferArray(scip, &treearray);

   return SCIP_OKAY;
}


/** inserts a new node into the pairing heap */
SCIP_RETCODE SCIPpairheapInsert(
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
      (*root) = SCIPpairheapAddtoheap(scip, (*root), node);
   }
   return SCIP_OKAY;
}

/** deletes the root of the paring heap, concomitantly storing its data and key in '*element' and '*key' respectively */
void SCIPpairheapDeletemin(
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
         pairheapCombineSiblings(scip, &newroot, --(*size));
      }

      SCIPfreeBuffer(scip, root);
      (*root) = newroot;
   }
}


/** links nodes 'root1' and 'root2' together, roots the resulting tree at root1 and sets root2 to NULL */
void SCIPpairheapMeldheaps(
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

   (*root1) = SCIPpairheapMergeheaps(scip, *root1, *root2);
   (*sizeroot1) += (*sizeroot2);
   (*root2) = NULL;
}


/** frees the paring heap with root 'p' */
void SCIPpairheapFree(SCIP* scip, PHNODE** root)
{
   if( (*root) == NULL )
   {
      return;
   }
   if( (*root)->sibling != NULL )
   {
      SCIPpairheapFree(scip, &((*root)->sibling));
   }
   if( (*root)->child != NULL )
   {
      SCIPpairheapFree(scip, &((*root)->child));
   }

   SCIPfreeBuffer(scip, root);
   (*root) = NULL;

}


/** internal method used by 'pairheap_buffarr' */
static
void pairheapRec(PHNODE* p, int** arr, int* n)
{
   if( p == NULL )
   {
      return;
   }
   (*arr)[(*n)++] = p->element;
   pairheapRec(p->sibling, arr, n);
   pairheapRec(p->child, arr, n);
}


/** stores all elements of the pairing heap in an array */
SCIP_RETCODE SCIPpairheapBuffarr(SCIP* scip, PHNODE* root, int size, int** elements)
{
   int n = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, elements, size) );
   pairheapRec(root, elements, &n);
   return SCIP_OKAY;
}


/*
 *Union-Find data structure
 */

/** initializes the union-find structure 'uf' with 'length' many components (of size one) */
SCIP_RETCODE SCIPunionfindInit(
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
int SCIPunionfindFind(
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
void SCIPunionfindUnion(
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
   idp = SCIPunionfindFind(uf, p);
   idq = SCIPunionfindFind(uf, q);

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
void SCIPunionfindFree(
   SCIP* scip,
   UF* uf
   )
{
   SCIPfreeMemoryArray(scip, &uf->parent);
   SCIPfreeMemoryArray(scip, &uf->size);
   uf = NULL;
}
