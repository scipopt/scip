/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc_stp.c
 * @brief  Miscellaneous methods used for solving Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file includes miscellaneous methods used for solving Steiner problems. For more details see \ref MISCSTP page.
 *
 * @page MISCSTP Miscellaneous methods used for Steiner tree problems
 *
 * This file implements an integer data linked list, a linear link-cut tree, a union-find data structure
 * and a pairing heap.
 *
 * A list of all interface methods can be found in misc_stp.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "heur_tm.h"
#include "probdata_stp.h"
#include "portab.h"
#include "scip/misc.h"


/** returns maximum of given SCIP_Real values */
SCIP_Real misc_stp_maxReal(
   SCIP_Real*            realarr,            /**< array of reals */
   unsigned              nreals              /**< size of array of reals */
  )
{
   SCIP_Real max = realarr[0];

   assert(nreals >= 1);

   for( unsigned i = 1; i < nreals; i++ )
      if( realarr[i] > max )
         max = realarr[i];

   return max;
}

/** compares distances of two GNODE structures */
SCIP_DECL_SORTPTRCOMP(GNODECmpByDist)
{
   SCIP_Real first = ((GNODE*)elem1)->dist;
   SCIP_Real second = ((GNODE*)elem2)->dist;
   if( LT(first,second) )
   {
      return -1;
   }
   else if( EQ(first, second) )  /* first == second */
   {
      return 0;
   }
   else
   {
      return 1;
   }
}

/** insert a new node */
SCIP_RETCODE SCIPintListNodeInsert(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX**                 node,               /**< pointer to the last list node */
   int                   nodeval             /**< data of the new node */
   )
{
   IDX* curr;
   curr = *node;

   SCIP_CALL( SCIPallocBlockMemory(scip, node) );
   (*node)->index = nodeval;
   (*node)->parent = (curr);

   return SCIP_OKAY;
}

/** append copy of list pertaining to node2 to node1 */
SCIP_RETCODE SCIPintListNodeAppendCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX**                 node1,              /**< pointer to the last node of list to be enlarged */
   IDX*                  node2,              /**< pointer to the last node of source list */
   SCIP_Bool*            conflict            /**< pointer to store whether a conflict has been detected by the method */
   )
{
   IDX* new;
   IDX* curr1;
   IDX* curr2;
   int curr2idx;
   SCIP_Bool checkconflict;

   assert(scip != NULL);

   if( node2 == NULL )
      return SCIP_OKAY;

   if( conflict != NULL )
   {
      *conflict = FALSE;
      checkconflict = TRUE;
   }
   else
   {
      checkconflict = FALSE;
   }

   new = NULL;
   curr1 = *node1;
   curr2 = node2;

   if( curr1 != NULL )
   {
      int* pointer;
      int* hasharr;
      int i;
      int nelems = 0;
      const SCIP_Bool ignoreconflicts = (SCIPprobdataGetType(scip) == STP_SPG);

      if( !ignoreconflicts )
      {
         /* todo fix this hack; don't check at all, but just add */
         const int maxlength = SCIPprobdataGetNorgEdges(scip);

         assert(maxlength > 0);

         SCIP_CALL(SCIPallocCleanBufferArray(scip, &hasharr, maxlength));
         SCIP_CALL(SCIPallocCleanBufferArray(scip, &pointer, maxlength));
         while( curr1 != NULL )
         {
            i = curr1->index;
            assert(i < maxlength && nelems < maxlength);
            pointer[nelems++] = i;
            hasharr[i] = 1;
            curr1 = curr1->parent;
         }
      }

      curr1 = *node1;
      while( curr2 != NULL )
      {
         curr2idx = curr2->index;

         if( ignoreconflicts || hasharr[curr2idx] == 0 )
         {
            SCIP_CALL( SCIPallocBlockMemory(scip, &new) );
            new->index = curr2idx;
            new->parent = curr1;
            curr1 = new;
         }
         else if( checkconflict )
         {
            assert(conflict != NULL);
            (*conflict) = TRUE;
         }

         curr2 = curr2->parent;
      }

      if( !ignoreconflicts )
      {
         for( i = 0; i < nelems; i++ )
         {
            hasharr[pointer[i]] = 0;
            pointer[i] = 0;
         }
         SCIPfreeCleanBufferArray(scip, &pointer);
         SCIPfreeCleanBufferArray(scip, &hasharr);
      }
   }
   else
   {
      while( curr2 != NULL )
      {
         SCIP_CALL( SCIPallocBlockMemory(scip, &new) );
         new->index = curr2->index;
         new->parent = curr1;
         curr1 = new;

         curr2 = curr2->parent;
      }
   }

   if( new != NULL )
      *node1 = new;

   return SCIP_OKAY;
}

/** free list */
void SCIPintListNodeFree(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX**                 node                /**< pointer to the last list node */
   )
{
   IDX* curr;
   curr = *node;

   while( curr != NULL )
   {
      *node = curr->parent;
      SCIPfreeBlockMemory(scip, &curr);
      curr = *node;
   }
   assert(*node == NULL);
}

/*
 * Linear Link Cut Tree
 */

/** inits a node, setting 'parent' and 'edge' to its default values */
void SCIPlinkcuttreeInit(
   NODE*                 v                   /**< pointer to node representing the tree */
   )
{
   v->parent = NULL;
   v->edge = -1;
}

/** renders v a child of w; v has to be the root of its tree */
void SCIPlinkcuttreeLink(
   NODE*                 v,                  /**< pointer to node representing the tree */
   NODE*                 w,                  /**< pointer to the child */
   int                   edge                /**< link edge */
   )
{
   assert(v->parent == NULL);
   assert(v->edge == -1);
   v->parent = w;
   v->edge = edge;
}

/** cut tree at given node */
void SCIPlinkcuttreeCut(
   NODE*                 v                   /**< node to cut at */
   )
{
   v->edge = -1;
   v->parent = NULL;
}

/** finds minimum weight chain between node 'start' and distinct root node **/
SCIP_Real SCIPlinkcuttreeFindMinChain(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      nodeweight,         /**< node weight array */
   const int*            head,               /**< head of an arc */
   const int*            stdeg,              /**< degree in Steiner tree */
   const NODE*           start,              /**< the node to start at */
   NODE**                first,              /**< first node of chain */
   NODE**                last                /**< last node of chain */
   )
{
   NODE* curr;
   NODE* tmpfirst;
   SCIP_Real min;
   SCIP_Real tmpmin;
   SCIP_Bool stopped;
   int node;

   assert(scip != NULL);
   assert(head != NULL);
   assert(nodeweight != NULL);
   assert(stdeg != NULL);
   assert(start != NULL);
   assert(first != NULL);
   assert(last != NULL);

   *last = NULL;
   *first = NULL;

   min = 0.0;
   tmpmin = 0.0;
   stopped = TRUE;

   /* while p is not root */
   for( curr = (NODE*) start; curr->parent != NULL; curr = curr->parent )
   {
      assert(curr->edge >= 0);

      node = head[curr->edge];

      if( stdeg[node] == 2 && !SCIPisPositive(scip, nodeweight[node]) && curr->parent->parent != NULL )
      {
         if( stopped )
         {
            stopped = FALSE;
            tmpmin = 0.0;
            tmpfirst = curr;
         }
         tmpmin += nodeweight[node];
      }
      else
      {
         /* better chain found? */
         if( !stopped && SCIPisLT(scip, tmpmin, min) )
         {
            assert(tmpfirst != NULL);
            min = tmpmin;
            *first = tmpfirst;
            *last = curr;
         }
         stopped = TRUE;
      }
   }

   return min;
}



/** finds the max edge value between node 'v' and the root of the tree **/
NODE* SCIPlinkcuttreeFindMax(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      cost,               /**< edge cost array */
   NODE*                 v                   /**< the node */
   )
{
   NODE* p = v;
   NODE* q = v;
   SCIP_Real max = -1;

   /* while p is not the root */
   while( p->parent != NULL )
   {
      assert(p->edge >= 0);
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
   NODE*                 v                   /**< the vertex to become the root */
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
      if( val != -1 )
         q->edge = flipedge(val);
      else
         q->edge = -1;
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
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE                *root1,             /**< pointer to root of first heap */
   PHNODE                *root2              /**< pointer to root of second heap */
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

/** add heap to heap */
PHNODE* SCIPpairheapAddtoheap(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE*               root1,              /**< pointer to root of first heap */
   PHNODE*               root2               /**< pointer to root of second heap */
   )
{
   assert(root2 != NULL);
   assert(root1 != NULL);

   if( root1->key <= root2->key )
   {
      /* attach root2 as (the leftmost) child of root1 */
      root2->prev = root1;
      root1->sibling = root2->sibling;
      /* @todo: should never happen */
      if( root1->sibling != NULL )
      {
         root1->sibling->prev = root1;
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
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              p,                  /**< the first sibling */
   int                   size                /**< the size of the heap */
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
      assert(size > nsiblings);
      treearray[nsiblings] = (*p);
      if( (*p)->prev != NULL )
         (*p)->prev->sibling = NULL;
      (*p) = (*p)->sibling;
   }
   assert(size > nsiblings);
   treearray[nsiblings] = NULL;

   /* combine the subtrees (two at a time) */
   for( i = 0; i < nsiblings - 1; i += 2 )
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
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              root,               /**< pointer to root of the heap */
   int                   element,            /**< data of new node */
   SCIP_Real             key,                /**< key of new node */
   int*                  size                /**< pointer to size of the heap */
   )
{
   if( (*root) == NULL )
   {
      (*size) = 1;
      SCIP_CALL( SCIPallocBlockMemory(scip, root) );
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
      SCIP_CALL( SCIPallocBlockMemory(scip, &node) );
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
SCIP_RETCODE SCIPpairheapDeletemin(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  element,            /**< data of the root */
   SCIP_Real*            key,                /**< key of the root */
   PHNODE**              root,               /**< pointer to root of the heap */
   int*                  size                /**< pointer to size of the heap */
   )
{
   assert(scip != NULL);
   if( (*root) == NULL )
   {
      *element = -1;
      return SCIP_OKAY;
   }
   else
   {
      PHNODE *newroot = NULL;

      assert(key != NULL);
      assert(size != NULL);

      *element = (*root)->element;
      *key = (*root)->key;
      if( (*root)->child != NULL )
      {
         newroot = (*root)->child;
         SCIP_CALL( pairheapCombineSiblings(scip, &newroot, (*size)--) );
      }

      SCIPfreeBlockMemory(scip, root);
      (*root) = newroot;
   }
   return SCIP_OKAY;
}


/** links nodes 'root1' and 'root2' together, roots the resulting tree at root1 and sets root2 to NULL */
void SCIPpairheapMeldheaps(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              root1,              /**< pointer to root of first heap */
   PHNODE**              root2,              /**< pointer to root of second heap */
   int*                  sizeroot1,          /**< pointer to size of first heap */
   int*                  sizeroot2           /**< pointer to size of second heap */
   )
{
   assert(scip != NULL);
   assert(root1 != NULL);
   assert(root2 != NULL);
   assert(sizeroot1 != NULL);
   assert(sizeroot2 != NULL);

   if( *root1 == NULL && *root2 == NULL )
   {
      assert(*sizeroot1 == 0);
      assert(*sizeroot2 == 0);
      return;
   }

   (*root1) = SCIPpairheapMergeheaps(scip, *root1, *root2);
   (*sizeroot1) += (*sizeroot2);
   (*root2) = NULL;
}


/** frees the paring heap with root 'p' */
void SCIPpairheapFree(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE**              root                /**< root of heap to be freed */
   )
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

   SCIPfreeBlockMemory(scip, root);
   (*root) = NULL;

}


/** internal method used by 'pairheap_buffarr' */
static
void pairheapRec(
   PHNODE* p,
   int** arr,
   int* n
   )
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
SCIP_RETCODE SCIPpairheapBuffarr(
   SCIP*                 scip,               /**< SCIP data structure */
   PHNODE*               root,               /**< root of the heap */
   int                   size,               /**< size of the array */
   int**                 elements            /**< pointer to array */
   )
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
SCIP_RETCODE SCIPStpunionfindInit(
   SCIP*                 scip,               /**< SCIP data structure */
   UF*                   uf,                 /**< union find data structure */
   int                   length              /**< number of components */
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

/** clears the union-find structure 'uf'*/
void SCIPStpunionfindClear(
   SCIP*                 scip,               /**< SCIP data structure */
   UF*                   uf,                 /**< union find data structure */
   int                   length              /**< number of components */
   )
{
   int i;
   uf->count = length;

   for( i = 0; i < length; i++ )
   {
      uf->parent[i] = i;
      uf->size[i] = 1;
   }

   return;
}


/** finds and returns the component identifier */
int SCIPStpunionfindFind(
   UF*                   uf,                 /**< union find data structure */
   int                   element             /**< element to be found */
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
void SCIPStpunionfindUnion(
   UF*                   uf,                 /**< union find data structure */
   int                   p,                  /**< first component */
   int                   q,                  /**< second component*/
   SCIP_Bool             compress            /**< compress union find structure? */
   )
{
   int idp;
   int idq;
   int* size = uf->size;
   int* parent = uf->parent;
   idp = SCIPStpunionfindFind(uf, p);
   idq = SCIPStpunionfindFind(uf, q);

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
void SCIPStpunionfindFree(
   SCIP*                 scip,               /**< SCIP data structure */
   UF*                   uf                  /**< union find data structure */
   )
{
   SCIPfreeMemoryArray(scip, &uf->parent);
   SCIPfreeMemoryArray(scip, &uf->size);
}
