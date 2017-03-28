/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sgtrie.h
 * @brief  interface for signature trie
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/sgtrie.h"
#include "scip/struct_sgtrie.h"

#define INIT_SIZE                            8
#define GROW_FAC                             2.0

#define NEW_STACK(type, name)               struct { int size; int nelems; type* elems; } name = {0, 0, NULL}
#define STACK_PUT(bufmem, stack, elem)  do \
{ \
   if( (stack).size == (stack).nelems ) \
   { \
      (stack).size = calcGrowSize(INIT_SIZE, GROW_FAC, (stack).nelems + 1); \
      BMSreallocBufferMemoryArray((bufmem), &(stack).elems, (stack).size); \
   } \
   (stack).elems[(stack).nelems++] = (elem); \
} while(0)

#define STACK_POP(stack)  ( (stack).elems[--(stack).nelems]   )
#define STACK_PEEK(stack) ( (stack).elems[(stack).nelems - 1] )
#define STACK_SIZE(stack) ( (stack).nelems )
#define DESTROY_STACK(bufmem, stack) BMSfreeBufferMemoryArrayNull((bufmem), &(stack).elems)

#define ISOLATE_LSB(x)   ((x) & (~(x)+1u))
#define IS_LEAF(node)  ((node)->mask & UINT64_C(0x1))

/** calculate memory size for dynamically allocated arrays (copied from scip/set.c) */
static
int calcGrowSize(
   int                   initsize,           /**< initial size of array */
   SCIP_Real             growfac,            /**< growing factor of array */
   int                   num                 /**< minimum number of entries to store */
   )
{
   int size;

   assert(initsize >= 0);
   assert(growfac >= 1.0);
   assert(num >= 0);

   if( growfac == 1.0 )
      size = MAX(initsize, num);
   else
   {
      int oldsize;

      /* calculate the size with this loop, such that the resulting numbers are always the same (-> block memory) */
      initsize = MAX(initsize, 4);
      size = initsize;
      oldsize = size - 1;

      /* second condition checks against overflow */
      while( size < num && size > oldsize )
      {
         oldsize = size;
         size = (int)(growfac * size + initsize);
      }

      /* if an overflow happened, set the correct value */
      if( size <= oldsize )
         size = num;
   }

   assert(size >= initsize);
   assert(size >= num);

   return size;
}


static
uint64_t getMaskFromDiff(
   uint64_t x
   )
{
   x |= x >> 32;
   x |= x >> 16;
   x |= x >> 8;
   x |= x >> 4;
   x |= x >> 2;
   x |= x >> 1;

   return x;
}

/** returns the number of bits set to one in the given integer */
static
int populationCount(
   uint64_t x
   )
{
#if defined(__GNUC__) && defined(__x86_64__)
   /* for gcc in 64bit mode use the intrinsic since it will compile to a single pop count instruction */
   return __builtin_popcountll(x);
#else
   /* otherwise use a software implementation (see https://en.wikipedia.org/wiki/Hamming_weight) */
   x -= (x >> 1) & UINT64_C(0x5555555555555555);
   x = (x & UINT64_C(0x3333333333333333)) + ((x >> 2) & UINT64_C(0x3333333333333333));
   x = (x + (x >> 4)) & UINT64_C(0x0f0f0f0f0f0f0f0f);
   x += x >> 8;
   x += x >> 16;
   x += x >> 32;
   return (int) (x & 0x7fu);
#endif
}

static
SCIP_RETCODE allocNode(
   SCIP_SGTRIE*          sgtrie,
   SCIP_SGTRIENODE**     node,
   uint64_t              signature,
   uint64_t              mask
   )
{
   SCIP_ALLOC( BMSallocBlockMemory(sgtrie->blkmem, node) );

   (*node)->prefix = signature;
   (*node)->mask = mask;

   return SCIP_OKAY;
}

static
void freeNode(
   SCIP_SGTRIE*          sgtrie,
   SCIP_SGTRIENODE**     node
   )
{
   BMSfreeBlockMemory(sgtrie->blkmem, node);

   assert(*node == NULL);
}

static
SCIP_Bool signatureIsSubset(
   uint64_t a,
   uint64_t b,
   uint64_t mask
   )
{
   return (mask & ( a & ~b )) == 0;
}

static
SCIP_Bool signatureIsEqual(
   uint64_t a,
   uint64_t b,
   uint64_t mask
   )
{
   return ((a ^ b) & mask) == 0;
}

int SCIPsgtrieGetNElems(
   SCIP_SGTRIE*          sgtrie
   )
{
   return sgtrie->nelements;
}

SCIP_RETCODE SCIPsgtrieCreate(
   SCIP_SGTRIE**         sgtrie,
   BMS_BLKMEM*           blkmem,
   BMS_BUFMEM*           bufmem,
   SCIP_DECL_GETSIGNATURE ((*getsignature)),
   SCIP_DECL_SETCMP      ((*setcmp))
   )
{
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sgtrie) );

   /* init other fields */
   (*sgtrie)->blkmem = blkmem;
   (*sgtrie)->bufmem = bufmem;
   (*sgtrie)->setcmp = setcmp;
   (*sgtrie)->getsignature = getsignature;
   (*sgtrie)->root = NULL;
   (*sgtrie)->nelements = 0;

   return SCIP_OKAY;
}

void SCIPsgtrieFree(
   SCIP_SGTRIE**         sgtrie
   )
{
   SCIP_SGTRIENODE* node;
   NEW_STACK(SCIP_SGTRIENODE*, nodes);

   node = (*sgtrie)->root;
   STACK_PUT((*sgtrie)->bufmem, nodes, NULL);

   while( node != NULL )
   {
      if( IS_LEAF(node) )
      {
         LEAFNODEDATA* data;
         LEAFNODEDATA* tmp;

         data = node->data.leaf.next;
         while( data != NULL )
         {
            tmp = data;
            data = data->next;
            BMSfreeBlockMemory((*sgtrie)->blkmem, &tmp);
         }
      }
      else
      {
         STACK_PUT((*sgtrie)->bufmem, nodes, node->data.inner.left);
         STACK_PUT((*sgtrie)->bufmem, nodes, node->data.inner.right);
      }
      freeNode(*sgtrie, &node);

      node = STACK_POP(nodes);
   }

   DESTROY_STACK((*sgtrie)->bufmem, nodes);
   BMSfreeBlockMemory((*sgtrie)->blkmem, sgtrie);
}

SCIP_RETCODE SCIPsgtrieInsert(
   SCIP_SGTRIE*           sgtrie,
   void*                  elem
   )
{
   uint64_t         signature;
   SCIP_SGTRIENODE* currnode;
   LEAFNODEDATA*    leafdata;
   LEAFNODEDATA*    prevleafdata;

   currnode = sgtrie->root;
   signature = sgtrie->getsignature(elem);
   if( currnode == NULL )
   {
      SCIP_CALL( allocNode(sgtrie, &currnode, signature, UINT64_C(0xFFFFFFFFFFFFFFFF)) );
      currnode->data.leaf.element = elem;
      currnode->data.leaf.next = NULL;
      sgtrie->root = currnode;
      return SCIP_OKAY;
   }

   while( TRUE )
   {
      uint64_t diff = currnode->mask & (currnode->prefix ^ signature);

      if( diff != 0 )
      {
         uint64_t prefixmask;
         uint64_t oldsuffixmask;
         uint64_t newsuffixmask;
         SCIP_SGTRIENODE* oldsuffix;
         SCIP_SGTRIENODE* newsuffix;

         /* compute the masks for the new nodes */
         newsuffixmask = getMaskFromDiff(diff);
         oldsuffixmask = currnode->mask & newsuffixmask;
         prefixmask = currnode->mask & ~newsuffixmask;

         /* allocate new nodes */
         SCIP_CALL( allocNode(sgtrie, &oldsuffix, currnode->prefix, oldsuffixmask) );
         SCIP_CALL( allocNode(sgtrie, &newsuffix, signature, newsuffixmask) );

         /* set data in the new node */
         newsuffix->data.leaf.element = elem;
         newsuffix->data.leaf.next = NULL;

         /* copy data from old node */
         oldsuffix->data = currnode->data;

         /* update mask and links for current node */
         currnode->mask = prefixmask;
         if( signature & (ISOLATE_LSB(prefixmask)>>1) )
         {
            currnode->data.inner.left = oldsuffix;
            currnode->data.inner.right = newsuffix;
         }
         else
         {
            currnode->data.inner.left = newsuffix;
            currnode->data.inner.right = oldsuffix;
         }
         ++sgtrie->nelements;

         return SCIP_OKAY;
      }
      /* if we reached a leaf stop looping */
      if( IS_LEAF(currnode) )
         break;

      /* continue search for the correct leaf in proper subtree */
      if( signature & (ISOLATE_LSB(currnode->mask)>>1) )
         currnode = currnode->data.inner.right;
      else
         currnode = currnode->data.inner.left;
   }

   leafdata = &currnode->data.leaf;
   prevleafdata = NULL;
   do
   {
      if( elem == leafdata->element || (sgtrie->setcmp != NULL && sgtrie->setcmp(SCIP_SGTRIE_SETEQ, elem, leafdata->element)) )
         return SCIP_KEYALREADYEXISTING;

      prevleafdata = leafdata;
      leafdata = leafdata->next;
   } while( leafdata != NULL );

   /* create new data element and link it to the end */
   BMSallocBlockMemory(sgtrie->blkmem, &prevleafdata->next);
   prevleafdata->next->element = elem;
   prevleafdata->next->next = NULL;
   ++sgtrie->nelements;

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPsgtrieRemove(
   SCIP_SGTRIE*           sgtrie,
   void*                  elem
   )
{
   uint64_t          signature;
   SCIP_SGTRIENODE** prevnode;
   SCIP_SGTRIENODE** currnode;
   LEAFNODEDATA*     leafdata;
   LEAFNODEDATA*     prevleafdata;

   prevnode = NULL;
   currnode = &sgtrie->root;
   signature = sgtrie->getsignature(elem);

   /* find the correct leaf and store the parent in prevnode */
   while( ! IS_LEAF(*currnode) )
   {
      prevnode = currnode;
      if( signature & (ISOLATE_LSB((*currnode)->mask)>>1) )
         currnode = &(*prevnode)->data.inner.right;
      else
         currnode = &(*prevnode)->data.inner.left;
   }

   /* find the leafnodedata that contains the element */
   leafdata = &(*currnode)->data.leaf;
   prevleafdata = NULL;
   while( leafdata != NULL )
   {
      if( leafdata->element == elem )
         break;
      prevleafdata = leafdata;
      leafdata = leafdata->next;
   }

   /* element is must be contained in the trie */
   assert(leafdata != NULL);

   /* element is contained and will now be removed */
   --sgtrie->nelements;
   /* if the node contains other elements we just remove the element from the linked list
    * otherwise the node must be deleted
    */
   if( prevleafdata != NULL )
   {
      /* the element was not the first one in the linked list so
       * we just set the next link of the previous one to the next one
       * of the removed element and free it
       */
      prevleafdata->next = leafdata->next;
      BMSfreeBlockMemory(sgtrie->blkmem, &leafdata);
   }
   else if( leafdata->next != NULL )
   {
      /* the element is the first one in the list and there are more elements after it
       * thus we move the content of the next leafnodedata into the leafnode data which
       * is embedded in the trie node itself and free the old data
       */
      prevleafdata = leafdata->next;
      *leafdata = *(leafdata->next);
      BMSfreeBlockMemory(sgtrie->blkmem, &prevleafdata);
   }
   else
   {
      /* the element was the only one in this leaf node, thus the whole trie node must be removed */
      SCIP_SGTRIENODE* sibling;

      if( prevnode == NULL )
      {
         /* if the pointer for prevnode is NULL this was the only node in the trie, i.e. we now delete
          * the root node
          */
         freeNode(sgtrie, currnode);
         return SCIP_OKAY;
      }

      /* get the sibling of the current node, since it
       * needs to be merged with the nodes parent
       */
      if( &(*prevnode)->data.inner.left == currnode )
      {
         sibling = (*prevnode)->data.inner.right;
      }
      else
      {
         assert(&(*prevnode)->data.inner.right == currnode);
         sibling = (*prevnode)->data.inner.left;
      }
      /* free the current node */
      freeNode(sgtrie, currnode);

      /* the sibling's mask is merged with the parent node's mask */
      sibling->mask |= (*prevnode)->mask;
      /* now the parent node can be replaced by the current node's sibling without any other changes */
      freeNode(sgtrie, prevnode);
      *prevnode = sibling;
   }

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPsgtrieFindSubsets(
   SCIP_SGTRIE*          sgtrie,
   void*                 elem,
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   )
{
   uint64_t signature;
   SCIP_SGTRIENODE* node;
   NEW_STACK(SCIP_SGTRIENODE*, candnodes);

   node = sgtrie->root;
   STACK_PUT(sgtrie->bufmem, candnodes, NULL);

   *nmatches = 0;
   signature = sgtrie->getsignature(elem);
   while( node != NULL )
   {
      if( signatureIsSubset(node->prefix, signature, node->mask) )
      {
         /* if the node is a leaf node we check all elements if they
          * are a subset and store the matches
          */
         if( IS_LEAF(node) )
         {
            LEAFNODEDATA* leafdata;

            /* iterate all elements in this leaf node */
            leafdata = &node->data.leaf;
            do {
               /* call user callback to check if element is a subset
                * and add it to the matches if it is
                */
               if( sgtrie->setcmp == NULL || sgtrie->setcmp(SCIP_SGTRIE_SUBSET, leafdata->element, elem) )
                  matches[(*nmatches)++] = leafdata->element;

               leafdata = leafdata->next;
            } while(leafdata != NULL);
         }
         else
         {
            /* the left subtree can always contain a subset regardless of the bit value
             * at the split index
             */
            STACK_PUT(sgtrie->bufmem, candnodes, node->data.inner.left);

            assert(((ISOLATE_LSB(node->mask)>>1) & node->mask) == 0);
            assert(populationCount((ISOLATE_LSB(node->mask)>>1)) == 1);
            assert((ISOLATE_LSB(node->mask)) & node->mask);

            /* if the signature has a one at the split index we must also check the right subtree */
            if( signature & (ISOLATE_LSB(node->mask)>>1) )
               STACK_PUT(sgtrie->bufmem, candnodes, node->data.inner.right);
         }
      }
      node = STACK_POP(candnodes);
   }

   DESTROY_STACK(sgtrie->bufmem, candnodes);

   return SCIP_OKAY;
}


typedef struct
{
   SCIP_SGTRIENODE*      node;
   int                   distance;
}  NODEDISTPAIR;

SCIP_RETCODE SCIPsgtrieFindSubsetsPlusOne(
   SCIP_SGTRIE*          sgtrie,
   void*                 elem,
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   )
{
   uint64_t signature;
   NODEDISTPAIR current;
   NEW_STACK(NODEDISTPAIR, stack);

   /* put null node on the stack to mark the end */
   current.node = NULL;
   STACK_PUT(sgtrie->bufmem, stack, current);

   /* start with root */
   current.node = sgtrie->root;
   current.distance = 0;

   *nmatches = 0;
   signature = sgtrie->getsignature(elem);
   while( current.node != NULL )
   {
      uint64_t subsetmask = current.node->mask & ( current.node->prefix & ~signature );
      current.distance += populationCount(subsetmask);

      if( current.distance <= 1 )
      {
         /* if the node is a leaf node all elements are candidates for being a subset plus one extra element
          */
         if( IS_LEAF(current.node) )
         {
            LEAFNODEDATA* leafdata;
            /* iterate all elements in this leaf node */
            leafdata = &current.node->data.leaf;
            do {
               /* call user callback to check if element is a subset
                * and add it to the matches if it is
                */

               if( sgtrie->setcmp == NULL || sgtrie->setcmp(SCIP_SGTRIE_SUBSETPLUSONE, leafdata->element, elem) )
                  matches[(*nmatches)++] = leafdata->element;

               leafdata = leafdata->next;
            } while(leafdata != NULL);
         }
         else
         {
            /* add the left and the right node to the stack with the current distance */
            SCIP_SGTRIENODE* node;

            node = current.node;

            /* the left subtree can always contain a subset regardless of the bit value
             * at the split index
             */
            current.node = node->data.inner.left;
            STACK_PUT(sgtrie->bufmem, stack, current);

            current.node = node->data.inner.right;
            STACK_PUT(sgtrie->bufmem, stack, current);
         }
      }

      current = STACK_POP(stack);
   }

   DESTROY_STACK(sgtrie->bufmem, stack);

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPsgtrieFindSupersets(
   SCIP_SGTRIE*          sgtrie,
   void*                 elem,
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   )
{
   uint64_t signature;
   SCIP_SGTRIENODE* node;
   NEW_STACK(SCIP_SGTRIENODE*, candnodes);

   node = sgtrie->root;
   STACK_PUT(sgtrie->bufmem, candnodes, NULL);

   *nmatches = 0;
   signature = sgtrie->getsignature(elem);
   while( node != NULL )
   {
      if( signatureIsSubset(signature, node->prefix, node->mask) )
      {
         /* if the node is a leaf node we check all elements if they
          * are a subset and store the matches
          */
         if( IS_LEAF(node) )
         {
            LEAFNODEDATA* leafdata;

            /* iterate all elements in this leaf node */
            leafdata = &node->data.leaf;
            do {
               /* call user callback to check if element is a subset
                * and add it to the matches if it is
                */
               if( sgtrie->setcmp == NULL || sgtrie->setcmp(SCIP_SGTRIE_SUBSET, elem, leafdata->element) )
                  matches[(*nmatches)++] = leafdata->element;

               leafdata = leafdata->next;
            } while(leafdata != NULL);
         }
         else
         {
            /* the right subtree can always contain a superset regardless of the bit value
             * at the split index
             */
            STACK_PUT(sgtrie->bufmem, candnodes, node->data.inner.right);

            assert(((ISOLATE_LSB(node->mask)>>1) & node->mask) == 0);
            assert(populationCount((ISOLATE_LSB(node->mask)>>1)) == 1);
            assert((ISOLATE_LSB(node->mask)) & node->mask);

            /* if the signature has a zero at the split index we must also check the left subtree */
            if( ! (signature & (ISOLATE_LSB(node->mask)>>1)) )
               STACK_PUT(sgtrie->bufmem, candnodes, node->data.inner.left);
         }
      }
      node = STACK_POP(candnodes);
   }

   DESTROY_STACK(sgtrie->bufmem, candnodes);

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPsgtrieFindSupersetsPlusOne(
   SCIP_SGTRIE*          sgtrie,
   void*                 elem,
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   )
{
   uint64_t signature;
   NODEDISTPAIR current;
   NEW_STACK(NODEDISTPAIR, stack);

   /* put null node on the stack to mark the end */
   current.node = NULL;
   STACK_PUT(sgtrie->bufmem, stack, current);

   /* start with root */
   current.node = sgtrie->root;
   current.distance = 0;

   *nmatches = 0;
   signature = sgtrie->getsignature(elem);
   while( current.node != NULL )
   {
      uint64_t subsetmask = current.node->mask & ( signature & ~current.node->prefix );
      current.distance += populationCount(subsetmask);

      if( current.distance <= 1 )
      {
         /* if the node is a leaf node all elements are candidates for being a subset plus one extra element
          */
         if( IS_LEAF(current.node) )
         {
            LEAFNODEDATA* leafdata;
            /* iterate all elements in this leaf node */
            leafdata = &current.node->data.leaf;
            do {
               /* call user callback to check if element is a subset
                * and add it to the matches if it is
                */

               if( sgtrie->setcmp == NULL || sgtrie->setcmp(SCIP_SGTRIE_SUBSETPLUSONE, elem, leafdata->element) )
                  matches[(*nmatches)++] = leafdata->element;

               leafdata = leafdata->next;
            } while(leafdata != NULL);
         }
         else
         {
            /* add the left and the right node to the stack with the current distance */
            SCIP_SGTRIENODE* node;

            node = current.node;

            /* the left subtree can always contain a subset regardless of the bit value
             * at the split index
             */
            current.node = node->data.inner.left;
            STACK_PUT(sgtrie->bufmem, stack, current);

            current.node = node->data.inner.right;
            STACK_PUT(sgtrie->bufmem, stack, current);
         }
      }

      current = STACK_POP(stack);
   }

   DESTROY_STACK(sgtrie->bufmem, stack);

   return SCIP_OKAY;
}

void* SCIPsgtrieFindEq(
   SCIP_SGTRIE*          sgtrie,
   void*                 elem
   )
{
   uint64_t signature;
   SCIP_SGTRIENODE* currnode;

   /* start with the root node and distance zero */
   currnode = sgtrie->root;
   if( currnode == NULL )
      return NULL;

   signature = sgtrie->getsignature(elem);

   while( TRUE )
   {
      if( !signatureIsEqual(signature, currnode->prefix, currnode->mask) )
         break;

      if( IS_LEAF(currnode) )
      {
         LEAFNODEDATA* leafdata;
         /* iterate all elements in this leaf node */
         leafdata = &currnode->data.leaf;
         do {
            /* call user callback to check if element is a subset
             * and add it to the matches if it is
             */
            if( elem == leafdata->element || (sgtrie->setcmp != NULL && sgtrie->setcmp(SCIP_SGTRIE_SETEQ, elem, leafdata->element)) )
               return leafdata->element;

            leafdata = leafdata->next;
         } while(leafdata != NULL);

         break;
      }

      if( signature & (ISOLATE_LSB(currnode->mask)>>1) )
         currnode = currnode->data.inner.right;
      else
         currnode = currnode->data.inner.left;
   }

   return NULL;
}
