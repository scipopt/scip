/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc.c
 * @brief  miscellaneous datastructures and methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/** priority queue data structure
 *  Elements are stored in an array, which grows dynamically in size as new elements are added to the queue.
 *  The ordering is done through a pointer comparison function.
 *  The array is organized as follows. The root element (that is the "best" element $r$ with $r <= x$ for all $x$)
 *  is stored in position 0. The children of an element at position $p$ are stored at positions $q_1 = 2*p+1$ and
 *  $q_2 = 2*p+2$. That means, the parent of the element at position $q$ is at position $p = \lfloor (q-1)/2 \rfloor$.
 *  At any time, the condition holds that $p <= q$ for each parent $p$ and its children $q$.
 *  Insertion and removal of single elements needs time $\mathcal{O}(\log n)$.
 */

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "memory.h"
#include "message.h"
#include "retcode.h"
#include "misc.h"


#if 0 /* PRIORITY QUEUE NOT NEEDED */

#define PQ_PARENT(q) (((q)-1)/2)
#define PQ_LEFTCHILD(p) (2*(p)+1)
#define PQ_RIGHTCHILD(p) (2*(p)+2)


/** priority queue data structure */
struct PQueue
{
   int              len;                /**< number of used element slots */
   int              size;               /**< total number of available element slots */
   Real             sizefac;            /**< memory growing factor */
   void**           slots;              /**< array of element slots */
   DECL_SORTPTRCOMP((*ptrcmp));         /**< compares two data elements */
};


/** resizes element memory to hold at least the given number of elements */
static
RETCODE pqueueResize(
   PQUEUE*          pqueue,             /**< pointer to a priority queue */
   int              minsize             /**< minimal number of storeable elements */
   )
{
   assert(pqueue != NULL);
   
   if( minsize <= pqueue->size )
      return SCIP_OKAY;

   pqueue->size = MAX(minsize, (int)(pqueue->size * pqueue->sizefac));
   ALLOC_OKAY( reallocMemoryArray(&pqueue->slots, pqueue->size) );

   return SCIP_OKAY;
}

/** initializes priority queue */
RETCODE SCIPpqueueInit(
   PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int              initsize,           /**< initial number of available element slots */
   Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   )
{
   assert(pqueue != NULL);
   assert(ptrcmp != NULL);

   initsize = MAX(1, initsize);
   sizefac = MAX(1.0, sizefac);

   ALLOC_OKAY( allocMemory(pqueue) );
   (*pqueue)->len = 0;
   (*pqueue)->size = 0;
   (*pqueue)->sizefac = sizefac;
   (*pqueue)->slots = NULL;
   (*pqueue)->ptrcmp = ptrcmp;
   CHECK_OKAY( pqueueResize(*pqueue, initsize) );

   return SCIP_OKAY;
}

/** frees priority queue, but not the data elements themselves */
void SCIPpqueueFree(
   PQUEUE**         pqueue              /**< pointer to a priority queue */
   )
{
   assert(pqueue != NULL);

   freeMemoryArray(&(*pqueue)->slots);
   freeMemory(pqueue);
}

/** inserts element into priority queue */
RETCODE SCIPpqueueInsert(
   PQUEUE*          pqueue,             /**< pointer to a priority queue */
   void*            elem                /**< element to be inserted */
   )
{
   int pos;

   assert(pqueue != NULL);
   assert(pqueue->len >= 0);
   assert(elem != NULL);

   CHECK_OKAY( pqueueResize(pqueue, pqueue->len+1) );

   /* insert element as leaf in the tree, move it towards the root as long it is better than its parent */
   pos = pqueue->len;
   pqueue->len++;
   while( pos > 0 && (*pqueue->ptrcmp)(elem, pqueue->slots[PQ_PARENT(pos)]) < 0 )
   {
      pqueue->slots[pos] = pqueue->slots[PQ_PARENT(pos)];
      pos = PQ_PARENT(pos);
   }
   pqueue->slots[pos] = elem;

   return SCIP_OKAY;
}

/** removes and returns best element from the priority queue */
void* SCIPpqueueRemove(
   PQUEUE*          pqueue              /**< pointer to a priority queue */
   )
{
   void* root;
   void* last;
   int pos;
   int childpos;
   int brotherpos;

   assert(pqueue != NULL);
   assert(pqueue->len >= 0);
   
   if( pqueue->len == 0 )
      return NULL;

   /* remove root element of the tree, move the better child to its parents position until the last element
    * of the queue could be placed in the empty slot */
   root = pqueue->slots[0];
   last = pqueue->slots[pqueue->len-1];
   pqueue->len--;
   pos = 0;
   while( pos < PQ_PARENT(pqueue->len-1) )
   {
      childpos = PQ_LEFTCHILD(pos);
      brotherpos = PQ_RIGHTCHILD(pos);
      if( brotherpos <= pqueue->len && (*pqueue->ptrcmp)(pqueue->slots[brotherpos], pqueue->slots[childpos]) < 0 )
         childpos = brotherpos;
      if( (*pqueue->ptrcmp)(last, pqueue->slots[childpos]) <= 0 )
         break;
      pqueue->slots[pos] = pqueue->slots[childpos];
      pos = childpos;
   }
   assert(pos <= pqueue->len);
   pqueue->slots[pos] = last;

   return root;
}

/** returns the best element of the queue without removing it */
void* SCIPpqueueFirst(
   const PQUEUE*    pqueue              /**< pointer to a priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   if( pqueue->len == 0 )
      return NULL;

   return pqueue->slots[0];
}

#endif




/*
 * Hash Table
 */

typedef struct HashList HASHLIST;       /**< element list to store in a hash table */

/** element list to store in a hash table */
struct HashList
{
   void*            element;            /**< this element */
   HASHLIST*        next;               /**< rest of the hash list */
};

/** hash table data structure */
struct HashTable
{
   HASHLIST**       lists;              /**< hash lists of the hash table */
   int              nlists;             /**< number of lists stored in the hash table */
   DECL_HASHGETKEY((*hashgetkey));      /**< gets the key of the given element */
   DECL_HASHKEYEQ ((*hashkeyeq));       /**< returns TRUE iff both keys are equal */
   DECL_HASHKEYVAL((*hashkeyval));      /**< returns the hash value of the key */
};


/** appends element to the hash list */
static
RETCODE hashlistAppend(
   HASHLIST**       hashlist,           /**< pointer to hash list to free */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to append to the list */
   )
{
   HASHLIST* newlist;

   assert(hashlist != NULL);
   assert(memhdr != NULL);
   assert(element != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, &newlist) );
   newlist->element = element;
   newlist->next = *hashlist;
   *hashlist = newlist;

   return SCIP_OKAY;
}

/** frees a hash list entry and all its successors */
static
void hashlistFree(
   HASHLIST**       hashlist,           /**< pointer to hash list to free */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   HASHLIST* actlist;
   HASHLIST* nextlist;

   assert(hashlist != NULL);
   assert(memhdr != NULL);
   
   actlist = *hashlist;
   while( actlist != NULL )
   {
      nextlist = actlist->next;
      freeBlockMemory(memhdr, &actlist);
      actlist = nextlist;
   }

   *hashlist = NULL;
}

/** retrieves element with given key from the hash list, or NULL */
static
void* hashlistRetrieve(
   HASHLIST*        hashlist,           /**< hash list */
   DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   unsigned int     keyval,             /**< hash value of key */
   void*            key                 /**< key to retrieve */
   )
{
   unsigned int actkeyval;
   void* actkey;

   assert(hashkeyeq != NULL);
   assert(key != NULL);

   while( hashlist != NULL )
   {
      actkey = hashgetkey(hashlist->element);
      actkeyval = hashkeyval(actkey);
      if( actkeyval == keyval && hashkeyeq(actkey, key) )
         return hashlist->element;
      hashlist = hashlist->next;
   }

   return NULL;
}

/** removes element from the hash list */
static
RETCODE hashlistRemove(
   HASHLIST**       hashlist,           /**< pointer to hash list to free */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to remove from the list */
   )
{
   HASHLIST* nextlist;

   assert(hashlist != NULL);
   assert(memhdr != NULL);
   assert(element != NULL);

   while( *hashlist != NULL && (*hashlist)->element != element )
   {
      hashlist = &(*hashlist)->next;
   }
   if( *hashlist != NULL )
   {
      nextlist = (*hashlist)->next;
      freeBlockMemory(memhdr, hashlist);
      *hashlist = nextlist;

      return SCIP_OKAY;
   }
   else
   {
      errorMessage("element not found in the hash table");
      return SCIP_INVALIDDATA;
   }
}

   
/** creates a hash table */
RETCODE SCIPhashtableCreate(
   HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   int              tablesize,          /**< size of the hash table */
   DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   DECL_HASHKEYVAL((*hashkeyval))       /**< returns the hash value of the key */
   )
{
   int i;

   assert(hashtable != NULL);
   assert(tablesize > 0);
   assert(hashgetkey != NULL);
   assert(hashkeyeq != NULL);
   assert(hashkeyval != NULL);

   ALLOC_OKAY( allocMemory(hashtable) );
   ALLOC_OKAY( allocMemoryArray(&(*hashtable)->lists, tablesize) );
   (*hashtable)->nlists = tablesize;
   (*hashtable)->hashgetkey = hashgetkey;
   (*hashtable)->hashkeyeq = hashkeyeq;
   (*hashtable)->hashkeyval = hashkeyval;

   /* initialize hash lists */
   for( i = 0; i < tablesize; ++i )
      (*hashtable)->lists[i] = NULL;

   return SCIP_OKAY;
}

/** frees the hash table */
void SCIPhashtableFree(
   HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   int i;

   assert(hashtable != NULL);
   assert(memhdr != NULL);

   /* free hash lists */
   for( i = 0; i < (*hashtable)->nlists; ++i )
      hashlistFree(&(*hashtable)->lists[i], memhdr);

   /* free main hast table data structure */
   freeMemoryArray(&(*hashtable)->lists);
   freeMemory(hashtable);
}

/** inserts element in hash table (multiple inserts of same element possible) */
RETCODE SCIPhashtableInsert(
   HASHTABLE*       hashtable,          /**< hash table */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to insert into the table */
   )
{
   void* key;
   unsigned int keyval;
   unsigned int hashval;

   assert(hashtable != NULL);
   assert(hashtable->lists != NULL);
   assert(hashtable->nlists > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(memhdr != NULL);
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(element);
   keyval = hashtable->hashkeyval(key);
   hashval = keyval % hashtable->nlists;

   /* append element to the list at the hash position */
   CHECK_OKAY( hashlistAppend(&hashtable->lists[hashval], memhdr, element) );
   
   return SCIP_OKAY;
}

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
RETCODE SCIPhashtableSafeInsert(
   HASHTABLE*       hashtable,          /**< hash table */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to insert into the table */
   )
{
   assert(hashtable != NULL);
   assert(hashtable->hashgetkey != NULL);

   /* check, if key is already existing */
   if( SCIPhashtableRetrieve(hashtable, hashtable->hashgetkey(element)) != NULL )
      return SCIP_KEYALREADYEXISTING;

   /* insert element in hash table */
   CHECK_OKAY( SCIPhashtableInsert(hashtable, memhdr, element) );
   
   return SCIP_OKAY;
}

/** retrieve element with key from hash table, returns NULL if not existing */
void* SCIPhashtableRetrieve(
   HASHTABLE*       hashtable,          /**< hash table */
   void*            key                 /**< key to retrieve */
   )
{
   unsigned int keyval;
   unsigned int hashval;

   assert(hashtable != NULL);
   assert(hashtable->lists != NULL);
   assert(hashtable->nlists > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(key != NULL);

   /* get the hash value of the key */
   keyval = hashtable->hashkeyval(key);
   hashval = keyval % hashtable->nlists;

   return hashlistRetrieve(hashtable->lists[hashval], hashtable->hashgetkey, hashtable->hashkeyeq, hashtable->hashkeyval,
      keyval, key);
}

/** removes existing element from the hash table */
RETCODE SCIPhashtableRemove(
   HASHTABLE*       hashtable,          /**< hash table */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to remove from the table */
   )
{
   void* key;
   unsigned int keyval;
   unsigned int hashval;

   assert(hashtable != NULL);
   assert(hashtable->lists != NULL);
   assert(hashtable->nlists > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(memhdr != NULL);
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(element);
   keyval = hashtable->hashkeyval(key);
   hashval = keyval % hashtable->nlists;

   /* append element to the list at the hash position */
   CHECK_OKAY( hashlistRemove(&hashtable->lists[hashval], memhdr, element) );
   
   return SCIP_OKAY;
}

/** prints statistics about hash table usage */
void SCIPhashtablePrintStatistics(
   HASHTABLE*       hashtable           /**< hash table */
   )
{
   HASHLIST* hashlist;
   int usedslots;
   int maxslotsize;
   int sumslotsize;
   int slotsize;
   int i;

   assert(hashtable != NULL);

   usedslots = 0;
   maxslotsize = 0;
   sumslotsize = 0;
   for( i = 0; i < hashtable->nlists; ++i )
   {
      hashlist = hashtable->lists[i];
      if( hashlist != NULL )
      {
         usedslots++;
         slotsize = 0;
         while( hashlist != NULL )
         {
            slotsize++;
            hashlist = hashlist->next;
         }
         maxslotsize = MAX(maxslotsize, slotsize);
         sumslotsize += slotsize;
      }
   }

   printf("%d hash entries, used %d/%d slots (%.1f%%)",
      sumslotsize, usedslots, hashtable->nlists, 100.0*(Real)usedslots/(Real)(hashtable->nlists));
   if( usedslots > 0 )
      printf(", avg. %.1f entries/used slot, max. %d entries in slot", (Real)sumslotsize/(Real)usedslots, maxslotsize);
   printf("\n");
}


/** returns TRUE iff both keys (i.e. strings) are equal */
DECL_HASHKEYEQ(SCIPhashKeyEqString)
{
   const char* string1 = (const char*)key1;
   const char* string2 = (const char*)key2;

   return (strcmp(string1, string2) == 0);
}

static
const int strhashshift[28] = { 
   0, 16,  8, 24,  4, 20, 12, /*28, */
   2, 18, 10, 26,  6, 22, 14, /*30, */
   1, 17,  9, 25,  5, 21, 13, /*29, */
   3, 19, 11, 27,  7, 23, 15  /*31  */
};

/** returns the hash value of the key (i.e. string) */
DECL_HASHKEYVAL(SCIPhashKeyValString)
{
   unsigned char* string = (unsigned char*)key;
   unsigned int sum;
   unsigned int val;
   int i;

   sum = 0;
   i = 0;
   while( *string != '\0' )
   {
      assert(0 <= i && i < 28);
      val = *string;
      val <<= strhashshift[i];
      sum += val;
      i++;
      i %= 28;
      string++;
   }

   return sum;
}



/*
 * Dynamic Arrays
 */

/** dynamic array for storing real values */
struct RealArray
{
   MEMHDR*          memhdr;             /**< block memory that stores the vals array */
   Real*            vals;               /**< array values */
   int              valssize;           /**< size of vals array */
   int              firstidx;           /**< index of first element in vals array */
   int              minusedidx;         /**< index of first non zero element in vals array */
   int              maxusedidx;         /**< index of last non zero element in vals array */
};

/** creates a dynamic array of real values */
RETCODE SCIPrealarrayCreate(
   REALARRAY**      realarray,          /**< pointer to store the real array */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(realarray != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, realarray) );
   (*realarray)->memhdr = memhdr;
   (*realarray)->vals = NULL;
   (*realarray)->valssize = 0;
   (*realarray)->firstidx = -1;
   (*realarray)->minusedidx = INT_MAX;
   (*realarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of real values */
RETCODE SCIPrealarrayCopy(
   REALARRAY**      realarray,          /**< pointer to store the copied real array */
   MEMHDR*          memhdr,             /**< block memory */
   REALARRAY*       sourcerealarray     /**< dynamic real array to copy */
   )
{
   assert(realarray != NULL);
   assert(sourcerealarray != NULL);

   CHECK_OKAY( SCIPrealarrayCreate(realarray, memhdr) );
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*realarray)->vals, sourcerealarray->vals, sourcerealarray->valssize) );
   (*realarray)->valssize = sourcerealarray->valssize;
   (*realarray)->firstidx = sourcerealarray->firstidx;
   (*realarray)->minusedidx = sourcerealarray->minusedidx;
   (*realarray)->maxusedidx = sourcerealarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
RETCODE SCIPrealarrayFree(
   REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   assert(realarray != NULL);
   assert(*realarray != NULL);

   freeBlockMemoryArrayNull((*realarray)->memhdr, &(*realarray)->vals, (*realarray)->valssize);
   freeBlockMemory((*realarray)->memhdr, realarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPrealarrayExtend(
   REALARRAY*       realarray,          /**< dynamic real array */
   const SET*       set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(realarray != NULL);
   assert(realarray->minusedidx == INT_MAX || realarray->firstidx >= 0);
   assert(realarray->maxusedidx == INT_MIN || realarray->firstidx >= 0);
   assert(realarray->minusedidx == INT_MAX || realarray->minusedidx >= realarray->firstidx);
   assert(realarray->maxusedidx == INT_MIN || realarray->maxusedidx < realarray->firstidx + realarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);
   
   minidx = MIN(minidx, realarray->minusedidx);
   maxidx = MAX(maxidx, realarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   debugMessage("extending realarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > realarray->valssize )
   {
      Real* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = SCIPsetCalcMemGrowSize(set, nused);
      ALLOC_OKAY( allocBlockMemoryArray(realarray->memhdr, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( realarray->firstidx != -1 )
      {
         for( i = 0; i < realarray->minusedidx - newfirstidx; ++i )
            newvals[i] = 0.0;
         copyMemoryArray(&newvals[realarray->minusedidx - newfirstidx],
            &realarray->vals[realarray->minusedidx - realarray->firstidx],
            realarray->maxusedidx - realarray->minusedidx + 1);
         for( i = realarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = 0.0;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = 0.0;
      }

      /* free old memory storage, and set the new array parameters */
      freeBlockMemoryArrayNull(realarray->memhdr, &realarray->vals, realarray->valssize);
      realarray->vals = newvals;
      realarray->valssize = newvalssize;
      realarray->firstidx = newfirstidx;
   }
   else if( realarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = realarray->valssize - nused;
      assert(nfree >= 0);
      realarray->firstidx = minidx - nfree/2;
      assert(realarray->firstidx <= minidx);
      assert(maxidx < realarray->firstidx + realarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < realarray->valssize; ++i )
         assert(realarray->vals[i] == 0.0);
#endif
   }
   else if( minidx < realarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = realarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + realarray->valssize);
      
      if( realarray->minusedidx <= realarray->maxusedidx )
      {
         int shift;

         assert(realarray->firstidx <= realarray->minusedidx);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);

         /* shift used part of array to the right */
         shift = realarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = realarray->maxusedidx - realarray->firstidx; i >= realarray->minusedidx - realarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < realarray->valssize);
            realarray->vals[i + shift] = realarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            realarray->vals[realarray->minusedidx - realarray->firstidx + i] = 0.0;
      }
      realarray->firstidx = newfirstidx;
   }
   else if( maxidx >= realarray->firstidx + realarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = realarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + realarray->valssize);
      
      if( realarray->minusedidx <= realarray->maxusedidx )
      {
         int shift;

         assert(realarray->firstidx <= realarray->minusedidx);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - realarray->firstidx;
         assert(shift > 0);
         for( i = realarray->minusedidx - realarray->firstidx; i <= realarray->maxusedidx - realarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < realarray->valssize);
            realarray->vals[i - shift] = realarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            realarray->vals[realarray->maxusedidx - realarray->firstidx - i] = 0.0;
      }
      realarray->firstidx = newfirstidx;
   }

   assert(minidx >= realarray->firstidx);
   assert(maxidx < realarray->firstidx + realarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic real array */
void SCIPrealarrayClear(
   REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   debugMessage("clearing realarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx);

   if( realarray->minusedidx <= realarray->maxusedidx )
   {
      int i;
   
      assert(realarray->firstidx <= realarray->minusedidx);
      assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
      assert(realarray->firstidx != -1);
      assert(realarray->valssize > 0);

      /* clear the used part of array */
      for( i = realarray->minusedidx - realarray->firstidx; i <= realarray->maxusedidx - realarray->firstidx; ++i )
         realarray->vals[i] = 0.0;

      /* mark the array cleared */
      realarray->minusedidx = INT_MAX;
      realarray->maxusedidx = INT_MIN;
   }
   assert(realarray->minusedidx == INT_MAX);
   assert(realarray->maxusedidx == INT_MIN);
}

/** gets value of entry in dynamic array */
Real SCIPrealarrayGetVal(
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx                 /**< array index to get value for */
   )
{
   assert(realarray != NULL);
   assert(idx >= 0);
   
   if( idx < realarray->minusedidx || idx > realarray->maxusedidx )
      return 0.0;
   else
   {
      assert(realarray->vals != NULL);
      assert(idx - realarray->firstidx >= 0);
      assert(idx - realarray->firstidx < realarray->valssize);

      return realarray->vals[idx - realarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
RETCODE SCIPrealarraySetVal(
   REALARRAY*       realarray,          /**< dynamic real array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   Real             val                 /**< value to set array index to */
   )
{
   assert(realarray != NULL);
   assert(idx >= 0);

   debugMessage("setting realarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %g\n", 
      realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx, idx, val);

   if( !SCIPsetIsZero(set, val) )
   {
      /* extend array to be able to store the index */
      CHECK_OKAY( SCIPrealarrayExtend(realarray, set, idx, idx) );
      assert(idx >= realarray->firstidx);
      assert(idx < realarray->firstidx + realarray->valssize);
      
      /* set the array value of the index */
      realarray->vals[idx - realarray->firstidx] = val;

      /* update min/maxusedidx */
      realarray->minusedidx = MIN(realarray->minusedidx, idx);
      realarray->maxusedidx = MAX(realarray->maxusedidx, idx);
   }
   else if( idx >= realarray->firstidx && idx < realarray->firstidx + realarray->valssize )
   {
      /* set the array value of the index to zero */
      realarray->vals[idx - realarray->firstidx] = 0.0;
      
      /* check, if we can tighten the min/maxusedidx */
      if( idx == realarray->minusedidx )
      {
         assert(realarray->maxusedidx >= 0);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
         do
         {
            realarray->minusedidx++;
         }
         while( realarray->minusedidx <= realarray->maxusedidx
            && SCIPsetIsZero(set, realarray->vals[realarray->minusedidx - realarray->firstidx]) );
         if( realarray->minusedidx > realarray->maxusedidx )
         {
            realarray->minusedidx = INT_MAX;
            realarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == realarray->maxusedidx )
      {
         assert(realarray->minusedidx >= 0);
         assert(realarray->minusedidx < realarray->maxusedidx);
         assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
         do
         {
            realarray->maxusedidx--;
            assert(realarray->minusedidx <= realarray->maxusedidx);
         }
         while( SCIPsetIsZero(set, realarray->vals[realarray->maxusedidx - realarray->firstidx]) );
      }      
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
RETCODE SCIPrealarrayIncVal(
   REALARRAY*       realarray,          /**< dynamic real array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to increase value for */
   Real             incval              /**< value to increase array index */
   )
{
   return SCIPrealarraySetVal(realarray, set, idx, SCIPrealarrayGetVal(realarray, idx) + incval);
}


/** dynamic array for storing int values */
struct IntArray
{
   MEMHDR*          memhdr;             /**< block memory that stores the vals array */
   int*             vals;               /**< array values */
   int              valssize;           /**< size of vals array */
   int              firstidx;           /**< index of first element in vals array */
   int              minusedidx;         /**< index of first non zero element in vals array */
   int              maxusedidx;         /**< index of last non zero element in vals array */
};

/** creates a dynamic array of int values */
RETCODE SCIPintarrayCreate(
   INTARRAY**       intarray,           /**< pointer to store the int array */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(intarray != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, intarray) );
   (*intarray)->memhdr = memhdr;
   (*intarray)->vals = NULL;
   (*intarray)->valssize = 0;
   (*intarray)->firstidx = -1;
   (*intarray)->minusedidx = INT_MAX;
   (*intarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of real values */
RETCODE SCIPintarrayCopy(
   INTARRAY**       intarray,           /**< pointer to store the copied real array */
   MEMHDR*          memhdr,             /**< block memory */
   INTARRAY*        sourceintarray      /**< dynamic real array to copy */
   )
{
   assert(intarray != NULL);
   assert(sourceintarray != NULL);

   CHECK_OKAY( SCIPintarrayCreate(intarray, memhdr) );
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*intarray)->vals, sourceintarray->vals, sourceintarray->valssize) );
   (*intarray)->valssize = sourceintarray->valssize;
   (*intarray)->firstidx = sourceintarray->firstidx;
   (*intarray)->minusedidx = sourceintarray->minusedidx;
   (*intarray)->maxusedidx = sourceintarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of int values */
RETCODE SCIPintarrayFree(
   INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   assert(intarray != NULL);
   assert(*intarray != NULL);

   freeBlockMemoryArrayNull((*intarray)->memhdr, &(*intarray)->vals, (*intarray)->valssize);
   freeBlockMemory((*intarray)->memhdr, intarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPintarrayExtend(
   INTARRAY*        intarray,           /**< dynamic int array */
   const SET*       set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(intarray != NULL);
   assert(intarray->minusedidx == INT_MAX || intarray->firstidx >= 0);
   assert(intarray->maxusedidx == INT_MIN || intarray->firstidx >= 0);
   assert(intarray->minusedidx == INT_MAX || intarray->minusedidx >= intarray->firstidx);
   assert(intarray->maxusedidx == INT_MIN || intarray->maxusedidx < intarray->firstidx + intarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);
   
   minidx = MIN(minidx, intarray->minusedidx);
   maxidx = MAX(maxidx, intarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   debugMessage("extending intarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > intarray->valssize )
   {
      int* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = SCIPsetCalcMemGrowSize(set, nused);
      ALLOC_OKAY( allocBlockMemoryArray(intarray->memhdr, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( intarray->firstidx != -1 )
      {
         for( i = 0; i < intarray->minusedidx - newfirstidx; ++i )
            newvals[i] = 0;
         copyMemoryArray(&newvals[intarray->minusedidx - newfirstidx],
            &intarray->vals[intarray->minusedidx - intarray->firstidx],
            intarray->maxusedidx - intarray->minusedidx + 1);
         for( i = intarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = 0;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = 0;
      }

      /* free old memory storage, and set the new array parameters */
      freeBlockMemoryArrayNull(intarray->memhdr, &intarray->vals, intarray->valssize);
      intarray->vals = newvals;
      intarray->valssize = newvalssize;
      intarray->firstidx = newfirstidx;
   }
   else if( intarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = intarray->valssize - nused;
      assert(nfree >= 0);
      intarray->firstidx = minidx - nfree/2;
      assert(intarray->firstidx <= minidx);
      assert(maxidx < intarray->firstidx + intarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < intarray->valssize; ++i )
         assert(intarray->vals[i] == 0);
#endif
   }
   else if( minidx < intarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = intarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + intarray->valssize);
      
      if( intarray->minusedidx <= intarray->maxusedidx )
      {
         int shift;

         assert(intarray->firstidx <= intarray->minusedidx);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);

         /* shift used part of array to the right */
         shift = intarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = intarray->maxusedidx - intarray->firstidx; i >= intarray->minusedidx - intarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < intarray->valssize);
            intarray->vals[i + shift] = intarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            intarray->vals[intarray->minusedidx - intarray->firstidx + i] = 0;
      }
      intarray->firstidx = newfirstidx;
   }
   else if( maxidx >= intarray->firstidx + intarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = intarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + intarray->valssize);
      
      if( intarray->minusedidx <= intarray->maxusedidx )
      {
         int shift;

         assert(intarray->firstidx <= intarray->minusedidx);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - intarray->firstidx;
         assert(shift > 0);
         for( i = intarray->minusedidx - intarray->firstidx; i <= intarray->maxusedidx - intarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < intarray->valssize);
            intarray->vals[i - shift] = intarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            intarray->vals[intarray->maxusedidx - intarray->firstidx - i] = 0;
      }
      intarray->firstidx = newfirstidx;
   }

   assert(minidx >= intarray->firstidx);
   assert(maxidx < intarray->firstidx + intarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic int array */
void SCIPintarrayClear(
   INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   debugMessage("clearing intarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx);

   if( intarray->minusedidx <= intarray->maxusedidx )
   {
      int i;
   
      assert(intarray->firstidx <= intarray->minusedidx);
      assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
      assert(intarray->firstidx != -1);
      assert(intarray->valssize > 0);

      /* clear the used part of array */
      for( i = intarray->minusedidx - intarray->firstidx; i <= intarray->maxusedidx - intarray->firstidx; ++i )
         intarray->vals[i] = 0;

      /* mark the array cleared */
      intarray->minusedidx = INT_MAX;
      intarray->maxusedidx = INT_MIN;
   }
   assert(intarray->minusedidx == INT_MAX);
   assert(intarray->maxusedidx == INT_MIN);
}

/** gets value of entry in dynamic array */
int SCIPintarrayGetVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx                 /**< array index to get value for */
   )
{
   assert(intarray != NULL);
   assert(idx >= 0);
   
   if( idx < intarray->minusedidx || idx > intarray->maxusedidx )
      return 0;
   else
   {
      assert(intarray->vals != NULL);
      assert(idx - intarray->firstidx >= 0);
      assert(idx - intarray->firstidx < intarray->valssize);

      return intarray->vals[idx - intarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
RETCODE SCIPintarraySetVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   int              val                 /**< value to set array index to */
   )
{
   assert(intarray != NULL);
   assert(idx >= 0);

   debugMessage("setting intarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %d\n", 
      intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx, idx, val);

   if( !SCIPsetIsZero(set, val) )
   {
      /* extend array to be able to store the index */
      CHECK_OKAY( SCIPintarrayExtend(intarray, set, idx, idx) );
      assert(idx >= intarray->firstidx);
      assert(idx < intarray->firstidx + intarray->valssize);
      
      /* set the array value of the index */
      intarray->vals[idx - intarray->firstidx] = val;

      /* update min/maxusedidx */
      intarray->minusedidx = MIN(intarray->minusedidx, idx);
      intarray->maxusedidx = MAX(intarray->maxusedidx, idx);
   }
   else if( idx >= intarray->firstidx && idx < intarray->firstidx + intarray->valssize )
   {
      /* set the array value of the index to zero */
      intarray->vals[idx - intarray->firstidx] = 0;
      
      /* check, if we can tighten the min/maxusedidx */
      if( idx == intarray->minusedidx )
      {
         assert(intarray->maxusedidx >= 0);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
         do
         {
            intarray->minusedidx++;
         }
         while( intarray->minusedidx <= intarray->maxusedidx
            && SCIPsetIsZero(set, intarray->vals[intarray->minusedidx - intarray->firstidx]) );
         if( intarray->minusedidx > intarray->maxusedidx )
         {
            intarray->minusedidx = INT_MAX;
            intarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == intarray->maxusedidx )
      {
         assert(intarray->minusedidx >= 0);
         assert(intarray->minusedidx < intarray->maxusedidx);
         assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
         do
         {
            intarray->maxusedidx--;
            assert(intarray->minusedidx <= intarray->maxusedidx);
         }
         while( SCIPsetIsZero(set, intarray->vals[intarray->maxusedidx - intarray->firstidx]) );
      }      
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
RETCODE SCIPintarrayIncVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to increase value for */
   int              incval              /**< value to increase array index */
   )
{
   return SCIPintarraySetVal(intarray, set, idx, SCIPintarrayGetVal(intarray, idx) + incval);
}


/** dynamic array for storing bool values */
struct BoolArray
{
   MEMHDR*          memhdr;             /**< block memory that stores the vals array */
   Bool*            vals;               /**< array values */
   int              valssize;           /**< size of vals array */
   int              firstidx;           /**< index of first element in vals array */
   int              minusedidx;         /**< index of first non zero element in vals array */
   int              maxusedidx;         /**< index of last non zero element in vals array */
};

/** creates a dynamic array of bool values */
RETCODE SCIPboolarrayCreate(
   BOOLARRAY**      boolarray,          /**< pointer to store the bool array */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(boolarray != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, boolarray) );
   (*boolarray)->memhdr = memhdr;
   (*boolarray)->vals = NULL;
   (*boolarray)->valssize = 0;
   (*boolarray)->firstidx = -1;
   (*boolarray)->minusedidx = INT_MAX;
   (*boolarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of real values */
RETCODE SCIPboolarrayCopy(
   BOOLARRAY**      boolarray,          /**< pointer to store the copied real array */
   MEMHDR*          memhdr,             /**< block memory */
   BOOLARRAY*       sourceboolarray     /**< dynamic real array to copy */
   )
{
   assert(boolarray != NULL);
   assert(sourceboolarray != NULL);

   CHECK_OKAY( SCIPboolarrayCreate(boolarray, memhdr) );
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*boolarray)->vals, sourceboolarray->vals, sourceboolarray->valssize) );
   (*boolarray)->valssize = sourceboolarray->valssize;
   (*boolarray)->firstidx = sourceboolarray->firstidx;
   (*boolarray)->minusedidx = sourceboolarray->minusedidx;
   (*boolarray)->maxusedidx = sourceboolarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values */
RETCODE SCIPboolarrayFree(
   BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   assert(boolarray != NULL);
   assert(*boolarray != NULL);

   freeBlockMemoryArrayNull((*boolarray)->memhdr, &(*boolarray)->vals, (*boolarray)->valssize);
   freeBlockMemory((*boolarray)->memhdr, boolarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPboolarrayExtend(
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   const SET*       set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(boolarray != NULL);
   assert(boolarray->minusedidx == INT_MAX || boolarray->firstidx >= 0);
   assert(boolarray->maxusedidx == INT_MIN || boolarray->firstidx >= 0);
   assert(boolarray->minusedidx == INT_MAX || boolarray->minusedidx >= boolarray->firstidx);
   assert(boolarray->maxusedidx == INT_MIN || boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);
   
   minidx = MIN(minidx, boolarray->minusedidx);
   maxidx = MAX(maxidx, boolarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   debugMessage("extending boolarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > boolarray->valssize )
   {
      Bool* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = SCIPsetCalcMemGrowSize(set, nused);
      ALLOC_OKAY( allocBlockMemoryArray(boolarray->memhdr, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( boolarray->firstidx != -1 )
      {
         for( i = 0; i < boolarray->minusedidx - newfirstidx; ++i )
            newvals[i] = FALSE;
         copyMemoryArray(&newvals[boolarray->minusedidx - newfirstidx],
            &boolarray->vals[boolarray->minusedidx - boolarray->firstidx],
            boolarray->maxusedidx - boolarray->minusedidx + 1);
         for( i = boolarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = FALSE;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = FALSE;
      }

      /* free old memory storage, and set the new array parameters */
      freeBlockMemoryArrayNull(boolarray->memhdr, &boolarray->vals, boolarray->valssize);
      boolarray->vals = newvals;
      boolarray->valssize = newvalssize;
      boolarray->firstidx = newfirstidx;
   }
   else if( boolarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = boolarray->valssize - nused;
      assert(nfree >= 0);
      boolarray->firstidx = minidx - nfree/2;
      assert(boolarray->firstidx <= minidx);
      assert(maxidx < boolarray->firstidx + boolarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < boolarray->valssize; ++i )
         assert(boolarray->vals[i] == FALSE);
#endif
   }
   else if( minidx < boolarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = boolarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + boolarray->valssize);
      
      if( boolarray->minusedidx <= boolarray->maxusedidx )
      {
         int shift;

         assert(boolarray->firstidx <= boolarray->minusedidx);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);

         /* shift used part of array to the right */
         shift = boolarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = boolarray->maxusedidx - boolarray->firstidx; i >= boolarray->minusedidx - boolarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < boolarray->valssize);
            boolarray->vals[i + shift] = boolarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            boolarray->vals[boolarray->minusedidx - boolarray->firstidx + i] = FALSE;
      }
      boolarray->firstidx = newfirstidx;
   }
   else if( maxidx >= boolarray->firstidx + boolarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = boolarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + boolarray->valssize);
      
      if( boolarray->minusedidx <= boolarray->maxusedidx )
      {
         int shift;

         assert(boolarray->firstidx <= boolarray->minusedidx);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - boolarray->firstidx;
         assert(shift > 0);
         for( i = boolarray->minusedidx - boolarray->firstidx; i <= boolarray->maxusedidx - boolarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < boolarray->valssize);
            boolarray->vals[i - shift] = boolarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            boolarray->vals[boolarray->maxusedidx - boolarray->firstidx - i] = FALSE;
      }
      boolarray->firstidx = newfirstidx;
   }

   assert(minidx >= boolarray->firstidx);
   assert(maxidx < boolarray->firstidx + boolarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic bool array */
void SCIPboolarrayClear(
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   debugMessage("clearing boolarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx);

   if( boolarray->minusedidx <= boolarray->maxusedidx )
   {
      int i;
   
      assert(boolarray->firstidx <= boolarray->minusedidx);
      assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
      assert(boolarray->firstidx != -1);
      assert(boolarray->valssize > 0);

      /* clear the used part of array */
      for( i = boolarray->minusedidx - boolarray->firstidx; i <= boolarray->maxusedidx - boolarray->firstidx; ++i )
         boolarray->vals[i] = FALSE;

      /* mark the array cleared */
      boolarray->minusedidx = INT_MAX;
      boolarray->maxusedidx = INT_MIN;
   }
   assert(boolarray->minusedidx == INT_MAX);
   assert(boolarray->maxusedidx == INT_MIN);
}

/** gets value of entry in dynamic array */
Bool SCIPboolarrayGetVal(
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx                 /**< array index to get value for */
   )
{
   assert(boolarray != NULL);
   assert(idx >= 0);
   
   if( idx < boolarray->minusedidx || idx > boolarray->maxusedidx )
      return FALSE;
   else
   {
      assert(boolarray->vals != NULL);
      assert(idx - boolarray->firstidx >= 0);
      assert(idx - boolarray->firstidx < boolarray->valssize);

      return boolarray->vals[idx - boolarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
RETCODE SCIPboolarraySetVal(
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   Bool             val                 /**< value to set array index to */
   )
{
   assert(boolarray != NULL);
   assert(idx >= 0);

   debugMessage("setting boolarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %d\n", 
      boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx, idx, val);

   if( !SCIPsetIsZero(set, val) )
   {
      /* extend array to be able to store the index */
      CHECK_OKAY( SCIPboolarrayExtend(boolarray, set, idx, idx) );
      assert(idx >= boolarray->firstidx);
      assert(idx < boolarray->firstidx + boolarray->valssize);
      
      /* set the array value of the index */
      boolarray->vals[idx - boolarray->firstidx] = val;

      /* update min/maxusedidx */
      boolarray->minusedidx = MIN(boolarray->minusedidx, idx);
      boolarray->maxusedidx = MAX(boolarray->maxusedidx, idx);
   }
   else if( idx >= boolarray->firstidx && idx < boolarray->firstidx + boolarray->valssize )
   {
      /* set the array value of the index to zero */
      boolarray->vals[idx - boolarray->firstidx] = FALSE;
      
      /* check, if we can tighten the min/maxusedidx */
      if( idx == boolarray->minusedidx )
      {
         assert(boolarray->maxusedidx >= 0);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
         do
         {
            boolarray->minusedidx++;
         }
         while( boolarray->minusedidx <= boolarray->maxusedidx
            && SCIPsetIsZero(set, boolarray->vals[boolarray->minusedidx - boolarray->firstidx]) );
         if( boolarray->minusedidx > boolarray->maxusedidx )
         {
            boolarray->minusedidx = INT_MAX;
            boolarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == boolarray->maxusedidx )
      {
         assert(boolarray->minusedidx >= 0);
         assert(boolarray->minusedidx < boolarray->maxusedidx);
         assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
         do
         {
            boolarray->maxusedidx--;
            assert(boolarray->minusedidx <= boolarray->maxusedidx);
         }
         while( SCIPsetIsZero(set, boolarray->vals[boolarray->maxusedidx - boolarray->firstidx]) );
      }      
   }

   return SCIP_OKAY;
}





/*
 * Sorting algorithms
 */

/** bubble sort of an array of pointers */
void SCIPbsortPtr(
   void**           ptrarray,           /**< pointer array to be sorted */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int actpos;
   int sortpos;
   void* tmpptr;

   assert(len == 0 || ptrarray != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      actpos = firstpos;
      sortpos = firstpos;
      while( actpos < lastpos )
      {
         while( actpos < lastpos && ptrcmp(ptrarray[actpos], ptrarray[actpos+1]) <= 0 )
            actpos++;
         if( actpos >= lastpos )
            break;
         assert( ptrcmp(ptrarray[actpos], ptrarray[actpos+1]) > 0 );
         tmpptr = ptrarray[actpos];
         do
         {
            ptrarray[actpos] = ptrarray[actpos+1];
            actpos++;
         }
         while( actpos < lastpos && ptrcmp(tmpptr, ptrarray[actpos+1]) > 0 );
         ptrarray[actpos] = tmpptr;
         sortpos = actpos;
         actpos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      actpos = lastpos;
      sortpos = lastpos;
      while( actpos > firstpos )
      {
         while( actpos > firstpos && ptrcmp(ptrarray[actpos-1], ptrarray[actpos]) <= 0 )
            actpos--;
         if( actpos <= firstpos )
            break;
         assert( ptrcmp(ptrarray[actpos-1], ptrarray[actpos]) > 0 );
         tmpptr = ptrarray[actpos];
         do
         {
            ptrarray[actpos] = ptrarray[actpos-1];
            actpos--;
         }
         while( actpos > firstpos && ptrcmp(ptrarray[actpos-1], tmpptr) > 0 );
         ptrarray[actpos] = tmpptr;
         sortpos = actpos;
         actpos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of two joint arrays of pointers/Reals, sorted by first array */
void SCIPbsortPtrDbl(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            dblarray,           /**< Real array to be permuted in the same way */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int actpos;
   int sortpos;
   void* tmpptr;
   Real tmpdbl;

   assert(len == 0 || ptrarray != NULL);
   assert(len == 0 || dblarray != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      actpos = firstpos;
      sortpos = firstpos;
      while( actpos < lastpos )
      {
         while( actpos < lastpos && ptrcmp(ptrarray[actpos], ptrarray[actpos+1]) <= 0 )
            actpos++;
         if( actpos >= lastpos )
            break;
         assert( ptrcmp(ptrarray[actpos], ptrarray[actpos+1]) > 0 );
         tmpptr = ptrarray[actpos];
         tmpdbl = dblarray[actpos];
         do
         {
            ptrarray[actpos] = ptrarray[actpos+1];
            dblarray[actpos] = dblarray[actpos+1];
            actpos++;
         }
         while( actpos < lastpos && ptrcmp(tmpptr, ptrarray[actpos+1]) > 0 );
         ptrarray[actpos] = tmpptr;
         dblarray[actpos] = tmpdbl;
         sortpos = actpos;
         actpos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      actpos = lastpos;
      sortpos = lastpos;
      while( actpos > firstpos )
      {
         while( actpos > firstpos && ptrcmp(ptrarray[actpos-1], ptrarray[actpos]) <= 0 )
            actpos--;
         if( actpos <= firstpos )
            break;
         assert( ptrcmp(ptrarray[actpos-1], ptrarray[actpos]) > 0 );
         tmpptr = ptrarray[actpos];
         tmpdbl = dblarray[actpos];
         do
         {
            ptrarray[actpos] = ptrarray[actpos-1];
            dblarray[actpos] = dblarray[actpos-1];
            actpos--;
         }
         while( actpos > firstpos && ptrcmp(ptrarray[actpos-1], tmpptr) > 0 );
         ptrarray[actpos] = tmpptr;
         dblarray[actpos] = tmpdbl;
         sortpos = actpos;
         actpos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of three joint arrays of pointers/Reals/Ints, sorted by first */
void SCIPbsortPtrDblInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            dblarray,           /**< Real array to be permuted in the same way */
   int*             intarray,           /**< int array to be permuted in the same way */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int actpos;
   int sortpos;
   void* tmpptr;
   Real tmpdbl;
   int tmpint;

   assert(len == 0 || ptrarray != NULL);
   assert(len == 0 || dblarray != NULL);
   assert(len == 0 || intarray != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      actpos = firstpos;
      sortpos = firstpos;
      while( actpos < lastpos )
      {
         while( actpos < lastpos && ptrcmp(ptrarray[actpos], ptrarray[actpos+1]) <= 0 )
            actpos++;
         if( actpos >= lastpos )
            break;
         assert( ptrcmp(ptrarray[actpos], ptrarray[actpos+1]) > 0 );
         tmpptr = ptrarray[actpos];
         tmpdbl = dblarray[actpos];
         tmpint = intarray[actpos];
         do
         {
            ptrarray[actpos] = ptrarray[actpos+1];
            dblarray[actpos] = dblarray[actpos+1];
            intarray[actpos] = intarray[actpos+1];
            actpos++;
         }
         while( actpos < lastpos && ptrcmp(tmpptr, ptrarray[actpos+1]) > 0 );
         ptrarray[actpos] = tmpptr;
         dblarray[actpos] = tmpdbl;
         intarray[actpos] = tmpint;
         sortpos = actpos;
         actpos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      actpos = lastpos;
      sortpos = lastpos;
      while( actpos > firstpos )
      {
         while( actpos > firstpos && ptrcmp(ptrarray[actpos-1], ptrarray[actpos]) <= 0 )
            actpos--;
         if( actpos <= firstpos )
            break;
         assert( ptrcmp(ptrarray[actpos-1], ptrarray[actpos]) > 0 );
         tmpptr = ptrarray[actpos];
         tmpdbl = dblarray[actpos];
         tmpint = intarray[actpos];
         do
         {
            ptrarray[actpos] = ptrarray[actpos-1];
            dblarray[actpos] = dblarray[actpos-1];
            intarray[actpos] = intarray[actpos-1];
            actpos--;
         }
         while( actpos > firstpos && ptrcmp(ptrarray[actpos-1], tmpptr) > 0 );
         ptrarray[actpos] = tmpptr;
         dblarray[actpos] = tmpdbl;
         intarray[actpos] = tmpint;
         sortpos = actpos;
         actpos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of four joint arrays of pointers/Reals/Ints/Ints, sorted by first */
void SCIPbsortPtrDblIntInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            dblarray,           /**< Real array to be permuted in the same way */
   int*             intarray1,          /**< first int array to be permuted in the same way */
   int*             intarray2,          /**< second int array to be permuted in the same way */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int actpos;
   int sortpos;
   void* tmpptr;
   Real tmpdbl;
   int tmpint1;
   int tmpint2;

   assert(len == 0 || ptrarray != NULL);
   assert(len == 0 || dblarray != NULL);
   assert(len == 0 || intarray1 != NULL);
   assert(len == 0 || intarray2 != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      actpos = firstpos;
      sortpos = firstpos;
      while( actpos < lastpos )
      {
         while( actpos < lastpos && ptrcmp(ptrarray[actpos], ptrarray[actpos+1]) <= 0 )
            actpos++;
         if( actpos >= lastpos )
            break;
         assert( ptrcmp(ptrarray[actpos], ptrarray[actpos+1]) > 0 );
         tmpptr = ptrarray[actpos];
         tmpdbl = dblarray[actpos];
         tmpint1 = intarray1[actpos];
         tmpint2 = intarray2[actpos];
         do
         {
            ptrarray[actpos] = ptrarray[actpos+1];
            dblarray[actpos] = dblarray[actpos+1];
            intarray1[actpos] = intarray1[actpos+1];
            intarray2[actpos] = intarray2[actpos+1];
            actpos++;
         }
         while( actpos < lastpos && ptrcmp(tmpptr, ptrarray[actpos+1]) > 0 );
         ptrarray[actpos] = tmpptr;
         dblarray[actpos] = tmpdbl;
         intarray1[actpos] = tmpint1;
         intarray2[actpos] = tmpint2;
         sortpos = actpos;
         actpos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      actpos = lastpos;
      sortpos = lastpos;
      while( actpos > firstpos )
      {
         while( actpos > firstpos && ptrcmp(ptrarray[actpos-1], ptrarray[actpos]) <= 0 )
            actpos--;
         if( actpos <= firstpos )
            break;
         assert( ptrcmp(ptrarray[actpos-1], ptrarray[actpos]) > 0 );
         tmpptr = ptrarray[actpos];
         tmpdbl = dblarray[actpos];
         tmpint1 = intarray1[actpos];
         tmpint2 = intarray2[actpos];
         do
         {
            ptrarray[actpos] = ptrarray[actpos-1];
            dblarray[actpos] = dblarray[actpos-1];
            intarray1[actpos] = intarray1[actpos-1];
            intarray2[actpos] = intarray2[actpos-1];
            actpos--;
         }
         while( actpos > firstpos && ptrcmp(ptrarray[actpos-1], tmpptr) > 0 );
         ptrarray[actpos] = tmpptr;
         dblarray[actpos] = tmpdbl;
         intarray1[actpos] = tmpint1;
         intarray2[actpos] = tmpint2;
         sortpos = actpos;
         actpos--;
      }
      firstpos = sortpos+1;
   }
}




/*
 * Numerical methods
 */

/** returns the machine epsilon: the smallest number eps > 0, for which 1.0 + eps > 1.0 */
Real SCIPcalcMachineEpsilon(
   void
   )
{
   Real eps;
   Real lasteps;
   Real one;
   Real onepluseps;

   one = 1.0;
   eps = 1.0;
   do
   {
      lasteps = eps;
      eps /= 2.0;
      onepluseps = one + eps;
   }
   while( onepluseps > one );

   return lasteps;
}

/** calculates the greatest common divisor of the two given values */
Longint SCIPcalcGreComDiv(
   Longint          val1,               /**< first value of greatest common devisor calculation */
   Longint          val2                /**< second value of greatest common devisor calculation */
   )
{
   Longint t;
   Longint gcd;

   assert(val1 >= 0);
   assert(val2 >= 0);

   /* extract all prime factors 2 */
   gcd = 1;
   while( !(val1 & 1) && !(val2 & 1) )
   {
      val1 /= 2;
      val2 /= 2;
      gcd *= 2;
   }

   t = val1 & 1 ? -val2 : val1;
   do
   {
      while( !(t & 1) )
	 t /= 2;

      if( t > 0 )
	 val1 = t;
      else
	 val2 = -t;

      t = val1 - val2;
   }
   while( t != 0 );
   gcd *= val1;

   return gcd;
}

/** calculates the smallest common multiple of the two given values */
Longint SCIPcalcSmaComMul(
   Longint          val1,               /**< first value of greatest common devisor calculation */
   Longint          val2                /**< second value of greatest common devisor calculation */
   )
{
   Longint gcd;

   assert(val1 >= 0);
   assert(val2 >= 0);

   gcd = SCIPcalcGreComDiv(val1, val2);
   
   return val1/gcd * val2;
}




/*
 * File methods
 */

/** returns, whether the given file exists */
Bool SCIPfileExists(
   const char*      filename            /**< file name */
   )
{
   FILE* f;

   f = fopen(filename, "r");
   if( f == NULL )
      return FALSE;

   fclose(f);

   return TRUE;
}
