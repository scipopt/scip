/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: misc.c,v 1.32 2004/12/06 14:11:23 bzfpfend Exp $"

/**@file   misc.c
 * @brief  miscellaneous methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "misc.h"

#include "struct_misc.h"



/*
 * Priority Queue
 */

#define PQ_PARENT(q) (((q)+1)/2-1)
#define PQ_LEFTCHILD(p) (2*(p)+1)
#define PQ_RIGHTCHILD(p) (2*(p)+2)


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

/** creates priority queue */
RETCODE SCIPpqueueCreate(
   PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int              initsize,           /**< initial number of available element slots */
   Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   assert(pqueue != NULL);
   assert(ptrcomp != NULL);

   initsize = MAX(1, initsize);
   sizefac = MAX(1.0, sizefac);

   ALLOC_OKAY( allocMemory(pqueue) );
   (*pqueue)->len = 0;
   (*pqueue)->size = 0;
   (*pqueue)->sizefac = sizefac;
   (*pqueue)->slots = NULL;
   (*pqueue)->ptrcomp = ptrcomp;
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

/** clears the priority queue, but doesn't free the data elements themselves */
void SCIPpqueueClear(
   PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);

   pqueue->len = 0;
}

/** inserts element into priority queue */
RETCODE SCIPpqueueInsert(
   PQUEUE*          pqueue,             /**< priority queue */
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
   while( pos > 0 && (*pqueue->ptrcomp)(elem, pqueue->slots[PQ_PARENT(pos)]) < 0 )
   {
      pqueue->slots[pos] = pqueue->slots[PQ_PARENT(pos)];
      pos = PQ_PARENT(pos);
   }
   pqueue->slots[pos] = elem;

   return SCIP_OKAY;
}

/** removes and returns best element from the priority queue */
void* SCIPpqueueRemove(
   PQUEUE*          pqueue              /**< priority queue */
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
    * of the queue could be placed in the empty slot
    */
   root = pqueue->slots[0];
   last = pqueue->slots[pqueue->len-1];
   pqueue->len--;
   pos = 0;
   while( pos <= PQ_PARENT(pqueue->len-1) )
   {
      childpos = PQ_LEFTCHILD(pos);
      brotherpos = PQ_RIGHTCHILD(pos);
      if( brotherpos <= pqueue->len && (*pqueue->ptrcomp)(pqueue->slots[brotherpos], pqueue->slots[childpos]) < 0 )
         childpos = brotherpos;
      if( (*pqueue->ptrcomp)(last, pqueue->slots[childpos]) <= 0 )
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
   PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   if( pqueue->len == 0 )
      return NULL;

   return pqueue->slots[0];
}

/** returns the number of elements in the queue */
int SCIPpqueueNElems(
   PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   return pqueue->len;
}

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
void** SCIPpqueueElems(
   PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   return pqueue->slots;
}




/*
 * Hash Table
 */

/** appends element to the hash list */
static
RETCODE hashtablelistAppend(
   HASHTABLELIST**  hashtablelist,      /**< pointer to hash list */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to append to the list */
   )
{
   HASHTABLELIST* newlist;

   assert(hashtablelist != NULL);
   assert(memhdr != NULL);
   assert(element != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, &newlist) );
   newlist->element = element;
   newlist->next = *hashtablelist;
   *hashtablelist = newlist;

   return SCIP_OKAY;
}

/** frees a hash list entry and all its successors */
static
void hashtablelistFree(
   HASHTABLELIST**  hashtablelist,      /**< pointer to hash list to free */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   HASHTABLELIST* list;
   HASHTABLELIST* nextlist;

   assert(hashtablelist != NULL);
   assert(memhdr != NULL);
   
   list = *hashtablelist;
   while( list != NULL )
   {
      nextlist = list->next;
      freeBlockMemory(memhdr, &list);
      list = nextlist;
   }

   *hashtablelist = NULL;
}

/** finds hash list entry pointing to element with given key in the hash list, returns NULL if not found */
static
HASHTABLELIST* hashtablelistFind(
   HASHTABLELIST*   hashtablelist,      /**< hash list */
   DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   unsigned int     keyval,             /**< hash value of key */
   void*            key                 /**< key to retrieve */
   )
{
   unsigned int currentkeyval;
   void* currentkey;

   assert(hashkeyeq != NULL);
   assert(key != NULL);

   while( hashtablelist != NULL )
   {
      currentkey = hashgetkey(hashtablelist->element);
      currentkeyval = hashkeyval(currentkey);
      if( currentkeyval == keyval && hashkeyeq(currentkey, key) )
         return hashtablelist;
      hashtablelist = hashtablelist->next;
   }

   return NULL;
}

/** retrieves element with given key from the hash list, or NULL */
static
void* hashtablelistRetrieve(
   HASHTABLELIST*   hashtablelist,      /**< hash list */
   DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   unsigned int     keyval,             /**< hash value of key */
   void*            key                 /**< key to retrieve */
   )
{
   HASHTABLELIST* h;

   /* find hash list entry */
   h = hashtablelistFind(hashtablelist, hashgetkey, hashkeyeq, hashkeyval, keyval, key);

   /* return element */
   if( h != NULL )
      return h->element;
   else
      return NULL;
}

/** removes element from the hash list */
static
RETCODE hashtablelistRemove(
   HASHTABLELIST**  hashtablelist,      /**< pointer to hash list */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to remove from the list */
   )
{
   HASHTABLELIST* nextlist;

   assert(hashtablelist != NULL);
   assert(memhdr != NULL);
   assert(element != NULL);

   while( *hashtablelist != NULL && (*hashtablelist)->element != element )
   {
      hashtablelist = &(*hashtablelist)->next;
   }
   if( *hashtablelist != NULL )
   {
      nextlist = (*hashtablelist)->next;
      freeBlockMemory(memhdr, hashtablelist);
      *hashtablelist = nextlist;
   }
   else
   {
      errorMessage("element not found in the hash table\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

   
/** creates a hash table */
RETCODE SCIPhashtableCreate(
   HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   MEMHDR*          memhdr,             /**< block memory used to store hash table entries */
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
   (*hashtable)->memhdr = memhdr;
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
   HASHTABLE**      hashtable           /**< pointer to the hash table */
   )
{
   int i;

   assert(hashtable != NULL);
   assert(*hashtable != NULL);

   /* free hash lists */
   for( i = 0; i < (*hashtable)->nlists; ++i )
      hashtablelistFree(&(*hashtable)->lists[i], (*hashtable)->memhdr);

   /* free main hast table data structure */
   freeMemoryArray(&(*hashtable)->lists);
   freeMemory(hashtable);
}

/** inserts element in hash table (multiple inserts of same element possible) */
RETCODE SCIPhashtableInsert(
   HASHTABLE*       hashtable,          /**< hash table */
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
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(element);
   keyval = hashtable->hashkeyval(key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   /* append element to the list at the hash position */
   CHECK_OKAY( hashtablelistAppend(&hashtable->lists[hashval], hashtable->memhdr, element) );
   
   return SCIP_OKAY;
}

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
RETCODE SCIPhashtableSafeInsert(
   HASHTABLE*       hashtable,          /**< hash table */
   void*            element             /**< element to insert into the table */
   )
{
   assert(hashtable != NULL);
   assert(hashtable->hashgetkey != NULL);

   /* check, if key is already existing */
   if( SCIPhashtableRetrieve(hashtable, hashtable->hashgetkey(element)) != NULL )
      return SCIP_KEYALREADYEXISTING;

   /* insert element in hash table */
   CHECK_OKAY( SCIPhashtableInsert(hashtable, element) );
   
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
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   return hashtablelistRetrieve(hashtable->lists[hashval], hashtable->hashgetkey, hashtable->hashkeyeq, 
      hashtable->hashkeyval, keyval, key);
}

/** removes existing element from the hash table */
RETCODE SCIPhashtableRemove(
   HASHTABLE*       hashtable,          /**< hash table */
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
   assert(element != NULL);

   /* get the hash key and its hash value */
   key = hashtable->hashgetkey(element);
   keyval = hashtable->hashkeyval(key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   /* append element to the list at the hash position */
   CHECK_OKAY( hashtablelistRemove(&hashtable->lists[hashval], hashtable->memhdr, element) );
   
   return SCIP_OKAY;
}

/** prints statistics about hash table usage */
void SCIPhashtablePrintStatistics(
   HASHTABLE*       hashtable           /**< hash table */
   )
{
   HASHTABLELIST* hashtablelist;
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
      hashtablelist = hashtable->lists[i];
      if( hashtablelist != NULL )
      {
         usedslots++;
         slotsize = 0;
         while( hashtablelist != NULL )
         {
            slotsize++;
            hashtablelist = hashtablelist->next;
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
{  /*lint --e{715}*/
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
{  /*lint --e{715}*/
   const char* string = (const char*)key;
   unsigned int sum;
   unsigned int val;
   int i;

   sum = 0;
   i = 0;
   while( *string != '\0' )
   {
      assert(0 <= i && i < 28);
      val = (unsigned int)(*string); /*lint !e571*/
      val <<= strhashshift[i];
      sum += val;
      i++;
      i %= 28;
      string++;
   }

   return sum;
}




/*
 * Hash Map
 */

/** appends origin->image pair to the hash list */
static
RETCODE hashmaplistAppend(
   HASHMAPLIST**    hashmaplist,        /**< pointer to hash list */
   MEMHDR*          memhdr,             /**< block memory */
   void*            origin,             /**< origin of the mapping origin -> image */
   void*            image               /**< image of the mapping origin -> image */
   )
{
   HASHMAPLIST* newlist;

   assert(hashmaplist != NULL);
   assert(memhdr != NULL);
   assert(origin != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, &newlist) );
   newlist->origin = origin;
   newlist->image = image;
   newlist->next = *hashmaplist;
   *hashmaplist = newlist;

   return SCIP_OKAY;
}

/** frees a hash list entry and all its successors */
static
void hashmaplistFree(
   HASHMAPLIST**    hashmaplist,        /**< pointer to hash list to free */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   HASHMAPLIST* list;
   HASHMAPLIST* nextlist;

   assert(hashmaplist != NULL);
   assert(memhdr != NULL);
   
   list = *hashmaplist;
   while( list != NULL )
   {
      nextlist = list->next;
      freeBlockMemory(memhdr, &list);
      list = nextlist;
   }

   *hashmaplist = NULL;
}

/** finds hash list entry pointing to given origin in the hash list, returns NULL if not found */
static
HASHMAPLIST* hashmaplistFind(
   HASHMAPLIST*     hashmaplist,        /**< hash list */
   void*            origin              /**< origin to find */
   )
{
   assert(origin != NULL);

   while( hashmaplist != NULL )
   {
      if( hashmaplist->origin == origin )
         return hashmaplist;
      hashmaplist = hashmaplist->next;
   }

   return NULL;
}

/** retrieves image of given origin from the hash list, or NULL */
static
void* hashmaplistGetImage(
   HASHMAPLIST*     hashmaplist,        /**< hash list */
   void*            origin              /**< origin to retrieve image for */
   )
{
   HASHMAPLIST* h;

   /* find hash list entry */
   h = hashmaplistFind(hashmaplist, origin);

   /* return image */
   if( h != NULL )
      return h->image;
   else
      return NULL;
}

/** sets image for given origin in the hash list, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
static
RETCODE hashmaplistSetImage(
   HASHMAPLIST**    hashmaplist,        /**< pointer to hash list */
   MEMHDR*          memhdr,             /**< block memory */
   void*            origin,             /**< origin to set image for */
   void*            image               /**< new image for origin */
   )
{
   HASHMAPLIST* h;

   /* find hash list entry */
   h = hashmaplistFind(*hashmaplist, origin);

   /* set image or add origin->image pair */
   if( h != NULL )
      h->image = image;
   else
   {
      CHECK_OKAY( hashmaplistAppend(hashmaplist, memhdr, origin, image) );
   }

   return SCIP_OKAY;
}

/** removes origin->image pair from the hash list */
static
RETCODE hashmaplistRemove(
   HASHMAPLIST**    hashmaplist,        /**< pointer to hash list */
   MEMHDR*          memhdr,             /**< block memory */
   void*            origin              /**< origin to remove from the list */
   )
{
   HASHMAPLIST* nextlist;

   assert(hashmaplist != NULL);
   assert(memhdr != NULL);
   assert(origin != NULL);

   while( *hashmaplist != NULL && (*hashmaplist)->origin != origin )
   {
      hashmaplist = &(*hashmaplist)->next;
   }
   if( *hashmaplist != NULL )
   {
      nextlist = (*hashmaplist)->next;
      freeBlockMemory(memhdr, hashmaplist);
      *hashmaplist = nextlist;
   }
   else
   {
      errorMessage("origin not found in the hash map\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}


/** creates a hash map mapping pointers to pointers */
RETCODE SCIPhashmapCreate(
   HASHMAP**        hashmap,            /**< pointer to store the created hash map */
   MEMHDR*          memhdr,             /**< block memory used to store hash map entries */
   int              mapsize             /**< size of the hash map */
   )
{
   int i;

   assert(hashmap != NULL);
   assert(mapsize > 0);

   ALLOC_OKAY( allocMemory(hashmap) );
   ALLOC_OKAY( allocMemoryArray(&(*hashmap)->lists, mapsize) );
   (*hashmap)->memhdr = memhdr;
   (*hashmap)->nlists = mapsize;

   /* initialize hash lists */
   for( i = 0; i < mapsize; ++i )
      (*hashmap)->lists[i] = NULL;

   return SCIP_OKAY;
}

/** frees the hash map */
void SCIPhashmapFree(
   HASHMAP**        hashmap             /**< pointer to the hash map */
   )
{
   int i;

   assert(hashmap != NULL);
   assert(*hashmap != NULL);

   /* free hash lists */
   for( i = 0; i < (*hashmap)->nlists; ++i )
      hashmaplistFree(&(*hashmap)->lists[i], (*hashmap)->memhdr);

   /* free main hast map data structure */
   freeMemoryArray(&(*hashmap)->lists);
   freeMemory(hashmap);
}

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
RETCODE SCIPhashmapInsert(
   HASHMAP*         hashmap,            /**< hash map */
   void*            origin,             /**< origin to set image for */
   void*            image               /**< new image for origin */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)origin % hashmap->nlists;

   /* append origin->image pair to the list at the hash position */
   CHECK_OKAY( hashmaplistAppend(&hashmap->lists[hashval], hashmap->memhdr, origin, image) );
   
   return SCIP_OKAY;
}

/** retrieves image of given origin from the hash map, or NULL if no image exists */
void* SCIPhashmapGetImage(
   HASHMAP*         hashmap,            /**< hash map */
   void*            origin              /**< origin to retrieve image for */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)origin % hashmap->nlists;

   /* get image for origin from hash list */
   return hashmaplistGetImage(hashmap->lists[hashval], origin);
}

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
RETCODE SCIPhashmapSetImage(
   HASHMAP*         hashmap,            /**< hash map */
   void*            origin,             /**< origin to set image for */
   void*            image               /**< new image for origin */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)origin % hashmap->nlists;

   /* set image for origin in hash list */
   CHECK_OKAY( hashmaplistSetImage(&hashmap->lists[hashval], hashmap->memhdr, origin, image) );
   
   return SCIP_OKAY;
}

/** removes existing origin->image pair from the hash map */
RETCODE SCIPhashmapRemove(
   HASHMAP*         hashmap,            /**< hash map */
   void*            origin              /**< origin to remove from the list */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)origin % hashmap->nlists;

   /* append element to the list at the hash position */
   CHECK_OKAY( hashmaplistRemove(&hashmap->lists[hashval], hashmap->memhdr, origin) );
   
   return SCIP_OKAY;
}

/** prints statistics about hash map usage */
void SCIPhashmapPrintStatistics(
   HASHMAP*         hashmap             /**< hash map */
   )
{
   HASHMAPLIST* hashmaplist;
   int usedslots;
   int maxslotsize;
   int sumslotsize;
   int slotsize;
   int i;

   assert(hashmap != NULL);

   usedslots = 0;
   maxslotsize = 0;
   sumslotsize = 0;
   for( i = 0; i < hashmap->nlists; ++i )
   {
      hashmaplist = hashmap->lists[i];
      if( hashmaplist != NULL )
      {
         usedslots++;
         slotsize = 0;
         while( hashmaplist != NULL )
         {
            slotsize++;
            hashmaplist = hashmaplist->next;
         }
         maxslotsize = MAX(maxslotsize, slotsize);
         sumslotsize += slotsize;
      }
   }

   printf("%d hash entries, used %d/%d slots (%.1f%%)",
      sumslotsize, usedslots, hashmap->nlists, 100.0*(Real)usedslots/(Real)(hashmap->nlists));
   if( usedslots > 0 )
      printf(", avg. %.1f entries/used slot, max. %d entries in slot", (Real)sumslotsize/(Real)usedslots, maxslotsize);
   printf("\n");
}




/*
 * Dynamic Arrays
 */

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
   if( sourcerealarray->valssize > 0 )
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*realarray)->vals, sourcerealarray->vals,
                     sourcerealarray->valssize) );
   }
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
   SET*             set,                /**< global SCIP settings */
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
RETCODE SCIPrealarrayClear(
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

   return SCIP_OKAY;
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
   SET*             set,                /**< global SCIP settings */
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
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to increase value for */
   Real             incval              /**< value to increase array index */
   )
{
   return SCIPrealarraySetVal(realarray, set, idx, SCIPrealarrayGetVal(realarray, idx) + incval);
}

/** returns the minimal index of all stored non-zero elements */
int SCIPrealarrayGetMinIdx(
   REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   return realarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPrealarrayGetMaxIdx(
   REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   return realarray->maxusedidx;
}

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

/** creates a copy of a dynamic array of int values */
RETCODE SCIPintarrayCopy(
   INTARRAY**       intarray,           /**< pointer to store the copied int array */
   MEMHDR*          memhdr,             /**< block memory */
   INTARRAY*        sourceintarray      /**< dynamic int array to copy */
   )
{
   assert(intarray != NULL);
   assert(sourceintarray != NULL);

   CHECK_OKAY( SCIPintarrayCreate(intarray, memhdr) );
   if( sourceintarray->valssize > 0 )
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*intarray)->vals, sourceintarray->vals, sourceintarray->valssize) );
   }
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
   SET*             set,                /**< global SCIP settings */
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
RETCODE SCIPintarrayClear(
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

   return SCIP_OKAY;
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
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   int              val                 /**< value to set array index to */
   )
{
   assert(intarray != NULL);
   assert(idx >= 0);

   debugMessage("setting intarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %d\n", 
      intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx, idx, val);

   if( val != 0 )
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
            && intarray->vals[intarray->minusedidx - intarray->firstidx] == 0 );
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
         while( intarray->vals[intarray->maxusedidx - intarray->firstidx] == 0 );
      }      
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
RETCODE SCIPintarrayIncVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to increase value for */
   int              incval              /**< value to increase array index */
   )
{
   return SCIPintarraySetVal(intarray, set, idx, SCIPintarrayGetVal(intarray, idx) + incval);
}

/** returns the minimal index of all stored non-zero elements */
int SCIPintarrayGetMinIdx(
   INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   return intarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPintarrayGetMaxIdx(
   INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   return intarray->maxusedidx;
}


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

/** creates a copy of a dynamic array of bool values */
RETCODE SCIPboolarrayCopy(
   BOOLARRAY**      boolarray,          /**< pointer to store the copied bool array */
   MEMHDR*          memhdr,             /**< block memory */
   BOOLARRAY*       sourceboolarray     /**< dynamic bool array to copy */
   )
{
   assert(boolarray != NULL);
   assert(sourceboolarray != NULL);

   CHECK_OKAY( SCIPboolarrayCreate(boolarray, memhdr) );
   if( sourceboolarray->valssize > 0 )
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*boolarray)->vals, sourceboolarray->vals, 
                     sourceboolarray->valssize) );
   }
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
   SET*             set,                /**< global SCIP settings */
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
RETCODE SCIPboolarrayClear(
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

   return SCIP_OKAY;
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
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   Bool             val                 /**< value to set array index to */
   )
{
   assert(boolarray != NULL);
   assert(idx >= 0);

   debugMessage("setting boolarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %d\n", 
      boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx, idx, val);

   if( val != FALSE )
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
            && boolarray->vals[boolarray->minusedidx - boolarray->firstidx] == FALSE );
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
         while( boolarray->vals[boolarray->maxusedidx - boolarray->firstidx] == FALSE );
      }      
   }

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPboolarrayGetMinIdx(
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   return boolarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPboolarrayGetMaxIdx(
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   return boolarray->maxusedidx;
}


/** creates a dynamic array of pointer values */
RETCODE SCIPptrarrayCreate(
   PTRARRAY**       ptrarray,           /**< pointer to store the ptr array */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(ptrarray != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, ptrarray) );
   (*ptrarray)->memhdr = memhdr;
   (*ptrarray)->vals = NULL;
   (*ptrarray)->valssize = 0;
   (*ptrarray)->firstidx = -1;
   (*ptrarray)->minusedidx = INT_MAX;
   (*ptrarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of pointer values */
RETCODE SCIPptrarrayCopy(
   PTRARRAY**       ptrarray,           /**< pointer to store the copied ptr array */
   MEMHDR*          memhdr,             /**< block memory */
   PTRARRAY*        sourceptrarray      /**< dynamic ptr array to copy */
   )
{
   assert(ptrarray != NULL);
   assert(sourceptrarray != NULL);

   CHECK_OKAY( SCIPptrarrayCreate(ptrarray, memhdr) );
   if( sourceptrarray->valssize > 0 )
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*ptrarray)->vals, sourceptrarray->vals, sourceptrarray->valssize) );
   }
   (*ptrarray)->valssize = sourceptrarray->valssize;
   (*ptrarray)->firstidx = sourceptrarray->firstidx;
   (*ptrarray)->minusedidx = sourceptrarray->minusedidx;
   (*ptrarray)->maxusedidx = sourceptrarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of pointer values */
RETCODE SCIPptrarrayFree(
   PTRARRAY**       ptrarray            /**< pointer to the ptr array */
   )
{
   assert(ptrarray != NULL);
   assert(*ptrarray != NULL);

   freeBlockMemoryArrayNull((*ptrarray)->memhdr, &(*ptrarray)->vals, (*ptrarray)->valssize);
   freeBlockMemory((*ptrarray)->memhdr, ptrarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPptrarrayExtend(
   PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   SET*             set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   int nused;
   int nfree;
   int newfirstidx;
   int i;

   assert(ptrarray != NULL);
   assert(ptrarray->minusedidx == INT_MAX || ptrarray->firstidx >= 0);
   assert(ptrarray->maxusedidx == INT_MIN || ptrarray->firstidx >= 0);
   assert(ptrarray->minusedidx == INT_MAX || ptrarray->minusedidx >= ptrarray->firstidx);
   assert(ptrarray->maxusedidx == INT_MIN || ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
   assert(0 <= minidx);
   assert(minidx <= maxidx);
   
   minidx = MIN(minidx, ptrarray->minusedidx);
   maxidx = MAX(maxidx, ptrarray->maxusedidx);
   assert(0 <= minidx);
   assert(minidx <= maxidx);

   debugMessage("extending ptrarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > ptrarray->valssize )
   {
      void** newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = SCIPsetCalcMemGrowSize(set, nused);
      ALLOC_OKAY( allocBlockMemoryArray(ptrarray->memhdr, &newvals, newvalssize) );
      nfree = newvalssize - nused;
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + newvalssize);

      /* initialize memory array by copying old values and setting new values to zero */
      if( ptrarray->firstidx != -1 )
      {
         for( i = 0; i < ptrarray->minusedidx - newfirstidx; ++i )
            newvals[i] = NULL;
         copyMemoryArray(&newvals[ptrarray->minusedidx - newfirstidx],
            &ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx],
            ptrarray->maxusedidx - ptrarray->minusedidx + 1);
         for( i = ptrarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = NULL;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = NULL;
      }

      /* free old memory storage, and set the new array parameters */
      freeBlockMemoryArrayNull(ptrarray->memhdr, &ptrarray->vals, ptrarray->valssize);
      ptrarray->vals = newvals;
      ptrarray->valssize = newvalssize;
      ptrarray->firstidx = newfirstidx;
   }
   else if( ptrarray->firstidx == -1 )
   {
      /* a sufficiently large memory storage exists, but it was cleared */
      nfree = ptrarray->valssize - nused;
      assert(nfree >= 0);
      ptrarray->firstidx = minidx - nfree/2;
      assert(ptrarray->firstidx <= minidx);
      assert(maxidx < ptrarray->firstidx + ptrarray->valssize);
#ifndef NDEBUG
      for( i = 0; i < ptrarray->valssize; ++i )
         assert(ptrarray->vals[i] == NULL);
#endif
   }
   else if( minidx < ptrarray->firstidx )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the right */
      nfree = ptrarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + ptrarray->valssize);
      
      if( ptrarray->minusedidx <= ptrarray->maxusedidx )
      {
         int shift;

         assert(ptrarray->firstidx <= ptrarray->minusedidx);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);

         /* shift used part of array to the right */
         shift = ptrarray->firstidx - newfirstidx;
         assert(shift > 0);
         for( i = ptrarray->maxusedidx - ptrarray->firstidx; i >= ptrarray->minusedidx - ptrarray->firstidx; --i )
         {
            assert(0 <= i + shift && i + shift < ptrarray->valssize);
            ptrarray->vals[i + shift] = ptrarray->vals[i];
         }
         /* clear the formerly used head of the array */
         for( i = 0; i < shift; ++i )
            ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx + i] = NULL;
      }
      ptrarray->firstidx = newfirstidx;
   }
   else if( maxidx >= ptrarray->firstidx + ptrarray->valssize )
   {
      /* a sufficiently large memory storage exists, but it has to be shifted to the left */
      nfree = ptrarray->valssize - nused;
      assert(nfree >= 0);
      newfirstidx = minidx - nfree/2;
      newfirstidx = MAX(newfirstidx, 0);
      assert(newfirstidx <= minidx);
      assert(maxidx < newfirstidx + ptrarray->valssize);
      
      if( ptrarray->minusedidx <= ptrarray->maxusedidx )
      {
         int shift;

         assert(ptrarray->firstidx <= ptrarray->minusedidx);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);

         /* shift used part of array to the left */
         shift = newfirstidx - ptrarray->firstidx;
         assert(shift > 0);
         for( i = ptrarray->minusedidx - ptrarray->firstidx; i <= ptrarray->maxusedidx - ptrarray->firstidx; ++i )
         {
            assert(0 <= i - shift && i - shift < ptrarray->valssize);
            ptrarray->vals[i - shift] = ptrarray->vals[i];
         }
         /* clear the formerly used tail of the array */
         for( i = 0; i < shift; ++i )
            ptrarray->vals[ptrarray->maxusedidx - ptrarray->firstidx - i] = NULL;
      }
      ptrarray->firstidx = newfirstidx;
   }

   assert(minidx >= ptrarray->firstidx);
   assert(maxidx < ptrarray->firstidx + ptrarray->valssize);

   return SCIP_OKAY;
}

/** clears a dynamic pointer array */
RETCODE SCIPptrarrayClear(
   PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   debugMessage("clearing ptrarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx);

   if( ptrarray->minusedidx <= ptrarray->maxusedidx )
   {
      int i;
   
      assert(ptrarray->firstidx <= ptrarray->minusedidx);
      assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
      assert(ptrarray->firstidx != -1);
      assert(ptrarray->valssize > 0);

      /* clear the used part of array */
      for( i = ptrarray->minusedidx - ptrarray->firstidx; i <= ptrarray->maxusedidx - ptrarray->firstidx; ++i )
         ptrarray->vals[i] = NULL;

      /* mark the array cleared */
      ptrarray->minusedidx = INT_MAX;
      ptrarray->maxusedidx = INT_MIN;
   }
   assert(ptrarray->minusedidx == INT_MAX);
   assert(ptrarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
void* SCIPptrarrayGetVal(
   PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int              idx                 /**< array index to get value for */
   )
{
   assert(ptrarray != NULL);
   assert(idx >= 0);
   
   if( idx < ptrarray->minusedidx || idx > ptrarray->maxusedidx )
      return NULL;
   else
   {
      assert(ptrarray->vals != NULL);
      assert(idx - ptrarray->firstidx >= 0);
      assert(idx - ptrarray->firstidx < ptrarray->valssize);

      return ptrarray->vals[idx - ptrarray->firstidx];
   }
}

/** sets value of entry in dynamic array */
RETCODE SCIPptrarraySetVal(
   PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   SET*             set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   void*            val                 /**< value to set array index to */
   )
{
   assert(ptrarray != NULL);
   assert(idx >= 0);

   debugMessage("setting ptrarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %p\n", 
      ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx, idx, val);

   if( val != NULL )
   {
      /* extend array to be able to store the index */
      CHECK_OKAY( SCIPptrarrayExtend(ptrarray, set, idx, idx) );
      assert(idx >= ptrarray->firstidx);
      assert(idx < ptrarray->firstidx + ptrarray->valssize);
      
      /* set the array value of the index */
      ptrarray->vals[idx - ptrarray->firstidx] = val;

      /* update min/maxusedidx */
      ptrarray->minusedidx = MIN(ptrarray->minusedidx, idx);
      ptrarray->maxusedidx = MAX(ptrarray->maxusedidx, idx);
   }
   else if( idx >= ptrarray->firstidx && idx < ptrarray->firstidx + ptrarray->valssize )
   {
      /* set the array value of the index to zero */
      ptrarray->vals[idx - ptrarray->firstidx] = NULL;
      
      /* check, if we can tighten the min/maxusedidx */
      if( idx == ptrarray->minusedidx )
      {
         assert(ptrarray->maxusedidx >= 0);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
         do
         {
            ptrarray->minusedidx++;
         }
         while( ptrarray->minusedidx <= ptrarray->maxusedidx
            && ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx] == NULL );
         if( ptrarray->minusedidx > ptrarray->maxusedidx )
         {
            ptrarray->minusedidx = INT_MAX;
            ptrarray->maxusedidx = INT_MIN;
         }
      }
      else if( idx == ptrarray->maxusedidx )
      {
         assert(ptrarray->minusedidx >= 0);
         assert(ptrarray->minusedidx < ptrarray->maxusedidx);
         assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
         do
         {
            ptrarray->maxusedidx--;
            assert(ptrarray->minusedidx <= ptrarray->maxusedidx);
         }
         while( ptrarray->vals[ptrarray->maxusedidx - ptrarray->firstidx] == NULL );
      }      
   }

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPptrarrayGetMinIdx(
   PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   return ptrarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPptrarrayGetMaxIdx(
   PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   return ptrarray->maxusedidx;
}




/*
 * Sorting algorithms
 */

/** bubble sort an indexed element set, resulting in a permutation index array */
void SCIPbsort(
   void*            dataptr,            /**< pointer to data field that is given to the external compare method */
   int              len,                /**< number of elements to be sorted (valid index range) */
   DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   int*             indarray            /**< pointer to store the sorted index array */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   int tmpind;

   assert(indcomp != NULL);
   assert(len == 0 || indarray != NULL);

   /* create identity permutation */
   for( pos = 0; pos < len; ++pos )
      indarray[pos] = pos;

   /* bubble sort index array */
   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && indcomp(dataptr, indarray[pos], indarray[pos+1]) <= 0 )
            pos++;
         if( pos >= lastpos )
            break;
         assert( indcomp(dataptr, indarray[pos], indarray[pos+1]) > 0 );
         tmpind = indarray[pos];
         do
         {
            indarray[pos] = indarray[pos+1];
            pos++;
         }
         while( pos < lastpos && indcomp(dataptr, tmpind, indarray[pos+1]) > 0 );
         indarray[pos] = tmpind;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && indcomp(dataptr, indarray[pos-1], indarray[pos]) <= 0 )
            pos--;
         if( pos <= firstpos )
            break;
         assert( indcomp(dataptr, indarray[pos-1], indarray[pos]) > 0 );
         tmpind = indarray[pos];
         do
         {
            indarray[pos] = indarray[pos-1];
            pos--;
         }
         while( pos > firstpos && indcomp(dataptr, indarray[pos-1], tmpind) > 0 );
         indarray[pos] = tmpind;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of an array of pointers */
void SCIPbsortPtr(
   void**           ptrarray,           /**< pointer array to be sorted */
   int              len,                /**< length of array */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   void* tmpptr;

   assert(len == 0 || ptrarray != NULL);
   assert(ptrcomp != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && ptrcomp(ptrarray[pos], ptrarray[pos+1]) <= 0 )
            pos++;
         if( pos >= lastpos )
            break;
         assert( ptrcomp(ptrarray[pos], ptrarray[pos+1]) > 0 );
         tmpptr = ptrarray[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos+1];
            pos++;
         }
         while( pos < lastpos && ptrcomp(tmpptr, ptrarray[pos+1]) > 0 );
         ptrarray[pos] = tmpptr;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], ptrarray[pos]) <= 0 )
            pos--;
         if( pos <= firstpos )
            break;
         assert( ptrcomp(ptrarray[pos-1], ptrarray[pos]) > 0 );
         tmpptr = ptrarray[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos-1];
            pos--;
         }
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], tmpptr) > 0 );
         ptrarray[pos] = tmpptr;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of two joint arrays of pointers/Reals, sorted by first array */
void SCIPbsortPtrReal(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            realarray,          /**< Real array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   void* tmpptr;
   Real tmpreal;

   assert(len == 0 || ptrarray != NULL);
   assert(len == 0 || realarray != NULL);
   assert(ptrcomp != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && ptrcomp(ptrarray[pos], ptrarray[pos+1]) <= 0 )
            pos++;
         if( pos >= lastpos )
            break;
         assert( ptrcomp(ptrarray[pos], ptrarray[pos+1]) > 0 );
         tmpptr = ptrarray[pos];
         tmpreal = realarray[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos+1];
            realarray[pos] = realarray[pos+1];
            pos++;
         }
         while( pos < lastpos && ptrcomp(tmpptr, ptrarray[pos+1]) > 0 );
         ptrarray[pos] = tmpptr;
         realarray[pos] = tmpreal;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], ptrarray[pos]) <= 0 )
            pos--;
         if( pos <= firstpos )
            break;
         assert( ptrcomp(ptrarray[pos-1], ptrarray[pos]) > 0 );
         tmpptr = ptrarray[pos];
         tmpreal = realarray[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos-1];
            realarray[pos] = realarray[pos-1];
            pos--;
         }
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], tmpptr) > 0 );
         ptrarray[pos] = tmpptr;
         realarray[pos] = tmpreal;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of two joint arrays of pointers/ints, sorted by first array */
void SCIPbsortPtrInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   int*             intarray,           /**< int array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   void* tmpptr;
   int tmpint;

   assert(len == 0 || ptrarray != NULL);
   assert(len == 0 || intarray != NULL);
   assert(ptrcomp != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && ptrcomp(ptrarray[pos], ptrarray[pos+1]) <= 0 )
            pos++;
         if( pos >= lastpos )
            break;
         assert( ptrcomp(ptrarray[pos], ptrarray[pos+1]) > 0 );
         tmpptr = ptrarray[pos];
         tmpint = intarray[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos+1];
            intarray[pos] = intarray[pos+1];
            pos++;
         }
         while( pos < lastpos && ptrcomp(tmpptr, ptrarray[pos+1]) > 0 );
         ptrarray[pos] = tmpptr;
         intarray[pos] = tmpint;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], ptrarray[pos]) <= 0 )
            pos--;
         if( pos <= firstpos )
            break;
         assert( ptrcomp(ptrarray[pos-1], ptrarray[pos]) > 0 );
         tmpptr = ptrarray[pos];
         tmpint = intarray[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos-1];
            intarray[pos] = intarray[pos-1];
            pos--;
         }
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], tmpptr) > 0 );
         ptrarray[pos] = tmpptr;
         intarray[pos] = tmpint;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of three joint arrays of pointers/ints/ints, sorted by first array */
void SCIPbsortPtrIntInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   int*             intarray1,          /**< first int array to be permuted in the same way */
   int*             intarray2,          /**< second int array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   void* tmpptr;
   int tmpint1;
   int tmpint2;

   assert(len == 0 || ptrarray != NULL);
   assert(len == 0 || intarray1 != NULL);
   assert(len == 0 || intarray2 != NULL);
   assert(ptrcomp != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && ptrcomp(ptrarray[pos], ptrarray[pos+1]) <= 0 )
            pos++;
         if( pos >= lastpos )
            break;
         assert( ptrcomp(ptrarray[pos], ptrarray[pos+1]) > 0 );
         tmpptr = ptrarray[pos];
         tmpint1 = intarray1[pos];
         tmpint2 = intarray2[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos+1];
            intarray1[pos] = intarray1[pos+1];
            intarray2[pos] = intarray2[pos+1];
            pos++;
         }
         while( pos < lastpos && ptrcomp(tmpptr, ptrarray[pos+1]) > 0 );
         ptrarray[pos] = tmpptr;
         intarray1[pos] = tmpint1;
         intarray2[pos] = tmpint2;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], ptrarray[pos]) <= 0 )
            pos--;
         if( pos <= firstpos )
            break;
         assert( ptrcomp(ptrarray[pos-1], ptrarray[pos]) > 0 );
         tmpptr = ptrarray[pos];
         tmpint1 = intarray1[pos];
         tmpint2 = intarray2[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos-1];
            intarray1[pos] = intarray1[pos-1];
            intarray2[pos] = intarray2[pos-1];
            pos--;
         }
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], tmpptr) > 0 );
         ptrarray[pos] = tmpptr;
         intarray1[pos] = tmpint1;
         intarray2[pos] = tmpint2;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of three joint arrays of pointers/Reals/ints, sorted by first */
void SCIPbsortPtrRealInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            realarray,          /**< Real array to be permuted in the same way */
   int*             intarray,           /**< int array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   void* tmpptr;
   Real tmpreal;
   int tmpint;

   assert(len == 0 || ptrarray != NULL);
   assert(len == 0 || realarray != NULL);
   assert(len == 0 || intarray != NULL);
   assert(ptrcomp != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && ptrcomp(ptrarray[pos], ptrarray[pos+1]) <= 0 )
            pos++;
         if( pos >= lastpos )
            break;
         assert( ptrcomp(ptrarray[pos], ptrarray[pos+1]) > 0 );
         tmpptr = ptrarray[pos];
         tmpreal = realarray[pos];
         tmpint = intarray[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos+1];
            realarray[pos] = realarray[pos+1];
            intarray[pos] = intarray[pos+1];
            pos++;
         }
         while( pos < lastpos && ptrcomp(tmpptr, ptrarray[pos+1]) > 0 );
         ptrarray[pos] = tmpptr;
         realarray[pos] = tmpreal;
         intarray[pos] = tmpint;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], ptrarray[pos]) <= 0 )
            pos--;
         if( pos <= firstpos )
            break;
         assert( ptrcomp(ptrarray[pos-1], ptrarray[pos]) > 0 );
         tmpptr = ptrarray[pos];
         tmpreal = realarray[pos];
         tmpint = intarray[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos-1];
            realarray[pos] = realarray[pos-1];
            intarray[pos] = intarray[pos-1];
            pos--;
         }
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], tmpptr) > 0 );
         ptrarray[pos] = tmpptr;
         realarray[pos] = tmpreal;
         intarray[pos] = tmpint;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of four joint arrays of pointers/Reals/ints/ints, sorted by first */
void SCIPbsortPtrRealIntInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            realarray,          /**< Real array to be permuted in the same way */
   int*             intarray1,          /**< first int array to be permuted in the same way */
   int*             intarray2,          /**< second int array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   void* tmpptr;
   Real tmpreal;
   int tmpint1;
   int tmpint2;

   assert(len == 0 || ptrarray != NULL);
   assert(len == 0 || realarray != NULL);
   assert(len == 0 || intarray1 != NULL);
   assert(len == 0 || intarray2 != NULL);
   assert(ptrcomp != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && ptrcomp(ptrarray[pos], ptrarray[pos+1]) <= 0 )
            pos++;
         if( pos >= lastpos )
            break;
         assert( ptrcomp(ptrarray[pos], ptrarray[pos+1]) > 0 );
         tmpptr = ptrarray[pos];
         tmpreal = realarray[pos];
         tmpint1 = intarray1[pos];
         tmpint2 = intarray2[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos+1];
            realarray[pos] = realarray[pos+1];
            intarray1[pos] = intarray1[pos+1];
            intarray2[pos] = intarray2[pos+1];
            pos++;
         }
         while( pos < lastpos && ptrcomp(tmpptr, ptrarray[pos+1]) > 0 );
         ptrarray[pos] = tmpptr;
         realarray[pos] = tmpreal;
         intarray1[pos] = tmpint1;
         intarray2[pos] = tmpint2;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], ptrarray[pos]) <= 0 )
            pos--;
         if( pos <= firstpos )
            break;
         assert( ptrcomp(ptrarray[pos-1], ptrarray[pos]) > 0 );
         tmpptr = ptrarray[pos];
         tmpreal = realarray[pos];
         tmpint1 = intarray1[pos];
         tmpint2 = intarray2[pos];
         do
         {
            ptrarray[pos] = ptrarray[pos-1];
            realarray[pos] = realarray[pos-1];
            intarray1[pos] = intarray1[pos-1];
            intarray2[pos] = intarray2[pos-1];
            pos--;
         }
         while( pos > firstpos && ptrcomp(ptrarray[pos-1], tmpptr) > 0 );
         ptrarray[pos] = tmpptr;
         realarray[pos] = tmpreal;
         intarray1[pos] = tmpint1;
         intarray2[pos] = tmpint2;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of two joint arrays of Reals/pointers, sorted s.t. the Real array is in non-decreasing order */
void SCIPbsortRealPtr(
   Real*            realarray,          /**< Real array to be permuted in the same way */
   void**           ptrarray,           /**< pointer array to be sorted */
   int              len                 /**< length of arrays */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   Real tmpreal;
   void* tmpptr;

   assert(len == 0 || realarray != NULL);
   assert(len == 0 || ptrarray != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && realarray[pos] <= realarray[pos+1] )
            pos++;
         if( pos >= lastpos )
            break;
         tmpreal = realarray[pos];
         tmpptr = ptrarray[pos];
         do
         {
            realarray[pos] = realarray[pos+1];
            ptrarray[pos] = ptrarray[pos+1];
            pos++;
         }
         while( pos < lastpos && tmpreal > realarray[pos+1] );
         realarray[pos] = tmpreal;
         ptrarray[pos] = tmpptr;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && realarray[pos-1] <= realarray[pos] )
            pos--;
         if( pos <= firstpos )
            break;
         tmpreal = realarray[pos];
         tmpptr = ptrarray[pos];
         do
         {
            realarray[pos] = realarray[pos-1];
            ptrarray[pos] = ptrarray[pos-1];
            pos--;
         }
         while( pos > firstpos && realarray[pos-1] > tmpreal );
         realarray[pos] = tmpreal;
         ptrarray[pos] = tmpptr;
         sortpos = pos;
         pos--;
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

   assert(val1 > 0);
   assert(val2 > 0);

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
   Longint          val1,               /**< first value of smallest common multiple calculation */
   Longint          val2                /**< second value of smallest common multiple calculation */
   )
{
   Longint gcd;

   assert(val1 > 0);
   assert(val2 > 0);

   gcd = SCIPcalcGreComDiv(val1, val2);
   
   return val1/gcd * val2;
}

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
Bool SCIPrealToRational(
   Real             val,                /**< real value r to convert into rational number */
   Real             mindelta,           /**< minimal allowed difference r - q of real r and rational q = n/d */
   Real             maxdelta,           /**< maximal allowed difference r - q of real r and rational q = n/d */
   Longint          maxdnom,            /**< maximal denominator allowed */
   Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   )
{
   Real a;
   Real b;
   Real g0;
   Real g1;
   Real gx;
   Real h0;
   Real h1;
   Real hx;
   Real delta0;
   Real delta1;
   Real epsilon;

   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(nominator != NULL);
   assert(denominator != NULL);

   epsilon = MIN(-mindelta, maxdelta);

   b = val;
   a = EPSFLOOR(b, epsilon);
   g0 = a;
   h0 = 1.0;
   g1 = 1.0;
   h1 = 0.0;
   delta0 = val - g0/h0;
   delta1 = (delta0 < 0.0 ? val - (g0-1.0)/h0 : val - (g0+1.0)/h0);
  
   while( (delta0 < mindelta || delta0 > maxdelta) && (delta1 < mindelta || delta1 > maxdelta) )
   {
      assert(EPSGT(b, a, epsilon));
      assert(h0 >= 0.0);
      assert(h1 >= 0.0);

      b = 1.0 / (b - a);
      a = EPSFLOOR(b, epsilon);

      assert(a >= 0.0);
      gx = g0;
      hx = h0;

      g0 = a * g0 + g1;
      h0 = a * h0 + h1;

      g1 = gx;
      h1 = hx;
      
      if( h0 > maxdnom )
         return FALSE;
      
      delta0 = val - g0/h0;
      delta1 = (delta0 < 0.0 ? val - (g0-1.0)/h0 : val - (g0+1.0)/h0);
   }

   if( REALABS(g0) > (LONGINT_MAX >> 4) || h0 > (LONGINT_MAX >> 4) )
      return FALSE;

   assert(h0 > 0.5);

   if( delta0 < mindelta )
   {
      assert(mindelta <= delta1 && delta1 <= maxdelta);
      *nominator = (Longint)(g0 - 1.0);
      *denominator = (Longint)h0;
   }
   else if( delta0 > maxdelta )
   {
      assert(mindelta <= delta1 && delta1 <= maxdelta);
      *nominator = (Longint)(g0 + 1.0);
      *denominator = (Longint)h0;
   }
   else
   {
      *nominator = (Longint)g0;
      *denominator = (Longint)h0;
   }
   assert(*denominator >= 1);
   assert(val - (Real)(*nominator)/(Real)(*denominator) >= mindelta);
   assert(val - (Real)(*nominator)/(Real)(*denominator) <= maxdelta);

   return TRUE;
}

/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */
static
Bool isIntegralScalar(
   Real             val,                /**< value that should be scaled to an integral value */
   Real             scalar,             /**< scalar that should be tried */
   Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   Real sval;
   Real abssval;
   Real downval;
   Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = EPSFLOOR(sval, 0.0);
   upval = EPSCEIL(sval, 0.0);

   return (SCIPrelDiff(sval, downval) <= maxdelta || SCIPrelDiff(sval, upval) >= mindelta);
}

/** additional scalars that are tried in integrality scaling */
static const Real scalars[] = {3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};
static const int nscalars = 9;

/** tries to find a value, such that all given values, if scaled with this value become integral */
RETCODE SCIPcalcIntegralScalar(
   Real*            vals,               /**< values to scale */
   int              nvals,              /**< number of values to scale */
   Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal allowed scalar */
   Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   Bool*            success             /**< stores whether returned value is valid */
   )
{
   Longint gcd;
   Longint scm;
   Longint nominator;
   Longint denominator;
   Real val;
   Real minval;
   Real usedtol;
   Real absval;
   Real scaleval;
   Real twomultval;
   Bool scalable;
   Bool twomult;
   Bool rational;
   int c;
   int s;

   assert(vals != NULL);
   assert(nvals > 0);
   assert(maxdnom >= 1);
   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(success != NULL);

   debugMessage("trying to find rational representation for given values\n");

   if( intscalar != NULL )
      *intscalar = SCIP_INVALID;
   *success = FALSE;

   /* get minimal absolute non-zero value */
   minval = REAL_MAX;
   for( c = 0; c < nvals; ++c )
   {
      val = vals[c];
      if( val < mindelta || val > maxdelta )
      {
         absval = REALABS(val);
         minval = MIN(minval, absval);
      }
   }

   if( minval == REAL_MAX )
   {
      /* all coefficients are zero (inside tolerances) */
      if( intscalar != NULL )
         *intscalar = 1.0;
      *success = TRUE;
      debugMessage(" -> all values are zero (inside tolerances)\n");

      return SCIP_OKAY;
   }
   assert(minval > MIN(-mindelta, maxdelta));

   /* try, if values can be made integral multiplying them with the reciprocal of the smallest value and a power of 2 */
   scalable = TRUE;
   scaleval = 1.0/minval;
   for( c = 0; c < nvals && scalable; ++c )
   {
      /* check, if the value can be scaled with a simple scalar */
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;
    
      absval = REALABS(val);
      while( scaleval <= maxscale
         && (absval * scaleval < 0.5 || !isIntegralScalar(val, scaleval, mindelta, maxdelta)) )
      {
         for( s = 0; s < nscalars; ++s )
         {
            if( isIntegralScalar(val, scaleval * scalars[s], mindelta, maxdelta) )
            {
               scaleval *= scalars[s];
               break;
            }
         }
         if( s >= nscalars )
            scaleval *= 2.0;
      }
      scalable = (scaleval <= maxscale);
      debugMessage(" -> val=%g, scaleval=%g, val*scaleval=%g, scalable=%d\n", 
         val, scaleval, val*scaleval, scalable);
   }
   if( scalable )
   {
      /* make values integral by dividing them by the smallest value (and multiplying them with a power of 2) */
      assert(scaleval <= maxscale);
      if( intscalar != NULL )
         *intscalar = scaleval;
      *success = TRUE;
      debugMessage(" -> integrality can be achieved by scaling with %g (minval=%g)\n", scaleval, minval);
      
      return SCIP_OKAY;
   }

   /* try, if values can be made integral by multiplying them by a power of 2 */
   twomult = TRUE;
   twomultval = 1.0;
   for( c = 0; c < nvals && twomult; ++c )
   {
      /* check, if the value can be scaled with a simple scalar */
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;
      
      absval = REALABS(val);
      while( twomultval <= maxscale
         && (absval * twomultval < 0.5 || !isIntegralScalar(val, twomultval, mindelta, maxdelta)) )
      {
         for( s = 0; s < nscalars; ++s )
         {
            if( isIntegralScalar(val, twomultval * scalars[s], mindelta, maxdelta) )
            {
               twomultval *= scalars[s];
               break;
            }
         }
         if( s >= nscalars )
            twomultval *= 2.0;
      }
      twomult = (twomultval <= maxscale);
      debugMessage(" -> val=%g, twomult=%g, val*twomult=%g, twomultable=%d\n",
         val, twomultval, val*twomultval, twomult);
   }
   if( twomult )
   {
      /* make values integral by multiplying them with a power of 2 */
      assert(twomultval <= maxscale);
      if( intscalar != NULL )
         *intscalar = twomultval;
      *success = TRUE;
      debugMessage(" -> integrality can be achieved by scaling with %g (power of 2)\n", twomultval);
      
      return SCIP_OKAY;
   }

   /* convert each value into a rational number, calculate the greatest common divisor of the nominators
    * and the smallest common multiple of the denominators
    */
   gcd = 1;
   scm = 1;
   rational = TRUE;

   /* first value (to initialize gcd) */
   for( c = 0; c < nvals && rational; ++c )
   {
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;

      rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
      if( rational && nominator != 0 )
      {
         assert(denominator > 0);
         gcd = ABS(nominator);
         scm = denominator;
         rational = ((Real)scm/(Real)gcd <= maxscale);
         debugMessage(" -> c=%d first rational: val: %g == %lld/%lld, gcd=%lld, scm=%lld, rational=%d\n",
            c, val, nominator, denominator, gcd, scm, rational);
         break;
      }
   }

   /* remaining values */
   for( ++c; c < nvals && rational; ++c )
   {
      val = vals[c];
      if( val == 0.0 ) /* zeros are allowed in the vals array */
         continue;

      rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
      if( rational && nominator != 0 )
      {
         assert(denominator > 0);
         gcd = SCIPcalcGreComDiv(gcd, ABS(nominator));
         scm *= denominator / SCIPcalcGreComDiv(scm, denominator);
         rational = ((Real)scm/(Real)gcd <= maxscale);
         debugMessage(" -> c=%d next rational : val: %g == %lld/%lld, gcd=%lld, scm=%lld, rational=%d\n",
            c, val, nominator, denominator, gcd, scm, rational);
      }
      else
      {
         debugMessage(" -> failed to convert %g into a rational representation\n", val);
      }
   }

   if( rational )
   {
      /* make values integral by multiplying them with the smallest common multiple of the denominators */
      assert((Real)scm/(Real)gcd <= maxscale);
      if( intscalar != NULL )
         *intscalar = (Real)scm/(Real)gcd;
      *success = TRUE;
      debugMessage(" -> integrality can be achieved by scaling with %g (rational:%lld/%lld)\n", 
         (Real)scm/(Real)gcd, scm, gcd);
   }

   return SCIP_OKAY;
}


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
Real SCIPrelDiff(
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real absval1;
   Real absval2;
   Real quot;

   absval1 = REALABS(val1);
   absval2 = REALABS(val2);
   quot = MAX3(1.0, absval1, absval2);
   
   return (val1-val2)/quot;
}

#endif




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

/** splits filename into path, name, and extension */
void SCIPsplitFilename(
   char*            filename,           /**< filename to split; is destroyed (but not freed) during process */
   char**           path,               /**< pointer to store path, or NULL if not needed */
   char**           name,               /**< pointer to store name, or NULL if not needed */
   char**           extension           /**< pointer to store extension, or NULL if not needed */
   )
{
   char* lastslash;
   char* lastdot;

   assert(filename != NULL);

   lastslash = strrchr(filename, '/');
   lastdot = strrchr(filename, '.');
   if( lastslash != NULL && lastdot != NULL && lastdot < lastslash ) /* is the last dot belonging to the path? */
      lastdot = NULL;

   if( lastslash == NULL )
   {
      if( path != NULL )
         *path = NULL;
      if( name != NULL )
         *name = filename;
   }
   else
   {
      if( path != NULL )
         *path = filename;
      if( name != NULL )
         *name = lastslash+1;
      *lastslash = '\0';
   }

   if( lastdot == NULL )
   {
      if( extension != NULL )
         *extension = NULL;
   }
   else
   {
      if( extension != NULL )
         *extension = lastdot+1;
      *lastdot = '\0';
   }
}

