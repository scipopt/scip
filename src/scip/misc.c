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

/**@file   misc.c
 * @brief  miscellaneous methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "scip/def.h"
#include "scip/pub_message.h"
#include "scip/misc.h"
#include "scip/intervalarith.h"

#ifndef NDEBUG
#include "scip/struct_misc.h"
#endif

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

/*
 * GML graphical printing methods
 * For a detailed format decription see http://docs.yworks.com/yfiles/doc/developers-guide/gml.html
 */

#define GMLNODEWIDTH 120.0
#define GMLNODEHEIGTH 30.0
#define GMLFONTSIZE 13
#define GMLNODETYPE "rectangle"
#define GMLNODEFILLCOLOR "#ff0000"
#define GMLEDGECOLOR "black"
#define GMLNODEBORDERCOLOR "#000000"


/** writes a node section to the given graph file */
void SCIPgmlWriteNode(
   FILE*                 file,               /**< file to write to */
   unsigned int          id,                 /**< id of the node */
   const char*           label,              /**< label of the node */
   const char*           nodetype,           /**< type of the node, or NULL */
   const char*           fillcolor,          /**< color of the node's interior, or NULL */
   const char*           bordercolor         /**< color of the node's border, or NULL */
   )
{
   assert(file != NULL);
   assert(label != NULL);

   fprintf(file, "  node\n");
   fprintf(file, "  [\n");
   fprintf(file, "    id      %u\n", id);
   fprintf(file, "    label   \"%s\"\n", label);
   fprintf(file, "    graphics\n");
   fprintf(file, "    [\n");
   fprintf(file, "      w       %g\n", GMLNODEWIDTH);
   fprintf(file, "      h       %g\n", GMLNODEHEIGTH);

   if( nodetype != NULL )
      fprintf(file, "      type    \"%s\"\n", nodetype);
   else
      fprintf(file, "      type    \"%s\"\n", GMLNODETYPE);

   if( fillcolor != NULL )
      fprintf(file, "      fill    \"%s\"\n", fillcolor);
   else
      fprintf(file, "      fill    \"%s\"\n", GMLNODEFILLCOLOR);

   if( bordercolor != NULL )
      fprintf(file, "      outline \"%s\"\n", bordercolor);
   else
      fprintf(file, "      outline \"%s\"\n", GMLNODEBORDERCOLOR);

   fprintf(file, "    ]\n");
   fprintf(file, "    LabelGraphics\n");
   fprintf(file, "    [\n");
   fprintf(file, "      text      \"%s\"\n", label);
   fprintf(file, "      fontSize  %d\n", GMLFONTSIZE);
   fprintf(file, "      fontName  \"Dialog\"\n");
   fprintf(file, "      anchor    \"c\"\n");
   fprintf(file, "    ]\n");
   fprintf(file, "  ]\n");
}

/** writes an edge section to the given graph file */
void SCIPgmlWriteEdge(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   )
{
   assert(file != NULL);

   fprintf(file, "  edge\n");
   fprintf(file, "  [\n");
   fprintf(file, "    source  %u\n", source);
   fprintf(file, "    target  %u\n", target);

   if( label != NULL)
      fprintf(file, "    label   \"%s\"\n", label);

   fprintf(file, "    graphics\n");
   fprintf(file, "    [\n");

   if( color != NULL )
      fprintf(file, "      fill    \"%s\"\n", color);
   else
      fprintf(file, "      fill    \"%s\"\n", GMLEDGECOLOR);

   /* fprintf(file, "      arrow     \"both\"\n"); */
   fprintf(file, "    ]\n");

   if( label != NULL)
   {
      fprintf(file, "    LabelGraphics\n");
      fprintf(file, "    [\n");
      fprintf(file, "      text      \"%s\"\n", label);
      fprintf(file, "      fontSize  %d\n", GMLFONTSIZE);
      fprintf(file, "      fontName  \"Dialog\"\n");
      fprintf(file, "      anchor    \"c\"\n");
      fprintf(file, "    ]\n");
   }

   fprintf(file, "  ]\n");
}

/** writes an arc section to the given graph file */
void SCIPgmlWriteArc(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   )
{
   assert(file != NULL);

   fprintf(file, "  edge\n");
   fprintf(file, "  [\n");
   fprintf(file, "    source  %u\n", source);
   fprintf(file, "    target  %u\n", target);

   if( label != NULL)
      fprintf(file, "    label   \"%s\"\n", label);

   fprintf(file, "    graphics\n");
   fprintf(file, "    [\n");

   if( color != NULL )
      fprintf(file, "      fill    \"%s\"\n", color);
   else
      fprintf(file, "      fill    \"%s\"\n", GMLEDGECOLOR);

   fprintf(file, "      targetArrow     \"standard\"\n");
   fprintf(file, "    ]\n");

   if( label != NULL)
   {
      fprintf(file, "    LabelGraphics\n");
      fprintf(file, "    [\n");
      fprintf(file, "      text      \"%s\"\n", label);
      fprintf(file, "      fontSize  %d\n", GMLFONTSIZE);
      fprintf(file, "      fontName  \"Dialog\"\n");
      fprintf(file, "      anchor    \"c\"\n");
      fprintf(file, "    ]\n");
   }

   fprintf(file, "  ]\n");
}

/** writes the starting line to a GML graph file, does not open a file */
void SCIPgmlWriteOpening(
   FILE*                 file,               /**< file to write to */
   SCIP_Bool             directed            /**< is the graph directed */
   )
{
   assert(file != NULL);

   fprintf(file, "graph\n");
   fprintf(file, "[\n");
   fprintf(file, "  hierarchic      1\n");

   if( directed )
      fprintf(file, "  directed        1\n");
}

/** writes the ending lines to a GML graph file, does not close a file */
void SCIPgmlWriteClosing(
   FILE*                 file                /**< file to close */
   )
{
   assert(file != NULL);

   fprintf(file, "]\n");
}


/*
 * Sparse solution
 */

/** creates a sparse solution */
SCIP_RETCODE SCIPsparseSolCreate(
   SCIP_SPARSESOL**      sparsesol,          /**< pointer to store the created sparse solution */
   SCIP_VAR**            vars,               /**< variables in the sparse solution, must not contain continuous
					      *   variables
					      */
   int                   nvars,              /**< number of variables to store, size of the lower and upper bound
					      *   arrays
					      */
   SCIP_Bool             cleared             /**< should the lower and upper bound arrays be cleared (entries set to
					      *	  0)
					      */
   )
{
   assert(sparsesol != NULL);
   assert(vars != NULL);
   assert(nvars > 0);

   SCIP_ALLOC( BMSallocMemory(sparsesol) );

#ifndef NDEBUG
   {
      int v;

      for( v = nvars - 1; v >= 0; --v )
      {
	 assert(vars[v] != NULL);
	 /* assert(SCIPvarGetType(vars[v]) != SCIP_VARTYPE_CONTINUOUS); */
      }
   }
#endif

   /* copy variables */
   SCIP_ALLOC( BMSduplicateMemoryArray(&((*sparsesol)->vars), vars, nvars) );

   /* create bound arrays */
   if( cleared )
   {
      SCIP_ALLOC( BMSallocClearMemoryArray(&((*sparsesol)->lbvalues), nvars) );
      SCIP_ALLOC( BMSallocClearMemoryArray(&((*sparsesol)->ubvalues), nvars) );
   }
   else
   {
      SCIP_ALLOC( BMSallocMemoryArray(&((*sparsesol)->lbvalues), nvars) );
      SCIP_ALLOC( BMSallocMemoryArray(&((*sparsesol)->ubvalues), nvars) );
   }

   (*sparsesol)->nvars = nvars;

   return SCIP_OKAY;
}

/** frees priority queue, but not the data elements themselves */
void SCIPsparseSolFree(
   SCIP_SPARSESOL**      sparsesol           /**< pointer to a sparse solution */
   )
{
   assert(sparsesol != NULL);
   assert(*sparsesol != NULL);

   BMSfreeMemoryArray(&((*sparsesol)->vars));
   BMSfreeMemoryArray(&((*sparsesol)->ubvalues));
   BMSfreeMemoryArray(&((*sparsesol)->lbvalues));
   BMSfreeMemory(sparsesol);
}

/** returns the variables stored in the given sparse solution */
SCIP_VAR** SCIPsparseSolGetVars(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   )
{
   assert(sparsesol != NULL);

   return sparsesol->vars;
}

/** returns the number of variables stored in the given sparse solution */
int SCIPsparseSolGetNVars(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   )
{
   assert(sparsesol != NULL);

   return sparsesol->nvars;
}

/** returns the lower bound array for all variables for a given sparse solution */
SCIP_Longint* SCIPsparseSolGetLbs(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   )
{
   assert(sparsesol != NULL);

   return sparsesol->lbvalues;
}

/** returns the upper bound array for all variables for a given sparse solution */
SCIP_Longint* SCIPsparseSolGetUbs(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   )
{
   assert(sparsesol != NULL);

   return sparsesol->ubvalues;
}


/*
 * Priority Queue
 */

#define PQ_PARENT(q) (((q)+1)/2-1)
#define PQ_LEFTCHILD(p) (2*(p)+1)
#define PQ_RIGHTCHILD(p) (2*(p)+2)


/** resizes element memory to hold at least the given number of elements */
static
SCIP_RETCODE pqueueResize(
   SCIP_PQUEUE*          pqueue,             /**< pointer to a priority queue */
   int                   minsize             /**< minimal number of storable elements */
   )
{
   assert(pqueue != NULL);
   
   if( minsize <= pqueue->size )
      return SCIP_OKAY;

   pqueue->size = MAX(minsize, (int)(pqueue->size * pqueue->sizefac));
   SCIP_ALLOC( BMSreallocMemoryArray(&pqueue->slots, pqueue->size) );

   return SCIP_OKAY;
}

/** creates priority queue */
SCIP_RETCODE SCIPpqueueCreate(
   SCIP_PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   )
{
   assert(pqueue != NULL);
   assert(ptrcomp != NULL);

   initsize = MAX(1, initsize);
   sizefac = MAX(1.0, sizefac);

   SCIP_ALLOC( BMSallocMemory(pqueue) );
   (*pqueue)->len = 0;
   (*pqueue)->size = 0;
   (*pqueue)->sizefac = sizefac;
   (*pqueue)->slots = NULL;
   (*pqueue)->ptrcomp = ptrcomp;
   SCIP_CALL( pqueueResize(*pqueue, initsize) );

   return SCIP_OKAY;
}

/** frees priority queue, but not the data elements themselves */
void SCIPpqueueFree(
   SCIP_PQUEUE**         pqueue              /**< pointer to a priority queue */
   )
{
   assert(pqueue != NULL);

   BMSfreeMemoryArray(&(*pqueue)->slots);
   BMSfreeMemory(pqueue);
}

/** clears the priority queue, but doesn't free the data elements themselves */
void SCIPpqueueClear(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);

   pqueue->len = 0;
}

/** inserts element into priority queue */
SCIP_RETCODE SCIPpqueueInsert(
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   void*                 elem                /**< element to be inserted */
   )
{
   int pos;

   assert(pqueue != NULL);
   assert(pqueue->len >= 0);
   assert(elem != NULL);

   SCIP_CALL( pqueueResize(pqueue, pqueue->len+1) );

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
   SCIP_PQUEUE*          pqueue              /**< priority queue */
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
   SCIP_PQUEUE*          pqueue              /**< priority queue */
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
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   return pqueue->len;
}

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
void** SCIPpqueueElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   )
{
   assert(pqueue != NULL);
   assert(pqueue->len >= 0);

   return pqueue->slots;
}




/*
 * Hash Table
 */

/** table of some prime numbers */
static int primetable[] = {
   2,
   7,
   19,
   31,
   59,
   227,
   617,
   1523,
   3547,
   8011,
   17707,
   38723,
   83833,
   180317,
   385897,
   821411,
   1742369,
   3680893,
   5693959,
   7753849,
   9849703,
   11973277,
   14121853,
   17643961,
   24273817,
   32452843,
   49979687,
   67867967,
   86028121,
   104395301,
   122949823,
   141650939,
   160481183,
   179424673,
   198491317,
   217645177,
   256203161,
   314606869,
   373587883,
   433024223,
   492876847,
   553105243,
   613651349,
   694847533,
   756065159,
   817504243,
   879190747,
   941083981,
   982451653,
   INT_MAX
};
static const int primetablesize = sizeof(primetable)/sizeof(int);

/** returns a reasonable hash table size (a prime number) that is at least as large as the specified value */
int SCIPcalcHashtableSize(
   int                   minsize             /**< minimal size of the hash table */
   )
{
   int pos;

   (void) SCIPsortedvecFindInt(primetable, minsize, primetablesize, &pos);
   assert(pos < primetablesize);

   return primetable[pos];
}

/** appends element to the hash list */
static
SCIP_RETCODE hashtablelistAppend(
   SCIP_HASHTABLELIST**  hashtablelist,      /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 element             /**< element to append to the list */
   )
{
   SCIP_HASHTABLELIST* newlist;

   assert(hashtablelist != NULL);
   assert(blkmem != NULL);
   assert(element != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &newlist) );
   newlist->element = element;
   newlist->next = *hashtablelist;
   *hashtablelist = newlist;

   return SCIP_OKAY;
}

/** frees a hash list entry and all its successors */
static
void hashtablelistFree(
   SCIP_HASHTABLELIST**  hashtablelist,      /**< pointer to hash list to free */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_HASHTABLELIST* list;
   SCIP_HASHTABLELIST* nextlist;

   assert(hashtablelist != NULL);
   assert(blkmem != NULL);
   
   list = *hashtablelist;
   while( list != NULL )
   {
      nextlist = list->next;
      BMSfreeBlockMemory(blkmem, &list);
      list = nextlist;
   }

   *hashtablelist = NULL;
}

/** finds hash list entry pointing to element with given key in the hash list, returns NULL if not found */
static
SCIP_HASHTABLELIST* hashtablelistFind(
   SCIP_HASHTABLELIST*   hashtablelist,      /**< hash list */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr,            /**< user pointer */
   unsigned int          keyval,             /**< hash value of key */
   void*                 key                 /**< key to retrieve */
   )
{
   unsigned int currentkeyval;
   void* currentkey;

   assert(hashkeyeq != NULL);
   assert(key != NULL);

   while( hashtablelist != NULL )
   {
      currentkey = hashgetkey(userptr, hashtablelist->element);
      currentkeyval = hashkeyval(userptr, currentkey);
      if( currentkeyval == keyval && hashkeyeq(userptr, currentkey, key) )
         return hashtablelist;
      
      hashtablelist = hashtablelist->next;
   }

   return NULL;
}

/** retrieves element with given key from the hash list, or NULL */
static
void* hashtablelistRetrieve(
   SCIP_HASHTABLELIST*   hashtablelist,      /**< hash list */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr,            /**< user pointer */
   unsigned int          keyval,             /**< hash value of key */
   void*                 key                 /**< key to retrieve */
   )
{
   SCIP_HASHTABLELIST* h;

   /* find hash list entry */
   h = hashtablelistFind(hashtablelist, hashgetkey, hashkeyeq, hashkeyval, userptr, keyval, key);

   /* return element */
   if( h != NULL )
   {
#ifndef NDEBUG
      SCIP_HASHTABLELIST* h2;

      h2 = hashtablelistFind(h->next, hashgetkey, hashkeyeq, hashkeyval, userptr, keyval, key);

      if( h2 != NULL )
      {
         void* key1;
         void* key2;

         key1 = hashgetkey(userptr, h->element);
         key2 = hashgetkey(userptr, h2->element);
         assert(hashkeyval(userptr, key1) == hashkeyval(userptr, key2));

         if( hashkeyeq(userptr, key1, key2) )
         {
            SCIPerrorMessage("WARNING: hashkey with same value exists multiple times (e.g. duplicate constraint/variable names), so the return value is maybe not correct\n");
         }
      }
#endif

      return h->element;
   }
   else
      return NULL;
}


/** retrieves element with given key from the hash list, or NULL
 * returns pointer to hash table list entry */
static
void* hashtablelistRetrieveNext(
   SCIP_HASHTABLELIST**  hashtablelist,      /**< on input: hash list to search; on exit: hash list entry corresponding to element after retrieved one, or NULL */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr,            /**< user pointer */
   unsigned int          keyval,             /**< hash value of key */
   void*                 key                 /**< key to retrieve */
   )
{
   SCIP_HASHTABLELIST* h;

   assert(hashtablelist != NULL);
   
   /* find hash list entry */
   h = hashtablelistFind(*hashtablelist, hashgetkey, hashkeyeq, hashkeyval, userptr, keyval, key);

   /* return element */
   if( h != NULL )
   {
      *hashtablelist = h->next;
      
      return h->element;
   }
   
   *hashtablelist = NULL;
   
   return NULL;
}

/** removes element from the hash list */
static
SCIP_RETCODE hashtablelistRemove(
   SCIP_HASHTABLELIST**  hashtablelist,      /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 element             /**< element to remove from the list */
   )
{
   SCIP_HASHTABLELIST* nextlist;

   assert(hashtablelist != NULL);
   assert(blkmem != NULL);
   assert(element != NULL);

   while( *hashtablelist != NULL && (*hashtablelist)->element != element )
      hashtablelist = &(*hashtablelist)->next;

   if( *hashtablelist != NULL )
   {
      nextlist = (*hashtablelist)->next;
      BMSfreeBlockMemory(blkmem, hashtablelist);
      *hashtablelist = nextlist;
   }

   return SCIP_OKAY;
}

/** creates a hash table */
SCIP_RETCODE SCIPhashtableCreate(
   SCIP_HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash table entries */
   int                   tablesize,          /**< size of the hash table */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr             /**< user pointer */
   )
{
   assert(hashtable != NULL);
   assert(tablesize > 0);
   assert(hashgetkey != NULL);
   assert(hashkeyeq != NULL);
   assert(hashkeyval != NULL);

   SCIP_ALLOC( BMSallocMemory(hashtable) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*hashtable)->lists, tablesize) );
   (*hashtable)->blkmem = blkmem;
   (*hashtable)->nlists = tablesize;
   (*hashtable)->hashgetkey = hashgetkey;
   (*hashtable)->hashkeyeq = hashkeyeq;
   (*hashtable)->hashkeyval = hashkeyval;
   (*hashtable)->userptr = userptr;

   BMSclearMemoryArray((*hashtable)->lists, tablesize);

   return SCIP_OKAY;
}

/** frees the hash table */
void SCIPhashtableFree(
   SCIP_HASHTABLE**      hashtable           /**< pointer to the hash table */
   )
{
   int i;
   SCIP_HASHTABLE* table;
   BMS_BLKMEM* blkmem;
   SCIP_HASHTABLELIST** lists;

   assert(hashtable != NULL);
   assert(*hashtable != NULL);

   table = (*hashtable);
   blkmem = table->blkmem;
   lists = table->lists;

   /* free hash lists */
   for( i = table->nlists - 1; i >= 0; --i )
      hashtablelistFree(&lists[i], blkmem);

   /* free main hast table data structure */
   BMSfreeMemoryArray(&table->lists);
   BMSfreeMemory(hashtable);
}

/** removes all elements of the hash table
 *
 *  @note From a performance point of view you should not fill and clear a hash table too often since the clearing can
 *        be expensive. Clearing is done by looping over all buckets and removing the hash table lists one-by-one.
 */
void SCIPhashtableClear(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   )
{
   int i;
   BMS_BLKMEM* blkmem;
   SCIP_HASHTABLELIST** lists;

   assert(hashtable != NULL);

   blkmem = hashtable->blkmem;
   lists = hashtable->lists;

   /* free hash lists */
   for( i = hashtable->nlists - 1; i >= 0; --i )
      hashtablelistFree(&lists[i], blkmem);
}

/** inserts element in hash table (multiple inserts of same element possible) */
SCIP_RETCODE SCIPhashtableInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
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
   key = hashtable->hashgetkey(hashtable->userptr, element);
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   /* append element to the list at the hash position */
   SCIP_CALL( hashtablelistAppend(&hashtable->lists[hashval], hashtable->blkmem, element) );
   
   return SCIP_OKAY;
}

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
SCIP_RETCODE SCIPhashtableSafeInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   )
{
   assert(hashtable != NULL);
   assert(hashtable->hashgetkey != NULL);

   /* check, if key is already existing */
   if( SCIPhashtableRetrieve(hashtable, hashtable->hashgetkey(hashtable->userptr, element)) != NULL )
      return SCIP_KEYALREADYEXISTING;

   /* insert element in hash table */
   SCIP_CALL( SCIPhashtableInsert(hashtable, element) );
   
   return SCIP_OKAY;
}

/** retrieve element with key from hash table, returns NULL if not existing */
void* SCIPhashtableRetrieve(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 key                 /**< key to retrieve */
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
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   return hashtablelistRetrieve(hashtable->lists[hashval], hashtable->hashgetkey, hashtable->hashkeyeq, 
      hashtable->hashkeyval, hashtable->userptr, keyval, key);
}

/** retrieve element with key from hash table, returns NULL if not existing
 * can be used to retrieve all entries with the same key (one-by-one) */
void* SCIPhashtableRetrieveNext(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_HASHTABLELIST**  hashtablelist,      /**< input: entry in hash table list from which to start searching, or NULL; output: entry in hash table list corresponding to element after retrieved one, or NULL */
   void*                 key                 /**< key to retrieve */
   )
{
   unsigned int keyval;

   assert(hashtable != NULL);
   assert(hashtable->lists != NULL);
   assert(hashtable->nlists > 0);
   assert(hashtable->hashgetkey != NULL);
   assert(hashtable->hashkeyeq != NULL);
   assert(hashtable->hashkeyval != NULL);
   assert(hashtablelist != NULL);
   assert(key != NULL);

   keyval = hashtable->hashkeyval(hashtable->userptr, key);

   if( *hashtablelist == NULL )
   {
      unsigned int hashval;
      
      /* get the hash value of the key */
      hashval = keyval % hashtable->nlists; /*lint !e573*/
      
      *hashtablelist = hashtable->lists[hashval];
   }

   return hashtablelistRetrieveNext(hashtablelist, hashtable->hashgetkey, hashtable->hashkeyeq, 
      hashtable->hashkeyval, hashtable->userptr, keyval, key);
}

/** returns whether the given element exists in the table */
SCIP_Bool SCIPhashtableExists(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to search in the table */
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
   key = hashtable->hashgetkey(hashtable->userptr, element);
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   return (hashtablelistFind(hashtable->lists[hashval], hashtable->hashgetkey, hashtable->hashkeyeq,
         hashtable->hashkeyval, hashtable->userptr, keyval, key) != NULL);
}

/** removes element from the hash table, if it exists */
SCIP_RETCODE SCIPhashtableRemove(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to remove from the table */
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
   key = hashtable->hashgetkey(hashtable->userptr, element);
   keyval = hashtable->hashkeyval(hashtable->userptr, key);
   hashval = keyval % hashtable->nlists; /*lint !e573*/

   /* remove element from the list at the hash position */
   SCIP_CALL( hashtablelistRemove(&hashtable->lists[hashval], hashtable->blkmem, element) );
   
   return SCIP_OKAY;
}

/** prints statistics about hash table usage */
void SCIPhashtablePrintStatistics(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_HASHTABLELIST* hashtablelist;
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

   SCIPmessagePrintInfo(messagehdlr, "%d hash entries, used %d/%d slots (%.1f%%)",
      sumslotsize, usedslots, hashtable->nlists, 100.0*(SCIP_Real)usedslots/(SCIP_Real)(hashtable->nlists));
   if( usedslots > 0 )
      SCIPmessagePrintInfo(messagehdlr, ", avg. %.1f entries/used slot, max. %d entries in slot",
         (SCIP_Real)sumslotsize/(SCIP_Real)usedslots, maxslotsize);
   SCIPmessagePrintInfo(messagehdlr, "\n");
}


/** returns TRUE iff both keys (i.e. strings) are equal */
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqString)
{  /*lint --e{715}*/
   const char* string1 = (const char*)key1;
   const char* string2 = (const char*)key2;

   return (strcmp(string1, string2) == 0);
}

/** returns the hash value of the key (i.e. string) */
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValString)
{  /*lint --e{715}*/
   const char* str;
   unsigned int sum;

   str = (const char*)key;
   sum = 0;

   while( *str != '\0' )
   {
      sum *= 31;
      sum += (unsigned int)(*str); /*lint !e571*/
      str++;
   }

   return sum;
}


/** gets the element as the key */
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyStandard)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys(pointer) are equal */
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqPtr)
{  /*lint --e{715}*/
   return (key1 == key2);
}

/** returns the hash value of the key */
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValPtr)
{  /*lint --e{715}*/
   /* the key is used as the keyvalue too */
   return (unsigned int)(size_t) key;
}



/*
 * Hash Map
 */

/** appends origin->image pair to the hash list */
static
SCIP_RETCODE hashmaplistAppend(
   SCIP_HASHMAPLIST**    hashmaplist,        /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL */
   void*                 origin,             /**< origin of the mapping origin -> image */
   void*                 image               /**< image of the mapping origin -> image */
   )
{
   SCIP_HASHMAPLIST* newlist;

   assert(hashmaplist != NULL);
   assert(origin != NULL);

   if( blkmem != NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &newlist) );
   }
   else
   {
      SCIP_ALLOC( BMSallocMemory(&newlist) );
   }

   newlist->origin = origin;
   newlist->image = image;
   newlist->next = *hashmaplist;
   *hashmaplist = newlist;

   return SCIP_OKAY;
}

/** frees a hash list entry and all its successors */
static
void hashmaplistFree(
   SCIP_HASHMAPLIST**    hashmaplist,        /**< pointer to hash list to free */
   BMS_BLKMEM*           blkmem              /**< block memory, or NULL */
   )
{
   SCIP_HASHMAPLIST* list;
   SCIP_HASHMAPLIST* nextlist;

   assert(hashmaplist != NULL);
   
   list = *hashmaplist;
   while( list != NULL )
   {
      nextlist = list->next;

      if( blkmem != NULL )
      {
         BMSfreeBlockMemory(blkmem, &list);
      }
      else
      {
         BMSfreeMemory(&list);
      }

      list = nextlist;
   }

   *hashmaplist = NULL;
}

/** finds hash list entry pointing to given origin in the hash list, returns NULL if not found */
static
SCIP_HASHMAPLIST* hashmaplistFind(
   SCIP_HASHMAPLIST*     hashmaplist,        /**< hash list */
   void*                 origin              /**< origin to find */
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
   SCIP_HASHMAPLIST*     hashmaplist,        /**< hash list */
   void*                 origin              /**< origin to retrieve image for */
   )
{
   SCIP_HASHMAPLIST* h;

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
SCIP_RETCODE hashmaplistSetImage(
   SCIP_HASHMAPLIST**    hashmaplist,        /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   )
{
   SCIP_HASHMAPLIST* h;

   /* find hash list entry */
   h = hashmaplistFind(*hashmaplist, origin);

   /* set image or add origin->image pair */
   if( h != NULL )
      h->image = image;
   else
   {
      SCIP_CALL( hashmaplistAppend(hashmaplist, blkmem, origin, image) );
   }

   return SCIP_OKAY;
}

/** removes origin->image pair from the hash list */
static
SCIP_RETCODE hashmaplistRemove(
   SCIP_HASHMAPLIST**    hashmaplist,        /**< pointer to hash list */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL */
   void*                 origin              /**< origin to remove from the list */
   )
{
   SCIP_HASHMAPLIST* nextlist;

   assert(hashmaplist != NULL);
   assert(origin != NULL);

   while( *hashmaplist != NULL && (*hashmaplist)->origin != origin )
   {
      hashmaplist = &(*hashmaplist)->next;
   }
   if( *hashmaplist != NULL )
   {
      nextlist = (*hashmaplist)->next;

      if( blkmem != NULL )
      {
         BMSfreeBlockMemory(blkmem, hashmaplist);
      }
      else
      {
         BMSfreeMemory(hashmaplist);
      }

      *hashmaplist = nextlist;
   }

   return SCIP_OKAY;
}


/** creates a hash map mapping pointers to pointers 
 *
 * @note if possible always use a blkmem pointer instead of NULL, otherwise it could slow down the map
 */
SCIP_RETCODE SCIPhashmapCreate(
   SCIP_HASHMAP**        hashmap,            /**< pointer to store the created hash map */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash map entries, or NULL */
   int                   mapsize             /**< size of the hash map */
   )
{
   int i;

   assert(hashmap != NULL);
   assert(mapsize > 0);

   SCIP_ALLOC( BMSallocMemory(hashmap) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*hashmap)->lists, mapsize) );
   (*hashmap)->blkmem = blkmem;
   (*hashmap)->nlists = mapsize;

   /* initialize hash lists */
   for( i = 0; i < mapsize; ++i )
      (*hashmap)->lists[i] = NULL;

   return SCIP_OKAY;
}

/** frees the hash map */
void SCIPhashmapFree(
   SCIP_HASHMAP**        hashmap             /**< pointer to the hash map */
   )
{
   int i;

   assert(hashmap != NULL);
   assert(*hashmap != NULL);

   /* free hash lists */
   for( i = 0; i < (*hashmap)->nlists; ++i )
      hashmaplistFree(&(*hashmap)->lists[i], (*hashmap)->blkmem);

   /* free main hash map data structure */
   BMSfreeMemoryArray(&(*hashmap)->lists);
   BMSfreeMemory(hashmap);
}

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
SCIP_RETCODE SCIPhashmapInsert(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   /* append origin->image pair to the list at the hash position */
   SCIP_CALL( hashmaplistAppend(&hashmap->lists[hashval], hashmap->blkmem, origin, image) );
   
   return SCIP_OKAY;
}

/** retrieves image of given origin from the hash map, or NULL if no image exists */
void* SCIPhashmapGetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   /* get image for origin from hash list */
   return hashmaplistGetImage(hashmap->lists[hashval], origin);
}

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
SCIP_RETCODE SCIPhashmapSetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   /* set image for origin in hash list */
   SCIP_CALL( hashmaplistSetImage(&hashmap->lists[hashval], hashmap->blkmem, origin, image) );
   
   return SCIP_OKAY;
}

/** checks whether an image to the given origin exists in the hash map */
SCIP_Bool SCIPhashmapExists(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to search for */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   return (hashmaplistFind(hashmap->lists[hashval], origin) != NULL);
}

/** removes origin->image pair from the hash map, if it exists */
SCIP_RETCODE SCIPhashmapRemove(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to remove from the list */
   )
{
   unsigned int hashval;

   assert(hashmap != NULL);
   assert(hashmap->lists != NULL);
   assert(hashmap->nlists > 0);
   assert(origin != NULL);

   /* get the hash value */
   hashval = (unsigned int)((size_t)origin % (unsigned int)hashmap->nlists);

   /* remove element from the list at the hash position */
   SCIP_CALL( hashmaplistRemove(&hashmap->lists[hashval], hashmap->blkmem, origin) );
   
   return SCIP_OKAY;
}

/** prints statistics about hash map usage */
void SCIPhashmapPrintStatistics(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_HASHMAPLIST* hashmaplist;
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

   SCIPmessagePrintInfo(messagehdlr, "%d hash entries, used %d/%d slots (%.1f%%)",
      sumslotsize, usedslots, hashmap->nlists, 100.0*(SCIP_Real)usedslots/(SCIP_Real)(hashmap->nlists));
   if( usedslots > 0 )
      SCIPmessagePrintInfo(messagehdlr, ", avg. %.1f entries/used slot, max. %d entries in slot", 
         (SCIP_Real)sumslotsize/(SCIP_Real)usedslots, maxslotsize);
   SCIPmessagePrintInfo(messagehdlr, "\n");
}

/** indicates whether a hash map has no entries */
SCIP_Bool SCIPhashmapIsEmpty(
   SCIP_HASHMAP*         hashmap             /**< hash map */
)
{
   int i;
   assert(hashmap != NULL);
   
   for( i = 0; i < hashmap->nlists; ++i )
      if( hashmap->lists[i] )
         return FALSE;
   
   return TRUE;
}

/** gives the number of entries in a hash map */ 
int SCIPhashmapGetNEntries(
   SCIP_HASHMAP*         hashmap             /**< hash map */
)
{
   int count = 0;
   int i;
   assert(hashmap != NULL);
   
   for( i = 0; i < hashmap->nlists; ++i )
      count += SCIPhashmapListGetNEntries(hashmap->lists[i]);

   return count;
}

/** gives the number of lists (buckets) in a hash map */ 
int SCIPhashmapGetNLists(
   SCIP_HASHMAP*         hashmap             /**< hash map */
)
{
   assert(hashmap != NULL);
   
   return hashmap->nlists;
}

/** gives a specific list (bucket) in a hash map */
SCIP_HASHMAPLIST* SCIPhashmapGetList(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   int                   listindex           /**< index of hash map list */
)
{
   assert(hashmap != NULL);
   assert(listindex >= 0);
   assert(listindex < hashmap->nlists);
   
   return hashmap->lists[listindex];
}

/** gives the number of entries in a list of a hash map */ 
int SCIPhashmapListGetNEntries(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list, can be NULL */
)
{
   int count = 0;
   
   for( ; hashmaplist; hashmaplist = hashmaplist->next )
      ++count;
   
   return count;
}

/** retrieves origin of given entry in a hash map */ 
void* SCIPhashmapListGetOrigin(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list */
)
{
   assert(hashmaplist != NULL);
   
   return hashmaplist->origin;
}

/** retrieves image of given entry in a hash map */ 
void* SCIPhashmapListGetImage(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list */
)
{
   assert(hashmaplist != NULL);
   
   return hashmaplist->image;
}

/** retrieves next entry from given entry in a hash map list, or NULL if at end of list. */ 
SCIP_HASHMAPLIST* SCIPhashmapListGetNext(
   SCIP_HASHMAPLIST*     hashmaplist         /**< hash map list */
)
{
   assert(hashmaplist != NULL);
   
   return hashmaplist->next;
}

/** removes all entries in a hash map. */ 
SCIP_RETCODE SCIPhashmapRemoveAll(
   SCIP_HASHMAP*         hashmap             /**< hash map */
)
{
   int listidx;

   assert(hashmap != NULL);
   
   /* free hash lists */
   for( listidx = 0; listidx < hashmap->nlists; ++listidx )
      hashmaplistFree(&hashmap->lists[listidx], hashmap->blkmem);

   return SCIP_OKAY;
}



/*
 * Dynamic Arrays
 */

/** creates a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayCreate(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(realarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, realarray) );
   (*realarray)->blkmem = blkmem;
   (*realarray)->vals = NULL;
   (*realarray)->valssize = 0;
   (*realarray)->firstidx = -1;
   (*realarray)->minusedidx = INT_MAX;
   (*realarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayCopy(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the copied real array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REALARRAY*       sourcerealarray     /**< dynamic real array to copy */
   )
{
   assert(realarray != NULL);
   assert(sourcerealarray != NULL);

   SCIP_CALL( SCIPrealarrayCreate(realarray, blkmem) );
   if( sourcerealarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*realarray)->vals, sourcerealarray->vals,
                     sourcerealarray->valssize) );
   }
   (*realarray)->valssize = sourcerealarray->valssize;
   (*realarray)->firstidx = sourcerealarray->firstidx;
   (*realarray)->minusedidx = sourcerealarray->minusedidx;
   (*realarray)->maxusedidx = sourcerealarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayFree(
   SCIP_REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   assert(realarray != NULL);
   assert(*realarray != NULL);

   BMSfreeBlockMemoryArrayNull((*realarray)->blkmem, &(*realarray)->vals, (*realarray)->valssize);
   BMSfreeBlockMemory((*realarray)->blkmem, realarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPrealarrayExtend(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
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

   SCIPdebugMessage("extending realarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n",
      (void*)realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > realarray->valssize )
   {
      SCIP_Real* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = calcGrowSize(arraygrowinit, arraygrowfac, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(realarray->blkmem, &newvals, newvalssize) );
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

         /* check for possible overflow or negative value */
         assert(realarray->maxusedidx - realarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[realarray->minusedidx - newfirstidx],
            &(realarray->vals[realarray->minusedidx - realarray->firstidx]),
            realarray->maxusedidx - realarray->minusedidx + 1); /*lint !e866*/
         for( i = realarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = 0.0;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = 0.0;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(realarray->blkmem, &realarray->vals, realarray->valssize);
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
SCIP_RETCODE SCIPrealarrayClear(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   SCIPdebugMessage("clearing realarray %p (firstidx=%d, size=%d, range=[%d,%d])\n",
      (void*)realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx);

   if( realarray->minusedidx <= realarray->maxusedidx )
   {
      assert(realarray->firstidx <= realarray->minusedidx);
      assert(realarray->maxusedidx < realarray->firstidx + realarray->valssize);
      assert(realarray->firstidx != -1);
      assert(realarray->valssize > 0);

      /* clear the used part of array */
      BMSclearMemoryArray(&realarray->vals[realarray->minusedidx - realarray->firstidx],
         realarray->maxusedidx - realarray->minusedidx + 1); /*lint !e866*/

      /* mark the array cleared */
      realarray->minusedidx = INT_MAX;
      realarray->maxusedidx = INT_MIN;
   }
   assert(realarray->minusedidx == INT_MAX);
   assert(realarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
SCIP_Real SCIPrealarrayGetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx                 /**< array index to get value for */
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
SCIP_RETCODE SCIPrealarraySetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   SCIP_Real             eps,                /**< epsilon value */
   int                   idx,                /**< array index to set value for */
   SCIP_Real             val                 /**< value to set array index to */
   )
{
   assert(realarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting realarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %g\n",
      (void*)realarray, realarray->firstidx, realarray->valssize, realarray->minusedidx, realarray->maxusedidx, idx, val);

   if( val != 0.0 )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPrealarrayExtend(realarray, arraygrowinit, arraygrowfac, idx, idx) );
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
            && realarray->vals[realarray->minusedidx - realarray->firstidx] == 0.0 );

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
         while( realarray->vals[realarray->maxusedidx - realarray->firstidx] == 0.0 );
      }
   }

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPrealarrayIncVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   SCIP_Real             eps,                /**< epsilon value */
   int                   idx,                /**< array index to increase value for */
   SCIP_Real             incval              /**< value to increase array index */
   )
{
   SCIP_Real oldval;

   oldval = SCIPrealarrayGetVal(realarray, idx);
   if( oldval != SCIP_INVALID ) /*lint !e777*/
      return SCIPrealarraySetVal(realarray, arraygrowinit, arraygrowfac, eps, idx, oldval + incval);
   else
      return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPrealarrayGetMinIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   return realarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPrealarrayGetMaxIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   assert(realarray != NULL);

   return realarray->maxusedidx;
}

/** creates a dynamic array of int values */
SCIP_RETCODE SCIPintarrayCreate(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the int array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(intarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, intarray) );
   (*intarray)->blkmem = blkmem;
   (*intarray)->vals = NULL;
   (*intarray)->valssize = 0;
   (*intarray)->firstidx = -1;
   (*intarray)->minusedidx = INT_MAX;
   (*intarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of int values */
SCIP_RETCODE SCIPintarrayCopy(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the copied int array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_INTARRAY*        sourceintarray      /**< dynamic int array to copy */
   )
{
   assert(intarray != NULL);
   assert(sourceintarray != NULL);

   SCIP_CALL( SCIPintarrayCreate(intarray, blkmem) );
   if( sourceintarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*intarray)->vals, sourceintarray->vals, sourceintarray->valssize) );
   }
   (*intarray)->valssize = sourceintarray->valssize;
   (*intarray)->firstidx = sourceintarray->firstidx;
   (*intarray)->minusedidx = sourceintarray->minusedidx;
   (*intarray)->maxusedidx = sourceintarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of int values */
SCIP_RETCODE SCIPintarrayFree(
   SCIP_INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   assert(intarray != NULL);
   assert(*intarray != NULL);

   BMSfreeBlockMemoryArrayNull((*intarray)->blkmem, &(*intarray)->vals, (*intarray)->valssize);
   BMSfreeBlockMemory((*intarray)->blkmem, intarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPintarrayExtend(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
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

   SCIPdebugMessage("extending intarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      (void*)intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > intarray->valssize )
   {
      int* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = calcGrowSize(arraygrowinit, arraygrowfac, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(intarray->blkmem, &newvals, newvalssize) );
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

         /* check for possible overflow or negative value */
         assert(intarray->maxusedidx - intarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[intarray->minusedidx - newfirstidx],
            &intarray->vals[intarray->minusedidx - intarray->firstidx],
            intarray->maxusedidx - intarray->minusedidx + 1); /*lint !e866*/
         for( i = intarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = 0;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = 0;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(intarray->blkmem, &intarray->vals, intarray->valssize);
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
SCIP_RETCODE SCIPintarrayClear(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   SCIPdebugMessage("clearing intarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      (void*)intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx);

   if( intarray->minusedidx <= intarray->maxusedidx )
   {
      assert(intarray->firstidx <= intarray->minusedidx);
      assert(intarray->maxusedidx < intarray->firstidx + intarray->valssize);
      assert(intarray->firstidx != -1);
      assert(intarray->valssize > 0);

      /* clear the used part of array */
      BMSclearMemoryArray(&intarray->vals[intarray->minusedidx - intarray->firstidx],
         intarray->maxusedidx - intarray->minusedidx + 1); /*lint !e866*/

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
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
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
SCIP_RETCODE SCIPintarraySetVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   int                   val                 /**< value to set array index to */
   )
{
   assert(intarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting intarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %d\n", 
      (void*)intarray, intarray->firstidx, intarray->valssize, intarray->minusedidx, intarray->maxusedidx, idx, val);

   if( val != 0 )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPintarrayExtend(intarray, arraygrowinit, arraygrowfac, idx, idx) );
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
SCIP_RETCODE SCIPintarrayIncVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
   int                   incval              /**< value to increase array index */
   )
{
   return SCIPintarraySetVal(intarray, arraygrowinit, arraygrowfac, idx, SCIPintarrayGetVal(intarray, idx) + incval);
}

/** returns the minimal index of all stored non-zero elements */
int SCIPintarrayGetMinIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   return intarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPintarrayGetMaxIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   assert(intarray != NULL);

   return intarray->maxusedidx;
}


/** creates a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayCreate(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the bool array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(boolarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, boolarray) );
   (*boolarray)->blkmem = blkmem;
   (*boolarray)->vals = NULL;
   (*boolarray)->valssize = 0;
   (*boolarray)->firstidx = -1;
   (*boolarray)->minusedidx = INT_MAX;
   (*boolarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayCopy(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the copied bool array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_BOOLARRAY*       sourceboolarray     /**< dynamic bool array to copy */
   )
{
   assert(boolarray != NULL);
   assert(sourceboolarray != NULL);

   SCIP_CALL( SCIPboolarrayCreate(boolarray, blkmem) );
   if( sourceboolarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*boolarray)->vals, sourceboolarray->vals, 
                     sourceboolarray->valssize) );
   }
   (*boolarray)->valssize = sourceboolarray->valssize;
   (*boolarray)->firstidx = sourceboolarray->firstidx;
   (*boolarray)->minusedidx = sourceboolarray->minusedidx;
   (*boolarray)->maxusedidx = sourceboolarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayFree(
   SCIP_BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   assert(boolarray != NULL);
   assert(*boolarray != NULL);

   BMSfreeBlockMemoryArrayNull((*boolarray)->blkmem, &(*boolarray)->vals, (*boolarray)->valssize);
   BMSfreeBlockMemory((*boolarray)->blkmem, boolarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPboolarrayExtend(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
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

   SCIPdebugMessage("extending boolarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      (void*)boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > boolarray->valssize )
   {
      SCIP_Bool* newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = calcGrowSize(arraygrowinit, arraygrowfac, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(boolarray->blkmem, &newvals, newvalssize) );
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

         /* check for possible overflow or negative value */
         assert(boolarray->maxusedidx - boolarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[boolarray->minusedidx - newfirstidx],
            &boolarray->vals[boolarray->minusedidx - boolarray->firstidx],
            boolarray->maxusedidx - boolarray->minusedidx + 1); /*lint !e866*/
         for( i = boolarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = FALSE;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = FALSE;
      }

      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(boolarray->blkmem, &boolarray->vals, boolarray->valssize);
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

	 assert(0 <= boolarray->minusedidx - boolarray->firstidx - shift);
         assert(boolarray->maxusedidx - boolarray->firstidx - shift < boolarray->valssize);
	 BMSmoveMemoryArray(&(boolarray->vals[boolarray->minusedidx - boolarray->firstidx - shift]),
            &(boolarray->vals[boolarray->minusedidx - boolarray->firstidx]),
            boolarray->maxusedidx - boolarray->minusedidx + 1); /*lint !e866*/

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
SCIP_RETCODE SCIPboolarrayClear(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   SCIPdebugMessage("clearing boolarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      (void*)boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx);

   if( boolarray->minusedidx <= boolarray->maxusedidx )
   {
      assert(boolarray->firstidx <= boolarray->minusedidx);
      assert(boolarray->maxusedidx < boolarray->firstidx + boolarray->valssize);
      assert(boolarray->firstidx != -1);
      assert(boolarray->valssize > 0);

      /* clear the used part of array */
      BMSclearMemoryArray(&boolarray->vals[boolarray->minusedidx - boolarray->firstidx],
         boolarray->maxusedidx - boolarray->minusedidx + 1); /*lint !e866*/

      /* mark the array cleared */
      boolarray->minusedidx = INT_MAX;
      boolarray->maxusedidx = INT_MIN;
   }
   assert(boolarray->minusedidx == INT_MAX);
   assert(boolarray->maxusedidx == INT_MIN);

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
SCIP_Bool SCIPboolarrayGetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx                 /**< array index to get value for */
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
SCIP_RETCODE SCIPboolarraySetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   SCIP_Bool             val                 /**< value to set array index to */
   )
{
   assert(boolarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting boolarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %u\n", 
      (void*)boolarray, boolarray->firstidx, boolarray->valssize, boolarray->minusedidx, boolarray->maxusedidx, idx, val);

   if( val != FALSE )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPboolarrayExtend(boolarray, arraygrowinit, arraygrowfac, idx, idx) );
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
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   return boolarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPboolarrayGetMaxIdx(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   assert(boolarray != NULL);

   return boolarray->maxusedidx;
}


/** creates a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayCreate(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the ptr array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(ptrarray != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, ptrarray) );
   (*ptrarray)->blkmem = blkmem;
   (*ptrarray)->vals = NULL;
   (*ptrarray)->valssize = 0;
   (*ptrarray)->firstidx = -1;
   (*ptrarray)->minusedidx = INT_MAX;
   (*ptrarray)->maxusedidx = INT_MIN;

   return SCIP_OKAY;
}

/** creates a copy of a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayCopy(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the copied ptr array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PTRARRAY*        sourceptrarray      /**< dynamic ptr array to copy */
   )
{
   assert(ptrarray != NULL);
   assert(sourceptrarray != NULL);

   SCIP_CALL( SCIPptrarrayCreate(ptrarray, blkmem) );
   if( sourceptrarray->valssize > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*ptrarray)->vals, sourceptrarray->vals, sourceptrarray->valssize) );
   }
   (*ptrarray)->valssize = sourceptrarray->valssize;
   (*ptrarray)->firstidx = sourceptrarray->firstidx;
   (*ptrarray)->minusedidx = sourceptrarray->minusedidx;
   (*ptrarray)->maxusedidx = sourceptrarray->maxusedidx;

   return SCIP_OKAY;
}

/** frees a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayFree(
   SCIP_PTRARRAY**       ptrarray            /**< pointer to the ptr array */
   )
{
   assert(ptrarray != NULL);
   assert(*ptrarray != NULL);

   BMSfreeBlockMemoryArrayNull((*ptrarray)->blkmem, &(*ptrarray)->vals, (*ptrarray)->valssize);
   BMSfreeBlockMemory((*ptrarray)->blkmem, ptrarray);

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPptrarrayExtend(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
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

   SCIPdebugMessage("extending ptrarray %p (firstidx=%d, size=%d, range=[%d,%d]) to range [%d,%d]\n", 
      (void*)ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx, minidx, maxidx);

   /* check, whether we have to allocate additional memory, or shift the array */
   nused = maxidx - minidx + 1;
   if( nused > ptrarray->valssize )
   {
      void** newvals;
      int newvalssize;

      /* allocate new memory storage */
      newvalssize = calcGrowSize(arraygrowinit, arraygrowfac, nused);
      SCIP_ALLOC( BMSallocBlockMemoryArray(ptrarray->blkmem, &newvals, newvalssize) );
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

         /* check for possible overflow or negative value */
         assert(ptrarray->maxusedidx - ptrarray->minusedidx + 1 > 0);

         BMScopyMemoryArray(&newvals[ptrarray->minusedidx - newfirstidx],
            &(ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx]),
            ptrarray->maxusedidx - ptrarray->minusedidx + 1); /*lint !e866*/
         for( i = ptrarray->maxusedidx - newfirstidx + 1; i < newvalssize; ++i )
            newvals[i] = NULL;
      }
      else
      {
         for( i = 0; i < newvalssize; ++i )
            newvals[i] = NULL;
      }
      
      /* free old memory storage, and set the new array parameters */
      BMSfreeBlockMemoryArrayNull(ptrarray->blkmem, &ptrarray->vals, ptrarray->valssize);
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
SCIP_RETCODE SCIPptrarrayClear(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   SCIPdebugMessage("clearing ptrarray %p (firstidx=%d, size=%d, range=[%d,%d])\n", 
      (void*)ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx);

   if( ptrarray->minusedidx <= ptrarray->maxusedidx )
   {
      assert(ptrarray->firstidx <= ptrarray->minusedidx);
      assert(ptrarray->maxusedidx < ptrarray->firstidx + ptrarray->valssize);
      assert(ptrarray->firstidx != -1);
      assert(ptrarray->valssize > 0);

      /* clear the used part of array */
      BMSclearMemoryArray(&ptrarray->vals[ptrarray->minusedidx - ptrarray->firstidx],
         ptrarray->maxusedidx - ptrarray->minusedidx + 1); /*lint !e866*/

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
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   idx                 /**< array index to get value for */
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
SCIP_RETCODE SCIPptrarraySetVal(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   void*                 val                 /**< value to set array index to */
   )
{
   assert(ptrarray != NULL);
   assert(idx >= 0);

   SCIPdebugMessage("setting ptrarray %p (firstidx=%d, size=%d, range=[%d,%d]) index %d to %p\n", 
      (void*)ptrarray, ptrarray->firstidx, ptrarray->valssize, ptrarray->minusedidx, ptrarray->maxusedidx, idx, val);

   if( val != NULL )
   {
      /* extend array to be able to store the index */
      SCIP_CALL( SCIPptrarrayExtend(ptrarray, arraygrowinit, arraygrowfac, idx, idx) );
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
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   return ptrarray->minusedidx;
}

/** returns the maximal index of all stored non-zero elements */
int SCIPptrarrayGetMaxIdx(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   assert(ptrarray != NULL);

   return ptrarray->maxusedidx;
}




/*
 * Sorting algorithms
 */

/** default comparer for integers */
SCIP_DECL_SORTPTRCOMP(SCIPsortCompInt)
{
   int value1;
   int value2;

   value1 = (int)(size_t)elem1;
   value2 = (int)(size_t)elem2;

   if( value1 < value2 )
      return -1;

   if( value2 < value1 )
      return 1;

   return 0;
}

/* first all upwards-sorting methods */

/** sort an indexed element set in non-decreasing order, resulting in a permutation index array */
void SCIPsort(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   )
{
   int pos;

   assert(indcomp != NULL);
   assert(len == 0 || perm != NULL);

   /* create identity permutation */
   for( pos = 0; pos < len; ++pos )
      perm[pos] = pos;

   SCIPsortInd(perm, indcomp, dataptr, len);
}

/* SCIPsortInd(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Ind
#define SORTTPL_KEYTYPE     int
#define SORTTPL_INDCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Ptr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrRealIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrLongInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrLongInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrLongIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrPtrLongIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Real
#define SORTTPL_KEYTYPE     SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealIntLong
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealIntPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealLongRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealLongRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Longint
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortRealRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealRealRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortRealRealRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     RealRealRealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Int
#define SORTTPL_KEYTYPE     int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntLong
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtrReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntIntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtrIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtrIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Real
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     Long
#define SORTTPL_KEYTYPE     SCIP_Longint
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtr
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrPtrIntInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortLongPtrPtrBoolInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     LongPtrPtrBoolInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  int
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrIntIntBoolBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     PtrIntIntBoolBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortIntPtrIntIntBoolBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     IntPtrIntIntBoolBool
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_FIELD5TYPE  SCIP_Bool
#include "scip/sorttpl.c" /*lint !e451*/


/* now all downwards-sorting methods */


/** sort an indexed element set in non-increasing order, resulting in a permutation index array */
void SCIPsortDown(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   )
{
   int pos;

   assert(indcomp != NULL);
   assert(len == 0 || perm != NULL);

   /* create identity permutation */
   for( pos = 0; pos < len; ++pos )
      perm[pos] = pos;
   
   SCIPsortDownInd(perm, indcomp, dataptr, len);
}


/* SCIPsortDownInd(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownInd
#define SORTTPL_KEYTYPE     int
#define SORTTPL_INDCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtr
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownPtrBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/

/* SCIPsortDownPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrReal
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrRealIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortPtrPtrRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrRealInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrLongInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrLongInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrPtrLongIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrPtrLongIntInt
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownReal
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Bool
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealIntLong
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealIntPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealPtrPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealPtrPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealLongRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealLongRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Longint
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownRealRealRealBoolPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownRealRealRealBoolPtr
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntLong
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  SCIP_Longint
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntIntIntPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntIntIntPtr
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntPtrIntReal(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntPtrIntReal
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Real
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLong(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLong
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtr(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtr
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrPtrInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrPtrInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrPtrIntInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrPtrIntInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownLongPtrPtrBoolInt(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownLongPtrPtrBoolInt
#define SORTTPL_KEYTYPE     SCIP_Longint
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  void*
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  int
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownPtrIntIntBoolBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownPtrIntIntBoolBool
#define SORTTPL_KEYTYPE     void*
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  SCIP_Bool
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_PTRCOMP
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/* SCIPsortDownIntPtrIntIntBoolBool(), SCIPsortedvecInsert...(), SCIPsortedvecDelPos...(), SCIPsortedvecFind...() via sort template */
#define SORTTPL_NAMEEXT     DownIntPtrIntIntBoolBool
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  void*
#define SORTTPL_FIELD2TYPE  int
#define SORTTPL_FIELD3TYPE  int
#define SORTTPL_FIELD4TYPE  SCIP_Bool
#define SORTTPL_FIELD5TYPE  SCIP_Bool
#define SORTTPL_BACKWARDS
#include "scip/sorttpl.c" /*lint !e451*/


/*
 * Resource Profile
 */

/** creates resource profile */
SCIP_RETCODE SCIPprofileCreate(
   SCIP_PROFILE**        profile,            /**< pointer to store the resource profile */
   int                   capacity            /**< resource capacity */
   )
{
   assert(profile != NULL);
   assert(capacity > 0);

   SCIP_ALLOC( BMSallocMemory(profile) );

   (*profile)->arraysize = 10;
   SCIP_ALLOC( BMSallocMemoryArray(&(*profile)->timepoints, (*profile)->arraysize) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*profile)->loads, (*profile)->arraysize) );

   /* setup resource profile for use */
   (*profile)->ntimepoints = 1;
   (*profile)->timepoints[0] = 0;
   (*profile)->loads[0] = 0;
   (*profile)->capacity = capacity;

   return SCIP_OKAY;
}

/** frees given resource profile */
void SCIPprofileFree(
   SCIP_PROFILE**        profile             /**< pointer to the resource profile */
   )
{
   assert(profile != NULL);
   assert(*profile != NULL);

   /* free main hash map data structure */
   BMSfreeMemoryArray(&(*profile)->loads);
   BMSfreeMemoryArray(&(*profile)->timepoints);
   BMSfreeMemory(profile);
}

/** output of the given resource profile */
void SCIPprofilePrint(
   SCIP_PROFILE*         profile,            /**< resource profile to output */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int t;

   SCIPmessageFPrintInfo(messagehdlr, file, "Profile <%p> (capacity %d) --> ", profile, profile->capacity);

   for( t = 0; t < profile->ntimepoints; ++t )
   {
      if( t == 0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "%d:(%d,%d)", t, profile->timepoints[t], profile->loads[t]);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, ", %d:(%d,%d)", t, profile->timepoints[t], profile->loads[t]);
   }

   SCIPmessageFPrintInfo(messagehdlr, file,"\n");
}

/** returns the capacity of the resource profile */
int SCIPprofileGetCapacity(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   )
{
   assert(profile != NULL);

   return profile->capacity;
}

/** returns the number time points of the resource profile */
int SCIPprofileGetNTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   )
{
   assert(profile != NULL);

   return profile->ntimepoints;
}

/** returns the time points of the resource profile */
int* SCIPprofileGetTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   )
{
   assert(profile != NULL);

   return profile->timepoints;
}

/** returns the loads of the resource profile */
int* SCIPprofileGetLoads(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   )
{
   assert(profile != NULL);

   return profile->loads;
}

/** returns the time point for given position of the resource profile */
int SCIPprofileGetTime(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   pos                 /**< position */
   )
{
   assert(profile != NULL);
   assert(pos >= 0 && pos < profile->ntimepoints);

   return profile->timepoints[pos];
}

/** returns the loads of the resource profile at the given position */
int SCIPprofileGetLoad(
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   pos                 /**< position */
   )
{
   assert(profile != NULL);
   assert(pos >= 0 && pos < profile->ntimepoints);

   return profile->loads[pos];
}

/** returns if the given time point exists in the resource profile and stores the position of the given time point if it
 *  exists; otherwise the position of the next smaller existing time point is stored
 */
SCIP_Bool SCIPprofileFindLeft(
   SCIP_PROFILE*         profile,            /**< resource profile to search */
   int                   timepoint,          /**< time point to search for */
   int*                  pos                 /**< pointer to store the position */
   )
{
   assert(profile != NULL);
   assert(timepoint >= 0);
   assert(profile->ntimepoints > 0);
   assert(profile->timepoints[0] == 0);

   /* find the position of time point in the time points array via binary search */
   if( SCIPsortedvecFindInt(profile->timepoints, timepoint, profile->ntimepoints, pos) )
      return TRUE;

   assert(*pos > 0);
   (*pos)--;

   return FALSE;
}

/* ensures that resource profile arrays is big enough */
static
SCIP_RETCODE ensureProfileSize(
   SCIP_PROFILE*         profile,            /**< resource profile to insert the time point */
   int                   neededsize          /**< needed size */
   )
{
   assert(profile->arraysize > 0);

   /* check whether the arrays are big enough */
   if( neededsize <= profile->arraysize )
      return SCIP_OKAY;

   profile->arraysize *= 2;

   SCIP_ALLOC( BMSreallocMemoryArray(&profile->timepoints, profile->arraysize) );
   SCIP_ALLOC( BMSreallocMemoryArray(&profile->loads, profile->arraysize) );

   return SCIP_OKAY;
}

/** inserts the given time point into the resource profile if it this time point does not exists yet; returns its
 *  position in the time point array
 */
static
SCIP_RETCODE profileInsertTimepoint(
   SCIP_PROFILE*         profile,            /**< resource profile to insert the time point */
   int                   timepoint,          /**< time point to insert */
   int*                  pos                 /**< pointer to store the insert position */
   )
{
   assert(profile != NULL);
   assert(timepoint >= 0);
   assert(profile->arraysize >= profile->ntimepoints);

   /* get the position of the given time point in the resource profile array if it exists; otherwise the position of the
    * next smaller existing time point
    */
   if( !SCIPprofileFindLeft(profile, timepoint, pos) )
   {
      assert(*pos >= 0 && *pos < profile->ntimepoints);
      assert(timepoint >= profile->timepoints[*pos]);

      /* ensure that the arrays are big enough */
      SCIP_CALL( ensureProfileSize(profile, profile->ntimepoints + 1) );
      assert(profile->arraysize > profile->ntimepoints);

      /* insert new time point into the (sorted) resource profile */
      SCIPsortedvecInsertIntInt(profile->timepoints, profile->loads, timepoint, profile->loads[*pos],
         &profile->ntimepoints, pos);
   }

#ifndef NDEBUG
   /* check if the time points are sorted */
   {
      int i;
      for( i = 1; i < profile->ntimepoints; ++i )
         assert(profile->timepoints[i-1] < profile->timepoints[i]);
   }
#endif

   return SCIP_OKAY;
}

/** updates the resource profile due to inserting of a core */
static
SCIP_RETCODE profileUpdate(
   SCIP_PROFILE*         profile,            /**< resource profile to update */
   int                   left,               /**< left side of core interval */
   int                   right,              /**< right side of core interval */
   int                   demand,             /**< demand of the core */
   int*                  pos,                /**< pointer to store the first position were it gets infeasible */
   SCIP_Bool*            infeasible          /**< pointer to store if the update is infeasible */
   )
{
   int startpos;
   int endpos;
   int i;

   assert(profile != NULL);
   assert(profile->arraysize >= profile->ntimepoints);
   assert(left >= 0);
   assert(left < right);
   assert(infeasible != NULL);

   (*infeasible) = FALSE;
   (*pos) = -1;

   /* get position of the starttime in profile */
   SCIP_CALL( profileInsertTimepoint(profile, left, &startpos) );
   assert(profile->timepoints[startpos] == left);

   /* get position of the endtime in profile */
   SCIP_CALL( profileInsertTimepoint(profile, right, &endpos) );
   assert(profile->timepoints[endpos] == right);

   assert(startpos < endpos);
   assert(profile->arraysize >= profile->ntimepoints);

   /* remove/add the given demand from the core */
   for( i = startpos; i < endpos; ++i )
   {
      profile->loads[i] += demand;

      /* check if the core fits */
      if( profile->loads[i] > profile->capacity )
      {
         SCIPdebugMessage("core insertion detected infeasibility (pos %d)\n", i);

         (*infeasible) = TRUE;
         (*pos) = i;

         /* remove the partly inserted core since it does fit completely */
         for( ; i >= startpos; --i ) /*lint !e445*/
            profile->loads[i] -= demand;

         break;
      }
   }

   return SCIP_OKAY;
}

/** insert a core into resource profile; if the core is non-empty the resource profile will be updated otherwise nothing
 *  happens
 */
SCIP_RETCODE SCIPprofileInsertCore(
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   demand,             /**< demand of the core */
   int*                  pos,                /**< pointer to store the first position were it gets infeasible */
   SCIP_Bool*            infeasible          /**< pointer to store if the core does not fit due to capacity */
   )
{
   assert(profile != NULL);
   assert(left < right);
   assert(demand >= 0);
   assert(infeasible != NULL);

   (*infeasible) = FALSE;
   (*pos) = -1;

   /* insert core into the resource profile */
   SCIPdebugMessage("insert core [%d,%d] with demand %d\n", left, right, demand);

   if( demand > 0 )
   {
      /* try to insert core into the resource profile */
      SCIP_CALL( profileUpdate(profile, left, right, demand, pos, infeasible) );
   }

   return SCIP_OKAY;
}

/** subtracts the demand from the resource profile during core time */
SCIP_RETCODE SCIPprofileDeleteCore(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   demand              /**< demand of the core */
   )
{
   SCIP_Bool infeasible;
   int pos;

   assert(left < right);
#ifndef NDEBUG
   {
      /* check if the left and right time points of the core correspond to a time point in the resource profile; this
       * should be the case since we added the core before to the resource profile
       */
      assert(SCIPprofileFindLeft(profile, left, &pos));
      assert(SCIPprofileFindLeft(profile, right, &pos));
   }
#endif

   /* remove the core from the resource profile */
   SCIPdebugMessage("delete core [%d,%d] with demand %d\n", left, right, demand);

   SCIP_CALL( profileUpdate(profile, left, right, -demand, &pos, &infeasible) );
   assert(!infeasible);

   return SCIP_OKAY;
}

/** returns TRUE if the core (given by its demand and during) can be inserted at the given time point; otherwise FALSE */
static
int profileFindFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   pos,                /**< pointer to store the position in the profile to start the serch */
   int                   lst,                /**< latest start time */
   int                   duration,           /**< duration of the core */
   int                   demand,             /**< demand of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the corer cannot be inserted */
   )
{
   int remainingduration;
   int startpos;

   assert(profile != NULL);
   assert(pos >= 0);
   assert(pos < profile->ntimepoints);
   assert(duration > 0);
   assert(demand > 0);
   assert(profile->loads[profile->ntimepoints-1] == 0);

   remainingduration = duration;
   startpos = pos;
   (*infeasible) = FALSE;

   if( profile->timepoints[startpos] > lst )
   {
      (*infeasible) = TRUE;
      return pos;
   }

   while( pos < profile->ntimepoints - 1 )
   {
      if( profile->loads[pos] + demand > profile->capacity )
      {
         SCIPdebugMessage("profile <%p>: core does not fit at time point %d (pos %d)\n", (void*)profile, profile->timepoints[pos], pos);
         startpos = pos + 1;
         remainingduration = duration;

         if( profile->timepoints[startpos] > lst )
         {
            (*infeasible) = TRUE;
            return pos;
         }
      }
      else
         remainingduration -= profile->timepoints[pos+1] - profile->timepoints[pos];

      if( remainingduration <= 0 )
         break;

      pos++;
   }

   return startpos;
}

/** return the earliest possible starting point within the time interval [lb,ub] for a given core (given by its demand
 *  and duration)
 */
int SCIPprofileGetEarliestFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   est,                /**< earliest starting time of the given core */
   int                   lst,                /**< latest starting time of the given core */
   int                   duration,           /**< duration of the core */
   int                   demand,             /**< demand of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the corer cannot be inserted */
   )
{
   SCIP_Bool found;
   int pos;

   assert(profile != NULL);
   assert(est >= 0);
   assert(est <= lst);
   assert(duration >= 0);
   assert(demand >= 0);
   assert(infeasible != NULL);
   assert(profile->ntimepoints > 0);
   assert(profile->loads[profile->ntimepoints-1] == 0);

   SCIPdebugMessage("profile <%p>: find earliest start time (demad %d, duration %d) [%d,%d]\n", (void*)profile, demand, duration, est, lst);

   if( duration == 0 || demand == 0 )
   {
      *infeasible = FALSE;
      return est;
   }

   found = SCIPprofileFindLeft(profile, est, &pos);
   SCIPdebugMessage("profile <%p>: earliest start time does %s exist as time point (pos %d)\n", (void*)profile, found ? "" : "not", pos);

   /* if the position is the last time point in the profile, the core can be inserted at its earliest start time */
   if( pos == profile->ntimepoints - 1 )
   {
      (*infeasible) = FALSE;
      return est;
   }

   if( found )
   {
      /* if the start time matches a time point in the profile we can just search */
      assert(profile->timepoints[pos] == est);
      pos = profileFindFeasibleStart(profile, pos, lst, duration, demand, infeasible);

      assert(pos < profile->ntimepoints);
      est = profile->timepoints[pos];
   }
   else if( profile->loads[pos] + demand > profile->capacity )
   {
      /* if the the time point left to the start time has not enough free capacity we can just search the profile
       * starting from the next time point
       */
      assert(profile->timepoints[pos] <= est);
      pos = profileFindFeasibleStart(profile, pos+1, lst, duration, demand, infeasible);

      assert(pos < profile->ntimepoints);
      est = profile->timepoints[pos];
   }
   else
   {
      int remainingduration;

      /* check if the core can be placed at its earliest start time */

      assert(pos < profile->ntimepoints - 1);

      remainingduration = duration - (profile->timepoints[pos+1] - est);
      SCIPdebugMessage("remaining duration %d\n", remainingduration);


      if( remainingduration <= 0 )
         (*infeasible) = FALSE;
      else
      {
         pos = profileFindFeasibleStart(profile, pos+1, profile->timepoints[pos+1], remainingduration, demand, infeasible);
         SCIPdebugMessage("remaining duration can%s be processed\n", *infeasible ? "not" : "");

         if( *infeasible )
         {
            pos = profileFindFeasibleStart(profile, pos+1, lst, duration, demand, infeasible);

            assert(pos < profile->ntimepoints);
            est = profile->timepoints[pos];
         }
      }
   }

   return est;
}

/** returns TRUE if the core (given by its demand and during) can be inserted at the given time point; otherwise FALSE */
static
int profileFindDownFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   pos,                /**< pointer to store the position in the profile to start the search */
   int                   ect,                /**< earliest completion time */
   int                   duration,           /**< duration of the core */
   int                   demand,             /**< demand of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the corer cannot be inserted */
   )
{
   int remainingduration;
   int endpos;

   assert(profile != NULL);
   assert(pos >= 0);
   assert(pos < profile->ntimepoints);
   assert(duration > 0);
   assert(demand > 0);
   assert(profile->ntimepoints > 0);
   assert(profile->loads[profile->ntimepoints-1] == 0);

   remainingduration = duration;
   endpos = pos;
   (*infeasible) = TRUE;

   if( profile->timepoints[endpos] < ect - duration )
      return pos;

   while( pos > 0 )
   {
      if( profile->loads[pos-1] + demand > profile->capacity )
      {
         SCIPdebugMessage("profile <%p>: core does not fit at time point %d (pos %d)\n", (void*)profile, profile->timepoints[pos-1], pos-1);

         endpos = pos - 1;
         remainingduration = duration;

         if( profile->timepoints[endpos] < ect - duration )
            return pos;
      }
      else
         remainingduration -= profile->timepoints[pos] - profile->timepoints[pos-1];

      if( remainingduration <= 0 )
      {
         *infeasible = FALSE;
         break;
      }

      pos--;
   }

   return endpos;
}

/** return the latest possible starting point within the time interval [lb,ub] for a given core (given by its demand and
 *  duration)
 */
int SCIPprofileGetLatestFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   est,                /**< earliest possible start point */
   int                   lst,                /**< latest possible start point */
   int                   duration,           /**< duration of the core */
   int                   demand,             /**< demand of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the core cannot be inserted */
   )
{
   SCIP_Bool found;
   int ect;
   int lct;
   int pos;

   assert(profile != NULL);
   assert(est >= 0);
   assert(est <= lst);
   assert(duration >= 0);
   assert(demand >= 0);
   assert(infeasible != NULL);
   assert(profile->ntimepoints > 0);
   assert(profile->loads[profile->ntimepoints-1] == 0);

   if( duration == 0 || demand == 0 )
   {
      *infeasible = FALSE;
      return lst;
   }

   ect = est + duration;
   lct = lst + duration;

   found = SCIPprofileFindLeft(profile, lct, &pos);
   SCIPdebugMessage("profile <%p>: latest completion time %d does %s exist as time point (pos %d)\n", (void*)profile, lct, found ? "" : "not", pos);

   if( found )
   {
      /* if the start time matches a time point in the profile we can just search */
      assert(profile->timepoints[pos] == lct);
      pos = profileFindDownFeasibleStart(profile, pos, ect, duration, demand, infeasible);

      assert(pos < profile->ntimepoints && pos >= 0);
      lct = profile->timepoints[pos];
   }
   else if( profile->loads[pos] + demand > profile->capacity )
   {
      /* if the time point left to the start time has not enough free capacity we can just search the profile starting
       * from the next time point
       */
      assert(profile->timepoints[pos] < lct);
      pos = profileFindDownFeasibleStart(profile, pos, ect, duration, demand, infeasible);

      assert(pos < profile->ntimepoints && pos >= 0);
      lct = profile->timepoints[pos];
   }
   else
   {
      int remainingduration;

      /* check if the core can be placed at its latest start time */
      assert(profile->timepoints[pos] < lct);

      remainingduration = duration - (lct - profile->timepoints[pos]);

      if( remainingduration <= 0 )
         (*infeasible) = FALSE;
      else
      {
         pos = profileFindDownFeasibleStart(profile, pos, profile->timepoints[pos], remainingduration, demand, infeasible);

         if( *infeasible )
         {
            pos = profileFindDownFeasibleStart(profile, pos, ect, duration, demand, infeasible);

            assert(pos < profile->ntimepoints && pos >= 0);
            lct = profile->timepoints[pos];
         }
      }
   }

   return lct - duration;
}

/*
 * Directed graph
 */

/** creates directed graph structure */
SCIP_RETCODE SCIPdigraphCreate(
   SCIP_DIGRAPH**        digraph,            /**< pointer to store the created directed graph */
   int                   nnodes              /**< number of nodes */
   )
{
   assert(digraph != NULL);
   assert(nnodes > 0);

   /* allocate memory for the graph and the arrays storing arcs and datas */
   SCIP_ALLOC( BMSallocMemory(digraph) );
   SCIP_ALLOC( BMSallocClearMemoryArray(&(*digraph)->successors, nnodes) );
   SCIP_ALLOC( BMSallocClearMemoryArray(&(*digraph)->arcdatas, nnodes) );
   SCIP_ALLOC( BMSallocClearMemoryArray(&(*digraph)->successorssize, nnodes) );
   SCIP_ALLOC( BMSallocClearMemoryArray(&(*digraph)->nsuccessors, nnodes) );

   /* store number of nodes */
   (*digraph)->nnodes = nnodes;

   /* at the beginning, no components are stored */
   (*digraph)->ncomponents = 0;
   (*digraph)->componentstartsize = 0;
   (*digraph)->components = NULL;
   (*digraph)->componentstarts = NULL;

   return SCIP_OKAY;
}

/** copies directed graph structure */
SCIP_RETCODE SCIPdigraphCopy(
   SCIP_DIGRAPH**        targetdigraph,      /**< pointer to store the copied directed graph */
   SCIP_DIGRAPH*         sourcedigraph       /**< source directed graph */
   )
{
   int ncomponents;
   int nnodes;
   int i;

   SCIP_ALLOC( BMSallocMemory(targetdigraph) );

   nnodes = sourcedigraph->nnodes;
   ncomponents = sourcedigraph->ncomponents;
   (*targetdigraph)->nnodes = nnodes;
   (*targetdigraph)->ncomponents = ncomponents;

   /* copy arcs and datas */
   SCIP_ALLOC( BMSallocClearMemoryArray(&(*targetdigraph)->successors, nnodes) );
   SCIP_ALLOC( BMSallocClearMemoryArray(&(*targetdigraph)->arcdatas, nnodes) );

   /* copy lists of successors and arc datas */
   for( i = 0; i < nnodes; ++i )
   {
      if( sourcedigraph->nsuccessors[i] > 0 )
      {
         assert(sourcedigraph->successors[i] != NULL);
         assert(sourcedigraph->arcdatas[i] != NULL);
         SCIP_ALLOC( BMSduplicateMemoryArray(&((*targetdigraph)->successors[i]),
               sourcedigraph->successors[i], sourcedigraph->nsuccessors[i]) ); /*lint !e866*/
         SCIP_ALLOC( BMSduplicateMemoryArray(&((*targetdigraph)->arcdatas[i]),
               sourcedigraph->arcdatas[i], sourcedigraph->nsuccessors[i]) ); /*lint !e866*/
      }
   }
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*targetdigraph)->successorssize, sourcedigraph->nsuccessors, nnodes) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*targetdigraph)->nsuccessors, sourcedigraph->nsuccessors, nnodes) );

   /* copy component data */
   if( ncomponents > 0 )
   {
      SCIP_ALLOC( BMSduplicateMemoryArray(&(*targetdigraph)->components, sourcedigraph->components,
            sourcedigraph->componentstarts[ncomponents]) );
      SCIP_ALLOC( BMSduplicateMemoryArray(&(*targetdigraph)->componentstarts,
            sourcedigraph->componentstarts,ncomponents + 1) );
      (*targetdigraph)->componentstartsize = ncomponents + 1;
   }
   else
   {
      (*targetdigraph)->components = NULL;
      (*targetdigraph)->componentstarts = NULL;
      (*targetdigraph)->componentstartsize = 0;
   }

   return SCIP_OKAY;
}

/** sets the sizes of the successor lists for the nodes in a directed graph and allocates memory for the lists */
SCIP_RETCODE SCIPdigraphSetSizes(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int*                  sizes               /**< sizes of the successor lists */
   )
{
   int i;

   assert(digraph != NULL);
   assert(digraph->nnodes > 0);

   for( i = 0; i < digraph->nnodes; ++i )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&digraph->successors[i], sizes[i]) ); /*lint !e866*/
      SCIP_ALLOC( BMSallocMemoryArray(&digraph->arcdatas[i], sizes[i]) ); /*lint !e866*/
      digraph->successorssize[i] = sizes[i];
      digraph->nsuccessors[i] = 0;
   }

   return SCIP_OKAY;
}

/** frees given directed graph structure */
void SCIPdigraphFree(
   SCIP_DIGRAPH**        digraph             /**< pointer to the directed graph */
   )
{
   int i;

   assert(digraph != NULL);
   assert(*digraph != NULL);

   /* free arrays storing the successor nodes and arc datas */
   for( i = (*digraph)->nnodes - 1; i >= 0; --i )
   {
      BMSfreeMemoryArrayNull(&(*digraph)->successors[i]);
      BMSfreeMemoryArrayNull(&(*digraph)->arcdatas[i]);
   }

   /* free components structure */
   SCIPdigraphFreeComponents(*digraph);
   assert((*digraph)->ncomponents == 0);
   assert((*digraph)->componentstartsize == 0);
   assert((*digraph)->components == NULL);
   assert((*digraph)->componentstarts == NULL);

   /* free directed graph data structure */
   BMSfreeMemoryArray(&(*digraph)->successorssize);
   BMSfreeMemoryArray(&(*digraph)->nsuccessors);
   BMSfreeMemoryArray(&(*digraph)->successors);
   BMSfreeMemoryArray(&(*digraph)->arcdatas);

   BMSfreeMemory(digraph);
}

#define STARTSUCCESSORSSIZE 5

/* ensures that successors array of one node in a directed graph is big enough */
static
SCIP_RETCODE ensureSuccessorsSize(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   idx,                /**< index for which the size is ensured */
   int                   newsize             /**< needed size */
   )
{
   assert(digraph != NULL);
   assert(idx >= 0);
   assert(idx < digraph->nnodes);
   assert(newsize > 0);

   /* check whether array is big enough, and realloc, if needed */
   if( newsize > digraph->successorssize[idx] )
   {
      if( digraph->successorssize[idx] == 0 )
      {
         digraph->successorssize[idx] = STARTSUCCESSORSSIZE;
         SCIP_ALLOC( BMSallocMemoryArray(&digraph->successors[idx], digraph->successorssize[idx]) ); /*lint !e866*/
         SCIP_ALLOC( BMSallocMemoryArray(&digraph->arcdatas[idx], digraph->successorssize[idx]) ); /*lint !e866*/
      }
      else
      {
         digraph->successorssize[idx] = 2 * digraph->successorssize[idx];
         SCIP_ALLOC( BMSreallocMemoryArray(&digraph->successors[idx], digraph->successorssize[idx]) ); /*lint !e866*/
         SCIP_ALLOC( BMSreallocMemoryArray(&digraph->arcdatas[idx], digraph->successorssize[idx]) ); /*lint !e866*/
      }
   }

   return SCIP_OKAY;
}

/** add (directed) arc and a related data to the directed graph structure
 *
 *  @note if the arc is already contained, it is added a second time
 */
SCIP_RETCODE SCIPdigraphAddArc(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the arc */
   int                   endnode,            /**< start node of the arc */
   void*                 data                /**< data that should be stored for the arc; or NULL */
   )
{
   assert(digraph != NULL);
   assert(startnode >= 0);
   assert(endnode >= 0);
   assert(startnode < digraph->nnodes);
   assert(endnode < digraph->nnodes);

   SCIP_CALL( ensureSuccessorsSize(digraph, startnode, digraph->nsuccessors[startnode] + 1) );

   /* add arc */
   digraph->successors[startnode][digraph->nsuccessors[startnode]] = endnode;
   digraph->arcdatas[startnode][digraph->nsuccessors[startnode]] = data;
   digraph->nsuccessors[startnode]++;

   return SCIP_OKAY;
}

/** add (directed) arc to the directed graph structure, if it is not contained, yet
 *
 * @note if there already exists an arc from startnode to endnode, the new arc is not added,
 *       even if its data is different
 */
SCIP_RETCODE SCIPdigraphAddArcSafe(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the arc */
   int                   endnode,            /**< start node of the arc */
   void*                 data                /**< data that should be stored for the arc; or NULL */
   )
{
   int nsuccessors;
   int i;

   assert(digraph != NULL);
   assert(startnode >= 0);
   assert(endnode >= 0);
   assert(startnode < digraph->nnodes);
   assert(endnode < digraph->nnodes);

   nsuccessors = digraph->nsuccessors[startnode];

   /* search for the arc in existing arcs */
   for( i = 0; i < nsuccessors; ++i )
      if( digraph->successors[startnode][i] == endnode )
         return SCIP_OKAY;

   SCIP_CALL( ensureSuccessorsSize(digraph, startnode, nsuccessors + 1) );

   /* add arc */
   digraph->successors[startnode][nsuccessors] = endnode;
   digraph->arcdatas[startnode][nsuccessors] = data;
   ++(digraph->nsuccessors[startnode]);

   return SCIP_OKAY;
}

/** returns the number of nodes of the given digraph */
int SCIPdigraphGetNNodes(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   assert(digraph != NULL);

   return digraph->nnodes;
}

/** returns the total number of arcs in the given digraph */
int SCIPdigraphGetNArcs(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   int i;
   int narcs;

   assert(digraph != NULL);

   /* count number of arcs */
   narcs = 0;
   for( i = 0; i < digraph->nnodes; ++i )
      narcs += digraph->nsuccessors[i];

   return narcs;
}

/** returns the number of successor nodes of the given node */
int SCIPdigraphGetNSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the number of outgoing arcs is returned */
   )
{
   assert(digraph != NULL);
   assert(node >= 0);
   assert(node < digraph->nnodes);
   assert(digraph->nsuccessors[node] >= 0);
   assert(digraph->nsuccessors[node] <= digraph->successorssize[node]);

   return digraph->nsuccessors[node];
}

/** returns the array of indices of the successor nodes; this array must not be changed from outside */
int* SCIPdigraphGetSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the array of outgoing arcs is returned */
   )
{
   assert(digraph != NULL);
   assert(node >= 0);
   assert(node < digraph->nnodes);
   assert(digraph->nsuccessors[node] >= 0);
   assert(digraph->nsuccessors[node] <= digraph->successorssize[node]);
   assert((digraph->nsuccessors[node] == 0) || (digraph->successors[node] != NULL));

   return digraph->successors[node];
}

/** returns the array of datas corresponding to the arcs originating at the given node, or NULL if no data exist; this
 *  array must not be changed from outside
 */
void** SCIPdigraphGetSuccessorsDatas(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the data corresponding to the outgoing arcs is returned */
   )
{
   assert(digraph != NULL);
   assert(node >= 0);
   assert(node < digraph->nnodes);
   assert(digraph->nsuccessors[node] >= 0);
   assert(digraph->nsuccessors[node] <= digraph->successorssize[node]);
   assert(digraph->arcdatas != NULL);

   return digraph->arcdatas[node];
}

/** performs depth-first-search in the given directed graph from the given start node */
static
void depthFirstSearch(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< node to start the depth-first-search */
   SCIP_Bool*            visited,            /**< array to store for each node, whether it was already visited */
   int*                  dfsstack,           /**< array of size number of nodes to store the stack;
                                              *   only needed for performance reasons */
   int*                  stackadjvisited,    /**< array of size number of nodes to store the number of adjacent nodes already visited
                                              *   for each node on the stack; only needed for performance reasons */
   int*                  dfsnodes,           /**< array of nodes that can be reached starting at startnode, in reverse dfs order */
   int*                  ndfsnodes           /**< pointer to store number of nodes that can be reached starting at startnode */
   )
{
   int stacksize;
   int currnode;

   assert(digraph != NULL);
   assert(startnode >= 0);
   assert(startnode < digraph->nnodes);
   assert(visited != NULL);
   assert(visited[startnode] == FALSE);
   assert(dfsstack != NULL);
   assert(dfsnodes != NULL);
   assert(ndfsnodes != NULL);

   /* put start node on the stack */
   dfsstack[0] = startnode;
   stackadjvisited[0] = 0;
   stacksize = 1;

   while( stacksize > 0 )
   {
      /* get next node from stack */
      currnode = dfsstack[stacksize - 1];

      /* mark current node as visited */
      assert(visited[currnode] == (stackadjvisited[stacksize - 1] > 0));
      visited[currnode] = TRUE;

      /* iterate through the successor list until we reach unhandled node */
      while( stackadjvisited[stacksize - 1] < digraph->nsuccessors[currnode]
         && visited[digraph->successors[currnode][stackadjvisited[stacksize - 1]]] )
      {
         stackadjvisited[stacksize - 1]++;
      }

      /* the current node was completely handled, remove it from stack */
      if( stackadjvisited[stacksize - 1] == digraph->nsuccessors[currnode] )
      {
         stacksize--;

         /* store node in the sorted nodes array */
         dfsnodes[(*ndfsnodes)] = currnode;
         (*ndfsnodes)++;
      }
      /* handle next unhandled successor node */
      else
      {
         assert(!visited[digraph->successors[currnode][stackadjvisited[stacksize - 1]]]);

         /* put the successor node onto the stack */
         dfsstack[stacksize] = digraph->successors[currnode][stackadjvisited[stacksize - 1]];
         stackadjvisited[stacksize] = 0;
         stackadjvisited[stacksize - 1]++;
         stacksize++;
         assert(stacksize <= digraph->nnodes);
      }
   }
}

/** Compute undirected connected components on the given graph.
 *
 *  @note For each arc, its reverse is added, so the graph does not need to be the directed representation of an
 *        undirected graph.
 */
SCIP_RETCODE SCIPdigraphComputeUndirectedComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   minsize,            /**< all components with less nodes are ignored */
   int*                  components,         /**< array with as many slots as there are nodes in the directed graph
                                              *   to store for each node the component to which it belongs
                                              *   (components are numbered 0 to ncomponents - 1); or NULL, if components
                                              *   are accessed one-by-one using SCIPdigraphGetComponent() */
   int*                  ncomponents         /**< pointer to store the number of components; or NULL, if the
                                              *   number of components is accessed by SCIPdigraphGetNComponents() */
   )
{
   SCIP_Bool* visited;
   int* ndirectedsuccessors;
   int* stackadjvisited;
   int* dfsstack;
   int ndfsnodes;
   int compstart;
   int v;
   int i;
   int j;

   assert(digraph != NULL);
   assert(digraph->nnodes > 0);

   digraph->ncomponents = 0;
   digraph->componentstartsize = 10;

   SCIP_ALLOC( BMSallocClearMemoryArray(&visited, digraph->nnodes) );
   SCIP_ALLOC( BMSallocMemoryArray(&digraph->components, digraph->nnodes) );
   SCIP_ALLOC( BMSallocMemoryArray(&digraph->componentstarts, digraph->componentstartsize) );
   SCIP_ALLOC( BMSallocMemoryArray(&dfsstack, digraph->nnodes) );
   SCIP_ALLOC( BMSallocMemoryArray(&stackadjvisited, digraph->nnodes) );
   SCIP_ALLOC( BMSallocMemoryArray(&ndirectedsuccessors, digraph->nnodes) );

   digraph->componentstarts[0] = 0;

   /* store the number of directed arcs per node */
   BMScopyMemoryArray(ndirectedsuccessors, digraph->nsuccessors, digraph->nnodes);

   /* add reverse arcs to the graph */
   for( i = digraph->nnodes - 1; i >= 0; --i )
   {
      for( j = 0; j < ndirectedsuccessors[i]; ++j )
      {
         SCIP_CALL( SCIPdigraphAddArc(digraph, digraph->successors[i][j], i, NULL) );
      }
   }

   for( v = 0; v < digraph->nnodes; ++v )
   {
      if( visited[v] )
         continue;

      compstart = digraph->componentstarts[digraph->ncomponents];
      ndfsnodes = 0;
      depthFirstSearch(digraph, v, visited, dfsstack, stackadjvisited,
         &digraph->components[compstart], &ndfsnodes);

      /* forget about this component if it is too small */
      if( ndfsnodes >= minsize )
      {
         digraph->ncomponents++;

         /* enlarge componentstartsize array, if needed */
         if( digraph->ncomponents >= digraph->componentstartsize )
         {
            digraph->componentstartsize = 2 * digraph->componentstartsize;
            assert(digraph->ncomponents < digraph->componentstartsize);

            SCIP_ALLOC( BMSreallocMemoryArray(&digraph->componentstarts, digraph->componentstartsize) );
         }
         digraph->componentstarts[digraph->ncomponents] = compstart + ndfsnodes;

         /* store component number for contained nodes if array was given */
         if( components != NULL )
         {
            for( i = digraph->componentstarts[digraph->ncomponents] - 1; i >= compstart; --i )
            {
               components[digraph->components[i]] = digraph->ncomponents - 1;
            }
         }
      }
   }

   /* restore the number of directed arcs per node */
   BMScopyMemoryArray(digraph->nsuccessors, ndirectedsuccessors, digraph->nnodes);
   BMSclearMemoryArray(visited, digraph->nnodes);

   /* return number of components, if the pointer was given */
   if( ncomponents != NULL )
      (*ncomponents) = digraph->ncomponents;

   BMSfreeMemoryArray(&ndirectedsuccessors);
   BMSfreeMemoryArray(&stackadjvisited);
   BMSfreeMemoryArray(&dfsstack);
   BMSfreeMemoryArray(&visited);

   return SCIP_OKAY;
}

/** Performes an (almost) topological sort on the undirected components of the given directed graph. The undirected
 *  components should be computed before using SCIPdigraphComputeUndirectedComponents().
 *
 *  @note In general a topological sort is not unique.  Note, that there might be directed cycles, that are randomly
 *        broken, which is the reason for having only almost topologically sorted arrays.
 */
SCIP_RETCODE SCIPdigraphTopoSortComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   SCIP_Bool* visited;
   int* comps;
   int* compstarts;
   int* stackadjvisited;
   int* dfsstack;
   int* dfsnodes;
   int ndfsnodes;
   int ncomps;
   int i;
   int j;
   int k;
   int endidx;

   assert(digraph != NULL);

   ncomps = digraph->ncomponents;
   comps = digraph->components;
   compstarts = digraph->componentstarts;

   SCIP_ALLOC( BMSallocClearMemoryArray(&visited, digraph->nnodes) );
   SCIP_ALLOC( BMSallocMemoryArray(&dfsnodes, digraph->nnodes) );
   SCIP_ALLOC( BMSallocMemoryArray(&dfsstack, digraph->nnodes) );
   SCIP_ALLOC( BMSallocMemoryArray(&stackadjvisited, digraph->nnodes) );

   /* sort the components (almost) topologically */
   for( i = 0; i < ncomps; ++i )
   {
      endidx = compstarts[i+1] - 1;
      ndfsnodes = 0;
      for( j = compstarts[i]; j < compstarts[i+1]; ++j )
      {
         if( visited[comps[j]] )
            continue;

         /* perform depth first search, nodes visited in this call are appended to the list dfsnodes in reverse
          * dfs order, after the nodes already contained;
          * so at every point in time, the nodes in dfsnode are in reverse (almost) topological order
          */
         depthFirstSearch(digraph, comps[j], visited, dfsstack, stackadjvisited, dfsnodes, &ndfsnodes);
      }
      assert(endidx - ndfsnodes == compstarts[i] - 1);

      /* copy reverse (almost) topologically sorted array of nodes reached by the dfs searches;
       * reverse their order to get an (almost) topologically sort
       */
      for( k = 0; k < ndfsnodes; ++k )
      {
         digraph->components[endidx - k] = dfsnodes[k];
      }
   }

   BMSfreeMemoryArray(&stackadjvisited);
   BMSfreeMemoryArray(&dfsstack);
   BMSfreeMemoryArray(&dfsnodes);
   BMSfreeMemoryArray(&visited);

   return SCIP_OKAY;
}

/** returns the number of previously computed undirected components for the given directed graph */
int SCIPdigraphGetNComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   assert(digraph != NULL);
   assert(digraph->componentstartsize > 0); /* components should have been computed */

   return digraph->ncomponents;
}

/** Returns the previously computed undirected component of the given number for the given directed graph.
 *  If the components were sorted using SCIPdigraphTopoSortComponents(), the component is (almost) topologically sorted.
 */
void SCIPdigraphGetComponent(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   compidx,            /**< number of the component to return */
   int**                 nodes,              /**< pointer to store the nodes in the component; or NULL, if not needed */
   int*                  nnodes              /**< pointer to store the number of nodes in the component;
                                              *   or NULL, if not needed */
   )
{
   assert(digraph != NULL);
   assert(compidx >= 0);
   assert(compidx < digraph->ncomponents);
   assert(nodes != NULL || nnodes != NULL);

   if( nodes != NULL )
      (*nodes) = &(digraph->components[digraph->componentstarts[compidx]]);
   if( nnodes != NULL )
      (*nnodes) = digraph->componentstarts[compidx + 1] - digraph->componentstarts[compidx];
}

/** frees the component information for the given directed graph */
void SCIPdigraphFreeComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   )
{
   assert(digraph != NULL);

   /* free components structure */
   if( digraph->componentstartsize > 0 )
   {
      BMSfreeMemoryArray(&digraph->componentstarts);
      BMSfreeMemoryArray(&digraph->components);
      digraph->components = NULL;
      digraph->componentstarts = NULL;
      digraph->ncomponents = 0;
      digraph->componentstartsize = 0;
   }
#ifndef NDEBUG
   else
   {
      assert(digraph->components == NULL);
      assert(digraph->componentstarts == NULL);
      assert(digraph->ncomponents == 0);
   }
#endif
}

/** output of the given directed graph via the given message handler */
void SCIPdigraphPrint(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int n;

   for( n = 0; n < digraph->nnodes; ++n )
   {
      int* successors;
      int nsuccessors;
      int m;

      nsuccessors = digraph->nsuccessors[n];
      successors = digraph->successors[n];

      SCIPmessageFPrintInfo(messagehdlr, file, "node %d --> ", n);

      for( m = 0; m < nsuccessors ; ++m )
      {
         if( m == 0 )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%d", successors[m]);
         }
         else
         {
            SCIPmessageFPrintInfo(messagehdlr, file, ", %d", successors[m]);
         }
      }
      SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }
}

/** prints the given directed graph structure in GML format into the given file */
void SCIPdigraphPrintGml(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   FILE*                 file                /**< file to write to */
   )
{
   int n;

   /* write GML format opening */
   SCIPgmlWriteOpening(file, TRUE);

   /* write all nodes of the graph */
   for( n = 0; n < digraph->nnodes; ++n )
   {
      char label[SCIP_MAXSTRLEN];

      (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%d", n);
      SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", NULL, NULL);
   }

   /* write all edges */
   for( n = 0; n < digraph->nnodes; ++n )
   {
      int* successors;
      int nsuccessors;
      int m;

      nsuccessors = digraph->nsuccessors[n];
      successors = digraph->successors[n];

      for( m = 0; m < nsuccessors; ++m )
      {
         SCIPgmlWriteArc(file, (unsigned int)n, (unsigned int)successors[m], NULL, NULL);
      }
   }
   /* write GML format closing */
   SCIPgmlWriteClosing(file);
}

/** output of the given directed graph via the given message handler */
void SCIPdigraphPrintComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int c;
   int i;

   for( c = 0; c < digraph->ncomponents; ++c )
   {
      int start = digraph->componentstarts[c];
      int end =  digraph->componentstarts[c+1];

      SCIPmessageFPrintInfo(messagehdlr, file, "Components %d --> ", c);

      for( i = start; i < end; ++i )
      {
         if( i == start )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%d", digraph->components[i]);
         }
         else
         {
            SCIPmessageFPrintInfo(messagehdlr, file, ", %d", digraph->components[i]);
         }
      }
      SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }
}

/*
 * Binary tree
 */

/** creates a node for a binary tree */
static
SCIP_RETCODE btnodeCreateEmpty(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node                /**< pointer to store the created node */
   )
{
   SCIP_ALLOC( BMSallocBlockMemory(tree->blkmem, node) );

   (*node)->parent = NULL;
   (*node)->left = NULL;
   (*node)->right = NULL;
   (*node)->dataptr = NULL;

   return SCIP_OKAY;
}

/** creates a tree node with (optinal) user data */
SCIP_RETCODE SCIPbtnodeCreate(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node,               /**< pointer to store the created node */
   void*                 dataptr             /**< user node data pointer, or NULL */
   )
{
   assert(tree != NULL);
   assert(node != NULL);

   SCIP_CALL( btnodeCreateEmpty(tree, node) );

   assert((*node)->parent == NULL);
   assert((*node)->left == NULL);
   assert((*node)->right == NULL);

   /* initialize user data */
   (*node)->dataptr = dataptr;

   return SCIP_OKAY;
}

/** frees a tree leaf */
static
void btnodeFreeLeaf(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node                /**< pointer to node which has to be freed */
   )
{
   assert(tree != NULL);
   assert(node != NULL);
   assert(*node != NULL);

   assert((*node)->left == NULL);
   assert((*node)->right == NULL);

#if 0
   /* remove reference from parent node */
   if( (*node)->parent != NULL )
   {
      assert(*node != NULL);

      assert((*node)->parent->left == *node || ((*node)->parent->right == *node));

      if( (*node)->parent->left == *node )
      {
         (*node)->parent->left = NULL;
      }
      else
      {
         assert((*node)->parent->right == *node);
         (*node)->parent->right = NULL;
      }
   }
#endif

   assert(*node != NULL);
   BMSfreeBlockMemory(tree->blkmem, node);
   assert(*node == NULL);
}

/** frees the node including the rooted subtree
 *
 *  @note The user pointer (object) is not freed. If needed, it has to be done by the user.
 */
void SCIPbtnodeFree(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node                /**< node to be freed */
   )
{
   assert(tree != NULL);
   assert(node != NULL);
   assert(*node != NULL);

   if( (*node)->left != NULL )
   {
      SCIPbtnodeFree(tree, &(*node)->left);
      assert((*node)->left == NULL);
   }

   if( (*node)->right != NULL )
   {
      SCIPbtnodeFree(tree, &(*node)->right);
      assert((*node)->right == NULL);
   }

   btnodeFreeLeaf(tree, node);
   assert(*node == NULL);
}

/* some simple variable functions implemented as defines */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPbtnodeGetData
#undef SCIPbtnodeGetKey
#undef SCIPbtnodeGetParent
#undef SCIPbtnodeGetLeftchild
#undef SCIPbtnodeGetRightchild
#undef SCIPbtnodeGetSibling
#undef SCIPbtnodeIsRoot
#undef SCIPbtnodeIsLeaf
#undef SCIPbtnodeIsLeftchild
#undef SCIPbtnodeIsRightchild

/** returns the user data pointer stored in that node */
void* SCIPbtnodeGetData(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return node->dataptr;
}

/** returns the parent which can be NULL if the given node is the root */
SCIP_BTNODE* SCIPbtnodeGetParent(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return node->parent;
}

/** returns left child which can be NULL if the given node is a leaf */
SCIP_BTNODE* SCIPbtnodeGetLeftchild(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return node->left;
}

/** returns right child which can be NULL if the given node is a leaf */
SCIP_BTNODE* SCIPbtnodeGetRightchild(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return node->right;
}

/** returns the sibling of the node or NULL if does not exist */
SCIP_BTNODE* SCIPbtnodeGetSibling(
   SCIP_BTNODE*          node                /**< node */
   )
{
   SCIP_BTNODE* parent;

   parent = SCIPbtnodeGetParent(node);

   if( parent == NULL )
      return NULL;

   if( SCIPbtnodeGetLeftchild(parent) == node )
      return SCIPbtnodeGetRightchild(parent);

   assert(SCIPbtnodeGetRightchild(parent) == node);

   return SCIPbtnodeGetLeftchild(parent);
}

/** returns whether the node is a root node */
SCIP_Bool SCIPbtnodeIsRoot(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return (node->parent == NULL);
}

/** returns whether the node is a leaf */
SCIP_Bool SCIPbtnodeIsLeaf(
   SCIP_BTNODE*          node                /**< node */
   )
{
   assert(node != NULL);

   return (node->left == NULL && node->right == NULL);
}

/** returns TRUE if the given node is left child */
SCIP_Bool SCIPbtnodeIsLeftchild(
   SCIP_BTNODE*          node                /**< node */
   )
{
   SCIP_BTNODE* parent;

   if( SCIPbtnodeIsRoot(node) )
      return FALSE;

   parent = SCIPbtnodeGetParent(node);

   if( SCIPbtnodeGetLeftchild(parent) == node )
      return TRUE;

   return FALSE;
}

/** returns TRUE if the given node is right child */
SCIP_Bool SCIPbtnodeIsRightchild(
   SCIP_BTNODE*          node                /**< node */
   )
{
   SCIP_BTNODE* parent;

   if( SCIPbtnodeIsRoot(node) )
      return FALSE;

   parent = SCIPbtnodeGetParent(node);

   if( SCIPbtnodeGetRightchild(parent) == node )
      return TRUE;

   return FALSE;
}

/** sets the give node data
 *
 *  @note The old user pointer is not freed.
 */
void SCIPbtnodeSetData(
   SCIP_BTNODE*          node,               /**< node */
   void*                 dataptr             /**< node user data pointer */
   )
{
   assert(node != NULL);

   node->dataptr = dataptr;
}

/** sets parent node
 *
 *  @note The old parent including the rooted subtree is not delete.
 */
void SCIPbtnodeSetParent(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          parent              /**< new parent node, or NULL */
   )
{
   assert(node != NULL);

   node->parent = parent;
}

/** sets left child
 *
 *  @note The old left child including the rooted subtree is not delete.
 */
void SCIPbtnodeSetLeftchild(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          left                /**< new left child, or NULL */
   )
{
   assert(node != NULL);

   node->left = left;
}

/** sets right child
 *
 *  @note The old right child including the rooted subtree is not delete.
 */
void SCIPbtnodeSetRightchild(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          right               /**< new right child, or NULL */
   )
{
   assert(node != NULL);

   node->right = right;
}

/** creates an binary tree */
SCIP_RETCODE SCIPbtCreate(
   SCIP_BT**             tree,               /**< pointer to store the created binary tree */
   BMS_BLKMEM*           blkmem              /**< block memory used to createnode */
   )
{
   assert(tree != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocMemory(tree) );
   (*tree)->blkmem = blkmem;
   (*tree)->root = NULL;

   return SCIP_OKAY;
}

/** frees binary tree
 *
 *  @note The user pointers (object) of the nodes are not freed. If needed, it has to be done by the user.
 */
void SCIPbtFree(
   SCIP_BT**             tree                /**< pointer to binary tree */
   )
{
   assert(tree != NULL);

   if( (*tree)->root != NULL )
   {
      SCIPbtnodeFree(*tree, &((*tree)->root));
   }

   BMSfreeMemory(tree);
}

/** prints the rooted subtree of the given binary tree node in GML format into the given file */
static
void btPrintSubtree(
   SCIP_BTNODE*          node,               /**< binary tree node */
   FILE*                 file,               /**< file to write to */
   int*                  nnodes              /**< pointer to count the number of nodes */
   )
{
   SCIP_BTNODE* left;
   SCIP_BTNODE* right;
   char label[SCIP_MAXSTRLEN];

   assert(node != NULL);

   (*nnodes)++;
   (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%d", *nnodes);

   SCIPgmlWriteNode(file, (unsigned int)(size_t)node, label, "circle", NULL, NULL);

   left = SCIPbtnodeGetLeftchild(node);
   right = SCIPbtnodeGetRightchild(node);

   if( left != NULL )
   {
      btPrintSubtree(left, file, nnodes);

      SCIPgmlWriteArc(file, (unsigned int)(size_t)node, (unsigned int)(size_t)left, NULL, NULL);
   }

   if( right != NULL )
   {
      btPrintSubtree(right, file, nnodes);

      SCIPgmlWriteArc(file, (unsigned int)(size_t)node, (unsigned int)(size_t)right, NULL, NULL);
   }
}

/** prints the binary tree in GML format into the given file */
void SCIPbtPrintGml(
   SCIP_BT*              tree,               /**< binary tree */
   FILE*                 file                /**< file to write to */
   )
{
   /* write GML opening */
   SCIPgmlWriteOpening(file, TRUE);

   if( !SCIPbtIsEmpty(tree) )
   {
      SCIP_BTNODE* root;
      int nnodes;

      root = SCIPbtGetRoot(tree);
      assert(root != NULL);

      nnodes = 0;

      btPrintSubtree(root, file, &nnodes);
   }

   /* write GML closing */
   SCIPgmlWriteClosing(file);
}

/* some simple variable functions implemented as defines */
#undef SCIPbtIsEmpty
#undef SCIPbtGetRoot

/** returns whether the binary tree is empty (has no nodes) */
SCIP_Bool SCIPbtIsEmpty(
   SCIP_BT*              tree                /**< binary tree */
   )
{
   assert(tree != NULL);

   return (tree->root == NULL);
}

/** returns the the root node of the binary or NULL if the binary tree is empty */
SCIP_BTNODE* SCIPbtGetRoot(
   SCIP_BT*              tree                /**< tree to be evaluated */
   )
{
   assert(tree != NULL);

   return tree->root;
}

/** sets root node
 *
 *  @note The old root including the rooted subtree is not delete.
 */
void SCIPbtSetRoot(
   SCIP_BT*              tree,               /**< tree to be evaluated */
   SCIP_BTNODE*          root                /**< new root, or NULL */
   )
{
   assert(tree != NULL);

   tree->root = root;
}


/*
 * Numerical methods
 */

/** returns the machine epsilon: the smallest number eps > 0, for which 1.0 + eps > 1.0 */
SCIP_Real SCIPcalcMachineEpsilon(
   void
   )
{
   SCIP_Real eps;
   SCIP_Real lasteps;
   SCIP_Real one;
   SCIP_Real onepluseps;

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
SCIP_Longint SCIPcalcGreComDiv(
   SCIP_Longint          val1,               /**< first value of greatest common devisor calculation */
   SCIP_Longint          val2                /**< second value of greatest common devisor calculation */
   )
{
   int t;

   assert(val1 > 0);
   assert(val2 > 0);
   
   t = 0;
   /* if val1 is even, divide it by 2 */
   while( !(val1 & 1) )
   {
      val1 >>= 1; /*lint !e704*/
      
      /* if val2 is even too, divide it by 2 and increase t(=number of e) */
      if( !(val2 & 1) )
      {
         val2 >>= 1; /*lint !e704*/
         ++t;
      }
      /* only val1 can be odd */
      else
      {
         /* while val1 is even, divide it by 2 */
         while( !(val1 & 1) )
            val1 >>= 1; /*lint !e704*/
         
         break;
      }
   }
   /* while val2 is even, divide it by 2 */
   while( !(val2 & 1) )
      val2 >>= 1; /*lint !e704*/
   
   /* val1 and val 2 are odd */
   while( val1 != val2 )
   {
      if( val1 > val2 )
      {
         val1 -= val2;
         /* val1 is now even, divide it by 2  */
         do 
         {
            val1 >>= 1;   /*lint !e704*/
         }
         while( !(val1 & 1) );
      }
      else 
      {
         val2 -= val1;
         /* val2 is now even, divide it by 2  */
         do 
         {
            val2 >>= 1;  /*lint !e704*/
         }
         while( !(val2 & 1) );
      }
   }

   return (val1 << t);  /*lint !e703*/
}

/** calculates the smallest common multiple of the two given values */
SCIP_Longint SCIPcalcSmaComMul(
   SCIP_Longint          val1,               /**< first value of smallest common multiple calculation */
   SCIP_Longint          val2                /**< second value of smallest common multiple calculation */
   )
{
   SCIP_Longint gcd;

   assert(val1 > 0);
   assert(val2 > 0);

   gcd = SCIPcalcGreComDiv(val1, val2);
   
   return val1/gcd * val2;
}

static const SCIP_Real simplednoms[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
                                   17.0, 18.0, 19.0, 25.0, -1.0};

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
SCIP_Bool SCIPrealToRational(
   SCIP_Real             val,                /**< real value r to convert into rational number */
   SCIP_Real             mindelta,           /**< minimal allowed difference r - q of real r and rational q = n/d */
   SCIP_Real             maxdelta,           /**< maximal allowed difference r - q of real r and rational q = n/d */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   )
{
   SCIP_Real a;
   SCIP_Real b;
   SCIP_Real g0;
   SCIP_Real g1;
   SCIP_Real gx;
   SCIP_Real h0;
   SCIP_Real h1;
   SCIP_Real hx;
   SCIP_Real delta0;
   SCIP_Real delta1;
   SCIP_Real epsilon;
   int i;

   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(nominator != NULL);
   assert(denominator != NULL);

   /* try the simple denominators first: each value of the simpledenoms table multiplied by powers of 10
    * is tried as denominator
    */
   for( i = 0; simplednoms[i] > 0.0; ++i )
   {
      SCIP_Real nom;
      SCIP_Real dnom;
      SCIP_Real ratval0;
      SCIP_Real ratval1;

      /* try powers of 10 (including 10^0) */
      dnom = simplednoms[i];
      while( dnom <= maxdnom )
      {
         nom = floor(val * dnom);
         ratval0 = nom/dnom;
         ratval1 = (nom+1.0)/dnom;
         if( mindelta <= val - ratval0 && val - ratval1 <= maxdelta )
         {
            if( val - ratval0 <= maxdelta )
            {
               *nominator = (SCIP_Longint)nom;
               *denominator = (SCIP_Longint)dnom;
               return TRUE;
            }
            if( mindelta <= val - ratval1 )
            {
               *nominator = (SCIP_Longint)(nom+1.0);
               *denominator = (SCIP_Longint)dnom;
               return TRUE;
            }
         }
         dnom *= 10.0;
      }
   }

   /* the simple denominators didn't work: calculate rational representation with arbitrary denominator */
   epsilon = MIN(-mindelta, maxdelta)/2.0;

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

   if( REALABS(g0) > (SCIP_LONGINT_MAX >> 4) || h0 > (SCIP_LONGINT_MAX >> 4) )
      return FALSE;

   assert(h0 > 0.5);

   if( delta0 < mindelta )
   {
      assert(mindelta <= delta1 && delta1 <= maxdelta);
      *nominator = (SCIP_Longint)(g0 - 1.0);
      *denominator = (SCIP_Longint)h0;
   }
   else if( delta0 > maxdelta )
   {
      assert(mindelta <= delta1 && delta1 <= maxdelta);
      *nominator = (SCIP_Longint)(g0 + 1.0);
      *denominator = (SCIP_Longint)h0;
   }
   else
   {
      *nominator = (SCIP_Longint)g0;
      *denominator = (SCIP_Longint)h0;
   }
   assert(*denominator >= 1);
   assert(val - (SCIP_Real)(*nominator)/(SCIP_Real)(*denominator) >= mindelta);
   assert(val - (SCIP_Real)(*nominator)/(SCIP_Real)(*denominator) <= maxdelta);

   return TRUE;
}

/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */
static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = floor(sval);
   upval = ceil(sval);

   return (SCIPrelDiff(sval, downval) <= maxdelta || SCIPrelDiff(sval, upval) >= mindelta);
}

/** additional scalars that are tried in integrality scaling */
static const SCIP_Real scalars[] = {3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};
static const int nscalars = 9;

/** tries to find a value, such that all given values, if scaled with this value become integral */
SCIP_RETCODE SCIPcalcIntegralScalar(
   SCIP_Real*            vals,               /**< values to scale */
   int                   nvals,              /**< number of values to scale */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   )
{
   SCIP_Real bestscalar;
   SCIP_Longint gcd;
   SCIP_Longint scm;
   SCIP_Longint nominator;
   SCIP_Longint denominator;
   SCIP_Real val;
   SCIP_Real minval;
   SCIP_Real absval;
   SCIP_Real scaleval;
   SCIP_Bool scalable;
   SCIP_Bool rational;
   int c;
   int s;
   int i;

   assert(vals != NULL);
   assert(nvals >= 0);
   assert(maxdnom >= 1);
   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(success != NULL);

   SCIPdebugMessage("trying to find rational representation for given values\n");

   if( intscalar != NULL )
      *intscalar = SCIP_INVALID;
   *success = FALSE;

   /* get minimal absolute non-zero value */
   minval = SCIP_REAL_MAX;
   for( c = 0; c < nvals; ++c )
   {
      val = vals[c];
      if( val < mindelta || val > maxdelta )
      {
         absval = REALABS(val);
         minval = MIN(minval, absval);
      }
   }

   if( minval == SCIP_REAL_MAX )
   {
      /* all coefficients are zero (inside tolerances) */
      if( intscalar != NULL )
         *intscalar = 1.0;
      *success = TRUE;
      SCIPdebugMessage(" -> all values are zero (inside tolerances)\n");

      return SCIP_OKAY;
   }
   assert(minval > MIN(-mindelta, maxdelta));

   bestscalar = SCIP_INVALID;

   for( i = 0; i < 2; ++i )
   {
      scalable = TRUE;

      /* try, if values can be made integral multiplying them with the reciprocal of the smallest value and a power of 2 */
      if( i == 0 )
	 scaleval = 1.0/minval;
      /* try, if values can be made integral by multiplying them by a power of 2 */
      else
	 scaleval = 1.0;

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
	 SCIPdebugMessage(" -> val=%g, scaleval=%g, val*scaleval=%g, scalable=%u\n",
	    val, scaleval, val*scaleval, scalable);
      }
      if( scalable )
      {
	 /* make values integral by dividing them by the smallest value (and multiplying them with a power of 2) */
	 assert(scaleval <= maxscale);

	 /* check if we found a better scaling value */
	 if( scaleval < bestscalar )
	    bestscalar = scaleval;

	 SCIPdebugMessage(" -> integrality could be achieved by scaling with %g\n", scaleval);

	 /* if the scalar is still the reciprocal of the minimal value, all coeffcients are the same and we do not get a better scalar */
	 if( i == 0 && EPSEQ(scaleval, 1.0/minval, SCIP_DEFAULT_EPSILON) )
	 {
	    if( intscalar != NULL )
	       *intscalar = bestscalar;
	    *success = TRUE;

	    return SCIP_OKAY;
	 }
      }
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
         rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
         SCIPdebugMessage(" -> c=%d first rational: val: %g == %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", gcd=%"SCIP_LONGINT_FORMAT", scm=%"SCIP_LONGINT_FORMAT", rational=%u\n",
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
         rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
         SCIPdebugMessage(" -> c=%d next rational : val: %g == %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", gcd=%"SCIP_LONGINT_FORMAT", scm=%"SCIP_LONGINT_FORMAT", rational=%u\n",
            c, val, nominator, denominator, gcd, scm, rational);
      }
      else
      {
         SCIPdebugMessage(" -> failed to convert %g into a rational representation\n", val);
      }
   }

   if( rational )
   {
      /* make values integral by multiplying them with the smallest common multiple of the denominators */
      assert((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);

      /* check if we found a better scaling value */
      if( (SCIP_Real)scm/(SCIP_Real)gcd < bestscalar )
	 bestscalar = (SCIP_Real)scm/(SCIP_Real)gcd;

      SCIPdebugMessage(" -> integrality could be achieved by scaling with %g (rational:%"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT")\n",
         (SCIP_Real)scm/(SCIP_Real)gcd, scm, gcd);
   }

   if( bestscalar < SCIP_INVALID )
   {
      if( intscalar != NULL )
         *intscalar = bestscalar;
      *success = TRUE;

      SCIPdebugMessage(" -> smallest value to achieve integrality is %g \n", bestscalar);
   }

   return SCIP_OKAY;
}

/** given a (usually very small) interval, tries to find a rational number with simple denominator (i.e. a small
 *  number, probably multiplied with powers of 10) out of this interval; returns TRUE iff a valid rational
 *  number inside the interval was found
 */
SCIP_Bool SCIPfindSimpleRational(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed for resulting rational number */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   )
{
   SCIP_Real center;
   SCIP_Real delta;

   assert(lb <= ub);

   center = 0.5*(lb+ub);

   /* in order to compute a rational number that is exactly within the bounds (as the user expects),
    * we computed the allowed delta with downward rounding, if available
    */
   if( SCIPintervalHasRoundingControl() )
   {
      SCIP_ROUNDMODE roundmode;

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      delta = 0.5*(ub-lb);

      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {
      delta = 0.5*(ub-lb);
   }

   return SCIPrealToRational(center, -delta, +delta, maxdnom, nominator, denominator);
}

/** given a (usually very small) interval, selects a value inside this interval; it is tried to select a rational number
 *  with simple denominator (i.e. a small number, probably multiplied with powers of 10);
 *  if no valid rational number inside the interval was found, selects the central value of the interval
 */
SCIP_Real SCIPselectSimpleValue(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom             /**< maximal denominator allowed for resulting rational number */
   )
{
   SCIP_Real val;

   val = 0.5*(lb+ub);
   if( lb < ub )
   {
      SCIP_Longint nominator;
      SCIP_Longint denominator;
      SCIP_Bool success;
      
      /* try to find a "simple" rational number inside the interval */
      SCIPdebugMessage("simple rational in [%.9f,%.9f]:", lb, ub);
      success = SCIPfindSimpleRational(lb, ub, maxdnom, &nominator, &denominator);
      if( success )
      {
         val = (SCIP_Real)nominator/(SCIP_Real)denominator;
         SCIPdebugPrintf(" %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT" == %.9f\n", nominator, denominator, val);

         if( val - lb < 0.0 || val - ub > 0.0 )
         {
            SCIPdebugPrintf(" value is out of interval bounds by %g -> failed\n", MAX(lb-val, val-ub));
            val = 0.5*(lb+ub);
         }
      }
      else
      {
         SCIPdebugPrintf(" failed\n");
      }
   }
   
   return val;
}




/*
 * Random Numbers
 */

#ifdef NO_RAND_R

#define SCIP_RAND_MAX 32767
/** returns a random number between 0 and SCIP_RAND_MAX */
static
int getRand(
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   SCIP_Longint nextseed;

   assert(seedp != NULL);

   nextseed = (*seedp) * 1103515245 + 12345;
   *seedp = (unsigned int)nextseed;

   return (int)((unsigned int)(nextseed/(2*(SCIP_RAND_MAX+1))) % (SCIP_RAND_MAX+1));
}

#else

#define SCIP_RAND_MAX RAND_MAX

/** returns a random number between 0 and SCIP_RAND_MAX */
static
int getRand(
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return rand_r(seedp);
}

#endif

/** returns a random integer between minrandval and maxrandval */
int SCIPgetRandomInt(
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return minrandval + (int) ((maxrandval - minrandval + 1)*(SCIP_Real)getRand(seedp)/(SCIP_RAND_MAX+1.0));
}

/** returns a random real between minrandval and maxrandval */
SCIP_Real SCIPgetRandomReal(
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return minrandval + (maxrandval - minrandval)*(SCIP_Real)getRand(seedp)/(SCIP_Real)SCIP_RAND_MAX;
}


/*
 * Additional math functions
 */

/** calculates a binomial coefficient n over m, choose m elements out of n, maximal value will be 33 over 16 (because
 *  the n=33 is the last line in the Pascal's triangle where each entry fits in a 4 byte value), an error occurs due to
 *  big numbers or an negative value m (and m < n) and -1 will be returned
 */
SCIP_Longint SCIPcalcBinomCoef(
   int                   n,                  /**< number of different elements */
   int                   m                   /**< number to choose out of the above */
   )
{
   if( m == 0 || m >= n )
      return 1;

   if( m < 0 )
      return -1;

   /* symmetry of the binomial coefficient, choose smaller m */
   if( m > n/2 )
      m = n - m;

   /* trivial case m == 1 */
   if( m == 1 )
      return n;

   /* simple case m == 2 */
   if( m == 2 )
   {
      if( ((SCIP_Real)SCIP_LONGINT_MAX) / n >= (n-1) * 2 ) /*lint !e790*/
	 return (n*(n-1)/2); /*lint !e647*/
      else
	 return -1;
   }

   /* abort on to big numbers */
   if( m > 16 || n > 33 )
      return -1;

   /* simple case m == 3 */
   if( m == 3 )
      return (n*(n-1)*(n-2)/6); /*lint !e647*/
   else
   {
      /* first half of Pascal's triangle numbers(without the symmetric part) backwards from (33,16) over (32,16),
       * (33,15), (32,15),(31,15, (30,15), (33,14) to (8,4) (rest is calculated directly)
       *
       * due to this order we can extract the right binomial coefficient by (16-m)^2+(16-m)+(33-n)
       */
      static const SCIP_Longint binoms[182] = {
         1166803110, 601080390, 1037158320, 565722720, 300540195, 155117520, 818809200, 471435600, 265182525, 145422675,
         77558760, 40116600, 573166440, 347373600, 206253075, 119759850, 67863915, 37442160, 20058300, 10400600,
         354817320, 225792840, 141120525, 86493225, 51895935, 30421755, 17383860, 9657700, 5200300, 2704156, 193536720,
         129024480, 84672315, 54627300, 34597290, 21474180, 13037895, 7726160, 4457400, 2496144, 1352078, 705432,
         92561040, 64512240, 44352165, 30045015, 20030010, 13123110, 8436285, 5311735, 3268760, 1961256, 1144066,
         646646, 352716, 184756, 38567100, 28048800, 20160075, 14307150, 10015005, 6906900, 4686825, 3124550, 2042975,
         1307504, 817190, 497420, 293930, 167960, 92378, 48620, 13884156, 10518300, 7888725, 5852925, 4292145, 3108105,
         2220075, 1562275, 1081575, 735471, 490314, 319770, 203490, 125970, 75582, 43758, 24310, 12870, 4272048, 3365856,
         2629575, 2035800, 1560780, 1184040, 888030, 657800, 480700, 346104, 245157, 170544, 116280, 77520, 50388, 31824,
         19448, 11440, 6435, 3432, 1107568, 906192, 736281, 593775, 475020, 376740, 296010, 230230, 177100, 134596,
         100947, 74613, 54264, 38760, 27132, 18564, 12376, 8008, 5005, 3003, 1716, 924, 237336, 201376, 169911, 142506,
         118755, 98280, 80730, 65780, 53130, 42504, 33649, 26334, 20349, 15504, 11628, 8568, 6188, 4368, 3003, 2002,
         1287, 792, 462, 252, 40920, 35960, 31465, 27405, 23751, 20475, 17550, 14950, 12650, 10626, 8855, 7315, 5985,
         4845, 3876, 3060, 2380, 1820, 1365, 1001, 715, 495, 330, 210, 126, 70};

      /* m can at most be 16 */
      const int t = 16-m;
      assert(t >= 0);
      assert(n <= 33);

      /* binoms array hast exactly 182 elements */
      assert(t*(t+1)+(33-n) < 182);

      return binoms[t*(t+1)+(33-n)]; /*lint !e662 !e661*/
   }
}

/** negates a number */
SCIP_Real SCIPnegateReal(
   SCIP_Real             x                   /**< value to negate */
   )
{
   return -x;
}

/*
 * Permutations / Shuffling
 */

/** swaps two ints */
void SCIPswapInts(
   int*                  value1,             /**< pointer to first integer */
   int*                  value2              /**< pointer ti second integer */
   )
{
   int tmp;

   tmp = *value1;
   *value1 = *value2;
   *value2 = tmp;
}

/** swaps the addresses of two pointers */
void SCIPswapPointers(
   void**                pointer1,           /**< first pointer */
   void**                pointer2            /**< second pointer */
   )
{
   void* tmp;

   tmp = *pointer1;
   *pointer1 = *pointer2;
   *pointer2 = tmp;
}

/** randomly shuffles parts of an integer array using the Fisher-Yates algorithm */
void SCIPpermuteIntArray(
   int*                  array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end,                /**< last index that should be subject to shuffling (array size for whole
					      *   array)
					      */
   unsigned int*         randseed            /**< seed value for the random generator */
   )
{
   int tmp;
   int i;

   /* loop backwards through all elements and always swap the current last element to a random position */
   while( end > begin+1 )
   {
      --end;

      /* get a random position into which the last entry should be shuffled */
      i = SCIPgetRandomInt(begin, end, randseed);

      /* swap the last element and the random element */
      tmp = array[i];
      array[i] = array[end];
      array[end] = tmp;
   }
}


/** randomly shuffles parts of an array using the Fisher-Yates algorithm */
void SCIPpermuteArray(
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end,                /**< last index that should be subject to shuffling (array size for whole
					      *   array)
					      */
   unsigned int*         randseed            /**< seed value for the random generator */
   )
{
   void* tmp;
   int i;

   /* loop backwards through all elements and always swap the current last element to a random position */
   while( end > begin+1 )
   {
      end--;

      /* get a random position into which the last entry should be shuffled */
      i = SCIPgetRandomInt(begin, end, randseed);

      /* swap the last element and the random element */
      tmp = array[i];
      array[i] = array[end];
      array[end] = tmp;
   }
}

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 */
SCIP_RETCODE SCIPgetRandomSubset(
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems,          /**< number of elements that should be drawn and stored */
   unsigned int          randseed            /**< seed value for random generator */
   )
{
   int i;
   int j;

   /* if both sets are of equal size, we just copy the array */
   if( nelems == nsubelems)
   {
      BMScopyMemoryArray(subset,set,nelems);
      return SCIP_OKAY;
   }

   /* abort, if size of subset is too big */
   if( nsubelems > nelems )
   {
      SCIPerrorMessage("Cannot create %d-elementary subset of %d-elementary set.\n", nsubelems, nelems);
      return SCIP_INVALIDDATA;
   }
#ifndef NDEBUG
   for( i = 0; i < nsubelems; i++ )
      for( j = 0; j < i; j++ )
         assert(set[i] != set[j]);
#endif

   /* draw each element individually */
   i = 0;
   while( i < nsubelems )
   {
      int r;

      r = SCIPgetRandomInt(0, nelems-1, &randseed);
      subset[i] = set[r];

      /* if we get an element that we already had, we will draw again */
      for( j = 0; j < i; j++ ) 
      {
         if( subset[i] == subset[j] ) 
         {
            --i;
            break;
         }
      }
      ++i;
   }
   return SCIP_OKAY;
}



/*
 * Strings
 */


/** copies characters from 'src' to 'dest', copying is stopped when either the 'stop' character is reached or after
 *  'cnt' characters have been copied, whichever comes first.
 *
 *  @note undefined behaviuor on overlapping arrays
 */
int SCIPmemccpy(
   char*                 dest,               /**< destination pointer to copy to */
   const char*           src,                /**< source pointer to copy to */
   char                  stop,               /**< character when found stop copying */
   unsigned int          cnt                 /**< maximal number of characters to copy too */
   )
{
   if( dest == NULL || src == NULL || cnt == 0 )
      return -1;
   else
   {
      char* destination = dest;

      while( cnt-- && (*destination++ = *src++) != stop ); /*lint !e722*/

      return (destination - dest);
   }
}

/** prints an error message containing of the given string followed by a string describing the current system error;
 *  prefers to use the strerror_r method, which is threadsafe; on systems where this method does not exist,
 *  NO_STRERROR_R should be defined (see INSTALL), in this case, strerror is used which is not guaranteed to be
 *  threadsafe (on SUN-systems, it actually is)
 */
void SCIPprintSysError(
   const char*           message             /**< first part of the error message, e.g. the filename */
   )
{
#ifdef NO_STRERROR_R
   char* buf;
   buf = strerror(errno);
#else
   char buf[SCIP_MAXSTRLEN];

#if defined(_WIN32) || defined(_WIN64)
   (void) strerror_s(buf, SCIP_MAXSTRLEN, errno);
#else
   (void) strerror_r(errno, buf, SCIP_MAXSTRLEN);
#endif

   buf[SCIP_MAXSTRLEN - 1] = '\0';
#endif
   SCIPmessagePrintError("%s: %s\n", message, buf);
}

/** extracts tokens from strings - wrapper method for strtok_r() */
char* SCIPstrtok(
   char*                 s,                  /**< string to parse */
   const char*           delim,              /**< delimiters for parsing */
   char**                ptrptr              /**< pointer to working char pointer - must stay the same while parsing */
   )
{
#ifdef NO_STRTOK_R
   return strtok(s, delim);
#else
   return strtok_r(s, delim, ptrptr);
#endif
}

/** translates the given string into a string where symbols ", ', and spaces are escaped with a \ prefix */
void SCIPescapeString(
   char*                 t,                  /**< target buffer to store escaped string */
   int                   bufsize,            /**< size of buffer t */
   const char*           s                   /**< string to transform into escaped string */
   )
{
   int len;
   int i;
   int p;

   assert(t != NULL);
   assert(bufsize > 0);

   len = (int)strlen(s);
   for( p = 0, i = 0; i <= len && p < bufsize; ++i, ++p )
   {
      if( s[i] == ' ' || s[i] == '"' || s[i] == '\'' )
      {
         t[p] = '\\';
         p++;
      }
      if( p < bufsize )
         t[p] = s[i];
   }
   t[bufsize-1] = '\0';
}

/* safe version of snprintf */
int SCIPsnprintf(
   char*                 t,                  /**< target string */
   int                   len,                /**< length of the string to copy */
   const char*           s,                  /**< source string */
   ...                                       /**< further parameters */
   )
{
   va_list ap;
   int n;

   assert(t != NULL);
   assert(len > 0);

   va_start(ap, s); /*lint !e826*/

#if defined(_WIN32) || defined(_WIN64)
   n = _vsnprintf(t, (size_t) len, s, ap);
#else
   n = vsnprintf(t, (size_t) len, s, ap);
#endif
   va_end(ap);

   if( n < 0 || n >= len )
   {
#ifndef NDEBUG
      if( n < 0 )
      {
         SCIPerrorMessage("vsnprintf returned %d\n",n);
      }
#endif
      t[len-1] = '\0';
      n = len-1;
   }
   return n;
}

/** extract the next token as a integer value if it is one; in case no value is parsed the endptr is set to @p str
 *
 *  @return Returns TRUE if a value could be extracted, otherwise FALSE
 */
SCIP_Bool SCIPstrToIntValue(
   const char*           str,                /**< string to search */
   int*                  value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   )
{
   assert(str != NULL);
   assert(value != NULL);
   assert(endptr != NULL);

   /* init errno to detect possible errors */
   errno = 0;

   *value = strtol(str, endptr, 10);

   if( *endptr != str && *endptr != NULL )
   {
      SCIPdebugMessage("parsed integer value <%d>\n", *value);
      return TRUE;
   }
   *endptr = (char*)str;

   SCIPdebugMessage("failed parsing integer value <%s>\n", str);

   return FALSE;
}

/** extract the next token as a double value if it is one; in case no value is parsed the endptr is set to @p str
 *
 *  @return Returns TRUE if a value could be extracted, otherwise FALSE
 */
SCIP_Bool SCIPstrToRealValue(
   const char*           str,                /**< string to search */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   )
{
   assert(str != NULL);
   assert(value != NULL);
   assert(endptr != NULL);

   /* init errno to detect possible errors */
   errno = 0;

   *value = strtod(str, endptr);

   if( *endptr != str && *endptr != NULL )
   {
      SCIPdebugMessage("parsed real value <%g>\n", *value);
      return TRUE;
   }
   *endptr = (char*)str;

   SCIPdebugMessage("failed parsing real value <%s>\n", str);

   return FALSE;
}

/** copies the first size characters between a start and end character of str into token, if no error occured endptr
 *  will point to the position after the read part, otherwise it will point to @p str
 */
void SCIPstrCopySection(
   const char*           str,                /**< string to search */
   char                  startchar,          /**< character which defines the beginning */
   char                  endchar,            /**< character which defines the ending */
   char*                 token,              /**< string to store the copy */
   int                   size,               /**< size of the token char array */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   )
{
   const char* copystr;
   int nchars;

   assert(str != NULL);
   assert(token != NULL);
   assert(size > 0);
   assert(endptr != NULL);

   nchars = 0;

   copystr = str;

   /* find starting character */
   while( *str != '\0' && *str != startchar )
      ++str;

   /* did not find start character */
   if( *str == '\0' )
   {
      *endptr = (char*)copystr;
      return;
   }

   /* skip start character */
   ++str;

   /* copy string */
   while( *str != '\0' && *str != endchar && nchars < size-1 )
   {
      assert(nchars < SCIP_MAXSTRLEN);
      token[nchars] = *str;
      nchars++;
      ++str;
   }

   /* add end to token */
   token[nchars] = '\0';

   /* if section was longer than size, we want to reach the end of the parsing section anyway */
   if( nchars == (size-1) )
      while( *str != '\0' && *str != endchar )
         ++str;

   /* did not find end character */
   if( *str == '\0' )
   {
      *endptr = (char*)copystr;
      return;
   }

   /* skip end character */
   ++str;

   SCIPdebugMessage("parsed section <%s>\n", token);

   *endptr = (char*) str;
}

/*
 * File methods
 */

/** returns, whether the given file exists */
SCIP_Bool SCIPfileExists(
   const char*           filename            /**< file name */
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
   char*                 filename,           /**< filename to split; is destroyed (but not freed) during process */
   char**                path,               /**< pointer to store path, or NULL if not needed */
   char**                name,               /**< pointer to store name, or NULL if not needed */
   char**                extension,          /**< pointer to store extension, or NULL if not needed */
   char**                compression         /**< pointer to store compression extension, or NULL if not needed */
   )
{
   char* lastslash;
   char* lastbackslash;
   char* lastdot;

   assert(filename != NULL);

   if( path != NULL )
      *path = NULL;
   if( name != NULL )
      *name = NULL;
   if( extension != NULL )
      *extension = NULL;
   if( compression != NULL )
      *compression = NULL;

   /* treat both slashes '/' and '\' as directory delimiters */
   lastslash = strrchr(filename, '/');
   lastbackslash = strrchr(filename, '\\');
   lastslash = MAX(lastslash, lastbackslash); /*lint !e613*/
   lastdot = strrchr(filename, '.');
   if( lastslash != NULL && lastdot != NULL && lastdot < lastslash ) /* is the last dot belonging to the path? */
      lastdot = NULL;

   /* detect known compression extensions */
#ifdef WITH_ZLIB
   if( lastdot != NULL )
   {
      char* compext;

      compext = lastdot+1;
      if( strcmp(compext, "gz") == 0
        || strcmp(compext, "z") == 0
        || strcmp(compext, "Z") == 0 )
      {
         if( compression != NULL )
            *compression = compext;
         *lastdot = '\0';
      }

      /* find again the last dot in the filename without compression extension */
      lastdot = strrchr(filename, '.');
      if( lastslash != NULL && lastdot != NULL && lastdot < lastslash ) /* is the last dot belonging to the path? */
         lastdot = NULL;
   }
#endif

   if( lastslash == NULL )
   {
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

   if( lastdot != NULL )
   {
      if( extension != NULL )
         *extension = lastdot+1;
      *lastdot = '\0';
   }
}




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPrelDiff

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
SCIP_Real SCIPrelDiff(
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real absval1;
   SCIP_Real absval2;
   SCIP_Real quot;

   absval1 = REALABS(val1);
   absval2 = REALABS(val2);
   quot = MAX3(1.0, absval1, absval2);
   
   return (val1-val2)/quot;
}
