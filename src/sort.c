/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sort.c
 * @brief  datastructures and algorithms for sorting and queueing data elements
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

#include "memory.h"
#include "retcode.h"
#include "sort.h"


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


static
RETCODE pqueueResize(                   /**< resizes element memory to hold at least the given number of elements */
   PQUEUE*          pqueue,             /**< pointer to a priority queue */
   int              minsize             /**< minimal number of storeable elements */
   )
{
   assert(pqueue != NULL);
   
   if( minsize <= pqueue->size )
      return SCIP_OKAY;

   pqueue->size = MAX(minsize, (int)(pqueue->size * pqueue->sizefac));
   ALLOC_OKAY( reallocMemoryArray(pqueue->slots, pqueue->size) );

   return SCIP_OKAY;
}

RETCODE SCIPpqueueInit(                 /**< initializes priority queue */
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

   ALLOC_OKAY( allocMemory(*pqueue) );
   (*pqueue)->len = 0;
   (*pqueue)->size = 0;
   (*pqueue)->sizefac = sizefac;
   (*pqueue)->slots = NULL;
   (*pqueue)->ptrcmp = ptrcmp;
   CHECK_OKAY( pqueueResize(*pqueue, initsize) );

   return SCIP_OKAY;
}

void SCIPpqueueFree(                    /**< frees priority queue, but not the data elements themselves */
   PQUEUE**         pqueue              /**< pointer to a priority queue */
   )
{
   assert(pqueue != NULL);

   freeMemoryArray((*pqueue)->slots);
   freeMemory(*pqueue);
}

RETCODE SCIPpqueueInsert(               /**< inserts element into priority queue */
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

void* SCIPpqueueRemove(                 /**< removes and returns best element from the priority queue */
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

void* SCIPpqueueFirst(                  /**< returns the best element of the queue without removing it */
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




void SCIPbsortPtr(                      /**< bubble sort of an array of pointers */
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

void SCIPbsortPtrDbl(                   /**< bubble sort of two joint arrays of pointers/Reals, sorted by first array */
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
