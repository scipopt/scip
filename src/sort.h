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

/**@file   sort.h
 * @brief  datastructures and algorithms for sorting and queueing data elements
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SORT_H__
#define __SORT_H__


typedef struct PQueue PQUEUE;           /**< priority queue */



#include "def.h"
#include "retcode.h"



/** compares two data element pointers
 *  result:
 *    < 0: elem1 comes before (is better than) elem2
 *    = 0: both elements have the same value
 *    > 0: elem2 comes after (is worse than) elem2
 */
#define DECL_SORTPTRCOMP(x) int x (void* elem1, void* elem2)


extern
RETCODE SCIPpqueueInit(                 /**< initializes priority queue */
   PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int              initsize,           /**< initial number of available element slots */
   Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );

extern
void SCIPpqueueFree(                    /**< frees priority queue, but not the data elements themselves */
   PQUEUE**         pqueue              /**< pointer to a priority queue */
   );

extern
RETCODE SCIPpqueueInsert(               /**< inserts element into priority queue */
   PQUEUE*          pqueue,             /**< pointer to a priority queue */
   void*            elem                /**< element to be inserted */
   );

extern
void* SCIPpqueueRemove(                 /**< removes and returns best element from the priority queue */
   PQUEUE*          pqueue              /**< pointer to a priority queue */
   );

extern
void* SCIPpqueueFirst(                  /**< returns the best element of the queue without removing it */
   const PQUEUE*    pqueue              /**< pointer to a priority queue */
   );

extern
void SCIPbsortPtrDbl(                   /**< bubble sort of two joint arrays of pointers/Reals, sorted by first array */
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            dblarray,           /**< Real array to be permuted in the same way */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );


#endif
