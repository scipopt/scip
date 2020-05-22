/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   graphheaps.h
 * @brief  includes (inline) graph heap methods used for Steiner tree problems
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_GRAPHHEAPS_H_
#define APPLICATIONS_STP_SRC_GRAPHHEAPS_H_

#include "scip/scip.h"
#include "graphdefs.h"
#include "portab.h"



/** gets closest node from heap */
inline static int nearestX(
   int* RESTRICT heap,
   int* RESTRICT state,
   int* RESTRICT count,    /* pointer to store the number of elements on the heap */
   const SCIP_Real* pathdist)
{
   int   k;
   int   t;
   int   c;
   int   j;
   int   dcount;

   /* Heap shift down
    */
   k              = heap[1];
   j              = 1;
   c              = 2;
   heap[1]        = heap[(*count)--];
   state[heap[1]] = 1;

   dcount = *count;

   if (dcount > 2)
      if (LT(pathdist[heap[3]], pathdist[heap[2]]))
         c++;

   while((c <= dcount) && GT(pathdist[heap[j]], pathdist[heap[c]]))
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c             += c;

      if ((c + 1) <= dcount)
         if (LT(pathdist[heap[c + 1]], pathdist[heap[c]]))
            c++;
   }
   return(k);
}

/** corrects the heap */
inline static void correctX(
   int* RESTRICT heap,
   int* RESTRICT state,
   int* RESTRICT count,    /* pointer to store the number of elements on the heap */
   SCIP_Real* RESTRICT pathdist,
   int* RESTRICT pathedge,
   int    l,
   int    k,
   int    e,
   SCIP_Real cost
   )
{
   int    t;
   int    c;
   int    j;

   pathdist[l] = (pathdist[k] + cost);

   if( pathedge != NULL )
      pathedge[l] = e;

   /* not yet in heap? */
   if( state[l] == UNKNOWN )
   {
      heap[++(*count)] = l;
      state[l]      = (*count);
   }

   /* Heap shift up
    */
   j = state[l];
   c = j / 2;

   while( (j > 1) && pathdist[heap[c]] > pathdist[heap[j]] )
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c              = j / 2;
   }
}


/** resets element from heap */
inline static void resetX(
   SCIP_Real* RESTRICT pathdist,
   int* RESTRICT heap,
   int* RESTRICT state,
   int* RESTRICT count,
   int node,
   SCIP_Real distnew
   )
{
   int c;
   int j;

   pathdist[node] = distnew;

   heap[++(*count)] = node;
   state[node] = (*count);

   /* heap shift up */
   j = state[node];
   c = j / 2;

   while( (j > 1) && pathdist[heap[c]] > pathdist[heap[j]] )
   {
      const int t = heap[c];
      heap[c] = heap[j];
      heap[j] = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j = c;
      c = j / 2;
   }
}



#endif /* APPLICATIONS_STP_SRC_GRAPHHEAPS_H_ */
