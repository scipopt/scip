/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: history.c,v 1.1 2004/01/07 13:14:13 bzfpfend Exp $"

/**@file   history.c
 * @brief  methods for branching history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "set.h"
#include "history.h"

#ifndef NDEBUG
#include "struct_history.h"
#endif




/*
 * methods for branching history
 */

/** creates an empty history entry */
RETCODE SCIPhistoryCreate(
   HISTORY**        history,            /**< pointer to store branching history */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(history != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, history) );

   SCIPhistoryReset(*history);

   return SCIP_OKAY;
}

/** frees a history entry */
void SCIPhistoryFree(
   HISTORY**        history,            /**< pointer to branching history */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(history != NULL);
   assert(*history != NULL);

   freeBlockMemory(memhdr, history);
}

/** resets history entry to zero */
void SCIPhistoryReset(
   HISTORY*         history             /**< branching history */
   )
{
   assert(history != NULL);

   history->count[0] = 0.0;
   history->count[1] = 0.0;
   history->sum[0] = 0.0;
   history->sum[1] = 0.0;
}

/** updates the history for a change of "solvaldelta" in the variable's LP solution value and a change of "objdelta"
 *  in the LP's objective value
 */
void SCIPhistoryUpdate(
   HISTORY*         history,            /**< branching history */
   const SET*       set,                /**< global SCIP settings */
   Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   Real             weight              /**< weight of this update in history sum (added to count) */
   )
{
   Real distance;
   int dir;

   assert(history != NULL);
   assert(set != NULL);
   assert(!SCIPsetIsInfinity(set, ABS(solvaldelta)));
   assert(!SCIPsetIsInfinity(set, objdelta));
   assert(!SCIPsetIsNegative(set, objdelta));
   assert(0.0 < weight && weight <= 1.0);
   
   if( SCIPsetIsPositive(set, solvaldelta) )
   {
      /* variable's solution value moved upwards */
      dir = 1;
      distance = solvaldelta;
   }
   else if( SCIPsetIsNegative(set, solvaldelta) )
   {
      /* variable's solution value moved downwards */
      dir = 0;
      distance = -solvaldelta;
   }
   else
   {
      /* the variable's solution value didn't change, and the history cannot be updated */
      return;
   }
   assert(dir == 0 || dir == 1);
   assert(SCIPsetIsPositive(set, distance));

   /* apply a lower limit on the distance to avoid numerical instabilities due to very large summands */
   distance = MAX(distance, set->historyeps);

   /* update the history values */
   history->count[dir] += weight;
   history->sum[dir] += weight * objdelta/distance;

   debugMessage("updated history %p: dir=%d, distance=%g, objdelta=%g, weight=%g  ->  %g/%g\n",
      history, dir, distance, objdelta, weight, history->count[dir], history->sum[dir]);
}


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the expected dual gain for moving the corresponding variable by "solvaldelta" */
Real SCIPhistoryGetValue(
   HISTORY*         history,            /**< branching history */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   assert(history != NULL);
   
   if( solvaldelta >= 0.0 )
      return solvaldelta * (history->count[1] > 0.0 ? history->sum[1] / history->count[1] : 1.0);
   else
      return -solvaldelta * (history->count[0] > 0.0 ? history->sum[0] / history->count[0] : 1.0);
}

/** returns the (possible fractional) number of (partial) history updates performed on this history entry in 
 *  the given direction
 */
Real SCIPhistoryGetCount(
   HISTORY*         history,            /**< branching history */
   int              dir                 /**< direction: downwards (0), or upwards (1) */
   )
{
   assert(history != NULL);
   assert(dir == 0 || dir == 1);

   return history->count[dir];
}

/** returns whether the history entry is empty in the given direction (whether no value was added since initialization) */
Bool SCIPhistoryIsEmpty(
   HISTORY*         history,            /**< branching history */
   int              dir                 /**< direction: downwards (0), or upwards (1) */
   )
{
   assert(history != NULL);
   assert(dir == 0 || dir == 1);
   
   return (history->count[dir] == 0.0);
}

#endif
