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
#pragma ident "@(#) $Id: history.c,v 1.7 2004/04/15 10:41:23 bzfpfend Exp $"

/**@file   history.c
 * @brief  methods for branching and inference history
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
 * methods for branching and inference history
 */

/** creates an empty history entry */
RETCODE SCIPhistoryCreate(
   HISTORY**        history,            /**< pointer to store branching and inference history */
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
   HISTORY**        history,            /**< pointer to branching and inference history */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(history != NULL);
   assert(*history != NULL);

   freeBlockMemory(memhdr, history);
}

/** resets history entry to zero */
void SCIPhistoryReset(
   HISTORY*         history             /**< branching and inference history */
   )
{
   assert(history != NULL);

   history->pscostcount[0] = 0.0;
   history->pscostcount[1] = 0.0;
   history->pscostsum[0] = 0.0;
   history->pscostsum[1] = 0.0;
   history->nbranchings[0] = 0;
   history->nbranchings[1] = 0;
   history->ninferences[0] = 0;
   history->ninferences[1] = 0;
   history->branchdepthsum[0] = 0;
   history->branchdepthsum[1] = 0;
}

/** updates the pseudo costs for a change of "solvaldelta" in the variable's LP solution value and a change of "objdelta"
 *  in the LP's objective value
 */
void SCIPhistoryUpdatePseudocost(
   HISTORY*         history,            /**< branching and inference history */
   const SET*       set,                /**< global SCIP settings */
   Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   Real             weight              /**< weight of this update in pseudo cost sum (added to pscostcount) */
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
      /* the variable's solution value didn't change, and the pseudo costs cannot be updated */
      return;
   }
   assert(dir == 0 || dir == 1);
   assert(SCIPsetIsPositive(set, distance));

   /* apply a lower limit on the distance to avoid numerical instabilities due to very large summands */
   distance = MAX(distance, set->pseudocosteps);

   /* slightly increase objective delta, s.t. pseudo cost values are not zero, and fractionalities are
    * always used at least a bit
    */
   objdelta += set->pseudocostdelta;

   /* update the pseudo cost values */
   history->pscostcount[dir] += weight;
   history->pscostsum[dir] += weight * objdelta/distance;

   debugMessage("updated pseudo costs of history %p: dir=%d, distance=%g, objdelta=%g, weight=%g  ->  %g/%g\n",
      history, dir, distance, objdelta, weight, history->pscostcount[dir], history->pscostsum[dir]);
}


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the expected dual gain for moving the corresponding variable by "solvaldelta" */
Real SCIPhistoryGetPseudocost(
   HISTORY*         history,            /**< branching and inference history */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   assert(history != NULL);
   
   if( solvaldelta >= 0.0 )
      return solvaldelta * (history->pscostcount[1] > 0.0 ? history->pscostsum[1] / history->pscostcount[1] : 1.0);
   else
      return -solvaldelta * (history->pscostcount[0] > 0.0 ? history->pscostsum[0] / history->pscostcount[0] : 1.0);
}

/** returns the (possible fractional) number of (partial) pseudo cost updates performed on this pseudo cost entry in 
 *  the given branching direction
 */
Real SCIPhistoryGetPseudocostCount(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == 0 || dir == 1);

   return history->pscostcount[dir];
}

/** returns whether the pseudo cost entry is empty in the given branching direction (whether no value was added yet) */
Bool SCIPhistoryIsPseudocostEmpty(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == 0 || dir == 1);
   
   return (history->pscostcount[dir] == 0.0);
}

/** increases the number of branchings counter */
void SCIPhistoryIncNBranchings(
   HISTORY*         history,            /**< branching and inference history */
   int              depth,              /**< depth at which the bound change took place */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(depth >= 1);
   assert(dir == 0 || dir == 1);

   history->nbranchings[dir]++;
   history->branchdepthsum[dir] += depth;
}

/** increases the number of inferences counter */
void SCIPhistoryIncNInferences(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == 0 || dir == 1);

   history->ninferences[dir]++;
}

/** get number of branchings counter */
Longint SCIPhistoryGetNBranchings(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == 0 || dir == 1);

   return history->nbranchings[dir];
}

/** get number of branchings counter */
Longint SCIPhistoryGetNInferences(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == 0 || dir == 1);

   return history->ninferences[dir];
}

/** returns the average number of inferences per branching */
Real SCIPhistoryGetAvgInferences(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == 0 || dir == 1);

   return history->nbranchings[dir] > 0 ? (Real)history->ninferences[dir]/(Real)history->nbranchings[dir] : 0;
}

/** returns the average depth of bound changes due to branching */
Real SCIPhistoryGetAvgBranchdepth(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == 0 || dir == 1);

   return history->nbranchings[dir] > 0 ? (Real)history->branchdepthsum[dir]/(Real)history->nbranchings[dir] : 0;
}

#endif
