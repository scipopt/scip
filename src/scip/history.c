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

/**@file   history.c
 * @brief  methods for branching and inference history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/history.h"

#ifndef NDEBUG
#include "scip/struct_history.h"
#endif

/*
 * methods for branching and inference history
 */

/** creates an empty history entry */
SCIP_RETCODE SCIPhistoryCreate(
   SCIP_HISTORY**        history,            /**< pointer to store branching and inference history */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(history != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, history) );

   SCIPhistoryReset(*history);

   return SCIP_OKAY;
}

/** frees a history entry */
void SCIPhistoryFree(
   SCIP_HISTORY**        history,            /**< pointer to branching and inference history */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(history != NULL);
   assert(*history != NULL);

   BMSfreeBlockMemory(blkmem, history);
}

/** resets history entry to zero */
void SCIPhistoryReset(
   SCIP_HISTORY*         history             /**< branching and inference history */
   )
{
   assert(history != NULL);

   history->pscostcount[0] = 0.0;
   history->pscostcount[1] = 0.0;
   history->pscostsum[0] = 0.0;
   history->pscostsum[1] = 0.0;
   history->vsids[0] = 0.0;
   history->vsids[1] = 0.0;
   history->conflengthsum[0] = 0.0;
   history->conflengthsum[1] = 0.0;
   history->inferencesum[0] = 0.0;
   history->inferencesum[1] = 0.0;
   history->cutoffsum[0] = 0.0;
   history->cutoffsum[1] = 0.0;
   history->nactiveconflicts[0] = 0;
   history->nactiveconflicts[1] = 0;
   history->nbranchings[0] = 0;
   history->nbranchings[1] = 0;
   history->branchdepthsum[0] = 0;
   history->branchdepthsum[1] = 0;
}

/** unites two history entries by adding the values of the second one to the first one */
void SCIPhistoryUnite(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_HISTORY*         addhistory,         /**< history values to add to history */
   SCIP_Bool             switcheddirs        /**< should the history entries be united with switched directories */
   )
{
   int d;

   assert(history != NULL);
   assert(addhistory != NULL);

   d = switcheddirs ? 1 : 0;

   history->pscostcount[0] += addhistory->pscostcount[d];
   history->pscostcount[1] += addhistory->pscostcount[1-d];
   history->pscostsum[0] += addhistory->pscostsum[d];
   history->pscostsum[1] += addhistory->pscostsum[1-d];
   history->vsids[0] += addhistory->vsids[d];
   history->vsids[1] += addhistory->vsids[1-d];
   history->conflengthsum[0] += addhistory->conflengthsum[d];
   history->conflengthsum[1] += addhistory->conflengthsum[1-d];
   history->inferencesum[0] += addhistory->inferencesum[d];
   history->inferencesum[1] += addhistory->inferencesum[1-d];
   history->cutoffsum[0] += addhistory->cutoffsum[d];
   history->cutoffsum[1] += addhistory->cutoffsum[1-d];
   history->nactiveconflicts[0] += addhistory->nactiveconflicts[d];
   history->nactiveconflicts[1] += addhistory->nactiveconflicts[1-d];
   history->nbranchings[0] += addhistory->nbranchings[d];
   history->nbranchings[1] += addhistory->nbranchings[1-d];
   history->branchdepthsum[0] += addhistory->branchdepthsum[d];
   history->branchdepthsum[1] += addhistory->branchdepthsum[1-d];
}

/** updates the pseudo costs for a change of "solvaldelta" in the variable's LP solution value and a change of "objdelta"
 *  in the LP's objective value
 */
void SCIPhistoryUpdatePseudocost(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   SCIP_Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   SCIP_Real             weight              /**< weight of this update in pseudo cost sum (added to pscostcount) */
   )
{
   SCIP_Real distance;
   SCIP_Real eps;
   int dir;

   assert(history != NULL);
   assert(set != NULL);
   assert(!SCIPsetIsInfinity(set, REALABS(solvaldelta)));
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
   eps = SCIPsetPseudocosteps(set);
   distance = MAX(distance, eps);

   /* slightly increase objective delta, s.t. pseudo cost values are not zero, and fractionalities are
    * always used at least a bit
    */
   objdelta += SCIPsetPseudocostdelta(set);

   /* update the pseudo cost values */
   history->pscostcount[dir] += weight;
   history->pscostsum[dir] += weight * objdelta/distance;

   SCIPdebugMessage("updated pseudo costs of history %p: dir=%d, distance=%g, objdelta=%g, weight=%g  ->  %g/%g\n",
      (void*)history, dir, distance, objdelta, weight, history->pscostcount[dir], history->pscostsum[dir]);
}




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPbranchdirOpposite
#undef SCIPhistoryGetPseudocost
#undef SCIPhistoryGetPseudocostCount
#undef SCIPhistoryIsPseudocostEmpty
#undef SCIPhistoryIncVSIDS
#undef SCIPhistoryScaleVSIDS
#undef SCIPhistoryGetVSIDS
#undef SCIPhistoryIncNActiveConflicts
#undef SCIPhistoryGetNActiveConflicts
#undef SCIPhistoryGetAvgConflictlength
#undef SCIPhistoryIncNBranchings
#undef SCIPhistoryIncInferenceSum
#undef SCIPhistoryIncCutoffSum
#undef SCIPhistoryGetNBranchings
#undef SCIPhistoryGetInferenceSum
#undef SCIPhistoryGetAvgInferences
#undef SCIPhistoryGetCutoffSum
#undef SCIPhistoryGetAvgCutoffs
#undef SCIPhistoryGetAvgBranchdepth

/** returns the opposite direction of the given branching direction */
SCIP_BRANCHDIR SCIPbranchdirOpposite(
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   )
{
   return (dir == SCIP_BRANCHDIR_DOWNWARDS ? SCIP_BRANCHDIR_UPWARDS
      : (dir == SCIP_BRANCHDIR_UPWARDS ? SCIP_BRANCHDIR_DOWNWARDS : SCIP_BRANCHDIR_AUTO));
}

/** returns the expected dual gain for moving the corresponding variable by "solvaldelta" */
SCIP_Real SCIPhistoryGetPseudocost(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
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
SCIP_Real SCIPhistoryGetPseudocostCount(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->pscostcount[dir];
}

/** returns whether the pseudo cost entry is empty in the given branching direction (whether no value was added yet) */
SCIP_Bool SCIPhistoryIsPseudocostEmpty(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);
   
   return (history->pscostcount[dir] == 0.0);
}

/** increases the conflict score of the history entry by the given weight */
void SCIPhistoryIncVSIDS(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_Real             weight              /**< weight of this update in conflict score */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   history->vsids[dir] += weight;
}

/** scales the conflict score values with the given scalar */
void SCIPhistoryScaleVSIDS(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_Real             scalar              /**< scalar to multiply the conflict scores with */
   )
{
   assert(history != NULL);
   
   history->vsids[0] *= scalar;
   history->vsids[1] *= scalar;
}

/** gets the conflict score of the history entry */
SCIP_Real SCIPhistoryGetVSIDS(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->vsids[dir];
}

/** increases the number of active conflicts by one and the overall length of the history entry by the given weight */
void SCIPhistoryIncNActiveConflicts(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_Real             length              /**< length of the conflict */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1); 
   assert(length >= 0.0);
   
   history->nactiveconflicts[dir]++;
   history->conflengthsum[dir] += length;
}

/** gets the number of active conflicts of the history entry */
SCIP_Longint SCIPhistoryGetNActiveConflicts(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->nactiveconflicts[dir];
}

/** gets the average conflict length of the history entry */
SCIP_Real SCIPhistoryGetAvgConflictlength(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->conflengthsum[dir] > 0.0 ? (SCIP_Real)history->nactiveconflicts[dir]/(SCIP_Real)history->conflengthsum[dir] : 0.0;
}

/** increases the number of branchings counter */
void SCIPhistoryIncNBranchings(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   int                   depth,              /**< depth at which the bound change took place */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(history != NULL);
   assert(depth >= 1);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   history->nbranchings[dir]++;
   history->branchdepthsum[dir] += depth;
}

/** increases the number of inferences counter by a certain value */
void SCIPhistoryIncInferenceSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             weight              /**< weight of this update in inference score */
   ) 
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);
   assert(history->nbranchings[dir] >= 1);
   assert(weight >= 0.0);

   history->inferencesum[dir] += weight;
} 

/** increases the number of cutoffs counter */
void SCIPhistoryIncCutoffSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             weight              /**< weight of this update in cutoff score */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);
   assert(history->nbranchings[dir] >= 1);
   assert(weight >= 0.0);

   history->cutoffsum[dir] += weight;
}

/** get number of branchings counter */
SCIP_Longint SCIPhistoryGetNBranchings(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->nbranchings[dir];
}

/** get number of inferences counter */
SCIP_Real SCIPhistoryGetInferenceSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->inferencesum[dir];
}

/** returns the average number of inferences per branching */
SCIP_Real SCIPhistoryGetAvgInferences(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->nbranchings[dir] > 0 ? (SCIP_Real)history->inferencesum[dir]/(SCIP_Real)history->nbranchings[dir] : 0.0;
}

/** get number of cutoffs counter */
SCIP_Real SCIPhistoryGetCutoffSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->cutoffsum[dir];
}

/** returns the average number of cutoffs per branching */
SCIP_Real SCIPhistoryGetAvgCutoffs(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->nbranchings[dir] > 0 ? (SCIP_Real)history->cutoffsum[dir]/(SCIP_Real)history->nbranchings[dir] : 0.0;
}

/** returns the average depth of bound changes due to branching */
SCIP_Real SCIPhistoryGetAvgBranchdepth(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(history != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert((int)dir == 0 || (int)dir == 1);

   return history->nbranchings[dir] > 0 ? (SCIP_Real)history->branchdepthsum[dir]/(SCIP_Real)history->nbranchings[dir] : 1.0;
}
