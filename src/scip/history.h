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
#pragma ident "@(#) $Id: history.h,v 1.7 2004/04/27 15:50:00 bzfpfend Exp $"

/**@file   history.h
 * @brief  internal methods for branching and inference history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __HISTORY_H__
#define __HISTORY_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_history.h"

#ifdef NDEBUG
#include "struct_history.h"
#endif



/** creates an empty history entry */
extern
RETCODE SCIPhistoryCreate(
   HISTORY**        history,            /**< pointer to store branching and inference history */
   MEMHDR*          memhdr              /**< block memory */
   );

/** frees a history entry */
extern
void SCIPhistoryFree(
   HISTORY**        history,            /**< pointer to branching and inference history */
   MEMHDR*          memhdr              /**< block memory */
   );

/** resets history entry to zero */
extern
void SCIPhistoryReset(
   HISTORY*         history             /**< branching and inference history */
   );

/** unites two history entries by adding the values of the second one to the first one */
extern
void SCIPhistoryUnite(
   HISTORY*         history,            /**< branching and inference history */
   HISTORY*         addhistory,         /**< history values to add to history */
   Bool             switcheddirs        /**< should the history entries be united with switched directories */
   );
   
/** updates the pseudo costs for a change of "solvaldelta" in the variable's LP solution value and a change of "objdelta"
 *  in the LP's objective value
 */
extern
void SCIPhistoryUpdatePseudocost(
   HISTORY*         history,            /**< branching and inference history */
   const SET*       set,                /**< global SCIP settings */
   Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   Real             weight              /**< weight of this update in pseudo cost sum (added to pscostcount) */
   );


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the opposite direction of the given branching direction */
extern
BRANCHDIR SCIPbranchdirOpposite(
   BRANCHDIR        dir                 /**< branching direction */
   );

/** returns the expected dual gain for moving the corresponding variable by "solvaldelta" */
extern
Real SCIPhistoryGetPseudocost(
   HISTORY*         history,            /**< branching and inference history */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** returns the (possible fractional) number of (partial) pseudo cost updates performed on this pseudo cost entry in 
 *  the given branching direction
 */
extern
Real SCIPhistoryGetPseudocostCount(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** returns whether the pseudo cost entry is empty in the given branching direction (whether no value was added yet) */
extern
Bool SCIPhistoryIsPseudocostEmpty(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** increases the number of branchings counter */
extern
void SCIPhistoryIncNBranchings(
   HISTORY*         history,            /**< branching and inference history */
   int              depth,              /**< depth at which the bound change took place */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** increases the number of inferences counter */
extern
void SCIPhistoryIncNInferences(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** get number of branchings counter */
extern
Longint SCIPhistoryGetNBranchings(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** get number of branchings counter */
extern
Longint SCIPhistoryGetNInferences(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** returns the average number of inferences per branching */
extern
Real SCIPhistoryGetAvgInferences(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** returns the average depth of bound changes due to branching */
extern
Real SCIPhistoryGetAvgBranchdepth(
   HISTORY*         history,            /**< branching and inference history */
   BRANCHDIR        dir                 /**< branching direction */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPbranchdirOpposite(dir)                 ((dir) == SCIP_BRANCHDIR_DOWNWARDS \
                                                   ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS)
#define SCIPhistoryGetPseudocost(history,solvaldelta)                                       \
   ( (solvaldelta) >= 0.0 ? (solvaldelta) * ((history)->pscostcount[1] > 0.0                \
                            ? (history)->pscostsum[1] / (history)->pscostcount[1] : 1.0)    \
                          : -(solvaldelta) * ((history)->pscostcount[0] > 0.0               \
                            ? (history)->pscostsum[0] / (history)->pscostcount[0] : 1.0) )
#define SCIPhistoryGetPseudocostCount(history,dir) ((history)->pscostcount[dir])
#define SCIPhistoryIsPseudocostEmpty(history,dir)  ((history)->pscostcount[dir] == 0.0)
#define SCIPhistoryIncNBranchings(history,depth,dir) { (history)->nbranchings[dir]++; \
                                                       (history)->branchdepthsum[dir] += depth; }
#define SCIPhistoryIncNInferences(history,dir)     (history)->ninferences[dir]++;
#define SCIPhistoryGetNBranchings(history,dir)     ((history)->nbranchings[dir])
#define SCIPhistoryGetNInferences(history,dir)     ((history)->ninferences[dir])
#define SCIPhistoryGetAvgInferences(history,dir)   ((history)->nbranchings[dir] > 0 \
                                                   ? (Real)(history)->ninferences[dir]/(Real)(history)->nbranchings[dir] \
                                                   : 0)
#define SCIPhistoryGetAvgBranchdepth(history,dir)  ((history)->nbranchings[dir] > 0 \
                                                   ? (Real)(history)->branchdepthsum[dir]/(Real)(history)->nbranchings[dir] \
                                                   : 0)

#endif


#endif
