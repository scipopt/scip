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
#pragma ident "@(#) $Id: history.h,v 1.3 2004/04/06 15:21:05 bzfpfend Exp $"

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

/** returns the expected dual gain for moving the corresponding variable by "solvaldelta" */
extern
Real SCIPhistoryGetPseudocost(
   HISTORY*         history,            /**< branching and inference history */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** returns the (possible fractional) number of (partial) pseudo cost updates performed on this pseudo cost entry in 
 *  the given direction
 */
extern
Real SCIPhistoryGetPseudocostCount(
   HISTORY*         history,            /**< branching and inference history */
   int              dir                 /**< direction: downwards (0), or upwards (1) */
   );

/** returns whether the pseudo cost entry is empty in the given direction (whether no value was added yet) */
extern
Bool SCIPhistoryIsPseudocostEmpty(
   HISTORY*         history,            /**< branching and inference history */
   int              dir                 /**< direction: downwards (0), or upwards (1) */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPhistoryGetPseudocost(history,solvaldelta)                                       \
   ( (solvaldelta) >= 0.0 ? (solvaldelta) * ((history)->pscostcount[1] > 0.0                \
                            ? (history)->pscostsum[1] / (history)->pscostcount[1] : 1.0)    \
                          : -(solvaldelta) * ((history)->pscostcount[0] > 0.0               \
                            ? (history)->pscostsum[0] / (history)->pscostcount[0] : 1.0) )
#define SCIPhistoryGetPseudocostCount(history,dir) ((history)->pscostcount[dir])
#define SCIPhistoryIsPseudocostEmpty(history,dir)  ((history)->pscostcount[dir] == 0.0)

#endif


#endif
