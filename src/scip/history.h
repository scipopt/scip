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
#pragma ident "@(#) $Id: history.h,v 1.2 2004/02/04 17:27:26 bzfpfend Exp $"

/**@file   history.h
 * @brief  internal methods for branching history
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
   HISTORY**        history,            /**< pointer to store branching history */
   MEMHDR*          memhdr              /**< block memory */
   );

/** frees a history entry */
extern
void SCIPhistoryFree(
   HISTORY**        history,            /**< pointer to branching history */
   MEMHDR*          memhdr              /**< block memory */
   );

/** resets history entry to zero */
extern
void SCIPhistoryReset(
   HISTORY*         history             /**< branching history */
   );

/** updates the history for a change of "solvaldelta" in the variable's LP solution value and a change of "objdelta"
 *  in the LP's objective value
 */
extern
void SCIPhistoryUpdate(
   HISTORY*         history,            /**< branching history */
   const SET*       set,                /**< global SCIP settings */
   Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   Real             weight              /**< weight of this update in history sum (added to count) */
   );


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the expected dual gain for moving the corresponding variable by "solvaldelta" */
extern
Real SCIPhistoryGetValue(
   HISTORY*         history,            /**< branching history */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** returns the (possible fractional) number of (partial) history updates performed on this history entry in 
 *  the given direction
 */
extern
Real SCIPhistoryGetCount(
   HISTORY*         history,            /**< branching history */
   int              dir                 /**< direction: downwards (0), or upwards (1) */
   );

/** returns whether the history entry is empty in the given direction (whether no value was added since initialization) */
extern
Bool SCIPhistoryIsEmpty(
   HISTORY*         history,            /**< branching history */
   int              dir                 /**< direction: downwards (0), or upwards (1) */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPhistoryGetValue(history,solvaldelta) \
   ( (solvaldelta) >= 0.0 ? (solvaldelta) * ((history)->count[1] > 0.0 ? (history)->sum[1] / (history)->count[1] : 1.0) \
                          : -(solvaldelta) * ((history)->count[0] > 0.0 ? (history)->sum[0] / (history)->count[0] : 1.0) )
#define SCIPhistoryGetCount(history,dir) ((history)->count[dir])
#define SCIPhistoryIsEmpty(history,dir)  ((history)->count[dir] == 0.0)

#endif


#endif
