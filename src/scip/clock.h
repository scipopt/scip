/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: clock.h,v 1.12 2005/01/21 09:16:48 bzfpfend Exp $"

/**@file   clock.h
 * @brief  internal methods for clocks and timing issues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CLOCK_H__
#define __CLOCK_H__


#include "def.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_clock.h"


/** creates a clock and initializes it */
extern
RETCODE SCIPclockCreate(
   CLOCK**          clck,               /**< pointer to clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   );

/** frees a clock */
extern
void SCIPclockFree(
   CLOCK**          clck                /**< pointer to clock timer */
   );

/** initializes and resets a clock */
extern
void SCIPclockInit(
   CLOCK*           clck,               /**< clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   );

/** completely stop the clock and reset the clock's counter to zero */
extern
void SCIPclockReset(
   CLOCK*           clck                /**< clock timer */
   );

/** enables the clock */
extern
void SCIPclockEnable(
   CLOCK*           clck                /**< clock timer */
   );

/** disables and resets the clock */
extern
void SCIPclockDisable(
   CLOCK*           clck                /**< clock timer */
   );

/** sets the type of the clock, overriding the default clock type, and resets the clock */
extern
void SCIPclockSetType(
   CLOCK*           clck,               /**< clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   );

/** starts measurement of time in the given clock, update the clock's type if it is bound to the default type */
extern
void SCIPclockStart(
   CLOCK*           clck,               /**< clock timer */
   SET*             set                 /**< global SCIP settings */
   );

/** stops measurement of time in the given clock */
extern
void SCIPclockStop(
   CLOCK*           clck,               /**< clock timer */
   SET*             set                 /**< global SCIP settings */
   );

/** returns whether the clock is currently running */
extern
Bool SCIPclockIsRunning(
   CLOCK*           clck                /**< clock timer */
   );

/** gets the used time of this clock in seconds */
extern
Real SCIPclockGetTime(
   CLOCK*           clck                /**< clock timer */
   );

/** sets the used time of this clock in seconds */
extern
void SCIPclockSetTime(
   CLOCK*           clck,               /**< clock timer */
   Real             sec                 /**< time in seconds to set the clock's timer to */
   );

/** gets current time of day in seconds (standard time zone) */
extern
Real SCIPclockGetTimeOfDay(
   void
   );


#endif
