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
#pragma ident "@(#) $Id: clock.h,v 1.6 2004/02/04 17:27:17 bzfpfend Exp $"

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
   CLOCK**          clock,              /**< pointer to clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   );

/** frees a clock */
extern
void SCIPclockFree(
   CLOCK**          clock               /**< pointer to clock timer */
   );

/** initializes and resets a clock */
extern
void SCIPclockInit(
   CLOCK*           clock,              /**< clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   );

/** completely stop the clock and reset the clock's counter to zero */
extern
void SCIPclockReset(
   CLOCK*           clock               /**< clock timer */
   );

/** enables the clock */
extern
void SCIPclockEnable(
   CLOCK*           clock               /**< clock timer */
   );

/** disables and resets the clock */
extern
void SCIPclockDisable(
   CLOCK*           clock               /**< clock timer */
   );

/** sets the type of the clock, overriding the default clock type, and resets the clock */
extern
void SCIPclockSetType(
   CLOCK*           clock,              /**< clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   );

/** starts measurement of time in the given clock, update the clock's type if it is bound to the default type */
extern
void SCIPclockStart(
   CLOCK*           clock,              /**< clock timer */
   const SET*       set                 /**< global SCIP settings */
   );

/** stops measurement of time in the given clock */
extern
void SCIPclockStop(
   CLOCK*           clock,              /**< clock timer */
   const SET*       set                 /**< global SCIP settings */
   );

/** gets the used time of this clock in seconds */
extern
Real SCIPclockGetTime(
   CLOCK*           clock               /**< clock timer */
   );

/** sets the used time of this clock in seconds */
extern
void SCIPclockSetTime(
   CLOCK*           clock,              /**< clock timer */
   Real             sec                 /**< time in seconds to set the clock's timer to */
   );


#endif
