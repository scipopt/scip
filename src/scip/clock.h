/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   clock.h
 * @brief  methods and datastructures for taking timings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CLOCK_H__
#define __CLOCK_H__

enum ClockType
{
   SCIP_CLOCKTYPE_DEFAULT = 0,          /**< use default clock type */
   SCIP_CLOCKTYPE_CPU     = 1,          /**< use CPU clock */
   SCIP_CLOCKTYPE_WALL    = 2           /**< use wall clock */
};
typedef enum ClockType CLOCKTYPE;       /**< clock type to use */

typedef struct Clock CLOCK;             /**< clock timer */
typedef struct CPUClock CPUCLOCK;       /**< CPU clock counter */
typedef struct WallClock WALLCLOCK;     /**< wall clock counter */

#include <sys/times.h>
#include <sys/time.h>
#include <time.h>

#include "def.h"
#include "retcode.h"
#include "set.h"


/** CPU clock counter */
struct CPUClock
{
   clock_t          user;               /**< clock ticks for user CPU time */
};

/** wall clock counter */
struct WallClock
{
   long             sec;                /**< seconds counter */
   long             usec;               /**< microseconds counter */
};

/** clock timer */
struct Clock
{
   union
   {
      CPUCLOCK      cpuclock;           /**< CPU clock counter */
      WALLCLOCK     wallclock;          /**< wall clock counter */
   } data;
   CLOCKTYPE        clocktype;          /**< current type of clock used */
   Bool             usedefault;         /**< should the clock's type be overruled by the default clock type? */
   Bool             enabled;            /**< should the clock be used? */
   int              nruns;              /**< number of SCIPclockStart() calls without SCIPclockStop() calls */
};


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
