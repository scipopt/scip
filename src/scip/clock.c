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
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: clock.c,v 1.14 2005/01/18 09:26:42 bzfpfend Exp $"

/**@file   clock.c
 * @brief  methods for clocks and timing issues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>

#include "def.h"
#include "message.h"
#include "memory.h"
#include "set.h"
#include "clock.h"

#include "struct_clock.h"


/*lint -esym(*,timeval)*/
/*lint -esym(*,gettimeofday)*/



/** sets the clock's type and converts the clock timer accordingly */
static
void clockSetType(
   CLOCK*           clock,              /**< clock timer */
   CLOCKTYPE        newtype             /**< new clock type */
   )
{
   assert(clock != NULL);
   assert(newtype != SCIP_CLOCKTYPE_DEFAULT);

   if( clock->clocktype != newtype )
   {
      if( clock->clocktype == SCIP_CLOCKTYPE_DEFAULT )
      {
         assert(clock->nruns == 0);
         clock->clocktype = newtype;
         SCIPclockReset(clock);
         debugMessage("switched clock type to %d\n", newtype);
      }
      else
      {
         Real sec;
         
         sec = SCIPclockGetTime(clock);
         clock->clocktype = newtype;
         SCIPclockSetTime(clock, sec);
         debugMessage("switched clock type to %d (%g seconds -> %g seconds)\n", newtype, sec, SCIPclockGetTime(clock));
      }
   }
}

/** if the clock uses the default clock type and the default changed, converts the clock timer to the new type */
static
void clockUpdateDefaultType(
   CLOCK*           clock,              /**< clock timer */
   CLOCKTYPE        defaultclocktype    /**< default type of clock to use */
   )
{
   assert(clock != NULL);
   assert(defaultclocktype != SCIP_CLOCKTYPE_DEFAULT);

   if( clock->usedefault && clock->clocktype != defaultclocktype )
      clockSetType(clock, defaultclocktype);
}

/** creates a clock and initializes it */
RETCODE SCIPclockCreate(
   CLOCK**          clock,              /**< pointer to clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert(clock != NULL);

   ALLOC_OKAY( allocMemory(clock) );

   SCIPclockInit(*clock, clocktype);

   return SCIP_OKAY;
}

/** frees a clock */
void SCIPclockFree(
   CLOCK**          clock               /**< pointer to clock timer */
   )
{
   assert(clock != NULL);
   
   freeMemory(clock);
}

/** initializes and resets a clock */
void SCIPclockInit(
   CLOCK*           clock,              /**< clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert(clock != NULL);

   debugMessage("initializing clock %p of type %d\n", clock, clocktype);
   clock->enabled = TRUE;
   SCIPclockSetType(clock, clocktype);
}

/** completely stop the clock and reset the clock's counter to zero */
void SCIPclockReset(
   CLOCK*           clock               /**< clock timer */
   )
{
   assert(clock != NULL);

   debugMessage("resetting clock %p of type %d (usedefault=%d)\n", clock, clock->clocktype, clock->usedefault);
   switch( clock->clocktype )
   {
   case SCIP_CLOCKTYPE_DEFAULT:
      break;
   case SCIP_CLOCKTYPE_CPU:
      clock->data.cpuclock.user = 0;
      break;
   case SCIP_CLOCKTYPE_WALL:
      clock->data.wallclock.sec = 0;
      clock->data.wallclock.usec = 0;
      break;
   default:
      errorMessage("invalid clock type\n");
      abort();
   }
   clock->nruns = 0;
}

/** enables the clock */
void SCIPclockEnable(
   CLOCK*           clock               /**< clock timer */
   )
{
   assert(clock != NULL);

   debugMessage("enabling clock %p of type %d (usedefault=%d)\n", clock, clock->clocktype, clock->usedefault);

   clock->enabled = TRUE;
}

/** disables and resets the clock */
void SCIPclockDisable(
   CLOCK*           clock               /**< clock timer */
   )
{
   assert(clock != NULL);

   debugMessage("disabling clock %p of type %d (usedefault=%d)\n", clock, clock->clocktype, clock->usedefault);

   clock->enabled = FALSE;
   SCIPclockReset(clock);
}

/** sets the type of the clock, overriding the default clock type, and resets the clock */
void SCIPclockSetType(
   CLOCK*           clock,              /**< clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert(clock != NULL);

   debugMessage("setting type of clock %p (type %d, usedefault=%d) to %d\n", 
      clock, clock->clocktype, clock->usedefault, clocktype);

   clock->clocktype = clocktype;
   clock->usedefault = (clocktype == SCIP_CLOCKTYPE_DEFAULT);
   SCIPclockReset(clock);
}

/** starts measurement of time in the given clock */
void SCIPclockStart(
   CLOCK*           clock,              /**< clock timer */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(clock != NULL);
   assert(set != NULL);

   if( set->time_enabled && clock->enabled )
   {
      clockUpdateDefaultType(clock, set->time_clocktype);

      if( clock->nruns == 0 )
      {
         struct timeval tp; /*lint !e86*/
         struct tms now;
         
         debugMessage("starting clock %p (type %d, usedefault=%d)\n", clock, clock->clocktype, clock->usedefault);

         switch( clock->clocktype )
         {
         case SCIP_CLOCKTYPE_CPU:
            (void)times(&now);
            clock->data.cpuclock.user -= now.tms_utime;
            break;

         case SCIP_CLOCKTYPE_WALL:            
            gettimeofday(&tp, NULL);
            if( tp.tv_usec > clock->data.wallclock.usec ) /*lint !e115 !e40*/
            {
               clock->data.wallclock.sec -= (tp.tv_sec + 1); /*lint !e115 !e40*/
               clock->data.wallclock.usec += (1000000 - tp.tv_usec); /*lint !e115 !e40*/
            }
            else
            {
               clock->data.wallclock.sec -= tp.tv_sec; /*lint !e115 !e40*/
               clock->data.wallclock.usec -= tp.tv_usec; /*lint !e115 !e40*/
            }
            break;

         case SCIP_CLOCKTYPE_DEFAULT:            
         default:
            errorMessage("invalid clock type\n");
            abort();
         }
      }
      clock->nruns++;
   }
}

/** stops measurement of time in the given clock */
void SCIPclockStop(
   CLOCK*           clock,              /**< clock timer */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(clock != NULL);
   assert(set != NULL);

   if( set->time_enabled && clock->enabled )
   {
      assert(clock->nruns >= 1);

      clock->nruns--;
      if( clock->nruns == 0 )
      {
         struct timeval tp; /*lint !e86*/
         struct tms now;
         
         debugMessage("stopping clock %p (type %d, usedefault=%d)\n", clock, clock->clocktype, clock->usedefault);

         switch( clock->clocktype )
         {
         case SCIP_CLOCKTYPE_CPU:
            (void)times(&now);
            clock->data.cpuclock.user += now.tms_utime;
            break;

         case SCIP_CLOCKTYPE_WALL:            
            gettimeofday(&tp, NULL);
            if( tp.tv_usec + clock->data.wallclock.usec > 1000000 ) /*lint !e115 !e40*/
            {
               clock->data.wallclock.sec += (tp.tv_sec + 1); /*lint !e115 !e40*/
               clock->data.wallclock.usec -= (1000000 - tp.tv_usec); /*lint !e115 !e40*/
            }
            else
            {
               clock->data.wallclock.sec += tp.tv_sec; /*lint !e115 !e40*/
               clock->data.wallclock.usec += tp.tv_usec; /*lint !e115 !e40*/
            }
            break;

         case SCIP_CLOCKTYPE_DEFAULT:
         default:
            errorMessage("invalid clock type\n");
            abort();
         }
      }
   }
}

/** returns whether the clock is currently running */
Bool SCIPclockIsRunning(
   CLOCK*           clock               /**< clock timer */
   )
{
   assert(clock != NULL);

   return (clock->nruns > 0);
}

/** converts CPU clock ticks into seconds */
static
Real cputime2sec(
   clock_t          cputime             /**< clock ticks for CPU time */
   )
{
   clock_t clocks_per_second;

#ifndef CLK_TCK
   clocks_per_second = sysconf(_SC_CLK_TCK);
#else
   clocks_per_second = CLK_TCK;
#endif

   return (Real)cputime / (Real)clocks_per_second;
}

/** converts wall clock time into seconds */
static
Real walltime2sec(
   long             sec,                /**< seconds counter */
   long             usec                /**< microseconds counter */
   )
{
   return (Real)sec + 0.000001 * (Real)usec;
}

/** converts seconds into CPU clock ticks */
static
void sec2cputime(
   Real             sec,                /**< seconds */
   clock_t*         cputime             /**< pointer to store clock ticks for CPU time */
   )
{
   clock_t clocks_per_second;

   assert(cputime != NULL);

#ifndef CLK_TCK
   clocks_per_second = sysconf(_SC_CLK_TCK);
#else
   clocks_per_second = CLK_TCK;
#endif
   *cputime = (clock_t)(sec * clocks_per_second);
}

/** converts wall clock time into seconds */
static
void sec2walltime(
   Real             sec,                /**< seconds */
   long*            wallsec,            /**< pointer to store seconds counter */
   long*            wallusec            /**< pointer to store microseconds counter */
   )
{
   assert(wallsec != NULL);
   assert(wallusec != NULL);

   *wallsec = (long)sec;
   *wallusec = (long)(sec * 1000000.0) - (*wallsec);
}

/** gets the used time of this clock in seconds */
Real SCIPclockGetTime(
   CLOCK*           clock               /**< clock timer */
   )
{
   assert(clock != NULL);

   debugMessage("getting time of clock %p (type %d, usedefault=%d, nruns=%d)\n",
      clock, clock->clocktype, clock->usedefault, clock->nruns);

   if( clock->nruns == 0 )
   {
      /* the clock is not running: convert the clocks timer into seconds */
      switch( clock->clocktype )
      {
      case SCIP_CLOCKTYPE_DEFAULT:
         return 0.0;

      case SCIP_CLOCKTYPE_CPU:
         return cputime2sec(clock->data.cpuclock.user);

      case SCIP_CLOCKTYPE_WALL:            
         return walltime2sec(clock->data.wallclock.sec, clock->data.wallclock.usec);

      default:
         errorMessage("invalid clock type\n");
         abort();
      }
   }
   else
   {
      struct timeval tp; /*lint !e86*/
      struct tms now;
         
      /* the clock is currently running: we have to add the current time to the clocks timer */
      switch( clock->clocktype )
      {
      case SCIP_CLOCKTYPE_CPU:
         (void)times(&now);
         return cputime2sec(clock->data.cpuclock.user + now.tms_utime);

      case SCIP_CLOCKTYPE_WALL:            
         gettimeofday(&tp, NULL);
         if( tp.tv_usec + clock->data.wallclock.usec > 1000000 ) /*lint !e115 !e40*/
            return walltime2sec(clock->data.wallclock.sec + tp.tv_sec + 1, /*lint !e115 !e40*/
               (clock->data.wallclock.usec - 1000000) + tp.tv_usec); /*lint !e115 !e40*/
         else
            return walltime2sec(clock->data.wallclock.sec + tp.tv_sec, /*lint !e115 !e40*/
               clock->data.wallclock.usec + tp.tv_usec); /*lint !e115 !e40*/

      case SCIP_CLOCKTYPE_DEFAULT:
      default:
         errorMessage("invalid clock type\n");
         abort();
      }
   }
}

/** sets the used time of this clock in seconds */
void SCIPclockSetTime(
   CLOCK*           clock,              /**< clock timer */
   Real             sec                 /**< time in seconds to set the clock's timer to */
   )
{
   assert(clock != NULL);

   debugMessage("setting time of clock %p (type %d, usedefault=%d, nruns=%d) to %g\n",
      clock, clock->clocktype, clock->usedefault, clock->nruns, sec);

   /* if the clock type is not yet set, set it to an arbitrary value to be able to store the number */
   if( clock->clocktype == SCIP_CLOCKTYPE_DEFAULT )
      clockSetType(clock, SCIP_CLOCKTYPE_WALL);

   switch( clock->clocktype )
   {
   case SCIP_CLOCKTYPE_CPU:
      sec2cputime(sec, &clock->data.cpuclock.user);
      break;
      
   case SCIP_CLOCKTYPE_WALL:            
      sec2walltime(sec, &clock->data.wallclock.sec, &clock->data.wallclock.usec);
      break;
      
   case SCIP_CLOCKTYPE_DEFAULT:
   default:
      errorMessage("invalid clock type\n");
      abort();
   }
   
   if( clock->nruns >= 1 )
   {
      struct timeval tp; /*lint !e86*/
      struct tms now;
         
      /* the clock is currently running: we have to subtract the current time from the new timer value */
      switch( clock->clocktype )
      {
      case SCIP_CLOCKTYPE_CPU:
         (void)times(&now);
         clock->data.cpuclock.user -= now.tms_utime;
         break;

      case SCIP_CLOCKTYPE_WALL:            
         gettimeofday(&tp, NULL);
         if( tp.tv_usec > clock->data.wallclock.usec ) /*lint !e115 !e40*/
         {
            clock->data.wallclock.sec -= (tp.tv_sec + 1); /*lint !e115 !e40*/
            clock->data.wallclock.usec += (1000000 - tp.tv_usec); /*lint !e115 !e40*/
         }
         else
         {
            clock->data.wallclock.sec -= tp.tv_sec; /*lint !e115 !e40*/
            clock->data.wallclock.usec -= tp.tv_usec; /*lint !e115 !e40*/
         }
         break;

      case SCIP_CLOCKTYPE_DEFAULT:
      default:
         errorMessage("invalid clock type\n");
         abort();
      }
   }
}

/** gets current time of day in seconds (standard time zone) */
Real SCIPclockGetTimeOfDay(
   void
   )
{
   struct timeval tp; /*lint !e86*/
   
   gettimeofday(&tp, NULL);

   return (Real)(tp.tv_sec % (24*3600)) + (Real)tp.tv_usec / 1e+6;
}
