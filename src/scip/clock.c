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
#pragma ident "@(#) $Id: clock.c,v 1.15 2005/01/18 14:34:27 bzfpfend Exp $"

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
   CLOCK*           clck,               /**< clock timer */
   CLOCKTYPE        newtype             /**< new clock type */
   )
{
   assert(clck != NULL);
   assert(newtype != SCIP_CLOCKTYPE_DEFAULT);

   if( clck->clocktype != newtype )
   {
      if( clck->clocktype == SCIP_CLOCKTYPE_DEFAULT )
      {
         assert(clck->nruns == 0);
         clck->clocktype = newtype;
         SCIPclockReset(clck);
         debugMessage("switched clock type to %d\n", newtype);
      }
      else
      {
         Real sec;
         
         sec = SCIPclockGetTime(clck);
         clck->clocktype = newtype;
         SCIPclockSetTime(clck, sec);
         debugMessage("switched clock type to %d (%g seconds -> %g seconds)\n", newtype, sec, SCIPclockGetTime(clck));
      }
   }
}

/** if the clock uses the default clock type and the default changed, converts the clock timer to the new type */
static
void clockUpdateDefaultType(
   CLOCK*           clck,               /**< clock timer */
   CLOCKTYPE        defaultclocktype    /**< default type of clock to use */
   )
{
   assert(clck != NULL);
   assert(defaultclocktype != SCIP_CLOCKTYPE_DEFAULT);

   if( clck->usedefault && clck->clocktype != defaultclocktype )
      clockSetType(clck, defaultclocktype);
}

/** creates a clock and initializes it */
RETCODE SCIPclockCreate(
   CLOCK**          clck,               /**< pointer to clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert(clck != NULL);

   ALLOC_OKAY( allocMemory(clck) );

   SCIPclockInit(*clck, clocktype);

   return SCIP_OKAY;
}

/** frees a clock */
void SCIPclockFree(
   CLOCK**          clck                /**< pointer to clock timer */
   )
{
   assert(clck != NULL);
   
   freeMemory(clck);
}

/** initializes and resets a clock */
void SCIPclockInit(
   CLOCK*           clck,               /**< clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert(clck != NULL);

   debugMessage("initializing clock %p of type %d\n", clck, clocktype);
   clck->enabled = TRUE;
   SCIPclockSetType(clck, clocktype);
}

/** completely stop the clock and reset the clock's counter to zero */
void SCIPclockReset(
   CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   debugMessage("resetting clock %p of type %d (usedefault=%d)\n", clck, clck->clocktype, clck->usedefault);
   switch( clck->clocktype )
   {
   case SCIP_CLOCKTYPE_DEFAULT:
      break;
   case SCIP_CLOCKTYPE_CPU:
      clck->data.cpuclock.user = 0;
      break;
   case SCIP_CLOCKTYPE_WALL:
      clck->data.wallclock.sec = 0;
      clck->data.wallclock.usec = 0;
      break;
   default:
      errorMessage("invalid clock type\n");
      abort();
   }
   clck->nruns = 0;
}

/** enables the clock */
void SCIPclockEnable(
   CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   debugMessage("enabling clock %p of type %d (usedefault=%d)\n", clck, clck->clocktype, clck->usedefault);

   clck->enabled = TRUE;
}

/** disables and resets the clock */
void SCIPclockDisable(
   CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   debugMessage("disabling clock %p of type %d (usedefault=%d)\n", clck, clck->clocktype, clck->usedefault);

   clck->enabled = FALSE;
   SCIPclockReset(clck);
}

/** sets the type of the clock, overriding the default clock type, and resets the clock */
void SCIPclockSetType(
   CLOCK*           clck,               /**< clock timer */
   CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert(clck != NULL);

   debugMessage("setting type of clock %p (type %d, usedefault=%d) to %d\n", 
      clck, clck->clocktype, clck->usedefault, clocktype);

   clck->clocktype = clocktype;
   clck->usedefault = (clocktype == SCIP_CLOCKTYPE_DEFAULT);
   SCIPclockReset(clck);
}

/** starts measurement of time in the given clock */
void SCIPclockStart(
   CLOCK*           clck,               /**< clock timer */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(clck != NULL);
   assert(set != NULL);

   if( set->time_enabled && clck->enabled )
   {
      clockUpdateDefaultType(clck, set->time_clocktype);

      if( clck->nruns == 0 )
      {
         struct timeval tp; /*lint !e86*/
         struct tms now;
         
         debugMessage("starting clock %p (type %d, usedefault=%d)\n", clck, clck->clocktype, clck->usedefault);

         switch( clck->clocktype )
         {
         case SCIP_CLOCKTYPE_CPU:
            (void)times(&now);
            clck->data.cpuclock.user -= now.tms_utime;
            break;

         case SCIP_CLOCKTYPE_WALL:            
            gettimeofday(&tp, NULL);
            if( tp.tv_usec > clck->data.wallclock.usec ) /*lint !e115 !e40*/
            {
               clck->data.wallclock.sec -= (tp.tv_sec + 1); /*lint !e115 !e40*/
               clck->data.wallclock.usec += (1000000 - tp.tv_usec); /*lint !e115 !e40*/
            }
            else
            {
               clck->data.wallclock.sec -= tp.tv_sec; /*lint !e115 !e40*/
               clck->data.wallclock.usec -= tp.tv_usec; /*lint !e115 !e40*/
            }
            break;

         case SCIP_CLOCKTYPE_DEFAULT:            
         default:
            errorMessage("invalid clock type\n");
            abort();
         }
      }
      clck->nruns++;
   }
}

/** stops measurement of time in the given clock */
void SCIPclockStop(
   CLOCK*           clck,               /**< clock timer */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(clck != NULL);
   assert(set != NULL);

   if( set->time_enabled && clck->enabled )
   {
      assert(clck->nruns >= 1);

      clck->nruns--;
      if( clck->nruns == 0 )
      {
         struct timeval tp; /*lint !e86*/
         struct tms now;
         
         debugMessage("stopping clock %p (type %d, usedefault=%d)\n", clck, clck->clocktype, clck->usedefault);

         switch( clck->clocktype )
         {
         case SCIP_CLOCKTYPE_CPU:
            (void)times(&now);
            clck->data.cpuclock.user += now.tms_utime;
            break;

         case SCIP_CLOCKTYPE_WALL:            
            gettimeofday(&tp, NULL);
            if( tp.tv_usec + clck->data.wallclock.usec > 1000000 ) /*lint !e115 !e40*/
            {
               clck->data.wallclock.sec += (tp.tv_sec + 1); /*lint !e115 !e40*/
               clck->data.wallclock.usec -= (1000000 - tp.tv_usec); /*lint !e115 !e40*/
            }
            else
            {
               clck->data.wallclock.sec += tp.tv_sec; /*lint !e115 !e40*/
               clck->data.wallclock.usec += tp.tv_usec; /*lint !e115 !e40*/
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
   CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   return (clck->nruns > 0);
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
   CLOCK*           clck                /**< clock timer */
   )
{
   assert(clck != NULL);

   debugMessage("getting time of clock %p (type %d, usedefault=%d, nruns=%d)\n",
      clck, clck->clocktype, clck->usedefault, clck->nruns);

   if( clck->nruns == 0 )
   {
      /* the clock is not running: convert the clocks timer into seconds */
      switch( clck->clocktype )
      {
      case SCIP_CLOCKTYPE_DEFAULT:
         return 0.0;

      case SCIP_CLOCKTYPE_CPU:
         return cputime2sec(clck->data.cpuclock.user);

      case SCIP_CLOCKTYPE_WALL:            
         return walltime2sec(clck->data.wallclock.sec, clck->data.wallclock.usec);

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
      switch( clck->clocktype )
      {
      case SCIP_CLOCKTYPE_CPU:
         (void)times(&now);
         return cputime2sec(clck->data.cpuclock.user + now.tms_utime);

      case SCIP_CLOCKTYPE_WALL:            
         gettimeofday(&tp, NULL);
         if( tp.tv_usec + clck->data.wallclock.usec > 1000000 ) /*lint !e115 !e40*/
            return walltime2sec(clck->data.wallclock.sec + tp.tv_sec + 1, /*lint !e115 !e40*/
               (clck->data.wallclock.usec - 1000000) + tp.tv_usec); /*lint !e115 !e40*/
         else
            return walltime2sec(clck->data.wallclock.sec + tp.tv_sec, /*lint !e115 !e40*/
               clck->data.wallclock.usec + tp.tv_usec); /*lint !e115 !e40*/

      case SCIP_CLOCKTYPE_DEFAULT:
      default:
         errorMessage("invalid clock type\n");
         abort();
      }
   }
}

/** sets the used time of this clock in seconds */
void SCIPclockSetTime(
   CLOCK*           clck,               /**< clock timer */
   Real             sec                 /**< time in seconds to set the clock's timer to */
   )
{
   assert(clck != NULL);

   debugMessage("setting time of clock %p (type %d, usedefault=%d, nruns=%d) to %g\n",
      clck, clck->clocktype, clck->usedefault, clck->nruns, sec);

   /* if the clock type is not yet set, set it to an arbitrary value to be able to store the number */
   if( clck->clocktype == SCIP_CLOCKTYPE_DEFAULT )
      clockSetType(clck, SCIP_CLOCKTYPE_WALL);

   switch( clck->clocktype )
   {
   case SCIP_CLOCKTYPE_CPU:
      sec2cputime(sec, &clck->data.cpuclock.user);
      break;
      
   case SCIP_CLOCKTYPE_WALL:            
      sec2walltime(sec, &clck->data.wallclock.sec, &clck->data.wallclock.usec);
      break;
      
   case SCIP_CLOCKTYPE_DEFAULT:
   default:
      errorMessage("invalid clock type\n");
      abort();
   }
   
   if( clck->nruns >= 1 )
   {
      struct timeval tp; /*lint !e86*/
      struct tms now;
         
      /* the clock is currently running: we have to subtract the current time from the new timer value */
      switch( clck->clocktype )
      {
      case SCIP_CLOCKTYPE_CPU:
         (void)times(&now);
         clck->data.cpuclock.user -= now.tms_utime;
         break;

      case SCIP_CLOCKTYPE_WALL:            
         gettimeofday(&tp, NULL);
         if( tp.tv_usec > clck->data.wallclock.usec ) /*lint !e115 !e40*/
         {
            clck->data.wallclock.sec -= (tp.tv_sec + 1); /*lint !e115 !e40*/
            clck->data.wallclock.usec += (1000000 - tp.tv_usec); /*lint !e115 !e40*/
         }
         else
         {
            clck->data.wallclock.sec -= tp.tv_sec; /*lint !e115 !e40*/
            clck->data.wallclock.usec -= tp.tv_usec; /*lint !e115 !e40*/
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
