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
#pragma ident "@(#) $Id: struct_clock.h,v 1.5 2005/01/21 09:17:07 bzfpfend Exp $"

/**@file   struct_clock.h
 * @brief  datastructures for clocks and timing issues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_CLOCK_H__
#define __STRUCT_CLOCK_H__


#include <sys/times.h>

#include "def.h"
#include "type_clock.h"



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
   int              nruns;              /**< number of SCIPclockStart() calls without SCIPclockStop() calls */
   CLOCKTYPE        clocktype;          /**< current type of clock used */
   Bool             usedefault;         /**< should the clock's type be overruled by the default clock type? */
   Bool             enabled;            /**< should the clock be used? */
};


#endif
