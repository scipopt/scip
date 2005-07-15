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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_clock.h,v 1.7 2005/07/15 17:20:22 bzfpfend Exp $"

/**@file   type_clock.h
 * @brief  type definitions for clocks and timing issues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_CLOCK_H__
#define __SCIP_TYPE_CLOCK_H__

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


#include "scip/def.h"


#endif
