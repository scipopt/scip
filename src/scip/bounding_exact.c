/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bounding_exact.h
 * @brief  safe exact rational bounding methods
 * @author Leon Eifler
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BOUNDING_EXACT_C__
#define __SCIP_BOUNDING_EXACT_C__


#include <stdio.h>
#include <assert.h>

#include "scip/bounding_exact.h"
#include "scip/struct_set.h"

static
SCIP_RETCODE boundShiftRational(
   void
)
{
   return SCIP_OKAY;
}

static
SCIP_RETCODE boundShiftInterval(
   void
)
{
   return SCIP_OKAY;
}

static
SCIP_RETCODE boundShiftFP(
   void
)
{
   return SCIP_OKAY;
}

static
SCIP_RETCODE projectShiftInterval(
   void
)
{
   return SCIP_OKAY;
}

static
SCIP_RETCODE projectShiftRational(
   void
)
{
   return SCIP_OKAY;
}

static
SCIP_RETCODE basisVerification(
   void
)
{
   return SCIP_OKAY;
}

static
SCIP_RETCODE solveLpExact(
   void
)
{
   return SCIP_OKAY;
}

SCIP_RETCODE computeSafeBound(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEX*            lpex,               /**< Exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Longint          itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Real*            safebound
   )
{
   /* If we are not in exact solving mode, just return */
   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   assert(set->misc_exactsolve);

   /* Not implemented yet */

   /* choose which bounding method should be calles and return a safe objective bound */
   return SCIP_OKAY;
}

#endif