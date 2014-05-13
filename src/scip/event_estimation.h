/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_estimation.h
 * @ingroup EVENTS
 * @brief  eventhdlr for estimation event
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_ESTIMATION_H__
#define __SCIP_EVENT_ESTIMATION_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

EXTERN
SCIP_Real SCIPgetRootLPSolPscostEstimate(
   SCIP*                 scip
   );

SCIP_RETCODE SCIPgetCorrectedEstimateData(
   SCIP*                 scip,
   SCIP_Real*            mincorrectedestimate,
   SCIP_Real*            rootcorrectedestim,
   SCIP_Real*            minestimate,
   int*                  nnodesbelowincumbentcorrected,
   int*                  nnodesbelowincumbent,
   SCIP_Bool             recalcestim
);

/** creates event handler for estimation event */
EXTERN
SCIP_RETCODE SCIPincludeEventHdlrEstimation(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
