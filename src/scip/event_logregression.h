/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_logregression.h
 * @ingroup EVENTS
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_LOGREGRESSION_H__
#define __SCIP_EVENT_LOGREGRESSION_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates event handler for solving stage event */
extern
SCIP_RETCODE SCIPincludeEventHdlrLogregression(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get axis intercept of current tangent to logarithmic regression curve */
extern
SCIP_Real getCurrentRegressionTangentAxisIntercept(
  SCIP*                 scip,
  SCIP_EVENTHDLR*       eventhdlr
);

#ifdef __cplusplus
}
#endif

#endif
