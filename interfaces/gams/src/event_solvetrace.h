/* Copyright (C) GAMS Development and others 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

/**@file   event_bbtrace.h
 * @brief  event handler to write GAMS solve trace file
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_SOLVETRACE_H__
#define __SCIP_EVENT_SOLVETRACE_H__

#include "scip/scip.h"
#include "gmomcc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for solve trace event */
extern
SCIP_RETCODE SCIPincludeEventHdlrSolveTrace(
   SCIP*                 scip,               /**< SCIP data structure */
   gmoHandle_t           gmo                 /**< GAMS model object */
   );

#ifdef __cplusplus
}
#endif

#endif
