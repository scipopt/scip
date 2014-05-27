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

/**@file   event_nodereopt.h
 * @ingroup EVENTS 
 * @brief  eventhdlr for nodereopt event
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_NODEREOPT_H__
#define __SCIP_EVENT_NODEREOPT_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for nodereopt event */
extern
SCIP_RETCODE SCIPincludeEventHdlrNodereopt(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
void SCIPeventhdlrNodereoptDisable(
   SCIP*                 scip
   );

#ifdef CCHECK
SCIP_Bool SCIPeventhdlrNodereoptCheckConsistency(
   SCIP*                scip,                /**< SCIP data structure */
   SCIP_Longint         ncreatednodes        /**< number of generated nodes */
   );
#endif

#ifdef __cplusplus
}
#endif

#endif
