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

/**@file   event_solustore.h
 * @ingroup EVENTS 
 * @brief  eventhdlr for solustore event
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_SOLUSTORE_H__
#define __SCIP_EVENT_SOLUSTORE_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for solustore event */
extern
SCIP_RETCODE SCIPincludeEventHdlrSolustore(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
void SCIPeventhdlrSolustoreDisable(
   SCIP*                 scip
   );

#ifdef __cplusplus
}
#endif

#endif
