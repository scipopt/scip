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

/**@file   event_treeinfos.h
 * @ingroup EVENTS
 * @brief  eventhdlr for treeinfos event
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_TREEINFOS_H__
#define __SCIP_EVENT_TREEINFOS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for treeinfos event */
EXTERN
SCIP_RETCODE SCIPincludeEventHdlrTreeinfos(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the current number of rank 1 nodes in the tree */
EXTERN
int SCIPgetNRank1Nodes(
   SCIP* scip
   );

#ifdef __cplusplus
}
#endif

#endif
