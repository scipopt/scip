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

/**@file   event_treeprofile.h
 * @ingroup EVENTS 
 * @brief  eventhdlr for treeprofile event
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_TREEPROFILE_H__
#define __SCIP_EVENT_TREEPROFILE_H__


#include "scip/type_scip.h"
#include "scip/type_event.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for treeprofile event */
EXTERN
SCIP_RETCODE SCIPincludeEventHdlrTreeprofile(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns a prediction of the total size of the final B&B tree, or -1 if no suitable prediction can be computed (yet) */
EXTERN
SCIP_Real SCIPpredictTotalSizeTreeprofile(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
