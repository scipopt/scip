/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_conflict.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for conflict analysis handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_CONFLICT_H__
#define __SCIP_PUB_CONFLICT_H__



#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_conflict.h"

#ifdef __cplusplus
extern "C" {
#endif

/** compares two conflict handlers w. r. to their priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPconflicthdlrComp);

/** gets user data of conflict handler */
extern
SCIP_CONFLICTHDLRDATA* SCIPconflicthdlrGetData(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** sets user data of conflict handler; user has to free old data in advance! */
extern
void SCIPconflicthdlrSetData(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< new conflict handler user data */
   );

/** gets name of conflict handler */
extern
const char* SCIPconflicthdlrGetName(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets description of conflict handler */
extern
const char* SCIPconflicthdlrGetDesc(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets priority of conflict handler */
extern
int SCIPconflicthdlrGetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** is conflict handler initialized? */
extern
SCIP_Bool SCIPconflicthdlrIsInitialized(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets time in seconds used in this conflict handler for setting up for next stages */
extern
SCIP_Real SCIPconflicthdlrGetSetupTime(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets time in seconds used in this conflict handler */
extern
SCIP_Real SCIPconflicthdlrGetTime(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

#ifdef __cplusplus
}
#endif

#endif
