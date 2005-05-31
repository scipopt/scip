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
#pragma ident "@(#) $Id: pub_conflict.h,v 1.7 2005/05/31 17:20:18 bzfpfend Exp $"

/**@file   pub_conflict.h
 * @brief  public methods for conflict analysis handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_CONFLICT_H__
#define __PUB_CONFLICT_H__



#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_conflict.h"



/** compares two conflict handlers w. r. to their priority */
extern
DECL_SORTPTRCOMP(SCIPconflicthdlrComp);

/** gets user data of conflict handler */
extern
CONFLICTHDLRDATA* SCIPconflicthdlrGetData(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** sets user data of conflict handler; user has to free old data in advance! */
extern
void SCIPconflicthdlrSetData(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata   /**< new conflict handler user data */
   );

/** gets name of conflict handler */
extern
const char* SCIPconflicthdlrGetName(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets description of conflict handler */
extern
const char* SCIPconflicthdlrGetDesc(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets priority of conflict handler */
extern
int SCIPconflicthdlrGetPriority(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** is conflict handler initialized? */
extern
Bool SCIPconflicthdlrIsInitialized(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );


#endif
