/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_conflict.h,v 1.1 2003/12/01 14:41:29 bzfpfend Exp $"

/**@file   pub_conflict.h
 * @brief  public methods for conflict analysis handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_CONFLICT_H__
#define __PUB_CONFLICT_H__



#include "def.h"
#include "type_conflict.h"



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
