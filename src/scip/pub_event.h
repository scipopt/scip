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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_event.h,v 1.6 2005/02/14 13:35:48 bzfpfend Exp $"

/**@file   pub_event.h
 * @brief  public methods for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_EVENT_H__
#define __PUB_EVENT_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_event.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"
#include "scip/type_tree.h"



/*
 * Event handler methods
 */

/** gets name of event handler */
extern
const char* SCIPeventhdlrGetName(
   EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** gets user data of event handler */
extern
EVENTHDLRDATA* SCIPeventhdlrGetData(
   EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** sets user data of event handler; user has to free old data in advance! */
extern
void SCIPeventhdlrSetData(
   EVENTHDLR*       eventhdlr,          /**< event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< new event handler user data */
   );

/** is event handler initialized? */
extern
Bool SCIPeventhdlrIsInitialized(
   EVENTHDLR*       eventhdlr           /**< event handler */
   );



/*
 * Event methods
 */

/** gets type of event */
extern
EVENTTYPE SCIPeventGetType(
   EVENT*           event               /**< event */
   );

/** gets variable for a variable event (var added, var fixed, objective value or domain change) */
extern
VAR* SCIPeventGetVar(
   EVENT*           event               /**< event */
   );

/** gets old objective value for an objective value change event */
extern
Real SCIPeventGetOldobj(
   EVENT*           event               /**< event */
   );

/** gets new objective value for an objective value change event */
extern
Real SCIPeventGetNewobj(
   EVENT*           event               /**< event */
   );

/** gets old bound for a bound change event */
extern
Real SCIPeventGetOldbound(
   EVENT*           event               /**< event */
   );

/** gets new bound for a bound change event */
extern
Real SCIPeventGetNewbound(
   EVENT*           event               /**< event */
   );

/** gets node for a node or LP event */
extern
NODE* SCIPeventGetNode(
   EVENT*           event               /**< event */
   );

/** gets solution for a primal solution event */
extern
SOL* SCIPeventGetSol(
   EVENT*           event               /**< event */
   );


#endif
