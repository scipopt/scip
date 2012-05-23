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

/**@file   pub_event.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_EVENT_H__
#define __SCIP_PUB_EVENT_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_event.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"
#include "scip/type_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Event handler methods
 */

/** gets name of event handler */
extern
const char* SCIPeventhdlrGetName(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** gets user data of event handler */
extern
SCIP_EVENTHDLRDATA* SCIPeventhdlrGetData(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** sets user data of event handler; user has to free old data in advance! */
extern
void SCIPeventhdlrSetData(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< new event handler user data */
   );

/** sets copy callback for all events of this event handler */
extern
void SCIPeventhdlrSetCopy(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTCOPY   ((*eventcopy))      /**< copy callback for events */
   );

/** sets destructor callback of this event handler */
extern
void SCIPeventhdlrSetFree(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTFREE   ((*eventfree))      /**< destructor callback of event handler */
   );

/** sets initialization callback of this event handler */
extern
void SCIPeventhdlrSetInit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit))      /**< initialization callback of event handler */
   );

/** sets deinitialization callback of this event handler */
extern
void SCIPeventhdlrSetExit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit))      /**< deinitialization callback of event handler */
   );

/** sets solving process initialization callback of this event handler */
extern
void SCIPeventhdlrSetInitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol))   /**< solving process initialization callback of event handler */
   );

/** sets solving process deinitialization callback of this event handler */
extern
void SCIPeventhdlrSetExitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol))   /**< solving process deinitialization callback of event handler */
   );

/** sets callback to free specific event data */
extern
void SCIPeventhdlrSetDelete(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete))    /**< callback to free specific event data */
   );

/** is event handler initialized? */
extern
SCIP_Bool SCIPeventhdlrIsInitialized(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** gets time in seconds used in this event handler for setting up for next stages */
extern
SCIP_Real SCIPeventhdlrGetSetupTime(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** gets time in seconds used in this event handler */
extern
SCIP_Real SCIPeventhdlrGetTime(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );


/*
 * Event methods
 */

/** gets type of event */
extern
SCIP_EVENTTYPE SCIPeventGetType(
   SCIP_EVENT*           event               /**< event */
   );

/** gets variable for a variable event (var added, var deleted, var fixed, 
 *  objective value or domain change, domain hole added or removed) */
extern
SCIP_VAR* SCIPeventGetVar(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old objective value for an objective value change event */
extern
SCIP_Real SCIPeventGetOldobj(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new objective value for an objective value change event */
extern
SCIP_Real SCIPeventGetNewobj(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old bound for a bound change event */
extern
SCIP_Real SCIPeventGetOldbound(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new bound for a bound change event */
extern
SCIP_Real SCIPeventGetNewbound(
   SCIP_EVENT*           event               /**< event */
   );

/** gets node for a node or LP event */
extern
SCIP_NODE* SCIPeventGetNode(
   SCIP_EVENT*           event               /**< event */
   );

/** gets solution for a primal solution event */
extern
SCIP_SOL* SCIPeventGetSol(
   SCIP_EVENT*           event               /**< event */
   );

/** gets the left bound of open interval in the hole */
extern
SCIP_Real SCIPeventGetHoleLeft(
   SCIP_EVENT*           event               /**< event */
   );

/** gets the right bound of open interval in the hole */
extern
SCIP_Real SCIPeventGetHoleRight(
   SCIP_EVENT*           event               /**< event */
   );

/** gets row for a row event */
extern
SCIP_ROW* SCIPeventGetRow(
   SCIP_EVENT*           event               /**< event */
   );

/** gets column for a row change coefficient event */
extern
SCIP_COL* SCIPeventGetRowCol(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old coefficient value for a row change coefficient event */
extern
SCIP_Real SCIPeventGetRowOldCoefVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new coefficient value for a row change coefficient event */
extern
SCIP_Real SCIPeventGetRowNewCoefVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old constant value for a row change constant event */
extern
SCIP_Real SCIPeventGetRowOldConstVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new constant value for a row change constant event */
extern
SCIP_Real SCIPeventGetRowNewConstVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets side for a row change side event */
extern
SCIP_SIDETYPE SCIPeventGetRowSide(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old side value for a row change side event */
extern
SCIP_Real SCIPeventGetRowOldSideVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new side value for a row change side event */
extern
SCIP_Real SCIPeventGetRowNewSideVal(
   SCIP_EVENT*           event               /**< event */
   );

#ifdef __cplusplus
}
#endif

#endif
