/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event.h
 * @brief  datastructures and methods for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __EVENT_H__
#define __EVENT_H__

enum Eventtype
{
   SCIP_EVENTTYPE_DISABLED       = 0x00000000, /**< the event was disabled and has no effect any longer */
   SCIP_EVENTTYPE_VARCREATED     = 0x00000001, /**< a variable has been created */
   SCIP_EVENTTYPE_VARFIXED       = 0x00000002, /**< a variable has been fixed, aggregated, or multiaggregated */
   SCIP_EVENTTYPE_LBTIGHTENED    = 0x00000004, /**< the lower bound of a variable has been increased */
   SCIP_EVENTTYPE_LBRELAXED      = 0x00000008, /**< the lower bound of a variable has been decreased */
   SCIP_EVENTTYPE_UBTIGHTENED    = 0x00000010, /**< the upper bound of a variable has been decreased */
   SCIP_EVENTTYPE_UBRELAXED      = 0x00000020, /**< the upper bound of a variable has been increased */
   SCIP_EVENTTYPE_HOLEADDED      = 0x00000040, /**< a hole has been added to the hole list of a variable's domain */
   SCIP_EVENTTYPE_HOLEREMOVED    = 0x00000080, /**< a hole has been removed from the hole list of a variable's domain */
   SCIP_EVENTTYPE_NODEACTIVATED  = 0x00000100, /**< a node has been activated and is now the current active node */
   SCIP_EVENTTYPE_NODESOLVED     = 0x00000200, /**< the active node has been solved by branching or bounding */
   SCIP_EVENTTYPE_NODEBRANCHED   = 0x00000400, /**< the active node has been solved by branching */
   SCIP_EVENTTYPE_NODEINFEASIBLE = 0x00000800, /**< the active node has been proven to be infeasible */
   SCIP_EVENTTYPE_NODEFEASIBLE   = 0x00001000, /**< the LP or pseudo solution of the node was feasible */
   SCIP_EVENTTYPE_FIRSTLPSOLVED  = 0x00002000, /**< the node's initial LP was solved */
   SCIP_EVENTTYPE_LPSOLVED       = 0x00004000, /**< the node's LP was completely solved with cut & price */
   SCIP_EVENTTYPE_LBCHANGED      = SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
   SCIP_EVENTTYPE_UBCHANGED      = SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
   SCIP_EVENTTYPE_BOUNDCHANGED   = SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED,
   SCIP_EVENTTYPE_HOLECHANGED    = SCIP_EVENTTYPE_HOLEADDED | SCIP_EVENTTYPE_HOLEREMOVED,
   SCIP_EVENTTYPE_DOMCHANGED     = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_HOLECHANGED
};
typedef enum Eventtype EVENTTYPE;       /**< type of event (bit field) */

typedef struct EventHdlr EVENTHDLR;     /**< event handler for a specific events */
typedef struct EventHdlrData EVENTHDLRDATA; /**< event handler data */
typedef struct Event EVENT;             /**< event data structure */
typedef struct EventData EVENTDATA;     /**< locally defined event specific data */
typedef struct EventFilter EVENTFILTER; /**< event filter to select events to be processed by an event handler */
typedef struct EventQueue EVENTQUEUE;   /**< event queue to cache events and process them later */



/** destructor of event handler to free user data (called when SCIP is exiting)
 *
 *  input:
 *    eventhdlr       : the event handler itself
 *    scip            : SCIP main data structure
 */
#define DECL_EVENTFREE(x) RETCODE x (EVENTHDLR* eventhdlr, SCIP* scip)

/** initialization method of event handler (called when problem solving starts)
 *
 *  input:
 *    eventhdlr       : the event handler itself
 *    scip            : SCIP main data structure
 */
#define DECL_EVENTINIT(x) RETCODE x (EVENTHDLR* eventhdlr, SCIP* scip)

/** deinitialization method of event handler (called when problem solving exits)
 *
 *  input:
 *    eventhdlr       : the event handler itself
 *    scip            : SCIP main data structure
 */
#define DECL_EVENTEXIT(x) RETCODE x (EVENTHDLR* eventhdlr, SCIP* scip)

/** frees specific event data
 *
 *  input:
 *    eventhdlr       : the event handler itself
 *    scip            : SCIP main data structure
 *    eventdata       : pointer to the event data to free
 */
#define DECL_EVENTDELE(x) RETCODE x (EVENTHDLR* eventhdlr, SCIP* scip, EVENTDATA** eventdata)

/** execution method of event handler
 *
 *  Processes the event. The method is called every time an event occurs, for which the event handler
 *  is responsible. Event handlers may declare themselves resposible for events by calling the
 *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
 *  given event handler and event data.
 *
 *  input:
 *    eventhdlr       : the event handler itself
 *    scip            : SCIP main data structure
 *    event           : event to process
 *    eventdata       : user data for the event
 */
#define DECL_EVENTEXEC(x) RETCODE x (EVENTHDLR* eventhdlr, SCIP* scip, EVENT* event, EVENTDATA* eventdata)




#include "scip.h"
#include "memory.h"
#include "retcode.h"
#include "var.h"
#include "tree.h"
#include "lp.h"



/** event filter to select events to be processed by an event handler */
struct EventFilter
{
   EVENTTYPE*       eventtypes;         /**< array with types of event to process */
   EVENTHDLR**      eventhdlrs;         /**< array with event handlers to process the event */
   EVENTDATA**      eventdatas;         /**< array with user data for the issued event */
   int              size;               /**< size of filter arrays (available slots in arrays) */
   int              len;                /**< number entries in filter arrays */
};



/*
 * Event handler methods
 */

extern
RETCODE SCIPeventhdlrCreate(            /**< creates an event handler */
   EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   const char*      name,               /**< name of event handler */
   const char*      desc,               /**< description of event handler */
   DECL_EVENTFREE((*eventfree)),        /**< destructor of event handler */
   DECL_EVENTINIT((*eventinit)),        /**< initialise event handler */
   DECL_EVENTEXIT((*eventexit)),        /**< deinitialise event handler */
   DECL_EVENTDELE((*eventdele)),        /**< free specific event data */
   DECL_EVENTEXEC((*eventexec)),        /**< execute event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   );

extern
RETCODE SCIPeventhdlrFree(              /**< calls destructor and frees memory of event handler */
   EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPeventhdlrInit(              /**< initializes event handler */
   EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPeventhdlrExit(              /**< calls exit method of event handler */
   EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPeventhdlrExec(              /**< calls execution method of event handler */
   EVENTHDLR*       eventhdlr,          /**< event handler */
   const SET*       set,                /**< global SCIP settings */
   EVENT*           event,              /**< event to call event handler with */
   EVENTDATA*       eventdata           /**< user data for the issued event */
   );

extern
const char* SCIPeventhdlrGetName(       /**< gets name of event handler */
   EVENTHDLR*       eventhdlr           /**< event handler */
   );

extern
EVENTHDLRDATA* SCIPeventhdlrGetData(    /**< gets user data of event handler */
   EVENTHDLR*       eventhdlr           /**< event handler */
   );

extern
void SCIPeventhdlrSetData(              /**< sets user data of event handler; user has to free old data in advance! */
   EVENTHDLR*       eventhdlr,          /**< event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< new event handler user data */
   );

extern
Bool SCIPeventhdlrIsInitialized(        /**< is event handler initialized? */
   EVENTHDLR*       eventhdlr           /**< event handler */
   );



/*
 * Event methods
 */

extern
RETCODE SCIPeventCreateLbChanged(       /**< creates an event for a change in the lower bound of a variable */
   EVENT**          event,              /**< pointer to store the event */
   MEMHDR*          memhdr,             /**< block memory */
   VAR*             var,                /**< variable whose bound changed */
   Real             oldbound,           /**< old bound before bound changed */
   Real             newbound            /**< new bound after bound changed */
   );

extern
RETCODE SCIPeventCreateUbChanged(       /**< creates an event for a change in the upper bound of a variable */
   EVENT**          event,              /**< pointer to store the event */
   MEMHDR*          memhdr,             /**< block memory */
   VAR*             var,                /**< variable whose bound changed */
   Real             oldbound,           /**< old bound before bound changed */
   Real             newbound            /**< new bound after bound changed */
   );

extern
RETCODE SCIPeventFree(                  /**< frees an event */
   EVENT**          event,              /**< event to free */
   MEMHDR*          memhdr              /**< block memory buffer */
   );

extern
RETCODE SCIPeventGetType(               /**< gets type of event */
   EVENT*           event,              /**< event */
   EVENTTYPE*       eventtype           /**< pointer to store the event type */
   );

extern
RETCODE SCIPeventGetVar(                /**< gets variable for a domain change event */
   EVENT*           event,              /**< event */
   VAR**            var                 /**< pointer to store the variable */
   );

extern
RETCODE SCIPeventGetOldbound(           /**< gets old bound for a bound change event */
   EVENT*           event,              /**< event */
   Real*            bound               /**< pointer to store the bound */
   );

extern
RETCODE SCIPeventGetNewbound(           /**< gets new bound for a bound change event */
   EVENT*           event,              /**< event */
   Real*            bound               /**< pointer to store the bound */
   );

extern
RETCODE SCIPeventProcess(               /**< processes event by calling the appropriate event handlers */
   EVENT*           event,              /**< event */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   );



/*
 * Event filter methods
 */

extern
RETCODE SCIPeventfilterCreate(          /**< creates an event filter */
   EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   MEMHDR*          memhdr              /**< block memory buffer */
   );

extern
RETCODE SCIPeventfilterFree(            /**< frees an event filter and the associated event data entries */
   EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPeventfilterAdd(             /**< adds element to event filter */
   EVENTFILTER*     eventfilter,        /**< event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   );

extern
RETCODE SCIPeventfilterDel(             /**< deletes element from event filter */
   EVENTFILTER*     eventfilter,        /**< event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   );

extern
RETCODE SCIPeventfilterProcess(         /**< processes the event with all event handlers with matching filter setting */
   EVENTFILTER*     eventfilter,        /**< event filter */
   const SET*       set,                /**< global SCIP settings */
   EVENT*           event               /**< event to process */
   );



/*
 * Event queue methods
 */

extern
RETCODE SCIPeventqueueCreate(           /**< creates an event queue */
   EVENTQUEUE**     eventqueue          /**< pointer to store the event queue */
   );

extern
RETCODE SCIPeventqueueFree(             /**< frees event queue; there must not be any unprocessed eventy in the queue! */
   EVENTQUEUE**     eventqueue          /**< pointer to the event queue */
   );

extern
RETCODE SCIPeventqueueAdd(              /**< processes event or adds event to the event queue */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENT**          event               /**< pointer to event to add to the queue; will be NULL after queue addition */
   );

extern
RETCODE SCIPeventqueueDelay(            /**< marks queue to delay incoming events until a call to SCIPeventqueueProcess() */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

extern
RETCODE SCIPeventqueueProcess(          /**< processes all events in the queue */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   );


#endif
