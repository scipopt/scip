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
 * @brief  methods and datastructures for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __EVENT_H__
#define __EVENT_H__


/*
 * event types
 */

#define SCIP_EVENTTYPE_DISABLED       0x00000000 /**< the event was disabled and has no effect any longer */

/* variable events */
#define SCIP_EVENTTYPE_VARCREATED     0x00000001 /**< ??? TODO: a variable has been created */
#define SCIP_EVENTTYPE_VARFIXED       0x00000002 /**< ??? TODO: a variable has been fixed, aggregated, or multiaggregated */
#define SCIP_EVENTTYPE_OBJCHANGED     0x00000004 /**< the objective value of a variable has been changed */
#define SCIP_EVENTTYPE_LBTIGHTENED    0x00000008 /**< the lower bound of a variable has been increased */
#define SCIP_EVENTTYPE_LBRELAXED      0x00000010 /**< the lower bound of a variable has been decreased */
#define SCIP_EVENTTYPE_UBTIGHTENED    0x00000020 /**< the upper bound of a variable has been decreased */
#define SCIP_EVENTTYPE_UBRELAXED      0x00000040 /**< the upper bound of a variable has been increased */
#define SCIP_EVENTTYPE_HOLEADDED      0x00000080 /**< ??? TODO: a hole has been added to the hole list of a variable's domain */
#define SCIP_EVENTTYPE_HOLEREMOVED    0x00000100 /**< ??? TODO: a hole has been removed from the hole list of a variable's domain */

/* node events */
#define SCIP_EVENTTYPE_NODEACTIVATED  0x00000200 /**< a node has been activated and is now the current active node */
#define SCIP_EVENTTYPE_NODEFEASIBLE   0x00000400 /**< the LP/pseudo solution of the node was feasible */
#define SCIP_EVENTTYPE_NODEINFEASIBLE 0x00000800 /**< the active node has been proven to be infeasible or was bounded */
#define SCIP_EVENTTYPE_NODEBRANCHED   0x00001000 /**< the active node has been solved by branching */

/* LP events */
#define SCIP_EVENTTYPE_FIRSTLPSOLVED  0x00002000 /**< the node's initial LP was solved */
#define SCIP_EVENTTYPE_LPSOLVED       0x00004000 /**< the node's LP was completely solved with cut & price */

/* primal solution events */
#define SCIP_EVENTTYPE_POORSOLFOUND   0x00008000 /**< a good enough primal feasible (but not new best) solution was found */
#define SCIP_EVENTTYPE_BESTSOLFOUND   0x00010000 /**< a new best primal feasible solution was found */

/* event masks for variable events */
#define SCIP_EVENTTYPE_LBCHANGED      (SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED)
#define SCIP_EVENTTYPE_UBCHANGED      (SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED)
#define SCIP_EVENTTYPE_BOUNDTIGHTENED (SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBTIGHTENED)
#define SCIP_EVENTTYPE_BOUNDRELAXED   (SCIP_EVENTTYPE_LBRELAXED | SCIP_EVENTTYPE_UBRELAXED)
#define SCIP_EVENTTYPE_BOUNDCHANGED   (SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED)
#define SCIP_EVENTTYPE_HOLECHANGED    (SCIP_EVENTTYPE_HOLEADDED | SCIP_EVENTTYPE_HOLEREMOVED)
#define SCIP_EVENTTYPE_DOMCHANGED     (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_HOLECHANGED)
#define SCIP_EVENTTYPE_VARCHANGED     (SCIP_EVENTTYPE_OBJCHANGED | SCIP_EVENTTYPE_DOMCHANGED)
#define SCIP_EVENTTYPE_VAREVENT       (SCIP_EVENTTYPE_VARCREATED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARCHANGED)
   
/* event masks for node events */
#define SCIP_EVENTTYPE_NODESOLVED     (SCIP_EVENTTYPE_NODEFEASIBLE | SCIP_EVENTTYPE_NODEINFEASIBLE \
                                       | SCIP_EVENTTYPE_NODEBRANCHED)
#define SCIP_EVENTTYPE_NODEEVENT      (SCIP_EVENTTYPE_NODEACTIVATED | SCIP_EVENTTYPE_NODESOLVED)

/* event masks for LP events */
#define SCIP_EVENTTYPE_LPEVENT        (SCIP_EVENTTYPE_FIRSTLPSOLVED | SCIP_EVENTTYPE_LPSOLVED)

/* event masks for primal solution events */
#define SCIP_EVENTTYPE_SOLFOUND       (SCIP_EVENTTYPE_POORSOLFOUND | SCIP_EVENTTYPE_BESTSOLFOUND)
#define SCIP_EVENTTYPE_SOLEVENT       (SCIP_EVENTTYPE_SOLFOUND)


typedef unsigned int EVENTTYPE;         /**< type of event (bit field) */

typedef struct EventHdlr EVENTHDLR;     /**< event handler for a specific events */
typedef struct EventHdlrData EVENTHDLRDATA; /**< event handler data */
typedef struct Event EVENT;             /**< event data structure */
typedef struct EventObjChg EVENTOBJCHG; /**< data for objective value change events */
typedef struct EventBdChg EVENTBDCHG;   /**< data for bound change events */
typedef struct EventData EVENTDATA;     /**< locally defined event specific data */
typedef struct EventFilter EVENTFILTER; /**< event filter to select events to be processed by an event handler */
typedef struct EventQueue EVENTQUEUE;   /**< event queue to cache events and process them later */



/** destructor of event handler to free user data (called when SCIP is exiting)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    eventhdlr       : the event handler itself
 */
#define DECL_EVENTFREE(x) RETCODE x (SCIP* scip, EVENTHDLR* eventhdlr)

/** initialization method of event handler (called when problem solving starts)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    eventhdlr       : the event handler itself
 */
#define DECL_EVENTINIT(x) RETCODE x (SCIP* scip, EVENTHDLR* eventhdlr)

/** deinitialization method of event handler (called when problem solving exits)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    eventhdlr       : the event handler itself
 */
#define DECL_EVENTEXIT(x) RETCODE x (SCIP* scip, EVENTHDLR* eventhdlr)

/** frees specific event data
 *
 *  input:
 *    scip            : SCIP main data structure
 *    eventhdlr       : the event handler itself
 *    eventdata       : pointer to the event data to free
 */
#define DECL_EVENTDELETE(x) RETCODE x (SCIP* scip, EVENTHDLR* eventhdlr, EVENTDATA** eventdata)

/** execution method of event handler
 *
 *  Processes the event. The method is called every time an event occurs, for which the event handler
 *  is responsible. Event handlers may declare themselves resposible for events by calling the
 *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
 *  given event handler and event data.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    eventhdlr       : the event handler itself
 *    event           : event to process
 *    eventdata       : user data for the event
 */
#define DECL_EVENTEXEC(x) RETCODE x (SCIP* scip, EVENTHDLR* eventhdlr, EVENT* event, EVENTDATA* eventdata)




#include "scip.h"
#include "memory.h"
#include "retcode.h"
#include "var.h"
#include "tree.h"
#include "lp.h"
#include "branch.h"



/** data for objective value change events */
struct EventObjChg
{
   VAR*             var;                /**< variable whose objective value changed */
   Real             oldobj;             /**< old objective value before value changed */
   Real             newobj;             /**< new objective value after value changed */
};

/** data for bound change events */
struct EventBdChg
{
   VAR*             var;                /**< variable whose bound changed */
   Real             oldbound;           /**< old bound before bound changed */
   Real             newbound;           /**< new bound after bound changed */
};

/** event data structure */
struct Event
{
   union
   {
      EVENTOBJCHG   eventobjchg;        /**< data for objective value change events */
      EVENTBDCHG    eventbdchg;         /**< data for bound change events */
      NODE*         node;               /**< data for node and LP events */
      SOL*          sol;                /**< data for primal solution events */
   } data;
   EVENTTYPE        eventtype;          /**< type of event */
};

/** event filter to select events to be processed by an event handler */
struct EventFilter
{
   EVENTTYPE*       eventtypes;         /**< array with types of event to process */
   EVENTHDLR**      eventhdlrs;         /**< array with event handlers to process the event */
   EVENTDATA**      eventdatas;         /**< array with user data for the issued event */
   int*             eventnuses;         /**< array with number of times, the eventhandler/data was added to the filter */
   int              size;               /**< size of filter arrays (available slots in arrays) */
   int              len;                /**< number entries in filter arrays */
   unsigned int     eventmask;          /**< mask for events that are handled by any event handler in the filter */
};



/*
 * Event handler methods
 */

/** creates an event handler */
extern
RETCODE SCIPeventhdlrCreate(
   EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   const char*      name,               /**< name of event handler */
   const char*      desc,               /**< description of event handler */
   DECL_EVENTFREE   ((*eventfree)),     /**< destructor of event handler */
   DECL_EVENTINIT   ((*eventinit)),     /**< initialize event handler */
   DECL_EVENTEXIT   ((*eventexit)),     /**< deinitialize event handler */
   DECL_EVENTDELETE ((*eventdelete)),   /**< free specific event data */
   DECL_EVENTEXEC   ((*eventexec)),     /**< execute event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   );

/** calls destructor and frees memory of event handler */
extern
RETCODE SCIPeventhdlrFree(
   EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes event handler */
extern
RETCODE SCIPeventhdlrInit(
   EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls exit method of event handler */
extern
RETCODE SCIPeventhdlrExit(
   EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls execution method of event handler */
extern
RETCODE SCIPeventhdlrExec(
   EVENTHDLR*       eventhdlr,          /**< event handler */
   const SET*       set,                /**< global SCIP settings */
   EVENT*           event,              /**< event to call event handler with */
   EVENTDATA*       eventdata           /**< user data for the issued event */
   );

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

/** creates an event for a change in the objective value of a variable */
extern
RETCODE SCIPeventCreateObjChanged(
   EVENT**          event,              /**< pointer to store the event */
   MEMHDR*          memhdr,             /**< block memory */
   VAR*             var,                /**< variable whose objective value changed */
   Real             oldobj,             /**< old objective value before value changed */
   Real             newobj              /**< new objective value after value changed */
   );

/** creates an event for a change in the lower bound of a variable */
extern
RETCODE SCIPeventCreateLbChanged(
   EVENT**          event,              /**< pointer to store the event */
   MEMHDR*          memhdr,             /**< block memory */
   VAR*             var,                /**< variable whose bound changed */
   Real             oldbound,           /**< old bound before bound changed */
   Real             newbound            /**< new bound after bound changed */
   );

/** creates an event for a change in the upper bound of a variable */
extern
RETCODE SCIPeventCreateUbChanged(
   EVENT**          event,              /**< pointer to store the event */
   MEMHDR*          memhdr,             /**< block memory */
   VAR*             var,                /**< variable whose bound changed */
   Real             oldbound,           /**< old bound before bound changed */
   Real             newbound            /**< new bound after bound changed */
   );

/** frees an event */
extern
RETCODE SCIPeventFree(
   EVENT**          event,              /**< event to free */
   MEMHDR*          memhdr              /**< block memory buffer */
   );

/** gets type of event */
extern
EVENTTYPE SCIPeventGetType(
   EVENT*           event               /**< event */
   );

/** sets type of event */
extern
RETCODE SCIPeventChgType(
   EVENT*           event,              /**< event */
   EVENTTYPE        eventtype           /**< new event type */
   );

/** gets variable for a variable change event (objective value or domain change) */
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

/** sets node for a node or LP event */
extern
RETCODE SCIPeventChgNode(
   EVENT*           event,              /**< event */
   NODE*            node                /**< new node */
   );

/** gets solution for a primal solution event */
extern
SOL* SCIPeventGetSol(
   EVENT*           event               /**< event */
   );

/** sets solution for a primal solution event */
extern
RETCODE SCIPeventChgSol(
   EVENT*           event,              /**< event */
   SOL*             sol                 /**< new primal solution */
   );

/** processes event by calling the appropriate event handlers */
extern
RETCODE SCIPeventProcess(
   EVENT*           event,              /**< event */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree; only needed for BOUNDCHANGED events */
   LP*              lp,                 /**< actual LP data; only needed for BOUNDCHANGED events */
   BRANCHCAND*      branchcand,         /**< branching candidate storage; only needed for BOUNDCHANGED events */
   EVENTFILTER*     eventfilter         /**< event filter for global events; not needed for BOUNDCHANGED events */
   );



/*
 * Event filter methods
 */

/** creates an event filter */
extern
RETCODE SCIPeventfilterCreate(
   EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   MEMHDR*          memhdr              /**< block memory buffer */
   );

/** frees an event filter and the associated event data entries */
extern
RETCODE SCIPeventfilterFree(
   EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   );

/** adds element to event filter */
extern
RETCODE SCIPeventfilterAdd(
   EVENTFILTER*     eventfilter,        /**< event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   );

/** deletes element from event filter */
extern
RETCODE SCIPeventfilterDel(
   EVENTFILTER*     eventfilter,        /**< event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   );

/** processes the event with all event handlers with matching filter setting */
extern
RETCODE SCIPeventfilterProcess(
   EVENTFILTER*     eventfilter,        /**< event filter */
   const SET*       set,                /**< global SCIP settings */
   EVENT*           event               /**< event to process */
   );



/*
 * Event queue methods
 */

/** creates an event queue */
extern
RETCODE SCIPeventqueueCreate(
   EVENTQUEUE**     eventqueue          /**< pointer to store the event queue */
   );

/** frees event queue; there must not be any unprocessed eventy in the queue! */
extern
RETCODE SCIPeventqueueFree(
   EVENTQUEUE**     eventqueue          /**< pointer to the event queue */
   );

/** processes event or adds event to the event queue */
extern
RETCODE SCIPeventqueueAdd(
   EVENTQUEUE*      eventqueue,         /**< event queue */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree; only needed for BOUNDCHANGED events */
   LP*              lp,                 /**< actual LP data; only needed for BOUNDCHANGED events */
   BRANCHCAND*      branchcand,         /**< branching candidate storage; only needed for BOUNDCHANGED events */
   EVENTFILTER*     eventfilter,        /**< event filter for global events; not needed for BOUNDCHANGED events */
   EVENT**          event               /**< pointer to event to add to the queue; will be NULL after queue addition */
   );

/** marks queue to delay incoming events until a call to SCIPeventqueueProcess() */
extern
RETCODE SCIPeventqueueDelay(
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** processes all events in the queue */
extern
RETCODE SCIPeventqueueProcess(
   EVENTQUEUE*      eventqueue,         /**< event queue */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTFILTER*     eventfilter         /**< event filter for global (not variable dependent) events */
   );

/** returns TRUE iff events of the queue are delayed until the next SCIPeventqueueProcess() call */
extern
Bool SCIPeventqueueIsDelayed(
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

#endif
