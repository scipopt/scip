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
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_event.h,v 1.11 2005/01/18 09:26:58 bzfpfend Exp $"

/**@file   type_event.h
 * @brief  type definitions for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_EVENT_H__
#define __TYPE_EVENT_H__


/*
 * event types
 */

#define SCIP_EVENTTYPE_DISABLED       0x00000000 /**< the event was disabled and has no effect any longer */

/* variable events */
#define SCIP_EVENTTYPE_VARADDED       0x00000001 /**< a variable has been added to the transformed problem */
#define SCIP_EVENTTYPE_VARFIXED       0x00000002 /**< a variable has been fixed, aggregated, or multiaggregated */
#define SCIP_EVENTTYPE_LOCKSCHANGED   0x00000004 /**< the number of rounding locks of a variable changed */
#define SCIP_EVENTTYPE_OBJCHANGED     0x00000008 /**< the objective value of a variable has been changed */
#define SCIP_EVENTTYPE_LBTIGHTENED    0x00000010 /**< the lower bound of a variable has been increased */
#define SCIP_EVENTTYPE_LBRELAXED      0x00000020 /**< the lower bound of a variable has been decreased */
#define SCIP_EVENTTYPE_UBTIGHTENED    0x00000040 /**< the upper bound of a variable has been decreased */
#define SCIP_EVENTTYPE_UBRELAXED      0x00000080 /**< the upper bound of a variable has been increased */
#define SCIP_EVENTTYPE_HOLEADDED      0x00000100 /**< ??? TODO: a hole has been added to the hole list of a variable's domain */
#define SCIP_EVENTTYPE_HOLEREMOVED    0x00000200 /**< ??? TODO: a hole has been removed from the hole list of a variable's domain */

/* node events */
#define SCIP_EVENTTYPE_NODEFOCUSED    0x00000400 /**< a node has been focused and is now the focus node */
#define SCIP_EVENTTYPE_NODEFEASIBLE   0x00000800 /**< the LP/pseudo solution of the node was feasible */
#define SCIP_EVENTTYPE_NODEINFEASIBLE 0x00001000 /**< the focus node has been proven to be infeasible or was bounded */
#define SCIP_EVENTTYPE_NODEBRANCHED   0x00002000 /**< the focus node has been solved by branching */

/* LP events */
#define SCIP_EVENTTYPE_FIRSTLPSOLVED  0x00004000 /**< the node's initial LP was solved */
#define SCIP_EVENTTYPE_LPSOLVED       0x00008000 /**< the node's LP was completely solved with cut & price */

/* primal solution events */
#define SCIP_EVENTTYPE_POORSOLFOUND   0x00010000 /**< a good enough primal feasible (but not new best) solution was found */
#define SCIP_EVENTTYPE_BESTSOLFOUND   0x00020000 /**< a new best primal feasible solution was found */

/* event masks for variable events */
#define SCIP_EVENTTYPE_LBCHANGED      (SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED)
#define SCIP_EVENTTYPE_UBCHANGED      (SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED)
#define SCIP_EVENTTYPE_BOUNDTIGHTENED (SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBTIGHTENED)
#define SCIP_EVENTTYPE_BOUNDRELAXED   (SCIP_EVENTTYPE_LBRELAXED | SCIP_EVENTTYPE_UBRELAXED)
#define SCIP_EVENTTYPE_BOUNDCHANGED   (SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED)
#define SCIP_EVENTTYPE_HOLECHANGED    (SCIP_EVENTTYPE_HOLEADDED | SCIP_EVENTTYPE_HOLEREMOVED)
#define SCIP_EVENTTYPE_DOMCHANGED     (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_HOLECHANGED)
#define SCIP_EVENTTYPE_VARCHANGED     (SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_OBJCHANGED | SCIP_EVENTTYPE_DOMCHANGED \
                                       | SCIP_EVENTTYPE_LOCKSCHANGED)
#define SCIP_EVENTTYPE_VAREVENT       (SCIP_EVENTTYPE_VARCREATED | SCIP_EVENTTYPE_VARCHANGED)
   
/* event masks for node events */
#define SCIP_EVENTTYPE_NODESOLVED     (SCIP_EVENTTYPE_NODEFEASIBLE | SCIP_EVENTTYPE_NODEINFEASIBLE \
                                       | SCIP_EVENTTYPE_NODEBRANCHED)
#define SCIP_EVENTTYPE_NODEEVENT      (SCIP_EVENTTYPE_NODEFOCUSED | SCIP_EVENTTYPE_NODESOLVED)

/* event masks for LP events */
#define SCIP_EVENTTYPE_LPEVENT        (SCIP_EVENTTYPE_FIRSTLPSOLVED | SCIP_EVENTTYPE_LPSOLVED)

/* event masks for primal solution events */
#define SCIP_EVENTTYPE_SOLFOUND       (SCIP_EVENTTYPE_POORSOLFOUND | SCIP_EVENTTYPE_BESTSOLFOUND)
#define SCIP_EVENTTYPE_SOLEVENT       (SCIP_EVENTTYPE_SOLFOUND)

typedef unsigned int EVENTTYPE;         /**< type of event (bit field) */


typedef struct Eventhdlr EVENTHDLR;     /**< event handler for a specific events */
typedef struct EventhdlrData EVENTHDLRDATA; /**< event handler data */
typedef struct Event EVENT;             /**< event data structure */
typedef struct EventVarAdded EVENTVARADDED; /**< data for variable addition events */
typedef struct EventVarFixed EVENTVARFIXED; /**< data for variable fixing events */
typedef struct EventLocksChg EVENTLOCKSCHG; /**< data for locks change events */
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

/** initialization method of event handler (called after problem was transformed)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    eventhdlr       : the event handler itself
 */
#define DECL_EVENTINIT(x) RETCODE x (SCIP* scip, EVENTHDLR* eventhdlr)

/** deinitialization method of event handler (called before transformed problem is freed)
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



#include "def.h"
#include "type_retcode.h"
#include "type_scip.h"


#endif
