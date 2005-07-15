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
#pragma ident "@(#) $Id: struct_event.h,v 1.16 2005/07/15 17:20:19 bzfpfend Exp $"

/**@file   struct_event.h
 * @brief  datastructures for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_EVENT_H__
#define __SCIP_STRUCT_EVENT_H__


#include "scip/def.h"
#include "scip/type_event.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"
#include "scip/type_tree.h"



/** data for variable addition events */
struct EventVarAdded
{
   VAR*             var;                /**< variable that was added to the problem */
};

/** data for variable deletion events */
struct EventVarDeleted
{
   VAR*             var;                /**< variable that will be deleted from the problem */
};

/** data for variable fixing events */
struct EventVarFixed
{
   VAR*             var;                /**< variable that was fixed */
};

/** data for locks change events */
struct EventLocksChg
{
   VAR*             var;                /**< variable for which the lock numbers were changed */
};

/** data for objective value change events */
struct EventObjChg
{
   Real             oldobj;             /**< old objective value before value changed */
   Real             newobj;             /**< new objective value after value changed */
   VAR*             var;                /**< variable whose objective value changed */
};

/** data for bound change events */
struct EventBdChg
{
   Real             oldbound;           /**< old bound before bound changed */
   Real             newbound;           /**< new bound after bound changed */
   VAR*             var;                /**< variable whose bound changed */
};

/** event data structure */
struct Event
{
   union
   {
      EVENTVARADDED eventvaradded;      /**< data for variable addition events */
      EVENTVARDELETED eventvardeleted;  /**< data for variable deletion events */
      EVENTVARFIXED eventvarfixed;      /**< data for variable fixing events */
      EVENTLOCKSCHG eventlockschg;      /**< data for locks change events */
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
   EVENTTYPE*       eventtypes;         /**< array with types of event to process; 0 marks a deleted event catch entry */
   EVENTHDLR**      eventhdlrs;         /**< array with event handlers to process the event */
   EVENTDATA**      eventdatas;         /**< array with user data for the issued event */
   int*             nextpos;            /**< linked lists for free, delayed added and delayed deleted slot positions */
   int              size;               /**< size of filter arrays (available slots in arrays) */
   int              len;                /**< number entries in filter arrays (used and deleted) */
   int              firstfreepos;       /**< first deleted slot; remaining slots are in poslist */
   int              firstdeletedpos;    /**< first delayed deleted slot; remaining slots are in poslist */
   unsigned int     eventmask;          /**< mask for events that are handled by any event handler in the filter */
   unsigned int     delayedeventmask;   /**< mask for delayed added events */
   Bool             delayupdates;       /**< should additions and deletions to the filter be delayed? */
};

/** event handler */
struct Eventhdlr
{
   char*            name;               /**< name of event handler */
   char*            desc;               /**< description of event handler */
   DECL_EVENTFREE   ((*eventfree));     /**< destructor of event handler */
   DECL_EVENTINIT   ((*eventinit));     /**< initialize event handler */
   DECL_EVENTEXIT   ((*eventexit));     /**< deinitialize event handler */
   DECL_EVENTINITSOL((*eventinitsol));  /**< solving process initialization method of event handler */
   DECL_EVENTEXITSOL((*eventexitsol));  /**< solving process deinitialization method of event handler */
   DECL_EVENTDELETE ((*eventdelete));   /**< free specific event data */
   DECL_EVENTEXEC   ((*eventexec));     /**< execute event handler */
   EVENTHDLRDATA*   eventhdlrdata;      /**< event handler data */
   Bool             initialized;        /**< is event handler initialized? */
};

/** event queue to cache events and process them later */
struct EventQueue
{
   EVENT**          events;             /**< array with queued events */
   int              eventssize;         /**< number of available slots in events array */
   int              nevents;            /**< number of events in queue (used slots if events array) */
   Bool             delayevents;        /**< should the events be delayed and processed later? */
};

#endif
