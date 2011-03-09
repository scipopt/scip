/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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

#ifdef __cplusplus
extern "C" {
#endif

/** data for variable addition events */
struct SCIP_EventVarAdded
{
   SCIP_VAR*             var;                /**< variable that was added to the problem */
};

/** data for variable deletion events */
struct SCIP_EventVarDeleted
{
   SCIP_VAR*             var;                /**< variable that will be deleted from the problem */
};

/** data for variable fixing events */
struct SCIP_EventVarFixed
{
   SCIP_VAR*             var;                /**< variable that was fixed */
};

/** data for locks change events */
struct SCIP_EventVarUnlocked
{
   SCIP_VAR*             var;                /**< variable for which the lock numbers were changed */
};

/** data for objective value change events */
struct SCIP_EventObjChg
{
   SCIP_Real             oldobj;             /**< old objective value before value changed */
   SCIP_Real             newobj;             /**< new objective value after value changed */
   SCIP_VAR*             var;                /**< variable whose objective value changed */
};

/** data for bound change events */
struct SCIP_EventBdChg
{
   SCIP_Real             oldbound;           /**< old bound before bound changed */
   SCIP_Real             newbound;           /**< new bound after bound changed */
   SCIP_VAR*             var;                /**< variable whose bound changed */
};

/** data for implication added events */
struct SCIP_EventImplAdd
{
   SCIP_VAR*             var;                /**< variable for which the lock numbers were changed */
};

/** event data structure */
struct SCIP_Event
{
   union
   {
      SCIP_EVENTVARADDED eventvaradded;      /**< data for variable addition events */
      SCIP_EVENTVARDELETED eventvardeleted;  /**< data for variable deletion events */
      SCIP_EVENTVARFIXED eventvarfixed;      /**< data for variable fixing events */
      SCIP_EVENTVARUNLOCKED eventvarunlocked;/**< data for locks change events */
      SCIP_EVENTOBJCHG   eventobjchg;        /**< data for objective value change events */
      SCIP_EVENTBDCHG    eventbdchg;         /**< data for bound change events */
      SCIP_EVENTIMPLADD  eventimpladd;       /**< data for implication added events */
      SCIP_NODE*         node;               /**< data for node and LP events */
      SCIP_SOL*          sol;                /**< data for primal solution events */
   } data;
   SCIP_EVENTTYPE        eventtype;          /**< type of event */
};

/** event filter to select events to be processed by an event handler */
struct SCIP_EventFilter
{
   SCIP_EVENTTYPE*       eventtypes;         /**< array with types of event to process; 0 marks a deleted event catch entry */
   SCIP_EVENTHDLR**      eventhdlrs;         /**< array with event handlers to process the event */
   SCIP_EVENTDATA**      eventdatas;         /**< array with user data for the issued event */
   int*                  nextpos;            /**< linked lists for free, delayed added and delayed deleted slot positions */
   int                   size;               /**< size of filter arrays (available slots in arrays) */
   int                   len;                /**< number entries in filter arrays (used and deleted) */
   int                   firstfreepos;       /**< first deleted slot; remaining slots are in poslist */
   int                   firstdeletedpos;    /**< first delayed deleted slot; remaining slots are in poslist */
   unsigned int          eventmask;          /**< mask for events that are handled by any event handler in the filter */
   unsigned int          delayedeventmask;   /**< mask for delayed added events */
   SCIP_Bool             delayupdates;       /**< should additions and deletions to the filter be delayed? */
};

/** event handler */
struct SCIP_Eventhdlr
{
   char*                 name;               /**< name of event handler */
   char*                 desc;               /**< description of event handler */
   SCIP_DECL_EVENTFREE   ((*eventfree));     /**< destructor of event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit));     /**< initialize event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit));     /**< deinitialize event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol));  /**< solving process initialization method of event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol));  /**< solving process deinitialization method of event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete));   /**< free specific event data */
   SCIP_DECL_EVENTEXEC   ((*eventexec));     /**< execute event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata;      /**< event handler data */
   SCIP_Bool             initialized;        /**< is event handler initialized? */
};

/** event queue to cache events and process them later */
struct SCIP_EventQueue
{
   SCIP_EVENT**          events;             /**< array with queued events */
   int                   eventssize;         /**< number of available slots in events array */
   int                   nevents;            /**< number of events in queue (used slots if events array) */
   SCIP_Bool             delayevents;        /**< should the events be delayed and processed later? */
};

#ifdef __cplusplus
}
#endif

#endif
