/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event.c
 * @brief  datastructures and methods for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "event.h"
#include "var.h"


/** event handler */
struct EventHdlr
{
   char*            name;               /**< name of event handler */
   char*            desc;               /**< description of event handler */
   DECL_EVENTFREE((*eventfree));        /**< destructor of event handler */
   DECL_EVENTINIT((*eventinit));        /**< initialise event handler */
   DECL_EVENTEXIT((*eventexit));        /**< deinitialise event handler */
   DECL_EVENTDELE((*eventdele));        /**< free specific event data */
   DECL_EVENTEXEC((*eventexec));        /**< execute event handler */
   EVENTHDLRDATA*   eventhdlrdata;      /**< event handler data */
   unsigned int     initialized:1;      /**< is event handler initialized? */
};

/** data for bound change events */
struct Boundchange
{
   VAR*             var;                /**< variable whose bound changed */
   Real             oldbound;           /**< old bound before bound changed */
   Real             newbound;           /**< new bound after bound changed */
};
typedef struct Boundchange BOUNDCHANGE; /**< data for bound change events */

/** event data structure */
struct Event
{
   union
   {
      BOUNDCHANGE   boundchange;        /**< data for bound change events */
   } data;
   EVENTTYPE        eventtype;          /**< type of event */
};

/** event queue to cache events and process them later */
struct EventQueue
{
   EVENT**          events;             /**< array with queued events */
   int              eventssize;         /**< number of available slots in events array */
   int              nevents;            /**< number of events in queue (used slots if events array) */
   unsigned int     delayevents:1;      /**< should the events be delayed and processed later? */
};



/*
 * Event handler methods
 */

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
   )
{
   assert(eventhdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(eventexec != NULL);

   ALLOC_OKAY( allocMemory(*eventhdlr) );
   ALLOC_OKAY( duplicateMemoryArray((*eventhdlr)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray((*eventhdlr)->desc, desc, strlen(desc)+1) );
   (*eventhdlr)->eventfree = eventfree;
   (*eventhdlr)->eventinit = eventinit;
   (*eventhdlr)->eventexit = eventexit;
   (*eventhdlr)->eventdele = eventdele;
   (*eventhdlr)->eventexec = eventexec;
   (*eventhdlr)->eventhdlrdata = eventhdlrdata;
   (*eventhdlr)->initialized = FALSE;

   return SCIP_OKAY;
}

RETCODE SCIPeventhdlrFree(              /**< calls destructor and frees memory of event handler */
   EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(eventhdlr != NULL);
   assert(*eventhdlr != NULL);
   assert(!(*eventhdlr)->initialized);
   assert(scip != NULL);

   /* call destructor of event handler */
   if( (*eventhdlr)->eventfree != NULL )
   {
      CHECK_OKAY( (*eventhdlr)->eventfree(*eventhdlr, scip) );
   }

   freeMemoryArray((*eventhdlr)->name);
   freeMemoryArray((*eventhdlr)->desc);
   freeMemory(*eventhdlr);

   return SCIP_OKAY;
}

RETCODE SCIPeventhdlrInit(              /**< initializes event handler */
   EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(eventhdlr != NULL);
   assert(scip != NULL);

   if( eventhdlr->initialized )
   {
      char s[255];
      sprintf(s, "Event handler <%s> already initialized", eventhdlr->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( eventhdlr->eventinit != NULL )
   {
      CHECK_OKAY( eventhdlr->eventinit(eventhdlr, scip) );
   }
   eventhdlr->initialized = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPeventhdlrExit(              /**< calls exit method of event handler */
   EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(eventhdlr != NULL);
   assert(scip != NULL);

   if( !eventhdlr->initialized )
   {
      char s[255];
      sprintf(s, "Event handler <%s> not initialized", eventhdlr->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( eventhdlr->eventexit != NULL )
   {
      CHECK_OKAY( eventhdlr->eventexit(eventhdlr, scip) );
   }
   eventhdlr->initialized = FALSE;

   return SCIP_OKAY;
}

RETCODE SCIPeventhdlrExec(              /**< calls execution method of event handler */
   EVENTHDLR*       eventhdlr,          /**< event handler */
   const SET*       set,                /**< global SCIP settings */
   EVENT*           event,              /**< event to call event handler with */
   EVENTDATA*       eventdata           /**< user data for the issued event */
   )
{
   assert(eventhdlr != NULL);
   assert(eventhdlr->eventexec != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(event != NULL);

   debugMessage("execute event of handler <%s> with event %p of type 0x%x\n", eventhdlr->name, event, event->eventtype);

   CHECK_OKAY( eventhdlr->eventexec(eventhdlr, set->scip, event, eventdata) );

   return SCIP_OKAY;
}

const char* SCIPeventhdlrGetName(       /**< gets name of event handler */
   EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return eventhdlr->name;
}

EVENTHDLRDATA* SCIPeventhdlrGetData(    /**< gets user data of event handler */
   EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return eventhdlr->eventhdlrdata;
}

void SCIPeventhdlrSetData(              /**< sets user data of event handler; user has to free old data in advance! */
   EVENTHDLR*       eventhdlr,          /**< event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< new event handler user data */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventhdlrdata = eventhdlrdata;
}

Bool SCIPeventhdlrIsInitialized(        /**< is event handler initialized? */
   EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return eventhdlr->initialized;
}




/*
 * Event methods
 */

RETCODE SCIPeventCreateLbChanged(       /**< creates an event for a change in the lower bound of a variable */
   EVENT**          event,              /**< pointer to store the event */
   MEMHDR*          memhdr,             /**< block memory */
   VAR*             var,                /**< variable whose bound changed */
   Real             oldbound,           /**< old bound before bound changed */
   Real             newbound            /**< new bound after bound changed */
   )
{
   assert(event != NULL);
   assert(memhdr != NULL);

   /* create event data */
   ALLOC_OKAY( allocBlockMemory(memhdr, *event) );
   if( newbound > oldbound )
      (*event)->eventtype = SCIP_EVENTTYPE_LBTIGHTENED;
   else
      (*event)->eventtype = SCIP_EVENTTYPE_LBRELAXED;
   (*event)->data.boundchange.var = var;
   (*event)->data.boundchange.oldbound = oldbound;
   (*event)->data.boundchange.newbound = newbound;

   return SCIP_OKAY;
}

RETCODE SCIPeventCreateUbChanged(       /**< creates an event for a change in the upper bound of a variable */
   EVENT**          event,              /**< pointer to store the event */
   MEMHDR*          memhdr,             /**< block memory */
   VAR*             var,                /**< variable whose bound changed */
   Real             oldbound,           /**< old bound before bound changed */
   Real             newbound            /**< new bound after bound changed */
   )
{
   assert(event != NULL);
   assert(memhdr != NULL);

   /* create event data */
   ALLOC_OKAY( allocBlockMemory(memhdr, *event) );
   if( newbound < oldbound )
      (*event)->eventtype = SCIP_EVENTTYPE_UBTIGHTENED;
   else
      (*event)->eventtype = SCIP_EVENTTYPE_UBRELAXED;
   (*event)->data.boundchange.var = var;
   (*event)->data.boundchange.oldbound = oldbound;
   (*event)->data.boundchange.newbound = newbound;

   return SCIP_OKAY;
}

RETCODE SCIPeventFree(                  /**< frees an event */
   EVENT**          event,              /**< event to free */
   MEMHDR*          memhdr              /**< block memory buffer */
   )
{
   assert(event != NULL);
   assert(memhdr != NULL);

   freeBlockMemory(memhdr, *event);

   return SCIP_OKAY;
}

static
void eventDisable(                      /**< disables an event */
   EVENT*           event               /**< event to disable */
   )
{
   assert(event != NULL);

   event->eventtype = SCIP_EVENTTYPE_DISABLED;
}

RETCODE SCIPeventGetType(               /**< gets type of event */
   EVENT*           event,              /**< event */
   EVENTTYPE*       eventtype           /**< pointer to store the event type */
   )
{
   assert(event != NULL);
   assert(eventtype != NULL);

   *eventtype = event->eventtype;

   return SCIP_OKAY;
}

RETCODE SCIPeventGetVar(                /**< gets variable for a domain change event */
   EVENT*           event,              /**< event */
   VAR**            var                 /**< pointer to store the variable */
   )
{
   assert(event != NULL);
   assert(var != NULL);

   switch( event->eventtype )
   {
   case SCIP_EVENTTYPE_VARCREATED:
      errorMessage("VARCREATED event not implemented yet");
      abort();
   case SCIP_EVENTTYPE_VARFIXED:
      errorMessage("VARFIXED event not implemented yet");
      abort();

   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      *var = event->data.boundchange.var;
      assert(*var != NULL);
      break;

   case SCIP_EVENTTYPE_HOLEADDED:
      errorMessage("HOLEADDED event not implemented yet");
      abort();
   case SCIP_EVENTTYPE_HOLEREMOVED:
      errorMessage("HOLEREMOVED event not implemented yet");
      abort();

   default:
      errorMessage("event does not belong to a variable");
      *var = NULL;
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

RETCODE SCIPeventGetOldbound(           /**< gets old bound for a bound change event */
   EVENT*           event,              /**< event */
   Real*            bound               /**< pointer to store the bound */
   )
{
   assert(event != NULL);
   assert(bound != NULL);

   switch( event->eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      *bound = event->data.boundchange.oldbound;
      break;

   default:
      errorMessage("event is not a bound change event");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

RETCODE SCIPeventGetNewbound(           /**< gets new bound for a bound change event */
   EVENT*           event,              /**< event */
   Real*            bound               /**< pointer to store the bound */
   )
{
   assert(event != NULL);
   assert(bound != NULL);

   switch( event->eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      *bound = event->data.boundchange.newbound;
      break;

   default:
      errorMessage("event is not a bound change event");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

RETCODE SCIPeventProcess(               /**< processes event by calling the appropriate event handlers */
   EVENT*           event,              /**< event */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   VAR* var;

   assert(event != NULL);

   debugMessage("processing event of type 0x%x\n", event->eventtype);
   
   switch( event->eventtype )
   {
   case SCIP_EVENTTYPE_DISABLED:
      break;
   case SCIP_EVENTTYPE_VARCREATED:
   case SCIP_EVENTTYPE_VARFIXED:
   case SCIP_EVENTTYPE_HOLEADDED:
   case SCIP_EVENTTYPE_HOLEREMOVED:
   case SCIP_EVENTTYPE_NODEACTIVATED:
   case SCIP_EVENTTYPE_NODESOLVED:
   case SCIP_EVENTTYPE_NODEBRANCHED:
   case SCIP_EVENTTYPE_NODEINFEASIBLE:
   case SCIP_EVENTTYPE_NODEFEASIBLE:
   case SCIP_EVENTTYPE_FIRSTLPSOLVED:
   case SCIP_EVENTTYPE_LPSOLVED:
      todoMessage("standard event filter not implemented yet");
      break;

   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
      var = event->data.boundchange.var;
      assert(var != NULL);
      assert(var->eventqueueindexlb == -1);

      /* inform LP and tree about bound change */
      if( var->varstatus == SCIP_VARSTATUS_COLUMN || var->varstatus == SCIP_VARSTATUS_LOOSE )
      {
         if( var->varstatus == SCIP_VARSTATUS_COLUMN )
         {
            CHECK_OKAY( SCIPcolBoundChanged(var->data.col, memhdr, set, lp, SCIP_BOUNDTYPE_LOWER,
                           event->data.boundchange.oldbound, event->data.boundchange.newbound) );
         }
         CHECK_OKAY( SCIPtreeBoundChanged(tree, memhdr, set, var, SCIP_BOUNDTYPE_LOWER,
                        event->data.boundchange.oldbound, event->data.boundchange.newbound) );
         CHECK_OKAY( SCIPbranchcandUpdateVar(branchcand, set, var) );
      }

      /* process variable's event filter */
      CHECK_OKAY( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      var = event->data.boundchange.var;
      assert(var != NULL);
      assert(var->eventqueueindexub == -1);

      /* inform LP and tree about bound change */
      if( var->varstatus == SCIP_VARSTATUS_COLUMN || var->varstatus == SCIP_VARSTATUS_LOOSE )
      {
         if( var->varstatus == SCIP_VARSTATUS_COLUMN )
         {
            CHECK_OKAY( SCIPcolBoundChanged(var->data.col, memhdr, set, lp, SCIP_BOUNDTYPE_UPPER,
                           event->data.boundchange.oldbound, event->data.boundchange.newbound) );
         }
         CHECK_OKAY( SCIPtreeBoundChanged(tree, memhdr, set, var, SCIP_BOUNDTYPE_UPPER,
                        event->data.boundchange.oldbound, event->data.boundchange.newbound) );
         CHECK_OKAY( SCIPbranchcandUpdateVar(branchcand, set, var) );
      }

      /* process variable's event filter */
      CHECK_OKAY( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   default:
      {
         char s[255];
         sprintf(s, "unknown event type <%d>", event->eventtype);
         errorMessage(s);
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}



/*
 * Event filter methods
 */

static
RETCODE eventfilterEnsureMem(           /**< resizes eventfilter arrays to be able to store at least num entries */
   EVENTFILTER*     eventfilter,        /**< event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(eventfilter != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);

   if( num > eventfilter->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, eventfilter->eventtypes, eventfilter->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, eventfilter->eventhdlrs, eventfilter->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, eventfilter->eventdatas, eventfilter->size, newsize) );
      eventfilter->size = newsize;
   }
   assert(num <= eventfilter->size);
   
   return SCIP_OKAY;
}

RETCODE SCIPeventfilterCreate(          /**< creates an event filter */
   EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   MEMHDR*          memhdr              /**< block memory buffer */
   )
{
   assert(eventfilter != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, *eventfilter) );
   (*eventfilter)->eventtypes = NULL;
   (*eventfilter)->eventhdlrs = NULL;
   (*eventfilter)->eventdatas = NULL;
   (*eventfilter)->size = 0;
   (*eventfilter)->len = 0;
   
   return SCIP_OKAY;
}

RETCODE SCIPeventfilterFree(            /**< frees an event filter and the associated event data entries */
   EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;

   assert(eventfilter != NULL);
   assert(*eventfilter != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   /* free event data */
   for( i = 0; i < (*eventfilter)->len; ++i )
   {
      assert((*eventfilter)->eventhdlrs[i] != NULL);
      if( (*eventfilter)->eventhdlrs[i]->eventdele != NULL )
      {
         CHECK_OKAY( (*eventfilter)->eventhdlrs[i]->eventdele((*eventfilter)->eventhdlrs[i], set->scip,
                        &(*eventfilter)->eventdatas[i]) );
      }
   }

   /* free event filter data */
   freeBlockMemoryArrayNull(memhdr, (*eventfilter)->eventtypes, (*eventfilter)->size);
   freeBlockMemoryArrayNull(memhdr, (*eventfilter)->eventhdlrs, (*eventfilter)->size);
   freeBlockMemoryArrayNull(memhdr, (*eventfilter)->eventdatas, (*eventfilter)->size);
   freeBlockMemory(memhdr, *eventfilter);

   return SCIP_OKAY;
}

RETCODE SCIPeventfilterAdd(             /**< adds element to event filter */
   EVENTFILTER*     eventfilter,        /**< event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   )
{
   int left;
   int middle;
   int right;
   int i;

   assert(eventfilter != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(eventhdlr != NULL);

   debugMessage("adding event handler %p with data %p for type mask 0x%x to event filter %p\n",
      eventhdlr, eventdata, eventtype, eventfilter);

   /* binary search the insert position, such that elements are sorted by eventhdlr and eventdata pointers */
   left = -1;
   right = eventfilter->len;
   while( left < right-1 )
   {
      middle = (left+right)/2;
      assert(left < middle && middle < right);
      assert(0 <= middle && middle < eventfilter->len);
      if( eventhdlr < eventfilter->eventhdlrs[middle]
         || (eventhdlr == eventfilter->eventhdlrs[middle] && eventdata < eventfilter->eventdatas[middle]) )
      {
         right = middle;
      }
      else
      {
         assert(eventhdlr > eventfilter->eventhdlrs[middle]
            || (eventhdlr == eventfilter->eventhdlrs[middle] && eventdata > eventfilter->eventdatas[middle]));
         left = middle;
      }
   }
   assert(left == right-1);
   assert(0 <= right && right <= eventfilter->len);

   /* insert event catch at position 'right' */
   CHECK_OKAY( eventfilterEnsureMem(eventfilter, memhdr, set, eventfilter->len+1) );
   eventfilter->len++;
   for( i = eventfilter->len-1; i > right; --i )
   {
      eventfilter->eventtypes[i] = eventfilter->eventtypes[i-1];
      eventfilter->eventhdlrs[i] = eventfilter->eventhdlrs[i-1];
      eventfilter->eventdatas[i] = eventfilter->eventdatas[i-1];
   }
   eventfilter->eventtypes[right] = eventtype;
   eventfilter->eventhdlrs[right] = eventhdlr;
   eventfilter->eventdatas[right] = eventdata;

   return SCIP_OKAY;
}

static
int eventfilterSearch(                  /**< binary search for the given event catch in event filter */
   EVENTFILTER*     eventfilter,        /**< event filter */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   )
{
   int left;
   int middle;
   int right;

   assert(eventfilter != NULL);
   assert(eventhdlr != NULL);

   /* binary search the position of the event catch */
   left = 0;
   right = eventfilter->len-1;
   while( left <= right )
   {
      middle = (left+right)/2;
      assert(left <= middle && middle <= right);
      assert(0 <= middle && middle < eventfilter->len);
      if( eventhdlr < eventfilter->eventhdlrs[middle]
         || (eventhdlr == eventfilter->eventhdlrs[middle] && eventdata < eventfilter->eventdatas[middle]) )
      {
         right = middle-1;
      }
      else if( eventhdlr > eventfilter->eventhdlrs[middle]
         || (eventhdlr == eventfilter->eventhdlrs[middle] && eventdata > eventfilter->eventdatas[middle]) )
      {
         left = middle+1;
      }
      else
      {
         return middle;
      }
   }

   return -1;
}

RETCODE SCIPeventfilterDel(             /**< deletes element from event filter */
   EVENTFILTER*     eventfilter,        /**< event filter */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   )
{
   int pos;
   int i;

   assert(eventfilter != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);
   assert(eventhdlr != NULL);

   debugMessage("deleting event handler %p with data %p from event filter %p\n",
      eventhdlr, eventdata, eventfilter);

   pos = eventfilterSearch(eventfilter, eventhdlr, eventdata);
   if( pos == -1 )
   {
      char s[255];
      sprintf(s, "no event for event handler %p with data %p found in event filter %p\n",
         eventhdlr, eventdata, eventfilter);
      errorMessage(s);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < eventfilter->len);
   assert(eventfilter->eventhdlrs[pos] == eventhdlr);
   assert(eventfilter->eventdatas[pos] == eventdata);

   /* move remaining events to fill the empty slot */
   eventfilter->len--;
   for( i = pos; i < eventfilter->len; ++i )
   {
      eventfilter->eventtypes[i] = eventfilter->eventtypes[i+1];
      eventfilter->eventhdlrs[i] = eventfilter->eventhdlrs[i+1];
      eventfilter->eventdatas[i] = eventfilter->eventdatas[i+1];
   }

   return SCIP_OKAY;
}

RETCODE SCIPeventfilterProcess(         /**< processes the event with all event handlers with matching filter setting */
   EVENTFILTER*     eventfilter,        /**< event filter */
   const SET*       set,                /**< global SCIP settings */
   EVENT*           event               /**< event to process */
   )
{
   int i;

   assert(eventfilter != NULL);
   assert(set != NULL);
   assert(event != NULL);

   debugMessage("processing event filter %p with event type 0x%x\n", eventfilter, event->eventtype);

   for( i = 0; i < eventfilter->len; ++i )
   {
      /* check, if event is applicable for the filter element */
      if( (event->eventtype & eventfilter->eventtypes[i]) != 0 )
      {
         /* call event handler */
         CHECK_OKAY( SCIPeventhdlrExec(eventfilter->eventhdlrs[i], set, event, eventfilter->eventdatas[i]) );         
      }
   }

   return SCIP_OKAY;
}



/*
 * Event queue methods
 */

static
RETCODE eventqueueEnsureEventsMem(      /**< resizes events array to be able to store at least num entries */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(eventqueue != NULL);
   assert(set != NULL);

   if( num > eventqueue->eventssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(eventqueue->events, newsize) );
      eventqueue->eventssize = newsize;
   }
   assert(num <= eventqueue->eventssize);
   
   return SCIP_OKAY;
}

RETCODE SCIPeventqueueCreate(           /**< creates an event queue */
   EVENTQUEUE**     eventqueue          /**< pointer to store the event queue */
   )
{
   assert(eventqueue != NULL);

   ALLOC_OKAY( allocMemory(*eventqueue) );
   (*eventqueue)->events = NULL;
   (*eventqueue)->eventssize = 0;
   (*eventqueue)->nevents = 0;
   (*eventqueue)->delayevents = FALSE;

   return SCIP_OKAY;
}

RETCODE SCIPeventqueueFree(             /**< frees event queue; there must not be any unprocessed eventy in the queue! */
   EVENTQUEUE**     eventqueue          /**< pointer to the event queue */
   )
{
   assert(eventqueue != NULL);
   assert(*eventqueue != NULL);
   assert((*eventqueue)->nevents == 0);

   freeMemoryArrayNull((*eventqueue)->events);
   freeMemory(*eventqueue);
   
   return SCIP_OKAY;
}

static
RETCODE eventqueueAppend(               /**< appends event to the event queue; sets event to NULL afterwards */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   const SET*       set,                /**< global SCIP settings */
   EVENT**          event               /**< pointer to event to append to the queue */
   )
{
   int insertpos;

   assert(eventqueue != NULL);
   assert(eventqueue->delayevents);
   assert(event != NULL);
   assert(*event != NULL);

   debugMessage("appending event %p of type 0x%x to event queue %p at position %d\n",
      *event, (*event)->eventtype, eventqueue, eventqueue->nevents);

   CHECK_OKAY( eventqueueEnsureEventsMem(eventqueue, set, eventqueue->nevents+1) );
   eventqueue->events[eventqueue->nevents] = *event;
   eventqueue->nevents++;

   *event = NULL;

   return SCIP_OKAY;
}

RETCODE SCIPeventqueueAdd(              /**< processes event or adds event to the event queue */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENT**          event               /**< pointer to event to add to the queue; will be NULL after queue addition */
   )
{
   VAR* var;
   EVENT* qevent;
   int pos;

   assert(eventqueue != NULL);
   assert(event != NULL);
   assert(*event != NULL);

   if( !eventqueue->delayevents )
   {
      /* immediately process event */
      CHECK_OKAY( SCIPeventProcess(*event, memhdr, set, tree, lp, branchcand) );
      CHECK_OKAY( SCIPeventFree(event, memhdr) );
   }
   else
   {
      /* delay processing of event by appending it to the event queue */
      debugMessage("adding event %p of type 0x%x to event queue %p\n", *event, (*event)->eventtype, eventqueue);

      switch( (*event)->eventtype )
      {
      case SCIP_EVENTTYPE_DISABLED:
         errorMessage("cannot add a disabled event to the event queue");
         return SCIP_INVALIDDATA;
      case SCIP_EVENTTYPE_VARCREATED:
      case SCIP_EVENTTYPE_VARFIXED:
      case SCIP_EVENTTYPE_HOLEADDED:
      case SCIP_EVENTTYPE_HOLEREMOVED:
      case SCIP_EVENTTYPE_NODEACTIVATED:
      case SCIP_EVENTTYPE_NODESOLVED:
      case SCIP_EVENTTYPE_NODEBRANCHED:
      case SCIP_EVENTTYPE_NODEINFEASIBLE:
      case SCIP_EVENTTYPE_NODEFEASIBLE:
      case SCIP_EVENTTYPE_FIRSTLPSOLVED:
      case SCIP_EVENTTYPE_LPSOLVED:
         /* these events cannot be merged; just add them to the queue */
         CHECK_OKAY( eventqueueAppend(eventqueue, set, event) );
         break;
      case SCIP_EVENTTYPE_LBTIGHTENED:
      case SCIP_EVENTTYPE_LBRELAXED:
         /* changes in lower bound may be merged with older changes in lower bound */
         var = (*event)->data.boundchange.var;
         assert(var != NULL);
         pos = var->eventqueueindexlb;
         if( pos >= 0 )
         {
            /* the lower bound change event already exists -> modifiy it accordingly */
            assert(pos < eventqueue->nevents);
            qevent = eventqueue->events[pos];
            assert(qevent != NULL);
            assert(qevent->eventtype == SCIP_EVENTTYPE_LBTIGHTENED || qevent->eventtype == SCIP_EVENTTYPE_LBRELAXED);
            assert(qevent->data.boundchange.var == var);
            assert(SCIPsetIsEQ(set, (*event)->data.boundchange.oldbound, qevent->data.boundchange.newbound));

            debugMessage(" -> merging LB event (<%s>,%g -> %g) with event at position %d (<%s>,%g -> %g)\n",
               (*event)->data.boundchange.var->name, (*event)->data.boundchange.oldbound,
               (*event)->data.boundchange.newbound,
               pos, qevent->data.boundchange.var->name, qevent->data.boundchange.oldbound, 
               qevent->data.boundchange.newbound);

            qevent->data.boundchange.newbound = (*event)->data.boundchange.newbound;
            if( SCIPsetIsL(set, qevent->data.boundchange.newbound, qevent->data.boundchange.oldbound) )
               qevent->eventtype = SCIP_EVENTTYPE_LBRELAXED;
            else if( SCIPsetIsG(set, qevent->data.boundchange.newbound, qevent->data.boundchange.oldbound) )
               qevent->eventtype = SCIP_EVENTTYPE_LBTIGHTENED;
            else
            {
               /* the queued bound change was reversed -> disable the event in the queue */
               assert(SCIPsetIsEQ(set, qevent->data.boundchange.newbound, qevent->data.boundchange.oldbound));
               eventDisable(qevent);
               var->eventqueueindexlb = -1;
               debugMessage(" -> event disabled\n");
            }

            /* free the event that is of no use any longer */
            CHECK_OKAY( SCIPeventFree(event, memhdr) );
         }
         else
         {
            /* the lower bound change event doesn't exist -> add it to the queue, and remember the array index */
            var->eventqueueindexlb = eventqueue->nevents;
            CHECK_OKAY( eventqueueAppend(eventqueue, set, event) );
         }
         break;
      case SCIP_EVENTTYPE_UBTIGHTENED:
      case SCIP_EVENTTYPE_UBRELAXED:
         /* changes in upper bound may be merged with older changes in upper bound */
         var = (*event)->data.boundchange.var;
         assert(var != NULL);
         pos = var->eventqueueindexub;
         if( pos >= 0 )
         {
            /* the upper bound change event already exists -> modifiy it accordingly */
            assert(pos < eventqueue->nevents);
            qevent = eventqueue->events[pos];
            assert(qevent != NULL);
            assert(qevent->eventtype == SCIP_EVENTTYPE_UBTIGHTENED || qevent->eventtype == SCIP_EVENTTYPE_UBRELAXED);
            assert(qevent->data.boundchange.var == var);
            assert(SCIPsetIsEQ(set, (*event)->data.boundchange.oldbound, qevent->data.boundchange.newbound));

            debugMessage(" -> merging UB event (<%s>,%g -> %g) with event at position %d (<%s>,%g -> %g)\n",
               (*event)->data.boundchange.var->name, (*event)->data.boundchange.oldbound,
               (*event)->data.boundchange.newbound,
               pos, qevent->data.boundchange.var->name, qevent->data.boundchange.oldbound,
               qevent->data.boundchange.newbound);

            qevent->data.boundchange.newbound = (*event)->data.boundchange.newbound;
            if( SCIPsetIsL(set, qevent->data.boundchange.newbound, qevent->data.boundchange.oldbound) )
               qevent->eventtype = SCIP_EVENTTYPE_UBTIGHTENED;
            else if( SCIPsetIsG(set, qevent->data.boundchange.newbound, qevent->data.boundchange.oldbound) )
               qevent->eventtype = SCIP_EVENTTYPE_UBRELAXED;
            else
            {
               /* the queued bound change was reversed -> disable the event in the queue */
               assert(SCIPsetIsEQ(set, qevent->data.boundchange.newbound, qevent->data.boundchange.oldbound));
               eventDisable(qevent);
               var->eventqueueindexub = -1;
            }

            /* free the event that is of no use any longer */
            CHECK_OKAY( SCIPeventFree(event, memhdr) );
         }
         else
         {
            /* the upper bound change event doesn't exist -> add it to the queue, and remember the array index */
            var->eventqueueindexub = eventqueue->nevents;
            CHECK_OKAY( eventqueueAppend(eventqueue, set, event) );
         }
         break;
      default:
         {
            char s[255];
            sprintf(s, "unknown event type <%d>", (*event)->eventtype);
            errorMessage(s);
            return SCIP_INVALIDDATA;
         }
      }
   }
   
   assert(*event == NULL);

   return SCIP_OKAY;
}

RETCODE SCIPeventqueueDelay(            /**< marks queue to delay incoming events until a call to SCIPeventqueueProcess() */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(eventqueue != NULL);
   assert(!eventqueue->delayevents);

   eventqueue->delayevents = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPeventqueueProcess(          /**< processes all delayed events, marks queue to process events immediately */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   EVENT* event;
   int i;

   assert(eventqueue != NULL);
   assert(eventqueue->delayevents);

   /* pass events to the responsible event filters
    * During event processing, new events may be raised. We have to loop to the mutable eventqueue->nevents.
    * A loop to something like "nevents = eventqueue->nevents; for(...; i < nevents; ...)" would miss the
    * newly created events. The same holds for eventqueue->events, which can be moved in memory due to
    * memory reallocation in eventqueueAppend().
    */
   for( i = 0; i < eventqueue->nevents; ++i )
   {
      event = eventqueue->events[i];
      assert(event != NULL);

      debugMessage("processing event %d of %d events in queue: eventtype=0x%x\n", i, eventqueue->nevents, event->eventtype);

      /* unmark the event queue index of a variable with changed bounds */
      if( (event->eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
      {
         assert(event->data.boundchange.var->eventqueueindexlb == i);
         event->data.boundchange.var->eventqueueindexlb = -1;
      }
      else if( (event->eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0 )
      {
         assert(event->data.boundchange.var->eventqueueindexub == i);
         event->data.boundchange.var->eventqueueindexub = -1;
      }

      CHECK_OKAY( SCIPeventProcess(event, memhdr, set, tree, lp, branchcand) );

      /* free the event immediately, because additionally raised events during event processing
       * can lead to a large event queue
       */
      CHECK_OKAY( SCIPeventFree(&eventqueue->events[i], memhdr) );
   }

   assert(i == eventqueue->nevents);
   eventqueue->nevents = 0;
   eventqueue->delayevents = FALSE;

   return SCIP_OKAY;
}
