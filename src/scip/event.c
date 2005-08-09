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
#pragma ident "@(#) $Id: event.c,v 1.49 2005/08/09 16:27:06 bzfpfend Exp $"

/**@file   event.c
 * @brief  methods and datastructures for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/primal.h"
#include "scip/branch.h"




/*
 * Event handler methods
 */

/** creates an event handler */
RETCODE SCIPeventhdlrCreate(
   EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   const char*      name,               /**< name of event handler */
   const char*      desc,               /**< description of event handler */
   DECL_EVENTFREE   ((*eventfree)),     /**< destructor of event handler */
   DECL_EVENTINIT   ((*eventinit)),     /**< initialize event handler */
   DECL_EVENTEXIT   ((*eventexit)),     /**< deinitialize event handler */
   DECL_EVENTINITSOL((*eventinitsol)),  /**< solving process initialization method of event handler */
   DECL_EVENTEXITSOL((*eventexitsol)),  /**< solving process deinitialization method of event handler */
   DECL_EVENTDELETE ((*eventdelete)),   /**< free specific event data */
   DECL_EVENTEXEC   ((*eventexec)),     /**< execute event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   assert(eventhdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(eventexec != NULL);

   ALLOC_OKAY( allocMemory(eventhdlr) );
   ALLOC_OKAY( duplicateMemoryArray(&(*eventhdlr)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*eventhdlr)->desc, desc, strlen(desc)+1) );
   (*eventhdlr)->eventfree = eventfree;
   (*eventhdlr)->eventinit = eventinit;
   (*eventhdlr)->eventexit = eventexit;
   (*eventhdlr)->eventinitsol = eventinitsol;
   (*eventhdlr)->eventexitsol = eventexitsol;
   (*eventhdlr)->eventdelete = eventdelete;
   (*eventhdlr)->eventexec = eventexec;
   (*eventhdlr)->eventhdlrdata = eventhdlrdata;
   (*eventhdlr)->initialized = FALSE;

   return SCIP_OKAY;
}

/** calls destructor and frees memory of event handler */
RETCODE SCIPeventhdlrFree(
   EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(*eventhdlr != NULL);
   assert(!(*eventhdlr)->initialized);
   assert(set != NULL);

   /* call destructor of event handler */
   if( (*eventhdlr)->eventfree != NULL )
   {
      CHECK_OKAY( (*eventhdlr)->eventfree(set->scip, *eventhdlr) );
   }

   freeMemoryArray(&(*eventhdlr)->name);
   freeMemoryArray(&(*eventhdlr)->desc);
   freeMemory(eventhdlr);

   return SCIP_OKAY;
}

/** initializes event handler */
RETCODE SCIPeventhdlrInit(
   EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(set != NULL);

   if( eventhdlr->initialized )
   {
      errorMessage("event handler <%s> already initialized\n", eventhdlr->name);
      return SCIP_INVALIDCALL;
   }

   if( eventhdlr->eventinit != NULL )
   {
      CHECK_OKAY( eventhdlr->eventinit(set->scip, eventhdlr) );
   }
   eventhdlr->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of event handler */
RETCODE SCIPeventhdlrExit(
   EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(set != NULL);

   if( !eventhdlr->initialized )
   {
      errorMessage("event handler <%s> not initialized\n", eventhdlr->name);
      return SCIP_INVALIDCALL;
   }

   if( eventhdlr->eventexit != NULL )
   {
      CHECK_OKAY( eventhdlr->eventexit(set->scip, eventhdlr) );
   }
   eventhdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs event handler that the branch and bound process is being started */
RETCODE SCIPeventhdlrInitsol(
   EVENTHDLR*       eventhdlr,          /**< event handler */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(set != NULL);

   /* call solving process initialization method of event handler */
   if( eventhdlr->eventinitsol != NULL )
   {
      CHECK_OKAY( eventhdlr->eventinitsol(set->scip, eventhdlr) );
   }

   return SCIP_OKAY;
}

/** informs event handler that the branch and bound process data is being freed */
RETCODE SCIPeventhdlrExitsol(
   EVENTHDLR*       eventhdlr,          /**< event handler */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of event handler */
   if( eventhdlr->eventexitsol != NULL )
   {
      CHECK_OKAY( eventhdlr->eventexitsol(set->scip, eventhdlr) );
   }

   return SCIP_OKAY;
}

/** calls execution method of event handler */
RETCODE SCIPeventhdlrExec(
   EVENTHDLR*       eventhdlr,          /**< event handler */
   SET*             set,                /**< global SCIP settings */
   EVENT*           event,              /**< event to call event handler with */
   EVENTDATA*       eventdata           /**< user data for the issued event */
   )
{
   assert(eventhdlr != NULL);
   assert(eventhdlr->eventexec != NULL);
   assert(set != NULL);
   assert(event != NULL);

   debugMessage("execute event of handler <%s> with event %p of type 0x%x\n", eventhdlr->name, event, event->eventtype);

   CHECK_OKAY( eventhdlr->eventexec(set->scip, eventhdlr, event, eventdata) );

   return SCIP_OKAY;
}

/** gets name of event handler */
const char* SCIPeventhdlrGetName(
   EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return eventhdlr->name;
}

/** gets user data of event handler */
EVENTHDLRDATA* SCIPeventhdlrGetData(
   EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return eventhdlr->eventhdlrdata;
}

/** sets user data of event handler; user has to free old data in advance! */
void SCIPeventhdlrSetData(
   EVENTHDLR*       eventhdlr,          /**< event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< new event handler user data */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventhdlrdata = eventhdlrdata;
}

/** is event handler initialized? */
Bool SCIPeventhdlrIsInitialized(
   EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return eventhdlr->initialized;
}




/*
 * Event methods
 */

/** creates an event for an addition of a variable to the problem */
RETCODE SCIPeventCreateVarAdded(
   EVENT**          event,              /**< pointer to store the event */
   BLKMEM*          blkmem,             /**< block memory */
   VAR*             var                 /**< variable that was added to the problem */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   /* create event data */
   ALLOC_OKAY( allocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_VARADDED;
   (*event)->data.eventvaradded.var = var;

   return SCIP_OKAY;
}

/** creates an event for a deletion of a variable from the problem */
RETCODE SCIPeventCreateVarDeleted(
   EVENT**          event,              /**< pointer to store the event */
   BLKMEM*          blkmem,             /**< block memory */
   VAR*             var                 /**< variable that is to be deleted from the problem */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   /* create event data */
   ALLOC_OKAY( allocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_VARDELETED;
   (*event)->data.eventvardeleted.var = var;

   return SCIP_OKAY;
}

/** creates an event for a fixing of a variable */
RETCODE SCIPeventCreateVarFixed(
   EVENT**          event,              /**< pointer to store the event */
   BLKMEM*          blkmem,             /**< block memory */
   VAR*             var                 /**< variable that was fixed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);

   /* create event data */
   ALLOC_OKAY( allocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_VARFIXED;
   (*event)->data.eventvarfixed.var = var;

   return SCIP_OKAY;
}

/** creates an event for a change in the number of locks of a variable */
RETCODE SCIPeventCreateLocksChanged(
   EVENT**          event,              /**< pointer to store the event */
   BLKMEM*          blkmem,             /**< block memory */
   VAR*             var                 /**< variable that was fixed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED);

   /* create event data */
   ALLOC_OKAY( allocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_LOCKSCHANGED;
   (*event)->data.eventvarfixed.var = var;

   return SCIP_OKAY;
}

/** creates an event for a change in the objective value of a variable */
RETCODE SCIPeventCreateObjChanged(
   EVENT**          event,              /**< pointer to store the event */
   BLKMEM*          blkmem,             /**< block memory */
   VAR*             var,                /**< variable whose objective value changed */
   Real             oldobj,             /**< old objective value before value changed */
   Real             newobj              /**< new objective value after value changed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(oldobj != newobj); /*lint !e777*/

   /* create event data */
   ALLOC_OKAY( allocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_OBJCHANGED;
   (*event)->data.eventobjchg.var = var;
   (*event)->data.eventobjchg.oldobj = oldobj;
   (*event)->data.eventobjchg.newobj = newobj;

   return SCIP_OKAY;
}

/** creates an event for a change in the lower bound of a variable */
RETCODE SCIPeventCreateLbChanged(
   EVENT**          event,              /**< pointer to store the event */
   BLKMEM*          blkmem,             /**< block memory */
   VAR*             var,                /**< variable whose bound changed */
   Real             oldbound,           /**< old bound before bound changed */
   Real             newbound            /**< new bound after bound changed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(oldbound != newbound); /*lint !e777*/

   /* create event data */
   ALLOC_OKAY( allocBlockMemory(blkmem, event) );
   if( newbound > oldbound )
      (*event)->eventtype = SCIP_EVENTTYPE_LBTIGHTENED;
   else
      (*event)->eventtype = SCIP_EVENTTYPE_LBRELAXED;
   (*event)->data.eventbdchg.var = var;
   (*event)->data.eventbdchg.oldbound = oldbound;
   (*event)->data.eventbdchg.newbound = newbound;

   return SCIP_OKAY;
}

/** creates an event for a change in the upper bound of a variable */
RETCODE SCIPeventCreateUbChanged(
   EVENT**          event,              /**< pointer to store the event */
   BLKMEM*          blkmem,             /**< block memory */
   VAR*             var,                /**< variable whose bound changed */
   Real             oldbound,           /**< old bound before bound changed */
   Real             newbound            /**< new bound after bound changed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(oldbound != newbound); /*lint !e777*/

   /* create event data */
   ALLOC_OKAY( allocBlockMemory(blkmem, event) );
   if( newbound < oldbound )
      (*event)->eventtype = SCIP_EVENTTYPE_UBTIGHTENED;
   else
      (*event)->eventtype = SCIP_EVENTTYPE_UBRELAXED;
   (*event)->data.eventbdchg.var = var;
   (*event)->data.eventbdchg.oldbound = oldbound;
   (*event)->data.eventbdchg.newbound = newbound;

   return SCIP_OKAY;
}

/** frees an event */
RETCODE SCIPeventFree(
   EVENT**          event,              /**< event to free */
   BLKMEM*          blkmem              /**< block memory buffer */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   freeBlockMemory(blkmem, event);

   return SCIP_OKAY;
}

/** disables an event */
static
void eventDisable(
   EVENT*           event               /**< event to disable */
   )
{
   assert(event != NULL);

   event->eventtype = SCIP_EVENTTYPE_DISABLED;
}

/** gets type of event */
EVENTTYPE SCIPeventGetType(
   EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   return event->eventtype;
}

/** sets type of event */
RETCODE SCIPeventChgType(
   EVENT*           event,              /**< event */
   EVENTTYPE        eventtype           /**< new event type */
   )
{
   assert(event != NULL);

   event->eventtype = eventtype;

   return SCIP_OKAY;
}

/** gets variable for a variable event (var added, var deleted, var fixed, objective value or domain change) */
VAR* SCIPeventGetVar(
   EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   switch( event->eventtype )
   {  
   case SCIP_EVENTTYPE_VARADDED:
      assert(event->data.eventvaradded.var != NULL);
      return event->data.eventvaradded.var;

   case SCIP_EVENTTYPE_VARDELETED:
      assert(event->data.eventvardeleted.var != NULL);
      return event->data.eventvardeleted.var;

   case SCIP_EVENTTYPE_VARFIXED:
      assert(event->data.eventvarfixed.var != NULL);
      return event->data.eventvarfixed.var;

   case SCIP_EVENTTYPE_LOCKSCHANGED:
      assert(event->data.eventlockschg.var != NULL);
      return event->data.eventlockschg.var;

   case SCIP_EVENTTYPE_OBJCHANGED:
      assert(event->data.eventobjchg.var != NULL);
      return event->data.eventobjchg.var;

   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      assert(event->data.eventbdchg.var != NULL);
      return event->data.eventbdchg.var;

   case SCIP_EVENTTYPE_HOLEADDED:
      errorMessage("HOLEADDED event not implemented yet\n");
      SCIPABORT();
      return NULL;

   case SCIP_EVENTTYPE_HOLEREMOVED:
      errorMessage("HOLEREMOVED event not implemented yet\n");
      SCIPABORT();
      return NULL;

   default:
      errorMessage("event does not belong to a variable\n");
      SCIPABORT();
      return NULL;
   }  /*lint !e788*/
}

/** gets old objective value for an objective value change event */
Real SCIPeventGetOldobj(
   EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( event->eventtype != SCIP_EVENTTYPE_OBJCHANGED )
   {
      errorMessage("event is not an objective value change event\n");
      return SCIP_INVALID;
   }

   return event->data.eventobjchg.oldobj;
}

/** gets new objective value for an objective value change event */
Real SCIPeventGetNewobj(
   EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( event->eventtype != SCIP_EVENTTYPE_OBJCHANGED )
   {
      errorMessage("event is not an objective value change event\n");
      return SCIP_INVALID;
   }

   return event->data.eventobjchg.newobj;
}

/** gets old bound for a bound change event */
Real SCIPeventGetOldbound(
   EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   switch( event->eventtype )
   {  
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      return event->data.eventbdchg.oldbound;

   default:
      errorMessage("event is not a bound change event\n");
      SCIPABORT();
      return 0.0;
   }  /*lint !e788*/
}

/** gets new bound for a bound change event */
Real SCIPeventGetNewbound(
   EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   switch( event->eventtype )
   {  
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      return event->data.eventbdchg.newbound;

   default:
      errorMessage("event is not a bound change event\n");
      SCIPABORT();
      return 0.0;
   }  /*lint !e788*/
}

/** gets node for a node or LP event */
NODE* SCIPeventGetNode(
   EVENT*           event               /**< event */
   )
{
   assert(event != NULL);
   
   if( (event->eventtype & (SCIP_EVENTTYPE_NODEEVENT | SCIP_EVENTTYPE_LPEVENT)) == 0 )
   {
      errorMessage("event is neither node nor LP event\n");
      return NULL;
   }

   return event->data.node;
}

/** sets node for a node or LP event */
RETCODE SCIPeventChgNode(
   EVENT*           event,              /**< event */
   NODE*            node                /**< new node */
   )
{
   assert(event != NULL);

   if( (event->eventtype & (SCIP_EVENTTYPE_NODEEVENT | SCIP_EVENTTYPE_LPEVENT)) == 0 )
   {
      errorMessage("event is neither node nor LP event\n");
      return SCIP_INVALIDDATA;
   }

   event->data.node = node;

   return SCIP_OKAY;
}

/** gets solution for a primal solution event */
SOL* SCIPeventGetSol(
   EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_SOLEVENT) == 0 )
   {
      errorMessage("event is not a primal solution event\n");
      return NULL;
   }

   return event->data.sol;
}

/** sets solution for a primal solution event */
RETCODE SCIPeventChgSol(
   EVENT*           event,              /**< event */
   SOL*             sol                 /**< new primal solution */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_SOLEVENT) == 0 )
   {
      errorMessage("event is not a primal solution event\n");
      return SCIP_INVALIDDATA;
   }

   event->data.sol = sol;

   return SCIP_OKAY;
}

/** processes event by calling the appropriate event handlers */
RETCODE SCIPeventProcess(
   EVENT*           event,              /**< event */
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data; only needed for objchanged events */
   LP*              lp,                 /**< current LP data; only needed for obj/boundchanged events */
   BRANCHCAND*      branchcand,         /**< branching candidate storage; only needed for boundchange events */
   EVENTFILTER*     eventfilter         /**< event filter for global events; not needed for variable specific events */
   )
{
   VAR* var;

   assert(event != NULL);
   
   debugMessage("processing event of type 0x%x\n", event->eventtype);
   
   switch( event->eventtype )
   {
   case SCIP_EVENTTYPE_DISABLED:
      break;

   case SCIP_EVENTTYPE_VARADDED:
   case SCIP_EVENTTYPE_VARDELETED:
   case SCIP_EVENTTYPE_PRESOLVEROUND:
   case SCIP_EVENTTYPE_NODEFOCUSED:
   case SCIP_EVENTTYPE_NODEFEASIBLE:
   case SCIP_EVENTTYPE_NODEINFEASIBLE:
   case SCIP_EVENTTYPE_NODEBRANCHED:
   case SCIP_EVENTTYPE_FIRSTLPSOLVED:
   case SCIP_EVENTTYPE_LPSOLVED:
   case SCIP_EVENTTYPE_POORSOLFOUND:
   case SCIP_EVENTTYPE_BESTSOLFOUND:
      CHECK_OKAY( SCIPeventfilterProcess(eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_VARFIXED:
      var = event->data.eventvarfixed.var;
      assert(var != NULL);

      /* process variable's event filter */
      CHECK_OKAY( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_LOCKSCHANGED:
      var = event->data.eventlockschg.var;
      assert(var != NULL);

      /* process variable's event filter */
      CHECK_OKAY( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_OBJCHANGED:
      var = event->data.eventobjchg.var;
      assert(var != NULL);
      assert(var->eventqueueindexobj == -1);
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

      /* inform LP about the objective change */
      if( SCIPvarGetProbindex(var) >= 0 )
      {
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
         {
            CHECK_OKAY( SCIPcolChgObj(SCIPvarGetCol(var), set, lp, event->data.eventobjchg.newobj) );
         }
         CHECK_OKAY( SCIPlpUpdateVarObj(lp, set, var, event->data.eventobjchg.oldobj, event->data.eventobjchg.newobj) );
      }

      /* inform all existing primal solutions about the objective change */
      SCIPprimalUpdateVarObj(primal, var, event->data.eventobjchg.oldobj, event->data.eventobjchg.newobj);

      /* process variable's event filter */
      CHECK_OKAY( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
      var = event->data.eventbdchg.var;
      assert(var != NULL);
      assert(var->eventqueueindexlb == -1);

      /* inform LP about bound change and update branching candidates */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         assert(SCIPvarGetProbindex(var) >= 0);
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
         {
            CHECK_OKAY( SCIPcolChgLb(SCIPvarGetCol(var), set, lp, event->data.eventbdchg.newbound) );
         }
         CHECK_OKAY( SCIPlpUpdateVarLb(lp, set, var, event->data.eventbdchg.oldbound,
                        event->data.eventbdchg.newbound) );
         CHECK_OKAY( SCIPbranchcandUpdateVar(branchcand, set, var) );
      }

      /* process variable's event filter */
      CHECK_OKAY( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      var = event->data.eventbdchg.var;
      assert(var != NULL);
      assert(var->eventqueueindexub == -1);

      /* inform LP about bound change and update branching candidates */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         assert(SCIPvarGetProbindex(var) >= 0);
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
         {
            CHECK_OKAY( SCIPcolChgUb(SCIPvarGetCol(var), set, lp, event->data.eventbdchg.newbound) );
         }
         CHECK_OKAY( SCIPlpUpdateVarUb(lp, set, var, event->data.eventbdchg.oldbound, 
                        event->data.eventbdchg.newbound) );
         CHECK_OKAY( SCIPbranchcandUpdateVar(branchcand, set, var) );
      }

      /* process variable's event filter */
      CHECK_OKAY( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_HOLEADDED:
      errorMessage("HOLEADDED event not implemented yet\n");
      SCIPABORT();

   case SCIP_EVENTTYPE_HOLEREMOVED:
      errorMessage("HOLEREMOVED event not implemented yet\n");
      SCIPABORT();

   default:
      errorMessage("unknown event type <%d>\n", event->eventtype);
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}



/*
 * Event filter methods
 */

/** resizes eventfilter arrays to be able to store at least num entries */
static
RETCODE eventfilterEnsureMem(
   EVENTFILTER*     eventfilter,        /**< event filter */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(eventfilter != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);

   if( num > eventfilter->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &eventfilter->eventtypes, eventfilter->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &eventfilter->eventhdlrs, eventfilter->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &eventfilter->eventdatas, eventfilter->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &eventfilter->nextpos, eventfilter->size, newsize) );
      eventfilter->size = newsize;
   }
   assert(num <= eventfilter->size);
   
   return SCIP_OKAY;
}

/** creates an event filter */
RETCODE SCIPeventfilterCreate(
   EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   BLKMEM*          blkmem              /**< block memory buffer */
   )
{
   assert(eventfilter != NULL);
   assert(blkmem != NULL);

   ALLOC_OKAY( allocBlockMemory(blkmem, eventfilter) );
   (*eventfilter)->eventtypes = NULL;
   (*eventfilter)->eventhdlrs = NULL;
   (*eventfilter)->eventdatas = NULL;
   (*eventfilter)->nextpos = NULL;
   (*eventfilter)->size = 0;
   (*eventfilter)->len = 0;
   (*eventfilter)->firstfreepos = -1;
   (*eventfilter)->firstdeletedpos = -1;
   (*eventfilter)->eventmask = SCIP_EVENTTYPE_DISABLED;
   (*eventfilter)->delayedeventmask = SCIP_EVENTTYPE_DISABLED;
   (*eventfilter)->delayupdates = FALSE;

   return SCIP_OKAY;
}

/** frees an event filter and the associated event data entries */
RETCODE SCIPeventfilterFree(
   EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(eventfilter != NULL);
   assert(*eventfilter != NULL);
   assert(!(*eventfilter)->delayupdates);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   /* free event data */
   for( i = 0; i < (*eventfilter)->len; ++i )
   {
      if( (*eventfilter)->eventtypes[i] != SCIP_EVENTTYPE_DISABLED )
      {
         assert((*eventfilter)->eventhdlrs[i] != NULL);
         if( (*eventfilter)->eventhdlrs[i]->eventdelete != NULL )
         {
            CHECK_OKAY( (*eventfilter)->eventhdlrs[i]->eventdelete(set->scip, (*eventfilter)->eventhdlrs[i],
                  &(*eventfilter)->eventdatas[i]) );
         }
      }
   }

   /* free event filter data */
   freeBlockMemoryArrayNull(blkmem, &(*eventfilter)->eventtypes, (*eventfilter)->size);
   freeBlockMemoryArrayNull(blkmem, &(*eventfilter)->eventhdlrs, (*eventfilter)->size);
   freeBlockMemoryArrayNull(blkmem, &(*eventfilter)->eventdatas, (*eventfilter)->size);
   freeBlockMemoryArrayNull(blkmem, &(*eventfilter)->nextpos, (*eventfilter)->size);
   freeBlockMemory(blkmem, eventfilter);

   return SCIP_OKAY;
}

/** adds element to event filter */
RETCODE SCIPeventfilterAdd(
   EVENTFILTER*     eventfilter,        /**< event filter */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int*             filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   int pos;

   assert(eventfilter != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(eventhdlr != NULL);

   if( eventfilter->delayupdates )
   {
      /* insert addition to the end of the arrays;
       * in delayed addition we have to add to the end of the arrays, in order to not destroy the validity of the
       * arrays we are currently iterating over
       */
      CHECK_OKAY( eventfilterEnsureMem(eventfilter, blkmem, set, eventfilter->len + 1) );
      pos = eventfilter->len;
      eventfilter->len++;

      /* update delayed event filter mask */
      eventfilter->delayedeventmask |= eventtype;
   }
   else
   {
      if( eventfilter->firstfreepos == -1 )
      {
         /* insert addition to the end of the arrays */
         CHECK_OKAY( eventfilterEnsureMem(eventfilter, blkmem, set, eventfilter->len + 1) );
         pos = eventfilter->len;
         eventfilter->len++;
      }
      else
      {
         /* use the first free slot to store the added event filter entry */
         pos = eventfilter->firstfreepos;
         assert(0 <= pos && pos < eventfilter->len);
         assert(eventfilter->eventtypes[pos] == SCIP_EVENTTYPE_DISABLED);
         eventfilter->firstfreepos = eventfilter->nextpos[pos];
         assert(-1 <= eventfilter->firstfreepos && eventfilter->firstfreepos < eventfilter->len);
      }

      /* update event filter mask */
      eventfilter->eventmask |= eventtype;
   }
   assert(0 <= pos && pos < eventfilter->len);

   eventfilter->eventtypes[pos] = eventtype;
   eventfilter->eventhdlrs[pos] = eventhdlr;
   eventfilter->eventdatas[pos] = eventdata;
   eventfilter->nextpos[pos] = -2;

   if( filterpos != NULL )
      *filterpos = pos;

   return SCIP_OKAY;
}

/** linear search for the given entry in event filter */
static
int eventfilterSearch(
   EVENTFILTER*     eventfilter,        /**< event filter */
   EVENTTYPE        eventtype,          /**< event type */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   )
{
   int i;

   assert(eventfilter != NULL);
   assert(eventtype != SCIP_EVENTTYPE_DISABLED);
   assert(eventhdlr != NULL);

   for( i = 0; i < eventfilter->len; ++i )
   {
      if( eventhdlr == eventfilter->eventhdlrs[i]
         && eventdata == eventfilter->eventdatas[i]
         && eventtype == eventfilter->eventtypes[i]
         && eventfilter->nextpos[i] == -2 )
         return i;
   }
   
   return -1;
}

/** deletes element from event filter */
RETCODE SCIPeventfilterDel(
   EVENTFILTER*     eventfilter,        /**< event filter */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int              filterpos           /**< position of event filter entry, or -1 if unknown */
   )
{
   assert(eventfilter != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(eventtype != SCIP_EVENTTYPE_DISABLED);
   assert(eventhdlr != NULL);
   assert(-1 <= filterpos && filterpos < eventfilter->len);

   /* search position of event filter entry, if not given by the user */
   if( filterpos == -1 )
      filterpos = eventfilterSearch(eventfilter, eventtype, eventhdlr, eventdata);
   if( filterpos == -1 )
   {
      errorMessage("no event for event handler %p with data %p and event mask 0x%x found in event filter %p\n",
         eventhdlr, eventdata, eventtype, eventfilter);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= filterpos && filterpos < eventfilter->len);
   assert(eventfilter->eventtypes[filterpos] == eventtype);
   assert(eventfilter->eventhdlrs[filterpos] == eventhdlr);
   assert(eventfilter->eventdatas[filterpos] == eventdata);
   assert(eventfilter->nextpos[filterpos] == -2);

   /* if updates are delayed, insert entry into the list of delayed deletions;
    * otherwise, delete the entry from the filter directly and add the slot to the free list
    */
   if( eventfilter->delayupdates )
   {
      /* append filterpos to the list of deleted entries */
      eventfilter->nextpos[filterpos] = eventfilter->firstdeletedpos;
      eventfilter->firstdeletedpos = filterpos;
   }
   else
   {
      /* disable the entry in the filter and add the slot to the free list */
      assert(eventfilter->nextpos[filterpos] == -2);
      eventfilter->eventtypes[filterpos] = SCIP_EVENTTYPE_DISABLED;
      eventfilter->nextpos[filterpos] = eventfilter->firstfreepos;
      eventfilter->firstfreepos = filterpos;
   }

   return SCIP_OKAY;
}

/** makes the event filter to delay and buffer all updates until eventfilterProcessUpdates() is called */
static
void eventfilterDelayUpdates(
   EVENTFILTER*     eventfilter         /**< event filter */
   )
{
   assert(eventfilter != NULL);
   assert(!eventfilter->delayupdates);
   assert(eventfilter->delayedeventmask == SCIP_EVENTTYPE_DISABLED);

   eventfilter->delayupdates = TRUE;
}

/** processes all delayed additions and deletions */
static
void eventfilterProcessUpdates(
   EVENTFILTER*     eventfilter         /**< event filter */
   )
{
   int pos;
   int nextpos;

   assert(eventfilter != NULL);
   assert(eventfilter->delayupdates);

   /* move deleted entries into the free list and disable them */
   pos = eventfilter->firstdeletedpos;
   while( pos != -1 )
   {
      assert(0 <= pos && pos < eventfilter->len);
      assert(eventfilter->eventtypes[pos] != SCIP_EVENTTYPE_DISABLED);
      assert(eventfilter->nextpos[pos] >= -1);

      nextpos = eventfilter->nextpos[pos];
      eventfilter->nextpos[pos] = eventfilter->firstfreepos;
      eventfilter->firstfreepos = pos;
      eventfilter->eventtypes[pos] = SCIP_EVENTTYPE_DISABLED;
      pos = nextpos;
   }
   eventfilter->firstdeletedpos = -1;

   /* update event mask */
   eventfilter->eventmask |= eventfilter->delayedeventmask;
   eventfilter->delayedeventmask = SCIP_EVENTTYPE_DISABLED;

   /* mark the event filter updates to be no longer delayed */
   eventfilter->delayupdates = FALSE;
}

/** processes the event with all event handlers with matching filter setting */
RETCODE SCIPeventfilterProcess(
   EVENTFILTER*     eventfilter,        /**< event filter */
   SET*             set,                /**< global SCIP settings */
   EVENT*           event               /**< event to process */
   )
{
   EVENTTYPE eventtype;
   EVENTTYPE* eventtypes;
   Bool processed;
   int len;
   int i;

   assert(eventfilter != NULL);
   assert(set != NULL);
   assert(event != NULL);

   debugMessage("processing event filter %p (len %d, mask 0x%x) with event type 0x%x\n",
      eventfilter, eventfilter->len, eventfilter->eventmask, event->eventtype);

   eventtype = event->eventtype;

   /* check, if there may be any event handler for specific event */
   if( (eventtype & eventfilter->eventmask) == 0 )
      return SCIP_OKAY;

   /* delay the updates on this eventfilter, such that changes during event processing to the event filter
    * don't destroy necessary information of the arrays we are currently using
    */
   eventfilterDelayUpdates(eventfilter);

   /* process the event by calling the event handlers */
   processed = FALSE;
   len = eventfilter->len;
   eventtypes = eventfilter->eventtypes;
   for( i = 0; i < len; ++i )
   {
      /* check, if event is applicable for the filter element */
      if( (eventtype & eventtypes[i]) != 0 )
      {
         /* call event handler */
         CHECK_OKAY( SCIPeventhdlrExec(eventfilter->eventhdlrs[i], set, event, eventfilter->eventdatas[i]) );
         processed = TRUE;
      }
   }
      
   /* update eventfilter mask, if event was not processed by any event handler */
   if( !processed )
   {
      eventfilter->eventmask &= ~event->eventtype;
      debugMessage(" -> event type 0x%x not processed. new mask of event filter %p: 0x%x\n",
         event->eventtype, eventfilter, eventfilter->eventmask);
   }

   /* process delayed events on this eventfilter */
   eventfilterProcessUpdates(eventfilter);

   return SCIP_OKAY;
}



/*
 * Event queue methods
 */

/** resizes events array to be able to store at least num entries */
static
RETCODE eventqueueEnsureEventsMem(
   EVENTQUEUE*      eventqueue,         /**< event queue */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(eventqueue != NULL);
   assert(set != NULL);

   if( num > eventqueue->eventssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&eventqueue->events, newsize) );
      eventqueue->eventssize = newsize;
   }
   assert(num <= eventqueue->eventssize);
   
   return SCIP_OKAY;
}

/** creates an event queue */
RETCODE SCIPeventqueueCreate(
   EVENTQUEUE**     eventqueue          /**< pointer to store the event queue */
   )
{
   assert(eventqueue != NULL);

   ALLOC_OKAY( allocMemory(eventqueue) );
   (*eventqueue)->events = NULL;
   (*eventqueue)->eventssize = 0;
   (*eventqueue)->nevents = 0;
   (*eventqueue)->delayevents = FALSE;

   return SCIP_OKAY;
}

/** frees event queue; there must not be any unprocessed eventy in the queue! */
RETCODE SCIPeventqueueFree(
   EVENTQUEUE**     eventqueue          /**< pointer to the event queue */
   )
{
   assert(eventqueue != NULL);
   assert(*eventqueue != NULL);
   assert((*eventqueue)->nevents == 0);

   freeMemoryArrayNull(&(*eventqueue)->events);
   freeMemory(eventqueue);
   
   return SCIP_OKAY;
}

/** appends event to the event queue; sets event to NULL afterwards */
static
RETCODE eventqueueAppend(
   EVENTQUEUE*      eventqueue,         /**< event queue */
   SET*             set,                /**< global SCIP settings */
   EVENT**          event               /**< pointer to event to append to the queue */
   )
{
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

/** processes event or adds event to the event queue */
RETCODE SCIPeventqueueAdd(
   EVENTQUEUE*      eventqueue,         /**< event queue */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data; only needed for objchanged events */
   LP*              lp,                 /**< current LP data; only needed for obj/boundchanged events */
   BRANCHCAND*      branchcand,         /**< branching candidate storage; only needed for boundchange events */
   EVENTFILTER*     eventfilter,        /**< event filter for global events; not needed for variable specific events */
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
      CHECK_OKAY( SCIPeventProcess(*event, set, primal, lp, branchcand, eventfilter) );
      CHECK_OKAY( SCIPeventFree(event, blkmem) );
   }
   else
   {
      /* delay processing of event by appending it to the event queue */
      debugMessage("adding event %p of type 0x%x to event queue %p\n", *event, (*event)->eventtype, eventqueue);

      switch( (*event)->eventtype )
      {
      case SCIP_EVENTTYPE_DISABLED:
         errorMessage("cannot add a disabled event to the event queue\n");
         return SCIP_INVALIDDATA;

      case SCIP_EVENTTYPE_VARADDED:
      case SCIP_EVENTTYPE_VARDELETED:
      case SCIP_EVENTTYPE_VARFIXED:
      case SCIP_EVENTTYPE_LOCKSCHANGED:
      case SCIP_EVENTTYPE_PRESOLVEROUND:
      case SCIP_EVENTTYPE_NODEFOCUSED:
      case SCIP_EVENTTYPE_NODEFEASIBLE:
      case SCIP_EVENTTYPE_NODEINFEASIBLE:
      case SCIP_EVENTTYPE_NODEBRANCHED:
      case SCIP_EVENTTYPE_FIRSTLPSOLVED:
      case SCIP_EVENTTYPE_LPSOLVED:
      case SCIP_EVENTTYPE_POORSOLFOUND:
      case SCIP_EVENTTYPE_BESTSOLFOUND:
         /* these events cannot (or need not) be merged; just add them to the queue */
         CHECK_OKAY( eventqueueAppend(eventqueue, set, event) );
         break;

      case SCIP_EVENTTYPE_OBJCHANGED:
         /* changes in objective value may be merged with older changes in objective value */
         var = (*event)->data.eventobjchg.var;
         assert(var != NULL);
         pos = var->eventqueueindexobj;
         if( pos >= 0 )
         {
            /* the objective value change event already exists -> modifiy it accordingly */
            assert(pos < eventqueue->nevents);
            qevent = eventqueue->events[pos];
            assert(qevent != NULL);
            assert(qevent->eventtype == SCIP_EVENTTYPE_OBJCHANGED);
            assert(qevent->data.eventobjchg.var == var);
            assert(SCIPsetIsEQ(set, (*event)->data.eventobjchg.oldobj, qevent->data.eventobjchg.newobj));

            debugMessage(" -> merging OBJ event (<%s>,%g -> %g) with event at position %d (<%s>,%g -> %g)\n",
               SCIPvarGetName((*event)->data.eventobjchg.var), (*event)->data.eventobjchg.oldobj,
               (*event)->data.eventobjchg.newobj,
               pos, SCIPvarGetName(qevent->data.eventobjchg.var), qevent->data.eventobjchg.oldobj, 
               qevent->data.eventobjchg.newobj);

            qevent->data.eventobjchg.newobj = (*event)->data.eventobjchg.newobj;
            if( qevent->data.eventobjchg.newobj == qevent->data.eventobjchg.oldobj ) /*lint !e777*/
            {
               /* the queued objective value change was reversed -> disable the event in the queue */
               eventDisable(qevent);
               var->eventqueueindexobj = -1;
               debugMessage(" -> event disabled\n");
            }

            /* free the event that is of no use any longer */
            CHECK_OKAY( SCIPeventFree(event, blkmem) );
         }
         else
         {
            /* the objective value change event doesn't exist -> add it to the queue, and remember the array index */
            var->eventqueueindexobj = eventqueue->nevents;
            CHECK_OKAY( eventqueueAppend(eventqueue, set, event) );
         }
         break;

      case SCIP_EVENTTYPE_LBTIGHTENED:
      case SCIP_EVENTTYPE_LBRELAXED:
         /* changes in lower bound may be merged with older changes in lower bound */
         var = (*event)->data.eventbdchg.var;
         assert(var != NULL);
         pos = var->eventqueueindexlb;
         if( pos >= 0 )
         {
            /* the lower bound change event already exists -> modifiy it accordingly */
            assert(pos < eventqueue->nevents);
            qevent = eventqueue->events[pos];
            assert(qevent != NULL);
            assert(qevent->eventtype == SCIP_EVENTTYPE_LBTIGHTENED || qevent->eventtype == SCIP_EVENTTYPE_LBRELAXED);
            assert(qevent->data.eventbdchg.var == var);
            assert(SCIPsetIsEQ(set, (*event)->data.eventbdchg.oldbound, qevent->data.eventbdchg.newbound));

            debugMessage(" -> merging LB event (<%s>,%g -> %g) with event at position %d (<%s>,%g -> %g)\n",
               SCIPvarGetName((*event)->data.eventbdchg.var), (*event)->data.eventbdchg.oldbound,
               (*event)->data.eventbdchg.newbound,
               pos, SCIPvarGetName(qevent->data.eventbdchg.var), qevent->data.eventbdchg.oldbound, 
               qevent->data.eventbdchg.newbound);

            qevent->data.eventbdchg.newbound = (*event)->data.eventbdchg.newbound;
            /*if( SCIPsetIsLT(set, qevent->data.eventbdchg.newbound, qevent->data.eventbdchg.oldbound) )*/
            if( qevent->data.eventbdchg.newbound < qevent->data.eventbdchg.oldbound )
               qevent->eventtype = SCIP_EVENTTYPE_LBRELAXED;
            /*else if( SCIPsetIsGT(set, qevent->data.eventbdchg.newbound, qevent->data.eventbdchg.oldbound) )*/
            else if( qevent->data.eventbdchg.newbound > qevent->data.eventbdchg.oldbound )
               qevent->eventtype = SCIP_EVENTTYPE_LBTIGHTENED;
            else
            {
               /* the queued bound change was reversed -> disable the event in the queue */
               assert(qevent->data.eventbdchg.newbound == qevent->data.eventbdchg.oldbound); /*lint !e777*/
               eventDisable(qevent);
               var->eventqueueindexlb = -1;
               debugMessage(" -> event disabled\n");
            }

            /* free the event that is of no use any longer */
            CHECK_OKAY( SCIPeventFree(event, blkmem) );
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
         var = (*event)->data.eventbdchg.var;
         assert(var != NULL);
         pos = var->eventqueueindexub;
         if( pos >= 0 )
         {
            /* the upper bound change event already exists -> modifiy it accordingly */
            assert(pos < eventqueue->nevents);
            qevent = eventqueue->events[pos];
            assert(qevent != NULL);
            assert(qevent->eventtype == SCIP_EVENTTYPE_UBTIGHTENED || qevent->eventtype == SCIP_EVENTTYPE_UBRELAXED);
            assert(qevent->data.eventbdchg.var == var);
            assert(SCIPsetIsEQ(set, (*event)->data.eventbdchg.oldbound, qevent->data.eventbdchg.newbound));

            debugMessage(" -> merging UB event (<%s>,%g -> %g) with event at position %d (<%s>,%g -> %g)\n",
               SCIPvarGetName((*event)->data.eventbdchg.var), (*event)->data.eventbdchg.oldbound,
               (*event)->data.eventbdchg.newbound,
               pos, SCIPvarGetName(qevent->data.eventbdchg.var), qevent->data.eventbdchg.oldbound,
               qevent->data.eventbdchg.newbound);

            qevent->data.eventbdchg.newbound = (*event)->data.eventbdchg.newbound;
            /*if( SCIPsetIsLT(set, qevent->data.eventbdchg.newbound, qevent->data.eventbdchg.oldbound) )*/
            if( qevent->data.eventbdchg.newbound < qevent->data.eventbdchg.oldbound )
               qevent->eventtype = SCIP_EVENTTYPE_UBTIGHTENED;
            /*else if( SCIPsetIsGT(set, qevent->data.eventbdchg.newbound, qevent->data.eventbdchg.oldbound) )*/
            else if( qevent->data.eventbdchg.newbound > qevent->data.eventbdchg.oldbound )
               qevent->eventtype = SCIP_EVENTTYPE_UBRELAXED;
            else
            {
               /* the queued bound change was reversed -> disable the event in the queue */
               assert(qevent->data.eventbdchg.newbound == qevent->data.eventbdchg.oldbound); /*lint !e777*/
               eventDisable(qevent);
               var->eventqueueindexub = -1;
               debugMessage(" -> event disabled\n");
            }

            /* free the event that is of no use any longer */
            CHECK_OKAY( SCIPeventFree(event, blkmem) );
         }
         else
         {
            /* the upper bound change event doesn't exist -> add it to the queue, and remember the array index */
            var->eventqueueindexub = eventqueue->nevents;
            CHECK_OKAY( eventqueueAppend(eventqueue, set, event) );
         }
         break;

      case SCIP_EVENTTYPE_HOLEADDED:
         errorMessage("HOLEADDED event not implemented yet\n");
         SCIPABORT();

      case SCIP_EVENTTYPE_HOLEREMOVED:
         errorMessage("HOLEREMOVED event not implemented yet\n");
         SCIPABORT();

      default:
         errorMessage("unknown event type <%d>\n", (*event)->eventtype);
         return SCIP_INVALIDDATA;
      }
   }
   
   assert(*event == NULL);

   return SCIP_OKAY;
}

/** marks queue to delay incoming events until a call to SCIPeventqueueProcess() */
RETCODE SCIPeventqueueDelay(
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(eventqueue != NULL);
   assert(!eventqueue->delayevents);

   debugMessage("event processing is delayed\n");

   eventqueue->delayevents = TRUE;

   return SCIP_OKAY;
}

/** processes all delayed events, marks queue to process events immediately */
RETCODE SCIPeventqueueProcess(
   EVENTQUEUE*      eventqueue,         /**< event queue */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTFILTER*     eventfilter         /**< event filter for global (not variable dependent) events */
   )
{
   EVENT* event;
   int i;

   assert(eventqueue != NULL);
   assert(eventqueue->delayevents);

   debugMessage("processing %d queued events\n", eventqueue->nevents);

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

      /* unmark the event queue index of a variable with changed objective value or bounds */
      if( (event->eventtype & SCIP_EVENTTYPE_OBJCHANGED) != 0 )
      {
         assert(event->data.eventobjchg.var->eventqueueindexobj == i);
         event->data.eventobjchg.var->eventqueueindexobj = -1;
      }
      else if( (event->eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
      {
         assert(event->data.eventbdchg.var->eventqueueindexlb == i);
         event->data.eventbdchg.var->eventqueueindexlb = -1;
      }
      else if( (event->eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0 )
      {
         assert(event->data.eventbdchg.var->eventqueueindexub == i);
         event->data.eventbdchg.var->eventqueueindexub = -1;
      }

      /* process event */
      CHECK_OKAY( SCIPeventProcess(event, set, primal, lp, branchcand, eventfilter) );

      /* free the event immediately, because additionally raised events during event processing
       * can lead to a large event queue
       */
      CHECK_OKAY( SCIPeventFree(&eventqueue->events[i], blkmem) );
   }

   assert(i == eventqueue->nevents);
   eventqueue->nevents = 0;
   eventqueue->delayevents = FALSE;

   return SCIP_OKAY;
}

/** returns TRUE iff events of the queue are delayed until the next SCIPeventqueueProcess() call */
Bool SCIPeventqueueIsDelayed(
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(eventqueue != NULL);

   return eventqueue->delayevents;
}
