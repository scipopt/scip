/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_nodereopt.c
 * @brief  eventhdlr for nodereopt event
 * @author Jakob Witzig
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include "scip/branch_nodereopt.h"
#include "scip/branch_pseudo.h"
#include "scip/event_nodereopt.h"
#include "scip/pub_event.h"
#include "scip/tree.h"
#include <string.h>

#define EVENTHDLR_NAME         "nodereopt"
#define EVENTHDLR_DESC         "event handler for nodereopt event"

//#define DEBUG_MODE

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Longint          highestseennodenr;       /** highest node number of a seen node in this round */
   int                   lastnodenr;              /** number of last seen node */
   int                   lastbranched;            /** number of last branched node */
   int                   allocseennodes;          /** site of allocated memory for seen nodes */
   int*                  seennodes;               /** array of already seen nodes and their event, 0: not seen yet, 1: branched, 2: infeasible, 3: feasible */
   SCIP_EVENTTYPE        lasteventtype;           /** last caught eventtype */
   SCIP_Bool             reopt;                   /** is reoptimization enabled */
   SCIP_Bool             strbrinleafs;            /** save strong branching information in leaf nodes. */
   SCIP_Bool             init;                    /** data is initialized */
};

void SCIPeventhdlrNodereoptDisable(
   SCIP*                 scip
)
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->init = FALSE;
   eventhdlrdata->reopt = FALSE;

   return;
}

/*
 * static methods
 */

/** Check the reason of the cutoff */
static SCIP_RETCODE
checkCutoffReason(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata,
   SCIP_EVENT*           event
)
{
   SCIP_NODE* eventnode;
   SCIP_LPSOLSTAT solstat;
   SCIP_EVENTTYPE eventtype;
   SCIP_Bool strongbranched;

   assert(scip != NULL );
   assert(eventhdlrdata->reopt);

   eventnode = SCIPeventGetNode(event);
   solstat = SCIPgetLPSolstat(scip);
   eventtype = SCIPeventGetType(event);
   strongbranched = FALSE;

   assert(eventnode != NULL );

   if (eventtype == SCIP_EVENTTYPE_NODEBRANCHED
    || eventtype == SCIP_EVENTTYPE_NODEFEASIBLE
    || eventtype == SCIP_EVENTTYPE_NODEINFEASIBLE)
   {
      if (SCIPgetEffectiveRootDepth(scip) == SCIPnodeGetDepth(eventnode))
      {
         strongbranched = SCIPbranchrulePseudoIsPseudoBranched(scip, eventnode);
      }
      else
      {
         strongbranched = SCIPnodeGetNPseudoBranchings(eventnode) > 0 ? TRUE : FALSE;
      }
   }

   SCIPdebugMessage("check the reason of cutoff for node %lld:\n", SCIPnodeGetNumber(eventnode));
   SCIPdebugMessage(" -> focusnode: %u\n", SCIPgetCurrentNode(scip) == eventnode);
   SCIPdebugMessage(" -> depth: %d, eff. root depth: %d\n", SCIPnodeGetDepth(eventnode), SCIPgetEffectiveRootDepth(scip));
   SCIPdebugMessage(" -> strong branched: %u\n", strongbranched);
   SCIPdebugMessage(" -> LP solstat     : %d\n", solstat);

   /** three different cases are possible: NODEBRANCHED, NODEFEASIBLE, NODEINFEASIBLE */
   switch (eventtype)
      {
   case SCIP_EVENTTYPE_NODEBRANCHED:
      /** current node has to be the eventnode */
      assert(SCIPgetCurrentNode(scip) == eventnode);

      eventhdlrdata->lastbranched = SCIPnodeGetNumber(SCIPeventGetNode(event));

      /**
       * we have to check the depth of the current node. if the depth is equal to the effective
       * root depth, then all information about pseudo-branchings already exists, else we have to
       * look at the domchg-data-structure to check if the node is pseudo-branched.
       */
      if (SCIPnodeGetDepth(eventnode) == SCIPgetEffectiveRootDepth(scip))
      {
//         /** if the node is pseudo-branched, finish the collecting of pseudo-data for this node and save them */
//         if (strongbranched || SCIPbranchruleNodereoptGetNAddedConss(scip, eventnode) > 0)
//         {
            /*
             * Save the node if there are added constraints, because this means the node a copy of pseudo-branched node
             * and contains a pseudo logic-or-constraint
             */
         if( strongbranched )
         {
            SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_STRBRANCHED);
            SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_STRBRANCHED);

            SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_STRBRANCHED) );
         }
         else if( SCIPbranchruleNodereoptGetNAddedConss(scip, eventnode) > 0 )
         {
            SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_LOGICORNODE);

            SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_LOGICORNODE) );
         }
         else
         {
            SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_TRANSIT);

            SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_TRANSIT) );
         }
      }
      else
      {
         /**
          * we only branch on binary variables and a -1 indicates memory allocation w/o saving information.
          *
          * we have to do this in the following order:
          * 1) all bound-changes are local, thats way we have to mark the node as pseudo-branched for branch_pseudo
          *    by adding a NULL variable
          * 2) save ancestor-branchings before finishing the node, because the check if a node is pseudo-branched in
          *    branch_pseudo only looks at the last seen node, NodeFinished() will clear the pointer to the last seen node.
          * 3) call NodeFinished(); reset all intern pointers, erase number of pseudo-branched nodes
          */
         if( strongbranched )
         {
            SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_STRBRANCHED);
            SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_STRBRANCHED);

            SCIP_CALL( SCIPbranchrulePseudoAddPseudoVar(scip, eventnode, NULL, SCIP_BOUNDTYPE_LOWER, -1) );
            SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_STRBRANCHED) );
         }
         else if( SCIPbranchruleNodereoptGetNAddedConss(scip, eventnode) > 0 )
         {
            SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_LOGICORNODE);

            SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_LOGICORNODE) );
         }
         else
         {
            SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_TRANSIT);

            SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_TRANSIT) );
         }
      }
      break;

   case SCIP_EVENTTYPE_NODEFEASIBLE:
      /** current node has to be the eventnode */
      assert(SCIPgetCurrentNode(scip) == eventnode);

      /** no after-branch heuristic should throw a NODEFEASIBLE event */
      if (SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == eventhdlrdata->lastbranched)
      {
         printf("new event %x for node %llu\n", eventhdlrdata->lasteventtype,
               SCIPnodeGetNumber(eventnode));
         assert(eventhdlrdata->lastbranched != SCIPnodeGetNumber(eventnode));
      }

      /** we can delete all pseudo information; because we have to revive this node and we don't wont to split up */
      if( strongbranched )
      {
         if(SCIPnodeGetDepth(eventnode) == SCIPgetEffectiveRootDepth(scip) )
         {
            if( !eventhdlrdata->strbrinleafs )
            {
               SCIP_CALL(SCIPbranchrulePseudoDeleteLastNodeInfo(scip, eventnode));

               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_FEASIBLE);
               SCIP_CALL(SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_FEASIBLE));
            }
            else
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_STRBRANCHED);
               SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_STRBRANCHED);

               SCIP_CALL(SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_STRBRANCHED));
            }
         }
         else if(SCIPnodeGetDepth(eventnode) > SCIPgetEffectiveRootDepth(scip) )
         {
            if( !eventhdlrdata->strbrinleafs )
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_FEASIBLE);
               SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_FEASIBLE) );
            }
            else
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_STRBRANCHED);
               SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_STRBRANCHED);

               SCIP_CALL( SCIPbranchrulePseudoAddPseudoVar(scip, eventnode, NULL, SCIP_BOUNDTYPE_LOWER, -1) );
               SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_STRBRANCHED) );
            }
         }
      }
      else
      {
         /* if the node was created by branch_nodereopt, only the LP basis will be refreshed (if enabled) */
         SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_FEASIBLE);
         SCIP_CALL(SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_FEASIBLE));
      }

      break;

   case SCIP_EVENTTYPE_NODEINFEASIBLE:
      /**
       * we have to check if the current node is the event node.
       * if the current node is not the event node, we have to save this node, else we have to
       * look at LP solstat and decide.
       */
      if (SCIPgetCurrentNode(scip) == eventnode)
      {
         /**
          * an after-branch heuristic says NODEINFEASIBLE, maybe the cutoff bound is reached.
          * because the node is already branched we have all children and can delete this node.
          */
         if (SCIPnodeGetNumber(eventnode) == eventhdlrdata->lastbranched)
            break;

         /*
          * if the node is strong branched we possible detects an infeasible subtree, if not,
          * the whole node is either infeasible or exceeds the cutoff bound.
          */
         if( strongbranched )
         {
            /*
             * 1. the LP is not solved or infeasible: the subnode is infeasible and can be discarded
             *    because either the LP proves infeasibility or a constraint handler.
             *    We have to store an infeasible subtree constraint
             * 2. the LP exceeds the objective limit, we have to store the node and can delete the
             *    strong branching information
             */
            if( solstat == SCIP_LPSOLSTAT_INFEASIBLE )
            {
               /* add a dummy variable, because the bound changes were not global in the
                * sense of effective root depth */
               if( SCIPnodeGetDepth(eventnode) > SCIPgetEffectiveRootDepth(scip) )
               {
                  SCIP_CALL( SCIPbranchrulePseudoAddPseudoVar(scip, eventnode, NULL, SCIP_BOUNDTYPE_LOWER, -1) );
               }

               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_INFSUBTREE);
               SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_INFSUBTREE);

               /* save the node as a strong branched node */
               SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_INFSUBTREE) );
            }
            else
            {
               assert(SCIP_LPSOLSTAT_OBJLIMIT || SCIP_LPSOLSTAT_OPTIMAL || SCIP_LPSOLSTAT_NOTSOLVED);

               /* delete strong branching information of some exists */
               if( SCIPnodeGetDepth(eventnode) <= SCIPgetEffectiveRootDepth(scip))
               {
                  if( !eventhdlrdata->strbrinleafs )
                  {
                     SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_PRUNED);

                     SCIP_CALL(SCIPbranchrulePseudoDeleteLastNodeInfo(scip, eventnode));
                     SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_PRUNED) );
                  }
                  else
                  {
                     SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_STRBRANCHED);
                     SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_STRBRANCHED);

                     SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_STRBRANCHED) );
                  }
               }
               else
               {
                  if( !eventhdlrdata->strbrinleafs )
                  {
                     SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_PRUNED);
                     SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_PRUNED) );
                  }
                  else
                  {
                     SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_STRBRANCHED);
                     SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_STRBRANCHED);

                     SCIP_CALL( SCIPbranchrulePseudoAddPseudoVar(scip, eventnode, NULL, SCIP_BOUNDTYPE_LOWER, -1) );
                     SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_STRBRANCHED) );
                  }
               }
            }
         }
         else
         {
            /*
             * 1. the LP is not solved or infeasible: the whole node is infeasible and can be discarded
             *    because either the LP proves infeasibility or a constraint handler.
             * 2. the LP exceeds the objective limit, we have to store the node and can delete the
             *    strong branching information
             */
            if( solstat == SCIP_LPSOLSTAT_INFEASIBLE )
            {
               /* save the information of an infeasible node */
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_INFEASIBLE);
               SCIP_CALL(SCIPbranchruleNodereoptInfNode(scip, eventnode));
            }
            else
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_PRUNED);

               /* store the node */
               SCIP_CALL(SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_PRUNED));
            }
         }
      }
      else
      {
         SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_PRUNED);

         /* if the node was created by branch_nodereopt, nothing happens */
         SCIP_CALL(SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_PRUNED));

      }

      break;

   default:
      SCIPerrorMessage( "the node event is neither NODEBRANCHED, NODEFEASIBLE nor NODEINFEASIBLE\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** reset the data structure */
static
SCIP_RETCODE Reset(
   SCIP*                 scip,                     /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata
)
{
   assert(scip != NULL );
   assert(eventhdlrdata != NULL );

   BMSclearMemoryArray(eventhdlrdata->seennodes, eventhdlrdata->allocseennodes);

   eventhdlrdata->highestseennodenr = -1;
   eventhdlrdata->lastnodenr = -1;
   eventhdlrdata->lastbranched = -1;
   eventhdlrdata->lasteventtype = 0;

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 */

#define eventCopyNodereopt NULL;
#define eventExitsolNodereopt NULL;
#define eventDeleteNodereopt NULL;

static
SCIP_DECL_EVENTINIT(eventInitNodereopt)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL );
   assert(eventhdlr != NULL );
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );

   /** HACK */
   /** check if all variable are binary, if not, disable reoptimization */
   if( !eventhdlrdata->init )
   {
      int maxsavednodes;

      SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/enable", &eventhdlrdata->reopt) );
      if( eventhdlrdata->reopt && SCIPgetNImplVars(scip) + SCIPgetNIntVars(scip) > 0 )
         eventhdlrdata->reopt = FALSE;

      SCIP_CALL( SCIPgetIntParam(scip, "reoptimization/maxsavednodes", &maxsavednodes) );

      if( maxsavednodes == 0 )
         eventhdlrdata->reopt = FALSE;
   }

   if( !eventhdlrdata->init && eventhdlrdata->reopt )
   {
      eventhdlrdata->allocseennodes = 500;
      eventhdlrdata->highestseennodenr = eventhdlrdata->allocseennodes;

      SCIP_CALL( SCIPallocMemoryArray(scip, &(eventhdlrdata->seennodes), eventhdlrdata->allocseennodes) );
      SCIP_CALL( Reset(scip, eventhdlrdata) );

      eventhdlrdata->init = TRUE;
   }

   return SCIP_OKAY;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeNodereopt)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL );
   assert(eventhdlr != NULL );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );

   if( eventhdlrdata->init )
   {
      SCIPfreeMemoryArray(scip, &eventhdlrdata->seennodes);
      eventhdlrdata->init = FALSE;
      eventhdlrdata->reopt = FALSE;
   }

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitNodereopt)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL );
   assert(eventhdlr != NULL );
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );

   if(eventhdlrdata->init && eventhdlrdata->reopt)
   {
      SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEINFEASIBLE | SCIP_EVENTTYPE_NODEFEASIBLE | SCIP_EVENTTYPE_NODEBRANCHED, eventhdlr, NULL, -1));
   }

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolNodereopt)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL );
   assert(eventhdlr != NULL );
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );

   if(eventhdlrdata->init && eventhdlrdata->reopt)
   {
      SCIP_CALL( Reset(scip, eventhdlrdata) );
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEINFEASIBLE | SCIP_EVENTTYPE_NODEFEASIBLE | SCIP_EVENTTYPE_NODEBRANCHED, eventhdlr, NULL, NULL) );
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecNodereopt)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_NODE* eventnode;

   assert(scip != NULL );
   assert(eventhdlr != NULL );
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );

   eventnode = SCIPeventGetNode(event);

   if (!eventhdlrdata->reopt)
      return SCIP_OKAY;

   /* return of the event node is a probing or refocus node */
   if( SCIPnodeGetType(eventnode) == SCIP_NODETYPE_PROBINGNODE || SCIPnodeGetType(eventnode) == SCIP_NODETYPE_REFOCUSNODE )
      return SCIP_OKAY;

   if( eventhdlrdata->lastnodenr == SCIPnodeGetNumber(eventnode) )
   {
      switch ( SCIPeventGetType(event) ) {
         case SCIP_EVENTTYPE_NODEFEASIBLE:
            if( eventhdlrdata->lasteventtype != SCIP_EVENTTYPE_NODEFEASIBLE )
            {
               SCIPdebugMessage("node %lld was already seen as <%s>, new event is <SCIP_EVENTTYPE_FEASIBLE>.\n",
                     SCIPnodeGetNumber(eventnode),
                     eventhdlrdata->lasteventtype == SCIP_EVENTTYPE_NODEBRANCHED ? "SCIP_EVENTTYPE_NODEBRANCHED" : "SCIP_EVENTTYPE_NODEINFEASIBLE");
               assert(eventhdlrdata->lasteventtype == SCIP_EVENTTYPE_NODEFEASIBLE);
            }
            break;

         case SCIP_EVENTTYPE_NODEBRANCHED:
            SCIPdebugMessage("node %lld was already seen as <%s>, new event is <SCIP_EVENTTYPE_BRANCHED>.\n",
                  SCIPnodeGetNumber(eventnode),
                     eventhdlrdata->lasteventtype == SCIP_EVENTTYPE_NODEFEASIBLE ? "SCIP_EVENTTYPE_NODEFEASIBLE" : "SCIP_EVENTTYPE_NODEINFEASIBLE");
            assert(SCIPeventGetType(event) != SCIP_EVENTTYPE_NODEFEASIBLE);
            break;

         case SCIP_EVENTTYPE_NODEINFEASIBLE:
            if( eventhdlrdata->lasteventtype == SCIP_EVENTTYPE_NODEBRANCHED)
            {
               SCIPdebugMessage("node %lld was already seen as <%s>, new event is <SCIP_EVENTTYPE_NODEBRANCHED>.\n",
                     SCIPnodeGetNumber(eventnode),
                     eventhdlrdata->lasteventtype == SCIP_EVENTTYPE_NODEINFEASIBLE ? "SCIP_EVENTTYPE_NODEINFEASIBLE" : "SCIP_EVENTTYPE_NODEBRANCHED");
               assert(eventhdlrdata->lasteventtype == SCIP_EVENTTYPE_NODEBRANCHED);
            }
            break;

         default:
            break;
      }

      return SCIP_OKAY;
   }

   SCIPdebugMessage("catch event %x for node %lld\n", SCIPeventGetType(event), SCIPnodeGetNumber(SCIPeventGetNode(event)));

   eventhdlrdata->lastnodenr = SCIPnodeGetNumber(eventnode);
   eventhdlrdata->lasteventtype = SCIPeventGetType(event);

   /** if the current node is the root, than clear the data structure, set the new event and return */
   if (SCIPgetRootNode(scip) == eventnode)
   {
      if (SCIPbranchrulePseudoIsPseudoBranched(scip, eventnode))
      {
         SCIP_CALL(checkCutoffReason(scip, eventhdlrdata, event));
      }
      else if (SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEBRANCHED)
      {
         SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_TRANSIT) );
      }
      else if (SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEFEASIBLE)
      {
         SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_FEASIBLE) );
      }
      else if (SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEINFEASIBLE
            && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE)
         assert(FALSE);
      else if (SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEINFEASIBLE)
      {
         SCIP_CALL( SCIPbranchruleNodereoptAddNode(scip, eventnode, SCIP_REOPTTYPE_PRUNED) );
      }

      return SCIP_OKAY;
   }


   /** check if we have already seen this node */
   if (SCIPgetRootNode(scip) == eventnode && SCIPgetNChildren(scip) > 0
    && SCIPeventGetType(event) != SCIP_EVENTTYPE_NODEBRANCHED)
   {
      /** the node was already seen as branched node and the only allowed event is NODEINFEASIBLE
       * we don't save this node because we have already saved the children.
       */
      if (SCIPnodeGetReopttype(eventnode) <= SCIP_REOPTTYPE_STRBRANCHED)
      {
         assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEINFEASIBLE);
         return SCIP_OKAY;
      }
      else
      {
         assert(SCIPnodeGetReopttype(eventnode) <= SCIP_REOPTTYPE_STRBRANCHED);
      }
   }

   checkCutoffReason(scip, eventhdlrdata, event);

   return SCIP_OKAY;
}

#ifdef DEBUG_MODE
/*
 * Check the consistency of reoptimization nodes. Each node has to be seen once.
 */
SCIP_Bool SCIPeventhdlrNodereoptCheckConsistency(
   SCIP*                 scip,
   SCIP_Longint          ncreatednodes
)
{
   switch (SCIPgetStatus(scip))
      {
   case SCIP_STATUS_USERINTERRUPT:
   case SCIP_STATUS_NODELIMIT:
   case SCIP_STATUS_TOTALNODELIMIT:
   case SCIP_STATUS_STALLNODELIMIT:
   case SCIP_STATUS_TIMELIMIT:
   case SCIP_STATUS_MEMLIMIT:
   case SCIP_STATUS_GAPLIMIT:
   case SCIP_STATUS_SOLLIMIT:
   case SCIP_STATUS_BESTSOLLIMIT:
      return TRUE;
      break;
   default:
      break;
      }

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;
   int i;

   assert(scip != NULL );

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );

   for (i = 0; i < ncreatednodes; ++i)
   {
      if (eventhdlrdata->seennodes[i] == 0)
      {
         printf("Node %u not seen!\n", i + 1);
         return FALSE;
      }
   }

   return TRUE;
}
#endif

/** creates event handler for nodereopt event */
SCIP_RETCODE
SCIPincludeEventHdlrNodereopt(
   SCIP*                 scip                     /**< SCIP data structure */
)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   eventhdlrdata = NULL;

   /** create nodereopt event handler data */
   SCIP_CALL(SCIPallocMemory(scip, &eventhdlrdata));

   eventhdlr = NULL;

   SCIP_CALL(SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecNodereopt, eventhdlrdata));

   assert(eventhdlr != NULL );

   /** set non fundamental callbacks via setter functions */
   SCIP_CALL(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitNodereopt));
   SCIP_CALL(SCIPsetEventhdlrExit(scip, eventhdlr, eventExitNodereopt));
   SCIP_CALL(SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolNodereopt));
   SCIP_CALL(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeNodereopt));

   eventhdlrdata->init = FALSE;


   SCIP_CALL( SCIPaddBoolParam(scip, "reoptimization/strbrinleafs", "store determined strong branching information in leaf nodes.",
         &eventhdlrdata->strbrinleafs, TRUE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}

