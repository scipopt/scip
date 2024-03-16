/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_shadowtree.c
 * @ingroup DEFPLUGINS_EVENT
 * @brief  event handler for maintaining the unmodified branch-and-bound tree
 * @author Jasper van Doornmalen
 *
 * It is possible that SCIP detects that variable bounds can be restricted globally further than formerly known.
 * In that case, it is decided to update the global bounds of these variables, and modify the history of the branching
 * decisions this way. This breaks methods that depend on the assumption that historic choices in the branch-and-bound
 * tree remain unmodified througout the search, e.g., dynamic symmetry handling constraints.
 *
 * This event handler registers decisions made by the branch-and-bound tree directly at the moment of branching, and
 * does not modify those at later stages of the solve.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/debug.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/struct_var.h"
#include "scip/type_var.h"
#include "scip/scip.h"
#include "scip/scip_branch.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include "scip/struct_scip.h"
#include "scip/struct_mem.h"
#include "scip/struct_tree.h"
#include "scip/symmetry.h"
#include <ctype.h>
#include <string.h>
#include <memory.h>
#include "scip/event_shadowtree.h"

#define EVENTHDLR_NAME         "event_shadowtree"
#define EVENTHDLR_DESC         "event handler for maintaining the unmodified branch-and-bound tree"
#define NODEMAP_MAX_INITIAL_SIZE 10000
#define NODEMAP_MAX_INITIAL_SIZE_2LOG 14


/*
 * Data structures
 */


/** wrapper for shadow tree eventhandler data */
struct SCIP_EventhdlrData
{
#ifndef NDEBUG
   SCIP*                 scip;               /**< SCIP data structure */
#endif
   SCIP_SHADOWTREE*      shadowtree;         /**< Shadow tree structure */
   SCIP_CLOCK*           clock;              /**< clock for measuring time in shadow tree events */
   SCIP_Bool             active;             /**< whether a shadow tree should be maintained */
};


/*
 * Local methods
 */

/** hash key for SCIP_SHADOWNODE */
static
SCIP_DECL_HASHGETKEY(hashGetKeyShadowNode)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the indices of both node numbers are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqShadowNode)
{  /*lint --e{715}*/
   return ((SCIP_SHADOWNODE*) key1)->nodeid == ((SCIP_SHADOWNODE*) key2)->nodeid;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValShadowNode)
{  /*lint --e{715}*/
   return (unsigned int) ((SCIP_SHADOWNODE*) key)->nodeid;
}


/** get the time spent in the shadow tree eventhdlr */
SCIP_Real SCIPgetShadowTreeEventHandlerExecutionTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = (SCIP_EVENTHDLRDATA*) SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != NULL );
   assert( eventhdlrdata->scip != NULL );
   assert( eventhdlrdata->scip == scip );
   assert( eventhdlrdata->clock != NULL );

   return SCIPgetClockTime(scip, eventhdlrdata->clock);
}


/** given a node number, returns the node in the shadow tree, or NULL if it doesn't exist */
SCIP_SHADOWNODE* SCIPshadowTreeGetShadowNodeFromNodeNumber(
   SCIP_SHADOWTREE*      shadowtree,         /**< pointer to the shadow tree */
   SCIP_Longint          nodeid              /**< index of the node, equivalent to the standard branch and bound tree */
)
{
   SCIP_SHADOWNODE tmpnode;

   assert( shadowtree != NULL );
   assert( nodeid >= 0 );

   tmpnode.nodeid = nodeid;

   /* the following line of code returns NULL if it cannot find the entry in the hashtable */
   return (SCIP_SHADOWNODE*) SCIPhashtableRetrieve(shadowtree->nodemap, (void*) &tmpnode);
}

/** given a node, returns the node in the shadowtree, or NULL if it doesn't exist */
SCIP_SHADOWNODE* SCIPshadowTreeGetShadowNode(
   SCIP_SHADOWTREE*      shadowtree,         /**< pointer to the shadow tree */
   SCIP_NODE*            node                /**< node from the actual branch-and-bound tree */
)
{
   assert( shadowtree != NULL );
   assert( node != NULL );

   return SCIPshadowTreeGetShadowNodeFromNodeNumber(shadowtree, SCIPnodeGetNumber(node));
}

/*
 * Callback methods of event handler
 */

/** event handler for branching event */
static
SCIP_DECL_EVENTEXEC(eventExecNodeBranched)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_SHADOWTREE* shadowtree;
   SCIP_SHADOWNODE* eventshadownode;
   SCIP_SHADOWNODE* childshadownode;
   SCIP_NODE* eventnode;
   SCIP_NODE** children;
   SCIP_NODE* childnode;
   SCIP_DOMCHG* domchg;
   SCIP_BOUNDCHG* boundchg;
   SCIP_SHADOWBOUNDUPDATE* branchingdecisions;
   SCIP_SHADOWBOUNDUPDATE* update;
   int maxnbranchingdecisions;
   int nbranchingdecisions;
   int nboundchgs;
   int nchildren;
   int i;
   int c;

   assert( scip != NULL );
   assert( eventhdlr != NULL );
   assert( event != NULL );
   assert( SCIPeventGetType(event) & SCIP_EVENTTYPE_NODEBRANCHED );

   /* no branching during probing */
   assert( !SCIPinProbing(scip) );

   eventnode = SCIPeventGetNode(event);
   assert( SCIPgetFocusNode(scip) == eventnode );
   assert( SCIPnodeGetType(eventnode) == SCIP_NODETYPE_FOCUSNODE );

   eventhdlrdata = (SCIP_EVENTHDLRDATA*) SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != NULL );
   assert( scip == eventhdlrdata->scip );

   shadowtree = eventhdlrdata->shadowtree;
   assert( shadowtree != NULL );

   eventshadownode = SCIPshadowTreeGetShadowNode(shadowtree, eventnode);

   /* only add children to the shadowtree if eventnode is in the shadowtree */
   if ( eventshadownode == NULL )
      return SCIP_OKAY;

   assert( eventshadownode->nchildren == 0 );
   assert( eventshadownode->children == NULL );

   SCIP_CALL( SCIPgetChildren(scip, &children, &nchildren) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &eventshadownode->children, nchildren) );
   eventshadownode->nchildren = nchildren;

   maxnbranchingdecisions = 1;  /* good guess that there's one branching variable, because that's likely the number */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchingdecisions, maxnbranchingdecisions) );

   /* get all variables branched upon and check all branches */
   for (c = 0; c < nchildren; ++c)
   {
      nbranchingdecisions = 0;

      childnode = children[c];
      domchg = SCIPnodeGetDomchg(childnode);

      /* loop through all bound changes */
      nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
      for (i = 0; i < nboundchgs; ++i)
      {
         /* get bound change info */
         boundchg = SCIPdomchgGetBoundchg(domchg, i);
         assert( boundchg != NULL );

         /* branching decisions have to be in the beginning of the bound change array */
         if ( SCIPboundchgGetBoundchgtype(boundchg) != SCIP_BOUNDCHGTYPE_BRANCHING )
            break;

         if ( nbranchingdecisions >= maxnbranchingdecisions )
         {
            assert( nbranchingdecisions == maxnbranchingdecisions );
            assert( maxnbranchingdecisions > 0 );
            maxnbranchingdecisions = SCIPcalcMemGrowSize(scip, maxnbranchingdecisions + 1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &branchingdecisions, maxnbranchingdecisions) );
         }
         assert( nbranchingdecisions < maxnbranchingdecisions );

         /* get corresponding branching step */
         update = &branchingdecisions[nbranchingdecisions++];
         update->var = SCIPboundchgGetVar(boundchg);
         update->boundchgtype = SCIPboundchgGetBoundtype(boundchg);
         update->newbound = SCIPboundchgGetNewbound(boundchg);
      }

      /* create the child in the shadow tree */
      SCIP_CALL( SCIPallocBlockMemory(scip, &childshadownode) );
      eventshadownode->children[c] = childshadownode;

      childshadownode->nodeid = SCIPnodeGetNumber(childnode);
      childshadownode->parent = eventshadownode;

      /* children are only set after this node is focused and branched on */
      childshadownode->children = NULL;
      childshadownode->nchildren = 0;

      if ( nbranchingdecisions <= 0 )
         childshadownode->branchingdecisions = NULL;
      else
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &childshadownode->branchingdecisions, nbranchingdecisions) );
         for (i = 0; i < nbranchingdecisions; ++i)
         {
            /* this copies the whole struct */
            childshadownode->branchingdecisions[i] = branchingdecisions[i];
         }
      }
      childshadownode->nbranchingdecisions = nbranchingdecisions;

      /* propagations are only set after this node is focused and branched on */
      childshadownode->propagations = NULL;
      childshadownode->npropagations = 0;

      /* add childshadownode to the nodemap as well
       *
       * The hashtable only checks by the 'nodeid' field, so we just check if there's none with this nodeid.
       */
      assert( !SCIPhashtableExists(shadowtree->nodemap, (void*) childshadownode));
      SCIP_CALL( SCIPhashtableInsert(shadowtree->nodemap, childshadownode) );
   }
   SCIPfreeBufferArray(scip, &branchingdecisions);

   /* also store the propagations in the eventnode (the node that got solved by branching) */
   domchg = SCIPnodeGetDomchg(eventnode);

   /* loop through all bound changes in the focus node */
   nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
   if ( nboundchgs <= 0 )
   {
      assert( nboundchgs == 0 );

      /* this is set to NULL at initialization of this shadownode, already */
      assert( eventshadownode->npropagations == 0 );
      assert( eventshadownode->branchingdecisions == NULL );
   }
   else
   {
      /* just include everything, even the branching decisions! */
      eventshadownode->npropagations = nboundchgs;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &eventshadownode->propagations, nboundchgs) );
      for (i = 0; i < nboundchgs; ++i)
      {
         boundchg = SCIPdomchgGetBoundchg(domchg, i);
         assert( boundchg != NULL );
         update = &(eventshadownode->propagations[i]);
         update->var = SCIPboundchgGetVar(boundchg);
         update->boundchgtype = SCIPboundchgGetBoundtype(boundchg);
         update->newbound = SCIPboundchgGetNewbound(boundchg);
      }
   }

   return SCIP_OKAY;
} /*lint !e715*/


/** event handler for node deletion event */
static
SCIP_DECL_EVENTEXEC(eventExecNodeDeleted)
{ /*lint !e396*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_SHADOWTREE* shadowtree;
   SCIP_NODE* deletednode;
   SCIP_SHADOWNODE* deletedshadownode;
   int c;
   SCIP_SHADOWNODE* childshadownode;

   assert( scip != NULL );
   assert( eventhdlr != NULL );
   assert( event != NULL );
   assert( SCIPeventGetType(event) & SCIP_EVENTTYPE_NODEDELETE );

   deletednode = SCIPeventGetNode(event);
   assert( deletednode != NULL );

   /* probing nodes are not stored */
   if( SCIPnodeGetType(deletednode) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   eventhdlrdata = (SCIP_EVENTHDLRDATA*) SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != NULL );
   assert( scip == eventhdlrdata->scip );

   shadowtree = eventhdlrdata->shadowtree;
   assert( shadowtree != NULL );

   deletedshadownode = SCIPshadowTreeGetShadowNode(shadowtree, deletednode);

   /* no need to delete if not included in the shadowtree */
   if ( deletedshadownode == NULL )
      return SCIP_OKAY;
   assert( deletedshadownode->nodeid == SCIPnodeGetNumber(deletednode) );

   /* It is possible that deletedshadownode has a non-deleted sibling.
    * If the branching variable of this sibling differs from deletedshadownode's,
    * then in the variable branching order also the branching variables of deletedshadownode must be included,
    * e.g., see `shadowtreeFillNodeDepthBranchIndices` in symmetry_lexred.c.
    * As such, we may not delete deletedshadownode just yet. However, we can delete its children.
    * So, mark deletedshadownode as 'ready to delete' by freeing its children, and setting nchildren to -1.
    * SCIP always deletes leaf nodes only, so if `deletedshadownode` is removed,
    * its children in the shadowtree (if they exist) in the 'ready to delete' state. */
   assert( deletedshadownode->nchildren >= 0 );
   assert( (deletedshadownode->nchildren == 0) == (deletedshadownode->children == NULL) );
   for (c = 0; c < deletedshadownode->nchildren; ++c)
   {
      childshadownode = deletedshadownode->children[c];

      /* remove from hashtable */
      SCIP_CALL( SCIPhashtableRemove(shadowtree->nodemap, (void*) childshadownode) );

      /* clean childshadownode */
      assert( childshadownode->npropagations >= 0 );
      assert( (childshadownode->npropagations > 0) != (childshadownode->propagations == NULL) );
      SCIPfreeBlockMemoryArrayNull(scip, &childshadownode->propagations, childshadownode->npropagations);

      assert( childshadownode->nbranchingdecisions >= 0 );
      assert( (childshadownode->nbranchingdecisions > 0) != (childshadownode->branchingdecisions == NULL) );
      SCIPfreeBlockMemoryArrayNull(scip, &childshadownode->branchingdecisions, childshadownode->nbranchingdecisions);

      /* childshadownode must be in the 'ready to delete'-state */
      assert( childshadownode->nchildren < 0 );

      SCIPfreeBlockMemory(scip, &childshadownode);
   }

   assert( (deletedshadownode->nchildren > 0) != (deletedshadownode->children == NULL) );
   if ( deletedshadownode->nchildren > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &deletedshadownode->children, deletedshadownode->nchildren);
   }

   /* mark deletedshadownode as 'ready to delete' */
   deletedshadownode->children = NULL;
   deletedshadownode->nchildren = -1;

   return SCIP_OKAY;
} /*lint !e715*/


/** execution method for all events handled by this eventhandler */
static
SCIP_DECL_EVENTEXEC(eventExec)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );

   eventhdlrdata = (SCIP_EVENTHDLRDATA*) SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != NULL );
   assert( scip == eventhdlrdata->scip );
   assert( eventhdlrdata->clock != NULL );

   SCIP_CALL( SCIPstartClock(scip, eventhdlrdata->clock) );

   switch (SCIPeventGetType(event))
   {
   case SCIP_EVENTTYPE_NODEBRANCHED:
      SCIP_CALL( eventExecNodeBranched(scip, eventhdlr, event, eventdata) );
      break;
   case SCIP_EVENTTYPE_NODEDELETE:
      SCIP_CALL( eventExecNodeDeleted(scip, eventhdlr, event, eventdata) );
      break;
   default:
      SCIPerrorMessage("unrecognized eventtype in shadowtree event handler\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPstopClock(scip, eventhdlrdata->clock) );

   return SCIP_OKAY;
}


/** frees shadow tree data structure */
static
SCIP_RETCODE freeShadowTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SHADOWTREE*      shadowtree          /**< pointer to shadow tree*/
)
{
   int i;
   int nentries;
   SCIP_SHADOWNODE* shadownode;

   assert( scip != NULL );
   assert( shadowtree != NULL );
   assert( shadowtree->nodemap != NULL );

   nentries = SCIPhashtableGetNEntries(shadowtree->nodemap);

   /* free all shadow tree nodes */
   for (i = 0; i < nentries; ++i)
   {
      shadownode = (SCIP_SHADOWNODE*) SCIPhashtableGetEntry(shadowtree->nodemap, i);
      if ( shadownode == NULL )
         continue;

      assert( shadownode != NULL );

      assert( shadownode->npropagations >= 0 );
      assert( (shadownode->npropagations > 0) != (shadownode->propagations == NULL) );
      SCIPfreeBlockMemoryArrayNull(scip, &shadownode->propagations, shadownode->npropagations);

      assert( shadownode->nbranchingdecisions >= 0 );
      assert( (shadownode->nbranchingdecisions > 0) != (shadownode->branchingdecisions == NULL) );
      SCIPfreeBlockMemoryArrayNull(scip, &shadownode->branchingdecisions, shadownode->nbranchingdecisions);

      assert( shadownode->nchildren >= -1 );
      assert( (shadownode->nchildren > 0) != (shadownode->children == NULL) );
      SCIPfreeBlockMemoryArrayNull(scip, &shadownode->children, shadownode->nchildren);

      SCIPfreeBlockMemory(scip, &shadownode);
   }
   SCIPhashtableFree(&(shadowtree->nodemap));

   return SCIP_OKAY;
}


/** destructor of event handler to free shadow tree data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeShadowTree)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert( scip != NULL );
   assert( eventhdlr != NULL );
   assert( SCIPgetStage(scip) != SCIP_STAGE_SOLVING );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != NULL );
   assert( eventhdlrdata->scip == scip );
   assert( eventhdlrdata->clock != NULL );

   SCIP_CALL( SCIPfreeClock(scip, &eventhdlrdata->clock) );

   if ( eventhdlrdata->shadowtree != NULL )
   {
      SCIP_CALL( freeShadowTree(scip, eventhdlrdata->shadowtree) );
      SCIPfreeBlockMemory(scip, &eventhdlrdata->shadowtree);
   }

   SCIPfreeBlockMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}


/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolShadowTree)
{
   int initialnodemapsize;

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_SHADOWTREE* shadowtree;
   SCIP_SHADOWNODE* rootnode;

   assert( scip != NULL );
   assert( eventhdlr != NULL );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != NULL );
   assert( eventhdlrdata->scip == scip );

   assert( eventhdlrdata->shadowtree == NULL );
   assert( SCIPisTransformed(scip) );

   /* early termination */
   if ( !eventhdlrdata->active )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventhdlrdata->shadowtree) );
   shadowtree = eventhdlrdata->shadowtree;

   /* prevent unnecessary reallocations by having a good initial guess for the tree size
    *
    * By default, we initialize NODEMAP_MAX_INITIAL_SIZE slots, unless reasonably fewer nodes suffice.
    * Knowing that a full enumeration tree on n binary variables has size 2^n, we base our guess on this number,
    * counting with the number of binary and integer variables in the problem.
    */
   initialnodemapsize = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) >= NODEMAP_MAX_INITIAL_SIZE_2LOG ?
      NODEMAP_MAX_INITIAL_SIZE :
      MIN(NODEMAP_MAX_INITIAL_SIZE, 1 << (SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip)));  /*lint !e666 !e701 !e747*/
   SCIP_CALL( SCIPhashtableCreate(&shadowtree->nodemap, scip->mem->probmem, initialnodemapsize,
      hashGetKeyShadowNode, hashKeyEqShadowNode, hashKeyValShadowNode, NULL) );

   /* the root node is the only branch-and-bound tree node not created by branching, so add. */
   SCIP_CALL( SCIPallocBlockMemory(scip, &rootnode) );
   rootnode->nodeid = 1ll;  /*lint !e620*/  /* root node has number 1 */
   rootnode->parent = NULL;
   rootnode->children = NULL;
   rootnode->nchildren = 0;
   rootnode->branchingdecisions = NULL;
   rootnode->nbranchingdecisions = 0;
   rootnode->propagations = NULL;
   rootnode->npropagations = 0;

   /* add to the nodemap structure */
   SCIP_CALL( SCIPhashtableInsert(shadowtree->nodemap, rootnode) );

   /* catch NODEBRANCHED and NODEDELETE events */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEBRANCHED | SCIP_EVENTTYPE_NODEDELETE, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolShadowTree)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert( scip != NULL );
   assert( eventhdlr != NULL );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != NULL );
   assert( eventhdlrdata->scip == scip );
   assert( SCIPisTransformed(scip) );

   /* early termination */
   if ( !eventhdlrdata->active )
   {
      assert( eventhdlrdata->shadowtree == NULL );
      return SCIP_OKAY;
   }

   assert( eventhdlrdata->shadowtree != NULL );

   SCIP_CALL( freeShadowTree(scip, eventhdlrdata->shadowtree) );
   SCIPfreeBlockMemory(scip, &eventhdlrdata->shadowtree);
   eventhdlrdata->shadowtree = NULL;

   /* do not listen for NODEBRANCHED events */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEBRANCHED | SCIP_EVENTTYPE_NODEDELETE, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}


/** gets the shadow tree */
SCIP_SHADOWTREE* SCIPgetShadowTree(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert( eventhdlr != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != NULL );

   return eventhdlrdata->shadowtree;
}


/** activates shadow tree eventhandler if it is not already activated (which keeps a copy of the tree) */
SCIP_RETCODE SCIPactivateShadowTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert( eventhdlr != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != NULL );
   assert( eventhdlrdata->scip == scip );
   assert( eventhdlrdata->shadowtree == NULL );

   /* active param may not be changed between (and including) the initsol and exitsol stages */
   SCIP_CALL( SCIPcheckStage(scip, "SCIPactivateShadowTree", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE) );

   eventhdlrdata->active = TRUE;

   return SCIP_OKAY;
}


/** creates event handler for event */
SCIP_RETCODE SCIPincludeEventHdlrShadowTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR**      eventhdlrptr        /**< pointer to store the event handler */
)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &eventhdlrdata) );

#ifndef NDEBUG
   /* only needed for assertions, to check whether we're working with the correct SCIP. */
   eventhdlrdata->scip = scip;
#endif

   /* shadow tree must be activated */
   eventhdlrdata->active = FALSE;

   /* do not start with a shadow tree by default. Initialize at initsol, remove at exitsol. */
   eventhdlrdata->shadowtree = NULL;
   eventhdlr = NULL;

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExec, eventhdlrdata) );
   assert(eventhdlr != NULL);
   *eventhdlrptr = eventhdlr;

   /* clock */
   SCIP_CALL( SCIPcreateClock(scip, &eventhdlrdata->clock) );

   /* set non fundamental callbacks via setter functions */

   /* frees the event handler */
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeShadowTree) );

   /* initialize the shadowtree data structure, initialize by setting the root node */
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolShadowTree) );

   /* free the shadowtree data structure */
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolShadowTree) );

   return SCIP_OKAY;
}
