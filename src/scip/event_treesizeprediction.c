/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_treesizeprediction
 * @brief  eventhdlr for tree-size prediction related events
 * @author Pierre Le Bodic
 */


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_DEBUG 1

#include "scip/event_treesizeprediction.h"
#include "scip/misc.h"
#include "scip/struct_tree.h"

#include "string.h"
#define EVENTHDLR_NAME         "treesizeprediction"
#define EVENTHDLR_DESC         "event handler for tree-size prediction related events"

#define DEFAULT_HASHMAP_SIZE   100000

/*
 * Data structures
 */

/*
 * Indicates for a given node if/how the size of its subtree is computed.
 * UNKNOWN: the node has children, both with UNKNOWN sizes. No tree-size estimate at this node, it is UNKNOWN.
 * ESTIMATED: the node has children, exactly one of them has UNKNOWN size. The tree-size at this node is ESTIMATED.
 * KNOWN: the node is a leaf or both its children have KNOWN size. The tree-size at this node is thus KNOWN.
 */

typedef enum {UNKNOWN, ESTIMATED, KNOWN} SizeStatus;

typedef struct TreeSizeEstimateTree TSEtree;
struct TreeSizeEstimateTree
{
   TSEtree *parent;
   TSEtree *leftchild;
   TSEtree *rightchild;

   SCIP_Bool prunedinPQ; /* Whether the node has been pruned while in the priority queue, and thus never focused */
   SCIP_Longint number; /* The number (id) of the node, as assigned by SCIP */
   SCIP_Real lowerbound; /* The lower bound at that node. TODO update this if a node gets an update */
   SizeStatus status; /* See enum SizeStatus documentation */
   /*SCIP_Longint treesize;*/ /* The computed tree-size */
};

/**
 * Estimates the tree-size of a tree, using the given upperbound to determine
 * if a node is counted as a leaf (independent of whether it has children).
 */
static
SizeStatus estimateTreeSize(SCIP* scip, TSEtree *node, SCIP_Longint* size, SCIP_Real upperbound)
{
   assert(node != NULL);
   assert(size != NULL);
   /* base cases: determine if the current node is a leaf */
   if( node->prunedinPQ == TRUE )
   {
      assert(node->leftchild == NULL);
      assert(node->rightchild == NULL);
      *size = 1;
      return KNOWN;
   }
   else if( SCIPisGE(scip, node->lowerbound, upperbound) )
   {
      *size = 1;
      return KNOWN;
   } 
   else if(node->leftchild == NULL) /* The node is not a leaf but still needs to be solved (and possibly branched on) */
   {
      assert(node->rightchild == NULL);
      *size = -1; /* TODO remove once code is finished & debugged */
      return UNKNOWN;
   }
   else /* The node has children */
   {
      SCIP_Longint leftsize;
      SCIP_Longint rightsize;
      SizeStatus leftstatus;
      SizeStatus rightstatus;

      assert(node->leftchild != NULL && node->rightchild != NULL);

      leftstatus = estimateTreeSize(scip, node->leftchild, &leftsize, upperbound);
      rightstatus = estimateTreeSize(scip, node->rightchild, &rightsize, upperbound);

      assert(leftsize > 0 || leftstatus == UNKNOWN);
      assert(rightsize > 0 || rightstatus == UNKNOWN);

      if(leftstatus == UNKNOWN && rightstatus == UNKNOWN) /* Neither child has information on tree-size*/
      {
         *size = -1; /* TODO remove once code is finished & debugged */
         return UNKNOWN;
      }
      else if ( leftstatus == KNOWN && rightstatus == KNOWN  ) /* If both left and right subtrees are known */
      {
         *size = 1 + leftsize + rightsize;
         return KNOWN;
      }
      else /* Exactly one subtree is UNKNOWN: we estimate */
      {
         assert((leftstatus == UNKNOWN) ^ (rightstatus == UNKNOWN));
         if( leftstatus == UNKNOWN )
            *size = 1 + 2*rightsize;
         else
            *size = 1 + 2*leftsize;
         return ESTIMATED;
      }
   }
}

/**
 * Recursively frees memory in the tree.
 */
static
void freeTreeMemory(SCIP *scip, TSEtree **tree)
{
   assert(tree != NULL);
   assert(*tree != NULL);
   /* postfix traversal */
   if( (*tree)->leftchild != NULL )
      freeTreeMemory(scip, &((*tree)->leftchild));
   if( (*tree)->rightchild != NULL )
      freeTreeMemory(scip, &((*tree)->rightchild));
   SCIPfreeMemory(scip, tree);
   *tree = NULL;
}

/** event handler data */
struct SCIP_EventhdlrData
{
   /* Parameters */
   int hashmapsize;

   /* Internal variables */
//   SCIP_Bool initialized;
   unsigned int nodesfound;
   TSEtree *tree; /* The representation of the B&B tree */
   SCIP_HASHMAP *opennodes; /* The open nodes (that have yet to be branched on). The key is the (scip) id/number of the SCIP_Node */
};

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolTreeSizePrediction)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   /* SCIPdebugMsg(scip, "initsol method of eventhdlr "EVENTHDLR_NAME"\n"); */
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

//   eventhdlrdata->initialized = TRUE;
   eventhdlrdata->nodesfound = 0;
   eventhdlrdata->tree = NULL;
   SCIP_CALL( SCIPhashmapCreate(&(eventhdlrdata->opennodes), SCIPblkmem(scip), eventhdlrdata->hashmapsize) );

   /* We catch node solved events */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, NULL) );

   /* We catch priority queue nodes being removed by bound */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_PQNODEINFEASIBLE, eventhdlr, NULL, NULL) );

   /* We catch updates to the primal bound */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolTreeSizePrediction)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   /* SCIPdebugMsg(scip, "exitsol method of eventhdlr "EVENTHDLR_NAME"\n"); */
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPdebugMessage("Found %u nodes in the B&B tree\n", eventhdlrdata->nodesfound);

   #ifdef SCIP_DEBUG
   {
      SCIP_Longint estimate;
      SizeStatus status = estimateTreeSize(scip, eventhdlrdata->tree, &estimate, SCIPgetUpperbound(scip));
      assert(status == KNOWN);/* TODO remove after debug */
      assert(estimate > 0);
      SCIPdebugMessage("Final estimate: %" SCIP_LONGINT_FORMAT " nodes in the B&B tree\n", estimate);
   }
   #endif

//   if( eventhdlrdata->initialized == TRUE )
//   {
      assert(SCIPhashmapIsEmpty(eventhdlrdata->opennodes) == TRUE);
      SCIPhashmapFree(&(eventhdlrdata->opennodes));
      if( eventhdlrdata->nodesfound > 0 )
      {
         assert(eventhdlrdata->tree != NULL);
         freeTreeMemory(scip, &(eventhdlrdata->tree));
      }
//   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecTreeSizePrediction)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_NODE* scipnode; /* The node found by the event (if any) */

   /* SCIPdebugMsg(scip, "exec method of eventhdlr "EVENTHDLR_NAME"\n"); */
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->opennodes != NULL);

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND )
   {
      /* When a new primal bound is found, then some of the leaves that were
       * previously infeasible could have an ancester that would have been
       * pruned by this new primal bound. We are going to trim our
       * representation of the tree so that previous branched nodes may
       * become leaves in our data structure.
       */
      //newlowerbound = SCIPgetLowerbound(scip);
      /* SCIPdebugMsg(scip, "New best solution found\n"); */
      scipnode = NULL;
   }
   else
   {
      SCIP_NODE *scipparent; /* The parent of the scipnode we found */
      scipnode = SCIPeventGetNode(event);
      assert(scipnode != NULL);
      scipparent = SCIPnodeGetParent(scipnode);
 
//      if( SCIPeventGetType(event) == SCIP_EVENTTYPE_PQNODEINFEASIBLE )
//      {
//         SCIPdebugMessage("Node event for PQ node #%" SCIP_LONGINT_FORMAT "\n", SCIPnodeGetNumber(scipnode) );
//      }
//      else
      {
         /* We create a new node in our tree structure */
         TSEtree *newnode;
         TSEtree *parentnode;
         eventhdlrdata->nodesfound += 1;

         SCIP_CALL( SCIPallocMemory(scip, &newnode) );
         newnode->lowerbound = SCIPnodeGetLowerbound(scipnode);
         newnode->leftchild = NULL;
         newnode->rightchild = NULL;
         newnode->number = SCIPnodeGetNumber(scipnode);
         /* SCIPdebugMessage("Node event for focus node #%" SCIP_LONGINT_FORMAT "\n", newnode->number); */
         /* TODO handle parent, which should be in the open nodes  */
         if( scipparent != NULL )
         {
            SCIP_Longint parentnumber = SCIPnodeGetNumber(scipparent);
            parentnode = (TSEtree*) SCIPhashmapGetImage(eventhdlrdata->opennodes, (void *)parentnumber);
            assert(parentnode != NULL);
            assert(parentnode->number == parentnumber);

            newnode->parent = parentnode;

            if(parentnode->leftchild == NULL)
               parentnode->leftchild = newnode;
            else
            {
               assert(parentnode->rightchild == NULL);
               parentnode->rightchild = newnode;
               /* We have seen all the children of this parent, thus we can remove it from open nodes */
               SCIP_CALL( SCIPhashmapRemove(eventhdlrdata->opennodes, (void *)parentnumber) );
            }
         }
         else
         {
            /* Then this should be the root node */
            parentnode = NULL;
            newnode->parent = parentnode;
            assert(eventhdlrdata->tree == NULL);
            eventhdlrdata->tree = newnode;
         }
    
         newnode->prunedinPQ = FALSE;
         switch(SCIPeventGetType(event)) {
            case SCIP_EVENTTYPE_PQNODEINFEASIBLE:
               newnode->prunedinPQ = TRUE;
            case SCIP_EVENTTYPE_NODEFEASIBLE:
            case SCIP_EVENTTYPE_NODEINFEASIBLE:
               /* When an (in)feasible node is found, this corresponds to a new sample
                * (in Knuth's algorithm). This may change the tree-size estimate. */
               newnode->status = KNOWN;
               break;
            case SCIP_EVENTTYPE_NODEBRANCHED:
               /* When a node is branched on, we need to add the corresponding nodes
                * to our own data structure */
               newnode->status = UNKNOWN;
               SCIP_CALL( SCIPhashmapInsert(eventhdlrdata->opennodes, (void *) (newnode->number), newnode) ); //TODO since the number is a unique indentifier, perhaps we could use an array?
               break;
            default:
               SCIPerrorMessage("Missing case in this switch.\n");
               SCIPABORT();
         }
      }
   }
   return SCIP_OKAY;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeTreeSizePrediction)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   /* SCIPdebugMsg(scip, "eventfree method of eventhdlr "EVENTHDLR_NAME"\n"); */
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);
   return SCIP_OKAY;
}

/** creates event handler for tree-size prediction event */
SCIP_RETCODE SCIPincludeEventHdlrTreeSizePrediction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create tree-size prediction event handler data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   /* TODO: (optional) create event handler specific data here */

   eventhdlr = NULL;

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecTreeSizePrediction, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, NULL) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeTreeSizePrediction) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, NULL) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, NULL) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolTreeSizePrediction) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolTreeSizePrediction) );
   SCIP_CALL( SCIPsetEventhdlrDelete(scip, eventhdlr, NULL) );

   /* add tree-size prediction event handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "estimates/hashmapsize", "Default hashmap size to store the open nodes of the B&B tree", &(eventhdlrdata->hashmapsize), TRUE, DEFAULT_HASHMAP_SIZE, 0, INT_MAX, NULL, NULL) );


   return SCIP_OKAY;
}


