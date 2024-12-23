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

#ifndef __SCIP_EVENT_SHADOWTREE_H__
#define __SCIP_EVENT_SHADOWTREE_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_event.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** bound change for branch-and-bound tree node in shadow tree */
struct SCIP_ShadowBoundUpdate
{
   SCIP_VAR*             var;                /**< changed variable */
   SCIP_Real             newbound;           /**< bound change */
   SCIP_BOUNDTYPE        boundchgtype;       /**< which bound of variable is changed (upper or lower) */
};
typedef struct SCIP_ShadowBoundUpdate SCIP_SHADOWBOUNDUPDATE;

/** branch-and-bound tree node for the shadowtree */
struct SCIP_ShadowNode
{
   SCIP_Longint          nodeid;             /**< ID of corresponding branch-and-bound tree node */
   struct SCIP_ShadowNode* parent;           /**< parent of this shadowtree node. NULL iff it is the root node */
   struct SCIP_ShadowNode** children;        /**< list of children of this shadowtree node. NULL iff it is a leaf */
   int                   nchildren;          /**< number of elements in children
                                              *   0 iff it is a leaf, -1 iff original node is DELETED */
   SCIP_SHADOWBOUNDUPDATE* branchingdecisions;/**< the variables branched on in this node.
                                               *  NULL iff nbranchingdecisions == 0 */
   int                   nbranchingdecisions;/**< the number of variables branched on in this node
                                              *   0 iff branchingdecisions == NULL */
   SCIP_SHADOWBOUNDUPDATE* propagations;     /**< the propagation (and branching decisions) updates in the node
                                              *   This is populated after branching with the propagations in that node.
                                              *   NULL iff npropagations == 0 */
   int                   npropagations;      /**< the number of propagations. 0 iff propagations == NULL */
};
typedef struct SCIP_ShadowNode SCIP_SHADOWNODE;

/** shadow tree data structure */
struct SCIP_ShadowTree
{
   SCIP_HASHTABLE*       nodemap;            /**< pointer to the hashmap containing all shadow tree nodes */
};
typedef struct SCIP_ShadowTree SCIP_SHADOWTREE;

/** get the time spent in the shadow tree eventhdlr */
SCIP_EXPORT
SCIP_Real SCIPgetShadowTreeEventHandlerExecutionTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
);

/** given a node number, returns the node in the shadow tree, or NULL if it doesn't exist */
SCIP_EXPORT
SCIP_SHADOWNODE* SCIPshadowTreeGetShadowNodeFromNodeNumber(
   SCIP_SHADOWTREE*      shadowtree,         /**< pointer to the shadow tree */
   SCIP_Longint          nodeno              /**< index of the node, equivalent to the standard branch-and-bound tree */
);

/** given a node, returns the node in the shadowtree, or NULL if it doesn't exist */
SCIP_EXPORT
SCIP_SHADOWNODE* SCIPshadowTreeGetShadowNode(
   SCIP_SHADOWTREE*      shadowtree,         /**< pointer to the shadow tree */
   SCIP_NODE*            node                /**< node from the actual branch-and-bound tree */
);

/** gets the shadow tree */
SCIP_EXPORT
SCIP_SHADOWTREE* SCIPgetShadowTree(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
);

/** activates shadow tree eventhandler if it is not already activated (which keeps a copy of the tree) */
SCIP_EXPORT
SCIP_RETCODE SCIPactivateShadowTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
);

/** creates event handler for event */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeEventHdlrShadowTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR**      eventhdlrptr        /**< pointer to store the event handler */
);

#ifdef __cplusplus
}
#endif

#endif
