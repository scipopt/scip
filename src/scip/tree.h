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

/**@file   tree.h
 * @brief  branch-and-bound tree datastructures and operations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TREE_H__
#define __TREE_H__

enum NodeType
{
   SCIP_NODETYPE_ACTNODE  = 0,          /**< the active node, whose data is stored in the dynamic data object */
   SCIP_NODETYPE_SIBLING  = 1,          /**< unsolved sibling of the active node */
   SCIP_NODETYPE_CHILD    = 2,          /**< unsolved child of the active node */
   SCIP_NODETYPE_LEAF     = 3,          /**< unsolved leaf of the tree, stored in the tree's queue */
   SCIP_NODETYPE_DEADEND  = 4,          /**< temporary type of active node, if it was solved completely */
   SCIP_NODETYPE_JUNCTION = 5,          /**< fork without LP solution */
   SCIP_NODETYPE_FORK     = 6,          /**< fork with solved LP and added rows and columns */
   SCIP_NODETYPE_SUBROOT  = 7           /**< fork with solved LP and arbitrarily changed rows and columns */
};
typedef enum NodeType NODETYPE;         /**< type of node */

typedef struct Actnode ACTNODE;         /**< data for the actual node */
typedef struct Child CHILD;             /**< data for child nodes */
typedef struct Sibling SIBLING;         /**< data for sibling nodes */
typedef struct Leaf LEAF;               /**< data for leaf nodes */
typedef struct Junction JUNCTION;       /**< data for junction nodes */
typedef struct Fork FORK;               /**< data for fork nodes */
typedef struct Subroot SUBROOT;         /**< data for subroot nodes */
typedef struct Node NODE;               /**< node data structure */
typedef struct Tree TREE;               /**< branch and bound tree */



#include "lp.h"
#include "lpi.h"
#include "cons.h"
#include "nodesel.h"
#include "prob.h"
#include "sol.h"
#include "branch.h"
#include "event.h"



/** child information (should not exceed the size of a pointer) */
struct Child
{
   int              arraypos;           /**< position of node in the children array */
};

/** sibling information (should not exceed the size of a pointer) */
struct Sibling
{
   int              arraypos;           /**< position of node in the siblings array */
};

/** leaf information (should not exceed the size of a pointer) */
struct Leaf
{
   NODE*            lpfork;             /**< fork/subroot node defining the LP state of the leaf */
};

/** fork without LP solution, where only bounds and constraints have been changed */
struct Junction
{
   int              nchildren;          /**< number of children of this parent node */
};

/** fork with solved LP, where bounds and constraints have been changed, and rows and columns were added */
struct Fork
{
   COL**            addedcols;          /**< array with pointers to new columns added at this node into the LP */
   ROW**            addedrows;          /**< array with pointers to new rows added at this node into the LP */
   LPISTATE*        lpistate;           /**< LP state information */
   int              naddedcols;         /**< number of columns added at this node */
   int              naddedrows;         /**< number of rows added at this node */
   int              nchildren;          /**< number of children of this parent node */
   int              nlpistateref;       /**< number of times, the LP state is needed */   
};

/** fork with solved LP, where bounds and constraints have been changed, and rows and columns were removed and added */
struct Subroot
{
   COL**            cols;               /**< array with pointers to the columns in the same order as in the LP */
   ROW**            rows;               /**< array with pointers to the rows in the same order as in the LP */
   LPISTATE*        lpistate;           /**< LP state information */
   int              ncols;              /**< number of columns in the LP */
   int              nrows;              /**< number of rows in the LP */
   int              nchildren;          /**< number of children of this parent node */
   int              nlpistateref;       /**< number of times, the LP state is needed */   
};

/** node data structure */
struct Node
{
   union
   {
      SIBLING       sibling;            /**< data for sibling nodes */
      CHILD         child;              /**< data for child nodes */
      LEAF          leaf;               /**< data for leaf nodes */
      JUNCTION      junction;           /**< data for junction nodes */
      FORK*         fork;               /**< data for fork nodes */
      SUBROOT*      subroot;            /**< data for subroot nodes */
   } data;
   NODE*            parent;             /**< parent node in the tree */
   CONSSETCHG*      conssetchg;         /**< constraint set changes at this node or NULL */
   DOMCHG*          domchg;             /**< domain changes at this node or NULL */
   Real             lowerbound;         /**< lower (dual) LP bound of subtree */
   unsigned int     depth:16;           /**< depth in the tree */
   unsigned int     nodetype:3;         /**< type of node */
   unsigned int     active:1;           /**< is node in the path to the actual active node? */
};

/** branch and bound tree */
struct Tree
{
   NODE*            root;               /**< root node of the tree */
   NODEPQ*          leaves;             /**< leaves of the tree */
   NODE**           path;               /**< array of fork/subtree nodes storing the active path from root to leaf */
   NODE*            actnode;            /**< active node */
   NODE*            actlpfork;          /**< fork/subroot node defining the LP state of the active node */
   NODE*            actsubroot;         /**< root of the active subtree */
   NODE**           children;           /**< array with children of the active node */
   NODE**           siblings;           /**< array with siblings of the active node */
   CONSSETCHGDYN*   actnodeconssetchg;  /**< constraint set changed of the active node */
   CONSSETCHGDYN**  childrenconssetchg; /**< constraint set changed of the child nodes */
   CONSSETCHGDYN**  siblingsconssetchg; /**< constraint set changed of the sibling nodes */
   DOMCHGDYN*       actnodedomchg;      /**< domain changes of the active node */
   DOMCHGDYN**      childrendomchg;     /**< domain changes of the child nodes */
   DOMCHGDYN**      siblingsdomchg;     /**< domain changes of the sibling nodes */
   int*             pathnlpcols;        /**< array with number of LP columns for each problem in active path */
   int*             pathnlprows;        /**< array with number of LP rows for each problem in active path */
   Real             actpseudoobjval;    /**< actual pseudo solution value with all variables set to their best bounds,
                                         *   ignoring variables, with infinite best bound */
   int              actpseudoobjvalinf; /**< number of variables with infinite best bound in actual pseudo solution */
   int              pathlen;            /**< length of the actual path (== depth of the current node + 1) */
   int              pathsize;           /**< number of available slots in path arrays */
   int              correctlpdepth;     /**< depth to which current LP data corresponds to LP data of active path */
   int              childrensize;       /**< available slots in children vector */
   int              nchildren;          /**< actual number of children (number of used slots in children vector) */
   int              siblingssize;       /**< available slots in siblings vector */
   int              nsiblings;          /**< actual number of siblings (number of used slots in siblings vector) */
   unsigned int     actnodehaslp:1;     /**< is LP being processed in the active node? */
   unsigned int     cutoffdelayed:1;    /**< the treeCutoff() call was delayed because of diving and has to be executed */
};



/*
 * Node methods
 */

/** node comparator for best lower bound */
extern
DECL_SORTPTRCOMP(SCIPnodeCmpLowerbound);

/** creates a child node of the active node */
extern
RETCODE SCIPnodeCreate(
   NODE**           node,               /**< pointer to node data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree                /**< branch-and-bound tree */
   );

/** frees node */
extern
RETCODE SCIPnodeFree(
   NODE**           node,               /**< node data */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

/** increases the reference counter of the LP state in the fork or subroot node */
extern
void SCIPnodeCaptureLPIState(
   NODE*            node,               /**< fork/subroot node */
   int              nuses               /**< number to add to the usage counter */
   );

/** decreases the reference counter of the LP state in the fork or subroot node */
extern
RETCODE SCIPnodeReleaseLPIState(
   NODE*            node,               /**< fork/subroot node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   );

/** activates a leaf node */
extern
RETCODE SCIPnodeActivate(
   NODE*            node,               /**< leaf node to activate (or NULL to deactivate all nodes) */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** adds local constraint to the node and captures it */
extern
RETCODE SCIPnodeAddCons(
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   CONS*            cons                /**< constraint to add */
   );

/** disables constraint's separation, enforcing, and propagation capabilities at the node, and captures constraint */
extern
RETCODE SCIPnodeDisableCons(
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   CONS*            cons                /**< constraint to disable */
   );

/** adds bound change to actual node, child or sibling of actual node */
extern
RETCODE SCIPnodeAddBoundchg(
   NODE*            node,               /**< node to add bound change to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** if given value is larger than the node's lower bound, sets the node's lower bound to the new value */
extern
void SCIPnodeUpdateLowerbound(
   NODE*            node,               /**< node to update lower bound for */
   Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets the type of the node */
extern
NODETYPE SCIPnodeGetType(
   NODE*            node                /**< node */
   );

/** gets the depth of the node */
extern
int SCIPnodeGetDepth(
   NODE*            node                /**< node */
   );

/** gets the lower bound of the node */
extern
Real SCIPnodeGetLowerbound(
   NODE*            node                /**< node */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPnodeGetType(node)           ( (node)->nodetype )
#define SCIPnodeGetDepth(node)          ( (node)->depth )
#define SCIPnodeGetLowerbound(node)     ( (node)->lowerbound )


#endif



/*
 * Tree methods
 */

/** creates an initialized tree data structure */
extern
RETCODE SCIPtreeCreate(
   TREE**           tree,               /**< pointer to tree data structure */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   );

/** frees tree data structure */
extern
RETCODE SCIPtreeFree(
   TREE**           tree,               /**< pointer to tree data structure */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

/** resorts the leave priority queue (necessary for changes in node selector) */
extern
RETCODE SCIPtreeResortLeaves(
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** cuts off nodes with lower bound not better than given upper bound */
extern
RETCODE SCIPtreeCutoff(
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             upperbound          /**< upper bound: all nodes with lowerbound >= upperbound are cut off */
   );

/** constructs the LP and loads LP state for fork/subroot of the active node */
extern
RETCODE SCIPtreeLoadLP(
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

/** branches on a variable; if solution value x' is fractional, two child nodes are created
 *  (x <= floor(x'), x >= ceil(x')), if solution value is integral, three child nodes are created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
extern
RETCODE SCIPtreeBranchVar(
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var                 /**< variable to branch on */
   );

/** updates actual pseudo objective value for a change in a variable's objective value or bounds */
extern
RETCODE SCIPtreeUpdateVar(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldobj,             /**< old objective value of variable */
   Real             oldlb,              /**< old objective value of variable */
   Real             oldub,              /**< old objective value of variable */
   Real             newobj,             /**< new objective value of variable */
   Real             newlb,              /**< new objective value of variable */
   Real             newub               /**< new objective value of variable */
   );

/** updates actual pseudo objective value for a change in a variable's objective value */
extern
RETCODE SCIPtreeUpdateVarObj(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldobj,             /**< old objective value of variable */
   Real             newobj              /**< new objective value of variable */
   );

/** updates actual pseudo objective value for a change in a variable's lower bound */
extern
RETCODE SCIPtreeUpdateVarLb(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldlb,              /**< old lower bound of variable */
   Real             newlb               /**< new lower bound of variable */
   );

/** updates actual pseudo objective value for a change in a variable's upper bound */
extern
RETCODE SCIPtreeUpdateVarUb(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   Real             oldub,              /**< old upper bound of variable */
   Real             newub               /**< new upper bound of variable */
   );

/** gets number of leaves */
extern
int SCIPtreeGetNLeaves(
   TREE*            tree                /**< branch-and-bound tree */
   );

/** gets number of nodes (children + siblings + leaves) */
extern   
int SCIPtreeGetNNodes(
   TREE*            tree                /**< branch-and-bound tree */
   );

/** gets the best child of the active node */
extern
NODE* SCIPtreeGetBestChild(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   );

/** gets the best sibling of the active node */
extern
NODE* SCIPtreeGetBestSibling(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   );

/** gets the best leaf from the node queue */
extern
NODE* SCIPtreeGetBestLeaf(
   TREE*            tree                /**< branch-and-bound tree */
   );

/** gets the best node from the tree (child, sibling, or leaf) */
extern
NODE* SCIPtreeGetBestNode(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   );

/** gets the minimal lower bound of all nodes in the tree */
extern
Real SCIPtreeGetLowerbound(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   );

/** gets the node with minimal lower bound of all nodes in the tree (child, sibling, or leaf) */
extern
NODE* SCIPtreeGetLowerboundNode(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   );

/** gets the lower bound of the active node */
extern
Real SCIPtreeGetActLowerbound(
   TREE*            tree                /**< branch-and-bound tree */
   );

/** gets the average lower bound of all nodes in the tree */
extern
Real SCIPtreeGetAvgLowerbound(
   TREE*            tree,               /**< branch-and-bound tree */
   Real             upperbound          /**< global upper bound */
   );

/** gets the pseudo objective value of the active node */
extern
Real SCIPtreeGetActPseudoobjval(
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   );

#endif
