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

/**@file   tree.h
 * @brief  branch-and-bound tree datastructures and operations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TREE_H__
#define __TREE_H__

enum Nodetype
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
typedef enum Nodetype NODETYPE;         /**< type of node */

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
      JUNCTION*     junction;           /**< data for junction nodes */
      FORK*         fork;               /**< data for fork nodes */
      SUBROOT*      subroot;            /**< data for subroot nodes */
   } data;
   NODE*            parent;             /**< parent node in the tree */
   CONSLIST*        conslist;           /**< list of constraints created at this node */
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
   DOMCHGDYN*       actnodedomchg;      /**< domain changes of the active node */
   DOMCHGDYN**      childrendomchg;     /**< domain changes of the child nodes */
   DOMCHGDYN**      siblingsdomchg;     /**< domain changes of the sibling nodes */
   SOL*             actpseudosol;       /**< actual pseudosolution with all variables set to their best bounds */
   int*             pathnlpcols;        /**< array with number of LP columns for each problem in active path */
   int*             pathnlprows;        /**< array with number of LP rows for each problem in active path */
   int              pathlen;            /**< length of the actual path (== depth of the current node + 1) */
   int              pathsize;           /**< number of available slots in path arrays */
   int              correctlpdepth;     /**< depth to which current LP data corresponds to LP data of active path */
   int              childrensize;       /**< available slots in children vector */
   int              nchildren;          /**< actual number of children (number of used slots in children vector) */
   int              siblingssize;       /**< available slots in siblings vector */
   int              nsiblings;          /**< actual number of siblings (number of used slots in siblings vector) */
};



/*
 * Node methods
 */

extern
DECL_SORTPTRCOMP(SCIPnodeCmpLowerbound);/**< node comparator for best lower bound */

extern
RETCODE SCIPnodeCreate(                 /**< creates a child node of the active node */
   NODE**           node,               /**< pointer to node data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
RETCODE SCIPnodeFree(                   /**< frees node */
   NODE**           node,               /**< node data */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIPnodeCaptureLPIState(           /**< increases the reference counter of the LP state in the fork or subroot node */
   NODE*            node,               /**< fork/subroot node */
   int              numuses             /**< number to add to the usage counter */
   );

extern
RETCODE SCIPnodeReleaseLPIState(        /**< decreases the reference counter of the LP state in the fork or subroot node */
   NODE*            node,               /**< fork/subroot node */
   MEMHDR*          memhdr,             /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPnodeActivate(               /**< activates a leaf node */
   NODE*            node,               /**< leaf node to activate (or NULL to deactivate all nodes) */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
RETCODE SCIPnodeToDeadend(              /**< converts the active node into a deadend node */
   NODE*            node,               /**< node to convert */
   MEMHDR*          memhdr,             /**< block memory buffers */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPnodeToJunction(             /**< converts the active node into a junction node */
   NODE*            node,               /**< node to convert */
   MEMHDR*          memhdr,             /**< block memory */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
RETCODE SCIPnodeToFork(                 /**< converts the active node into a fork node */
   NODE*            node,               /**< node to convert */
   MEMHDR*          memhdr,             /**< block memory */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPnodeToSubroot(              /**< converts the active node into a subroot node */
   NODE*            node,               /**< node to convert */
   MEMHDR*          memhdr,             /**< block memory */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPnodeAddCons(                /**< adds local constraint to the node and captures it */
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   );

extern
RETCODE SCIPnodeAddBoundchg(            /**< adds bound change to actual node, child or sibling of actual node */
   NODE*            node,               /**< node to add bound change to */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

extern
NODETYPE SCIPnodeGetType(               /**< gets the type of the node */
   NODE*            node                /**< node */
   );

extern
int SCIPnodeGetDepth(                   /**< gets the depth of the node */
   NODE*            node                /**< node */
   );

extern
Real SCIPnodeGetLowerBound(             /**< gets the lower bound of the node */
   NODE*            node                /**< node */
   );


/*
 * Tree methods
 */

extern
RETCODE SCIPtreeCreate(                 /**< creates an initialized tree data structure */
   TREE**           tree,               /**< pointer to tree data structure */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   PROB*            prob                /**< problem data */
   );

extern
RETCODE SCIPtreeFree(                   /**< frees tree data structure */
   TREE**           tree,               /**< pointer to tree data structure */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPtreeLoadLP(                 /**< constructs the LP and loads LP state for fork/subroot of the active node */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPtreeAddLocalCons(           /**< adds local constraint to the active node and captures it */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   );

extern
RETCODE SCIPtreeAddGlobalCons(          /**< adds global constraint to the problem and captures it */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   );

extern
RETCODE SCIPtreeBoundChanged(           /**< notifies tree, that a bound of a variable changed */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable that changed */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

extern
int SCIPtreeGetNLeaves(                 /**< gets number of leaves */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern   
int SCIPtreeGetNNodes(                  /**< gets number of nodes (children + siblings + leaves) */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
NODE* SCIPtreeGetBestLeaf(              /**< gets the best leaf from the node queue */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
NODE* SCIPtreeGetBestNode(              /**< gets the best node from the tree (child, sibling, or leaf) */
   TREE*            tree,               /**< branch-and-bound tree */
   const SET*       set                 /**< global SCIP settings */
   );

#endif
