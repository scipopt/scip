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
   SCIP_NODETYPE_JUNCTION = 4,          /**< fork without LP solution */
   SCIP_NODETYPE_FORK     = 5,          /**< fork with solved LP and added rows and columns */
   SCIP_NODETYPE_SUBROOT  = 6           /**< fork with solved LP and arbitrarily changed rows and columns */
};
typedef enum Nodetype NODETYPE;         /**< type of node */

typedef struct Actnode ACTNODE;         /**< data for the actual node */
typedef struct Child CHILD;             /**< data for child nodes */
typedef struct Sibling SIBLING;         /**< data for sibling nodes */
typedef struct Junction JUNCTION;       /**< data for junction nodes */
typedef struct Fork FORK;               /**< data for fork nodes */
typedef struct Subroot SUBROOT;         /**< data for subroot nodes */
typedef struct Node NODE;               /**< node data structure */
typedef struct Tree TREE;               /**< branch and bound tree */



#include "lp.h"
#include "cons.h"



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
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPnodeActivate(               /**< activates a leaf node */
   NODE*            node,               /**< leaf node to activate */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
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
RETCODE SCIPnodeAddCons(                /**< adds local constraint to the node */
   NODE*            node,               /**< node to add constraint to */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
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
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPtreeFree(                   /**< frees tree data structure */
   TREE**           tree                /**< pointer to tree data structure */
   );

extern
RETCODE SCIPtreeLoadLP(                 /**< constructs the LP and loads LP state for fork/subroot of the active node */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPtreeAddLocalCons(           /**< adds local constraint to the active node */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   );

extern
RETCODE SCIPtreeAddGlobalCons(          /**< adds global constraint to the problem */
   TREE*            tree,               /**< branch-and-bound tree */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   );

extern
NODE** SCIPtreeGetChildren(             /**< gets children array of actual node */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
int SCIPtreeGetNChildren(               /**< gets number of children of actual node */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
NODE** SCIPtreeGetSiblings(             /**< gets siblings array of actual node */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
int SCIPtreeGetNSiblings(               /**< gets number of siblings of actual node */
   TREE*            tree                /**< branch-and-bound tree */
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
