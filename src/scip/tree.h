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

#include "lp.h"


enum Nodetype
{
   SCIP_NODETYPE_LEAF    = 0,           /**< unsolved leaf of the tree */
   SCIP_NODETYPE_ACTNODE = 1,           /**< the active node, whose data is stored in the dynamic data object */
   SCIP_NODETYPE_FORK    = 2,           /**< solved fork, where rows were only added to the father */
   SCIP_NODETYPE_SUBROOT = 3            /**< solved fork, where rows were added, deleted, or rearranged */
};

typedef enum Nodetype NODETYPE;         /**< type of node */
typedef struct Actnode ACTNODE;         /**< data for the actual node */
typedef struct Leaf LEAF;               /**< data for leaf nodes */
typedef struct Fork FORK;               /**< data for fork nodes */
typedef struct Subroot SUBROOT;         /**< data for subroot nodes */
typedef struct Node NODE;               /**< node data structure */
typedef struct Tree TREE;               /**< branch and bound tree */



/*
 * Node methods
 */

extern
DECL_SORTPTRCOMP(SCIPnodeCmpLowerbound);/**< node comparator for best lower bound */

extern
NODE* SCIPnodeCreate(                   /**< creates a leaf node */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   NODE*            parent,             /**< parent node in the tree */
   LPSTATE*         lpstate             /**< pointer to LP state information */
   );

extern
void SCIPnodeFree(                      /**< frees node */
   NODE**           node,               /**< node data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPnodeActivate(               /**< activates a leaf node */
   NODE*            node,               /**< leaf node to activate */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPactnodeToFork(              /**< converts the active node into a fork node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPactnodeToSubroot(           /**< converts the active node into a subroot node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );


/*
 * Tree methods
 */

extern
TREE* SCIPtreeCreate(                   /**< creates an initialized tree data structure */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPtreeExtendPath(             /**< appends node to the path of active nodes and marks node active */
   TREE*            tree,               /**< branch-and-bound tree */
   MEM*             mem,                /**< block memory buffers */
   NODE*            node                /**< node to append */
   );

extern
void SCIPtreeShrinkPath(                /**< cuts off path of active nodes after given node and marks them inactive */
   TREE*            tree,               /**< branch-and-bound tree */
   NODE*            node                /**< last node in the path */
   );

#endif
