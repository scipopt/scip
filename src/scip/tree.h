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
   SCIP_NODETYPE_LEAF     = 0,          /**< unsolved leaf of the tree */
   SCIP_NODETYPE_ACTNODE  = 1,          /**< the active node, whose data is stored in the dynamic data object */
   SCIP_NODETYPE_JUNCTION = 2,          /**< fork without LP solution */
   SCIP_NODETYPE_FORK     = 3,          /**< fork with solved LP and added rows and columns */
   SCIP_NODETYPE_SUBROOT  = 4           /**< fork with solved LP and arbitrarily changed rows and columns */
};

typedef enum Nodetype NODETYPE;         /**< type of node */
typedef struct Actnode ACTNODE;         /**< data for the actual node */
typedef struct Junction JUNCTION;       /**< data for junction nodes */
typedef struct Fork FORK;               /**< data for fork nodes */
typedef struct Subroot SUBROOT;         /**< data for subroot nodes */
typedef struct Node NODE;               /**< node data structure */
typedef struct Tree TREE;               /**< branch and bound tree */



#include "lp.h"



/*
 * Node methods
 */

extern
DECL_SORTPTRCOMP(SCIPnodeCmpLowerbound);/**< node comparator for best lower bound */

extern
RETCODE SCIPnodeCreate(                 /**< creates a child node of the active node */
   NODE**           node,               /**< pointer to node data structure */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
RETCODE SCIPnodeFree(                   /**< frees node */
   NODE**           node,               /**< node data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPnodeActivate(               /**< activates a leaf node */
   NODE*            node,               /**< leaf node to activate */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
RETCODE SCIPnodeToJunction(             /**< converts the active node into a junction node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
RETCODE SCIPnodeToFork(                 /**< converts the active node into a fork node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPnodeToSubroot(              /**< converts the active node into a subroot node */
   NODE*            node,               /**< node to convert */
   MEM*             mem,                /**< block memory buffers */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   );



/*
 * Tree methods
 */

extern
RETCODE SCIPtreeCreate(                 /**< creates an initialized tree data structure */
   TREE**           tree,               /**< pointer to tree data structure */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPtreeLoadLP(                 /**< constructs the LP and loads LP state for fork/subroot of the active node */
   TREE*            tree,               /**< branch-and-bound tree */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );


#endif
