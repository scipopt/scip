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
   SCIP_NODETYPE_LEAF    = 0,           /**< unsolved child node without children */
   SCIP_NODETYPE_FORK    = 1            /**< solved parent node with children */
};

typedef enum Nodetype NODETYPE;         /**< type of node */
typedef struct Leaf LEAF;               /**< data for leaf nodes */
typedef struct Fork FORK;               /**< data for fork nodes */
typedef struct Node NODE;               /**< node data structure */
typedef struct Tree TREE;               /**< branch and bound tree */


extern
NODE* SCIPcreateLeaf(                   /**< creates a leaf node */
   NODE*            parent,             /**< parent node in the tree */
   BASIS*           basis               /**< pointer to LP basis information */
   );


#endif
