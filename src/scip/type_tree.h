/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_tree.h,v 1.6 2004/09/07 18:22:21 bzfpfend Exp $"

/**@file   type_tree.h
 * @brief  type definitions for branch and bound tree
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_TREE_H__
#define __TYPE_TREE_H__


enum NodeType
{
   SCIP_NODETYPE_ACTNODE     = 0,       /**< the active node, whose data is stored in the dynamic data object */
   SCIP_NODETYPE_PROBINGNODE = 1,       /**< temporary child node of the active node used for probing */
   SCIP_NODETYPE_SIBLING     = 2,       /**< unsolved sibling of the active node */
   SCIP_NODETYPE_CHILD       = 3,       /**< unsolved child of the active node */
   SCIP_NODETYPE_LEAF        = 4,       /**< unsolved leaf of the tree, stored in the tree's queue */
   SCIP_NODETYPE_DEADEND     = 5,       /**< temporary type of active node, if it was solved completely */
   SCIP_NODETYPE_JUNCTION    = 6,       /**< fork without LP solution */
   SCIP_NODETYPE_FORK        = 7,       /**< fork with solved LP and added rows and columns */
   SCIP_NODETYPE_SUBROOT     = 8        /**< fork with solved LP and arbitrarily changed rows and columns */
};
typedef enum NodeType NODETYPE;         /**< type of node */

typedef struct Child CHILD;             /**< data for child nodes */
typedef struct Sibling SIBLING;         /**< data for sibling nodes */
typedef struct Leaf LEAF;               /**< data for leaf nodes */
typedef struct Junction JUNCTION;       /**< data for junction nodes */
typedef struct Fork FORK;               /**< data for fork nodes */
typedef struct Subroot SUBROOT;         /**< data for subroot nodes */
typedef struct Node NODE;               /**< node data structure */
typedef struct Tree TREE;               /**< branch and bound tree */


#endif
