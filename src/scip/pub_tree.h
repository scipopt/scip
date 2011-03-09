/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_tree.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for branch and bound tree
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_TREE_H__
#define __SCIP_PUB_TREE_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_tree.h"

#ifdef NDEBUG
#include "scip/struct_tree.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Node methods
 */

/** node comparator for best lower bound */
extern
SCIP_DECL_SORTPTRCOMP(SCIPnodeCompLowerbound);

/** returns the set of variable branchings that were performed in the parent node to create this node */
extern
void SCIPnodeGetParentBranchings(
   SCIP_NODE*            node,                /**< node data */
   SCIP_VAR**            branchvars,          /**< array of variables on which the branching has been performed in the parent node */
   SCIP_Real*            branchbounds,        /**< array of bounds which the branching in the parent node set */
   SCIP_BOUNDTYPE*       boundtypes,          /**< array of boundtypes which the branching in the parent node set */
   int*                  nbranchvars,         /**< number of variables on which branching has been performed in the parent node 
                                               *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int                   branchvarssize       /**< available slots in arrays */
   );

/** returns the set of variable branchings that were performed in all ancestor nodes (nodes on the path to the root) to create this node */
extern
void SCIPnodeGetAncestorBranchings(
   SCIP_NODE*            node,                /**< node data */
   SCIP_VAR**            branchvars,          /**< array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real*            branchbounds,        /**< array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes,          /**< array of boundtypes which the branchings in all ancestors set */
   int*                  nbranchvars,         /**< number of variables on which branchings have been performed in all ancestors 
                                               *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int                   branchvarssize       /**< available slots in arrays */
   );

/*  returns the set of variable branchings that were performed in all ancestor nodes (nodes on the path to the root) to create this node 
 *  sorted by the nodes, starting from the current node going up to the root */
extern
void SCIPnodeGetAncestorBranchingPath(
   SCIP_NODE*            node,                /**< node data */
   SCIP_VAR**            branchvars,          /**< array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real*            branchbounds,        /**< array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes,          /**< array of boundtypes which the branchings in all ancestors set */
   int*                  nbranchvars,         /**< number of variables on which branchings have been performed in all ancestors 
                                               *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int                   branchvarssize,      /**< available slots in arrays */   
   int*                  nodeswitches,        /**< marks, where in the arrays the branching decisions of the next node on the path start 
                                               * branchings performed at the parent of node always start at position 0. For single variable branching,
                                               * nodeswitches[i] = i holds
                                               */
   int*                  nnodes,              /**< number of nodes in the nodeswitch array */
   int                   nodeswitchsize       /**< available slots in node switch array */   
   );


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets the type of the node */
extern
SCIP_NODETYPE SCIPnodeGetType(
   SCIP_NODE*            node                /**< node */
   );

/** gets successively assigned number of the node */
extern
SCIP_Longint SCIPnodeGetNumber(
   SCIP_NODE*            node                /**< node */
   );

/** gets the depth of the node */
extern
int SCIPnodeGetDepth(
   SCIP_NODE*            node                /**< node */
   );

/** gets the lower bound of the node */
extern
SCIP_Real SCIPnodeGetLowerbound(
   SCIP_NODE*            node                /**< node */
   );

/** gets the estimated value of the best feasible solution in subtree of the node */
extern
SCIP_Real SCIPnodeGetEstimate(
   SCIP_NODE*            node                /**< node */
   );

/** gets the domain change information of the node, i.e., the information about the differences in the
 *  variables domains to the parent node
 */
extern
SCIP_DOMCHG* SCIPnodeGetDomchg(
   SCIP_NODE*            node                /**< node */
   );


/** returns whether node is in the path to the current node */
extern
SCIP_Bool SCIPnodeIsActive(
   SCIP_NODE*            node                /**< node */
   );

/** returns whether the node is marked to be propagated again */
extern
SCIP_Bool SCIPnodeIsPropagatedAgain(
   SCIP_NODE*            node                /**< node data */
   );



#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPnodeGetType(node)           ((SCIP_NODETYPE)(node)->nodetype)
#define SCIPnodeGetNumber(node)         ((node)->number)
#define SCIPnodeGetDepth(node)          ((node)->depth)
#define SCIPnodeGetLowerbound(node)     ((node)->lowerbound)
#define SCIPnodeGetEstimate(node)       ((node)->estimate)
#define SCIPnodeGetDomchg(node)         ((node)->domchg)
#define SCIPnodeIsActive(node)          ((node)->active)
#define SCIPnodeIsPropagatedAgain(node) ((node)->reprop)

#endif

#ifdef __cplusplus
}
#endif

#endif
