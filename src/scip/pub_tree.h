/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_tree.h,v 1.14 2005/07/15 17:20:16 bzfpfend Exp $"

/**@file   pub_tree.h
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



/*
 * Node methods
 */

/** node comparator for best lower bound */
extern
DECL_SORTPTRCOMP(SCIPnodeCompLowerbound);

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets the type of the node */
extern
NODETYPE SCIPnodeGetType(
   NODE*            node                /**< node */
   );

/** gets successively assigned number of the node */
extern
Longint SCIPnodeGetNumber(
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

/** gets the node selection priority of the node assigned by the branching rule */
extern
Real SCIPnodeGetPriority(
   NODE*            node                /**< node */
   );

/** returns whether node is in the path to the current node */
extern
Bool SCIPnodeIsActive(
   NODE*            node                /**< node */
   );

/** returns whether the node is marked to be propagated again */
extern
Bool SCIPnodeIsPropagatedAgain(
   NODE*            node                /**< node data */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPnodeGetType(node)           ((node)->nodetype)
#define SCIPnodeGetNumber(node)         ((node)->number)
#define SCIPnodeGetDepth(node)          ((node)->depth)
#define SCIPnodeGetLowerbound(node)     ((node)->lowerbound)
#define SCIPnodeGetPriority(node)       ((node)->priority)
#define SCIPnodeIsActive(node)          ((node)->active)
#define SCIPnodeIsPropagatedAgain(node) ((node)->reprop)

#endif



#endif
