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
#pragma ident "@(#) $Id: pub_tree.h,v 1.3 2004/02/04 17:27:38 bzfpfend Exp $"

/**@file   pub_tree.h
 * @brief  public methods for branch-and-bound tree
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_TREE_H__
#define __PUB_TREE_H__


#include "def.h"
#include "type_misc.h"
#include "type_tree.h"

#ifdef NDEBUG
#include "struct_tree.h"
#endif



/*
 * Node methods
 */

/** node comparator for best lower bound */
extern
DECL_SORTPTRCOMP(SCIPnodeCmpLowerbound);

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



#endif
