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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: vbc.h,v 1.5 2005/01/21 09:17:14 bzfpfend Exp $"

/**@file   vbc.h
 * @brief  methods for VBC Tool output
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __VBC_H__
#define __VBC_H__


#include "def.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_tree.h"
#include "type_vbc.h"



/** creates VBCTool data structure */
extern
RETCODE SCIPvbcCreate(
   VBC**            vbc                 /**< pointer to store the VBC information */
   );

/** frees VBC Tool data structure */
extern
void SCIPvbcFree(
   VBC**            vbc                 /**< pointer to store the VBC information */
   );

/** initializes VBC information and creates a file for VBC output */
extern
RETCODE SCIPvbcInit(
   VBC*             vbc,                /**< VBC information */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   );

/** closes the VBC output file */
extern
void SCIPvbcExit(
   VBC*             vbc,                /**< VBC information */
   SET*             set                 /**< global SCIP settings */
   );

/** creates a new node entry in the VBC output file */
extern
RETCODE SCIPvbcNewChild(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of solved nodes */
extern
void SCIPvbcSolvedNode(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of cutoff nodes */
extern
void SCIPvbcCutoffNode(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of nodes where a conflict clause was found */
extern
void SCIPvbcFoundConflict(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of nodes that were marked to be repropagated */
extern
void SCIPvbcMarkedRepropagateNode(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of repropagated nodes */
extern
void SCIPvbcRepropagatedNode(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of nodes with a primal solution */
extern
void SCIPvbcFoundSolution(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   );

/** outputs a new global lower bound to the VBC output file */
extern
void SCIPvbcLowerbound(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   Real             lowerbound          /**< new lower bound */
   );

/** outputs a new global upper bound to the VBC output file */
extern
void SCIPvbcUpperbound(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   Real             upperbound          /**< new upper bound */
   );


#endif
