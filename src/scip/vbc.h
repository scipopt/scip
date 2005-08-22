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
#pragma ident "@(#) $Id: vbc.h,v 1.10 2005/08/22 18:35:55 bzfpfend Exp $"

/**@file   vbc.h
 * @brief  methods for SCIP_VBC Tool output
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VBC_H__
#define __SCIP_VBC_H__


#include "scip/def.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"
#include "scip/type_vbc.h"



/** creates VBCTool data structure */
extern
SCIP_RETCODE SCIPvbcCreate(
   SCIP_VBC**            vbc                 /**< pointer to store the SCIP_VBC information */
   );

/** frees SCIP_VBC Tool data structure */
extern
void SCIPvbcFree(
   SCIP_VBC**            vbc                 /**< pointer to store the SCIP_VBC information */
   );

/** initializes SCIP_VBC information and creates a file for SCIP_VBC output */
extern
SCIP_RETCODE SCIPvbcInit(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** closes the SCIP_VBC output file */
extern
void SCIPvbcExit(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** creates a new node entry in the SCIP_VBC output file */
extern
SCIP_RETCODE SCIPvbcNewChild(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of solved nodes */
extern
void SCIPvbcSolvedNode(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of cutoff nodes */
extern
void SCIPvbcCutoffNode(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of nodes where a conflict clause was found */
extern
void SCIPvbcFoundConflict(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of nodes that were marked to be repropagated */
extern
void SCIPvbcMarkedRepropagateNode(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of repropagated nodes */
extern
void SCIPvbcRepropagatedNode(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** changes the color of the node to the color of nodes with a primal solution */
extern
void SCIPvbcFoundSolution(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** outputs a new global lower bound to the SCIP_VBC output file */
extern
void SCIPvbcLowerbound(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             lowerbound          /**< new lower bound */
   );

/** outputs a new global upper bound to the SCIP_VBC output file */
extern
void SCIPvbcUpperbound(
   SCIP_VBC*             vbc,                /**< SCIP_VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             upperbound          /**< new upper bound */
   );


#endif
