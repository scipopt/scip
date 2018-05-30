/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file   branch_stp.h
 * @brief  Steiner vertex branching rule
 * @author Daniel Rehfeldt
 *
 * The Steiner branching rule implemented in this file is described in
 * "A Generic Approach to Solving the Steiner Tree Problem and Variants" by Daniel Rehfeldt.
 * It removes includes and exludes Steiner vertices during branching.
 *
*/
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_STP_H__
#define __SCIP_BRANCH_STP_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BRANCH_STP_VERTEX_KILLED      -1
#define BRANCH_STP_VERTEX_NONTERM      0
#define BRANCH_STP_VERTEX_TERM         1


/** parse constraint name and apply changes to graph or array */
EXTERN
SCIP_RETCODE STPStpBranchruleParseConsname(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  vertexchgs,         /**< array to store changes or NULL */
   GRAPH*                graph,              /**< graph to modify or NULL */
   const char*           consname,           /**< constraint name */
   SCIP_Bool             deletehistory       /**< delete history of graph? */
   );

/** applies vertex changes caused by this branching rule, either on a graph or on an array */
EXTERN
SCIP_RETCODE SCIPStpBranchruleApplyVertexChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  vertexchgs,         /**< array to store changes or NULL */
   GRAPH*                graph               /**< graph to apply changes on or NULL */
   );

/** applies vertex changes caused by this branching rule, either on a graph or on an array */
EXTERN
void SCIPStpBranchruleInitNodeState(
   const GRAPH*          g,                  /**< graph data structure */
   int*                  nodestate           /**< node state array */
   );

/** creates the stp branching rule and includes it to SCIP */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleStp(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
