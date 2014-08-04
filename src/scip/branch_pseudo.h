/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_Pseudo.h
 * @ingroup BRANCHINGRULES
 * @brief  Pseudo branching rule
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_PSEUDE_H__
#define __SCIP_BRANCH_PSEUDO_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct LogicOrData LOGICORDATA;

enum Reopt_ConsType
{
   REOPT_CONSTYPE_SEPASOLUTION = 0,
   REOPT_CONSTYPE_INFSUBTREE   = 1,
   REOPT_CONSTYPE_STRBRANCHED  = 2
};
typedef enum Reopt_ConsType REOPT_CONSTYPE;

/** creates the Pseudo branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchrulePseudo(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoAddPseudoVar(
   SCIP*                 scip,
   SCIP_NODE*            node,
   SCIP_VAR*             var,
   SCIP_BOUNDTYPE        boundtype,
   int                   newbound
   );

extern
int SCIPbranchrulePseudoGetNPseudoVars(
   SCIP*                 scip,
   int                   externID
   );

extern
void SCIPbranchrulePseudoDisable(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoSetParams(
   SCIP*                 scip,
   SCIP_Bool             reopt
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoReset(
   SCIP*                 scip,
   SCIP_Bool             restart
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoNodeFinished(
   SCIP*                 scip,
   SCIP_NODE*            node,
   REOPT_CONSTYPE        constype
   );

/*
 * delete branching information for a given nodeID
 */
extern
SCIP_RETCODE SCIPbranchrulePseudoDelInformation(
   SCIP*                 scip,
   int                   nodeID
   );

extern
SCIP_Bool SCIPbranchrulePseudoIsPseudoBranched(
   SCIP*                 scip,
   SCIP_NODE*            node
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoDeleteLastNodeInfo(
   SCIP*                 scip,
   SCIP_NODE*            node
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoGetNextPseudoNode(
   SCIP*                 scip,
   SCIP_VAR**            var,
   SCIP_Real*            bound,
   SCIP_BOUNDTYPE*       boundtype,
   int*                  nbranchvars,
   int                   nbranchvarsize
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoGetNextPseudoCons(
   SCIP*                 scip,
   SCIP_Bool*            consadded
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoGetNextPseudoFixing(
   SCIP*                 scip,
   SCIP_VAR*             var,
   double*               bound
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoFeasLeaf(
   SCIP*                 scip,
   SCIP_Bool             feasible
   );

extern
int SCIPbranchrulePseudoGetNPseudovarsNextPseudoNode(
   SCIP*                 scip
   );

int SCIPbranchrulePseudoGetNPseudoNodes(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoAddGlobalCons(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoGenerateCons(
   SCIP*                 scip,
   LOGICORDATA*          consdata,
   int*                  nvars,
   int                   nallocvars,
   int                   ID,
   SCIP_Bool             local,
   SCIP_Bool             cleardata
   );

extern
SCIP_Bool SCIPbranchrulePseudoIsRootPseudo(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoSaveParentID(
   SCIP*                 scip,
   SCIP_NODE*            node
   );

extern
int SCIPbranchrulePseudoGetID(
   SCIP*                 scip,
   int                   nodenr
   );

extern
SCIP_RETCODE SCIPbranchrulePseudoLinkIDs(
   SCIP*                 scip,
   int                   ID
   );

#ifdef __cplusplus
}
#endif

#endif
