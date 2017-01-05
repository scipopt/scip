/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_nodereopt.h
 * @ingroup BRANCHINGRULES
 * @brief  nodereopt branching rule
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_NODEREOPT_H__
#define __SCIP_BRANCH_NODEREOPT_H__


#include "scip/scip.h"
#include "scip/type_branch.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the nodereopt branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleNodereopt(
   SCIP*                   scip                /**< SCIP data structure */
   );

extern
SCIP_RETCODE SCIPbranchruleNodereoptInfNode(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_NODE*              cnode            /**< current node */
   );

/*
 * Remove node from data structure
 */
extern
SCIP_RETCODE SCIPbranchruleNodereoptRemoveNode(
   SCIP*                 scip,
   SCIP_NODE*            node,
   SCIP_Bool             branched,
   SCIP_Bool             infeasible
   );

extern
SCIP_RETCODE SCIPbranchruleNodereoptAddConfCons(
   SCIP*                 scip              /**< SCIP data structure */
   );

extern
SCIP_RETCODE SCIPbranchruleNodereoptForceRestart(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchnodereoptCheckFeasibility(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchruleNodereoptSetRootLPI(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchruleNodereoptAddNode(
   SCIP*                 scip,
   SCIP_NODE*            node,               /** current node */
   SCIP_REOPTTYPE        reopttype,
   SCIP_Bool             saveafterdual
   );

extern
SCIP_RETCODE SCIPbranchruleNodereoptAddGlobalCons(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchruleNodereoptRestartCheck(
   SCIP*                 scip
   );

extern
/* get all statistic information */
SCIP_RETCODE SCIPbranchruleNodereoptGetStatistic(
   SCIP*                 scip,
   int*                  nfeasnodes,
   int*                  nfeasnodesround,
   int*                  ninfeasnodes,
   int*                  ninfeasnodesround,
   int*                  nprunednodes,
   int*                  nprunednodesround,
   int*                  nrediednodes,
   int*                  nrediednodesround,
   int*                  nruns,
   int*                  nrestarts,
   int*                  firstrestart,
   int*                  lastrestart,
   int*                  nrestartsround,
   int*                  ninfsubtrees,
   int*                  mfscalls,
   int*                  mfssucces,
   int*                  lisk,
   SCIP_Real*            mfstime,
   int*                  rtfcalls,
   int*                  rtfsucces,
   int*                  lck,
   SCIP_Real*            rtftime,
   SCIP_Real*            lptime,
   SCIP_Real*            lptime_maxthread
   );

extern
/* start clock for update solutions */
SCIP_RETCODE SCIPbranchruleNodereoptStartUpdatesoluTime(
   SCIP*                 scip
   );

extern
/* stop clock for update solutions */
SCIP_RETCODE SCIPbranchruleNodereoptStopUpdatesoluTime(
   SCIP*                 scip
   );

extern
SCIP_Real SCIPbranchruleNodereoptGetCutoffbound(
   SCIP*                 scip,
   int                   nodeID
   );

extern
SCIP_RETCODE SCIPbranchruledataNodereoptLPTimeStart(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchruledataNodereoptLPTimeStop(
   SCIP*                 scip
   );

extern
SCIP_RETCODE SCIPbranchruleNodereoptSolveLP(
   SCIP*                 scip,
   SCIP_NODE*            node,
   SCIP_Bool*            solvelp
   );

/*
 * Save global constraint to separate solution
 */
extern
SCIP_RETCODE SCIPbranchruleNodereoptSaveGlobaleCons(
   SCIP*                 scip,
   SCIP_SOL*             sol,
   SCIP_SET*             set,
   SCIP_STAT*            stat
   );

extern
int SCIPbranchruleNodereoptGetNAddedConss(
   SCIP*                 scip,
   SCIP_NODE*            node
   );

/*
 * save unexplored nodes
 */
extern
SCIP_RETCODE SCIPbranchruleNodereoptSaveOpenNodes(
   SCIP*                 scip
   );

#ifdef __cplusplus
}
#endif

#endif
