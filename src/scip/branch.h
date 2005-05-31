/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch.h,v 1.40 2005/05/31 17:20:09 bzfpfend Exp $"

/**@file   branch.h
 * @brief  internal methods for branching rules and branching candidate storage
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BRANCH_H__
#define __BRANCH_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_misc.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_sepastore.h"
#include "scip/type_branch.h"
#include "scip/pub_branch.h"




/*
 * branching candidate storage methods
 */

/** creates a branching candidate storage */
extern
RETCODE SCIPbranchcandCreate(
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   );

/** frees branching candidate storage */
extern
RETCODE SCIPbranchcandFree(
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   );

/** gets branching candidates for LP solution branching (fractional variables) */
extern
RETCODE SCIPbranchcandGetLPCands(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands,           /**< pointer to store the number of LP branching candidates, or NULL */
   int*             npriolpcands        /**< pointer to store the number of candidates with maximal priority, or NULL */
   );

/** gets branching candidates for pseudo solution branching (nonfixed variables) */
extern
RETCODE SCIPbranchcandGetPseudoCands(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands,       /**< pointer to store the number of pseudo branching candidates, or NULL */
   int*             npriopseudocands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   );

/** gets number of branching candidates for pseudo solution branching (nonfixed variables) */
extern
int SCIPbranchcandGetNPseudoCands(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPbranchcandGetNPrioPseudoCands(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPbranchcandGetNPrioPseudoBins(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPbranchcandGetNPrioPseudoInts(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPbranchcandGetNPrioPseudoImpls(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** removes variable from branching candidate list */
extern
RETCODE SCIPbranchcandRemoveVar(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var                 /**< variable that changed its bounds */
   );

/** updates branching candidate list for a given variable */
extern
RETCODE SCIPbranchcandUpdateVar(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< variable that changed its bounds */
   );




/*
 * branching rules
 */

/** creates a branching rule */
extern
RETCODE SCIPbranchruleCreate(
   BRANCHRULE**     branchrule,         /**< pointer to store branching rule */
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   int              maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                         *   compared to best node's dual bound for applying branching rule
                                         *   (0.0: only on current best node, 1.0: on all nodes) */
   DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   DECL_BRANCHINITSOL((*branchinitsol)),/**< solving process initialization method of branching rule */
   DECL_BRANCHEXITSOL((*branchexitsol)),/**< solving process deinitialization method of branching rule */
   DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   );

/** frees memory of branching rule */
extern
RETCODE SCIPbranchruleFree(
   BRANCHRULE**     branchrule,         /**< pointer to branching rule data structure */
   SET*             set                 /**< global SCIP settings */
   );

/** initializes branching rule */
extern
RETCODE SCIPbranchruleInit(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set                 /**< global SCIP settings */
   );

/** deinitializes branching rule */
extern
RETCODE SCIPbranchruleExit(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set                 /**< global SCIP settings */
   );

/** informs branching rule that the branch and bound process is being started */
extern
RETCODE SCIPbranchruleInitsol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set                 /**< global SCIP settings */
   );

/** informs branching rule that the branch and bound process data is being freed */
extern
RETCODE SCIPbranchruleExitsol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set                 /**< global SCIP settings */
   );

/** executes branching rule for fractional LP solution */
extern
RETCODE SCIPbranchruleExecLPSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   SEPASTORE*       sepastore,          /**< separation storage */
   Real             cutoffbound,        /**< global upper cutoff bound */
   Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** executes branching rule for not completely fixed pseudo solution */
extern
RETCODE SCIPbranchruleExecPseudoSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   Real             cutoffbound,        /**< global upper cutoff bound */
   Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of branching rule */
extern
void SCIPbranchruleSetPriority(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the branching rule */
   );

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
extern
void SCIPbranchruleSetMaxdepth(
   BRANCHRULE*      branchrule,         /**< branching rule */
   int              maxdepth            /**< new maxdepth of the branching rule */
   );

/** sets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
extern
void SCIPbranchruleSetMaxbounddist(
   BRANCHRULE*      branchrule,         /**< branching rule */
   Real             maxbounddist        /**< new maxbounddist of the branching rule */
   );



/*
 * branching methods
 */

/** calculates the branching score out of the gain predictions for a binary branching */
extern
Real SCIPbranchGetScore(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   Real             downgain,           /**< prediction of objective gain for rounding downwards */
   Real             upgain              /**< prediction of objective gain for rounding upwards */
   );

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children */
extern
Real SCIPbranchGetScoreMultiple(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   int              nchildren,          /**< number of children that the branching will create */
   Real*            gains               /**< prediction of objective gain for each child */
   );

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 */
extern
RETCODE SCIPbranchExecLP(
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             cutoffbound,        /**< global upper cutoff bound */
   Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN */
extern
RETCODE SCIPbranchExecPseudo(
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             cutoffbound,        /**< global upper cutoff bound */
   Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );


#endif
