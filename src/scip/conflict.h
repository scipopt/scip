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
#pragma ident "@(#) $Id: conflict.h,v 1.29 2005/02/14 13:35:40 bzfpfend Exp $"

/**@file   conflict.h
 * @brief  internal methods for conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONFLICT_H__
#define __CONFLICT_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_conflict.h"
#include "scip/pub_conflict.h"



/*
 * Conflict Handler
 */

/** creates a conflict handler */
extern
RETCODE SCIPconflicthdlrCreate(
   CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   const char*      name,               /**< name of conflict handler */
   const char*      desc,               /**< description of conflict handler */
   int              priority,           /**< priority of the conflict handler */
   DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   DECL_CONFLICTINITSOL((*conflictinitsol)),/**< solving process initialization method of conflict handler */
   DECL_CONFLICTEXITSOL((*conflictexitsol)),/**< solving process deinitialization method of conflict handler */
   DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   );

/** calls destructor and frees memory of conflict handler */
extern
RETCODE SCIPconflicthdlrFree(
   CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SET*             set                 /**< global SCIP settings */
   );

/** calls init method of conflict handler */
extern
RETCODE SCIPconflicthdlrInit(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of conflict handler */
extern
RETCODE SCIPconflicthdlrExit(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SET*             set                 /**< global SCIP settings */
   );

/** informs conflict handler that the branch and bound process is being started */
extern
RETCODE SCIPconflicthdlrInitsol(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SET*             set                 /**< global SCIP settings */
   );

/** informs conflict handler that the branch and bound process data is being freed */
extern
RETCODE SCIPconflicthdlrExitsol(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of conflict handler */
extern
RETCODE SCIPconflicthdlrExec(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SET*             set,                /**< global SCIP settings */
   NODE*            node,               /**< node to add conflict constraint to */
   VAR**            conflictvars,       /**< variables of the conflict set */
   int              nconflictvars,      /**< number of variables in the conflict set */
   Bool             local,              /**< is the conflict set only valid locally? */
   Bool             resolved,           /**< is the conflict set already used to create a constraint? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of conflict handler */
extern
void SCIPconflicthdlrSetPriority(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the conflict handler */
   );




/*
 * Conflict Analysis
 */

/** creates conflict analysis data for propagation conflicts */
extern
RETCODE SCIPconflictCreate(
   CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   SET*             set                 /**< global SCIP settings */
   );

/** frees conflict analysis data for propagation conflicts */
extern
RETCODE SCIPconflictFree(
   CONFLICT**       conflict            /**< pointer to conflict analysis data */
   );

/** initializes the propagation conflict analysis by clearing the conflict candidate queue */
extern
RETCODE SCIPconflictInit(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** adds variable's bound to conflict candidate queue */
extern
RETCODE SCIPconflictAddBound(
   CONFLICT*        conflict,           /**< conflict analysis data */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   );

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set;
 *  updates statistics for propagation conflict analysis
 */
extern
RETCODE SCIPconflictAnalyze(
   CONFLICT*        conflict,           /**< conflict analysis data */
   BLKMEM*          blkmem,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   int              validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/** gets time in seconds used for analyzing propagation conflicts */
extern
Real SCIPconflictGetPropTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to propagation conflict analysis */
extern
Longint SCIPconflictGetNPropCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of conflict clauses detected in propagation conflict analysis */
extern
Longint SCIPconflictGetNPropConflictClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in conflict clauses created in propagation conflict analysis */
extern
Longint SCIPconflictGetNPropConflictLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of reconvergence clauses detected in propagation conflict analysis */
extern
Longint SCIPconflictGetNPropReconvergenceClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in reconvergence clauses created in propagation conflict analysis */
extern
Longint SCIPconflictGetNPropReconvergenceLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   );




/*
 * Infeasible LP Conflict Analysis
 */

/** analyzes an infeasible LP to find out the bound changes on binary variables that were responsible for the 
 *  infeasibility;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible LP conflict analysis
 */
extern
RETCODE SCIPconflictAnalyzeLP(
   CONFLICT*        conflict,           /**< conflict analysis data */
   BLKMEM*          blkmem,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/** gets time in seconds used for analyzing infeasible LP conflicts */
extern
Real SCIPconflictGetLPTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to infeasible LP conflict analysis */
extern
Longint SCIPconflictGetNLPCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of conflict clauses detected in infeasible LP conflict analysis */
extern
Longint SCIPconflictGetNLPConflictClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in conflict clauses created in infeasible LP conflict analysis */
extern
Longint SCIPconflictGetNLPConflictLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of reconvergence clauses detected in infeasible LP conflict analysis */
extern
Longint SCIPconflictGetNLPReconvergenceClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in reconvergence clauses created in infeasible LP conflict analysis */
extern
Longint SCIPconflictGetNLPReconvergenceLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of LP iterations in infeasible LP conflict analysis */
extern
Longint SCIPconflictGetNLPIterations(
   CONFLICT*        conflict            /**< conflict analysis data */
   );




/*
 * infeasible strong branching conflict analysis
 */

/** analyses infeasible strong branching sub problems for conflicts */
extern
RETCODE SCIPconflictAnalyzeStrongbranch(
   CONFLICT*        conflict,           /**< conflict analysis data */
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   COL*             col,                /**< LP column with at least one infeasible strong branching subproblem */
   Bool*            downconflict,       /**< pointer to store whether a conflict clause was created for an infeasible
                                         *   downwards branch, or NULL */
   Bool*            upconflict          /**< pointer to store whether a conflict clause was created for an infeasible
                                         *   upwards branch, or NULL */
   );

/** gets time in seconds used for analyzing infeasible strong branching conflicts */
extern
Real SCIPconflictGetStrongbranchTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to infeasible strong branching conflict analysis */
extern
Longint SCIPconflictGetNStrongbranchCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of conflict clauses detected in infeasible strong branching conflict analysis */
extern
Longint SCIPconflictGetNStrongbranchConflictClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in conflict clauses created in infeasible strong branching conflict analysis */
extern
Longint SCIPconflictGetNStrongbranchConflictLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of reconvergence clauses detected in infeasible strong branching conflict analysis */
extern
Longint SCIPconflictGetNStrongbranchReconvergenceClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in reconvergence clauses created in infeasible strong branching conflict analysis */
extern
Longint SCIPconflictGetNStrongbranchReconvergenceLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of LP iterations in infeasible strong branching conflict analysis */
extern
Longint SCIPconflictGetNStrongbranchIterations(
   CONFLICT*        conflict            /**< conflict analysis data */
   );




/*
 * pseudo solution conflict analysis
 */

/** analyzes a pseudo solution with objective value exceeding the current cutoff to find out the bound changes on
 *  variables that were responsible for the objective value degradation;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for pseudo solution conflict analysis
 */
extern
RETCODE SCIPconflictAnalyzePseudo(
   CONFLICT*        conflict,           /**< conflict analysis data */
   BLKMEM*          blkmem,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/** gets time in seconds used for analyzing pseudo solution conflicts */
extern
Real SCIPconflictGetPseudoTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to pseudo solution conflict analysis */
extern
Longint SCIPconflictGetNPseudoCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of conflict clauses detected in pseudo solution conflict analysis */
extern
Longint SCIPconflictGetNPseudoConflictClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in conflict clauses created in pseudo solution conflict analysis */
extern
Longint SCIPconflictGetNPseudoConflictLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of reconvergence clauses detected in pseudo solution conflict analysis */
extern
Longint SCIPconflictGetNPseudoReconvergenceClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in reconvergence clauses created in pseudo solution conflict analysis */
extern
Longint SCIPconflictGetNPseudoReconvergenceLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   );


#endif
