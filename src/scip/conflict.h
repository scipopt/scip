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
#pragma ident "@(#) $Id: conflict.h,v 1.15 2004/04/27 15:49:57 bzfpfend Exp $"

/**@file   conflict.h
 * @brief  internal methods for conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONFLICT_H__
#define __CONFLICT_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_prob.h"
#include "type_tree.h"
#include "type_conflict.h"
#include "type_scip.h"
#include "pub_conflict.h"



/*
 * Conflict Handler
 */

/** creates a conflict handler */
extern
RETCODE SCIPconflicthdlrCreate(
   CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of conflict handler */
   const char*      desc,               /**< description of conflict handler */
   int              priority,           /**< priority of the conflict handler */
   DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   );

/** calls destructor and frees memory of conflict handler */
extern
RETCODE SCIPconflicthdlrFree(
   CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls init method of conflict handler */
extern
RETCODE SCIPconflicthdlrInit(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls exit method of conflict handler */
extern
RETCODE SCIPconflicthdlrExit(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls execution method of conflict handler */
extern
RETCODE SCIPconflicthdlrExec(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP*            scip,               /**< SCIP data structure */   
   NODE*            node,               /**< node to add conflict constraint to */
   VAR**            conflictvars,       /**< variables of the conflict set */
   int              nconflictvars,      /**< number of variables in the conflict set */
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
   const SET*       set                 /**< global SCIP settings */
   );

/** frees conflict analysis data for propagation conflicts */
extern
RETCODE SCIPconflictFree(
   CONFLICT**       conflict            /**< pointer to conflict analysis data */
   );

/** initializes the propagation conflict analysis by clearing the conflict variable candidate queue */
extern
RETCODE SCIPconflictInit(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** adds variable to conflict variable candidates */
extern
RETCODE SCIPconflictAddVar(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR*             var                 /**< problem variable */
   );

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set;
 *  updates statistics for propagation conflict analysis
 */
extern
RETCODE SCIPconflictAnalyze(
   CONFLICT*        conflict,           /**< conflict analysis data */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
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

/** gets number of valid conflicts detected in propagation conflict analysis */
extern
Longint SCIPconflictGetNPropConflicts(
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
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
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

/** gets number of valid conflicts detected in infeasible LP conflict analysis */
extern
Longint SCIPconflictGetNLPConflicts(
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
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   COL*             col,                /**< LP column with at least one infeasible strong branching subproblem */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
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

/** gets number of valid conflicts detected in infeasible strong branching conflict analysis */
extern
Longint SCIPconflictGetNStrongbranchConflicts(
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

/** analyzes a pseudo solution with objective value exceeding the current cutoff to find out the bound changes on binary
 *  variables that were responsible for the objective value degradation;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for pseudo solution conflict analysis
 */
extern
RETCODE SCIPconflictAnalyzePseudo(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
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

/** gets number of valid conflicts detected in pseudo solution conflict analysis */
extern
Longint SCIPconflictGetNPseudoConflicts(
   CONFLICT*        conflict            /**< conflict analysis data */
   );



#endif
