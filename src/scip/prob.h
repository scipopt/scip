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
#pragma ident "@(#) $Id: prob.h,v 1.46 2005/05/31 17:20:18 bzfpfend Exp $"

/**@file   prob.h
 * @brief  internal methods for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PROB_H__
#define __PROB_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_branch.h"
#include "scip/type_cons.h"

#include "scip/struct_prob.h"



/*
 * problem creation
 */

/** creates problem data structure
 *  If the problem type requires the use of variable pricers, these pricers should be activated with calls
 *  to SCIPactivatePricer(). These pricers are automatically deactivated, when the problem is freed.
 */
extern
RETCODE SCIPprobCreate(
   PROB**           prob,               /**< pointer to problem data structure */
   BLKMEM*          blkmem,             /**< block memory */
   const char*      name,               /**< problem name */
   DECL_PROBDELORIG ((*probdelorig)),   /**< frees user data of original problem */
   DECL_PROBTRANS   ((*probtrans)),     /**< creates user data of transformed problem by transforming original user data */
   DECL_PROBDELTRANS((*probdeltrans)),  /**< frees user data of transformed problem */
   DECL_PROBINITSOL ((*probinitsol)),   /**< solving process initialization method of transformed data */
   DECL_PROBEXITSOL ((*probexitsol)),   /**< solving process deinitialization method of transformed data */
   PROBDATA*        probdata,           /**< user problem data set by the reader */
   Bool             transformed         /**< is this the transformed problem? */
   );

/** frees problem data structure */
extern
RETCODE SCIPprobFree(
   PROB**           prob,               /**< pointer to problem data structure */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   LP*              lp                  /**< current LP data (or NULL, if it's the original problem) */
   );

/** transform problem data into normalized form */
extern
RETCODE SCIPprobTransform(
   PROB*            source,             /**< problem to transform */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   PROB**           target              /**< pointer to target problem data structure */
   );

/** resets the global and local bounds of original variables in original problem to their original values */
extern
RETCODE SCIPprobResetBounds(
   PROB*            prob,               /**< original problem data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   );




/*
 * problem modification
 */

/** sets user problem data */
extern
void SCIPprobSetData(
   PROB*            prob,               /**< problem */
   PROBDATA*        probdata            /**< user problem data to use */
   );

/** adds variable to the problem and captures it */
extern
RETCODE SCIPprobAddVar(
   PROB*            prob,               /**< problem data */
   BLKMEM*          blkmem,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var                 /**< variable to add */
   );

/** marks variable to be removed from the problem; however, the variable is NOT removed from the constraints */
extern
RETCODE SCIPprobDelVar(
   PROB*            prob,               /**< problem data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var                 /**< problem variable */
   );

/** actually removes the deleted variables from the problem and releases them */
extern
RETCODE SCIPprobPerformVarDeletions(
   PROB*            prob,               /**< problem data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data (may be NULL) */
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** changes the type of a variable in the problem */
extern
RETCODE SCIPprobChgVarType(
   PROB*            prob,               /**< problem data */
   SET*             set,                /**< global SCIP settings */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var,                /**< variable to add */
   VARTYPE          vartype             /**< new type of variable */
   );

/** informs problem, that the given loose problem variable changed its status */
extern
RETCODE SCIPprobVarChangedStatus(
   PROB*            prob,               /**< problem data */
   SET*             set,                /**< global SCIP settings */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var                 /**< problem variable */
   );

/** adds constraint to the problem and captures it;
 *  a local constraint is automatically upgraded into a global constraint
 */
extern
RETCODE SCIPprobAddCons(
   PROB*            prob,               /**< problem data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   CONS*            cons                /**< constraint to add */
   );

/** releases and removes constraint from the problem; if the user has not captured the constraint for his own use, the
 *  constraint may be invalid after the call
 */
extern
RETCODE SCIPprobDelCons(
   PROB*            prob,               /**< problem data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   CONS*            cons                /**< constraint to remove */
   );

/** remembers the current number of constraints in the problem's internal data structure
 *  - resets maximum number of constraints to current number of constraints
 *  - remembers current number of constraints as starting number of constraints
 */
extern
void SCIPprobMarkNConss(
   PROB*            prob                /**< problem data */
   );

/** sets objective sense: minimization or maximization */
extern
void SCIPprobSetObjsense(
   PROB*            prob,               /**< problem data */
   OBJSENSE         objsense            /**< new objective sense */
   );

/** adds value to objective offset */
extern
void SCIPprobAddObjoffset(
   PROB*            prob,               /**< problem data */
   Real             addval              /**< value to add to objective offset */
   );

/** sets limit on objective function, such that only solutions better than this limit are accepted */
extern
void SCIPprobSetObjlim(
   PROB*            prob,               /**< problem data */
   Real             objlim              /**< external objective limit */
   );

/** informs the problem, that its objective value is always integral in every feasible solution */
extern
void SCIPprobSetObjIntegral(
   PROB*            prob                /**< problem data */
   );

/** sets integral objective value flag, if all variables with non-zero objective values are integral and have 
 *  integral objective value
 */
extern
void SCIPprobCheckObjIntegral(
   PROB*            prob,               /**< problem data */
   SET*             set                 /**< global SCIP settings */
   );

/** remembers the current solution as root solution in the problem variables */
extern
void SCIPprobStoreRootSol(
   PROB*            prob,               /**< problem data */
   Bool             roothaslp           /**< is the root solution from LP? */
   );

/** informs problem, that the presolving process was finished, and updates all internal data structures */
extern
RETCODE SCIPprobExitPresolve(
   PROB*            prob,               /**< problem data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            infeasible          /**< pointer to store TRUE, if an infeasibility was detected */
   );

/** initializes problem for branch and bound process */
extern
RETCODE SCIPprobInitSolve(
   PROB*            prob,               /**< problem data */
   SET*             set                 /**< global SCIP settings */
   );

/** deinitializes problem after branch and bound process, and converts all COLUMN variables back into LOOSE variables */
extern
RETCODE SCIPprobExitSolve(
   PROB*            prob,               /**< problem data */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   );




/*
 * problem information
 */

/** gets problem name */
extern
const char* SCIPprobGetName(
   PROB*            prob                /**< problem data */
   );

/** gets user problem data */
extern
PROBDATA* SCIPprobGetData(
   PROB*            prob                /**< problem */
   );

/** returns the external value of the given internal objective value */
extern
Real SCIPprobExternObjval(
   PROB*            prob,               /**< problem data */
   SET*             set,                /**< global SCIP settings */
   Real             objval              /**< internal objective value */
   );

/** returns the internal value of the given external objective value */
extern
Real SCIPprobInternObjval(
   PROB*            prob,               /**< problem data */
   SET*             set,                /**< global SCIP settings */
   Real             objval              /**< external objective value */
   );

/** gets limit on objective function in external space */
extern
Real SCIPprobGetObjlim(
   PROB*            prob,               /**< problem data */
   SET*             set                 /**< global SCIP settings */
   );

/** returns whether the objective value is known to be integral in every feasible solution */
extern
Bool SCIPprobIsObjIntegral(
   PROB*            prob                /**< problem data */
   );

/** returns variable of the problem with given name */
extern
VAR* SCIPprobFindVar(
   PROB*            prob,               /**< problem data */
   const char*      name                /**< name of variable to find */
   );

/** returns constraint of the problem with given name */
extern
CONS* SCIPprobFindCons(
   PROB*            prob,               /**< problem data */
   const char*      name                /**< name of variable to find */
   );

/** returns TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing
 */
extern
Bool SCIPprobAllColsInLP(
   PROB*            prob,               /**< problem data */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data */
   );

/** displays current pseudo solution */
extern
void SCIPprobPrintPseudoSol(
   PROB*            prob,               /**< problem data */
   SET*             set                 /**< global SCIP settings */
   );

/** outputs problem statistics */
extern
void SCIPprobPrintStatistics(
   PROB*            prob,               /**< problem data */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** outputs problem to file stream */
extern
RETCODE SCIPprobPrint(
   PROB*            prob,               /**< problem data */
   SET*             set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   );


#endif
