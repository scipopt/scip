/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: prob.h,v 1.26 2003/12/15 17:45:33 bzfpfend Exp $"

/**@file   prob.h
 * @brief  internal methods for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PROB_H__
#define __PROB_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_prob.h"
#include "type_tree.h"
#include "type_branch.h"
#include "type_cons.h"

#include "struct_prob.h"



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
   const char*      name,               /**< problem name */
   DECL_PROBDELORIG ((*probdelorig)),   /**< frees user data of original problem */
   DECL_PROBTRANS   ((*probtrans)),     /**< creates user data of transformed problem by transforming original user data */
   DECL_PROBDELTRANS((*probdeltrans)),  /**< frees user data of transformed problem */
   PROBDATA*        probdata,           /**< user problem data set by the reader */
   Bool             transformed         /**< is this the transformed problem? */
   );

/** frees problem data structure */
extern
RETCODE SCIPprobFree(
   PROB**           prob,               /**< pointer to problem data structure */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's the original problem) */
   );

/** transform problem data into normalized form */
extern
RETCODE SCIPprobTransform(
   PROB*            source,             /**< problem to transform */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   PROB**           target              /**< pointer to target problem data structure */
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
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var                 /**< variable to add */
   );

/** changes the type of a variable in the problem */
extern
RETCODE SCIPprobChgVarType(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var,                /**< variable to add */
   VARTYPE          vartype             /**< new type of variable */
   );

/** informs problem, that the given loose problem variable changed its status */
extern
RETCODE SCIPprobVarChangedStatus(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var                 /**< problem variable */
   );

/** adds constraint to the problem and captures it; a local constraint is automatically upgraded into a global constraint */
extern
RETCODE SCIPprobAddCons(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   );

/** releases and removes constraint from the problem; if the user has not captured the constraint for his own use, the
 *  constraint may be invalid after the call
 */
extern
RETCODE SCIPprobDelCons(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to remove */
   );

/** resets maximum number of constraints to current number of constraints, remembers actual number of constraints
 *  as starting number of constraints
 */
extern
void SCIPprobSolvingStarts(
   PROB*            prob                /**< problem data */
   );

/** sets objective sense: minimization or maximization */
extern
void SCIPprobSetObjsense(
   PROB*            prob,               /**< problem data */
   OBJSENSE         objsense            /**< new objective sense */
   );

/** increases objective offset */
extern
void SCIPprobIncObjoffset(
   PROB*            prob,               /**< problem data */
   Real             incval              /**< value to add to objective offset */
   );

/** sets limit on objective function, such that only solutions better than this limit are accepted */
extern
void SCIPprobSetExternObjlim(
   PROB*            prob,               /**< problem data */
   Real             objlim              /**< external objective limit */
   );

/** sets limit on objective function as transformed internal objective value */
extern
void SCIPprobSetInternObjlim(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   Real             objlim              /**< transformed internal objective limit */
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
   const SET*       set,                /**< global SCIP settings */
   Real             objval              /**< internal objective value */
   );

/** returns the internal value of the given external objective value */
extern
Real SCIPprobInternObjval(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   Real             objval              /**< external objective value */
   );

/** gets limit on objective function in external space */
extern
Real SCIPprobGetExternObjlim(
   PROB*            prob                /**< problem data */
   );

/** gets limit on objective function as transformed internal objective value */
extern
Real SCIPprobGetInternObjlim(
   PROB*            prob,               /**< problem data */
   const SET*       set                 /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

/** displays actual pseudo solution */
extern
void SCIPprobPrintPseudoSol(
   PROB*            prob,               /**< problem data */
   const SET*       set                 /**< global SCIP settings */
   );

/** outputs problem statistics */
extern
void SCIPprobPrintStatistics(
   PROB*            prob,               /**< problem data */
   FILE*            file                /**< output file (or NULL for standard output) */
   );



#endif
