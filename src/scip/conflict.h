/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   conflict.h
 * @brief  methods and datastructures for conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONFLICT_H__
#define __CONFLICT_H__


typedef struct Conflict CONFLICT;       /**< conflict analysis data structure */


#include "def.h"
#include "retcode.h"
#include "memory.h"
#include "set.h"
#include "stat.h"
#include "var.h"



/** creates conflict analysis data */
extern
RETCODE SCIPconflictCreate(
   CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   const SET*       set                 /**< global SCIP settings */
   );

/** frees conflict analysis data */
extern
RETCODE SCIPconflictFree(
   CONFLICT**       conflict            /**< pointer to conflict analysis data */
   );

/** initializes the conflict analysis by clearing the conflict variable candidate queue */
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

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and returns a conflict set, that
 *  can be used to create a conflict constraint
 */
extern
RETCODE SCIPconflictAnalyze(
   CONFLICT*        conflict,           /**< conflict analysis data */
   const SET*       set,                /**< global SCIP settings */
   int              maxsize,            /**< maximal size of the conflict set or -1 for no restriction */
   VAR***           conflictvars,       /**< pointer to store the conflict set (user must not change this array) */
   int*             nconflictvars,      /**< pointer to store the number of conflict variables */
   Bool*            success             /**< pointer to store whether the conflict set is valid */
   );

/** gets time in seconds used for analyzing conflicts */
extern
Real SCIPconflictGetTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to conflict analysis */
extern
Longint SCIPconflictGetNCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of valid conflicts detected in conflict analysis */
extern
Longint SCIPconflictGetNConflicts(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

#endif
