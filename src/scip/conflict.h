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


typedef struct ConflictHdlr CONFLICTHDLR; /**< conflict handler to process conflict sets */
typedef struct ConflictHdlrData CONFLICTHDLRDATA; /**< conflict handler data */
typedef struct Conflict CONFLICT;       /**< conflict analysis data structure for propagation conflicts */
typedef struct LPConflict LPCONFLICT;   /**< conflict analysis data structure for infeasible LP conflicts */


/** destructor of conflict handler to free conflict handler data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conflicthdlr    : the conflict handler itself
 */
#define DECL_CONFLICTFREE(x) RETCODE x (SCIP* scip, CONFLICTHDLR* conflicthdlr)

/** initialization method of conflict handler (called when problem solving starts)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conflicthdlr    : the conflict handler itself
 */
#define DECL_CONFLICTINIT(x) RETCODE x (SCIP* scip, CONFLICTHDLR* conflicthdlr)

/** deinitialization method of conflict handler (called when problem solving exits)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conflicthdlr    : the conflict handler itself
 */
#define DECL_CONFLICTEXIT(x) RETCODE x (SCIP* scip, CONFLICTHDLR* conflicthdlr)

/** conflict processing method of conflict handler (called when conflict was found)
 *
 *  This method is called, when the conflict analysis found a conflict on binary variable assignments.
 *  The conflict handler may update its data accordingly and create a constraint out of the conflict.
 *  If the parameter "resolved" is set, the conflict handler should not create a constraint, because
 *  a different conflict handler with higher priority already created a constraint.
 *  The variables in the conflict set lead to a conflict (i.e. an infeasibility) when all set to FALSE.
 *  Thus, a feasible conflict constraint must demand, that at least one of the variables in the conflict
 *  set is set to TRUE.
 *  The given "conflictvars" array representing the conflict set is only a reference to an internal
 *  buffer, that may be modified at any time by SCIP. The user must copy the needed information from the
 *  "conflictvars" array to its own data structures, if he wants to use the information later.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conflicthdlr    : the conflict handler itself
 *  - conflictvars    : array with binary variables in the conflict set
 *  - nconflictvars   : number of binary variables in the conflict set
 *  - resolved        : is the conflict set already used to create a constraint?
 *  - result          : pointer to store the result of the conflict processing call
 *
 *  possible return values for *result:
 *  - SCIP_CONSADDED  : the conflict handler created a constraint out of the conflict set
 *  - SCIP_DIDNOTFIND : the conflict handler could not create a constraint out of the conflict set
 *  - SCIP_DIDNOTRUN  : the conflict handler was skipped
 */
#define DECL_CONFLICTEXEC(x) RETCODE x (SCIP* scip, CONFLICTHDLR* conflicthdlr, VAR** conflictvars, int nconflictvars, \
                                        Bool resolved, RESULT* result)



#include "def.h"
#include "retcode.h"
#include "memory.h"
#include "set.h"
#include "stat.h"
#include "var.h"



/*
 * Conflict Handler
 */

/** creates a conflict handler */
extern
RETCODE SCIPconflicthdlrCreate(
   CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SET*             set,                /**< global SCIP settings */
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
   VAR**            conflictvars,       /**< variables of the conflict set */
   int              nconflictvars,      /**< number of variables in the conflict set */
   Bool             resolved,           /**< is the conflict set already used to create a constraint? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** gets name of conflict handler */
extern
const char* SCIPconflicthdlrGetName(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets user data of conflict handler */
extern
CONFLICTHDLRDATA* SCIPconflicthdlrGetData(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** sets user data of conflict handler; user has to free old data in advance! */
extern
void SCIPconflicthdlrSetData(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata   /**< new conflict handler user data */
   );

/** gets priority of conflict handler */
extern
int SCIPconflicthdlrGetPriority(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** is conflict handler initialized? */
extern
Bool SCIPconflicthdlrIsInitialized(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
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

/** analyses conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set
 */
extern
RETCODE SCIPconflictAnalyse(
   CONFLICT*        conflict,           /**< conflict analysis data */
   const SET*       set,                /**< global SCIP settings */
   int              maxsize,            /**< maximal size of the conflict set or -1 for no restriction */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/** gets time in seconds used for analyzing propagation conflicts */
extern
Real SCIPconflictGetTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to propagation conflict analysis */
extern
Longint SCIPconflictGetNCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of valid conflicts detected in propagation conflict analysis */
extern
Longint SCIPconflictGetNConflicts(
   CONFLICT*        conflict            /**< conflict analysis data */
   );





/*
 * Infeasible LP Conflict Analysis
 */

/** creates conflict analysis data for infeasible LP conflicts */
extern
RETCODE SCIPlpconflictCreate(
   LPCONFLICT**     lpconflict,         /**< pointer to LP conflict analysis data */
   const SET*       set                 /**< global SCIP settings */
   );

/** frees conflict analysis data for infeasible LP conflicts */
extern
RETCODE SCIPlpconflictFree(
   LPCONFLICT**     lpconflict          /**< pointer to LP conflict analysis data */
   );

/** analyses conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set
 */
extern
RETCODE SCIPlpconflictAnalyse(
   LPCONFLICT*      lpconflict,         /**< LP conflict analysis data */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   int              maxsize,            /**< maximal size of the conflict set or -1 for no restriction */
   Bool*            success             /**< pointer to store whether a conflict constraint was created */
   );

/** gets time in seconds used for analyzing infeasible LP conflicts */
extern
Real SCIPlpconflictGetTime(
   LPCONFLICT*      lpconflict          /**< LP conflict analysis data */
   );

/** gets number of calls to infeasible LP conflict analysis */
extern
Longint SCIPlpconflictGetNCalls(
   LPCONFLICT*      lpconflict          /**< LP conflict analysis data */
   );

/** gets number of valid conflicts detected in infeasible LP conflict analysis */
extern
Longint SCIPlpconflictGetNConflicts(
   LPCONFLICT*      lpconflict          /**< LP conflict analysis data */
   );

#endif
