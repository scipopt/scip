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
#pragma ident "@(#) $Id: type_conflict.h,v 1.5 2004/04/27 15:50:06 bzfpfend Exp $"

/**@file   type_conflict.h
 * @brief  type definitions for conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_CONFLICT_H__
#define __TYPE_CONFLICT_H__


typedef struct Conflicthdlr CONFLICTHDLR; /**< conflict handler to process conflict sets */
typedef struct ConflicthdlrData CONFLICTHDLRDATA; /**< conflict handler data */
typedef struct Conflict CONFLICT;       /**< conflict analysis data structure */


/** destructor of conflict handler to free conflict handler data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conflicthdlr    : the conflict handler itself
 */
#define DECL_CONFLICTFREE(x) RETCODE x (SCIP* scip, CONFLICTHDLR* conflicthdlr)

/** initialization method of conflict handler (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conflicthdlr    : the conflict handler itself
 */
#define DECL_CONFLICTINIT(x) RETCODE x (SCIP* scip, CONFLICTHDLR* conflicthdlr)

/** deinitialization method of conflict handler (called before transformed problem is freed)
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
 *  - node            : node to add resulting conflict clause to (with SCIPaddConsNode())
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
#define DECL_CONFLICTEXEC(x) RETCODE x (SCIP* scip, CONFLICTHDLR* conflicthdlr, NODE* node, \
                                        VAR** conflictvars, int nconflictvars, Bool resolved, RESULT* result)



#include "def.h"
#include "type_retcode.h"
#include "type_result.h"
#include "type_var.h"



#endif
