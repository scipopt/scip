/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_branch.h,v 1.9 2005/01/21 09:17:10 bzfpfend Exp $"

/**@file   type_branch.h
 * @brief  type definitions for branching rules
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_BRANCH_H__
#define __TYPE_BRANCH_H__


typedef struct BranchCand BRANCHCAND;   /**< branching candidate storage */
typedef struct Branchrule BRANCHRULE;   /**< branching method data structure */
typedef struct BranchruleData BRANCHRULEDATA; /**< branching method specific data */


/** destructor of branching method to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 */
#define DECL_BRANCHFREE(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** initialization method of branching rule (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 */
#define DECL_BRANCHINIT(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** deinitialization method of branching rule (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 */
#define DECL_BRANCHEXIT(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** solving process initialization method of branching rule (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The branching rule may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 */
#define DECL_BRANCHINITSOL(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** solving process deinitialization method of branching rule (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The branching rule should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 */
#define DECL_BRANCHEXITSOL(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** branching execution method for fractional LP solutions
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 *  - allowaddcons    : is the branching rule allowed to add constraints to the current node in order to cut off the
 *                      current solution instead of creating a branching?
 *  - result          : pointer to store the result of the branching call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the current node was detected to be infeasible
 *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict clause) was generated; this result code must not be
 *                      returned, if allowaddcons is FALSE
 *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current LP solution infeasible
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_BRANCHED   : branching was applied
 *  - SCIP_DIDNOTRUN  : the branching rule was skipped
 */
#define DECL_BRANCHEXECLP(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule, Bool allowaddcons, RESULT* result)

/** branching execution method for not completely fixed pseudo solutions
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 *  - allowaddcons    : is the branching rule allowed to add constraints to the current node in order to cut off the
 *                      current solution instead of creating a branching?
 *  - result          : pointer to store the result of the branching call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the current node was detected to be infeasible
 *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict clause) was generated; this result code must not be
 *                      returned, if allowaddcons is FALSE
 *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current pseudo solution infeasible
 *  - SCIP_BRANCHED   : branching was applied
 *  - SCIP_DIDNOTRUN  : the branching rule was skipped
 */
#define DECL_BRANCHEXECPS(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule, Bool allowaddcons, RESULT* result)



#include "def.h"
#include "type_result.h"
#include "type_scip.h"



#endif
