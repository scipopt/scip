/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_branch.h,v 1.3 2004/02/05 14:12:44 bzfpfend Exp $"

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

/** initialization method of branching rule (called when problem solving starts)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 */
#define DECL_BRANCHINIT(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** deinitialization method of branching rule (called when problem solving exits)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 */
#define DECL_BRANCHEXIT(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** branching execution method for fractional LP solutions
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 *  - result          : pointer to store the result of the branching call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : the current node was detected to be infeasible
 *  - SCIP_BRANCHED   : branching was applied
 *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current LP solution infeasible
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_DIDNOTRUN  : the branching rule was skipped
 */
#define DECL_BRANCHEXECLP(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule, RESULT* result)

/** branching execution method for not completely fixed pseudo solutions
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - branchrule      : the branching rule itself
 *  - result          : pointer to store the result of the branching call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : the current node was detected to be infeasible
 *  - SCIP_BRANCHED   : branching was applied
 *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current pseudo solution infeasible
 *  - SCIP_DIDNOTRUN  : the branching rule was skipped
 */
#define DECL_BRANCHEXECPS(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule, RESULT* result)



#include "def.h"
#include "type_result.h"
#include "type_scip.h"



#endif
