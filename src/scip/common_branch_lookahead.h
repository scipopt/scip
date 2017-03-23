/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   common_branch_Lookahead.c
 * @brief  Common functions used by the lookahead branching rules
 * @author Christoph Schubert
 */
#ifndef __SCIP_COMMON_BRANCH_LOOKAHEADABBREVIATED_H__
#define __SCIP_COMMON_BRANCH_LOOKAHEADABBREVIATED_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

struct BranchingDecision
{
   SCIP_VAR*             bestvar;
   SCIP_Real             bestval;
   SCIP_Real             bestdown;
   SCIP_Bool             bestdownvalid;
   SCIP_Real             bestup;
   SCIP_Bool             bestupvalid;
   SCIP_Real             provedbound;
};
typedef struct BranchingDecision BRANCHINGDECISION;

SCIP_RETCODE allocateBranchingDecision(
   SCIP*                 scip,
   BRANCHINGDECISION**   decision,
   SCIP_Real             lpobjval
);

SCIP_Bool isBranchingDecisionValid(
   SCIP*                 scip,
   BRANCHINGDECISION*    decision
);

void copyBranchingDecision(
   BRANCHINGDECISION*    sourcedecision,
   BRANCHINGDECISION*    targetdecision
);

void freeBranchingDecision(
   SCIP*                 scip,
   BRANCHINGDECISION**   decision
);


struct BranchRuleResult
{
   BRANCHINGDECISION*    decision;
   SCIP_Real*            candscores;
   SCIP_Real*            candslpvalue;
   SCIP_VAR**            candswithscore;
   int                   ncandscores;
};
typedef struct BranchRuleResult BRANCHRULERESULT;

SCIP_RETCODE allocateBranchRuleResultFull(
   SCIP*                 scip,
   BRANCHRULERESULT**    branchruleresult,
   SCIP_Real             lpobjval,
   int                   ncands
);

SCIP_RETCODE allocateBranchRuleResultReduced(
   SCIP*                 scip,
   BRANCHRULERESULT**    branchruleresult,
   SCIP_Real             lpobjval
);

void freeBranchRuleResultFull(
   SCIP*                 scip,
   BRANCHRULERESULT**    branchruleresult
);

void freeBranchRuleResultReduced(
   SCIP*                 scip,
   BRANCHRULERESULT**    branchruleresult
);

/**
 * Executes the branching on a given variable with a given value.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE branchOnVar(
   SCIP*                 scip,               /**< SCIP data structure */
   BRANCHINGDECISION*    decision
);

/** get a copy of the fractional candidates we can branch on
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE copyLPBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           lpcands,            /**< a pointer to store the variables */
   SCIP_Real**           lpcandssol,         /**< a pointer to store the solution values of the vars */
   SCIP_Real**           lpcandsfrac,
   int*                  nlpcands            /**< a pointer to store the number of candidates */
);

void freeLPBranchCands(
   SCIP*                 scip,
   SCIP_VAR***           lpcands,
   SCIP_Real**           lpcandssol,
   SCIP_Real**           lpcandsfrac
);

EXTERN
const char* getStatusString(
   SCIP_RESULT           result
);

#ifdef __cplusplus
}
#endif

#endif
