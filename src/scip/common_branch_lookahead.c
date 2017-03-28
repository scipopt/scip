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
/*#define SCIP_DEBUG*/
/**@file   common_branch_Lookahead.c
 * @brief  LookaheadAbbreviated branching rule
 * @author Christoph Schubert
 */
/*
#define SCIP_DEBUG
 */
#define SCIP_STATISTIC

#include "scip/common_branch_lookahead.h"

SCIP_RETCODE branchOnVar(
   SCIP*                 scip                /**< SCIP data structure */,
   BRANCHINGDECISION*    decision
)
{
   SCIP_VAR* bestvar = decision->bestvar;
   SCIP_Real bestval = decision->bestval;

   SCIP_NODE* downchild = NULL;
   SCIP_NODE* upchild = NULL;

   assert(!SCIPisIntegral(scip, bestval));

   SCIPdebugMessage("Effective branching on var <%s> with value <%g>. Old domain: [%g..%g].\n",
      SCIPvarGetName(bestvar), bestval, SCIPvarGetLbLocal(bestvar), SCIPvarGetUbLocal(bestvar));

   SCIP_CALL( SCIPbranchVarVal(scip, bestvar, bestval, &downchild, NULL, &upchild) );

   assert(downchild != NULL);
   assert(upchild != NULL);

   /* update the lower bounds in the children; we must not do this if columns are missing in the LP
    * (e.g., because we are doing branch-and-price) or the problem should be solved exactly
    */
   if( SCIPallColsInLP(scip) && !SCIPisExactSolve(scip) )
   {
      SCIP_Real bestdown = decision->bestdown;
      SCIP_Bool bestdownvalid = decision->bestdownvalid;
      SCIP_Real bestup = decision->bestup;
      SCIP_Bool bestupvalid = decision->bestupvalid;
      SCIP_Real provedbound = decision->provedbound;

      SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdownvalid ? MAX(bestdown, provedbound) : provedbound) );
      SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestupvalid ? MAX(bestup, provedbound) : provedbound) );
   }
   SCIPdebugMessage(" -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
   SCIPdebugMessage(" -> up child's lowerbound: %g\n", SCIPnodeGetLowerbound(upchild));

   return SCIP_OKAY;
}

SCIP_RETCODE allocateBranchingDecision(
   SCIP*                 scip,
   BRANCHINGDECISION**   decision,
   SCIP_Real             lpobjval
   )
{
   SCIPallocBuffer(scip, decision);
   (*decision)->bestvar = NULL;
   (*decision)->bestdownvalid = FALSE;
   (*decision)->bestupvalid = FALSE;
   (*decision)->provedbound = lpobjval;

   return SCIP_OKAY;
}

void copyBranchingDecision(
   BRANCHINGDECISION*    sourcedecision,
   BRANCHINGDECISION*    targetdecision
   )
{
   targetdecision->bestvar = sourcedecision->bestvar;
   targetdecision->bestval = sourcedecision->bestval;
   targetdecision->bestdown = sourcedecision->bestdown;
   targetdecision->bestdownvalid = sourcedecision->bestdownvalid;
   targetdecision->bestup = sourcedecision->bestup;
   targetdecision->bestupvalid = sourcedecision->bestupvalid;
   targetdecision->provedbound = sourcedecision->provedbound;
}

SCIP_Bool isBranchingDecisionValid(
   SCIP*                 scip,
   BRANCHINGDECISION*    decision
   )
{
   return decision->bestvar != NULL;
}

void freeBranchingDecision(
   SCIP*                 scip,
   BRANCHINGDECISION**   decision
   )
{
   SCIPfreeBuffer(scip, decision);
}

SCIP_RETCODE allocateBranchRuleResultFull(
   SCIP*                 scip,
   BRANCHRULERESULT**    branchruleresult,
   SCIP_Real             lpobjval,
   int                   ncands
   )
{
   int i;

   SCIP_CALL( allocateBranchRuleResultReduced(scip, branchruleresult, lpobjval) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->lpinorms, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->lpistates, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->candscores, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->candswithscore, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->candslpvalue, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*branchruleresult)->candsvalfrac, ncands) );
   (*branchruleresult)->ncandscores = ncands;

   for( i = 0; i < ncands; i++ )
   {
      (*branchruleresult)->candscores[i] = -1;
      (*branchruleresult)->candswithscore[i] = NULL;
   }

   return SCIP_OKAY;
}

SCIP_RETCODE allocateBranchRuleResultReduced(
   SCIP*                 scip,
   BRANCHRULERESULT**    branchruleresult,
   SCIP_Real             lpobjval
)
{
   SCIP_CALL( SCIPallocBuffer(scip, branchruleresult) );
   SCIP_CALL( allocateBranchingDecision(scip, &(*branchruleresult)->decision, lpobjval) );
   (*branchruleresult)->candscores = NULL;
   (*branchruleresult)->candswithscore = NULL;
   (*branchruleresult)->candslpvalue = NULL;
   (*branchruleresult)->candsvalfrac = NULL;
   (*branchruleresult)->lpistates = NULL;
   (*branchruleresult)->lpinorms = NULL;
   (*branchruleresult)->ncandscores = -1;

   return SCIP_OKAY;
}

void freeBranchRuleResultFull(
   SCIP*                 scip,
   BRANCHRULERESULT**    branchruleresult
   )
{
   assert(&(*branchruleresult)->candslpvalue != NULL);
   assert(&(*branchruleresult)->candswithscore != NULL);
   assert(&(*branchruleresult)->candscores != NULL);

   SCIPfreeBufferArray(scip, &(*branchruleresult)->candsvalfrac);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->candslpvalue);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->candswithscore);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->candscores);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->lpistates);
   SCIPfreeBufferArray(scip, &(*branchruleresult)->lpinorms);
   freeBranchRuleResultReduced(scip, branchruleresult);
}

void freeBranchRuleResultReduced(
   SCIP*                 scip,
   BRANCHRULERESULT**    branchruleresult
   )
{
   freeBranchingDecision(scip, &(*branchruleresult)->decision);
   SCIPfreeBuffer(scip, branchruleresult);
}

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
   )
{
   SCIP_VAR** tmplpcands;
   SCIP_Real* tmplpcandssol;
   SCIP_Real* tmplpcandsfrac;

   /* get branching candidates and their solution values (integer variables with fractional value in the current LP) */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, &tmplpcandsfrac, nlpcands, NULL, NULL) );

   /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
    * solution during the second level branchings */
   SCIP_CALL( SCIPduplicateBufferArray(scip, lpcands, tmplpcands, *nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, lpcandssol, tmplpcandssol, *nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, lpcandsfrac, tmplpcandsfrac, *nlpcands) );

   return SCIP_OKAY;
}

void freeLPBranchCands(
   SCIP*                 scip,
   SCIP_VAR***           lpcands,
   SCIP_Real**           lpcandssol,
   SCIP_Real**           lpcandsfrac
   )
{
   SCIPfreeBuffer(scip, lpcandsfrac);
   SCIPfreeBuffer(scip, lpcandssol);
   SCIPfreeBuffer(scip, lpcands);
}

static const char* names[18] = { "", "SCIP_DIDNOTRUN", "SCIP_DELAYED", "SCIP_DIDNOTFIND", "SCIP_FEASIBLE", "SCIP_INFEASIBLE",
   "SCIP_UNBOUNDED", "SCIP_CUTOFF", "SCIP_SEPARATED", "SCIP_NEWROUND", "SCIP_REDUCEDDOM", "SCIP_CONSADDED",
   "SCIP_CONSCHANGED", "SCIP_BRANCHED", "SCIP_SOLVELP", "SCIP_FOUNDSOL", "SCIP_SUSPENDED", "SCIP_SUCCESS" };

const char* getStatusString(
   SCIP_RESULT           result
)
{
   return names[result];
}
