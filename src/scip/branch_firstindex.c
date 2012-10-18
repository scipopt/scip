/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_firstindex.c
 * @brief  pseudo costs branching rule
 * @author Ambros Gleixner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_firstindex.h"

#define BRANCHRULE_NAME                      "firstindex"
#define BRANCHRULE_DESC                      "branching on pseudo cost values"
#define BRANCHRULE_PRIORITY           -10000
#define BRANCHRULE_MAXDEPTH               -1
#define BRANCHRULE_MAXBOUNDDIST          1.0


/*
 * Local methods
 */

static
SCIP_RETCODE performBranchingFirstindex(
   SCIP*                 scip,
   SCIP_VAR**            cands,
   SCIP_Real*            candssol,
   int                   ncands,
   SCIP_RESULT*          result
   )
{
   SCIP_VAR* bestcand;
   SCIP_Real bestcandsol;
   int bestindex;
   int nchildren;
   int c;

   assert(scip != NULL);
   assert(cands != NULL);
   assert(ncands > 0);
   assert(result != NULL);

   /* select first candidate */
   bestcand = cands[0];
   bestcandsol = (candssol == NULL ? SCIPvarGetSol(bestcand, FALSE) : candssol[0]);
   bestindex = SCIPvarGetProbindex(bestcand);

   /* go through remaining candidates and find the one with smallest problem index */
   for( c = 1; c < ncands; c++ )
   {
      if( SCIPvarGetProbindex(cands[c]) < bestindex )
      {
         bestcand = cands[c];
         bestcandsol = (candssol == NULL ? SCIPvarGetSol(bestcand, FALSE) : candssol[c]);
         bestindex = SCIPvarGetProbindex(bestcand);
      }
   }
   assert(bestcand != NULL);

   /* perform the branching */
   SCIPdebugMessage(" -> %d cands, selected variable <%s> (type=%u, solval=%g, index=%d)\n",
      ncands, SCIPvarGetName(bestcand), SCIPvarGetType(bestcand), bestcandsol, bestindex);

   if( SCIPvarIsIntegral(bestcand) )
   {
      SCIP_CALL( SCIPbranchVar(scip, bestcand, NULL, NULL, NULL) );
      *result = SCIP_BRANCHED;
   }
   else
   {
      SCIP_CALL( SCIPbranchVarValNary(scip, bestcand, SCIPgetBranchingPoint(scip, bestcand, bestcandsol), 2, 0.0, 1.0, &nchildren) );

      if( nchildren > 1 )
         *result = SCIP_BRANCHED;
      else
      {
         /* if there are no children, then the variable should have been fixed by SCIPbranchVarValNary() */
         assert(SCIPisEQ(scip, SCIPvarGetLbLocal(bestcand), SCIPvarGetUbLocal(bestcand)));
         *result = SCIP_REDUCEDDOM;
      }
   }

   return SCIP_OKAY;
}



/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyFirstindex)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleFirstindex(scip) );

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpFirstindex)
{  /*lint --e{715}*/
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   int nlpcands;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of firstindex branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, NULL, &nlpcands) );
   assert(nlpcands > 0);

   /* perform the branching */
   SCIP_CALL( performBranchingFirstindex(scip, lpcands, lpcandssol, nlpcands, result) );

   return SCIP_OKAY;
}


/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextFirstindex)
{  /*lint --e{715}*/
   SCIP_VAR** externcands;
   SCIP_Real* externcandssol;
   int nprioexterncands;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execext method of firstindex branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &externcands, &externcandssol, NULL, NULL, &nprioexterncands, NULL, NULL, NULL) );
   assert(nprioexterncands > 0);

   /* perform the branching */
   SCIP_CALL( performBranchingFirstindex(scip, externcands, externcandssol, nprioexterncands, result) );

   return SCIP_OKAY;
}


/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsFirstindex)
{  /*lint --e{715}*/
   SCIP_VAR** pscands;
   int npriopscands;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of firstindex branching\n");

   /* get pseudo candidates (non-fixed integer variables) */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &pscands, &npriopscands, NULL) );
   assert(npriopscands > 0);

   /* perform the branching */
   SCIP_CALL( performBranchingFirstindex(scip, pscands, NULL, npriopscands, result) );

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the pseudo cost branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleFirstindex(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULE* branchrule;

   /* include allfullstrong branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, NULL) );

   assert(branchrule != NULL);

   /* set non-fundamental callbacks via specific setter functions*/
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyFirstindex) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpFirstindex) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextFirstindex) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsFirstindex) );

   return SCIP_OKAY;
}
