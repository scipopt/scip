/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_leastinf.c,v 1.38 2010/09/26 11:31:35 bzfviger Exp $"

/**@file   branch_leastinf.c
 * @ingroup BRANCHINGRULES
 * @brief  least infeasible LP branching rule
 * @author Tobias Achterberg
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_leastinf.h"

#define BRANCHRULE_NAME          "leastinf"
#define BRANCHRULE_DESC          "least infeasible branching"
#define BRANCHRULE_PRIORITY      50
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

/*
 * Local methods
 */

/** compares the so far best branching candidate with a new candidate and updates best candidate, if new candidate is better */
static
void updateBestCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            bestvar,            /**< best branching candidate */
   SCIP_Real*            bestscore,          /**< score of best branching candidate */
   SCIP_Real*            bestobj,            /**< absolute objective value of best branching candidate */
   SCIP_Real*            bestsol,            /**< proposed branching point of best branching candidate */
   SCIP_VAR*             cand,               /**< branching candidate to consider */
   SCIP_Real             candscore,          /**< scoring of branching candidate */
   SCIP_Real             candsol             /**< proposed branching point of branching candidate */
   )
{
   SCIP_Real obj;
   
   assert(scip != NULL);
   assert(bestvar != NULL);
   assert(bestscore != NULL);
   assert(bestobj != NULL);
   assert(*bestobj >= 0.0);
   assert(cand != NULL);
   
   /* a branching variable candidate should either be an active problem variable or a multiaggregated variable */
   assert(SCIPvarIsActive(SCIPvarGetProbvar(cand)) ||
      SCIPvarGetStatus(SCIPvarGetProbvar(cand)) == SCIP_VARSTATUS_MULTAGGR);
   
   if( SCIPvarGetStatus(SCIPvarGetProbvar(cand)) == SCIP_VARSTATUS_MULTAGGR )
   {
      /* for a multiaggregated variable, we call updateBestCandidate function recursively with all variables in the multiaggregation */
      int i;
      
      cand = SCIPvarGetProbvar(cand);
      
      for( i = 0; i < SCIPvarGetMultaggrNVars(cand); ++i )
      {
         /* skip fixed variables */
         if( SCIPrelDiff(SCIPvarGetUbLocal(SCIPvarGetMultaggrVars(cand)[i]), SCIPvarGetLbLocal(SCIPvarGetMultaggrVars(cand)[i])) <= 2.0*SCIPepsilon(scip) )
            continue;
         
         updateBestCandidate(scip, bestvar, bestscore, bestobj, bestsol,
            SCIPvarGetMultaggrVars(cand)[i], candscore, SCIP_INVALID);
      }
      assert(*bestvar != NULL); /* if all variables were fixed, something is strange */
      
      return;
   }
   
   candscore *= SCIPvarGetBranchFactor(cand);
   obj = SCIPvarGetObj(cand);
   obj = REALABS(obj);
   if( SCIPisInfinity(scip, *bestscore)
      || (!SCIPisInfinity(scip, candscore) && 
          (SCIPisLT(scip, candscore, *bestscore) || (SCIPisLE(scip, candscore, *bestscore) && obj > *bestobj))) )
   {
      *bestvar = cand;
      *bestscore = candscore;
      *bestobj = obj;
      *bestsol = candsol;
   }
}

/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyLeastinf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleLeastinf(scip) );

   *valid = TRUE;
   
   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeLeastinf NULL


/** initialization method of branching rule (called after problem was transformed) */
#define branchInitLeastinf NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitLeastinf NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolLeastinf NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolLeastinf NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpLeastinf)
{  /*lint --e{715}*/
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   SCIP_Real infeasibility;
   SCIP_Real score;
   SCIP_Real obj;
   SCIP_Real bestscore;
   SCIP_Real bestobj;
   int bestcand;
   int i;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of leastinf branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands) );
   assert(nlpcands > 0);

   /* search the least infeasible candidate */
   bestscore = SCIP_REAL_MIN;
   bestobj = 0.0;
   bestcand = -1;
   for( i = 0; i < nlpcands; ++i )
   {
      assert(lpcands[i] != NULL);

      infeasibility = lpcandsfrac[i];
      infeasibility = MIN(infeasibility, 1.0-infeasibility);
      score = 1.0 - infeasibility;
      score *= SCIPvarGetBranchFactor(lpcands[i]);
      obj = SCIPvarGetObj(lpcands[i]);
      obj = REALABS(obj);
      if( SCIPisGT(scip, score, bestscore)
         || (SCIPisGE(scip, score, bestscore) && obj > bestobj) )
      {
         bestscore = score;
         bestobj = obj;
         bestcand = i;
      }
   }
   assert(bestcand >= 0);

   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (frac=%g, obj=%g, factor=%g, score=%g)\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandsfrac[bestcand], bestobj,
      SCIPvarGetBranchFactor(lpcands[bestcand]), bestscore);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECREL(branchExecrelLeastinf)
{  /*lint --e{715}*/
   SCIP_VAR** relaxcands;
   SCIP_Real* relaxcandssol;
   SCIP_Real* relaxcandsscore;
   int nrelaxcands;
   SCIP_VAR* bestcand;
   SCIP_Real bestscore;
   SCIP_Real bestobj;
   SCIP_Real bestsol;
   SCIP_Real brpoint;
   int i;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execrel method of leastinf branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetRelaxBranchCands(scip, &relaxcands, &relaxcandssol, &relaxcandsscore, NULL, &nrelaxcands, NULL, NULL, NULL) );
   assert(nrelaxcands > 0);

   /* search the least infeasible candidate */
   bestscore = SCIPinfinity(scip);
   bestobj  = 0.0;
   bestcand = NULL;
   bestsol = SCIP_INVALID;
   for( i = 0; i < nrelaxcands; ++i )
   {
      updateBestCandidate(scip, &bestcand, &bestscore, &bestobj, &bestsol, relaxcands[i], relaxcandsscore[i], relaxcandssol[i]);
   }

   if( bestcand == NULL )
   {
      SCIPerrorMessage("branchExecrelLeastinf failed to select a branching variable from %d candidates\n", nrelaxcands);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   brpoint = SCIPgetBranchingPoint(scip, bestcand, bestsol);

   SCIPdebugMessage(" -> %d candidates, selected variable <%s> (infeas=%g, obj=%g, factor=%g, score=%g), branching point=%g\n",
      nrelaxcands, SCIPvarGetName(bestcand), bestsol, bestobj,
      SCIPvarGetBranchFactor(bestcand), bestscore, brpoint);
   
   /* perform the branching */
   SCIP_CALL( SCIPbranchVarVal(scip, bestcand, brpoint, NULL, NULL, NULL) );
   
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsLeastinf NULL




/*
 * branching specific interface methods
 */

/** creates the least infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleLeastinf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchCopyLeastinf,
         branchFreeLeastinf, branchInitLeastinf, branchExitLeastinf, branchInitsolLeastinf, branchExitsolLeastinf, 
         branchExeclpLeastinf, branchExecrelLeastinf, branchExecpsLeastinf,
         NULL) );

   return SCIP_OKAY;
}
