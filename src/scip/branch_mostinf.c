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
#pragma ident "@(#) $Id: branch_mostinf.c,v 1.41 2010/11/23 19:47:58 bzfviger Exp $"

/**@file   branch_mostinf.c
 * @ingroup BRANCHINGRULES
 * @brief  most infeasible LP branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_mostinf.h"


#define BRANCHRULE_NAME          "mostinf"
#define BRANCHRULE_DESC          "most infeasible branching"
#define BRANCHRULE_PRIORITY      100
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
         if( SCIPisEQ(scip, SCIPvarGetLbLocal(SCIPvarGetMultaggrVars(cand)[i]), SCIPvarGetUbLocal(SCIPvarGetMultaggrVars(cand)[i])) )
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
   if( SCIPisInfinity(scip, candscore)
      || (!SCIPisInfinity(scip, *bestscore) && 
          (SCIPisGT(scip, candscore, *bestscore) || (SCIPisGE(scip, candscore, *bestscore) && obj > *bestobj))) )
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
SCIP_DECL_BRANCHCOPY(branchCopyMostinf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleMostinf(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeMostinf NULL


/** initialization method of branching rule (called after problem was transformed) */
#define branchInitMostinf NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitMostinf NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolMostinf NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolMostinf NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMostinf)
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

   SCIPdebugMessage("Execlp method of mostinf branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands) );
   assert(nlpcands > 0);

   /* search the most infeasible candidate */
   bestscore = SCIP_REAL_MIN;
   bestobj = 0.0;
   bestcand = -1;
   for( i = 0; i < nlpcands; ++i )
   {
      assert(lpcands[i] != NULL);

      infeasibility = lpcandsfrac[i];
      infeasibility = MIN(infeasibility, 1.0-infeasibility);
      score = infeasibility;
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

/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextMostinf)
{  /*lint --e{715}*/
   SCIP_VAR** externcands;
   SCIP_Real* externcandssol;
   SCIP_Real* externcandsscore;
   int nexterncands;
   SCIP_VAR* bestcand;
   SCIP_Real bestscore;
   SCIP_Real bestobj;
   SCIP_Real bestsol;
   SCIP_Real brpoint;
   int i;
   SCIP_NODE* downchild;
   SCIP_NODE* eqchild;
   SCIP_NODE* upchild;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execext method of mostinf branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &externcands, &externcandssol, &externcandsscore, NULL, &nexterncands, NULL, NULL, NULL) );
   assert(nexterncands > 0);

   /* search the most infeasible candidate */
   bestscore = SCIP_REAL_MIN;
   bestobj = 0.0;
   bestcand = NULL;
   bestsol = SCIP_INVALID;
   for( i = 0; i < nexterncands; ++i )
   {
      updateBestCandidate(scip, &bestcand, &bestscore, &bestobj, &bestsol, externcands[i], externcandsscore[i], externcandssol[i]);
   }
   
   if( bestcand == NULL )
   {
      SCIPerrorMessage("branchExecextMostinf failed to select a branching variable from %d candidates\n", nexterncands);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   brpoint = SCIPgetBranchingPoint(scip, bestcand, bestsol);

   SCIPdebugMessage(" -> %d candidates, selected variable <%s> (infeas=%g, obj=%g, factor=%g, score=%g), branching point=%g\n",
      nexterncands, SCIPvarGetName(bestcand), bestsol, bestobj,
      SCIPvarGetBranchFactor(bestcand), bestscore, brpoint);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVarVal(scip, bestcand, brpoint, &downchild, &eqchild, &upchild) );

   if( downchild != NULL || eqchild != NULL || upchild != NULL )
   {
      *result = SCIP_BRANCHED;
   }
   else
   {
      /* if there are no children, then variable should have been fixed by SCIPbranchVarVal */
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(bestcand), SCIPvarGetUbLocal(bestcand)));
      *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsMostinf NULL




/*
 * branching specific interface methods
 */

/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMostinf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchCopyMostinf,
         branchFreeMostinf, branchInitMostinf, branchExitMostinf, branchInitsolMostinf, branchExitsolMostinf, 
         branchExeclpMostinf, branchExecextMostinf, branchExecpsMostinf,
         NULL) );

   return SCIP_OKAY;
}
