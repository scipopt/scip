/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_leastinf.c,v 1.17 2005/02/14 13:35:39 bzfpfend Exp $"

/**@file   branch_leastinf.c
 * @brief  least infeasible LP branching rule
 * @author Tobias Achterberg
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
 * Callback methods
 */

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
DECL_BRANCHEXECLP(branchExeclpLeastinf)
{  /*lint --e{715}*/
   VAR** lpcands;
   Real* lpcandsfrac;
   int nlpcands;
   Real infeasibility;
   Real score;
   Real obj;
   Real bestscore;
   Real bestobj;
   int bestcand;
   int i;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Execlp method of leastinf branching\n");

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands) );
   assert(nlpcands > 0);

   /* search the least infeasible candidate */
   bestscore = REAL_MIN;
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

   debugMessage(" -> %d candidates, selected candidate %d: variable <%s> (frac=%g, obj=%g, factor=%g, score=%g)\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandsfrac[bestcand], bestobj,
      SCIPvarGetBranchFactor(lpcands[bestcand]), bestscore);

   /* perform the branching */
   CHECK_OKAY( SCIPbranchVar(scip, lpcands[bestcand]) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsLeastinf NULL




/*
 * branching specific interface methods
 */

/** creates the least infeasible LP braching rule and includes it in SCIP */
RETCODE SCIPincludeBranchruleLeastinf(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create inference branching rule data */
   branchruledata = NULL;

   /* include branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreeLeastinf, branchInitLeastinf, branchExitLeastinf, branchInitsolLeastinf, branchExitsolLeastinf, 
         branchExeclpLeastinf, branchExecpsLeastinf,
         branchruledata) );

   return SCIP_OKAY;
}
