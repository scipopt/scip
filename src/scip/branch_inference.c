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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_inference.c,v 1.17 2005/11/02 11:14:43 bzfpfend Exp $"

/**@file   branch_inference.c
 * @brief  inference history branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_inference.h"


#define BRANCHRULE_NAME          "inference"
#define BRANCHRULE_DESC          "inference history branching"
#define BRANCHRULE_PRIORITY      1000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_CONFLICTWEIGHT  1000.0  /**< factor to weigh conflict score against inference score */
#define DEFAULT_CUTOFFWEIGHT       1.0  /**< factor to weigh average number of cutoffs in branching score */
#define DEFAULT_FRACTIONALS        TRUE /**< should branching on LP solution be restricted to the fractional variables? */


/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Real             conflictweight;     /**< factor to weigh conflict score against inference score */
   SCIP_Real             cutoffweight;       /**< factor to weigh average number of cutoffs in branching score */
   SCIP_Bool             fractionals;        /**< should branching on LP solution be restricted to the fractional variables? */
};



/** selects a variable out of the given candidate array and performs the branching */
static
SCIP_RETCODE performBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            cands,              /**< candidate array */
   int                   ncands,             /**< number of candidates */
   SCIP_Real             conflictweight,     /**< factor to weigh conflict score against inference score */
   SCIP_Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   )
{
   SCIP_VAR* var;
   SCIP_BRANCHDIR branchdir;
   SCIP_Real bestscore;
   SCIP_Real score;
   int bestcand;
   int c;
   SCIP_Real downscore;
   SCIP_Real upscore;

   /* search for variable with best score w.r.t. average inferences per branching */
   bestscore = 0.0;
   bestcand = -1;
   for( c = 0; c < ncands; ++c )
   {
      score = conflictweight * SCIPgetVarConflictScore(scip, cands[c])
         + SCIPgetVarAvgInferenceCutoffScore(scip, cands[c], cutoffweight);
      if( score > bestscore )
      {
         bestscore = score;
         bestcand = c;
      }
      SCIPdebugMessage(" -> cand <%s>: prio=%d, score=%g\n", SCIPvarGetName(cands[c]), SCIPvarGetBranchPriority(cands[c]),
         score);
   }
   assert(0 <= bestcand && bestcand < ncands);
   var = cands[bestcand];

   /* select the preferred branching direction: try to find a conflict as fast as possible */
   downscore = SCIPgetVarAvgInferences(scip, var, SCIP_BRANCHDIR_DOWNWARDS);
   upscore = SCIPgetVarAvgInferences(scip, var, SCIP_BRANCHDIR_UPWARDS);
   if( downscore > upscore + 1.0 )
      branchdir = SCIP_BRANCHDIR_DOWNWARDS;
   else if( downscore < upscore - 1.0 )
      branchdir = SCIP_BRANCHDIR_UPWARDS;
   else
      branchdir = SCIP_BRANCHDIR_AUTO;

   /* perform the branching */
   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (prio=%d, solval=%.12f, score=%g)\n",
      ncands, bestcand, SCIPvarGetName(var), SCIPvarGetBranchPriority(var), SCIPgetVarSol(scip, var), bestscore);
   SCIP_CALL( SCIPbranchVar(scip, var, branchdir) );

   return SCIP_OKAY;
}




/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeInference)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
#define branchInitInference NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitInference NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolInference NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolInference NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpInference)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** cands;
   int ncands;

   SCIPdebugMessage("Execlp method of inference branching\n");
   
   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( branchruledata->fractionals )
   {
      /* get LP candidates (fractional integer variables) */
      SCIP_CALL( SCIPgetLPBranchCands(scip, &cands, NULL, NULL, NULL, &ncands) );
   }
   else
   {
      /* get pseudo candidates (non-fixed integer variables) */
      SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands) );
   }

   /* perform the branching */
   SCIP_CALL( performBranching(scip, cands, ncands, branchruledata->conflictweight, branchruledata->cutoffweight) );

   *result = SCIP_BRANCHED;
   
   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsInference)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** cands;
   int ncands;

   SCIPdebugMessage("Execps method of inference branching\n");
   
   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get pseudo candidates (non-fixed integer variables) */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands) );

   /* perform the branching */
   SCIP_CALL( performBranching(scip, cands, ncands, branchruledata->conflictweight, branchruledata->cutoffweight) );

   *result = SCIP_BRANCHED;
   
   return SCIP_OKAY;
}




/*
 * branching specific interface methods
 */

/** creates the inference history braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleInference(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create inference branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreeInference, branchInitInference, branchExitInference, branchInitsolInference, branchExitsolInference, 
         branchExeclpInference, branchExecpsInference,
         branchruledata) );

   /* inference branching rule parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/inference/conflictweight", 
         "factor to weigh conflict score against inference score",
         &branchruledata->conflictweight, DEFAULT_CONFLICTWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/inference/cutoffweight", 
         "factor to weigh average number of cutoffs in branching score",
         &branchruledata->cutoffweight, DEFAULT_CUTOFFWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/inference/fractionals", 
         "should branching on LP solution be restricted to the fractional variables?",
         &branchruledata->fractionals, DEFAULT_FRACTIONALS, NULL, NULL) );

   return SCIP_OKAY;
}
