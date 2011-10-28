/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_inference.c
 * @brief  inference history branching rule
 * @author Tobias Achterberg
 * @author Timo Berthold
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

#define DEFAULT_CONFLICTWEIGHT  1000.0  /**< weight in score calculations for conflict score */
#define DEFAULT_CUTOFFWEIGHT       1.0  /**< weight in score calculations for cutoff score */
#define DEFAULT_INFERENCEWEIGHT    1.0  /**< weight in score calculations for inference score */
#define DEFAULT_FRACTIONALS        TRUE /**< should branching on LP solution be restricted to the fractional variables? */
#define DEFAULT_USEWEIGHTEDSUM     TRUE /**< should a weighted sum of inference, conflict and cutoff weights be used? */



/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Real             conflictweight;     /**< weight in score calculations for conflict score */
   SCIP_Real             cutoffweight;       /**< weight in score calculations for cutoff score */
   SCIP_Real             inferenceweight;    /**< weight in score calculations for inference score */
   SCIP_Bool             fractionals;        /**< should branching on LP solution be restricted to the fractional variables? */
   SCIP_Bool             useweightedsum;     /**< should a weighted sum of inference, conflict and cutoff weights be used? */
};



/** selects a variable out of the given candidate array and performs the branching */
static
SCIP_RETCODE performBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            cands,              /**< candidate array */
   int                   ncands,             /**< number of candidates */
   SCIP_Real             conflictweight,     /**< weight in score calculations for conflict score */
   SCIP_Real             inferenceweight,    /**< weight in score calculations for inference score */
   SCIP_Real             cutoffweight,       /**< weight in score calculations for cutoff score */
   SCIP_Bool             useweightedsum      /**< should a weighted sum of inference, conflict and cutoff weights be used? */
   )
{
   SCIP_VAR* var;
   SCIP_Real bestscore;
   SCIP_Real score;
   int bestcand;
   int c;

   /* search for variable with best score w.r.t. average inferences per branching */
   bestscore = 0.0;
   bestcand = -1;
   if( useweightedsum )
   { 
      for( c = 0; c < ncands; ++c )
      {
         score = conflictweight * SCIPgetVarConflictScore(scip, cands[c])
            + inferenceweight * SCIPgetVarAvgInferenceCutoffScore(scip, cands[c], cutoffweight);
         if( score > bestscore )
         {
            bestscore = score;
            bestcand = c;
         }
         SCIPdebugMessage(" -> cand <%s>: prio=%d, solval=%g, score=%g\n", SCIPvarGetName(cands[c]), SCIPvarGetBranchPriority(cands[c]),
            SCIPgetVarSol(scip, cands[c]), score);
      }
   }
   else
   {
      for( c = 0; c < ncands; ++c )
      {
         score = SCIPgetVarAvgInferenceScore(scip, cands[c]);
         if( score > bestscore )
         {
            bestscore = score;
            bestcand = c;
         }
         SCIPdebugMessage(" -> cand <%s>: prio=%d, solval=%g, score=%g\n", SCIPvarGetName(cands[c]), SCIPvarGetBranchPriority(cands[c]),
            SCIPgetVarSol(scip, cands[c]), score);
      }
   }
   assert(0 <= bestcand && bestcand < ncands);
   var = cands[bestcand];

   /* perform the branching */
   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (prio=%d, solval=%.12f, score=%g)\n",
      ncands, bestcand, SCIPvarGetName(var), SCIPvarGetBranchPriority(var), SCIPgetVarSol(scip, var), bestscore);
   SCIP_CALL( SCIPbranchVar(scip, var, NULL, NULL, NULL) );

   return SCIP_OKAY;
}




/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyInference)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleInference(scip) );
   
   return SCIP_OKAY;
}

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
   SCIP_CALL( performBranching(scip, cands, ncands, branchruledata->conflictweight, 
         branchruledata->inferenceweight, branchruledata->cutoffweight, branchruledata->useweightedsum) );

   *result = SCIP_BRANCHED;
   
   return SCIP_OKAY;
}


/** branching execution method for relaxation solutions */
#define branchExecextInference NULL


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
   SCIP_CALL( performBranching(scip, cands, ncands, branchruledata->conflictweight, 
         branchruledata->inferenceweight, branchruledata->cutoffweight, branchruledata->useweightedsum) );

   *result = SCIP_BRANCHED;
   
   return SCIP_OKAY;
}




/*
 * branching specific interface methods
 */

/** creates the inference history branching rule and includes it in SCIP */
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
         branchCopyInference,
         branchFreeInference, branchInitInference, branchExitInference, branchInitsolInference, branchExitsolInference, 
         branchExeclpInference, branchExecextInference, branchExecpsInference,
         branchruledata) );

   /* inference branching rule parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/inference/conflictweight", 
         "weight in score calculations for conflict score",
         &branchruledata->conflictweight, TRUE, DEFAULT_CONFLICTWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/inference/inferenceweight", 
         "weight in score calculations for inference score",
         &branchruledata->inferenceweight, TRUE, DEFAULT_INFERENCEWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "branching/inference/cutoffweight", 
         "weight in score calculations for cutoff score",
         &branchruledata->cutoffweight, TRUE, DEFAULT_CUTOFFWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/inference/fractionals", 
         "should branching on LP solution be restricted to the fractional variables?",
         &branchruledata->fractionals, TRUE, DEFAULT_FRACTIONALS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/inference/useweightedsum", 
         "should a weighted sum of inference, conflict and cutoff weights be used?",
         &branchruledata->useweightedsum, FALSE, DEFAULT_USEWEIGHTEDSUM, NULL, NULL) );

   return SCIP_OKAY;
}
