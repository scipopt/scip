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

/**@file   branch_inference.c
 * @brief  inference history branching rule
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Stefan Heinz
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

/** evaluate the given candidate with the given score against the currently best know candidate */
static
void evaluateCand(
   SCIP_VAR*            cand,                /**< candidate to be checked */
   SCIP_Real            score,               /**< score of the candidate */
   SCIP_VAR**           bestcand,            /**< pointer to the currently best candidate */
   SCIP_Real*           bestscore            /**< pointer to the score of the currently best candidate */
   )
{

   /* evaluate the candidate against the currently best candidate */
   if( (*bestscore) < score )
   {
      /* the score of the candidate is better than the currently best know candidate */
      (*bestscore) = score;
      (*bestcand) = cand;
   }
   else if( (*bestscore) == score ) /*lint !e777*/
   {
      SCIP_Real bestobj;
      SCIP_Real candobj;

      bestobj = REALABS(SCIPvarGetObj(*bestcand));
      candobj = REALABS(SCIPvarGetObj(cand));

      /* the candidate has the same score as the best known candidate; therefore we use a second and third
       * criteria to detect a unique best candidate;
       *
       * - the second criteria prefers the candidate with a larger absolute value of its objective coefficient
       *   since branching on that variable might trigger further propagation w.r.t. objective function
       * - if the absolute values of the objective coefficient are equal the variable index is used to define a
       *   unique best candidate
       *
       * @note It is very important to select a unique best candidate. Otherwise the solver might vary w.r.t. the
       *       performance to much since the candidate array which is used here (SCIPgetPseudoBranchCands() or
       *       SCIPgetLPBranchCands()) gets dynamically changed during the solution process. In particular,
       *       starting a probing mode might already change the order of these arrays. To be independent of that
       *       the selection should be unique. Otherwise, to selection process can get influenced by starting a
       *       probing or not.
       */
      if( bestobj < candobj || (bestobj == candobj && SCIPvarGetIndex(*bestcand) < SCIPvarGetIndex(cand)) ) /*lint !e777*/
      {
         (*bestcand) = cand;
      }
   }
}

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
   SCIP_VAR* bestcand;
   SCIP_VAR* cand;
   SCIP_Real bestscore;
   SCIP_Real score;
   int c;

   assert(ncands > 0);

   /* check if the weighted sum between the average inferences and conflict score should be used */
   if( useweightedsum )
   {
      bestcand = cands[0];
      bestscore = conflictweight * SCIPgetVarConflictScore(scip, cands[0])
         + inferenceweight * SCIPgetVarAvgInferenceCutoffScore(scip, cands[0], cutoffweight);

      for( c = 1; c < ncands; ++c )
      {
         cand = cands[c];
         assert(cand != NULL);

         /* compute weighted score for the candidate */
         score = conflictweight * SCIPgetVarConflictScore(scip, cand)
            + inferenceweight * SCIPgetVarAvgInferenceCutoffScore(scip, cand, cutoffweight);

         SCIPdebugMessage(" -> cand <%s>: prio=%d, solval=%g, score=%g\n", SCIPvarGetName(cand), SCIPvarGetBranchPriority(cand),
            SCIPgetVarSol(scip, cand), score);

         /* evaluate the candidate against the currently best candidate */
         evaluateCand(cand, score, &bestcand, &bestscore);
      }
   }
   else
   {
      bestcand = cands[0];
      bestscore = SCIPgetVarAvgInferenceScore(scip, cands[0]);

      /* search for variable with best score w.r.t. average inferences per branching */
      for( c = 1; c < ncands; ++c )
      {
         cand = cands[c];
         assert(cand != NULL);

         score = SCIPgetVarAvgInferenceScore(scip, cand);

         SCIPdebugMessage(" -> cand <%s>: prio=%d, solval=%g, score=%g\n", SCIPvarGetName(cand), SCIPvarGetBranchPriority(cand),
            SCIPgetVarSol(scip, cand), score);

         /* evaluate the candidate against the currently best candidate */
         evaluateCand(cand, score, &bestcand, &bestscore);
      }
   }

   assert(bestcand != NULL);

   SCIPdebugMessage(" -> %d candidates, selected variable <%s> (prio=%d, solval=%.12f, score=%g)\n",
      ncands, SCIPvarGetName(bestcand), SCIPvarGetBranchPriority(bestcand), SCIPgetVarSol(scip, bestcand), bestscore);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, bestcand, NULL, NULL, NULL) );

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
