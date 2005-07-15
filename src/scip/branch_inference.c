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
#pragma ident "@(#) $Id: branch_inference.c,v 1.12 2005/07/15 17:20:04 bzfpfend Exp $"

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

#define DEFAULT_CUTOFFWEIGHT     1.0    /**< factor to weigh average number of cutoffs in branching score */
#define DEFAULT_FRACTIONALS      TRUE   /**< should branching on LP solution be restricted to the fractional variables? */


/** branching rule data */
struct BranchruleData
{
   Real             cutoffweight;       /**< factor to weigh average number of cutoffs in branching score */
   Bool             fractionals;        /**< should branching on LP solution be restricted to the fractional variables? */
};



/** selects a variable out of the given candidate array and performs the branching */
static
RETCODE performBranching(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            cands,              /**< candidate array */
   int              ncands,             /**< number of candidates */
   Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   )
{
   VAR* var;
   BRANCHDIR branchdir;
   Real bestscore;
   Real score;
   int bestcand;
   int c;
   Real downinfs;
   Real upinfs;

   /* search for variable with best score w.r.t. average inferences per branching */
   bestscore = 0.0;
   bestcand = -1;
   for( c = 0; c < ncands; ++c )
   {
      score = SCIPgetVarAvgInferenceCutoffScore(scip, cands[c], cutoffweight);
      if( score > bestscore )
      {
         bestscore = score;
         bestcand = c;
      }
   }
   assert(0 <= bestcand && bestcand < ncands);
   var = cands[bestcand];

   /* select the preferred branching direction: try to find a conflict as fast as possible */
   downinfs = SCIPgetVarAvgInferences(scip, var, SCIP_BRANCHDIR_DOWNWARDS);
   upinfs = SCIPgetVarAvgInferences(scip, var, SCIP_BRANCHDIR_UPWARDS);
   if( downinfs > upinfs + 1.0 )
      branchdir = SCIP_BRANCHDIR_DOWNWARDS;
   else if( downinfs < upinfs - 1.0 )
      branchdir = SCIP_BRANCHDIR_UPWARDS;
   else
      branchdir = SCIP_BRANCHDIR_AUTO;

   /* perform the branching */
   debugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%.12f, score=%g)\n",
      ncands, bestcand, SCIPvarGetName(var), SCIPgetVarSol(scip, var), bestscore);
   CHECK_OKAY( SCIPbranchVar(scip, var, branchdir) );

   return SCIP_OKAY;
}




/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
DECL_BRANCHFREE(branchFreeInference)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;

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
DECL_BRANCHEXECLP(branchExeclpInference)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;
   VAR** cands;
   int ncands;

   debugMessage("Execlp method of inference branching\n");
   
   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( branchruledata->fractionals )
   {
      /* get LP candidates (fractional integer variables) */
      CHECK_OKAY( SCIPgetLPBranchCands(scip, &cands, NULL, NULL, NULL, &ncands) );
   }
   else
   {
      /* get pseudo candidates (non-fixed integer variables) */
      CHECK_OKAY( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands) );
   }

   /* perform the branching */
   CHECK_OKAY( performBranching(scip, cands, ncands, branchruledata->cutoffweight) );

   *result = SCIP_BRANCHED;
   
   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
static
DECL_BRANCHEXECPS(branchExecpsInference)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;
   VAR** cands;
   int ncands;

   debugMessage("Execps method of inference branching\n");
   
   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get pseudo candidates (non-fixed integer variables) */
   CHECK_OKAY( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands) );

   /* perform the branching */
   CHECK_OKAY( performBranching(scip, cands, ncands, branchruledata->cutoffweight) );

   *result = SCIP_BRANCHED;
   
   return SCIP_OKAY;
}




/*
 * branching specific interface methods
 */

/** creates the inference history braching rule and includes it in SCIP */
RETCODE SCIPincludeBranchruleInference(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create inference branching rule data */
   CHECK_OKAY( SCIPallocMemory(scip, &branchruledata) );
   
   /* include branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreeInference, branchInitInference, branchExitInference, branchInitsolInference, branchExitsolInference, 
         branchExeclpInference, branchExecpsInference,
         branchruledata) );

   /* inference branching rule parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
         "branching/inference/cutoffweight", 
         "factor to weigh average number of cutoffs in branching score",
         &branchruledata->cutoffweight, DEFAULT_CUTOFFWEIGHT, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddBoolParam(scip,
         "branching/inference/fractionals", 
         "should branching on LP solution be restricted to the fractional variables?",
         &branchruledata->fractionals, DEFAULT_FRACTIONALS, NULL, NULL) );

   return SCIP_OKAY;
}
