/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_inference.c,v 1.1 2004/04/15 10:41:21 bzfpfend Exp $"

/**@file   branch_inference.c
 * @brief  inference history branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_inference.h"


#define BRANCHRULE_NAME          "inference"
#define BRANCHRULE_DESC          "inference history branching"
#define BRANCHRULE_PRIORITY      1000
#define BRANCHRULE_MAXDEPTH      -1




/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeInference NULL


/** initialization method of branching rule (called when problem solving starts) */
#define branchInitInference NULL


/** deinitialization method of branching rule (called when problem solving exits) */
#define branchExitInference NULL


/** branching execution method for fractional LP solutions */
#define branchExeclpInference NULL


/** branching execution method for not completely fixed pseudo solutions */
static
DECL_BRANCHEXECPS(branchExecpsInference)
{
   NODE* node;
   VAR** pseudocands;
   Real bestscore;
   Real score;
   int npseudocands;
   int bestcand;
   int c;

   debugMessage("Execps method of relpscost branching\n");
   
   /* get pseudo candidates (non-fixed integer variables) */
   CHECK_OKAY( SCIPgetPseudoBranchCands(scip, &pseudocands, NULL, &npseudocands) );

   /* search for variable with best score w.r.t. average inferences per branching */
   bestscore = 0.0;
   bestcand = -1;
   for( c = 0; c < npseudocands; ++c )
   {
      score = SCIPgetVarAvgInferenceScore(scip, pseudocands[c]);
      if( score > bestscore )
      {
         bestscore = score;
         bestcand = c;
      }
   }
   assert(0 <= bestcand && bestcand < npseudocands);

   /* perform the branching */
   debugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%.12f, score=%g)\n",
      npseudocands, bestcand, SCIPvarGetName(pseudocands[bestcand]), 
      SCIPvarGetPseudoSol(pseudocands[bestcand]), bestscore);
   CHECK_OKAY( SCIPbranchVar(scip, pseudocands[bestcand]) );

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
   branchruledata = NULL;
   
   /* include branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                  branchFreeInference, branchInitInference, branchExitInference, 
                  branchExeclpInference, branchExecpsInference,
                  branchruledata) );

   return SCIP_OKAY;
}
