/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_mostinf.c
 * @brief  most infeasible LP branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_mostinf.h"


#define BRANCHRULE_NAME          "mostinf"
#define BRANCHRULE_DESC          "most infeasible branching"
#define BRANCHRULE_PRIORITY      100



/*
 * Callback methods
 */

static
DECL_BRANCHEXLP(branchExlpMostinf)
{
   VAR** lpcands;
   Real* lpcandsfrac;
   int nlpcands;
   Real infeasibility;
   Real maxinfeasibility;
   int bestcand;
   int i;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Exlp method of mostinf branching\n");

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandsfrac, &nlpcands) );
   assert(nlpcands > 0);

   /* search the most infeasible candidate */
   maxinfeasibility = 0.0;
   bestcand = -1;
   for( i = 0; i < nlpcands; ++i )
   {
      assert(lpcands[i] != NULL);

      infeasibility = lpcandsfrac[i];
      infeasibility = MIN(infeasibility, 1.0-infeasibility);
      if( infeasibility > maxinfeasibility )
      {
         maxinfeasibility = infeasibility;
         bestcand = i;
      }
   }
   assert(bestcand >= 0);

   debugMessage(" -> %d candidates, selected candidate %d: variable <%s> (frac=%g, infeasibility=%g)\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandsfrac[bestcand], maxinfeasibility);

   /* perform the branching */
   CHECK_OKAY( SCIPbranchVar(scip, lpcands[bestcand]) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}





/*
 * branching specific interface methods
 */

RETCODE SCIPincludeBranchruleMostinf(   /**< creates the most infeasible LP braching rule and includes it in SCIP */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
                  NULL, NULL, NULL, branchExlpMostinf, NULL,
                  NULL) );

   return SCIP_OKAY;
}
