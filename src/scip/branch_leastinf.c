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

/**@file   branch_leastinf.c
 * @brief  least infeasible LP branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_leastinf.h"


#define BRANCHRULE_NAME          "leastinf"
#define BRANCHRULE_DESC          "least infeasible branching"
#define BRANCHRULE_PRIORITY      50



/*
 * Callback methods
 */

static
DECL_BRANCHEXLP(branchExlpLeastinf)
{
   VAR** lpcands;
   Real* lpcandsfrac;
   int nlpcands;
   Real infeasibility;
   Real mininfeasibility;
   int bestcand;
   int i;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Exlp method of leastinf branching\n");

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, &nlpcands) );
   assert(nlpcands > 0);

   /* search the least infeasible candidate */
   mininfeasibility = 1.0;
   bestcand = -1;
   for( i = 0; i < nlpcands; ++i )
   {
      assert(lpcands[i] != NULL);

      infeasibility = lpcandsfrac[i];
      infeasibility = MIN(infeasibility, 1.0-infeasibility);
      if( infeasibility < mininfeasibility )
      {
         mininfeasibility = infeasibility;
         bestcand = i;
      }
   }
   assert(bestcand >= 0);

   debugMessage(" -> %d candidates, selected candidate %d: variable <%s> (frac=%g, infeasibility=%g)\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandsfrac[bestcand], mininfeasibility);

   /* perform the branching */
   CHECK_OKAY( SCIPbranchVar(scip, lpcands[bestcand]) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}





/*
 * branching specific interface methods
 */

RETCODE SCIPincludeBranchruleLeastinf(  /**< creates the least infeasible LP braching rule and includes it in SCIP */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
                  NULL, NULL, NULL, branchExlpLeastinf, NULL,
                  NULL) );

   return SCIP_OKAY;
}
