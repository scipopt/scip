/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_fullstrong.c,v 1.11 2004/01/13 11:58:29 bzfpfend Exp $"

/**@file   branch_fullstrong.c
 * @brief  full strong LP branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_fullstrong.h"


#define BRANCHRULE_NAME          "fullstrong"
#define BRANCHRULE_DESC          "full strong branching"
#define BRANCHRULE_PRIORITY      1000



/*
 * Callback methods
 */

static
DECL_BRANCHEXECLP(branchExeclpFullstrong)
{  /*lint --e{715}*/
   VAR** lpcands;
   Real* lpcandssol;
   Real upperbound;
   Real lowerbound;
   Real down;
   Real up;
   Real score;
   Real bestdown;
   Real bestup;
   Real bestscore;
   Bool allcolsinlp;
   int nlpcands;
   int bestcand;
   int i;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Execlp method of fullstrong branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get current lower objective bound of the local sub problem and global upper bound */
   lowerbound = SCIPgetLocalLowerbound(scip);
   upperbound = SCIPgetUpperbound(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands) );
   assert(nlpcands > 0);

   /* search the full strong candidate */
   bestscore = -SCIPinfinity(scip);
   bestcand = -1;
   bestdown = 0.0;
   bestup = 0.0;
   for( i = 0; i < nlpcands && *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM; ++i )
   {
      assert(lpcands[i] != NULL);

      CHECK_OKAY( SCIPgetVarStrongbranch(scip, lpcands[i], INT_MAX, &down, &up) );
      down = MAX(down, lowerbound);
      up = MAX(up, lowerbound);

      if( allcolsinlp )
      {
         Bool downinf;
         Bool upinf;

         /* because all existing columns are in LP, the strong branching bounds are feasible lower bounds */
         downinf = SCIPisGE(scip, down, upperbound);
         upinf = SCIPisGE(scip, up, upperbound);

         if( downinf && upinf )
         {
            /* both roundings are infeasible -> node is infeasible */
            *result = SCIP_CUTOFF;
            debugMessage(" -> variable <%s> is infeasible in both directions\n", SCIPvarGetName(lpcands[i]));
         }
         else if( downinf )
         {
            /* downwards rounding is infeasible -> change lower bound of variable to upward rounding */
            CHECK_OKAY( SCIPchgVarLb(scip, lpcands[i], SCIPceil(scip, lpcandssol[i])) );
            *result = SCIP_REDUCEDDOM;
            debugMessage(" -> variable <%s> is infeasible in downward branch\n", SCIPvarGetName(lpcands[i]));
         }
         else if( upinf )
         {
            /* upwards rounding is infeasible -> change upper bound of variable to downward rounding */
            CHECK_OKAY( SCIPchgVarUb(scip, lpcands[i], SCIPfloor(scip, lpcandssol[i])) );
            *result = SCIP_REDUCEDDOM;
            debugMessage(" -> variable <%s> is infeasible in upward branch\n", SCIPvarGetName(lpcands[i]));
         }
      }

      score = SCIPgetBranchScore(scip, down - lowerbound, up - lowerbound) + 1e-6; /* no gain -> use fractionalities */
      score *= SCIPvarGetBranchingPriority(lpcands[i]);
      debugMessage("   -> var <%s> (solval=%g, down=%g, up=%g, prio=%g, score=%g)\n",
         SCIPvarGetName(lpcands[i]), lpcandssol[i], down, up, SCIPvarGetBranchingPriority(lpcands[i]), score);

      if( score > bestscore )
      {
         bestcand = i;
         bestdown = down;
         bestup = up;
         bestscore = score;
      }
   }
   assert(bestcand >= 0);

   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM )
   {
      NODE* node;

      assert(*result == SCIP_DIDNOTRUN);

      /* perform the branching */
      debugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g, down=%g, up=%g, prio=%g, score=%g)\n",
         nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandssol[bestcand], bestdown, bestup, 
         SCIPvarGetBranchingPriority(lpcands[bestcand]), bestscore);

      /* create child node with x <= floor(x') */
      debugMessage(" -> creating child: <%s> <= %g\n",
         SCIPvarGetName(lpcands[bestcand]), SCIPfloor(scip, lpcandssol[bestcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      CHECK_OKAY( SCIPchgVarUbNode(scip, node, lpcands[bestcand], SCIPfloor(scip, lpcandssol[bestcand])) );
      if( allcolsinlp )
      {
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, bestdown) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));
      
      /* create child node with x >= ceil(x') */
      debugMessage(" -> creating child: <%s> >= %g\n", 
         SCIPvarGetName(lpcands[bestcand]), SCIPceil(scip, lpcandssol[bestcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      CHECK_OKAY( SCIPchgVarLbNode(scip, node, lpcands[bestcand], SCIPceil(scip, lpcandssol[bestcand])) );
      if( allcolsinlp )
      {
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, bestup) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}





/*
 * branching specific interface methods
 */

/** creates the full strong LP braching rule and includes it in SCIP */
RETCODE SCIPincludeBranchruleFullstrong(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
                  NULL, NULL, NULL, branchExeclpFullstrong, NULL,
                  NULL) );

   return SCIP_OKAY;
}
