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
#pragma ident "@(#) $Id: branch_fullstrong.c,v 1.14 2004/01/22 14:42:26 bzfpfend Exp $"

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
#define BRANCHRULE_PRIORITY      0



/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeFullstrong NULL


/** initialization method of branching rule (called when problem solving starts) */
#define branchInitFullstrong NULL


/** deinitialization method of branching rule (called when problem solving exits) */
#define branchExitFullstrong NULL


/** branching execution method for fractional LP solutions */
static
DECL_BRANCHEXECLP(branchExeclpFullstrong)
{  /*lint --e{715}*/
   VAR** lpcands;
   Real* lpcandssol;
   Real* lpcandsfrac;
   Real cutoffbound;
   Real lowerbound;
   Real bestdown;
   Real bestup;
   Real bestscore;
   Bool allcolsinlp;
   int nlpcands;
   int bestlpcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Execlp method of fullstrong branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get current lower objective bound of the local sub problem and global cutoff bound */
   lowerbound = SCIPgetLocalLowerbound(scip);
   cutoffbound = SCIPgetCutoffbound(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );
   assert(nlpcands > 0);

   if( nlpcands == 1 )
   {
      /* only one candidate: nothing has to be done */
      bestlpcand = 0;
      bestdown = lowerbound;
      bestup = lowerbound;
   }
   else
   {
      Real down;
      Real up;
      Real downgain;
      Real upgain;
      Real score;
      int c;

      /* search the full strong candidate */
      bestscore = -SCIPinfinity(scip);
      bestlpcand = -1;
      bestdown = 0.0;
      bestup = 0.0;
      for( c = 0; c < nlpcands; ++c )
      {
         assert(lpcands[c] != NULL);

         CHECK_OKAY( SCIPgetVarStrongbranch(scip, lpcands[c], INT_MAX, &down, &up) );
         down = MAX(down, lowerbound);
         up = MAX(up, lowerbound);
         downgain = down - lowerbound;
         upgain = up - lowerbound;

         if( allcolsinlp )
         {
            Bool downinf;
            Bool upinf;

            /* because all existing columns are in LP, the strong branching bounds are feasible lower bounds */
            downinf = SCIPisGE(scip, down, cutoffbound);
            upinf = SCIPisGE(scip, up, cutoffbound);

            if( downinf && upinf )
            {
               /* both roundings are infeasible -> node is infeasible */
               *result = SCIP_CUTOFF;
               debugMessage(" -> variable <%s> is infeasible in both directions\n", SCIPvarGetName(lpcands[c]));
               break; /* terminate initialization loop, because node is infeasible */
            }
            else if( downinf )
            {
               /* downwards rounding is infeasible -> change lower bound of variable to upward rounding */
               CHECK_OKAY( SCIPchgVarLb(scip, lpcands[c], SCIPceil(scip, lpcandssol[c])) );
               *result = SCIP_REDUCEDDOM;
               debugMessage(" -> variable <%s> is infeasible in downward branch\n", SCIPvarGetName(lpcands[c]));
               break; /* terminate initialization loop, because LP was changed */
            }
            else if( upinf )
            {
               /* upwards rounding is infeasible -> change upper bound of variable to downward rounding */
               CHECK_OKAY( SCIPchgVarUb(scip, lpcands[c], SCIPfloor(scip, lpcandssol[c])) );
               *result = SCIP_REDUCEDDOM;
               debugMessage(" -> variable <%s> is infeasible in upward branch\n", SCIPvarGetName(lpcands[c]));
               break; /* terminate initialization loop, because LP was changed */
            }
         }

         /* check for a better score */
         score = SCIPgetBranchScore(scip, downgain, upgain) + 1e-4; /* no gain -> use fractionalities */
         score *= SCIPvarGetBranchingPriority(lpcands[c]);
         if( score > bestscore )
         {
            bestlpcand = c;
            bestdown = down;
            bestup = up;
            bestscore = score;
         }

         /* update history values */
         CHECK_OKAY( SCIPupdateVarLPHistory(scip, lpcands[c], 0.0-lpcandsfrac[c], downgain, 1.0) );
         CHECK_OKAY( SCIPupdateVarLPHistory(scip, lpcands[c], 1.0-lpcandsfrac[c], upgain, 1.0) );

         debugMessage(" -> var <%s> (solval=%g, downgain=%g, upgain=%g, prio=%g, score=%g) -- best: <%s> (%g)\n",
            SCIPvarGetName(lpcands[c]), lpcandssol[c], downgain, upgain, SCIPvarGetBranchingPriority(lpcands[c]), score,
            SCIPvarGetName(lpcands[bestlpcand]), bestscore);
      }
   }

   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM )
   {
      NODE* node;

      assert(*result == SCIP_DIDNOTRUN);
      assert(bestlpcand >= 0);

      /* perform the branching */
      debugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g, down=%g, up=%g, prio=%g, score=%g)\n",
         nlpcands, bestlpcand, SCIPvarGetName(lpcands[bestlpcand]), lpcandssol[bestlpcand], bestdown, bestup, 
         SCIPvarGetBranchingPriority(lpcands[bestlpcand]), bestscore);

      /* create child node with x <= floor(x') */
      debugMessage(" -> creating child: <%s> <= %g\n",
         SCIPvarGetName(lpcands[bestlpcand]), SCIPfloor(scip, lpcandssol[bestlpcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      CHECK_OKAY( SCIPchgVarUbNode(scip, node, lpcands[bestlpcand], SCIPfloor(scip, lpcandssol[bestlpcand])) );
      if( allcolsinlp )
      {
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, bestdown) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));
      
      /* create child node with x >= ceil(x') */
      debugMessage(" -> creating child: <%s> >= %g\n", 
         SCIPvarGetName(lpcands[bestlpcand]), SCIPceil(scip, lpcandssol[bestlpcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      CHECK_OKAY( SCIPchgVarLbNode(scip, node, lpcands[bestlpcand], SCIPceil(scip, lpcandssol[bestlpcand])) );
      if( allcolsinlp )
      {
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, bestup) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsFullstrong NULL




/*
 * branching specific interface methods
 */

/** creates the full strong LP braching rule and includes it in SCIP */
RETCODE SCIPincludeBranchruleFullstrong(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create fullstrong branching rule data */
   branchruledata = NULL;

   /* include fullstrong branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
                  branchFreeFullstrong, branchInitFullstrong, branchExitFullstrong, 
                  branchExeclpFullstrong, branchExecpsFullstrong,
                  branchruledata) );

   return SCIP_OKAY;
}
