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
#pragma ident "@(#) $Id: branch_conffullstrong.c,v 1.3 2004/03/30 12:51:40 bzfpfend Exp $"

/**@file   branch_conffullstrong.c
 * @brief  full strong LP branching rule, that creates infeasible children to give input to conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_conffullstrong.h"


#define BRANCHRULE_NAME          "conffullstrong"
#define BRANCHRULE_DESC          "full strong branching that creates infeasible children (input for conflict analysis)"
#define BRANCHRULE_PRIORITY      -100
#define BRANCHRULE_MAXDEPTH      -1




/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeConffullstrong NULL


/** initialization method of branching rule (called when problem solving starts) */
#define branchInitConffullstrong NULL


/** deinitialization method of branching rule (called when problem solving exits) */
#define branchExitConffullstrong NULL


/** branching execution method for fractional LP solutions */
static
DECL_BRANCHEXECLP(branchExeclpConffullstrong)
{  /*lint --e{715}*/
   NODE* node;
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

   debugMessage("Execlp method of conffullstrong branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get current lower objective bound of the local sub problem and global cutoff bound */
   lowerbound = SCIPgetLocalLowerbound(scip);
   cutoffbound = SCIPgetCutoffbound(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );
   assert(nlpcands > 0);

   /* if only one candidate exists, choose this one without applying strong branching */
   bestlpcand = 0;
   bestdown = lowerbound;
   bestup = lowerbound;
   bestscore = -SCIPinfinity(scip);
   if( nlpcands > 1 )
   {
      Real down;
      Real up;
      Real downgain;
      Real upgain;
      Real score;
      Bool lperror;
      int c;

      /* search the full strong candidate */
      for( c = 0; c < nlpcands; ++c )
      {
         assert(lpcands[c] != NULL);

         CHECK_OKAY( SCIPgetVarStrongbranch(scip, lpcands[c], INT_MAX, &down, &up, &lperror) );

         /* check for an error in strong branching */
         if( lperror )
         {
            SCIPmessage(scip, SCIP_VERBLEVEL_HIGH,
               "(node %lld) error in strong branching call for variable <%s> with solution %g\n", 
               SCIPgetNodenum(scip), SCIPvarGetName(lpcands[c]), lpcandssol[c]);
            break;
         }

         /* evaluate strong branching */
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

            if( downinf || upinf )
            {
               /* we found at least one infeasible child: generate this branching immediately */
               bestlpcand = c;
               bestdown = lowerbound; /* make sure, the conflict analysis gets the child nodes soon */
               bestup = lowerbound;
               bestscore = SCIPinfinity(scip);
               break;
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

   assert(*result == SCIP_DIDNOTRUN);
   assert(0 <= bestlpcand && bestlpcand < nlpcands);
   
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

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsConffullstrong NULL




/*
 * branching specific interface methods
 */

/** creates the full strong LP braching rule and includes it in SCIP */
RETCODE SCIPincludeBranchruleConffullstrong(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create conffullstrong branching rule data */
   branchruledata = NULL;

   /* include conffullstrong branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                  branchFreeConffullstrong, branchInitConffullstrong, branchExitConffullstrong, 
                  branchExeclpConffullstrong, branchExecpsConffullstrong,
                  branchruledata) );

   return SCIP_OKAY;
}
