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
#pragma ident "@(#) $Id: branch_pscost.c,v 1.6 2004/11/29 12:17:14 bzfpfend Exp $"

/**@file   branch_pscost.c
 * @brief  pseudo costs branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_pscost.h"


#define BRANCHRULE_NAME          "pscost"
#define BRANCHRULE_DESC          "branching on pseudo cost values"
#define BRANCHRULE_PRIORITY      2000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0




/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreePscost NULL


/** initialization method of branching rule (called after problem was transformed) */
#define branchInitPscost NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitPscost NULL


/** branching execution method for fractional LP solutions */
static
DECL_BRANCHEXECLP(branchExeclpPscost)
{  /*lint --e{715}*/
   NODE* node;
   VAR** lpcands;
   Real* lpcandssol;
   Real* lpcandsfrac;
   Real bestscore;
   Real bestrootdiff;
   Real downprio;
   int nlpcands;
   int bestcand;
   int c;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Execlp method of pscost branching\n");

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, NULL, &nlpcands) );
   assert(nlpcands > 0);

   bestcand = -1;
   bestscore = -SCIPinfinity(scip);
   bestrootdiff = 0.0;
   for( c = 0; c < nlpcands; ++c )
   {
      Real score;
      Real rootsolval;
      Real rootdiff;

      score = SCIPgetVarPseudocostScore(scip, lpcands[c], lpcandssol[c]);
      rootsolval = SCIPvarGetRootSol(lpcands[c]);
      rootdiff = REALABS(lpcandssol[c] - rootsolval);
      if( SCIPisSumGT(scip, score, bestscore) || (SCIPisSumEQ(scip, score, bestscore) && rootdiff > bestrootdiff) )
      {
         bestcand = c;
         bestscore = score;
         bestrootdiff = rootdiff;
      }
   }
   assert(0 <= bestcand && bestcand < nlpcands);
   assert(!SCIPisFeasIntegral(scip, lpcandssol[bestcand]));

   /* perform the branching */
   debugMessage(" -> %d cands, selected cand %d: variable <%s> (solval=%g, score=%g)\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandssol[bestcand], bestscore);

   /* choose preferred branching direction */
   switch( SCIPvarGetBranchDirection(lpcands[bestcand]) )
   {
   case SCIP_BRANCHDIR_DOWNWARDS:
      downprio = 1.0;
      break;
   case SCIP_BRANCHDIR_UPWARDS:
      downprio = -1.0;
      break;
   case SCIP_BRANCHDIR_AUTO:
      downprio = SCIPvarGetRootSol(lpcands[bestcand]) - lpcandssol[bestcand];
      break;
   default:
      errorMessage("invalid preferred branching direction <%d> of variable <%s>\n", 
         SCIPvarGetBranchDirection(lpcands[bestcand]), SCIPvarGetName(lpcands[bestcand]));
      return SCIP_INVALIDDATA;
   }

   /* create child node with x <= floor(x') */
   debugMessage(" -> creating child: <%s> <= %g\n",
      SCIPvarGetName(lpcands[bestcand]), SCIPfeasFloor(scip, lpcandssol[bestcand]));
   CHECK_OKAY( SCIPcreateChild(scip, &node, downprio) );
   CHECK_OKAY( SCIPchgVarUbNode(scip, node, lpcands[bestcand], SCIPfeasFloor(scip, lpcandssol[bestcand])) );
      
   /* create child node with x >= ceil(x') */
   debugMessage(" -> creating child: <%s> >= %g\n", 
      SCIPvarGetName(lpcands[bestcand]), SCIPfeasCeil(scip, lpcandssol[bestcand]));
   CHECK_OKAY( SCIPcreateChild(scip, &node, -downprio) );
   CHECK_OKAY( SCIPchgVarLbNode(scip, node, lpcands[bestcand], SCIPfeasCeil(scip, lpcandssol[bestcand])) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsPscost NULL




/*
 * branching specific interface methods
 */

/** creates the pseudo cost braching rule and includes it in SCIP */
RETCODE SCIPincludeBranchrulePscost(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create pscost branching rule data */
   branchruledata = NULL;
   
   /* include branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreePscost, branchInitPscost, branchExitPscost, 
         branchExeclpPscost, branchExecpsPscost,
         branchruledata) );

   return SCIP_OKAY;
}
