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
#pragma ident "@(#) $Id: branch_history.c,v 1.4 2004/01/19 14:10:03 bzfpfend Exp $"

/**@file   branch_history.c
 * @brief  history branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_history.h"


#define BRANCHRULE_NAME          "history"
#define BRANCHRULE_DESC          "history branching"
#define BRANCHRULE_PRIORITY      10000

#define DEFAULT_RELIABLE         4.0    /**< minimum history size to regard history value as reliable */
#define DEFAULT_MAXLOOKAHEAD     8      /**< maximal number of further variables evaluated without better score */
#define DEFAULT_STRONGBRANCHCAND 100    /**< maximal number of candidates initialized with strong branching per node */
#define DEFAULT_STRONGBRANCHITER 0      /**< iteration limit for strong branching init of history entries (0: auto) */


/** branching rule data */
struct BranchruleData
{
   Real             reliable;           /**< minimum history size to regard history value as reliable */
   int              maxlookahead;       /**< maximal number of further variables evaluated without better score */
   int              strongbranchcand;   /**< maximal number of candidates initialized with strong branching per node */
   int              strongbranchiter;   /**< iteration limit for strong branching init of history entries (0: auto) */
};



/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
DECL_BRANCHFREE(branchFreeHistory)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called when problem solving starts) */
#define branchInitHistory NULL


/** deinitialization method of branching rule (called when problem solving exits) */
#define branchExitHistory NULL


/** branching execution method for fractional LP solutions */
static
DECL_BRANCHEXECLP(branchExeclpHistory)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;
   VAR** lpcands;
   Real* lpcandssol;
   Real* lpcandsfrac;
   Real lowerbound;
   Real cutoffbound;
   Real bestdown;
   Real bestup;
   Real bestscore;
   Bool bestisstrongbranch;
   Bool allcolsinlp;
   int nlpcands;
   int bestlpcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Execlp method of history branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get current lower objective bound of the local sub problem and global cutoff bound */
   lowerbound = SCIPgetLocalLowerbound(scip);
   cutoffbound = SCIPgetCutoffbound(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );
   assert(nlpcands > 0);

   bestlpcand = -1;
   bestisstrongbranch = FALSE;
   bestdown = SCIP_INVALID;
   bestup = SCIP_INVALID;

   if( nlpcands == 1 )
   {
      /* only one candidate: nothing has to be done */
      bestlpcand = 0;
   }
   else
   {
      int* sbcands;
      Real* sbcandscores;
      int nsbcands;
      int maxnsbcands;
      Real downsize;
      Real upsize;
      Real size;
      Real score;
      Real down;
      Real up;
      Real downgain;
      Real upgain;
      Longint actnodenum;
      int maxlookahead;
      int lookahead;
      int sbiter;
      int i;
      int j;
      int c;

      /* get buffer for storing the unreliable candidates */
      maxnsbcands = MIN(nlpcands, branchruledata->strongbranchcand);
      CHECK_OKAY( SCIPallocBufferArray(scip, &sbcands, maxnsbcands+1) ); /* allocate one additional slot for convenience */
      CHECK_OKAY( SCIPallocBufferArray(scip, &sbcandscores, maxnsbcands+1) );
      nsbcands = 0;

      /* get current node number */
      actnodenum = SCIPgetNodenum(scip);

      /* search for the best history candidate, while remembering unreliable candidates in a sorted buffer */
      bestscore = -SCIPinfinity(scip);
      for( c = 0; c < nlpcands; ++c )
      {
         assert(lpcands[c] != NULL);
         assert(!SCIPisIntegral(scip, lpcandssol[c]));

         /* if strong branching was already applied to the variable at the current node, use the old strong branching
          * values and don't use strong branching again
          */
         if( SCIPgetVarStrongbranchNode(scip, lpcands[c]) == actnodenum )
         {
            Real solval;
            Real frac;

            /* get the old strong branching values */
            CHECK_OKAY( SCIPgetVarStrongbranchLast(scip, lpcands[c], &down, &up, &solval) );
            frac = solval - SCIPfloor(scip, solval);
            assert(SCIPisPositive(scip, frac));
            downgain = MAX(down - lowerbound, 0.0) * lpcandsfrac[c]/frac;
            upgain = MAX(up - lowerbound, 0.0) * (1.0-lpcandsfrac[c])/(1.0-frac);

            /* calculate score out of old strong branching values */
            score = SCIPgetBranchScore(scip, downgain, upgain);
            score *= SCIPvarGetBranchingPriority(lpcands[c]);
         }
         else
         {
            /* get history score of candidate */
            score = SCIPgetVarLPHistoryScore(scip, lpcands[c], lpcandssol[c]);
            score *= SCIPvarGetBranchingPriority(lpcands[c]);

            /* check, if the history score of the variable is reliable */
            downsize = SCIPgetVarLPHistoryCount(scip, lpcands[c], 0);
            upsize = SCIPgetVarLPHistoryCount(scip, lpcands[c], 1);
            size = MIN(downsize, upsize);
            if( size < branchruledata->reliable )
            {
               /* history of variable is not reliable: insert candidate in sbcands buffer */
               for( j = nsbcands; j > 0 && score > sbcandscores[j-1]; --j )
               {
                  sbcands[j] = sbcands[j-1];
                  sbcandscores[j] = sbcandscores[j-1];
               }
               sbcands[j] = c;
               sbcandscores[j] = score;
               nsbcands++;
               nsbcands = MIN(nsbcands, maxnsbcands);
            }
         }

         /* check for better score of candidate */
         if( score > bestscore )
         {
            bestscore = score;
            bestlpcand = c;
         }
      }
      assert(bestlpcand >= 0);

      /* initialize unreliable candidates with strong branching until maxlookahead is reached */
      maxlookahead = branchruledata->maxlookahead;
      sbiter = branchruledata->strongbranchiter;
      if( sbiter == 0 )
      {
         int nlps = SCIPgetNLPs(scip);
         sbiter = 2*(SCIPgetNLPIterations(scip)) / MAX(1, nlps);
         sbiter = MAX(sbiter, 10);
      }
      lookahead = 0;
      for( i = 0; i < nsbcands && lookahead < maxlookahead; ++i )
      {
         /* get candidate number to initialize */
         c = sbcands[i];
         assert(!SCIPisIntegral(scip, lpcandssol[c]));

         debugMessage("initializing history (%g/%g) of variable <%s> at %g with strong branching (%d iterations)\n",
            SCIPgetVarLPHistoryCount(scip, lpcands[c], 0), SCIPgetVarLPHistoryCount(scip, lpcands[c], 1), 
            SCIPvarGetName(lpcands[c]), lpcandssol[c], sbiter);

         /* use strong branching on candidate */
         CHECK_OKAY( SCIPgetVarStrongbranch(scip, lpcands[c], sbiter, &down, &up) );
         down = MAX(down, lowerbound);
         up = MAX(up, lowerbound);
         downgain = down - lowerbound;
         upgain = up - lowerbound;

         /* check for possible fixings */
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
            lookahead = 0;
            bestisstrongbranch = TRUE;
         }
         else if( SCIPisLT(scip, score, bestscore) )
            lookahead++;
      
         /* update history values */
         CHECK_OKAY( SCIPupdateVarLPHistory(scip, lpcands[c], 0.0-lpcandsfrac[c], downgain, 1.0) );
         CHECK_OKAY( SCIPupdateVarLPHistory(scip, lpcands[c], 1.0-lpcandsfrac[c], upgain, 1.0) );
      
         debugMessage(" -> var <%s> (solval=%g, downgain=%g, upgain=%g, prio=%g, score=%g) -- best: <%s> (%g), lookahead=%d/%d\n",
            SCIPvarGetName(lpcands[c]), lpcandssol[c], downgain, upgain, SCIPvarGetBranchingPriority(lpcands[c]), score,
            SCIPvarGetName(lpcands[bestlpcand]), bestscore, lookahead, maxlookahead);
      }

      /* free buffer for the unreliable candidates */
      CHECK_OKAY( SCIPfreeBufferArray(scip, &sbcandscores) );
      CHECK_OKAY( SCIPfreeBufferArray(scip, &sbcands) );
   }

   /* if no domain could be reduced, create the branching */
   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM )
   {
      NODE* node;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestlpcand && bestlpcand < nlpcands);
      assert(!SCIPisIntegral(scip, lpcandssol[bestlpcand]));

      /* perform the branching */
      debugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%.12f, down=%g, up=%g, prio=%g, score=%g)\n",
         nlpcands, bestlpcand, SCIPvarGetName(lpcands[bestlpcand]), lpcandssol[bestlpcand], bestdown, bestup, 
         SCIPvarGetBranchingPriority(lpcands[bestlpcand]), bestscore);

      /* create child node with x <= floor(x') */
      debugMessage(" -> creating child: <%s> <= %g\n",
         SCIPvarGetName(lpcands[bestlpcand]), SCIPfloor(scip, lpcandssol[bestlpcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      CHECK_OKAY( SCIPchgVarUbNode(scip, node, lpcands[bestlpcand], SCIPfloor(scip, lpcandssol[bestlpcand])) );
      if( allcolsinlp && bestisstrongbranch )
      {
         assert(SCIPisLT(scip, bestdown, cutoffbound));
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, bestdown) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));
      
      /* create child node with x >= ceil(x') */
      debugMessage(" -> creating child: <%s> >= %g\n", 
         SCIPvarGetName(lpcands[bestlpcand]), SCIPceil(scip, lpcandssol[bestlpcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      CHECK_OKAY( SCIPchgVarLbNode(scip, node, lpcands[bestlpcand], SCIPceil(scip, lpcandssol[bestlpcand])) );
      if( allcolsinlp && bestisstrongbranch )
      {
         assert(SCIPisLT(scip, bestup, cutoffbound));
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, bestup) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsHistory NULL




/*
 * branching specific interface methods
 */

/** creates the history braching rule and includes it in SCIP */
RETCODE SCIPincludeBranchruleHistory(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create history branching rule data */
   CHECK_OKAY( SCIPallocMemory(scip, &branchruledata) );
   
   /* include branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
                  branchFreeHistory, branchInitHistory, branchExitHistory, branchExeclpHistory, branchExecpsHistory,
                  branchruledata) );

   /* history branching rule parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "branching/history/reliable", 
                  "minimum history size to regard history value as reliable",
                  &branchruledata->reliable, DEFAULT_RELIABLE, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "branching/history/maxlookahead", 
                  "maximal number of further variables evaluated without better score",
                  &branchruledata->maxlookahead, DEFAULT_MAXLOOKAHEAD, 1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "branching/history/strongbranchcand", 
                  "maximal number of candidates initialized with strong branching per node",
                  &branchruledata->strongbranchcand, DEFAULT_STRONGBRANCHCAND, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "branching/history/strongbranchiter", 
                  "iteration limit for strong branching initializations of history entries (0: auto)",
                  &branchruledata->strongbranchiter, DEFAULT_STRONGBRANCHITER, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
