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
#pragma ident "@(#) $Id: branch_relpscost.c,v 1.8 2004/08/03 16:02:50 bzfpfend Exp $"

/**@file   branch_relpscost.c
 * @brief  reliable pseudo costs branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_relpscost.h"


#define BRANCHRULE_NAME          "relpscost"
#define BRANCHRULE_DESC          "reliability branching on pseudo cost values"
#define BRANCHRULE_PRIORITY      10000
#define BRANCHRULE_MAXDEPTH      -1

#define DEFAULT_MINRELIABLE      8.0    /**< minimal value for minimum pseudo cost size to regard pseudo cost value as reliable */
#define DEFAULT_MAXRELIABLE      8.0    /**< maximal value for minimum pseudo cost size to regard pseudo cost value as reliable */
#define DEFAULT_SBITERQUOT       0.5    /**< maximal fraction of strong branching LP iterations compared to normal iters */
#define DEFAULT_SBITEROFS   100000      /**< additional number of allowed strong branching LP iterations */
#define DEFAULT_MAXLOOKAHEAD     4      /**< maximal number of further variables evaluated without better score */
#define DEFAULT_INITCAND       100      /**< maximal number of candidates initialized with strong branching per node */
#define DEFAULT_INITITER         0      /**< iteration limit for strong branching init of pseudo cost entries (0: auto) */


/** branching rule data */
struct BranchruleData
{
   Real             minreliable;        /**< minimal value for minimum pseudo cost size to regard pseudo cost value as reliable */
   Real             maxreliable;        /**< maximal value for minimum pseudo cost size to regard pseudo cost value as reliable */
   Real             sbiterquot;         /**< maximal fraction of strong branching LP iterations compared to normal iters */
   int              sbiterofs;          /**< additional number of allowed strong branching LP iterations */
   int              maxlookahead;       /**< maximal number of further variables evaluated without better score */
   int              initcand;           /**< maximal number of candidates initialized with strong branching per node */
   int              inititer;           /**< iteration limit for strong branching init of pseudo cost entries (0: auto) */
};



/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
DECL_BRANCHFREE(branchFreeRelpscost)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
#define branchInitRelpscost NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitRelpscost NULL


#define MINMAXDEPTH   20
#define MAXMAXDEPTH  100
#define MAXSIZE     5000
#define FRACSCORE   1e-3
/** branching execution method for fractional LP solutions */
static
DECL_BRANCHEXECLP(branchExeclpRelpscost)
{  /*lint --e{715}*/
   BRANCHRULEDATA* branchruledata;
   VAR** lpcands;
   Real* lpcandssol;
   Real* lpcandsfrac;
   Real lpobjval;
   Real cutoffbound;
   Real bestsbdown;
   Real bestsbup;
   Real provedbound;
   Bool bestisstrongbranch;
   Bool allcolsinlp;
   Bool exactsolve;
   int nlpcands;
   int bestcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   debugMessage("Execlp method of relpscost branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get current LP objective bound of the local sub problem and global cutoff bound */
   lpobjval = SCIPgetLPObjval(scip);
   cutoffbound = SCIPgetCutoffbound(scip);

   /* check, if we want to solve the problem exactly, meaning that strong branching information is not useful
    * for cutting off sub problems and improving lower bounds of children
    */
   exactsolve = SCIPisExactSolve(scip);

   /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
   allcolsinlp = SCIPallColsInLP(scip);

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, NULL, &nlpcands) );
   assert(nlpcands > 0);

   bestcand = -1;
   bestisstrongbranch = FALSE;
   bestsbdown = SCIP_INVALID;
   bestsbup = SCIP_INVALID;
   provedbound = lpobjval;

   if( nlpcands == 1 )
   {
      /* only one candidate: nothing has to be done */
      bestcand = 0;
   }
   else
   {
      int* initcands;
      Real* initcandscores;
      int ninitcands;
      int maxninitcands;
      Real downsize;
      Real upsize;
      Real size;
      Real score;
      Real down;
      Real up;
      Real downgain;
      Real upgain;
      Real fracscore;
      Real domainscore;
      Real bestpsscore;
      Real bestpsfracscore;
      Real bestpsdomainscore;
      Real bestsbscore;
      Real bestuninitsbscore;
      Real bestsbfracscore;
      Real bestsbdomainscore;
      Real depthfac;
      Real sizefac;
      Real prio;
      Longint nodenum;
      Longint maxnsblpiterations;
      Bool usesb;
      int depth;
      int maxdepth;
      int nintvars;
      int reliable;
      int bestpscand;
      int bestsbcand;
      int maxlookahead;
      int lookahead;
      int inititer;
      int i;
      int j;
      int c;

      /* get buffer for storing the unreliable candidates */
      maxninitcands = MIN(nlpcands, branchruledata->initcand);
      CHECK_OKAY( SCIPallocBufferArray(scip, &initcands, maxninitcands+1) ); /* allocate one additional slot for convenience */
      CHECK_OKAY( SCIPallocBufferArray(scip, &initcandscores, maxninitcands+1) );
      ninitcands = 0;

      /* get current node number, depth, maximal depth, and number of binary/integer variables */
      nodenum = SCIPgetNNodes(scip);
      depth = SCIPgetDepth(scip);
      maxdepth = SCIPgetMaxDepth(scip);
      maxdepth = MAX(maxdepth, MINMAXDEPTH);
      maxdepth = MIN(maxdepth, MAXMAXDEPTH);
      nintvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);

      /* calculate value used as reliability */
      depthfac = 1.1 - (Real)depth/(Real)maxdepth;
      depthfac = MAX(depthfac, 0.0);
      depthfac = MIN(depthfac, 1.0);
      sizefac = 1.2 - sqrt((Real)nintvars/(Real)MAXSIZE);
      sizefac = MAX(sizefac, 0.0);
      sizefac = MIN(sizefac, 1.0);
      prio = depthfac * sizefac;
      reliable = (1.0-prio) * branchruledata->minreliable + prio * branchruledata->maxreliable;

      /* search for the best pseudo cost candidate, while remembering unreliable candidates in a sorted buffer */
      bestpscand = -1;
      bestpsscore = -SCIPinfinity(scip);
      bestpsfracscore = -SCIPinfinity(scip);
      bestpsdomainscore = -SCIPinfinity(scip);
      for( c = 0; c < nlpcands; ++c )
      {
         assert(lpcands[c] != NULL);
         assert(!SCIPisIntegral(scip, lpcandssol[c]));

         /* if strong branching was already applied to the variable at the current node, use the old strong branching
          * values and don't use strong branching again
          */
         if( SCIPgetVarStrongbranchNode(scip, lpcands[c]) == nodenum )
         {
            Real solval;
            Real frac;
            Real oldlpobjval;

            /* get the old strong branching values */
            CHECK_OKAY( SCIPgetVarStrongbranchLast(scip, lpcands[c], &down, &up, &solval, &oldlpobjval) );
            frac = solval - SCIPfloor(scip, solval);
            assert(SCIPisPositive(scip, frac));
            downgain = MAX(down - oldlpobjval, 0.0) * lpcandsfrac[c]/frac;
            upgain = MAX(up - oldlpobjval, 0.0) * (1.0-lpcandsfrac[c])/(1.0-frac);

            /* calculate score out of old strong branching values */
            score = SCIPgetBranchScore(scip, lpcands[c], downgain, upgain);

            /* don't use strong branching on variables that have already been initialized at the current node */
            usesb = FALSE;
         }
         else
         {
            /* get pseudo cost score of candidate */
            score = SCIPgetVarPseudocostScore(scip, lpcands[c], lpcandssol[c]);

            /* check, if the pseudo cost score of the variable is reliable */
            downsize = SCIPgetVarPseudocostCountCurrentRun(scip, lpcands[c], 0);
            upsize = SCIPgetVarPseudocostCountCurrentRun(scip, lpcands[c], 1);
            size = MIN(downsize, upsize);

            /* use strong branching on variables with unreliable pseudo cost scores */
            usesb = (size < reliable);
         }

         if( usesb )
         {
            /* pseudo cost of variable is not reliable: insert candidate in initcands buffer */
            for( j = ninitcands; j > 0 && score > initcandscores[j-1]; --j )
            {
               initcands[j] = initcands[j-1];
               initcandscores[j] = initcandscores[j-1];
            }
            initcands[j] = c;
            initcandscores[j] = score;
            ninitcands++;
            ninitcands = MIN(ninitcands, maxninitcands);
         }
         else
         {
            /* variable will keep it's pseudo cost value: check for better pseudo cost score of candidate */
            if( SCIPisSumGE(scip, score, bestpsscore) )
            {
               fracscore = MIN(lpcandsfrac[c], 1.0 - lpcandsfrac[c]);
               domainscore = -(SCIPvarGetUbLocal(lpcands[c]) - SCIPvarGetLbLocal(lpcands[c]));
               if( SCIPisSumGT(scip, score, bestpsscore)
                  || SCIPisSumGT(scip, fracscore, bestpsfracscore)
                  || (SCIPisSumGE(scip, fracscore, bestpsfracscore) && domainscore > bestpsdomainscore) )
               {
                  bestpscand = c;
                  bestpsscore = score;
                  bestpsfracscore = fracscore;
                  bestpsdomainscore = domainscore;
               }
            }
         }
      }

      /* calculate maximal number of strong branching LP iterations */
      maxnsblpiterations = branchruledata->sbiterquot * SCIPgetNNodeLPIterations(scip) + branchruledata->sbiterofs;

      /* initialize unreliable candidates with strong branching until maxlookahead is reached,
       * search best strong branching candidate
       */
      maxlookahead = branchruledata->maxlookahead;
      inititer = branchruledata->inititer;
      if( inititer == 0 )
      {
         int nlps = SCIPgetNLPs(scip);
         inititer = 2*SCIPgetNLPIterations(scip) / MAX(1, nlps);
         inititer = MAX(inititer, 10);
         inititer = MIN(inititer, 1000);
      }
      bestsbcand = -1;
      bestsbscore = -SCIPinfinity(scip);
      bestsbfracscore = -SCIPinfinity(scip);
      bestsbdomainscore = -SCIPinfinity(scip);
      lookahead = 0;
      for( i = 0; i < ninitcands && lookahead < maxlookahead
              && SCIPgetNStrongbranchLPIterations(scip) < maxnsblpiterations; ++i )
      {
         Bool lperror;
         Bool downinf;
         Bool upinf;

         /* get candidate number to initialize */
         c = initcands[i];
         assert(!SCIPisIntegral(scip, lpcandssol[c]));

         debugMessage("init pseudo cost (%g/%g) of <%s> at %g with strong branching (%d iters) -- %lld/%lld iterations\n",
            SCIPgetVarPseudocostCountCurrentRun(scip, lpcands[c], 0), 
            SCIPgetVarPseudocostCountCurrentRun(scip, lpcands[c], 1), 
            SCIPvarGetName(lpcands[c]), lpcandssol[c], inititer,
            SCIPgetNStrongbranchLPIterations(scip), maxnsblpiterations);

         /* use strong branching on candidate */
         CHECK_OKAY( SCIPgetVarStrongbranch(scip, lpcands[c], inititer, &down, &up, &lperror) );

         /* check for an error in strong branching */
         if( lperror )
         {
            SCIPmessage(scip, SCIP_VERBLEVEL_HIGH,
               "(node %lld) error in strong branching call for variable <%s> with solution %g\n", 
               SCIPgetNNodes(scip), SCIPvarGetName(lpcands[c]), lpcandssol[c]);
            break;
         }

         /* evaluate strong branching */
         down = MAX(down, lpobjval);
         up = MAX(up, lpobjval);
         downinf = SCIPisGE(scip, down, cutoffbound);
         upinf = SCIPisGE(scip, up, cutoffbound);
         downgain = down - lpobjval;
         upgain = up - lpobjval;

         /* check for possible fixings */
         if( allcolsinlp && !exactsolve )
         {
            /* because all existing columns are in LP, the strong branching bounds are feasible lower bounds */
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
            else
            {
               Real minbound;
               
               /* the minimal lower bound of both children is a proved lower bound of the current subtree */
               minbound = MIN(down, up);
               provedbound = MAX(provedbound, minbound);
            }
         }

         /* check for a better score */
         score = SCIPgetBranchScore(scip, lpcands[c], downgain, upgain);
         if( SCIPisSumGE(scip, score, bestsbscore) )
         {
            fracscore = MIN(lpcandsfrac[c], 1.0 - lpcandsfrac[c]);
            domainscore = -(SCIPvarGetUbLocal(lpcands[c]) - SCIPvarGetLbLocal(lpcands[c]));
            if( SCIPisSumGT(scip, score, bestsbscore)
               || SCIPisSumGT(scip, fracscore, bestsbfracscore)
               || (SCIPisSumGE(scip, fracscore, bestsbfracscore) && domainscore > bestsbdomainscore) )
            {
               bestsbcand = c;
               bestsbscore = score;
               bestsbdown = down;
               bestsbup = up;
               bestsbfracscore = fracscore;
               bestsbdomainscore = domainscore;
               lookahead = 0;
            }
         }
         else
            lookahead++;
      
         /* update pseudo cost values */
         if( !downinf )
         {
            CHECK_OKAY( SCIPupdateVarPseudocost(scip, lpcands[c], 0.0-lpcandsfrac[c], downgain, 1.0) );
         }
         if( !upinf )
         {
            CHECK_OKAY( SCIPupdateVarPseudocost(scip, lpcands[c], 1.0-lpcandsfrac[c], upgain, 1.0) );
         }
      
         debugMessage(" -> var <%s> (solval=%g, down=%g (%+g), up=%g (%+g), score=%g) -- best: <%s> (%g), lookahead=%d/%d\n",
            SCIPvarGetName(lpcands[c]), lpcandssol[c], down, downgain, up, upgain, score, 
            SCIPvarGetName(lpcands[bestsbcand]), bestsbscore, lookahead, maxlookahead);
      }

      /* get the score of the best uninitialized strong branching candidate */
      if( i < ninitcands )
         bestuninitsbscore = initcandscores[i];
      else
         bestuninitsbscore = -SCIPinfinity(scip);

      /* if the best pseudo cost candidate is better than the best uninitialized strong branching candidate,
       * compare it to the best initialized strong branching candidate
       */
      if( bestpsscore > bestuninitsbscore && SCIPisGT(scip, bestpsscore, bestsbscore) )
      {
         bestcand = bestpscand;
         bestisstrongbranch = FALSE;
      }
      else if( bestsbcand >= 0 )
      {
         bestcand = bestsbcand;
         bestisstrongbranch = TRUE;
      }
      else
      {
         /* no candidate was initialized, and the best score is the one of the first candidate in the initialization
          * queue
          */
         assert(ninitcands >= 1);
         bestcand = initcands[0];
         bestisstrongbranch = FALSE;
      }

      /* free buffer for the unreliable candidates */
      SCIPfreeBufferArray(scip, &initcandscores);
      SCIPfreeBufferArray(scip, &initcands);
   }

   /* if no domain could be reduced, create the branching */
   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM )
   {
      NODE* node;
      Real rootsolval;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestcand && bestcand < nlpcands);
      assert(!SCIPisIntegral(scip, lpcandssol[bestcand]));
      assert(SCIPisLT(scip, provedbound, cutoffbound));

      /* perform the branching */
      debugMessage(" -> %d cands, selected cand %d: variable <%s> (solval=%.12f, down=%g (%+g), up=%g (%+g), sb=%d)\n",
         nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandssol[bestcand], bestsbdown, bestsbdown - lpobjval,
         bestsbup, bestsbup - lpobjval, bestisstrongbranch);

      rootsolval = SCIPvarGetRootSol(lpcands[bestcand]);

      /* create child node with x <= floor(x') */
      debugMessage(" -> creating child: <%s> <= %g\n",
         SCIPvarGetName(lpcands[bestcand]), SCIPfloor(scip, lpcandssol[bestcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node, rootsolval - lpcandssol[bestcand]) );
      CHECK_OKAY( SCIPchgVarUbNode(scip, node, lpcands[bestcand], SCIPfloor(scip, lpcandssol[bestcand])) );
      if( allcolsinlp && !exactsolve && bestisstrongbranch )
      {
         assert(SCIPisLT(scip, bestsbdown, cutoffbound));
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, MAX(provedbound, bestsbdown)) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));
      
      /* create child node with x >= ceil(x') */
      debugMessage(" -> creating child: <%s> >= %g\n", 
         SCIPvarGetName(lpcands[bestcand]), SCIPceil(scip, lpcandssol[bestcand]));
      CHECK_OKAY( SCIPcreateChild(scip, &node, lpcandssol[bestcand] - rootsolval) );
      CHECK_OKAY( SCIPchgVarLbNode(scip, node, lpcands[bestcand], SCIPceil(scip, lpcandssol[bestcand])) );
      if( allcolsinlp && !exactsolve && bestisstrongbranch )
      {
         assert(SCIPisLT(scip, bestsbup, cutoffbound));
         CHECK_OKAY( SCIPupdateNodeLowerbound(scip, node, MAX(provedbound, bestsbup)) );
      }
      debugMessage(" -> child's lowerbound: %g\n", SCIPnodeGetLowerbound(node));

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsRelpscost NULL




/*
 * branching specific interface methods
 */

/** creates the reliable pseudo cost braching rule and includes it in SCIP */
RETCODE SCIPincludeBranchruleRelpscost(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   BRANCHRULEDATA* branchruledata;

   /* create relpscost branching rule data */
   CHECK_OKAY( SCIPallocMemory(scip, &branchruledata) );
   
   /* include branching rule */
   CHECK_OKAY( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                  branchFreeRelpscost, branchInitRelpscost, branchExitRelpscost, 
                  branchExeclpRelpscost, branchExecpsRelpscost,
                  branchruledata) );

   /* relpscost branching rule parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "branching/relpscost/minreliable", 
                  "minimal value for minimum pseudo cost size to regard pseudo cost value as reliable",
                  &branchruledata->minreliable, DEFAULT_MINRELIABLE, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "branching/relpscost/maxreliable", 
                  "maximal value for minimum pseudo cost size to regard pseudo cost value as reliable",
                  &branchruledata->maxreliable, DEFAULT_MAXRELIABLE, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "branching/relpscost/sbiterquot", 
                  "maximal fraction of strong branching LP iterations compared to node relaxation LP iterations",
                  &branchruledata->sbiterquot, DEFAULT_SBITERQUOT, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "branching/relpscost/sbiterofs", 
                  "additional number of allowed strong branching LP iterations",
                  &branchruledata->sbiterofs, DEFAULT_SBITEROFS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "branching/relpscost/maxlookahead", 
                  "maximal number of further variables evaluated without better score",
                  &branchruledata->maxlookahead, DEFAULT_MAXLOOKAHEAD, 1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "branching/relpscost/initcand", 
                  "maximal number of candidates initialized with strong branching per node",
                  &branchruledata->initcand, DEFAULT_INITCAND, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "branching/relpscost/inititer", 
                  "iteration limit for strong branching initializations of pseudo cost entries (0: auto)",
                  &branchruledata->inititer, DEFAULT_INITITER, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
