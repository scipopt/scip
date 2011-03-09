/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_pscost.c
 * @ingroup BRANCHINGRULES
 * @brief  pseudo costs branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_pscost.h"
#include "scip/cons_linear.h"


#define BRANCHRULE_NAME          "pscost"
#define BRANCHRULE_DESC          "branching on pseudo cost values"
#define BRANCHRULE_PRIORITY      2000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Real             mindistbrpointtobound; /**< minimal (relative) distance of branching point to its bounds (for continuous variables) */
   char                  strategy;           /**< strategy for computing score of relaxation branching candidates */
};

/*
 * Local methods
 */

/** determines branching point for a variable */
static
SCIP_RETCODE selectBranchingPoint(
   SCIP*             scip,           /**< SCIP data structure */
   SCIP_VAR*         var,            /**< branching variable */
   SCIP_Real         suggestion,     /**< suggestion for branching point, or SCIP_INVALID if no suggestion */
   SCIP_Real         mindistbrpointtobound, /**< minimal relative distance of branching point from variable bounds */
   SCIP_Real*        leftub,         /**< buffer to store new upper bound of variable in left  branch */
   SCIP_Real*        rightlb         /**< buffer to store new lower bound of variable in right branch */
   )
{
   SCIP_Real branchpoint;
   SCIP_Real lb, ub;

   assert(scip != NULL);
   assert(var  != NULL);
   assert(leftub  != NULL);
   assert(rightlb != NULL);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   if( suggestion != SCIP_INVALID && !SCIPisInfinity(scip, ABS(suggestion)) )
   { /* user suggested branching point */
      /* first, project it onto the current domain */
      branchpoint = MAX(lb, MIN(suggestion, ub));
      if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      { /* if it is a discrete variable, then round it down and up and accept this choice */
         *leftub  = SCIPfloor(scip, branchpoint);
         *rightlb = MAX(SCIPceil(scip, branchpoint), *leftub + 1);
         return SCIP_OKAY;
      }
      else if( (SCIPisInfinity(scip, -lb) || SCIPisGT(scip, branchpoint, lb)) && (SCIPisInfinity(scip, ub) || SCIPisLT(scip, branchpoint, ub)) )
      { /* if it is continuous and inside the box, then accept it */ 
         *leftub = *rightlb = branchpoint;
         return SCIP_OKAY;
      }
   }
   else
   { /* try the LP or pseudo LP solution */
      branchpoint = SCIPgetVarSol(scip, var);
   }

   /* if value is at +/- infty, then choose some value a bit off from bounds or 0.0 */
   if( SCIPisInfinity(scip, branchpoint) )
   { /* if value is at +infty, then the upper bound should be at infinity; choose 0.0 or something above lower bound if lower bound > 0 */
      assert(SCIPisInfinity(scip, ub));
      if( SCIPisPositive(scip, lb) )
         branchpoint = lb + 1000.0;
      else
         branchpoint = 0.0;
   }
   else if( SCIPisInfinity(scip, -branchpoint) )
   { /* if value is at -infty, then the lower bound should be at -infinity; choose 0.0 or something below upper bound if upper bound < 0 */
      assert(SCIPisInfinity(scip, -lb));
      if( SCIPisNegative(scip, ub) )
         branchpoint = ub - 1000.0;
      else
         branchpoint = 0.0;
   }

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
      { /* if branching point is too close to the bounds, move more into the middle of the interval */
         if ( branchpoint < (1.0 - mindistbrpointtobound) * lb + mindistbrpointtobound * ub )
            branchpoint = (1.0 - mindistbrpointtobound) * lb + mindistbrpointtobound * ub;
         else if( branchpoint > mindistbrpointtobound * lb + (1.0 - mindistbrpointtobound) * ub )
            branchpoint = mindistbrpointtobound * lb + (1.0 - mindistbrpointtobound) * ub;

         /* for very tiny intervals we set it exactly into the middle */
         if( SCIPisEQ(scip, lb/2.0, ub/2.0) )
            branchpoint = (lb+ub)/2.0;
      }
      else if( !SCIPisLT(scip, lb, branchpoint) )
      { /* if branching point is too close to the lower bound and there is no upper bound, then move it to somewhere away from the lower bound */
         assert(SCIPisInfinity(scip,  ub));
         branchpoint = lb + MAX(0.5*ABS(lb), 1000);
      }
      else if( !SCIPisGT(scip, ub, branchpoint) )
      { /* if branching point is too close to the upper bound and there is no lower bound, then move it to somewhere away from the upper bound */
         assert(SCIPisInfinity(scip, -lb));
         branchpoint = ub - MAX(0.5*ABS(ub), 1000);
      }

      *leftub = *rightlb = branchpoint;
   }
   else
   { /* integer variables */
      assert(SCIPisInfinity(scip,  ub) || SCIPisLE(scip, branchpoint, ub));
      assert(SCIPisInfinity(scip, -lb) || SCIPisGE(scip, branchpoint, lb));
      if( SCIPisEQ(scip, branchpoint, lb) )
      { /* if branchpoint is on lower bound, create one branch with x = lb and one with x >= lb+1 */
         *leftub  = lb;
         *rightlb = lb + 1.0;
      }
      else if( SCIPisEQ(scip, branchpoint, ub) )
      { /* if branchpoint is on upper bound, create one branch with x = ub and one with x <= ub-1 */
         *leftub  = ub - 1.0;
         *rightlb = ub;
      }
      else if( SCIPisIntegral(scip, branchpoint) )
      { /* if branchpoint is integral, create one branch with x <= x'-1 and one with x >= x'
           TODO: could in the same way be x <= x' and x >= x'+1; is there some easy way to know which is better? */
         *leftub  = branchpoint - 1.0;
         *rightlb = branchpoint;
      }
      else
      { /* branchpoint is somewhere between bounds and not fractional, so just round down and up */
         *leftub  = SCIPfloor(scip, branchpoint);
         *rightlb = SCIPceil(scip, branchpoint);
      }
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_SORTPTRCOMP(compptr)
{
   if( elem1 < elem2 )
      return -1;
   else if( elem1 == elem2 )
      return 0;
   return 1;
}

/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreePscost)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
#define branchInitPscost NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitPscost NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolPscost NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolPscost NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpPscost)
{  /*lint --e{715}*/
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real bestscore;
   SCIP_Real bestrootdiff;
   int nlpcands;
   int bestcand;
   int c;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of pscost branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, NULL, &nlpcands) );
   assert(nlpcands > 0);

   bestcand = -1;
   bestscore = -SCIPinfinity(scip);
   bestrootdiff = 0.0;
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_Real score;
      SCIP_Real rootsolval;
      SCIP_Real rootdiff;

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
   SCIPdebugMessage(" -> %d cands, selected cand %d: variable <%s> (solval=%g, score=%g)\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandssol[bestcand], bestscore);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECREL(branchExecrelPscost)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** relaxcands;
   SCIP_Real* relaxcandssol;
   SCIP_Real* relaxcandsscore;
   int npriorelaxcands;
   SCIP_VAR* cand;
   SCIP_Real candsol;
   SCIP_Real candleftub, candrightlb;
   SCIP_VAR* brvar;
   SCIP_Real leftub, rightlb;
   SCIP_Real brscore;
   SCIP_Real bestbrscore;
   SCIP_Real scoremin, scoresum, scoremax;
   SCIP_Real deltaminus, deltaplus;
   SCIP_Real pscostdown, pscostup;
   int i, j;

   SCIP_VAR** candssorted;
   int*       candsorigidx;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execrel method of random branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetRelaxBranchCands(scip, &relaxcands, &relaxcandssol, &relaxcandsscore, NULL, &npriorelaxcands, NULL, NULL, NULL) );
   assert(npriorelaxcands > 0);

   /* sort branching candidates (in a copy of relaxcands), so that same variables are on consecutive positions */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &candssorted, relaxcands, npriorelaxcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candsorigidx, npriorelaxcands) );
   for( i = 0; i < npriorelaxcands; ++i )
      candsorigidx[i] = i;

   SCIPsortPtrInt((void*)candssorted, candsorigidx, compptr, npriorelaxcands);

   bestbrscore = -1.0;
   brvar = NULL;
   leftub = 0.0; rightlb = 0.0; /* to quiet the compiler */
   for( i = 0; i < npriorelaxcands; ++i )
   {
      cand = candssorted[i];

      /* compute min, sum, and max of all registered scores for this variables
       * set candsol to a valid value, if someone registered one */
      scoremin = relaxcandsscore[candsorigidx[i]];
      scoresum = scoremin;
      scoremax = scoremin;
      candsol  = relaxcandssol[candsorigidx[i]];
      for( j = i+1 ; j < npriorelaxcands && candssorted[j] == cand; ++j )
      {
         assert(relaxcandsscore[candsorigidx[j]] >= 0.0);
         scoresum += relaxcandsscore[candsorigidx[j]];
         if( relaxcandsscore[candsorigidx[j]] < scoremin )
            scoremin = relaxcandsscore[candsorigidx[j]];
         else if( relaxcandsscore[candsorigidx[j]] > scoremax )
            scoremax = relaxcandsscore[candsorigidx[j]];

         /* TODO if there are two valid relaxcandssol available, should we take the one closer to the middle of the domain? */
         if( candsol == SCIP_INVALID || SCIPisInfinity(scip, ABS(candsol)) )
            candsol = relaxcandssol[candsorigidx[j]];
      }
      /* set i to last occurence of cand in candssorted (instead of first one as before), so in next round we look at another variable */
      i = j-1;
      assert(candssorted[i] == cand);

      switch( branchruledata->strategy )
      {
         case 'b':
            SCIP_CALL( selectBranchingPoint(scip, cand, candsol, branchruledata->mindistbrpointtobound, &candleftub, &candrightlb) );
            assert(candleftub == candrightlb || SCIPvarGetType(cand) <= SCIP_VARTYPE_IMPLINT);

            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) )
               deltaminus = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
            else
               deltaminus = candleftub - SCIPvarGetLbLocal(cand);

            if( SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)) )
               deltaplus  = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
            else
               deltaplus  = SCIPvarGetUbLocal(cand) - candrightlb;

            break;

         case 'r':
            SCIP_CALL( selectBranchingPoint(scip, cand, candsol, branchruledata->mindistbrpointtobound, &candleftub, &candrightlb) );
            assert(candleftub == candrightlb || SCIPvarGetType(cand) <= SCIP_VARTYPE_IMPLINT);

            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) )
               deltaplus  = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
            else
               deltaplus  = candleftub - SCIPvarGetLbLocal(cand);

            if( SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)) )
               deltaminus = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
            else
               deltaminus = SCIPvarGetUbLocal(cand) - candrightlb;

            break;

         case 'i':
            deltaminus = deltaplus = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
            break;

         default :
            SCIPerrorMessage("branching strategy %c unknown\n", branchruledata->strategy);
            return SCIP_ERROR;
      }

      if( SCIPisInfinity(scip, deltaminus) || SCIPisInfinity(scip, deltaplus) )
      {
         brscore = SCIPinfinity(scip);
      }
      else
      {
         pscostdown = SCIPgetVarPseudocost(scip, cand, -deltaminus);
         pscostup   = SCIPgetVarPseudocost(scip, cand,  deltaplus);
         brscore    = SCIPgetBranchScore(scip, cand, pscostdown, pscostup);
      }
      SCIPdebugMessage("branching score variable <%s> = %g; \tscore = %g; \ttype=%d  bestbrscore=%g\n", SCIPvarGetName(cand), brscore, 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax, SCIPvarGetType(cand), bestbrscore);

      if( SCIPisInfinity(scip, brscore) )
         brscore = 0.9* SCIPinfinity(scip);

      if( SCIPisSumGT(scip, brscore, bestbrscore) )
      {
         bestbrscore = brscore;
         brvar       = cand;
         leftub      = candleftub;
         rightlb     = candrightlb;
      }
      else if( SCIPisSumEQ(scip, brscore, bestbrscore) && !(SCIPisInfinity(scip, -SCIPvarGetLbLocal(brvar)) && SCIPisInfinity(scip, SCIPvarGetUbLocal(brvar))) )
      { /* if best candidate so far is bounded or unbounded at atmost one side, maybe take new candidate */
         if( (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(cand))) &&
            (SCIPisInfinity(scip, -SCIPvarGetLbLocal(brvar)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(brvar))) )
         { /* if both variables are unbounded but one of them is bounded on one side, take the one with the larger bound on this side (hope that this avoids branching on always the same variable) */
            if( SCIPvarGetUbLocal(cand) > SCIPvarGetUbLocal(brvar) || SCIPvarGetLbLocal(cand) < SCIPvarGetLbLocal(brvar) )
            {
               brvar   = cand;
               leftub  = candleftub;
               rightlb = candrightlb;
            }
         }
         else if( SCIPvarGetType(brvar) == SCIPvarGetType(cand) )
         { /* if both have the same type, take the one with larger diameter */
            if( SCIPisLT(scip, SCIPvarGetUbLocal(brvar) - SCIPvarGetLbLocal(brvar), SCIPvarGetUbLocal(cand) - SCIPvarGetLbLocal(cand)) )
            {
               brvar   = cand;
               leftub  = candleftub;
               rightlb = candrightlb;
            }
         }
         else if( SCIPvarGetType(brvar) > SCIPvarGetType(cand) )
         { /* take the one with better type ("more discrete") */
            brvar   = cand;
            leftub  = candleftub;
            rightlb = candrightlb;
         }
      }
   }

   SCIPfreeBufferArray(scip, &candssorted);
   SCIPfreeBufferArray(scip, &candsorigidx);

   if( brvar == NULL )
   {
      SCIPwarningMessage("branching variable selection failed to select a variable\n");
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   /* perform the branching */
   if( SCIPvarGetStatus(brvar) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_NODE* node;
      SCIP_CONS* cons;
      SCIP_Real  one;
      SCIP_Real  priority;
      SCIP_Real  estimate;

      one = 1.0;

      SCIPdebugMessage("branching on multiaggregated variable <%s>: new intervals: [%g, %g] and [%g, %g]\n",
         SCIPvarGetName(brvar), SCIPvarGetLbLocal(brvar), leftub, rightlb, SCIPvarGetUbLocal(brvar));

      priority = SCIPcalcNodeselPriority(scip, brvar, leftub);
      estimate = SCIPcalcChildEstimate  (scip, brvar, leftub);
      if( SCIPisInfinity(scip, estimate) )
         estimate = SCIPinfinity(scip) / 5.0;

      SCIPdebugMessage(" -> creating child: <%s> <= %g (priority: %g, estimate: %g)\n",
         SCIPvarGetName(brvar), leftub, priority, estimate);

      SCIP_CALL( SCIPcreateChild(scip, &node, priority, estimate) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &brvar, &one, SCIPvarGetLbLocal(brvar),
            leftub, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      if( leftub != rightlb )
      {
         priority = SCIPcalcNodeselPriority(scip, brvar, rightlb);
         estimate = SCIPcalcChildEstimate  (scip, brvar, rightlb);
         if( SCIPisInfinity(scip, estimate) )
            estimate = SCIPinfinity(scip) / 5.0;
      }

      SCIPdebugMessage(" -> creating child: <%s> >= %g (priority: %g, estimate: %g)\n",
         SCIPvarGetName(brvar), rightlb, priority, estimate);

      SCIP_CALL( SCIPcreateChild(scip, &node, priority, estimate) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &brvar, &one, rightlb, SCIPvarGetUbLocal(brvar),
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   else
   {
      SCIPdebugMessage("branching on variable <%s>: new intervals: [%g, %g] and [%g, %g]\n",
         SCIPvarGetName(brvar), SCIPvarGetLbLocal(brvar), leftub, rightlb, SCIPvarGetUbLocal(brvar));

      SCIP_CALL( SCIPbranchVarVal(scip, brvar, (leftub+rightlb)/2.0, NULL, NULL, NULL) );
   }

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsPscost NULL




/*
 * branching specific interface methods
 */

/** creates the pseudo cost braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchrulePscost(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create pscost branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreePscost, branchInitPscost, branchExitPscost, branchInitsolPscost, branchExitsolPscost, 
         branchExeclpPscost, branchExecrelPscost, branchExecpsPscost,
         branchruledata) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/mindistbrpointtobound",
         "minimal fractional distance of branching point to a continuous variable' bounds; a value of 0.5 leads to branching always in the middle of a bounded domain",
         &branchruledata->mindistbrpointtobound, FALSE, 0.2, 0.0001, 0.5, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "branching/"BRANCHRULE_NAME"/strategy",
         "strategy for computing score of relaxation branching candidates: b: rb-int-br, r: rb-int-br-rev, i: rb-inf",
         &branchruledata->strategy, FALSE, 'r', "bri", NULL, NULL) );

   return SCIP_OKAY;
}
