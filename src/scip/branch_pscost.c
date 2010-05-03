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
#pragma ident "@(#) $Id: branch_pscost.c,v 1.30 2010/05/03 16:32:45 bzfviger Exp $"

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
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< branching variable */
   SCIP_Real             suggestion,         /**< suggestion for branching point, or SCIP_INVALID if no suggestion */
   SCIP_Real             mindistbrpointtobound, /**< minimal relative distance of branching point from variable bounds */
   SCIP_Real*            leftub,             /**< buffer to store new upper bound of variable in left  branch */
   SCIP_Real*            rightlb             /**< buffer to store new lower bound of variable in right branch */
   )
{
   SCIP_Real branchpoint;
   SCIP_Real lb;
   SCIP_Real ub;
   
   assert(scip != NULL);
   assert(var  != NULL);
   assert(leftub  != NULL);
   assert(rightlb != NULL);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   
   if( !SCIPisInfinity(scip, REALABS(suggestion)) )
   { 
      /* use user suggested branching point */

      /* first, project it onto the current domain */
      branchpoint = MAX(lb, MIN(suggestion, ub));
      
      if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      { 
         /* if it is a discrete variable, then round it down and up and accept this choice */
         if( SCIPisEQ(scip, branchpoint, ub) )
         {
            (*leftub)  = SCIPfloor(scip, branchpoint) - 1.0;
            (*rightlb) = (*leftub) + 1.0;
         }
         else
         {
            (*leftub)  = SCIPfloor(scip, branchpoint);
            (*rightlb) = (*leftub) + 1.0;
         }
         return SCIP_OKAY;
      }
      else if( (SCIPisInfinity(scip, -lb) || SCIPisGT(scip, branchpoint, lb)) && (SCIPisInfinity(scip, ub) || SCIPisLT(scip, branchpoint, ub)) )
      {
         /* if it is continuous and inside the box, then accept it */ 
         (*leftub)  = branchpoint;
         (*rightlb) = branchpoint;
         return SCIP_OKAY;
      }
   }
   else
   { 
      /* try the LP or pseudo LP solution */
      branchpoint = SCIPgetVarSol(scip, var);
   }

   /* if value is at +/- infty, then choose some value a bit off from bounds or 0.0 */
   if( SCIPisInfinity(scip, branchpoint) )
   { 
      /* if value is at +infty, then the upper bound should be at infinity; choose 0.0 or something above lower bound if lower bound > 0 */
      assert(SCIPisInfinity(scip, ub));

      if( SCIPisPositive(scip, lb) )
         branchpoint = lb + 1000.0;
      else
         branchpoint = 0.0;
   }
   else if( SCIPisInfinity(scip, -branchpoint) )
   { 
      /* if value is at -infty, then the lower bound should be at -infinity; choose 0.0 or something below upper bound if upper bound < 0 */
      assert(SCIPisInfinity(scip, -lb));
      if( SCIPisNegative(scip, ub) )
         branchpoint = ub - 1000.0;
      else
         branchpoint = 0.0;
   }

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
      { 
         /* if branching point is too close to the bounds, move more into the middle of the interval */
         if( ub - lb < 2.02*SCIPepsilon(scip) )
         { 
            /* for very tiny intervals we set it exactly into the middle */
            branchpoint = (lb+ub)/2.0;
         }
         else
         { 
            /* otherwise we project it away from the bounds */
            SCIP_Real minbrpoint;
            SCIP_Real maxbrpoint;

            minbrpoint = (1.0 - mindistbrpointtobound) * lb + mindistbrpointtobound * ub;
            minbrpoint = MAX(lb + 1.01*SCIPepsilon(scip), minbrpoint);

            maxbrpoint = mindistbrpointtobound * lb + (1.0 - mindistbrpointtobound) * ub;
            maxbrpoint = MIN(ub - 1.01*SCIPepsilon(scip), maxbrpoint);

            branchpoint = MAX(minbrpoint, MIN(branchpoint, maxbrpoint));
         }
         assert(SCIPisLT(scip, lb, branchpoint));
         assert(SCIPisLT(scip, branchpoint, ub));
      }
      else if( !SCIPisLT(scip, lb, branchpoint) )
      { 
         /* if branching point is too close to the lower bound and there is no upper bound, then move it to somewhere away from the lower bound */
         assert(SCIPisInfinity(scip,  ub));
         branchpoint = lb + MAX(0.5*REALABS(lb), 1000);
      }
      else if( !SCIPisGT(scip, ub, branchpoint) )
      { 
         /* if branching point is too close to the upper bound and there is no lower bound, then move it to somewhere away from the upper bound */
         assert(SCIPisInfinity(scip, -lb));
         branchpoint = ub - MAX(0.5*REALABS(ub), 1000);
      }

      *leftub = *rightlb = branchpoint;
   }
   else
   { 
      /* integer variables */
      assert(SCIPisInfinity(scip,  ub) || SCIPisLE(scip, branchpoint, ub));
      assert(SCIPisInfinity(scip, -lb) || SCIPisGE(scip, branchpoint, lb));
      if( SCIPisEQ(scip, branchpoint, lb) )
      { 
         /* if branchpoint is on lower bound, create one branch with x = lb and one with x >= lb+1 */
         (*leftub)  = lb;
         (*rightlb) = lb + 1.0;
      }
      else if( SCIPisEQ(scip, branchpoint, ub) )
      { 
         /* if branchpoint is on upper bound, create one branch with x = ub and one with x <= ub-1 */
         (*leftub)  = ub - 1.0;
         (*rightlb) = ub;
      }
      else if( SCIPisIntegral(scip, branchpoint) )
      { 
         /* if branchpoint is integral, create one branch with x <= x'-1 and one with x >= x'
          * @todo could in the same way be x <= x' and x >= x'+1; is there some easy way to know which is better? */
         (*leftub)  = branchpoint - 1.0;
         (*rightlb) = branchpoint;
      }
      else
      { 
         /* branchpoint is somewhere between bounds and not fractional, so just round down and up */
         (*leftub) = SCIPfloor(scip, branchpoint);
         (*rightlb) = SCIPceil(scip, branchpoint);
      }
   }
   
   return SCIP_OKAY;
}

/** selects the branching variable from given candidate array */
static
SCIP_RETCODE selectBranchVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_VAR**            cands,              /**< array of branching candidates */
   SCIP_Real*            candssol,           /**< array of candidate solution values */
   SCIP_Real*            candsscore,         /**< array of candidate scores */
   int                   ncands,             /**< the number of candidates */
   SCIP_VAR**            brvar,              /**< pointer to store the selected branching candidate or NULL if none */
   SCIP_Real*            rightlb,            /**< lower bound for branching variable for the right child */
   SCIP_Real*            leftub              /**< upper bound for branching variable for the left child */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   SCIP_VAR* cand;
   SCIP_Real candsol;
   SCIP_Real candleftub;
   SCIP_Real candrightlb;

   SCIP_Real branchscore;
   SCIP_Real bestbranchscore;

   SCIP_Real scoremin;
   SCIP_Real scoresum;
   SCIP_Real scoremax;

   SCIP_Real deltaminus;
   SCIP_Real deltaplus;

   SCIP_Real pscostdown;
   SCIP_Real pscostup;


   SCIP_VAR** candssorted;
   int* candsorigidx;
   
   int i;
   int j;
   
   (*brvar) = NULL;
   (*rightlb) = 0.0; 
   (*leftub) = 0.0; 

   if( ncands == 0 )
      return SCIP_OKAY;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   
   /* sort branching candidates (in a copy), such that same variables are on consecutive positions */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &candssorted, cands, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candsorigidx, ncands) );
   for( i = 0; i < ncands; ++i )
      candsorigidx[i] = i;
   
   SCIPsortPtrInt((void*)candssorted, candsorigidx, SCIPvarComp, ncands);

   candrightlb = 0.0;
   candleftub = 0.0;
   bestbranchscore = -1.0;

   for( i = 0; i < ncands; ++i )
   {
      cand = candssorted[i];

      /* compute min, sum, and max of all registered scores for this variables
       * set candsol to a valid value, if someone registered one */
      scoremin = candsscore[candsorigidx[i]];
      scoresum = scoremin;
      scoremax = scoremin;
      candsol  = candssol[candsorigidx[i]];
      for( j = i+1 ; j < ncands && SCIPvarCompare(candssorted[j], cand) == 0; ++j )
      {
         assert(candsscore[candsorigidx[j]] >= 0.0);
         scoresum += candsscore[candsorigidx[j]];
         if( candsscore[candsorigidx[j]] < scoremin )
            scoremin = candsscore[candsorigidx[j]];
         else if( candsscore[candsorigidx[j]] > scoremax )
            scoremax = candsscore[candsorigidx[j]];

         /* @todo if there are two valid relaxcandssol available, should we take the one closer to the middle of the domain? */
         if( SCIPisInfinity(scip, REALABS(candsol)) )
            candsol = candssol[candsorigidx[j]];
      }
      /* set i to last occurence of cand in candssorted (instead of first one as before), so in next round we look at another variable */
      i = j-1;
      assert(candssorted[i] == cand);

      switch( branchruledata->strategy )
      {
      case 'b':
         SCIP_CALL( selectBranchingPoint(scip, cand, candsol, branchruledata->mindistbrpointtobound, &candleftub, &candrightlb) );
         assert(SCIPisEQ(scip, candleftub, candrightlb) || SCIPvarGetType(cand) <= SCIP_VARTYPE_IMPLINT);
         
         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) )
            deltaminus = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
         else
            deltaminus = candleftub - SCIPvarGetLbLocal(cand);

         if( SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)) )
            deltaplus = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
         else
            deltaplus = SCIPvarGetUbLocal(cand) - candrightlb;
            
         break;
      case 'r':
         SCIP_CALL( selectBranchingPoint(scip, cand, candsol, branchruledata->mindistbrpointtobound, &candleftub, &candrightlb) );
         assert(SCIPisEQ(scip, candleftub, candrightlb) || SCIPvarGetType(cand) <= SCIP_VARTYPE_IMPLINT);

         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) )
            deltaplus = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
         else
            deltaplus = candleftub - SCIPvarGetLbLocal(cand);

         if( SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)) )
            deltaminus = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
         else
            deltaminus = SCIPvarGetUbLocal(cand) - candrightlb;

         break;

      case 'i':
         deltaplus = SCIPisInfinity(scip, scoremax) ? SCIPinfinity(scip) : 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax;
         deltaminus = deltaplus;
         break;

      default :
         SCIPerrorMessage("branching strategy %c unknown\n", branchruledata->strategy);
         return SCIP_ERROR;
      }

      if( SCIPisInfinity(scip, deltaminus) || SCIPisInfinity(scip, deltaplus) )
         branchscore = SCIPinfinity(scip);
      else
      {
         pscostdown = SCIPgetVarPseudocostVal(scip, cand, -deltaminus);
         pscostup = SCIPgetVarPseudocostVal(scip, cand,  deltaplus);
         branchscore = SCIPgetBranchScore(scip, cand, pscostdown, pscostup);
      }
      SCIPdebugMessage("branching score variable <%s> = %g; score = %g; type=%d bestbrscore=%g\n", 
         SCIPvarGetName(cand), branchscore, 0.1 * scoresum + 0.8 * scoremin + 1.3 * scoremax, 
         SCIPvarGetType(cand), bestbranchscore);

      if( SCIPisInfinity(scip, branchscore) )
         branchscore = 0.9*SCIPinfinity(scip);
      
      if( SCIPisSumGT(scip, branchscore, bestbranchscore) )
      {
         bestbranchscore = branchscore;
         (*brvar) = cand;
         (*leftub) = candleftub;
         (*rightlb) = candrightlb;
      }
      else if( SCIPisSumEQ(scip, branchscore, bestbranchscore) 
         && !(SCIPisInfinity(scip, -SCIPvarGetLbLocal(*brvar)) && SCIPisInfinity(scip, SCIPvarGetUbLocal(*brvar))) )
      { 
         /* if best candidate so far is bounded or unbounded at atmost one side, maybe take new candidate */
         if( (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(cand))) &&
            (SCIPisInfinity(scip, -SCIPvarGetLbLocal(*brvar)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(*brvar))) )
         { 
            /* if both variables are unbounded but one of them is bounded on one side, take the one with the larger bound on this side (hope that this avoids branching on always the same variable) */
            if( SCIPvarGetUbLocal(cand) > SCIPvarGetUbLocal(*brvar) || SCIPvarGetLbLocal(cand) < SCIPvarGetLbLocal(*brvar) )
            {
               (*brvar) = cand;
               (*leftub) = candleftub;
               (*rightlb) = candrightlb;
            }
         }
         else if( SCIPvarGetType(*brvar) == SCIPvarGetType(cand) )
         { 
            /* if both have the same type, take the one with larger diameter */
            if( SCIPisLT(scip, SCIPvarGetUbLocal(*brvar) - SCIPvarGetLbLocal(*brvar), SCIPvarGetUbLocal(cand) - SCIPvarGetLbLocal(cand)) )
            {
               (*brvar) = cand;
               (*leftub) = candleftub;
               (*rightlb) = candrightlb;
            }
         }
         else if( SCIPvarGetType(*brvar) > SCIPvarGetType(cand) )
         { 
            /* take the one with better type ("more discrete") */
            (*brvar) = cand;
            (*leftub) = candleftub;
            (*rightlb) = candrightlb;
         }
      }
   }
   
   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &candssorted);
   SCIPfreeBufferArray(scip, &candsorigidx);

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyPscost)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchrulePscost(scip) );
 
   return SCIP_OKAY;
}

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
   SCIP_VAR** relaxcands;
   SCIP_Real* relaxcandssol;
   SCIP_Real* relaxcandsscore;
   int npriorelaxcands;
   SCIP_VAR* brvar;
   SCIP_Real leftub;
   SCIP_Real rightlb;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   
   SCIPdebugMessage("Execrel method of pscost branching\n");
   
   /* get branching candidates */
   SCIP_CALL( SCIPgetRelaxBranchCands(scip, &relaxcands, &relaxcandssol, &relaxcandsscore, NULL, &npriorelaxcands, NULL, NULL, NULL) );
   assert(npriorelaxcands > 0);
   
   /* select braching variable */
   SCIP_CALL( selectBranchVar(scip, branchrule, relaxcands, relaxcandssol, relaxcandsscore, npriorelaxcands, &brvar, &rightlb, &leftub) );
   
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
      
      if( !SCIPisEQ(scip, leftub, rightlb) )
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
         branchCopyPscost,
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

/** selects a branching variable, due to pseudo cost, from the given candidate array and returns this variable together
 *  with a lower and upper bound for the right and left branch, respectively */
SCIP_RETCODE SCIPselectBranchVarPscost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            branchcands,        /**< brancing candidates */
   SCIP_Real*            branchcandssol,     /**< solution value for the branching candidates */
   SCIP_Real*            branchcandsscore,   /**< array of candidate scores */
   int                   nbranchcands,       /**< number of branching candidates */
   SCIP_VAR**            var,                /**< pointer to store the variable to branch on, or NULL if none */
   SCIP_Real*            rightlb,            /**< pointer to store the lower bound for the branching variable in the right brach */
   SCIP_Real*            leftub              /**< pointer to store thr upper bound for the branching variable in the left brach */
   )
{
   SCIP_BRANCHRULE* branchrule;
   
   assert(scip != NULL);
   
   /* find branching rule */
   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);
   
   /* select braching variable */
   SCIP_CALL( selectBranchVar(scip, branchrule, branchcands, branchcandssol, branchcandsscore, nbranchcands, var, rightlb, leftub) );

   return SCIP_OKAY;
}
