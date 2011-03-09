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

/**@file   branch_mostinf.c
 * @ingroup BRANCHINGRULES
 * @brief  most infeasible LP branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_mostinf.h"
#include "scip/cons_linear.h"


#define BRANCHRULE_NAME          "mostinf"
#define BRANCHRULE_DESC          "most infeasible branching"
#define BRANCHRULE_PRIORITY      100
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Real             mindistbrpointtobound; /**< minimal (relative) distance of branching point to its bounds (for continuous variables) */
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

/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeMostinf)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
#define branchInitMostinf NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitMostinf NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolMostinf NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolMostinf NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMostinf)
{  /*lint --e{715}*/
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   SCIP_Real infeasibility;
   SCIP_Real score;
   SCIP_Real obj;
   SCIP_Real bestscore;
   SCIP_Real bestobj;
   int bestcand;
   int i;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of mostinf branching\n");

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands) );
   assert(nlpcands > 0);

   /* search the most infeasible candidate */
   bestscore = SCIP_REAL_MIN;
   bestobj = 0.0;
   bestcand = -1;
   for( i = 0; i < nlpcands; ++i )
   {
      assert(lpcands[i] != NULL);

      infeasibility = lpcandsfrac[i];
      infeasibility = MIN(infeasibility, 1.0-infeasibility);
      score = infeasibility;
      score *= SCIPvarGetBranchFactor(lpcands[i]);
      obj = SCIPvarGetObj(lpcands[i]);
      obj = REALABS(obj);
      if( SCIPisGT(scip, score, bestscore)
         || (SCIPisGE(scip, score, bestscore) && obj > bestobj) )
      {
         bestscore = score;
         bestobj = obj;
         bestcand = i;
      }
   }
   assert(bestcand >= 0);

   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (frac=%g, obj=%g, factor=%g, score=%g)\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandsfrac[bestcand], bestobj,
      SCIPvarGetBranchFactor(lpcands[bestcand]), bestscore);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECREL(branchExecrelMostinf)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** relaxcands;
   SCIP_Real* relaxcandssol;
   SCIP_Real* relaxcandsscore;
   int nrelaxcands;
   SCIP_Real score;
   SCIP_Real obj;
   SCIP_Real bestscore;
   SCIP_Real bestobj;
   int bestcand;
   int i;
   SCIP_VAR* brvar;
   SCIP_Real leftub, rightlb;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execrel method of mostinf branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetRelaxBranchCands(scip, &relaxcands, &relaxcandssol, &relaxcandsscore, NULL, &nrelaxcands, NULL, NULL, NULL) );
   assert(nrelaxcands > 0);

   /* search the least infeasible candidate */
   bestscore = SCIP_REAL_MIN;
   bestobj = 0.0;
   bestcand = -1;
   for( i = 0; i < nrelaxcands; ++i )
   {
      assert(relaxcands[i] != NULL);

      score = relaxcandsscore[i];
      score *= SCIPvarGetBranchFactor(relaxcands[i]);
      obj = SCIPvarGetObj(relaxcands[i]);
      obj = REALABS(obj);
      if( SCIPisInfinity(scip, score)
         || (!SCIPisInfinity(scip, bestscore) && 
             (SCIPisGT(scip, score, bestscore) || (SCIPisGE(scip, score, bestscore) && obj > bestobj))) )
      {
         bestscore = score;
         bestobj = obj;
         bestcand = i;
      }
   }
   assert(bestcand >= 0);

   brvar = relaxcands[bestcand];

   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (infeas=%g, obj=%g, factor=%g, score=%g)\n",
      nrelaxcands, bestcand, SCIPvarGetName(brvar), relaxcandsscore[bestcand], bestobj,
      SCIPvarGetBranchFactor(brvar), bestscore);

   SCIP_CALL( selectBranchingPoint(scip, brvar, relaxcandssol[bestcand], branchruledata->mindistbrpointtobound, &leftub, &rightlb) );
   
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
      SCIP_CALL( SCIPbranchVarVal(scip, brvar, (leftub + rightlb) / 2.0, NULL, NULL, NULL) );
   }
   
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
#define branchExecpsMostinf NULL




/*
 * branching specific interface methods
 */

/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMostinf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create inference branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreeMostinf, branchInitMostinf, branchExitMostinf, branchInitsolMostinf, branchExitsolMostinf, 
         branchExeclpMostinf, branchExecrelMostinf, branchExecpsMostinf,
         branchruledata) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/mindistbrpointtobound",
         "minimal fractional distance of branching point to a continuous variable' bounds; a value of 0.5 leads to branching always in the middle of a bounded domain",
         &branchruledata->mindistbrpointtobound, FALSE, 0.2, 0.0001, 0.5, NULL, NULL) );

   return SCIP_OKAY;
}
