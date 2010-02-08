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
#pragma ident "@(#) $Id: branch_random.c,v 1.8 2010/02/08 20:07:50 bzfviger Exp $"

/**@file   branch_random.c
 * @ingroup BRANCHINGRULES
 * @brief  random variable branching rule
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_random.h"
#include "scip/cons_linear.h"


#define BRANCHRULE_NAME          "random"
#define BRANCHRULE_DESC          "random variable branching"
#define BRANCHRULE_PRIORITY      -100000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0



/** branching rule data */
struct SCIP_BranchruleData
{
   unsigned int          randseed;           /**< seed value for random number generator */
   SCIP_Real             mindistbrpointtobound; /**< minimal (relative) distance of branching point to its bounds (for continuous variables) */
};



/*
 * Callback methods
 */

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   branchruledata->randseed = 0;

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitRandom NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolRandom NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolRandom NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** lpcands;
   int nlpcands;
   int bestcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of random branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, NULL, NULL, &nlpcands) );
   assert(nlpcands > 0);

   /* get random branching candidate */
   bestcand = SCIPgetRandomInt(0, nlpcands-1, &branchruledata->randseed);
   assert(bestcand >= 0);

   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s>\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]));

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECREL(branchExecrelRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** relaxcands;
   SCIP_Real* relaxcandssol;
   int npriorelaxcands;
   int bestcand;
   SCIP_VAR* brvar;
   SCIP_Real brval;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execrel method of random branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetRelaxBranchCands(scip, &relaxcands, &relaxcandssol, NULL, NULL, &npriorelaxcands, NULL, NULL, NULL) );
   assert(npriorelaxcands > 0);

   /* get random branching candidate */
   bestcand = SCIPgetRandomInt(0, npriorelaxcands-1, &branchruledata->randseed);
   assert(bestcand >= 0);
   
   brvar = relaxcands[bestcand];
   brval = relaxcandssol[bestcand];
   assert(brvar != NULL);
   
   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> with solution value %g\n",
      npriorelaxcands, bestcand, SCIPvarGetName(brvar), brval);

   if( brval == SCIP_INVALID || SCIPisInfinity(scip, ABS(brval)) )
   { /* branching point not specified or infinity, so select one */
      SCIP_Real lb;
      SCIP_Real ub;

      brval = SCIPgetSolVal(scip, NULL, brvar);
      lb = SCIPvarGetLbLocal(brvar);
      ub = SCIPvarGetUbLocal(brvar);

      if( SCIPisInfinity(scip, brval) )
      {
         if( SCIPisPositive(scip, lb) )
            brval = lb + 1000;
         else
            brval = 0.0;
      }
      else if( SCIPisInfinity(scip, -brval) )
      {
         if( SCIPisNegative(scip, ub) )
            brval = ub - 1000;
         else
            brval = 0.0;
      }

      if( SCIPvarGetType(brvar) == SCIP_VARTYPE_CONTINUOUS )
      {
         if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
         {
            /* branch on value of LP solution
             * if it is too close to the bounds, move more into the middle of the interval */
            if( brval < (1.0 - branchruledata->mindistbrpointtobound) * lb + branchruledata->mindistbrpointtobound * ub )
               brval = (1.0 - branchruledata->mindistbrpointtobound) * lb + branchruledata->mindistbrpointtobound * ub;
            else if( brval > branchruledata->mindistbrpointtobound * lb + (1.0 - branchruledata->mindistbrpointtobound) * ub )
               brval = branchruledata->mindistbrpointtobound * lb + (1.0 - branchruledata->mindistbrpointtobound) * ub;

            /* for very tiny intervals we set it into the middle */
            if( SCIPisEQ(scip, lb/2.0, ub/2.0) )
               brval = (lb+ub)/2.0;
         }
         else if( !SCIPisLT(scip, lb, brval) )
         {
            assert(SCIPisInfinity(scip, ub));
            brval = lb + MAX(0.5*ABS(lb), 1000);
         }
         else if( !SCIPisGT(scip, ub, brval) )
         {
            assert(SCIPisInfinity(scip, -lb));
            brval = ub - MAX(0.5*ABS(ub), 1000);
         }
      }
      else
      {
         assert(ub - lb > 0.9);
         if( brval > ub )
            brval = ub - 0.5;
         else if( brval < lb )
            brval = lb + 0.5;
      }
      
      SCIPdebugMessage(" -> selected branching point %g  (bounds = [%g, %g], solval = %g)\n", brval, lb, ub, SCIPgetSolVal(scip, NULL, brvar));
   }

   /* perform the branching */
   if( SCIPvarGetStatus(brvar) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_NODE* node;
      SCIP_CONS* cons;
      SCIP_Real  one;
      SCIP_Real  priority;
      SCIP_Real  estimate;
      SCIP_Real  leftub;
      SCIP_Real  rightlb;

      one = 1.0;

      if( SCIPvarGetType(brvar) < SCIP_VARTYPE_CONTINUOUS )
      {
         leftub  = SCIPfloor(scip, brval);
         rightlb = SCIPceil (scip, brval);
         assert(SCIPisGE(scip, leftub,  SCIPvarGetLbLocal(brvar)));
         assert(SCIPisLE(scip, rightlb, SCIPvarGetUbLocal(brvar)));
      }
      else
      {
         leftub  = brval;
         rightlb = brval;
      }

      SCIPdebugMessage("branching on multiaggregated variable %s: new intervals: [%g, %g] and [%g, %g]\n", 
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

      SCIP_CALL( SCIPcreateChild(scip, &node, priority, estimate) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &brvar, &one, rightlb, SCIPvarGetUbLocal(brvar), 
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   else
   {
      SCIP_CALL( SCIPbranchVarVal(scip, brvar, brval, NULL, NULL, NULL) );
   }
   
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** pseudocands;
   int npseudocands;
   int bestcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of random branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &pseudocands, NULL, &npseudocands) );
   assert(npseudocands > 0);

   /* get random branching candidate */
   bestcand = SCIPgetRandomInt(0, npseudocands-1, &branchruledata->randseed);
   assert(bestcand >= 0);

   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s>\n",
      npseudocands, bestcand, SCIPvarGetName(pseudocands[bestcand]));

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, pseudocands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}





/*
 * branching specific interface methods
 */

/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleRandom(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create inference branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   branchruledata->randseed = 0;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchFreeRandom, branchInitRandom, branchExitRandom, branchInitsolRandom, branchExitsolRandom, 
         branchExeclpRandom, branchExecrelRandom, branchExecpsRandom,
         branchruledata) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/"BRANCHRULE_NAME"/mindistbrpointtobound",
         "minimal fractional distance of branching point to a continuous variable' bounds; a value of 0.5 leads to branching always in the middle of a bounded domain",
         &branchruledata->mindistbrpointtobound, FALSE, 0.2, 0.0001, 0.5, NULL, NULL) );

   return SCIP_OKAY;
}
