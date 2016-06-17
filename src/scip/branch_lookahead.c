/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*#define SCIP_DEBUG*/
/**@file   branch_lookahead.c
 * @brief  lookahead branching rule
 * @author Christoph Schubert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_lookahead.h"
#include "scip/branch_fullstrong.h"
#include "scip/var.h"

#define BRANCHRULE_NAME            "lookahead"
#define BRANCHRULE_DESC            "fullstrong branching with depth of 2" /* TODO CS: expand description */
#define BRANCHRULE_PRIORITY        536870911
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool somerandomfield;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static SCIP_RETCODE executeBranchingOnUpperBound(
   SCIP*                 scip,
   SCIP_VAR*             branchingvar,
   SCIP_Real             branchingvarsolvalue,
   SCIP_Real*            objval,
   SCIP_Bool*            cutoff
   )
{
   SCIP_Bool lperror;
   SCIP_LPSOLSTAT solstat;

   assert( scip != NULL);
   assert( branchingvar != NULL);
   assert( objval != NULL);
   assert( cutoff != NULL);

   SCIPdebugMessage("Started branching on upper bound.\n");

   SCIP_CALL( SCIPnewProbingNode(scip) );

   SCIP_CALL( SCIPchgVarUbProbing(scip, branchingvar, SCIPfeasFloor(scip, branchingvarsolvalue)) );

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, cutoff) );
   solstat = SCIPgetLPSolstat(scip);

   lperror = lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && *cutoff == 0) ||
         (solstat == SCIP_LPSOLSTAT_ITERLIMIT) || (solstat == SCIP_LPSOLSTAT_TIMELIMIT);
   assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

   if( !lperror )
   {
      *objval = SCIPgetLPObjval(scip);
      *cutoff = *cutoff || SCIPisGE(scip, *objval, SCIPgetCutoffbound(scip));
      assert(((solstat != SCIP_LPSOLSTAT_INFEASIBLE) && (solstat != SCIP_LPSOLSTAT_OBJLIMIT)) || *cutoff);
   }

   SCIPdebugMessage("Finished branching on upper bound.\n");

   return SCIP_OKAY;
}

static SCIP_Bool executeBranchingOnLowerBound(
   SCIP*                 scip,
   SCIP_VAR*             fixedvar,
   SCIP_Real             fixedvarsol,
   SCIP_Real*            upobjval,
   SCIP_Bool*            cutoff
   )
{
   SCIP_Bool lperror;
   SCIP_LPSOLSTAT solstat;

   assert( scip != NULL );
   assert( fixedvar != NULL );
   assert( upobjval != NULL );
   assert( cutoff != NULL );

   SCIPdebugMessage("Started branching on lower bound.\n");

   SCIP_CALL( SCIPnewProbingNode(scip) );
   SCIP_CALL( SCIPchgVarLbProbing(scip, fixedvar, SCIPfeasCeil(scip, fixedvarsol)) );

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, cutoff) );
   solstat = SCIPgetLPSolstat(scip);

   lperror = lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && *cutoff == 0) ||
         (solstat == SCIP_LPSOLSTAT_ITERLIMIT) || (solstat == SCIP_LPSOLSTAT_TIMELIMIT);
   assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

   if( !lperror )
   {
      *upobjval = SCIPgetLPObjval(scip);
      *cutoff = *cutoff || SCIPisGE(scip, *upobjval, SCIPgetCutoffbound(scip));
      assert(((solstat != SCIP_LPSOLSTAT_INFEASIBLE) && (solstat != SCIP_LPSOLSTAT_OBJLIMIT)) || *cutoff);
   }

   SCIPdebugMessage("Finished branching on lower bound.\n");

   return SCIP_OKAY;
}

static
SCIP_RETCODE calculateWeight(
   SCIP*                 scip,
   SCIP_Real             lowerbounddiff,
   SCIP_Real             upperbounddiff,
   SCIP_Real*            result
)
{
   SCIP_Real min;
   SCIP_Real max;
   SCIP_Real minweight = 4;
   SCIP_Real maxweight = 1;

   assert(scip != NULL);
   assert(result != NULL);

   if( SCIPisFeasGE(scip, lowerbounddiff, upperbounddiff) )
   {
      min = upperbounddiff;
      max = lowerbounddiff;
   }
   else
   {
      min = lowerbounddiff;
      max = upperbounddiff;
   }

   *result = minweight * min + maxweight * max;

   SCIPdebugMessage("The calculated weight of <%g> and <%g> is <%g>.\n", lowerbounddiff, upperbounddiff, *result);

   return SCIP_OKAY;
}

static SCIP_RETCODE executeDeepBranchingOnVar(
   SCIP*                 scip,
   SCIP_Real             lpobjval,
   SCIP_VAR*             deepbranchvar,
   SCIP_Real             deepbranchvarsolval,
   SCIP_Real*            highestweight,
   SCIP_Real*            sumweights,
   int*                  nweights,
   int*                  ncutoffs
)
{
   SCIP_Bool deepdowncutoff;
   SCIP_Bool deepupcutoff;
   SCIP_Real deepdownobjval;
   SCIP_Real deepupobjval;
   SCIP_Real upperbounddiff;
   SCIP_Real lowerbounddiff;
   SCIP_Real currentweight;

   assert(scip != NULL);
   assert(deepbranchvar != NULL);
   assert(highestweight != NULL);
   assert(sumweights != NULL);
   assert(nweights != NULL);
   assert(ncutoffs != NULL);

   SCIPdebugMessage("Second level branching on variable <%s>\n", SCIPvarGetName(deepbranchvar));

   /* NOTE CS: Start of the probe branching on x <= floor(x') and y <= floor(y') */
   SCIP_CALL(executeBranchingOnUpperBound(scip, deepbranchvar, deepbranchvarsolval, &deepdownobjval, &deepdowncutoff));

   SCIPdebugMessage("Going back to layer 1.\n");
   /* go back one layer (we are currently in depth 2) */
   SCIP_CALL(SCIPbacktrackProbing(scip, 1));

   /* NOTE CS: Start of the probe branching on x <= floor(x') and y >= ceil(y') */
   SCIP_CALL(executeBranchingOnLowerBound(scip, deepbranchvar, deepbranchvarsolval, &deepupobjval, &deepupcutoff));

   SCIPdebugMessage("Going back to layer 1.\n");
   /* go back one layer (we are currently in depth 2) */
   SCIP_CALL(SCIPbacktrackProbing(scip, 1));

   if( !deepdowncutoff && !deepupcutoff )
   {
      upperbounddiff = deepdownobjval - lpobjval;
      lowerbounddiff = deepupobjval - lpobjval;

      SCIPdebugMessage("The difference between the objective values of the base lp and the upper bounded lp is <%g>\n",
         upperbounddiff);
      SCIPdebugMessage("The difference between the objective values of the base lp and the lower bounded lp is <%g>\n",
         lowerbounddiff);

      assert(!SCIPisFeasNegative(scip, upperbounddiff));
      assert(!SCIPisFeasNegative(scip, lowerbounddiff));

      SCIP_CALL(calculateWeight(scip, lowerbounddiff, upperbounddiff, &currentweight));
      if( SCIPisFeasGE(scip, currentweight, *highestweight) )
      {
         *highestweight = currentweight;
      }
      *sumweights = *sumweights + currentweight;
      *nweights = *nweights + 1;

      SCIPdebugMessage("The sum of weights is <%g>.\n", *sumweights);
      SCIPdebugMessage("The number of weights is <%i>.\n", *nweights);
   }
   if( deepdowncutoff )
   {
      *ncutoffs = *ncutoffs + 1;
   }
   if( deepdowncutoff )
   {
      *ncutoffs = *ncutoffs + 1;
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE executeDeepBranching(
   SCIP*                 scip,
   SCIP_Real             lpobjval,
   SCIP_Real*            highestweight,
   SCIP_Real*            sumweights,
   int*                  nweights,
   int*                  ncutoffs
)
{
   SCIP_VAR**  lpcands;
   SCIP_VAR**  tmplpcands;
   SCIP_Real*  lpcandssol;
   SCIP_Real*  tmplpcandssol;
   int         nlpcands;
   int         j;

   assert(scip != NULL);
   assert(highestweight != NULL);
   assert(sumweights != NULL);
   assert(nweights != NULL);
   assert(ncutoffs != NULL);

   SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, NULL, &nlpcands, NULL, NULL) );

   /* copy LP branching candidates and solution values, because they will be updated w.r.t. the strong branching LP
    * solution
    */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcands, tmplpcands, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandssol, tmplpcandssol, nlpcands) );

   SCIPdebugMessage("The deeper lp has <%i> variables with fractional value.\n", nlpcands);

   for( j = 0; j < nlpcands; j++ )
   {
      SCIP_VAR* deepbranchvar = lpcands[j];
      SCIP_Real deepbranchvarsolval = lpcandssol[j];

      SCIP_CALL( executeDeepBranchingOnVar(scip, lpobjval, deepbranchvar, deepbranchvarsolval,
         highestweight, sumweights, nweights, ncutoffs) );
   }

   SCIPfreeBufferArray(scip, &lpcandssol);
   SCIPfreeBufferArray(scip, &lpcands);

   return SCIP_OKAY;
}

static
SCIP_RETCODE calculateCurrentWeight(
   SCIP*                 scip,
   int                   currentvarindex,
   SCIP_Real             highestweightupperbound,
   SCIP_Real             sumweightsupperbound,
   int                   nweightsupperbound,
   SCIP_Real             highestweightlowerbound,
   SCIP_Real             sumweightslowerbound,
   int                   nweightslowerbound,
   int                   ncutoffs,
   SCIP_Real*            highestweight,
   int*                  highestweightindex
)
{
   SCIP_Real lambdaupperbound = 0;
   SCIP_Real lambdalowerbound = 0;
   SCIP_Real lambda;
   SCIP_Real totalweight;

   assert(scip != NULL);
   assert(!SCIPisFeasNegative(scip, highestweightupperbound));
   assert(!SCIPisFeasNegative(scip, sumweightsupperbound));
   assert(nweightsupperbound >= 0);
   assert(!SCIPisFeasNegative(scip, highestweightlowerbound));
   assert(!SCIPisFeasNegative(scip, sumweightslowerbound));
   assert(nweightslowerbound >= 0);
   assert(!SCIPisFeasNegative(scip, ncutoffs));
   assert(highestweight != NULL);
   assert(highestweightindex != NULL);

   if( nweightsupperbound > 0 )
   {
      lambdaupperbound = (1 / nweightsupperbound) * sumweightsupperbound;
   }
   if( nweightslowerbound )
   {
      lambdalowerbound = (1 / nweightslowerbound) * sumweightslowerbound;
   }
   lambda = lambdaupperbound + lambdalowerbound;

   assert(!SCIPisFeasNegative(scip, lambda));

   SCIPdebugMessage("The lambda value is <%g>.\n", lambda);

   totalweight = highestweightlowerbound + highestweightupperbound + lambda * ncutoffs;
   if( SCIPisFeasGT(scip, totalweight, *highestweight) )
   {
      *highestweight = totalweight;
      *highestweightindex = currentvarindex;
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE selectVarLookaheadBranching(
   SCIP*                 scip,    /**< original SCIP data structure */
   SCIP_VAR**            lpcands,
   SCIP_Real*            lpcandssol,
   int                   nlpcands,
   int*                  bestcand,
   SCIP_RESULT*          result   /**< pointer to store results of branching */){

   assert(scip != NULL);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(bestcand != NULL);
   assert(result != NULL);

   if( SCIPgetDepthLimit(scip) <= (SCIPgetDepth(scip) + 2) )
   {
      SCIPdebugMessage("cannot perform probing in selectVarLookaheadBranching, depth limit reached.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( nlpcands != 0 )
   {
      SCIP_Real lpobjval;
      SCIP_Real downobjval;
      SCIP_Bool downcutoff;
      SCIP_Real upobjval;
      SCIP_Bool upcutoff;
      SCIP_Real highestweight = 0;
      int highestweightindex = -1;
      int i;

      lpobjval = SCIPgetLPObjval(scip);

      SCIPdebugMessage("The objective value of the base lp is <%g>.\n", lpobjval);

      SCIP_CALL( SCIPstartProbing(scip) );
      SCIPdebugMessage("Start Probing Mode\n");

      for( i = 0; i < nlpcands; i++ )
      {
         SCIP_Real sumweightsupperbound = 0;
         SCIP_Real sumweightslowerbound = 0;
         SCIP_Real highestweightupperbound = 0;
         SCIP_Real highestweightlowerbound = 0;
         int nweightsupperbound = 0;
         int nweightslowerbound = 0;
         int ncutoffs = 0;

         assert(lpcands[i] != NULL);

         SCIPdebugMessage("First level branching on variable <%s>\n", SCIPvarGetName(lpcands[i]));

         SCIP_CALL( executeBranchingOnUpperBound(scip, lpcands[i], lpcandssol[i], &downobjval, &downcutoff) );

         SCIPdebugMessage("objective gain: %g\n", downobjval - lpobjval);

         if( !downcutoff )
         {
            SCIP_CALL( executeDeepBranching(scip, lpobjval,
               &highestweightupperbound, &sumweightsupperbound, &nweightsupperbound, &ncutoffs) );
         }
         else
         {
            /* Both deeper branches are cutoff */
            ncutoffs = ncutoffs + 2;
         }

         SCIPdebugMessage("Going back to layer 0.\n");
         SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

         SCIP_CALL( executeBranchingOnLowerBound(scip, lpcands[i], lpcandssol[i], &upobjval, &upcutoff) );

         if( !upcutoff )
         {
            SCIP_CALL( executeDeepBranching(scip, lpobjval,
               &highestweightlowerbound, &sumweightslowerbound, &nweightslowerbound, &ncutoffs) );
         }
         else
         {
            /* Both deeper branches are cutoff */
            ncutoffs = ncutoffs + 2;
         }

         SCIPdebugMessage("Going back to layer 0.\n");
         SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

         SCIP_CALL( calculateCurrentWeight(scip, i, highestweightupperbound, sumweightsupperbound, nweightsupperbound,
            highestweightlowerbound, sumweightslowerbound, nweightslowerbound, ncutoffs,
            &highestweight, &highestweightindex) );
      }

      SCIPdebugMessage("End Probing Mode\n");
      SCIP_CALL( SCIPendProbing(scip) );

      if( highestweightindex != -1 )
      {
         *bestcand = highestweightindex;
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyLookahead)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   SCIP_CALL( SCIPincludeBranchruleLookahead(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeLookahead)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitLookahead)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitLookahead)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpLookahead)
{  /*lint --e{715}*/

   SCIP_VAR** tmplpcands;
   SCIP_VAR** lpcands;
   SCIP_Real* tmplpcandssol;
   SCIP_Real* lpcandssol;
   SCIP_Real* tmplpcandsfrac;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   int npriolpcands;
   int bestcand = -1;

   SCIPdebugMessage("Entering branchExeclpLookahead.");

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*SCIPdebugMessage("Execlp method of lookahead branching\n");*/
   *result = SCIP_DIDNOTRUN;

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, &tmplpcandsfrac, &nlpcands, &npriolpcands, NULL) );
   assert(nlpcands > 0);
   assert(npriolpcands > 0);

   /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
    * solution
    */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcands, tmplpcands, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandssol, tmplpcandssol, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandsfrac, tmplpcandsfrac, nlpcands) );

   SCIPdebugMessage("The base lp has <%i> variables with fractional value.\n", nlpcands);

   SCIP_CALL( selectVarLookaheadBranching(scip, lpcands, lpcandssol, nlpcands, &bestcand, result) );

   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED
      && 0 <= bestcand && bestcand < nlpcands )
   {
      SCIP_NODE* downchild = NULL;
      SCIP_NODE* upchild = NULL;
      SCIP_VAR* var;
      SCIP_Real val;

      assert(*result == SCIP_DIDNOTRUN);

      var = lpcands[bestcand];
      val = lpcandssol[bestcand];

      SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g)\n",
         nlpcands, bestcand, SCIPvarGetName(var), val);
      SCIP_CALL( SCIPbranchVarVal(scip, var, val, &downchild, NULL, &upchild) );

      assert(downchild != NULL);
      assert(upchild != NULL);

      SCIPdebugMessage("Branched on variable <%s>\n", SCIPvarGetName(var));
      *result = SCIP_BRANCHED;
   }
   else
   {
      SCIPdebugMessage("Could not find any variable to branch\n");
   }

   SCIPfreeBufferArray(scip, &lpcandsfrac);
   SCIPfreeBufferArray(scip, &lpcandssol);
   SCIPfreeBufferArray(scip, &lpcands);

   SCIPdebugMessage("Exiting branchExeclpLookahead.\n");

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the lookahead branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleLookahead(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create lookahead branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   /* TODO: (optional) create branching rule specific data here */

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyLookahead) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeLookahead) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitLookahead) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitLookahead) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpLookahead) );

   /* add lookahead branching rule parameters */

   return SCIP_OKAY;
}
