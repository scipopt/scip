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
#define SCIP_DEBUG
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

/* TODO CS: reduce to the needed checks and "return" (via parameter?) the objective value and the pruning status */
static SCIP_RETCODE executeBranchingOnUpperBound(
   SCIP* scip,
   SCIP_VAR* branchingvar,
   SCIP_Real branchingvarsolvalue,
   SCIP_Real* objval,
   SCIP_Bool* cutoff
   )
{
   SCIP_Bool lperror;
   SCIP_LPSOLSTAT solstat;

   assert( scip != NULL);
   assert( branchingvar != NULL);

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

/* TODO CS: reduce to the needed checks and "return" (via parameter?) the objective value and the pruning status */
static SCIP_Bool executeBranchingOnLowerBound(
   SCIP* scip,
   SCIP_VAR* fixedvar,
   SCIP_Real fixedvarsol,
   SCIP_Real* upobjval,
   SCIP_Bool* cutoff
   )
{
   SCIP_Bool lperror;
   SCIP_LPSOLSTAT solstat;

   assert( scip != NULL);
   assert( fixedvar != NULL);

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
SCIP_RETCODE calculateWeight(SCIP* scip, SCIP_Real lowerbounddiff, SCIP_Real upperbounddiff, SCIP_Real* result)
{
   SCIP_Real min;
   SCIP_Real max;
   SCIP_Real minweight = 4;
   SCIP_Real maxweight = 1;

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

   return SCIP_OKAY;
}

static SCIP_RETCODE executeDeeperBranching(
   SCIP* scip,
   SCIP_Real lpobjval,
   SCIP_VAR* deepbranchvar,
   SCIP_Real deepbranchvarsolval,
   SCIP_Real* highestweight,
   SCIP_Real* sumweights,
   SCIP_Real* nweights,
   int* ncutoffs)
{
   SCIP_Bool deepdowncutoff;
   SCIP_Bool deepupcutoff;
   SCIP_Real deepdownobjval;
   SCIP_Real deepupobjval;
   SCIP_Real upperbounddiff;
   SCIP_Real lowerbounddiff;
   SCIP_Real currentweight;

   assert(deepbranchvar != NULL);

   SCIPdebugMessage("Going to branch on variable <%s>\n", SCIPvarGetName(deepbranchvar));

   /* NOTE CS: Start of the probe branching on x <= floor(x') and y <= floor(y') */
   SCIP_CALL(executeBranchingOnUpperBound(scip, deepbranchvar, deepbranchvarsolval, &deepdownobjval, &deepdowncutoff));

   /* go back one layer (we are currently in depth 2) */
   SCIP_CALL(SCIPbacktrackProbing(scip, 1));

   /* NOTE CS: Start of the probe branching on x <= floor(x') and y >= ceil(y') */
   SCIP_CALL(executeBranchingOnLowerBound(scip, deepbranchvar, deepbranchvarsolval, &deepupobjval, &deepupcutoff));

   if( !deepdowncutoff && !deepupcutoff )
   {
      upperbounddiff = lpobjval - deepdownobjval;
      lowerbounddiff = lpobjval - deepupobjval;

      assert(SCIPisFeasPositive(scip, upperbounddiff));
      assert(SCIPisFeasPositive(scip, lowerbounddiff));

      SCIP_CALL(calculateWeight(scip, lowerbounddiff, upperbounddiff, &currentweight));
      if( SCIPisFeasGE(scip, currentweight, *highestweight) )
      {
         *highestweight = currentweight;
      }
      *sumweights = *sumweights + currentweight;
      *nweights = *nweights + 1;
   }
   if( deepdowncutoff )
   {
      *ncutoffs = *ncutoffs + 1;
   }
   if( deepdowncutoff )
   {
      *ncutoffs = *ncutoffs + 1;
   }

   /* go back one layer (we are currently in depth 2) */
   SCIP_CALL(SCIPbacktrackProbing(scip, 1));

   return SCIP_OKAY;
}

static
SCIP_RETCODE selectVarLookaheadBranching(
   SCIP*          scip,    /**< original SCIP data structure */
   SCIP_VAR**     lpcands,
   SCIP_Real*     lpcandssol,
   int            nlpcands,
   int*           bestcand,
   SCIP_RESULT*   result   /**< pointer to store results of branching */){


   assert(scip != NULL);

   if( SCIPgetDepthLimit(scip) <= (SCIPgetDepth(scip) + 2) )
   {
      SCIPdebugMessage("cannot perform probing in selectVarLookaheadBranching, depth limit reached.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( nlpcands != 0)
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

      SCIP_CALL( SCIPstartProbing(scip) );
      SCIPdebugMessage("Start Probing Mode\n");

      for( i = 0; i < nlpcands; i++ )
      {
         SCIP_Real sumweightupperbound = 0;
         SCIP_Real sumweightsupperbound = 0;
         SCIP_Real sumweightlowerbound = 0;
         SCIP_Real sumweightslowerbound = 0;
         SCIP_Real highestweightupperbound = 0;
         SCIP_Real highestweightlowerbound = 0;
         SCIP_Real lambda;
         SCIP_Real totalweight;
         int ncutoffs = 0;

         assert(lpcands[i] != NULL);
         assert(lpcandssol[i] != NULL);

         /* NOTE CS: Start of the probe branching on x <= floor(x') */
         SCIP_CALL( executeBranchingOnUpperBound(scip, lpcands[i], lpcandssol[i], &downobjval, &downcutoff) );

         if( !downcutoff )
         {
            SCIP_VAR**  ublpcands;
            SCIP_VAR**  ubtmplpcands;
            SCIP_Real*  ublpcandssol;
            SCIP_Real*  ubtmplpcandssol;
            int         ubnlpcands;

            SCIP_CALL( SCIPgetLPBranchCands(scip, &ubtmplpcands, &ubtmplpcandssol, NULL, &ubnlpcands, NULL, NULL) );
            assert(ubnlpcands > 0);

            /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
             * solution
             */
            SCIP_CALL( SCIPduplicateBufferArray(scip, &ublpcands, ubtmplpcands, nlpcands) );
            SCIP_CALL( SCIPduplicateBufferArray(scip, &ublpcandssol, ubtmplpcandssol, nlpcands) );

            if( ubnlpcands != 0 )
            {
               int j;

               for( j = 0; j < ubnlpcands; j++ )
               {

                  SCIP_VAR* deepbranchvar = ublpcands[j];
                  SCIP_Real deepbranchvarsolval = ublpcandssol[j];

                  SCIP_CALL( executeDeeperBranching(scip, lpobjval, deepbranchvar, deepbranchvarsolval,
                     &highestweightupperbound, &sumweightupperbound, &sumweightsupperbound, &ncutoffs) );
               }
            }
            SCIPfreeBufferArray(scip, &ublpcandssol);
            SCIPfreeBufferArray(scip, &ublpcands);
         }

         SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

         SCIP_CALL( executeBranchingOnLowerBound(scip, lpcands[i], lpcandssol[i], &upobjval, &upcutoff) );

         if( !upcutoff )
         {
            SCIP_VAR**  lblpcands;
            SCIP_VAR**  lbtmplpcands;
            SCIP_Real*  lblpcandssol;
            SCIP_Real*  lbtmplpcandssol;
            int         lbnlpcands;

            SCIP_CALL( SCIPgetLPBranchCands(scip, &lbtmplpcands, &lbtmplpcandssol, NULL, &lbnlpcands, NULL, NULL) );
            assert(lbnlpcands > 0);

            /* copy LP banching candidates and solution values, because they will be updated w.r.t. the strong branching LP
             * solution
             */
            SCIP_CALL( SCIPduplicateBufferArray(scip, &lblpcands, lbtmplpcands, nlpcands) );
            SCIP_CALL( SCIPduplicateBufferArray(scip, &lblpcandssol, lbtmplpcandssol, nlpcands) );

            if( lbnlpcands != 0 )
            {
               int k;

               for( k = 0; k < lbnlpcands; k++ )
               {
                  SCIP_VAR* deepbranchvar = lblpcands[k];
                  SCIP_Real deepbranchvarsolval = lblpcandssol[k];

                  SCIP_CALL( executeDeeperBranching(scip, lpobjval, deepbranchvar, deepbranchvarsolval,
                     &highestweightlowerbound, &sumweightlowerbound, &sumweightslowerbound, &ncutoffs) );
               }
            }
         }
         SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

         lambda = (1/sumweightsupperbound)*sumweightupperbound + (1/sumweightslowerbound)*sumweightlowerbound;
         totalweight = highestweightlowerbound + highestweightupperbound + lambda*ncutoffs;
         if( SCIPisFeasGT(scip, totalweight, highestweight) )
         {
            highestweight = totalweight;
            highestweightindex = i;
         }
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

   SCIP_CALL( selectVarLookaheadBranching(scip, lpcands, lpcandssol, nlpcands, &bestcand, result) );

   if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED )
   {
      SCIP_NODE* downchild = NULL;
      SCIP_NODE* upchild = NULL;
      SCIP_VAR* var;
      SCIP_Real val;

      assert(*result == SCIP_DIDNOTRUN);
      assert(0 <= bestcand && bestcand < nlpcands);

      var = lpcands[bestcand];
      val = lpcandssol[bestcand];

      SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s> (solval=%g)\n",
         nlpcands, bestcand, SCIPvarGetName(var), lpcandssol[bestcand]);
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
