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

/**@file   branch_LookaheadAbbreviated.c
 * @brief  LookaheadAbbreviated branching rule
 * @author Christoph Schubert
 */
#define SCIP_DEBUG
#define SCIP_STATISTIC
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_lookahead.h"
#include "scip/branch_lookaheadabbreviated.h"
#include "scip/common_branch_lookahead.h"
#include "scip/pub_misc_sort.h"


#define BRANCHRULE_NAME            "lookahead-abbreviated"
#define BRANCHRULE_DESC            "branching rule template" /* TODO CS*/
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_NSBCANDIDATES                0.5
#define DEFAULT_NSBLPITERATION               10
#define DEFAULT_NALABCANDIDATES              3
#define DEFAULT_NALABLPITERATIONS            10


/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Real             nsbcandidates;      /* alpha in the paper */
   int                   nsblpiterations;    /* beta in the paper */
   int                   nalabcandidates;    /* gamma in the paper */
   int                   nalablpiterations;  /* delta in the paper */
};

typedef struct
{
   SCIP_VAR**            vars;
   SCIP_Real*            fractions;
   SCIP_Real*            scores;
   int*                  sortedindices;
} SCORING;

typedef struct
{
   SCIP_VAR*             branchvar;
   SCIP_Real             branchval;
   /* TODO: add the lp basis etc. here for re-usage? */
} CANDIDATE;

typedef struct
{
   CANDIDATE**           candidates;
   int                   ncandidates;
} CANDIDATELIST;

static
SCIP_RETCODE allocCandidate(
   SCIP*                 scip,
   CANDIDATE**           candidate
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, candidate) );

   return SCIP_OKAY;
}

static
void freeCandidate(
   SCIP*                 scip,
   CANDIDATE**           candidate
   )
{
   SCIPfreeBuffer(scip, candidate);
}

static
SCIP_RETCODE allocCandidateList(
   SCIP*                 scip,
   CANDIDATELIST**       list,
   int                   ncandidates
   )
{
   int i;

   SCIP_CALL( SCIPallocBuffer(scip, list) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*list)->candidates, ncandidates) );
   (*list)->ncandidates = ncandidates;

   for( i = 0; i < ncandidates; i++ )
   {
      SCIP_CALL( allocCandidate(scip, &(*list)->candidates[i]) );
   }

   return SCIP_OKAY;
}

static
void freeCandidateList(
   SCIP*                 scip,
   CANDIDATELIST**       list
   )
{
   int i;

   for( i = (*list)->ncandidates-1; i >= 0; i-- )
   {
      freeCandidate(scip, &(*list)->candidates[i]);
   }

   SCIPfreeBufferArray(scip, &(*list)->candidates);
   SCIPfreeBuffer(scip, list);
}

/*
 * Local methods
 */

static
SCIP_RETCODE getFSBResult(
   SCIP*                 scip,
   SCIP_Real             lpobjval,
   BRANCHRULERESULT*     branchruleresult
   )
{
   CONFIGURATION* config;
   STATUS* status;
   STATISTICS* statistics;
   LOCALSTATISTICS* localstats;

   SCIP_CALL( allocConfiguration(scip, &config) );
   SCIP_CALL( allocateStatus(scip, &status) );
   SCIP_CALL( allocStatistics(scip, &statistics, 1) );
   SCIP_CALL( allocateLocalStatistics(scip, &localstats) );

   config->usebincons = FALSE;
   config->usedomainreduction = FALSE;
   config->recursiondepth = 1;

   SCIP_CALL( selectVarStart(scip, config, NULL, status, branchruleresult, lpobjval, statistics, localstats) );

   freeLocalStatistics(scip, &localstats);
   freeStatistics(scip, &statistics);
   freeStatus(scip, &status);
   freeConfiguration(scip, &config);

   return SCIP_OKAY;
}

static
SCIP_DECL_SORTINDCOMP(branchRuleScoreComp)
{  /*lint --e{715}*/
   BRANCHRULERESULT* branchruleresult = (BRANCHRULERESULT*)dataptr;
   SCIP_Real score1;
   SCIP_Real score2;

   assert(branchruleresult != NULL);
   assert(0 <= ind1 && ind1 < branchruleresult->ncandscores);
   assert(0 <= ind2 && ind2 < branchruleresult->ncandscores);

   score1 = branchruleresult->candscores[ind1];
   score2 = branchruleresult->candscores[ind2];

   /* TODO: replace with the scip internal comparisons (containing an eps)*/
   if( score1 == score2 )
   {
      return 0;
   }
   else if( score1 < score2 )
   {
      return -1;
   }
   else
   {
      return 1;
   }
}

static
SCIP_RETCODE getNBestCandidates(
   SCIP*                 scip,
   BRANCHRULERESULT*     branchruleresult,
   CANDIDATELIST*        candidates
   )
{
   int* permutation;
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &permutation, branchruleresult->ncandscores) );

   SCIPsortDown(permutation, branchRuleScoreComp, branchruleresult, branchruleresult->ncandscores);

   SCIPdebugMessage("After sort:\n");

   SCIPdebug(
      for( i = 0; i < branchruleresult->ncandscores; i++ )
      {
         int sortedindex = permutation[i];
         SCIP_VAR* var = branchruleresult->candswithscore[sortedindex];
         SCIP_Real score = branchruleresult->candscores[sortedindex];

         SCIPdebugMessage("Index %2i: Var %s Score %g\n", i, SCIPvarGetName(var), score);
      }
   )

   SCIPdebugMessage("Selection:\n");

   for( i = 0; i < candidates->ncandidates; i++)
   {
      int sortedindex = permutation[i];
      CANDIDATE* candidate = candidates->candidates[i];

      candidate->branchvar = branchruleresult->candswithscore[sortedindex];
      candidate->branchval = branchruleresult->candslpvalue[sortedindex];

      SCIPdebugMessage("Index %i: Var %s Val %g\n", i, SCIPvarGetName(candidate->branchvar), candidate->branchval);
   }

   SCIPfreeBufferArray(scip, &permutation);

   return SCIP_OKAY;
}
/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyLookaheadAbbreviated)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   SCIP_CALL( SCIPincludeBranchruleLookaheadAbbreviated(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeLookaheadAbbreviated)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitLookaheadAbbreviated)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitLookaheadAbbreviated)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpLookaheadAbbreviated)
{  /*lint --e{715}*/
   BRANCHRULERESULT* branchruleresult;
   CANDIDATELIST* candidates;
   SCIP_Real lpobjval;
   int ncands;
   int i;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( result != NULL );

   SCIPdebugMessage("--- Starting Abbreviated Lookahead Branching\n");

   SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &ncands, NULL, NULL) );
   lpobjval = SCIPgetLPObjval(scip);

   SCIP_CALL( allocateBranchRuleResultFull(scip, &branchruleresult, lpobjval, ncands) );

   SCIP_CALL( getFSBResult(scip, lpobjval, branchruleresult) );

   SCIP_CALL( allocCandidateList(scip, &candidates, 4));

   getNBestCandidates(scip, branchruleresult, candidates);

   SCIPstartProbing(scip);

   for( i = 0; i < candidates->ncandidates; i++ )
   {
      /*SCIP_VAR* branchvar = candidates->candidates[i]->branchvar;
      SCIP_Real branchval = candidates->candidates[i]->branchval;*/

      /* TODO: down branching */

      /* TODO: up branching */

      /* TODO: evaluate the branching results */
   }

   SCIPendProbing(scip);

   freeCandidateList(scip, &candidates);
   freeBranchRuleResultFull(scip, &branchruleresult);

   SCIPdebugMessage("--- Finished Abbreviated Lookahead Branching\n");

   SCIPABORT();

   return SCIP_OKAY;
}


/*
 * branching rule specific interface methods
 */

/** creates the LookaheadAbbreviated branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleLookaheadAbbreviated(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create LookaheadAbbreviated branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   /* TODO: (optional) create branching rule specific data here */

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitLookaheadAbbreviated) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpLookaheadAbbreviated) );

   /* add LookaheadAbbreviated branching rule parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "branching/lookahead-abbreviated/nsbcandidates",
         "share of candidates that should be given to the internal strong branching", &branchruledata->nsbcandidates, TRUE,
         DEFAULT_NSBCANDIDATES, 0, 1, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead-abbreviated/nsblpiterations",
         "number of iterations that are executed for each strong branching lp", &branchruledata->nsblpiterations, TRUE,
         DEFAULT_NSBLPITERATION, 0, 1000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead-abbreviated/sbcandidates",
         "number of candidates that should be given to the lookahead branching branching", &branchruledata->nalabcandidates,
         TRUE, DEFAULT_NALABCANDIDATES, 0, 1000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead-abbreviated/nlablpiterations",
         "number of iterations that are executed for each lookahead branching lp", &branchruledata->nalablpiterations, TRUE,
         DEFAULT_NSBLPITERATION, 0, 1000, NULL, NULL) );

   return SCIP_OKAY;
}
