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
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_lookaheadabbreviated.h"
#include "scip/pub_misc.h"
#include "scip/branch_fullstrong.h"


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

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Real             nsbcandidates;      /* alpha in the paper */
   SCIP_Real             nsblpiterations;    /* beta in the paper */
   SCIP_Real             nalabcandidates;    /* gamma in the paper */
   SCIP_Real             nalablpiterations;  /* delta in the paper */
};

typedef struct
{
   SCIP_VAR**            vars;
   SCIP_Real*            vals;
   int                   ncandidates;
   int                   memsize;
} Candidates;

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static
SCIP_RETCODE allocCandidates(
   SCIP*                 scip,
   Candidates**          candidates,
   int                   initialsize
)
{
   SCIP_CALL( SCIPallocBuffer(scip, candidates) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*candidates)->vars, initialsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*candidates)->vals, initialsize) );

   candidates->ncandidates = 0;
   candidates->memsize = initialsize;

   return SCIP_OKAY;
}

static
SCIP_RETCODE addCandidate(
   SCIP_VAR*             var,
   SCIP_Real             val,
   Candidates*           candidates
)
{
   int emptyindex = candidates->ncandidates;

   if( candidates->memsize == emptyindex )
   {
      int newmemsize = SCIPcalcMemGrowSize(scip, emptyindex + 1);
      SCIP_CALL( SCIPreallocBufferArray(scip, &candidates->vars, newmemsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &candidates->vals, newmemsize) );
      candidates->memsize = newmemsize;
   }
   candidates->vars[emptyindex] = var;
   candidates->vals[emptyindex] = val;
   candidates->ncandidates++;

   return SCIP_OKAY;
}

static
void freeCandidates(
   SCIP*                 scip,
   Candidates**          candidates
)
{
   SCIPfreeBufferArray(scip, &(*candidates)->vals);
   SCIPfreeBufferArray(scip, &(*candidates)->vars);
   SCIPfreeBuffer(scip, candidates);
}

/*
 * Other callback methods
 */


/*static
SCIP_DECL_SORTINDCOMP(branchCandFractionalityComparator)
{  *//*lint --e{715}*/
/*}*/

/* put your local methods here, and declare them static */

static
int getNumberOfCands(
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nlpcands
)
{
   int minnumber;
   SCIP_Real candidatesshare;

   assert(nlpcands > 0);

   minnumber = MIN(10, nlpcands); /* minimal number of candidates to consider */
   candidatesshare = branchruledata->nsbcandidates; /* share of candidates to consider */

   return MAX(candidatesshare*nlpcands, minnumber);
}

static
SCIP_RETCODE getReducedNumberOfCandidates(
   SCIP*                 scip,
   SCIP_VAR**            cands,
   SCIP_Real*            candssol,
   SCIP_Real*            candsfrac,
   int                   ncands
)
{
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   SCIP_VAR** tmplpcands;
   SCIP_Real* tmplpcandssol;
   SCIP_Real* tmplpcandsfrac;
   int nlpcands;
   int i;

   /* get branching candidates and their solution values (integer variables with fractional value in the current LP) */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &tmplpcands, &tmplpcandssol, &tmplpcandsfrac, &nlpcands, NULL, NULL) );

   assert(ncands <= nlpcands);

   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcands, tmplpcands, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandssol, tmplpcandssol, nlpcands) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &lpcandsfrac, tmplpcandsfrac, nlpcands) );

   /* TODO: maybe sort the arrays according to some (arbitrary?) rule here, otherwise remove the duplicate calls above */

   for( i = 0; i < ncands; i++ )
   {
      cands[i] = lpcands[i];
      candssol[i] = lpcandssol[i];
      candsfrac[i] = lpcandsfrac[i];
   }

   SCIPfreeBufferArray(scip, &lpcandsfrac);
   SCIPfreeBufferArray(scip, &lpcandssol);
   SCIPfreeBufferArray(scip, &lpcands);

   return SCIP_OKAY;
}

static
SCIP_RETCODE calculateStrongBranchingScores(
   SCIP*                 scip,
   SCIP_VAR**            cands,
   SCIP_Real*            candssol,
   SCIP_Real*            candsfrac,
   SCIP_Real*            scores,
   int                   ncands,
   SCIP_RESULT*          result

)
{
   SCIP_Bool* skipup;
   SCIP_Bool* skipdown;
   SCIP_Bool allowaddcons = TRUE;
   SCIP_Bool probingbounds = TRUE;
   SCIP_Bool forcestrongbranch = FALSE;
   SCIP_Bool bestdownvalid;
   SCIP_Bool bestupvalid;
   SCIP_Real provedbound;
   SCIP_Real bestdown;
   SCIP_Real bestup;
   SCIP_Real bestscore;
   int npriolpcands = ncands; /* TODO: what exactly are those priocands? */
   int lastcand = 0;
   int maxproprounds = -2;
   int bestcand;
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &skipdown, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &skipup, ncands) );

   for( i = 0; i < ncands; i++ )
   {
      skipup[i] = FALSE;
      skipdown[i] = FALSE;
   }

   SCIP_CALL( SCIPselectVarStrongBranchingRanking(scip, cands, candssol, candsfrac, skipdown, skipup, scores, ncands,
      npriolpcands, ncands, &lastcand, allowaddcons, maxproprounds, probingbounds, forcestrongbranch, &bestcand, &bestdown,
      &bestup, &bestscore, &bestdownvalid, &bestupvalid, &provedbound, result) );

   SCIPfreeBufferArray(scip, &skipup);
   SCIPfreeBufferArray(scip, &skipdown);

   return SCIP_OKAY;
}

static
SCIP_RETCODE getLookaheadBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA   branchruledata,
   Candidates*           candidates,         /**<  */
   int                   ncanditatestoadd,
   SCIP_RESULT*          result
)
{
   SCIP_VAR** cands;
   SCIP_Real* candssol;
   SCIP_Real* candsfrac;
   SCIP_Real* scores;
   int nlpcands;
   int ncands;
   int i;

   SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &nlpcands, NULL, NULL) );
   ncands = getNumberOfCands(branchruledata, nlpcands);

   SCIPdebugMessage("Considering <%i> candidates instead of the all <%i>\n", ncands, nlpcands);

   SCIP_CALL( SCIPallocBufferArray(scip, &cands, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candssol, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candsfrac, ncands) );

   SCIP_CALL( getReducedNumberOfCandidates(scip, cands, candssol, candsfrac, ncands) );

   SCIP_CALL( SCIPallocBufferArray(scip, &scores, ncands) );
   for( i = 0; i < ncands; i++ )
   {
      scores[i] = SCIPnegativeInfinity(scip);
   }

   SCIP_CALL( calculateStrongBranchingScores(scip, cands, candssol, candsfrac, scores, ncands, result) );

   /* sort the candidates descending according to their score */
   SCIPsortDownRealRealPtr(scores, candssol, (void**)cands, ncands);

   for( i = 0; i < ncands && candidates->ncandidates < ncanditatestoadd; i++ )
   {
      if( !SCIPisNegativeInfinity(scip, scores[i]) )
      {
         addCandidate(cands[i], candssol[i], candidates);

         SCIPdebugMessage("Var: <%s>, Score: <%g>, Val: <%g>\n", SCIPvarGetName(cands[i]), scores[i], candssol[i]);
      }
   }

   SCIPfreeBufferArray(scip, &scores);

   SCIPfreeBufferArray(scip, &candsfrac);
   SCIPfreeBufferArray(scip, &candssol);
   SCIPfreeBufferArray(scip, &cands);

   return SCIP_OKAY;
}

static
SCIP_RETCODE executeAbbreviatedLookaheadBranchingOnVars(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_VAR*             firstvar,
   SCIP_Real             firstval,
   SCIP_VAR*             secondvar,
   SCIP_Real             secondval
)
{
   /* TODO: implement*/


   return SCIP_OKAY;
}

static
SCIP_RETCODE executeAbbreviatedLookaheadBranching(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   Candidates*           candidates,
   SCIP_RESULT*          tmpresult
)
{
   int i;
   int j;

   SCIPstartProbing(scip);

   for( i = 0; i < candidates->ncandidates; i++ )
   {
      SCIP_VAR* firstvar = candidates->vars[i];
      SCIP_Real* firstval = candidates->vals[i];
      for( j = i+1; j < candidates->ncandidates; j++)
      {
         SCIP_VAR* secondvar = candidates->vars[j];
         SCIP_Real* secondval = candidates->vals[j];

         SCIP_CALL( executeAbbreviatedLookaheadBranchingOnVars(scip, branchruledata, firstvar, firstval, secondvar, secondval) );
      }
   }

   SCIPendProbing(scip);

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
   SCIP_RESULT tmpresult = SCIP_DIDNOTRUN;
   SCIP_BRANCHRULEDATA branchruledata;
   Candidates* candidates;
   int nlookaheadcandidates;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( result != NULL );

   SCIPdebugMessage("--- Starting Abbreviated Lookahead Branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   nlookaheadcandidates = branchruledata->nalabcandidates;

   allocCandidates(scip, &candidates, nlookaheadcandidates);

   SCIP_CALL( getLookaheadBranchingCandidates(scip, branchruledata, candidates, nlookaheadcandidates, &tmpresult) );

   if( tmpresult == SCIP_CUTOFF || tmpresult == SCIP_REDUCEDDOM || tmpresult == SCIP_CONSADDED )
   {
      *result = tmpresult;
   }
   else
   {
      SCIP_CALL( executeAbbreviatedLookaheadBranching(scip, branchruledata, candidates, &tmpresult) );
      *result = SCIP_DIDNOTRUN;
   }

   freeCandidates(scip, &candidates);

   SCIPdebugMessage("--- Finished Abbreviated Lookahead Branching\n");

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
   SCIP_CALL( SCIPaddRealParam(scip, "branching/lookahead-abbreviated/sbcandidates",
      "number of candidates that should be given to the lookahead branching branching", &branchruledata->nalabcandidates,
      TRUE, DEFAULT_NALABCANDIDATES, 0, 1, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead-abbreviated/nlablpiterations",
      "number of iterations that are executed for each lookahead branching lp", &branchruledata->nalablpiterations, TRUE,
      DEFAULT_NSBLPITERATION, 0, 1000, NULL, NULL) );

   return SCIP_OKAY;
}
