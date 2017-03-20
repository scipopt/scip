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

#include "scip/branch_fullstrong.h"
#include "scip/branch_lookaheadabbreviated.h"
#include "scip/common_branch_lookahead.h"
#include "scip/pub_misc.h"


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

/*
 * Local methods
 */

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
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** branchvars;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   int* basisind;
   int nrows;
   int i;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( result != NULL );

   SCIPdebugMessage("--- Starting Abbreviated Lookahead Branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIP_CALL( copyLPBranchCands(scip, &branchvars, &lpcandssol, &lpcandsfrac, &nlpcands) );
   SCIP_CALL( SCIPgetLPRowsData(scip, NULL, &nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );

   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );

   for( i = 0; i < nrows; i++)
   {
      SCIPdebugMessage("Element %i = %i\n", i, basisind[i]);
   }

   SCIPfreeBufferArray(scip, &basisind);
   freeLPBranchCands(scip, &branchvars, &lpcandssol, &lpcandsfrac);

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
