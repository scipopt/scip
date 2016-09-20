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


#define BRANCHRULE_NAME            "lookahead-abbreviated"
#define BRANCHRULE_DESC            "branching rule template" /* TODO CS*/
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
   int                   randomfield;
};

typedef struct
{
   SCIP_VAR*             var;
   SCIP_Real             val;
} Candidate;

typedef struct
{
   Candidate**           candidates;
   int                   ncandidates;
} Candidates;

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/*static
SCIP_RETCODE allocCandidates(
   SCIP*                 scip,
   Candidates**          candidates,
   int                   ncandidates
)
{
   SCIP_CALL( SCIPallocBuffer(scip, candidates) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*candidates)->candidates, ncandidates) );

   return SCIP_OKAY;
}

static
void initCandidates(
   Candidates*           candidates
)
{
   candidates->ncandidates = 0;
}

static
void addCandidate(
   Candidate*            candidate,
   Candidates*           candidates
)
{
   candidates->candidates[candidates->ncandidates] = candidate;
   candidates->ncandidates++;
}

static
void freeCandidates(
   SCIP*                 scip,
   Candidates**          candidates
)
{
   SCIPfreeBufferArray(scip, &(*candidates)->candidates);
   SCIPfreeBuffer(scip, candidates);
}*/

/*
 * Other callback methods
 */

static
SCIP_DECL_SORTINDCOMP(branchCandFractionalityComparator)
{  /*lint --e{715}*/
   SCIP_Real* valuesptr = (SCIP_Real*)dataptr;
   SCIP_Real value1;
   SCIP_Real value2;

   assert(valuesptr != NULL);

   value1 = valuesptr[ind1];
   value2 = valuesptr[ind2];

   /* TODO: implement real fractionality comp */
   return value2 - value1;/*SCIPvarCompare(valuesptr[ind1], valuesptr[ind2]);*/
}


/* put your local methods here, and declare them static */

static
SCIP_RETCODE getLookaheadBranchingCandidates(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   Candidates* allcandidates;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   int* sortedindices;
   int nlpcands;
   int i;

   /* get branching candidates and their solution values (integer variables with fractional value in the current LP) */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );

   SCIP_CALL( SCIPallocBufferArray(scip, &sortedindices, nlpcands) );

   SCIPsort(sortedindices, branchCandFractionalityComparator, (void*)lpcandssol, nlpcands);

   for( i = 0; i < nlpcands; i++ )
   {
      int sortedindex = sortedindices[i];
      SCIP_VAR* sortedvar = lpcands[sortedindex];
      SCIP_Real sortedval = lpcandssol[sortedindex];

      SCIPdebugMessage("Index: <%i>, VarIndex: <%i>, Var: <%s>, Val: <%g>\n", i, sortedindex, SCIPvarGetName(sortedvar), sortedval);
   }

   SCIPfreeBufferArray(scip, &sortedindices);


   /*SCIP_CALL( allocCandidates(scip, &allcandidates, nlpcands) );
   initCandidates(allcandidates);

   for( i = 0; i < nlpcands; i++ )
   {
      Candidate* candidate;

      SCIP_CALL( SCIPallocBuffer(scip, &candidate) );
      candidate->var = lpcands[i];
      candidate->val = lpcandssol[i];

      addCandidate(candidate, allcandidates);
   }

   for( i = 0; i < allcandidates->ncandidates; i++ )
   {
      Candidate* candidate = allcandidates->candidates[i];
      SCIPdebugMessage("Added candidate <%s> with value <%g>\n", SCIPvarGetName(candidate->var), candidate->val);
   }

   for( i = allcandidates->ncandidates-1; i >= 0; i-- )
   {
      SCIPfreeBuffer(scip, &allcandidates->candidates[i]);
   }

   freeCandidates(scip, &allcandidates);*/
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
   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( result != NULL );

   SCIP_CALL( getLookaheadBranchingCandidates(scip) );

   *result = SCIP_DIDNOTRUN;

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
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
