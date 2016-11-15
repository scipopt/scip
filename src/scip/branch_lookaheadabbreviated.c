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
   SCIP_Real*            vals;
   int                   ncandidates;
   int                   memsize;
} CANDIDATES;

typedef struct
{
   SCIP_Real*            scoresum;
   int*                  nscores;
   SCIP_Real*            maxscore;
   SCIP_Real*            minscore;
   int*                  ncutoffs;
} BRANCHINGRESULTS;

typedef struct
{
   SCIP_Bool             cutoff;
   SCIP_Real             objval;
} LPResult;

typedef struct
{
   SCIP_Bool             lperror;
   SCIP_Bool             domredcutoff;
   SCIP_Bool             domred;
} STATUS;

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static
SCIP_RETCODE allocCandidates(
   SCIP*                 scip,
   CANDIDATES**          candidates,
   int                   initialsize
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, candidates) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*candidates)->vars, initialsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*candidates)->vals, initialsize) );

   (*candidates)->ncandidates = 0;
   (*candidates)->memsize = initialsize;

   return SCIP_OKAY;
}

static
SCIP_RETCODE addCandidate(
   SCIP*                 scip,
   SCIP_VAR*             var,
   SCIP_Real             val,
   CANDIDATES*           candidates
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
   CANDIDATES**          candidates
   )
{
   SCIPfreeBufferArray(scip, &(*candidates)->vals);
   SCIPfreeBufferArray(scip, &(*candidates)->vars);
   SCIPfreeBuffer(scip, candidates);
}

static
SCIP_RETCODE allocBranchingResults(
   SCIP*                 scip,
   BRANCHINGRESULTS**    results,
   int                   nentries
   )
{
   int i;

   SCIP_CALL( SCIPallocBuffer(scip, results) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*results)->scoresum, nentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*results)->maxscore, nentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*results)->minscore, nentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*results)->ncutoffs, nentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*results)->nscores, nentries) );

   for(i = 0; i < nentries; i++)
   {
      (*results)->scoresum[i] = 0;
      (*results)->maxscore[i] = 0;
      (*results)->minscore[i] = SCIPinfinity(scip);
      (*results)->ncutoffs[i] = 0;
      (*results)->nscores[i] = 0;
   }
   return SCIP_OKAY;
}

static
void freeBranchingResults(
   SCIP*                 scip,
   BRANCHINGRESULTS**    results
   )
{
   SCIPfreeBufferArray(scip, &(*results)->nscores);
   SCIPfreeBufferArray(scip, &(*results)->ncutoffs);
   SCIPfreeBufferArray(scip, &(*results)->minscore);
   SCIPfreeBufferArray(scip, &(*results)->maxscore);
   SCIPfreeBufferArray(scip, &(*results)->scoresum);
   SCIPfreeBuffer(scip, results);
}

static
SCIP_RETCODE allocLPResult(
   SCIP*                 scip,
   LPResult**            result
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, result) );
   (*result)->cutoff = FALSE;
   return SCIP_OKAY;
}

static
void freeLPResult(
   SCIP*                 scip,
   LPResult**            result
   )
{
   SCIPfreeBuffer(scip, result);
}

static
SCIP_RETCODE allocStatus(
   SCIP*                 scip,
   STATUS**              status
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, status) );
   (*status)->lperror = FALSE;
   (*status)->domred = FALSE;
   (*status)->domredcutoff = FALSE;
   return SCIP_OKAY;
}

static
void freeStatus(
   SCIP*                 scip,
   STATUS**              status
   )
{
   SCIPfreeBuffer(scip, status);
}

/*
 * Other callback methods
 */


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

   /* minimal number of candidates to consider */
   minnumber = MIN(10, nlpcands);
   /* share of candidates to consider */
   candidatesshare = branchruledata->nsbcandidates;

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
   SCIP_BRANCHRULEDATA*  branchruledata,
   CANDIDATES*           candidates,         /**<  */
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

   SCIPdebugMessage("Calculating SB score for <%i> candidates instead of all <%i>\n", ncands, nlpcands);

   SCIP_CALL( SCIPallocBufferArray(scip, &cands, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candssol, ncands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candsfrac, ncands) );

   SCIP_CALL( getReducedNumberOfCandidates(scip, cands, candssol, candsfrac, ncands) );

   SCIP_CALL( SCIPallocBufferArray(scip, &scores, ncands) );
   for( i = 0; i < ncands; i++ )
   {
      scores[i] = -SCIPinfinity(scip);
   }

   SCIP_CALL( calculateStrongBranchingScores(scip, cands, candssol, candsfrac, scores, ncands, result) );


   SCIPdebugMessage("Sorting the candidates descending according to their SB score\n");
   /* sort the candidates descending according to their score */
   SCIPsortDownRealRealPtr(scores, candssol, (void**)cands, ncands);

   SCIPdebugMessage("Considering (at most) the <%i> candidates with the best (non-infinite) score\n", ncanditatestoadd);
   for( i = 0; i < ncands && candidates->ncandidates < ncanditatestoadd; i++ )
   {
      if( !SCIPisInfinity(scip, -scores[i]) )
      {
         addCandidate(scip, cands[i], candssol[i], candidates);

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
SCIP_RETCODE setUpperBound(
   SCIP*                 scip,
   SCIP_VAR*             var,
   SCIP_Real             val
   )
{
   SCIP_Real varlowerbound;
   SCIP_Real varupperbound;
   SCIP_Real newupperbound;

   varlowerbound = SCIPvarGetLbLocal(var);
   varupperbound = SCIPvarGetUbLocal(var);
   newupperbound = SCIPfeasFloor(scip, val);

   SCIPdebugMessage("Changing bounds for var <%s> from [<%g>..<%g>] to [<%g>..<%g>]\n", SCIPvarGetName(var),
      varlowerbound, varupperbound, varlowerbound, newupperbound);

   if( SCIPisFeasLT(scip, newupperbound, varupperbound) )
   {
      SCIP_CALL( SCIPchgVarUbProbing(scip, var, newupperbound) );
   }
   else
   {
      SCIPdebugMessage("The new upper bound <%g> is not less than the old upper bound <%g>\n", newupperbound,
         varupperbound);
      return SCIP_PARAMETERWRONGVAL;
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE setLowerBound(
   SCIP*                 scip,
   SCIP_VAR*             var,
   SCIP_Real             val
   )
{
   SCIP_Real varlowerbound;
   SCIP_Real varupperbound;
   SCIP_Real newlowerbound;

   varlowerbound = SCIPvarGetLbLocal(var);
   varupperbound = SCIPvarGetUbLocal(var);
   newlowerbound = SCIPfeasCeil(scip, val);

   SCIPdebugMessage("Changing bounds for var <%s> from [<%g>..<%g>] to [<%g>..<%g>]\n", SCIPvarGetName(var),
      varlowerbound, varupperbound, newlowerbound, varupperbound);

   if( SCIPisFeasGT(scip, newlowerbound, varlowerbound) )
   {
      SCIP_CALL( SCIPchgVarLbProbing(scip, var, newlowerbound) );
   }
   else
   {
      SCIPdebugMessage("The new lower bound <%g> is not greater than the old lower bound <%g>\n", newlowerbound,
         varlowerbound);
      return SCIP_PARAMETERWRONGVAL;
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE solveLP(
   SCIP*                 scip,
   SCIP_Bool*            lperror,
   LPResult*             result
   )
{
   SCIP_LPSOLSTAT solstat;

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, lperror, &result->cutoff) );

   solstat = SCIPgetLPSolstat(scip);

   assert(solstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

   /* for us an error occurred, if an error during the solving occurred, or the lp could not be solved but was not
    * cutoff, or if the iter or time limit was reached. */
   *lperror = *lperror || (solstat == SCIP_LPSOLSTAT_NOTSOLVED && result->cutoff == FALSE)
      || (solstat == SCIP_LPSOLSTAT_ITERLIMIT) || (solstat == SCIP_LPSOLSTAT_TIMELIMIT);

   if( !*lperror )
   {
      /* if we have no error, we save the new objective value and the cutoff decision in the resultdata */
      result->objval = SCIPgetLPObjval(scip);
      result->cutoff = result->cutoff || SCIPisGE(scip, result->objval, SCIPgetCutoffbound(scip));
      assert(((solstat != SCIP_LPSOLSTAT_INFEASIBLE) && (solstat != SCIP_LPSOLSTAT_OBJLIMIT)) || result->cutoff);
   }

   return SCIP_OKAY;
}

static
SCIP_Real getGain(
   SCIP_Real             baselpsol,
   LPResult*             result
   )
{
   if( !result->cutoff )
   {
      return MAX(0, result->objval - baselpsol);
   }
   else
   {
      /* TODO: what to return here? */
      return 0;
   }
}

static
SCIP_Real calculateScores(
   SCIP*                 scip,
   SCIP_Real             baselpsol,
   LPResult*             uuresult,
   LPResult*             ulresult,
   LPResult*             luresult,
   LPResult*             llresult
   )
{
   SCIP_Real uscore;
   SCIP_Real uugain;
   SCIP_Real ulgain;
   SCIP_Real lscore;
   SCIP_Real lugain;
   SCIP_Real llgain;

   uugain = getGain(baselpsol, uuresult);
   ulgain = getGain(baselpsol, ulresult);
   lugain = getGain(baselpsol, luresult);
   llgain = getGain(baselpsol, luresult);

   uscore = SCIPgetBranchScore(scip, NULL, uugain, ulgain);
   lscore = SCIPgetBranchScore(scip, NULL, lugain, llgain);

   return SCIPgetBranchScore(scip, NULL, uscore, lscore);
}

static
int countCutoffs(
   LPResult*             uuresult,
   LPResult*             ulresult,
   LPResult*             luresult,
   LPResult*             llresult
   )
{
   int count = 0;
   if( uuresult->cutoff )
   {
      count++;
   }
   if( ulresult->cutoff )
   {
      count++;
   }
   if( luresult->cutoff )
   {
      count++;
   }
   if( llresult->cutoff )
   {
      count++;
   }
   return count;
}

static
SCIP_RETCODE executeAbbreviatedLookaheadBranchingOnVars(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   STATUS*               status,
   VALIDDOMREDDATA*      domainreductions,
   SCIP_VAR*             firstvar,
   SCIP_Real             firstval,
   SCIP_VAR*             secondvar,
   SCIP_Real             secondval,
   SCIP_Real             baselpsol,
   SCIP_Real*            score,
   int*                  ncutoffs
   )
{
   SCIP_Real firstvallowerbound;
   SCIP_Real firstvalupperbound;
   SCIP_Real secondvallowerbound;
   SCIP_Real secondvalupperbound;
   LPResult* uuresult;
   LPResult* ulresult;
   LPResult* luresult;
   LPResult* llresult;

   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(firstvar != NULL);
   assert(secondvar != NULL);

   firstvallowerbound = SCIPvarGetLbLocal(firstvar);
   firstvalupperbound = SCIPvarGetUbLocal(firstvar);
   secondvallowerbound = SCIPvarGetLbLocal(secondvar);
   secondvalupperbound = SCIPvarGetUbLocal(secondvar);

   SCIPdebugMessage("Execute branchings on var <%s> with val <%g> in [<%g>..<%g>]\n", SCIPvarGetName(firstvar),
      firstval, firstvallowerbound, firstvalupperbound);
   SCIPdebugMessage("\tand on var <%s> with val <%g> in [<%g>..<%g>]\n", SCIPvarGetName(secondvar),
      secondval, secondvallowerbound, secondvalupperbound);

   allocLPResult(scip, &uuresult);
   allocLPResult(scip, &ulresult);
   allocLPResult(scip, &luresult);
   allocLPResult(scip, &llresult);

   if( !status->lperror )
   {
      SCIPdebugMessage("1. Adding upper bounds for both variables.\n");
      SCIP_CALL( SCIPnewProbingNode(scip) );
      SCIP_CALL( setUpperBound(scip, firstvar, firstval) );
      SCIP_CALL( setUpperBound(scip, secondvar, secondval) );
      SCIP_CALL( solveLP(scip, &status->lperror, uuresult) );

      SCIP_CALL( SCIPbacktrackProbing(scip, 0) );
   }

   if( !status->lperror )
   {
      SCIPdebugMessage("2. Adding upper bound for var <%s> and lower bound for <%s>.\n", SCIPvarGetName(firstvar),
         SCIPvarGetName(secondvar));
      SCIP_CALL( SCIPnewProbingNode(scip) );
      SCIP_CALL( setUpperBound(scip, firstvar, firstval) );
      SCIP_CALL( setLowerBound(scip, secondvar, secondval) );
      SCIP_CALL( solveLP(scip, &status->lperror, ulresult) );

      SCIP_CALL( SCIPbacktrackProbing(scip, 0) );
   }

   if( !status->lperror )
   {
      SCIPdebugMessage("3. Adding lower bound for var <%s> and upper bound for <%s>.\n", SCIPvarGetName(firstvar),
         SCIPvarGetName(secondvar));
      SCIP_CALL( SCIPnewProbingNode(scip) );
      SCIP_CALL( setLowerBound(scip, firstvar, firstval) );
      SCIP_CALL( setUpperBound(scip, secondvar, secondval) );
      SCIP_CALL( solveLP(scip, &status->lperror, luresult) );

      SCIP_CALL( SCIPbacktrackProbing(scip, 0) );
   }

   if( !status->lperror )
   {
      SCIPdebugMessage("4. Adding lower bounds for both variables.\n");
      SCIP_CALL( SCIPnewProbingNode(scip) );
      SCIP_CALL( setLowerBound(scip, firstvar, firstval) );
      SCIP_CALL( setLowerBound(scip, secondvar, secondval) );
      SCIP_CALL( solveLP(scip, &status->lperror, llresult) );

      SCIP_CALL( SCIPbacktrackProbing(scip, 0) );
   }

   if( !status->lperror )
   {
      *score = calculateScores(scip, baselpsol, uuresult, ulresult, luresult, llresult);
      *ncutoffs = countCutoffs(uuresult, ulresult, luresult, llresult);

      if( uuresult->cutoff && ulresult->cutoff )
      {
         addValidLowerBound(scip, baselpsol, firstvar, firstval, domainreductions);
      }

      if( llresult->cutoff && luresult->cutoff )
      {
         addValidUpperBound(scip, baselpsol, firstvar, firstval, domainreductions);
      }

      if( uuresult->cutoff && luresult->cutoff )
      {
         addValidLowerBound(scip, baselpsol, secondvar, secondval, domainreductions);
      }

      if( ulresult->cutoff && llresult->cutoff )
      {
         addValidUpperBound(scip, baselpsol, secondvar, secondval, domainreductions);
      }
   }


   freeLPResult(scip, &llresult);
   freeLPResult(scip, &luresult);
   freeLPResult(scip, &ulresult);
   freeLPResult(scip, &uuresult);

   return SCIP_OKAY;
}

static
void updateResult(
   BRANCHINGRESULTS*     results,
   int                   varindex,
   SCIP_Real             score,
   int                   ncutoffs
   )
{

   results->ncutoffs[varindex] = results->ncutoffs[varindex] + ncutoffs;
   results->maxscore[varindex] = MAX(results->maxscore[varindex], score);
   results->minscore[varindex] = MIN(results->maxscore[varindex], score);
   results->nscores[varindex] = results->nscores[varindex] + 1;
   results->scoresum[varindex] = results->scoresum[varindex] + score;
}

static
SCIP_RETCODE executeAbbreviatedLookaheadBranching(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   STATUS*               status,
   VALIDDOMREDDATA*      domainreductions,
   CANDIDATES*           candidates
   )
{
   BRANCHINGRESULTS* results;
   SCIP_Real baselpsol;
   int i;
   int j;

   allocBranchingResults(scip, &results, candidates->ncandidates);
   baselpsol = SCIPgetLPObjval(scip);

   SCIPstartProbing(scip);

   for( i = 0; i < candidates->ncandidates && !status->lperror; i++ )
   {
      SCIP_VAR* firstvar = candidates->vars[i];
      SCIP_Real firstval = candidates->vals[i];
      for( j = i+1; j < candidates->ncandidates && !status->lperror; j++)
      {
         SCIP_Real score;
         int ncutoffs;
         SCIP_VAR* secondvar = candidates->vars[j];
         SCIP_Real secondval = candidates->vals[j];

         SCIP_CALL( executeAbbreviatedLookaheadBranchingOnVars(scip, branchruledata, status, domainreductions, firstvar,
               firstval, secondvar, secondval, baselpsol, &score, &ncutoffs) );

         if( !status->lperror )
         {
            updateResult(results, i, score, ncutoffs);
            updateResult(results, j, score, ncutoffs);
         }
      }
   }

   SCIPendProbing(scip);

   SCIPsortDownRealRealPtr(results->maxscore, candidates->vals, (void**)candidates->vars, candidates->ncandidates);

   freeBranchingResults(scip, &results);

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
   SCIP_BRANCHRULEDATA* branchruledata;
   CANDIDATES* candidates;
   STATUS* status;
   int nlookaheadcandidates;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( result != NULL );

   SCIPdebugMessage("--- Starting Abbreviated Lookahead Branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   nlookaheadcandidates = branchruledata->nalabcandidates;

   allocCandidates(scip, &candidates, nlookaheadcandidates);
   allocStatus(scip, &status);

   SCIP_CALL( getLookaheadBranchingCandidates(scip, branchruledata, candidates, nlookaheadcandidates, result) );

   if( *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM || *result == SCIP_CONSADDED )
   {
      SCIPdebugMessage("Strong Branching finished with status <%s>\n", getStatusString(*result));
   }
   else if( candidates->ncandidates >= 1)
   {
      VALIDDOMREDDATA* domainreductions;

      allocValidBoundData(scip, &domainreductions);

      SCIPdebugMessage("Found <%i> candidates with non-infinite score for branching\n", candidates->ncandidates);

      SCIP_CALL( executeAbbreviatedLookaheadBranching(scip, branchruledata, status, domainreductions, candidates) );
      if( !status->lperror )
      {
         if( domainreductions->nboundedvars > 0 )
         {
            /* if we have no other result status set and found (potential) implied domain reductions, we add those here */
            SCIP_CALL( addDomainReductions(scip, domainreductions, &status->domredcutoff, &status->domred) );
         }

         if( status->domredcutoff )
         {
            SCIPdebugMessage("Found a domain reduction that resulted in a cutoff of the base node\n");
            *result = SCIP_CUTOFF;
         }
         else if( status->domred )
         {
            SCIPdebugMessage("Found a domain reduction\n");
            *result = SCIP_REDUCEDDOM;
         }
         else
         {
            SCIPerrorMessage("The branching has to be reimplemented, due to the update of the dualbounds.\n");
            SCIPABORT();
            /*branchOnVar(scip, candidates->vars[0], candidates->vals[0]);
            *result = SCIP_BRANCHED;*/
         }
      }
      else
      {
         *result = SCIP_DIDNOTFIND;
      }

      freeValidBoundData(scip, &domainreductions);
   }
   else
   {
      SCIPdebugMessage("Found no candidates with non-infinite score for branching\n");
   }


   freeStatus(scip, &status);
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
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead-abbreviated/sbcandidates",
         "number of candidates that should be given to the lookahead branching branching", &branchruledata->nalabcandidates,
         TRUE, DEFAULT_NALABCANDIDATES, 0, 1000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/lookahead-abbreviated/nlablpiterations",
         "number of iterations that are executed for each lookahead branching lp", &branchruledata->nalablpiterations, TRUE,
         DEFAULT_NSBLPITERATION, 0, 1000, NULL, NULL) );

   return SCIP_OKAY;
}
