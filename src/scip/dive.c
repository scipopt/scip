/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip.c
 * @brief  library methods for diving heuristics
 * @author Gregor Hendel
 *
 * @todo check all checkStage() calls, use bit flags instead of the SCIP_Bool parameters
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/* the indicator constraint handler is included for the diving algorithm SCIPperformGenericDivingAlgorithm() */
#include "scip/pub_dive.h"
#include "pub_heur.h"

/* the indicator constraint handler is included for the diving algorithm SCIPperformGenericDivingAlgorithm() */
#include "scip/cons_indicator.h"

/** get candidates for diving from indicator constraints */
static
SCIP_RETCODE getIndCandVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           indconss,           /**< indicator constraints */
   int                   nindconss,          /**< number of indicator constraints */
   SCIP_VAR**            indcands,           /**< indicator candidate variables */
   SCIP_Real*            indcandssol,        /**< solution values of candidates */
   SCIP_Real*            indcandfrac,        /**< fractionalities of candidates */
   int*                  nindcands           /**< number of candidates */
   )
{
   SCIP_VAR* binvar;
   SCIP_Real val;
   int c;

   assert(scip != NULL);
   assert(indconss != NULL);
   assert(indcands != NULL);
   assert(nindcands != NULL);
   assert(indcandssol != NULL);
   assert(indcandfrac != NULL);

   *nindcands = 0;
   for( c = 0; c < nindconss; ++c )
   {
      /* check whether constraint is violated */
      if( SCIPisViolatedIndicator(scip, indconss[c], NULL) )
      {
         binvar = SCIPgetBinaryVarIndicator(indconss[c]);
         val = SCIPgetSolVal(scip, NULL, binvar);

         /* fractional indicator variables are treated by lpcands */
         if( SCIPisFeasIntegral(scip, val) )
         {
            indcands[*nindcands] = binvar;
            indcandssol[*nindcands] = val;
            indcandfrac[*nindcands] = SCIPfrac(scip, val);
            ++(*nindcands);
         }
      }
   }

   return SCIP_OKAY;
}
/** assign the current candidates for diving, together with their LP solution values and fractionalities */
static
SCIP_RETCODE getDivingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< settings for diving */
   SCIP_CONSHDLR*        indconshdlr,        /**< indicator constraint handler */
   SCIP_VAR***           lpcands,            /**< pointer to store the array of LP candidates */
   SCIP_Real**           lpcandssol,         /**< pointer to store the array of LP solution values */
   SCIP_Real**           lpcandsfrac,        /**< pointer to store the array of LP fractionalities*/
   int*                  nlpcands,           /**< pointer to store the number of LP branching candidates */
   SCIP_VAR**            indcands,           /**< pointer to store the indicator candidate variables, or NULL */
   SCIP_Real*            indcandssol,        /**< pointer to store the indicator LP solution values, or NULL */
   SCIP_Real*            indcandsfrac,       /**< pointer to store the indicator candidate fractionalities, or NULL */
   int*                  nindcands           /**< pointer to store the number of indicator candidate variables, or NULL */
   )
{

   int nindconss;
   assert(diveset != NULL);
   assert(nlpcands != NULL);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(lpcandsfrac != NULL);

   assert(lpcands != NULL);
   assert(lpcandssol != NULL);
   assert(lpcandsfrac != NULL);

   assert((indcands == NULL) == (indcandssol == NULL));
   assert((indcandssol == NULL) == (indcandsfrac == NULL));

   /* store the LP candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, lpcands, lpcandssol, lpcandsfrac, nlpcands, NULL, NULL) );

   nindconss = 0;
   if( indconshdlr != NULL )
      nindconss = SCIPconshdlrGetNConss(indconshdlr);

   /* get indicator candidates */
   if( nindconss > 0 )
   {
      SCIP_CONS** indconss;

      indconss = SCIPconshdlrGetConss(indconshdlr);

      SCIP_CALL( getIndCandVars(scip, indconss, nindconss, indcands, indcandssol, indcandsfrac, nindcands) );
   }
   else if( nindcands != NULL )
      *nindcands = 0;

   return SCIP_OKAY;
}

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */

/** performs a diving within the limits of the diveset parameters
 *
 *  This method performs a diving according to the settings defined by the diving settings @p diveset; Contrary to the
 *  name, SCIP enters probing mode (not diving mode) and dives along a path into the tree. Domain propagation
 *  is applied at every node in the tree, whereas probing LPs might be solved less frequently.
 *
 *  Starting from the current LP candidates, the algorithm determines a fraction of the candidates that should be
 *  branched on; if a single candidate should be fixed, the algorithm selects a candidate which minimizes the
 *  score defined by the @p diveset.
 *  If more than one candidate should be selected, the candidates are sorted in non-decreasing order
 *  of their score.
 *
 *  The algorithm iteratively selects the the next (unfixed) candidate in the list, until the
 *  targeted depth is reached, or the last node is proven to be infeasible. It optionally backtracks and tries the
 *  other branching direction.
 *
 *  After the set of remaining candidates is empty or the targeted depth is reached, the node LP is
 *  solved, and the old candidates are replaced by the new LP candidates.
 *
 *  @see heur_guideddiving.c for an example implementation of a dive set controlling the diving algorithm.
 *
 *  @see the parameter @p heuristics/startdivefrac to determine the fraction of candidates that should be dived on at the
 *       beginning. Setting this parameter to 0.0 will result in an LP solved after every candidate selection.
 *
 *  @note the fraction of candidate variables is subject to change during solving. It is decreased by a factor of
 *        2 every time the algorithm could not dive half as deep as desired. However, if it succeeded, the fraction
 *        is multiplied by a factor of 1.1.
 *
 *  @note the node from where the algorithm is called is checked for a basic LP solution. If the solution
 *        is non-basic, e.g., when barrier without crossover is used, the method returns without performing a dive.
 *
 *  @note currently, when multiple diving heuristics call this method and solve an LP at the same node, only the first
 *        call will be executed, @see SCIPgetLastDiveNode()
 *
 *  @todo generalize method to work correctly with pseudo or external branching/diving candidates
 */
#include "scip/struct_heur.h"
#include "scip/struct_stat.h"

SCIP_RETCODE SCIPperformGenericDivingAlgorithm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< settings for diving */
   SCIP_SOL*             worksol,            /**< non-NULL working solution */
   SCIP_HEUR*            heur,               /**< the calling primal heuristic */
   SCIP_RESULT*          result,             /**< SCIP result pointer */
   SCIP_Bool             nodeinfeasible      /**< is the current node known to be infeasible? */
   )
{
   SCIP_CONSHDLR* indconshdlr;               /* constraint handler for indicator constraints */
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real objval;
   SCIP_Real nextcandsol;
   SCIP_Real ubquot;
   SCIP_Real avgquot;
   SCIP_Longint ncalls;
   SCIP_Longint oldsolsuccess;
   SCIP_Longint nlpiterations;
   SCIP_Longint maxnlpiterations;
   int ndivecands;
   int startndivecands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int nextcand;
   int targetdepth;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_Bool backtracked;
   SCIP_Bool backtrack;

   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   int oldnsolsfound;
   int oldnbestsolsfound;

   SCIP_VAR** indcands = NULL;
   SCIP_Real* indcandssol = NULL;
   SCIP_Real* indcandsfrac = NULL;
   int nindconss = 0;
   int nindcands;

   assert(scip != NULL);
   assert(result != NULL);
   assert(worksol != NULL);

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < SCIPdivesetGetMinRelDepth(diveset) * maxdepth || depth > SCIPdivesetGetMaxRelDepth(diveset) * maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNNodeLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   oldsolsuccess = SCIPdivesetGetSolSuccess(diveset);

   /*todo another factor of 10, REALLY? */
   maxnlpiterations = (SCIP_Longint)((1.0 + 10*(oldsolsuccess+1.0)/(ncalls+1.0)) * SCIPdivesetGetMaxLPIterQuot(diveset) * nlpiterations);
   maxnlpiterations += SCIPdivesetGetMaxLPIterOffset(diveset);


   /* don't try to dive, if we took too many LP iterations during diving */
   if( SCIPdivesetGetNLPIterations(diveset) >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   if( SCIPdivesetGetNLPIterations(diveset) + MINLPITER > maxnlpiterations )
      maxnlpiterations = SCIPdivesetGetNLPIterations(diveset) + MINLPITER;

   /* if indicator variables are present, add them to the set of diving candidates */
   /* todo maybe store those constraints once and not every time */
   indconshdlr = SCIPfindConshdlr(scip, "indicator");
   if( indconshdlr != NULL )
   {
      nindconss = SCIPconshdlrGetNConss(indconshdlr);
      if( nindconss > 0 )
      {
         /* get storage for candidate variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &indcands, nindconss) );
         SCIP_CALL( SCIPallocBufferArray(scip, &indcandssol, nindconss) );
         SCIP_CALL( SCIPallocBufferArray(scip, &indcandsfrac, nindconss) );
      }
   }

   /* get all current diving candidates */
   SCIP_CALL( getDivingCandidates(scip, diveset, indconshdlr, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands,
         indcands, indcandssol, indcandsfrac, &nindcands) );
   ndivecands = nlpcands + nindcands;

   /* don't try to dive, if there are no diving candidates */
   if( ndivecands == 0 )
   {
      if( nindconss > 0 )
      {
         /* free storage for indicator variables */
         SCIPfreeBufferArray(scip, &indcands);
         SCIPfreeBufferArray(scip, &indcandssol);
         SCIPfreeBufferArray(scip, &indcandsfrac);
      }

      return SCIP_OKAY;
   }

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      ubquot = SCIPdivesetGetUbQuotNoSol(diveset);
      avgquot = SCIPdivesetGetAvgQuotNoSol(diveset);
   }
   else
   {
      ubquot = SCIPdivesetGetUbQuot(diveset);
      avgquot = SCIPdivesetGetAvgQuot(diveset);
   }

   if( ubquot > 0.0 )
      searchubbound = SCIPgetLowerbound(scip) + ubquot * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
   else
      searchubbound = SCIPinfinity(scip);

   if( avgquot > 0.0 )
      searchavgbound = SCIPgetLowerbound(scip) + avgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
   else
      searchavgbound = SCIPinfinity(scip);

   searchbound = MIN(searchubbound, searchavgbound);

   if( SCIPisObjIntegral(scip) )
      searchbound = SCIPceil(scip, searchbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;

   *result = SCIP_DIDNOTFIND;

   oldnsolsfound = SCIPgetNSolsFound(scip);
   oldnbestsolsfound = SCIPgetNBestSolsFound(scip);

   /* start probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   /* get LP objective value */
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetLPObjval(scip);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing %s heuristic: depth=%d, %d fractionals, dualbound=%g, avgbound=%g, cutoffbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPheurGetName(heur), SCIPgetDepth(scip), ndivecands, SCIPgetDualbound(scip), SCIPgetAvgDualbound(scip),
      SCIPretransformObj(scip, SCIPgetCutoffbound(scip)), SCIPretransformObj(scip, searchbound));

   lperror = FALSE;
   cutoff = FALSE;
   divedepth = 0;
   startndivecands = ndivecands;

   /* LP loop; every time a new LP was solved, conditions are checked
    * dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   while( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && ndivecands > 0
      && (divedepth < 10
         || ndivecands <= startndivecands - divedepth / 2
         || (divedepth < maxdivedepth && SCIPdivesetGetNLPIterations(diveset) < maxnlpiterations && objval < searchbound))
         && !SCIPisStopped(scip) )
   {
      SCIP_VAR** divecands;
      SCIP_Real* divecandssol;
      SCIP_Real* divecandsfrac;
      SCIP_Bool* candsroundup;  /* stores for every candidate if it should be rounded up or down */
      SCIP_Bool nextcandroundup;
      SCIP_VAR* nextcandvar;
      int nbacktracks;
      int startdepth;
      int ncandstofix;
      int maxnbacktracks;
      SCIP_Bool allroundable;
      int c;

      /* determine the target depth (depth where the next LP should be solved) */
      startdepth = divedepth;
#if 0
      ncandstofix = MIN(ndivecands, maxdivedepth - divedepth);
      ncandstofix = (int)SCIPceil(scip, ncandstofix * SCIPdivesetGetTargetdepthfrac(diveset));
      ncandstofix = MAX(ncandstofix, 1);
#endif
      if( ndivecands <= SCIPgetDiveLPSolveFreq(scip) || SCIPgetDiveLPSolveFreq(scip) == 0 )
         ncandstofix = ndivecands;
      else
         ncandstofix = SCIPgetDiveLPSolveFreq(scip);

      targetdepth = divedepth + ncandstofix;

      SCIPdebugMessage("%s heuristic continues diving at depth %d, %d candidates left, %d candidates to fix\n",
         SCIPheurGetName(heur), startdepth, ndivecands, ncandstofix);

      /* loop over candidates and determine if they are roundable */
      allroundable = TRUE;
      c = 0;
      while( allroundable && c < nlpcands )
      {
         if( SCIPvarMayRoundDown(lpcands[c]) || SCIPvarMayRoundUp(lpcands[c]) )
            allroundable = TRUE;
         else
            allroundable = FALSE;
         ++c;
      }

      /* if all candidates are roundable, try to round the solution */
      if( allroundable )
      {
         SCIP_Bool success = FALSE;

         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIPlinkLPSol(scip, worksol) );
         SCIP_CALL( SCIProundSol(scip, worksol, &success) );

         /* succesfully rounded solutions are tried for primal feasibility */
         if( success )
         {
            SCIP_Bool changed = FALSE;
            SCIPdebugMessage("%s found roundable primal solution: obj=%g\n", SCIPheurGetName(heur), SCIPgetSolOrigObj(scip, worksol));

            /* adjust indicator constraints */
            if( indconshdlr != NULL )
            {
               SCIP_CALL( SCIPmakeIndicatorsFeasible(scip, indconshdlr, worksol, &changed) );
            }

            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, &success) );

            /* check, if solution was feasible and good enough */
            if( success )
            {
               SCIPdebugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      nextcand = -1;
      divecands = NULL;
      divecandssol = NULL;
      divecandsfrac = NULL;
      candsroundup = NULL;
      nextcandvar = NULL;
      nextcandsol = SCIP_INVALID;
      nextcandroundup = FALSE;

      /* sort the candidates in nondecreasing score order */
      if( ncandstofix > 1 )
      {
         SCIP_Real* candscores;
         int s;

         SCIP_CALL( SCIPallocBufferArray(scip, &divecands, ndivecands) );
         SCIP_CALL( SCIPallocBufferArray(scip, &divecandssol, ndivecands) );
         SCIP_CALL( SCIPallocBufferArray(scip, &divecandsfrac, ndivecands) );

         /* copy LP branching candidates */
         if( nlpcands > 0 )
         {
            BMScopyMemoryArray(divecands, lpcands, nlpcands);
            BMScopyMemoryArray(divecandssol, lpcandssol, nlpcands);
            BMScopyMemoryArray(divecandsfrac, lpcandsfrac, nlpcands);
         }

         /* copy indicator variables */
         if( nindcands > 0 )
         {
            BMScopyMemoryArray(&(divecands[nlpcands]), indcands, nindcands); /*lint !e866*/
            BMScopyMemoryArray(&(divecandssol[nlpcands]), indcandssol, nindcands); /*lint !e866*/
            BMScopyMemoryArray(&(divecandsfrac[nlpcands]), indcandsfrac, nindcands); /*lint !e866*/
         }

         SCIP_CALL( SCIPallocBufferArray(scip, &candscores, ndivecands) );
         SCIP_CALL( SCIPallocBufferArray(scip, &candsroundup, ndivecands) );

         /* store variable scores in temporary array */
         /* todo maybe move the for loop into the callback of the heuristics */
         for( s = 0; s < ndivecands; ++s )
         {
            SCIP_CALL( SCIPgetDivesetScore(scip, diveset, divecands[s], divecandssol[s], divecandsfrac[s],
                  &(candscores[s]), &(candsroundup[s])) );
         }

         SCIPsortRealRealRealBoolPtr(candscores, divecandssol, divecandsfrac, candsroundup, (void **)divecands, ndivecands);

         SCIPfreeBufferArray(scip, &candscores);

         nextcand = 0;
      }
      else
      {
         /* search for the candidate which minimizes the score */
         SCIP_Real minscore;

         /* find LP candidate with minimum score */
         minscore = SCIPinfinity(scip);
         for( c = 0; c < nlpcands; ++c )
         {
            SCIP_Real score;
            SCIP_Bool roundup;

            SCIP_CALL( SCIPgetDivesetScore(scip, diveset, lpcands[c], lpcandssol[c], lpcandsfrac[c], &score, &roundup) );

            /* new minimum found (or no variable has been found previously and score is infinity) */
            if( score < minscore || nextcandvar == NULL )
            {
               assert( score < minscore || SCIPisInfinity(scip, score) );
               minscore = score;
               nextcandvar = lpcands[c];
               nextcandsol = lpcandssol[c];
               nextcandroundup = roundup;
            }
         }

         /* search indicator candidates for even better scores */
         for( c = 0; c < nindcands; ++c )
         {
            SCIP_Real score;
            SCIP_Bool roundup;

            assert( indcands != NULL && indcandssol != NULL && indcandsfrac != NULL );  /* for lint */
            SCIP_CALL( SCIPgetDivesetScore(scip, diveset, indcands[c], indcandssol[c], indcandsfrac[c], &score, &roundup) );

            /* new minimum found (or no variable has been found previously and score is infinity) */
            if( score < minscore || nextcandvar == NULL )
            {
               assert( score < minscore || SCIPisInfinity(scip, score) );
               minscore = score;
               nextcandvar = indcands[c];
               nextcandsol = indcandssol[c];
               nextcandroundup = roundup;
            }
         }
      }

      /* limit the number of allowed backtracks */
      nbacktracks = 0;
      maxnbacktracks = (1 + (ncandstofix / 2));

      /* start propagating candidate variables
       *   - until the desired targetdepth is reached,
       *   - or there is no further candidate variable left because of intermediate bound changes,
       *   - or a cutoff is detected
       */
      do
      {
         assert(nextcand >= 0 || ncandstofix == 1);
         assert(nextcand >= 0 || nextcandvar != NULL);

         /* a next cand of -1 means that the best candidate was already selected prior to the loop */
         if( nextcand >= 0 )
         {
            assert( divecands != NULL && divecandssol != NULL && candsroundup != NULL ); /* for lint */
            nextcandvar = divecands[nextcand];
            nextcandsol = divecandssol[nextcand];
            nextcandroundup = candsroundup[nextcand];
         }

         /* dive deeper into the tree */
         SCIP_CALL( SCIPnewProbingNode(scip) );
         divedepth++;

         backtracked = FALSE;
         do
         {
            backtrack = FALSE;

            /* if the variable is already fixed or if the solution value is outside the domain, numerical troubles may have
             * occured or variable was fixed by propagation while backtracking => Abort diving!
             */
            if( SCIPvarGetLbLocal(nextcandvar) >= SCIPvarGetUbLocal(nextcandvar) - 0.5 )
            {
               SCIPdebugMessage("Selected variable <%s> already fixed to [%g,%g] (solval: %.9f), diving aborted \n",
                     SCIPvarGetName(nextcandvar), SCIPvarGetLbLocal(nextcandvar), SCIPvarGetUbLocal(nextcandvar), nextcandsol);
               cutoff = TRUE;
               break;
            }
            if( SCIPisFeasLT(scip, nextcandsol, SCIPvarGetLbLocal(nextcandvar)) || SCIPisFeasGT(scip, nextcandsol, SCIPvarGetUbLocal(nextcandvar)) )
            {
               SCIPdebugMessage("selected variable's <%s> solution value is outside the domain [%g,%g] (solval: %.9f), diving aborted\n",
                     SCIPvarGetName(nextcandvar), SCIPvarGetLbLocal(nextcandvar), SCIPvarGetUbLocal(nextcandvar), nextcandsol);
               cutoff = TRUE;
               break;
            }

            /* apply rounding of best candidate */
            if( nextcandroundup == !backtracked )
            {
               SCIP_Real value = SCIPfeasCeil(scip, nextcandsol);
               if( SCIPisFeasIntegral(scip, nextcandsol) )
               {
                  /* only indicator variables can have integral solution value */
                  assert(SCIPvarGetType(nextcandvar) == SCIP_VARTYPE_BINARY);
                  value = 1.0;
               }

               /* round variable up */
               SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
                  divedepth, maxdivedepth, SCIPdivesetGetNLPIterations(diveset), maxnlpiterations,
                  SCIPvarGetName(nextcandvar),
                  nextcandsol, SCIPvarGetLbLocal(nextcandvar), SCIPvarGetUbLocal(nextcandvar),
                  value, SCIPvarGetUbLocal(nextcandvar));

               SCIP_CALL( SCIPchgVarLbProbing(scip, nextcandvar, value) );
            }
            else
            {
               SCIP_Real value = SCIPfeasFloor(scip, nextcandsol);

               if( SCIPisFeasIntegral(scip, nextcandsol) )
               {
                  /* only indicator variables can have integral solution value */
                  assert(SCIPvarGetType(nextcandvar) == SCIP_VARTYPE_BINARY);
                  value = 0.0;
               }
               /* round variable down */
               SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
                  divedepth, maxdivedepth, SCIPdivesetGetNLPIterations(diveset), maxnlpiterations,
                  SCIPvarGetName(nextcandvar),
                  nextcandsol, SCIPvarGetLbLocal(nextcandvar), SCIPvarGetUbLocal(nextcandvar),
                  SCIPvarGetLbLocal(nextcandvar), value);

               SCIP_CALL( SCIPchgVarUbProbing(scip, nextcandvar, value) );
            }

            /* apply domain propagation */
            SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );

            /* perform backtracking if a cutoff was detected */
            if( cutoff && !backtracked && SCIPdivesetUseBacktrack(diveset) && nbacktracks < maxnbacktracks )
            {
               SCIPdebugMessage("  *** cutoff detected at level %d - backtracking\n", SCIPgetProbingDepth(scip));
               SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );
               SCIP_CALL( SCIPnewProbingNode(scip) );
               backtracked = TRUE;
               backtrack = TRUE;
               cutoff = FALSE;
               ++nbacktracks;
            }
            else
               backtrack = FALSE;
         }
         while( backtrack );

         if( !cutoff && targetdepth > divedepth )
         {
            assert(nextcand >= 0 && nextcand < ndivecands);
            assert(divecands != NULL);
            assert(divecandssol != NULL);

            /* we need to search for the next candidate in our list which was not previously fixed or whose LP solution
             * is not already infeasible
             */
            while( nextcand < ndivecands )
            {
               ++nextcand;
               if( nextcand == ndivecands || (SCIPvarGetLbLocal(divecands[nextcand]) >= SCIPvarGetUbLocal(divecands[nextcand]) - 0.5)
                  || SCIPisLT(scip, divecandssol[nextcand], SCIPvarGetLbLocal(divecands[nextcand]))
                  || SCIPisGT(scip, divecandssol[nextcand], SCIPvarGetUbLocal(divecands[nextcand])) )
               {
                  SCIPdebugMessage(" <%s> solution value is fixed/outside the domain [%g,%g] (solval: %.9f), variable is skipped\n",
                     SCIPvarGetName(divecands[nextcand]), SCIPvarGetLbLocal(divecands[nextcand]), SCIPvarGetUbLocal(divecands[nextcand]), divecandssol[nextcand]);
               }
               else
               {
                  /* proceed with candidate variable whose LP solution value was not already cut off by domain propagation */
                  break;
               }
            }
         }
      }
      while( !cutoff && targetdepth > divedepth && nextcand < ndivecands );

      assert(cutoff || nextcand == ndivecands || targetdepth == divedepth);

      /* resolve the diving LP */
      if( !cutoff )
      {
         int lpiterationlimit;
         SCIP_RETCODE retstat;

         nlpiterations = SCIPgetNLPIterations(scip);

         /* allow at least MINLPITER more iterations */
         lpiterationlimit = (int)(maxnlpiterations - SCIPdivesetGetNLPIterations(diveset));
         lpiterationlimit = MAX(lpiterationlimit, MINLPITER);

         retstat = SCIPsolveProbingLP(scip, lpiterationlimit, &lperror, &cutoff);
         /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
          * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
          */
#ifdef NDEBUG
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while solving LP in %s heuristic; LP solve terminated with code <%d>\n", SCIPheurGetName(heur), retstat);
         }
#else
         SCIP_CALL( retstat );
#endif

         if( lperror )
         {
            /* free stored buffer arrays */
            if( ncandstofix > 1 )
            {
               SCIPfreeBufferArray(scip, &candsroundup);
               SCIPfreeBufferArray(scip, &divecandsfrac);
               SCIPfreeBufferArray(scip, &divecandssol);
               SCIPfreeBufferArray(scip, &divecands);
            }

            break;
         }

         /* update iteration count */
         SCIPupdateDivesetLPStats(scip, diveset, SCIPgetNLPIterations(scip) - nlpiterations);

         /* get LP solution status, objective value, and fractional variables, that should be integral */
         objval = SCIPgetLPObjval(scip);
         lpsolstat = SCIPgetLPSolstat(scip);
         assert(cutoff || (lpsolstat != SCIP_LPSOLSTAT_OBJLIMIT && lpsolstat != SCIP_LPSOLSTAT_INFEASIBLE &&
               (lpsolstat != SCIP_LPSOLSTAT_OPTIMAL || SCIPisLT(scip, objval, SCIPgetCutoffbound(scip)))));
      }

      SCIPdebugMessage("   -> lpsolstat=%d, objval=%g/%g, nfrac=%d\n", lpsolstat, objval, searchbound, ndivecands);

      /* reward the diving setting by increasing the dive depth quotient for future purpose */
      if( !cutoff && (divedepth == targetdepth || nextcand == ndivecands ) )
      {
         SCIPdivesetSetTargetdepthfrac(diveset, 1.1 * SCIPdivesetGetTargetdepthfrac(diveset));
      }
      else if( cutoff )
      {
         /* if the last node was cut off, we terminate the heuristic */
         int reacheddepth = divedepth - startdepth - 1;
         assert(reacheddepth >= 0);

         /* evaluate how deep we went this time and be more conservative in the future, if possible
          * if not even half the targeted depth was reached, decrease the depth quotient
          */
         if( reacheddepth < (targetdepth - startdepth) / 2 )
         {
            SCIPdivesetSetTargetdepthfrac(diveset, 0.5 * SCIPdivesetGetTargetdepthfrac(diveset));
         }
      }

      /* free stored buffer arrays */
      if( ncandstofix > 1 )
      {
         SCIPfreeBufferArray(scip, &candsroundup);
         SCIPfreeBufferArray(scip, &divecandsfrac);
         SCIPfreeBufferArray(scip, &divecandssol);
         SCIPfreeBufferArray(scip, &divecands);
      }

      /* todo if only one candidate was fixed since last LP, use the LP Objective gain to update pseudo cost information */
      if( !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new diving candidate variables */
         SCIP_CALL( getDivingCandidates(scip, diveset, indconshdlr,
               &lpcands, &lpcandssol, &lpcandsfrac,
               &nlpcands, indcands, indcandssol, indcandsfrac, &nindcands) );

         ndivecands = nlpcands + nindcands;
      }
      else
         ndivecands = 0;
   }

   /* check if a solution has been found */
   if( ndivecands == 0 && !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Bool success;

      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, worksol) );
      SCIPdebugMessage("%s found primal solution: obj=%g\n", SCIPheurGetName(heur), SCIPgetSolOrigObj(scip, worksol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* free storage for indicator variables */
   if( nindconss > 0 )
   {
      SCIPfreeBufferArray(scip, &indcands);
      SCIPfreeBufferArray(scip, &indcandssol);
      SCIPfreeBufferArray(scip, &indcandsfrac);
   }

   SCIPupdateDivesetStats(scip, diveset, SCIPgetDepth(scip), 10 * (SCIPgetNBestSolsFound(scip) - oldnbestsolsfound) + SCIPgetNSolsFound(scip) - oldnsolsfound);

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );


   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") finished %s heuristic: %d fractionals, dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", objval=%g/%g, lpsolstat=%d, cutoff=%u\n",
      SCIPgetNNodes(scip), SCIPheurGetName(heur), ndivecands, divedepth, maxdivedepth, SCIPdivesetGetNLPIterations(diveset), maxnlpiterations,
      SCIPretransformObj(scip, objval), SCIPretransformObj(scip, searchbound), lpsolstat, cutoff);

   return SCIP_OKAY;
}
