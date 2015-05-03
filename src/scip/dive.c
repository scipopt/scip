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
/*#define SCIP_DEBUG*/
/**@file   scip.c
 * @brief  library methods for diving heuristics
 * @author Gregor Hendel
 *
 * @todo check all checkStage() calls, use bit flags instead of the SCIP_Bool parameters
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

#include "scip/pub_dive.h"
#include "pub_heur.h"
#include "scip/struct_heur.h"
#include "scip/struct_stat.h"

/* the indicator constraint handler is included for the diving algorithm SCIPperformGenericDivingAlgorithm() */
#include "scip/cons_indicator.h"

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */


/** solve probing LP */
static
SCIP_RETCODE solveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_Longint          maxnlpiterations,   /**< maximum number of allowed LP iterations */
   SCIP_Bool*            lperror,            /**< pointer to store if an internal LP error occurred */
   SCIP_Bool*            cutoff              /**< pointer to store whether the LP was infeasible */
   )
{
   int lpiterationlimit;
   SCIP_RETCODE retstat;
   SCIP_Longint nlpiterations;

   assert(lperror != NULL);
   assert(cutoff != NULL);

   nlpiterations = SCIPgetNLPIterations(scip);

   /* allow at least MINLPITER more iterations */
   lpiterationlimit = (int)(maxnlpiterations - SCIPdivesetGetNLPIterations(diveset));
   lpiterationlimit = MAX(lpiterationlimit, MINLPITER);

   retstat = SCIPsolveProbingLP(scip, lpiterationlimit, lperror, cutoff);

   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
#ifdef NDEBUG
   if( retstat != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while solving LP in %s diving heuristic; LP solve terminated with code <%d>.\n", SCIPdivesetGetName(diveset), retstat);
   }
#else
   SCIP_CALL( retstat );
#endif

   /* update iteration count */
   SCIPupdateDivesetLPStats(scip, diveset, SCIPgetNLPIterations(scip) - nlpiterations);

   return SCIP_OKAY;
}



/** performs a diving within the limits of the diveset parameters
 *
 *  This method performs a diving according to the settings defined by the diving settings @p diveset; Contrary to the
 *  name, SCIP enters probing mode (not diving mode) and dives along a path into the tree. Domain propagation
 *  is applied at every node in the tree, whereas probing LPs might be solved less frequently.
 *
 *  Starting from the current LP candidates, the algorithm determines a fraction of the candidates that should be
 *  branched on; if a single candidate should be fixed, the algorithm selects a candidate which minimizes the score
 *  defined by the @p diveset.  If more than one candidate should be selected, the candidates are sorted in
 *  non-decreasing order of their score.
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
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;

   SCIP_Real* vals;
   SCIP_VAR** previouscands;
   SCIP_Real* previousvals;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real nextcandsol;
   SCIP_Real ubquot;
   SCIP_Real avgquot;
   SCIP_Longint ncalls;
   SCIP_Longint oldsolsuccess;
   SCIP_Longint nlpiterations;
   SCIP_Longint maxnlpiterations;
   SCIP_Longint domreds;
   int startndivecands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int totalnbacktracks;
   int totalnprobingnodes;
   int lastlpdepth;
   int previouscandssize;

   SCIP_Bool success;
   SCIP_Bool enfosuccess;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_Bool backtracked;
   SCIP_Bool backtrack;

   int nlpcands;

   assert(scip != NULL);
   assert(result != NULL);
   assert(worksol != NULL);

   *result = SCIP_DELAYED;

   /* do not call heuristic in node that was already detected to be infeasible */
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
   ncalls = SCIPdivesetGetNCalls(diveset);
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

   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );

   /* don't try to dive, if there are no diving candidates */
   if( nlpcands == 0 )
      return SCIP_OKAY;

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

   /* start probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing %s heuristic: depth=%d, %d fractionals, dualbound=%g, avgbound=%g, cutoffbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPheurGetName(heur), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPgetAvgDualbound(scip),
      SCIPretransformObj(scip, SCIPgetCutoffbound(scip)), SCIPretransformObj(scip, searchbound));


   /* allocate buffer storage for previous candidates and their branching values for pseudo cost updates */
   previouscandssize = MAX(1, SCIPgetDiveLPSolveFreq(scip));
   SCIP_CALL( SCIPallocBufferArray(scip, &previouscands, previouscandssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &previousvals, previouscandssize) );

   /* keep some statistics */
   lperror = FALSE;
   cutoff = FALSE;
   lastlpdepth = -1;
   startndivecands = nlpcands;
   totalnbacktracks = 0;
   totalnprobingnodes = 0;

   /* link the working solution to the dive set */
   SCIPdivesetSetWorkSolution(diveset, worksol);

   /* allocate buffer array to store enforcement decision */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );
   enfosuccess = TRUE;

   /* LP loop; every time a new LP was solved, conditions are checked
    * dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   while( !lperror && !cutoff && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL && enfosuccess
      && (SCIPgetProbingDepth(scip) < 10
         || nlpcands <= startndivecands - SCIPgetProbingDepth(scip) / 2
         || (SCIPgetProbingDepth(scip) < maxdivedepth && SCIPdivesetGetNLPIterations(diveset) < maxnlpiterations && SCIPgetLPObjval(scip) < searchbound))
         && !SCIPisStopped(scip) )
   {
      SCIP_VAR* nextcandvar;
      SCIP_Real lastlpobjval;
      SCIP_Bool nextcandroundup;
      SCIP_Bool allroundable;
      SCIP_Bool infeasible;
      int c;

      /* remember the last LP depth  */
      assert(lastlpdepth < SCIPgetProbingDepth(scip));
      lastlpdepth = SCIPgetProbingDepth(scip);
      domreds = 0;

      SCIPdebugMessage("%s heuristic continues diving at depth %d, %d candidates left\n",
         SCIPdivesetGetName(diveset), lastlpdepth, nlpcands);


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
         success = FALSE;

         /* working solution must be linked to LP solution */
         SCIP_CALL( SCIPlinkLPSol(scip, worksol) );
         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIProundSol(scip, worksol, &success) );

         /* succesfully rounded solutions are tried for primal feasibility */
         if( success )
         {
            SCIP_Bool changed = FALSE;
            SCIPdebugMessage("%s found roundable primal solution: obj=%g\n", SCIPdivesetGetName(diveset), SCIPgetSolOrigObj(scip, worksol));

            /* adjust indicator constraints */
            if( indconshdlr != NULL )
            {
               SCIP_CALL( SCIPmakeIndicatorsFeasible(scip, indconshdlr, worksol, &changed) );
            }

            success = FALSE;
            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, &success) );

            /* check, if solution was feasible and good enough */
            if( success )
            {
               SCIPdebugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;

               /* the rounded solution found above led to a cutoff of the node LP solution */
               if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OBJLIMIT )
               {
                  cutoff = TRUE;
                  break;
               }
            }
         }
      }

      /* working solution must be linked to LP solution */
      assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);
      lastlpobjval = SCIPgetLPObjval(scip);
      SCIP_CALL( SCIPlinkLPSol(scip, worksol) );


      nextcandvar = NULL;
      nextcandsol = SCIP_INVALID;
      nextcandroundup = FALSE;
      /* choose next candidate variable, usually done at the end of the previous iteration */
      enfosuccess = FALSE;
      vals[0] = vals[1] = SCIP_INVALID;

      SCIP_CALL( SCIPenforceDiveSolution(scip, diveset, worksol, &nextcandvar, vals, &enfosuccess, &infeasible) );

      /* if we did not succeed finding an enforcement, the solution is potentially feasible and we break immediately */
      if( ! enfosuccess )
         break;

      /* start propagating candidate variables
       *   - until the desired targetdepth is reached,
       *   - or there is no further candidate variable left because of intermediate bound changes,
       *   - or a cutoff is detected
       */
      do
      {
         SCIP_Longint localdomreds;

         /* ensure that a new candidate was successfully determined (usually at the end of the previous loop iteration) */
         assert(enfosuccess);
         assert(nextcandvar != NULL);
         assert(vals[1] != SCIP_INVALID);
         assert(vals[0] != SCIP_INVALID);
         nextcandsol = SCIPgetSolVal(scip, worksol, nextcandvar);

         /* treat indicator variables specially, they might have integral solution values */
         nextcandroundup = ((SCIPvarIsBinary(nextcandvar) && vals[0] > 0.5) || vals[0] > nextcandsol);

         backtracked = FALSE;
         do
         {
            assert(nextcandvar != NULL);

            /* dive deeper into the tree */
            SCIP_CALL( SCIPnewProbingNode(scip) );
            ++totalnprobingnodes;

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
               SCIP_Real value = vals[backtracked ? 1 : 0];
               if( SCIPisFeasIntegral(scip, nextcandsol) )
               {
                  /* only indicator variables can have integral solution value */
                  assert(SCIPvarGetType(nextcandvar) == SCIP_VARTYPE_BINARY);
                  value = 1.0;
               }

               /* round variable up */
               SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
                  SCIPgetProbingDepth(scip), maxdivedepth, SCIPdivesetGetNLPIterations(diveset), maxnlpiterations,
                  SCIPvarGetName(nextcandvar),
                  nextcandsol, SCIPvarGetLbLocal(nextcandvar), SCIPvarGetUbLocal(nextcandvar),
                  value, SCIPvarGetUbLocal(nextcandvar));

               SCIP_CALL( SCIPchgVarLbProbing(scip, nextcandvar, value) );
            }
            else
            {
               SCIP_Real value = vals[backtracked ? 1 : 0];

               if( SCIPisFeasIntegral(scip, nextcandsol) )
               {
                  /* only indicator variables can have integral solution value */
                  assert(SCIPvarGetType(nextcandvar) == SCIP_VARTYPE_BINARY);
                  value = 0.0;
               }
               /* round variable down */
               SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
                  SCIPgetProbingDepth(scip), maxdivedepth, SCIPdivesetGetNLPIterations(diveset), maxnlpiterations,
                  SCIPvarGetName(nextcandvar),
                  nextcandsol, SCIPvarGetLbLocal(nextcandvar), SCIPvarGetUbLocal(nextcandvar),
                  SCIPvarGetLbLocal(nextcandvar), value);

               SCIP_CALL( SCIPchgVarUbProbing(scip, nextcandvar, value) );
            }

            /* apply domain propagation */
            localdomreds = 0;
            SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, &localdomreds) );

            /* resolve the diving LP if the diving resolve frequency is reached or a sufficient number of intermediate bound changes
             * was reached
             */
            if( !cutoff && ((SCIPgetDiveLPSolveFreq(scip) > 0 && ((SCIPgetProbingDepth(scip) - lastlpdepth) % SCIPgetDiveLPSolveFreq(scip)) == 0)
                  || (domreds + localdomreds > SCIPgetDiveLPResolveDomChgQuot(scip) * SCIPgetNVars(scip))) )
            {
               SCIP_CALL( solveLP(scip, diveset, maxnlpiterations, &lperror, &cutoff) );

               /* lp errors lead to early termination */
               if( lperror )
               {
                  cutoff = TRUE;
                  break;
               }
            }

            /* perform backtracking if a cutoff was detected */
            if( cutoff && !backtracked && SCIPdivesetUseBacktrack(diveset) )
            {
               SCIPdebugMessage("  *** cutoff detected at level %d - backtracking\n", SCIPgetProbingDepth(scip));
               SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );
               ++totalnbacktracks;
               backtracked = TRUE;
               backtrack = TRUE;
               cutoff = FALSE;
            }
            else
               backtrack = FALSE;
         }
         while( backtrack );

         /* we add the domain reductions from the last evaluated node */
         domreds += localdomreds;

         /* store candidate for pseudo cost update and choose next candidate only if no cutoff was detected */
         if( ! cutoff )
         {
            int insertidx = SCIPgetProbingDepth(scip) - lastlpdepth - 1;
            assert(SCIPgetProbingDepth(scip) > 0);

            /* extend array in case of a dynamic, domain change based LP resolve strategy */
            if( insertidx >= previouscandssize )
            {
               previouscandssize *= 2;
               SCIP_CALL( SCIPreallocBufferArray(scip, &previouscands, previouscandssize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &previousvals, previouscandssize) );
            }
            assert(previouscandssize > insertidx);

            /* store candidate for pseudo cost update */
            previouscands[insertidx] = nextcandvar;
            previousvals[insertidx] = vals[backtracked ? 1 : 0];

            /* choose next candidate variable and resolve the LP if none is found. */
            if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_NOTSOLVED )
            {
               assert(SCIPgetProbingDepth(scip) > lastlpdepth);
               enfosuccess = FALSE;
               nextcandvar = NULL;
               vals[0] = vals[1] = SCIP_INVALID;

               SCIP_CALL( SCIPenforceDiveSolution(scip, diveset, worksol, &nextcandvar, vals, &enfosuccess, &infeasible) );

               /* in case of an unsuccesful candidate search, we solve the node LP */
               if( !enfosuccess )
               {
                  SCIP_CALL( solveLP(scip, diveset, maxnlpiterations, &lperror, &cutoff) );

                  /* check for an LP error and terminate in this case, cutoffs lead to termination anyway */
                  if( lperror )
                     cutoff = TRUE;

                  /* enfosuccess must be set to TRUE for entering the main LP loop again */
                  enfosuccess = TRUE;
               }
            }
         }
      }
      while( !cutoff && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_NOTSOLVED );

      assert(cutoff || lperror || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_NOTSOLVED);

      assert(cutoff || (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OBJLIMIT && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_INFEASIBLE &&
            (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL || SCIPisLT(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))));


      /* check new LP candidates and use the LP Objective gain to update pseudo cost information */
      if( ! cutoff && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         int v;
         SCIP_Real gain;

         SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );

         /* distribute the gain equally over all variables that we rounded since the last LP */
         gain = SCIPgetLPObjval(scip) - lastlpobjval;
         gain = MAX(gain, 0.0);
         gain /= (1.0 * (SCIPgetProbingDepth(scip) - lastlpdepth));

         /* loop over previously fixed candidates and share gain improvement */
         for( v = 0; v < (SCIPgetProbingDepth(scip) - lastlpdepth); ++v )
         {
            SCIP_VAR* cand = previouscands[v];
            SCIP_Real val = previousvals[v];
            SCIP_Real solval = SCIPgetSolVal(scip, worksol, cand);

            /* it may happen that a variable had an integral solution value beforehand, e.g., for indicator variables */
            if( !SCIPisZero(scip, val - solval) )
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, cand, val - solval, gain, 1.0) );
            }
         }
      }
      else
         nlpcands = 0;
      SCIPdebugMessage("   -> lpsolstat=%d, objval=%g/%g, nfrac=%d\n", SCIPgetLPSolstat(scip), SCIPgetLPObjval(scip), searchbound, nlpcands);
   }

   success = FALSE;
   /* check if a solution has been found */
   if( !enfosuccess && !lperror && !cutoff && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, worksol) );
      SCIPdebugMessage("%s found primal solution: obj=%g\n", SCIPdivesetGetName(diveset), SCIPgetSolOrigObj(scip, worksol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   SCIPupdateDivesetStats(scip, diveset, totalnprobingnodes, totalnbacktracks, success);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") finished %s heuristic: %d fractionals, dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", objval=%g/%g, lpsolstat=%d, cutoff=%u\n",
      SCIPgetNNodes(scip), SCIPdivesetGetName(diveset), nlpcands, SCIPgetProbingDepth(scip), maxdivedepth, SCIPdivesetGetNLPIterations(diveset), maxnlpiterations,
      SCIPretransformObj(scip, SCIPgetLPObjval(scip)), SCIPretransformObj(scip, searchbound), SCIPgetLPSolstat(scip), cutoff);

   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   SCIPdivesetSetWorkSolution(diveset, NULL);

   SCIPfreeBufferArray(scip, &previousvals);
   SCIPfreeBufferArray(scip, &previouscands);
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}
