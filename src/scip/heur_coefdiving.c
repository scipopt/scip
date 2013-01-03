/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_coefdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the matrix coefficients
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 *
 * Indicator constraints are taken into account if present.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_coefdiving.h"
#include "scip/cons_indicator.h"


#define HEUR_NAME             "coefdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the matrix coefficients"
#define HEUR_DISPCHAR         'c'
#define HEUR_PRIORITY         -1001000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          1
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.05 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int                   maxlpiterofs;       /**< additional number of allowed LP iterations */
   SCIP_Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Bool             backtrack;          /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Longint          nlpiterations;      /**< LP iterations used in this heuristic */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */
   SCIP_CONSHDLR*        indconshdlr;        /**< indicator constraint handler (or NULL) */
};


/*
 * local methods
 */


/** get indicator candidate variables */
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

   assert( scip != NULL );
   assert( indconss != NULL );
   assert( indcands != NULL );
   assert( nindcands != NULL );
   assert( indcandssol != NULL );
   assert( indcandfrac != NULL );

   *nindcands = 0;
   for (c = 0; c < nindconss; ++c)
   {
      /* check whether constraint is violated */
      if ( SCIPisViolatedIndicator(scip, indconss[c], NULL) )
      {
         binvar = SCIPgetBinaryVarIndicator(indconss[c]);
         val = SCIPgetSolVal(scip, NULL, binvar);

         /* fractional indicator variables are treated by lpcands */
         if ( SCIPisFeasIntegral(scip, val) )
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

/** choose best candidate variable */
static
SCIP_RETCODE getBestCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            cands,              /**< candidate variables */
   SCIP_Real*            candssol,           /**< solution values of candidates */
   SCIP_Real*            candsfrac,          /**< fractional solution values of candidates */
   int                   ncands,             /**< number of candidates */
   int*                  bestcand,           /**< bestcandidate */
   int*                  bestnviolrows,      /**< number of violated rows for best candidate */
   SCIP_Real*            bestcandsol,        /**< solution of best candidate */
   SCIP_Real*            bestcandfrac,       /**< fractionality of best candidate */
   SCIP_Bool*            bestcandmayrounddown,/**< whether best candidate may be rounded down */
   SCIP_Bool*            bestcandmayroundup, /**< whether best candidate may be rounded down */
   SCIP_Bool*            bestcandroundup     /**< whether the best candidate should be rounded up */
   )
{
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   SCIP_Bool roundup;
   SCIP_Real frac;
   SCIP_VAR* var;
   int nlocksdown;
   int nlocksup;
   int nviolrows;
   int c;

   assert( cands != NULL );
   assert( candsfrac != NULL );
   assert( candssol != NULL );
   assert( bestcand != NULL );
   assert( bestnviolrows != NULL );
   assert( bestcandfrac != NULL );
   assert( bestcandsol != NULL );
   assert( bestcandfrac != NULL );
   assert( bestcandmayroundup != NULL );
   assert( bestcandmayrounddown != NULL );
   assert( bestcandroundup != NULL );

   for( c = 0; c < ncands; ++c )
   {
      var = cands[c];
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      frac = candsfrac[c];
      if( mayrounddown || mayroundup )
      {
         /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
         if( *bestcandmayrounddown || *bestcandmayroundup )
         {
            /* choose rounding direction:
             * - if variable may be rounded in both directions, round corresponding to the fractionality
             * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
             *   the current fractional solution
             */
            if( mayrounddown && mayroundup )
               roundup = (frac > 0.5);
            else
               roundup = mayrounddown;

            if( roundup )
            {
               frac = 1.0 - frac;
               nviolrows = SCIPvarGetNLocksUp(var);
            }
            else
               nviolrows = SCIPvarGetNLocksDown(var);

            /* penalize too small fractions */
            if( frac < 0.01 )
               nviolrows *= 100;

            /* prefer decisions on binary variables */
            if( !SCIPvarIsBinary(var) )
               nviolrows *= 100;

            /* check, if candidate is new best candidate */
            assert( (0.0 < frac && frac < 1.0) || SCIPvarIsBinary(var) );
            if( nviolrows + frac < *bestnviolrows + *bestcandfrac )
            {
               *bestcand = c;
               *bestnviolrows = nviolrows;
               *bestcandsol = candssol[c];
               *bestcandfrac = frac;
               *bestcandmayrounddown = mayrounddown;
               *bestcandmayroundup = mayroundup;
               *bestcandroundup = roundup;
            }
         }
      }
      else
      {
         /* the candidate may not be rounded */
         nlocksdown = SCIPvarGetNLocksDown(var);
         nlocksup = SCIPvarGetNLocksUp(var);
         roundup = (nlocksdown > nlocksup || (nlocksdown == nlocksup && frac > 0.5));
         if( roundup )
         {
            nviolrows = nlocksup;
            frac = 1.0 - frac;
         }
         else
            nviolrows = nlocksdown;

         /* penalize too small fractions */
         if( frac < 0.01 )
            nviolrows *= 100;

         /* prefer decisions on binary variables */
         if( !SCIPvarIsBinary(var) )
            nviolrows *= 100;

         /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
         assert( (0.0 < frac && frac < 1.0) || SCIPvarIsBinary(var) );
         if( *bestcandmayrounddown || *bestcandmayroundup || nviolrows + frac < *bestnviolrows + *bestcandfrac )
         {
            *bestcand = c;
            *bestnviolrows = nviolrows;
            *bestcandsol = candssol[c];
            *bestcandfrac = frac;
            *bestcandmayrounddown = FALSE;
            *bestcandmayroundup = FALSE;
            *bestcandroundup = roundup;
         }
         assert( *bestcandfrac < SCIP_INVALID );
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyCoefdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeHeurCoefdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeCoefdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitCoefdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize data */
   heurdata->nlpiterations = 0;
   heurdata->nsuccess = 0;

   /* get indicator constraint hanlder */
   heurdata->indconshdlr = SCIPfindConshdlr(scip, "indicator");

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitCoefdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecCoefdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_CONS** indconss;
   SCIP_VAR** indcands;
   SCIP_VAR** lpcands;
   SCIP_VAR* bestcandvar;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   SCIP_Real* indcandssol;
   SCIP_Real* indcandfrac;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real objval;
   SCIP_Real oldobjval;
   SCIP_Real bestcandsol;
   SCIP_Real bestcandfrac;   
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   SCIP_Bool bestcandroundup;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_Bool backtracked;
   SCIP_Longint ncalls;
   SCIP_Longint nsolsfound;
   SCIP_Longint nlpiterations;
   SCIP_Longint maxnlpiterations;
   int nindconss;
   int nlpcands;
   int nindcands;
   int startnlpcands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int bestnviolrows;
   int bestindcand;
   int bestlpcand;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < heurdata->minreldepth*maxdepth || depth > heurdata->maxreldepth*maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNNodeLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = 10*SCIPheurGetNBestSolsFound(heur) + heurdata->nsuccess;
   maxnlpiterations = (SCIP_Longint)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * nlpiterations);
   maxnlpiterations += heurdata->maxlpiterofs;

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + MINLPITER);

   /* get fractional variables that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );

   /* get indicator variable candidates */
   nindconss = 0;
   indconss = NULL;
   nindcands = 0;
   indcands = NULL;
   indcandssol = NULL;
   indcandfrac = NULL;
   if ( heurdata->indconshdlr != NULL )
   {
      indconss = SCIPconshdlrGetConss(heurdata->indconshdlr);
      nindconss = SCIPconshdlrGetNConss(heurdata->indconshdlr);

      if ( nindconss > 0 )
      {
         /* get storage for candidate variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &indcands, nindconss) );
         SCIP_CALL( SCIPallocBufferArray(scip, &indcandssol, nindconss) );
         SCIP_CALL( SCIPallocBufferArray(scip, &indcandfrac, nindconss) );

         /* get indicator canditates */
         SCIP_CALL( getIndCandVars(scip, indconss, nindconss, indcands, indcandssol, indcandfrac, &nindcands) );
      }
   }

   /* don't try to dive, if there are no fractional variables and no indicator candidates */
   if( nlpcands == 0 && nindcands == 0 )
   {
      SCIPfreeBufferArrayNull(scip, &indcandfrac);
      SCIPfreeBufferArrayNull(scip, &indcandssol);
      SCIPfreeBufferArrayNull(scip, &indcands);
      return SCIP_OKAY;
   }

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      if( heurdata->maxdiveubquotnosol > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquotnosol * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquotnosol > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquotnosol * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   else
   {
      if( heurdata->maxdiveubquot > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquot * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquot > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   searchbound = MIN(searchubbound, searchavgbound);
   if( SCIPisObjIntegral(scip) )
      searchbound = SCIPceil(scip, searchbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;

   *result = SCIP_DIDNOTFIND;

   /* start diving */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   /* get LP objective value */
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetLPObjval(scip);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing coefdiving heuristic: depth=%d, %d fractionals, dualbound=%g, avgbound=%g, cutoffbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPgetAvgDualbound(scip),
      SCIPretransformObj(scip, SCIPgetCutoffbound(scip)), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   lperror = FALSE;
   cutoff = FALSE;
   divedepth = 0;
   bestcandmayrounddown = FALSE;
   bestcandmayroundup = FALSE;
   startnlpcands = nlpcands + nindcands;
   while( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations && objval < searchbound))
      && !SCIPisStopped(scip) )
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      divedepth++;

      /* choose variable fixing:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, round variable with least number of locks in corresponding direction
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - round variable with least number of locks in opposite of its feasible rounding direction
       */
      bestlpcand = -1;
      bestindcand = -1;
      bestcandvar = NULL;
      bestnviolrows = INT_MAX;
      bestcandsol = SCIP_INVALID;
      bestcandfrac = SCIP_INVALID;
      bestcandmayrounddown = TRUE;
      bestcandmayroundup = TRUE;
      bestcandroundup = FALSE;

      /* get best lp candidate */
      if ( nlpcands > 0 )
      {
         SCIP_CALL( getBestCandidate(scip, lpcands, lpcandssol, lpcandsfrac, nlpcands, &bestlpcand, &bestnviolrows, &bestcandsol, &bestcandfrac,
               &bestcandmayrounddown, &bestcandmayroundup, &bestcandroundup) );
         bestcandvar = lpcands[bestlpcand];
         assert( bestlpcand >= 0 );
      }

      /* get best indicator candidate */
      if ( nindconss > 0 )
      {
         assert( indcands != NULL );
         assert( indcandssol != NULL );
         assert( indcandfrac != NULL );
         SCIP_CALL( getBestCandidate(scip, indcands, indcandssol, indcandfrac, nindcands, &bestindcand, &bestnviolrows, &bestcandsol, &bestcandfrac, 
               &bestcandmayrounddown, &bestcandmayroundup, &bestcandroundup) );
         if ( bestindcand >= 0 )
            bestcandvar = indcands[bestindcand];
      }

      /* if all candidates are roundable, try to round the solution */
      if( bestcandmayrounddown || bestcandmayroundup )
      {
         SCIP_Bool success;

         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
         SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

         if( success )
         {
            SCIPdebugMessage("coefdiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

            /* try to correct indicator constraints */
            if ( nindconss > 0 )
            {
               assert( heurdata->indconshdlr != NULL );
               SCIP_CALL( SCIPmakeIndicatorsFeasible(scip, heurdata->indconshdlr, heurdata->sol, &success) );
            }

            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, &success) );

            /* check, if solution was feasible and good enough */
            if( success )
            {
               SCIPdebugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      backtracked = FALSE;
      do
      {
         /* if the variable is already fixed or if the solution value is outside the domain, numerical troubles may have
          * occured or variable was fixed by propagation while backtracking => Abort diving!
          */
         if( SCIPvarGetLbLocal(bestcandvar) >= SCIPvarGetUbLocal(bestcandvar) - 0.5 )
         {
            SCIPdebugMessage("Selected variable <%s> already fixed to [%g,%g] (solval: %.9f), diving aborted \n",
               SCIPvarGetName(bestcandvar), SCIPvarGetLbLocal(bestcandvar), SCIPvarGetUbLocal(bestcandvar), bestcandsol);
            cutoff = TRUE;
            break;
         }
         if( SCIPisFeasLT(scip, bestcandsol, SCIPvarGetLbLocal(bestcandvar)) || SCIPisFeasGT(scip, bestcandsol, SCIPvarGetUbLocal(bestcandvar)) )
         {
            SCIPdebugMessage("selected variable's <%s> solution value is outside the domain [%g,%g] (solval: %.9f), diving aborted\n",
               SCIPvarGetName(bestcandvar), SCIPvarGetLbLocal(bestcandvar), SCIPvarGetUbLocal(bestcandvar), bestcandsol);
            assert(backtracked);
            break;
         }

         /* apply rounding of best candidate */
         if( bestcandroundup == !backtracked )
         {
	    SCIP_Real value = SCIPfeasCeil(scip, bestcandsol);
            
	    if ( SCIPisFeasIntegral(scip, bestcandsol) )
	    {
               /* only indicator variables can have integral solution value */
	       assert( SCIPvarGetType(bestcandvar) == SCIP_VARTYPE_BINARY );
	       value = 1.0;
	    }

            /* round variable up */
            SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
               SCIPvarGetName(bestcandvar), bestcandmayrounddown, bestcandmayroundup,
               bestcandsol, SCIPvarGetLbLocal(bestcandvar), SCIPvarGetUbLocal(bestcandvar),
               value, SCIPvarGetUbLocal(bestcandvar));

            SCIP_CALL( SCIPchgVarLbProbing(scip, bestcandvar, value) );
         }
         else
         {
	    SCIP_Real value = SCIPfeasFloor(scip, bestcandsol);

	    if ( SCIPisFeasIntegral(scip, bestcandsol) )
	    {
               /* only indicator variables can have integral solution value */
	       assert( SCIPvarGetType(bestcandvar) == SCIP_VARTYPE_BINARY );
	       value = 0.0;
	    }

            /* round variable down */
            SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
               SCIPvarGetName(bestcandvar), bestcandmayrounddown, bestcandmayroundup,
               bestcandsol, SCIPvarGetLbLocal(bestcandvar), SCIPvarGetUbLocal(bestcandvar),
               SCIPvarGetLbLocal(bestcandvar), value);

            SCIP_CALL( SCIPchgVarUbProbing(scip, bestcandvar, value) );
         }

         /* apply domain propagation */
         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
         if( !cutoff )
         {
            /* resolve the diving LP */
            /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
             * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
             */
#ifdef NDEBUG
            SCIP_RETCODE retstat;
            nlpiterations = SCIPgetNLPIterations(scip);
            retstat = SCIPsolveProbingLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror);
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving LP in Coefdiving heuristic; LP solve terminated with code <%d>\n",retstat);
            }
#else
            nlpiterations = SCIPgetNLPIterations(scip);
            SCIP_CALL( SCIPsolveProbingLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror) );
#endif

            if( lperror )
               break;

            /* update iteration count */
            heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;

            /* get LP solution status, objective value, and fractional variables, that should be integral */
            lpsolstat = SCIPgetLPSolstat(scip);
            cutoff = (lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
         }

         /* perform backtracking if a cutoff was detected */
         if( cutoff && !backtracked && heurdata->backtrack )
         {
            SCIPdebugMessage("  *** cutoff detected at level %d - backtracking\n", SCIPgetProbingDepth(scip));
            SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
            SCIP_CALL( SCIPnewProbingNode(scip) );
            backtracked = TRUE;
         }
         else
            backtracked = FALSE;
      }
      while( backtracked );

      if( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new objective value */
         oldobjval = objval;
         objval = SCIPgetLPObjval(scip);

         /* update pseudo cost values */
         if( SCIPisGT(scip, objval, oldobjval) )
         {
            if( bestcandroundup )
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, bestcandvar, 1.0 - bestcandfrac, objval - oldobjval, 1.0) );
            }
            else
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, bestcandvar, 0.0 - bestcandfrac, objval - oldobjval, 1.0) );
            }
         }

         /* get new fractional variables */
         SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );

         /* get indicator canditates */
         if ( nindconss > 0 )
         {
            SCIP_CALL( getIndCandVars(scip, indconss, nindconss, indcands, indcandssol, indcandfrac, &nindcands) );
         }
      }
      SCIPdebugMessage("   -> lpsolstat=%d, objval=%g/%g, nfrac=%d\n", lpsolstat, objval, searchbound, nlpcands);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Bool success;

      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
      SCIPdebugMessage("coefdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* free storage */
   if ( nindconss > 0 )
   {
      assert( indconss != NULL );
      assert( indcands != NULL );
      assert( indcandssol != NULL );
      assert( indcandfrac != NULL );

      SCIPfreeBufferArray(scip, &indcandfrac);
      SCIPfreeBufferArray(scip, &indcandssol);
      SCIPfreeBufferArray(scip, &indcands);
   }

   /* end diving */
   SCIP_CALL( SCIPendProbing(scip) );

   if( *result == SCIP_FOUNDSOL )
      heurdata->nsuccess++;

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") finished coefdiving heuristic: %d fractionals, dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", objval=%g/%g, lpsolstat=%d, cutoff=%u\n",
      SCIPgetNNodes(scip), nlpcands, divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
      SCIPretransformObj(scip, objval), SCIPretransformObj(scip, searchbound), lpsolstat, cutoff);

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the coefdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurCoefdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create coefdiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecCoefdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyCoefdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeCoefdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitCoefdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitCoefdiving) );

   /* coefdiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/coefdiving/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/coefdiving/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/coefdiving/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/coefdiving/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/coefdiving/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveubquot, TRUE, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/coefdiving/maxdiveavgquot",
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveavgquot, TRUE, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/coefdiving/maxdiveubquotnosol",
         "maximal UBQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveubquotnosol, TRUE, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/coefdiving/maxdiveavgquotnosol",
         "maximal AVGQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveavgquotnosol, TRUE, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/coefdiving/backtrack",
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->backtrack, FALSE, DEFAULT_BACKTRACK, NULL, NULL) );

   return SCIP_OKAY;
}

