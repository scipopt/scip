/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_linesearchdiving.c,v 1.4 2004/11/29 12:17:15 bzfpfend Exp $"

/**@file   heur_linesearchdiving.c
 * @brief  linesearchdiving primal heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "heur_linesearchdiving.h"


#define HEUR_NAME             "linesearchdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings following the line from root solution to current solution"
#define HEUR_DISPCHAR         'l'
#define HEUR_PRIORITY         -1006000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          6
#define HEUR_MAXDEPTH         -1
#define HEUR_PSEUDONODES      FALSE     /* call heuristic at nodes where only a pseudo solution exist? */
#define HEUR_DURINGPLUNGING   FALSE     /* call heuristic during plunging? (should be FALSE for diving heuristics!) */




/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0  /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0  /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT       0.02 /**< maximal fraction of diving LP iterations compared to total iteration number */
#define DEFAULT_MAXDIVEUBQUOT       0.8  /**< maximal quotient (curlowerbound - lowerbound)/(upperbound - lowerbound)
                                          *   where diving is performed */
#define DEFAULT_MAXDIVEAVGQUOT      4.0  /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                          *   where diving is performed */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1  /**< maximal UBQUOT when no solution was found yet */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 8.0  /**< maximal AVGQUOT when no solution was found yet */



/*
 * Data structures
 */

/** primal heuristic data */
struct HeurData
{
   SOL*             sol;                /**< working solution */
   Real             minreldepth;        /**< minimal relative depth to start diving */
   Real             maxreldepth;        /**< maximal relative depth to start diving */
   Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to total iteration number */
   Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(upperbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet */
   Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet */
   Longint          nlpiterations;      /**< LP iterations used in this heuristic */
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
DECL_HEURFREE(heurFreeLinesearchdiving)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
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
DECL_HEURINIT(heurInitLinesearchdiving)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   CHECK_OKAY( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize data */
   heurdata->nlpiterations = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
DECL_HEUREXIT(heurExitLinesearchdiving)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   CHECK_OKAY( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
DECL_HEUREXEC(heurExecLinesearchdiving)
{  /*lint --e{715}*/
   HEURDATA* heurdata;
   LPSOLSTAT lpsolstat;
   VAR* var;
   VAR** lpcands;
   Real* lpcandssol;
   Real* lpcandsfrac;
   Real searchubbound;
   Real searchavgbound;
   Real searchbound;
   Real objval;
   Real oldobjval;
   Real solval;
   Real rootsolval;
   Real distquot;
   Real bestdistquot;
   Bool roundup;
   Bool bestcandroundup;
   Bool lperror;
   Bool success;
   Longint ncalls;
   Longint nsolsfound;
   Longint nlpiterations;
   Longint maxnlpiterations;
   int nlpcands;
   int startnlpcands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int bestcand;
   int c;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) )
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
   nlpiterations = SCIPgetNLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = SCIPheurGetNSolsFound(heur);
   maxnlpiterations = (1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * (nlpiterations + 10000);

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + 10000);

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      searchubbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveubquotnosol * (SCIPgetUpperbound(scip) - SCIPgetLowerbound(scip));
      searchavgbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveavgquotnosol * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
   }
   else
   {
      searchubbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveubquot * (SCIPgetUpperbound(scip) - SCIPgetLowerbound(scip));
      searchavgbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveavgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
   }
   searchbound = MIN(searchubbound, searchavgbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;


   *result = SCIP_DIDNOTFIND;

   /* start diving */
   CHECK_OKAY( SCIPstartDive(scip) );

   /* get LP objective value, and fractional variables, that should be integral */
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetLPObjval(scip);
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );

   debugMessage("(node %lld) executing linesearchdiving heuristic: depth=%d, %d fractionals, dualbound=%g, searchbound=%g\n", 
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if the last rounding was in a direction, that never destroys feasibility, we continue in any case
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   lperror = FALSE;
   divedepth = 0;
   startnlpcands = nlpcands;
   while( !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations && objval < searchbound)) )
   {
      divedepth++;

      /* choose variable fixing:
       * - in the projected space of fractional variables, extend the line segment connecting the root solution and
       *   the current LP solution up to the point, where one of the fractional variables becomes integral
       * - round this variable to the integral value
       */
      bestcand = 0;
      bestdistquot = SCIPinfinity(scip);
      bestcandroundup = FALSE;
      for( c = 0; c < nlpcands; ++c )
      {
         var = lpcands[c];
         solval = lpcandssol[c];
         rootsolval = SCIPvarGetRootSol(var);

         /* calculate distance to integral value divided by distance to root solution */
         if( SCIPisLT(scip, solval, rootsolval) )
         {
            roundup = FALSE;
            distquot = (solval - SCIPfeasFloor(scip, solval)) / (rootsolval - solval);

            /* avoid roundable candidates */
            if( SCIPvarMayRoundDown(var) )
               distquot *= 1000.0;
         }
         else if( SCIPisGT(scip, solval, rootsolval) )
         {
            roundup = TRUE;
            distquot = (SCIPfeasCeil(scip, solval) - solval) / (solval - rootsolval);

            /* avoid roundable candidates */
            if( SCIPvarMayRoundUp(var) )
               distquot *= 1000.0;
         }
         else
         {
            roundup = FALSE;
            distquot = SCIPinfinity(scip);
         }

         /* check, if candidate is new best candidate */
         if( distquot < bestdistquot )
         {
            bestcand = c;
            bestcandroundup = roundup;
            bestdistquot = distquot;
         }
      }

      /* create solution from diving LP and try to round it */
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      CHECK_OKAY( SCIProundSol(scip, heurdata->sol, &success) );
      if( success )
      {
         debugMessage("linesearchdiving found roundable primal solution: obj=%g\n",
            SCIPgetSolOrigObj(scip, heurdata->sol));
         
         /* try to add solution to SCIP */
         CHECK_OKAY( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, &success) );
         
         /* check, if solution was feasible and good enough */
         if( success )
         {
            debugMessage(" -> solution was feasible and good enough\n");
            *result = SCIP_FOUNDSOL;
         }
      }

      /* apply rounding of best candidate */
      var = lpcands[bestcand];

      /* if the variable is already fixed, abort diving due to numerical troubles */
      if( SCIPgetVarLbDive(scip, var) >= SCIPgetVarUbDive(scip, var) - 0.5 )
         break;

      if( bestcandroundup )
      {
         /* round variable up */
         debugMessage("  dive %d/%d, LP iter %lld/%lld: var <%s>, sol=%g, root=%g, [%g,%g] -> [%g,%g]\n",
            divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
            SCIPvarGetName(var), lpcandssol[bestcand], SCIPvarGetRootSol(var),
            SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var),
            SCIPfeasCeil(scip, lpcandssol[bestcand]), SCIPgetVarUbDive(scip, var));
         CHECK_OKAY( SCIPchgVarLbDive(scip, var, SCIPfeasCeil(scip, lpcandssol[bestcand])) );
      }
      else
      {
         /* round variable down */
         debugMessage("  dive %d/%d, LP iter %lld/%lld: var <%s>, sol=%g, root=%g, [%g,%g] -> [%g,%g]\n",
            divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
            SCIPvarGetName(var), lpcandssol[bestcand], SCIPvarGetRootSol(var),
            SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var),
            SCIPgetVarLbDive(scip, var), SCIPfeasFloor(scip, lpcandssol[bestcand]));
         CHECK_OKAY( SCIPchgVarUbDive(scip, var, SCIPfeasFloor(scip, lpcandssol[bestcand])) );
      }

      /* resolve the diving LP */
      CHECK_OKAY( SCIPsolveDiveLP(scip, maxnlpiterations, &lperror) );
      if( lperror )
         break;

      /* update iteration count */
      heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;
      nlpiterations = SCIPgetNLPIterations(scip);

      /* get LP solution status, objective value, and fractional variables, that should be integral */
      lpsolstat = SCIPgetLPSolstat(scip);
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new objective value */
         oldobjval = objval;
         objval = SCIPgetLPObjval(scip);

         /* update pseudo cost values */
         if( SCIPisGT(scip, objval, oldobjval) )
         {
            if( bestcandroundup )
            {
               CHECK_OKAY( SCIPupdateVarPseudocost(scip, lpcands[bestcand], 1.0-lpcandsfrac[bestcand], 
                     objval - oldobjval, 1.0) );
            }
            else
            {
               CHECK_OKAY( SCIPupdateVarPseudocost(scip, lpcands[bestcand], 0.0-lpcandsfrac[bestcand], 
                     objval - oldobjval, 1.0) );
            }
         }

         /* get new fractional variables */
         CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );
      }
      debugMessage("   -> lpsolstat=%d, objval=%g, nfrac=%d\n", lpsolstat, objval, nlpcands);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      Bool success;

      /* create solution from diving LP */
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      debugMessage("linesearchdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      CHECK_OKAY( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         debugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   CHECK_OKAY( SCIPendDive(scip) );

   debugMessage("linesearchdiving heuristic finished\n");

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the linesearchdiving primal heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurLinesearchdiving(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create linesearchdiving primal heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING,
         heurFreeLinesearchdiving, heurInitLinesearchdiving, heurExitLinesearchdiving, heurExecLinesearchdiving,
         heurdata) );

   /* add linesearchdiving primal heuristic parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/linesearchdiving/minreldepth", 
         "minimal relative depth to start diving",
         &heurdata->minreldepth, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/linesearchdiving/maxreldepth", 
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/linesearchdiving/maxlpiterquot", 
         "maximal fraction of diving LP iterations compared to total iteration number",
         &heurdata->maxlpiterquot, DEFAULT_MAXLPITERQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/linesearchdiving/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(upperbound - lowerbound) where diving is performed",
         &heurdata->maxdiveubquot, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/linesearchdiving/maxdiveavgquot", 
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed",
         &heurdata->maxdiveavgquot, DEFAULT_MAXDIVEAVGQUOT, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/linesearchdiving/maxdiveubquotnosol", 
         "maximal UBQUOT when no solution was found yet",
         &heurdata->maxdiveubquotnosol, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/linesearchdiving/maxdiveavgquotnosol", 
         "maximal AVGQUOT when no solution was found yet",
         &heurdata->maxdiveavgquotnosol, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
