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
#pragma ident "@(#) $Id: heur_feaspump.c,v 1.6 2004/11/29 12:17:15 bzfpfend Exp $"

/**@file   heur_feaspump.c
 * @brief  feasibility pump heuristic by Fischetti, Glover and Lodi 
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_feaspump.h"


#define HEUR_NAME         "feaspump"
#define HEUR_DESC         "feasibility pump heuristic by Fischetti, Glover and Lodi"
#define HEUR_DISPCHAR     'F'
#define HEUR_PRIORITY     -1003000
#define HEUR_FREQ         10
#define HEUR_FREQOFS       3
#define HEUR_MAXDEPTH     -1
#define HEUR_PSEUDONODES  FALSE         /** call heuristic at nodes where only a pseudo solution exist? */
#define HEUR_DURINGPLUNGING   FALSE     /** call heuristic during plunging? (should be FALSE for diving heuristics!) */



/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH        0.0  /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH        1.0  /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.01 /**< maximal fraction of diving LP iterations compared to total iteration number */
#define DEFAULT_DEPTHFAC           0.5  /**< maximal diving depth: number of binary/integer variables times depthfac */
#define DEFAULT_DEPTHFACNOSOL      2.0  /**< maximal diving depth factor if no feasible solution was found yet */


/* locally defined heuristic data */
struct HeurData
{
   SOL*             sol;                /**< working solution */
   Real             minreldepth;        /**< minimal relative depth to start diving */
   Real             maxreldepth;        /**< maximal relative depth to start diving */
   Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to total iteration number */
   Real             depthfac;           /**< maximal diving depth: number of binary/integer variables times depthfac */
   Real             depthfacnosol;      /**< maximal diving depth factor if no feasible solution was found yet */
   Longint          nlpiterations;      /**< LP iterations used in this heuristic */
};




/*
 * Callback methods
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
DECL_HEURFREE(heurFreeFeaspump) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

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
DECL_HEURINIT(heurInitFeaspump) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

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
DECL_HEUREXIT(heurExitFeaspump) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   CHECK_OKAY( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
DECL_HEUREXEC(heurExecFeaspump) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata; 

   VAR** vars;
   int nvars;           
   int nbinvars;
   int nintvars;
   int i;
   
   LPSOLSTAT lpsolstat; 
   Real alpha;
   Bool lperror;          
   int nlpcands;
   
   Real orgobj;
   Real newobj;

   Longint nsolsfound;
   Longint ncalls;
   Longint iter;
   Longint nlpiterations; 
   Longint maxnlpiterations;
   int depth;       
   int maxdepth;    
   int maxdivedepth;
   int divedepth;   
   int startnlpcands;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
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


   *result = SCIP_DIDNOTFIND;

   /* calculate the maximal diving depth */
   nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   if( SCIPgetNSolsFound(scip) == 0 )
      maxdivedepth = heurdata->depthfacnosol * nvars;
   else
      maxdivedepth = heurdata->depthfac * nvars;

   /* start diving */
   CHECK_OKAY( SCIPstartDive(scip) );

   /* get all variables of LP and number of fractional variables in LP solution that should be integral */
   CHECK_OKAY( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nlpcands = SCIPgetNLPBranchCands(scip);
   assert(0 <= nlpcands && nlpcands <= nbinvars + nintvars);

   debugMessage("(node %lld) executing feaspump heuristic: depth=%d, %d fractionals, dualbound=%g, maxnlpiterations=%lld, maxdivedepth=%d\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), maxnlpiterations, maxdivedepth);

   lperror = FALSE;
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   divedepth = 0;
   alpha = 1.0;
   startnlpcands = nlpcands;
   while( !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations)) )
   {
      Bool success;
         
      /* create solution from diving LP and try to round it */
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      CHECK_OKAY( SCIProundSol(scip, heurdata->sol, &success) );

      if( success )
      {
         debugMessage("feaspump found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));
         
         /* try to add solution to SCIP */
         CHECK_OKAY( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, &success) );
            
         /* check, if solution was feasible and good enough */
         if( success )
         {
            debugMessage(" -> solution was feasible and good enough\n");
            *result = SCIP_FOUNDSOL;
         }
      }

      divedepth++;
      alpha *= 0.9;
      if( alpha < 0.05 )
         alpha = 0.0;

      /* round solution x* from diving LP: 
       *     x~_j = down(x*_j)    if x*_j is integer or binary variable and frac(x*_j) <= 0.5
       *     x~_j = up(x*_j)      if x*_j is integer or binary variable and frac(x*_j)  > 0.5
       *     x~_j = x*_j          if x*_j is continuous variable
       * change objective function in diving LP:
       *     (1 - alpha) * SUM abs(x~_j - x_j) + alpha * objective function of original LP
       */
      for( i = 0; i < nvars; i++ )
      {
         VAR* var;

         var = vars[i];
         orgobj = SCIPvarGetObj(var);

         if( i < nbinvars )
         {
            Real solval;

            /* binary variable: 
             *   abs(x~_j - x_j) = x_j - 0   if x*_j was rounded down to 0 
             *   abs(x~_j - x_j) = 1 - x_j   if x*_j was rounded up to 1  
             */
            solval =  SCIPvarGetLPSol(var);
            if( solval <= 0.5 )
               newobj = (1.0 - alpha) + alpha * orgobj;
            else
               newobj = - (1.0 - alpha) + alpha * orgobj;

            debugMessage("dive %d/%d, LP iter %lld/%lld: i=%d  var <%s>, sol=%g, orgobj=%g, newobj=%g\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, i,
               SCIPvarGetName(var), solval, orgobj, newobj);
         }
         else if( i < nintvars )
         {
            Real solval;
            Real frac;

            /* integer variable (introduction of x_j+ and x_j- to get abs(x~_j - x_j) not done here;
             * instead, we try to push the variable to the nearest integer, if it is fractional):
             *   abs(x~_j - x_j) = x_j - l     if x*_j was rounded down to l (doesn't have to be lower bound)
             *   abs(x~_j - x_j) = u - x_j     if x*_j was rounded up to u (doesn't have to be upper bound)
             */
            solval =  SCIPvarGetLPSol(var);
            frac = SCIPfeasFrac(scip, solval);
          
            if ( SCIPisFeasZero(scip, frac) ) 
                newobj = alpha * orgobj;
            else if( frac <= 0.5 )
               newobj = (1.0 - alpha) + alpha * orgobj;
            else
               newobj = - (1.0 - alpha) + alpha * orgobj;
         }
         else
         {
            /* continuous variable */
            newobj = alpha * orgobj;
         }
         
         CHECK_OKAY( SCIPchgVarObjDive(scip, var, newobj) ); 
      }

      /* resolve the diving LP */
      CHECK_OKAY( SCIPsolveDiveLP(scip, maxnlpiterations, &lperror) );
      if( lperror )
         break;

      /* update iteration count */
      iter = SCIPgetNLPIterations(scip);
      heurdata->nlpiterations += iter - nlpiterations;
      nlpiterations = SCIPgetNLPIterations(scip);

      /* if we stayed at the same LP solution, rapidly decrease alpha */
      if( iter == nlpiterations )
      {
         /* if alpha == 0.0 was tried, abort */
         if( alpha == 0.0 )
            break;
         alpha *= 0.5;
      }

      /* get LP solution status and number of fractional variables, that should be integral */
      lpsolstat = SCIPgetLPSolstat(scip);
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         nlpcands = SCIPgetNLPBranchCands(scip);
      }
      debugMessage("   -> lpsolstat=%d, nfrac=%d\n", lpsolstat, nlpcands);
   }

   debugMessage("---> diving finished: lpsolstat = %d, depth %d/%d, LP iter %lld/%lld\n", 
      lpsolstat, divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations);

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      Bool success;

      /* create solution from diving LP */
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      debugMessage("feaspump found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

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

   debugMessage("feaspump heuristic finished\n");

   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the feaspump heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurFeaspump(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
                  HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING,
                  heurFreeFeaspump, heurInitFeaspump, heurExitFeaspump, heurExecFeaspump,
                  heurdata) );

   /* feaspump heuristic parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/feaspump/minreldepth", 
                  "minimal relative depth to start diving",
                  &heurdata->minreldepth, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/feaspump/maxreldepth", 
                  "maximal relative depth to start diving",
                  &heurdata->maxreldepth, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/feaspump/maxlpiterquot", 
                  "maximal fraction of diving LP iterations compared to total iteration number",
                  &heurdata->maxlpiterquot, DEFAULT_MAXLPITERQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/feaspump/depthfac",
                  "maximal diving depth: number of binary/integer variables times depthfac",
                  &heurdata->depthfac, DEFAULT_DEPTHFAC, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/feaspump/depthfacnosol",
                  "maximal diving depth factor if no feasible solution was found yet",
                  &heurdata->depthfacnosol, DEFAULT_DEPTHFACNOSOL, 0.0, REAL_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}

