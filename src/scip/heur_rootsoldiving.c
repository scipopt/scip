/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_rootsoldiving.c,v 1.9 2005/01/21 09:16:54 bzfpfend Exp $"

/**@file   heur_rootsoldiving.c
 * @brief  LP diving heuristic that changes variable's objective values using root LP solution as guide
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_rootsoldiving.h"


#define HEUR_NAME         "rootsoldiving"
#define HEUR_DESC         "LP diving heuristic that changes variable's objective values using root LP solution as guide"
#define HEUR_DISPCHAR     's'
#define HEUR_PRIORITY     -1005000
#define HEUR_FREQ         10
#define HEUR_FREQOFS       5
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
DECL_HEURFREE(heurFreeRootsoldiving) /*lint --e{715}*/
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
DECL_HEURINIT(heurInitRootsoldiving) /*lint --e{715}*/
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
DECL_HEUREXIT(heurExitRootsoldiving) /*lint --e{715}*/
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
DECL_HEUREXEC(heurExecRootsoldiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata; 
   VAR** vars;
   Real* rootsol;
   int* softroundings;
   int nvars;           
   int nbinvars;
   int nintvars;
   int nlpcands;
   LPSOLSTAT lpsolstat; 
   Real alpha;
   Real oldobj;
   Real newobj;
   Bool lperror;          
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
   int ncycles;
   int i;

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

   /* calculate the maximal diving depth */
   nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   if( SCIPgetNSolsFound(scip) == 0 )
      maxdivedepth = heurdata->depthfacnosol * nvars;
   else
      maxdivedepth = heurdata->depthfac * nvars;
   maxdivedepth = MAX(maxdivedepth, 10);


   *result = SCIP_DIDNOTFIND;

   /* get all varibales of LP */
   CHECK_OKAY( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* get root solution value of all binary and integer variables */
   CHECK_OKAY( SCIPallocBufferArray(scip, &rootsol, nbinvars + nintvars) );
   for( i = 0; i < nbinvars + nintvars; i++ )
      rootsol[i] = SCIPvarGetRootSol(vars[i]);

   /* initialize array storing the preferred soft rounding direction for the variables */
   CHECK_OKAY( SCIPallocBufferArray(scip, &softroundings, nbinvars + nintvars) );
   clearMemoryArray(softroundings, nbinvars + nintvars);

   /* start diving */
   CHECK_OKAY( SCIPstartDive(scip) );

   /* get number of fractional variables, that should be integral */
   nlpcands = SCIPgetNLPBranchCands(scip);

   debugMessage("(node %lld) executing rootsoldiving heuristic: depth=%d, %d fractionals, dualbound=%g, maxnlpiterations=%lld, maxdivedepth=%d\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), maxnlpiterations, maxdivedepth);

   lperror = FALSE;
   divedepth = 0;
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   alpha = 0.9; /**@todo make thís constant a parameter setting */
   ncycles = 0;
   startnlpcands = nlpcands;
   while( !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0 && ncycles < 10
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
         debugMessage("rootsoldiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));
         
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

      /* round solution x* from diving LP: 
       *   - x~_j = down(x*_j)    if x*_j is integer or binary variable and x*_j <= root solution_j
       *   - x~_j = up(x*_j)      if x*_j is integer or binary variable and x*_j  > root solution_j
       *   - x~_j = x*_j          if x*_j is continuous variable
       * change objective function in diving LP:
       *   - if x*_j is integral, or j is a continuous variable, set obj'_j = alpha * obj_j
       *   - otherwise, set obj'_j = alpha * obj_j + sign(x*_j - x~_j)
       */
      for( i = 0; i < nbinvars + nintvars; i++ )
      {
         VAR* var;
         Real solval;

         var = vars[i];
         oldobj = SCIPgetVarObjDive(scip, var);
         newobj = oldobj;
         
         solval =  SCIPvarGetLPSol(var);
         if( SCIPisIntegral(scip, solval) )
         {
            /* if the variable became integral after a soft rounding, fix it to this value;
             * otherwise, fade out the objective value
             */
            if( softroundings[i] != 0 )
            {
               solval = SCIPfeasFloor(scip, solval);
               CHECK_OKAY( SCIPchgVarLbDive(scip, var, solval) );
               CHECK_OKAY( SCIPchgVarUbDive(scip, var, solval) );
            }
            else
               newobj = alpha * oldobj;
         }
         else if( solval <= rootsol[i] )
         {
            /* if the variable was soft rounded most of the time downwards, round it downwards by changing the bounds;
             * otherwise, apply soft rounding by changing the objective value
             */
            softroundings[i]--;
            if( softroundings[i] <= -10 )
            {
               CHECK_OKAY( SCIPchgVarUbDive(scip, var, SCIPfeasFloor(scip, solval)) );
            }
            else
               newobj = alpha * oldobj + 1.0;
         }
         else
         {
            /* if the variable was soft rounded most of the time upwards, round it upwards by changing the bounds;
             * otherwise, apply soft rounding by changing the objective value
             */
            softroundings[i]++;
            if( softroundings[i] >= +10 )
            {
               CHECK_OKAY( SCIPchgVarLbDive(scip, var, SCIPfeasCeil(scip, solval)) );
            }
            else
               newobj = alpha * oldobj - 1.0;
         }
         
         debugMessage("dive %d/%d, LP iter %lld/%lld: i=%d  var <%s>, solval= %g, rootsol=%g, oldobj=%g, newobj=%g\n",
            divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, i,
            SCIPvarGetName(var), solval, rootsol[i], oldobj, newobj);
         
         CHECK_OKAY( SCIPchgVarObjDive(scip, var, newobj) ); 
      }
      
      /* fade out the objective values of the continuous variables */
      for( i = nbinvars + nintvars; i < nvars; i++ )
      {
         VAR* var;

         var = vars[i];
         oldobj = SCIPvarGetObj(var);
         newobj = alpha * oldobj;
         
         debugMessage("dive %d/%d, LP iter %lld/%lld: i=%d  var <%s>, oldobj=%g, newobj=%g\n",
            divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, i,
            SCIPvarGetName(var), oldobj, newobj);
         
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

      /* if no LP iterations were performed, we stayed at the same solution -> count this cycling */
      if( iter == nlpiterations )
         ncycles++;
      else
         ncycles = 0;

      /* get LP solution status and number of fractional variables, that should be integral */
      lpsolstat = SCIPgetLPSolstat(scip);
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
         nlpcands = SCIPgetNLPBranchCands(scip);
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
      debugMessage("rootsoldiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

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

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &softroundings);
   SCIPfreeBufferArray(scip, &rootsol);

   debugMessage("rootsoldiving heuristic finished\n");

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the rootsoldiving heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurRootsoldiving(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
                  HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING,
                  heurFreeRootsoldiving, heurInitRootsoldiving, heurExitRootsoldiving, heurExecRootsoldiving,
                  heurdata) );

   /* rootsoldiving heuristic parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/rootsoldiving/minreldepth", 
                  "minimal relative depth to start diving",
                  &heurdata->minreldepth, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/rootsoldiving/maxreldepth", 
                  "maximal relative depth to start diving",
                  &heurdata->maxreldepth, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/rootsoldiving/maxlpiterquot", 
                  "maximal fraction of diving LP iterations compared to total iteration number",
                  &heurdata->maxlpiterquot, DEFAULT_MAXLPITERQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/rootsoldiving/depthfac",
                  "maximal diving depth: number of binary/integer variables times depthfac",
                  &heurdata->depthfac, DEFAULT_DEPTHFAC, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/rootsoldiving/depthfacnosol",
                  "maximal diving depth factor if no feasible solution was found yet",
                  &heurdata->depthfacnosol, DEFAULT_DEPTHFACNOSOL, 0.0, REAL_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}

