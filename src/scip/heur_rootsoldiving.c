/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_rootsoldiving.c,v 1.31 2006/03/23 14:04:48 bzfpfend Exp $"

/**@file   heur_rootsoldiving.c
 * @brief  LP diving heuristic that changes variable's objective values using root LP solution as guide
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_rootsoldiving.h"


#define HEUR_NAME         "rootsoldiving"
#define HEUR_DESC         "LP diving heuristic that changes variable's objective values using root LP solution as guide"
#define HEUR_DISPCHAR     'S'
#define HEUR_PRIORITY     -1005000
#define HEUR_FREQ         20
#define HEUR_FREQOFS       5
#define HEUR_MAXDEPTH     -1
#define HEUR_PSEUDONODES      FALSE     /* call heuristic at nodes where only a pseudo solution exist? */
#define HEUR_DURINGPLUNGING   FALSE     /* call heuristic during plunging? (should be FALSE for diving heuristics!) */
#define HEUR_DURINGLPLOOP     FALSE     /* call heuristic during the LP price-and-cut loop? */
#define HEUR_AFTERNODE        TRUE      /* call heuristic after or before the current node was solved? */



/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.01 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXSOLS              -1 /**< total number of feasible solutions found up to which heuristic is called
                                              *   (-1: no limit) */
#define DEFAULT_DEPTHFAC            0.5 /**< maximal diving depth: number of binary/integer variables times depthfac */
#define DEFAULT_DEPTHFACNOSOL       2.0 /**< maximal diving depth factor if no feasible solution was found yet */

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */



/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int                   maxlpiterofs;       /**< additional number of allowed LP iterations */
   int                   maxsols;            /**< total number of feasible solutions found up to which heuristic is called
                                              *   (-1: no limit) */
   SCIP_Real             depthfac;           /**< maximal diving depth: number of binary/integer variables times depthfac */
   SCIP_Real             depthfacnosol;      /**< maximal diving depth factor if no feasible solution was found yet */
   SCIP_Longint          nlpiterations;      /**< LP iterations used in this heuristic */
};

/*
 * Callback methods
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRootsoldiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitRootsoldiving) /*lint --e{715}*/
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

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRootsoldiving) /*lint --e{715}*/
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


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolRootsoldiving NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolRootsoldiving NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRootsoldiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata; 
   SCIP_VAR** vars;
   SCIP_Real* rootsol;
   int* softroundings;
   int nvars;           
   int nbinvars;
   int nintvars;
   int nlpcands;
   SCIP_LPSOLSTAT lpsolstat; 
   SCIP_Real absstartobjval;
   SCIP_Real objstep;
   SCIP_Real alpha;
   SCIP_Real oldobj;
   SCIP_Real newobj;
   SCIP_Bool lperror;          
   SCIP_Longint nsolsfound;
   SCIP_Longint ncalls;
   SCIP_Longint nlpiterations; 
   SCIP_Longint maxnlpiterations;
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

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only apply heuristic, if only a few solutions have been found */
   if( heurdata->maxsols >= 0 && SCIPgetNSolsFound(scip) >= heurdata->maxsols )
      return SCIP_OKAY;

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < heurdata->minreldepth*maxdepth || depth > heurdata->maxreldepth*maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNNodeLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = SCIPheurGetNSolsFound(heur);
   maxnlpiterations = (SCIP_Longint)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * nlpiterations);
   maxnlpiterations += heurdata->maxlpiterofs;

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + MINLPITER);

   /* calculate the maximal diving depth */
   nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   if( SCIPgetNSolsFound(scip) == 0 )
      maxdivedepth = (int)(heurdata->depthfacnosol * nvars);
   else
      maxdivedepth = (int)(heurdata->depthfac * nvars);
   maxdivedepth = MAX(maxdivedepth, 10);


   *result = SCIP_DIDNOTFIND;

   /* get all varibales of LP */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* get root solution value of all binary and integer variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &rootsol, nbinvars + nintvars) );
   for( i = 0; i < nbinvars + nintvars; i++ )
      rootsol[i] = SCIPvarGetRootSol(vars[i]);

   /* get current LP objective value, and calculate length of a single step in an objective coefficient */
   absstartobjval = SCIPgetLPObjval(scip);
   absstartobjval = ABS(absstartobjval);
   absstartobjval = MAX(absstartobjval, 1.0);
   objstep = absstartobjval / 10.0;

   /* initialize array storing the preferred soft rounding direction for the variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &softroundings, nbinvars + nintvars) );
   BMSclearMemoryArray(softroundings, nbinvars + nintvars);

   /* start diving */
   SCIP_CALL( SCIPstartDive(scip) );

   /* get number of fractional variables, that should be integral */
   nlpcands = SCIPgetNLPBranchCands(scip);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing rootsoldiving heuristic: depth=%d, %d fractionals, dualbound=%g, maxnlpiterations=%"SCIP_LONGINT_FORMAT", maxdivedepth=%d, LPobj=%g, objstep=%g\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), maxnlpiterations, maxdivedepth,
      SCIPgetLPObjval(scip), objstep);

   lperror = FALSE;
   divedepth = 0;
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   alpha = 0.9; /**@todo make thís constant a parameter setting */
   ncycles = 0;
   startnlpcands = nlpcands;
   while( !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0 && ncycles < 10
      && (FALSE /*??????????????divedepth < 10*/
#if 1 /*?????????????????????*/
         || divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations)
#else
         || divedepth < maxdepth
         || nlpcands <= startnlpcands - divedepth/2
         || (nlpcands <= startnlpcands - divedepth/5 && divedepth < maxdivedepth
            && heurdata->nlpiterations < maxnlpiterations)
#endif
          ) )
   {
      SCIP_Bool success;
      int hardroundingidx;
      int hardroundingdir;
      SCIP_Real hardroundingoldbd;
      SCIP_Real hardroundingnewbd;

      /* create solution from diving LP and try to round it */
      SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
      SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

      if( success )
      {
         SCIPdebugMessage("rootsoldiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));
         
         /* try to add solution to SCIP */
         SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, &success) );
            
         /* check, if solution was feasible and good enough */
         if( success )
         {
            SCIPdebugMessage(" -> solution was feasible and good enough\n");
            *result = SCIP_FOUNDSOL;
         }
      }

      divedepth++;
      hardroundingidx = -1;
      hardroundingdir = 0;
      hardroundingoldbd = 0.0;
      hardroundingnewbd = 0.0;

      SCIPdebugMessage("dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT":\n", divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations);
         
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
         SCIP_VAR* var;
         SCIP_Real solval;

         var = vars[i];
         oldobj = SCIPgetVarObjDive(scip, var);
         newobj = oldobj;
         
         solval =  SCIPvarGetLPSol(var);
         if( SCIPisFeasIntegral(scip, solval) )
         {
            /* if the variable became integral after a soft rounding, fix it to this value;
             * otherwise, fade out the objective value
             */
            if( softroundings[i] != 0 )
            {
               solval = SCIPfeasFloor(scip, solval);
               SCIP_CALL( SCIPchgVarLbDive(scip, var, solval) );
               SCIP_CALL( SCIPchgVarUbDive(scip, var, solval) );
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
            if( softroundings[i] <= -10 && hardroundingidx == -1 )
            {
               SCIPdebugMessage(" -> hard rounding <%s>[%g] <= %g\n", 
                  SCIPvarGetName(var), solval, SCIPfeasFloor(scip, solval));
               hardroundingidx = i;
               hardroundingdir = -1;
               hardroundingoldbd = SCIPgetVarUbDive(scip, var);
               hardroundingnewbd = SCIPfeasFloor(scip, solval);
               SCIP_CALL( SCIPchgVarUbDive(scip, var, hardroundingnewbd) );
            }
            else
               newobj = alpha * oldobj + objstep;
         }
         else
         {
            /* if the variable was soft rounded most of the time upwards, round it upwards by changing the bounds;
             * otherwise, apply soft rounding by changing the objective value
             */
            softroundings[i]++;
            if( softroundings[i] >= +10 && hardroundingidx == -1 )
            {
               SCIPdebugMessage(" -> hard rounding <%s>[%g] >= %g\n", 
                  SCIPvarGetName(var), solval, SCIPfeasCeil(scip, solval));
               hardroundingidx = i;
               hardroundingdir = +1;
               hardroundingoldbd = SCIPgetVarLbDive(scip, var);
               hardroundingnewbd = SCIPfeasCeil(scip, solval);
               SCIP_CALL( SCIPchgVarLbDive(scip, var, hardroundingnewbd) );
            }
            else
               newobj = alpha * oldobj - objstep;
         }
         
         SCIPdebugMessage(" -> i=%d  var <%s>, solval=%g, rootsol=%g, oldobj=%g, newobj=%g\n",
            i, SCIPvarGetName(var), solval, rootsol[i], oldobj, newobj);
         
         SCIP_CALL( SCIPchgVarObjDive(scip, var, newobj) ); 
      }
      
      /* fade out the objective values of the continuous variables */
      for( i = nbinvars + nintvars; i < nvars; i++ )
      {
         SCIP_VAR* var;

         var = vars[i];
         oldobj = SCIPgetVarObjDive(scip, var);
         newobj = alpha * oldobj;
         
         SCIPdebugMessage(" -> i=%d  var <%s>, solval=%g, oldobj=%g, newobj=%g\n", 
            i, SCIPvarGetName(var), SCIPvarGetLPSol(var), oldobj, newobj);
         
         SCIP_CALL( SCIPchgVarObjDive(scip, var, newobj) ); 
      }

   SOLVEAGAIN:
      /* resolve the diving LP */
      nlpiterations = SCIPgetNLPIterations(scip);
      SCIP_CALL( SCIPsolveDiveLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror) );
      if( lperror )
         break;

      /* update iteration count */
      heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;

      /* if no LP iterations were performed, we stayed at the same solution -> count this cycling */
      if( SCIPgetNLPIterations(scip) == nlpiterations )
         ncycles++;
      else
         ncycles = 0;

      /* get LP solution status and number of fractional variables, that should be integral */
      lpsolstat = SCIPgetLPSolstat(scip);
      if( lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE && hardroundingidx != -1 )
      {
         SCIP_VAR* var;

         var = vars[hardroundingidx];

         /* round the hard rounded variable to the opposite direction and resolve the LP */
         if( hardroundingdir == -1 )
         {
            SCIPdebugMessage(" -> opposite hard rounding <%s> >= %g\n", SCIPvarGetName(var), hardroundingnewbd + 1.0);
            SCIP_CALL( SCIPchgVarUbDive(scip, var, hardroundingoldbd) );
            SCIP_CALL( SCIPchgVarLbDive(scip, var, hardroundingnewbd + 1.0) );
         }
         else
         {
            SCIPdebugMessage(" -> opposite hard rounding <%s> <= %g\n", SCIPvarGetName(var), hardroundingnewbd - 1.0);
            SCIP_CALL( SCIPchgVarLbDive(scip, var, hardroundingoldbd) );
            SCIP_CALL( SCIPchgVarUbDive(scip, var, hardroundingnewbd - 1.0) );
         }
         hardroundingidx = -1;
         goto SOLVEAGAIN;
      }
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
         nlpcands = SCIPgetNLPBranchCands(scip);
      SCIPdebugMessage("   -> lpsolstat=%d, nfrac=%d\n", lpsolstat, nlpcands);
   }

   SCIPdebugMessage("---> diving finished: lpsolstat = %d, depth %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT"\n", 
      lpsolstat, divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations);

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Bool success;

      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
      SCIPdebugMessage("rootsoldiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   SCIP_CALL( SCIPendDive(scip) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &softroundings);
   SCIPfreeBufferArray(scip, &rootsol);

   SCIPdebugMessage("rootsoldiving heuristic finished\n");

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the rootsoldiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRootsoldiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING, HEUR_DURINGLPLOOP, HEUR_AFTERNODE,
         heurFreeRootsoldiving, heurInitRootsoldiving, heurExitRootsoldiving, 
         heurInitsolRootsoldiving, heurExitsolRootsoldiving, heurExecRootsoldiving,
         heurdata) );

   /* rootsoldiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/minreldepth", 
         "minimal relative depth to start diving",
         &heurdata->minreldepth, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/maxreldepth", 
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/maxlpiterquot", 
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/rootsoldiving/maxlpiterofs", 
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/rootsoldiving/maxsols", 
         "total number of feasible solutions found up to which heuristic is called (-1: no limit)",
         &heurdata->maxsols, DEFAULT_MAXSOLS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/depthfac",
         "maximal diving depth: number of binary/integer variables times depthfac",
         &heurdata->depthfac, DEFAULT_DEPTHFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/depthfacnosol",
         "maximal diving depth factor if no feasible solution was found yet",
         &heurdata->depthfacnosol, DEFAULT_DEPTHFACNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}

