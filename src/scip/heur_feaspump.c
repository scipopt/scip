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
#pragma ident "@(#) $Id: heur_feaspump.c,v 1.1 2004/07/07 09:52:41 bzfwolte Exp $"

/**@file   heur_feaspump.c
 * @brief  feasibility pump heuristic from Fischetti, Glover and Lodi 
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_feaspump.h"


#define HEUR_NAME         "feaspump"
#define HEUR_DESC         "feasibility pump heuristic from Fischetti, Glover and Lodi"
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

#define DEFAULT_DIVESTARTDEPTH      0.5 /**< minimal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT       0.1 /**< maximal fraction of diving LP iterations compared to total iteration number */
#define DEFAULT_DEPTHFAC            0.5 /**< maximal diving depth: number of binary/integer variables times depthfac */
#define DEFAULT_DEPTHFACNOSOL       2.0 /**< maximal diving depth factor if no feasible solution was found yet */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(upperbound - lowerbound)
                                         *   where diving is performed */
#define DEFAULT_MAXDIVEAVGQUOT      4.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 8.0 /**< maximal AVGQUOT when no solution was found yet */


/* locally defined heuristic data */
struct HeurData
{
   SOL*             sol;                /**< working solution */
   Real             divestartdepth;     /**< minimal relative depth to start diving */
   Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to total iteration number */
   Real             depthfac;           /**< maximal diving depth: number of binary/integer variables times depthfac */
   Real             depthfacnosol;      /**< maximal diving depth factor if no feasible solution was found yet */
   Longint          nlpiterations;      /**< LP iterations used in this heuristic */
   Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(upperbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet */
   Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet */

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
   
   Real alpha;
   LPSOLSTAT lpsolstat; 
   Bool lperror;          
   int nlpcands;
   
   Real orgobj;
   Real newobj;

   Real objval;
   Longint nsolsfound;
  
   Longint ncalls;
   Longint nlpiterations; 
   Longint maxnlpiterations;
   int depth;       
   int maxdepth;    
   int maxdivedepth;
   int divedepth;   
   
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasActNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* don't try to dive, if we are in the higher fraction of the tree, given by divestartdepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   if( depth < heurdata->divestartdepth*maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = SCIPheurGetNSolsFound(heur);
   maxnlpiterations = (1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * MAX(nlpiterations, 1000);

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + 10000);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
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
   objval = SCIPgetLPObjval(scip);
   alpha = 1;
   while( !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations )
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

         /* binary variable: 
          *   abs(x~_j - x_j) = x_j - 0   if x*_j was rounded down to 0 
          *   abs(x~_j - x_j) = 1 - x_j   if x*_j was rounded up to 1  
          */
         if( i < nbinvars )
         {
            Real solval;
            solval =  SCIPvarGetLPSol(var);
            
            /* x~_j = 0 */
            if( solval <= 0.5 )
               newobj = (1 - alpha) + alpha * orgobj;
            /* x~_j = 1 */
            else
               newobj = - (1 - alpha) + alpha * orgobj;

            if( !SCIPisFeasZero(scip, solval) )
            {
               debugMessage("dive %d/%d, LP iter %lld/%lld: i=%d  var <%s>, sol=%g, orgobj=%g, newobj=%g\n",
                  divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, i,
                  SCIPvarGetName(var), solval, orgobj, newobj);
            }
            
         }
         /* integer variable (introduction of x_j+ and x_j- to get abs(x~_j - x_j) not done here):
          *   abs(x~_j - x_j) = x_j - l     if x*_j was rounded down to l (doesn't have to be lower bound)
          *   abs(x~_j - x_j) = u - x_j     if x*_j was rounded up to u (doesn't have to be upper bound)
          */
         else if (i < nintvars )
         {
            Real solval;
            Real frac;
            solval =  SCIPvarGetLPSol(var);
            frac = SCIPfrac(scip, solval);
          
            if ( SCIPisFeasZero(scip, frac) ) 
                newobj = alpha * orgobj;
            else
            {
               /* x~_j = down(x*_j) */
               if( frac <= 0.5 )
                  newobj = (1 - alpha) + alpha * orgobj;
               /* x~_j = up(x*_j) */
               else
                  newobj = - (1 - alpha) + alpha * orgobj;
            }
         }
         /* continuous variable */
         else
            newobj = alpha * orgobj;
         
         CHECK_OKAY( SCIPchgVarObjDive(scip, var, newobj) ); 
      }

      /* resolve the diving LP */
      CHECK_OKAY( SCIPsolveDiveLP(scip, maxnlpiterations, &lperror) );
      if( lperror )
         break;

      /* update iteration count */
      heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;
      nlpiterations = SCIPgetNLPIterations(scip);

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
                  "heuristics/feaspump/divestartdepth", 
                  "minimal relative depth to start diving",
                  &heurdata->divestartdepth, DEFAULT_DIVESTARTDEPTH, 0.0, 1.0, NULL, NULL) );
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
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/feaspump/maxdiveubquot",
                  "maximal quotient (curlowerbound - lowerbound)/(upperbound - lowerbound) where diving is performed",
                  &heurdata->maxdiveubquot, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/feaspump/maxdiveavgquot", 
                  "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed",
                  &heurdata->maxdiveavgquot, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_INVALID, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/feaspump/maxdiveubquotnosol", 
                  "maximal UBQUOT when no solution was found yet",
                  &heurdata->maxdiveubquotnosol, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/feaspump/maxdiveavgquotnosol", 
                  "maximal AVGQUOT when no solution was found yet",
                  &heurdata->maxdiveavgquotnosol, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_INVALID, NULL, NULL) );
   
   return SCIP_OKAY;
}

