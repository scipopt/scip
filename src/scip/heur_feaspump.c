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
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_feaspump.c,v 1.25 2005/03/21 11:37:31 bzfpfend Exp $"

/**@file   heur_feaspump.c
 * @brief  feasibility pump primal heuristic
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "scip/heur_feaspump.h"


#define HEUR_NAME             "feaspump"
#define HEUR_DESC             "feasibility pump heuristic by Fischetti, Glover and Lodi"
#define HEUR_DISPCHAR         'F'
#define HEUR_PRIORITY         -1000000
#define HEUR_FREQ             20
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_PSEUDONODES      FALSE      /* call heuristic at nodes where only a pseudo solution exist? */
#define HEUR_DURINGPLUNGING   FALSE      /* call heuristic during plunging? (should be FALSE for diving heuristics!) */
#define HEUR_AFTERNODE        TRUE      /* call heuristic after or before the current node was solved? */

#define DEFAULT_MAXLPITERQUOT    0.01   /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS    10000   /**< additional number of allowed LP iterations */
#define DEFAULT_MAXLOOPS        10000   /**< maximal number of pumping rounds (-1: no limit) */
#define DEFAULT_MINFLIPS           10   /**< minimum number of random variables to flip, if a 1-cycle is encountered */
#define DEFAULT_CYCLELENGTH         3   /**< maximum length of cycles to be checked explicitly in each round */
#define DEFAULT_PERTURBFREQ       100   /**< number of iterations until a random perturbation is forced */
#define DEFAULT_OBJFACTOR         1.0   /**< factor by which the regard of the objective is decreased in each round, 
                                         * 1.0 for dynamic, depending on solutions already found */


/** primal heuristic data */
struct HeurData
{
   SOL*             sol;                /**< working solution */
   SOL*             roundedsol;         /**< rounded solution */ 
   Longint          nlpiterations;      /**< number of LP iterations used in this heuristic */
   Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int              maxlpiterofs;       /**< additional number of allowed LP iterations */
   Real             objfactor;          /**< factor by which the regard of the objective is decreased in each round, 
                                         *   1.0 for dynamic, depending on solutions already found */
   int              maxloops;           /**< maximum number of loops (-1: no limit) */ 
   int              minflips;           /**< minimum number of random variables to flip, if a 1-cycle is encountered */
   int              cyclelength;        /**< maximum length of cycles to be checked explicitly in each round */
   int              perturbfreq;        /**< number of iterations until a random perturbation is forced */
   unsigned int     randseed;           /**< seed value for random number generator */
};

/** checks whether a variable is one of the currently most fractional ones */
static
void insertFlipCand( 
   VAR**            mostfracvars,       /**< sorted array of the currently most fractional variables */
   Real*            mostfracvals,       /**< array of their fractionality, decreasingly sorted */
   int*             nflipcands,         /**< number of fractional variables already labeled to be flipped*/
   int              maxnflipcands,      /**< typically randomized number of maximum amount of variables to flip */
   VAR*             var,                /**< variable to be checked */
   Real             frac                /**< fractional value of the variable */
   )
{
   int i;

   assert(mostfracvars != NULL);
   assert(mostfracvals != NULL);
   assert(nflipcands != NULL);

   /* instead of the fractional value use the fractionality */
   if( frac > 0.5 )
      frac = 1 - frac;

   /* if there are already enough candidates and the variable is less fractional, return, else reserve the last entry */
   if( *nflipcands >= maxnflipcands )
   {
      if( frac <= mostfracvals[*nflipcands-1] )
         return;
      else
         (*nflipcands)--;
   }

   /* shifting var and frac through the (sorted) arrays */
   for( i = *nflipcands; i > 0 && mostfracvals[i-1] < frac; i-- )
   {
      mostfracvars[i] = mostfracvars[i-1];
      mostfracvals[i] = mostfracvals[i-1];
   }
   assert(0 <= i && i <= *nflipcands && *nflipcands < maxnflipcands);

   /* insert the variable and its fractionality */
   mostfracvars[i] = var;
   mostfracvals[i] = frac;
   
   /* we've found another candidate */
   (*nflipcands)++;
}

/** flips the roundings of the most fractional variables, if a 1-cycle was found */
static
RETCODE handle1Cycle(
   SCIP*            scip,               /**< SCIP data structure  */
   HEURDATA*        heurdata,           /**< data of this special heuristic */
   VAR**            mostfracvars,       /**< sorted array of the currently most fractional variables */
   int              nflipcands,          /**< number of variables to flip */
   Real             alpha               /**< factor how much the original objective is regarded */
   )
{
   VAR* var;
   Real solval;
   Real frac;
   Real newobjcoeff;  
   Real orgobjcoeff;
   int  i;

   /* just flipping the objective coefficients from +1 to -1 and the rounding from floor to ceil */
   for( i = 0; i < nflipcands; i++ )
   {
      var = mostfracvars[i];
      solval = SCIPvarGetLPSol(var);   
      orgobjcoeff = SCIPvarGetObj(var);
      frac = SCIPfeasFrac(scip, solval);
      if( frac > 0.5 )
      {
         newobjcoeff = (1.0 - alpha) + alpha * orgobjcoeff;
         solval = SCIPfeasFloor(scip, solval);
      }         
      else
      {
         newobjcoeff = - (1.0 - alpha) + alpha * orgobjcoeff;
         solval = SCIPfeasCeil(scip, solval);
      }
      /* updating the rounded solution and the objective */
      CHECK_OKAY( SCIPsetSolVal(scip, heurdata->roundedsol, var, solval) );
      CHECK_OKAY( SCIPchgVarObjDive(scip, var, newobjcoeff) );
   }
   return SCIP_OKAY;
}

/** flips the roundings of the most fractional variables, if a 1-cycle was found */
static
RETCODE handleCycle(
   SCIP*            scip,               /**< SCIP data structure  */
   HEURDATA*        heurdata,           /**< data of this special heuristic */
   VAR**            vars,               /**< array of all variables */
   int              nbinandintvars,      /**< number of general integer and 0-1 variables */
   Real             alpha               /**< factor how much the original objective is regarded */
   )
{
   VAR* var;
   Real solval;
   Real frac;
   Real newobjcoeff;  
   Real orgobjcoeff;
   Real flipprob;
   int i;

   /* just flipping the objective coefficients from +1 to -1 and the rounding from floor to ceil */
   for( i = 0; i < nbinandintvars; i++ )
   {
      /* decide arbitraryly whether the variable will be flipped or not */
      var = vars[i];
      solval = SCIPvarGetLPSol(var);   
      orgobjcoeff = SCIPvarGetObj(var);
      frac = SCIPfeasFrac(scip, solval);
      flipprob = -0.3 + ((Real)(rand_r(&heurdata->randseed) % 1000000))/1000000.0 ;

      /* flip, iff the sum of the randomized number and the fractionality is big enough */
      if( MIN(frac,1.0-frac)+MAX(flipprob,0.0) > 0.5 )
      {
         if( frac > 0.5 )
         {
            newobjcoeff = (1.0 - alpha) + alpha * orgobjcoeff;
            solval = SCIPfeasFloor(scip, solval);
         }         
         else
         {
            newobjcoeff = - (1.0 - alpha) + alpha * orgobjcoeff;
            solval = SCIPfeasCeil(scip, solval);
         } 
         CHECK_OKAY( SCIPsetSolVal(scip, heurdata->roundedsol, var, solval) );
         CHECK_OKAY( SCIPchgVarObjDive(scip, var, newobjcoeff) );
      }
   }

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
DECL_HEURFREE(heurFreeFeaspump)
{   /*lint --e{715}*/
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
DECL_HEURINIT(heurInitFeaspump)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   CHECK_OKAY( SCIPcreateSol(scip, &heurdata->sol, heur) );
   CHECK_OKAY( SCIPcreateSol(scip, &heurdata->roundedsol, heur) );

   /* initialize data */
   heurdata->nlpiterations = 0;
   heurdata->randseed = 0;

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
DECL_HEUREXIT(heurExitFeaspump)
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   CHECK_OKAY( SCIPfreeSol(scip, &heurdata->sol) );
   CHECK_OKAY( SCIPfreeSol(scip, &heurdata->roundedsol) );

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolFeaspump NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolFeaspump NULL


/** execution method of primal heuristic */
static
DECL_HEUREXEC(heurExecFeaspump)
{  /*lint --e{715}*/
   HEURDATA* heurdata; 
   SOL* tmpsol;          /* only used for swapping */
   SOL** lastroundedsols;/* solutions of the last pumping rounds (depending on heurdata->cyclelength) */
   LPSOLSTAT lpsolstat;  /* status of the LP solution */

   VAR** vars;
   VAR* var;
   VAR** mostfracvars;   /* the 30 most fractional variables, needed to avoid 1-cycles */

   Real* mostfracvals;   /* the values of the variables above */
   Real newobjcoeff;     /* used for changing the objective */
   Real orgobjcoeff;     /* used for regarding the original objective */
   Real oldsolval;       /* one value of the last solution */ 
   Real solval;          /* one value of the actual solution */ 
   Real frac;            /* the fractional part of the value above */  
   Real objfactor;       /* factor by which the regard of the objective is decreased in each round, in [0,0.99] */
   Real alpha;           /* factor how the original objective is regarded, used for convex combination of two functions */
   Real objnorm;         /* Euclidean norm of the objective function, used for scaling */

   int nvars;            /* number of variables  */
   int nbinvars;         /* number of 0-1-variables */
   int nintvars;         /* number of integer variables */
   int nstartfracs;      /* number of fractional variables that should be integer at the beginning of the heuristic */
   int nfracs;           /* number of fractional variables updated after each pumping round*/
   int i;
   int j;
   int nflipcands;       /* how many flipcands (most frac. var.) have been found */
   int maxnflipcands;    /* maximal number of candidates to flip in the current pumping round */
   int nloops;           /* how many pumping rounds have been made */
   int maxflips;         /* maximum number of flips, if a 1-cycle is found (depending on heurdata->minflips) */ 
   int maxloops;         /* maximum number of pumping rounds */

   Longint nlpiterations;   /* number of LP iterations done during one pumping round */
   Longint maxnlpiterations; /* maximum number of LP iterations fpr this heuristic */
   Longint nsolsfound;   /* number of solutions found by this heuristic */
   Longint ncalls;       /* number of calls of this heuristic */  

   Bool success;         
   Bool lperror; 
   Bool* cycles;           /* are there short cycles */
   
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* get all variables of LP and number of fractional variables in LP solution that should be integral */
   CHECK_OKAY( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nfracs = SCIPgetNLPBranchCands(scip);
   assert(0 <= nfracs && nfracs <= nbinvars + nintvars);
   if( nfracs == 0 )
      return SCIP_OKAY;
   
   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNNodeLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = SCIPheurGetNSolsFound(heur);
   maxnlpiterations = (1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * nlpiterations;
   maxnlpiterations += heurdata->maxlpiterofs;
  
   /* initialize some heuristic data */
   maxflips = 3*heurdata->minflips;
   maxloops = heurdata->maxloops;
   if( maxloops == -1 )
      maxloops = INT_MAX;

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;
   
   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + 10000);
   
   debugMessage("executing feasibility pump heuristic, maxnlpit:%lld, maxflips:%d \n", maxnlpiterations, maxflips);

   *result = SCIP_DIDNOTFIND;

   /* memory allocation */
   CHECK_OKAY( SCIPallocBufferArray(scip, &mostfracvars, maxflips) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &mostfracvals, maxflips) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &lastroundedsols, heurdata->cyclelength) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &cycles, heurdata->cyclelength) );

   for( j = 0; j < heurdata->cyclelength; j++ )
   {
      CHECK_OKAY( SCIPcreateSol(scip, &lastroundedsols[j], heur) ); 
   }

   /* start diving */
   CHECK_OKAY( SCIPstartDive(scip) );
 
   /* lp was solved optimal */
   lperror = FALSE;
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;

   /* pumping rounds */
   objnorm = SCIPgetObjNorm(scip);
   objnorm = MAX(objnorm, 1.0);
   nsolsfound = SCIPgetNSolsFound(scip);
   if( heurdata->objfactor == 1.0 )
      objfactor = MIN(1.0 - 0.5 / (Real)(1 + nsolsfound), 0.99);
   else  
      objfactor = heurdata->objfactor;
   alpha = 1.0;
   nstartfracs = nfracs;
   nloops = 0;
   while( nfracs > 0 &&  heurdata->nlpiterations < maxnlpiterations && nloops < maxloops )
   {
      nloops++;
      alpha *= objfactor;

      debugMessage("feasibility pump loop %d: %d fractional variables\n", nloops, nfracs);

      /* create solution from diving LP and try to round it */
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      CHECK_OKAY( SCIProundSol(scip, heurdata->sol, &success) );  
         
      /* if the rounded solution is feasible and better, add it to SCIP */ 
      if( success )
      {
         CHECK_OKAY( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, &success) );
         if( success )
            *result = SCIP_FOUNDSOL; 
      }
      
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->roundedsol) );

      /* randomly choose maximum number of variables to flip in current pumping round in case of a 1-cycle */
      maxnflipcands = heurdata->minflips + rand_r(&heurdata->randseed) % (maxflips - heurdata->minflips + 1);
      nflipcands = 0;

      /* check, whether there is the possibility of j-cycling */
      for( j = 0; j <  heurdata->cyclelength; j++ )
         cycles[j] = (nloops > j+1);
         
      /* change objective function to Manhattan-distance of the integer variables to the LP and get the rounded solution */
      for( i = 0; i < nvars; i++ )
      {
         var = vars[i];
         solval = SCIPvarGetLPSol(var);
         orgobjcoeff = SCIPvarGetObj(var) * SQRT(nbinvars + nintvars) / objnorm;

         /* handle all integer variables*/
         if( i < nbinvars + nintvars )
         {  
            frac = SCIPfeasFrac(scip, solval);
            /* variables which are already integral, are treated separately */
            if( SCIPisFeasZero(scip,frac) )
            {
               /* variables at their bounds should be kept there */
               if( SCIPisFeasEQ(scip, solval, SCIPvarGetLbLocal(var)) )
                  newobjcoeff = (1.0 - alpha) + alpha * orgobjcoeff;
               else if( SCIPisFeasEQ(scip, solval, SCIPvarGetUbLocal(var)) )
                  newobjcoeff = - (1.0 - alpha) + alpha * orgobjcoeff;
               else
                  newobjcoeff = alpha * orgobjcoeff;
            }
            else 
            {
               /* check whether the variable is one of the most fractionals and label if so */
               if( cycles[0] )
                  insertFlipCand(mostfracvars, mostfracvals, &nflipcands, maxnflipcands, var, frac);
               if( frac > 0.5 )
               {
                  newobjcoeff = - (1.0 - alpha) + alpha * orgobjcoeff;
                  solval = SCIPfeasCeil(scip, solval);
               }            
               else
               {
                  newobjcoeff = (1.0 - alpha) + alpha * orgobjcoeff;
                  solval = SCIPfeasFloor(scip, solval);
               }

               /* update the rounded solution */
               CHECK_OKAY( SCIPsetSolVal(scip, heurdata->roundedsol, var, solval) );
            }
         }
         else
            newobjcoeff = alpha * orgobjcoeff;
         
         /* change one coefficient of the objective */
         CHECK_OKAY( SCIPchgVarObjDive(scip, var, newobjcoeff) );
         
         /* check, whether there is still the possibility of j-cycles */
         for( j = 0; j < MIN(heurdata->cyclelength, nloops-1); j++ ) 
         {
            /* cycles exist, iff all solution values are equal */
            if( cycles[j] )
            {
               oldsolval = SCIPgetSolVal(scip, lastroundedsols[j], var);
               cycles[j] = SCIPisFeasEQ(scip, solval, oldsolval);
            }
         }
      }
  
      /* force to flip variables at random after a couple of pumping rounds */
      assert(heurdata->perturbfreq > 0);
      if( nloops % heurdata->perturbfreq == 0 )
      {
         debugMessage(" -> random perturbation\n");
         CHECK_OKAY( handleCycle(scip, heurdata, vars, nintvars+nbinvars, alpha) );
      }
      else 
      {
         for( j = 0; j < MIN(heurdata->cyclelength, nloops-1); j++ ) 
         {
            /* if we got the same rounded solution as in some step before, we have to flip some variables */
            if( cycles[j] )
            {
               /* 1-cycles have a special flipping rule (flip most fractional variables) */
               if( j == 0 )
               {
                  debugMessage(" -> avoiding 1-cycle: flipping %d candidates\n", nflipcands);
                  CHECK_OKAY( handle1Cycle(scip, heurdata, mostfracvars, nflipcands, alpha) );
               }
               else 
               {
                  debugMessage(" -> avoiding %d-cycle by random flip\n", j+1);
                  CHECK_OKAY( handleCycle(scip, heurdata, vars, nintvars+nbinvars, alpha) );
               }
               break;
            }
         }
      }
      
      /* the LP with the new (distance) objective is solved */
      nlpiterations = SCIPgetNLPIterations(scip);
      CHECK_OKAY( SCIPsolveDiveLP(scip, maxnlpiterations , &lperror) );
      lpsolstat = SCIPgetLPSolstat(scip);

      /* check whether LP was solved optimal */
      if( lperror || lpsolstat != SCIP_LPSOLSTAT_OPTIMAL )
         break; 
      
      /* update iteration count */
      heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;
      nfracs = SCIPgetNLPBranchCands(scip);
      debugMessage(" -> number of iterations: %lld/%lld\n", heurdata->nlpiterations, maxnlpiterations);

      /* swap the last solutions */
      tmpsol = lastroundedsols[heurdata->cyclelength-1];
      for( j = heurdata->cyclelength-1; j > 0; j-- )
         lastroundedsols[j] = lastroundedsols[j-1]; 
      lastroundedsols[0] = heurdata->roundedsol;
      heurdata->roundedsol = tmpsol;
   }

   /* try final solution, if no more fractional variables are left */
   if( nfracs == 0 && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      CHECK_OKAY( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, &success) );
      if( success )
         *result = SCIP_FOUNDSOL;  
   }

   /* end diving */
   CHECK_OKAY( SCIPendDive(scip) );

   /* free memory */
   for( j = 0; j < heurdata->cyclelength; j++ )
   {
      CHECK_OKAY( SCIPfreeSol(scip, &lastroundedsols[j]) );
   }
   SCIPfreeBufferArray(scip, &cycles);
   SCIPfreeBufferArray(scip, &lastroundedsols);
   SCIPfreeBufferArray(scip, &mostfracvals);
   SCIPfreeBufferArray(scip, &mostfracvars);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the feaspump primal heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurFeaspump(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create feaspump primal heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING, HEUR_AFTERNODE,
         heurFreeFeaspump, heurInitFeaspump, heurExitFeaspump, 
         heurInitsolFeaspump, heurExitsolFeaspump, heurExecFeaspump,
         heurdata) );

   /* add feaspump primal heuristic parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/feaspump/maxlpiterquot", 
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, DEFAULT_MAXLPITERQUOT, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "heuristics/feaspump/maxlpiterofs", 
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/feaspump/objfactor", 
         "factor by which the regard of the objective is decreased in each round, 1.0 for dynamic",
         &heurdata->objfactor, DEFAULT_OBJFACTOR, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "heuristics/feaspump/maxloops",
         "maximal number of pumping loops (-1: no limit)",
         &heurdata->maxloops, DEFAULT_MAXLOOPS, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "heuristics/feaspump/minflips", 
         "minimum number of random variables to flip, if a 1-cycle is encountered",
         &heurdata->minflips, DEFAULT_MINFLIPS, 1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "heuristics/feaspump/cyclelength", 
         "maximum length of cycles to be checked explicitly in each round",
         &heurdata->cyclelength, DEFAULT_CYCLELENGTH, 1, 100, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "heuristics/feaspump/perturbfreq", 
         "number of iterations until a random perturbation is forced",
         &heurdata->perturbfreq, DEFAULT_PERTURBFREQ, 1, INT_MAX, NULL, NULL) );
   return SCIP_OKAY;
}
