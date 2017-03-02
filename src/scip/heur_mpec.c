/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_mpec.c
 * @brief  mpec primal heuristic
 * @author Felipe Serrano
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_mpec.h"
#include "scip/heur_subnlp.h"
#include "nlpi/nlpi.h"


#define HEUR_NAME             "mpec"
#define HEUR_DESC             "regularization heuristic for convex and nonconvex MINLPs"
#define HEUR_DISPCHAR         'W'
#define HEUR_PRIORITY         -2050000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_NLPI* nlpi;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_HASHMAP* var2idx;
};


/*
 * Local methods
 */

static
SCIP_RETCODE createNLP(
   SCIP* scip,
   SCIP_HEURDATA* heurdata
   )
{
   assert(heurdata != NULL);

   if( heurdata->nlpiprob != NULL )
      return SCIP_OKAY;

   assert(SCIPgetNNlpis(scip) > 0);
   heurdata->nlpi = SCIPgetNlpis(scip)[0];
   assert(heurdata->nlpi != NULL);

   SCIP_CALL( SCIPnlpiCreateProblem(heurdata->nlpi, &heurdata->nlpiprob, "MPEC-nlp") );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->var2idx, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   SCIP_CALL( SCIPcreateNlpiProb(scip, heurdata->nlpi, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip),
         heurdata->nlpiprob, heurdata->var2idx, NULL, SCIPinfinity(scip), TRUE, FALSE) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE freeNLP(
   SCIP* scip,
   SCIP_HEURDATA* heurdata
   )
{
   assert(heurdata != NULL);

   if( heurdata->nlpiprob == NULL )
      return SCIP_OKAY;

   assert(heurdata->var2idx != NULL);

   SCIPhashmapFree(&heurdata->var2idx);
   SCIP_CALL( SCIPnlpiFreeProblem(heurdata->nlpi, &heurdata->nlpiprob) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE addRegularScholtes(
   SCIP* scip,
   SCIP_HEURDATA* heurdata,
   SCIP_VAR** binvars,
   int nbinvars,
   SCIP_Real theta,
   SCIP_Bool update
   )
{
   int i;

   assert(binvars != NULL);
   assert(nbinvars > 0);

   /* add or update regularization for each non-fixed binary variables */
   if( !update )
   {
      SCIP_QUADELEM* quadelems;
      SCIP_Real* linvals;
      int* lininds;

      SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds, 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals, 1) );

      for( i = 0; i < nbinvars; ++i )
      {
         SCIP_VAR* var = binvars[i];
         SCIP_Real lhs = -SCIPinfinity(scip);
         SCIP_Real rhs = theta;
         int nlininds = 1;
         int nquadelems = 1;
         int idx;

         assert(var != NULL);
         assert(heurdata->var2idx != NULL);
         assert(SCIPhashmapExists(heurdata->var2idx, (void*)var));
         idx = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)var);

         lininds[0] = idx;
         linvals[0] = 1.0;
         quadelems->idx1 = lininds[0];
         quadelems->idx2 = lininds[0];
         quadelems->coef = -1.0;

         SCIP_CALL( SCIPnlpiAddConstraints(heurdata->nlpi, heurdata->nlpiprob, 1, &lhs, &rhs, &nlininds,
               &lininds, &linvals, &nquadelems, &quadelems, NULL, NULL, NULL) );
      }

      SCIPfreeBufferArray(scip, &linvals);
      SCIPfreeBufferArray(scip, &lininds);
      SCIPfreeBufferArray(scip, &quadelems);
   }
   else
   {
      int startidx = SCIPgetNNLPNlRows(scip) + 1; /* the cutoff is a separate constraint */
      SCIP_Real* lhss;
      SCIP_Real* rhss;
      int* indices;

      SCIP_CALL( SCIPallocBufferArray(scip, &lhss, nbinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nbinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &indices, nbinvars) );

      for( i = 0; i < nbinvars; ++i )
      {
         lhss[i] = -SCIPinfinity(scip);
         rhss[i] = theta;
         indices[i] = startidx + i;
      }

      SCIP_CALL( SCIPnlpiChgConsSides(heurdata->nlpi, heurdata->nlpiprob, nbinvars, indices, lhss, rhss) );

      SCIPfreeBufferArray(scip, &indices);
      SCIPfreeBufferArray(scip, &rhss);
      SCIPfreeBufferArray(scip, &lhss);
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE heurExec(
   SCIP* scip,
   SCIP_HEUR* heur,
   SCIP_HEURDATA* heurdata,
   SCIP_RESULT* result
   )
{
   SCIP_NLPSTATISTICS* nlpstatistics = NULL;
   SCIP_VAR** binvars = NULL;
   SCIP_Real* initguess = NULL;
   SCIP_Real* ubs = NULL;
   SCIP_Real* lbs = NULL;
   int* indices = NULL;
   SCIP_Real timelim;
   SCIP_Bool reinit = TRUE;
   SCIP_Bool fixed = FALSE;
   SCIP_Bool subnlpcalled = FALSE;
   int nbinvars = 0;
   int i;

   SCIP_Real theta = 1.0 / 8.0;
   SCIP_Real sigma = 0.70;

   assert(heurdata != NULL);
   assert(heurdata->nlpiprob != NULL);
   assert(heurdata->var2idx != NULL);
   assert(heurdata->nlpi != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, SCIPgetNBinVars(scip)) );
   SCIP_CALL( SCIPnlpStatisticsCreate(&nlpstatistics) );

   /* collect all non-fixed binary variables */
   for( i = 0; i < SCIPgetNBinVars(scip); ++i )
   {
      SCIP_VAR* var = SCIPgetVars(scip)[i];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

      if( !SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         binvars[nbinvars++] = var;
   }

   /* all binary variables are fixed */
   if( nbinvars == 0 )
      goto TERMINATE;

   SCIP_CALL( SCIPallocBufferArray(scip, &initguess, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, nbinvars) );

   /* set initial guess */
   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      SCIP_VAR* var = SCIPgetVars(scip)[i];
      initguess[i] = SCIPgetSolVal(scip, NULL, var);
      SCIPdebugMsg(scip, "set initial value for %s to %g\n", SCIPvarGetName(var), initguess[i]);
   }
   SCIP_CALL( SCIPnlpiSetInitialGuess(heurdata->nlpi, heurdata->nlpiprob, initguess, NULL, NULL, NULL) );

   /* set parameters of NLP solver */
   SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_FEASTOL,
         SCIPfeastol(scip) / 10.0) );
   SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_RELOBJTOL,
         SCIPdualfeastol(scip) / 10.0) );
   SCIP_CALL( SCIPnlpiSetIntPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_VERBLEVEL, 0) );
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelim) );

   /* main loop */
   for( i = 0; i < 100 && *result != SCIP_FOUNDSOL; ++i )
   {
      SCIP_Real* primal = NULL;
      SCIP_Real timeleft = SCIPinfinity(scip);
      SCIP_Bool binaryfeasible;
      SCIP_Bool regularfeasible;
      SCIP_NLPSOLSTAT solstat;
      SCIP_Real maxviolbin = 0.0;
      SCIP_Real maxviol = 0.0;
      int iterleft = 500;
      int j;

      /* add or update regularization */
      SCIP_CALL( addRegularScholtes(scip, heurdata, binvars, nbinvars, theta, i > 0) );

      /* set working limits */
      if( !SCIPisInfinity(scip, timelim) )
      {
         timeleft = timelim - SCIPgetSolvingTime(scip);
         if( timeleft <= 0.0 )
         {
            SCIPdebugMsg(scip, "skip NLP solve; no time left\n");
            break;
         }
      }
      SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_TILIM, timeleft) );
      SCIP_CALL( SCIPnlpiSetIntPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_ITLIM, iterleft) );

      /* solve NLP */
      SCIP_CALL( SCIPnlpiSolve(heurdata->nlpi, heurdata->nlpiprob) );
      SCIP_CALL( SCIPnlpiGetStatistics(heurdata->nlpi, heurdata->nlpiprob, nlpstatistics) );
      solstat = SCIPnlpiGetSolstat(heurdata->nlpi, heurdata->nlpiprob);

      /* give up if an error occurred or no primal values are accessible */
      if( solstat > SCIP_NLPSOLSTAT_LOCINFEASIBLE )
      {
         SCIPdebugMsg(scip, "error occured during NLP solve -> stop!\n");
         break;
      }

      SCIP_CALL( SCIPnlpiGetSolution(heurdata->nlpi, heurdata->nlpiprob, &primal, NULL, NULL, NULL) );
      assert(primal != NULL);

      for( j = 0; j < SCIPgetNVars(scip); ++j )
      {
         SCIP_VAR* var = SCIPgetVars(scip)[j];
         SCIPdebugMsg(scip, "NLP sol for %s = %g\n", SCIPvarGetName(var), primal[j]);
      }

      /* check for binary feasibility */
      binaryfeasible = TRUE;
      regularfeasible = TRUE;
      for( j = 0; j < nbinvars; ++j )
      {
         int idx = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)binvars[j]);
         binaryfeasible = binaryfeasible && SCIPisFeasIntegral(scip, primal[idx]);
         regularfeasible = regularfeasible && SCIPisLE(scip, primal[idx] - SQR(primal[idx]), theta);

         maxviol = MAX(maxviol, primal[idx] - SQR(primal[idx]) - theta);
         maxviolbin = MAX(maxviolbin, MIN(primal[idx], 1.0-primal[idx]));
      }
      SCIPdebugMsg(scip, "maxviol-regularization %g maxviol-integrality %g\n", maxviol, maxviolbin);

      /* call sub-NLP heursitic when the maximum binary infeasibility is small enough */
      if( !subnlpcalled && SCIPisLE(scip, maxviolbin, 1e-3) )
      {
         SCIP_SOL* refpoint;
         SCIP_RESULT subnlpresult;

         SCIPdebugMsg(scip, "call sub-NLP heuristic because binary infeasibility is small enough\n");
         SCIP_CALL( SCIPcreateSol(scip, &refpoint, heur) );

         for( j = 0; j < SCIPgetNVars(scip); ++j )
         {
            SCIP_VAR* var = SCIPgetVars(scip)[j];
            SCIP_Real val = SCIPvarIsBinary(var) ? SCIPfeasRound(scip, primal[j]) : primal[j];
            SCIP_CALL( SCIPsetSolVal(scip, refpoint, var, val) );
         }

         SCIP_CALL( SCIPapplyHeurSubNlp(scip, SCIPfindHeur(scip, "subnlp"), &subnlpresult, refpoint, -1LL, timeleft,
               0.0, NULL, NULL) );
         SCIP_CALL( SCIPfreeSol(scip, &refpoint) );
         SCIPdebugMsg(scip, "result of sub-NLP call: %d\n", subnlpresult);

         if( subnlpresult == SCIP_FOUNDSOL )
         {
            SCIPdebugMsg(scip, "sub-NLP found a feasible solution -> stop!\n");
            break;
         }

         subnlpcalled = TRUE;
      }

      /* NLP feasible + binary feasible -> add solution and stop */
      if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE && binaryfeasible )
      {
         SCIP_SOL* sol;
         SCIP_Bool stored;

         SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

         for( j = 0; j < SCIPgetNVars(scip); ++j )
         {
            SCIP_VAR* var = SCIPgetVars(scip)[j];
            assert(j == (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)var));
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, primal[j]) );
         }

         SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, FALSE, &stored) );
         SCIPdebugMsg(scip, "found a solution (stored = %u)\n", stored);

         if( stored )
            *result = SCIP_FOUNDSOL;
         break;
      }

      /* NLP feasible + binary infeasible -> reduce theta */
      else if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE && !binaryfeasible )
      {
         BMScopyMemoryArray(initguess, primal, SCIPgetNVars(scip));
         SCIP_CALL( SCIPnlpiSetInitialGuess(heurdata->nlpi, heurdata->nlpiprob, primal, NULL, NULL, NULL) );
         SCIPdebugMsg(scip, "update theta from %g -> %g\n", theta, theta*sigma);

         if( !reinit )
         {
            SCIPdebugMsg(scip, "reinit fixed the infeasibility\n");
            reinit = TRUE;
         }

         theta *= sigma;

         /* unfix binary variables */
         if( fixed )
         {
            SCIPdebugMsg(scip, "unfixing binary variables\n");
            for( j = 0; j < nbinvars; ++j )
            {
               lbs[j] = 0.0;
               ubs[j] = 1.0;
               indices[j] = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)binvars[j]);
            }
            SCIP_CALL( SCIPnlpiChgVarBounds(heurdata->nlpi, heurdata->nlpiprob, nbinvars, indices, lbs, ubs) );
            fixed = FALSE;
         }
      }

      /* NLP infeasible + regularization feasible -> stop (give up) */
      else if( solstat > SCIP_NLPSOLSTAT_FEASIBLE && regularfeasible )
      {
         SCIPdebugMsg(scip, "NLP is infeasible but regularization constraints are satisfied -> stop!\n");
         break;
      }

      /* NLP infeasible + binary infeasible -> set initial point / fix binary variables */
      else
      {
         assert(solstat > SCIP_NLPSOLSTAT_FEASIBLE && !regularfeasible);

         SCIPdebugMsg(scip, "NLP solution is not feasible for the NLP and the binary variables\n");

         /* fix variables if reinit is FALSE; otherwise set another initial point */
         if( !reinit )
         {
            /* fix binary variables */
            for( j = 0; j < nbinvars; ++j )
            {
               int idx = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)binvars[j]);
               indices[j] = idx;

               if( SCIPisFeasLE(scip, primal[idx] - SQR(primal[idx]), theta) )
               {
                  lbs[j] = 0.0;
                  ubs[j] = 1.0;
               }
               else
               {
                  lbs[j] = primal[idx] >= 0.5 ? 0.0 : 1.0;
                  ubs[j] = primal[idx] >= 0.5 ? 0.0 : 1.0;
                  SCIPdebugMsg(scip, "fix binary variable %s = %g\n", SCIPvarGetName(binvars[j]), ubs[j]);
               }
            }
            SCIP_CALL( SCIPnlpiChgVarBounds(heurdata->nlpi, heurdata->nlpiprob, nbinvars, indices, lbs, ubs) );
            fixed = TRUE;
         }
         else
         {
            /* set initial point */
            for( j = 0; j < nbinvars; ++j )
            {
               int idx = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)binvars[j]);
               initguess[idx] = primal[idx] >= 0.5 ? 0.0 : 1.0;
               SCIPdebugMsg(scip, "update init guess for %s to %g\n", SCIPvarGetName(binvars[j]), initguess[idx]);
            }
            SCIP_CALL( SCIPnlpiSetInitialGuess(heurdata->nlpi, heurdata->nlpiprob, initguess, NULL, NULL, NULL) );
            reinit = FALSE;
         }
      }
   }

TERMINATE:
   SCIPfreeBufferArrayNull(scip, &indices);
   SCIPfreeBufferArrayNull(scip, &ubs);
   SCIPfreeBufferArrayNull(scip, &lbs);
   SCIPfreeBufferArrayNull(scip, &initguess);
   SCIPnlpStatisticsFree(&nlpstatistics);
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyMpec)
{  /*lint --e{715}*/
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurMpec(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolMpec)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolMpec)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNIntVars(scip) > 0 || SCIPgetNBinVars(scip) == 0 )
      /* || SCIPgetNNlpis(scip) == 0 || !SCIPisNLPConstructed(scip) ) */
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* create NLP */
   SCIP_CALL( createNLP(scip, heurdata) );

   SCIP_CALL( heurExec(scip, heur, heurdata, result) );

   /* free NLP */
   SCIP_CALL( freeNLP(scip, heurdata) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the mpec primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMpec(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata = NULL;
   SCIP_HEUR* heur = NULL;

   /* create mpec primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecMpec, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyMpec) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeMpec) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolMpec) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolMpec) );

   /* add mpec primal heuristic parameters */

   return SCIP_OKAY;
}
