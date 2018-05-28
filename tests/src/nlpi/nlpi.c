/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlpi.c
 * @brief  unit test for nonlinear interface methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/nlpi_ipopt.h"
#include "nlpi/nlpi_worhp.h"
#include "nlpi/nlpi_all.h"
#include "nlpi/nlpioracle.h"
#include "nlpi/nlpi.h"

#include "include/scip_test.h"

#define INF 1e+20

static SCIP* scip = NULL;
static SCIP_NLPI* ipopt = NULL;
static SCIP_NLPI* worhpip = NULL;
static SCIP_NLPI* worhpsqp = NULL;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* check whether WORHP is available */
   if( SCIPisWorhpAvailableWorhp() )
   {
      /* WORHP-IP */
      SCIP_CALL( SCIPcreateNlpSolverWorhp(SCIPblkmem(scip), &worhpip, TRUE) );
      cr_assert(worhpip != NULL);
      SCIP_CALL( SCIPincludeNlpi(scip, worhpip) );

      /* WORHP-SQP */
      SCIP_CALL( SCIPcreateNlpSolverWorhp(SCIPblkmem(scip), &worhpsqp, FALSE) );
      cr_assert(worhpsqp != NULL);
      SCIP_CALL( SCIPincludeNlpi(scip, worhpsqp) );

      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameWorhp(), SCIPgetSolverDescWorhp()) );
   }

   /* check whether Ipopt is available */
   if( SCIPisIpoptAvailableIpopt() )
   {
      SCIP_CALL( SCIPcreateNlpSolverIpopt(SCIPblkmem(scip), &ipopt) );
      cr_assert(ipopt != NULL);
      SCIP_CALL( SCIPincludeNlpi(scip, ipopt) );
      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameIpopt(), SCIPgetSolverDescIpopt()) );
   }

   if( SCIPgetNNlpis(scip) > 1 )
   {
      SCIP_NLPI* nlpi;

      SCIP_CALL( SCIPcreateNlpSolverAll(SCIPblkmem(scip), &nlpi, SCIPgetNlpis(scip), SCIPgetNNlpis(scip)) );
      SCIP_CALL( SCIPincludeNlpi(scip, nlpi) );
      cr_assert(nlpi != NULL);
   }
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

/* helper function to test NLPI */
static
SCIP_RETCODE testNlpi(SCIP_NLPI* nlpi)
{
   SCIP_NLPSTATISTICS* statistics;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_Real* dualcons;
   SCIP_Real* duallb;
   SCIP_Real* dualub;
   SCIP_Real* primal;

   /* variables */
   const char* varnames[4] = {"x0", "x1", "x2", "x3"};
   SCIP_Real lbs[4] = {-0.5, -2, 0, -2};
   SCIP_Real ubs[4] = {INF, INF, 2, 2};

   /* objective */
   SCIP_Real objvals[1] = {-1};
   int objinds[1] = {2};

   /* constraints */
   const char* consnames[3] = {"c1", "c2", "c3"};
   SCIP_Real lhss[3] = {1, -INF, 2.5};
   SCIP_Real rhss[3] = {1, -1, 5};
   SCIP_Real initguess[4];
   SCIP_QUADELEM* quadelems;
   SCIP_Real* linvals;
   int* lininds;
   int nquadelems;
   int nlinds;

   /*-----------------------------------------------------------------------
    *
    * min             f(x) = x0^2 + 2 x1^2 - x2
    *
    * subject to      -0.5 <= x0 <= INF
    *                   -2 <= x1 <= INF
    *                    0 <= x2 <= 2
    *                   -2 <= x3 <= 2
    *                    1 <= x0^2 + x2^2 + x0x2 <= 1
    *                 -INF <= x2 - x3 <= -1
    *                  2.5 <= x1 + x3 <= 5
    *
    * optimal solution
    *                     x*  = (0, 0.5, 1, 2)
    *                   f(x*) = -0.5
    *
    *-----------------------------------------------------------------------*/

   /* create NLPI problem */
   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, 2) );
   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, &nlpiprob, "convex_NLP") );
   SCIP_CALL( SCIPnlpiAddVars(nlpi, nlpiprob, 4, lbs, ubs, varnames) );

   /* set objective */
   quadelems[0].idx1 = 0; quadelems[0].idx2 = 0; quadelems[0].coef = 1;
   quadelems[1].idx1 = 1; quadelems[1].idx2 = 1; quadelems[1].coef = 2;
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 1, objinds, objvals, 2, quadelems, NULL, NULL, 0.0) );

   /* constraints */
   nquadelems = 3;
   quadelems[0].idx1 = 0; quadelems[0].idx2 = 0; quadelems[0].coef = 1;
   quadelems[1].idx1 = 2; quadelems[1].idx2 = 2; quadelems[1].coef = 1;
   quadelems[2].idx1 = 0; quadelems[2].idx2 = 2; quadelems[2].coef = 1;
   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, 1, &lhss[0], &rhss[0], NULL, NULL, NULL, &nquadelems, &quadelems,
         NULL, NULL, &consnames[0]) );

   nlinds = 2;
   lininds[0] = 2;
   linvals[0] = 1;
   lininds[1] = 3;
   linvals[1] = -1;
   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, 1, &lhss[1], &rhss[1], &nlinds, &lininds, &linvals, NULL, NULL,
         NULL, NULL, &consnames[1]) );

   nlinds = 2;
   lininds[0] = 1;
   linvals[0] = 1;
   lininds[1] = 3;
   linvals[1] = 1;
   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, 1, &lhss[2], &rhss[2], &nlinds, &lininds, &linvals, NULL, NULL,
         NULL, NULL, &consnames[2]) );

   /* solve NLP */
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_FEASTOL, 1e-9) );
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   cr_expect(SCIPnlpiGetTermstat(nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPnlpiGetSolstat(nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);

   /* collect statistics */
   SCIP_CALL( SCIPnlpStatisticsCreate(SCIPblkmem(scip), &statistics) );
   SCIP_CALL( SCIPnlpiGetStatistics(nlpi, nlpiprob, statistics) );
   cr_expect(SCIPnlpStatisticsGetNIterations(statistics) > 0);
   cr_expect(SCIPnlpStatisticsGetTotalTime(statistics) >= 0.0);
   SCIPnlpStatisticsFree(SCIPblkmem(scip), &statistics);

   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, &dualcons, &duallb, &dualub, NULL) );

   /* check primal solution */
   cr_expect(SCIPisFeasEQ(scip, primal[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, primal[1], 0.5));
   cr_expect(SCIPisFeasEQ(scip, primal[2], 1.0));
   cr_expect(SCIPisFeasEQ(scip, primal[3], 2.0));

   /* check dual solution */
   cr_expect(SCIPisFeasEQ(scip, dualcons[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, dualcons[1], 1.0));
   cr_expect(SCIPisFeasEQ(scip, dualcons[2], -2.0));
   cr_expect(SCIPisFeasEQ(scip, duallb[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, duallb[1], 0.0));
   cr_expect(SCIPisFeasEQ(scip, duallb[2], 0.0));
   cr_expect(SCIPisFeasEQ(scip, duallb[3], 0.0));
   cr_expect(SCIPisFeasEQ(scip, dualub[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, dualub[1], 0.0));
   cr_expect(SCIPisFeasEQ(scip, dualub[2], 0.0));
   cr_expect(SCIPisFeasEQ(scip, dualub[3], 3.0));

   /* change upper bound of x2 */
   lininds[0] = 2;
   lbs[0] = 0.0;
   ubs[0] = 0.5;
   SCIP_CALL( SCIPnlpiChgVarBounds(nlpi, nlpiprob, 1, lininds, lbs, ubs) );

   /* set the initial guess to the previous solution */
   initguess[0] = 0.6;
   initguess[1] = 0.5;
   initguess[2] = 0.4;
   initguess[3] = 2.0;
   SCIP_CALL( SCIPnlpiSetInitialGuess(nlpi, nlpiprob, initguess, NULL, NULL, NULL) );

   /* solve NLP */
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_FEASTOL, 1e-9) );
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   cr_expect(SCIPnlpiGetTermstat(nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPnlpiGetSolstat(nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);

   /* collect statistics */
   SCIP_CALL( SCIPnlpStatisticsCreate(SCIPblkmem(scip), &statistics) );
   SCIP_CALL( SCIPnlpiGetStatistics(nlpi, nlpiprob, statistics) );
   cr_expect(SCIPnlpStatisticsGetNIterations(statistics) > 0);
   cr_expect(SCIPnlpStatisticsGetTotalTime(statistics) >= 0.0);
   SCIPnlpStatisticsFree(SCIPblkmem(scip), &statistics);

   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, &dualcons, &duallb, &dualub, NULL) );

   /* check primal solution */
   cr_expect(SCIPisFeasEQ(scip, primal[0], 0.6513878189));
   cr_expect(SCIPisFeasEQ(scip, primal[1], 0.5));
   cr_expect(SCIPisFeasEQ(scip, primal[2], 0.5));
   cr_expect(SCIPisFeasEQ(scip, primal[3], 2.0));

   /* check dual solution */
   cr_expect(SCIPisFeasEQ(scip, dualcons[0], -0.7226499019));
   cr_expect(SCIPisFeasEQ(scip, dualcons[1], 0.0));
   cr_expect(SCIPisFeasEQ(scip, dualcons[2], -2.0));
   cr_expect(SCIPisFeasEQ(scip, duallb[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, duallb[1], 0.0));
   cr_expect(SCIPisFeasEQ(scip, duallb[2], 0.0));
   cr_expect(SCIPisFeasEQ(scip, duallb[3], 0.0));
   cr_expect(SCIPisFeasEQ(scip, dualub[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, dualub[1], 0.0));
   cr_expect(SCIPisFeasEQ(scip, dualub[2], 2.1933752453));
   cr_expect(SCIPisFeasEQ(scip, dualub[3], 2.0));

   /* free memory */
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &quadelems);
   SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &nlpiprob) );

   return SCIP_OKAY;
}

/* helper function to solve a convex QP */
static
SCIP_RETCODE solveQP(
   SCIP_NLPI* nlpi,
   unsigned int rndseed,
   int         n,
   SCIP_Real   minlb,
   SCIP_Real   maxub,
   SCIP_Real   timelim,
   int         iterlim,
   SCIP_Real*  solval,
   SCIP_NLPSOLSTAT* nlpsolstat,
   SCIP_NLPTERMSTAT* nlptermstat
   )
{
   SCIP_NLPSTATISTICS* statistics;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_QUADELEM* quadelems;
   SCIP_Real* primal;
   SCIP_Real* vals;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real objval;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int* inds;
   int nquadelems;
   int nlininds;
   int objind;
   int i;
   int k;


   *solval = SCIP_INVALID;
   *nlpsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
   *nlptermstat = SCIP_NLPTERMSTAT_OTHER;

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, rndseed, TRUE) );

   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, n+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, n+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, n*n) );

   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, &nlpiprob, "QP") );

   /* create variables */
   for( i = 0; i < n; ++i )
   {
      lbs[i] = SCIPrandomGetReal(randnumgen, minlb, maxub);
      ubs[i] = SCIPrandomGetReal(randnumgen, lbs[i], maxub);
   }
   lbs[n] = -SCIPinfinity(scip);
   ubs[n] = SCIPinfinity(scip);
   SCIP_CALL( SCIPnlpiAddVars(nlpi, nlpiprob, n+1, lbs, ubs, NULL) );

   /* set objective */
   objind = n;
   objval = 1.0;
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 1, &objind, &objval, 0, NULL, NULL, NULL, 0.0) );

   /* create constraint */
   for( i = 0; i < n; ++i )
      vals[i] = SCIPrandomGetReal(randnumgen, -10.0, 10.0);

   k = 0;
   for( i = 0; i < n; ++i )
   {
      int j;

      for( j = 0; j <= i; ++j )
      {
         quadelems[k].coef = (j < i) ? 2.0 * vals[i] * vals[j] : vals[i] * vals[j];
         quadelems[k].idx1 = j;
         quadelems[k].idx2 = i;
         ++k;
      }
   }
   assert(k <= n*n);

   nquadelems = k;
   nlininds = 1;
   inds[0] = n;
   vals[0] = -1.0;
   lhs = 0.0;
   rhs = 0.0;
   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, 1, &lhs, &rhs, &nlininds, &inds, &vals, &nquadelems, &quadelems,
         NULL, NULL, NULL) );

   /* set parameters */
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_FEASTOL, 1e-9) );
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_TILIM, timelim) );
   SCIP_CALL( SCIPnlpiSetIntPar(nlpi, nlpiprob, SCIP_NLPPAR_ITLIM, iterlim) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPnlpiSetIntPar(nlpi, nlpiprob, SCIP_NLPPAR_VERBLEVEL, 1) );
#endif

   /* solve NLP */
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );

   *solval = primal[n];
   *nlpsolstat = SCIPnlpiGetSolstat(nlpi, nlpiprob);
   *nlptermstat = SCIPnlpiGetTermstat(nlpi, nlpiprob);

   /* collect statistics */
   SCIP_CALL( SCIPnlpStatisticsCreate(SCIPblkmem(scip), &statistics) );
   SCIP_CALL( SCIPnlpiGetStatistics(nlpi, nlpiprob, statistics) );
   SCIPnlpStatisticsFree(SCIPblkmem(scip), &statistics);

   /* free memory */
   SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &nlpiprob) );
   SCIPfreeBufferArray(scip, &quadelems);
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);
   SCIPfreeRandom(scip, &randnumgen);

   return SCIP_OKAY;
}

Test(nlpi, interface, .init = setup, .fini = teardown,
   .description = "checks fundamental NLPI functions"
   )
{
   int i;

   for( i = 0; i < SCIPgetNNlpis(scip); ++i )
   {
      SCIP_CALL( testNlpi(SCIPgetNlpis(scip)[i]) );
   }
}

Test(nlpi, solveQP, .init = setup, .fini = teardown,
   .description = "solves convex QP with different NLPIs"
   )
{
   SCIP_NLPTERMSTAT ipopttermstat;
   SCIP_NLPTERMSTAT worhpiptermstat;
   SCIP_NLPTERMSTAT worhpsqptermstat;
   SCIP_NLPSOLSTAT ipoptsolstat;
   SCIP_NLPSOLSTAT worhpipsolstat;
   SCIP_NLPSOLSTAT worhpsqpsolstat;
   SCIP_Real ipoptval;
   SCIP_Real worhpipval;
   SCIP_Real worhpsqpval;
   int i;

   for( i = 0; i < 10; ++i )
   {
      /* solve QP with Ipopt */
      if( ipopt != NULL )
      {
         SCIP_CALL( solveQP(ipopt, i+1, 100, -100.0, 100.0, SCIPinfinity(scip), INT_MAX, &ipoptval,
               &ipoptsolstat, &ipopttermstat) );
      }

      /* solve QP with WORHP-IP */
      if( worhpip != NULL )
      {
         SCIP_CALL( solveQP(worhpip, i+1, 100, -100.0, 100.0, SCIPinfinity(scip), INT_MAX, &worhpipval,
               &worhpipsolstat, &worhpiptermstat) );
      }

      /* solve QP with WORHP-SQP */
      if( worhpsqp != NULL )
      {
         SCIP_CALL( solveQP(worhpsqp, i+1, 100, -100.0, 100.0, SCIPinfinity(scip), INT_MAX, &worhpsqpval,
               &worhpsqpsolstat, &worhpsqptermstat) );
      }

      /* compare the solution values of WORHP and Ipopt */
      if( ipopt != NULL && worhpip != NULL &&
         ipoptsolstat == SCIP_NLPSOLSTAT_LOCOPT && worhpipsolstat == SCIP_NLPSOLSTAT_LOCOPT )
      {
         cr_assert(SCIPisFeasLE(scip, worhpipval, ipoptval));
      }

      /* compare the solution values of WORHP and Ipopt */
      if( ipopt != NULL && worhpsqp != NULL &&
         ipoptsolstat == SCIP_NLPSOLSTAT_LOCOPT && worhpsqpsolstat == SCIP_NLPSOLSTAT_LOCOPT )
      {
         cr_assert(SCIPisFeasLE(scip, worhpsqpval, ipoptval));
      }
   }
}

Test(nlpi, workinglimits, .init = setup, .fini = teardown,
   .description = "solves convex QP with a small iteration limit"
   )
{
   SCIP_NLPTERMSTAT termstat;
   SCIP_NLPSOLSTAT solstat;
   SCIP_Real solval;

   if( worhpip != NULL )
   {
      /* set a small iteration limit */
      SCIP_CALL( solveQP(worhpip, 1, 100, -100.0, 100.0, SCIPinfinity(scip), 5, &solval, &solstat, &termstat) );
      cr_expect(termstat == SCIP_NLPTERMSTAT_ITLIM);

      /* set a small time limit */
      SCIP_CALL( solveQP(worhpip, 1, 500, -100.0, 100.0, 1.0, INT_MAX, &solval, &solstat, &termstat) );
      cr_expect(termstat == SCIP_NLPTERMSTAT_TILIM);
   }

   if( worhpsqp != NULL )
   {
      /* set a small iteration limit */
      SCIP_CALL( solveQP(worhpsqp, 1, 100, -100.0, 100.0, SCIPinfinity(scip), 5, &solval, &solstat, &termstat) );
      cr_expect(termstat == SCIP_NLPTERMSTAT_ITLIM);

      /* set a small time limit */
      SCIP_CALL( solveQP(worhpsqp, 1, 500, -100.0, 100.0, 1.0, INT_MAX, &solval, &solstat, &termstat) );
      cr_expect(termstat == SCIP_NLPTERMSTAT_TILIM);
   }
}

static
SCIP_RETCODE resolveAfterFixingVars(
   SCIP_NLPI* nlpi
   )
{
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_Real* primal;
   SCIP_Real* vals;
   int* inds;
   SCIP_Real lbs[2] = {0.0, -2.0};
   SCIP_Real ubs[2] = {1.0, +2.0};
   const char* varnames[2] = {"x0", "x1"};

   SCIP_CALL( SCIPallocBufferArray(scip, &inds, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, &nlpiprob, "QP") );
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_FEASTOL, 1e-6) );

   /* add variables */
   SCIP_CALL( SCIPnlpiAddVars(nlpi, nlpiprob, 2, lbs, ubs, varnames) );

   /* add objective */
   inds[0] = 0;
   inds[1] = 1;
   vals[0] = -1.0;
   vals[1] = +1.0;
   SCIP_CALL( SCIPnlpiSetObjective(nlpi, nlpiprob, 2, inds, vals, 0, NULL, NULL, NULL, 0.0) );

   /* first solve */
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   cr_expect(SCIPnlpiGetTermstat(nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPnlpiGetSolstat(nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_expect(SCIPisFeasEQ(scip, primal[0], 1.0));
   cr_expect(SCIPisFeasEQ(scip, primal[1], -2.0));

   /* fix x0 and resolve */
   lbs[0] = 0.0;
   ubs[0] = 0.0;
   SCIP_CALL( SCIPnlpiChgVarBounds(nlpi, nlpiprob, 1, inds, lbs, ubs) );
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );
   cr_expect(SCIPnlpiGetTermstat(nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPnlpiGetSolstat(nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_expect(SCIPisFeasEQ(scip, primal[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, primal[1], -2.0));

   /* free memory */
   SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &nlpiprob) );
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}

Test(nlpi, fixvars, .init = setup, .fini = teardown,
   .description = "resolves a QP after fixing some variables"
   )
{
   int i;

   for( i = 0; i < SCIPgetNNlpis(scip); ++i )
   {
      SCIP_NLPI* nlpi = SCIPgetNlpis(scip)[i];
      assert(nlpi != NULL);

      SCIP_CALL( resolveAfterFixingVars(nlpi) );
   }
}
