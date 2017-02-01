/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlpi.c
 * @brief  unit test for nonlinear interface methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/nlpi_ipopt.h"
#include "nlpi/nlpi_worhp.h"
#include "nlpi/nlpioracle.h"
#include "nlpi/nlpi.h"

#include "include/scip_test.h"

#define INF 1e+20

/* helper function to test NLPI */
static
SCIP_RETCODE testNlpi(SCIP* scip, SCIP_NLPI* nlpi)
{
   SCIP_NLPSTATISTICS* statistics;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_Real* primal;
   SCIP_Real* dualcons;
   SCIP_Real* duallb;
   SCIP_Real* dualub;

   /* variables */
   const char* varnames[4] = {"x0", "x1", "x2", "x3"};
   SCIP_Real lbs[4] = {-0.5, -2, 0, -2};
   SCIP_Real ubs[4] = {INF, INF, 2, 2};

   /* objective */
   int objinds[1] = {2};
   SCIP_Real objvals[1] = {-1};

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
   int i;

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
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_FEASTOL, 1e-6) );
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );

   cr_expect(SCIPnlpiGetTermstat(nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPnlpiGetSolstat(nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);

   /* print statistics */
   SCIP_CALL( SCIPnlpStatisticsCreate(&statistics) );
   SCIP_CALL( SCIPnlpiGetStatistics(nlpi, nlpiprob, statistics) );
   printf("TIME = %f ITER = %d SOLSTAT = %d\n", SCIPnlpStatisticsGetTotalTime(statistics), SCIPnlpStatisticsGetNIterations(statistics), SCIPnlpiGetSolstat(nlpi, nlpiprob));
   SCIPnlpStatisticsFree(&statistics);

   /* check primal solution */
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, &dualcons, &duallb, &dualub) );


   for( i = 0; i < 3; ++i )
   {
      printf("dualcons[%d] = %g\n", i, dualcons[i]);
   }

   for( i = 0; i < 4; ++i )
   {
      printf("duallb[%d] = %g\n", i, duallb[i]);
      printf("dualub[%d] = %g\n", i, dualub[i]);
   }

   cr_expect(SCIPisFeasEQ(scip, primal[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, primal[1], 0.5));
   cr_expect(SCIPisFeasEQ(scip, primal[2], 1.0));
   cr_expect(SCIPisFeasEQ(scip, primal[3], 2.0));

   /* change bounds of x2 */
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
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );

   cr_expect(SCIPnlpiGetTermstat(nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPnlpiGetSolstat(nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);

   /* print statistics */
   SCIP_CALL( SCIPnlpStatisticsCreate(&statistics) );
   SCIP_CALL( SCIPnlpiGetStatistics(nlpi, nlpiprob, statistics) );
   printf("TIME = %f ITER = %d SOLSTAT = %d\n", SCIPnlpStatisticsGetTotalTime(statistics), SCIPnlpStatisticsGetNIterations(statistics), SCIPnlpiGetSolstat(nlpi, nlpiprob));
   SCIPnlpStatisticsFree(&statistics);

   /* check primal solution */
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL) );

   cr_expect(SCIPisFeasEQ(scip, primal[0], 0.651388));
   cr_expect(SCIPisFeasEQ(scip, primal[1], 0.5));
   cr_expect(SCIPisFeasEQ(scip, primal[2], 0.5));
   cr_expect(SCIPisFeasEQ(scip, primal[3], 2.0));

   /* for( i = 0; i < 4; ++i ) */
   /*    printf("x[%d] = %g\n", i, primal[i]); */

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
   SCIP* scip,
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
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_NLPSTATISTICS* statistics;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_Real* primal;
   int i;
   int k;

   SCIP_QUADELEM* quadelems;
   SCIP_Real* vals;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real objval;
   int objind;
   int* inds;
   int nlininds;
   int nquadelems;

   *solval = SCIP_INVALID;
   *nlpsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
   *nlptermstat = SCIP_NLPTERMSTAT_OTHER;

   SCIP_CALL( SCIPrandomCreate(&randnumgen, SCIPblkmem(scip), rndseed) );

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
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL) );

   /* for( i = 0; i < n+1; ++i ) */
   /*    printf("x[%d] = %.10f\n", i, primal[i]); */

   *solval = primal[n];
   *nlpsolstat = SCIPnlpiGetSolstat(nlpi, nlpiprob);
   *nlptermstat = SCIPnlpiGetTermstat(nlpi, nlpiprob);

   /* print statistics */
   SCIP_CALL( SCIPnlpStatisticsCreate(&statistics) );
   SCIP_CALL( SCIPnlpiGetStatistics(nlpi, nlpiprob, statistics) );
   printf("TIME = %f ITER = %d SOLSTAT = %d\n", SCIPnlpStatisticsGetTotalTime(statistics), SCIPnlpStatisticsGetNIterations(statistics), SCIPnlpiGetSolstat(nlpi, nlpiprob));
   SCIPnlpStatisticsFree(&statistics);

   /* free memory */
   SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &nlpiprob) );
   SCIPfreeBufferArray(scip, &quadelems);
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);
   SCIPrandomFree(&randnumgen);

   return SCIP_OKAY;
}

Test(nlpi, worhp, .description = "checks the NLPI for WORHP"
   )
{
   SCIP* scip;
   SCIP_NLPI* nlpi;

   if( !SCIPisWorhpAvailableWorhp() )
      return;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateNlpSolverWorhp(SCIPblkmem(scip), &nlpi) );
   cr_assert(nlpi != NULL);

   SCIP_CALL( SCIPincludeNlpi(scip, nlpi) );
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameWorhp(), SCIPgetSolverDescWorhp()) );

   SCIP_CALL( testNlpi(scip, nlpi) );
   SCIP_CALL( SCIPfree(&scip) );
}

Test(nlpi, ipopt, .description = "checks the NLPI for Ipopt")
{
   SCIP* scip;
   SCIP_NLPI* nlpi;

   if( !SCIPisIpoptAvailableIpopt() )
      return;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateNlpSolverIpopt(SCIPblkmem(scip), &nlpi) );
   cr_assert(nlpi != NULL);

   SCIP_CALL( SCIPincludeNlpi(scip, nlpi) );
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameIpopt(), SCIPgetSolverDescIpopt()) );

   SCIP_CALL( testNlpi(scip, nlpi) );
   SCIP_CALL( SCIPfree(&scip) );
}

Test(nlpi, solveQP, .description = "solves convex QP with different NLPIs")
{
   SCIP* scipworhp;
   SCIP* scipipopt;
   SCIP_NLPI* worhp;
   SCIP_NLPI* ipopt;
   int i;

   if( !SCIPisIpoptAvailableIpopt() || !SCIPisWorhpAvailableWorhp() )
      return;

   SCIP_CALL( SCIPcreate(&scipworhp) );
   SCIP_CALL( SCIPcreateNlpSolverWorhp(SCIPblkmem(scipworhp), &worhp) );
   cr_assert(worhp != NULL);
   SCIP_CALL( SCIPincludeNlpi(scipworhp, worhp) );
   SCIP_CALL( SCIPincludeExternalCodeInformation(scipworhp, SCIPgetSolverNameWorhp(), SCIPgetSolverDescWorhp()) );

   SCIP_CALL( SCIPcreate(&scipipopt) );
   SCIP_CALL( SCIPcreateNlpSolverIpopt(SCIPblkmem(scipipopt), &ipopt) );
   cr_assert(ipopt != NULL);
   SCIP_CALL( SCIPincludeNlpi(scipipopt, ipopt) );
   SCIP_CALL( SCIPincludeExternalCodeInformation(scipipopt, SCIPgetSolverNameIpopt(), SCIPgetSolverDescIpopt()) );

   for( i = 0; i < 5; ++i )
   {
      SCIP_NLPTERMSTAT ipopttermstat;
      SCIP_NLPTERMSTAT worhptermstat;
      SCIP_NLPSOLSTAT ipoptsolstat;
      SCIP_NLPSOLSTAT worhpsolstat;
      SCIP_Real ipoptval;
      SCIP_Real worhpval;

      SCIP_CALL( solveQP(scipipopt, ipopt, i+1, 100, -100.0, 100.0, SCIPinfinity(scipipopt), INT_MAX, &ipoptval, &ipoptsolstat, &ipopttermstat) );
      SCIP_CALL( solveQP(scipworhp, worhp, i+1, 100, -100.0, 100.0, SCIPinfinity(scipworhp), INT_MAX, &worhpval, &worhpsolstat, &worhptermstat) );

      cr_assert(ipoptsolstat != SCIP_NLPSOLSTAT_UNKNOWN);
      cr_assert(ipopttermstat != SCIP_NLPTERMSTAT_OTHER);
      cr_assert(worhpsolstat != SCIP_NLPSOLSTAT_UNKNOWN);
      cr_assert(worhptermstat != SCIP_NLPTERMSTAT_OTHER);

      if( ipoptsolstat == SCIP_NLPSOLSTAT_LOCOPT && worhpsolstat == SCIP_NLPSOLSTAT_LOCOPT )
      {
         cr_assert(SCIPisFeasLE(scipipopt, worhpval, ipoptval));
      }
   }

   SCIP_CALL( SCIPfree(&scipipopt) );
   SCIP_CALL( SCIPfree(&scipworhp) );
}

Test(nlpi, workinglimits, .description = "solves convex QP with a small iteration limit")
{
   SCIP* scip;
   SCIP_NLPI* worhp;
   SCIP_NLPTERMSTAT termstat;
   SCIP_NLPSOLSTAT solstat;
   SCIP_Real solval;

   if( !SCIPisWorhpAvailableWorhp() )
      return;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateNlpSolverWorhp(SCIPblkmem(scip), &worhp) );
   cr_assert(worhp != NULL);
   SCIP_CALL( SCIPincludeNlpi(scip, worhp) );
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameWorhp(), SCIPgetSolverDescWorhp()) );

   /* set a small iteration limit */
   SCIP_CALL( solveQP(scip, worhp, 1, 100, -100.0, 100.0, SCIPinfinity(scip), 5, &solval, &solstat, &termstat) );
   cr_expect(termstat == SCIP_NLPTERMSTAT_ITLIM);

   /* set a small time limit */
   SCIP_CALL( solveQP(scip, worhp, 1, 500, -100.0, 100.0, 1.0, INT_MAX, &solval, &solstat, &termstat) );
   cr_expect(termstat == SCIP_NLPTERMSTAT_TILIM);

   SCIP_CALL( SCIPfree(&scip) );
}
