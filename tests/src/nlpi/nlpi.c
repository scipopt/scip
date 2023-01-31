/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlpi.c
 * @brief  unit test for nonlinear interface methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "scip/nlpioracle.h"
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
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* need a problem to have stat created, which is used by expr iterators */
   SCIP_CALL( SCIPcreateProbBasic(scip, "dummy") );

   ipopt = SCIPfindNlpi(scip, "ipopt");
   worhpip = SCIPfindNlpi(scip, "worhp-ip");
   worhpsqp = SCIPfindNlpi(scip, "worhp-sqp");
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
   SCIP_NLPSTATISTICS statistics;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_Real* dualcons;
   SCIP_Real* duallb;
   SCIP_Real* dualub;
   SCIP_Real* primal;

   /* variables */
   const char* varnames[4] = {"x0", "x1", "x2", "x3"};
   SCIP_EXPR* varexprs[4];
   SCIP_Real lbs[4] = {-0.5, -2, 0, -2};
   SCIP_Real ubs[4] = {INF, INF, 2, 2};
   SCIP_Real start[4] = {0.1,0.6,1.1,2.9};

   /* objective */
   SCIP_Real objvals[1] = {-1};
   int objinds[1] = {2};

   /* constraints */
   const char* consnames[3] = {"c1", "c2", "c3"};
   SCIP_Real lhss[3] = {1, -INF, 2.5};
   SCIP_Real rhss[3] = {1, -1, 5};
   SCIP_Real initguess[4];
   SCIP_Real* linvals;
   int* lininds;
   int nlinds;

   SCIP_EXPR* x0sqr;
   SCIP_EXPR* x1sqr;
   SCIP_EXPR* x2sqr;
   SCIP_EXPR* x0x2;
   SCIP_EXPR* expr;

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
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, 2) );
   SCIP_CALL( SCIPcreateNlpiProblem(scip, nlpi, &nlpiprob, "convex_NLP") );
   SCIP_CALL( SCIPaddNlpiVars(scip, nlpi, nlpiprob, 4, lbs, ubs, varnames) );

   /* create expressions */
   SCIP_CALL( SCIPcreateExprVaridx(scip, &varexprs[0], 0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVaridx(scip, &varexprs[1], 1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVaridx(scip, &varexprs[2], 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVaridx(scip, &varexprs[3], 3, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &x0sqr, varexprs[0], 2.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &x1sqr, varexprs[1], 2.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &x2sqr, varexprs[2], 2.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &x0x2, 1, &varexprs[0], 1.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, x0x2, varexprs[2]) );

   /* set objective */
   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 1, &x0sqr, NULL, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr, x1sqr, 2.0) );
   SCIP_CALL( SCIPsetNlpiObjective(scip, nlpi, nlpiprob, 1, objinds, objvals, expr, 0.0) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* constraints */
   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 1, &x0sqr, NULL, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr, x2sqr, 1.0) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr, x0x2, 1.0) );
   SCIP_CALL( SCIPaddNlpiConstraints(scip, nlpi, nlpiprob, 1, &lhss[0], &rhss[0], NULL, NULL, NULL, &expr, &consnames[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   nlinds = 2;
   lininds[0] = 2;
   linvals[0] = 1;
   lininds[1] = 3;
   linvals[1] = -1;
   SCIP_CALL( SCIPaddNlpiConstraints(scip, nlpi, nlpiprob, 1, &lhss[1], &rhss[1], &nlinds, &lininds, &linvals, NULL, &consnames[1]) );

   nlinds = 2;
   lininds[0] = 1;
   linvals[0] = 1;
   lininds[1] = 3;
   linvals[1] = 1;
   SCIP_CALL( SCIPaddNlpiConstraints(scip, nlpi, nlpiprob, 1, &lhss[2], &rhss[2], &nlinds, &lininds, &linvals, NULL, &consnames[2]) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetNlpiIntPar(scip, nlpi, nlpiprob, SCIP_NLPPAR_VERBLEVEL, 1) );
#endif

   /* set a starting point close to solution to improve likelihood to converge to expected solution on this nonconvex NLP */
   SCIP_CALL( SCIPsetNlpiInitialGuess(scip, nlpi, nlpiprob, start, NULL, NULL, NULL) );

   /* solve NLP */
   SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiprob, .feastol = 1e-9) );
   cr_expect(SCIPgetNlpiTermstat(scip, nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPgetNlpiSolstat(scip, nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);

   /* collect statistics */
   SCIP_CALL( SCIPgetNlpiStatistics(scip, nlpi, nlpiprob, &statistics) );
   cr_expect(statistics.niterations > 0);
   cr_expect(statistics.totaltime >= 0.0);

   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, &dualcons, &duallb, &dualub, NULL) );

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
   SCIP_CALL( SCIPchgNlpiVarBounds(scip, nlpi, nlpiprob, 1, lininds, lbs, ubs) );

   /* set the initial guess to the previous solution */
   initguess[0] = 0.6;
   initguess[1] = 0.5;
   initguess[2] = 0.4;
   initguess[3] = 2.0;
   SCIP_CALL( SCIPsetNlpiInitialGuess(scip, nlpi, nlpiprob, initguess, NULL, NULL, NULL) );

   /* solve NLP */
   SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiprob, .feastol = 1e-9) );
   cr_expect(SCIPgetNlpiTermstat(scip, nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPgetNlpiSolstat(scip, nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);

   /* collect statistics */
   SCIP_CALL( SCIPgetNlpiStatistics(scip, nlpi, nlpiprob, &statistics) );
   cr_expect(statistics.niterations > 0);
   cr_expect(statistics.totaltime >= 0.0);

   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, &dualcons, &duallb, &dualub, NULL) );

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
   SCIP_CALL( SCIPfreeNlpiProblem(scip, nlpi, &nlpiprob) );

   SCIP_CALL( SCIPreleaseExpr(scip, &x0x2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &x2sqr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &x1sqr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &x0sqr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[3]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[2]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[0]) );

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
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_NLPPARAM nlpparam = SCIP_NLPPARAM_DEFAULT(scip);
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_EXPR** varexprs;
   SCIP_EXPR* sumexpr;
   SCIP_EXPR* prodexpr;
   SCIP_Real* primal;
   SCIP_Real* vals;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real objval;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int* inds;
   int nlininds;
   int objind;
   int i;

   *solval = SCIP_INVALID;
   *nlpsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
   *nlptermstat = SCIP_NLPTERMSTAT_OTHER;

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, rndseed, TRUE) );

   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, n+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, n+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, n) );

   SCIP_CALL( SCIPcreateNlpiProblem(scip, nlpi, &nlpiprob, "QP") );

   /* create variables */
   for( i = 0; i < n; ++i )
   {
      lbs[i] = SCIPrandomGetReal(randnumgen, minlb, maxub);
      ubs[i] = SCIPrandomGetReal(randnumgen, lbs[i], maxub);
      SCIP_CALL( SCIPcreateExprVaridx(scip, &varexprs[i], i, NULL, NULL) );
   }
   lbs[n] = -SCIPinfinity(scip);
   ubs[n] = SCIPinfinity(scip);
   SCIP_CALL( SCIPaddNlpiVars(scip, nlpi, nlpiprob, n+1, lbs, ubs, NULL) );

   /* set objective */
   objind = n;
   objval = 1.0;
   SCIP_CALL( SCIPsetNlpiObjective(scip, nlpi, nlpiprob, 1, &objind, &objval, NULL, 0.0) );

   /* create constraint */
   for( i = 0; i < n; ++i )
      vals[i] = SCIPrandomGetReal(randnumgen, -10.0, 10.0);

   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 0, NULL, NULL, 0.0, NULL, NULL) );
   for( i = 0; i < n; ++i )
   {
      int j;
      for( j = 0; j <= i; ++j )
      {
         SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr, 2, (SCIP_EXPR*[2]){ varexprs[j], varexprs[i] }, 1.0, NULL, NULL) );
         SCIP_CALL( SCIPappendExprSumExpr(scip, sumexpr, prodexpr, (j < i) ? 2.0 * vals[i] * vals[j] : vals[i] * vals[j]) );
         SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );
      }
   }

   nlininds = 1;
   inds[0] = n;
   vals[0] = -1.0;
   lhs = 0.0;
   rhs = 0.0;
   SCIP_CALL( SCIPaddNlpiConstraints(scip, nlpi, nlpiprob, 1, &lhs, &rhs, &nlininds, &inds, &vals, &sumexpr, NULL) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );

   /* set parameters */
   nlpparam.timelimit = timelim;
   nlpparam.iterlimit = iterlim;
#ifdef SCIP_DEBUG
   nlpparam.verblevel = 1;
#endif

   /* solve NLP */
   SCIP_CALL( SCIPsolveNlpiParam(scip, nlpi, nlpiprob, nlpparam) );
   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );

   *solval = primal[n];
   *nlpsolstat = SCIPgetNlpiSolstat(scip, nlpi, nlpiprob);
   *nlptermstat = SCIPgetNlpiTermstat(scip, nlpi, nlpiprob);

   /* free memory */
   SCIP_CALL( SCIPfreeNlpiProblem(scip, nlpi, &nlpiprob) );

   for( i = n-1; i >= 0; --i )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[i]) );
   }
   SCIPfreeBufferArray(scip, &varexprs);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &inds);
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

#define DIM 50

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

   for( i = 0; i < 3; ++i )
   {
      /* solve QP with Ipopt */
      if( ipopt != NULL )
      {
         SCIP_CALL( solveQP(ipopt, i+1, DIM, -100.0, 100.0, SCIPinfinity(scip), INT_MAX, &ipoptval,
               &ipoptsolstat, &ipopttermstat) );
      }

      /* solve QP with WORHP-IP */
      if( worhpip != NULL )
      {
         SCIP_CALL( solveQP(worhpip, i+1, DIM, -100.0, 100.0, SCIPinfinity(scip), INT_MAX, &worhpipval,
               &worhpipsolstat, &worhpiptermstat) );
      }

      /* solve QP with WORHP-SQP */
      if( worhpsqp != NULL )
      {
         SCIP_CALL( solveQP(worhpsqp, i+1, DIM, -100.0, 100.0, SCIPinfinity(scip), INT_MAX, &worhpsqpval,
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
      cr_expect(termstat == SCIP_NLPTERMSTAT_ITERLIMIT);

      /* set a small time limit */
      SCIP_CALL( solveQP(worhpip, 1, 500, -100.0, 100.0, 1.0, INT_MAX, &solval, &solstat, &termstat) );
      cr_expect(termstat == SCIP_NLPTERMSTAT_TIMELIMIT);
   }

   if( worhpsqp != NULL )
   {
      /* set a small iteration limit */
      SCIP_CALL( solveQP(worhpsqp, 1, 100, -100.0, 100.0, SCIPinfinity(scip), 5, &solval, &solstat, &termstat) );
      cr_expect(termstat == SCIP_NLPTERMSTAT_ITERLIMIT);

      /* set a small time limit */
      SCIP_CALL( solveQP(worhpsqp, 1, 500, -100.0, 100.0, 1.0, INT_MAX, &solval, &solstat, &termstat) );
      cr_expect(termstat == SCIP_NLPTERMSTAT_TIMELIMIT);
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

   SCIP_CALL( SCIPcreateNlpiProblem(scip, nlpi, &nlpiprob, "LP") );

   /* add variables */
   SCIP_CALL( SCIPaddNlpiVars(scip, nlpi, nlpiprob, 2, lbs, ubs, varnames) );

   /* add objective */
   inds[0] = 0;
   inds[1] = 1;
   vals[0] = -1.0;
   vals[1] = +1.0;
   SCIP_CALL( SCIPsetNlpiObjective(scip, nlpi, nlpiprob, 2, inds, vals, NULL, 0.0) );

   /* first solve */
   SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiprob) );
   cr_expect(SCIPgetNlpiTermstat(scip, nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPgetNlpiSolstat(scip, nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);
   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_expect(SCIPisFeasEQ(scip, primal[0], 1.0));
   cr_expect(SCIPisFeasEQ(scip, primal[1], -2.0));

   /* fix x0 and resolve */
   lbs[0] = 0.0;
   ubs[0] = 0.0;
   SCIP_CALL( SCIPchgNlpiVarBounds(scip, nlpi, nlpiprob, 1, inds, lbs, ubs) );
   SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiprob) );
   cr_expect(SCIPgetNlpiTermstat(scip, nlpi, nlpiprob) == SCIP_NLPTERMSTAT_OKAY);
   cr_expect(SCIPgetNlpiSolstat(scip, nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT);
   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, NULL, NULL, NULL, NULL) );
   cr_expect(SCIPisFeasEQ(scip, primal[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, primal[1], -2.0));

   /* free memory */
   SCIP_CALL( SCIPfreeNlpiProblem(scip, nlpi, &nlpiprob) );
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
