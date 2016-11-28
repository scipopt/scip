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

static SCIP* scip;
static SCIP_NLPI* nlpi;

/* creates SCIP, problem, and includes NLP */
static
void setup(void)
{
   /* skip the test if no NLP solver is available */
   if( !SCIPisIpoptAvailableIpopt() && !SCIPisWorhpAvailableWorhp() )
      return;

   SCIP_CALL( SCIPcreate(&scip) );

}

/* frees memory allocated in setup() */
static
void teardown(void)
{
   /* skip the test if no NLP solver is available */
   if( !SCIPisIpoptAvailableIpopt() && !SCIPisWorhpAvailableWorhp() )
      return;

   SCIP_CALL( SCIPfree(&scip) );
}

#define INF 1e+20

/* helper function to test NLPI */
static
SCIP_RETCODE testNlpi()
{
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_Real* primal;

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
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );

   /* check primal solution */
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL) );
   cr_expect(SCIPisFeasEQ(scip, primal[0], 0.0));
   cr_expect(SCIPisFeasEQ(scip, primal[1], 0.5));
   cr_expect(SCIPisFeasEQ(scip, primal[2], 1.0));
   cr_expect(SCIPisFeasEQ(scip, primal[3], 2.0));

   /* free memory */
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &quadelems);
   SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &nlpiprob) );
}

Test(nlpi, worhp, .init = setup, .fini = teardown,
   .description = "checks the NLPI for WORHP"
   )
{
   if( !SCIPisWorhpAvailableWorhp() )
      return;

   SCIP_CALL( SCIPcreateNlpSolverWorhp(SCIPblkmem(scip), &nlpi) );
   cr_assert(nlpi != NULL);

   SCIP_CALL( SCIPincludeNlpi(scip, nlpi) );
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameWorhp(), SCIPgetSolverDescWorhp()) );

   SCIP_CALL( testNlpi() );
}

Test(nlpi, Ipopt, .init = setup, .fini = teardown,
   .description = "checks the NLPI for Ipopt"
   )
{
   if( !SCIPisIpoptAvailableIpopt() )
      return;

   SCIP_CALL( SCIPcreateNlpSolverIpopt(SCIPblkmem(scip), &nlpi) );
   cr_assert(nlpi != NULL);

   SCIP_CALL( SCIPincludeNlpi(scip, nlpi) );
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameIpopt(), SCIPgetSolverDescIpopt()) );

   SCIP_CALL( testNlpi() );
}
