/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  unit test for checking bound changes
 * @author Marc Pfetsch
 *
 * We perform two tests:
 * - We change the bounds by a very small value and check whether this has an effect on the LP solver interface.
 * - We fix a variable to infinity.
 *
 * In both cases it is unclear what happens. Also LP-solvers might react differently. CPLEX simply ignores small bound changes.
 *
 * For fixing variables to infinity the following variants happen:
 * - SoPlex allows fixing of variables to infinity.
 * - Gurobi returns an error.
 * - CPLEX simply ignores such changes.
 *
 * These tests can be used for debugging or checking the behavior of LP-solvers.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <scip/lpi.h>


/** run unittest */
static
SCIP_RETCODE runUnittest(void)
{
   SCIP_LPI* lpi = NULL;
   SCIP_Real obj = 1.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 1.0;
   SCIP_Real lbnew = 0.0;
   SCIP_Real ubnew = 1.0;
   int ind = 0;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-boundchg ===========\n");
   printf("=opt=  unittest-boundchg 0\n\n");

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   /* add one column */
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* change bound to small value */
   lb = 1e-11;
   ub = 1.0 - 1e-11;
   SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );

   /* get bounds and compare */
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, 0, &lbnew, &ubnew) );

   if ( fabs(lb - lbnew) > 1e-6 )
   {
      SCIPerrorMessage("Violation of lower bounds: %g != %g\n", lb, lbnew);
      return SCIP_ERROR;
   }

   if ( fabs(ub - ubnew) > 1e-6 )
   {
      SCIPerrorMessage("Violation of upper bounds: %g != %g\n", ub, ubnew);
      return SCIP_ERROR;
   }

   /* fix variables to infinity */
   lb = SCIPlpiInfinity(lpi);
   ub = SCIPlpiInfinity(lpi);
   SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );

   /* get bounds and compare */
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, 0, &lbnew, &ubnew) );

   if ( ! SCIPlpiIsInfinity(lpi, ubnew) )
   {
      SCIPerrorMessage("Could not set upper bound to infinity.\n");
      return SCIP_ERROR;
   }

   if ( ! SCIPlpiIsInfinity(lpi, lbnew) )
   {
      SCIPerrorMessage("Could not set lower bound to infinity.\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPlpiFree(&lpi) );

   /* for automatic testing output the following */
   printf("Test passed.\n\n");
   printf("Ignore the following:\n");
   printf("SCIP Status        : problem is solved [optimal solution found]\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return SCIP_OKAY;
}


/** main function */
int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runUnittest();

   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
