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

/**@file   unittest-lpibases
 * @brief  unit test for checking the settings of slack variables in a basis of the lpi
 * @author Marc Pfetsch
 *
 * The behavior of different LP solvers w.r.t. the slack variables should not differ, if interfaced by LPI.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <lpi/lpi.h>


/** run unittest */
static
SCIP_RETCODE runUnittest(void)
{
   SCIP_LPI* lpi = NULL;
   SCIP_Real obj = 1.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 3.0;
   SCIP_Real vals[2];
   SCIP_Real lhs = 1.0;
   SCIP_Real rhs = 2.0;
   SCIP_Real val = 1.0;
   SCIP_Real objval;
   SCIP_Real binvrow[3];
   SCIP_Real coef[3];
   int inds[2];
   int cstats[3];
   int rstats[3];
   int basinds[3];
   int ind = 0;
   int beg = 0;
   int nrows;
   int ncols;
   int cstat;
   int rstat;
   int i;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-boundchg ===========\n");
   printf("=opt=  unittest-lpibases 0\n\n");

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   /* --- Test 1 ---------------------------------------------------- */

   /* use the following LP:
    *   max x
    *       1 <= x <= 2  (linear constraint)
    *       0 <= x <= 3  (bounds)
    */

   /* add one column */
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* add one row */
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val) );

   /* check size */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   assert( nrows == 1 );
   assert( ncols == 1 );

   /* debugging output */
   /* SCIP_CALL( SCIPlpiWriteLP(lpi, "debug.lp") ); */

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* get basis */
   SCIP_CALL( SCIPlpiGetBase(lpi, &cstat, &rstat) );

   /* the variable should be basic and the slack variable at the upper bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_UPPER );

   printf("Test 1 passed.\n");

   /* --- Test 2 ---------------------------------------------------- */

   /* use the following LP:
    *   min x
    *       1 <= x <= 2  (linear constraint)
    *       0 <= x <= 3  (bounds)
    */

   /* change sense */
   SCIP_CALL( SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MINIMIZE) );

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* get basis */
   SCIP_CALL( SCIPlpiGetBase(lpi, &cstat, &rstat) );

   /* the variable should be basic and the slack variable at the lower bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_LOWER );

   printf("Test 2 passed.\n");

   /* --- Test 3 ---------------------------------------------------- */

   /* use the following LP:
    *   min x
    *       1 <= x       (linear constraint)
    *       0 <= x <= 3  (bounds)
    */

   rhs = SCIPlpiInfinity(lpi);
   SCIP_CALL( SCIPlpiChgSides(lpi, 1, &ind, &lhs, &rhs) );

   /* debugging output */
   /* SCIP_CALL( SCIPlpiWriteLP(lpi, "debug.lp") ); */

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* get basis */
   SCIP_CALL( SCIPlpiGetBase(lpi, &cstat, &rstat) );

   /* the variable should be basic and the slack variable at the lower bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_LOWER );

   printf("Test 3 passed.\n");

   /* --- Test 4 ---------------------------------------------------- */

   /* use the following LP:
    *   max x
    *       x \leq 1         (linear constraint)
    *       0 \leq x \leq 3  (bounds)
    */

   /* change sense */
   SCIP_CALL( SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MAXIMIZE) );

   lhs = -SCIPlpiInfinity(lpi);
   rhs = 1.0;
   SCIP_CALL( SCIPlpiChgSides(lpi, 1, &ind, &lhs, &rhs) );

   /* debugging output */
   /* SCIP_CALL( SCIPlpiWriteLP(lpi, "debug.lp") ); */

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* get basis */
   SCIP_CALL( SCIPlpiGetBase(lpi, &cstat, &rstat) );

   /* the variable should be basic and the slack variable at the upper bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_UPPER );

   printf("Test 4 passed.\n");

   SCIP_CALL( SCIPlpiFree(&lpi) );

   /* --- Test 5 ---------------------------------------------------- */

   /* use the following LP:
    * max 1 x1 + 1 x2 + 1 x3
    *       -8 <= -x1 -          x3 <= -1
    *       -7 <= -x1 -   x2        <= -1
    *              x1 + 2 x2        <= 12
    *              x1,    x2,    x3 >= 0
    */

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   /* add columns */
   lb = 0.0;
   ub = SCIPlpiInfinity(lpi);
   obj = 1.0;

   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* add rows */
   lhs = -8.0;
   rhs = -1.0;
   inds[0] = 0;
   inds[1] = 2;
   vals[0] = -1.0;
   vals[1] = -1.0;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, inds, vals) );

   lhs = -7.0;
   rhs = -1.0;
   inds[0] = 0;
   inds[1] = 1;
   vals[0] = -1.0;
   vals[1] = -1.0;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, inds, vals) );

   lhs = -SCIPlpiInfinity(lpi);
   rhs = 12.0;
   inds[0] = 0;
   inds[1] = 1;
   vals[0] = 1.0;
   vals[1] = 2.0;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, inds, vals) );

   /* check size */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   assert( nrows == 3 );
   assert( ncols == 3 );

   /* debugging output */
   /* SCIP_CALL( SCIPlpiWriteLP(lpi, "debug.lp") ); */

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
   assert( fabs(objval - 14) < 1e-6 );

   /* the optimal basis should be: {x2, x3, slack for second row} */
   SCIP_CALL( SCIPlpiGetBase(lpi, cstats, rstats) );
   assert( cstats[0] == SCIP_BASESTAT_LOWER );
   assert( cstats[1] == SCIP_BASESTAT_BASIC );
   assert( cstats[2] == SCIP_BASESTAT_BASIC );

   assert( rstats[0] == SCIP_BASESTAT_LOWER );
   assert( rstats[1] == SCIP_BASESTAT_BASIC );
   assert( rstats[2] == SCIP_BASESTAT_UPPER );

   /* get basis indices */
   SCIP_CALL( SCIPlpiGetBasisInd(lpi, basinds) );

   /* search for slack variable in basis */
   for (i = 0; i < nrows; ++i)
   {
      if ( basinds[i] < 0 )
         break;
   }
   assert( i < nrows );

   /* check basis inverse */
   SCIP_CALL( SCIPlpiGetBInvRow(lpi, i, binvrow, NULL, NULL) );

   /* row of basis inverse should be (0, 1, 0.5) */
   assert( fabs(binvrow[0]) < 10e-6 );
   assert( fabs(binvrow[1] - 1.0) < 10e-6 );
   assert( fabs(binvrow[2] - 0.5) < 10e-6 );

   /* check basis inverse times nonbasic matrix */
   SCIP_CALL( SCIPlpiGetBInvARow(lpi, i, binvrow, coef, NULL, NULL) );

   /* row of basis inverse times nonbasic matrix should be (0.5, 0, 0) */
   assert( fabs(coef[0] + 0.5) < 10e-6 );
   assert( fabs(coef[1]) < 10e-6 );
   assert( fabs(coef[2]) < 10e-6 );

   printf("Test 5 passed.\n");

   /* --------------------------------------------------------------- */

   SCIP_CALL( SCIPlpiFree(&lpi) );

   /* for automatic testing output the following */
   printf("\n");
   printf("SCIP Status        : test passed\n");
   printf("Ignore the following:\n");
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
