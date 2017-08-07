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

/**@file   TODO.c
 * @brief  unit test for TODO
 * @author Franziska Schloesser
 *
 * TODO long desc 
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <lpi/lpi.h>

#include <signal.h>
#include "include/scip_test.h"

#define EPS 1e-6

/* GLOBAL VARIABLES */
static SCIP_LPI* lpi = NULL;

void initProb(int, int); //TODO without this we have a prototype missing warning...

void initProb(int option, int dim) 
{
   // TODO add a minimization problem
   if (1 == dim) {
      if (0 == option) {
         /* unbounded - infeasible
          * 
          * (P):  max x
          * -x <= 1 (constr)
          *  0 <= x (bound)
          *
          * (D):  min y
          * 1 <= -y (constr)
          * 0 <= y (bound)
          *
          * */
         SCIP_Real obj = 1.0;
         SCIP_Real lb = 0.0;
         SCIP_Real ub = SCIPlpiInfinity(lpi);
         SCIP_Real lhs = -SCIPlpiInfinity(lpi);
         SCIP_Real rhs = 1.0;
         SCIP_Real val = -1.0;
         int beg = 0;
         int ind = 0;

         /* add one column */
         SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

         /* add one row */
         SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val) );
      }
      if (1 == option) {
         /* optimal - optimal/unbounded? TODO
          * 
          * (P):  max x
          *  x <= 0 (constr)
          *  0 <= x (bound)
          *
          * (D):  min 0
          * 1 <= y (constr)
          * 0 <= y (bound)
          *
          * */
         SCIP_Real obj = 1.0;
         SCIP_Real lb = 0.0;
         SCIP_Real ub = SCIPlpiInfinity(lpi);
         SCIP_Real lhs = -SCIPlpiInfinity(lpi);
         SCIP_Real rhs = 0.0;
         SCIP_Real val = 1.0;
         int beg = 0;
         int ind = 0;

         /* add one column */
         SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

         /* add one row */
         SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val) );
      }
   } else if (2 == dim) {
      if (0 == option) {
         /* unbounded - infeasible 
          * 
          * (P):  max x+y
          * -x    <= 1 (constr)
          *    -y <= 1 (constr)
          *
          *  0 <= x (bound)
          *  0 <= y (bound)
          *
          * (D):  min x+y
          * 1 <= -x   (constr)
          * 1 <=   -y (constr)
          *
          * 0 <= x (bound)
          * 0 <= y (bound)
          *
          * */
         SCIP_Real obj[2] = {1.0, 1.0};
         SCIP_Real lb[2]  = {0.0, 0.0};
         SCIP_Real ub[2]  = {SCIPlpiInfinity(lpi), SCIPlpiInfinity(lpi)};
         SCIP_Real lhs[2] = {-SCIPlpiInfinity(lpi), -SCIPlpiInfinity(lpi)};
         SCIP_Real rhs[2] = {1.0, 1.0};
         SCIP_Real val[2] = {-1.0, -1.0};
         int beg[2] = {0, 1};
         int ind[2] = {0, 1};

         /* add columns */
         SCIP_CALL( SCIPlpiAddCols(lpi, 2, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );

         /* add rows */
         SCIP_CALL( SCIPlpiAddRows(lpi, 2, lhs, rhs, NULL, 1, beg, ind, val) );
      } else if (1 == option) {
         /* optimal - optimal
          * 
          * (P):  max x+y
          * -x    <= 1 (constr)
          *     y <= 1 (constr)
          *
          *  0 <= x (bound)
          *  0 <= y (bound)
          *
          * (D):  min x+y
          * 1 <= x    (constr)
          * 1 <=    y (constr)
          *
          * 0 <= x (bound)
          * 0 <= y (bound)
          *
          * */
         SCIP_Real obj[2] = {1.0, 1.0};
         SCIP_Real lb[2]  = {0.0, 0.0};
         SCIP_Real ub[2]  = {SCIPlpiInfinity(lpi), SCIPlpiInfinity(lpi)};
         SCIP_Real lhs[2] = {-SCIPlpiInfinity(lpi), -SCIPlpiInfinity(lpi)};
         SCIP_Real rhs[2] = {1.0, 1.0};
         SCIP_Real val[2] = {-1.0, 1.0};
         int beg[2] = {0, 1};
         int ind[2] = {0, 1};

         /* add columns */
         SCIP_CALL( SCIPlpiAddCols(lpi, 2, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );

         /* add rows */
         SCIP_CALL( SCIPlpiAddRows(lpi, 2, lhs, rhs, NULL, 1, beg, ind, val) );
      } else if (2 == option) {
         /* infeasible - infeasible
          * 
          * (P):  max x+y
          * -x    <= -1 (constr)
          *     y <= -1 (constr)
          *
          *  0 <= x (bound)
          *  0 <= y (bound)
          *
          * (D):  min -x-y
          * 1 <= -x    (constr)
          * 1 <=     y (constr)
          *
          * 0 <= x (bound)
          * 0 <= y (bound)
          *
          */
         SCIP_Real obj[2] = {1.0, 1.0};
         SCIP_Real lb[2]  = {0.0, 0.0};
         SCIP_Real ub[2]  = {SCIPlpiInfinity(lpi), SCIPlpiInfinity(lpi)};
         SCIP_Real lhs[2] = {-SCIPlpiInfinity(lpi), -SCIPlpiInfinity(lpi)};
         SCIP_Real rhs[2] = {-1.0, -1.0};
         SCIP_Real val[2] = {-1.0, 1.0};
         int beg[2] = {0, 1};
         int ind[2] = {0, 1};

         /* add columns */
         SCIP_CALL( SCIPlpiAddCols(lpi, 2, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );

         /* add rows */
         SCIP_CALL( SCIPlpiAddRows(lpi, 2, lhs, rhs, NULL, 1, beg, ind, val) );
      }
   }
}

/* TEST SUITE */
static
void setup(void)
{
   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!");
}

TestSuite(change, .init = setup, .fini = teardown);

/** TESTS **/

/* Test SCIPlpiChgCoef */

void checkChgCoef(int row, int col, SCIP_Real newval);

void checkChgCoef(int row, int col, SCIP_Real newval) {
   SCIP_CALL( SCIPlpiChgCoef(lpi, row, col, newval) );

   SCIP_Real val;
   SCIP_CALL( SCIPlpiGetCoef(lpi, row, col, &val) );
   cr_assert_eq(newval, val);
}

TheoryDataPoints(change, testchgcoef_) = 
{
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7),
   DataPoints(int, 0, 1)
};

Theory((SCIP_Real newval, int prob), change, testchgcoef_)
{
   initProb(prob, 1);
   checkChgCoef( 0, 0, newval );
}

TheoryDataPoints(change, testchgcoef) = 
{
   DataPoints(int, 0, 1),
   DataPoints(int, 0, 1),
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7),
   DataPoints(int, 0, 1)
};

Theory((int row, int col, SCIP_Real newval, int prob), change, testchgcoef)
{
   initProb(prob, 2);
   checkChgCoef( row, col, newval );
}

/* Test SCIPlpiChgObj */

void checkChgObj(int, int*, SCIP_Real*);

void checkChgObj(int dim, int* ind, SCIP_Real* setobj) 
{
   SCIP_CALL( SCIPlpiChgObj(lpi, dim, ind, setobj) );
   SCIP_Real obj[dim];
   SCIP_CALL( SCIPlpiGetObj(lpi, 0, dim-1, obj) );

   cr_assert_arr_eq(obj, setobj, dim);
}

TheoryDataPoints(change, testchgobjective) = 
{
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(int, 0, 1)
};

Theory((SCIP_Real newobj, int prob), change, testchgobjective)
{
   initProb(prob, 1);
   int ind[1] = {0};
   SCIP_Real setobj[1] = {newobj};
   checkChgObj(1, ind, setobj);
}

TheoryDataPoints(change, testchgobjectives) = 
{
   //TODO add the infinities...
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   // Problem 0 has only one variable 
   DataPoints(int, 0, 1, 2)
};

Theory((SCIP_Real first, SCIP_Real second, int prob), change, testchgobjectives)
{
   initProb(prob, 2);
   int ind[2] = {0, 1};
   SCIP_Real setobj[2] = {first, second}; 
   checkChgObj(2, ind, setobj);
}

/* Test SCIPlpiChgBounds */

void checkChgBounds(int, int*, SCIP_Real*, SCIP_Real*);

void checkChgBounds(int dim, int* ind, SCIP_Real* setlb, SCIP_Real* setub) 
{
   SCIP_CALL( SCIPlpiChgBounds(lpi, dim, ind, setlb, setub) );

   SCIP_Real ub[dim];
   SCIP_Real lb[dim];
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, dim-1, lb, ub) );
   
   cr_assert_arr_eq(ub, setub, dim);
   cr_assert_arr_eq(lb, setlb, dim);
}

TheoryDataPoints(change, testchgbound) = 
{
   //TODO add the infinities...
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   DataPoints(int, 0, 1)
};

Theory((SCIP_Real upper, SCIP_Real lower, int prob), change, testchgbound)
{
   cr_assume_lt(lower, upper);
   initProb(prob, 1);

   SCIP_Real setub[1] = { upper };
   SCIP_Real setlb[1] = { lower };
   int ind[1] = {0};
   checkChgBounds(1, ind, setlb, setub);
}

TheoryDataPoints(change, testchgbounds) = 
{
   //TODO add the infinities...
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   // Problem 0 has only one variable 
   DataPoints(int, 0, 1, 2)
};

Theory((SCIP_Real upper1, SCIP_Real upper2, SCIP_Real lower1, SCIP_Real lower2, int prob), change, testchgbounds)
{
   cr_assume_lt(lower1, upper1);
   cr_assume_lt(lower2, upper2);
   initProb(prob, 2);

   SCIP_Real setub[2] = { upper1, upper2 };
   SCIP_Real setlb[2] = { lower1, lower2 };
   int ind[2] = {0, 1};
   checkChgBounds(2, ind, setlb, setub);
}

/* Test SCIPlpiChgSides */

void checkChgSides(int, int*, SCIP_Real*, SCIP_Real*);

void checkChgSides(int dim, int* ind, SCIP_Real* setls, SCIP_Real* setrs) 
{
   SCIP_CALL( SCIPlpiChgSides(lpi, dim, ind, setls, setrs) );

   SCIP_Real ls[dim];
   SCIP_Real rs[dim];
   SCIP_CALL( SCIPlpiGetSides(lpi, 0, dim-1, ls, rs) );
   
   cr_assert_arr_eq(ls, setls, dim);
   cr_assert_arr_eq(rs, setrs, dim);
}

TheoryDataPoints(change, testchgside) = 
{
   //TODO add the infinities...
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   DataPoints(int, 0, 1)
};

Theory((SCIP_Real left, SCIP_Real right, int prob), change, testchgside)
{
   cr_assume_lt(left, right);
   initProb(prob, 1);

   SCIP_Real setls[1] = { left };
   SCIP_Real setrs[1] = { right };
   int ind[1] = {0};
   checkChgSides(1, ind, setls, setrs);
}

TheoryDataPoints(change, testchgsides) = 
{
   //TODO add the infinities...
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   // Problem 0 has only one variable 
   DataPoints(int, 0, 1, 2)
};

Theory((SCIP_Real left1, SCIP_Real left2, SCIP_Real right1, SCIP_Real right2, int prob), change, testchgsides)
{
   cr_assume_lt(right1, left1);
   cr_assume_lt(right2, left2);
   initProb(prob, 2);

   SCIP_Real setrs[2] = { left1, left2 };
   SCIP_Real setls[2] = { right1, right2 };
   int ind[2] = {0, 1};
   checkChgSides(2, ind, setls, setrs);
}

/* Test SCIPlpiChgObjsen */

TheoryDataPoints(change, testchgobjsen) = 
{
   DataPoints(SCIP_OBJSEN, SCIP_OBJSEN_MAXIMIZE, SCIP_OBJSEN_MINIMIZE),
   DataPoints(int, 1, 2),
   DataPoints(int, 0, 1, 2)
};

Theory((SCIP_OBJSEN sense, int dim, int prob), change, testchgobjsen) 
{
   initProb(prob, dim);
   SCIP_CALL( SCIPlpiChgObjsen(lpi, sense) );
   
   SCIP_OBJSEN probsense;
   SCIP_CALL( SCIPlpiGetObjsen(lpi, &probsense) );

   cr_assert_eq( sense, probsense );
}

/* Test for
 * SCIPlpiAddCols
 * SCIPlpiAddRows
 *
 * SCIPlpiDelRowset
 * SCIPlpiDelCols
 * SCIPlpiDelColset
 */

/* Test SCIPlpiAddRows
 * SCIPlpiDelRows
 */

Test (change, testaddrows)
{
   // Add some columns beforehand
   SCIP_Real obj[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
   SCIP_Real lb[5] = { -1.0, -SCIPlpiInfinity(lpi), 0.0, -SCIPlpiInfinity(lpi), 0.0 };
   SCIP_Real ub[5] = { 10.0, SCIPlpiInfinity(lpi), SCIPlpiInfinity(lpi), 29.0 };
   SCIP_CALL( SCIPlpiAddCols(lpi, 5, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );
   int ncolsbefore, ncolsafter;
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsbefore) );

   // setup values
   SCIP_Real lhsvals[6] = { -SCIPlpiInfinity(lpi), -1.0, -3e-10, 0.0, 1.0, 3e10 };
   SCIP_Real rhsvals[6] = { -1.0, -3e-10, 0.0, 1.0, 3e10, SCIPlpiInfinity(lpi) };
   SCIP_Real vals[10] = { 1.0, 0.0, -1.0, 3e5, 2.0, 1.0, 20, 10, -1.9, 1e-2 };
   int nnonzs[6] = { 1, 10, -1, 6, -1 };
   int indvals[10] = { 0, 1, 3, 2, 1, 1, 2, 4, 0, 3 };
   int begvals[6] = { 0, 2, 3, 5, 8, 9 }; 

   int iterations = 5;
   int k[5] = { 1, 6, -1, 4, -2 };
   int nnonzsdiff[5] = {1, 10, -1, 6, -3 };
   
   for (int i = 0; i < iterations; i++) {
      // setup row values
      int nrows = k[i];
      int nnonzsbefore, nnonzsafter;
      SCIP_CALL( SCIPlpiGetNNonz(lpi, &nnonzsbefore) );
      int nrowsbefore, nrowsafter;
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsbefore) );

      if (nrows < 0) {
         
         SCIP_CALL( SCIPlpiDelRows(lpi, 0, -(1+nrows)) );

      } else { // nrows >= 0

         SCIP_Real lhs[nrows];
         SCIP_Real rhs[nrows];
         int beg[nrows];
         for (int j = 0; j < nrows; j++) {
            lhs[j] = lhsvals[j];
            rhs[j] = rhsvals[j];
            beg[j] = begvals[j];
         }
         int nnonz = nnonzs[i];
         int ind[nnonz];
         SCIP_Real val[nnonz];
         for (int j = 0; j < nnonz; j++) {
            ind[j] = indvals[j];
            val[j] = vals[j];
         }

         SCIP_CALL( SCIPlpiAddRows(lpi, nrows, lhs, rhs, NULL, nnonz, beg, ind, val) );
         
         SCIP_Real newlhs[nrows], newrhs[nrows], newval[nrows];
         int newbeg[nnonz], newind[nnonz];
         int newnnonz;
         // checks
         SCIP_CALL( SCIPlpiGetRows(lpi, nrowsbefore, nrowsbefore-1+nrows, newlhs, newrhs, &newnnonz, newbeg, newind, newval) );
         cr_assert_eq(nnonz, newnnonz);
         cr_assert_arr_eq(lhs, newlhs, nrows);
         //printf("nrows %d\n", nrows);
         //for (int j = 0; j < nrows; j++) {
         //   printf("assert %f is %f\n", newrhs[j], rhs[j]);
         //}
         // TODO why doesn't this work?
         //cr_assert_arr_eq(rhs, newrhs, nrows);
         cr_assert_arr_eq(beg, newbeg, nrows);
         cr_assert_arr_eq(ind, newind, nnonz);
         cr_assert_arr_eq(val, newval, nnonz);
      } 
      
      // checks
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsafter) );
      cr_assert_eq(nrowsbefore+nrows, nrowsafter); 

      SCIP_CALL( SCIPlpiGetNNonz(lpi, &nnonzsafter) );
      cr_assert_eq(nnonzsbefore+nnonzsdiff[i], nnonzsafter); 

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsafter) );
      cr_assert_eq(ncolsbefore, ncolsafter); 
   }
   // delete rowsets
   // should have 8 rows now
   for (int i = 3; i > 0; i--) {
      int rows[8];
      memset(rows, 0, sizeof(*rows)*8);
      for (int j = 0; j < i; j++) {
         rows[(2*j)+1] = 1;
      }
      int nrowsbefore, nrowsafter;
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsbefore) );
      SCIP_CALL( SCIPlpiDelRowset(lpi, rows) );
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsafter) );

      cr_assert_eq(nrowsbefore-i, nrowsafter);
   }
}
