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

void initProb(int); //TODO without this we have a prototype missing warning...

void initProb(int option) 
{
   // TODO add a minimization problem
   if(0 == option) {
      /* unbounded - infeasible
       * 
       * (P):  max x
       * -x <= 1 (constr)
       *
       *  0 <= x (bound)
       *
       * (D):  min y
       * 1 <= -y (constr)
       *
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
   } else if(1 == option) {
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
      SCIP_Real val[2] = {1.0, 1.0};
      int beg[2] = {0, 1};
      int ind[2] = {0, 1};

      /* add columns */
      SCIP_CALL( SCIPlpiAddCols(lpi, 2, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );

      /* add rows */
      SCIP_CALL( SCIPlpiAddRows(lpi, 2, lhs, rhs, NULL, 1, beg, ind, val) );
   } else if(2 == option) {
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

TheoryDataPoints(change, testchgbound) = 
{
   //TODO add the infinities...
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   DataPoints(int, 0, 1, 2)
};

Theory((SCIP_Real upper, SCIP_Real lower, int prob), change, testchgbound)
{
   cr_assume_lt(lower, upper);
   initProb(prob);

   SCIP_Real setub[1] = { upper };
   SCIP_Real setlb[1] = { lower };
   int ind[1] = {0};
   SCIP_CALL( SCIPlpiChgBounds(lpi, 1, ind, setlb, setub) );

   SCIP_Real ub[1];
   SCIP_Real lb[1];
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, 0, lb, ub) );
   
   cr_assert_arr_eq(ub, setub, 1);
   cr_assert_arr_eq(lb, setlb, 1);
}

TheoryDataPoints(change, testchgbounds) = 
{
   //TODO add the infinities...
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   // Problem 0 has only one variable 
   DataPoints(int, 1, 2)
};

Theory((SCIP_Real upper1, SCIP_Real upper2, SCIP_Real lower1, SCIP_Real lower2, int prob), change, testchgbounds)
{
   cr_assume_lt(lower1, upper1);
   cr_assume_lt(lower2, upper2);
   initProb(prob);

   SCIP_Real setub[2] = { upper1, upper2 };
   SCIP_Real setlb[2] = { lower1, lower2 };
   int ind[2] = {0, 1};
   SCIP_CALL( SCIPlpiChgBounds(lpi, 2, ind, setlb, setub) );

   SCIP_Real ub[2];
   SCIP_Real lb[2];
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, 1, lb, ub) );
   
   cr_assert_arr_eq(ub, setub, 2);
   cr_assert_arr_eq(lb, setlb, 2);
}

TheoryDataPoints(change, testchgside) = 
{
   //TODO add the infinities...
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   DataPoints(int, 0, 1, 2)
};

Theory((SCIP_Real left, SCIP_Real right, int prob), change, testchgside)
{
   cr_assume_lt(left, right);
   initProb(prob);

   SCIP_Real setls[1] = { left };
   SCIP_Real setrs[1] = { right };
   int ind[1] = {0};
   SCIP_CALL( SCIPlpiChgBounds(lpi, 1, ind, setls, setrs) );

   SCIP_Real ls[1];
   SCIP_Real rs[1];
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, 0, ls, rs) );
   
   cr_assert_arr_eq(ls, setls, 1);
   cr_assert_arr_eq(rs, setrs, 1);
}

TheoryDataPoints(change, testchgsides) = 
{
   //TODO add the infinities...
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -3, 0, 5034.2, 2e15, 8e-7), 
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   DataPoints(SCIP_Real, -10, 0, -342.3, 2e13, -2e-9),
   // Problem 0 has only one variable 
   DataPoints(int, 1, 2)
};

Theory((SCIP_Real left1, SCIP_Real left2, SCIP_Real right1, SCIP_Real right2, int prob), change, testchgsides)
{
   cr_assume_lt(right1, left1);
   cr_assume_lt(right2, left2);
   initProb(prob);

   SCIP_Real setrs[2] = { left1, left2 };
   SCIP_Real setls[2] = { right1, right2 };
   int ind[2] = {0, 1};
   SCIP_CALL( SCIPlpiChgBounds(lpi, 2, ind, setls, setrs) );

   SCIP_Real rs[2];
   SCIP_Real ls[2];
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, 1, ls, rs) );
   
   cr_assert_arr_eq(rs, setrs, 2);
   cr_assert_arr_eq(ls, setls, 2);
}

TheoryDataPoints(change, testchgobjsen) = 
{
   DataPoints(SCIP_OBJSEN, SCIP_OBJSEN_MAXIMIZE, SCIP_OBJSEN_MINIMIZE),
   DataPoints(int, 0, 1, 2)
};

Theory((SCIP_OBJSEN sense, int prob), change, testchgobjsen) 
{
   initProb(prob);
   SCIP_CALL( SCIPlpiChgObjsen(lpi, sense) );
   
   SCIP_OBJSEN probsense;
   SCIP_CALL( SCIPlpiGetObjsen(lpi, &probsense) );

   cr_assert_eq( sense, probsense );
}
