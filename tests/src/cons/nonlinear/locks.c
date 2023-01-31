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

/**@file   locks.c
 * @brief  tests locking mechanism for nonlinear constraints
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_CONS* cons;
static SCIP_Bool success;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   cons = NULL;
   success = FALSE;
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/* macro for doing a cr_expect and return FALSE if condition is FALSE
 * note: x is evaluated twice
 */
#define EXPECTANDRETURN(x) do { cr_expect(x); if (!(x)) return FALSE; } while(FALSE)

/* auxiliary function to check locks of variables and their corresponding expressions */
static
SCIP_Bool checkVarLocks(
   int                  nxlockspos,         /**< target number of positive locks */
   int                  nxlocksneg,         /**< target number of negative locks */
   int                  nylockspos,         /**< target number of positive locks */
   int                  nylocksneg          /**< target number of negative locks */
   )
{
#if SCIP_DISABLED_CODE
   /* get expressions for x and y as used in conshdlr: these will have locks updated */
   SCIP_HASHMAP* var2expr = (SCIP_HASHMAP*)SCIPgetVarExprHashmapNonlinear(scip, conshdlr);
   SCIP_EXPR* curxexpr = (SCIP_EXPR*)SCIPhashmapGetImage(var2expr, x);
   SCIP_EXPR* curyexpr = (SCIP_EXPR*)SCIPhashmapGetImage(var2expr, y);

   /* if no locks at all, then there may not even be an expression for the var */
   EXPECTANDRETURN(nxlockspos != 0 || nxlocksneg != 0 || curxexpr == NULL || (SCIPgetExprNLocksPosNonlinear(curxexpr) == 0 && SCIPgetExprNLocksNegNonlinear(curxexpr) == 0));
   EXPECTANDRETURN((nxlockspos == 0 && nxlocksneg == 0) || SCIPgetExprNLocksPosNonlinear(curxexpr) == nxlockspos);
   EXPECTANDRETURN((nxlockspos == 0 && nxlocksneg == 0) || SCIPgetExprNLocksNegNonlinear(curxexpr) == nxlocksneg);

   EXPECTANDRETURN(nylockspos != 0 || nylocksneg != 0 || curyexpr == NULL || (SCIPgetExprNLocksPosNonlinear(curyexpr) == 0 && SCIPgetExprNLocksNegNonlinear(curyexpr) == 0));
   EXPECTANDRETURN((nylockspos == 0 && nylocksneg == 0) || SCIPgetExprNLocksPosNonlinear(curyexpr) == nylockspos);
   EXPECTANDRETURN((nylockspos == 0 && nylocksneg == 0) || SCIPgetExprNLocksNegNonlinear(curyexpr) == nylocksneg);
#endif
   EXPECTANDRETURN(SCIPvarGetNLocksDown(x) == nxlocksneg);
   EXPECTANDRETURN(SCIPvarGetNLocksUp(x) == nxlockspos);
   EXPECTANDRETURN(SCIPvarGetNLocksDown(y) == nylocksneg);
   EXPECTANDRETURN(SCIPvarGetNLocksUp(y) == nylockspos);

   return TRUE;
}

/* auxiliary function to set bounds of variables */
static
SCIP_RETCODE chgBounds(
   SCIP_Real lbx,          /**< new lower bound of x */
   SCIP_Real ubx,          /**< new upper bound of x */
   SCIP_Real lby,          /**< new lower bound of y */
   SCIP_Real uby           /**< new upper bound of y */
   )
{
   SCIP_CALL( SCIPchgVarLb(scip, x, lbx) );
   SCIP_CALL( SCIPchgVarUb(scip, x, ubx) );
   SCIP_CALL( SCIPchgVarLb(scip, y, lby) );
   SCIP_CALL( SCIPchgVarUb(scip, y, uby) );

   SCIPincrementCurBoundsTagNonlinear(SCIPfindConshdlr(scip, "nonlinear"), TRUE);

   return SCIP_OKAY;
}

/* define the test suite */
TestSuite(locks, .init = setup, .fini = teardown);

/*
 * tests
 */

Test(locks, sum)
{
   const char* input = "[nonlinear] <test> : <x> - <y> <= 1.0";

   SCIP_CALL( SCIPparseCons(scip, &cons, input, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success);
   cr_assert_not_null(cons);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );
   cr_expect( checkVarLocks(1, 0, 0, 1) );

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   cr_expect( checkVarLocks(0, 0, 0, 0) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* test for changing bounds between locking of a single constraint */
Test(locks, chg_bounds)
{
   const char* input = "[nonlinear] <test> : (<x>)^2 <= 1.0";

   SCIP_CALL( SCIPparseCons(scip, &cons, input, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success);
   cr_assert_not_null(cons);

   /*
    * x^2 is not monotone for [-1,1]
    */
   SCIP_CALL( chgBounds(-1.0, 1.0, 0.0, 0.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );
   cr_expect( checkVarLocks(1, 1, 0, 0) );

   /*
    * x^2 is monotone decreasing for [-1,0], but locking should still use old monotonicity
    */
   SCIP_CALL( chgBounds(-1.0, 0.0, 0.0, 0.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   cr_expect( checkVarLocks(0, 0, 0, 0) );

   /* locking x^2 again should result in only a single up-lock */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );
   cr_expect( checkVarLocks(0, 1, 0, 0) );

   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   cr_expect( checkVarLocks(0, 0, 0, 0) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* test locks for an non-monotone expression */
Test(locks, non_monotone)
{
   const char* input = "[nonlinear] <test> : sin(<x>^2) - cos(<y>^0.5) >= -.10";

   SCIP_CALL( SCIPparseCons(scip, &cons, input, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success);
   cr_assert_not_null(cons);

   /* add locks */
   SCIP_CALL( chgBounds(-10.0, 10.0, 1.0, 10.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );
   cr_expect( checkVarLocks(2, 2, 1, 1) );

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   cr_expect( checkVarLocks(0, 0, 0, 0) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* tests locks for a complex expression */
Test(locks, complex)
{
   const char* input = "[nonlinear] <test> : exp(<x>^2 + <x>*<y> - log(abs(<y>^3))) <= 0.0";

   SCIP_CALL( SCIPparseCons(scip, &cons, input, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success);
   cr_assert_not_null(cons);

   /* add locks */
   SCIP_CALL( chgBounds(-1.0, 1.0, -3.0, -2.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );
   cr_expect( checkVarLocks(1, 2, 2, 1) );

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   cr_expect( checkVarLocks(0, 0, 0, 0) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
