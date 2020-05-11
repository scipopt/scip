/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   locks.c
 * @brief  tests locking mechanism for expression constraints
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_CONSEXPR_EXPR* xexpr;
static SCIP_CONSEXPR_EXPR* yexpr;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );

   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/* auxiliary function to check locks of variables and their corresponding expressions */
static
SCIP_Bool checkVarLocks(
   int                  nxlockspos,         /**< target number of positive locks */
   int                  nxlocksneg,         /**< target number of negative locks */
   int                  nylockspos,         /**< target number of positive locks */
   int                  nylocksneg          /**< target number of negative locks */
   )
{
   return SCIPgetConsExprExprNLocksPos(xexpr) == nxlockspos
            && SCIPgetConsExprExprNLocksNeg(xexpr) == nxlocksneg
            && SCIPgetConsExprExprNLocksPos(yexpr) == nylockspos
            && SCIPgetConsExprExprNLocksNeg(yexpr) == nylocksneg
            && SCIPvarGetNLocksDown(x) == nxlocksneg
            && SCIPvarGetNLocksUp(x) == nxlockspos
            && SCIPvarGetNLocksDown(y) == nylocksneg
            && SCIPvarGetNLocksUp(y) == nylockspos;
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

   SCIPincrementConsExprCurBoundsTag(conshdlr, TRUE);

   return SCIP_OKAY;
}

/* define the test suite */
TestSuite(locks, .init = setup, .fini = teardown);

/*
 * tests
 */

Test(locks, sum)
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR* expr;

   const char* input = "<x> - <y>";

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -SCIPinfinity(scip), 1.0) );

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );
   cr_expect( checkVarLocks(1, 0, 0, 1) );

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   cr_expect( checkVarLocks(0, 0, 0, 0) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* test for changing bounds between locking of a single constraint */
Test(locks, chg_bounds)
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR* expr;

   const char* input = "(<x>)^2";

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -SCIPinfinity(scip), 1.0) );

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
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* test locks for an non-monotone expression */
Test(locks, non_monotone)
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR* expr;

   const char* input = "sin(<x>^2) - cos(<y>^0.5)";

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -1.0, SCIPinfinity(scip)) );

   /* add locks */
   SCIP_CALL( chgBounds(-10.0, 10.0, 1.0, 10.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );
   cr_expect( checkVarLocks(2, 2, 1, 1) );

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   cr_expect( checkVarLocks(0, 0, 0, 0) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* tests locks for a complex expression */
Test(locks, complex)
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR* expr;

   const char* input = "exp(<x>^2 + <x>*<y> - log(abs(<y>^3)))";

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, -SCIPinfinity(scip), 0.0) );

   /* add locks */
   SCIP_CALL( chgBounds(-1.0, 1.0, -3.0, -2.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );
   cr_expect( checkVarLocks(1, 2, 2, 1) );

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   cr_expect( checkVarLocks(0, 0, 0, 0) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
