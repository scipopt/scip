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

/**@file   unittest-cons_epxr.c
 * @brief  unit test for cons_epxr methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

/* creates scip, problem, includes expression constraint handler, creates  and adds variables */
static
void setup(void)
{
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 5.0, 2.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 5.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();
}

Test(conshdlr, conscheck, .init = setup, .fini = teardown,
   .description = "test feasibility check of the cons_expr constraint handler."
   )
{
   SCIP_SOL* sol;
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   const char* input = "[expr] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   /* create an infeasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, FALSE, FALSE, FALSE, &success) );
   if( success )
   {
      cr_assert(FALSE, "an ifeasible solution has been accepted");
   }

   /* create a feasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, FALSE, FALSE, FALSE, &success) );
   if( !success )
   {
      cr_assert(FALSE, "a feasible solution has been declined");
   }

   /* create an undefined solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, FALSE, FALSE, FALSE, &success) );
   if( success )
   {
      cr_assert(FALSE, "undefined solution has been accepted");
   }

   /* release solution */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}

