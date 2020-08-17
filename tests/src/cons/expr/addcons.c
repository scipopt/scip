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

/**@file   addcons.c
 * @brief  tests adding expression constraints during the solving process
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* test that creates and adds a nonchecked expression constraints during SCIP_STAGE_SOLVING */
Test(addcons, nonchecked, .init = setup, .fini = teardown)
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "[expr] <c1>: <t_x> + (<t_x>) == 2";
   SCIP_Bool success;
   SCIP_Bool cutoff;

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* create constraint from input string */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, input,
      TRUE,  /* initial */
      TRUE,  /* separate */
      FALSE, /* enforce */
      FALSE, /* check */
      TRUE,  /* propagate */
      FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add constraint */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   expr = SCIPgetExprConsExpr(scip, cons);
   cr_assert(expr != NULL);

   /* check locks */
   cr_expect(SCIPgetConsExprExprNLocksNeg(expr) == 1);
   cr_expect(SCIPgetConsExprExprNLocksPos(expr) == 1);

   /* expression should have been simplified in SCIP_DECL_CONSACTIVE */
   cr_expect(SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) );
   cr_expect(SCIPgetConsExprExprNChildren(expr) == 1);

   /* call SCIPconstructLP to trigger an INITLP call */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   cr_expect(!cutoff);

   /* check whether expression has been detected by at least one nonlinear handler (happens now already in addCons during solve */
   cr_expect(SCIPgetConsExprExprAuxVar(expr) != NULL);

   /* release the constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
