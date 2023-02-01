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

/**@file   addcons.c
 * @brief  tests adding nonlinear constraints during the solving process
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );

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

/* test that creates and adds a nonchecked nonlinear constraints during SCIP_STAGE_SOLVING */
Test(addcons, nonchecked, .init = setup, .fini = teardown)
{
   SCIP_CONS* cons;
   SCIP_CONS* t_cons;
   SCIP_EXPR* expr;
   const char* input = "[nonlinear] <c1>: <x> + (<x>) == 2";
   SCIP_Bool success;
   SCIP_Bool cutoff;

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

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPgetTransformedCons(scip, cons, &t_cons) );
   expr = SCIPgetExprNonlinear(t_cons);
   cr_assert(expr != NULL);

   /* check locks */
   cr_expect(SCIPgetExprNLocksNegNonlinear(expr) == 1);
   cr_expect(SCIPgetExprNLocksPosNonlinear(expr) == 1);

   /* expression should have been simplified in SCIP_DECL_CONSACTIVE */
   cr_expect(SCIPisExprSum(scip, expr) );
   cr_expect(SCIPexprGetNChildren(expr) == 1);

   /* call SCIPconstructLP to trigger an INITLP call */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   cr_expect(!cutoff);

   /* check whether expression has been detected by at least one nonlinear handler (happens now already in addCons during solve */
   cr_expect(SCIPgetExprAuxVarNonlinear(expr) != NULL);

   /* release the constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
