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

/**@file   interminor.c
 * @brief  unit test for sepa_interminor.c methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/nlpi_ipopt.h"

#include "scip/sepa_interminor.c"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* w;

/** helper method for setup */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.0, 1.0, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -1.0, 1.0,  1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -1.0, 1.0, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -1.0, 1.0, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPaddVar(scip, w) );
}

/** helper method for teardown */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/** tests the detection of minors; the artificial problem contains five minors: one principle minor for (x,y) and the
 * one corresponding to (x^2, xy, xz, yz), (x^2, xz, xy, yz), (xy, xz, y^2, yz), (xy, y^2, xz, yz)
 */
Test(interminor, detect, .init = setup, .fini = teardown)
{
   #define NCONSS 3
   const char* inputs[NCONSS] = {"[nonlinear] <c1>: 1<= <x> * <x> + <y> * <y> <= 2",
      "[nonlinear] <c2>: -0.5 <= <x> * <y> + <y> * <z> <= 0.5", "[nonlinear] <c3>: -0.5 <= <x> * <z> <= 0.5"};
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int c;

   /* add two expression constraints (in transformed space) */
   for( c = 0; c < 3; ++c )
   {
      SCIP_CALL( SCIPparseCons(scip, &cons, inputs[c], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
         &success) );
      cr_assert(success);

      /* add and release constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* go to solving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );
   cr_assert(SCIPgetNConss(scip) == NCONSS);
   cr_assert(SCIPconshdlrGetNConss(conshdlr) == NCONSS);

   /* make sure INITLP has been run to get auxiliary variables */
   SCIP_CALL( SCIPconstructLP(scip, &infeasible) );
   cr_assert(!infeasible);

   /* get separator data */
   sepa = SCIPfindSepa(scip, SEPA_NAME);
   cr_assert(sepa != NULL);
   sepadata = SCIPsepaGetData(sepa);
   cr_assert(sepadata != NULL);

   /* call minor detection */
   cr_expect(!sepadata->detectedminors);
   cr_expect(sepadata->nminors == 0);
   SCIP_CALL( detectMinors(scip, sepadata) );
   cr_expect(sepadata->detectedminors);
   cr_expect(sepadata->nminors == 5, "nminors = %d (expected 5)", sepadata->nminors);
}
