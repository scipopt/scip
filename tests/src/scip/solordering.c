/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   solordering.c
 * @brief  regression tests for ordering solutions with indistinguishable floating-point objectives
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

#ifdef SCIP_WITH_EXACTSOLVE

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPenableExactSolving(scip, TRUE) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "solution_order") );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 2.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVarExactData(scip, x, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVarExactData(scip, y, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0);
}

static
SCIP_SOL* makeSolution(
   const char*           objective,          /**< value of x, and hence of the objective */
   SCIP_Bool             yvalue              /**< distinguish the points independently of the objective */
   )
{
   SCIP_SOL* sol;
   SCIP_RATIONAL* value;

   SCIP_CALL( SCIPrationalCreate(&value) );
   SCIP_CALL( SCIPcreateSolExact(scip, &sol, NULL) );
   SCIPrationalSetString(value, objective);
   SCIP_CALL( SCIPsetSolValExact(scip, sol, x, value) );
   SCIPrationalSetReal(value, (SCIP_Real)yvalue);
   SCIP_CALL( SCIPsetSolValExact(scip, sol, y, value) );
   SCIPrationalFree(&value);
   return sol;
}

TestSuite(solordering, .init = setup, .fini = teardown);

/* Both exact objectives have the same floating approximation, also with upward rounding. */
#define BETTER "100000000000000000001/100000000000000000000"
#define WORSE  "100000000000000000002/100000000000000000000"

static
void checkBestObjective(void)
{
   SCIP_RATIONAL* expected;
   SCIP_RATIONAL* actual;

   SCIP_CALL( SCIPrationalCreate(&expected) );
   SCIP_CALL( SCIPrationalCreate(&actual) );
   SCIPrationalSetString(expected, BETTER);
   SCIPgetSolOrigObjExact(scip, SCIPgetBestSol(scip), actual);
   cr_expect(SCIPrationalIsEQ(actual, expected));
   SCIPrationalFree(&actual);
   SCIPrationalFree(&expected);
}

static
void presolve(void)
{
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPpresolve(scip) );
   cr_assert_eq(SCIPgetStage(scip), SCIP_STAGE_PRESOLVED);
}

Test(solordering, original_pool)
{
   SCIP_SOL* better;
   SCIP_SOL* worse;
   SCIP_Bool stored;

   worse = makeSolution(WORSE, FALSE);
   better = makeSolution(BETTER, TRUE);
   cr_assert_eq(SCIPgetSolOrigObj(scip, worse), SCIPgetSolOrigObj(scip, better));
   SCIP_CALL( SCIPaddSol(scip, worse, &stored) );
   cr_assert(stored);
   SCIP_CALL( SCIPaddSol(scip, better, &stored) );
   cr_assert(stored);
   checkBestObjective();
   SCIP_CALL( SCIPfreeSol(scip, &better) );
   SCIP_CALL( SCIPfreeSol(scip, &worse) );
}

Test(solordering, transformed_pool)
{
   SCIP_SOL* better;
   SCIP_SOL* worse;
   SCIP_Bool stored;

   presolve();
   worse = makeSolution(WORSE, FALSE);
   better = makeSolution(BETTER, TRUE);
   SCIP_CALL( SCIPaddSol(scip, worse, &stored) );
   cr_assert(stored);
   SCIP_CALL( SCIPaddSol(scip, better, &stored) );
   cr_assert(stored);
   checkBestObjective();
   SCIP_CALL( SCIPfreeSol(scip, &better) );
   SCIP_CALL( SCIPfreeSol(scip, &worse) );
}

Test(solordering, original_pool_survives_presolve)
{
   SCIP_SOL* sols[2];
   SCIP_Bool stored;

   sols[0] = makeSolution(WORSE, FALSE);
   sols[1] = makeSolution(BETTER, TRUE);
   SCIP_CALL( SCIPaddSol(scip, sols[0], &stored) );
   cr_assert(stored);
   SCIP_CALL( SCIPaddSol(scip, sols[1], &stored) );
   cr_assert(stored);
   SCIP_CALL( SCIPfreeSol(scip, &sols[1]) );
   SCIP_CALL( SCIPfreeSol(scip, &sols[0]) );
   presolve();
   checkBestObjective();
}

Test(solordering, transformed_pool_limit_one)
{
   SCIP_SOL* sol;
   SCIP_Bool stored;

   SCIP_CALL( SCIPsetIntParam(scip, "limits/maxsol", 1) );
   presolve();
   sol = makeSolution(WORSE, FALSE);
   SCIP_CALL( SCIPaddSol(scip, sol, &stored) );
   cr_assert(stored);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   sol = makeSolution(BETTER, TRUE);
   SCIP_CALL( SCIPaddSol(scip, sol, &stored) );
   cr_assert(stored, "An exactly better solution must not be discarded when the pool is full.");
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   cr_expect_eq(SCIPgetNSols(scip), 1);
   checkBestObjective();
}

Test(solordering, transformed_preference_requires_exact_tie)
{
   SCIP_SOL* original;
   SCIP_SOL* transformed;
   SCIP_Bool stored;

   presolve();
   original = makeSolution(BETTER, TRUE);
   SCIP_CALL( SCIPretransformSolExact(scip, original) );
   SCIP_CALL( SCIPsetBoolParam(scip, "exact/improvingsols", FALSE) );
   SCIP_CALL( SCIPaddSol(scip, original, &stored) );
   cr_assert(stored);
   transformed = makeSolution(WORSE, FALSE);
   SCIP_CALL( SCIPaddSol(scip, transformed, &stored) );
   cr_assert(stored);
   checkBestObjective();
   SCIP_CALL( SCIPfreeSol(scip, &transformed) );
   SCIP_CALL( SCIPfreeSol(scip, &original) );
}

Test(solordering, transformed_preference_on_exact_tie)
{
   SCIP_SOL* original;
   SCIP_SOL* transformed;
   SCIP_Bool stored;

   presolve();
   SCIP_CALL( SCIPsetBoolParam(scip, "exact/improvingsols", FALSE) );
   original = makeSolution(BETTER, TRUE);
   SCIP_CALL( SCIPretransformSolExact(scip, original) );
   SCIP_CALL( SCIPaddSol(scip, original, &stored) );
   cr_assert(stored);
   transformed = makeSolution(BETTER, FALSE);
   SCIP_CALL( SCIPaddSol(scip, transformed, &stored) );
   cr_assert(stored);
   checkBestObjective();
   cr_expect_not(SCIPsolIsOriginal(SCIPgetBestSol(scip)));
   SCIP_CALL( SCIPfreeSol(scip, &transformed) );
   SCIP_CALL( SCIPfreeSol(scip, &original) );
}

#endif
