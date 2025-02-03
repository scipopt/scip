/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   fixedvar.c
 * @brief  unit test that checks that the fixedvar constraint handler is doing something
 * @author Stefan Vigerske
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/struct_cons.h"
#include "include/scip_test.h"

/* TESTS  */

Test(fixedvar, check)
{
   SCIP* scip;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_SOL* sol;
   SCIP_Bool infeas;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;
   SCIP_Bool feasible;

   /* create */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "fixedvar_check") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -1.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 1) ); /* to be able to stop in presolve, it needs to run */

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );

   /* aggregate x = 1.0 + 2e5 y */
   SCIP_CALL( SCIPaggregateVars(scip, x, y, 1.0, -0.2 / SCIPfeastol(scip), 1.0, &infeas, &redundant, &aggregated) );
   cr_expect(!infeas);
   cr_expect(redundant);
   cr_expect(aggregated);

   /* SCIP_CALL( SCIPprintTransProblem(scip, NULL, NULL, FALSE) ); */

   /* solution y = 1e-6 -> x = 1.2 violates its upper bound 1 */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, SCIPfeastol(scip)) );

   /* without cons_fixedvar, this solution is feasible */
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/fixedvar/enabled", FALSE) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, TRUE, TRUE, TRUE, &feasible) );
   cr_expect(feasible);

   /* with cons_fixedvar, it is not feasible */
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/fixedvar/enabled", TRUE) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, TRUE, TRUE, TRUE, &feasible) );
   cr_expect(!feasible);

   /* free */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   SCIP_CALL( SCIPreleaseVar(scip , &y) );
   SCIP_CALL( SCIPreleaseVar(scip , &x) );

   SCIP_CALL( SCIPfree(&scip) );
}

Test(fixedvar, enforce)
{
   SCIP* scip;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_ROW* cut;
   SCIP_Real val = 1.0;
   SCIP_RESULT result;
   SCIP_Bool infeas;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;
   SCIP_Bool feasible;

   /* create */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "fixedvar_enforce") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -1.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* dummy variable and constraint to avoid pseudosolution to be feasible */
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "dummy", 1, &z, &val, 1, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 1) ); /* to be able to stop in presolve, it needs to run */
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );
   /* SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", SCIP_VERBLEVEL_FULL) ); */

   /* aggregate x = 1.0 + 2e5 y */
   SCIP_CALL( SCIPaggregateVars(scip, x, y, 1.0, -0.2 / SCIPfeastol(scip), 1.0, &infeas, &redundant, &aggregated) );
   cr_expect(!infeas);
   cr_expect(redundant);
   cr_expect(aggregated);

   /* SCIP_CALL( SCIPprintTransProblem(scip, NULL, NULL, FALSE) ); */

   /* go to solving stage, but stop there */
   SCIP_CALL( SCIPpresolve(scip) );

   SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1L) );
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );
   SCIP_CALL( SCIPsolve(scip) );

   /* solution y = 1e-6 -> x = 1.2 violates its upper bound 1 */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, SCIPfeastol(scip)) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1.0) );

   /* due to cons_fixedvar, this is not feasible */
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, TRUE, TRUE, TRUE, &feasible) );
   cr_expect(!feasible);
   cr_expect(SCIPgetNCuts(scip) == 0);

   /* enforce solution via the enforelax callback of cons_fixedvar */
   conshdlr = SCIPfindConshdlr(scip, "fixedvar");
   cr_assert_not_null(conshdlr);
   SCIP_CALL( conshdlr->consenforelax(scip, sol, conshdlr, NULL, 0, 0, FALSE, &result) );

   /* we should now have 1 cut to enforce bounds [0,1] on x: 0 <= 2e5 y + 1 <= 1 */
   cr_expect_eq(result, SCIP_SEPARATED);
   cr_assert(SCIPgetNCuts(scip) == 1);
   cut = SCIPgetCuts(scip)[0];
   SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
   cr_expect_eq(SCIProwGetLhs(cut), 0.0);
   cr_expect_eq(SCIProwGetConstant(cut), 1.0);
   cr_expect_eq(SCIProwGetRhs(cut), 1.0);
   cr_assert(SCIProwGetNNonz(cut) == 1);
   cr_expect_float_eq(SCIProwGetVals(cut)[0], 0.2 / SCIPfeastol(scip), SCIPepsilon(scip));
   cr_expect_eq(SCIProwGetCols(cut)[0], SCIPvarGetCol(SCIPvarGetTransVar(y)));

   /* free */
   SCIP_CALL( SCIPclearCuts(scip) );

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   SCIP_CALL( SCIPreleaseVar(scip , &z) );
   SCIP_CALL( SCIPreleaseVar(scip , &y) );
   SCIP_CALL( SCIPreleaseVar(scip , &x) );

   SCIP_CALL( SCIPfree(&scip) );
}
