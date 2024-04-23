/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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
#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include <stdio.h>

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

   /* aggregate x = 1.0 + 2e6 y */
   SCIP_CALL( SCIPaggregateVars(scip, x, y, 1.0, -2.0/SCIPfeastol(scip), 1.0, &infeas, &redundant, &aggregated) );
   cr_expect(!infeas);
   cr_expect(redundant);
   cr_expect(aggregated);

   /* SCIP_CALL( SCIPprintTransProblem(scip, NULL, NULL, FALSE) ); */

   /* solution y = 1e-6 -> x = 3 violates its upper bound 1 */
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

   SCIP_CALL( SCIPfree(&scip ) );
}
