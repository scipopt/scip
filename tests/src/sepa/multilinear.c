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

/**@file   multilinear.c
 * @brief  unit tests for multilinear cut separator
 * @author Matthias Walter
 */

#include <scip/scip.h>
#include <include/scip_test.h>
#include <scip/scipdefplugins.h>

/* global SCIP instance */
static SCIP* scip;

/** setup: create SCIP */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* turn off presolving, restarts, heuristics and other separators in order to avoid too clever simplifications .*/
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "branching/mostinf/priority", 100000) );
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/scaleobj", FALSE) );
   SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 0) );
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/and/sepafreq", 1) );
}

/** teardown: free SCIP */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!");
}

/** simple example with 4 variables and 2 linear constraints */
static
void hypergraph1flower(
   SCIP_Bool             separate            /**< Whether to generate 1-flowers. */
   )
{
   SCIP_VAR* x1;
   SCIP_VAR* x2;
   SCIP_VAR* x3;
   SCIP_VAR* x4;
   SCIP_VAR* y123;
   SCIP_VAR* y234;
   SCIP_VAR* y13;
   SCIP_VAR* y14;
   SCIP_VAR* y23;
   SCIP_VAR* y24;
   SCIP_VAR* y34;
   SCIP_VAR* z1;
   SCIP_VAR* z2;
   SCIP_VAR* z3;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];

   SCIP_CONS* cons;


   /* setup problem:
    * min y123 + 2y234 - y13 + 2y14 - y23 + y24 - 2y34 - x1 - 2x2 + x3 - 2x4  - z1 - z2 - z3
    *     y123 = x1 AND x2 AND X3
    *     y234 = x2 AND x3 AND x4
    *     y13 = x1 AND x3
    *     y14 = x1 AND x4
    *     y23 = x2 AND x3
    *     y24 = x2 AND x4
    *     y34 = x3 AND x4
    *     z1 + z2 <= 1
    *     z1 + z3 <= 1
    *     z2 + z3 <= 1
    *     all binary
    *
    * The only purpose of the z-variables is to avoid that the problem is solved "in the root node" by branching.
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "1flower"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", 0.0, 1.0, -2.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, x3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4, "x4", 0.0, 1.0, -2.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, x4) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &y123, "y123", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, y123) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y234, "y234", 0.0, 1.0, 2.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, y234) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y13, "y13", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, y13) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y14, "y14", 0.0, 1.0, 2.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, y14) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y23, "y23", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, y23) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y24, "y24", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, y24) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y34, "y34", 0.0, 1.0, -2.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, y34) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &z1, "z1", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, z1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z2, "z2", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, z2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z3, "z3", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, z3) );

   vars[0] = x1;
   vars[1] = x2;
   vars[2] = x3;
   SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, "and123", y123, 3, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = x2;
   vars[1] = x3;
   vars[2] = x4;
   SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, "and234", y234, 3, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = x1;
   vars[1] = x3;
   SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, "and13", y13, 2, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = x1;
   vars[1] = x4;
   SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, "and14", y14, 2, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = x2;
   vars[1] = x3;
   SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, "and23", y23, 2, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = x2;
   vars[1] = x4;
   SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, "and24", y24, 2, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = x3;
   vars[1] = x4;
   SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, "and34", y34, 2, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vals[0] = 1.0;
   vals[1] = 1.0;

   vars[0] = z1;
   vars[1] = z2;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "stable12", 2, vars, vals, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = z1;
   vars[1] = z3;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "stable13", 2, vars, vals, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = z2;
   vars[1] = z3;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "stable23", 2, vars, vals, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   if( separate )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "separating/multilinear/freq", 1) );
      SCIP_CALL( SCIPsetIntParam(scip, "separating/multilinear/maxtwoflower", 0) );
   }

   /* solve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPsolve(scip) );

#if SCIP_DEBUG
   SCIPprintSeparatorStatistics(scip, stdout);
   SCIPprintRootStatistics(scip, stdout);
#endif /* SCIP_DEBUG */

   SCIP_Real rootbound = SCIPgetDualboundRoot(scip);
   if( separate )
      cr_assert_float_eq(rootbound, -5.5, 1.0e-8);
   else
      cr_assert_float_eq(rootbound, -6.0, 1.0e-8);

   SCIP_CALL( SCIPreleaseVar(scip, &z3) );
   SCIP_CALL( SCIPreleaseVar(scip, &z2) );
   SCIP_CALL( SCIPreleaseVar(scip, &z1) );
   SCIP_CALL( SCIPreleaseVar(scip, &y34) );
   SCIP_CALL( SCIPreleaseVar(scip, &y24) );
   SCIP_CALL( SCIPreleaseVar(scip, &y23) );
   SCIP_CALL( SCIPreleaseVar(scip, &y14) );
   SCIP_CALL( SCIPreleaseVar(scip, &y13) );
   SCIP_CALL( SCIPreleaseVar(scip, &y234) );
   SCIP_CALL( SCIPreleaseVar(scip, &y123) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
}


/* TEST SUITE */
TestSuite(test_sepa_multilinear, .init = setup, .fini = teardown);

/* TEST 1 */
Test(test_sepa_multilinear, 1flower, .description = "trigger separation for a simple hypergraph")
{
   hypergraph1flower(false);

   hypergraph1flower(true);
}

