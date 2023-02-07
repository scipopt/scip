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

/**@file   compute.c
 * @brief  unit tests for computing symmetry
 * @author Marc Pfetsch
 * @author Fabian Wegscheider
 */

#include <scip/scip.h>
#include <include/scip_test.h>
#include <scip/prop_symmetry.c>
#include <scip/symmetry.h>
#include <symmetry/compute_symmetry.h>
#include <scip/scipdefplugins.h>

/* global SCIP instance */
static SCIP* scip;


/** check whether two int arrays are equal */
static
void checkIntArraysEqual(
   int*                  expected,           /**< array of expected values */
   int*                  candidate,          /**< array of values to be checked */
   int                   length,             /**< length of arrays */
   const char*           name                /**< name to be printed */
   )
{
   int i;

   for( i = 0; i < length; ++i )
      cr_expect(expected[i] == candidate[i], "%s[%d]: expected %d, but got %d\n", name, i, expected[i], candidate[i]);
}

/** setup: create SCIP */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* turn on symmetry computation */
   SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 1) );

#ifdef SCIP_DEBUG
   /* output external codes in order to see which external symmetry computation code is used */
   SCIPprintExternalCodes(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
#endif
}

/** teardown: free SCIP */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!");
}

/* TEST SUITE */
TestSuite(test_compute_symmetry, .init = setup, .fini = teardown);

/* TEST 1 */
Test(test_compute_symmetry, basic1, .description = "compute symmetry for a simple example with 4 variables and 2 linear constraints")
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_VAR** permvars;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int* components;
   int* componentbegins;
   int* vartocomponent;
   int ncomponents;
   int norbits;
   int npermvars;
   int nperms;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 + x4
    *     x1 + x2           = 1
    *               x3 + x4 = 1
    *     x1, ..., x4 binary
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "basic1"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var4) );

   vars[0] = var1;
   vars[1] = var2;
   vals[0] = 1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e1", 2, vars, vals, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var3;
   vars[1] = var4;
   vals[0] = 1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e2", 2, vars, vals, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* turn off subgroup detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/detectsubgroups", FALSE) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );
   cr_assert( nperms == 3 );
   cr_assert( ncomponents == 1 );
   cr_assert( componentbegins[0] == 0 );
   cr_assert( componentbegins[1] == 3 );
   cr_assert( vartocomponent[0] == 0 );
   cr_assert( vartocomponent[1] == 0 );
   cr_assert( vartocomponent[2] == 0 );
   cr_assert( vartocomponent[3] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   cr_assert( norbits == 1 );
   cr_assert( orbitbegins[0] == 0 );
   cr_assert( orbitbegins[1] == 4 );
   cr_assert( orbits[0] == 0 );
   cr_assert( orbits[1] == 1 );
   cr_assert( orbits[2] == 2 );
   cr_assert( orbits[3] == 3 );
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

#ifdef SCIP_DEBUG
   {
      int i;
      int j;

      for (i = 0; i < nperms; ++i)
      {
         SCIPinfoMessage(scip, NULL, "Permutation %d: (", i);
         for (j = 0; j < npermvars; ++j)
         {
            if ( j == 0 )
               SCIPinfoMessage(scip, NULL, "%d", perms[i][j]);
            else
               SCIPinfoMessage(scip, NULL, " %d", perms[i][j]);
         }
         SCIPinfoMessage(scip, NULL, ")\n");
      }
   }
#endif

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}


/* TEST 2 */
Test(test_compute_symmetry, basic2, .description = "compute symmetry for a simple example with 4 variables and 4 linear constraints - before presolving")
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_VAR** permvars;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int* components;
   int* componentbegins;
   int* vartocomponent;
   int ncomponents;
   int norbits;
   int npermvars;
   int nperms;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 + x4
    *     x1 + x2           =  1
    *               x3 + x4 =  1
    *    2x1 +           x4 <= 2
    *         2x2 + x3      <= 2
    *     x1, ..., x4 binary
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "basic2"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var4) );

   vars[0] = var1;
   vars[1] = var2;
   vals[0] = 1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e1", 2, vars, vals, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var3;
   vars[1] = var4;
   vals[0] = 1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e2", 2, vars, vals, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var1;
   vars[1] = var4;
   vals[0] = 2.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "i1", 2, vars, vals, -SCIPinfinity(scip), 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var2;
   vars[1] = var3;
   vals[0] = 2.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "i2", 2, vars, vals, -SCIPinfinity(scip), 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* turn off subgroup detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/detectsubgroups", FALSE) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );
   cr_assert( nperms == 1 );
   cr_assert( ncomponents == 1 );
   cr_assert( componentbegins[0] == 0 );
   cr_assert( componentbegins[1] == 1 );
   cr_assert( vartocomponent[0] == 0 );
   cr_assert( vartocomponent[1] == 0 );
   cr_assert( vartocomponent[2] == 0 );
   cr_assert( vartocomponent[3] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   cr_assert( norbits == 2 );
   cr_assert( orbitbegins[0] == 0 );
   cr_assert( orbitbegins[1] == 2 );
   cr_assert( orbits[0] == 0 );
   cr_assert( orbits[1] == 1 );
   cr_assert( orbits[2] == 2 );
   cr_assert( orbits[3] == 3 );
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

#ifdef SCIP_DEBUG
   {
      int i;
      int j;

      for (i = 0; i < nperms; ++i)
      {
         SCIPinfoMessage(scip, NULL, "Permutation %d: (", i);
         for (j = 0; j < npermvars; ++j)
         {
            if ( j == 0 )
               SCIPinfoMessage(scip, NULL, "%d", perms[i][j]);
            else
               SCIPinfoMessage(scip, NULL, " %d", perms[i][j]);
         }
         SCIPinfoMessage(scip, NULL, ")\n");
      }
   }
#endif

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}


/* TEST 3 */
Test(test_compute_symmetry, basic3, .description = "compute symmetry for a simple example with 4 variables and 4 linear constraints - after presolving")
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_VAR** permvars;
   int** perms;
   int* components;
   int* componentbegins;
   int* vartocomponent;
   int ncomponents;
   int nperms;
   int npermvars;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 + x4
    *     x1 + x2           =  1
    *               x3 + x4 =  1
    *    2x1 +           x4 <= 2
    *         2x2 + x3      <= 2
    *     x1, ..., x4 binary
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "basic3"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var4) );

   vars[0] = var1;
   vars[1] = var2;
   vals[0] = 1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e1", 2, vars, vals, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var3;
   vars[1] = var4;
   vals[0] = 1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e2", 2, vars, vals, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var1;
   vars[1] = var4;
   vals[0] = 2.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "i1", 2, vars, vals, -SCIPinfinity(scip), 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var2;
   vars[1] = var3;
   vals[0] = 2.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "i2", 2, vars, vals, -SCIPinfinity(scip), 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* turn off subgroup detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/detectsubgroups", FALSE) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );
   cr_assert( nperms == -1 );  /* problem should be empty */
   cr_assert( ncomponents == -1 );
   cr_assert( components == NULL );
   cr_assert( componentbegins == NULL );
   cr_assert( vartocomponent == NULL );

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}


/* TEST 4 */
Test(test_compute_symmetry, basic4, .description = "compute symmetry for a simple example with 5 variables and 2 linear constraints")
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_VAR* var5;
   SCIP_CONS* cons;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];
   SCIP_VAR** permvars;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int* components;
   int* componentbegins;
   int* vartocomponent;
   int ncomponents;
   int norbits;
   int npermvars;
   int nperms;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 + x4
    *     x1 + x2           + x5 = 1
    *               x3 + x4 + x5 = 2
    *     x1, ..., x4, x5  binary
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "basic4"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var4) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var5, "x5", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var5) );

   vars[0] = var1;
   vars[1] = var2;
   vars[2] = var5;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e1", 3, vars, vals, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var3;
   vars[1] = var4;
   vars[2] = var5;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e2", 3, vars, vals, 2.0, 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* turn off subgroup detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/detectsubgroups", FALSE) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );
   cr_assert( nperms == 2 );
   cr_assert( ncomponents == 2 );
   cr_assert( vartocomponent[0] == vartocomponent[1] );
   cr_assert( vartocomponent[2] == vartocomponent[3] );
   cr_assert( vartocomponent[0] != vartocomponent[2] );
   cr_assert( vartocomponent[1] != vartocomponent[3] );
   cr_assert( vartocomponent[4] == -1 );
   cr_assert( componentbegins[0] == 0 );
   cr_assert( componentbegins[1] == 1 );
   cr_assert( componentbegins[2] == 2 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   cr_assert( norbits == 2 );
   cr_assert( orbitbegins[0] == 0 );
   cr_assert( orbitbegins[1] == 2 );
   cr_assert( orbitbegins[2] == 4 );
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

#ifdef SCIP_DEBUG
   {
      int i;
      int j;

      for (i = 0; i < nperms; ++i)
      {
         SCIPinfoMessage(scip, NULL, "Permutation %d: (", i);
         for (j = 0; j < npermvars; ++j)
         {
            if ( j == 0 )
               SCIPinfoMessage(scip, NULL, "%d", perms[i][j]);
            else
               SCIPinfoMessage(scip, NULL, " %d", perms[i][j]);
         }
         SCIPinfoMessage(scip, NULL, ")\n");
      }
   }
#endif

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
   SCIP_CALL( SCIPreleaseVar(scip, &var5) );
}


/* TEST 5 */
Test(test_compute_symmetry, basic5, .description = "compute symmetry for a simple example with bounddisjunctions")
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_BOUNDTYPE boundtypes[2];
   SCIP_Real bounds[2];
   SCIP_VAR** permvars;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int npermvars;
   int nperms;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 + x4
    *     BD(x1 >= 1, x2 >= 1)
    *     BD(x3 >= 1, x4 >= 1)
    *     x1, ..., x4 binary
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "basic1"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var4) );

   vars[0] = var1;
   vars[1] = var2;
   boundtypes[0] = SCIP_BOUNDTYPE_LOWER;
   boundtypes[1] = SCIP_BOUNDTYPE_LOWER;
   bounds[0] = 1.0;
   bounds[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicBounddisjunction(scip, &cons, "d1", 2, vars, boundtypes, bounds) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var3;
   vars[1] = var4;
   boundtypes[0] = SCIP_BOUNDTYPE_LOWER;
   boundtypes[1] = SCIP_BOUNDTYPE_LOWER;
   bounds[0] = 1.0;
   bounds[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicBounddisjunction(scip, &cons, "d2", 2, vars, boundtypes, bounds) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* turn off subgroup detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/detectsubgroups", FALSE) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
#endif

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL) );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   cr_assert( norbits == 1 );
   cr_assert( orbitbegins[0] == 0 );
   cr_assert( orbitbegins[1] == 4 );
   cr_assert( orbits[0] == 0 );
   cr_assert( orbits[1] == 1 );
   cr_assert( orbits[2] == 2 );
   cr_assert( orbits[3] == 3 );
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

#ifdef SCIP_DEBUG
   {
      int i;
      int j;

      for (i = 0; i < nperms; ++i)
      {
         SCIPinfoMessage(scip, NULL, "Permutation %d: (", i);
         for (j = 0; j < npermvars; ++j)
         {
            if ( j == 0 )
               SCIPinfoMessage(scip, NULL, "%d", perms[i][j]);
            else
               SCIPinfoMessage(scip, NULL, " %d", perms[i][j]);
         }
         SCIPinfoMessage(scip, NULL, ")\n");
      }
   }
#endif

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}


/* TEST 6 */
Test(test_compute_symmetry, basic6, .description = "compute symmetry for a simple example with additional bounddisjunctions")
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_VAR* var5;
   SCIP_CONS* cons;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];
   SCIP_BOUNDTYPE boundtypes[2];
   SCIP_Real bounds[2];
   SCIP_VAR** permvars;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int npermvars;
   int nperms;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 + x4
    *     x1 + x2           + x5 =  3
    *               x3 + x4 + x5 =  3
    *    2x1 +           x4      <= 2
    *         2x2 + x3           <= 2
    *    BD(x5 <= 1, x5 >= 3)
    *     x1, ..., x4 binary, x5 continuous
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "basic5"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var4) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var5, "x5", 0.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var5) );

   vars[0] = var1;
   vars[1] = var2;
   vars[2] = var5;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e1", 3, vars, vals, 3.0, 3.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var3;
   vars[1] = var4;
   vars[2] = var5;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "e2", 3, vars, vals, 3.0, 3.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var1;
   vars[1] = var4;
   vals[0] = 2.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "i1", 2, vars, vals, -SCIPinfinity(scip), 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var2;
   vars[1] = var3;
   vals[0] = 2.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "i2", 2, vars, vals, -SCIPinfinity(scip), 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var5;
   vars[1] = var5;
   boundtypes[0] = SCIP_BOUNDTYPE_UPPER;
   boundtypes[1] = SCIP_BOUNDTYPE_LOWER;
   bounds[0] = 1.0;
   bounds[1] = 3.0;
   SCIP_CALL( SCIPcreateConsBasicBounddisjunction(scip, &cons, "d1", 2, vars, boundtypes, bounds) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* turn off subgroup detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/detectsubgroups", FALSE) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
#endif

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL) );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   cr_assert( nperms == 1 );
   cr_assert( norbits == 2 );
   cr_assert( orbitbegins[0] == 0 );
   cr_assert( orbitbegins[1] == 2 );
   cr_assert( orbits[0] == 0 );
   cr_assert( orbits[1] == 1 );
   cr_assert( orbits[2] == 2 );
   cr_assert( orbits[3] == 3 );
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

#ifdef SCIP_DEBUG
   {
      int i;
      int j;

      for (i = 0; i < nperms; ++i)
      {
         SCIPinfoMessage(scip, NULL, "Permutation %d: (", i);
         for (j = 0; j < npermvars; ++j)
         {
            if ( j == 0 )
               SCIPinfoMessage(scip, NULL, "%d", perms[i][j]);
            else
               SCIPinfoMessage(scip, NULL, " %d", perms[i][j]);
         }
         SCIPinfoMessage(scip, NULL, ")\n");
      }
   }
#endif

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
   SCIP_CALL( SCIPreleaseVar(scip, &var5) );
}

/* TEST 7 (subgroups) */
Test(test_compute_symmetry, subgroups1, .description = "detect symmetric subgroups for artificial propdata")
{
   SCIP_PROPDATA propdata;
   SCIP_VAR* dummyvar;
   SCIP_VAR* permvars[10];
   unsigned componentblocked = FALSE;
   int* graphcomponents;
   int* graphcompbegins;
   int* compcolorbegins;
   int* usedperms;
   int* perms[6];
   SCIP_Shortbool permused[6] = {FALSE};
   int permorder1[6] = {0,1,2,3,4,5};
   int permorder2[6] = {2,3,4,5,0,1};
   int permorder3[6] = {5,0,1,2,3,4};
   int perm1[10] = {1,0,2,4,3,5,6,7,8,9};
   int perm2[10] = {0,2,1,3,5,4,6,7,8,9};
   int perm3[10] = {0,1,2,3,4,5,7,6,8,9};
   int perm4[10] = {0,1,2,3,4,5,6,8,7,9};
   int perm5[10] = {0,1,2,3,4,5,6,7,9,8};
   int perm6[10] = {6,7,2,8,9,5,0,1,3,4};
   int components[6] = {0,1,2,3,4,5};
   int componentbegins[2] = {0,6};
   int expectedcomps[10] = {0,1,2,3,4,5,6,7,8,9};
   int expectedcompbegins[4] = {0,3,6,10};
   int expectedcolbegins[3] = {0,2,3};
   int expectedcomps1[10] = {0,1,2,3,4,5,6,7,8,9};
   int expectedcompbegins1[4] = {0,3,6,10};
   int expectedcolbegins1[3] = {0,2,3};
   int expectedcomps2[10] = {0,6,2,1,7,3,8,4,5,9};
   int expectedcompbegins2[5] = {0,2,5,7,10};
   int expectedcolbegins2[2] = {0,4};
   int ngraphcomponents;
   int ncompcolors;
   int nusedperms;
   int i;

   SCIP_CALL( SCIPcreateProbBasic(scip, "subgroup1"));

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   SCIP_CALL( SCIPcreateVarBasic(scip, &dummyvar, "dummyvar", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   for (i = 0; i < 10; ++i)
      permvars[i] = dummyvar;

   perms[0] = perm1;
   perms[1] = perm2;
   perms[2] = perm3;
   perms[3] = perm4;
   perms[4] = perm5;
   perms[5] = perm6;

   propdata.npermvars = 10;
   propdata.nperms = 6;
   propdata.perms = perms;
   propdata.ncomponents = 1;
   propdata.components = components;
   propdata.componentbegins = componentbegins;
   propdata.componentblocked = &componentblocked;
   propdata.computedsymmetry = TRUE;
   propdata.permvars = permvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &usedperms, 6) );

   /* check canonical order */
   for (i = 0; i < 6; ++i)
      permused[i] = FALSE;
   SCIP_CALL( buildSubgroupGraph(scip, &propdata, permorder1, 6, 0, &graphcomponents, &graphcompbegins,
         &compcolorbegins, &ngraphcomponents, &ncompcolors, &usedperms, &nusedperms, 6, permused) );

   cr_assert( graphcomponents != NULL );
   cr_assert( graphcompbegins != NULL );
   cr_assert( compcolorbegins != NULL );
   cr_assert( nusedperms == 5, "expected 5 used permutations, but got %d\n", nusedperms );
   cr_assert( ngraphcomponents == 3, "expected 3 graph components, but got %d\n", ngraphcomponents );
   cr_assert( ncompcolors == 2, "expected 2 component colors, but got %d\n", ncompcolors );

   checkIntArraysEqual(expectedcomps, graphcomponents, 10, "components");
   checkIntArraysEqual(expectedcompbegins, graphcompbegins, ngraphcomponents, "compbegins");
   checkIntArraysEqual(expectedcolbegins, compcolorbegins, ncompcolors, "colorbegins");

   SCIPfreeBlockMemoryArray(scip, &compcolorbegins, ncompcolors + 1);
   SCIPfreeBlockMemoryArray(scip, &graphcompbegins, ngraphcomponents + 1);
   SCIPfreeBlockMemoryArray(scip, &graphcomponents, 10);

   /* check different order */
   for (i = 0; i < 6; ++i)
      permused[i] = FALSE;
   SCIP_CALL( buildSubgroupGraph(scip, &propdata, permorder2, 6, 0, &graphcomponents, &graphcompbegins,
         &compcolorbegins, &ngraphcomponents, &ncompcolors, &usedperms, &nusedperms, 6, permused) );

   cr_assert( graphcomponents != NULL );
   cr_assert( graphcompbegins != NULL );
   cr_assert( compcolorbegins != NULL );
   cr_assert( nusedperms == 5, "expected 5 used permutations, but got %d\n", nusedperms );
   cr_assert( ngraphcomponents == 3, "expected 3 graph components, but got %d\n", ngraphcomponents );
   cr_assert( ncompcolors == 2, "expected 2 component colors, but got %d\n", ncompcolors );

   checkIntArraysEqual(expectedcomps1, graphcomponents, 10, "components");
   checkIntArraysEqual(expectedcompbegins1, graphcompbegins, ngraphcomponents, "compbegins");
   checkIntArraysEqual(expectedcolbegins1, compcolorbegins, ncompcolors, "colorbegins");

   SCIPfreeBlockMemoryArray(scip, &compcolorbegins, ncompcolors + 1);
   SCIPfreeBlockMemoryArray(scip, &graphcompbegins, ngraphcomponents + 1);
   SCIPfreeBlockMemoryArray(scip, &graphcomponents, 10);

   /* check order that leads to trivial subgroup */
   for (i = 0; i < 6; ++i)
      permused[i] = FALSE;
   SCIP_CALL( buildSubgroupGraph(scip, &propdata, permorder3, 6, 0, &graphcomponents, &graphcompbegins,
         &compcolorbegins, &ngraphcomponents, &ncompcolors, &usedperms, &nusedperms, 6, permused) );

   cr_assert( graphcomponents != NULL );
   cr_assert( graphcompbegins != NULL );
   cr_assert( compcolorbegins != NULL );
   cr_assert( nusedperms == 2, "expected 2 used permutations, but got %d\n", nusedperms );
   cr_assert( ngraphcomponents == 4, "expected 4 graph components, but got %d\n", ngraphcomponents );
   cr_assert( ncompcolors == 1, "expected 1 component colors, but got %d\n", ncompcolors );

   checkIntArraysEqual(expectedcomps2, graphcomponents, 10, "components");
   checkIntArraysEqual(expectedcompbegins2, graphcompbegins, ngraphcomponents, "compbegins");
   checkIntArraysEqual(expectedcolbegins2, compcolorbegins, ncompcolors, "colorbegins");

   SCIPfreeBlockMemoryArray(scip, &compcolorbegins, ncompcolors + 1);
   SCIPfreeBlockMemoryArray(scip, &graphcompbegins, ngraphcomponents + 1);
   SCIPfreeBlockMemoryArray(scip, &graphcomponents, 10);
   SCIPfreeBufferArray(scip, &usedperms);
   SCIP_CALL( SCIPreleaseVar(scip, &dummyvar) );
}

/* TEST 8 (subgroups) */
Test(test_compute_symmetry, subgroups2, .description = "detect symmetric subgroups for artificial propdata and different order")
{
   SCIP_PROPDATA propdata;
   SCIP_VAR* dummyvar;
   SCIP_VAR* permvars[10];
   unsigned componentblocked = FALSE;
   int* permorder;
   int* graphcomponents;
   int* graphcompbegins;
   int* compcolorbegins;
   int* usedperms;
   int* perms[6];
   SCIP_Shortbool permused[6] = {FALSE};
   int perm1[10] = {0,2,1,3,5,4,6,7,8,9};
   int perm2[10] = {0,1,2,3,4,5,7,6,8,9};
   int perm3[10] = {0,1,2,3,4,5,6,8,7,9};
   int perm4[10] = {6,7,0,8,9,5,2,1,3,4};
   int perm5[10] = {1,0,2,4,3,5,6,7,8,9};
   int perm6[10] = {0,1,2,3,4,5,6,7,9,8};
   int components[6] = {0,1,2,3,4,5};
   int componentbegins[2] = {0,6};
   int expectedpermorder[6] = {5,1,2,4,0,3};
   int expectedcomps[10] = {0,1,2,3,4,5,6,7,8,9};
   int expectedcompbegins[4] = {0,3,6,10};
   int expectedcolbegins[3] = {0,2,3};
   int ngraphcomponents;
   int ncompcolors;
   int ntwocycleperms;
   int nusedperms;
   int i;

   SCIP_CALL( SCIPcreateProbBasic(scip, "subgroup2"));

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   perms[0] = perm1;
   perms[1] = perm2;
   perms[2] = perm3;
   perms[3] = perm4;
   perms[4] = perm5;
   perms[5] = perm6;

   SCIP_CALL( SCIPcreateVarBasic(scip, &dummyvar, "dummyvar", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   for (i = 0; i < 10; ++i)
      permvars[i] = dummyvar;

   propdata.npermvars = 10;
   propdata.nperms = 6;
   propdata.perms = perms;
   propdata.ncomponents = 1;
   propdata.components = components;
   propdata.componentbegins = componentbegins;
   propdata.componentblocked = &componentblocked;
   propdata.computedsymmetry = TRUE;
   propdata.permvars = permvars;
   propdata.preferlessrows = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &usedperms, 6) );

   /* check sorted order */
   SCIP_CALL( SCIPallocBufferArray(scip, &permorder, 6) );

   for (i = 0; i < 6; ++i)
      permorder[i] = i;

   SCIP_CALL( chooseOrderOfGenerators(scip, &propdata, 0, &permorder, &ntwocycleperms) );
   cr_assert( ntwocycleperms == 5 );

   checkIntArraysEqual(expectedpermorder, permorder, 6, "permorder");

   SCIP_CALL( buildSubgroupGraph(scip, &propdata, permorder, ntwocycleperms, 0, &graphcomponents,
         &graphcompbegins, &compcolorbegins, &ngraphcomponents, &ncompcolors, &usedperms, &nusedperms, 6, permused) );

   cr_assert( graphcomponents != NULL );
   cr_assert( graphcompbegins != NULL );
   cr_assert( compcolorbegins != NULL );
   cr_assert( nusedperms == 5, "expected 5 used permutations, but got %d\n", nusedperms );
   cr_assert( ngraphcomponents == 3, "expected 3 graph components, but got %d\n", ngraphcomponents );
   cr_assert( ncompcolors == 2, "expected 2 component colors, but got %d\n", ncompcolors );

   checkIntArraysEqual(expectedcomps, graphcomponents, 10, "components");
   checkIntArraysEqual(expectedcompbegins, graphcompbegins, ngraphcomponents, "compbegins");
   checkIntArraysEqual(expectedcolbegins, compcolorbegins, ncompcolors, "colorbegins");

   SCIPfreeBlockMemoryArray(scip, &compcolorbegins, ncompcolors + 1);
   SCIPfreeBlockMemoryArray(scip, &graphcompbegins, ngraphcomponents + 1);
   SCIPfreeBlockMemoryArray(scip, &graphcomponents, 10);
   SCIPfreeBufferArray(scip, &permorder);
   SCIPfreeBufferArray(scip, &usedperms);
   SCIP_CALL( SCIPreleaseVar(scip, &dummyvar) );
}

/* TEST 4 */
Test(test_compute_symmetry, expr1, .description = "compute symmetry for a simple example with 4 variables and 2 nonlinear constraints - before presolving")
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_VAR* w;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** permvars;
   SCIP_EXPR* xexpr;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* zexpr;
   SCIP_EXPR* wexpr;
   SCIP_EXPR** powexprs1;
   SCIP_EXPR** powexprs2;
   SCIP_EXPR* sumexpr1;
   SCIP_EXPR* sumexpr2;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int npermvars;
   int nperms;

   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs1, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs2, 2) );

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert(conshdlr != NULL);

   /* skip test if no symmetry can be computed */
   if( !SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 + x4
    *     x1^2 + x2^2               =  1
    *                   x3^2 + x4^2 =  1
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr1"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, w) );

   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &yexpr, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &zexpr, z, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &wexpr, w, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs1[0], xexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs1[1], yexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs2[0], zexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs2[1], wexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr1, 2, powexprs1, NULL, 0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr2, 2, powexprs2, NULL, 0, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e1", sumexpr1, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e2", sumexpr2, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* make sure that symmetry is computed for all variable types */
   SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 4) );
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/onlybinarysymmetry", FALSE) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL) );

   cr_assert( nperms == 3, "number of permutations: %d, but expected: 3", nperms );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   cr_assert( norbits == 1 );
   cr_assert( orbitbegins[0] == 0 );
   cr_assert( orbitbegins[1] == 4 );
   cr_assert( orbits[0] == 0 );
   cr_assert( orbits[1] == 1 );
   cr_assert( orbits[2] == 2 );
   cr_assert( orbits[3] == 3 );
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

#ifdef SCIP_DEBUG
   {
         int i;
         int j;

         for( i = 0; i < nperms; ++i )
         {
            SCIPinfoMessage(scip, NULL, "Permutation %d: (", i);
            for( j = 0; j < npermvars; ++j )
            {
               if( j == 0 )
                  SCIPinfoMessage(scip, NULL, "%d", perms[i][j]);
               else
                  SCIPinfoMessage(scip, NULL, " %d", perms[i][j]);
            }
            SCIPinfoMessage(scip, NULL, ")\n");
         }
      }
#endif

   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs2[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs2[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs1[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs1[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIPfreeBufferArray(scip, &powexprs1);
   SCIPfreeBufferArray(scip, &powexprs2);
}

/* TEST 5 */
Test(test_compute_symmetry, expr2, .description = "compute symmetry for a more complex example with 5 variables and 3 nonlinear constraints - before presolving")
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_VAR* w;
   SCIP_VAR* v;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** permvars;
   SCIP_EXPR* xexpr;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* zexpr;
   SCIP_EXPR* wexpr;
   SCIP_EXPR* vexpr;
   SCIP_EXPR** powexprs1;
   SCIP_EXPR** powexprs2;
   SCIP_EXPR* sumexpr1;
   SCIP_EXPR* sumexpr2;
   SCIP_EXPR* sumexpr3;
   SCIP_EXPR* expexpr1;
   SCIP_EXPR* expexpr2;
   SCIP_EXPR* expexpr3;
   SCIP_Real vals1[5] = {1,1,1,1,-1};
   SCIP_Real vals2[5] = {5,5,5,5,5};
   int** perms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int npermvars;
   int nperms;

   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs1, 5) );
   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs2, 5) );

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert(conshdlr != NULL);

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
   * min x1 + x2 + 2x3 + x4 + x5
   *  0 <=  exp(  x1^2 +  x2^3 +  x3^2 +  x4^2 +  x5^2 )  <=  1
   *  0 <=  exp(  x1^2 +  x2^2 +  x3^2 +  x4^2 -  x5^2 )  <=  1
   *  0 <=  exp( 5x1^2 + 5x2^2 + 5x3^2 + 5x4^2 + 5x5^2 )  <=  1
   */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr2"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &v, "v", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, v) );

   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &yexpr, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &zexpr, z, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &wexpr, w, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &vexpr, v, NULL, NULL) );

   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs1[0], xexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs1[1], yexpr, 3, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs1[2], zexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs1[3], wexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs1[4], vexpr, 2, NULL, NULL) );

   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs2[0], xexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs2[1], yexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs2[2], zexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs2[3], wexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs2[4], vexpr, 2, NULL, NULL) );

   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr1, 5, powexprs1, NULL, 0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr2, 5, powexprs2, vals1, 0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr3, 5, powexprs2, vals2, 0, NULL, NULL) );

   SCIP_CALL( SCIPcreateExprExp(scip, &expexpr1, sumexpr1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprExp(scip, &expexpr2, sumexpr2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprExp(scip, &expexpr3, sumexpr3, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e1", expexpr1, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e2", expexpr2, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e3", expexpr3, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* make sure that symmetry is computed for all variable types */
   SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 4) );
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/onlybinarysymmetry", FALSE) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL) );

   cr_assert( nperms == 1, "number of permutations: %d, but expected: 1", nperms );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   cr_assert( norbits == 1 );
   cr_assert( orbitbegins[0] == 0 );
   cr_assert( orbitbegins[1] == 2 );
   cr_assert( orbits[0] == 0 );
   cr_assert( orbits[1] == 3 );
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

#ifdef SCIP_DEBUG
   {
      int i;
      int j;

      for (i = 0; i < nperms; ++i)
      {
         SCIPinfoMessage(scip, NULL, "Permutation %d: (", i);
         for (j = 0; j < npermvars; ++j)
         {
            if ( j == 0 )
               SCIPinfoMessage(scip, NULL, "%d", perms[i][j]);
            else
               SCIPinfoMessage(scip, NULL, " %d", perms[i][j]);
         }
         SCIPinfoMessage(scip, NULL, ")\n");
      }
   }
#endif

   SCIP_CALL( SCIPreleaseExpr(scip, &expexpr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs2[4]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs2[3]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs2[2]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs2[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs2[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs1[4]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs1[3]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs1[2]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs1[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs1[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &vexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &v) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIPfreeBufferArray(scip, &powexprs1);
   SCIPfreeBufferArray(scip, &powexprs2);
}

/* TEST 6 */
Test(test_compute_symmetry, expr3, .description = "compute symmetry for a simple example with 4 variables, 2 nonlinear constraints and 1 linear constaint - before presolving")
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_VAR* w;
   SCIP_VAR* vars[4];
   SCIP_Real vals[4];
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** permvars;
   SCIP_EXPR* xexpr;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* zexpr;
   SCIP_EXPR* wexpr;
   SCIP_EXPR** powexprs1;
   SCIP_EXPR** powexprs2;
   SCIP_EXPR* sumexpr1;
   SCIP_EXPR* sumexpr2;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int npermvars;
   int nperms;

   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs1, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs2, 2) );

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert(conshdlr != NULL);

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 + x4
    *           x1^2 + x2^2               =  1
    *                        x3^2 + x4^2  =  1
    * -inf <=   x1 + 2x2 + x3 + 2x4      <=  2
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr1"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, w) );

   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &yexpr, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &zexpr, z, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &wexpr, w, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs1[0], xexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs1[1], yexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs2[0], zexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexprs2[1], wexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr1, 2, powexprs1, NULL, 0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr2, 2, powexprs2, NULL, 0, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e1", sumexpr1, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e2", sumexpr2, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = x;
   vars[1] = y;
   vars[2] = z;
   vars[3] = w;
   vals[0] = 1.0;
   vals[1] = 2.0;
   vals[2] = 1.0;
   vals[3] = 2.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "i1", 4, vars, vals, -SCIPinfinity(scip), 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* make sure that symmetry is computed for all variable types */
   SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 4) );
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/onlybinarysymmetry", FALSE) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL) );

   cr_assert( nperms == 1, "number of permutations: %d, but expected: 1", nperms );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   cr_assert( norbits == 2 );
   cr_assert( orbitbegins[0] == 0 );
   cr_assert( orbitbegins[1] == 2 );
   cr_assert( orbitbegins[2] == 4 );
   cr_assert( orbits[0] == 0 );
   cr_assert( orbits[1] == 2 );
   cr_assert( orbits[2] == 1 );
   cr_assert( orbits[3] == 3 );
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   #ifdef SCIP_DEBUG
   {
         int i;
         int j;

         for (i = 0; i < nperms; ++i)
         {
            SCIPinfoMessage(scip, NULL, "Permutation %d: (", i);
            for (j = 0; j < npermvars; ++j)
            {
               if ( j == 0 )
                  SCIPinfoMessage(scip, NULL, "%d", perms[i][j]);
               else
                  SCIPinfoMessage(scip, NULL, " %d", perms[i][j]);
            }
            SCIPinfoMessage(scip, NULL, ")\n");
         }
      }
   #endif

   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs2[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs2[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs1[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexprs1[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIPfreeBufferArray(scip, &powexprs1);
   SCIPfreeBufferArray(scip, &powexprs2);
}

/* TEST 6 */
Test(test_compute_symmetry, expr4, .description = "compute symmetry for a simple example with 4 variables and 3 nonlinear constraints - before presolving")
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_VAR* w;
   SCIP_Real vals[2];
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** permvars;
   SCIP_EXPR** summands;
   SCIP_EXPR* xexpr;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* zexpr;
   SCIP_EXPR* wexpr;
   SCIP_EXPR* powexpr1;
   SCIP_EXPR* powexpr2;
   SCIP_EXPR* prodexpr;
   SCIP_EXPR* sumexpr1;
   SCIP_EXPR* sumexpr2;
   SCIP_EXPR* sumexpr3;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int npermvars;
   int nperms;
   SCIP_Bool binvarsaffected;

   SCIP_CALL( SCIPallocBufferArray(scip, &summands, 3) );

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert(conshdlr != NULL);

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2
    *           x1^2 - 4x3                <=  0
    *                        x2^2 - 4x4   <=  0
    * x1*x2 + x1 + x2                     <=  1
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr1"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 2.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 2.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, w) );

   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &yexpr, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &zexpr, z, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &wexpr, w, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr1, xexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr2, yexpr, 2, NULL, NULL) );
   summands[0] = xexpr;
   summands[1] = yexpr;
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr, 2, summands, 1, NULL, NULL) );

   summands[0] = powexpr1;
   summands[1] = zexpr;
   vals[0] = 1.0;
   vals[1] = -4.0;
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr1, 2, summands, vals, 0, NULL, NULL) );

   summands[0] = powexpr2;
   summands[1] = wexpr;
   vals[0] = 1.0;
   vals[1] = -4.0;
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr2, 2, summands, vals, 0, NULL, NULL) );

   summands[0] = xexpr;
   summands[1] = yexpr;
   summands[2] = prodexpr;
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr3, 3, summands, NULL, 0, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e1", sumexpr1, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e2", sumexpr2, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "e3", sumexpr3, -SCIPinfinity(scip), 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

   /* make sure that symmetry is computed for all variable types */
   SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 4) );
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/onlybinarysymmetry", FALSE) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL, NULL, &binvarsaffected,
         NULL, NULL, NULL, NULL) );

   cr_assert(!binvarsaffected);
   cr_assert( nperms == 1, "number of permutations: %d, but expected: 1", nperms );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, permvars, npermvars, perms, nperms, orbits, orbitbegins, &norbits) );
   cr_assert( norbits == 2 );
   cr_assert( orbitbegins[0] == 0 );
   cr_assert( orbitbegins[1] == 2 );
   cr_assert( orbitbegins[2] == 4 );
   cr_assert( orbits[0] == 0 );
   cr_assert( orbits[1] == 1 );
   cr_assert( orbits[2] == 2 );
   cr_assert( orbits[3] == 3 );
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

#ifdef SCIP_DEBUG
   {
      int i;
      int j;

      for (i = 0; i < nperms; ++i)
      {
         SCIPinfoMessage(scip, NULL, "Permutation %d: (", i);
         for (j = 0; j < npermvars; ++j)
         {
            if ( j == 0 )
               SCIPinfoMessage(scip, NULL, "%d", perms[i][j]);
            else
               SCIPinfoMessage(scip, NULL, " %d", perms[i][j]);
         }
         SCIPinfoMessage(scip, NULL, ")\n");
      }
   }
#endif

   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIPfreeBufferArray(scip, &summands);
}
