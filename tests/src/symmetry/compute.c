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

/**@file   compute.c
 * @brief  unit tests for computing symmetry
 * @author Marc Pfetsch
 * @author Fabian Wegscheider
 * @author Christopher Hojny
 */

#include <scip/scip.h>
#include <include/scip_test.h>
#include <scip/scip_sym.h>
#include <scip/symmetry.h>
#include <scip/symmetry_graph.h>
#include <symmetry/compute_symmetry.h>
#include <scip/scipdefplugins.h>

/* global SCIP instance */
static SCIP* scip;



/** sort orbits to compare */
static
void sortOrbits(
   int                   norbits,            /**< number of orbits */
   int*                  orbits,             /**< array that stores the orbits */
   int*                  orbitbegins         /**< array that marks the beginning of the orbits */
   )
{
   int i;

   for (i = 0; i < norbits; ++i)
   {
      SCIPsortInt(&orbits[orbitbegins[i]], orbitbegins[i+1] - orbitbegins[i]);
   }

}

/** setup: create SCIP */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* turn on symmetry computation */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

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
   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "Memory leak!");
}

/** simple example with 4 variables and 2 linear constraints */
static
void simpleExample1(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
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
   int permlen;

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

   /* determine symmetry type */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 8;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 4;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( componentbegins[0] == 0 );
   TEST_ASSERT( componentbegins[1] == nperms );
   TEST_ASSERT( vartocomponent[0] == 0 );
   TEST_ASSERT( vartocomponent[1] == 0 );
   TEST_ASSERT( vartocomponent[2] == 0 );
   TEST_ASSERT( vartocomponent[3] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 4 );
      TEST_ASSERT( orbitbegins[2] == 8 );
      TEST_ASSERT( orbits[0] == 0 );
      TEST_ASSERT( orbits[1] == 1 );
      TEST_ASSERT( orbits[2] == 2 );
      TEST_ASSERT( orbits[3] == 3 );
      TEST_ASSERT( orbits[4] == 4 );
      TEST_ASSERT( orbits[5] == 5 );
      TEST_ASSERT( orbits[6] == 6 );
      TEST_ASSERT( orbits[7] == 7 );
   }
   else
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 4 );
      TEST_ASSERT( orbits[0] == 0 );
      TEST_ASSERT( orbits[1] == 1 );
      TEST_ASSERT( orbits[2] == 2 );
      TEST_ASSERT( orbits[3] == 3 );
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}


/** simple example with 4 variables and 4 linear constraints */
static
void simpleExample2(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
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
   int permlen;

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

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 8;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 4;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   TEST_ASSERT( nperms == 1 );
   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( componentbegins[0] == 0 );
   TEST_ASSERT( componentbegins[1] == 1 );
   TEST_ASSERT( vartocomponent[0] == 0 );
   TEST_ASSERT( vartocomponent[1] == 0 );
   TEST_ASSERT( vartocomponent[2] == 0 );
   TEST_ASSERT( vartocomponent[3] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 4 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
      TEST_ASSERT( orbitbegins[3] == 6 );
      TEST_ASSERT( orbits[0] == 0 );
      TEST_ASSERT( orbits[1] == 1 );
      TEST_ASSERT( orbits[2] == 2 );
      TEST_ASSERT( orbits[3] == 3 );
      TEST_ASSERT( orbits[4] == 4 );
      TEST_ASSERT( orbits[5] == 5 );
      TEST_ASSERT( orbits[6] == 6 );
      TEST_ASSERT( orbits[7] == 7 );
   }
   else
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbits[0] == 0 );
      TEST_ASSERT( orbits[1] == 1 );
      TEST_ASSERT( orbits[2] == 2 );
      TEST_ASSERT( orbits[3] == 3 );
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}

/** simple example with 4 variables and 4 linear constraints */
static
void simpleExample3(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 + x4 + x5
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

   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 10;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 5;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );
   TEST_ASSERT( nperms == 2 );
   TEST_ASSERT( ncomponents == 2 );
   TEST_ASSERT( vartocomponent[0] == vartocomponent[1] );
   TEST_ASSERT( vartocomponent[2] == vartocomponent[3] );
   TEST_ASSERT( vartocomponent[0] != vartocomponent[2] );
   TEST_ASSERT( vartocomponent[1] != vartocomponent[3] );
   TEST_ASSERT( vartocomponent[4] == -1 );
   TEST_ASSERT( componentbegins[0] == 0 );
   TEST_ASSERT( componentbegins[1] == 1 );
   TEST_ASSERT( componentbegins[2] == 2 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 4 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
      TEST_ASSERT( orbitbegins[3] == 6 );
      TEST_ASSERT( orbitbegins[4] == 8 );
   }
   else
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
   SCIP_CALL( SCIPreleaseVar(scip, &var5) );
}

/** simple example with 6 variables and 3 bounddisjunction constraints */
static
void exampleBounddisjunction(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_VAR* var5;
   SCIP_VAR* var6;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_BOUNDTYPE btypes[2];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 - x2 + x3 - x4 + x5 - x6
    *     BD(x1 <= -1, x2 >= 1)
    *     BD(x3 <= 7, x4 >= 9)
    *     BD(x5 <= -1, x6 >= 1)
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "BD"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -10, 10, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -10, 10, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -10, 10, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -10, 10, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var5, "x5", -10, 10, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var5) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var6, "x6", -10, 10, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var6) );

   vars[0] = var1;
   vars[1] = var2;
   vals[0] = -1.0;
   vals[1] = 1.0;
   btypes[0] = SCIP_BOUNDTYPE_UPPER;
   btypes[1] = SCIP_BOUNDTYPE_LOWER;
   SCIP_CALL( SCIPcreateConsBasicBounddisjunction(scip, &cons, "c1", 2, vars, btypes, vals) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var3;
   vars[1] = var4;
   vals[0] = 7.0;
   vals[1] = 9.0;
   SCIP_CALL( SCIPcreateConsBasicBounddisjunction(scip, &cons, "c2", 2, vars, btypes, vals) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var5;
   vars[1] = var6;
   vals[0] = -1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicBounddisjunction(scip, &cons, "c3", 2, vars, btypes, vals) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 12;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 6;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   if ( detectsignedperms )
   {
      TEST_ASSERT( nperms == 2 || nperms == 3 );
      TEST_ASSERT( ncomponents == 1 );
      TEST_ASSERT( vartocomponent[0] == vartocomponent[1] );
      TEST_ASSERT( vartocomponent[1] == vartocomponent[4] );
      TEST_ASSERT( vartocomponent[4] == vartocomponent[5] );
      TEST_ASSERT( vartocomponent[2] == -1 );
      TEST_ASSERT( vartocomponent[3] == -1 );
      TEST_ASSERT( componentbegins[0] == 0 );
      TEST_ASSERT( componentbegins[1] == nperms );
   }
   else
   {
      TEST_ASSERT( nperms == 1 );
      TEST_ASSERT( ncomponents == 1 );
      TEST_ASSERT( vartocomponent[0] == vartocomponent[1] );
      TEST_ASSERT( vartocomponent[1] == vartocomponent[4] );
      TEST_ASSERT( vartocomponent[4] == vartocomponent[5] );
      TEST_ASSERT( vartocomponent[2] == -1 );
      TEST_ASSERT( vartocomponent[3] == -1 );
      TEST_ASSERT( componentbegins[0] == 0 );
      TEST_ASSERT( componentbegins[1] == 1 );
   }

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 4 );
      TEST_ASSERT( orbitbegins[2] == 8 );
   }
   else
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
   SCIP_CALL( SCIPreleaseVar(scip, &var5) );
   SCIP_CALL( SCIPreleaseVar(scip, &var6) );
}

/** simple example with 4 variables and a cardinality constraint */
static
void exampleCardinality(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_VAR* ind1;
   SCIP_VAR* ind2;
   SCIP_VAR* ind3;
   SCIP_VAR* ind4;
   SCIP_CONS* cons;
   SCIP_VAR* vars[4];
   SCIP_VAR* inds[4];
   SCIP_Real vals[4];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 - x2 + x3 + x4
    *     x1 - x2 + x3 + x4 >= 2
    *     CARD(x1, x2, x3, x4) <= 3
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "Card"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ind1, "ind1", 0, 1, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, ind1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ind2, "ind2", 0, 1, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, ind2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ind3, "ind3", 0, 1, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, ind3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ind4, "ind4", 0, 1, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, ind4) );

   vars[0] = var1;
   vars[1] = var2;
   vars[2] = var3;
   vars[3] = var4;
   vals[0] = 1.0;
   vals[1] = -1.0;
   vals[2] = 1.0;
   vals[3] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 4, vars, vals, 2.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   inds[0] = ind1;
   inds[1] = ind2;
   inds[2] = ind3;
   inds[3] = ind4;
   SCIP_CALL( SCIPcreateConsBasicCardinality(scip, &cons, "c2", 4, vars, 3, inds, NULL) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 16;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 8;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   if ( detectsignedperms )
   {
      TEST_ASSERT( nperms == 3 || nperms == 6 ); /* if more involutions are generated from existing one, then 6 */
      TEST_ASSERT( ncomponents == 1 );
      TEST_ASSERT( vartocomponent[0] == 0 );
      TEST_ASSERT( vartocomponent[1] == 0 );
      TEST_ASSERT( vartocomponent[2] == 0 );
      TEST_ASSERT( vartocomponent[3] == 0 );
      TEST_ASSERT( vartocomponent[4] == 0 );
      TEST_ASSERT( vartocomponent[5] == 0 );
      TEST_ASSERT( vartocomponent[6] == 0 );
      TEST_ASSERT( vartocomponent[7] == 0 );
   }
   else
   {
      TEST_ASSERT( nperms == 2 );
      TEST_ASSERT( ncomponents == 1 );
      TEST_ASSERT( vartocomponent[0] == vartocomponent[2] );
      TEST_ASSERT( vartocomponent[2] == vartocomponent[3] );
      TEST_ASSERT( vartocomponent[1] == -1 );
      TEST_ASSERT( componentbegins[0] == 0 );
      TEST_ASSERT( componentbegins[1] == 2 );
   }

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 4 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 4 );
      TEST_ASSERT( orbitbegins[2] == 8 );
      TEST_ASSERT( orbitbegins[3] == 12 );
      TEST_ASSERT( orbitbegins[4] == 16 );
   }
   else
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 3 );
      TEST_ASSERT( orbitbegins[2] == 6 );
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
   SCIP_CALL( SCIPreleaseVar(scip, &ind1) );
   SCIP_CALL( SCIPreleaseVar(scip, &ind2) );
   SCIP_CALL( SCIPreleaseVar(scip, &ind3) );
   SCIP_CALL( SCIPreleaseVar(scip, &ind4) );
}

/** simple example with 6 variables and indicator constraints */
static
void exampleIndicator(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_VAR* bin1;
   SCIP_VAR* bin2;
   SCIP_CONS* cons;
   SCIP_VAR* vars[4];
   SCIP_Real vals[4];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 - x2 + x3 - x4
    *     b1 = 1 --> x1 - x2 <= 2
    *     b2 = 1 --> x3 - x4 <= 2
    *     x1 - x2 + x3 - x4 >= 0
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "Indicator"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -SCIPinfinity(scip), SCIPinfinity(scip), -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &bin1, "bin1", 0, 1, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, bin1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &bin2, "bin2", 0, 1, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, bin2) );

   vars[0] = var1;
   vars[1] = var2;
   vars[2] = var3;
   vars[3] = var4;
   vals[0] = 1.0;
   vals[1] = -1.0;
   vals[2] = 1.0;
   vals[3] = -1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 4, vars, vals, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicIndicator(scip, &cons, "c2", bin1, 2, vars, vals, 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var3;
   vars[1] = var4;
   SCIP_CALL( SCIPcreateConsBasicIndicator(scip, &cons, "c3", bin2, 2, vars, vals, 2.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 16;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 8;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( vartocomponent[0] == 0 );
   TEST_ASSERT( vartocomponent[1] == 0 );
   TEST_ASSERT( vartocomponent[2] == 0 );
   TEST_ASSERT( vartocomponent[3] == 0 );
   TEST_ASSERT( vartocomponent[4] == 0 );
   TEST_ASSERT( vartocomponent[5] == 0 );
   TEST_ASSERT( vartocomponent[6] == 0 );
   TEST_ASSERT( vartocomponent[7] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      int orbitlens[6];
      int i;

      TEST_ASSERT( norbits == 6 );

      for (i = 0; i < 6; ++i)
         orbitlens[i] = orbitbegins[i+1] - orbitbegins[i];
      SCIPsortInt(orbitlens, 6);

      TEST_ASSERT( orbitlens[0] == 2 );
      TEST_ASSERT( orbitlens[1] == 2 );
      TEST_ASSERT( orbitlens[2] == 2 );
      TEST_ASSERT( orbitlens[3] == 2 );
      TEST_ASSERT( orbitlens[4] == 4 );
      TEST_ASSERT( orbitlens[5] == 4 );
   }
   else
   {
      TEST_ASSERT( norbits == 4 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
      TEST_ASSERT( orbitbegins[3] == 6 );
      TEST_ASSERT( orbitbegins[4] == 8 );
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
   SCIP_CALL( SCIPreleaseVar(scip, &bin1) );
   SCIP_CALL( SCIPreleaseVar(scip, &bin2) );
}

/** simple example with 4 variables and SOS1 constraints */
static
void exampleSOS1(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_VAR* vars[4];
   SCIP_Real vals[4];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 - x2 + x3 - x4
    *     SOS1(x1,x2)
    *     SOS1(X3,x4)
    *     x1 - x2 + x3 - x4 >= 0
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "SOS1"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -SCIPinfinity(scip), SCIPinfinity(scip), -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );

   vars[0] = var1;
   vars[1] = var2;
   vars[2] = var3;
   vars[3] = var4;
   vals[0] = 1.0;
   vals[1] = -1.0;
   vals[2] = 1.0;
   vals[3] = -1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 4, vars, vals, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &cons, "c2", 2, vars, NULL) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var3;
   vars[1] = var4;
   SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &cons, "c3", 2, vars, NULL) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 8;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 4;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   if( detectsignedperms )
   {
      TEST_ASSERT( nperms == 2 || nperms == 3 );
   }
   else
   {
      TEST_ASSERT( nperms == 1 );
   }
   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( vartocomponent[0] == 0 );
   TEST_ASSERT( vartocomponent[1] == 0 );
   TEST_ASSERT( vartocomponent[2] == 0 );
   TEST_ASSERT( vartocomponent[3] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 4 );
      TEST_ASSERT( orbitbegins[2] == 8 );
   }
   else
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}

/** simple example with 6 variables and SOS2 constraints */
static
void exampleSOS2(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_VAR* var5;
   SCIP_VAR* var6;
   SCIP_CONS* cons;
   SCIP_VAR* vars[6];
   SCIP_Real vals[6];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3 - x4 - x5 - x6
    *     SOS2(x1,x2,x3)
    *     SOS2(X4,x5,x6)
    *     x1 + x2 + x3 - x4 - x5 - x6 >= 0
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "SOS2"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -SCIPinfinity(scip), SCIPinfinity(scip), -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var5, "x5", -SCIPinfinity(scip), SCIPinfinity(scip), -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var5) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var6, "x6", -SCIPinfinity(scip), SCIPinfinity(scip), -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var6) );

   vars[0] = var1;
   vars[1] = var2;
   vars[2] = var3;
   vars[3] = var4;
   vars[4] = var5;
   vars[5] = var6;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;
   vals[3] = -1.0;
   vals[4] = -1.0;
   vals[5] = -1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 6, vars, vals, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicSOS2(scip, &cons, "c2", 3, vars, NULL) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   vars[0] = var4;
   vars[1] = var5;
   vars[2] = var6;
   SCIP_CALL( SCIPcreateConsBasicSOS2(scip, &cons, "c3", 3, vars, NULL) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 12;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 6;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   if( detectsignedperms )
   {
      TEST_ASSERT( ncomponents == 1 );
      TEST_ASSERT( vartocomponent[0] == 0 );
      TEST_ASSERT( vartocomponent[1] == 0 );
      TEST_ASSERT( vartocomponent[2] == 0 );
      TEST_ASSERT( vartocomponent[3] == 0 );
      TEST_ASSERT( vartocomponent[4] == 0 );
      TEST_ASSERT( vartocomponent[5] == 0 );
   }
   else
   {
      TEST_ASSERT( ncomponents == 2 );
      TEST_ASSERT( vartocomponent[0] == 0 );
      TEST_ASSERT( vartocomponent[1] == -1 );
      TEST_ASSERT( vartocomponent[2] == 0 );
      TEST_ASSERT( vartocomponent[3] == 1 );
      TEST_ASSERT( vartocomponent[4] == -1 );
      TEST_ASSERT( vartocomponent[5] == 1 );
   }

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      int orbitlens[4];
      int i;

      TEST_ASSERT( norbits == 4 );

      for (i = 0; i < 4; ++i)
         orbitlens[i] = orbitbegins[i+1] - orbitbegins[i];
      SCIPsortInt(orbitlens, 4);

      TEST_ASSERT( orbitlens[0] == 2 );
      TEST_ASSERT( orbitlens[1] == 2 );
      TEST_ASSERT( orbitlens[2] == 4 );
      TEST_ASSERT( orbitlens[3] == 4 );
   }
   else
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
   SCIP_CALL( SCIPreleaseVar(scip, &var5) );
   SCIP_CALL( SCIPreleaseVar(scip, &var6) );
}

/** simple example with 3 variables and a pseudoboolean constraint */
static
void examplePB(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR** terms[3];
   SCIP_VAR* term1[2];
   SCIP_VAR* term2[2];
   SCIP_VAR* term3[2];
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_Real vals[3];
   SCIP_CONS* cons;
   SCIP_VAR** permvars;
   int** perms;
   int nterms[3];
   int* components;
   int* componentbegins;
   int* vartocomponent;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int permlen;
   int ncomponents;
   int npermvars;
   int nperms;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    * min x1 + x2 + x3
    *     -2.0 <= x1 x2 + x1 x3 - x2 x3 <= 2.0
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "PB"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var3) );

   term1[0] = var1;
   term1[1] = var2;
   term2[0] = var1;
   term2[1] = var3;
   term3[0] = var2;
   term3[1] = var3;
   terms[0] = term1;
   terms[1] = term2;
   terms[2] = term3;
   nterms[0] = 2;
   nterms[1] = 2;
   nterms[2] = 2;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = -1.0;

   SCIP_CALL( SCIPcreateConsPseudoboolean(scip, &cons, "c1", NULL, 0, NULL, terms, 3, nterms, vals, NULL, 0.0, FALSE,
         -2.0, 2.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 12;             /* 3 binary variables + 3 artificial variables for AND-resultants + reflection */
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 6;              /* 3 binary variables + 3 artificial variables for AND-resultants */
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   TEST_ASSERT( nperms == 1 );
   TEST_ASSERT( ncomponents == 1 );
   SCIPsortInt(vartocomponent, 6);
   TEST_ASSERT( vartocomponent[0] == -1 );
   TEST_ASSERT( vartocomponent[2] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 4 );
      TEST_ASSERT( orbitbegins[1] - orbitbegins[0] == 2 );
      TEST_ASSERT( orbitbegins[2] - orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[3] - orbitbegins[2] == 2 );
      TEST_ASSERT( orbitbegins[4] - orbitbegins[3] == 2 );
   }
   else
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[1] - orbitbegins[0] == 2 );
      TEST_ASSERT( orbitbegins[2] - orbitbegins[1] == 2 );
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
}

/** simple example with 3 variables and nonlinear constraints */
static
void exampleExpr1(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPR* varexpr1;
   SCIP_EXPR* varexpr2;
   SCIP_EXPR* varexpr3;
   SCIP_EXPR* powexpr;
   SCIP_EXPR* prodexpr;
   SCIP_EXPR* exprs[3];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   TEST_ASSERT(conshdlr != NULL);

   /* setup problem:
    * min x1 + x2 + x3
    *     x1 + x2 + x3   >= 2
    *     x1^3 * x2 * x3 == 0
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr1"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );

   /* create linear constraint */
   vars[0] = var1;
   vars[1] = var2;
   vars[2] = var3;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 3, vars, vals, 2.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create nonlinear constraint */
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr1, var1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr2, var2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr3, var3, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr, varexpr1, 3, NULL, NULL) );
   exprs[0] = powexpr;
   exprs[1] = varexpr2;
   exprs[2] = varexpr3;
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr, 3, exprs, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c2", prodexpr, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 6;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 3;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   TEST_ASSERT( nperms == 1 );
   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( vartocomponent[0] == -1 );
   TEST_ASSERT( vartocomponent[1] == 0 );
   TEST_ASSERT( vartocomponent[2] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
   }
   else
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
}

/** simple example with 5 variables and nonlinear constraints */
static
void exampleExpr2(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_VAR* var5;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPR* varexpr1;
   SCIP_EXPR* varexpr2;
   SCIP_EXPR* varexpr3;
   SCIP_EXPR* varexpr4;
   SCIP_EXPR* varexpr5;
   SCIP_EXPR* powexpr1;
   SCIP_EXPR* powexpr2;
   SCIP_EXPR* prodexpr1;
   SCIP_EXPR* prodexpr2;
   SCIP_EXPR* exprs[3];
   SCIP_VAR* vars[5];
   SCIP_Real vals[5];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   TEST_ASSERT(conshdlr != NULL);

   /* setup problem:
    * min x1 + x2 + x3 + x4 + x5
    *     x1 + x2 + x3 + x4 + x5  >= 2
    *     x1^3 * x2 * x3          == 0
    *     x4^3 * x2 * x5          == 0
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr2"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var5, "x5", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var5) );

   /* create linear constraint */
   vars[0] = var1;
   vars[1] = var2;
   vars[2] = var3;
   vars[3] = var4;
   vars[4] = var5;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;
   vals[3] = 1.0;
   vals[4] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 5, vars, vals, 2.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create nonlinear constraints */
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr1, var1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr2, var2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr3, var3, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr4, var4, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr5, var5, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr1, varexpr1, 3, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr2, varexpr4, 3, NULL, NULL) );
   exprs[0] = powexpr1;
   exprs[1] = varexpr2;
   exprs[2] = varexpr3;
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr1, 3, exprs, 1.0, NULL, NULL) );
   exprs[0] = powexpr2;
   exprs[1] = varexpr2;
   exprs[2] = varexpr5;
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr2, 3, exprs, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c2", prodexpr1, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c3", prodexpr2, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 10;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 5;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   TEST_ASSERT( nperms == 1 );
   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( vartocomponent[0] == 0 );
   TEST_ASSERT( vartocomponent[1] == -1 );
   TEST_ASSERT( vartocomponent[2] == 0 );
   TEST_ASSERT( vartocomponent[3] == 0 );
   TEST_ASSERT( vartocomponent[4] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 4 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
      TEST_ASSERT( orbitbegins[3] == 6 );
      TEST_ASSERT( orbitbegins[4] == 8 );
   }
   else
   {
      TEST_ASSERT( norbits == 2 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
      TEST_ASSERT( orbitbegins[2] == 4 );
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr5) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr4) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
   SCIP_CALL( SCIPreleaseVar(scip, &var5) );
}

/** simple example with 4 variables and nonlinear constraints */
static
void exampleExpr3(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPR* varexpr1;
   SCIP_EXPR* varexpr2;
   SCIP_EXPR* powexpr1;
   SCIP_EXPR* powexpr2;
   SCIP_EXPR* sumexpr;
   SCIP_EXPR* exprs[2];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   TEST_ASSERT(conshdlr != NULL);

   /* setup problem:
    * min x3 + 2*x4
    *     x3 + x4     >= 2
    *     x1^2 + x2^2 == 1
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr3"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -SCIPinfinity(scip), SCIPinfinity(scip), 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );

   /* create linear constraint */
   vars[0] = var3;
   vars[1] = var4;
   vals[0] = 1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 2, vars, vals, 2.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create nonlinear constraints */
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr1, var1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr2, var2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr1, varexpr1, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr2, varexpr2, 2, NULL, NULL) );
   exprs[0] = powexpr1;
   exprs[1] = powexpr2;
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 2, exprs, vals, 0.0, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c2", sumexpr, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 8;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 4;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   if ( detectsignedperms )
   {
      TEST_ASSERT( nperms == 2 || nperms == 3 );
   }
   else
   {
      TEST_ASSERT( nperms == 1 );
   }
   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( vartocomponent[0] == 0 );
   TEST_ASSERT( vartocomponent[1] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 4 );
   }
   else
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}

/** simple example with 2 variables and nonlinear constraints */
static
void exampleExpr4(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPR* varexpr1;
   SCIP_EXPR* varexpr2;
   SCIP_EXPR* powexpr1;
   SCIP_EXPR* powexpr2;
   SCIP_EXPR* sumexpr;
   SCIP_EXPR* prodexpr;
   SCIP_EXPR* exprs[2];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   TEST_ASSERT(conshdlr != NULL);

   /* setup problem:
    * min x3 + 2*x4
    *     x3 + x4     >= 2
    *     x1^2 + x2^2 == 1
    *     x1 * x2     == 0
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr4"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -SCIPinfinity(scip), SCIPinfinity(scip), 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );

   /* create linear constraint */
   vars[0] = var3;
   vars[1] = var4;
   vals[0] = 1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 2, vars, vals, 2.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create nonlinear constraints */
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr1, var1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr2, var2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr1, varexpr1, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr2, varexpr2, 2, NULL, NULL) );
   exprs[0] = powexpr1;
   exprs[1] = powexpr2;
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 2, exprs, vals, 0.0, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c2", sumexpr, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   exprs[0] = varexpr1;
   exprs[1] = varexpr2;
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr, 2, exprs, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c3", prodexpr, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 8;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 4;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   if ( detectsignedperms )
   {
      TEST_ASSERT( nperms == 2 );
   }
   else
   {
      TEST_ASSERT( nperms == 1 );
   }
   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( vartocomponent[0] == 0 );
   TEST_ASSERT( vartocomponent[1] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 4 );
   }
   else
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}

/** simple example with 4 variables and nonlinear constraints */
static
void exampleExpr5(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPR* varexpr1;
   SCIP_EXPR* varexpr2;
   SCIP_EXPR* cosexpr1;
   SCIP_EXPR* cosexpr2;
   SCIP_EXPR* sumexpr;
   SCIP_EXPR* exprs[2];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   TEST_ASSERT(conshdlr != NULL);

   /* setup problem:
    * min x3 + 2*x4
    *     x3 + x4           >= 2
    *     cos(x1) + cos(x2) == 1
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr5"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -SCIPinfinity(scip), SCIPinfinity(scip), 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );

   /* create linear constraint */
   vars[0] = var3;
   vars[1] = var4;
   vals[0] = 1.0;
   vals[1] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 2, vars, vals, 2.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create nonlinear constraints */
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr1, var1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr2, var2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprCos(scip, &cosexpr1, varexpr1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprCos(scip, &cosexpr2, varexpr2, NULL, NULL) );
   exprs[0] = cosexpr1;
   exprs[1] = cosexpr2;
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 2, exprs, vals, 0.0, NULL, NULL) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c2", sumexpr, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 8;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 4;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( vartocomponent[0] == 0 );
   TEST_ASSERT( vartocomponent[1] == 0 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 4 );
   }
   else
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 2 );
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &cosexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &cosexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}

/** simple example with 5 variables and nonlinear constraints */
static
void exampleExpr6(
   SCIP_Bool             detectsignedperms   /**< whether signed permutations shall be detected */
   )
{
   SYM_SYMTYPE symtype;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_VAR* var5;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPR* varexpr1;
   SCIP_EXPR* varexpr2;
   SCIP_EXPR* varexpr3;
   SCIP_EXPR* varexpr4;
   SCIP_EXPR* varexpr5;
   SCIP_EXPR* powexpr1;
   SCIP_EXPR* powexpr2;
   SCIP_EXPR* powexpr3;
   SCIP_EXPR* powexpr4;
   SCIP_EXPR* prodexpr1;
   SCIP_EXPR* prodexpr2;
   SCIP_EXPR* sumexpr1;
   SCIP_EXPR* sumexpr2;
   SCIP_EXPR* exprs[4];
   SCIP_Real vals[4];
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
   int permlen;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   TEST_ASSERT(conshdlr != NULL);

   /* setup problem:
    * min -x5
    *     x1^2 -2 * x1 * x2 + x2^2 >= x5
    *     x3^2 -2 * x3 * x4 + x4^2 >= x5
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "expr6"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", -1, 1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", -1, 1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", -1, 1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", -1, 1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var4) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var5, "x5", 0, 4, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var5) );

   /* create nonlinear constraints */
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr1, var1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr2, var2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr3, var3, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr4, var4, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr5, var5, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr1, varexpr1, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr2, varexpr2, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr3, varexpr3, 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &powexpr4, varexpr4, 2, NULL, NULL) );

   exprs[0] = varexpr1;
   exprs[1] = varexpr2;
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr1, 2, exprs, 1.0, NULL, NULL) );

   exprs[0] = varexpr3;
   exprs[1] = varexpr4;
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr2, 2, exprs, 1.0, NULL, NULL) );

   exprs[0] = powexpr1;
   exprs[1] = powexpr2;
   exprs[2] = prodexpr1;
   exprs[3] = varexpr5;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = -2.0;
   vals[3] = -1.0;
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr1, 4, exprs, vals, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c1", sumexpr1, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   exprs[0] = powexpr3;
   exprs[1] = powexpr4;
   exprs[2] = prodexpr2;
   exprs[3] = varexpr5;
   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = -2.0;
   vals[3] = -1.0;
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr2, 4, exprs, vals, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c2", sumexpr2, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* general symmetry detection */
   SCIP_CALL( SCIPsetBoolParam(scip, "symmetries/enabled", TRUE) );

   /* note that indicator constraints introduce a slack variable, i.e., we have 8 variables */
   if ( detectsignedperms )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 1) );
      permlen = 10;
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(scip, "symmetries/symtype", 0) );
      permlen = 5;
   }

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetry(scip, &symtype,
         &npermvars, &permvars, NULL, &nperms, &perms, NULL,
         &components, &componentbegins, &vartocomponent, &ncomponents) );

   TEST_ASSERT( ncomponents == 1 );
   TEST_ASSERT( vartocomponent[0] == 0 );
   TEST_ASSERT( vartocomponent[1] == 0 );
   TEST_ASSERT( vartocomponent[2] == 0 );
   TEST_ASSERT( vartocomponent[3] == 0 );
   TEST_ASSERT( vartocomponent[4] == -1 );

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, permlen) );
   SCIP_CALL( SCIPcomputeOrbitsSym(scip, detectsignedperms, permvars, npermvars,
         perms, nperms, orbits, orbitbegins, &norbits) );

   /* make sure orbits are sorted for comparison below */
   sortOrbits(norbits, orbits, orbitbegins);

   if ( detectsignedperms )
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 8 );
   }
   else
   {
      TEST_ASSERT( norbits == 1 );
      TEST_ASSERT( orbitbegins[0] == 0 );
      TEST_ASSERT( orbitbegins[1] == 4 );
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr4) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &powexpr1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr5) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr4) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
   SCIP_CALL( SCIPreleaseVar(scip, &var5) );
}

/* TEST SUITE */

/* TEST 1 */
/** @brief compute permutation symmetries for a simple example with 4 variables and 2 linear constraints */
void test_test_compute_symmetry_basic1(void)
{
   simpleExample1(FALSE);
}

/* TEST 2 */
/** @brief compute signed symmetries for a simple example with 4 variables and 2 linear constraints */
void test_test_compute_symmetry_basic2(void)
{
   simpleExample1(TRUE);
}

/* TEST 3 */
/** @brief compute permutation symmetry for a simple example with 4 variables and 4 linear constraints */
void test_test_compute_symmetry_basic3(void)
{
   simpleExample2(FALSE);
}

/* TEST 4 */
/** @brief compute signed permutation symmetry for a simple example with 4 variables and 4 linear constraints */
void test_test_compute_symmetry_basic4(void)
{
   simpleExample2(TRUE);
}

/* TEST 5 */
/** @brief compute permutation symmetries for a simple example with 5 variables and 2 linear constraints */
void test_test_compute_symmetry_basic5(void)
{
   simpleExample3(FALSE);
}

/* TEST 6 */
/** @brief compute signed permutation symmetries for a simple example with 5 variables and 2 linear constraints */
void test_test_compute_symmetry_basic6(void)
{
   simpleExample3(TRUE);
}

/* TEST 7 */
/** @brief compute permutation symmetries for an example containing bounddisjunction constraints */
void test_test_compute_symmetry_special1(void)
{
   exampleBounddisjunction(FALSE);
}

/* TEST 8 */
/** @brief compute signed permutation symmetries for an example containing bounddisjunction constraints */
void test_test_compute_symmetry_special2(void)
{
   exampleBounddisjunction(TRUE);
}

/* TEST 9 */
/** @brief compute permutation symmetries for an example containing cardinality constraints */
void test_test_compute_symmetry_special3(void)
{
   exampleCardinality(FALSE);
}

/* TEST 10 */
/** @brief compute signed permutation symmetries for an example containing cardinality constraints */
void test_test_compute_symmetry_special4(void)
{
   exampleCardinality(TRUE);
}

/* TEST 11 */
/** @brief compute permutation symmetries for an example containing indicator constraints */
void test_test_compute_symmetry_special5(void)
{
   exampleIndicator(FALSE);
}

/* TEST 12 */
/** @brief compute signed permutation symmetries for an example containing indicator constraints */
void test_test_compute_symmetry_special6(void)
{
   exampleIndicator(TRUE);
}

/* TEST 13 */
/** @brief compute permutation symmetries for an example containing SOS1 constraints */
void test_test_compute_symmetry_special7(void)
{
   exampleSOS1(FALSE);
}

/* TEST 14 */
/** @brief compute signed permutation symmetries for an example containing SOS1 constraints */
void test_test_compute_symmetry_special8(void)
{
   exampleSOS1(TRUE);
}

/* TEST 15 */
/** @brief compute permutation symmetries for an example containing SOS2 constraints */
void test_test_compute_symmetry_special9(void)
{
   exampleSOS2(FALSE);
}

/* TEST 16 */
/** @brief compute signed permutation symmetries for an example containing SOS2 constraints */
void test_test_compute_symmetry_special10(void)
{
   exampleSOS2(TRUE);
}

/* TEST 17 */
/** @brief compute signed permutation symmetries for an example containing PB constraints */
void test_test_compute_symmetry_special11(void)
{
   examplePB(FALSE);
}

/* TEST 18 */
/** @brief compute signed permutation symmetries for an example containing PB constraints */
void test_test_compute_symmetry_special12(void)
{
   examplePB(TRUE);
}

/* TEST 19 */
/** @brief compute permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr1(void)
{
   exampleExpr1(FALSE);
}

/* TEST 20 */
/** @brief compute signed permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr2(void)
{
   exampleExpr1(TRUE);
}

/* TEST 21 */
/** @brief compute permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr3(void)
{
   exampleExpr2(FALSE);
}

/* TEST 22 */
/** @brief compute signed permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr4(void)
{
   exampleExpr2(TRUE);
}

/* TEST 23 */
/** @brief compute permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr5(void)
{
   exampleExpr3(FALSE);
}

/* TEST 24 */
/** @brief compute signed permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr6(void)
{
   exampleExpr3(TRUE);
}

/* TEST 25 */
/** @brief compute permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr7(void)
{
   exampleExpr4(FALSE);
}

/* TEST 26 */
/** @brief compute signed permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr8(void)
{
   exampleExpr4(TRUE);
}

/* TEST 27 */
/** @brief compute permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr9(void)
{
   exampleExpr5(FALSE);
}

/* TEST 28 */
/** @brief compute signed permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr10(void)
{
   exampleExpr5(TRUE);
}

/* TEST 29 */
/** @brief compute permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr11(void)
{
   exampleExpr6(FALSE);
}

/* TEST 30 */
/** @brief compute signed permutation symmetries for an example containing nonlinear constraints */
void test_test_compute_symmetry_expr12(void)
{
   exampleExpr6(TRUE);
}

/* TEST 31 (doublelex matrices) */
/** @brief detect action corresponding to double lex matrices */
void test_test_compute_symmetry_doublelex(void)
{
   int perm1[20] = {1,0,2,3,5,4,6,7,9,8,10,11,13,12,14,15,17,16,18,19};
   int perm2[20] = {0,1,3,2,4,5,7,6,8,9,11,10,12,13,15,14,16,17,19,18};
   int perm3[20] = {4,5,6,7,0,1,2,3,8,9,10,11,12,13,14,15,16,17,18,19};
   int perm4[20] = {0,1,2,3,8,9,10,11,4,5,6,7,12,13,14,15,16,17,18,19};
   int perm5[20] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
   int* perms[5];
   SCIP_Bool success;
   SCIP_Bool isorbitope;
   int** lexmatrix = NULL;
   int* lexrowsbegin = NULL;
   int* lexcolsbegin = NULL;
   int nrows = -1;
   int ncols = -1;
   int nrowmatrices = -1;
   int ncolmatrices = -1;
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

   SCIP_CALL( SCIPdetectSingleOrDoubleLexMatrices(scip, FALSE, perms, 5, 20,
         &success, &isorbitope, &lexmatrix, &nrows, &ncols,
         &lexrowsbegin, &lexcolsbegin, &nrowmatrices, &ncolmatrices) );

   TEST_ASSERT( success );
   TEST_ASSERT( lexmatrix != NULL );
   TEST_ASSERT( lexrowsbegin != NULL );
   TEST_ASSERT( lexcolsbegin != NULL );
   TEST_ASSERT( nrows == 5 );
   TEST_ASSERT( ncols == 4 );
   TEST_ASSERT( nrowmatrices == 2 );
   TEST_ASSERT( ncolmatrices == 2 );

   SCIPfreeBlockMemoryArray(scip, &lexcolsbegin, ncolmatrices + 1);
   SCIPfreeBlockMemoryArray(scip, &lexrowsbegin, nrowmatrices + 1);
   for (i = 0; i < nrows; ++i)
   {
      SCIPfreeBlockMemoryArray(scip, &lexmatrix[i], ncols);
   }
   SCIPfreeBlockMemoryArray(scip, &lexmatrix, nrows);
}

/* TEST 32 symmetry computation of SDG for permutation symmetries */
/** @brief detect symmetries of full SDG */
void test_test_compute_symmetry_symsdgnodes(void)
{
   SCIP_VAR* vars[4];
   SYM_GRAPH* graph;
   int opnode1;
   int opnode2;
   int opnode3;
   int nperms;
   int nmaxperms;
   int** perms;
   SCIP_Real log10groupsize;
   SCIP_Real symcodetime;
   int i;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   SCIP_CALL( SCIPcreateProbBasic(scip, "basic1"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[0], "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, vars[0]) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[1], "x2", 0.0, 1.0, 2.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, vars[1]) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[2], "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, vars[2]) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[3], "x4", 0.0, 1.0, 2.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, vars[3]) );

   /* create some symmetry detection graph  */
   SCIP_CALL( SCIPcreateSymgraph(scip, SYM_SYMTYPE_PERM, &graph, vars, 4, 3, 0, 0, 6, 0.0) );

   SCIP_CALL( SCIPaddSymgraphOpnode(scip, graph, 0, &opnode1) );
   SCIP_CALL( SCIPaddSymgraphOpnode(scip, graph, 1, &opnode2) );
   SCIP_CALL( SCIPaddSymgraphOpnode(scip, graph, 1, &opnode3) );

   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode1, opnode2, FALSE, 0.0) );
   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode2, SCIPgetSymgraphVarnodeidx(scip, graph, vars[0]), FALSE, 0.0) );
   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode2, SCIPgetSymgraphVarnodeidx(scip, graph, vars[1]), FALSE, 0.0) );

   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode1, opnode3, FALSE, 0.0) );
   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode3, SCIPgetSymgraphVarnodeidx(scip, graph, vars[2]), FALSE, 0.0) );
   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode3, SCIPgetSymgraphVarnodeidx(scip, graph, vars[3]), FALSE, 0.0) );

   /* compute its colors and compute symmetries */
   SCIP_CALL( SCIPcomputeSymgraphColors(scip, graph, 0) );
   SCIP_CALL( SYMcomputeSymmetryGeneratorsNode(scip, 0, graph, &nperms, &nmaxperms, &perms, &log10groupsize, &symcodetime) );

   TEST_ASSERT( nperms == 1 );
   /* check images of operator nodes */
   TEST_ASSERT( perms[0][opnode1] == opnode1 );
   TEST_ASSERT( perms[0][opnode2] == opnode3 );
   TEST_ASSERT( perms[0][opnode3] == opnode2 );
   /* check images of variables (1 -> 3, 2 -> 4)*/
   TEST_ASSERT( perms[0][opnode3 + 1] == opnode3 + 3 );
   TEST_ASSERT( perms[0][opnode3 + 2] == opnode3 + 4 );
   TEST_ASSERT( perms[0][opnode3 + 3] == opnode3 + 1 );
   TEST_ASSERT( perms[0][opnode3 + 4] == opnode3 + 2 );

   SCIPfreeBlockMemoryArray(scip, &perms[0], 7);
   SCIPfreeBlockMemoryArray(scip, &perms, nmaxperms);

   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   }
}

/* TEST 33 symmetry computation of SDG for reflection symmetries */
/** @brief detect symmetries of full SDG */
void test_test_compute_symmetry_symsdgnodes2(void)
{
   SCIP_VAR* vars[2];
   SYM_GRAPH* graph;
   int opnode1;
   int opnode2;
   int opnode3;
   int nperms;
   int nmaxperms;
   int** perms;
   SCIP_Real log10groupsize;
   SCIP_Real symcodetime;
   int i;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   SCIP_CALL( SCIPcreateProbBasic(scip, "basic1"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[0], "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, vars[0]) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[1], "x2", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, vars[1]) );

   /* create some symmetry detection graph  */
   SCIP_CALL( SCIPcreateSymgraph(scip, SYM_SYMTYPE_SIGNPERM, &graph, vars, 2, 3, 0, 0, 6, 0.0) );

   SCIP_CALL( SCIPaddSymgraphOpnode(scip, graph, 0, &opnode1) );
   SCIP_CALL( SCIPaddSymgraphOpnode(scip, graph, 1, &opnode2) );
   SCIP_CALL( SCIPaddSymgraphOpnode(scip, graph, 1, &opnode3) );

   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode1, opnode2, FALSE, 0.0) );
   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode2, SCIPgetSymgraphVarnodeidx(scip, graph, vars[0]), TRUE, 1.0) );
   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode2, SCIPgetSymgraphNegatedVarnodeidx(scip, graph, vars[0]),
         TRUE, -1.0) );

   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode1, opnode3, FALSE, 0.0) );
   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode3, SCIPgetSymgraphVarnodeidx(scip, graph, vars[1]), TRUE, -1.0) );
   SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, opnode3, SCIPgetSymgraphNegatedVarnodeidx(scip, graph, vars[1]),
         TRUE, 1.0) );

   /* compute its colors and compute symmetries */
   SCIP_CALL( SCIPcomputeSymgraphColors(scip, graph, 0) );
   SCIP_CALL( SYMcomputeSymmetryGeneratorsNode(scip, 0, graph, &nperms, &nmaxperms, &perms, &log10groupsize, &symcodetime) );

   TEST_ASSERT( nperms == 1 );
   /* check images of operator nodes */
   TEST_ASSERT( perms[0][opnode1] == opnode1 );
   TEST_ASSERT( perms[0][opnode2] == opnode3 );
   TEST_ASSERT( perms[0][opnode3] == opnode2 );
   /* check images of variables (1 -> -2, 2 -> -1) */
   TEST_ASSERT( perms[0][opnode3 + 1] == opnode3 + 4 );
   TEST_ASSERT( perms[0][opnode3 + 2] == opnode3 + 3 );
   TEST_ASSERT( perms[0][opnode3 + 3] == opnode3 + 2 );
   TEST_ASSERT( perms[0][opnode3 + 4] == opnode3 + 1 );

   SCIPfreeBlockMemoryArray(scip, &perms[0], 7);
   SCIPfreeBlockMemoryArray(scip, &perms, nmaxperms);

   for( i = 0; i < 2; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   }
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_test_compute_symmetry_basic1);
   RUN_TEST(test_test_compute_symmetry_basic2);
   RUN_TEST(test_test_compute_symmetry_basic3);
   RUN_TEST(test_test_compute_symmetry_basic4);
   RUN_TEST(test_test_compute_symmetry_basic5);
   RUN_TEST(test_test_compute_symmetry_basic6);
   RUN_TEST(test_test_compute_symmetry_doublelex);
   RUN_TEST(test_test_compute_symmetry_expr1);
   RUN_TEST(test_test_compute_symmetry_expr10);
   RUN_TEST(test_test_compute_symmetry_expr11);
   RUN_TEST(test_test_compute_symmetry_expr12);
   RUN_TEST(test_test_compute_symmetry_expr2);
   RUN_TEST(test_test_compute_symmetry_expr3);
   RUN_TEST(test_test_compute_symmetry_expr4);
   RUN_TEST(test_test_compute_symmetry_expr5);
   RUN_TEST(test_test_compute_symmetry_expr6);
   RUN_TEST(test_test_compute_symmetry_expr7);
   RUN_TEST(test_test_compute_symmetry_expr8);
   RUN_TEST(test_test_compute_symmetry_expr9);
   RUN_TEST(test_test_compute_symmetry_special1);
   RUN_TEST(test_test_compute_symmetry_special10);
   RUN_TEST(test_test_compute_symmetry_special11);
   RUN_TEST(test_test_compute_symmetry_special12);
   RUN_TEST(test_test_compute_symmetry_special2);
   RUN_TEST(test_test_compute_symmetry_special3);
   RUN_TEST(test_test_compute_symmetry_special4);
   RUN_TEST(test_test_compute_symmetry_special5);
   RUN_TEST(test_test_compute_symmetry_special6);
   RUN_TEST(test_test_compute_symmetry_special7);
   RUN_TEST(test_test_compute_symmetry_special8);
   RUN_TEST(test_test_compute_symmetry_special9);
   RUN_TEST(test_test_compute_symmetry_symsdgnodes);
   RUN_TEST(test_test_compute_symmetry_symsdgnodes2);
   return UNITY_END();
}
