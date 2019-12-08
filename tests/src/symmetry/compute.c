/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   compute.c
 * @brief  unit tests for computing symmetry
 * @author Marc Pfetsch
 * @author Fabian Wegscheider
 */

#include <scip/scip.h>
#include <scip/cons_expr.h>
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_var.h"
#include <scip/cons_expr_pow.h>
#include <include/scip_test.h>
#include <scip/symmetry.h>
#include <scip/prop_symmetry.h>
#include <symmetry/compute_symmetry.h>
#include <scip/scipdefplugins.h>
#include <scip/cons_expr_exp.h>
#include <scip/cons_expr_product.h>

static SCIP* scip;

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

/* TEST 4 */
Test(test_compute_symmetry, expr1, .description = "compute symmetry for a simple example with 4 variables and 2 expr constraints - before presolving")
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_VAR* w;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** permvars;
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_CONSEXPR_EXPR* yexpr;
   SCIP_CONSEXPR_EXPR* zexpr;
   SCIP_CONSEXPR_EXPR* wexpr;
   SCIP_CONSEXPR_EXPR** powexprs1;
   SCIP_CONSEXPR_EXPR** powexprs2;
   SCIP_CONSEXPR_EXPR* sumexpr1;
   SCIP_CONSEXPR_EXPR* sumexpr2;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int npermvars;
   int nperms;

   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs1, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs2, 2) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
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

   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &zexpr, z) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &wexpr, w) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs1[0], xexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs1[1], yexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs2[0], zexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs2[1], wexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr1, 2, powexprs1, NULL, 0) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr2, 2, powexprs2, NULL, 0) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e1", sumexpr1, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e2", sumexpr2, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

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

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs2[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs2[0]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs1[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs1[0]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIPfreeBufferArray(scip, &powexprs1);
   SCIPfreeBufferArray(scip, &powexprs2);
}

/* TEST 5 */
Test(test_compute_symmetry, expr2, .description = "compute symmetry for a more complex example with 5 variables and 3 expr constraints - before presolving")
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_VAR* w;
   SCIP_VAR* v;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** permvars;
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_CONSEXPR_EXPR* yexpr;
   SCIP_CONSEXPR_EXPR* zexpr;
   SCIP_CONSEXPR_EXPR* wexpr;
   SCIP_CONSEXPR_EXPR* vexpr;
   SCIP_CONSEXPR_EXPR** powexprs1;
   SCIP_CONSEXPR_EXPR** powexprs2;
   SCIP_CONSEXPR_EXPR* sumexpr1;
   SCIP_CONSEXPR_EXPR* sumexpr2;
   SCIP_CONSEXPR_EXPR* sumexpr3;
   SCIP_CONSEXPR_EXPR* expexpr1;
   SCIP_CONSEXPR_EXPR* expexpr2;
   SCIP_CONSEXPR_EXPR* expexpr3;
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

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
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

   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &zexpr, z) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &wexpr, w) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &vexpr, v) );

   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs1[0], xexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs1[1], yexpr, 3) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs1[2], zexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs1[3], wexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs1[4], vexpr, 2) );

   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs2[0], xexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs2[1], yexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs2[2], zexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs2[3], wexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs2[4], vexpr, 2) );

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr1, 5, powexprs1, NULL, 0) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr2, 5, powexprs2, vals1, 0) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr3, 5, powexprs2, vals2, 0) );

   SCIP_CALL( SCIPcreateConsExprExprExp(scip, conshdlr, &expexpr1, sumexpr1) );
   SCIP_CALL( SCIPcreateConsExprExprExp(scip, conshdlr, &expexpr2, sumexpr2) );
   SCIP_CALL( SCIPcreateConsExprExprExp(scip, conshdlr, &expexpr3, sumexpr3) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e1", expexpr1, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e2", expexpr2, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e3", expexpr3, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

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

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expexpr3) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expexpr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expexpr1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr3) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs2[4]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs2[3]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs2[2]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs2[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs2[0]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs1[4]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs1[3]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs1[2]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs1[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs1[0]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &vexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &v) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIPfreeBufferArray(scip, &powexprs1);
   SCIPfreeBufferArray(scip, &powexprs2);
}

/* TEST 6 */
Test(test_compute_symmetry, expr3, .description = "compute symmetry for a simple example with 4 variables, 2 expr constraints and 1 linear constaint - before presolving")
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
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_CONSEXPR_EXPR* yexpr;
   SCIP_CONSEXPR_EXPR* zexpr;
   SCIP_CONSEXPR_EXPR* wexpr;
   SCIP_CONSEXPR_EXPR** powexprs1;
   SCIP_CONSEXPR_EXPR** powexprs2;
   SCIP_CONSEXPR_EXPR* sumexpr1;
   SCIP_CONSEXPR_EXPR* sumexpr2;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int npermvars;
   int nperms;

   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs1, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &powexprs2, 2) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
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

   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &zexpr, z) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &wexpr, w) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs1[0], xexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs1[1], yexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs2[0], zexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexprs2[1], wexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr1, 2, powexprs1, NULL, 0) );
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr2, 2, powexprs2, NULL, 0) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e1", sumexpr1, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e2", sumexpr2, 1.0, 1.0) );
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

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs2[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs2[0]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs1[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexprs1[0]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIPfreeBufferArray(scip, &powexprs1);
   SCIPfreeBufferArray(scip, &powexprs2);
}

/* TEST 6 */
Test(test_compute_symmetry, expr4, .description = "compute symmetry for a simple example with 4 variables and 3 expr constraint - before presolving")
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* z;
   SCIP_VAR* w;
   SCIP_Real vals[2];
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** permvars;
   SCIP_CONSEXPR_EXPR** summands;
   SCIP_CONSEXPR_EXPR* xexpr;
   SCIP_CONSEXPR_EXPR* yexpr;
   SCIP_CONSEXPR_EXPR* zexpr;
   SCIP_CONSEXPR_EXPR* wexpr;
   SCIP_CONSEXPR_EXPR* powexpr1;
   SCIP_CONSEXPR_EXPR* powexpr2;
   SCIP_CONSEXPR_EXPR* prodexpr;
   SCIP_CONSEXPR_EXPR* sumexpr1;
   SCIP_CONSEXPR_EXPR* sumexpr2;
   SCIP_CONSEXPR_EXPR* sumexpr3;
   int** perms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int npermvars;
   int nperms;
   SCIP_Bool binvarsaffected;

   SCIP_CALL( SCIPallocBufferArray(scip, &summands, 3) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
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

   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &zexpr, z) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &wexpr, w) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexpr1, xexpr, 2) );
   SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &powexpr2, yexpr, 2) );
   summands[0] = xexpr;
   summands[1] = yexpr;
   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &prodexpr, 2, summands, 1) );

   summands[0] = powexpr1;
   summands[1] = zexpr;
   vals[0] = 1.0;
   vals[1] = -4.0;
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr1, 2, summands, vals, 0) );

   summands[0] = powexpr2;
   summands[1] = wexpr;
   vals[0] = 1.0;
   vals[1] = -4.0;
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr2, 2, summands, vals, 0) );

   summands[0] = xexpr;
   summands[1] = yexpr;
   summands[2] = prodexpr;
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr3, 3, summands, NULL, 0) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e1", sumexpr1, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e2", sumexpr2, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "e3", sumexpr3, -SCIPinfinity(scip), 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "propagating/symmetry/checksymmetries", TRUE) );

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

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr3) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   SCIPfreeBufferArray(scip, &summands);
}