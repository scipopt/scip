/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   compute.c
 * @brief  unit tests for computing symmetry
 * @author Marc Pfetsch
 */

#include <scip/scip.h>
#include <include/scip_test.h>
#include <scip/presol_symmetry.h>
#include <symmetry/compute_symmetry.h>
#include <scip/scipdefplugins.h>

static SCIP* scip;

/** setup: create SCIP */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* output external codes in order to see which esternal symmetry computation code is used */
   SCIPprintExternalCodes(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
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
   SCIP_PRESOL* presol;
   int npermvars;
   SCIP_VAR** permvars;
   int nperms;
   int** perms;
   int i;
   int j;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    *
    * x1 + x2          = 1
    *          x3 + x4 = 1
    * x1, ..., x4 binary
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

   /* get constraint handler */
   presol = SCIPfindPresol(scip, "symmetry");
   cr_assert( presol != NULL );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/symmetry/checksymmetries", TRUE) );

   /* mark that we want to have symmetry */
   SYMsetSpecRequirement(presol, SYM_SPEC_BINARY);

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetryGenerators(scip, presol, &npermvars, &permvars, &nperms, &perms) );
   cr_assert( nperms == 3 );

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

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}


/* TEST 2 */
Test(test_compute_symmetry, basic2, .description = "compute symmetry for a simple example with 4 variables and 2 linear constraints with fixed types")
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_PRESOL* presol;
   int npermvars;
   SCIP_VAR** permvars;
   int nperms;
   int** perms;
   int j;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    *
    * x1 + x2          = 1
    *          x3 + x4 = 1
    * x1, x2 binary; x3, x4 continuous
    */
   SCIP_CALL( SCIPcreateProbBasic(scip, "basic2"));

   SCIP_CALL( SCIPcreateVarBasic(scip, &var1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var1) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var2) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var3) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &var4, "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
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

   /* get constraint handler */
   presol = SCIPfindPresol(scip, "symmetry");
   cr_assert( presol != NULL );

   /* turn off presolving in order to avoid having trivial problem afterwards */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/symmetry/checksymmetries", TRUE) );

   /* mark that we want to have symmetry */
   SYMsetSpecRequirement(presol, SYM_SPEC_BINARY | SYM_SPEC_REAL);

   /* mark that we want to fix all binary variables */
   SYMsetSpecRequirementFixed(presol, SYM_SPEC_BINARY);

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetryGenerators(scip, presol, &npermvars, &permvars, &nperms, &perms) );
   cr_assert( nperms == 1 );

   SCIPinfoMessage(scip, NULL, "Permutation: (");
   for (j = 0; j < npermvars; ++j)
   {
      if ( j == 0 )
         SCIPinfoMessage(scip, NULL, "%d", perms[0][j]);
      else
         SCIPinfoMessage(scip, NULL, " %d", perms[0][j]);
   }
   SCIPinfoMessage(scip, NULL, ")\n");


   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}


/* TEST 3 */
Test(test_compute_symmetry, basic3, .description = "compute symmetry for a simple example with 4 variables and 4 linear constraints - before presolving")
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_PRESOL* presol;
   int npermvars;
   SCIP_VAR** permvars;
   int nperms;
   int** perms;
   int i;
   int j;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    *
    * x1 + x2          =  1
    *          x3 + x4 =  1
    * 2x1 +         x4 <= 2
    *     2x2 + x3     <= 2
    * x1, ..., x4 binary
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

   /* get constraint handler */
   presol = SCIPfindPresol(scip, "symmetry");
   cr_assert( presol != NULL );

   /* compute symmetry before presolving */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/symmetry/computepresolved", FALSE) );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/symmetry/checksymmetries", TRUE) );

   /* mark that we want to have symmetry */
   SYMsetSpecRequirement(presol, SYM_SPEC_BINARY);

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetryGenerators(scip, presol, &npermvars, &permvars, &nperms, &perms) );
   cr_assert( nperms == 1 );

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

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}


/* TEST 4 */
Test(test_compute_symmetry, basic4, .description = "compute symmetry for a simple example with 4 variables and 4 linear constraints - after presolving")
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* var4;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_PRESOL* presol;
   int npermvars;
   SCIP_VAR** permvars;
   int nperms;
   int** perms;
   int i;
   int j;

   /* skip test if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return;

   /* setup problem:
    *
    * x1 + x2          =  1
    *          x3 + x4 =  1
    * 2x1 +         x4 <= 2
    *     2x2 + x3     <= 2
    * x1, ..., x4 binary
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

   /* get constraint handler */
   presol = SCIPfindPresol(scip, "symmetry");
   cr_assert( presol != NULL );

   /* turn on checking of symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/symmetry/checksymmetries", TRUE) );

   /* mark that we want to have symmetry */
   SYMsetSpecRequirement(presol, SYM_SPEC_BINARY);

   /* presolve problem (symmetry will be available afterwards) */
   SCIP_CALL( SCIPpresolve(scip) );

   /* get symmetry */
   SCIP_CALL( SCIPgetSymmetryGenerators(scip, presol, &npermvars, &permvars, &nperms, &perms) );
   cr_assert( nperms == 1 );

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

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );
   SCIP_CALL( SCIPreleaseVar(scip, &var4) );
}
