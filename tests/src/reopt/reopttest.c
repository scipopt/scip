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

/**@file   reopttest.c
 * @brief  reoptimization unit test
 * @author Marc Pfetsch
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static const char* testfilename1 = "../check/instances/MIP/flugpl.mps";
static const char* testfilename2 = "../check/instances/MIP/MANN_a9.clq.lp";

static SCIP* scip;
static SCIP* reoptscip;
static SCIP_RANDNUMGEN* randgen;
static unsigned int randomseed = 42;

/* setup for tests */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 2) );
#else
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
#endif

   SCIP_CALL( SCIPcreate(&reoptscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(reoptscip) );
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(reoptscip, "display/verblevel", 2) );
#else
   SCIP_CALL( SCIPsetIntParam(reoptscip, "display/verblevel", 0) );
#endif
   SCIP_CALL( SCIPenableReoptimization(reoptscip, TRUE) );

   SCIPcreateRandom(scip, &randgen, randomseed, TRUE);
}

/* teardown for tests */
static
void teardown(void)
{
   SCIPfreeRandom(scip, &randgen);
   SCIPfree(&scip);
   SCIPfree(&reoptscip);

   BMScheckEmptyMemory();
}

TestSuite(reopt, .init = setup, .fini = teardown);

/* TEST 1: changing objective */
Test(reopt, objective)
{
   SCIP_VAR** reoptvars;
   SCIP_VAR** vars;
   SCIP_Real newobj;
   SCIP_STATUS status;
   SCIP_STATUS reoptstatus;
   SCIP_Real optval;
   SCIP_Real reoptval;
   SCIP_Real* objvals;
   SCIP_Real maxobj = 0.0;
   int* idx;
   int nvars;
   int nruns = 5;
   int i;
   int j;

   /* read problem */
   SCIP_CALL( SCIPreadProb(scip, testfilename1, NULL) );
   SCIP_CALL( SCIPreadProb(reoptscip, testfilename1, NULL) );

   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);

   assert( nvars == SCIPgetNOrigVars(reoptscip) );
   reoptvars = SCIPgetOrigVars(reoptscip);

   /* create arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idx, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &objvals, nvars) );
   for (j = 0; j < nvars; ++j)
   {
      objvals[j] = SCIPvarGetObj(vars[j]);
      idx[j] = j;
      if ( objvals[j] > maxobj )
         maxobj = objvals[j];
   }
   SCIPrandomPermuteIntArray(randgen, idx, 0, nvars);

   for (i = 0; i < nruns; ++i)
   {
      /* randomly change objective function of random variable */
      newobj = (SCIP_Real) SCIPrandomGetInt(randgen, -maxobj, maxobj);
      SCIP_CALL( SCIPchgVarObj(scip, vars[idx[i]], newobj) );

      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "Round %d: new objective value %g for variable <%s> (old: %g).\n", i+1, newobj, SCIPvarGetName(vars[idx[i]]), objvals[idx[i]]);
      objvals[idx[i]] = newobj;
      SCIP_CALL( SCIPchgReoptObjective(reoptscip, SCIPgetObjsense(scip), reoptvars, objvals, nvars) );

      /* solve both problems */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "         solving original problem ...\n");
      SCIP_CALL( SCIPsolve(scip) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "         solving reoptimization problem ...\n");
      SCIP_CALL( SCIPsolve(reoptscip) );

      status = SCIPgetStatus(scip);
      reoptstatus = SCIPgetStatus(reoptscip);

      cr_assert( status == reoptstatus );

      if ( status == SCIP_STATUS_OPTIMAL )
      {
         optval = SCIPgetPrimalbound(scip);
         reoptval = SCIPgetPrimalbound(reoptscip);

         SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "         objective = %g, reopt objecive = %g.\n", optval, reoptval);
         cr_assert_float_eq( optval, reoptval, 1e-6, "Optimal values are not equal: %g vs. %g.\n", optval, reoptval);
      }

      SCIP_CALL( SCIPfreeTransform(scip) );
      SCIP_CALL( SCIPfreeReoptSolve(reoptscip) );
   }

   SCIPfreeBlockMemoryArray(scip, &objvals, nvars);
   SCIPfreeBlockMemoryArray(scip, &idx, nvars);
}

/* TEST 2: changing objective for small instance */
Test(reopt, objectivesmall)
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_VAR* var3;
   SCIP_VAR* consvars[3];
   SCIP_Real consvals[3];
   SCIP_CONS* cons;
   SCIP_VAR** reoptvars;
   SCIP_VAR** vars;
   SCIP_Real newobj;
   SCIP_STATUS status;
   SCIP_STATUS reoptstatus;
   SCIP_Real optval;
   SCIP_Real reoptval;
   SCIP_Real* objvals;
   SCIP_Real maxobj = 0.0;
   int* idx;
   int nvars;
   int nruns = 3;
   int i;
   int j;

   /* create problems */
   SCIP_CALL( SCIPcreateProb(scip, "small1", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateProb(reoptscip, "reoptsmall", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* original */
   SCIP_CALL( SCIPcreateVar(scip, &var1, "x1", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, var1) );

   SCIP_CALL( SCIPcreateVar(scip, &var2, "x2", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, var2) );

   SCIP_CALL( SCIPcreateVar(scip, &var3, "x3", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, var3) );

   consvars[0] = var1;
   consvars[1] = var2;
   consvars[2] = var3;
   consvals[0] = 1.0;
   consvals[1] = 1.0;
   consvals[2] = 1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "pack", 3, consvars, consvals, -SCIPinfinity(scip), 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   consvars[0] = var1;
   consvars[1] = var2;
   consvals[0] = -1.0;
   consvals[1] = -1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "cover", 2, consvars, consvals, 1.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPreleaseVar(scip, &var1) );
   SCIP_CALL( SCIPreleaseVar(scip, &var2) );
   SCIP_CALL( SCIPreleaseVar(scip, &var3) );

   /* reopt */
   SCIP_CALL( SCIPcreateVar(reoptscip, &var1, "x1", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(reoptscip, var1) );

   SCIP_CALL( SCIPcreateVar(reoptscip, &var2, "x2", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(reoptscip, var2) );

   SCIP_CALL( SCIPcreateVar(reoptscip, &var3, "x3", 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(reoptscip, var3) );

   consvars[0] = var1;
   consvars[1] = var2;
   consvars[2] = var3;
   consvals[0] = 1.0;
   consvals[1] = 1.0;
   consvals[2] = 1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(reoptscip, &cons, "pack", 3, consvars, consvals, -SCIPinfinity(reoptscip), 1.0) );
   SCIP_CALL( SCIPaddCons(reoptscip, cons) );
   SCIP_CALL( SCIPreleaseCons(reoptscip, &cons) );

   consvars[0] = var1;
   consvars[1] = var2;
   consvals[0] = -1.0;
   consvals[1] = -1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(reoptscip, &cons, "cover", 2, consvars, consvals, 1.0, SCIPinfinity(reoptscip)) );
   SCIP_CALL( SCIPaddCons(reoptscip, cons) );
   SCIP_CALL( SCIPreleaseCons(reoptscip, &cons) );

   SCIP_CALL( SCIPreleaseVar(reoptscip, &var1) );
   SCIP_CALL( SCIPreleaseVar(reoptscip, &var2) );
   SCIP_CALL( SCIPreleaseVar(reoptscip, &var3) );

   /* loop preparation */
   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);

   assert( nvars == SCIPgetNOrigVars(reoptscip) );
   reoptvars = SCIPgetOrigVars(reoptscip);

   /* create arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idx, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &objvals, nvars) );
   for (j = 0; j < nvars; ++j)
   {
      objvals[j] = SCIPvarGetObj(vars[j]);
      idx[j] = j;
      if ( REALABS(objvals[j]) > maxobj )
         maxobj = REALABS(objvals[j]);
   }
   SCIPrandomPermuteIntArray(randgen, idx, 0, nvars);

   for (i = 0; i < nruns; ++i)
   {
      /* randomly change objective function of random variable */
      newobj = (SCIP_Real) SCIPrandomGetInt(randgen, -2 * maxobj, 2 * maxobj);
      SCIP_CALL( SCIPchgVarObj(scip, vars[idx[i]], newobj) );

      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "Round %d: new objective value %g for variable <%s> (old: %g).\n", i+1, newobj, SCIPvarGetName(vars[idx[i]]), objvals[idx[i]]);
      objvals[idx[i]] = newobj;
      SCIP_CALL( SCIPchgReoptObjective(reoptscip, SCIPgetObjsense(scip), reoptvars, objvals, nvars) );

      /* solve both problems */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "         solving original problem ...\n");
      SCIP_CALL( SCIPsolve(scip) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "         solving reoptimization problem ...\n");
      SCIP_CALL( SCIPsolve(reoptscip) );

      status = SCIPgetStatus(scip);
      reoptstatus = SCIPgetStatus(reoptscip);

      cr_assert( status == reoptstatus );

      if ( status == SCIP_STATUS_OPTIMAL )
      {
         optval = SCIPgetPrimalbound(scip);
         reoptval = SCIPgetPrimalbound(reoptscip);

         SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "         objective = %g, reopt objecive = %g.\n", optval, reoptval);
         cr_assert_float_eq( optval, reoptval, 1e-6, "Optimal values are not equal: %g vs. %g.\n", optval, reoptval);
      }

      SCIP_CALL( SCIPfreeTransform(scip) );
      SCIP_CALL( SCIPfreeReoptSolve(reoptscip) );
   }

   SCIPfreeBlockMemoryArray(scip, &objvals, nvars);
   SCIPfreeBlockMemoryArray(scip, &idx, nvars);
}

/* TEST 3: adding inequalities */
Test(reopt, conss)
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR** reoptvars;
   SCIP_VAR** vars;
   SCIP_STATUS status;
   SCIP_STATUS reoptstatus;
   SCIP_Real optval;
   SCIP_Real reoptval;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_SOL* bestsol;
   int nvars;
   int nruns = 5;
   int i;
   int j;

   /* read problem */
   SCIP_CALL( SCIPreadProb(scip, testfilename2, NULL) );
   SCIP_CALL( SCIPreadProb(reoptscip, testfilename2, NULL) );

   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);

   assert( nvars == SCIPgetNOrigVars(reoptscip) );
   reoptvars = SCIPgetOrigVars(reoptscip);

   /* create arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consvars, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consvals, nvars) );

   /* iteratively solve problem and forbid previous solution */
   for (i = 0; i < nruns; ++i)
   {
      /* solve both problems */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "Round %d: solving original problem ...\n", i+1);
      SCIP_CALL( SCIPsolve(scip) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "         solving reoptimization problem ...\n");
      SCIP_CALL( SCIPsolve(reoptscip) );

      status = SCIPgetStatus(scip);
      reoptstatus = SCIPgetStatus(reoptscip);

      cr_assert( status == reoptstatus );

      if ( status == SCIP_STATUS_OPTIMAL )
      {
         SCIP_CONS* cons;
         SCIP_Real solval;
         int nones = 0;

         optval = SCIPgetPrimalbound(scip);
         reoptval = SCIPgetPrimalbound(reoptscip);

         SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "         objective = %g, reopt objecive = %g.\n", optval, reoptval);
         cr_assert_float_eq( optval, reoptval, 1e-6, "Optimal values are not equal: %g vs. %g.\n", optval, reoptval);

         /* create no-good cut for previous solution */
         bestsol = SCIPgetBestSol(scip);
         for (j = 0; j < nvars; ++j)
         {
            /* assert( SCIPvarGetType(vars[j]) == SCIP_VARTYPE_BINARY ); */
            solval = SCIPgetSolVal(scip, bestsol, vars[j]);
            cr_assert( SCIPisFeasIntegral(scip, solval) );
            if ( solval < 0.5 )
               consvals[j] = 1.0;
            else
            {
               ++nones;
               consvals[j] = -1.0;
            }
         }

         SCIP_CALL( SCIPfreeTransform(scip) );
         SCIP_CALL( SCIPfreeReoptSolve(reoptscip) );

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nogood#%d", i);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, nvars, vars, consvals, 1.0 - (SCIP_Real) nones, SCIPinfinity(scip)) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

         SCIP_CALL( SCIPcreateConsBasicLinear(reoptscip, &cons, name, nvars, reoptvars, consvals, 1.0 - (SCIP_Real) nones, SCIPinfinity(scip)) );
         SCIP_CALL( SCIPaddCons(reoptscip, cons) );
         SCIP_CALL( SCIPreleaseCons(reoptscip, &cons) );
      }
      else
         break;
   }

   SCIPfreeBlockMemoryArray(scip, &consvars, nvars);
   SCIPfreeBlockMemoryArray(scip, &consvals, nvars);
}
