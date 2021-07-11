/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
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

/* TEST 1: changing objective  */
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

/* TEST 2: adding inequalities  */
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
   int nruns = 100;
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
