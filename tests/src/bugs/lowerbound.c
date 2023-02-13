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

/**@file   lowerbound.c
 * @brief  display strange behavior of lower bound
 * @author Marc Pfetsch
 *
 * This unit test is supposed to show the following behavior: By some error, it may happen that the lower bound is
 * larger than the optimal value. Previous SCIP versions accepted this behavior without warning/assert.
 *
 * The unit test builds up a problem in which the optimal value is SIZE/2. It adds a (wrong) relaxator that at the
 * beginning returns SIZE/2 + 1 and later SIZE/2 and a corresponding (integral) solution.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <signal.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/debug.h"
#include "include/scip_test.h"

/* UNIT TEST */
#define SIZE 10

/** relaxator data */
struct SCIP_RelaxData
{
   int ncalls;
};

/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecConstant)
{
   assert( lowerbound != NULL );
   assert( result != NULL );

   /* in the first nodes return a constant value that is larger than the optimal solution */
   if ( SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) <= 5 )
   {
      *lowerbound = SIZE/2 + 1;
   }
   else
   {
      int nvars;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      int j;

      /* in the later nodes return SIZE/2, which is the optimal value and an integral solution */
      *lowerbound = SIZE/2;

      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
      for (j = 0; j < nvars; ++j)
      {
         if ( j <= 1 )
            vals[j] = 0.0;
         else if ( j < nvars/2 + 2 )
            vals[j] = 1.0;
         else
            vals[j] = 0.0;
      }

      SCIP_CALL( SCIPsetRelaxSolVals(scip, relax, nvars, vars, vals, TRUE) );
      SCIP_CALL( SCIPmarkRelaxSolValid(scip, relax, TRUE) );

      SCIPfreeBufferArray(scip, &vals);
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}



/* TESTS */

/* This test constructs an instance on which the assert for testing the lower bound should fail. */
Test(lowerbound, constant, .signal = SIGABRT)
{
   char name[SCIP_MAXSTRLEN];
   SCIP* scip = NULL;
   SCIP_RELAX* relax = NULL;
   SCIP_VAR* vars[SIZE];
   SCIP_Real vals[SIZE];
   SCIP_CONS* cons;
   SCIP_Real rhs;
   int i;

   /* this test can only work with DEBUGSOL, so we make it pass otherwise */
   if ( ! SCIPwithDebugSol() )
      abort(); /* return SIGABORT */

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include constant relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, "constant", "constant relaxator", 101, 1, relaxExecConstant, NULL) );
   assert(relax != NULL);

   /* create a problem and disable the presolver */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   for (i = 0; i < SIZE; ++i)
   {
      SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &vars[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
      if ( i <= 2 )
         vals[i] = 1.0;
      else
         vals[i] = 2.0;
      SCIP_CALL( SCIPaddVar(scip, vars[i]) );
   }

   rhs = SIZE - 1;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "cons", SIZE, vars, vals, rhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) ); */

   /* turn off presolving */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

   /* turn off LP solving to use relaxator */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );

   /* turn off restarts */
   SCIP_CALL( SCIPsetIntParam(scip, "limits/restarts", 0) );

   /* turn off symmetry */
   SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 0) );

   SCIP_CALL( SCIPsolve(scip) );

   /* SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) ); */

   for (i = 0; i < SIZE; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   }

   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}
