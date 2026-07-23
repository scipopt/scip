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

/**@file   objoffset.c
 * @brief  Unittest that the transformed-problem writer preserves the objective offset
 * @author Mathieu Dutour Sikiric
 *
 * Regression test for a bug where SCIPwriteTransProblem (MPS and LP writers) wrote
 * only the transformed problem's objective offset and dropped the objective constant
 * of the original problem (origprob->objoffset). Writing, reading back and solving
 * then returned an optimum that was off by exactly the dropped constant.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "include/scip_test.h"

/* problem: min x s.t. x >= 3, 0 <= x <= 10, with an objective constant of 100.
 * The optimum is therefore 3 + 100 = 103. */
#define OBJOFFSET 100.0
#define OPTIMUM   103.0

static SCIP* scip;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
}

TestSuite(readerobjoffset, .init = setup, .fini = teardown);

/** write the transformed problem, read it back into a fresh SCIP and solve it;
 *  the recovered optimum must still include the original objective constant.
 */
static
void checkFormat(const char* extension, const char* filename)
{
   SCIP_VAR* x;
   SCIP_CONS* cons;
   SCIP_Real vals[1];
   SCIP_VAR* vars[1];
   SCIP* scip2;

   /* build the original problem with a nonzero objective constant */
   SCIP_CALL( SCIPcreateProbBasic(scip, "objoffset") );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 10.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   vars[0] = x;
   vals[0] = 1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c", 1, vars, vals, 3.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPaddOrigObjoffset(scip, OBJOFFSET) );

   /* transform and write the transformed problem (this is where the offset was lost) */
   SCIP_CALL( SCIPtransformProb(scip) );
   SCIP_CALL( SCIPwriteTransProblem(scip, filename, extension, FALSE) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   /* read the written transformed problem back into an independent SCIP and solve */
   SCIP_CALL( SCIPcreate(&scip2) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip2) );
   SCIP_CALL( SCIPsetIntParam(scip2, "display/verblevel", 0) );
   SCIP_CALL( SCIPreadProb(scip2, filename, NULL) );
   SCIP_CALL( SCIPsolve(scip2) );

   /* with the buggy writer the constant was dropped and the optimum came out as 3.0 */
   cr_expect_float_eq(SCIPgetPrimalbound(scip2), OPTIMUM, 1e-6,
      "%s writer dropped the objective offset: optimum %g, expected %g",
      extension, SCIPgetPrimalbound(scip2), OPTIMUM);

   SCIP_CALL( SCIPfree(&scip2) );
   remove(filename);
}

Test(readerobjoffset, mps, .description = "transformed MPS write must keep the original objective offset")
{
   checkFormat("mps", "objoffset_trans.mps");
}

Test(readerobjoffset, lp, .description = "transformed LP write must keep the original objective offset")
{
   checkFormat("lp", "objoffset_trans.lp");
}
