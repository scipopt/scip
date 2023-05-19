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

/**@file   fixedvar.c
 * @brief  unit test that checks that adding a linear constraint with a fixed variable works during solve
 * @author Stefan Vigerske
 */

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include <stdio.h>

/* TESTS  */
Test(fixedvar, addconswithfixedvar)
{
   SCIP* scip;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   char testfile[SCIP_MAXSTRLEN];
   int i;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* make things run faster */
   SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_EASYCIP, TRUE) );
   /* more output */
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", 1) );
   /* we want that the constraint that we add is propagated after it adds a row to the LP
    * to achieve this, we add the constraint in depth 1 and delay propagations until depth 3
    */
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/linear/propfreq", 3) );

   /* we take bell5.mps as our initial instance */
   strcpy(testfile, __FILE__);
   testfile[strlen(testfile) - 10] = '\0';  /* cutoff "fixedvar.c" */
   strcat(testfile, "../../../../check/instances/MIP/bell5.mps");
   SCIP_CALL( SCIPreadProb(scip, testfile, NULL) );

   /* the variable we want to SCIP to fix */
   SCIP_CALL( SCIPcreateVarBasic(scip, &var, "myvar", 0.0, 0.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var) );

   /* solve the root node */
   SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );
   SCIP_CALL( SCIPsolve(scip) );
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   /* now add a constraint that involves var and some more */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "mycons", 0, NULL, NULL, 0.0, SCIPinfinity(scip)) );
   for( i = 0; i < SCIPgetNOrigVars(scip) && i < 5; ++i )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, SCIPgetOrigVars(scip)[i], i + 1.0) );
   }
   SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* now continue solving until node 5 */
   SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 5) );
   SCIP_CALL( SCIPsolve(scip) );

   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   /** free SCIP */
   SCIP_CALL( SCIPfree(&scip) );
}
