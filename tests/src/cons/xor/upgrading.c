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

/**@file   upgrading.c
 * @brief  unit tests for checking the upgrading of XOR constraints
 * @author Marc Pfetsch
 */

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/cons_xor.h"
#include "scip/scipdefplugins.h"

#define MAXNVARS 3

/** GLOBAL VARIABLES **/
static SCIP* scip = NULL;

SCIP_VAR* xvars[MAXNVARS];
SCIP_VAR* ivar;
SCIP_VAR* consvars[MAXNVARS+2];
SCIP_Real consvals[MAXNVARS+2];

/** check whether two lists contain the same variables */
static
SCIP_Bool compareVariableLists(
   SCIP_VAR**            vars1,              /**< variable list 1 */
   SCIP_VAR**            vars2,              /**< variable list 1 */
   int                   nvars               /**< number of variables */
   )
{
   int i;
   int j;

   for (i = 0; i < nvars; ++i)
   {
      for (j = 0; j < nvars; ++j)
      {
         if ( vars1[i] == vars2[j] )
            break;
      }
      if ( j >= nvars )
         return FALSE;
   }
   return TRUE;
}

/* TEST SUITE */

/** create SCIP instance */
static
void setup(void)
{
   char name[SCIP_MAXSTRLEN];
   int i;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "upgradexor") );

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create variables */
   for (i = 0; i < MAXNVARS; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &xvars[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, xvars[i]) );
   }
   SCIP_CALL( SCIPcreateVarBasic(scip, &ivar, "i", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, ivar) );
}

/** free SCIP instance */
static
void teardown(void)
{
   int i;

   SCIP_CALL( SCIPreleaseVar(scip, &ivar) );
   for (i = 0; i < MAXNVARS; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &xvars[i]) );
   }

   SCIPfree(&scip);
}

TestSuite(upgradexor, .init = setup, .fini = teardown);

/* TESTS  */
Test(upgradexor, intvariables1, .description = "test upgrading method of XOR constraints with two variables")
{
   SCIP_CONS* cons;
   SCIP_CONS* xorcons;
   SCIP_VAR* xvars_t[MAXNVARS];

   assert( MAXNVARS >= 2 );

   /* first constraint */
   consvars[0] = xvars[0];
   consvars[1] = xvars[1];
   consvars[2] = ivar;

   consvals[0] = 1.0;
   consvals[1] = 1.0;
   consvals[2] = 2.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 3, consvars, consvals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );
   SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &xorcons) );
   cr_assert( xorcons != NULL );

   /* check xor constraint */
   cr_assert( SCIPgetNVarsXor(scip, xorcons) == 2 );
   cr_assert( SCIPgetRhsXor(scip, xorcons) == FALSE );

   SCIP_CALL( SCIPgetTransformedVar(scip, xvars[0], &xvars_t[0]) );
   SCIP_CALL( SCIPgetTransformedVar(scip, xvars[1], &xvars_t[1]) );

   cr_assert( compareVariableLists(xvars_t, SCIPgetVarsXor(scip, xorcons), 2) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &xorcons) );
}

Test(upgradexor, intvariables2, .description = "test upgrading method of XOR constraints with three variables")
{
   SCIP_CONS* cons;
   SCIP_CONS* xorcons;
   SCIP_VAR* xvars_t[MAXNVARS];
   SCIP_VAR* ivar_t;

   assert( MAXNVARS >= 3 );

   /* first constraint */
   consvars[0] = xvars[0];
   consvars[1] = xvars[1];
   consvars[2] = xvars[2];
   consvars[3] = ivar;

   consvals[0] = 1.0;
   consvals[1] = 1.0;
   consvals[2] = 1.0;
   consvals[3] = 2.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 4, consvars, consvals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );
   SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &xorcons) );
   cr_assert( xorcons != NULL );

   /* check xor constraint */
   cr_assert( SCIPgetNVarsXor(scip, xorcons) == 3 );
   cr_assert( SCIPgetRhsXor(scip, xorcons) == FALSE );

   SCIP_CALL( SCIPgetTransformedVar(scip, xvars[0], &xvars_t[0]) );
   SCIP_CALL( SCIPgetTransformedVar(scip, xvars[1], &xvars_t[1]) );
   SCIP_CALL( SCIPgetTransformedVar(scip, xvars[2], &xvars_t[2]) );
   SCIP_CALL( SCIPgetTransformedVar(scip, ivar, &ivar_t) );

   cr_assert( compareVariableLists(xvars_t, SCIPgetVarsXor(scip, xorcons), 3) );
   cr_assert( SCIPvarGetAggrVar(ivar_t) == SCIPgetIntVarXor(scip, xorcons) ); /* the integer variable should be aggregated */

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &xorcons) );
}

Test(upgradexor, intvariables3, .description = "test upgrading method of XOR constraints with three variables and -2.0 coefficient ")
{
   SCIP_CONS* cons;
   SCIP_CONS* xorcons;
   SCIP_VAR* xvars_t[MAXNVARS];
   SCIP_VAR* ivar_t;

   assert( MAXNVARS >= 3 );

   /* first constraint */
   consvars[0] = xvars[0];
   consvars[1] = xvars[1];
   consvars[2] = xvars[2];
   consvars[3] = ivar;

   consvals[0] = 1.0;
   consvals[1] = 1.0;
   consvals[2] = 1.0;
   consvals[3] = -2.0;  /* use -2.0 */

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 4, consvars, consvals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );
   SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &xorcons) );
   cr_assert( xorcons != NULL );

   /* check xor constraint */
   cr_assert( SCIPgetNVarsXor(scip, xorcons) == 3 );
   cr_assert( SCIPgetRhsXor(scip, xorcons) == FALSE );

   SCIP_CALL( SCIPgetTransformedVar(scip, xvars[0], &xvars_t[0]) );
   SCIP_CALL( SCIPgetTransformedVar(scip, xvars[1], &xvars_t[1]) );
   SCIP_CALL( SCIPgetTransformedVar(scip, xvars[2], &xvars_t[2]) );
   SCIP_CALL( SCIPgetTransformedVar(scip, ivar, &ivar_t) );

   cr_assert( compareVariableLists(xvars_t, SCIPgetVarsXor(scip, xorcons), 3) );
   cr_assert( SCIPgetIntVarXor(scip, xorcons) == ivar_t ); /* the integer variable should remain unchanged */

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &xorcons) );
}

Test(upgradexor, intvariables4, .description = "test solving of upgraded XOR constraints")
{
   SCIP_VAR* zvars[4];
   SCIP_VAR* ivars[4];
   SCIP_CONS* cons;

   assert( MAXNVARS >= 3 );

   /* setup F7 dual binary matroid */
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvars[0], "z1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvars[1], "z2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvars[2], "z3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvars[3], "z4", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &ivars[0], "i1", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ivars[1], "i2", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ivars[2], "i3", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ivars[3], "i4", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );

   SCIP_CALL( SCIPaddVar(scip, zvars[0]) );
   SCIP_CALL( SCIPaddVar(scip, zvars[1]) );
   SCIP_CALL( SCIPaddVar(scip, zvars[2]) );
   SCIP_CALL( SCIPaddVar(scip, zvars[3]) );

   SCIP_CALL( SCIPaddVar(scip, ivars[0]) );
   SCIP_CALL( SCIPaddVar(scip, ivars[1]) );
   SCIP_CALL( SCIPaddVar(scip, ivars[2]) );
   SCIP_CALL( SCIPaddVar(scip, ivars[3]) );

   /* first constraint - use global xvars variables */
   consvars[0] = xvars[1];
   consvars[1] = xvars[2];
   consvars[2] = zvars[0];
   consvars[3] = ivars[0];

   consvals[0] = 1.0;
   consvals[1] = 1.0;
   consvals[2] = 1.0;
   consvals[3] = 2.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 4, consvars, consvals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* second constraint */
   consvars[0] = xvars[0];
   consvars[1] = xvars[2];
   consvars[2] = zvars[1];
   consvars[3] = ivars[1];

   consvals[0] = 1.0;
   consvals[1] = 1.0;
   consvals[2] = 1.0;
   consvals[3] = -2.0;   /* use negative */

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c2", 4, consvars, consvals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* third constraint */
   consvars[0] = xvars[0];
   consvars[1] = xvars[1];
   consvars[2] = zvars[2];
   consvars[3] = ivars[2];

   consvals[0] = 1.0;
   consvals[1] = 1.0;
   consvals[2] = 1.0;
   consvals[3] = 2.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c3", 4, consvars, consvals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* fourth constraint */
   consvars[0] = xvars[0];
   consvars[1] = xvars[1];
   consvars[2] = xvars[2];
   consvars[3] = zvars[3];
   consvars[4] = ivars[3];

   consvals[0] = 1.0;
   consvals[1] = 1.0;
   consvals[2] = 1.0;
   consvals[3] = 1.0;
   consvals[4] = 2.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c4", 5, consvars, consvals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "lp", FALSE) );
   SCIP_CALL( SCIPpresolve(scip) );
   SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
#endif

   /* solve */
   SCIP_CALL( SCIPsolve(scip) );

   cr_assert( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL );
   cr_assert_float_eq(SCIPgetPrimalbound(scip), 4.0, 1e-6, "Wrong optimal value.\n");

   SCIP_CALL( SCIPreleaseVar(scip, &zvars[0]) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvars[1]) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvars[2]) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvars[3]) );
   SCIP_CALL( SCIPreleaseVar(scip, &ivars[0]) );
   SCIP_CALL( SCIPreleaseVar(scip, &ivars[1]) );
   SCIP_CALL( SCIPreleaseVar(scip, &ivars[2]) );
   SCIP_CALL( SCIPreleaseVar(scip, &ivars[3]) );
}
