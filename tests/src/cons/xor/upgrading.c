/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
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

/** GLOBAL VARIABLES **/
static SCIP* scip = NULL;

/* TEST SUITE */

/** create SCIP instance */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL (SCIPcreateProbBasic(scip, "upgradexor"));

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
}

/** free SCIP instance */
static
void teardown(void)
{
   SCIPfree(&scip);
}

TestSuite(upgradexor, .init = setup, .fini = teardown);

/* TESTS  */
Test(upgradexor, intvariables, .description = "test upgrading to XOR constraints using integer variables")
{
   SCIP_VAR* xvars[3];
   SCIP_VAR* zvars[4];
   SCIP_VAR* ivars[4];
   SCIP_VAR* consvars[5];
   SCIP_Real consvals[5];
   SCIP_CONS* cons;

   /* setup F7 dual binary matroid */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvars[0], "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvars[1], "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvars[2], "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &zvars[0], "z1", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvars[1], "z2", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvars[2], "z3", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvars[3], "z4", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &ivars[0], "i1", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ivars[1], "i2", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ivars[2], "i3", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ivars[3], "i4", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );

   SCIP_CALL( SCIPaddVar(scip, xvars[0]) );
   SCIP_CALL( SCIPaddVar(scip, xvars[1]) );
   SCIP_CALL( SCIPaddVar(scip, xvars[2]) );

   SCIP_CALL( SCIPaddVar(scip, zvars[0]) );
   SCIP_CALL( SCIPaddVar(scip, zvars[1]) );
   SCIP_CALL( SCIPaddVar(scip, zvars[2]) );
   SCIP_CALL( SCIPaddVar(scip, zvars[3]) );

   SCIP_CALL( SCIPaddVar(scip, ivars[0]) );
   SCIP_CALL( SCIPaddVar(scip, ivars[1]) );
   SCIP_CALL( SCIPaddVar(scip, ivars[2]) );
   SCIP_CALL( SCIPaddVar(scip, ivars[3]) );

   /* first constraint */
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

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "lp", FALSE) );

   /* solve */
   SCIP_CALL( SCIPsolve(scip) );

   cr_assert( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL );
   cr_assert_float_eq(SCIPgetPrimalbound(scip), 4.0, 1e-6, "Wrong optimal value.\n");

   SCIP_CALL( SCIPreleaseVar(scip, &xvars[0]) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvars[1]) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvars[2]) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvars[0]) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvars[1]) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvars[2]) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvars[3]) );
   SCIP_CALL( SCIPreleaseVar(scip, &ivars[0]) );
   SCIP_CALL( SCIPreleaseVar(scip, &ivars[1]) );
   SCIP_CALL( SCIPreleaseVar(scip, &ivars[2]) );
   SCIP_CALL( SCIPreleaseVar(scip, &ivars[3]) );
}
