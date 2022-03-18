/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   curvature.c
 * @brief  tests curvature for nonlinear constraint and row
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/nlpi_ipopt.h" /* to check whether LAPACK is around */
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 1.0, 2.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -4.0, -3.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(curvature, .init = setup, .fini = teardown);

/* check curvature in the constraint data and in the nonlinear rows */
Test(curvature, cons_and_nlrows)
{
   const char* inputs[3] = {
      "[nonlinear] <c1>: (<y>[C] + <z>[C])^2 <= 12;",
      "[nonlinear] <c2>: -(<x>[C] + <y>[C])^2 + <x>[C] + <y>[C] + <z>[C] == -5;",
      "[nonlinear] <c3>: <x>[C] * <y>[C] * <z>[C] - <x>[C] <= 1;"
      };
   SCIP_EXPRCURV targetcurvs[3] = {SCIP_EXPRCURV_CONVEX, SCIP_EXPRCURV_CONCAVE, SCIP_EXPRCURV_UNKNOWN};
   int ninputs = 3;
   int i;

   if( !SCIPisIpoptAvailableIpopt() )
   {
      targetcurvs[0] = SCIP_EXPRCURV_UNKNOWN;
      targetcurvs[1] = SCIP_EXPRCURV_UNKNOWN;
   }

   /* disable presolving */
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* create, add, and release nonlinear constraints */
   for( i = 0; i < ninputs; ++i )
   {
      SCIP_CONS* cons;
      SCIP_Bool success;

      SCIP_CALL( SCIPparseCons(scip, &cons, inputs[i],
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      cr_expect(success);
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* go to the solving stage; this should have triggered CONSINITSOL */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );
   cr_assert(SCIPgetNConss(scip) == ninputs);
   cr_assert(SCIPgetNNLPNlRows(scip) == ninputs);

   for( i = 0; i < ninputs; ++i )
   {
      SCIP_CONS* cons;
      SCIP_NLROW* nlrow;

      cons = SCIPgetConss(scip)[i];
      assert(cons != NULL);

      /* check curvature that is stored in the constraint data */
      cr_expect(SCIPgetCurvatureNonlinear(cons) == targetcurvs[i], "for cons %d (%s): expected %d got %d", i,
         SCIPconsGetName(cons), targetcurvs[i], SCIPgetCurvatureNonlinear(cons));

      /* check curvature that is stored in the nonlinear row */
      nlrow = SCIPgetNLPNlRows(scip)[i];
      assert(nlrow != NULL);
      cr_expect(SCIPnlrowGetCurvature(nlrow) == targetcurvs[i]);
   }
}

/* assume convenient curvature in the constraint data and in the nonlinear rows */
Test(curvature, assumeconvex)
{
   const char* inputs[3] = {
      "[nonlinear] <c1>: (<y>[C] + <z>[C])^2 <= 12;",
      "[nonlinear] <c2>: -(<x>[C] + <y>[C])^2 + <x>[C] + <y>[C] + <z>[C] == -5;",   /* this one will print a warning */
      "[nonlinear] <c3>: <x>[C] * <y>[C] * <z>[C] + <x>[C] >= 1;"
      };
   SCIP_EXPRCURV targetcurvs[3] = {SCIP_EXPRCURV_CONVEX, SCIP_EXPRCURV_LINEAR, SCIP_EXPRCURV_CONCAVE};
   int ninputs = 3;
   int i;

   /* disable presolving */
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* assume that constraints are convex */
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/nonlinear/assumeconvex", TRUE) );

   /* create, add, and release nonlinear constraints */
   for( i = 0; i < ninputs; ++i )
   {
      SCIP_CONS* cons;
      SCIP_Bool success;

      SCIP_CALL( SCIPparseCons(scip, &cons, inputs[i],
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      cr_expect(success);
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* go to the solving stage; this should have triggered CONSINITSOL */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );
   cr_assert(SCIPgetNConss(scip) == ninputs);
   cr_assert(SCIPgetNNLPNlRows(scip) == ninputs);

   for( i = 0; i < ninputs; ++i )
   {
      SCIP_CONS* cons;
      SCIP_NLROW* nlrow;

      cons = SCIPgetConss(scip)[i];
      assert(cons != NULL);

      /* check curvature that is stored in the constraint data */
      cr_expect(SCIPgetCurvatureNonlinear(cons) == targetcurvs[i], "for cons %d (%s): expected %d got %d", i,
         SCIPconsGetName(cons), targetcurvs[i], SCIPgetCurvatureNonlinear(cons));

      /* check curvature that is stored in the nonlinear row */
      nlrow = SCIPgetNLPNlRows(scip)[i];
      assert(nlrow != NULL);
      cr_expect(SCIPnlrowGetCurvature(nlrow) == targetcurvs[i]);
   }
}
