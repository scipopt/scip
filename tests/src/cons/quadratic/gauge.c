/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   unittest-cons_quadratic.c
 * @brief  unit test for cons_quadratic methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_quadratic.c"

#include "include/scip_test.h"

/*
 *****************************************
 *     TESTS FOR SEPARATION METHODS      *
 *****************************************
*/

#define EPS 1e-6
Test(separation, gauge, .description = "check computation of gauge function of a convex quadratic set")
{
   /* in a unittest we somehow assume that all the "setup" works.. I guess */
   SCIP* scip;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_SOL* point;
   SCIP_Bool success;
   SCIP_Real val;

   SCIP_CALL( SCIPcreate(&scip) );

   printf(" nnlpis %d\n", SCIPgetNNlpis(scip));
   /* include quadratic conshdlr (need to include nonlinear) */
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );
   printf(" nnlpis %d\n", SCIPgetNNlpis(scip));
   SCIP_CALL( SCIPincludeConshdlrQuadratic(scip) );
   printf(" nnlpis %d\n", SCIPgetNNlpis(scip));


   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
   printf(" nnlpis %d\n", SCIPgetNNlpis(scip));

   /* create convex quadratic constraint */
   {
      SCIP_VAR* vars[2];
      SCIP_Real vals[2];

      /* create variables */
      SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(scip, xvar) );
      SCIP_CALL( SCIPaddVar(scip, yvar) );

      /* create quadratic constraint: x^2 + y^2 <= 1 */
      vars[0] = xvar;
      vars[1] = yvar;

      vals[0] = 1.0;
      vals[1] = 1.0;

      SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "test", 0, NULL, NULL, 2, vars, vars, vals, -SCIPinfinity(scip), 1.0) );

      /* tell scip constraint is convex */
      consdata = SCIPconsGetData(cons);
      consdata->isconvex = TRUE;

      SCIP_CALL( SCIPaddCons(scip, cons) );
   }

   /* testing */
   {
      SCIP_Real expectedarray[2] = {0.0, 0.0};

      /* compute gauge */
      conshdlr = SCIPconsGetHdlr(cons);
      printf(" nnlpis %d\n", SCIPgetNNlpis(scip));
      SCIP_CALL( computeGauge(scip, conshdlr, cons) );

      /* this could be a single test probably */
      cr_assert(consdata->isgaugeavailable, "gauge unavailable, pointless to continue");
      cr_assert_arr_eq(consdata->interiorpoint, expectedarray, 2, "wrong interior point");
      cr_assert_float_eq(consdata->interiorpointval, 0.0, EPS, "wrong interior point value");
      cr_assert_arr_eq(consdata->gaugecoefs, expectedarray, 2, "wrong interior point");
      cr_assert_float_eq(consdata->gaugeconst, 0.0, EPS, "wrong gauge constant");

      /* create sol where to evaluate the gauge: one could parametrize this
       * or even create a theory out this, because in our case gauge(X) = X/norm(X)^2 */
      SCIP_CALL( SCIPcreateSol(scip, &point, NULL) );
      SCIP_CALL( SCIPsetSolVal(scip, point, xvar, 1.0) );
      SCIP_CALL( SCIPsetSolVal(scip, point, yvar, 1.0) );

      /* test evaluation */
      SCIP_CALL( evaluateGauge(scip, conshdlr, cons, point, &val, &success) );
      cr_assert(success, "unsuccessful evaluation of gauge");

      cr_assert_float_eq(val, 1.41421356237309504880, EPS, "wrong gauge evaluation");
   }

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );

}
