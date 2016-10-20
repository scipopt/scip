/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   unittest-qpkktref.c
 * @brief  Unittest for adding the KKT conditions to a quadratic program
 * @author Tobias Fischer
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "scip/dialog_default.h"


/** create quadratic problem with one linear constraint, set objective sense to minimize */
static
SCIP_RETCODE createProbLinear1(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-linear1") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create linear inequality */
   vars[0] = xvar;
   vars[1] = yvar;
   vars[2] = zvar;

   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "lower", 3, vars, vals, 0.25, 0.75) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one linear constraint, set objective sense to maximize */
static
SCIP_RETCODE createProbLinear2(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-linear2") );

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create linear inequality */
   vars[0] = xvar;
   vars[1] = yvar;
   vars[2] = zvar;

   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "lower", 3, vars, vals, 0.25, 0.75) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one knapsack constraint */
static
SCIP_RETCODE createProbKnapsack(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[2];
   SCIP_Longint vals[2];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-knapsack") );

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create knapsack inequality */
   vars[0] = xvar;
   vars[1] = yvar;

   vals[0] = 1;
   vals[1] = 1;

   SCIP_CALL( SCIPcreateConsBasicKnapsack(scip, &cons, "lower", 2, vars, vals, 1) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one set packing constraint */
static
SCIP_RETCODE createProbSetppc(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[3];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-setppc") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create set packing inequality */
   vars[0] = xvar;
   vars[1] = yvar;
   vars[2] = zvar;

   SCIP_CALL( SCIPcreateConsBasicSetpack(scip, &cons, "lower", 3, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one logicor constraint */
static
SCIP_RETCODE createProbLogicor(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[3];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-logicor") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, -1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -1.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, -2.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, 2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create logicor inequality */
   vars[0] = xvar;
   vars[1] = yvar;
   vars[2] = zvar;

   SCIP_CALL( SCIPcreateConsBasicLogicor(scip, &cons, "lower", 3, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one varbound */
static
SCIP_RETCODE createProbVarbound(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-varbound") );

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create varbound inequality */
   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, "lower", zvar, xvar, 1.0, -1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem */
static
SCIP_RETCODE createProb(
   SCIP*                 scip,               /**< SCIP instance */
   int                   instance            /**< instance to create */
   )
{
   switch( instance )
   {
      case 0:
         SCIP_CALL( createProbLinear1(scip) );
         break;
      case 1:
         SCIP_CALL( createProbLinear2(scip) );
         break;
      case 2:
         SCIP_CALL( createProbKnapsack(scip) );
         break;
      case 3:
         SCIP_CALL( createProbSetppc(scip) );
        break;
      case 4:
         SCIP_CALL( createProbLogicor(scip) );
        break;
      case 5:
         SCIP_CALL( createProbVarbound(scip) );
        break;
      default:
         SCIPerrorMessage("unknown instance number\n");
         return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** include QP settings (KKT conditions are not added, all presolvers are turned off) */
static
SCIP_RETCODE includeQPSettings(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   assert( scip != NULL );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* turn off presolving */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/trivial/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/inttobinary/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/gateextraction/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/dualcomp/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/stuffing/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/domcol/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/implics/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/components/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/qpkktref/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/dualfix/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/genvbounds/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/obbt/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/probing/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/pseudoobj/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/redcost/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/rootredcost/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/vbounds/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/soc/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/SOS1/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/SOS2/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/varbound/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/knapsack/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/setppc/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/linking/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/or/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/and/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/xor/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/conjunction/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/disjunction/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/linear/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/indicator/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/orbitope/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/logicor/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/bounddisjunction/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/cumulative/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/abspower/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/bivariate/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/quadratic/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/nonlinear/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/pseudoboolean/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/superindicator/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );

   return SCIP_OKAY;
}


/** include KKT settings (KKT conditions are added, all presolvers are turned off) */
static
SCIP_RETCODE includeKKTSettings(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   assert( scip != NULL );

   SCIP_CALL( includeQPSettings(scip) );

   SCIP_CALL( SCIPsetIntParam(scip, "presolving/qpkktref/maxrounds", -1) );

   /* allow variables to be unbounded */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/qpkktref/updatequadbounded", FALSE) );

   /* allow binary variables */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/qpkktref/addkktbinary", TRUE) );

   return SCIP_OKAY;
}


/** run unittest */
static
SCIP_RETCODE runUnittest(void)
{
   SCIP* scip1 = NULL;
   SCIP* scip2 = NULL;
   int j;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-qpkktref ===========\n");
   printf("=opt=  unittest-qpkktref 0\n\n");

   for (j = 0; j < 6; ++j)
   {
      SCIP_Bool equal = FALSE;

      /* initialize SCIP */
      SCIP_CALL( SCIPcreate(&scip1) );
      SCIP_CALL( SCIPcreate(&scip2) );

      /* include QP settings to scip1 (KKT conditions are not added, all presolvers are turned off) */
      SCIP_CALL( includeQPSettings(scip1) );

      /* include KKT settings to scip2 (KKT conditions are added, all presolvers are turned off) */
      SCIP_CALL( includeKKTSettings(scip2) );

      /* create problem */
      SCIP_CALL( createProb(scip1, j) );
      SCIP_CALL( createProb(scip2, j) );

      /* solve */
      SCIP_CALL( SCIPsolve(scip1) );
      SCIP_CALL( SCIPsolve(scip2) );
      /*
      SCIP_CALL( SCIPprintBestSol(scip1, NULL, FALSE) );
      SCIP_CALL( SCIPprintStatistics(scip1, NULL) );
      */
      /*SCIP_CALL( SCIPwriteTransProblem(scip2, "trafounittestQP.lp", NULL, FALSE ) );*/

      if ( SCIPisFeasEQ(scip1, SCIPgetPrimalbound(scip1), SCIPgetPrimalbound(scip2) ) )
         equal = TRUE;

      /* free transformed problem */
      SCIP_CALL( SCIPfreeTransform(scip1) );
      SCIP_CALL( SCIPfreeTransform(scip2) );

      /* free SCIP */
      SCIP_CALL( SCIPfree(&scip1) );
      SCIP_CALL( SCIPfree(&scip2) );

      /* check for memory leaks */
      BMScheckEmptyMemory();

      if ( ! equal )
      {
         SCIPerrorMessage("Optimal solution of original problem is not equal to optimal solution of reformulated problem.\n");
         return SCIP_ERROR;
      }
   }

   printf("Unit test for KKT-reformulation passed\n");

   /* for automatic testing output the following */
   printf("SCIP Status        : all tests passed\n");
   printf("Ignore the following:\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return SCIP_OKAY;
}

/** main function */
int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runUnittest();

   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
