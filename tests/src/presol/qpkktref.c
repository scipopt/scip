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

/**@file   qpkktref.c
 * @brief  Unittest for adding the KKT conditions to a quadratic program
 * @author Tobias Fischer
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/dialog_default.h"

#include "scip/presol_qpkktref.c"

#include "include/scip_test.h"


/** check number of added KKT constraints */
static
SCIP_RETCODE checkNAddConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side of linear constraint */
   SCIP_Real             rhs,                /**< right hand side of linear constraint */
   int                   naddconss           /**< number of added constraints */
   )
{
   int nconss = 0;
   int v;

   for (v = 0; v < nvars; ++v)
   {
      SCIP_VAR* var = vars[v];

      if ( SCIPvarIsBinary(var) )
         nconss += 4;   /* two linear constraints (for slack variables) and two SOS1 constraints */
      else
      {
         SCIP_Real ub;
         SCIP_Real lb;

         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);

         if ( ! SCIPisFeasZero(scip, lb) && ! SCIPisInfinity(scip, -lb) )
            nconss += 2;   /* one linear (for slack variable) and one SOS1 constraint */

         if ( ! SCIPisFeasZero(scip, ub) && ! SCIPisInfinity(scip, ub) )
            nconss += 2;   /* one linear (for slack variable) and one SOS1 constraint */
      }
   }

   if ( ! SCIPisFeasEQ(scip, lhs, rhs) )
   {
      if ( ! SCIPisInfinity(scip, -lhs) )
      {
         nconss += 2; /* one linear (for slack variable) and one SOS1 constraint */
      }

      if ( ! SCIPisInfinity(scip, rhs) )
      {
         nconss += 2; /* one linear (for slack variable) and one SOS1 constraint */
      }
   }

   if ( naddconss != nconss )
   {
      SCIPerrorMessage("failed test for number of added SOS1 constraints and constraints for slack variables.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/** check number of dual constraints */
static
SCIP_RETCODE checkNDualConss(
   SCIP*                 scip,               /**< SCIP instance */
   int                   nvars,              /**< number of variables */
   int                   ndualconss          /**< number of dual constraints */
   )
{
   if ( ndualconss != nvars )
   {
      SCIPerrorMessage("failed test for number of dual constraints.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/** checks linear constraint */
static
SCIP_RETCODE checkConsLinear(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_HASHMAP* varhash; /* hash map from variable to index of dual constraint */
   SCIP_CONS* objcons;
   SCIP_CONS** dualconss;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int ndualconss = 0;
   int naddconss = 0;
   int nvars;
   int j;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &objcons, "objcons", 0, NULL, NULL, 0.0, 0.0) );

   /* set up hash map */
   SCIP_CALL( SCIPhashmapCreate(&varhash, SCIPblkmem(scip), SCIPgetNVars(scip) + SCIPgetNFixedVars(scip)) );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &dualconss, 2 * SCIPgetNVars(scip) + 2 * SCIPgetNFixedVars(scip)) ); /*lint !e647*/

   lhs = SCIPgetLhsLinear(scip, cons);
   rhs = SCIPgetRhsLinear(scip, cons);
   nvars = SCIPgetNVarsLinear(scip, cons);
   vars = SCIPgetVarsLinear(scip, cons);
   vals = SCIPgetValsLinear(scip, cons);

   /* handle linear constraint for quadratic constraint update */
   SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPconsGetName(cons),
        vars, vals, lhs, rhs, nvars, varhash, dualconss, &ndualconss, &naddconss) );

   /* check number of added constraints */
   SCIP_CALL( checkNDualConss(scip, nvars, ndualconss) );
   SCIP_CALL( checkNAddConss(scip, vars, nvars, lhs, rhs, naddconss) );

   /* free buffer array */
   for (j = 0; j < ndualconss; ++j)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &dualconss[j]) );
   }
   SCIPfreeBufferArray(scip, &dualconss);

   /* free hash map */
   SCIPhashmapFree(&varhash);

   SCIP_CALL( SCIPreleaseCons(scip, &objcons) );

   return SCIP_OKAY;
}


/** checks knapsack constraint */
static
SCIP_RETCODE checkConsKnapsack(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONS*            cons                /**< knapsack constraint */
   )
{
   SCIP_HASHMAP* varhash; /* hash map from variable to index of dual constraint */
   SCIP_CONS* objcons;
   SCIP_CONS** dualconss;
   SCIP_VAR** vars;
   SCIP_Longint* weights;
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int ndualconss = 0;
   int naddconss = 0;
   int nvars;
   int v;
   int j;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &objcons, "objcons", 0, NULL, NULL, 0.0, 0.0) );

   /* set up hash map */
   SCIP_CALL( SCIPhashmapCreate(&varhash, SCIPblkmem(scip), SCIPgetNVars(scip) + SCIPgetNFixedVars(scip)) );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &dualconss, 2 * SCIPgetNVars(scip) + 2 * SCIPgetNFixedVars(scip)) ); /*lint !e647*/

   lhs = -SCIPinfinity(scip);
   rhs = (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons);
   nvars = SCIPgetNVarsKnapsack(scip, cons);
   vars = SCIPgetVarsKnapsack(scip, cons);
   weights = SCIPgetWeightsKnapsack(scip, cons);

   /* set coefficients of variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
   for (v = 0; v < nvars; ++v)
      vals[v] = (SCIP_Real) weights[v];

   /* handle linear constraint for quadratic constraint update */
   SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPconsGetName(cons),
         vars, vals, lhs, rhs, nvars, varhash, dualconss, &ndualconss, &naddconss) );

   /* check number of added constraints */
   SCIP_CALL( checkNDualConss(scip, nvars, ndualconss) );
   SCIP_CALL( checkNAddConss(scip, vars, nvars, lhs, rhs, naddconss) );

   /* free buffer arrays */
   for (j = 0; j < ndualconss; ++j)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &dualconss[j]) );
   }
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &dualconss);

   /* free hash map */
   SCIPhashmapFree(&varhash);

   SCIP_CALL( SCIPreleaseCons(scip, &objcons) );

   return SCIP_OKAY;
}


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

   /* check constraint */
   SCIP_CALL( checkConsLinear(scip, cons) );

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

   /* check constraint */
   SCIP_CALL( checkConsLinear(scip, cons) );

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

   /* check knapsack constraint */
   SCIP_CALL( checkConsKnapsack(scip, cons) );

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

   SCIP_CALL( SCIPsetIntParam(scip, "presolving/qpkktref/maxrounds", 0) );

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

   /* turn presolving of quadratic constraints off */
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/quadratic/maxprerounds", 0) );

   return SCIP_OKAY;
}

/* test suite */
TestSuite(qpkktref);

/** run unittest */
Test(qpkktref, runUnittest)
{
   SCIP* scip1 = NULL;
   SCIP* scip2 = NULL;
   int j;

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

      cr_assert( equal );

      /* check for memory leaks */
      BMScheckEmptyMemory();
   }
}
