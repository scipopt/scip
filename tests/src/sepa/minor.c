/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   minor.c
 * @brief  unit test for sepa_minor.c methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "nlpi/nlpi_ipopt.h"

#include "scip/sepa_minor.c"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

/** helper method for setup */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* include minor separator */
   SCIP_CALL( SCIPincludeSepaMinor(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.0, 1.0, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -1.0, 1.0,  1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -1.0, 1.0, -1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

/** helper method for teardown */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/** tests the detection of principal minors; the artificial problem contains two principal minors: for (x,y) and (y,z) */
Test(minor, detect, .init = setup, .fini = teardown)
{
   #define NCONSS 3
   const char* inputs[NCONSS] = {"[expr] <c1>: 1<= <x> * <x> + <y> * <y> <= 2",
      "[expr] <c2>: -0.5 <= <x> * <y> + <y> * <z> <= 0.5", "[expr] <c3>: -0.5 <= <z> * <z> <= 0.5"};
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int c;

   /* add two expression constraints (in transformed space) */
   for( c = 0; c < 3; ++c )
   {
      SCIP_CALL( SCIPparseCons(scip, &cons, inputs[c], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
         &success) );
      cr_assert(success);

      /* add and release constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* go to solving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );
   cr_assert(SCIPgetNConss(scip) == NCONSS);

   /* call the detection of the nonlinear handlers to create auxiliary variables */
   cr_assert(SCIPconshdlrGetNConss(conshdlr) == NCONSS);
   SCIP_CALL( SCIPdetectConsExprNlhdlrs(scip, conshdlr, SCIPconshdlrGetConss(conshdlr),
      SCIPconshdlrGetNConss(conshdlr), &infeasible) );
   cr_assert(!infeasible);

   /* get separator data */
   sepa = SCIPfindSepa(scip, SEPA_NAME);
   cr_assert(sepa != NULL);
   sepadata = SCIPsepaGetData(sepa);
   cr_assert(sepadata != NULL);

   /* call minor detection */
   cr_expect(!sepadata->detectedminors);
   cr_expect(sepadata->nminors == 0);
   SCIP_CALL( detectMinors(scip, sepadata) );
   cr_expect(sepadata->detectedminors);
   cr_expect(sepadata->nminors == 2, "nminors = %d (expected 2)", sepadata->nminors);
}

/** tests the eigenvalue and eigenvector computation */
Test(minor, eigenvals, .init = setup, .fini = teardown)
{
   SCIP_Real eigenvals[3];
   SCIP_Real eigenvecs[9];
   SCIP_Real xval = 1.5;
   SCIP_Real yval = 2.0;
   SCIP_Real xxval = -3.0;
   SCIP_Real yyval = 5.0;
   SCIP_Real xyval = -1.0;
   SCIP_Bool success;
   int i;

   /* LAPACK not available -> skip test */
   if( !SCIPisIpoptAvailableIpopt() )
      return;

   /* compute eigenvalues and eigenvectors */
   SCIP_CALL( getEigenValues(scip, xval, yval, xxval, yyval, xyval, eigenvals, eigenvecs, &success) );
   cr_assert(success);

   /* check whether A v_i = lambda_i v_i holds */
   for( i = 0; i < 1; ++i )
   {
      cr_assert(SCIPisRelEQ(scip,  1.0 * eigenvecs[3*i] +  xval * eigenvecs[3*i + 1] +  yval * eigenvecs[3*i + 2], eigenvals[i] * eigenvecs[3*i]));
      cr_assert(SCIPisRelEQ(scip, xval * eigenvecs[3*i] + xxval * eigenvecs[3*i + 1] + xyval * eigenvecs[3*i + 2], eigenvals[i] * eigenvecs[3*i + 1]));
      cr_assert(SCIPisRelEQ(scip, yval * eigenvecs[3*i] + xyval * eigenvecs[3*i + 1] + yyval * eigenvecs[3*i + 2], eigenvals[i] * eigenvecs[3*i + 2]));
   }
}