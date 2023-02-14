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

/**@file   convexproj.c
 * @brief  unit test for sepa_convexproj.c methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_nonlinear.h"
#include "scip/nlpi_ipopt.h"

#include "scip/sepa_convexproj.c"

#include "include/scip_test.h"

#define EPS 1e-5

static SCIP* scip = NULL;
static SCIP_SEPA* sepa = NULL;
static SCIP_NLPI* nlpi = NULL;
static SCIP_NLROW* nlrow1 = NULL;
static SCIP_NLROW* nlrow2 = NULL;
static SCIP_NLROW* nlrow3 = NULL;
static SCIP_VAR* x;
static SCIP_VAR* y;

/* This test computes the projection of x = 0.75, y = 0.25 onto 
 * log(exp(x) + exp(y)) <= 1
 * x^2 <= y
 * 1.1*x+2.4*x^2 + 0.01*x*y + 0.3*y^2 + 0.2*log(0.5*exp(0.12*x+0.1)+2*exp(0.1*y)+0.7)  <= 0.5
 * x, y \in [-2.5, 3.1]
 * and then computes the gradient cut at the projected point.
 * It does so, by considering each possible combination of representing the constraints,
 * log(exp(x) + exp(y)) <= 1 or -log(exp(x) + exp(y)) >= -1,
 * etc.
 *
 * In the test we use the following methods:
 *    storeNonlinearConvexNlrows
 *    setQuadraticObj
 */

/* create log(exp(x) + exp(y)) <= 1 */
static
void createNlRow1(CONVEXSIDE convexside)
{
   SCIP_EXPR* expr;
   SCIP_EXPRCURV curvature;
   SCIP_Real lhs;
   SCIP_Real rhs;

   /* decide curvature */
   if( convexside == RHS )
   {
      SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"log(exp(<x>) + exp(<y>))", NULL, NULL, NULL) );

      curvature = SCIP_EXPRCURV_CONVEX;
      lhs = -SCIPinfinity(scip);
      rhs = 1.0;
   }
   else
   {
      SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"-(log(exp(<x>) + exp(<y>)))", NULL, NULL, NULL) );

      curvature = SCIP_EXPRCURV_CONCAVE;
      lhs = -1.0;
      rhs = SCIPinfinity(scip);
   }

   /* create nonlinear row for expr <= 1 */
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow1, "lg(sum exp)", 0.0, 0, NULL, NULL, expr, lhs, rhs, curvature) );

   /* release */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* create x^2 <= y */
static
void createNlRow2(CONVEXSIDE convexside)
{
   SCIP_EXPR* expr;
   SCIP_EXPRCURV curvature;
   SCIP_Real lhs;
   SCIP_Real rhs;

   /* decide curvature */
   if( convexside == RHS )
   {
      SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x>^2 - <y>", NULL, NULL, NULL) );

      curvature = SCIP_EXPRCURV_CONVEX;
      lhs = -SCIPinfinity(scip);
      rhs = 0.0;
   }
   else
   {
      SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<y> - <x>^2", NULL, NULL, NULL) );

      curvature = SCIP_EXPRCURV_CONCAVE;
      lhs = 0.0;
      rhs = SCIPinfinity(scip);
   }

   /* create nonlinear row for expr <= 1 */
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow2, "quad", 0.0, 0, NULL, NULL, expr, lhs, rhs, curvature) );

   /* release */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* 1.1*x+2.4*x^2 + 0.01*x*y + 0.3*y^2 + 0.2*log(0.5*exp(0.12*x+0.1)+2*exp(0.1*y)+0.7)  <= 0.5 */
static
void createNlRow3(CONVEXSIDE convexside)
{
   SCIP_EXPR* expr;
   SCIP_EXPRCURV curvature;
   SCIP_Real lhs;
   SCIP_Real rhs;

   /* decide curvature */
   if( convexside == RHS )
   {
      SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"1.1*<x>+2.4*<x>^2 + 0.01*<x>*<y> + 0.3*<y>^2 + 0.2*log(0.5*exp(0.12*<x>+0.1)+2*exp(0.1*<y>)+0.7)", NULL, NULL, NULL) );

      curvature = SCIP_EXPRCURV_CONVEX;
      lhs = -SCIPinfinity(scip);
      rhs = 0.5;
   }
   else
   {
      SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"-(1.1*<x>+2.4*<x>^2 + 0.01*<x>*<y> + 0.3*<y>^2 + 0.2*log(0.5*exp(0.12*<x>+0.1)+2*exp(0.1*<y>)+0.7))", NULL, NULL, NULL) );

      curvature = SCIP_EXPRCURV_CONCAVE;
      lhs = -0.5;
      rhs = SCIPinfinity(scip);
   }


   /* create nonlinear row for expr <= 1 */
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow3, "complicated", 0.0, 0, NULL, NULL, expr, lhs, rhs, curvature) );

   /* release */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* setup the sepadata */
static
void setup_sepadata(void)
{
   SCIP_SEPADATA* sepadata;
   SCIP_NLROW* nlrows[3] = {nlrow1, nlrow2, nlrow3};

   sepadata = SCIPsepaGetData(sepa);
   cr_assert(sepadata != NULL);

   SCIP_CALL( storeNonlinearConvexNlrows(scip, sepadata, nlrows, 3) );
   cr_assert_eq(sepadata->nnlrows, 3, "error: received %d nlrows", sepadata->nnlrows);

   /* initialize some of the sepadata */
   sepadata->nlpinvars = SCIPgetNVars(scip);
   cr_assert_eq(sepadata->nlpinvars, 2, "error: received %d vars", sepadata->nlpinvars);

   /* create nlpi problem */
   sepadata->nlpi = SCIPgetNlpis(scip)[0];
   cr_assert_not_null(sepadata->nlpi);

   SCIP_CALL( SCIPhashmapCreate(&sepadata->var2nlpiidx, SCIPblkmem(scip), sepadata->nlpinvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &sepadata->nlpivars, SCIPgetVars(scip), sepadata->nlpinvars) );

   /* I shouldn't care about the cutoff, just assert that the lp solution satisfies the cutoff bound */
   SCIP_CALL( SCIPcreateNlpiProblemFromNlRows(scip, sepadata->nlpi, &sepadata->nlpiprob, "convexproj-nlp-unittest", nlrows, 3,
            sepadata->var2nlpiidx, NULL, NULL, SCIPgetCutoffbound(scip), FALSE, TRUE) );

   /* set quadratic part of objective function */
   SCIP_CALL( setQuadraticObj(scip, sepadata) );
}

/* setup of test */
static
void test_setup(void)
{
   SCIP_VAR* auxx;

   SCIP_CALL( SCIPcreate(&scip) );

   /* if no IPOPT available, don't run test */
   if( ! SCIPisIpoptAvailableIpopt() )
      return;

   /* include NLPI's */
   SCIP_CALL( SCIPincludeNlpSolverIpopt(scip) );


   /* include convexproj separator and get it */
   SCIP_CALL( SCIPincludeSepaConvexproj(scip) );
   sepa = SCIPfindSepa(scip, "convexproj");
   cr_assert(sepa != NULL);

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &auxx, "x", -2.5, 3.1, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, auxx) );
   /* change SCIP's stage to be able to create nlrows and rows */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   x = SCIPvarGetTransVar(auxx);
   SCIP_CALL( SCIPreleaseVar(scip, &auxx) );
   /* create and add variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -2.5, 3.1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   if( nlrow1 != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow1) );
   }
   if( nlrow2 != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow2) );
   }
   if( nlrow3 != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow3) );
   }

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

static
void project(SCIP_Bool* isrhsconvex)
{
   SCIP_SEPADATA* sepadata;
   SCIP_SOL* toseparate_sol;
   SCIP_ROW* gradcut;
   SCIP_Real* coefs;
   SCIP_Real maxvio;
   SCIP_RESULT result;

   test_setup();
   /* if no IPOPT available, don't run test */
   if( nlpi == NULL )
      return;

   /* create the nl rows */
   createNlRow1(isrhsconvex[0]);
   createNlRow2(isrhsconvex[1]);
   createNlRow3(isrhsconvex[2]);

   /* setup nlpi and other data */
   setup_sepadata();
   sepadata = SCIPsepaGetData(sepa);

   /** point active only on nlrow3 **/
   /* set point to separate */
   SCIP_CALL( SCIPcreateSol(scip, &toseparate_sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, toseparate_sol, x, 0.75) );
   SCIP_CALL( SCIPsetSolVal(scip, toseparate_sol, y, 0.25) );

   /* compute violations */
   SCIP_CALL( computeMaxViolation(scip, sepadata, toseparate_sol, &maxvio) );
   cr_expect_float_eq(0.224076984, sepadata->constraintviolation[0], EPS, "got %f\n", sepadata->constraintviolation[0]);
   cr_expect_float_eq(0.3125, sepadata->constraintviolation[1], EPS, "got %f\n", sepadata->constraintviolation[1]);
   cr_expect_float_eq(1.937730557, sepadata->constraintviolation[2], EPS, "got %f\n", sepadata->constraintviolation[2]);

   /* project and compute cuts */
   SCIP_CALL( separateCuts(scip, sepa, toseparate_sol, &result) );

   /* check cut */
   cr_assert_eq(result, SCIP_SEPARATED, "result is %d instead of SEPARATED", result);
   cr_assert_eq(SCIPgetNCuts(scip), 1, "got %d cuts", SCIPgetNCuts(scip));
   gradcut = SCIPgetCuts(scip)[0];
   coefs = SCIProwGetVals(gradcut);
   if( isrhsconvex[2] )
   {
      cr_expect_float_eq(1.900159674905818, coefs[0], EPS, "for x got %f\n", coefs[0]);
      cr_expect_float_eq(0.138455761838885, coefs[1], EPS, "for y got %f\n", coefs[1]);
      cr_expect_float_eq(0.343035343633182, SCIProwGetRhs(gradcut), EPS, "wrong rhs got %f\n", SCIProwGetRhs(gradcut));
   }
   else
   {
      cr_expect_float_eq(-1.900159674905818, coefs[0], EPS, "for x got %f\n", coefs[0]);
      cr_expect_float_eq(-0.138455761838885, coefs[1], EPS, "for y got %f\n", coefs[1]);
      cr_expect_float_eq(-0.343035343633182, SCIProwGetLhs(gradcut), EPS, "wrong rhs got %f\n", SCIProwGetLhs(gradcut));
   }

   /* to remove the added cuts */
   SCIP_CALL( SCIPremoveInefficaciousCuts(scip) );
   SCIP_CALL( SCIPfreeSol(scip, &toseparate_sol) );
   teardown();
}

TheoryDataPoints(evaluation, convex_is_minus_concave) =
{
   DataPoints(SCIP_Bool, TRUE, FALSE),
   DataPoints(SCIP_Bool, TRUE, FALSE),
   DataPoints(SCIP_Bool, TRUE, FALSE)
};

Theory((SCIP_Bool is1convex, SCIP_Bool is2convex, SCIP_Bool is3convex), evaluation, convex_is_minus_concave)
{
   fprintf(stderr, "calling projection with %d %d %d\n",is1convex, is2convex, is3convex);
   project((SCIP_Bool[]){is1convex, is2convex, is3convex});
}
