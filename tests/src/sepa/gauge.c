/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   gauge.c
 * @brief  unit test for sepa_gauge methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/nlpi_ipopt.h"
#include "scip/expr_varidx.h"
#include "scip/expr_abs.h"
#include "scip/expr_exp.h"
#include "scip/expr_log.h"
#include "scip/expr_pow.h"
#include "scip/expr_product.h"
#include "scip/expr_sum.h"
#include "scip/expr_trig.h"
#include "scip/expr_var.h"
#include "scip/expr_value.h"

#include "scip/sepa_gauge.c"

#include "include/scip_test.h"

#define EPS 1e-6

static SCIP* scip = NULL;
static SCIP_Bool haveipopt = FALSE;
static SCIP_NLROW* nlrow1 = NULL;
static SCIP_NLROW* nlrow2 = NULL;
static SCIP_NLROW* nlrow3 = NULL;
static SCIP_VAR* x = NULL;
static SCIP_VAR* y = NULL;
static SCIP_EXPRITER* exprit = NULL;

/* methods to test:
 * computeInteriorPoint: uses SCIP's NLP to build a convex NLP relaxation, solves it and store the solution in the sepadata.
 * findBoundaryPoint: receives the nlrows, on which sides they are convex, indices of nlrows to separate, number of indices, the interior solution, the solution to separate and the final solution
 * generateCut: receives point where to compute gradient, the nlrow, the side which is convex and where to store the cut
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

static
void evaluation_setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include some expr handlers */
   SCIP_CALL( SCIPincludeExprhdlrAbs(scip) );
   SCIP_CALL( SCIPincludeExprhdlrExp(scip) );
   SCIP_CALL( SCIPincludeExprhdlrLog(scip) );
   SCIP_CALL( SCIPincludeExprhdlrVar(scip) );
   SCIP_CALL( SCIPincludeExprhdlrValue(scip) );
   SCIP_CALL( SCIPincludeExprhdlrSum(scip) );
   SCIP_CALL( SCIPincludeExprhdlrPow(scip) );
   SCIP_CALL( SCIPincludeExprhdlrProduct(scip) );
   SCIP_CALL( SCIPincludeExprhdlrSin(scip) );
   SCIP_CALL( SCIPincludeExprhdlrCos(scip) );

   /* if no IPOPT available, don't run test */
   if( ! SCIPisIpoptAvailableIpopt() )
      return;

   SCIPincludeNlpSolverIpopt(scip);
   haveipopt = TRUE;

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* change SCIP's stage to be able to create nlrows and rows */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* create and add variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -2.5, 3.1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -2.5, 3.1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create expression iterator */
   SCIP_CALL( SCIPcreateExpriter(scip, &exprit) );
}

static
void teardown(void)
{
   if( exprit != NULL )
      SCIPfreeExpriter(&exprit);

   if( x != NULL )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &x) );
   }
   if( y != NULL )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &y) );
   }
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
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

static
void evaluate_gauge(CONVEXSIDE* convexsides)
{
   SCIP_SOL* interior_sol;
   SCIP_SOL* toseparate_sol;
   SCIP_SOL* boundary_sol;
   POSITION position;
   SCIP_Real feas;
   int nlrowsidx[] = {0, 1};
   int nnlrowsidx = 2;
   SCIP_NLROW* nlrows[2];

   /* if no IPOPT available, don't run test */
   if( !haveipopt )
      return;

   /* create the nl rows */
   createNlRow1(convexsides[0]);
   createNlRow2(convexsides[1]);
   nlrows[0] = nlrow1;
   nlrows[1] = nlrow2;


   /** first point active only on nlrow1 **/
   /* set interior point */
   SCIP_CALL( SCIPcreateSol(scip, &interior_sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, interior_sol, x, -0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, interior_sol, y,  0.5) );
   /* set point to separate */
   SCIP_CALL( SCIPcreateSol(scip, &toseparate_sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, toseparate_sol, x, 0.75) );
   SCIP_CALL( SCIPsetSolVal(scip, toseparate_sol, y, 0.25) );

   /* find boundary point */
   position = EXTERIOR;
   SCIP_CALL( SCIPcreateSol(scip, &boundary_sol, NULL) );
   SCIP_CALL( findBoundaryPoint(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, interior_sol, toseparate_sol, boundary_sol, &position) );
   cr_expect_eq(position, BOUNDARY, "position of boundary sol is %d, expected boundary", position);

   cr_expect_float_eq(SCIPgetSolVal(scip, boundary_sol, x), 0.265033442069232, EPS, "for x got %f\n", SCIPgetSolVal(scip, boundary_sol, x));
   cr_expect_float_eq(SCIPgetSolVal(scip, boundary_sol, y), 0.346993311586154, EPS, "for y got %f\n", SCIPgetSolVal(scip, boundary_sol, y));

   SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrow1, boundary_sol, &feas) );
   cr_expect_float_eq(feas, 0.0, EPS, "nlrow1 is not active: feas %f\n", feas);
   SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrow2, boundary_sol, &feas) );
   cr_expect_float_neq(feas, 0.0, EPS, "nlrow2 is active: feas %f\n", feas);

   /** second point active on both **/
   /* set interior point */
   SCIP_CALL( SCIPsetSolVal(scip, interior_sol, x, -0.8) );
   SCIP_CALL( SCIPsetSolVal(scip, interior_sol, y,  0.75) );
   /* set point to separate */
   SCIP_CALL( SCIPsetSolVal(scip, toseparate_sol, x, 0.80011239) );
   SCIP_CALL( SCIPsetSolVal(scip, toseparate_sol, y, 0.0) );

   /* find boundary point */
   position = EXTERIOR;
   SCIP_CALL( findBoundaryPoint(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, interior_sol, toseparate_sol, boundary_sol, &position) );
   cr_expect_eq(position, BOUNDARY, "position of boundary sol is %d, expected boundary", position);

   cr_expect_float_eq(SCIPgetSolVal(scip, boundary_sol, x), 0.4213473910790077, EPS, "for x got %f\n", SCIPgetSolVal(scip, boundary_sol, x));
   cr_expect_float_eq(SCIPgetSolVal(scip, boundary_sol, y), 0.1775336239690863, EPS, "for y got %f\n", SCIPgetSolVal(scip, boundary_sol, y));

   SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrow1, boundary_sol, &feas) );
   cr_expect_float_eq(feas, 0.0, EPS, "nlrow1 is not active: feas %f\n", feas);
   SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrow2, boundary_sol, &feas) );
   cr_expect_float_eq(feas, 0.0, EPS, "nlrow2 is not active: feas %f\n", feas);

   /** third point active on nlrow2 **/
   /* set interior point */
   SCIP_CALL( SCIPsetSolVal(scip, interior_sol, x, -0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, interior_sol, y,  0.5) );
   /* set point to separate */
   SCIP_CALL( SCIPsetSolVal(scip, toseparate_sol, x, -1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, toseparate_sol, y,  0.0) );

   /* find boundary point */
   position = EXTERIOR;
   SCIP_CALL( findBoundaryPoint(scip, nlrows, nlrowsidx, nnlrowsidx, convexsides, interior_sol, toseparate_sol, boundary_sol, &position) );
   cr_expect_eq(position, BOUNDARY, "position of boundary sol is %d, expected boundary", position);

   cr_expect_float_eq(SCIPgetSolVal(scip, boundary_sol, x), -0.618033988749895, EPS, "for x got %f\n", SCIPgetSolVal(scip, boundary_sol, x));
   cr_expect_float_eq(SCIPgetSolVal(scip, boundary_sol, y),  0.381966011250105, EPS, "for y got %f\n", SCIPgetSolVal(scip, boundary_sol, y));

   SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrow1, boundary_sol, &feas) );
   cr_expect_float_neq(feas, 0.0, EPS, "nlrow1 is active: feas %f\n", feas);
   SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrow2, boundary_sol, &feas) );
   cr_expect_float_eq(feas, 0.0, EPS, "nlrow2 is not active: feas %f\n", feas);

   SCIP_CALL( SCIPfreeSol(scip, &interior_sol) );
   SCIP_CALL( SCIPfreeSol(scip, &toseparate_sol) );
   SCIP_CALL( SCIPfreeSol(scip, &boundary_sol) );
}

/* TEST SUITE */
TestSuite(evaluation, .init = evaluation_setup, .fini = teardown);

/* Test: with nlrows having different convexity */
Test(evaluation, convex_convex)
{
   evaluate_gauge((CONVEXSIDE[]){RHS, RHS});
}
Test(evaluation, concave_convex)
{
   evaluate_gauge((CONVEXSIDE[]){LHS, RHS});
}
Test(evaluation, concave_concave)
{
   evaluate_gauge((CONVEXSIDE[]){LHS, LHS});
}
Test(evaluation, convex_concave)
{
   evaluate_gauge((CONVEXSIDE[]){RHS, LHS});
}

/* test that we get the correct gradient cut for convex regions, defined as both, convex <= rhs and concave >= lhs */
Test(evaluation, gradient_cut_convex)
{
   SCIP_SOL* x0;
   SCIP_ROW* gradcut;
   SCIP_Bool convex = TRUE;
   SCIP_Real* coefs;
   SCIP_Bool success = FALSE;

   /* if no IPOPT available, don't run test */
   if( !haveipopt )
      return;

   /* create the nl row */
   createNlRow1(convex);

   /** point where to compute gradient cut **/
   SCIP_CALL( SCIPcreateSol(scip, &x0, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, x0, x, -0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, x0, y,  0.5) );

   /* compute gradient cut */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &gradcut, "gradcut", -SCIPinfinity(scip), SCIPinfinity(scip),
            FALSE, FALSE , FALSE) );
   SCIP_CALL( generateCut(scip, x0, nlrow1, convex, exprit, gradcut, &success) );
   cr_assert(success);

   /* check coefficients */
   coefs = SCIProwGetVals(gradcut);
   cr_expect_eq(SCIProwGetNNonz(gradcut), 2);
   cr_expect_float_eq(0.268941421369995, coefs[0], EPS, "for x got %f\n", coefs[0]);
   cr_expect_float_eq(0.731058578630005, coefs[1], EPS, "for y got %f\n", coefs[1]);
   cr_expect_float_eq(0.417796891111782, SCIProwGetRhs(gradcut), EPS, "wrong rhs\n");

   SCIP_CALL( SCIPfreeSol(scip, &x0) );
   SCIP_CALL( SCIPreleaseRow(scip, &gradcut) );
}

Test(evaluation, gradient_cut_concave)
{
   SCIP_SOL* x0;
   SCIP_ROW* gradcut;
   SCIP_Bool convex = FALSE;
   SCIP_Real* coefs;
   SCIP_Bool success = FALSE;

   /* if no IPOPT available, don't run test */
   if( !haveipopt )
      return;

   /* create the nl row */
   createNlRow1(convex);

   /** point where to compute gradient cut **/
   SCIP_CALL( SCIPcreateSol(scip, &x0, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, x0, x, -0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, x0, y,  0.5) );

   /* compute gradient cut */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &gradcut, "gradcut", -SCIPinfinity(scip), SCIPinfinity(scip),
            FALSE, FALSE , FALSE) );
   SCIP_CALL( generateCut(scip, x0, nlrow1, convex, exprit, gradcut, &success) );
   cr_assert(success);

   /* check coefficients */
   coefs = SCIProwGetVals(gradcut);
   cr_expect_eq(SCIProwGetNNonz(gradcut), 2);
   cr_expect_float_eq(-0.268941421369995, coefs[0], EPS, "for x got %f\n", coefs[0]);
   cr_expect_float_eq(-0.731058578630005, coefs[1], EPS, "for y got %f\n", coefs[1]);
   cr_expect_float_eq(-0.417796891111782, SCIProwGetLhs(gradcut), EPS, "wrong rhs\n");

   SCIP_CALL( SCIPfreeSol(scip, &x0) );
   SCIP_CALL( SCIPreleaseRow(scip, &gradcut) );
}

Test(evaluation, gradient_complicated_convex)
{
   SCIP_SOL* x0;
   SCIP_ROW* gradcut;
   SCIP_Bool convex = TRUE;
   SCIP_Real* coefs;
   SCIP_Bool success = FALSE;

   /* if no IPOPT available, don't run test */
   if( !haveipopt )
      return;

   /* create compicated nlrow */
   createNlRow3(convex);

   /** point where to compute gradient cut **/
   SCIP_CALL( SCIPcreateSol(scip, &x0, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, x0, x, -0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, x0, y,  0.5) );

   /* compute gradient cut */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &gradcut, "gradcut", -SCIPinfinity(scip), SCIPinfinity(scip),
            FALSE, FALSE , FALSE) );
   SCIP_CALL( generateCut(scip, x0, nlrow3, convex, exprit, gradcut, &success) );
   cr_assert(success);

   /* check coefficients */
   coefs = SCIProwGetVals(gradcut);
   cr_expect_eq(SCIProwGetNNonz(gradcut), 2);
   cr_expect_float_eq(-1.29124137035262, coefs[0], EPS, "for x got %f\n", coefs[0]);
   cr_expect_float_eq(0.307654681677814, coefs[1], EPS, "for y got %f\n", coefs[1]);
   cr_expect_float_eq(0.936777583155184, SCIProwGetRhs(gradcut), EPS, "wrong rhs\n");

   SCIP_CALL( SCIPfreeSol(scip, &x0) );
   SCIP_CALL( SCIPreleaseRow(scip, &gradcut) );
}

Test(evaluation, gradient_complicated_concave)
{
   SCIP_SOL* x0;
   SCIP_ROW* gradcut;
   SCIP_Bool convex = FALSE;
   SCIP_Real* coefs;
   SCIP_Bool success = FALSE;

   /* if no IPOPT available, don't run test */
   if( !haveipopt )
      return;

   /* create compicated nlrow */
   createNlRow3(convex);

   /** point where to compute gradient cut **/
   SCIP_CALL( SCIPcreateSol(scip, &x0, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, x0, x, -0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, x0, y,  0.5) );

   /* compute gradient cut */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &gradcut, "gradcut", -SCIPinfinity(scip), SCIPinfinity(scip),
            FALSE, FALSE , FALSE) );
   SCIP_CALL( generateCut(scip, x0, nlrow3, convex, exprit, gradcut, &success) );
   cr_assert(success);

   /* check coefficients */
   coefs = SCIProwGetVals(gradcut);
   cr_expect_eq(SCIProwGetNNonz(gradcut), 2);
   cr_expect_float_eq(1.29124137035262, coefs[0], EPS, "for x got %f\n", coefs[0]);
   cr_expect_float_eq(-0.307654681677814, coefs[1], EPS, "for y got %f\n", coefs[1]);
   cr_expect_float_eq(-0.936777583155184, SCIProwGetLhs(gradcut), EPS, "wrong rhs\n");

   SCIP_CALL( SCIPfreeSol(scip, &x0) );
   SCIP_CALL( SCIPreleaseRow(scip, &gradcut) );
}

/* test that check interior point, notice that it doesn't belong to the test suite evaluation evaluation */
Test(interior_point, compute_interior_point)
{
   SCIP_SEPADATA* sepadata;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeExprhdlrAbs(scip) );
   SCIP_CALL( SCIPincludeExprhdlrExp(scip) );
   SCIP_CALL( SCIPincludeExprhdlrLog(scip) );
   SCIP_CALL( SCIPincludeExprhdlrValue(scip) );
   SCIP_CALL( SCIPincludeExprhdlrSum(scip) );
   SCIP_CALL( SCIPincludeExprhdlrPow(scip) );
   SCIP_CALL( SCIPincludeExprhdlrProduct(scip) );
   SCIP_CALL( SCIPincludeExprhdlrSin(scip) );
   SCIP_CALL( SCIPincludeExprhdlrCos(scip) );
   SCIP_CALL( SCIPincludeExprhdlrCos(scip) );
   SCIP_CALL( SCIPincludeExprhdlrVar(scip) );
   SCIP_CALL( SCIPincludeExprhdlrVaridx(scip) );

   /* if no IPOPT available, don't run test */
   if( ! SCIPisIpoptAvailableIpopt() )
      return;

   /* include NLPI's */
   SCIPincludeNlpSolverIpopt(scip);

   /* include gauge separator */
   SCIP_CALL( SCIPincludeSepaGauge(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* change SCIP's stage to be able to create nlrows and rows; ask SCIP to generate an NLP */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   /* create and add variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -2.5, 3.1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -2.5, 3.1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create nlrows, add them to SCIP and release them */
   createNlRow1(TRUE);
   createNlRow2(TRUE);
   createNlRow3(TRUE);

   SCIP_CALL( SCIPaddNlRow(scip, nlrow1) );
   SCIP_CALL( SCIPaddNlRow(scip, nlrow2) );
   SCIP_CALL( SCIPaddNlRow(scip, nlrow3) );

   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow1) );
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow2) );
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow3) );

   /* get sepadata */
   sepadata = SCIPsepaGetData(SCIPfindSepa(scip, "gauge"));

   /* compute interior point */
   SCIP_CALL( computeInteriorPoint(scip, sepadata) );

   /* check sepadata stuff changed in call */
   cr_assert_not(sepadata->skipsepa);
   cr_assert(sepadata->isintsolavailable);
   cr_assert_not_null(sepadata->intsol);

   /* check interior solution */
   cr_expect_float_eq(-0.275971168224138, SCIPgetSolVal(scip, sepadata->intsol, x), 1e-5, "received %.10f instead", SCIPgetSolVal(scip, sepadata->intsol, x));
   cr_expect_float_eq( 0.318323856389092, SCIPgetSolVal(scip, sepadata->intsol, y), 1e-5, "received %g instead", SCIPgetSolVal(scip, sepadata->intsol, y));

   /* free memory */
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}
