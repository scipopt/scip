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

/**@file   gauge.c
 * @brief  unit test for cons_quadratic methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scip/nodesel_bfs.h"
#include "scip/heur_multistart.c"
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

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_HEUR* heursubnlp;
static SCIP_HEUR* heurmultistart;
static SCIP_HASHMAP* varindex;
static SCIP_SOL* sol;
static SCIP_RANDNUMGEN* randumgen;

#define EPS    1e-5

static
void setup(void)
{
   SCIP_VAR* origx;
   SCIP_VAR* origy;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &origx, "x", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &origy, "y", -10.0, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, origx) );
   SCIP_CALL( SCIPaddVar(scip, origy) );

   SCIP_CALL( SCIPincludeHeurMultistart(scip) );
   SCIP_CALL( SCIPincludeHeurSubNlp(scip) );
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

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   x = SCIPvarGetTransVar(origx);
   y = SCIPvarGetTransVar(origy);

   /* create mapping between variables and 0,..,SCIPgetNVars(scip)-1 */
   SCIP_CALL( SCIPhashmapCreate(&varindex, SCIPblkmem(scip), 2) );
   SCIP_CALL( SCIPhashmapInsertInt(varindex, (void*)x, 0) );
   SCIP_CALL( SCIPhashmapInsertInt(varindex, (void*)y, 1) );

   SCIP_CALL( SCIPreleaseVar(scip, &origx) );
   SCIP_CALL( SCIPreleaseVar(scip, &origy) );

   /* create useless variables to get around assert(SCIPgetNTotalVars(scip) >= *nvarexprs + 1); in SCIPgetExprVarExprs */
   SCIP_CALL( SCIPcreateVarBasic(scip, &origx, "bla", -10.0, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPreleaseVar(scip, &origx) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &origx, "ble", -10.0, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPreleaseVar(scip, &origx) );

   /* solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   heurmultistart = SCIPfindHeur(scip, "multistart");
   heursubnlp = SCIPfindHeur(scip, "subnlp");
   cr_assert( heurmultistart != NULL );
   cr_assert( heursubnlp != NULL );

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randumgen, 777, TRUE) );
}

static
void teardown(void)
{
   SCIPfreeRandom(scip, &randumgen);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIPhashmapFree(&varindex);
   SCIP_CALL( SCIPfree(&scip) );
}


Test(heuristic, sampleRandomPoints, .init = setup, .fini = teardown,
   .description = "check sampleRandomPoints() subroutine of the multi-start heuristic"
   )
{
   SCIP_SOL** rndpoints;
   int nrndpoints;

   SCIP_CALL( SCIPallocBufferArray(scip, &rndpoints, 1) );

   /* compute a single random point */
   rndpoints[0] = NULL;
   SCIP_CALL( sampleRandomPoints(scip, rndpoints, 1, 1.0, randumgen, SCIPinfinity(scip), &nrndpoints) );
   cr_assert( nrndpoints == 1 );

   cr_assert( rndpoints[0] != NULL );
   cr_assert( SCIPgetSolVal(scip, rndpoints[0], x) <= SCIPvarGetUbLocal(x) );
   cr_assert( SCIPgetSolVal(scip, rndpoints[0], x) >= SCIPvarGetLbLocal(x) );
   cr_assert( SCIPgetSolVal(scip, rndpoints[0], y) <= SCIPvarGetUbLocal(y) );
   cr_assert( SCIPgetSolVal(scip, rndpoints[0], y) >= SCIPvarGetLbLocal(y) );

   SCIP_CALL( SCIPfreeSol(scip, &rndpoints[0]) );

   SCIPfreeBufferArray(scip, &rndpoints);
}

Test(heuristic, computeGradient, .init = setup, .fini = teardown,
   .description = "check computeGradient subroutine of the multi-start heuristic"
   )
{
   SCIP_NLROW* nlrow;
   SCIP_VAR* linvars[2];
   SCIP_Real lincoefs[2];
   SCIP_Real grad[2];
   SCIP_EXPR* expr;
   SCIP_Real norm;
   SCIP_EXPRITER* exprit;

   linvars[0] = x;
   linvars[1] = y;
   lincoefs[0] = 2.3;
   lincoefs[1] = -3.1;

   SCIP_CALL( SCIPcreateExpriter(scip, &exprit) );

   /* compute the gradient for 2.3*x - 3.1*y at point (x,y) = (-2,3) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 3.0) );
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, "nlrow", 5.0, 2, linvars, lincoefs, NULL, 1.0, 1.0,
         SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( computeGradient(scip, nlrow, sol, varindex, exprit, grad, &norm) );
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );

   cr_assert( SCIPisEQ(scip, grad[0], 2.3) );
   cr_assert( SCIPisEQ(scip, grad[1], -3.1) );

   /* compute the gradient for 2.3*x - 3.1*y + 2*x^2 -4*y^2 + 5xy at point (x,y) = (1,-7) */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"2*<t_x>^2 + 5.0*<t_x>*<t_y> -4*<t_y>^2", NULL, NULL, NULL) );

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -7.0) );
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, "nlrow", 5.0, 2, linvars, lincoefs, expr, 1.0, 1.0,
         SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( computeGradient(scip, nlrow, sol, varindex, exprit, grad, &norm) );
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );

   cr_assert( SCIPisEQ(scip, grad[0], 2.3 + 4 * 1 + 5 * (-7)), "expecting %g, got %g\n", 2.3 + 4 * 1 + 5 * (-7), grad[0]);
   cr_assert( SCIPisEQ(scip, grad[1], -3.1 - 8 * (-7) + 5 * 1) );

   /* create expression tree for 2.3*x - 3.1*y + 2*x^2 -4*y^2 + 5xy + x*e^y at point (x,y) = (3,3) */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"2*<t_x>^2 -4*<t_y>^2 + 5*<t_x>*<t_y> + <t_x>*exp(<t_y>)", NULL, NULL, NULL) );

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 3.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 3.0) );

   /* create expression tree for 2.3*x - 3.1*y + 2*x^2 -4*y^2 + 5xy + x*e^y at point (x,y) = (3,3) */
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, "nlrow", 5.0, 2, linvars, lincoefs, expr, 1.0, 1.0,
         SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( computeGradient(scip, nlrow, sol, varindex, exprit, grad, &norm) );

   cr_assert( SCIPisEQ(scip, grad[0], 2.3 + 4 * 3 + 5 * 3 + exp(3)) );
   cr_assert( SCIPisEQ(scip, grad[1], -3.1 - 8 * 3 + 5 * 3 + 3*exp(3)) );

   SCIPfreeExpriter(&exprit);
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(heuristic, improvePoint, .init = setup, .fini = teardown,
   .description = "check improvePoint subroutine of the multi-start heuristic"
   )
{
   SCIP_VAR* vars[2];
   SCIP_Real lincoefs[2];
   SCIP_NLROW* nlrows[2];
   SCIP_Real nlrowgradcosts[2];
   SCIP_Real gradcosts;
   SCIP_Real minfeas;
   SCIP_EXPR* expr;

   nlrowgradcosts[0] = 1.0;
   nlrowgradcosts[1] = 1.0;

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 5.0) );

   vars[0] = x;
   vars[1] = y;

   /* for one linear constraint the method should stop after one iteration (if the projection on the linear constraint
    * is not outside the domain)
    */
   lincoefs[0] = 5.0;
   lincoefs[1] = 1.0;
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[0], "nlrow1", 0.0, 2, vars, lincoefs, NULL, 1.0, 1.0,
         SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( improvePoint(scip, nlrows, 1, varindex, sol, 1, 0.0, INT_MAX, &minfeas, nlrowgradcosts, &gradcosts) );
   cr_expect(SCIPisFeasEQ(scip, minfeas, 0.0));
   cr_expect(gradcosts > 0.0);

   /* consider one quadratic constraint of the form sqrt(x^2 + y^2) = 0.5 */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<t_x>^2 + <t_y>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrows[1], "nlrow2", 0.0, 0, NULL, NULL, expr, 0.25, 0.25,
         SCIP_EXPRCURV_UNKNOWN) );

   /* start inside the ball */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 0.2) );
   SCIP_CALL( improvePoint(scip, &nlrows[1], 1, varindex, sol, 10, 0.0, INT_MAX, &minfeas, nlrowgradcosts, &gradcosts) );
   cr_expect(SCIPisFeasGE(scip, minfeas, -EPS), "expecting minfeas %g > -eps (%g)\n", minfeas, -EPS);
   cr_expect(gradcosts > 0.0);

   /* start outside the ball */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 5.0) );
   SCIP_CALL( improvePoint(scip, &nlrows[1], 1, varindex, sol, 100, 0.0, INT_MAX, &minfeas, nlrowgradcosts, &gradcosts) );
   cr_expect(SCIPisFeasGE(scip, minfeas, -EPS), "expecting minfeas %g > -eps (%g)\n", minfeas, -EPS);
   cr_expect(gradcosts > 0.0);

   /* consider linear and quadratic constraint */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -10.0) );
   SCIP_CALL( improvePoint(scip, nlrows, 2, varindex, sol, 100, 0.0, INT_MAX, &minfeas, nlrowgradcosts, &gradcosts) );
   cr_expect(SCIPisFeasGE(scip, minfeas, -EPS), "expecting minfeas %g > -eps (%g)\n", minfeas, -EPS);
   cr_expect(gradcosts > 0.0);

   /* free nlrows */
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrows[1]) );
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrows[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(heuristic, filterPoints, .init = setup, .fini = teardown,
   .description = "check filterPoints subroutine of the multi-start heuristic"
   )
{
   SCIP_SOL* points[10];
   SCIP_Real feasibilities[10];
   int nusefulpoints;
   int i;

   for( i = 0; i < 10; ++i )
      feasibilities[i] = i - 10;

   /* mean violation of the points should be -6.47 which means that the points with a violation of -7, -8. -9, or -10
    * are excluded
    */
   SCIP_CALL( filterPoints(scip, points, feasibilities, 10, &nusefulpoints) );
   cr_assert(nusefulpoints == 6);

   /* check is points are ordered decreasingly w.r.t their feasibilities */
   for( i = 0; i < 9; ++i )
      cr_assert(feasibilities[i] >= feasibilities[i+1]);
}

#define NPOINTS 100
Test(heuristic, clusterPointsGreedy, .init = setup, .fini = teardown,
   .description = "check clusterPointsGreedy subroutine of the multi-start heuristic"
   )
{
   SCIP_SOL* points[NPOINTS];
   SCIP_Real maxreldist;
   int clusteridx[NPOINTS];
   int nclusters;
   int i;

   maxreldist = 0.05;

   /* create solutions with random values */
   for( i = 0; i < NPOINTS; ++i )
   {
      SCIP_CALL( SCIPcreateSol(scip, &points[i], NULL) );
      SCIP_CALL( SCIPsetSolVal(scip, points[i], x, SCIPrandomGetReal(randumgen, -1.0, 1.0)) );
      SCIP_CALL( SCIPsetSolVal(scip, points[i], y, SCIPrandomGetReal(randumgen, -10.0, 10.0)) );
   }

   SCIP_CALL( clusterPointsGreedy(scip, points, NPOINTS, clusteridx, &nclusters, 1e+04, maxreldist, INT_MAX) );
   cr_assert(nclusters > 0 && nclusters <= NPOINTS);

   for( i = 0; i < NPOINTS; ++i )
   {
      int j;

      cr_assert(clusteridx[i] >= 0);
      cr_assert(clusteridx[i] < nclusters);

      for( j = i + 1; j < NPOINTS; ++j )
      {
         SCIP_Real dist = getRelDistance(scip, points[i], points[j], 1e+04);

         /* distance between i and j in the same cluster should not be twice as large as maxreldist */
         cr_assert(clusteridx[i] != clusteridx[j] || dist <= 2 * maxreldist);
      }
   }

   /* free solutions */
   for( i = NPOINTS - 1; i >= 0; --i )
   {
      SCIP_CALL( SCIPfreeSol(scip, &points[i]) );
   }
}
