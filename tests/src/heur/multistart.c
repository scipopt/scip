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

/**@file   gauge.c
 * @brief  unit test for cons_quadratic methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scip/nodesel_bfs.h"
#include "scip/heur_multistart.c"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_HEUR* heursubnlp;
static SCIP_HEUR* heurmultistart;

static
void setup(void)
{
   SCIP_VAR* origx;
   SCIP_VAR* origy;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &origx, "x", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &origy, "y", 0.0, 100.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, origx) );
   SCIP_CALL( SCIPaddVar(scip, origy) );

   SCIP_CALL( SCIPincludeHeurMultistart(scip) );
   SCIP_CALL( SCIPincludeHeurSubNlp(scip) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING) );

   x = SCIPvarGetTransVar(origx);
   y = SCIPvarGetTransVar(origy);

   SCIP_CALL( SCIPreleaseVar(scip, &origx) );
   SCIP_CALL( SCIPreleaseVar(scip, &origy) );

   heurmultistart = SCIPfindHeur(scip, "multistart");
   heursubnlp = SCIPfindHeur(scip, "subnlp");
   cr_assert( heurmultistart != NULL );
   cr_assert( heursubnlp != NULL );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
}


Test(heuristic, sampleRandomPoints, .init = setup, .fini = teardown,
   .description = "check sampleRandomPoints() subroutine of the multi-start heuristic"
   )
{
   SCIP_SOL** rndpoints;
   unsigned int rndseed;

   rndseed = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &rndpoints, 1) );

   /* compute a single random point */
   rndpoints[0] = NULL;
   SCIP_CALL( sampleRandomPoints(scip, rndpoints, 1, &rndseed, 1.0) );

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
   SCIP_VAR* quadvars[2];
   SCIP_QUADELEM quadelems[3];
   SCIP_Real grad[2];
   SCIP_EXPRINT* exprint;
   SCIP_HASHMAP* varindex;
   SCIP_SOL* sol;
   SCIP_EXPRTREE* tree;
   SCIP_EXPR* expexpr;
   SCIP_EXPR* xexpr;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* expr;
   SCIP_EXPR* exprs[2];
   SCIP_Real norm;

   linvars[0] = x;
   linvars[1] = y;
   lincoefs[0] = 2.3;
   lincoefs[1] = -3.1;
   quadvars[0] = x;
   quadvars[1] = y;
   quadelems[0].idx1 = 0;
   quadelems[0].idx2 = 0;
   quadelems[0].coef = 2.0;
   quadelems[1].idx1 = 1;
   quadelems[1].idx2 = 1;
   quadelems[1].coef = -4.0;
   quadelems[2].idx1 = 0;
   quadelems[2].idx2 = 1;
   quadelems[2].coef = 5.0;

   /* initialize data structures */
   SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &exprint) );
   SCIP_CALL( SCIPhashmapCreate(&varindex, SCIPblkmem(scip), SCIPcalcHashtableSize(2)) );
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPhashmapInsert(varindex, (void*)x, (void*)(size_t)0) );
   SCIP_CALL( SCIPhashmapInsert(varindex, (void*)y, (void*)(size_t)1) );

   /* compute the gradient for 2.3*x - 3.1*y at point (x,y) = (-2,3) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 3.0) );
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, "nlrow", 5.0, 2, linvars, lincoefs, 0, NULL, 0, NULL, NULL, 1.0, 1.0) );
   SCIP_CALL( computeGradient(scip, nlrow, exprint, sol, varindex, grad, &norm) );
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );

   cr_assert( SCIPisEQ(scip, grad[0], 2.3) );
   cr_assert( SCIPisEQ(scip, grad[1], -3.1) );

   /* compute the gradient for 2*x^2 -4*y^2 at point (x,y) = (-2,3) */
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, "nlrow", 5.0, 0, NULL, NULL, 2, quadvars, 2, quadelems, NULL, 1.0, 1.0) );
   SCIP_CALL( computeGradient(scip, nlrow, exprint, sol, varindex, grad, &norm) );
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );

   cr_assert( SCIPisEQ(scip, grad[0], 4 * (-2)) );
   cr_assert( SCIPisEQ(scip, grad[1], -8 * 3) );

   /* compute the gradient for 2.3*x - 3.1*y + 2*x^2 -4*y^2 + 5xy at point (x,y) = (1,-7) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -7.0) );
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, "nlrow", 5.0, 2, linvars, lincoefs, 2, quadvars, 3, quadelems, NULL, 1.0, 1.0) );
   SCIP_CALL( computeGradient(scip, nlrow, exprint, sol, varindex, grad, &norm) );
   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );

   cr_assert( SCIPisEQ(scip, grad[0], 2.3 + 4 * 1 + 5 * (-7)) );
   cr_assert( SCIPisEQ(scip, grad[1], -3.1 - 8 * (-7) + 5 * 1) );

   /* create expression tree for 2.3*x - 3.1*y + 2*x^2 -4*y^2 + 5xy + x*e^y at point (x,y) = (3,3) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 3.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 3.0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &xexpr, SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &yexpr, SCIP_EXPR_VARIDX, 1) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expexpr, SCIP_EXPR_EXP, yexpr) );
   exprs[0] = xexpr;
   exprs[1] = expexpr;
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_PRODUCT, 2, exprs) );
   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &tree, expr, 2, 0, NULL) );
   SCIP_CALL( SCIPexprtreeSetVars(tree, 2, linvars) );
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, "nlrow", 5.0, 2, linvars, lincoefs, 2, quadvars, 3, quadelems, tree, 1.0, 1.0) );
   SCIP_CALL( computeGradient(scip, nlrow, exprint, sol, varindex, grad, &norm) );

   cr_assert( SCIPisEQ(scip, grad[0], 2.3 + 4 * 3 + 5 * 3 + exp(3)) );
   cr_assert( SCIPisEQ(scip, grad[1], -3.1 - 8 * 3 + 5 * 3 + 3*exp(3)) );

   SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   SCIP_CALL( SCIPexprtreeFree(&tree) );

   /* free memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIPhashmapFree(&varindex);
   SCIP_CALL( SCIPexprintFree(&exprint) );
}

Test(heuristic, improvePoint, .init = setup, .fini = teardown,
   .description = "check improvePoint subroutine of the multi-start heuristic"
   )
{

}

Test(heuristic, filterPoints, .init = setup, .fini = teardown,
   .description = "check filterPoints subroutine of the multi-start heuristic"
   )
{

}

Test(heuristic, clusterPointsGreedy, .init = setup, .fini = teardown,
   .description = "check clusterPointsGreedy subroutine of the multi-start heuristic"
   )
{

}
