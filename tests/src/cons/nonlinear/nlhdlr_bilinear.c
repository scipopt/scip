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

/**@file   nlhdlr_bilinear.c
 * @brief  tests bilinear nonlinear handler methods
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/scipdefplugins.h"
#include "scip/scip.h"
#include "scip/nlhdlr_bilinear.c"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_NLHDLR* nlhdlr;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* disable heuristics so problems don't get solved; presolving and propagation so problem doesn't get modified */
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert_not_null(conshdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, "bilinear");
   cr_assert_not_null(nlhdlr);
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

/* define the test suite */
TestSuite(nlhdlrbilinear, .init = setup, .fini = teardown);

/** creates a constraint and calls detection methods of nonlinear handlers */
static
SCIP_RETCODE createAndDetect(
   SCIP_CONS**           cons,              /**< pointer to store the constraint */
   SCIP_EXPR*            rootexpr           /**< root expression of the constraint */
   )
{
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, cons, "cons", rootexpr, -100.0, 100.0) );

   /* detect gets called when adding a cons in solving stage */
   SCIP_CALL( SCIPaddCons(scip, *cons) );

   return SCIP_OKAY;
}

/*
 * tests
 */

/* creates and adds two nonlinear constraints and check output of SCIPgetBilinearExprsExpr */
Test(nlhdlrbilinear, collect_product_expressions)
{
   SCIP_CONS* conss[2];
   SCIP_EXPR** exprs;
   SCIP_EXPR* expr;
   const char* inputs[2] = {"<x> * <y> + <y> * <z>", "exp(<x>) * sin(<y>)"};
   SCIP_Bool infeasible;
   int nexprs;
   int c;

   /* create constraints */
   for( c = 0; c < 2; ++c )
   {
      /* when parsing we do not need to pass our own owner data, since cons nonlinear is going to own the expression in
       * SCIPcreateConsBasicNonlinear
       */
      SCIP_CALL( SCIPparseExpr(scip, &expr, inputs[c], NULL, NULL, NULL) );
      SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &conss[c], "cons", expr, -100.0, 100.0) );
      /* add and release constraint */
      SCIP_CALL( SCIPaddCons(scip, conss[c]) );
      SCIP_CALL( SCIPreleaseCons(scip, &conss[c]) );
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }

   /* transform problem, initialize solve, initlp */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   SCIP_CALL( SCIPconstructLP(scip, &infeasible) );

   /* check whether all product expressions could be found */
   exprs = SCIPgetExprsBilinear(nlhdlr);
   cr_expect(exprs != NULL);
   nexprs = SCIPgetNExprsBilinear(nlhdlr);
   cr_expect(nexprs == 3);

   for( c = 0; c < nexprs; ++c )
   {
      SCIP_EXPR* child1;
      SCIP_EXPR* child2;

      cr_assert(exprs[c] != NULL);
      cr_expect(SCIPexprGetNChildren(exprs[c]) == 2);
      cr_expect(SCIPgetExprAuxVarNonlinear(exprs[c]) != NULL);
      cr_expect(SCIPisExprProduct(scip, exprs[c]));

      child1 = SCIPexprGetChildren(exprs[c])[0];
      cr_assert(child1 != NULL);
      cr_expect(SCIPgetExprAuxVarNonlinear(child1) != NULL);
      child2 = SCIPexprGetChildren(exprs[c])[1];
      cr_assert(child2 != NULL);
      cr_expect(SCIPgetExprAuxVarNonlinear(child2) != NULL);

      /* check for x*y */
      if( SCIPgetExprAuxVarNonlinear(child1) == SCIPvarGetTransVar(x) )
         cr_expect(SCIPgetExprAuxVarNonlinear(child2) == SCIPvarGetTransVar(y));

      /* check for y*z */
      else if( SCIPgetExprAuxVarNonlinear(child1) == SCIPvarGetTransVar(y) )
         cr_expect(SCIPgetExprAuxVarNonlinear(child2) == SCIPvarGetTransVar(z));

      /* check for exp(x) * sin(y) */
      else
      {
         cr_expect(SCIPisExprExp(scip, child1));
         cr_expect(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(child2)), "sin") == 0);
      }
   }
}

/* adds valid inequalities to a product expression */
Test(nlhdlrbilinear, add_inequality)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 4.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x> * <t_y>", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   /*
    * add inequalities to product expression
    */

   /* separating inequality */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.0, -2.0, 4.5, &success) );
   cr_expect(success);

   /* non-separating inequality */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.0, -1.0, 5.0, &success) );
   cr_expect(!success);

   /* duplicates should not be accepted */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.0, -2.0, 4.5, &success) );
   cr_expect(!success);

   /* first inequality is still dominating this inequality -> do not accept it */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 0.9, -3.0, 5.5, &success) );
   cr_expect(!success);

   /* separates bottom-left corner */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, -1.0, -1.1, -2.0, &success) );
   cr_expect(success);

   /* dominates first inequality */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.1, 0.5, -3.0, &success) );
   cr_expect(success);

   /* worse than inequality for bottom-left corner */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, -1.0, -1.1, -2.2, &success) );
   cr_expect(!success);

   /* free memory */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/** auxiliary function to check coefficients of a cut */
static
SCIP_RETCODE checkCut(
   SCIP_ROW*            row,                /**< row */
   SCIP_Real            coefx,              /**< target coefficient for x */
   SCIP_Real            coefy,              /**< target coefficient for y */
   SCIP_Real            coefz,              /**< target coefficient for the auxiliary variable */
   SCIP_Real            constant            /**< target constant */
   )
{
   SCIP_VAR* tx = SCIPvarGetTransVar(x);
   SCIP_VAR* ty = SCIPvarGetTransVar(y);
   int i;

   assert(row != NULL);

   for( i = 0; i < SCIProwGetNNonz(row); ++i )
   {
      SCIP_VAR* var = SCIPcolGetVar(SCIProwGetCols(row)[i]);
      SCIP_Real coef = SCIProwGetVals(row)[i];

      cr_assert(var != NULL);

      if( var == tx )
         assert(SCIPisEQ(scip, coefx, coef));
      else if( var == ty )
         cr_expect(SCIPisEQ(scip, coefy, coef));
      else
         cr_expect(SCIPisEQ(scip, coefz, coef));

      if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
         cr_expect(SCIPisEQ(scip, -SCIProwGetLhs(row), constant));
      else
         cr_expect(SCIPisEQ(scip, -SCIProwGetRhs(row), constant));
   }

   return SCIP_OKAY;
}

static
SCIP_INTERVAL computeAndGetActivity(
   SCIP_EXPR*   expr                /**< expression where to free enfo data */
   )
{
   /* force recomputation of activity (SCIPincrementCurBoundsTagNonlinear would normally be called by prop_obbt) */
   SCIPincrementCurBoundsTagNonlinear(conshdlr, FALSE);

   /* this calls the conshdlr activity computation via the owner data inside expr */
   SCIP_CALL( SCIPevalExprActivity(scip, expr) );

   return SCIPexprGetActivity(expr);
}

/* adds a linear inequality to the product expression and computes a cut at a given reference point */
Test(nlhdlrbilinear, separation_single)
{
   SCIP_ROWPREP* rowprep;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_ROW* row;
   SCIP_Bool success;
   SCIP_Bool overestimate;
   SCIP_Bool dummy;
   SCIP_EXPR* expr;
   SCIP_PTRARRAY* rowpreps;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x> * <t_y>", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   SCIP_CALL( SCIPconstructLP(scip, &dummy) );

   /* INITLP should have added an auxiliary variable to the product expression (tight might change in the future) */
   cr_assert( SCIPgetExprAuxVarNonlinear(expr) != NULL);

   /* create a solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x), 1.1214285714285714) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(y), -0.1785714285714286) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr), -1.25) );

    /* add inequality for underestimation to product expression */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.0, 1.0, 1.3, &success) );
   cr_expect(success);

   /* auxvar = -1.25 => McCormick relaxation is violated so the resulting cut should be one of the McCormick inequalities */
   overestimate = FALSE;
   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );
   SCIP_CALL( nlhdlrEstimateBilinear(scip, conshdlr, nlhdlr, expr, SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr), sol,
            0.0, overestimate, 0.0, FALSE, rowpreps, &success, &dummy) );
   cr_expect(!success);
   cr_expect(SCIPgetPtrarrayMinIdx(scip, rowpreps) > SCIPgetPtrarrayMaxIdx(scip, rowpreps));
   SCIP_CALL( SCIPclearPtrarray(scip, rowpreps) );

   /* auxvar = -1.23 => point is in the convex hull of McCormick and thus we should use the bilinear inequality */
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr), -1.23) );
   overestimate = FALSE;
   SCIP_CALL( nlhdlrEstimateBilinear(scip, conshdlr, nlhdlr, expr, SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr), sol,
            0.0, overestimate, 0.0, FALSE, rowpreps, &success, &dummy) );
   cr_expect(success);
   cr_expect(SCIPgetPtrarrayMinIdx(scip, rowpreps) == 0);
   cr_expect(SCIPgetPtrarrayMaxIdx(scip, rowpreps) == 0);
   rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, 0);
   SCIP_CALL( SCIPgetRowprepRowConshdlr(scip, &row, rowprep, conshdlr) );
   SCIP_CALL( checkCut(row, 0.5790816326530611, 0.3637755102040817, -1.0, -0.7846938775510204) );
   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   SCIPfreeRowprep(scip, &rowprep);

   /* free memory */
   SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
}

/* adds two linear inequality to the product expression and computes a cut at a given reference point */
Test(nlhdlrbilinear, separation_two)
{
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_ROWPREP* rowprep;
   SCIP_ROW* row;
   SCIP_Bool success;
   SCIP_Bool overestimate;
   SCIP_Bool dummy;
   SCIP_EXPR* expr;
   SCIP_PTRARRAY* rowpreps;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 3.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x> * <t_y>", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   SCIP_CALL( SCIPconstructLP(scip, &dummy) );

   /* INITLP should have added an auxiliary variable to the product expression (tight might change in the future) */
   cr_assert( SCIPgetExprAuxVarNonlinear(expr) != NULL);

   /* create a solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x), 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(y), 1.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr), 4.4) );

   /* add two inequalities to the product expression */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, -1.2, -0.6, 0.9, &success) );
   cr_expect(success);
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 0.5, 0.3, 0.6, &success) );
   cr_expect(success);

   /* auxvar = -0.1 => McCormick should do the job */
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr), -0.1) );
   overestimate = FALSE;
   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );
   SCIP_CALL( nlhdlrEstimateBilinear(scip, conshdlr, nlhdlr, expr, SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr), sol, 0.0, overestimate, 0.0, FALSE, rowpreps, &success, &dummy) );
   cr_expect(!success);
   cr_expect(SCIPgetPtrarrayMinIdx(scip, rowpreps) > SCIPgetPtrarrayMaxIdx(scip, rowpreps));
   SCIP_CALL( SCIPclearPtrarray(scip, rowpreps) );

   /* auxvar = 0.1 => both inequalities needs to be used */
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr), 0.1) );
   overestimate = FALSE;
   SCIP_CALL( nlhdlrEstimateBilinear(scip, conshdlr, nlhdlr, expr, SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr), sol, 0.0, overestimate, 0.0, FALSE, rowpreps, &success, &dummy) );
   cr_expect(success);
   cr_expect(SCIPgetPtrarrayMinIdx(scip, rowpreps) == 0);
   cr_expect(SCIPgetPtrarrayMaxIdx(scip, rowpreps) == 0);
   rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, 0);
   SCIP_CALL( SCIPgetRowprepRowConshdlr(scip, &row, rowprep, conshdlr) );
   SCIP_CALL( checkCut(row, 1.45445115010332, 0.9772255750516624, -1.0, -1.9213268615442665) );
   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   SCIPfreeRowprep(scip, &rowprep);

   /* free memory */
   SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
}

/* tests the interval evaluation callback when linear inequalities are available */
Test(nlhdlrbilinear, inteval_corner)
{
   SCIP_CONS* cons;
   SCIP_EXPR* expr;
   SCIP_INTERVAL interval;
   SCIP_Bool success;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -3.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 5.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x> * <t_y>", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   /* INITSOL should have added a request for an auxiliary variable to the product expression (tight might change in the future) */
   cr_assert( SCIPgetExprNAuxvarUsesNonlinear(expr) > 0);

   /* add linear inequality x <= -y + 7 => maximum of xy is attained at (3,4) */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.0, -1.0, 7.0, &success) );
   cr_expect(success);

   interval = computeAndGetActivity(expr);

   /* compute interval */
   cr_expect(SCIPisEQ(scip, SCIPintervalGetInf(interval), -10.0), "expect -10.0 got %g\n", SCIPintervalGetInf(interval));
   cr_expect(SCIPisEQ(scip, SCIPintervalGetSup(interval), 12.0), "expect 12.0 got %g\n", SCIPintervalGetSup(interval));

   /* free memory */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* tests the interval evaluation callback when linear inequalities are available */
Test(nlhdlrbilinear, inteval_single_line)
{
   SCIP_CONS* cons;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -3.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 5.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x> * <t_y>", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   /* add linear inequality -x <= -y + 3.5 => maximum of xy is attained at (-1.75,1.75) */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, -1.0, -1.0, 3.5, &success) );
   cr_expect(success);

   interval = computeAndGetActivity(expr);

   /* compute interval */
   cr_expect(SCIPisEQ(scip, SCIPintervalGetInf(interval), -3.0625), "expect -3.0625 got %g\n", SCIPintervalGetInf(interval));
   cr_expect(SCIPisEQ(scip, SCIPintervalGetSup(interval), 15.0), "expect 15.0 got %g\n", SCIPintervalGetSup(interval));

   /* free memory */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* tests the interval evaluation callback when three intersecting linear inequalities are available */
Test(nlhdlrbilinear, inteval_three_lines)
{
   SCIP_EXPR* varexprs[2];
   SCIP_CONS* cons;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -4.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -2.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -3.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 5.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create nonlinear constraint */
   SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[0], SCIPvarGetTransVar(x), NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[1], SCIPvarGetTransVar(y), NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr, 2, varexprs, -2.0, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );

   /* simplify replaced expr by a sum of 1 term with coef -2, and term being the product we want */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   expr = SCIPexprGetChildren(SCIPgetExprNonlinear(cons))[0];
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* add -x <= -y + 3.5 and x <= -0.6 y -2.25 and -x <= 0.1 y + 3.5 */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, -1.0, -1.0, 3.5, &success) );
   cr_expect(success);

   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.0, -0.6, -2.25, &success) );
   cr_expect(success);

   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, -1.0, 0.1, 3.5, &success) );
   cr_expect(success);

   interval = computeAndGetActivity(SCIPgetExprNonlinear(cons));

   /* compute interval */
   cr_expect(SCIPisEQ(scip, SCIPintervalGetInf(interval), -19.2), "expect -19.2 got %g\n", SCIPintervalGetInf(interval));
   cr_expect(SCIPisEQ(scip, SCIPintervalGetSup(interval), 4.248046875), "expect 4.248046875 got %g\n", SCIPintervalGetSup(interval));

   /* free memory */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[0]) );
}

/* tests the interval evaluation callback when linear inequalities are available */
Test(nlhdlrbilinear, inteval_parallel)
{
   SCIP_CONS* cons;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x> * <t_y>", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   /* add inequalities with the same slope: x <= y + 1 and x >= y - 1 */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.0, 1.0, 1.0, &success) );
   cr_expect(success);
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, -1.0, -1.0, 1.0, &success) );
   cr_expect(success);

   interval = computeAndGetActivity(expr);

   /* compute interval */
   cr_expect(SCIPisEQ(scip, SCIPintervalGetInf(interval), -0.25), "expect -0.25 got %g\n", SCIPintervalGetInf(interval));
   cr_expect(SCIPisEQ(scip, SCIPintervalGetSup(interval), 1.0), "expect 1.0 got %g\n", SCIPintervalGetSup(interval));

   /* free memory */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* tests the reverse propagation callback when linear inequalities are available */
Test(nlhdlrbilinear, reverseprop_single)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata;
   SCIP_CONS* cons;
   SCIP_INTERVAL activity;
   SCIP_INTERVAL intervalx;
   SCIP_INTERVAL intervaly;
   SCIP_Bool success;
   SCIP_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x> * <t_y>", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   /* add inequality x <= y - 1 */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.0, 1.0, -1.0, &success) );
   cr_expect(success);

   activity = computeAndGetActivity(expr);

   cr_expect(SCIPisEQ(scip, activity.inf, -1.0));
   cr_expect(SCIPisEQ(scip, activity.sup, 0.0));

   /* compute intervals; inequality cuts off two corner points */
   nlhdlrexprdata = SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr);

   reversePropBilinear(scip, conshdlr, expr, activity, nlhdlrexprdata->underineqs, nlhdlrexprdata->nunderineqs, nlhdlrexprdata->overineqs, nlhdlrexprdata->noverineqs, &intervalx, &intervaly);
   cr_expect(SCIPisEQ(scip, intervalx.inf, -1.0));
   cr_expect(SCIPisEQ(scip, intervalx.sup, 0.0));
   cr_expect(SCIPisEQ(scip, intervaly.inf, 0.0));
   cr_expect(SCIPisEQ(scip, intervaly.sup, 1.0));

   /* free memory */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* tests the reverse propagation callback when a linear inequality together with the level set of a bilinear term
 * results in stronger variable bounds
 */
Test(nlhdlrbilinear, reverseprop_levelset)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata;
   SCIP_CONS* cons;
   SCIP_INTERVAL intervalx;
   SCIP_INTERVAL intervaly;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_Bool cutoff = FALSE;
   int ntightenings = 0;
   SCIP_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x> * <t_y>", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   /* add inequality x <= 0.7 y + 0.1 */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, 1.0, 0.7, 0.1, &success) );
   cr_expect(success);

   interval = computeAndGetActivity(expr);

   cr_expect(SCIPisEQ(scip, interval.inf, -1.0));
   cr_expect(SCIPisEQ(scip, interval.sup, 1.0));

   /* tighten the bounds of the expression */
   SCIPintervalSetBounds(&interval, -0.8, 0.5);
   SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, expr, interval, &cutoff, &ntightenings) );
   interval = SCIPgetExprBoundsNonlinear(scip, expr);
   cr_expect(!cutoff);
   cr_expect(SCIPisEQ(scip, interval.inf, -0.8));
   cr_expect(SCIPisEQ(scip, interval.sup, 0.5));

   /* compute intervals; inequality cuts off two corner points */
   nlhdlrexprdata = SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr);
   reversePropBilinear(scip, conshdlr, expr, interval, nlhdlrexprdata->underineqs, nlhdlrexprdata->nunderineqs, nlhdlrexprdata->overineqs, nlhdlrexprdata->noverineqs, &intervalx, &intervaly);
   cr_expect(SCIPisEQ(scip, intervalx.inf, -1.0));
   cr_expect_float_eq(intervalx.sup, 0.64371709509590, 1e-8);
   cr_expect(intervalx.sup >= 0.64371709509590); /* intervalarith.c gives us safe bounds */
   cr_expect_float_eq(intervaly.inf, -0.91959585365193, 1e-8);
   cr_expect(intervaly.inf <= -0.91959585365193); /* intervalarith.c gives us safe bounds */
   cr_expect(SCIPisEQ(scip, intervaly.sup, 1.0));

   /* free memory */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* tests the reverse propagation callback when a linear inequality does not intersect the level set */
Test(nlhdlrbilinear, reverseprop_levelset_nointersection)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata;
   SCIP_CONS* cons;
   SCIP_INTERVAL intervalx;
   SCIP_INTERVAL intervaly;
   SCIP_INTERVAL interval;
   SCIP_Bool success;
   SCIP_EXPR* expr;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<t_x> * <t_y>", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   /* add inequality -x <= y */
   SCIP_CALL( SCIPaddIneqBilinear(scip, nlhdlr, expr, -1.0, 1.0, 0.0, &success) );
   cr_expect(success);

   interval = computeAndGetActivity(expr);

   cr_expect(SCIPisEQ(scip, interval.inf, -1.0));
   cr_expect(SCIPisEQ(scip, interval.sup, 1.0));

   /* compute intervals; inequality does not intersect the level set in the interior of the domain -> no tightening */
   nlhdlrexprdata = SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr);
   SCIPintervalSetBounds(&interval, -1.0, 0.5);
   reversePropBilinear(scip, conshdlr, expr, interval, nlhdlrexprdata->underineqs, nlhdlrexprdata->nunderineqs, nlhdlrexprdata->overineqs, nlhdlrexprdata->noverineqs, &intervalx, &intervaly);
   cr_expect(SCIPisEQ(scip, intervalx.inf, -1.0));
   cr_expect(SCIPisEQ(scip, intervalx.sup, 1.0));
   cr_expect(SCIPisEQ(scip, intervaly.inf, -1.0));
   cr_expect(SCIPisEQ(scip, intervaly.sup, 1.0));

   /* free memory */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
