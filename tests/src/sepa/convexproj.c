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

/**@file   convexproj.c
 * @brief  unit test for sepa_convexproj.c methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_nonlinear.h"
#include "nlpi/nlpi_ipopt.h"

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
void createNlRow1(SCIP_Bool isrhsconvex)
{
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPRCURV curvature;
   SCIP_EXPR* xexpr;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* exp_x;
   SCIP_EXPR* exp_y;
   SCIP_EXPR* sum_exp;
   SCIP_EXPR* lgexp;
   SCIP_EXPR* finalexpr;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_VAR* vars[2];

   /* setup expression for x and y */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &xexpr, SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &yexpr, SCIP_EXPR_VARIDX, 1) );

   /* setup expression for exp(x) and exp(y) */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exp_x, SCIP_EXPR_EXP, xexpr) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exp_y, SCIP_EXPR_EXP, yexpr) );

   /* setup expression for exp(x) + exp(y) */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &sum_exp, SCIP_EXPR_PLUS, exp_x, exp_y) );

   /* setup expression for lg(sum) */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &lgexp, SCIP_EXPR_LOG, sum_exp) );

   /* decide curvature */
   if( isrhsconvex )
   {
      curvature = SCIP_EXPRCURV_CONVEX;
      lhs = -SCIPinfinity(scip);
      rhs = 1.0;
      finalexpr = lgexp;
   }
   else
   {
      SCIP_EXPR* minusone;
      curvature = SCIP_EXPRCURV_CONCAVE;
      lhs = -1.0;
      rhs = SCIPinfinity(scip);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &minusone, SCIP_EXPR_CONST, -1.0) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &finalexpr, SCIP_EXPR_MUL, minusone, lgexp) );
   }

   /* setup expression tree with finalexpr as root expression */
   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, finalexpr, 2, 0, NULL) );

   /* assign SCIP variables to tree */
   vars[0] = x;
   vars[1] = y;
   SCIP_CALL( SCIPexprtreeSetVars(exprtree, 2, vars) );

   /* create nonlinear row for exprtree <= 1 */
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow1, "lg(sum exp)", 0.0, 0, NULL, NULL, 0, NULL, 0, NULL, exprtree, lhs, rhs, curvature) );

   /* release */
   SCIP_CALL( SCIPexprtreeFree(&exprtree) );
}

/* create x^2 <= y */
static
void createNlRow2(SCIP_Bool isrhsconvex)
{
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPRCURV curvature;
   SCIP_EXPR* xexpr;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* sqr_x;
   SCIP_EXPR* diff;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_VAR* vars[2];

   /* setup expression for x and y */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &xexpr, SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &yexpr, SCIP_EXPR_VARIDX, 1) );

   /* setup expression for x^2 */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &sqr_x, SCIP_EXPR_SQUARE, xexpr) );

   /* decide curvature */
   if( isrhsconvex )
   {
      curvature = SCIP_EXPRCURV_CONVEX;
      lhs = -SCIPinfinity(scip);
      rhs = 0.0;
      /* setup expression for x^2 - y */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &diff, SCIP_EXPR_MINUS, sqr_x, yexpr) );
   }
   else
   {
      curvature = SCIP_EXPRCURV_CONCAVE;
      lhs = 0.0;
      rhs = SCIPinfinity(scip);
      /* setup expression for - x^2 + y */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &diff, SCIP_EXPR_MINUS, yexpr, sqr_x) );
   }

   /* setup expression trees with exprs as root expression */
   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, diff, 2, 0, NULL) );

   /* assign SCIP variables to tree */
   vars[0] = x;
   vars[1] = y;
   SCIP_CALL( SCIPexprtreeSetVars(exprtree, 2, vars) );

   /* create nonlinear row for exprtree <= 0 */
   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow2, "x^2 ? y", 0.0, 0, NULL, NULL, 0, NULL, 0, NULL, exprtree, lhs, rhs, curvature) );

   /* release */
   SCIP_CALL( SCIPexprtreeFree(&exprtree) );
}

/* 1.1*x+2.4*x^2 + 0.01*x*y + 0.3*y^2 + 0.2*log(0.5*exp(0.12*x+0.1)+2*exp(0.1*y)+0.7)  <= 0.5 */
static
void createNlRow3(SCIP_Bool isrhsconvex)
{
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPRCURV curvature;
   SCIP_EXPR* xexpr;
   SCIP_EXPR* arg_exp_x;
   SCIP_EXPR* arg_exp_y;
   SCIP_EXPR* arg_log;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* exp_x;
   SCIP_EXPR* exp_y;
   SCIP_EXPR* lgexp;
   SCIP_EXPR* finalexpr;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real side;
   SCIP_VAR* vars[2];

   /* expression for x and y */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &xexpr, SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &yexpr, SCIP_EXPR_VARIDX, 1) );

   /* expression for 0.12x+0.1 and 0.1*y */
   SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &arg_exp_x, 1, &xexpr, (SCIP_Real[]){0.12}, 0.1) );
   SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &arg_exp_y, 1, &yexpr, (SCIP_Real[]){0.1}, 0.0) );

   /* expression for exp(0.12*x+0.1) and exp(0.1*y) */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exp_x, SCIP_EXPR_EXP, arg_exp_x) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exp_y, SCIP_EXPR_EXP, arg_exp_y) );

   /* expression for 0.5*exp(0.12*x+0.1) + 2*exp(0.1*y) + 0.7 */
   SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &arg_log, 2,
            (SCIP_EXPR*[]){exp_x, exp_y}, (SCIP_Real[]){0.5, 2.0}, 0.7) );

   /* expression for lg(0.5*exp(0.12*x+0.1) + 2*exp(0.1*y) + 0.7) */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &lgexp, SCIP_EXPR_LOG, arg_log) );

   /* decide curvature */
   if( isrhsconvex )
   {
      curvature = SCIP_EXPRCURV_CONVEX;
      side = 1.0;
      lhs = -SCIPinfinity(scip);
      rhs = 0.5;
   }
   else
   {
      curvature = SCIP_EXPRCURV_CONCAVE;
      side = -1.0;
      lhs = -0.5;
      rhs = SCIPinfinity(scip);
   }

   SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &finalexpr, 1,
            &lgexp, (SCIP_Real[]){side * 0.2}, 0.0) );

   /* build expression tree with log as root expression */
   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, finalexpr, 2, 0, NULL) );

   /* assign SCIP variables to tree */
   vars[0] = x;
   vars[1] = y;
   SCIP_CALL( SCIPexprtreeSetVars(exprtree, 2, vars) );

   SCIP_CALL( SCIPcreateNlRow(scip, &nlrow3, "complicated",
            0.0, 1, &x, (SCIP_Real[]){side*1.1},
            2, vars,
            3,
            (SCIP_QUADELEM[]){
            {0, 0, side*2.4}, // 2.4 * x^2
            {0, 1, side*0.01}, // 0.01*x*y
            {1, 1, side*0.3}
            },
            exprtree, lhs, rhs, curvature) );



   /* release */
   SCIP_CALL( SCIPexprtreeFree(&exprtree) );
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
   SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &sepadata->exprinterpreter) );
   sepadata->nlpinvars = SCIPgetNVars(scip);
   cr_assert_eq(sepadata->nlpinvars, 2, "error: received %d vars", sepadata->nlpinvars);

   /* create nlpi problem */
   sepadata->nlpi = SCIPgetNlpis(scip)[0];
   cr_assert_not_null(sepadata->nlpi);

   SCIP_CALL( SCIPnlpiCreateProblem(sepadata->nlpi, &sepadata->nlpiprob, "convexproj-nlp-unittest") );
   SCIP_CALL( SCIPhashmapCreate(&sepadata->var2nlpiidx, SCIPblkmem(scip), sepadata->nlpinvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &sepadata->nlpivars, SCIPgetVars(scip), sepadata->nlpinvars) );

   /* I shouldn't care about the cutoff, just assert that the lp solution satisfies the cutoff bound */
   SCIP_CALL( SCIPcreateNlpiProb(scip, sepadata->nlpi, nlrows, 3,
            sepadata->nlpiprob, sepadata->var2nlpiidx, NULL, SCIPgetCutoffbound(scip), FALSE, TRUE) );

   /* set quadratic part of objective function */
   SCIP_CALL( setQuadraticObj(scip, sepadata) );
}

/* setup of test */
static
void test_setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_VAR* auxx;

   /* include NLPI's */
   SCIP_CALL( SCIPcreateNlpSolverIpopt(SCIPblkmem(scip), &nlpi) );
   /* if no IPOPT available, don't run test */
   if( nlpi == NULL )
      return;

   SCIP_CALL( SCIPincludeNlpi(scip, nlpi) );
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameIpopt(), SCIPgetSolverDescIpopt()) );

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
   //SCIP_CALL( SCIPreleaseVar(scip, &x) );
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
