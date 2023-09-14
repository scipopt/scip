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

/**@file   nlhdlr_concave.c
 * @brief  tests concave nonlinear handler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/scipdefplugins.h"
#include "scip/nlpi_ipopt.h" /* to check whether LAPACK is around */

#define NLHDLR_CONVEX_UNITTEST
#include "scip/nlhdlr_convex.c"


/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x_1;
static SCIP_VAR* x_2;
static SCIP_VAR* x_3;

static SCIP_CONSHDLR* conshdlr;
static SCIP_NLHDLR* nlhdlr = NULL;

/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert_not_null(conshdlr);

   /* get nlhdlr */
   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, CONCAVE_NLHDLR_NAME);
   cr_assert_not_null(nlhdlr);

   /* enable quadratic convexity check */
   SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/concave/cvxquadratic", TRUE) );

   SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE);
   SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x_1, "x1", 0.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x_2, "x2", 0.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x_3, "x3", 1.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x_1) );
   SCIP_CALL( SCIPaddVar(scip, x_2) );
   SCIP_CALL( SCIPaddVar(scip, x_3) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x_2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x_3) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/** given a string for f(x) and its curvature, run nlhdlr_concave detect on f(x) = 0 and see whether that gives correct flags */
static
SCIP_RETCODE detect(
   const char*           exprstr,
   SCIP_EXPRCURV         exprrootcurv,
   SCIP_Bool             simplify
)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* oexpr;
   SCIP_EXPR* expr;
   SCIP_Bool changed;
   SCIP_Bool infeas;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_CONS* cons;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseExpr(scip, &oexpr, (char*)exprstr, NULL, NULL, NULL) );
   if( simplify )
   {
      SCIP_CALL( SCIPsimplifyExpr(scip, oexpr, &expr, &changed, &infeas, NULL, NULL) );
      cr_expect(!infeas);
      SCIP_CALL( SCIPreleaseExpr(scip, &oexpr) );
   }
   else
      expr = oexpr;
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, (char*)"nlin", expr, 0.0, 0.0)  );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   SCIPprintCons(scip, cons, NULL);
   SCIPinfoMessage(scip, NULL, " and %s\n", SCIPexprcurvGetName(exprrootcurv));

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectConcave(scip, conshdlr, nlhdlr, expr, cons, &enforcing, &participating, &nlhdlrexprdata) );

   cr_expect_eq(enforcing, participating);
   if( (exprrootcurv & SCIP_EXPRCURV_CONVEX) != 0 )
   {
      cr_expect((enforcing & SCIP_NLHDLR_METHOD_SEPAABOVE) != 0);
   }

   if( (exprrootcurv & SCIP_EXPRCURV_CONCAVE) != 0 )
   {
      cr_expect((enforcing & SCIP_NLHDLR_METHOD_SEPABELOW) != 0);
   }

   if( participating != SCIP_NLHDLR_METHOD_NONE )
   {
      cr_assert_not_null(nlhdlrexprdata);
      SCIP_CALL( nlhdlrfreeExprDataConvexConcave(scip, nlhdlr, expr, &nlhdlrexprdata) );
   }

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/* tests detection of convex/concave subexpressions */
Test(nlhdlrconcave, detect, .init = setup, .fini = teardown)
{
   detect("exp(exp(<x1>))", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("exp(exp(log(<x1>)))", SCIP_EXPRCURV_CONVEX, FALSE);

   /* use signomial rule */
   detect("<x1>^0.5*<x2>^0.5)", SCIP_EXPRCURV_CONCAVE, FALSE);
   detect("<x1>^0.5*(exp(<x2>))^0.5", SCIP_EXPRCURV_CONCAVE, FALSE);  /* here an auxvar will be introduced for exp() */
   detect("<x1>^0.5*(log(<x2>))^0.5", SCIP_EXPRCURV_CONCAVE, FALSE);  /* here no auxvar will be introduced for exp() */
   detect("<x1>^(-0.5)*<x2>^(-0.5))", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("<x1>^3*<x2>^(-1)*<x3>^(-1)", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("<x1>^3*<x2>^(-1)*log(<x3>)^(-1)", SCIP_EXPRCURV_CONVEX, FALSE);

   /* product-composition */
   detect("<x1>*exp(<x1>)", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("<x1>*exp(<x1>+1.0)", SCIP_EXPRCURV_CONVEX, TRUE);
   detect("<x1>*exp(2*<x1>)", SCIP_EXPRCURV_CONVEX, TRUE);
   detect("(<x1>+<x2>)*exp(<x1>+<x2>)", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("exp(<x1>)*<x1>", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("exp(<x1>^2)*<x1>^2", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("exp(2*<x1>^2)*<x1>^2", SCIP_EXPRCURV_CONVEX, TRUE);
   detect("log(4-<x1>)*<x1>", SCIP_EXPRCURV_CONCAVE, TRUE);   /* similar to arki0017 */

   /* quadratic */
   detect("log(-<x1>^2-<x2>^2-2*<x1>*<x2>-<x3>^2)", SCIP_EXPRCURV_CONCAVE, FALSE);
   if( SCIPisIpoptAvailableIpopt() )
      detect("-2*<x1>^2+<x1>*<x2>-2*<x2>^2", SCIP_EXPRCURV_CONCAVE, TRUE);
}
