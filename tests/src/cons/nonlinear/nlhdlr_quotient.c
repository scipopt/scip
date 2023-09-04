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

/**@file   nlhdlr_quotient.c
 * @brief  tests quotient nonlinear handler methods
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include "scip/nlhdlr.c"
#include "scip/cons_nonlinear.c"
#include "scip/nlhdlr_quotient.c"


/*
 * TEST
 */

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* w;
static SCIP_NLHDLR* nlhdlr;
static SCIP_CONSHDLR* conshdlr;

/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, "quotient");
   cr_assert_not_null(nlhdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* go to SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 1.5, 5.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -4.0, 0.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 1.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -4.0, -1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPaddVar(scip, w) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   BMScheckEmptyMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* test suite */
TestSuite(nlhdlrquotient, .init = setup, .fini = teardown);

/** checks whether the values in nlhdlrexprdata are as expected */
static
void checkData(
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< the nlhdlr expression data */
   SCIP_VAR*             numvar,             /**< expected auxiliary variable of numerator expression */
   SCIP_Real             numcoef,            /**< expected numerator coefficient */
   SCIP_Real             numconst,           /**< expected numerator constant */
   SCIP_VAR*             denomvar,           /**< expected auxiliary variable of denominator expression */
   SCIP_Real             denomcoef,          /**< expected denominator coefficient */
   SCIP_Real             denomconst,         /**< expected denominator constant */
   SCIP_Real             constant            /**< expected constant */
   )
{
   cr_expect_not_null(nlhdlrexprdata->numexpr);
   cr_expect_not_null(nlhdlrexprdata->denomexpr);
   cr_expect(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->numexpr) == numvar);
   cr_expect(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->denomexpr) == denomvar);
   cr_expect(SCIPisEQ(scip, numcoef, nlhdlrexprdata->numcoef));
   cr_expect(SCIPisEQ(scip, numconst, nlhdlrexprdata->numconst));
   cr_expect(SCIPisEQ(scip, denomcoef, nlhdlrexprdata->denomcoef));
   cr_expect(SCIPisEQ(scip, denomconst, nlhdlrexprdata->denomconst));
   cr_expect(SCIPisEQ(scip, constant, nlhdlrexprdata->constant));
}

/* detects x / y */
Test(nlhdlrquotient, detectandfree1, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_Bool success;
   int i;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[nonlinear] <test>: <x> / <y> <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1) );

   expr = SCIPgetExprNonlinear(cons);
   ownerdata = SCIPexprGetOwnerData(expr);

   /* find the nlhdlr expr data */
   for( i = 0; i < SCIPgetExprNEnfosNonlinear(expr); ++i )
   {
      if( ownerdata->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = ownerdata->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 1.0, 0.0, y, 1.0, 0.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects (4x + 1) / (-3x - 3) */
Test(nlhdlrquotient, detectandfree2, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_Bool success;
   int i;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[nonlinear] <test>: (4*<x> + 1) / (-3*<x> - 3) <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1) );

   expr = SCIPgetExprNonlinear(cons);
   ownerdata = SCIPexprGetOwnerData(expr);

   /* find the nlhdlr expr data */
   for( i = 0; i < SCIPgetExprNEnfosNonlinear(expr); ++i )
   {
      if( ownerdata->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = ownerdata->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 4.0, 1.0, x, -3.0, -3.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects log((4x + 3) / (x + 1)) */
Test(nlhdlrquotient, detectandfree3, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_Bool success;
   int i;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[nonlinear] <test>: log((4*<x> + 3) / (<x> + 1)) <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1) );

   expr = SCIPexprGetChildren(SCIPgetExprNonlinear(cons))[0];
   ownerdata = SCIPexprGetOwnerData(expr);

   /* find the nlhdlr expr data */
   for( i = 0; i < SCIPgetExprNEnfosNonlinear(expr); ++i )
   {
      if( ownerdata->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = ownerdata->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 4.0, 3.0, x, 1.0, 1.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects that (4x + 2y + 3) / (x + 1)) is invalid */
Test(nlhdlrquotient, detectandfree4, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_Bool success;
   int i;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[nonlinear] <test>: (4*<x> + 2*<y> + 3) / (<x> + 1) <= 10",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1) );

   expr = SCIPgetExprNonlinear(cons);
   ownerdata = SCIPexprGetOwnerData(expr);

   /* find the nlhdlr expr data */
   for( i = 0; i < SCIPgetExprNEnfosNonlinear(expr); ++i )
   {
      if( ownerdata->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = ownerdata->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_null(nlhdlrexprdata);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects log(x) / log(y)) */
Test(nlhdlrquotient, detectandfree5, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_VAR* auxvarabs;
   SCIP_VAR* auxvarlog;
   SCIP_Bool success;
   int i;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[nonlinear] <test>: log(<x>) / abs(<y>) <= 10",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1) );

   expr = SCIPgetExprNonlinear(cons);
   ownerdata = SCIPexprGetOwnerData(expr);

   /* find the nlhdlr expr data */
   for( i = 0; i < SCIPgetExprNEnfosNonlinear(expr); ++i )
   {
      if( ownerdata->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = ownerdata->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   auxvarabs = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(SCIPexprGetChildren(expr)[0])[0]);
   auxvarlog = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[1]);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, auxvarlog, 1.0, 0.0, auxvarabs, 1.0, 0.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects (4x + 1) / (-3x + 3) + 2 after simplification */
Test(nlhdlrquotient, detectandfree6, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[nonlinear] <test>: ((4*<x> + 1) / (-3*<x> - 3) + 2) <= 3",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_expect_not(infeasible);

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1) );

   expr = SCIPgetExprNonlinear(cons);
   ownerdata = SCIPexprGetOwnerData(expr);

   /* find the nlhdlr expr data */
   for( i = 0; i < SCIPgetExprNEnfosNonlinear(expr); ++i )
   {
      if( ownerdata->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = ownerdata->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 4.0, 1.0, x, -3.0, -3.0, 2.0);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableNonlinear(scip, conshdlr, cons) );
   ownerdata->nconss = 0; /* TODO should consDisableNonlinear take care of this instead? */

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* tests interval evaluation for ((+/-)4x + 1) / (-3x + 3) - 2 */
Test(nlhdlrquotient, inteval, .description = "tests interval evaluation of simple quotient expression")
{
   SCIP_INTERVAL varbnds;
   SCIP_INTERVAL result;

   /* test interval including 0 in denominator*/

   varbnds.inf = 0.0;
   varbnds.sup = 2.0;

   result = intEvalQuotient(scip, varbnds, 4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, result));

   /* test positive denominator part for monotone increasing expression */

   varbnds.inf = 2.0;
   varbnds.sup = 9.0;

   result = intEvalQuotient(scip, varbnds, 4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, -5.0));
   cr_expect(SCIPisEQ(scip, result.sup, -37.0 / 24.0 - 2.0), "expected %f, but got %f\n",
      -37.0 / 24.0 - 2.0, result.sup);

   /* test negative denominator part for monotone increasing expression */

   varbnds.inf = -1.0;
   varbnds.sup = 0.0;

   result = intEvalQuotient(scip, varbnds, 4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, -2.5));
   cr_expect(SCIPisEQ(scip, result.sup, 1.0 / 3.0 - 2.0));


   /* test positive denominator part for monotone decreasing expression */

   varbnds.inf = 2.0;
   varbnds.sup = 9.0;

   result = intEvalQuotient(scip, varbnds, -4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, 35.0 / 24.0 - 2.0));
   cr_expect(SCIPisEQ(scip, result.sup, 7.0 / 3.0 - 2.0));

   /* test negative denominator part for monotone decreasing expression */

   varbnds.inf = -1.0;
   varbnds.sup = 0.0;

   result = intEvalQuotient(scip, varbnds, -4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, 1.0 / 3.0 - 2.0));
   cr_expect(SCIPisEQ(scip, result.sup, 5.0 / 6.0 - 2.0));
}

/* tests reverse propagation for univariate quotients */
Test(nlhdlrquotient, reverseprop, .description = "tests reverse propagation simple univariate quotient expressions")
{
   SCIP_INTERVAL bnds;
   SCIP_INTERVAL result;

   /* x / (x + 1) in [-3,-1] => x in [-0.75,0.5]*/
   SCIPintervalSetBounds(&bnds, -3.0, -1.0);
   result = reversepropQuotient(bnds, 1.0, 0.0, 1.0, 1.0, 0.0);
   cr_expect(SCIPisEQ(scip, result.inf, -0.75));
   cr_expect(SCIPisEQ(scip, result.sup, -0.5));

   /* x / (x + 1) in [-2,0.9] => x in [-2/3,9]*/
   SCIPintervalSetBounds(&bnds, -2.0, 0.9);
   result = reversepropQuotient(bnds, 1.0, 0.0, 1.0, 1.0, 0.0);
   cr_expect(SCIPisEQ(scip, result.inf, -2.0 / 3.0));
   cr_expect(SCIPisEQ(scip, result.sup, 9.0));

   /* (-5x + 2) / (3*x + 3) + 6 in [3,5] => x in [-inf,+inf]*/
   SCIPintervalSetBounds(&bnds, 3.0, 5.0);
   result = reversepropQuotient(bnds, -5.0, 2.0, 3.0, 3.0, 6.0);
   cr_expect(SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, result));

   /* (-5x + 2) / (3*x + 3) + 6 in [-2,-1] => x in [-23/16,-26/19]*/
   SCIPintervalSetBounds(&bnds, -2.0, -1.0);
   result = reversepropQuotient(bnds, -5.0, 2.0, 3.0, 3.0, 6.0);
   cr_expect(SCIPisEQ(scip, result.inf, -23.0/16.0));
   cr_expect(SCIPisEQ(scip, result.sup, -26.0/19.0));
}

/* estimates x = 2 for (4x + 1) / (-3x + 3) - 2 and x in [1.5,5] */
Test(nlhdlrquotient, estimation1, .description = "estimates simple univariate quotient expression")
{
   SCIP_Real constant;
   SCIP_Real coef;
   SCIP_Bool success;
   SCIP_Bool local;
   SCIP_Bool branchinguseful;
   SCIP_Real lbx = 1.5;
   SCIP_Real ubx = 5.0;
   SCIP_Real gllbx = 1.5;
   SCIP_Real glubx = 6.0;

   /*
    * tests overestimation
    */
   SCIP_CALL( estimateUnivariate(scip, lbx, ubx, gllbx, glubx, 2.0, 4.0, 1.0, -3.0, 3.0, -2.0, &coef, &constant, TRUE,
      &local, &branchinguseful, &success) );
   cr_expect(success);
   cr_expect(!branchinguseful);
   cr_expect(!local);
   cr_expect(SCIPisEQ(scip, coef, 5.0 / 3.0), "got %g expected %g", coef, 5.0 / 3.0);
   cr_expect(SCIPisEQ(scip, constant, -25.0 / 3.0), "got %g expected %g", constant, -25.0 / 3.0);

   /*
    * tests underestimation
    */
   SCIP_CALL( estimateUnivariate(scip, lbx, ubx, gllbx, glubx, 2.0, 4.0, 1.0, -3.0, 3.0, -2.0, &coef, &constant, FALSE,
      &local, &branchinguseful, &success) );
   cr_expect(success);
   cr_expect(branchinguseful);
   cr_expect(local);
   cr_expect(SCIPisEQ(scip, coef, 5.0 / 6.0), "got %g expected %g", coef, 5.0 / 6.0);
   cr_expect(SCIPisEQ(scip, constant, -95.0 / 12.0), "got %g expected %g", constant, -95.0 / 12.0);
}

/* estimates x = -1 for (4x + 1) / (-3x + 3) - 2 and x in [-4,0] */
Test(nlhdlrquotient, estimation2, .description = "estimates simple univariate quotient expression")
{
   SCIP_Real constant;
   SCIP_Real coef;
   SCIP_Bool success;
   SCIP_Bool local;
   SCIP_Bool branchinguseful;
   SCIP_Real lbx = -4.0;
   SCIP_Real ubx = 0.0;
   SCIP_Real gllbx = -5.0;
   SCIP_Real glubx = 0.0;

   /*
    * tests overestimation
    */
   SCIP_CALL( estimateUnivariate(scip, lbx, ubx, gllbx, glubx, -1.0, 4.0, 1.0, -3.0, 3.0, -2.0, &coef, &constant, TRUE,
      &local, &branchinguseful, &success) );
   cr_expect(success);
   cr_expect(branchinguseful);
   cr_expect(local);
   cr_expect(SCIPisEQ(scip, coef, 1.0 / 3.0), "got %g expected %g", coef, 1.0 / 3.0);
   cr_expect(SCIPisEQ(scip, constant, -5.0 / 3.0), "got %g expected %g", constant, -5.0 / 3.0);

   /*
    * tests underestimation
    */
   SCIP_CALL( estimateUnivariate(scip, lbx, ubx, gllbx, glubx, -1.0, 4.0, 1.0, -3.0, 3.0, -2.0, &coef, &constant, FALSE,
      &branchinguseful, &local, &success) );
   cr_expect(success);
   cr_expect(!branchinguseful);
   cr_expect(!local);
   cr_expect(SCIPisEQ(scip, coef, 5.0 / 12.0), "got %g expected %g", coef, 5.0 / 12.0);
   cr_expect(SCIPisEQ(scip, constant, -25.0 / 12.0), "got %g expected %g", constant, -25.0 / 12.0);
}

/* estimates (x,y) = (3,2) for x/y for x in [1,4] and y in [1.5,5] */
Test(nlhdlrquotient, estimation3, .description = "estimates simple bivariate quotient expression")
{
   SCIP_Real vals[3];
   SCIP_Bool success;
   SCIP_Bool branchingusefulx;
   SCIP_Bool branchingusefuly;

   /*
    * test overestimation
    */

   SCIP_CALL( estimateBivariate(scip, 1.0, 4.0, 1.5, 5.0, -SCIPinfinity(scip), SCIPinfinity(scip), 3.0, 2.0, 0.0, TRUE,
      &vals[0], &vals[1], &vals[2], &branchingusefulx, &branchingusefuly, &success) );
   cr_expect(branchingusefulx);
   cr_expect(branchingusefuly);
   cr_expect(success);
   cr_expect(SCIPisEQ(scip, vals[0], 2.0 / 3.0), "got %g expected %g", vals[0], 2.0 / 3.0);
   cr_expect(SCIPisEQ(scip, vals[1], -1.0 / 7.5), "got %g expected %g", vals[1], -1.0 / 7.5);
   cr_expect(SCIPisEQ(scip, vals[2], 1.0 / 5.0), "got %g expected %g", vals[2], 1.0 / 5.0);

   /*
    * test underestimation
    */

   SCIP_CALL( estimateBivariate(scip, 1.0, 4.0, 1.5, 5.0, -SCIPinfinity(scip), SCIPinfinity(scip), 3.0, 2.0, 0.0, FALSE,
      &vals[0], &vals[1], &vals[2], &branchingusefulx, &branchingusefuly, &success) );
   cr_expect(branchingusefulx);
   cr_expect(!branchingusefuly);
   cr_expect(success);
   cr_expect(SCIPisEQ(scip, vals[0], 5.0 / 9.0), "got %g expected %g", vals[0], 5.0 / 9.0);
   cr_expect(SCIPisEQ(scip, vals[1], -25.0 / 36.0), "got %g expected %g", vals[1], -25.0 / 36.0);
   cr_expect(SCIPisEQ(scip, vals[2], 10.0 / 9.0), "got %g expected %g", vals[2], 10.0 / 9.0);
}

/* estimates (x,y) = (-3,2) for x/y for x in [-4,-1] and y in [1.5,5] */
Test(nlhdlrquotient, estimation4, .description = "estimates simple bivariate quotient expression")
{
   SCIP_Real vals[3];
   SCIP_Bool success;
   SCIP_Bool branchingusefulx;
   SCIP_Bool branchingusefuly;

   /*
    * test overestimation
    */

   SCIP_CALL( estimateBivariate(scip, -4.0, -1.0, 1.5, 5.0, -SCIPinfinity(scip), SCIPinfinity(scip), -3.0, 2.0, 0.0, TRUE,
      &vals[0], &vals[1], &vals[2], &branchingusefulx, &branchingusefuly, &success) );
   cr_expect(success);
   cr_expect(branchingusefulx);
   cr_expect(!branchingusefuly);
   cr_expect(SCIPisEQ(scip, vals[0], 5.0 / 9.0), "got %g expected %g", vals[0], 5.0 / 9.0);
   cr_expect(SCIPisEQ(scip, vals[1], 25.0 / 36.0), "got %g expected %g", vals[1], 25.0 / 36.0);
   cr_expect(SCIPisEQ(scip, vals[2], -10.0 / 9.0), "got %g expected %g", vals[2], -10.0 / 9.0);

   /*
    * test underestimation
    */

   SCIP_CALL( estimateBivariate(scip, -4.0, -1.0, 1.5, 5.0, -SCIPinfinity(scip), SCIPinfinity(scip), -3.0, 2.0, 0.0, FALSE,
      &vals[0], &vals[1], &vals[2], &branchingusefulx, &branchingusefuly, &success) );
   cr_expect(success);
   cr_expect(branchingusefulx);
   cr_expect(branchingusefuly);
   cr_expect(SCIPisEQ(scip, vals[0], 2.0 / 3.0), "got %g expected %g", vals[0], 2.0 / 3.0);
   cr_expect(SCIPisEQ(scip, vals[1], 2.0 / 15.0), "got %g expected %g", vals[1], 2.0 / 15.0);
   cr_expect(SCIPisEQ(scip, vals[2], -1.0 / 5.0), "got %g expected %g", vals[2], -1.0 / 5.0);
}

/* TODO add a test for the case that 0 is in the domain of the numerator */
