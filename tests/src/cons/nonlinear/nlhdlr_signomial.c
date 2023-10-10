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

/**@file   nlhdlr_signomial.c
 * @brief  tests signomial nonlinear handler methods
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/scipdefplugins.h"
#include "scip/scip.h"
#include "scip/nlhdlr_signomial.c"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* x3;
static SCIP_VAR* x4;
static SCIP_VAR* x5;
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

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", 1e-2, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", 1e-2, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3, "x3", 1e-2, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4, "x4", 1e-2, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x5, "x5", 1e-2, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, x3) );
   SCIP_CALL( SCIPaddVar(scip, x4) );
   SCIP_CALL( SCIPaddVar(scip, x5) );

   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, "signomial");
   cr_assert_not_null(nlhdlr);
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4) );
   SCIP_CALL( SCIPreleaseVar(scip, &x5) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/* define the test suite */
TestSuite(nlhdlrsignomial, .init = setup, .fini = teardown);

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

/* creates and adds two nonlinear constraints and check output of SCIPgetSignomialExprsExpr */
Test(nlhdlrsignomial, collect_product_expressions)
{
   SCIP_CONS* conss[3];
   SCIP_EXPR** exprs;
   SCIP_EXPR* expr;
   const char* inputs[3] = {
    "<x1> * power(<x2>, 1.0) * power(<x4>, 1.0) + power(<x3>, -2.1) * power(<x4>, 1.1)",
    "power(<x1>, 0.8) * power(<x3>, 1.5) * power(<x4>, -1.5) * power(<x5>, -2) + power(<x2>, -2.1) * power(<x3>, 1.2) +  power(<x4>, 1.8)",
    "power(<x1>, 0.6) * power(<x2>, -1.5) * power(<x3>, -1) * power(<x4>, 2) * power(<x5>, 0.8) + power(<x1>, -1.1) * power(<x2>, 1.2)  * power(<x3>, 0.4)"};
   SCIP_Bool infeasible;
   int nexprs;
   int c;
   int t;

   /* create constraints */
   for( c = 0; c < 3; ++c )
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
   exprs = SCIPgetExprsSignomial(nlhdlr);
   cr_expect(exprs != NULL);
   nexprs = SCIPgetNExprsSignomial(nlhdlr);
   cr_expect(nexprs == 5);

   for( c = 0; c < nexprs; ++c )
   {
      SCIP_EXPR* child;

      cr_assert(exprs[c] != NULL);
      cr_expect(SCIPexprGetNChildren(exprs[c]) >= 2);
      cr_expect(SCIPgetExprAuxVarNonlinear(exprs[c]) != NULL);
      cr_expect(SCIPisExprProduct(scip, exprs[c]));

      for ( t = 0; t < SCIPgetExprAuxVarNonlinear(exprs[c]); ++t)
      {
         child = SCIPexprGetChildren(exprs[c])[t];
         cr_assert(child != NULL);
         cr_expect(SCIPgetExprAuxVarNonlinear(child) != NULL);
      }
   }
}



/** a probabilistic way to valid a rowprep by random sampling */
static 
SCIP_RETCODE validCutProb(
   SCIP*                 scip,               /**< scip pointer */ 
   SCIP_EXPR*            expr.               /**< expr pointer */ 
   SCIP_ROWPREP*         rowprep,            /**< rowprep pointer */ 
   int                   nsample,            /**< number of samples */ 
   SCIP_Bool*            isvalid             /**< whether the rowprep is valid */ 
   )
{
   int s;
   int i;
   int j;
   SCIP_Real cutval;
   SCIP_Real xval;
   SCIP_Real yval;
   SCIP_Real ycoef;
   SCIP_Real coef;
   SCIP_Real exponent;
   int nvars = SCIProwprepGetNVars(rowprep);
   SCIP_Real* coefs =  SCIProwprepGetCoefs(rowprep);
   SCIP_Real side = SCIProwprepGetSide(rowprep);
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata = SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr);    
   SCIP_VAR** vars = SCIProwprepGetVars(rowprep);
   SCIP_VAR* xvar;
   unsigned int seedp  = 2132;
   SCIP_RANDNUMGEN* randnumgen;


   SCIPcreateRandom(scip, &randnumgen, seedp, TRUE);
   assert(randnumgen);

   *isvalid = TRUE;

   for( s = 0; s < nsample; s++ )
   {
      cutval = 0;
      yval = nlhdlrexprdata->coef;
      for( i = 0; i < nlhdlrexprdata->nfactors; i++ )
      {
         xvar = nlhdlrexprdata->vars[i];
         coef = 0;
         for( j = 0; j < nvars; j++ )
         {
            if( vars[j] == xvar ){
               coef = coefs[j];
               break;
            }
         }
         xval = SCIPrandomGetReal(randnumgen, nlhdlrexprdata->intervals[i].inf, nlhdlrexprdata->intervals[i].sup);
         cutval += xval * coef;
         exponent = nlhdlrexprdata->exponents[i];
         yval *= pow(xval, exponent);
      }
      ycoef = 0;
      for( j = 0; j < nvars; j++ )
      {
         if( vars[j] == nlhdlrexprdata->vars[nlhdlrexprdata->nfactors] ){
            ycoef = coefs[j];
            break;
         }
      }     
      cutval += yval * ycoef;
      if( SCIProwprepGetSidetype(rowprep) == SCIP_SIDETYPE_LEFT)
      {
         if( cutval < side ){
            //SCIPdebugMsg(scip, "%d/%d, (lhs) %f <= %f \n", s, nsample, side, cutval);
            *isvalid = FALSE;
            break;
         }
      }
      else
      {
         if( cutval > side ){
            //SCIPdebugMsg(scip, "%d/%d, %f <= %f (rhs) \n", s, nsample, cutval, side);
            *isvalid = FALSE;
            break;
         }
      }
   }

   SCIPfreeRandom(scip, &randnumgen);
   return SCIP_OKAY;
}


/* adds a linear inequality to the product expression and computes a cut at a given reference point */
Test(nlhdlrsignomial, separation_single)
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
   SCIP_Real targetval;
   int nsample;
   SCIP_Bool isvalid;

   /* set variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, x1, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x1, 11.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x2, 1.1) );
   SCIP_CALL( SCIPchgVarUb(scip, x2, 12.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x3, 1.2) );
   SCIP_CALL( SCIPchgVarUb(scip, x3, 10.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x4, 1.3) );
   SCIP_CALL( SCIPchgVarUb(scip, x4, 9.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x5, 1.5) );
   SCIP_CALL( SCIPchgVarUb(scip, x5, 12.0) );

   /* transform problem */
   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE);

   /* create constraint containing a single product expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "power(<t_x1>, 0.6) * power(<t_x2>, -1.5) * power(<t_x3>, -1) * power(<t_x4>, 2) * power(<t_x5>, 0.8)", NULL, NULL, NULL) );
   SCIP_CALL( createAndDetect(&cons, expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   SCIP_CALL( SCIPconstructLP(scip, &dummy) );

   /* INITLP should have added an auxiliary variable to the product expression (tight might change in the future) */
   cr_assert( SCIPgetExprAuxVarNonlinear(expr) != NULL);

   /* create a solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x1), 2.11) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x2), 2.22) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x3), 3.11) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x4), 4.22) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x4), 5.22) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr), -1.23) );

   SCIP_CALL( nlhdlrEvalauxSignomial(scip, nlhdlr, expr, SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr), &targetval, sol));
   overestimate = FALSE;
   SCIP_CALL( nlhdlrEstimateSignomial(scip, conshdlr, nlhdlr, expr, SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr), sol,
            -1.23, overestimate, targetval, FALSE, rowpreps, &success, &dummy) );
   cr_expect(success);

   nsample = 10000;
   rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, 0);
   SCIP_CALL( validCutProb(scip, expr, rowprep, nsample, &isvalid));
   SCIPfreeRowprep(scip, &rowprep);

   /* free memory */
   SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   cr_expect(isvalid);
}



