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
#define SCIP_DEBUG
#include "scip/scipdefplugins.h"
#include "scip/scip.h"
#include "scip/nlhdlr.c"
#include "scip/cons_nonlinear.c"
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
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, "signomial");
   cr_assert_not_null(nlhdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );

   /* go to SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

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


   BMSdisplayMemory();
   BMScheckEmptyMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
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

/** a probabilistic way to valid a rowprep by random sampling */
static 
SCIP_RETCODE validCutProb(
   SCIP_EXPR*            expr,               /**< expr pointer */ 
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
            //SCIPdebugMsg("%d\n",s);
            *isvalid = FALSE;
            break;
         }
      }
      else
      {
         if( cutval > side ){
            //SCIPdebugMsg("%d\n",s);
            *isvalid = FALSE;
            break;
         }
      }
   }

   SCIPfreeRandom(scip, &randnumgen);
   return SCIP_OKAY;
}

/*
 * tests
 */


/* detects <x1> * (<x2>)^(1.0) * (<x4>)^(1.0) */
Test(nlhdlrsignomial, detectandfree1, .description = "detects signomial terms 1")
{
   SCIP_CONS* cons;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_Bool success;
   int i;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[nonlinear] <test>: <x1> * (<x2>)^(1.0) * (<x4>)^(1.0)  <= 100",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/*  detects  (<x3>)^(-1.5) */
Test(nlhdlrsignomial, detectandfree2, .description = "detects signomial terms 2")
{
   SCIP_CONS* cons;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_Bool success;
   int i;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[nonlinear] <test>:  (<x3>)^(-1.5)  <= 100",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}


/*  detects (<x3>)^(1.5) * (<x4>)^(-1.5) * (<x5>)^(-2)) */
Test(nlhdlrsignomial, detectandfree3, .description = "detects signomial terms 3")
{
   SCIP_CONS* cons;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_Bool success;
   int i;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create nonlinear constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[nonlinear] <test>:  (<x3>)^(1.5) * (<x4>)^(-1.5) * (<x5>)^(-2)  <= 100",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   cr_assert(nlhdlrexprdata->nvars == 4);
   cr_assert(nlhdlrexprdata->nposvars == 1);
   cr_assert(nlhdlrexprdata->nnegvars == 3);

   cr_assert(SCIPisEQ(scip, nlhdlrexprdata->exponents[0], 1.5));
   cr_assert(SCIPisEQ(scip, nlhdlrexprdata->refexponents[0], 1.5 / (1.5 + 2 + 1) ));
   cr_assert(SCIPisEQ(scip, nlhdlrexprdata->exponents[1], -1.5));
   cr_assert(SCIPisEQ(scip, nlhdlrexprdata->refexponents[1], 1.5 / (1.5 + 2 + 1) ));
   cr_assert(SCIPisEQ(scip, nlhdlrexprdata->exponents[2], -2));
   cr_assert(SCIPisEQ(scip, nlhdlrexprdata->refexponents[2], 2 / (1.5 + 2 + 1) ));
   cr_assert(SCIPisEQ(scip, nlhdlrexprdata->refexponents[3], 1 / (1.5 + 2 + 1) ));

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIPdebugMsg(scip, "1\n");
}






