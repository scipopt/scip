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

/**@file   prodyct_detection.c
 * @brief  tests hidden bilinear product detection for rlt cuts
 * @author Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/sepa_rlt.c"
#include "scip/cons_nonlinear.h"
#include "scip/cons_linear.h"
#include "include/scip_test.h"

/* in all comments it is assumed that the linear relations used for product detections have the form:
 * a1x + b1w + c1y <= d1 (considered to be active when x == 1) and
 * a2x + b2w + c2y <= d2 (considered to be active when x == 0)
 */

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SEPA* sepa;
static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* x3;
static SCIP_VAR* x4;
static SCIP_VAR* bvar1;
static SCIP_VAR* bvar2;

/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_VAR* x1o;
   SCIP_VAR* x2o;
   SCIP_VAR* x3o;
   SCIP_VAR* x4o;
   SCIP_VAR* bvar1o;
   SCIP_VAR* bvar2o;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include nonlinear and linear constraint handlers and rlt separator */
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) );
   SCIP_CALL( SCIPincludeSepaRlt(scip) );

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   assert(conshdlr != NULL);

   /* get rlt separator */
   sepa = SCIPfindSepa(scip, "rlt");
   assert(sepa != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1o, "x1", -1.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2o, "x2", -6.0, -3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3o, "x3", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4o, "x4", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &bvar1o, "bvar1", 0, 1, 1, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &bvar2o, "bvar2", 0, 1, 1, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, x1o) );
   SCIP_CALL( SCIPaddVar(scip, x2o) );
   SCIP_CALL( SCIPaddVar(scip, x3o) );
   SCIP_CALL( SCIPaddVar(scip, x4o) );
   SCIP_CALL( SCIPaddVar(scip, bvar1o) );
   SCIP_CALL( SCIPaddVar(scip, bvar2o) );

   /* get SCIP into SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* get transformed vars and release vars */
   SCIP_CALL( SCIPgetTransformedVar(scip, x1o, &x1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x2o, &x2) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x3o, &x3) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x4o, &x4) );
   SCIP_CALL( SCIPgetTransformedVar(scip, bvar1o, &bvar1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, bvar2o, &bvar2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4o) );
   SCIP_CALL( SCIPreleaseVar(scip, &bvar1o) );
   SCIP_CALL( SCIPreleaseVar(scip, &bvar2o) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* check an auxiliary expression by comparing simplified expressions */
static
void checkAuxExpr(
   SCIP_CONSNONLINEAR_BILINTERM term,
   int                   auxidx,             /* index of the term */
   SCIP_VAR*             auxvar,             /* expected auxiliary variable */
   SCIP_Real*            vals,               /* expected coefficients in the order w, x, y */
   SCIP_Real             cst,                /* expected constant */
   SCIP_Bool             underestimate,      /* should the auxexpr underestimate the product? */
   SCIP_Bool             overestimate        /* should the auxexpr overestimate the product? */
   )
{
   int i;

   cr_expect_eq(auxvar, term.aux.exprs[auxidx]->auxvar);

   for( i = 0; i < 3; ++i )
   {
      cr_expect(SCIPisEQ(scip, vals[i], term.aux.exprs[auxidx]->coefs[i]), "\ni = %d: expected != given: %g != %g", i,
                                                                           vals[i],
                                                                           term.aux.exprs[auxidx]->coefs[i]);
   }
   cr_expect(SCIPisEQ(scip, cst, term.aux.exprs[auxidx]->cst), "\nexpected != given: %g != %g", cst,
                                                               term.aux.exprs[auxidx]->cst);

   cr_expect(term.aux.exprs[auxidx]->underestimate == underestimate,
             "auxiliary expression [%d] should %sunderestimate", auxidx, underestimate ? "" : "NOT ");
   cr_expect(term.aux.exprs[auxidx]->overestimate == overestimate,
             "auxiliary expression [%d] should %soverestimate", auxidx, overestimate ? "" : "NOT ");
}

Test(product_detection, implrels, .init = setup, .fini = teardown, .description = "test extracting products from two implied relations")
{
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   SCIP_SEPADATA* sepadata;
   SCIP_VAR* vars[3];
   SCIP_Real coefs1[3];
   SCIP_Real coefs2[3];
   SCIP_Bool cutoff;
   SCIP_CONSNONLINEAR_BILINTERM* terms;

   vars[0] = x1;
   vars[1] = x2;
   vars[2] = bvar1;

   coefs1[0] = 2.0;
   coefs1[1] = 1.5;
   coefs1[2] = 3.0;

   coefs2[0] = 4.0;
   coefs2[1] = 2.0;
   coefs2[2] = -1.0;

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = conshdlr;
   sepadata->isinitialround = TRUE;
   sepadata->detecthidden = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear constraints */

   /* 0 <= 2x1 + 1.5x2 + 3bvar1 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons1, "c1", 3, vars, coefs1, 0.0, 1.0) );

   /* 4x1 + 2x2 - bvar1 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons2, "c2", 3, vars, coefs2, -SCIPinfinity(scip), 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons1) );
   SCIP_CALL( SCIPaddCons(scip, cons2) );

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* create separator data - this will also detect the products */
   SCIP_CALL( createSepaData(scip, sepadata) );

   /* the product should be: xy >=/<= (1/(b1c2 - c1b2))*(b1b2w + (b2(a1 - d1) + b1d2)x + b1c2y - b1d2)
    *
    * the first product (terms[0].auxexprs[0]) is derived from 2x1 + 1.5x2 + 3bvar1 <= 1 and
    * 4x1 + 2x2 - bvar1 <= 1 (x = bvar1, w = x1, y = x2):
    * xy <= (1/(2*2 - 1.5*4))*(2*4w + (4(3 - 1) + 2*1)x + 2*2y - 2*1)
    * <=>  xy <= -4w - 5x - 2y + 1, which in problem variables is: bvar1*x2 <= -4x1 - 5bvar1 - 2x2 + 1
    *
    * the second product (terms[1].auxexprs[3]) is derived from same relations (x = bvar1, w = x2, y = x1):
    * xy >= (1/(1.5*4 - 2*2))*(1.5*2w + (2(3 - 1) + 1.5*1)x + 1.5*4y - 1.5*1)
    * <=>  xy >= 1.5w + 2.75x + 3y - 0.75, which in problem variables is: bvar1*x1 >= 1.5x2 + 2.75bvar1 + 3x1 - 0.75
    */

   /* check the numbers */
   cr_expect_eq(sepadata->nbilinvars, 3, "\nExpected 3 bilinear vars, got %d", sepadata->nbilinvars);
   cr_expect_eq(SCIPgetNBilinTermsNonlinear(conshdlr), 2, "\nExpected 2 bilinear terms, got %d",
                                                          SCIPgetNBilinTermsNonlinear(conshdlr));

   terms = SCIPgetBilinTermsNonlinear(conshdlr);
   cr_assert(terms != NULL);

   cr_expect_eq(terms[0].nauxexprs, 4, "\nExpected 4 auxiliary expressions for product 0, got %d", terms[0].nauxexprs);
   cr_expect_eq(terms[1].nauxexprs, 4, "\nExpected 4 auxiliary expressions for product 1, got %d", terms[1].nauxexprs);

   /* check the product expressions */
   cr_expect_eq(terms[0].x, bvar1, "x var of product 0 should be bvar1, got %s", SCIPvarGetName(terms[0].x));
   cr_expect_eq(terms[0].y, x2, "y var of product 0 should be x2, got %s", SCIPvarGetName(terms[0].y));
   cr_expect_eq(terms[1].x, bvar1, "x var of product 1 should be bvar1, got %s", SCIPvarGetName(terms[1].x));
   cr_expect_eq(terms[1].y, x1, "y var of product 1 should be x1, got %s", SCIPvarGetName(terms[1].y));

   /* check the (sorted) auxiliary expressions and sides */
   cr_assert(terms[0].aux.exprs[0] != NULL);
   cr_assert(terms[1].aux.exprs[0] != NULL);

   /* first product: bvar1 * x2 */
   /* check 4 expressions with binary variable = bvar1:
    * 1) from two implied relations,
    * 2) from 1st implied relation (rhs) and x1 upper bound,
    * 3) from 1st implied relation (lhs) and x1 lower bound,
    * 4) from 2nd implied relation (rhs) and x1 upper bound
    */
   checkAuxExpr(terms[0], 0, x1, (SCIP_Real[3]) {-4.0, -5.0, -2.0}, 1.0, FALSE, TRUE);
   checkAuxExpr(terms[0], 1, x1, (SCIP_Real[3]) {-4.0/3.0, -8.0, 0.0}, 20.0/3.0, FALSE, TRUE);
   checkAuxExpr(terms[0], 2, x1, (SCIP_Real[3]) {4.0/3.0, 4.0/3.0, 1.0}, 0.0, FALSE, TRUE);
   checkAuxExpr(terms[0], 3, x1, (SCIP_Real[3]) {2.0, -9.5, 1.0}, -0.5, TRUE, FALSE);

   /* second product: bvar1 * x1, derived from two implied relations */
   checkAuxExpr(terms[1], 3, x2, (SCIP_Real[3]) {1.5, 2.75, 3.0}, -0.75, TRUE, FALSE);

   /* free memory */
   SCIP_CALL( freeSepaData(scip, sepadata) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(product_detection, implrelbnd, .init = setup, .fini = teardown, .description = "test extracting products from an implied relation and an implied bound")
{
   SCIP_CONS* cons;
   SCIP_SEPADATA* sepadata;
   SCIP_VAR* vars[3];
   SCIP_Real coefs[3];
   int nbdchgs;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_CONSNONLINEAR_BILINTERM* terms;

   vars[0] = x1;
   vars[1] = x2;
   vars[2] = bvar1;

   coefs[0] = 2.0;
   coefs[1] = 2.0;
   coefs[2] = 3.0;

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = conshdlr;
   sepadata->isinitialround = TRUE;
   sepadata->detecthidden = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear relations */

   /* 0 <= 2x1 + 2x2 + 3bvar1 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 3, vars, coefs, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* bvar1 == 0  =>  x1 <= 3 */
   /* x1 - Mbvar1 <= 3 */
   SCIP_CALL( SCIPaddVarImplication(scip, bvar1, FALSE, x1, SCIP_BOUNDTYPE_UPPER, 3.0, &infeasible, &nbdchgs) );

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* create separator data - this will also detect the products */
   SCIP_CALL( createSepaData(scip, sepadata) );

   /* the product should be: xy >=/<= (1/(b1c2 - c1b2))*(b1b2w + (b2(a1 - d1) + b1d2)x + b1c2y - b1d2)
    *
    * the product we check here (terms[0].auxexprs[1]) is derived from 2x1 + 2x2 + 3bvar1 <= 1 and the
    * implied bound (x = bvar1, w = x1, y = x2):
    * xy <= (1/(2*0 - 2*1))*(2*1w + (1*(3 - 1) + 2*3)x + 2*0y - 2*3)
    * <=> xy <= -w - 4x + 3, which in problem variables is: bvar1*x2 <= -x1 - 4bvar1 + 3
    */

   /* check the numbers */
   cr_expect_eq(sepadata->nbilinvars, 3, "\nExpected 3 bilinear vars, got %d", sepadata->nbilinvars);
   cr_expect_eq(SCIPgetNBilinTermsNonlinear(conshdlr), 2, "\nExpected 2 bilinear terms, got %d",
                                                          SCIPgetNBilinTermsNonlinear(conshdlr));

   terms = SCIPgetBilinTermsNonlinear(conshdlr);

   /* check the product expressions */
   cr_expect_eq(terms[0].x, bvar1, "x var of product 0 should be bvar1, got %s", SCIPvarGetName(terms[0].x));
   cr_expect_eq(terms[0].y, x2, "y var of product 0 should be x2, got %s", SCIPvarGetName(terms[0].y));
   cr_expect_eq(terms[1].x, bvar1, "x var of product 1 should be bvar1, got %s", SCIPvarGetName(terms[1].x));
   cr_expect_eq(terms[1].y, x1, "y var of product 1 should be x1, got %s", SCIPvarGetName(terms[1].y));

   cr_expect_eq(terms[0].nauxexprs, 3, "\nExpected 3 auxiliary expressions for product bvar1*x2, got %d", terms[0].nauxexprs);
   cr_expect_eq(terms[1].nauxexprs, 2, "\nExpected 2 auxiliary expressions for product bvar1*x1, got %d", terms[1].nauxexprs);

   /* product relation derived from an implied relation and an implied bound */
   checkAuxExpr(terms[0], 1, x1, (SCIP_Real[3]) {-1.0, -4.0, 0.0}, 3.0, FALSE, TRUE);

   /* free memory */
   SCIP_CALL( freeSepaData(scip, sepadata) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(product_detection, implrelclique, .init = setup, .fini = teardown, .description = "test extracting products from an implied relation and a clique")
{
   SCIP_CONS* cons1;
   SCIP_SEPADATA* sepadata;
   SCIP_VAR* vars[3];
   SCIP_Real coefs1[3];
   int nbdchgs;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_VAR* clique_vars[2];
   SCIP_Bool clique_vals[2];
   SCIP_CONSNONLINEAR_BILINTERM* terms;

   vars[0] = x1;
   vars[1] = bvar1;
   vars[2] = bvar2;

   coefs1[0] = 2.0;
   coefs1[1] = 2.0;
   coefs1[2] = 3.0;

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = conshdlr;
   sepadata->isinitialround = TRUE;
   sepadata->detecthidden = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear relations */

   /* 0 <= 2x1 + 2bvar1 + 3bvar2 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons1, "c1", 3, vars, coefs1, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons1) );

   /* add a clique bvar1 + (1-bvar2) <= 1 */
   clique_vars[0] = bvar1;
   clique_vars[1] = bvar2;
   clique_vals[0] = 1;
   clique_vals[1] = 0;
   scip->stat->nnz = 3; /* set this manually so that the clique is not rejected */
   SCIP_CALL( SCIPaddClique(scip, clique_vars, clique_vals, 2, FALSE, &infeasible, &nbdchgs) );

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* create separator data - this will also detect the products */
   SCIP_CALL( createSepaData(scip, sepadata) );

   /* the product should be: xy >=/<= (1/(b1c2 - c1b2))*(b1b2w + (b2(a1 - d1) + b1d2)x + b1c2y - b1d2)
    *
    * the first product relation (terms[0].auxexprs[1]) is derived from clique 0 <= -bvar1 + bvar2
    * and inequality 0 <= 2bvar1 + 3bvar2 + 2x1 (x = bvar1, w = bvar2, y = x1):
    * xy <= (1/(1*2 - 0*3))*(1*3w + (3(-1 - 0) + 1*0)x + 1*2y - 1*0)
    * <=> xy <= 1.5w - 1.5x + y, which in problem variables is bvar1*x1 <= 1.5bvar2 - 1.5bvar1 + x1
    *
    * the second product relation (terms[2].auxexprs[1]) is derived from inequality 3bvar2 + 2bvar1 + 2x1 <= 1
    * and clique -bvar2 + bvar1 <= 0 (x = bvar2, w = bvar1, y = x1)
    * xy <= (1/(2*0 - 2*1))*(2*1w + (1*(3 - 1) + 2*0)x + 2*0y - 2*0)
    * <=> xy <= -w - x, which in problem variables is bvar2*x1 <= -bvar1 - bvar2
    */

   terms = SCIPgetBilinTermsNonlinear(conshdlr);

   /* check the numbers */
   cr_expect_eq(sepadata->nbilinvars, 3, "\nExpected 3 bilinear vars, got %d", sepadata->nbilinvars);
   cr_expect_eq(SCIPgetNBilinTermsNonlinear(conshdlr), 3, "\nExpected 3 bilinear terms, got %d",
                                                          SCIPgetNBilinTermsNonlinear(conshdlr));
   cr_expect_eq(terms[0].nauxexprs, 3, "\nExpected 3 auxiliary expressions for product 0, got %d", terms[0].nauxexprs);
   cr_expect_eq(terms[1].nauxexprs, 4, "\nExpected 4 auxiliary expressions for product 1, got %d", terms[1].nauxexprs);
   cr_expect_eq(terms[2].nauxexprs, 3, "\nExpected 3 auxiliary expressions for product 2, got %d", terms[2].nauxexprs);

   /* check the product expressions */
   cr_expect_eq(terms[0].x, bvar1, "Var 0 of product 0 should be bvar1, got %s", SCIPvarGetName(terms[0].x));
   cr_expect_eq(terms[0].y, x1, "Var 1 of product 0 should be x1, got %s", SCIPvarGetName(terms[0].y));
   cr_expect_eq(terms[1].x, bvar1, "Var 0 of product 1 should be bvar1, got %s", SCIPvarGetName(terms[1].x));
   cr_expect_eq(terms[1].y, bvar2, "Var 1 of product 1 should be bvar2, got %s", SCIPvarGetName(terms[1].y));
   cr_expect_eq(terms[2].x, bvar2, "Var 0 of product 2 should be bvar2, got %s", SCIPvarGetName(terms[2].x));
   cr_expect_eq(terms[2].y, x1, "Var 1 of product 2 should be x1, got %s", SCIPvarGetName(terms[2].y));

   /* bvar1*x1 <= 1.5bvar2 - 1.5bvar1 + x1 (from constraint and clique, (x,w,y) = (bvar1,bvar2,x1)) */
   checkAuxExpr(terms[0], 1, bvar2, (SCIP_Real[3]) {1.5, -1.5, 1.0}, 0.0, FALSE, TRUE);

   /* bvar2*x1 <= -bvar1 - bvar2 (from constraint and clique, (x,w,y) = (bvar2,bvar1,x1)) */
   checkAuxExpr(terms[2], 1, bvar1, (SCIP_Real[3]) {-1.0, -1.0, 0.0}, 0.0, FALSE, TRUE);

   /* free memory */
   SCIP_CALL( freeSepaData(scip, sepadata) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(product_detection, implbnd, .init = setup, .fini = teardown, .description = "test extracting products from an implied bound and an unconditional relation")
{
   SCIP_CONS* cons;
   SCIP_SEPADATA* sepadata;
   SCIP_VAR* vars[2];
   SCIP_Real coefs[2];
   int nbdchgs;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_CONSNONLINEAR_BILINTERM* terms;

   vars[0] = x1;
   vars[1] = x2;

   coefs[0] = 2.0;
   coefs[1] = 2.0;

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = conshdlr;
   sepadata->isinitialround = TRUE;
   sepadata->detecthidden = TRUE;
   cr_assert(sepadata->conshdlr != NULL);

   /* create linear relations */

   /* 0 <= 2x1 + 2x2 <= 1 */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c1", 2, vars, coefs, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* bvar1 == 0  =>  x1 <= 1 */
   /* x1 <= -2bvar1 + 1 */
   SCIP_CALL( SCIPaddVarVub(scip, x1, bvar1, -2.0, 1.0, &infeasible, &nbdchgs) );

   /* x1 <= bvar2 + 0.0 */
   SCIP_CALL( SCIPaddVarVub(scip, x1, bvar2, 1.0, 0.0, &infeasible, &nbdchgs) );

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* create separator data - this will also detect the products */
   SCIP_CALL( createSepaData(scip, sepadata) );

   /* the product should be: xy >=/<= (1/(b1c2 - c1b2))*(b1b2w + (b2(a1 - d1) + b1d2)x + b1c2y - b1d2)
    *
    * the first product relation (terms[0].auxexprs[0]) is derived from two implied bounds on x1:
    * 2bvar1 + x1 <= 1 and x1 - bvar2 <= 0 (x = bvar1, w = x1, y = bvar2):
    * xy <= (1/(1*-1 - 0*1))*(1*1w + (1*(2 - 1) + 1*0)x + 1*-1y - 1*0)
    * <=> xy <= -w - x + y, which in problem variables is bvar1*bvar2 <= -x1 - bvar1 + bvar2
    *
    * the second product relation (terms[0].auxexprs[1]) is derived from two implied bounds on x1:
    * x1 + 2bvar1 <= 1 and -bvar2 + x1 <= 0 (x = bvar2, w = x1, y = bvar1):
    * xy <= (1/(1*0 - 2*1))*(1w + (1*(0 - 1) + 1*0)x + 1*0y - 1*0)
    * <=> xy <= -0.5w + 0.5x, which in problem variables is bvar2*bvar1 <= -0.5x1 + 0.5bvar2
    *
    * the third product relation (terms[1].auxexprs[0]) is derived from an implied bound
    * 2bvar1 + x1 <= 1 and inequality 2x1 + 2x2 <= 1 (x = bvar1, w = x1, y = x2):
    * xy >= (1/(b1c2 - c1b2))*(b1b2w + (b2(a1 - d1) + b1d2)x + b1c2y - b1d2)
    * xy >= (1/(1*2 - 0*2))*(1*2w + (2(2 - 1) + 1*1)x + 1*2y - 1*1)
    * <=> xy >= w + 1.5x + y - 0.5, which in problem variables is bvar1*x2 >= x1 + 1.5bvar1 + x2 - 0.5
    *
    * the fourth product relation (terms[2].auxexprs[0]) is derived from inqeuality
    * 2x1 + 2x2 <= 1 and implied bound -bvar2 + x1 <= 0 (x = bvar2, w = x1, y = x2):
    * xy <= (1/(b1c2 - c1b2))*(b1b2w + (b2(a1 - d1) + b1d2)x + b1c2y - b1d2)
    * xy <= (1/(2*0 - 2*1))*(2*1w + (1(0 - 1) + 2*0)x + 2*0y - 2*0)
    * <=> xy <= -w + 0.5x, which in problem variables is bvar2*x2 <= -x1 + 0.5bvar2
    */

   terms = SCIPgetBilinTermsNonlinear(conshdlr);

   /* check the numbers */
   cr_expect_eq(sepadata->nbilinvars, 3, "\nExpected 3 bilinear vars, got %d", sepadata->nbilinvars);
   cr_expect_eq(SCIPgetNBilinTermsNonlinear(conshdlr), 3, "\nExpected 3 bilinear terms, got %d",
                                                          SCIPgetNBilinTermsNonlinear(conshdlr));
   cr_expect_eq(terms[0].nauxexprs, 2, "\nExpected 2 auxiliary expressions for product 0, got %d", terms[0].nauxexprs);
   cr_expect_eq(terms[1].nauxexprs, 1, "\nExpected 1 auxiliary expressions for product 1, got %d", terms[1].nauxexprs);
   cr_expect_eq(terms[2].nauxexprs, 1, "\nExpected 1 auxiliary expressions for product 2, got %d", terms[2].nauxexprs);

   /* check the product expressions */
   cr_expect_eq(terms[0].x, bvar1, "x var of product 0 should be bvar1, got %s", SCIPvarGetName(terms[0].x));
   cr_expect_eq(terms[0].y, bvar2, "y var of product 0 should be bvar2, got %s", SCIPvarGetName(terms[0].y));
   cr_expect_eq(terms[1].x, bvar1, "x var of product 1 should be bvar1, got %s", SCIPvarGetName(terms[1].x));
   cr_expect_eq(terms[1].y, x2, "y var of product 1 should be x2, got %s", SCIPvarGetName(terms[1].y));
   cr_expect_eq(terms[2].x, bvar2, "x var of product 1 should be bvar2, got %s", SCIPvarGetName(terms[2].x));
   cr_expect_eq(terms[2].y, x2, "y var of product 1 should be x2, got %s", SCIPvarGetName(terms[2].y));

   /* check the auxiliary expressions obtained from the implied relation and the implied bound */

   /* expression from two implied bounds on x1: should be bvar1*bvar2 <= -x1 - bvar1 + bvar2 */
   checkAuxExpr(terms[0], 0, x1, (SCIP_Real[3]) {-1.0, -1.0, 1.0}, 0.0, FALSE, TRUE);

   /* from two implied bounds on x1: should be bvar2*bvar1 <= -0.5x1 + 0.5bvar2 */
   checkAuxExpr(terms[0], 1, x1, (SCIP_Real[3]) {-0.5, 0.0, 0.5}, 0.0, FALSE, TRUE);

   /* from implied bound on x1 (with bvar1) and unconditional: should be bvar1*x2 >= x1 + 1.5bvar1 + x2 - 0.5 */
   checkAuxExpr(terms[1], 0, x1, (SCIP_Real[3]) {1.0, 1.5, 1.0}, -0.5, TRUE, FALSE);

   /* from implied bound on x1 (with bvar2) and unconditional: should be bvar2*x2 <= -x1 + 0.5bvar2 */
   checkAuxExpr(terms[2], 0, x1, (SCIP_Real[3]) {-1.0, 0.5, 0.0}, 0.0, FALSE, TRUE);

   /* free memory */
   SCIP_CALL( freeSepaData(scip, sepadata) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIPfreeBuffer(scip, &sepadata);
}

Test(product_detection, reltables, .init = setup, .fini = teardown, .description = "test creating relation tables")
{
   SCIP_CONS* conss[5];
   SCIP_VAR* vars[3];
   SCIP_Real* coefs;
   SCIP_Bool cutoff;
   int nrows;
   int i;
   SCIP_ROW** prob_rows;
   int* row_list;
   SCIP_HASHTABLE* hashtable2;
   SCIP_HASHTABLE* hashtable3;
   SCIP_HASHMAP* vars_in_2rels;
   HASHDATA* foundhashdata;
   ADJACENTVARDATA* adjvardata;

   /* create linear relations */

   vars[0] = x1;
   vars[1] = x2;
   vars[2] = bvar1;

   /* 0 <= 2x1 + 2x2 <= 1 */
   coefs = (SCIP_Real[2]) {2.0, 2.0};
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &conss[0], "c1", 2, vars, coefs, 0.0, 1.0) );

   /* 0 <= 2x1 + 1.5x2 + 3bvar1 <= 1 */
   coefs = (SCIP_Real[3]) {2.0, 1.5, 3.0};
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &conss[1], "c2", 3, vars, coefs, 0.0, 1.0) );

   /* 4x1 + 2x2 - bvar1 <= 1 */
   coefs = (SCIP_Real[3]) {4.0, 2.0, -1.0};
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &conss[2], "c3", 3, vars, coefs, -SCIPinfinity(scip), 1.0) );

   /* x1 + 3x2 + bvar1 <= 2 */
   coefs = (SCIP_Real[3]) {1.0, 3.0, 1.0};
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &conss[3], "c4", 3, vars, coefs, -SCIPinfinity(scip), 2.0) );

   /* x1 - 2x2 <= 1 */
   coefs = (SCIP_Real[2]) {1.0, -2.0};
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &conss[4], "c5", 2, vars, coefs, -SCIPinfinity(scip), 1.0) );

   for( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPaddCons(scip, conss[i]) );
   }

   /* construct the LP */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   /* get the rows */
   nrows = 5;
   SCIP_CALL( getOriginalRows(scip, &prob_rows, &nrows) );

   /* create tables of implied and unconditional relations */
   SCIP_CALL( SCIPhashtableCreate(&hashtable3, SCIPblkmem(scip), nrows, SCIPhashGetKeyStandard,
      hashdataKeyEqConss, hashdataKeyValConss, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&hashtable2, SCIPblkmem(scip), nrows, SCIPhashGetKeyStandard,
      hashdataKeyEqConss, hashdataKeyValConss, (void*) scip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &row_list, nrows) );

   /* allocate the data structure for variables that appear in 2-var relations */
   SCIP_CALL( SCIPhashmapCreate(&vars_in_2rels, SCIPblkmem(scip), 10) );

   /* fill in the relevant data structures */
   SCIP_CALL( fillRelationTables(scip, prob_rows, nrows, hashtable2, hashtable3, vars_in_2rels, row_list) );

   /* check the linked list of rows (should correspond to 3 -> 2 -> 1 and 4 -> 0) */
   cr_expect_eq(row_list[0], -1);
   cr_expect_eq(row_list[1], -1);
   cr_expect_eq(row_list[2], 1);
   cr_expect_eq(row_list[3], 2);
   cr_expect_eq(row_list[4], 0);

   /* check the hashtable for 3-variable relations */
   cr_expect_eq(SCIPhashtableGetNElements(hashtable3), 1, "expected 1 var triple, got %lld",
         SCIPhashtableGetNElements(hashtable3));

   for( i = 0; i < SCIPhashtableGetNEntries(hashtable3); ++i )
   {
      foundhashdata = (HASHDATA*)SCIPhashtableGetEntry(hashtable3, i);
      if( foundhashdata == NULL )
         continue;

      cr_expect_eq(foundhashdata->nvars, 3, "expected 3 vars in an element of hashtable3, got %d", foundhashdata->nvars);
      cr_expect_eq(foundhashdata->vars[0], bvar1, "expected first var in the triple to be bvar1, got %s",
            SCIPvarGetName(foundhashdata->vars[0]));
      cr_expect_eq(foundhashdata->vars[1], x1, "expected second var in the triple to be x1, got %s",
            SCIPvarGetName(foundhashdata->vars[1]));
      cr_expect_eq(foundhashdata->vars[2], x2, "expected third var in the triple to be x2, got %s",
            SCIPvarGetName(foundhashdata->vars[2]));

      cr_expect_eq(foundhashdata->firstrow, 3);
   }

   /* check the hashtable for 2-variable relations */
   cr_expect_eq(SCIPhashtableGetNElements(hashtable2), 1, "expected 1 var pair, got %lld",
         SCIPhashtableGetNElements(hashtable2));

   for( i = 0; i < SCIPhashtableGetNEntries(hashtable2); ++i )
   {
      foundhashdata = (HASHDATA*)SCIPhashtableGetEntry(hashtable2, i);
      if( foundhashdata == NULL )
         continue;

      cr_expect_eq(foundhashdata->nvars, 2, "expected 2 vars in an element of hashtable2, got %d", foundhashdata->nvars);
      cr_expect_eq(foundhashdata->vars[0], x1, "expected first var in the triple to be x1, got %s",
            SCIPvarGetName(foundhashdata->vars[0]));
      cr_expect_eq(foundhashdata->vars[1], x2, "expected second var in the triple to be x2, got %s",
            SCIPvarGetName(foundhashdata->vars[1]));

      SCIPdebugMsg(scip, "(%s, %s): ", SCIPvarGetName(foundhashdata->vars[0]),
      SCIPvarGetName(foundhashdata->vars[1]));

      cr_expect_eq(foundhashdata->firstrow, 4);
   }

   /* check the data structure storing variables participating together in 2-variable relations */
   adjvardata = (ADJACENTVARDATA*) SCIPhashmapGetImage(vars_in_2rels, (void*)(size_t) SCIPvarGetIndex(x1));
   cr_assert(adjvardata != NULL);
   cr_expect_eq(adjvardata->nadjacentvars, 1, "expected 1 vars adjacent to x1, got %d", adjvardata->nadjacentvars);
   cr_expect_eq(adjvardata->adjacentvars[0], x2, "expected x2 to be adjacent to x1, got %s",
        SCIPvarGetName(adjvardata->adjacentvars[0]));

   adjvardata = (ADJACENTVARDATA*) SCIPhashmapGetImage(vars_in_2rels, (void*)(size_t) SCIPvarGetIndex(x2));
   cr_assert(adjvardata != NULL);
   cr_expect_eq(adjvardata->nadjacentvars, 1, "expected 1 vars adjacent to x2, got %d", adjvardata->nadjacentvars);
   cr_expect_eq(adjvardata->adjacentvars[0], x1, "expected x1 to be adjacent to x2, got %s",
        SCIPvarGetName(adjvardata->adjacentvars[0]));

   /* free memory */
   clearVarAdjacency(scip, vars_in_2rels);
   SCIPhashmapFree(&vars_in_2rels);

   SCIPfreeBufferArray(scip, &row_list);

   for( i = 0; i < SCIPhashtableGetNEntries(hashtable2); ++i )
   {
      foundhashdata = (HASHDATA*)SCIPhashtableGetEntry(hashtable2, i);
      if( foundhashdata == NULL )
         continue;

      SCIPfreeBuffer(scip, &foundhashdata);
   }

   for( i = 0; i < SCIPhashtableGetNEntries(hashtable3); ++i )
   {
      foundhashdata = (HASHDATA*)SCIPhashtableGetEntry(hashtable3, i);
      if( foundhashdata == NULL )
         continue;

      SCIPfreeBuffer(scip, &foundhashdata);
   }

   SCIPhashtableFree(&hashtable2);
   SCIPhashtableFree(&hashtable3);

   SCIPfreeBufferArray(scip, &prob_rows);

   for( i = 4; i >= 0; --i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conss[i]) );
   }
}
