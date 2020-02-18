/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bilinhash.c
 * @brief  tests functionalities to access all bilinear terms
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   /* create SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include expression constraint handler */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* store expression constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 2.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

static
void teardown(void)
{
   /* release variables */
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/* define test test suite */
TestSuite(bilinhash, .init = setup, .fini = teardown);

/* tests the creating and release of the hash table using non-API methods from cons_expr.c */
Test(bilinhash, createInsertFree)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create constraint handler data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &conshdlrdata) );

   /* inserts two bilinear terms into the hash table */
   SCIP_CALL( bilinearTermsInsert(scip, conshdlrdata, x, y, NULL) );
   SCIP_CALL( bilinearTermsInsert(scip, conshdlrdata, y, z, NULL) );
   cr_expect(conshdlrdata->nbilinterms == 2);
   cr_expect(conshdlrdata->bilinterms[0].x == x);
   cr_expect(conshdlrdata->bilinterms[0].y == y);
   cr_expect(conshdlrdata->bilinterms[1].x == y);
   cr_expect(conshdlrdata->bilinterms[1].y == z);

   /* free hash table */
   SCIP_CALL( bilinearTermsFree(scip, conshdlrdata) );

   /* free constraint handler data */
   SCIPfreeBlockMemory(scip, &conshdlrdata);
}

/* tests API methods for a simple problem containing two expression constraints */
Test(bilinhash, api_methods)
{
   const char* inputs[2] = {"[expr] <c1>: (<x>[C])^2 + <x>[C] * <y>[C] <= 4;",
      "[expr] <c2>: abs(<y>[C] * <z>[C] + <x>[C] * <y>[C]) * (log(<x>[C] + <z>[C]))^2 <= 1;"};
   SCIP_CONSEXPR_BILINTERM* bilinterms;
   SCIP_VAR* xs[3];
   SCIP_VAR* ys[3];
   SCIP_VAR* auxvars[3];
   SCIP_VAR* auxvar;
   SCIP_VAR* tx;
   SCIP_VAR* ty;
   SCIP_VAR* tz;
   int i;

   /* create, add, and release expression constraints */
   for( i = 0; i < 2; ++i )
   {
      SCIP_CONS* cons;
      SCIP_Bool success;

      SCIP_CALL( SCIPparseCons(scip, &cons, inputs[i],
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      cr_expect(success);
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* go to the solving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* collect transformed variables */
   tx = SCIPvarGetTransVar(x);
   ty = SCIPvarGetTransVar(y);
   tz = SCIPvarGetTransVar(z);

   /* collect all bilinear terms manually because CONSINITLP has not been called yet */
   SCIP_CALL( bilinearTermsInsertAll(scip, conshdlr, SCIPgetConss(scip), SCIPgetNConss(scip)) );

   /*
    * because no auxiliary variables are present, there are only three bilinear terms: xx, xy, yz
    */
   cr_expect(SCIPgetConsExprNBilinTerms(conshdlr) == 3);

   bilinterms = SCIPgetConsExprBilinTerms(conshdlr);
   cr_assert(bilinterms != NULL);
   cr_expect(bilinterms[0].x == tx && bilinterms[0].y == tx && bilinterms[0].auxvar == NULL);
   cr_expect(bilinterms[1].x == tx && bilinterms[1].y == ty && bilinterms[1].auxvar == NULL);
   cr_expect(bilinterms[2].x == ty && bilinterms[2].y == tz && bilinterms[2].auxvar == NULL);

   /* xx exists */
   cr_expect(SCIPgetConsExprBilinTerm(conshdlr, tx, tx) != NULL);

   /* xy exists */
   cr_expect(SCIPgetConsExprBilinTerm(conshdlr, tx, ty) != NULL);

   /* yx = xy exists */
   cr_expect(SCIPgetConsExprBilinTerm(conshdlr, ty, tx) != NULL);

   /* yz exists */
   cr_expect(SCIPgetConsExprBilinTerm(conshdlr, ty, tz) != NULL);

   /* xz does not exist */
   cr_expect(SCIPgetConsExprBilinTerm(conshdlr, tx, tz) == NULL);

   /* zz does not exist */
   cr_expect(SCIPgetConsExprBilinTerm(conshdlr, tz, tz) == NULL);
}