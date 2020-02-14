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

   /* create hash table */
   SCIP_CALL( bilinearHashTableCreate(scip, conshdlrdata) );

   /* inserts two bilinear terms into the hash table */
   SCIP_CALL( bilinearHashInsert(scip, conshdlrdata, x, y, NULL) );
   SCIP_CALL( bilinearHashInsert(scip, conshdlrdata, y, z, NULL) );
   cr_expect(conshdlrdata->nbilinentries == 2);
   cr_expect(conshdlrdata->bilinentries[0]->x == x);
   cr_expect(conshdlrdata->bilinentries[0]->y == y);
   cr_expect(conshdlrdata->bilinentries[1]->x == y);
   cr_expect(conshdlrdata->bilinentries[1]->y == z);

   /* free hash table */
   SCIP_CALL( bilinearHashTableFree(scip, conshdlrdata) );

   /* free constraint handler data */
   SCIPfreeBlockMemory(scip, &conshdlrdata);
}

/* tests API methods for a simple problem containing two expression constraints */
Test(bilinhash, api_methods)
{
   
}