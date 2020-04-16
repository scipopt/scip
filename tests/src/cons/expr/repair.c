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

/**@file   conscheck.c
 * @brief  tests repair mechanism in CONSCHECK
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"

/* we include the source file here to call the trysol heuristic manually */
#include "scip/heur_trysol.c"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );
   SCIP_CALL( SCIPincludeHeurTrySol(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 5.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 5.0, 3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

Test(repair, linvars1, .init = setup, .fini = teardown)
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   SCIP_VAR* var;
   SCIP_Real coef;
   const char* input = "[expr] <test>: 1.1*<x>*<y> + 3.2*<x>^2*<y>^(-5) + 0.5*<z> <= 2;";

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   cr_assert(SCIPgetNConss(scip) == 1);
   consexpr = SCIPgetConss(scip)[0];
   cr_assert(consexpr != NULL);

   SCIP_CALL( SCIPgetLinvarMayDecreaseExpr(scip, conshdlr, consexpr, &var, &coef) );
   cr_expect(var == SCIPvarGetTransVar(z));
   cr_expect(coef == 0.5);
   SCIP_CALL( SCIPgetLinvarMayIncreaseExpr(scip, conshdlr, consexpr, &var, &coef) );
   cr_expect(var == SCIPvarGetTransVar(z));
   cr_expect(coef == 0.5);
}

Test(repair, linvars2, .init = setup, .fini = teardown)
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   SCIP_VAR* var;
   SCIP_Real coef;
   int i;

   const char* input[2] = {
         "[expr] <test>: 1.1*<x>*<y> + 3.2*<x>^2*<y>^(-5) - 0.5*<z> <= 2;",
         "[expr] <test>: 1.1*<x>*<y> + 3.2*<x>^2*<y>^(-5) - 1.5*<z> <= 2;"
   };

   /* parse constraints */
   for( i = 0; i < 2; ++i )
   {
      success = FALSE;
      SCIP_CALL( SCIPparseCons(scip, &consexpr, input[i],
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      cr_assert(success);
      SCIP_CALL( SCIPaddCons(scip, consexpr) );
      SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
   }

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );
   cr_assert(SCIPgetNConss(scip) == 2);

   for( i = 0; i < 2; ++i )
   {
      consexpr = SCIPgetConss(scip)[i];
      cr_assert(consexpr != NULL);

      SCIP_CALL( SCIPgetLinvarMayDecreaseExpr(scip, conshdlr, consexpr, &var, &coef) );
      cr_expect(var == NULL);
      cr_expect(coef == 0.0);

      SCIP_CALL( SCIPgetLinvarMayIncreaseExpr(scip, conshdlr, consexpr, &var, &coef) );
      cr_expect(var == SCIPvarGetTransVar(z));
      cr_expect(i == 0 ? coef == -0.5 : coef == -1.5);
   }
}

Test(repair, sol, .init = setup, .fini = teardown)
{
   SCIP_CONS* consexpr;
   SCIP_SOL* sol;
   SCIP_Bool success;

   const char* input = "[expr] <test>: 1.1*<x>*<y> + <x>^2*<y> + <z> == 2;";
   SCIP_RESULT result;

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );

   /* go to solving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(x), 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(y), 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, SCIPvarGetTransVar(z), 0.0) );

   SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_expect(!success);

   /* call the execution method of trysol manually */
   SCIP_CALL( heurExecTrySol(scip, SCIPfindHeur(scip, "trysol"), SCIP_HEURTIMING_DURINGLPLOOP, FALSE, &result) );
   cr_assert(result == SCIP_FOUNDSOL);
   cr_assert(SCIPgetNSols(scip) == 1);

   sol = SCIPgetBestSol(scip);
   cr_assert(sol != NULL);
   cr_expect(SCIPgetSolVal(scip, sol, SCIPvarGetTransVar(x)) == 0.0);
   cr_expect(SCIPgetSolVal(scip, sol, SCIPvarGetTransVar(y)) == 0.0);
   cr_expect(SCIPgetSolVal(scip, sol, SCIPvarGetTransVar(z)) == 2.0);
}
