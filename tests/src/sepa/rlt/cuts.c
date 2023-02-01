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

/**@file   cuts.c
 * @brief  tests rlt cut functionalities
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_nonlinear.h"
#include "scip/sepa_rlt.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SEPA* sepa;
static SCIP_SEPADATA* sepadata;
static SCIP_EXPR* expr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* xx;
static SCIP_VAR* xy;
static SCIP_VAR* xz;
static SCIP_VAR* prodvar;
static SCIP_VAR* absvar;
static SCIP_VAR* powvar;
static SCIP_VAR* logvar;

static
void setup(void)
{
   SCIP_CONS* conss[2];
   SCIP_VAR* xo;
   SCIP_VAR* yo;
   SCIP_VAR* zo;
   SCIP_Bool success;
   SCIP_Bool infeasible;
   const char* input1 = "[nonlinear] <test1>: (<x>[C])^2 + <x>[C] * <y>[C] <= 4;";
   const char* input2 = "[nonlinear] <test2>: abs(<y>[C] * <x>[C]) * (log(<x>[C] * <z>[C]))^2 <= 1;";

   SCIP_CALL( SCIPcreate(&scip) );

   /* includes expression handlers as well as nonlinear constraint handler and rlt separator */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* get nonlinear conshdlr*/
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert(conshdlr != NULL);

   /* get separator */
   sepa = SCIPfindSepa(scip, "rlt");
   cr_assert(sepa != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &xo, "x", 0.0, 2.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yo, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zo, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, xo) );
   SCIP_CALL( SCIPaddVar(scip, yo) );
   SCIP_CALL( SCIPaddVar(scip, zo) );

   /* add nonlinear constraints */
   SCIP_CALL( SCIPparseCons(scip, &conss[0], input1,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success);
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &conss[1], input2,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success);
   SCIP_CALL( SCIPaddCons(scip, conss[0]) );
   SCIP_CALL( SCIPaddCons(scip, conss[1]) );
   SCIP_CALL( SCIPreleaseCons(scip, &conss[1]) );
   SCIP_CALL( SCIPreleaseCons(scip, &conss[0]) );

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* go to the solving stage (calls detectnlhdlrs) */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* initialize LP (creates auxvars and collects bilinear terms) */
   SCIP_CALL( SCIPconstructLP(scip, &infeasible) );
   cr_assert(!infeasible);

   /* create sepadata */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );

   /* setup some sepadata */
   sepadata->conshdlr = conshdlr;
   sepadata->maxusedvars = 3;
   sepadata->maxunknownterms = 1;
   sepadata->onlyoriginal = FALSE;
   sepadata->onlycontrows = TRUE;
   sepadata->onlyeqrows = FALSE;

   /* collect the data for RLT */
   SCIP_CALL( createSepaData(scip, sepadata) );

   /* get transformed vars and release vars */
   SCIP_CALL( SCIPgetTransformedVar(scip, xo, &x) );
   SCIP_CALL( SCIPgetTransformedVar(scip, yo, &y) );
   SCIP_CALL( SCIPgetTransformedVar(scip, zo, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &xo) );
   SCIP_CALL( SCIPreleaseVar(scip, &yo) );
   SCIP_CALL( SCIPreleaseVar(scip, &zo) );

   /* collect auxvars */
   expr = SCIPgetExprNonlinear(SCIPconshdlrGetConss(conshdlr)[0]);
   xx = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[0]);
   xy = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[1]);
   expr = SCIPgetExprNonlinear(SCIPconshdlrGetConss(conshdlr)[1]);
   prodvar = SCIPgetExprAuxVarNonlinear(expr);
   absvar = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[0]);
   powvar = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[1]);
   logvar = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(SCIPexprGetChildren(expr)[1])[0]);
   xz = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(SCIPexprGetChildren(SCIPexprGetChildren(expr)[1])[0])[0]);
}

static
void teardown(void)
{
   /* release sepadata */
   SCIP_CALL( freeSepaData(scip, sepadata) );
   SCIPfreeBlockMemory(scip, &sepadata);

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(cuts, .init = setup, .fini = teardown);

static
void checkCut(SCIP_ROW* cut, SCIP_VAR** vars, SCIP_Real* vals, int nvars, SCIP_Real lhs, SCIP_Real rhs)
{
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_Bool found;
   int i;
   int j;

   cr_assert(cut != NULL);
   cr_expect_eq(SCIProwGetNNonz(cut), nvars, "\nExpected %d nonz, got %d", nvars, SCIProwGetNNonz(cut));
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(cut), lhs));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(cut), rhs));

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];
      found = FALSE;

      for( j = 0; j < nvars; ++j )
      {
         if( var == vars[j] )
         {
            cr_expect(SCIPisEQ(scip, coef, vals[j]));
            found = TRUE;
         }
      }

      if( !found )
         cr_expect(FALSE, "found an unknown variable");
   }
}

/* checks the correspondence between variables and expressions */
Test(cuts, collect)
{
   /* check original variables */
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, x, x)->aux.var, xx);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, x, y)->aux.var, xy);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, x, z)->aux.var, xz);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, y, x)->aux.var, xy);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, z, x)->aux.var, xz);

   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, y, z), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, z, y), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, y, y), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, z, z), NULL);

   /* check auxiliary variables for second constraint */
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, logvar, logvar)->aux.var, powvar);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, absvar, powvar)->aux.var, prodvar);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, prodvar, prodvar), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, prodvar, absvar), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, prodvar, powvar), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, prodvar, logvar), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, absvar, absvar), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, absvar, logvar), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, powvar, powvar), NULL);
   cr_expect_eq(SCIPgetBilinTermNonlinear(conshdlr, powvar, logvar), NULL);
}

/* computes and checks cuts */
Test(cuts, separation)
{
   SCIP_ROW* row1;
   SCIP_ROW* cutlhs;
   SCIP_ROW* cutrhs;
   SCIP_VAR** cutvars;
   SCIP_Real* cutvals;
   SCIP_Bool result;
   SCIP_Bool success;
   int currentnunknown;
   int* bestunder;
   int* bestover;
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &bestunder, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestover, 3) );

   for( i = 0; i < 3; ++i )
   {
      bestunder[i] = -1;
      bestover[i] = -1;
   }

   /* create test row1: -10 <= 4x - 7y + z <= 5 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row1, "test_row", -10.0, 5.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, x, 4.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, y, -7.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, z, 1.0) );

   cutlhs = NULL;
   cutrhs = NULL;

   success = TRUE;
   /*
    * cut for row1 and (x-0)
    */
   SCIP_CALL( isAcceptableRow(sepadata, row1, x, &currentnunknown, &result) );
   cr_expect(result);
   cr_expect_eq(computeRltCut(scip, sepa, sepadata, &cutlhs, row1, NULL, NULL, bestunder, bestover, x, &success, TRUE, TRUE,
         TRUE, FALSE, FALSE), SCIP_OKAY);
   cr_assert(success);
   cr_expect_eq(computeRltCut(scip, sepa, sepadata, &cutrhs, row1, NULL, NULL, bestunder, bestover, x, &success, TRUE, FALSE,
         TRUE, FALSE, FALSE), SCIP_OKAY);
   cr_assert(success);
   cr_assert(cutlhs != NULL);
   cr_assert(cutrhs != NULL);

   /* check lhs cut */
   cutvars = (SCIP_VAR*[4]) {xx, xy, xz, x};
   cutvals = (SCIP_Real[4]) {4.0, -7.0, 1.0, 10.0};
   checkCut(cutlhs, cutvars, cutvals, 4, 0.0, SCIPinfinity(scip));

   /* check rhs cut */
   cutvals = (SCIP_Real[4]) {4.0, -7.0, 1.0, -5.0};
   checkCut(cutrhs, cutvars, cutvals, 4, -SCIPinfinity(scip), 0.0);

   SCIP_CALL( SCIPreleaseRow(scip, &cutlhs) );
   SCIP_CALL( SCIPreleaseRow(scip, &cutrhs) );

   /*
    * cut for row1 and (2-x)
    */
   cr_expect_eq(computeRltCut(scip, sepa, sepadata, &cutlhs, row1, NULL, NULL, bestunder, bestover, x, &success, FALSE,
         TRUE, TRUE, FALSE, FALSE), SCIP_OKAY);
   cr_assert(success);
   cr_expect_eq(computeRltCut(scip, sepa, sepadata, &cutrhs, row1, NULL, NULL, bestunder, bestover, x, &success, FALSE,
         FALSE, TRUE, FALSE, FALSE), SCIP_OKAY);
   cr_assert(success);
   cr_assert(cutlhs != NULL);
   cr_assert(cutrhs != NULL);

   /* check lhs cut */
   cutvars = (SCIP_VAR*[6]) {xx, xy, xz, x, y, z};
   cutvals = (SCIP_Real[6]) {-4.0, 7.0, -1.0, -2.0, -14.0, 2.0};
   checkCut(cutlhs, cutvars, cutvals, 6, -20.0, SCIPinfinity(scip));

   /* check rhs cut */
   cutvals = (SCIP_Real[6]) {-4.0, 7.0, -1.0, 13.0, -14.0, 2.0};
   checkCut(cutrhs, cutvars, cutvals, 6, -SCIPinfinity(scip), 10.0);

   SCIP_CALL( SCIPreleaseRow(scip, &cutlhs) );
   SCIP_CALL( SCIPreleaseRow(scip, &cutrhs) );

   /* check for not acceptable row */
   SCIP_CALL( isAcceptableRow(sepadata, row1, y, &currentnunknown, &result) );
   cr_expect(!result);
   SCIP_CALL( SCIPreleaseRow(scip, &row1) );

   SCIPfreeBufferArray(scip, &bestover);
   SCIPfreeBufferArray(scip, &bestunder);
}
