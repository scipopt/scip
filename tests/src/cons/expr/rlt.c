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

/**@file   rlt.c
 * @brief  tests rlt functionalities
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr.c"
#include "scip/sepa_rlt.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SEPA* sepa;
static SCIP_SEPADATA* sepadata;
static SCIP_CONSEXPR_EXPR* expr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* xx;
static SCIP_VAR* xy;
static SCIP_VAR* xz;
static SCIP_VAR* sumvar;
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
   const char* input1 = "[expr] <test1>: (<x>[C])^2 + <x>[C] * <y>[C] <= 4;";
   const char* input2 = "[expr] <test2>: abs(<y>[C] * <x>[C]) * (log(<x>[C] * <z>[C]))^2 <= 1;";

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr*/
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

   /* include separator and get it */
   SCIP_CALL( SCIPincludeSepaRlt(scip) );
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

   /* go to the presolving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* create auxiliary variables for all expressions */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, SCIPconshdlrGetConss(conshdlr), SCIPconshdlrGetNConss(conshdlr), &infeasible) );
   assert(!infeasible);

   /* store all bilinear terms in the data of the expression constraint handler */
   SCIP_CALL( SCIPcollectConsExprBilinTerms(scip, conshdlr, SCIPconshdlrGetConss(conshdlr), SCIPconshdlrGetNConss(conshdlr)) );

   /* create sepadata */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );

   /* setup some sepadata */
   sepadata->conshdlr = conshdlr;
   sepadata->maxusedvars = 3;
   sepadata->maxunknownterms = 1;
   sepadata->onlyinitial = FALSE;
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
   expr = SCIPconsGetData(SCIPconshdlrGetConss(conshdlr)[0])->expr;
   sumvar = SCIPgetConsExprExprAuxVar(expr);
   xx = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[0]);
   xy = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[1]);
   expr = SCIPconsGetData(SCIPconshdlrGetConss(conshdlr)[1])->expr;
   prodvar = SCIPgetConsExprExprAuxVar(expr);
   absvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[0]);
   powvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[1]);
   logvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(SCIPgetConsExprExprChildren(expr)[1])[0]);
   xz = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(SCIPgetConsExprExprChildren(SCIPgetConsExprExprChildren(expr)[1])[0])[0]);
}

static
void teardown(void)
{
   /* release sepadata */
   SCIP_CALL( freeSepaData(scip, sepadata) );
   SCIPfreeBlockMemory(scip, &sepadata);

   /* free rlt stuff */
   SCIP_CALL( freeAuxVars(scip, conshdlr, SCIPconshdlrGetConss(conshdlr), SCIPconshdlrGetNConss(conshdlr)) );

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(rlt, .init = setup, .fini = teardown);

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

/* helper method to check whether a bilinear term appears in the problem */
static
SCIP_VAR* getBilinVar(
   SCIP_VAR*             x_,                 /**< first variable */
   SCIP_VAR*             y_                  /**< second variable */
   )
{
   SCIP_CONSEXPR_BILINTERM* bilinterm;

   bilinterm = SCIPgetConsExprBilinTerm(conshdlr, x_, y_);

   return bilinterm == NULL ? NULL : bilinterm->auxvar;
}

Test(rlt, collect)
{
   /* check original variables */
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, x, x)->auxvar, xx);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, x, y)->auxvar, xy);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, x, z)->auxvar, xz);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, y, x)->auxvar, xy);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, z, x)->auxvar, xz);

   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, y, z), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, z, y), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, y, y), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, z, z), NULL);

   /* check auxiliary variables for second constraint */
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, logvar, logvar)->auxvar, powvar);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, absvar, powvar)->auxvar, prodvar);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, prodvar, prodvar), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, prodvar, absvar), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, prodvar, powvar), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, prodvar, logvar), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, absvar, absvar), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, absvar, logvar), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, powvar, powvar), NULL);
   cr_expect_eq(SCIPgetConsExprBilinTerm(conshdlr, powvar, logvar), NULL);
}

Test(rlt, separation)
{
   SCIP_ROW* row1;
   SCIP_ROW* cutlhs;
   SCIP_ROW* cutrhs;
   SCIP_VAR** cutvars;
   SCIP_Real* cutvals;
   SCIP_Bool result;
   SCIP_Bool success;
   int currentnunknown;

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
   SCIP_CALL( isAcceptableRow(scip, sepadata, row1, x, &currentnunknown, &result, NULL) );
   cr_expect(result);
   cr_expect_eq(computeRltCuts(scip, sepa, sepadata, &cutlhs, row1, NULL, x, &success, TRUE, TRUE, TRUE, FALSE), SCIP_OKAY);
   cr_assert(success);
   cr_expect_eq(computeRltCuts(scip, sepa, sepadata, &cutrhs, row1, NULL, x, &success, TRUE, FALSE, TRUE, FALSE), SCIP_OKAY);
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
   cr_expect_eq(computeRltCuts(scip, sepa, sepadata, &cutlhs, row1, NULL, x, &success, FALSE, TRUE, TRUE, FALSE), SCIP_OKAY);
   cr_assert(success);
   cr_expect_eq(computeRltCuts(scip, sepa, sepadata, &cutrhs, row1, NULL, x, &success, FALSE, FALSE, TRUE, FALSE), SCIP_OKAY);
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
   SCIP_CALL( isAcceptableRow(scip, sepadata, row1, y, &currentnunknown, &result, NULL) );
   cr_expect(!result);
   SCIP_CALL( SCIPreleaseRow(scip, &row1) );

}