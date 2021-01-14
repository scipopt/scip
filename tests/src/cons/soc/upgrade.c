/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   upgrade.c
 * @brief  unit test for upgrade to SOC constraint
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;

/* creates scip, problem, includes nonlinear and quadratic cons handlers, and includes NLP */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
}

static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

Test(upgrade, linearbinary, .init = setup, .fini = teardown,
   .description = "check upgrade from quadratic constraint with linear binary variables to SOC"
   )
{
   SCIP_VAR* b1;
   SCIP_VAR* b2;
   SCIP_VAR* b3;
   SCIP_VAR* x;
   SCIP_CONS* e1;
   SCIP_CONS* e1soc;
   int i;
   /*
      minimize
        obj: b1 - b2 - 12 b3 + x
      Subject to
        e1: 0.5 b1 + 2 b2 + b3 + [- x^2] <= 0
      Binary
        b1 b2 b3
      End
   */
   SCIP_CALL( SCIPcreateVarBasic(scip, &b1, "b1", 0.0, 1.0,   1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &b2, "b2", 0.0, 1.0,  -1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &b3, "b3", 0.0, 1.0, -12.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, b1) );
   SCIP_CALL( SCIPaddVar(scip, b2) );
   SCIP_CALL( SCIPaddVar(scip, b3) );
   SCIP_CALL( SCIPaddVar(scip, x) );

   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &e1, "e1", 0, NULL,NULL, 0, NULL, NULL, NULL, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, e1, b1, 0.5) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, e1, b2, 2.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, e1, b3, 1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, e1, x, -1.0) );
   SCIP_CALL( SCIPaddCons(scip, e1) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPprintTransProblem(scip, NULL, NULL, FALSE) );

   /*  [soc] <e1>: sqrt( + (1*(<t_b3>[B]+0))^2 + (1.4142135623731*(<t_b2>[B]+0))^2 ) <= 1*(<t_x>[C]-0); */

   cr_expect(SCIPgetNActiveConss(scip) == 1);
   e1soc = SCIPgetConss(scip)[0];
   cr_expect(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(e1soc)), "soc") == 0);

   cr_expect(SCIPgetRhsVarSOC(scip, e1soc) == SCIPvarGetTransVar(x));
   cr_expect(SCIPgetRhsCoefSOC(scip, e1soc) == 1.0);
   cr_expect(SCIPgetRhsOffsetSOC(scip, e1soc) == 0.0);

   for( i = 0; i < SCIPgetNLhsVarsSOC(scip, e1soc); ++i )
   {
      if( SCIPgetLhsVarsSOC(scip, e1soc)[i] == SCIPvarGetTransVar(b1) )
         cr_expect_float_eq(SCIPgetLhsCoefsSOC(scip, e1soc)[i], sqrt(0.5), 1e-10);
      if( SCIPgetLhsVarsSOC(scip, e1soc)[i] == SCIPvarGetTransVar(b2) )
         cr_expect_eq(SCIPgetLhsCoefsSOC(scip, e1soc)[i], sqrt(2.0));
      if( SCIPgetLhsVarsSOC(scip, e1soc)[i] == SCIPvarGetTransVar(b3) )
         cr_expect_eq(SCIPgetLhsCoefsSOC(scip, e1soc)[i], 1.0);
   }

   SCIP_CALL( SCIPreleaseCons(scip, &e1) );
   SCIP_CALL( SCIPreleaseVar(scip, &b1) );
   SCIP_CALL( SCIPreleaseVar(scip, &b2) );
   SCIP_CALL( SCIPreleaseVar(scip, &b3) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

}
