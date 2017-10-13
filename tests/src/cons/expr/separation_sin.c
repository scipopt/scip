/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   separation_sin.c
 * @brief  tests separation of sin()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_sin.c"
#include "separation.h"

Test(separation, sinus_x, .init = setup, .fini = teardown,
   .description = "test separation for a sinus expression in large range"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* secant;
   SCIP_ROW* ltangent;
   SCIP_ROW* rtangent;
   SCIP_ROW* lmidtangent;
   SCIP_ROW* rmidtangent;
   SCIP_ROW* soltangent;
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_Real newtonpoint;
   int i;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;
   soltangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr, xexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /*
    * test initial separation
    */

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
         TRUE, TRUE) );

   /* check infeasible cuts */
   cr_expect(secant == NULL);
   cr_expect(ltangent == NULL);

   /* check rtangent, should be underestimating */
   cr_assert(rtangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(rtangent), 2);
   cr_expect_eq(SCIProwGetLhs(rtangent), -SCIPinfinity(scip));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(rtangent), 5 * COS(5) - SIN(5)));

   for( i = 0; i < SCIProwGetNNonz(rtangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(rtangent)[i]);
      coef = SCIProwGetVals(rtangent)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_expect(SCIPisEQ(scip, coef, COS(5)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* check lmidtangent, should be overestimating */
   cr_assert(lmidtangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(lmidtangent), 2);
   newtonpoint = 0.4936608602;
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(lmidtangent), -COS(newtonpoint) - SIN(-1)));
   cr_expect_eq(SCIProwGetRhs(lmidtangent), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(lmidtangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(lmidtangent)[i]);
      coef = SCIProwGetVals(lmidtangent)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_expect(SCIPisEQ(scip, coef, COS(newtonpoint)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* check rmidtangent, should be overestimating */
   cr_assert(rmidtangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(rmidtangent), 2);
   cr_expect_eq(SCIProwGetLhs(rmidtangent), -SCIPinfinity(scip));
   newtonpoint = 2.2544608804;
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(rmidtangent), 5 * COS(newtonpoint) - SIN(5)));

   for( i = 0; i < SCIProwGetNNonz(rmidtangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(rmidtangent)[i]);
      coef = SCIProwGetVals(rmidtangent)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_expect(SCIPisEQ(scip, coef, COS(newtonpoint)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &rtangent) );
   SCIP_CALL( SCIPreleaseRow(scip, &lmidtangent) );
   SCIP_CALL( SCIPreleaseRow(scip, &rmidtangent) );

   /*
    * test overestimation
    */

   /* create solution that induces overestimating cut */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      TRUE, FALSE) );

   /* check infeasible cuts */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   cr_assert(soltangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(soltangent), 2);
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(soltangent), 1.5 * COS(1.5) - SIN(1.5)));
   cr_expect_eq(SCIProwGetRhs(soltangent), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(soltangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(soltangent)[i]);
      coef = SCIProwGetVals(soltangent)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_expect(SCIPisEQ(scip, coef, COS(1.5)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &soltangent) );

   /*
    * test underestimation
    */

   /* create solution that induces underestimating cut */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 4.8) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      FALSE, TRUE) );

   /* check infeasible cuts */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent */
   cr_assert(soltangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(soltangent), 2);
   cr_expect_eq(SCIProwGetLhs(soltangent), -SCIPinfinity(scip));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(soltangent), 4.8 * COS(4.8) - SIN(4.8)));

   for( i = 0; i < SCIProwGetNNonz(soltangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(soltangent)[i]);
      coef = SCIProwGetVals(soltangent)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_expect(SCIPisEQ(scip, coef, COS(4.8)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &soltangent) );

   /*
    * test point where solution tangent is not feasible
    */

   /* create solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 4.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      FALSE, TRUE) );

   /* check infeasible cuts */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);
   cr_expect(soltangent == NULL);


   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(separation, sinus_y, .init = setup, .fini = teardown,
   .description = "test separation for a sinus expression in mid size range"
)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* secant;
   SCIP_ROW* ltangent;
   SCIP_ROW* rtangent;
   SCIP_ROW* lmidtangent;
   SCIP_ROW* rmidtangent;
   SCIP_ROW* soltangent;
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_Real newtonpoint;
   int i;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;
   soltangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr, yexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /*
    * test initial separation
    */

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      TRUE, TRUE) );

   /* check infeasible cuts */
   cr_expect(secant == NULL);
   cr_expect(rtangent == NULL);

   /* check ltangent, should be overestimating */
   cr_assert(ltangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(ltangent), 2);
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(ltangent), (-6) * COS(-6) - SIN(-6)));
   cr_expect_eq(SCIProwGetRhs(ltangent), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(ltangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(ltangent)[i]);
      coef = SCIProwGetVals(ltangent)[i];

      if( var == SCIPvarGetTransVar(y) )
         cr_expect(SCIPisEQ(scip, coef, COS(-6)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* check lmidtangent, should be underestimating */
   cr_assert(lmidtangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(lmidtangent), 2);
   cr_expect_eq(SCIProwGetLhs(lmidtangent), -SCIPinfinity(scip));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(lmidtangent), 2.0 * SIN(-3) - SIN(-6)));

   for( i = 0; i < SCIProwGetNNonz(lmidtangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(lmidtangent)[i]);
      coef = SCIProwGetVals(lmidtangent)[i];

      if( var == SCIPvarGetTransVar(y) )
         cr_expect(SCIPisEQ(scip, coef, (SIN(-3) - SIN(-6)) / 3.0));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* check rmidtangent, should be overestimating */
   cr_assert(rmidtangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(rmidtangent), 2);
   newtonpoint = -3.2123712333;
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(rmidtangent), (-3) * COS(newtonpoint) - SIN(-3)));
   cr_expect_eq(SCIProwGetRhs(rmidtangent), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(rmidtangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(rmidtangent)[i]);
      coef = SCIProwGetVals(rmidtangent)[i];

      if( var == SCIPvarGetTransVar(y) )
         cr_expect(SCIPisEQ(scip, coef, COS(newtonpoint)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &ltangent) );
   SCIP_CALL( SCIPreleaseRow(scip, &lmidtangent) );
   SCIP_CALL( SCIPreleaseRow(scip, &rmidtangent) );

   /*
    * test overestimation
    */

   /* create solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -4.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      TRUE, FALSE) );

   /* check infeasible cuts */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent, should be overestimating */
   cr_assert(soltangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(soltangent), 2);
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(soltangent), (-4) * COS(-4) - SIN(-4)));
   cr_expect_eq(SCIProwGetRhs(soltangent), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(soltangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(soltangent)[i]);
      coef = SCIProwGetVals(soltangent)[i];

      if( var == SCIPvarGetTransVar(y) )
         cr_expect(SCIPisEQ(scip, coef, COS(-4)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &soltangent) );

   /*
    * test underestimation (not possible in this case)
    */

   /* create solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -3.1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      FALSE, TRUE) );

   /* check infeasible cuts */
   cr_expect(secant == NULL);
   cr_expect(rmidtangent == NULL);
   cr_expect(soltangent == NULL);

   /* check lmidtangent, should be underestimating */
   cr_assert(lmidtangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(lmidtangent), 2);
   cr_expect_eq(SCIProwGetLhs(lmidtangent), -SCIPinfinity(scip));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(lmidtangent), 2.0 * SIN(-3) - SIN(-6)));

   for( i = 0; i < SCIProwGetNNonz(lmidtangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(lmidtangent)[i]);
      coef = SCIProwGetVals(lmidtangent)[i];

      if( var == SCIPvarGetTransVar(y) )
         cr_expect(SCIPisEQ(scip, coef, (SIN(-3) - SIN(-6)) / 3.0));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &lmidtangent) );

   /*
    * test point where solution tangent is infeasible
    */

   /* create solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -3.2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      TRUE, FALSE) );

   /* check infeasible cuts */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(soltangent == NULL);

   /* check rmidtangent, should be overestimating */
   cr_assert(rmidtangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(rmidtangent), 2);
   newtonpoint = -3.2123712333;
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(rmidtangent), (-3) * COS(newtonpoint) - SIN(-3)));
   cr_expect_eq(SCIProwGetRhs(rmidtangent), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(rmidtangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(rmidtangent)[i]);
      coef = SCIProwGetVals(rmidtangent)[i];

      if( var == SCIPvarGetTransVar(y) )
         cr_expect(SCIPisEQ(scip, coef, COS(newtonpoint)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &rmidtangent) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

}

Test(separation, sinus_z, .init = setup, .fini = teardown,
   .description = "test separation for a sinus expression in short range"
)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* secant;
   SCIP_ROW* ltangent;
   SCIP_ROW* rtangent;
   SCIP_ROW* lmidtangent;
   SCIP_ROW* rmidtangent;
   SCIP_ROW* soltangent;
   SCIP_VAR* var;
   SCIP_Real coef;
   int i;

   secant = NULL;
   ltangent = NULL;
   rtangent = NULL;
   lmidtangent = NULL;
   rmidtangent = NULL;
   soltangent = NULL;

   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr, zexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /*
    * test initial separation
    */

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, &ltangent, &rtangent, &lmidtangent, &rmidtangent, NULL,
      TRUE, TRUE) );

   /* check infeasible cuts */
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check secant, should be underestimating */
   cr_assert(secant != NULL);
   cr_expect_eq(SCIProwGetNNonz(secant), 2);
   cr_expect_eq(SCIProwGetLhs(secant), -SCIPinfinity(scip));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(secant), 0.5 * SIN(3) - 1.5 * SIN(1)));

   for( i = 0; i < SCIProwGetNNonz(secant); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(secant)[i]);
      coef = SCIProwGetVals(secant)[i];

      if( var == SCIPvarGetTransVar(z) )
         cr_expect(SCIPisEQ(scip, coef, 0.5 * (SIN(3) - SIN(1))));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* check ltangent, should be overestimating */
   cr_assert(ltangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(ltangent), 2);
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(ltangent), COS(1) - SIN(1)));
   cr_expect_eq(SCIProwGetRhs(ltangent), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(ltangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(ltangent)[i]);
      coef = SCIProwGetVals(ltangent)[i];

      if( var == SCIPvarGetTransVar(z) )
         cr_expect(SCIPisEQ(scip, coef, COS(1)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* check rtangent, should be overestimating */
   cr_assert(rtangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(rtangent), 2);
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(rtangent), 3 * COS(3) - SIN(3)));
   cr_expect_eq(SCIProwGetRhs(rtangent), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(rtangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(rtangent)[i]);
      coef = SCIProwGetVals(rtangent)[i];

      if( var == SCIPvarGetTransVar(z) )
         cr_expect(SCIPisEQ(scip, coef, COS(3)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &secant) );
   SCIP_CALL( SCIPreleaseRow(scip, &ltangent) );
   SCIP_CALL( SCIPreleaseRow(scip, &rtangent) );

   /*
    * test overestimation
    */

   /* create solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 2.0) );

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      TRUE, FALSE) );

   /* check infeasible cuts */
   cr_expect(secant == NULL);
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);

   /* check soltangent, should be overestimating */
   cr_assert(soltangent != NULL);
   cr_expect_eq(SCIProwGetNNonz(soltangent), 2);
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(soltangent), 2 * COS(2) - SIN(2)));
   cr_expect_eq(SCIProwGetRhs(soltangent), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(soltangent); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(soltangent)[i]);
      coef = SCIProwGetVals(soltangent)[i];

      if( var == SCIPvarGetTransVar(z) )
         cr_expect(SCIPisEQ(scip, coef, COS(2)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &soltangent) );

   /* note: underestimation doesn't make sense in [1,3] */

   /*
    * test point where solution tangent is infeasible (underestimation)
    */

   /* create solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -2.0) );

   SCIP_CALL( computeCutsSin(scip, conshdlr, expr, sol, &secant, NULL, NULL, &lmidtangent, &rmidtangent, &soltangent,
      FALSE, TRUE) );

   /* check infeasible cuts */
   cr_expect(lmidtangent == NULL);
   cr_expect(rmidtangent == NULL);
   cr_expect(soltangent == NULL);

   /* check secant, should be underestimating */
   cr_assert(secant != NULL);
   cr_expect_eq(SCIProwGetNNonz(secant), 2);
   cr_expect_eq(SCIProwGetLhs(secant), -SCIPinfinity(scip));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(secant), 0.5 * SIN(3) - 1.5 * SIN(1)));

   for( i = 0; i < SCIProwGetNNonz(secant); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(secant)[i]);
      coef = SCIProwGetVals(secant)[i];

      if( var == SCIPvarGetTransVar(z) )
         cr_expect(SCIPisEQ(scip, coef, 0.5 * (SIN(3) - SIN(1))));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_expect_eq(coef, -1.0);
      else
         cr_expect(FALSE, "found an unknown variable");
   }

   /* release cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &secant) );
}