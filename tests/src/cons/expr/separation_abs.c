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

/**@file   unittest-cons_epxr.c
 * @brief  unit test for cons_epxr methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_abs.c"
#include "separation.h"


Test(separation, absolute, .init = setup, .fini = teardown,
   .description = "test separation for an absolute expression"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* rowneg;
   SCIP_ROW* rowpos;
   SCIP_ROW* secant;
   SCIP_VAR* var;
   SCIP_Real coef;
   int i;

   rowneg = NULL;
   rowpos = NULL;
   secant = NULL;

   SCIP_CALL( SCIPcreateConsExprExprAbs(scip, conshdlr, &expr, xexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /* compute all possible cuts */
   SCIP_CALL( computeCutsAbs(scip, conshdlr, expr, &rowneg, &rowpos, &secant) );

   /* check left tangent */
   cr_assert(rowneg != NULL);
   cr_assert_eq(SCIProwGetNNonz(rowneg), 2);
   cr_assert_eq(SCIProwGetLhs(rowneg), -SCIPinfinity(scip));
   cr_assert_eq(SCIProwGetRhs(rowneg), 0.0);

   for( i = 0; i < SCIProwGetNNonz(rowneg); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(rowneg)[i]);
      coef = SCIProwGetVals(rowneg)[i];

      if( var == SCIPvarGetTransVar(x) || var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }

   /* check right tangent */
   cr_assert(rowpos != NULL);
   cr_assert_eq(SCIProwGetNNonz(rowpos), 2);
   cr_assert_eq(SCIProwGetLhs(rowpos), 0.0);
   cr_assert_eq(SCIProwGetRhs(rowpos), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(rowpos); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(rowpos)[i]);
      coef = SCIProwGetVals(rowpos)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_assert_eq(coef, 1.0);
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }

   /* check secant */
   cr_assert(secant != NULL);
   cr_assert_eq(SCIProwGetNNonz(secant), 2);
   cr_assert(SCIPisEQ(scip, SCIProwGetLhs(secant), -5.0/3.0));
   cr_assert_eq(SCIProwGetRhs(secant), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(secant); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(secant)[0]);
      coef = SCIProwGetVals(secant)[0];

      if( var == SCIPvarGetTransVar(x) )
         cr_assert(SCIPisEQ(scip, coef, 2.0 / 3.0));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }

   /* release all cuts */
   SCIP_CALL( SCIPreleaseRow(scip, &rowneg) );
   SCIP_CALL( SCIPreleaseRow(scip, &rowpos) );
   SCIP_CALL( SCIPreleaseRow(scip, &secant) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
