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

/**@file   separation_prod.c
 * @brief  tests separation of products
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_product.c"
#include "separation.h"

Test(separation, bilinear, .init = setup, .fini = teardown,
   .description = "test separation for a bilinear expression"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* cut;
   int i;

   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr, 0, NULL, 1.5) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, xexpr) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, yexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /*
    * compute a cut for w = 1.5*x*y with x* = 0, y* = -4, w* = 1
    * together with the bounds this should result in a cut of the form
    *    w <= -4.5x - 1.5y - 4.5
    */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -4.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 1.0) );

   cut = NULL;
   SCIP_CALL( separatePointProduct(scip,  conshdlr, expr, sol, &cut) );

   cr_assert(cut != NULL);
   cr_assert_eq(SCIProwGetNNonz(cut), 3);
   cr_assert_eq(SCIProwGetLhs(cut), 4.5);
   cr_assert_eq(SCIProwGetRhs(cut), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_assert_eq(coef, -4.5);
      else if( var == SCIPvarGetTransVar(y) )
         cr_assert_eq(coef, -1.5);
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   /*
    * compute a cut for w = 1.5*x*y with x* = 0, y* = -4, w* = -1
    * together with the bounds this should result in a cut of the form
    *    w >= -9x - 1.5y - 9
    */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -4.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -1.0) );

   cut = NULL;
   SCIP_CALL( separatePointProduct(scip,  conshdlr, expr, sol, &cut) );

   cr_assert(cut != NULL);
   cr_assert_eq(SCIProwGetNNonz(cut), 3);
   cr_assert_eq(SCIProwGetLhs(cut), -SCIPinfinity(scip));
   cr_assert_eq(SCIProwGetRhs(cut), 9.0);

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_assert_eq(coef, -9.0);
      else if( var == SCIPvarGetTransVar(y) )
         cr_assert_eq(coef, -1.5);
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
