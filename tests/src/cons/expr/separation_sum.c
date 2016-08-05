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

/**@file   separation_sum.c
 * @brief  tests separation of sums
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_sumprod.c"
#include "separation.h"

Test(separation, sum, .init = setup, .fini = teardown,
   .description = "test separation for a sum expression"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* cut;
   int i;

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr, 0, NULL, NULL, 1.5) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, xexpr, 2.3) );
   SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, yexpr, -5.1) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /* compute cut */
   cut = NULL;
   SCIP_CALL( separatePointSum(scip,  conshdlr, expr, &cut) );

   assert(cut != NULL);
   cr_assert_eq(SCIProwGetNNonz(cut), 3);
   cr_assert_eq(SCIProwGetLhs(cut), -1.5);
   cr_assert_eq(SCIProwGetRhs(cut), -1.5);

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_assert_eq(coef, 2.3);
      else if( var == SCIPvarGetTransVar(y) )
         cr_assert_eq(coef, -5.1);
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
