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
   SCIP_ROW* cut;
   int i;

   SCIP_CALL( SCIPcreateConsExprExprAbs(scip, conshdlr, &expr, xexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /* compute a cut for which we need an overestimation (secant) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 3.0) );

   cut = NULL;
   SCIP_CALL( separatePointAbs(scip, conshdlr, expr, sol, &cut) );

   printf("%e %e\n", SCIProwGetLhs(cut), SCIProwGetRhs(cut));

   cr_assert(cut != NULL);
   cr_assert_eq(SCIProwGetNNonz(cut), 2);
   cr_assert(SCIPisEQ(scip, SCIProwGetLhs(cut), -5.0 / 3.0));
   cr_assert_eq(SCIProwGetRhs(cut), SCIPinfinity(scip));


   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_assert(SCIPisEQ(scip, coef, 2.0 / 3.0));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   /* compute a cut for which we need an underestimation (linearization) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -10.0) );

   cut = NULL;
   SCIP_CALL( separatePointAbs(scip, conshdlr, expr, sol, &cut) );

   cr_assert(cut != NULL);
   cr_assert_eq(SCIProwGetNNonz(cut), 2);
   cr_assert_eq(SCIProwGetLhs(cut), -SCIPinfinity(scip));
   cr_assert(SCIPisEQ(scip, SCIProwGetRhs(cut), 0.0));

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(x) )
         cr_assert_eq(coef, -1.0);
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
