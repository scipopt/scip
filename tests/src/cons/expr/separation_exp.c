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

/**@file   separation_exp.c
 * @brief  tests estimators of exp()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_exp.c"
#include "separation.h"

Test(separation, exponential, .init = setup, .fini = teardown,
   .description = "test separation for an exponential expression"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Real coef;
   SCIP_Real constant;
   SCIP_Bool local;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreateConsExprExprExp(scip, conshdlr, &expr, xexpr) );

   /* compute an overestimation (secant) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );

   SCIP_CALL( estimateExp(scip, conshdlr, expr, sol, TRUE, -SCIPinfinity(scip), &coef, &constant, &local, &success, NULL) );

   cr_assert(success);
   cr_assert_float_eq(constant, (exp(5) + 5 * exp(-1)) / 6.0, SCIPepsilon(scip));
   cr_assert_float_eq(coef, (exp(5) - exp(-1)) / 6.0, SCIPepsilon(scip));
   cr_assert(local);

   /* compute an underestimation (linearization) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );

   SCIP_CALL( estimateExp(scip, conshdlr, expr, sol, FALSE, SCIPinfinity(scip), &coef, &constant, &local, &success, NULL) );

   cr_assert(success);
   cr_assert_float_eq(constant, -exp(2), SCIPepsilon(scip));
   cr_assert_float_eq(coef, exp(2), SCIPepsilon(scip));
   cr_assert(!local);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
