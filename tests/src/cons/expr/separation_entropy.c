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

/**@file   separation_entropy.c
 * @brief  tests separation of entropy()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_entropy.c"
#include "separation.h"

Test(separation, entropy, .init = setup, .fini = teardown,
   .description = "test separation for an exponential expression"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Real coef;
   SCIP_Real constant;
   SCIP_Bool success;
   SCIP_Bool local;

   SCIP_CALL( SCIPcreateConsExprExprEntropy(scip, conshdlr, &expr, zexpr) );

   /* compute an overestimation (linearization) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 0.0) );

   SCIP_CALL( estimateEntropy(scip, conshdlr, expr, sol, TRUE, -SCIPinfinity(scip), &coef, &constant, &local, &success, NULL) );

   cr_assert(success);
   cr_assert_float_eq(constant, 2.0, SCIPepsilon(scip));
   cr_assert_float_eq(coef, -log(2.0) - 1.0, SCIPepsilon(scip));
   cr_assert(!local);

   /* compute an underestimation (secant) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -10.0) );

   SCIP_CALL( estimateEntropy(scip,  conshdlr, expr, sol, FALSE, SCIPinfinity(scip), &coef, &constant, &local, &success, NULL) );

   cr_assert(success);
   cr_assert_float_eq(constant, 1.5 * log(3.0) - 1.5 * log(1.0), SCIPepsilon(scip));
   cr_assert_float_eq(coef, 0.5 * (-3.0 * log(3.0) + log(1.0)), SCIPepsilon(scip));
   cr_assert(local);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
