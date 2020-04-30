/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   separation_log.c
 * @brief  tests separation of log()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_log.c"
#include "separation.h"

Test(separation, logarithmic, .init = setup, .fini = teardown,
   .description = "test separation for a logarithmic expression"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Real coef;
   SCIP_Real constant;
   SCIP_Bool islocal;
   SCIP_Bool branchcand;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreateConsExprExprLog(scip, conshdlr, &expr, zexpr) );

   /* compute an overestimation (linearization) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );

   branchcand = TRUE;
   SCIP_CALL( estimateLog(scip, conshdlr, expr, sol, TRUE, SCIPinfinity(scip), &coef, &constant, &islocal, &success, &branchcand) );

   cr_assert(success);
   cr_assert_float_eq(constant, -1.0 + log(2.0), SCIPepsilon(scip));
   cr_assert_float_eq(coef, 0.5, SCIPepsilon(scip));
   cr_assert(!islocal);
   cr_assert(!branchcand);

   /* compute an underestimation (secant) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );

   branchcand = TRUE;
   SCIP_CALL( estimateLog(scip, conshdlr, expr, sol, FALSE, -SCIPinfinity(scip), &coef, &constant, &islocal, &success, &branchcand) );

   cr_assert(success);
   cr_assert_float_eq(constant, -log(3.0)/2.0, SCIPepsilon(scip));
   cr_assert_float_eq(coef, log(3.0)/2.0, SCIPepsilon(scip));
   cr_assert(islocal);
   cr_assert(branchcand);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
