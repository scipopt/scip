/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   estimation.c
 * @brief  tests estimators of exp()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_exp.c"
#include "../estimation.h"

Test(estimation, exponential, .init = setup, .fini = teardown,
   .description = "test estimation for an exponential expression"
   )
{
   SCIP_EXPR* expr;
   SCIP_Real coef;
   SCIP_Real constant;
   SCIP_Bool local;
   SCIP_Bool success;
   SCIP_Bool branchcand = TRUE;
   SCIP_Real xval;
   SCIP_INTERVAL bnd;

   SCIP_CALL( SCIPcreateExprExp(scip, &expr, xexpr, NULL, NULL) );

   SCIP_CALL( SCIPevalExprActivity(scip, expr) );
   bnd = SCIPexprGetActivity(xexpr);

   /* compute an overestimation (secant) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   xval = 0.0;

   branchcand = TRUE;
   SCIP_CALL( estimateExp(scip, expr, &bnd, &bnd, &xval, TRUE, -SCIPinfinity(scip), &coef, &constant, &local, &success, &branchcand) );

   cr_assert(success);
   cr_assert_float_eq(constant, (exp(5) + 5 * exp(-1)) / 6.0, SCIPepsilon(scip));
   cr_assert_float_eq(coef, (exp(5) - exp(-1)) / 6.0, SCIPepsilon(scip));
   cr_assert(local);
   cr_assert(branchcand);

   /* compute an underestimation (linearization) */
   xval = 2.0;

   branchcand = TRUE;
   SCIP_CALL( estimateExp(scip, expr, &bnd, &bnd, &xval, FALSE, SCIPinfinity(scip), &coef, &constant, &local, &success, &branchcand) );

   cr_assert(success);
   cr_assert_float_eq(constant, -exp(2), SCIPepsilon(scip));
   cr_assert_float_eq(coef, exp(2), SCIPepsilon(scip));
   cr_assert(!local);
   cr_assert(!branchcand);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
