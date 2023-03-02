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
 * @brief  tests estimation of log()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_log.c"
#include "../estimation.h"

Test(separation, logarithmic, .init = setup, .fini = teardown,
   .description = "test estimation for a logarithmic expression"
   )
{
   SCIP_EXPR* expr;
   SCIP_Real coef;
   SCIP_Real constant;
   SCIP_Bool islocal;
   SCIP_Bool branchcand;
   SCIP_Bool success;
   SCIP_Real refpoint;
   SCIP_INTERVAL localbounds;
   SCIP_INTERVAL globalbounds;

   SCIP_CALL( SCIPcreateExprLog(scip, &expr, zexpr, NULL, NULL) );

   SCIPintervalSetBounds(&localbounds, SCIPvarGetLbLocal(z), SCIPvarGetUbLocal(z));
   SCIPintervalSetBounds(&globalbounds, SCIPvarGetLbGlobal(z), SCIPvarGetUbGlobal(z));

   /* compute an overestimation (linearization) */
   refpoint = 2.0;

   branchcand = TRUE;
   SCIP_CALL( estimateLog(scip, expr, &localbounds, &globalbounds, &refpoint, TRUE, SCIPinfinity(scip), &coef,
         &constant, &islocal, &success, &branchcand) );

   cr_assert(success);
   cr_assert_float_eq(constant, -1.0 + log(2.0), SCIPepsilon(scip));
   cr_assert_float_eq(coef, 0.5, SCIPepsilon(scip));
   cr_assert(!islocal);
   cr_assert(!branchcand);

   /* compute an underestimation (secant) */
   refpoint = 2.0;

   branchcand = TRUE;
   SCIP_CALL( estimateLog(scip, expr, &localbounds, &globalbounds, &refpoint, FALSE, -SCIPinfinity(scip), &coef,
         &constant, &islocal, &success, &branchcand) );

   cr_assert(success);
   cr_assert_float_eq(constant, -log(3.0)/2.0, SCIPepsilon(scip));
   cr_assert_float_eq(coef, log(3.0)/2.0, SCIPepsilon(scip));
   cr_assert(islocal);
   cr_assert(branchcand);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
