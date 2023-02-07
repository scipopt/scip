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
 * @brief  tests separation of abs()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_abs.c"
#include "../estimation.h"

static SCIP_Real coefs[SCIP_EXPR_MAXINITESTIMATES];
static SCIP_Real* coefsp[SCIP_EXPR_MAXINITESTIMATES];
static SCIP_Real constants[SCIP_EXPR_MAXINITESTIMATES];
static int nreturned;

Test(separation, absolute, .init = setup, .fini = teardown, .description = "test separation for an absolute expression")
{
   int i;
   for( i = 0; i < SCIP_EXPR_MAXINITESTIMATES; ++i )
      coefsp[i] = &coefs[i];

   SCIP_CALL( SCIPevalExprActivity(scip, xexpr) );

   /* compute overestimating secant */
   SCIP_CALL( computeCutsAbs(scip, SCIPexprGetActivity(xexpr), TRUE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 1);

   /* check secant */
   EXPECTFEQ( coefs[0], 2.0 / 3.0 );
   EXPECTFEQ( constants[0], 5.0 / 3.0 );

   /* compute underestimators */
   SCIP_CALL( computeCutsAbs(scip, SCIPexprGetActivity(xexpr), FALSE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 2);

   /* check left tangent */
   EXPECTFEQ( coefs[0], -1.0 );
   EXPECTFEQ( constants[0], 0.0 );

   /* check right tangent */
   EXPECTFEQ( coefs[1], 1.0 );
   EXPECTFEQ( constants[1], 0.0 );

}
