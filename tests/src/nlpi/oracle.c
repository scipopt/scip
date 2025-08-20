/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   oracle.c
 * @brief  unit test for nonlinear oracle methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include "scip/scipdefplugins.h"
#include "scip/nlpioracle.h"
#include "scip/nlpioracle.c"
#include "include/scip_test.h"

#define INF 1e+20

#define expecti(x, y) cr_expect_eq(x, y, "%s expected to be %d, got %d", #x, y, x)

static SCIP* scip = NULL;
static SCIP_NLPIORACLE* oracle = NULL;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* need a problem to have stat created, which is used by expr iterators */
   SCIP_CALL( SCIPcreateProbBasic(scip, "dummy") );

   /* create the oracle */
   SCIP_CALL( SCIPnlpiOracleCreate(scip, &oracle) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPnlpiOracleFree(scip, &oracle) );
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

Test(oracle, jacsparsity, .init = setup, .fini = teardown,
   .description = "checks Jacobian sparsity structure"
   )
{
   const char* varnames[4] = {"x0", "x1", "x2", "x3"};
   SCIP_EXPR* varexprs[4];
   SCIP_Real lbs[4] = {-0.5, -2, 0, -2};
   SCIP_Real ubs[4] = {INF, INF, 2, 2};

   /* constraints */
   const char* consnames[3] = {"c1", "c2", "c3"};
   SCIP_Real lhss[3] = {1, -INF, 2.5};
   SCIP_Real rhss[3] = {1, -1, 5};
   SCIP_Real* linvals;
   int* lininds;
   int nlinds;

   SCIP_EXPR* x0sqr;
   SCIP_EXPR* x1sqr;
   SCIP_EXPR* x2sqr;
   SCIP_EXPR* x0x2;
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, 2) );
   SCIP_CALL( SCIPnlpiOracleAddVars(scip, oracle, 4, lbs, ubs, varnames) );

   /* create expressions */
   SCIP_CALL( SCIPcreateExprVaridx(scip, &varexprs[0], 0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVaridx(scip, &varexprs[1], 1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVaridx(scip, &varexprs[2], 2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVaridx(scip, &varexprs[3], 3, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &x0sqr, varexprs[0], 2.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &x1sqr, varexprs[1], 2.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &x2sqr, varexprs[2], 2.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &x0x2, 1, &varexprs[0], 1.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, x0x2, varexprs[2]) );

   /* constraints */
   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 1, &x0sqr, NULL, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr, x2sqr, 1.0) );
   SCIP_CALL( SCIPappendExprSumExpr(scip, expr, x0x2, 1.0) );
   SCIP_CALL( SCIPnlpiOracleAddConstraints(scip, oracle, 1, &lhss[0], &rhss[0], NULL, NULL, NULL, &expr, &consnames[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   nlinds = 2;
   lininds[0] = 2;
   linvals[0] = 1;
   lininds[1] = 3;
   linvals[1] = -1;
   SCIP_CALL( SCIPnlpiOracleAddConstraints(scip, oracle, 1, &lhss[1], &rhss[1], &nlinds, &lininds, &linvals, NULL, &consnames[1]) );

   nlinds = 2;
   lininds[0] = 1;
   linvals[0] = 1;
   lininds[1] = 3;
   linvals[1] = 1;
   SCIP_CALL( SCIPnlpiOracleAddConstraints(scip, oracle, 1, &lhss[2], &rhss[2], &nlinds, &lininds, &linvals, NULL, &consnames[2]) );

   oracle->jaccoloffsets = NULL;
   oracle->jacrowoffsets = NULL;

   const int* jacrowoffsets = NULL;
   const int* jaccols = NULL;
   const int* jaccoloffsets = NULL;
   const int* jacrows = NULL;
   const SCIP_Bool* jaccolnlflags = NULL;
   const SCIP_Bool* jacrownlflags = NULL;
   int njacnlnz;

   /* rowwise sparsity */
   SCIP_CALL( SCIPnlpiOracleGetJacobianRowSparsity(scip, oracle, &jacrowoffsets, &jaccols, &jaccolnlflags, &njacnlnz) );

   expecti(njacnlnz, 2);

   expecti(jacrowoffsets[0], 0);
   expecti(jacrowoffsets[1], 2);
   expecti(jacrowoffsets[2], 4);
   expecti(jacrowoffsets[3], 6);

   expecti(jaccols[0], 0);
   expecti(jaccols[1], 2);
   expecti(jaccols[2], 2);
   expecti(jaccols[3], 3);
   expecti(jaccols[4], 1);
   expecti(jaccols[5], 3);

   cr_assert(jaccolnlflags != NULL);

   expecti(jaccolnlflags[0], TRUE);
   expecti(jaccolnlflags[1], TRUE);
   expecti(jaccolnlflags[2], FALSE);
   expecti(jaccolnlflags[3], FALSE);
   expecti(jaccolnlflags[4], FALSE);
   expecti(jaccolnlflags[5], FALSE);

   invalidateJacobiSparsity(scip, oracle);

   /* columnwise sparsity */
   SCIP_CALL( SCIPnlpiOracleGetJacobianColSparsity(scip, oracle, &jaccoloffsets, &jacrows, &jacrownlflags, &njacnlnz) );

   expecti(njacnlnz, 2);

   expecti(jaccoloffsets[0], 0);
   expecti(jaccoloffsets[1], 1);
   expecti(jaccoloffsets[2], 2);
   expecti(jaccoloffsets[3], 4);
   expecti(jaccoloffsets[4], 6);

   expecti(jacrows[0], 0);
   expecti(jacrows[1], 2);
   expecti(jacrows[2], 0);
   expecti(jacrows[3], 1);
   expecti(jacrows[4], 1);
   expecti(jacrows[5], 2);

   cr_assert(jacrownlflags != NULL);

   expecti(jacrownlflags[0], TRUE);
   expecti(jacrownlflags[1], FALSE);
   expecti(jacrownlflags[2], TRUE);
   expecti(jacrownlflags[3], FALSE);
   expecti(jacrownlflags[4], FALSE);
   expecti(jacrownlflags[5], FALSE);

   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &lininds);

   SCIP_CALL( SCIPreleaseExpr(scip, &x0x2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &x2sqr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &x1sqr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &x0sqr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[3]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[2]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[0]) );
}
