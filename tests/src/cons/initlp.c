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

/**@file   initlp.c
 * @brief  unit test for checking behaviour of the initlp callback
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include "scip/scipdefplugins.h"

/* dummy check that always succeeds unless scip has constructed its LP */
static
SCIP_DECL_CONSCHECK(consCheckTest)
{
   printf("checking test cons\n");
   *result = SCIP_INFEASIBLE;

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPisLPConstructed(scip) )
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}
/* dummy initlp that returns infeasible */
static
SCIP_DECL_CONSINITLP(consInitlpTest)
{
   printf("initing lp test cons\n");
   *infeasible = TRUE;
   return SCIP_OKAY;
}
/* dummy enfolp that always succeeds */
static
SCIP_DECL_CONSENFOLP(consEnfolpTest)
{
   printf("enforcing lp test cons\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}
/* dummy enfops that always succeeds */
static
SCIP_DECL_CONSENFOPS(consEnfopsTest)
{
   printf("enforcing ps test cons\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}
/* dummy locks */
static
SCIP_DECL_CONSLOCK(consLockTest)
{
   printf("locking test cons\n");
   return SCIP_OKAY;
}

/** TEST  **/
#include "include/scip_test.h"

Test(cons, initlp)
{
   SCIP* scip;
   SCIP_CONSHDLR* conshdlr;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* include test constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, "name", "test initlp",
            1, 1, -1, FALSE, consEnfolpTest, consEnfopsTest, consCheckTest,
            consLockTest, NULL) );

   /* SCIP can go to SOLVED after presolving if there are no vars, cons, nor pricer; we add a pricer to avoid this;
    * this pricer is defined in scip_test
    */
   SCIP_CALL( SCIPincludePricerBasic(scip, NULL, "pricerTest", "pricer to avoid SCIP skipping SOLVING", 0, FALSE,
            pricerRedcostTest, NULL, NULL) );
   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "pricerTest")) );

   /* set initlp callback */
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpTest) );

   /* solve problem and expect it to be infeasible */
   SCIP_CALL( SCIPsolve(scip) );
   cr_expect(SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE);

   /* free memory and assert no memory leak */
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!");
}
