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

   /* include test constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, "name", "test initlp",
            1, 1, -1, FALSE, consEnfolpTest, consEnfopsTest, consCheckTest,
            consLockTest, NULL) );

   /* set initlp callback */
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpTest) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* solve problem and expect it to be infeasible */
   SCIP_CALL( SCIPsolve(scip) );
   cr_expect(SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE);

   /* free memory and assert no memory leak */
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!");
}
