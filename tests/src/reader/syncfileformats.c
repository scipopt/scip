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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   syncfileformats.c
 * @brief  Unittest for checking whether a file written in different formats returns the same value
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "include/scip_test.h"

static SCIP* scip;

#define EPS 1e-4

/** create SCIP instance and problem */
static
void setup(void)
{
   /* create SCIP instance */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
}

/** free SCIP instance */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
}

/* TEST SUITE */
TestSuite(syncfileformats, .init = setup, .fini = teardown);

/** check writing MIPs */
Test(syncfileformats, mip, .description = "check writing MIPs")
{
   const char* filenamein = "../check/instances/MIP/p0033.osil";
   const char* filesout[] = {"p0033.cip", "p0033.lp", "p0033.mps", "p0033.opb", "p0033.pip", "p0033.fzn"};
   const SCIP_Real optval = 3089;
   int nfilesout;
   int i;

   nfilesout = (int) sizeof(filesout) / sizeof(filesout[0]);
   for (i = 0; i < nfilesout; ++i)
   {
      SCIPdebugMsg(scip, "Considering <%s> ...\n", filesout[i]);

      /* read problem */
      SCIP_CALL( SCIPreadProb(scip, filenamein, NULL) );

      /* write in different format, use generic variable names */
      SCIP_CALL( SCIPwriteOrigProblem(scip, filesout[i], NULL, TRUE) );

      /* free problem */
      SCIP_CALL( SCIPfreeProb(scip) );

      /* read problem */
      SCIP_CALL( SCIPreadProb(scip, filesout[i], NULL) );

      /* solve */
      SCIP_CALL( SCIPsolve(scip) );

      cr_assert( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL );
      cr_assert_float_eq(SCIPgetDualbound(scip), optval, EPS, "Problem <%s> produces different optimal value: %g != %g\n", filesout[i], SCIPgetDualbound(scip), optval);

      /* free problem */
      SCIP_CALL( SCIPfreeProb(scip) );

      (void)remove(filesout[i]);
   }
}


/** check writing MINLPs */
Test(syncfileformats, minlp, .description = "check writing quadratic MINLPs")
{
   const char* filenamein = "../check/instances/MINLP/ex1266.mps";
   const char* filesout[] = {"ex1266.cip", "ex1266.pip", "ex1266.lp"};   /* osil is not available for writing */
   const SCIP_Real optval = 16.3;
   int nfilesout;
   int i;

   nfilesout = (int) sizeof(filesout) / sizeof(filesout[0]);
   for (i = 0; i < nfilesout; ++i)
   {
      SCIPdebugMsg(scip, "Considering <%s> ...\n", filesout[i]);

      /* read problem */
      SCIP_CALL( SCIPreadProb(scip, filenamein, NULL) );

      /* write in different format */
      SCIP_CALL( SCIPwriteOrigProblem(scip, filesout[i], NULL, FALSE) );

      /* free problem */
      SCIP_CALL( SCIPfreeProb(scip) );

      /* read problem */
      SCIP_CALL( SCIPreadProb(scip, filesout[i], NULL) );

      /* solve */
      SCIP_CALL( SCIPsolve(scip) );

      cr_assert( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL );
      cr_assert_float_eq(SCIPgetDualbound(scip), optval, EPS, "Problem <%s> produces different optimal value: %g != %g\n", filesout[i], SCIPgetDualbound(scip), optval);

      /* free problem */
      SCIP_CALL( SCIPfreeProb(scip) );

      (void)remove(filesout[i]);
   }
}
