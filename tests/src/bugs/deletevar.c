/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   deletevar.c
 * @brief  Test for deleting a variable involved in a linear constraint at the SCIP_STAGE_PROBLEM
 * @author Mathieu Besançon
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/dialog_default.h"

#include "include/scip_test.h"

#undef SCIP_CALL
#define SCIP_CALL(x)   do                                                                 \
                       {                                                                  \
                          SCIP_RETCODE _restat_;                                          \
                          if( (_restat_ = (x)) != SCIP_OKAY )                             \
                          {                                                               \
                             cr_assert(FALSE, "Error <%d> in function call\n", _restat_); \
                          }                                                               \
                       }                                                                  \
                       while( FALSE )

SCIP* scip;

/** run unittest */
Test(deletevar, basictest)
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_Bool feasible;
   SCIP_RETCODE retcode;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPcreateProbBasic(scip, "unittest-deletevar") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "", 0.0, 1.0, 0.0,
                              SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, xvar) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "", 0.0, 1.0, 0.0,
                              SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "", 0, NULL, NULL, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, cons, xvar, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, cons, yvar, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPdelVar(scip, xvar, &feasible) );

   retcode = SCIPsolve(scip);
   cr_assert_eq(retcode, SCIP_INVALIDDATA);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
}
