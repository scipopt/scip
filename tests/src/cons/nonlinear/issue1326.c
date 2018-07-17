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

/**@file   issue1326.c
 * @brief  unit tests to check that bug1326 doesn't appear
 * @author Felipe Serrano
 */

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include <stdio.h>

Test(issue1326, read_solve_prob)
{
   SCIP* scip;
   char filename[SCIP_MAXSTRLEN];

   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* trick to get cip file on same directory without knowing path */
   (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%sip", __FILE__);
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );

   SCIP_CALL( SCIPsolve(scip) );

   SCIP_CALL( SCIPfree(&scip) );
}
