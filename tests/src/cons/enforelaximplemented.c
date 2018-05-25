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

/**@file   enforelaximplemented.c
 * @brief  unit test to check if all core constraint handlers implement enforelax callback
 * @author Tristan Gally
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons.h"
#include "include/scip_test.h"

/** All conshdlrs should implement enforelax callback, otherwise relaxation solution cannot be enforced. The callback
 *  is non-mandatory only to not break the code of users who are not interested in relaxations anyways.
 */
Test(cons, enforelaximplemented, .description = "tests if all constraint handlers implement enforelax callback")
{
   SCIP* scip;
   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;
   int h;

   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* iterate over all conshdlrs and check if they implement the enforelax callback */
   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);
   for( h = 0; h < nconshdlrs; h++ )
   {
      cr_expect_not_null(conshdlrs[h]->consenforelax, "conshdlr %s should implement enforelax callback",
            SCIPconshdlrGetName(conshdlrs[h]));
   }

   SCIP_CALL( SCIPfree(&scip) );
}
