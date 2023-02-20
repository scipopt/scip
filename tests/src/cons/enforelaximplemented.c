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
