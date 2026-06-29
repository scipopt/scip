/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   concurrent.c
 * @brief  unit test for concurrent solving combined with reoptimization
 * @author Mohammed Ghannam
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "tpi/tpi.h"
#include "include/scip_test.h"

/** concurrent solving does not support reoptimization: SCIPsolveConcurrent() must return SCIP_NOTIMPLEMENTED instead of
 *  running into a segmentation fault (see PySCIPOpt#884).
 */
Test(concurrent, reoptimization)
{
   SCIP* scip;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_CONS* cons;
   SCIP_VAR* consvars[2];
   SCIP_Real consvals[2];

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );

   SCIP_CALL( SCIPcreateProbBasic(scip, "concurrent-reopt") );
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 10.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 10.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* y - 2x >= 0 */
   consvars[0] = y;
   consvars[1] = x;
   consvals[0] = 1.0;
   consvals[1] = -2.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "c", 2, consvars, consvals, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPenableReoptimization(scip, TRUE) );

   /* without the task processing interface SCIPsolveConcurrent() bails out earlier with SCIP_PLUGINNOTFOUND, so only
    * check the reoptimization guard when concurrent solving is actually available
    */
   if( SCIPtpiIsAvailable() )
   {
      cr_assert_eq(SCIPsolveConcurrent(scip), SCIP_NOTIMPLEMENTED,
         "SCIPsolveConcurrent() must reject reoptimization with SCIP_NOTIMPLEMENTED.");
   }

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );
}
