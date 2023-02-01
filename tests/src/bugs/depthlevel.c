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

/**@file   depthlevel.c
 * @brief  Test for freeing SCIP if depth level has been reached
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include <signal.h>
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

/** create unbounded infeasible problem:
 * min   0
 * s.t.  1/4 <= x - y <= 3/4
 *       x,y \in Z
 */
static
void setup(void)
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPcreateProbBasic(scip, "unittest-depthlevel") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );

   /* create inequalities */
   vars[0] = xvar;
   vars[1] = yvar;

   vals[0] = 1.0;
   vals[1] = -1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "lower", 2, vars, vals, 0.25, 0.75) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );

   /*SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) ); */
}


/** run unittest */
Test(depthlevel, hit_depth_limit, .description = "show problem when hitting depth level", .init = setup, .signal = SIGABRT)
{
   SCIP_RETCODE retcode;

   /* this test can only work in debug mode, so we make it pass in opt mode */
#ifdef NDEBUG
   abort(); /* return SIGABORT */
#endif

   /* turn off presolving and rapid learning to avoid detecting infeasibility */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/rapidlearning/freq", -1) );

   /* use DFS */
   SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/dfs/stdpriority", 10000000) );

   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", 10000) );

   /* solve */
   retcode = SCIPsolve(scip);
   cr_expect(retcode == SCIP_OKAY || retcode == SCIP_MAXDEPTHLEVEL);

   /* print statistics */
   /* SCIP_CALL( SCIPprintStatistics(scip, NULL) ); */

   /* free transformed problem */
   SCIP_CALL( SCIPfreeTransform(scip) );

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}
