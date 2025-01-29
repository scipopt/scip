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

/**@file   copy.c
 * @brief  tests for default copy callbacks
 * @author Dominik Kamp
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"


/** GLOBAL CONSTANTS */

/* TODO: implement default copy */
const char* const missing[] =
{
   "benders",
   "concurrent",
   "display",
   "table",
   "compression",
   "solvingphases",
   "propagating/symmetry"
};
const int nmissing = sizeof(missing) / sizeof(missing[0]);

/** GLOBAL VARIABLES **/

static SCIP* sourcescip = NULL;
static SCIP* targetscip = NULL;

/* TEST SUITE */

/** create the necessary source SCIP with default plugins */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&sourcescip) );
   SCIP_CALL( SCIPcreate(&targetscip) );
}

/** free all allocated memory */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&targetscip) );
   SCIP_CALL( SCIPfree(&sourcescip) );
}

TestSuite(copy, .init = setup, .fini = teardown);

/* TESTS  */

Test(copy, MIP, .description="tests MIP copy")
{
   /* initialize validity flag */
   SCIP_PARAM** sourceparams;
   SCIP_PARAM** targetparams;
   SCIP_Bool valid = FALSE;
   const char* sourcename;
   const char* targetname;
   int nsourceparams;
   int ntargetparams;
   int i;
   int j;
   int k;

   /* setup source scip */
   SCIP_CALL( SCIPincludeDefaultPlugins(sourcescip) );
   /* TODO: remove symmetry restriction */
   SCIP_CALL( SCIPsetIntParam(sourcescip, "misc/usesymmetry", 0) );
   SCIP_CALL( SCIPreadProb(sourcescip, "../check/instances/MIP/misc03.mps", NULL) );
   sourceparams = SCIPgetParams(sourcescip);
   nsourceparams = SCIPgetNParams(sourcescip);

   /* copy target scip */
   SCIP_CALL( SCIPcopy(sourcescip, targetscip, NULL, NULL, "", TRUE, TRUE, TRUE, TRUE, &valid) );
   cr_assert(valid);
   targetparams = SCIPgetParams(targetscip);
   ntargetparams = SCIPgetNParams(targetscip);

   /* check copy dimension */
   /* TODO: complete the copy */
   cr_expect(ntargetparams <= nsourceparams);
   cr_expect(SCIPgetNVars(targetscip) == SCIPgetNVars(sourcescip));
   cr_expect(SCIPgetNConss(targetscip) == SCIPgetNConss(sourcescip));

   /* check matching parameters */
   i = 0;
   j = 0;
   k = 0;

   while( i < nsourceparams )
   {
      int l;

      sourcename = SCIPparamGetName(sourceparams[i]);

      for( l = 0; l < nmissing; ++l )
      {
         if( !strncmp(sourcename, missing[k], strlen(missing[k])) )
            break;

         ++k;

         if( k == nmissing )
            k = 0;
      }

      ++i;

      if( l < nmissing )
         continue;

      valid = FALSE;

      while( j < ntargetparams )
      {
         targetname = SCIPparamGetName(targetparams[j]);

         for( l = 0; l < nmissing; ++l )
         {
            if( !strncmp(targetname, missing[k], strlen(missing[k])) )
               break;

            ++k;

            if( k == nmissing )
               k = 0;
         }

         ++j;

         if( l < nmissing )
            continue;

         if( valid || strcmp(targetname, sourcename) )
         {
            --j;
            break;
         }

         valid = TRUE;
      }

      /* check target parameters */
      cr_expect(valid, "Parameter %s not in target.", sourcename);
   }

   /* check source parameters */
   cr_expect(j == ntargetparams, "Parameter %s not in source.", targetname);

   /* solve both scips */
   SCIP_CALL( SCIPsolve(sourcescip) );
   SCIP_CALL( SCIPsolve(targetscip) );

   /* compare optimal bounds */
   cr_expect(SCIPgetDualbound(targetscip) == SCIPgetPrimalbound(targetscip));
   cr_expect(SCIPgetDualbound(targetscip) == SCIPgetDualbound(sourcescip));
   cr_expect(SCIPgetPrimalbound(targetscip) == SCIPgetPrimalbound(sourcescip));

   /* compare solving performance */
   cr_expect(SCIPgetNLPIterations(targetscip) == SCIPgetNLPIterations(sourcescip));
   cr_expect(SCIPgetNTotalNodes(targetscip) == SCIPgetNTotalNodes(sourcescip));
   cr_expect(SCIPgetNSolsFound(targetscip) == SCIPgetNSolsFound(sourcescip));
}
